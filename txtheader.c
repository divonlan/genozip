// ------------------------------------------------------------------
//   txtheader.c
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "buffer.h"
#include "digest.h"
#include "flags.h"
#include "data_types.h"
#include "file.h"
#include "sections.h"
#include "txtfile.h"
#include "vblock.h"
#include "zfile.h"
#include "crypt.h"
#include "bgzf.h"
#include "writer.h"
#include "strings.h"
#include "endianness.h"
#include "contigs.h"

static bool is_first_txt = true; 

//----------
// ZIP stuff
//----------

// ZIP: reads txt header and writes its compressed form to the GENOZIP file
bool txtheader_zip_read_and_compress (uint64_t *txt_header_size)
{    
    Digest header_digest = DIGEST_NONE;
    digest_initialize(); 

    z_file->disk_at_beginning_of_this_txt_file = z_file->disk_so_far;

    if (DTPT(txt_header_required) == HDR_MUST || DTPT (txt_header_required) == HDR_OK)
        header_digest = txtfile_read_header (is_first_txt); // reads into evb->txt_data and evb->lines.len

    // for VCF, we need to check if the samples are the same before approving binding (other data types can bind without restriction)
    //          also: header is modified if --chain or compressing a Luft file
    // for SAM, we check that the contigs specified in the header are consistent with the reference given in --reference/--REFERENCE
    *txt_header_size = evb->txt_data.len;
    if (!(DT_FUNC_OPTIONAL (txt_file, inspect_txt_header, true)(evb, &evb->txt_data, (struct FlagsTxtHeader){}))) { 
        // this is the second+ file in a bind list, but its samples are incompatible
        buf_free (&evb->txt_data);
        return false;
    }

    // in the 2nd+ binding file, we require the dual-components status to be consistent with earlier files
    if (z_file) { 
        if (z_file->num_components) {
            ASSERT (!(txt_file->coords && !z_dual_coords), 
                    "All binding files must be either dual-coordinates or not. However %s is a dual-coordinates file while previous file(s) are not", txt_name);

            ASSERT (!(!txt_file->coords && z_dual_coords), 
                    "All binding files must be either dual-coordinates or not. However %s is not a dual-coordinates file while previous file(s) are", txt_name);
        }
        else  // first component 
            z_has_gencomp |= (txt_file->coords != DC_NONE);  // note: in case of --chain, z_flags.dual_coords is set in file_open_z
    }

    if (z_file && !flag.seg_only)       
        // we always write the txt_header section, even if we don't actually have a header, because the section
        // header contains the data about the file
        zfile_write_txt_header (&evb->txt_data, *txt_header_size, header_digest, is_first_txt); // we write all headers in bound mode too, to support genounzip

    // for stats: combined length of txt headers in this bound file, or only one file if not bound
    if (!flag.bind) z_file->txt_txtheader_so_far_bind=0;
    
    // VCF note: we don't account for DVCF rejects files as txt_len - the variant lines are already accounted for in the main file, and the added header lines are duplicates of the main header
    // SAM/BAM note: we don't account for PRIM/DEPN txt headers generated in gencomp_initialize_file
    if (!flag.gencomp_num) 
        z_file->txt_txtheader_so_far_bind += evb->txt_data.len; 

    z_file->num_txt_components_so_far++; // when compressing

    buf_free (&evb->txt_data);
    
    is_first_txt = false;

    return true; // everything's good
}

//----------
// PIZ stuff
//----------

// PIZ: reads the txt header from the genozip file and outputs it to the reconstructed txt file
Coords txtheader_piz_read_and_reconstruct (uint32_t component_i, Section sl)
{
    z_file->disk_at_beginning_of_this_txt_file = z_file->disk_so_far;

    VBlock *comp_vb = vb_get_vb ("piz", 0);

    zfile_read_section (z_file, comp_vb, 0, &comp_vb->z_data, "z_data", SEC_TXT_HEADER, sl);

    SectionHeaderTxtHeader *header = (SectionHeaderTxtHeader *)comp_vb->z_data.data;
    ASSERT0 (header, "Incorrectly skipped SEC_TXT_HEADER - check skip function");

    // 1. if flag.unbind (genounzip) - we open the output txt file of the component
    // 2. if flag.one_component (genocat) - output to stdout or --output
    // 3. when reading an auxiliary file or no_writer- we create txt_file here (but don't actually open the physical file)
    if (!txt_file) { 
        const char *filename = flag.unbind ? txtfile_piz_get_filename (header->txt_filename, flag.unbind, false) 
                                           : flag.out_filename;
        txt_file = file_open (filename, WRITE, TXT_FILE, flag.out_dt != DT_NONE ? flag.out_dt : z_file->data_type);
        if (flag.unbind) FREE (filename); // file_open copies the names
    }

    // initialize if needed - but only once per outputted txt file 
    //i.e. if we have rejects+normal, or concatenated, we will only init in the first)
    if (!txt_file->piz_header_init_has_run && DTPZ(piz_header_init))
        txt_file->piz_header_init_has_run = true;

    txt_file->max_lines_per_vb     = BGEN32 (header->max_lines_per_vb);
    txt_file->txt_flags            = header->h.flags.txt_header;

    if (txt_file->codec == CODEC_BGZF)
        memcpy (txt_file->bgzf_signature, header->codec_info, 3);
        
    // note: in case of txt files containing two components - --luft and --interleave - this will be nonsense, 
    // but we don't test digest anyway as they are both flag.data_modified 
    txt_file->digest = header->digest_single; 

    // case: we need to reconstruct (or not) the BGZF following the instructions from the z_file
    if (flag.bgzf == FLAG_BGZF_BY_ZFILE) {

        // load the source file isize if we have it and we are attempting to reconstruct an unmodifed file identical to the source
        bool loaded = false;
        if (!flag.data_modified    
        && (z_file->num_components == 1 || flag.unbind))  // not concatenating multiple files
            loaded = bgzf_load_isizes (sl); // also sets txt_file->bgzf_flags

        // case: user wants to see this section header, despite not needing BGZF data
        else if (flag.only_headers == SEC_BGZF+1 || flag.only_headers == -1) {
            bgzf_load_isizes (sl); 
            buf_free (&txt_file->bgzf_isizes);
        }

        // case: we need to reconstruct back to BGZF, but we don't have a SEC_BGZF to guide us - we'll creating our own BGZF blocks
        if (!loaded && z_file->z_flags.bgzf)
            txt_file->bgzf_flags = (struct FlagsBgzf){ // case: we're creating our own BGZF blocks
                .has_eof_block = true, // add an EOF block at the end
                .library       = BGZF_LIBDEFLATE, // default - libdeflate level 6
                .level         = BGZF_COMP_LEVEL_DEFAULT 
            };

        header = (SectionHeaderTxtHeader *)comp_vb->z_data.data; // re-assign after possible realloc of z_data in bgzf_load_isizes
    }

    // case: the user wants us to reconstruct (or not) the BGZF blocks in a particular way, this overrides the z_file instructions 
    else 
        txt_file->bgzf_flags = (struct FlagsBgzf){ // case: we're creating our own BGZF blocks
            .has_eof_block = true, // add an EOF block at the end
            .library       = BGZF_LIBDEFLATE, 
            .level         = flag.bgzf 
        };

    // sanity        
    ASSERT (txt_file->bgzf_flags.level >= 0 && txt_file->bgzf_flags.level <= 12, "txt_file->bgzf_flags.level=%u out of range [0,12]", 
            txt_file->bgzf_flags.level);

    ASSERT (txt_file->bgzf_flags.library >= 0 && txt_file->bgzf_flags.library < NUM_BGZF_LIBRARIES, "txt_file->bgzf_flags.library=%u out of range [0,%u]", 
            txt_file->bgzf_flags.level, NUM_BGZF_LIBRARIES-1);

    // now get the text of the txt header itself. note: we decompress even if --no-header, bc we need to inspect
    zfile_uncompress_section (comp_vb, header, &comp_vb->txt_data, "txt_data", 0, SEC_TXT_HEADER);

    // count header-lines (for --lines etc): before data-modifying inspect_txt_header
    if (writer_is_txtheader_in_plan (component_i)) {
        if (flag.header_one && Z_DT(DT_VCF))
            txt_file->num_lines += 1;
        else if (!DTPT (is_binary))
            txt_file->num_lines += str_count_char (STRb(comp_vb->txt_data), '\n'); // number of source-file lines
    }

    if (comp_vb->txt_data.len)
        DT_FUNC_OPTIONAL (z_file, inspect_txt_header, true)(comp_vb, &comp_vb->txt_data, header->h.flags.txt_header); // ignore return value

    // hand-over txt header if it is needed (it won't be if flag.no_header)
    if (writer_is_txtheader_in_plan (component_i)) {

        // if we're translating from one data type to another (SAM->BAM, BAM->FASTQ, ME23->VCF etc) translate the txt header 
        // note: in a header-less SAM, after translating to BAM, we will have a header
        DtTranslation trans = dt_get_translation(NULL);
        if (trans.txtheader_translator && !flag.no_header) trans.txtheader_translator (comp_vb, &comp_vb->txt_data); 

        if (comp_vb->txt_data.len) {

            bool test_digest = !digest_is_zero (header->digest_header) && // in v8 without --md5, we had no digest
                               !flag.data_modified; // no point calculating digest if we know already the file will be different

            if (test_digest) {
                if (flag.log_digest) digest_start_log (&txt_file->digest_ctx_bound); 
                digest_update (&txt_file->digest_ctx_bound, &comp_vb->txt_data, "txt_header:digest_ctx_bound");
            }

            if (txt_file->codec == CODEC_BGZF) {
                // inherit BGZF blocks from source file, if isizes was loaded (i.e. not flag.data_modified) - 
                // into comp_vb->bgzf_blocks
                bgzf_calculate_blocks_one_vb (comp_vb, comp_vb->txt_data.len); 

                // compress unless flag.maybe_vb_modified_by_writer (we compress in writer_flush_vb instead)
                if (!flag.maybe_vb_modified_by_writer)
                    bgzf_compress_vb (comp_vb); 
            }

            if (test_digest && z_file->genozip_version >= 9) {  // backward compatability with v8: we don't test against v8 MD5 for the header, as we had a bug in v8 in which we included a junk MD5 if they user didn't --md5 or --test. any file integrity problem will be discovered though on the whole-file MD5 so no harm in skipping this.
                Digest reconstructed_header_digest = digest_do (STRb(comp_vb->txt_data));
                
                TEMP_FLAG (quiet, flag.quiet && !flag.show_digest);

                ASSERTW (digest_is_equal (header->digest_header, DIGEST_NONE) || 
                         digest_recon_is_equal (reconstructed_header_digest, header->digest_header) ||
                         flag.data_modified,
                         "%s of reconstructed %s header (%s) differs from original file (%s)\n"
                         "Bad reconstructed header has been dumped to: %s\n", digest_name(),
                         dt_name (z_file->data_type), digest_display (reconstructed_header_digest).s, digest_display (header->digest_header).s,
                         txtfile_dump_vb (comp_vb, z_name));

                RESTORE_FLAG (quiet);
            }
        }

        // if entire header is excluded by --lines, don't write it
        if (txt_file->num_lines && flag.lines_first >= txt_file->num_lines) 
            comp_vb->txt_data.len = 0;

        writer_handover_txtheader (&comp_vb, component_i); // handover data to writer thread (even if the header is empty, as the writer thread is waiting for it)

        // accounting for data as in original source file - affects vb->vb_position_txt_file of next VB
        txt_file->txt_data_so_far_single_0 = z_file->genozip_version < 12 ? BGEN32 (header->h.data_uncompressed_len) : BGEN64 (header->txt_header_size); 
    }

    // case: component is not in plan - discard the VB
    else
        vb_release_vb (&comp_vb);    
    
    if (!flag.reading_chain && !flag.reading_reference)
        is_first_txt = false;

    return header->h.flags.txt_header.gencomp_num;
}
