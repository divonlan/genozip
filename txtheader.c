// ------------------------------------------------------------------
//   txtheader.c
//   Copyright (C) 2019-2022 Genozip Limited
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
#include "piz.h"
#include "gencomp.h"

static bool is_first_txt = true; 

//----------
// ZIP stuff
//----------

// ZIP: reads txt header and writes its compressed form to the GENOZIP file
int64_t txtheader_zip_read_and_compress (int64_t *txt_header_offset, CompIType comp_i) // out (-1 if SEC_TXT_HEADER not written)
{    
    Digest header_digest = DIGEST_NONE;
    evb->comp_i = comp_i; // used by def_is_header_done
    
    z_file->disk_at_beginning_of_this_txt_file = z_file->disk_so_far;

    TxtHeaderRequirement req = DTPT (txt_header_required); 
    if (req == HDR_MUST || req == HDR_OK || (comp_i==0 && (req == HDR_MUST_0 || req == HDR_OK_0))) {
        txtfile_read_header (is_first_txt); // reads into evb->txt_data and evb->lines.len

        // get header digest if needed
        if (!flag.data_modified) 
            header_digest = digest_txt_header (&evb->txt_data, DIGEST_NONE);
    }

    // for VCF, we need to check if the samples are the same before approving binding (other data types can bind without restriction)
    //          also: header is modified if --chain or compressing a Luft file
    // for SAM, we check that the contigs specified in the header are consistent with the reference given in --reference/--REFERENCE
    uint64_t txt_header_size = evb->txt_data.len;
    if (!(DT_FUNC_OPTIONAL (txt_file, inspect_txt_header, true)(evb, &evb->txt_data, (struct FlagsTxtHeader){}))) { 
        buf_free (evb->txt_data);
        return -1;
    }

    // note: we always write the txt_header for comp_i=0 even if we don't actually have a header, because the
    // section header contains the data about the file. Special case: we don't write headers of SAM DEPN
    if (z_file && !flag.seg_only && !(z_sam_gencomp && comp_i)) {
        *txt_header_offset = z_file->disk_so_far;     
        zfile_write_txt_header (&evb->txt_data, txt_header_size, header_digest, is_first_txt, comp_i); // we write all headers in bound mode too, to support genounzip
    }
    else 
        *txt_header_offset = -1; // no SEC_TXT_HEADER section written

    if (flag.show_lines)
        iprintf ("txtheader bytes=%"PRIu64"\n", txt_header_size);
    
    // DVCF note: we don't account for rejects files as txt_len - the variant lines are already accounted for in the main file, and the added header lines are duplicates of the main header
    // SAM/BAM note: we don't account for PRIM/DEPN txt headers generated in gencomp_initialize
    if (!comp_i) 
        z_file->header_size = txt_file->header_size; 

    z_file->num_txts_so_far++; // when compressing

    buf_free (evb->txt_data);
    
    is_first_txt = false;

    return txt_header_size; // all good
}

//----------
// PIZ stuff
//----------

// PIZ main thread: reads the txt header from the genozip file and outputs it to the reconstructed txt file
void txtheader_piz_read_and_reconstruct (Section sec)
{
    z_file->disk_at_beginning_of_this_txt_file = z_file->disk_so_far;

    VBlockP txt_header_vb = vb_get_vb (POOL_MAIN, PIZ_TASK_NAME, 0, sec->comp_i);

    zfile_read_section (z_file, txt_header_vb, 0, &txt_header_vb->z_data, "z_data", SEC_TXT_HEADER, sec);

    SectionHeaderTxtHeader *header = (SectionHeaderTxtHeader *)txt_header_vb->z_data.data;
    ASSERT0 (header, "Incorrectly skipped SEC_TXT_HEADER - check skip function");

    if (sec->comp_i == COMP_MAIN)
        z_file->txt_header_single = *header;
        
    // 1. if flag.unbind (genounzip) - we open the output txt file of the component
    // 2. if flag.one_component (genocat) - output to stdout or --output
    // 3. when reading an auxiliary file or no_writer- we create txt_file here (but don't actually open the physical file)
    if (!txt_file) { 
        rom filename = flag.unbind ? txtfile_piz_get_filename (header->txt_filename, flag.unbind, false) 
                                   : flag.out_filename;
        txt_file = file_open (filename, WRITE, TXT_FILE, flag.out_dt != DT_NONE ? flag.out_dt : z_file->data_type);
        if (flag.unbind) FREE (filename); // file_open copies the names
        
        z_file->digest_ctx = DIGEST_CONTEXT_NONE; // reset digest
    }

    // initialize if needed - but only once per outputted txt file 
    //i.e. if we have rejects+normal, or concatenated, we will only init in the first)
    if (!txt_file->piz_header_init_has_run && DTPZ(piz_header_init))
        txt_file->piz_header_init_has_run = true;

    evb->comp_i                = !VER(14) ? header->h.flags.txt_header.v13_dvcf_comp_i/*v12,13*/ : sec->comp_i/*since v14*/;
    txt_file->max_lines_per_vb = BGEN32 (header->max_lines_per_vb);
    txt_file->txt_flags        = header->h.flags.txt_header;
    txt_file->num_vbs          = sections_count_sections_until (SEC_VB_HEADER, sec, SEC_TXT_HEADER);
    
    if (txt_file->codec == CODEC_BGZF)
        memcpy (txt_file->bgzf_signature, header->codec_info, 3);
        
    // SAM: only the main component has a TxtHeader section. 
    // DVCF: always zero. 
    // FASTQ: flag.data_modified when interleaving, or just one file if --R1/2. 
    // V8: zero if not compressed with --md5
    if (!digest_is_zero(header->digest) && !flag.data_modified) 
        z_file->digest = header->digest; 

    ASSINP (!flag.test || !digest_is_zero(header->digest) || z_file->z_flags.adler, 
            "--test cannot be used with %s, as it was compressed without a digest. See " WEBSITE_DIGEST, z_name);

    // case: we need to reconstruct (or not) the BGZF following the instructions from the z_file
    if (flag.bgzf == BGZF_BY_ZFILE) {

        // load the source file isize if we have it and we are attempting to reconstruct an unmodifed file identical to the source
        bool loaded = false;
        if (!flag.data_modified    
        && (z_file->num_components == 1 || flag.unbind))  // not concatenating multiple files
            loaded = bgzf_load_isizes (sec); // also sets txt_file->bgzf_flags

        // case: user wants to see this section header, despite not needing BGZF data
        else if (flag.only_headers == SEC_BGZF+1 || flag.only_headers == -1) {
            bgzf_load_isizes (sec); 
            buf_free (txt_file->bgzf_isizes);
        }

        // case: we need to reconstruct back to BGZF, but we don't have a SEC_BGZF to guide us - we'll creating our own BGZF blocks
        if (!loaded && (txt_file->codec == CODEC_BGZF))
            txt_file->bgzf_flags = (struct FlagsBgzf){ // case: we're creating our own BGZF blocks
                .has_eof_block = true, // add an EOF block at the end
                .library       = BGZF_LIBDEFLATE, // default - libdeflate level 6
                .level         = BGZF_COMP_LEVEL_DEFAULT 
            };

        header = (SectionHeaderTxtHeader *)txt_header_vb->z_data.data; // re-assign after possible realloc of z_data in bgzf_load_isizes
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
    zfile_uncompress_section (txt_header_vb, header, &txt_header_vb->txt_data, "txt_data", 0, SEC_TXT_HEADER);

    // count header-lines (for --lines etc): before data-modifying inspect_txt_header
    if (writer_does_txtheader_need_write (sec)) {
        if (flag.header_one && Z_DT(DT_VCF))
            txt_file->num_lines += 1;
        else if (!DTPT (is_binary))
            txt_file->num_lines += str_count_char (STRb(txt_header_vb->txt_data), '\n'); // number of source-file lines
    }

    if (txt_header_vb->txt_data.len)
        DT_FUNC_OPTIONAL (z_file, inspect_txt_header, true)(txt_header_vb, &txt_header_vb->txt_data, header->h.flags.txt_header); // ignore return value

    // hand-over txt header if it is needed (it won't be if flag.no_header)
    if (writer_does_txtheader_need_write (sec)) {

        // if we're translating from one data type to another (SAM->BAM, BAM->FASTQ, ME23->VCF etc) translate the txt header 
        // note: in a header-less SAM, after translating to BAM, we will have a header
        DtTranslation trans = dt_get_translation(NULL);
        if (trans.txtheader_translator && !flag.no_header) trans.txtheader_translator (txt_header_vb, &txt_header_vb->txt_data); 

        // count textual lines in header, used for line= reporting in ASSPIZ
        if (txt_header_vb->txt_data.len && !DTPT (is_binary)) {
            uint32_t num_textual_lines = str_count_char (STRb(txt_header_vb->txt_data), '\n');
            writer_set_num_txtheader_lines (sec->comp_i, num_textual_lines);
        }

        if (piz_need_digest) 
            digest_txt_header (&txt_header_vb->txt_data, header->digest_header);

        writer_handover_txtheader (&txt_header_vb); // handover data to writer thread (even if the header is empty, as the writer thread is waiting for it)

        // accounting for data as in original source file - affects vb->vb_position_txt_file of next VB
        txt_file->txt_data_so_far_single_0 = !VER(12) ? BGEN32 (header->h.data_uncompressed_len) : BGEN64 (header->txt_header_size); 
    }

    // case: component is not in plan - discard the VB
    else
        vb_release_vb (&txt_header_vb, PIZ_TASK_NAME);    
    
    if (!flag.reading_chain && !flag.reading_reference)
        is_first_txt = false;
}
