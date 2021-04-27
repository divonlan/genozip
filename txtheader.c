// ------------------------------------------------------------------
//   txtheader.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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

// contigs loaded from the txt header (eg SAM or BAM, VCF header)
static bool has_contigs     = false;        // for BAM, this will be true even for contig-less header, indicting the file has no contigs (as opposed to contigs just not defined in the header)
static Buffer contigs       = EMPTY_BUFFER; // array of RefContig
static Buffer ocontigs      = EMPTY_BUFFER;
static Buffer contigs_dict  = EMPTY_BUFFER;
static Buffer ocontigs_dict = EMPTY_BUFFER;

static bool is_first_txt = true; 
static uint32_t total_bound_txt_headers_len = 0; 

//----------
// ZIP stuff
//----------

uint32_t txtheader_get_bound_headers_len(void) { return total_bound_txt_headers_len; }

ConstBufferP txtheader_get_contigs (void)
{
    return has_contigs ? &contigs : NULL; 
}

const char *txtheader_get_contig_name (uint32_t index, uint32_t *snip_len)
{
    RefContig *rc = ENT (RefContig, contigs, index);
    *snip_len = rc->snip_len;
    return ENT (char, contigs_dict, rc->char_index);
}

void txtheader_finalize (void)
{
    if (!flag.bind) { // in ZIP with bind we keep the header contigs
        buf_free (&contigs);
        buf_free (&contigs_dict);
        buf_free (&ocontigs);
        buf_free (&ocontigs_dict);
        has_contigs = false;
    }
}

void txtheader_alloc_contigs (uint32_t more_contigs, uint32_t more_dict_len, bool liftover)
{
    has_contigs = true;
    
    buf_alloc (evb, liftover ? &ocontigs : &contigs, more_contigs, 100, RefContig, 2, "txtheader:contigs"); 
    buf_alloc (evb, liftover ? &ocontigs_dict : &contigs_dict, more_dict_len, 2000, char, 2, "txtheader:contigs_dict"); 
}

void txtheader_add_contig (const char *chrom_name, unsigned chrom_name_len, PosType last_pos, void *liftover_)
{
    bool liftover = !!liftover_;
    WordIndex ref_chrom=WORD_INDEX_NONE;

    if (!liftover && (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE))
        ref_chrom = ref_contigs_ref_chrom_from_header_chrom (chrom_name, chrom_name_len, &last_pos, contigs.len);

    Buffer *this_contigs = liftover ? &ocontigs      : &contigs;
    Buffer *this_dict    = liftover ? &ocontigs_dict : &contigs_dict;

    // add to contigs
    NEXTENT (RefContig, *this_contigs) = (RefContig){ 
        .max_pos     = last_pos, 
        .char_index  = this_dict->len, 
        .snip_len    = chrom_name_len,
        .chrom_index = ref_chrom
    };

    if (flag.show_txt_contigs) 
        iprintf ("%sindex=%u \"%.*s\" LN=%"PRId64" ref_chrom_index=%u snip_len=%u\n", 
                 liftover ? "liftover: " : "", (unsigned)this_contigs->len-1, chrom_name_len, chrom_name, last_pos, ref_chrom, chrom_name_len);

    // add to contigs_dict
    buf_add (this_dict, chrom_name, chrom_name_len);
    NEXTENT (char, *this_dict) = 0; // nul-termiante
}

#define next_contig param // we use contigs.param as "next_contig"

void txtheader_verify_contig_init (void) { contigs.next_contig = ocontigs.next_contig = 0; }

void txtheader_verify_contig (const char *chrom_name, unsigned chrom_name_len, PosType last_pos, void *liftover_)
{
    bool liftover = !!liftover_;

    Buffer *this_contigs = liftover ? &ocontigs      : &contigs;
    Buffer *this_dict    = liftover ? &ocontigs_dict : &contigs_dict;

    ASSINP (this_contigs->next_contig < this_contigs->len, 
            "Error: txt header: contigs mismatch between files: first file has %u %scontigs, but %s has more",
             (unsigned)this_contigs->len, liftover ? "liftover " : "", txt_name);

    RefContig *rc = ENT (RefContig, *this_contigs, this_contigs->next_contig++);
    const char *rc_chrom_name = ENT (char, *this_dict, rc->char_index);
    
    ASSINP (chrom_name_len == rc->snip_len && !memcmp (chrom_name, rc_chrom_name, chrom_name_len),
            "Error: SAM header contig=%u: contig name mismatch between files: in first file: \"%s\", in %s: \"%.*s\"",
            (unsigned)this_contigs->next_contig, rc_chrom_name, txt_name, chrom_name_len, chrom_name);
            
    ASSINP (last_pos == rc->max_pos, "Error: SAM header in \"%s\": contig length mismatch between files: in first file: LN:%"PRId64", in %s: LN:%"PRId64,
            rc_chrom_name, rc->max_pos, txt_name, last_pos);
}

// ZIP: reads txt header and writes its compressed form to the GENOZIP file
bool txtheader_zip_read_and_compress (uint32_t *txt_line_i)
{    
    Digest header_digest = DIGEST_NONE;
    digest_initialize(); 

    z_file->disk_at_beginning_of_this_txt_file = z_file->disk_so_far;

    if (DTPT(txt_header_required) == HDR_MUST || DTPT (txt_header_required) == HDR_OK)
        header_digest = txtfile_read_header (is_first_txt); // reads into evb->txt_data and evb->lines.len
    
    *txt_line_i += (uint32_t)evb->lines.len;

    // for VCF, we need to check if the samples are the same before approving binding (other data types can bind without restriction)
    //          also: header is modified if --chain or compressing a Luft file
    // for SAM, we check that the contigs specified in the header are consistent with the reference given in --reference/--REFERENCE
    uint64_t unmodified_txt_header_len = evb->txt_data.len;
    if (!(DT_FUNC_OPTIONAL (txt_file, inspect_txt_header, true)(evb, &evb->txt_data))) { 
        // this is the second+ file in a bind list, but its samples are incompatible
        buf_free (&evb->txt_data);
        return false;
    }

    // in the 2nd+ binding file, we require the dual-components status to be consistent with earlier files
    if (z_file) { 
        if (z_file->num_components) {
            ASSERT (!(txt_file->dual_coords && !z_file->z_flags.dual_coords), 
                    "All binding files must be either dual-coordinates or not. However %s is a dual-coordinates file while previous file(s) are not", txt_name);

            ASSERT (!(!txt_file->dual_coords && z_file->z_flags.dual_coords), 
                    "All binding files must be either dual-coordinates or not. However %s is not a dual-coordinates file while previous file(s) are", txt_name);
        }
        else  // first component 
              // note: in case of --chain, z_flags.dual_coords is set in file_open_z
              // note: z_flags.dual_coords is bool and txt_file->dual_coords is an enum
            z_file->z_flags.dual_coords = z_file->z_flags.dual_coords || !!txt_file->dual_coords; 
    }

    if (z_file && !flag.seg_only)       
        // we always write the txt_header section, even if we don't actually have a header, because the section
        // header contains the data about the file
        zfile_write_txt_header (&evb->txt_data, unmodified_txt_header_len, header_digest, is_first_txt); // we write all headers in bound mode too, to support genounzip

    // for stats: combined length of txt headers in this bound file, or only one file if not bound
    if (!flag.bind) total_bound_txt_headers_len=0;
    total_bound_txt_headers_len += evb->txt_data.len; 

    z_file->num_txt_components_so_far++; // when compressing

    buf_free (&evb->txt_data);
    
    is_first_txt = false;

    return true; // everything's good
}

// copy contigs from reference or header to CHROM (and if relevant - also oCHROM, RNEXT for SAM/BAM)
void txtheader_zip_prepopulate_contig_ctxs (void)
{
    ConstBufferP this_contigs=NULL, this_contigs_dict=NULL;

    if (has_contigs) { // note: always true for BAM. for SAM and VCF - true if there are contigs in the header.
        this_contigs = &contigs;
        this_contigs_dict = &contigs_dict;
    }
    else if (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE) 
        ref_contigs_get (&this_contigs_dict, &this_contigs);

    if (this_contigs && this_contigs->len) {
        ctx_build_zf_ctx_from_contigs (CHROM, this_contigs, this_contigs_dict); 

        if (z_file->data_type == DT_SAM || z_file->data_type == DT_BAM)
            ctx_build_zf_ctx_from_contigs (SAM_RNEXT, this_contigs, this_contigs_dict);
    }

    if (ocontigs.len)
        ctx_build_zf_ctx_from_contigs (DTFT (ochrom), &ocontigs, &ocontigs_dict);
}

//----------
// PIZ stuff
//----------

// PIZ: reads the txt header from the genozip file and outputs it to the reconstructed txt file
void txtheader_piz_read_and_reconstruct (uint32_t component_i, Section sl)
{
    bool show_headers_only = (flag.show_headers && exe_type == EXE_GENOCAT);

    z_file->disk_at_beginning_of_this_txt_file = z_file->disk_so_far;

    VBlock *comp_vb = vb_get_vb ("piz", 0);

    zfile_read_section (z_file, comp_vb, 0, &comp_vb->z_data, "header_section", SEC_TXT_HEADER, sl);

    // handle the GENOZIP header of the txt header section
    SectionHeaderTxtHeader *header = (SectionHeaderTxtHeader *)comp_vb->z_data.data;

    flag.processing_rejects = flag.luft && header->h.flags.txt_header.liftover_rejects;

    // 1. if flag.unbind (genounzip) - we open the output txt file of the component
    // 2. if flag.one_component (genocat) - output to stdout or --output
    // 3. when reading an auxiliary file or no_writer- we create txt_file here (but don't actually open the physical file)
    if (!txt_file) { 
        const char *filename = flag.unbind ? txtfile_piz_get_filename (header->txt_filename, flag.unbind, false) 
                                           : flag.out_filename;
        txt_file = file_open (filename, WRITE, TXT_FILE, flag.out_dt != DT_NONE ? flag.out_dt : z_file->data_type);
        FREE (filename); // file_open copies the names
    }

    // initialize if needed - but only once per outputted txt file 
    //i.e. if we have rejects+normal, or concatenated, we will only init in the first)
    if (!txt_file->piz_header_init_has_run && DTPZ(piz_header_init))
        txt_file->piz_header_init_has_run = true;

    txt_file->txt_data_size_single = BGEN64 (header->txt_data_size); 
    txt_file->max_lines_per_vb     = BGEN32 (header->max_lines_per_vb);
    
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
        else if (exe_type == EXE_GENOCAT && (flag.show_headers == SEC_BGZF+1 || flag.show_headers == -1)) {
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

    // now get the text of the txt header itself
    if (!show_headers_only)
        zfile_uncompress_section (comp_vb, header, &comp_vb->txt_data, "txt_data", 0, SEC_TXT_HEADER);

    if (comp_vb->txt_data.len)
        DT_FUNC_OPTIONAL (z_file, inspect_txt_header, true)(comp_vb, &comp_vb->txt_data); // ignore return value

    // if we're translating from one data type to another (SAM->BAM, BAM->FASTQ, ME23->VCF etc) translate the txt header 
    // note: in a header-less SAM, after translating to BAM, we will have a header
    DtTranslation trans = dt_get_translation();
    if (trans.txtheader_translator && !show_headers_only) trans.txtheader_translator (comp_vb, &comp_vb->txt_data); 

    // hand-over txt header if it is needed:
    if (writer_is_txtheader_in_plan (component_i)) {

        if (comp_vb->txt_data.len) {

            if (!DTPT (is_binary))
                txt_file->num_lines += str_count_char (comp_vb->txt_data.data, comp_vb->txt_data.len, '\n');
                
            bool test_digest = !digest_is_zero (header->digest_header) && // in v8 without --md5, we had no digest
                                !flag.data_modified; // no point calculating digest if we know already the file will be different

            if (test_digest) digest_update (&txt_file->digest_ctx_bound, &comp_vb->txt_data, "txt_header:digest_ctx_bound");

            if (txt_file->codec == CODEC_BGZF) {
                // inherit BGZF blocks from source file, if isizes was loaded (i.e. not flag.data_modified) - 
                // into comp_vb->bgzf_blocks
                bgzf_calculate_blocks_one_vb (comp_vb, comp_vb->txt_data.len); 

                // compress unless flag.maybe_vb_modified_by_writer (we compress in writer_flush_vb instead)
                if (!flag.maybe_vb_modified_by_writer)
                    bgzf_compress_vb (comp_vb); 
            }

            if (test_digest && z_file->genozip_version >= 9) {  // backward compatability with v8: we don't test against v8 MD5 for the header, as we had a bug in v8 in which we included a junk MD5 if they user didn't --md5 or --test. any file integrity problem will be discovered though on the whole-file MD5 so no harm in skipping this.
                Digest reconstructed_header_digest = digest_do (comp_vb->txt_data.data, comp_vb->txt_data.len);
                
                TEMP_FLAG (quiet, flag.quiet && !flag.show_digest);

                ASSERTW (digest_is_equal (header->digest_header, DIGEST_NONE) || 
                        digest_is_equal (reconstructed_header_digest, header->digest_header) ||
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
    }

    // case: component is not in plan - discard the VB
    else
        vb_release_vb (&comp_vb);    
    
    if (!flag.reading_chain && !flag.reading_reference)
        is_first_txt = false;
}