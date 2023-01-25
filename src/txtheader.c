// ------------------------------------------------------------------
//   txtheader.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "genozip.h"
#include "buffer.h"
#include "digest.h"
#include "flags.h"
#include "data_types.h"
#include "filename.h"
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
#include "compressor.h"
#include "dispatcher.h"

static bool is_first_txt = true; 

//----------
// ZIP stuff
//----------

static SectionHeaderTxtHeader section_header = {};
static BufferP txt_header_buf = NULL;
static CompIType txt_header_comp_i = 0;

#define TXT_HEADER_VB_SIZE (1<<24) // 16MB

// compress the Txt Header fragment - to allow multi-threaded compression and decompression
static void txtheader_prepare_for_compress (VBlockP vb)
{
    vb->fragment_start = txt_header_buf->len ? Bc (*txt_header_buf, txt_header_buf->next) : NULL;
    vb->fragment_len   = MIN_(TXT_HEADER_VB_SIZE, txt_header_buf->len - txt_header_buf->next);
    vb->comp_i         = txt_header_comp_i; // goes into SectionEntFileFormat.comp_i
    
    txt_header_buf->next += vb->fragment_len;

    // we always compress vb_1 even if its empty - so that the TxtHeader section gets written even if there's no header (eg headerless SAM)
    if (vb->fragment_len || vb->vblock_i == 1) vb->dispatch = READY_TO_COMPUTE;
}

static void txtheader_compress_one_fragment (VBlockP vb)
{
    START_TIMER;

    SectionHeaderTxtHeader my_header = section_header; // duplicating the ~300 bytes header in each fragment, but this negligible vs the fragment length
    my_header.data_uncompressed_len = BGEN32 (vb->fragment_len);
    my_header.vblock_i              = BGEN32 (vb->vblock_i); // up to 14.0.8 (when TXT_HEADER fragmantization was introduced), this was always 0 for SEC_TXT_HEADER

    comp_compress (vb, NULL, &vb->z_data, (SectionHeader*)&my_header, vb->fragment_start, NO_CALLBACK, "SEC_TXT_HEADER");

    // copy first fragment header, for updating in zfile_update_txt_header_section_header
    if (vb->vblock_i == 1)
        memcpy (&z_file->txt_header_single, &my_header, sizeof (SectionHeaderTxtHeader));

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.

    COPY_TIMER (txtheader_compress_one_fragment)    
}

void txtheader_compress (BufferP txt_header, 
                         uint64_t unmodified_txt_header_len, // length of header before modifications, eg due to --chain or compressing a Luft file
                         Digest header_md5, bool is_first_txt, CompIType comp_i)
{
    START_TIMER;

    Codec codec = codec_assign_best_codec (evb, NULL, txt_header, SEC_TXT_HEADER);

    section_header = (SectionHeaderTxtHeader){
        .magic             = BGEN32 (GENOZIP_MAGIC),
        .section_type      = SEC_TXT_HEADER,
        .compressed_offset = BGEN32 (sizeof (SectionHeaderTxtHeader)),
        .codec             = (codec == CODEC_UNKNOWN) ? CODEC_NONE : codec,
        .src_codec         = txt_file->codec, 
        .digest_header     = flag.data_modified ? DIGEST_NONE : header_md5,
        .txt_header_size   = BGEN64 (unmodified_txt_header_len),
    };

    // data type specific fields
    if (DTPZ(zip_set_txt_header_specific)) DTPZ(zip_set_txt_header_specific)(&section_header);

    // In BGZF, we store the 3 least significant bytes of the file size, so check if the reconstructed BGZF file is likely the same
    if (txt_file->codec == CODEC_BGZF) 
        bgzf_sign (txt_file->disk_size, section_header.codec_info);
        
    filename_base (txt_file->name, false, FILENAME_STDIN, section_header.txt_filename, TXT_FILENAME_LEN);
    file_remove_codec_ext (section_header.txt_filename, txt_file->type); // eg "xx.fastq.gz -> xx.fastq"

    txt_header_buf = txt_header;
    txt_header_buf->next = 0;
    txt_header_comp_i = comp_i;

    dispatcher_fan_out_task ("compress_txt_header", NULL, 0, "Writing txt header...", false, false, false, 0, 20000,
                             txtheader_prepare_for_compress, 
                             txtheader_compress_one_fragment, 
                             zfile_output_processed_vb);

    // VCF note: we don't account for DVCF rejects components - the added header lines are duplicates of the main header
    if (!(Z_DT(VCF) || Z_DT(BCF)) || comp_i == VCF_COMP_MAIN) {        
        z_file->txt_data_so_far_single   += txt_header->len; // length of txt header as it would be reconstructed (possibly afer modifications)
        z_file->txt_data_so_far_bind     += txt_header->len;
        z_file->txt_data_so_far_single_0 += unmodified_txt_header_len; // length of the original txt header as read from the file
        z_file->txt_data_so_far_bind_0   += unmodified_txt_header_len;
    }

    z_file->txt_data_so_far_bind_comp[comp_i] += txt_header->len;
    z_file->txt_data_so_far_bind_0_comp[comp_i] += unmodified_txt_header_len;

    COPY_TIMER_VB (evb, txtheader_compress);
}

// ZIP: reads txt header and writes its compressed form to the GENOZIP file
int64_t txtheader_zip_read_and_compress (int64_t *txt_header_offset, CompIType comp_i) // out (-1 if SEC_TXT_HEADER not written)
{    
    START_TIMER;

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
    if (z_file && !flag.seg_only && !flag.make_reference && !(z_sam_gencomp && (comp_i == SAM_COMP_PRIM || comp_i == SAM_COMP_DEPN))) {
        *txt_header_offset = z_file->disk_so_far; // offset of first (vb=1) TXT_HEADER fragment
        txtheader_compress (&evb->txt_data, txt_header_size, header_digest, is_first_txt, comp_i); 
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

    COPY_TIMER_VB (evb, txtheader_zip_read_and_compress);

    return txt_header_size; // all good
}

//----------
// PIZ stuff
//----------

static Section txtheader_sec = NULL; 
static VBlockP txt_header_vb = NULL;

static void txtheader_read_one_vb (VBlockP vb)
{
    bool is_new_header = (txtheader_sec->st == SEC_TXT_HEADER) && (txtheader_sec->vblock_i == vb->vblock_i); // 14.0.9 or later
    bool is_old_header = (txtheader_sec->st == SEC_TXT_HEADER) && (txtheader_sec->vblock_i == 0) && (vb->vblock_i == 1); // up tp 14.0.8 - TXT_HEADER is always a single section with vblock_i=0
    
    if (!is_new_header && !is_old_header) return; // this is not a TXT_HEADER section at all, or not one that belongs to this component
    
    zfile_read_section (z_file, vb, txtheader_sec->vblock_i, &vb->z_data, "z_data", SEC_TXT_HEADER, txtheader_sec);    
    SectionHeaderTxtHeader *header = B1ST (SectionHeaderTxtHeader, vb->z_data);

    vb->fragment_len   = BGEN32 (header->data_uncompressed_len);
    vb->fragment_start = Bc (txt_header_vb->txt_data, txt_header_vb->txt_data.next);

    ASSERT (txt_header_vb->txt_data.next + vb->fragment_len <= txt_header_vb->txt_data.len, "TxtHeader fragments exceed length=%"PRIu64, txt_header_vb->txt_data.len); 

    // increment for next fragment
    txt_header_vb->txt_data.next += vb->fragment_len;
    txtheader_sec++;

    vb->dispatch = READY_TO_COMPUTE;
}

// entry point of compute thread of dictionary decompression
static void txtheader_uncompress_one_vb (VBlockP vb)
{
    SectionHeaderTxtHeader *header = B1ST (SectionHeaderTxtHeader, vb->z_data);
    zfile_uncompress_section_into_buf (vb, header, BGEN32 (header->vblock_i), SEC_TXT_HEADER, &txt_header_vb->txt_data, vb->fragment_start);

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}

// PIZ main thread: reads the txt header from the genozip file and outputs it to the reconstructed txt file
void txtheader_piz_read_and_reconstruct (Section sec)
{
    START_TIMER;

    z_file->disk_at_beginning_of_this_txt_file = z_file->disk_so_far;

    txtheader_sec = sec;

    txt_header_vb = vb_get_vb (POOL_MAIN, PIZ_TASK_NAME, 0, sec->comp_i);

    // allocate memory for uncompressed text header - the size is the sum of the fragment sizes
    Section my_sec = sec;
    SectionHeaderTxtHeader header = {};
    for (VBIType vb_i=1; my_sec->st == SEC_TXT_HEADER && (my_sec->vblock_i==vb_i || (!my_sec->vblock_i && vb_i==1 /* up to 14.0.8 */)); vb_i++, my_sec++) {
        SectionHeaderTxtHeader frag_header = zfile_read_section_header (txt_header_vb, my_sec->offset, my_sec->vblock_i, SEC_TXT_HEADER).txt_header;
        txt_header_vb->txt_data.len += BGEN32 (frag_header.data_uncompressed_len);

        if (vb_i == 1 && sec->comp_i == COMP_MAIN) 
            z_file->txt_header_single = frag_header;

        if (vb_i == 1)
            header = frag_header;
    }

    ASSERT0 (header.section_type == SEC_TXT_HEADER, "SEC_TXT_HEADER section not found");

    buf_alloc (txt_header_vb, &txt_header_vb->txt_data, 0, txt_header_vb->txt_data.len, char, 0, "txt_data");
    txt_header_vb->txt_data.next = 0;

    dispatcher_fan_out_task ("read_txt_header", NULL, 0, 0, true, flag.test, false, 0, 1000,
                             txtheader_read_one_vb, 
                             txtheader_uncompress_one_vb,
                             NO_CALLBACK);
        
    // 1. if flag.unbind (genounzip) - we open the output txt file of the component
    // 2. if flag.one_component (genocat) - output to stdout or --output
    // 3. when reading an auxiliary file or no_writer- we create txt_file here (but don't actually open the physical file)
    if (!txt_file) { 
        rom filename = flag.unbind ? txtfile_piz_get_filename (header.txt_filename, flag.unbind, false) 
                                   : flag.out_filename;
        txt_file = file_open (filename, WRITE, TXT_FILE, flag.out_dt != DT_NONE ? flag.out_dt : z_file->data_type);
        if (flag.unbind) FREE (filename); // file_open copies the names
        
        z_file->digest_ctx = DIGEST_CONTEXT_NONE; // reset digest
    }

    // initialize if needed - but only once per outputted txt file 
    //i.e. if we have rejects+normal, or concatenated, we will only init in the first)
    if (!txt_file->piz_header_init_has_run && DTPZ(piz_header_init))
        txt_file->piz_header_init_has_run = true;

    evb->comp_i                = !VER(14) ? header.flags.txt_header.v13_dvcf_comp_i/*v12,13*/ : sec->comp_i/*since v14*/;
    txt_file->max_lines_per_vb = BGEN32 (header.max_lines_per_vb);
    txt_file->txt_flags        = header.flags.txt_header;
    txt_file->num_vbs          = sections_count_sections_until (SEC_VB_HEADER, sec, SEC_TXT_HEADER);
    
    if (txt_file->codec == CODEC_BGZF)
        memcpy (txt_file->bgzf_signature, header.codec_info, 3);
        
    // SAM: only the main component has a TxtHeader section. 
    // DVCF: always zero. 
    // FASTQ: flag.data_modified when interleaving, or just one file if --R1/2. 
    // V8: zero if not compressed with --md5
    if (!digest_is_zero(header.digest) && !flag.data_modified) 
        z_file->digest = header.digest; 

    ASSINP (!flag.test || !digest_is_zero(header.digest) || z_file->z_flags.adler, 
            "--test cannot be used with %s, as it was compressed without a digest. See " WEBSITE_DIGEST, z_name);

    // case: we need to reconstruct (or not) the BGZF following the instructions from the z_file
    if (flag.bgzf == BGZF_BY_ZFILE) {

        // load the source file isize if we have it and we are attempting to reconstruct an unmodifed file identical to the source
        bool loaded = false;
        if (!flag.data_modified)     
//xxx        && (z_file->num_components == 1 || flag.unbind))  // not concatenating multiple files
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

    // count header-lines (for --lines etc): before data-modifying inspect_txt_header
    if (writer_does_txtheader_need_write (sec)) {
        if (flag.header_one && (Z_DT(VCF) || Z_DT(BCF)))
            txt_file->num_lines += 1;
        else if (!DTPT (is_binary))
            txt_file->num_lines += str_count_char (STRb(txt_header_vb->txt_data), '\n'); // number of source-file lines
    }

    if (txt_header_vb->txt_data.len)
        DT_FUNC_OPTIONAL (z_file, inspect_txt_header, true)(txt_header_vb, &txt_header_vb->txt_data, header.flags.txt_header); // ignore return value

    // hand-over txt header if it is needed (it won't be if flag.no_header)
    if (writer_does_txtheader_need_write (sec)) {

        // if we're translating from one data type to another (SAM->BAM, BAM->FASTQ, ME23->VCF etc) translate the txt header 
        // note: in a header-less SAM, after translating to BAM, we will have a header
        DtTranslation trans = dt_get_translation (NULL);
        if (trans.txtheader_translator && !flag.no_header) trans.txtheader_translator (txt_header_vb, &txt_header_vb->txt_data); 

        // count textual lines in header, used for line= reporting in ASSPIZ
        if (txt_header_vb->txt_data.len && !DTPT (is_binary)) {
            uint32_t num_textual_lines = str_count_char (STRb(txt_header_vb->txt_data), '\n');
            writer_set_num_txtheader_lines (sec->comp_i, num_textual_lines);
        }

        if (piz_need_digest) 
            digest_txt_header (&txt_header_vb->txt_data, header.digest_header);

        writer_handover_txtheader (&txt_header_vb); // handover data to writer thread (even if the header is empty, as the writer thread is waiting for it)

        // accounting for data as in original source file - affects vb->vb_position_txt_file of next VB
        txt_file->txt_data_so_far_single_0 = !VER(12) ? BGEN32 (header.data_uncompressed_len) : BGEN64 (header.txt_header_size); 
    }

    // case: component is not in plan - discard the VB
    else
        vb_release_vb (&txt_header_vb, PIZ_TASK_NAME);    
    
    if (!flag.reading_chain && !flag.reading_reference)
        is_first_txt = false;

    COPY_TIMER_VB (evb, txtheader_piz_read_and_reconstruct);
}
