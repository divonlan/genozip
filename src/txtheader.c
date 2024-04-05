// ------------------------------------------------------------------
//   txtheader.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
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

    comp_compress (vb, NULL, &vb->z_data, &my_header, vb->fragment_start, NO_CALLBACK, "SEC_TXT_HEADER");

    // copy first fragment header, for updating in zfile_update_txt_header_section_header
    if (vb->vblock_i == 1)
        z_file->txt_header_hdr = my_header;

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.

    COPY_TIMER (txtheader_compress_one_fragment)    
}

void txtheader_compress (BufferP txt_header, 
                         uint64_t unmodified_txt_header_len, // length of header before modifications, eg due to --optimize
                         Digest header_md5, bool is_first_txt, CompIType comp_i)
{
    START_TIMER;

    Codec codec = codec_assign_best_codec (evb, NULL, txt_header, SEC_TXT_HEADER);

    section_header = (SectionHeaderTxtHeader){
        .magic             = BGEN32 (GENOZIP_MAGIC),
        .section_type      = SEC_TXT_HEADER,
        .codec             = (codec == CODEC_UNKNOWN) ? CODEC_NONE : codec,
        .src_codec         = txt_file->source_codec, 
        .digest_header     = flag.zip_txt_modified ? DIGEST_NONE : header_md5,
        .txt_header_size   = BGEN64 (unmodified_txt_header_len), // length before zip-side modifications
    };

    // data type (of txt file!) specific fields
    if (DTPT(zip_set_txt_header_flags)) DTPT(zip_set_txt_header_flags)(&section_header.flags.txt_header);

    // true if filename is compressed with gz/bgzf but does not have a .gz/.bgz extension (e.g. usually true for BAM files)
    section_header.flags.txt_header.no_gz_ext = 
        (SOURCE_CODEC(GZ) || SOURCE_CODEC(BGZF) || SOURCE_CODEC(BAM) || SOURCE_CODEC(CRAM)/*reconstructed as BAM*/) &&
        txt_file->basename && !filename_has_ext (txt_file->basename, ".gz") && !filename_has_ext (txt_file->basename, ".bgz");

    // In BGZF, we store the 3 least significant bytes of the file size, so check if the reconstructed BGZF file is likely the same
    if (txt_file->codec == CODEC_BGZF) 
        bgzf_sign (txt_file->disk_size, section_header.codec_info);
        
    filename_base (txt_file->name, false, FILENAME_STDIN, section_header.txt_filename, TXT_FILENAME_LEN);
    filename_remove_codec_ext (section_header.txt_filename, txt_file->type); // eg "xx.fastq.gz -> xx.fastq"

    txt_header_buf = txt_header;
    txt_header_buf->next = 0;
    txt_header_comp_i = comp_i;

    dispatcher_fan_out_task ("compress_txt_header", NULL, 0, "Writing txt header...", false, false, false, 0, 20000, true,
                             txtheader_prepare_for_compress, 
                             txtheader_compress_one_fragment, 
                             zfile_output_processed_vb);

    z_file->txt_data_so_far_single   += txt_header->len; // length of txt header as it would be reconstructed (possibly after modifications)
    z_file->txt_data_so_far_bind     += txt_header->len;
    z_file->txt_data_so_far_single_0 += unmodified_txt_header_len; // length of txt header without applying piz-side modifications
    z_file->txt_data_so_far_bind_0   += unmodified_txt_header_len;

    z_file->txt_data_so_far_bind_comp[comp_i] += txt_header->len;
    z_file->txt_data_so_far_bind_0_comp[comp_i] += unmodified_txt_header_len;

    COPY_TIMER_EVB (txtheader_compress);
}

// ZIP: reads txt header and writes its compressed form to the GENOZIP file
int64_t txtheader_zip_read_and_compress (int64_t *txt_header_offset, CompIType comp_i) // out (-1 if SEC_TXT_HEADER not written)
{    
    START_TIMER;

    Digest header_digest = DIGEST_NONE;
    evb->comp_i = comp_i; // used by def_is_header_done
    
    TxtHeaderRequirement req = DTPT (txt_header_required); 
    if (req == HDR_MUST || req == HDR_OK || (comp_i==0 && (req == HDR_MUST_0 || req == HDR_OK_0))) 
        txtfile_read_header (is_first_txt); // reads into evb->txt_data and evb->lines.len

    // get header digest and initialize component digest
    if (zip_need_digest) 
        header_digest = digest_txt_header (&evb->txt_data, DIGEST_NONE, comp_i);

    // for VCF, we need to check if the samples are the same before approving binding (other data types can bind without restriction)
    //          also: header is modified if --chain or compressing a Luft file
    // for SAM, we check that the contigs specified in the header are consistent with the reference given in --reference/--REFERENCE
    uint64_t txt_header_size = evb->txt_data.len; // note: 0 if BAM with --show-bam
    if (txt_header_size && !(DT_FUNC_OPTIONAL (txt_file, inspect_txt_header, true)(evb, &evb->txt_data, (struct FlagsTxtHeader){}))) { 
        buf_free (evb->txt_data);
        return -1;
    }

    // note: we always write the txt_header for comp_i=0 even if we don't actually have a header, because the
    // section header contains the data about the file. Special case: we don't write headers of SAM PRIM/DEPN
    if (z_file && !flag.zip_no_z_file && !flag.make_reference && !(z_sam_gencomp && (comp_i == SAM_COMP_PRIM || comp_i == SAM_COMP_DEPN))) {
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
        z_file->header_size += txt_file->header_size; // note: header_size already contains the length difference if --match-chrom

    z_file->num_txts_so_far++; // when compressing

    buf_free (evb->txt_data);
    
    is_first_txt = false;

    COPY_TIMER_EVB (txtheader_zip_read_and_compress);

    return txt_header_size; // all good
}

//----------
// PIZ stuff
//----------

static Section txtheader_sec = NULL; 
static VBlockP txt_header_vb = NULL;
static uint64_t sum_fragment_len = 0;

static void txtheader_read_one_vb (VBlockP vb)
{
    bool is_new_header = (txtheader_sec->st == SEC_TXT_HEADER) && (txtheader_sec->vblock_i == vb->vblock_i); // 14.0.9 or later
    bool is_old_header = (txtheader_sec->st == SEC_TXT_HEADER) && (txtheader_sec->vblock_i == 0) && (vb->vblock_i == 1); // up tp 14.0.8 - TXT_HEADER is always a single section with vblock_i=0
    
    if (!is_new_header && !is_old_header) return; // this is not a TXT_HEADER section at all, or not one that belongs to this component
    
    zfile_read_section (z_file, vb, txtheader_sec->vblock_i, &vb->z_data, "z_data", SEC_TXT_HEADER, txtheader_sec);    
    SectionHeaderTxtHeaderP header = B1ST (SectionHeaderTxtHeader, vb->z_data);

    vb->fragment_len   = BGEN32 (header->data_uncompressed_len);
    vb->fragment_start = Bc (txt_header_vb->txt_data, txt_header_vb->txt_data.next);

    // accounting for data after zip-side modifications but before piz-side modifications - 
    // affects vb->vb_position_txt_file of next VB
    sum_fragment_len += vb->fragment_len; 

    ASSERT (txt_header_vb->txt_data.next + vb->fragment_len <= txt_header_vb->txt_data.len, "TxtHeader fragments exceed length=%"PRIu64, txt_header_vb->txt_data.len); 

    // increment for next fragment
    txt_header_vb->txt_data.next += vb->fragment_len;
    txtheader_sec++;

    vb->dispatch = READY_TO_COMPUTE;
}

// entry point of compute thread of dictionary decompression
static void txtheader_uncompress_one_vb (VBlockP vb)
{
    SectionHeaderTxtHeaderP header = B1ST (SectionHeaderTxtHeader, vb->z_data);
    zfile_uncompress_section_into_buf (vb, header, BGEN32 (header->vblock_i), SEC_TXT_HEADER, &txt_header_vb->txt_data, vb->fragment_start);

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}

// get filename of output txt file in genounzip if user didn't specific it with --output
// case 1: outputing a single file - generate txt_filename based on the z_file's name
// case 2: unbinding a genozip into multiple txt files - generate txt_filename of a component file from the
//         component name in SEC_TXT_HEADER 
static rom txtheader_piz_get_filename (rom orig_name, rom prefix, bool is_orig_name_genozip, bool with_bgzf, bool no_gz_ext)
{
    unsigned fn_len = strlen (orig_name);
    unsigned dn_len = flag.out_dirname ? strlen (flag.out_dirname) : 0;
    unsigned px_len = prefix ? strlen (prefix) : 0;
    unsigned genozip_ext_len = is_orig_name_genozip ? (sizeof GENOZIP_EXT - 1) : 0;
    int txt_filename_size = fn_len + dn_len + px_len + 10;
    char *txt_filename = (char *)CALLOC(txt_filename_size);

    #define EXT2_MATCHES_TRANSLATE(from,to,ext)  \
        ((z_file->data_type==(from) && flag.out_dt==(to) && \
         fn_len >= genozip_ext_len+strlen(ext) && \
         !strcmp (&txt_filename[fn_len-genozip_ext_len- (sizeof ext - 1)], (ext))) ? (int)(sizeof ext - 1) : 0) 

    // length of extension to remove if translating, eg remove ".sam" if .sam.genozip->.bam */
    int old_ext_removed_len = EXT2_MATCHES_TRANSLATE (DT_SAM,   DT_BAM,    ".sam")
                            + EXT2_MATCHES_TRANSLATE (DT_SAM,   DT_SAM,    ".bam")
                            + EXT2_MATCHES_TRANSLATE (DT_SAM,   DT_FASTQ,  ".sam")
                            + EXT2_MATCHES_TRANSLATE (DT_SAM,   DT_FASTQ,  ".bam")
                            + EXT2_MATCHES_TRANSLATE (DT_VCF,   DT_BCF,    ".vcf")
                            + EXT2_MATCHES_TRANSLATE (DT_ME23,  DT_VCF,    ".txt");

    // case: new directory - take only the basename
    if (dn_len) orig_name = filename_base (orig_name, 0,0,0,0); 
    
    // cases in which BGZF compression does not result in a ".gz" file name extension:
    // if original gz/bgzf-compressed file did not have a .gz/.bgz extension (since 15.0.23) or when outputing BAM or BCF
    no_gz_ext |= OUT_DT(BCF) || OUT_DT(BAM); 

    snprintf ((char *)txt_filename, txt_filename_size, "%s%s%s%.*s%s%s", prefix ? prefix : "",
             (dn_len ? flag.out_dirname : ""), (dn_len ? "/" : ""), 
             fn_len - genozip_ext_len - old_ext_removed_len, orig_name,
             old_ext_removed_len ? file_plain_ext_by_dt (flag.out_dt) : "", // add translated extension if needed
             (with_bgzf && !no_gz_ext) ? ".gz" : ""); // add .gz if needed

    if (dn_len) FREE (orig_name); // allocated by filename_base

    return txt_filename;
}

// get filename, even if txt_file has might not been open. 
StrTextLong txtheader_get_txt_filename_from_section (void)
{
    Section sec;

    if (flag.one_component && !flag.deep)
        sec = sections_get_comp_txt_header_sec (flag.one_component - 1);
    else 
        sec = sections_get_comp_txt_header_sec (COMP_MAIN);

    SectionHeaderTxtHeader header = zfile_read_section_header (evb, sec, SEC_TXT_HEADER).txt_header;

    // note: for bz2, xz, and zip - we reconstruct as gz too. better choice than plain.
    #define C(cdc) (header.src_codec == CODEC_##cdc)
    bool dot_gz = flag.bgzf == BGZF_NOT_INITIALIZED ? ((C(BGZF) || C(GZ) || C(BZ2) || C(XZ) || C(ZIP))) // note: similar logic to in bgzf_piz_calculate_bgzf_flags
                :                                     (flag.bgzf != 0);
    #undef C

    TEMP_FLAG(out_dirname, NULL); // only basename in progress string
    rom filename = txtheader_piz_get_filename (header.txt_filename, flag.unbind, false, dot_gz, header.flags.txt_header.no_gz_ext);
    RESTORE_FLAG(out_dirname);

    StrTextLong name = {};
    strncpy (name.s, filename, sizeof (StrTextLong)-1);
    FREE (filename);

    return name;
}          

// PIZ main thread: reads the txt header from the genozip file and outputs it to the reconstructed txt file
void txtheader_piz_read_and_reconstruct (Section sec)
{
    START_TIMER;

    txtheader_sec = sec;

    txt_header_vb = vb_get_vb (POOL_MAIN, PIZ_TASK_NAME, 0, sec->comp_i);

    // allocate memory for uncompressed text header - the size is the sum of the fragment sizes
    Section my_sec = sec;
    SectionHeaderTxtHeader header = {};
    for (VBIType vb_i=1; my_sec->st == SEC_TXT_HEADER && (my_sec->vblock_i==vb_i || (!my_sec->vblock_i && vb_i==1 /* up to 14.0.8 */)); vb_i++, my_sec++) {
        SectionHeaderTxtHeader frag_header = zfile_read_section_header (txt_header_vb, my_sec, SEC_TXT_HEADER).txt_header;
        txt_header_vb->txt_data.len += BGEN32 (frag_header.data_uncompressed_len);

        if (vb_i == 1 && sec->comp_i == COMP_MAIN) 
            memcpy (z_file->txt_filename, frag_header.txt_filename, TXT_FILENAME_LEN);

        if (vb_i == 1)
            header = frag_header;
    }

    ASSERT0 (header.section_type == SEC_TXT_HEADER, "SEC_TXT_HEADER section not found");

    buf_alloc (txt_header_vb, &txt_header_vb->txt_data, 0, txt_header_vb->txt_data.len, char, 0, "txt_data");
    txt_header_vb->txt_data.next = 0;

    // read actual txt header data only if we need to reconstruct the header
    bool needs_recon = writer_does_txtheader_need_recon (sec->comp_i);
    FlagsBgzf bgzf_flags = {};
    rom filename = NULL;
    sum_fragment_len = 0;

    if (needs_recon) {
        dispatcher_fan_out_task ("read_txt_header", NULL, 0, 0, true, flag.test, false, 0, 1000, true,
                                 txtheader_read_one_vb, 
                                 txtheader_uncompress_one_vb,
                                 NO_CALLBACK);

        bgzf_flags = bgzf_piz_calculate_bgzf_flags (sec->comp_i, header.src_codec);

        filename = flag.to_stdout    ? NULL 
                 : flag.out_filename ? flag.out_filename
                 : flag.unbind       ? txtheader_piz_get_filename (header.txt_filename, flag.unbind, false, bgzf_flags.level >= 1, header.flags.txt_header.no_gz_ext)
                 :                     txtheader_piz_get_filename (z_name,              "",          true,  bgzf_flags.level >= 1, header.flags.txt_header.no_gz_ext); // use genozip filename as a base regardless of original name
    }

    // note: when reading an auxiliary file or no_writer - we still create txt_file (but don't actually open the physical file)
    // note: if there are several components contributing to a single txt_file (DVCF, SAM w/gencomp) - we only open it once
    if (!txt_file)
        txt_file = file_open_txt_write (filename, flag_loading_auxiliary ? z_file->data_type : flag.out_dt, bgzf_flags.level);
    
    if (!flag.to_stdout && !flag.out_filename) FREE (filename); // file_open_z copies the names

    // set BGZF info in txt_file - either that originates from SEC_BGZF, or constructed based on bgzf_flags
    if (needs_recon && txt_file->codec == CODEC_BGZF)
        bgzf_piz_set_txt_file_bgzf_info (bgzf_flags, header.codec_info);

    // note: this is reset for each component:
    // since v14 it is used for the commulative component-scope MD5 used for both VBs and txt file verification
    // up to v13 it is for commulative digest of both MD5 and Adler, but used only for txt file verification,
    // while VB verification relies on v13_commulative_digest_ctx which is commulative across components.
    z_file->digest_ctx = DIGEST_CONTEXT_NONE; // reset digest

    if (!VER(14) && header.flags.txt_header.v13_dvcf_comp_i)
        ABORT0 ("DVCF files can be accessed with Genozip up to version 15.0.41");

    evb->comp_i                = VER(14) ? sec->comp_i/*since v14*/ : 0;
    txt_file->max_lines_per_vb = BGEN32 (header.max_lines_per_vb);
    txt_file->txt_flags        = header.flags.txt_header;
    txt_file->num_vbs          = sections_count_sections_until (SEC_VB_HEADER, sec, SEC_TXT_HEADER);    
    txt_file->txt_data_size    = BGEN64 (header.txt_data_size);
    txt_file->txt_data_so_far_single_0 = sum_fragment_len;

    if (VER(15)) 
        for (QType q=0; q < NUM_QTYPES; q++)
            segconf.flav_prop[q] = header.flav_prop[q];

    // initialize if needed - but only once per outputted txt file 
    // i.e. if we have rejects+normal, or concatenated, we will only init in the first)
    if (!txt_file->piz_header_init_has_run && DTPZ(piz_header_init)) {
        DTPZ(piz_header_init)();
        txt_file->piz_header_init_has_run = true;
    }
    
    bool needs_write = writer_does_txtheader_need_write (sec->comp_i);

    // count header-lines (for --lines etc): before data-modifying inspect_txt_header
    if (needs_write) {
        if (flag.header_one && (Z_DT(VCF) || Z_DT(BCF)))
            txt_file->num_lines += 1;
        else if (!DTPT (is_binary))
            txt_file->num_lines += str_count_char (STRb(txt_header_vb->txt_data), '\n'); // number of source-file lines
    }

    if (txt_header_vb->txt_data.len)
        DT_FUNC_OPTIONAL (z_file, inspect_txt_header, true)(txt_header_vb, &txt_header_vb->txt_data, header.flags.txt_header); // ignore return value

    // hand-over txt header if it is needed (it won't be if flag.no_header)
    if (needs_write) {

        // if we're translating from one data type to another (SAM->BAM, BAM->FASTQ, ME23->VCF etc) translate the txt header 
        // note: in a header-less SAM, after translating to BAM, we will have a header
        DtTranslation trans = dt_get_translation (NULL);
        if (trans.txtheader_translator && !flag.no_header) 
            trans.txtheader_translator (txt_header_vb, &txt_header_vb->txt_data); 

        // count textual lines in header, used for line= reporting in ASSPIZ
        if (txt_header_vb->txt_data.len && !DTPT (is_binary)) {
            uint32_t num_textual_lines = str_count_char (STRb(txt_header_vb->txt_data), '\n');
            writer_set_num_txtheader_lines (sec->comp_i, num_textual_lines);
        }
    }
    
    if (piz_need_digest) {
        // store txt-file-wide digest. If its 0, then the file-wide digest is not coming from this component
        // e.g. non-main component in SAM. If is also 0 if file not digested (--optimize, DVCF, v8 without --md5/--test...) or, since v14, if Adler32
        if (!digest_is_zero(header.digest)) 
            z_file->digest = header.digest; 

        digest_txt_header (&txt_header_vb->txt_data, header.digest_header, sec->comp_i); // verify txt header digest
    }

    if (!writer_handover_txtheader (&txt_header_vb)) {  // handover data to writer thread (even if the header is empty, as the writer thread is waiting for it)
        txt_file->txt_data_so_far_single += txt_header_vb->txt_data.len; // if writing, this is done in writer_write, caputring the processing in writer too

        vb_release_vb (&txt_header_vb, PIZ_TASK_NAME); // not handing over, so release here      
    }

    if (!flag.reading_reference)
        is_first_txt = false;

    COPY_TIMER_EVB (txtheader_piz_read_and_reconstruct);
}
