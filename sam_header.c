// ------------------------------------------------------------------
//   sam_header.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

// ----------------------
// SAM / BAM header stuff
// ----------------------

#include "libdeflate/libdeflate.h"
#include <time.h>
#if defined __APPLE__ 
#include "compatibility/mac_gettime.h"
#endif
#include "sam_private.h"
#include "file.h"
#include "random_access.h"
#include "reference.h"
#include "zfile.h"
#include "version.h"
#include "endianness.h"
#include "website.h"
#include "segconf.h"
#include "contigs.h"
#include "flags.h"
#include "buffer.h"
#include "strings.h"

// private globals
ContigPkgP sam_hdr_contigs = NULL; // points to contigs for: BAM (always) and headerful SAM and VCF

#define foreach_SQ_line(cb,cb_param) (IS_BAM ? bam_foreach_SQ_line (txt_header->data, txt_header->len, cb, cb_param) : sam_foreach_SQ_line(txt_header->data, cb, cb_param))

//---------------
// Header contigs
//---------------

// call a callback for each SQ line (contig). Note: callback function is the same as foreach_contig
void sam_foreach_SQ_line (const char *txt_header, // nul-terminated string
                          ContigsIteratorCallback callback, void *callback_param)
{
    const char *line = txt_header;

    while (1) {
        line = strchr (line, '@');
        if (!line) break;

        if (line[1] == 'S' && line[2] == 'Q') { // this test will always be in the buffer - possible in the overflow area
            char *newline = strchr (line, '\n');
            *newline = 0;

            // line looks something like: @SQ\tSN:chr10\tLN:135534747 (no strict order between SN and LN and there could be more parameters)
            const char *chrom_name   = strstr (line, "SN:"); 
            const char *last_pos_str = strstr (line, "LN:"); 

            if (chrom_name && last_pos_str) {
                unsigned chrom_name_len = strcspn (&chrom_name[3], "\t\n\r");
                PosType last_pos = (PosType)strtoull (&last_pos_str[3], NULL, 10);

                ASSINP (last_pos <= MAX_POS_SAM, "Error: @SQ record in header contains LN:%"PRId64" which is beyond the maximum permitted by the SAM spec of %"PRId64,
                        last_pos, MAX_POS_SAM);

                callback (&chrom_name[3], chrom_name_len, last_pos, callback_param);
            }

            line = newline+1;

            *newline = '\n'; // restore
        }
        else 
            line++;
    }
}

// call a callback for each SQ line (contig). Note: callback function is the same as foreach_contig
static void bam_foreach_SQ_line (const char *txt_header, uint32_t txt_header_len, // binary BAM header
                                 ContigsIteratorCallback callback, void *callback_param)
{
    #define HDRSKIP(n) if (txt_header_len < next + n) goto incomplete_header; next += n
    #define HDR32 (next + 4 <= txt_header_len ? GET_UINT32 (&txt_header[next]) : 0) ; if (txt_header_len < next + 4) goto incomplete_header; next += 4;

    uint32_t next=0;

    // skip magic, l_text, text
    HDRSKIP(4);
    uint32_t l_text = HDR32;
    HDRSKIP(l_text);

    uint32_t n_ref = HDR32;

    for (uint32_t ref_i=0; ref_i < n_ref; ref_i++) {
        uint32_t l_name = HDR32;
        const char *name = &txt_header[next];
        HDRSKIP (l_name);
        uint32_t l_ref = HDR32;

        callback (name, l_name-1, l_ref, callback_param);
    }
    return;

incomplete_header:
    ABORT ("Error in bam_foreach_SQ_line: incomplete BAM header (next=%u)", next);

    #undef HDRSKIP
    #undef HDR32
}   

typedef struct { uint32_t num_contigs, dict_len; } ContigsCbParam;
static ContigP sam_header_count_contigs_cb (const char *chrom_name, unsigned chrom_name_len, PosType last_pos, void *cb_param)
{
    // note: for simplicity, we consider contigs to be in the range [0, last_pos] even though POS=0 doesn't exist - last_pos+1 loci
    ((ContigsCbParam*)cb_param)->num_contigs++;
    ((ContigsCbParam*)cb_param)->dict_len += chrom_name_len + 1;

    return NULL;
}

static void sam_header_alloc_contigs (BufferP txt_header)
{
    sam_header_finalize(); // free previous header

    ContigsCbParam contigs_count = {}; 
    foreach_SQ_line (sam_header_count_contigs_cb, &contigs_count);

    // note: In BAM, a contig-less header is considered an unaligned BAM in which all RNAME=*
    // In SAM, a contig-less header just means we don't know the contigs and they are defined in RNAME 
    if (!IS_BAM && !contigs_count.num_contigs) return; // "headerless" SAM (not really headerless, just no contigs in the header)

    sam_hdr_contigs = CALLOC (sizeof (ContigPkg));
    sam_hdr_contigs->name = "sam_hdr_contigs";
    
    buf_alloc (evb, &sam_hdr_contigs->contigs, 0, contigs_count.num_contigs, Contig, 1, "JJContigPkg->contigs"); 
    buf_alloc (evb, &sam_hdr_contigs->dict,    0, contigs_count.dict_len, char, 1, "JJContigPkg->dict"); 
}

void sam_header_finalize (void)
{
    if (!flag.bind) { // in ZIP with bind we keep the header contigs
        contigs_free (sam_hdr_contigs);
        FREE (sam_hdr_contigs);
    }    
}

static void sam_header_match_contig_to_reference ()
{
/*    // case match_chrom_to_reference: update contig_name to the matching name in the reference
    // Note: we match reference contigs by name only, not searching LN, we just verify the LN after the fact.
    // This is so the process is the same when matching CHROM fields while segging, which will result in the header and lines
    // behaving the same
    WordIndex ref_index = ref_contigs_get_by_name (gref , STRa(contig_name), true, true);
    if (ref_index != WORD_INDEX_NONE) {
        contig_name = ref_contigs_get_name (gref, ref_index, &contig_name_len); // this contig as its called in the reference
        length = ref_contigs_get_contig_length (gref, ref_index, 0, 0, true);   // update - we check later that it is consistent
    }

    bufprintf (evb, new_txt_header, "%.*s<ID=%.*s,length=%"PRIu64">\n", key_len, line, STRf(contig_name), length);
    printed = true;
*/}

static ContigP sam_header_add_contig (STRp (hdr_contig_name), PosType LN, void *unused)
{
    if (flag.match_chrom_to_reference)
        sam_header_match_contig_to_reference();

    // case: we have a reference, we use the reference chrom_index
    WordIndex ref_index;
    if (flag.reference & REF_ZIP_LOADED)
        ref_index = ref_contigs_ref_chrom_from_header_chrom (gref, STRa(hdr_contig_name), &LN);

    // case: no reference, chrom_index is the index of the CHROM context, where these contigs will be copied by zip_prepopulate_contig_ctxs
    // note: if is_luft_contig, then we're compressing a file that's already DVCF, so no --chain and we don't have a second reference
    else
        ref_index = WORD_INDEX_NONE; // this_contigs->len;

    // add to contigs. note: index is sam_hdr_contigs is by order of appearance in header, not the same as the reference
    NEXTENT (Contig, sam_hdr_contigs->contigs) = (Contig){ 
        .max_pos    = LN,        // if last_pos was provided, it is unchanged, and verified = the pos in ref. if not provided, it is what the reference says.
        .char_index = sam_hdr_contigs->dict.len, 
        .snip_len   = hdr_contig_name_len,
        .ref_index  = ref_index  // matching reference contig (with possibly a different variation of the contig name)
    };

    if (flag.show_txt_contigs) 
        iprintf ("index=%u \"%.*s\" LN=%"PRId64" ref_index=%d snip_len=%u\n", 
                 (unsigned)sam_hdr_contigs->contigs.len-1, STRf(hdr_contig_name), LN, ref_index, hdr_contig_name_len);

    // add to contigs_dict
    buf_add (&sam_hdr_contigs->dict, hdr_contig_name, hdr_contig_name_len);
    NEXTENT (char, sam_hdr_contigs->dict) = 0; // nul-termiante

    return LASTENT (Contig, sam_hdr_contigs->contigs);
}

// verify that this contig in this 2nd+ component, also exists (with the same LN) in the contig list created in the 1st component
// Note: we don't allow components to add new contigs, as currently BAM expects header contigs to be first in the RNAME, and there might already
// be additional contigs after (eg reference contigs added by zip_prepopulate_contig_ctxs and/or "*"). To do: fix this
static ContigP sam_header_verify_contig (STRp (hdr_contig), PosType hdr_contig_LN, void *unused)
{
    ASSINP (sam_hdr_contigs->cp_next < sam_hdr_contigs->cp_num_contigs, 
            "Error: txt header: contigs mismatch between files: first file has %u contigs, but %s has more",
             (unsigned)sam_hdr_contigs->cp_num_contigs, txt_name);

    if (flag.match_chrom_to_reference)
        sam_header_match_contig_to_reference();

    STR(sam_hdr_contig_name);
    sam_hdr_contig_name = contigs_get_name (sam_hdr_contigs, sam_hdr_contigs->cp_next, &sam_hdr_contig_name_len);
    PosType sam_hdr_contig_LN   = contigs_get_LN (sam_hdr_contigs, sam_hdr_contigs->cp_next);
    
    ASSINP (str_issame (hdr_contig, sam_hdr_contig_name),
            "Error: %s header contig=%u: contig name mismatch between files: in first file: \"%s\", in %s: \"%.*s\"",
            dt_name(txt_file->data_type), (unsigned)sam_hdr_contigs->cp_next, sam_hdr_contig_name, txt_name, STRf(hdr_contig));

    ASSINP (hdr_contig_LN == sam_hdr_contig_LN, "Error: %s header in \"%s\": contig length mismatch between files: in first file: LN:%"PRId64", in %s: LN:%"PRId64,
            dt_name(txt_file->data_type), sam_hdr_contig_name, sam_hdr_contig_LN, txt_name, hdr_contig_LN);

    return contigs_get_by_index (sam_hdr_contigs, sam_hdr_contigs->cp_next++);
}

// constructs hdr_contigs, and in ZIP, also initialzes refererence ranges and random_access
bool sam_header_inspect (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags)
{    
    if (flag.show_txt_contigs && exe_type == EXE_GENOCAT) exit_ok();

    if (command != ZIP) return true; // nothing to inspect in in PIZ, all good

    if (!IS_BAM) *AFTERENT (char, *txt_header) = 0; // nul-terminate as required by sam_foreach_SQ_line

    // if there is no external reference provided, then we create our internal one, and store it
    // (if an external reference IS provided, the user can decide whether to store it or not, with --reference vs --REFERENCE)
    if (flag.reference == REF_NONE) flag.reference = REF_INTERNAL;

    // in case of internal reference, we need to initialize. in case of --reference, it was initialized by ref_load_external_reference()
    if (!flag.reference || flag.reference == REF_INTERNAL) 
        ref_initialize_ranges (gref, RT_DENOVO); // it will be REF_INTERNAL if this is the 2nd+ non-conatenated file

    // evb buffers must be alloced by main threads, since other threads cannot modify evb's buf_list
    random_access_alloc_ra_buf (evb, DC_PRIMARY, 0);

    // if its the 2nd+ file while binding in ZIP - we just check that the contigs are the same
    // (or a subset) of previous file's contigs (to do: merge headers, bug 327) 
    if (sam_hdr_contigs && flag.bind && z_file->num_txt_components_so_far /* 2nd+ file */) {
        sam_hdr_contigs->cp_next = 0; 
        foreach_SQ_line (sam_header_verify_contig, NULL);
    }
    
    // first bound file, or non-binding ; PIZ - 1st header when concatenating or all headers when unbinding 
    else {
        sam_header_alloc_contigs (txt_header); 
        if (sam_hdr_contigs)
            foreach_SQ_line (sam_header_add_contig, NULL);
    }

    return true;
}

// ZIP: called from txtfile_read_header
// returns header length if header read is complete + sets lines.len, -1 not complete yet 
// note: usually a BAM header fits into a single 512KB READ BUFFER, so this function is called only twice (without and then with data).
// callback from DataTypeProperties.is_header_done
int32_t bam_is_header_done (bool is_eof)
{
    #define HDRSKIP(n) if (evb->txt_data.len < next + n) goto incomplete_header; next += n
    #define HDR32 (next + 4 <= evb->txt_data.len ? GET_UINT32 (&evb->txt_data.data[next]) : 0) ; if (evb->txt_data.len < next + 4) goto incomplete_header; next += 4;

    uint32_t next=0;

    HDRSKIP(4); // magic
    ASSINP (!memcmp (evb->txt_data.data, BAM_MAGIC, 4), // magic
            "%s doesn't have a BAM magic - it doesn't seem to be a BAM file", txt_name);

    // sam header text
    uint32_t l_text = HDR32;
    const char *text = ENT (const char, evb->txt_data, next);
    HDRSKIP(l_text);

    uint32_t header_n_ref = HDR32;

    for (uint32_t ref_i=0; ref_i < header_n_ref; ref_i++) {
        uint32_t l_name = HDR32;
        HDRSKIP (l_name + 4); // skip name and l_ref
    }

    // we have the entire header - count the text lines in the SAM header
    evb->lines.len = 0;
    for (unsigned i=0; i < l_text; i++)
        if (text[i] == '\n') evb->lines.len++;

    return next; // return BAM header length

incomplete_header:
    return -1;

    #undef HDRSKIP
    #undef HDR32
}   

//------------------------------------------------------------------
// Translator functions for reconstructing SAM<-->BAM header formats
//------------------------------------------------------------------

// add @PG line
static inline void sam_header_add_PG (Buffer *txtheader_buf)
{
    if (flag.no_pg) return; // user specified --no-PG

    TimeSpecType tb;
    clock_gettime (CLOCK_REALTIME, &tb);

    uint32_t unique_id = libdeflate_adler32 (1, &tb, sizeof (tb));

    // the command line length is unbound, careful not to put it in a bufprintf
    bufprintf (txtheader_buf->vb, txtheader_buf, "@PG\tID:genozip-%u\tPN:genozip\tDS:%s\tVN:%s\tCL:", 
               unique_id, GENOZIP_URL, GENOZIP_CODE_VERSION);
    buf_add_string (txtheader_buf->vb, txtheader_buf, flags_command_line());
    buf_add_string (txtheader_buf->vb, txtheader_buf, "\n");
}

// PIZ main thread: make the txt header either SAM or BAM according to flag.out_dt, and regardless of the source file
TXTHEADER_TRANSLATOR (sam_header_bam2sam)
{
    if (flag.no_header) {
        txtheader_buf->len = 0; // remove the BAM header data
        return;
    }

    ASSERT0 (buf_is_alloc (txtheader_buf), "txtheader_buf not allocated");

    uint32_t l_text = GET_UINT32 (ENT (char, *txtheader_buf, 4));
    memmove (txtheader_buf->data, ENT (char, *txtheader_buf, 8), l_text);
    txtheader_buf->len = l_text;
    
    sam_header_add_PG (txtheader_buf);
}

static ContigP sam_header_sam2bam_count_sq (const char *chrom_name, unsigned chrom_name_len, PosType last_pos, void *callback_param)
{
    (*(uint32_t *)callback_param)++;
    return NULL;
}

static ContigP sam_header_sam2bam_ref_info (STRp (ref_contig_name), PosType last_pos, void *callback_param)
{
    Buffer *txtheader_buf = (Buffer *)callback_param;

    buf_alloc (txtheader_buf->vb, txtheader_buf, ref_contig_name_len+9, 0, char, 1, 0);

    // l_name
    ref_contig_name_len++; // inc. nul terminator
    *(uint32_t *)AFTERENT (char, *txtheader_buf) = LTEN32 (ref_contig_name_len); 
    txtheader_buf->len += sizeof (uint32_t);

    ASSINP (ref_contig_name_len <= INT32_MAX, "Error: cannot convert to BAM because l_name=%u exceeds BAM format maximum of %u", ref_contig_name_len, INT32_MAX);

    // name
    buf_add (txtheader_buf, ref_contig_name, ref_contig_name_len-1);                  
    NEXTENT (char, *txtheader_buf) = 0;

    // l_ref
    uint32_t last_pos32 = (uint32_t)last_pos;
    *(uint32_t *)AFTERENT (char, *txtheader_buf) = LTEN32 (last_pos32); 
    txtheader_buf->len += sizeof (uint32_t);

    ASSINP (last_pos <= INT32_MAX, "Error: cannot convert to BAM because contig %.*s has length=%"PRId64" that exceeds BAM format maximum of %u", 
            ref_contig_name_len-1, ref_contig_name, last_pos, INT32_MAX);

    return NULL;
}

// prepare BAM header from SAM header, according to https://samtools.github.io/hts-specs/SAMv1.pdf section 4.2
TXTHEADER_TRANSLATOR (sam_header_sam2bam)
{
    // grow buffer to accommodate the BAM header fixed size and text (inc. nul terminator)
    buf_alloc (comp_vb, txtheader_buf, 12 + 1, 0, char, 1, "txt_data");

    // nul-terminate text - required by sam_foreach_SQ_line - but without enlengthening buffer
    *AFTERENT (char, *txtheader_buf) = 0;

    // n_ref = count SQ lines in text
    uint32_t n_ref=0;
    sam_foreach_SQ_line (txtheader_buf->data, sam_header_sam2bam_count_sq, &n_ref);

    // we can't convert to BAM if its a SAM file with aligned reads but without SQ records, compressed with REF_INTERNAL - as using the REF_INTERNAL
    // contigs would produce lengths that don't match actual reference files - rendering the BAM file useless for downstream
    // analysis. Better give an error here than create confusion downstream. 
    ASSINP (n_ref ||  // has SQ records
            !z_file->z_flags.dts_ref_internal || // has external reference
            ZCTX(SAM_POS)->word_list.len==1, // has only one POS word = "Delta 0" = unaligned SAM that doesn't need contigs 
            "Failed to convert %s from SAM to BAM: genounzip requires that either the SAM header has SQ records (see https://samtools.github.io/hts-specs/SAMv1.pdf section 1.3), or the file was genozipped with --reference or --REFERENCE", z_name);

    // if no SQ lines - get lines from loaded contig (will be available only if file was compressed with --reference or --REFERENCE)
    bool from_SQ = !!n_ref;
    if (!from_SQ) n_ref = ref_num_contigs (gref);

    // grow buffer to accommodate estimated reference size (we will more in sam_header_sam2bam_ref_info if not enough)
    buf_alloc (comp_vb, txtheader_buf, n_ref * 100 + (50 + strlen (flags_command_line())), 0, char, 1, 0);

    // add PG
    sam_header_add_PG (txtheader_buf);

    // construct magic, l_text, text and n_ref fields of BAM header
    char *text = txtheader_buf->data + 8;
    uint32_t l_text = txtheader_buf->len;
    ASSINP (txtheader_buf->len <= INT32_MAX, "Cannot convert to BAM because SAM header length (%"PRIu64" bytes) exceeds BAM format maximum of %u", txtheader_buf->len, INT32_MAX);

    memmove (text, txtheader_buf->data, l_text);                          // text
    memcpy (FIRSTENT (char, *txtheader_buf), BAM_MAGIC, 4);               // magic
    *ENT (uint32_t, *txtheader_buf, 1) = LTEN32 (l_text);                 // l_text
    *(uint32_t *)ENT (char, *txtheader_buf, l_text + 8) = LTEN32 (n_ref); // n_ref
    txtheader_buf->len += 12; // BAM header length so far, before adding reference information

    // option 1: copy reference information from SAM SQ lines, if available
    if (from_SQ) 
        sam_foreach_SQ_line (text, sam_header_sam2bam_ref_info, txtheader_buf);

    // option 2: copy reference information from ref_contigs if available - i.e. if pizzed with external or stored reference
    else if (ref_num_contigs (gref)) 
        foreach_contig (ref_get_ctgs (gref), sam_header_sam2bam_ref_info, txtheader_buf);
}
