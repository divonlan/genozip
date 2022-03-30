// ------------------------------------------------------------------
//   sam_header.c
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
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

//---------------
// Header contigs
//---------------

ContigPkgP sam_hdr_contigs = NULL; // If no contigs in header: BAM: empty structure (RNAME must be * for all lines); SAM: NULL (RNAME in lines may have contigs) 

#define foreach_SQ_line(cb,cb_param,new_txt_header) (IS_BAM ? bam_foreach_SQ_line (txt_header->data, txt_header->len, cb, cb_param, new_txt_header) \
                                                            : sam_foreach_SQ_line (txt_header->data, cb, cb_param, new_txt_header))

static void sam_header_get_ref_index (STRp (contig_name), PosType LN, void *ref_index)
{
    *(WordIndex *)ref_index = ref_contigs_ref_chrom_from_header_chrom (gref, STRa(contig_name), &LN); // also verifies LN
}

static void sam_header_add_contig (STRp (contig_name), PosType LN, void *out_ref_index)
{
    // case: we have a reference, we use the reference chrom_index
    WordIndex ref_index = WORD_INDEX_NONE;
    if (flag.reference & REF_ZIP_LOADED) {
        ref_index = ref_contigs_ref_chrom_from_header_chrom (gref, STRa(contig_name), &LN); // also verifies LN

        // case --match-chrom-to-reference
        if (flag.match_chrom_to_reference) {
            if (ref_index != WORD_INDEX_NONE) // udpate contig name, if this contig is in the reference
                contig_name = ref_contigs_get_name (gref, ref_index, &contig_name_len);

            *(WordIndex*)out_ref_index = ref_index;
        }
    }

    // In case of REF_INTERNAL compression, and we have a header contigs - we pre-allocate the range space
    PosType gpos = ref_samheader_denovo_get_header_contig_gpos (sam_hdr_contigs->contigs.len ? BLST (Contig, sam_hdr_contigs->contigs) : NULL);

    // add to contigs. note: index is sam_hdr_contigs is by order of appearance in header, not the same as the reference
    BNXT (Contig, sam_hdr_contigs->contigs) = (Contig){ 
        .max_pos    = LN,        // if last_pos was provided, it is unchanged, and verified = the pos in ref. if not provided, it is what the reference says.
        .char_index = sam_hdr_contigs->dict.len, 
        .snip_len   = contig_name_len,
        .ref_index  = ref_index, // matching reference contig (with possibly a different variation of the contig name)
        .gpos       = gpos
    };

    if (flag.show_txt_contigs) 
        iprintf ("index=%u \"%.*s\" LN=%"PRId64" ref_index=%d snip_len=%u gpos=%"PRId64"\n", 
                 (unsigned)sam_hdr_contigs->contigs.len-1, STRf(contig_name), LN, ref_index, contig_name_len, gpos);

    // add to contigs_dict
    buf_add_more (NULL, &sam_hdr_contigs->dict, contig_name, contig_name_len, NULL);
    BNXTc (sam_hdr_contigs->dict) = 0; // nul-termiante
}

// verify that this contig in this 2nd+ component, also exists (with the same LN) in the contig list created in the 1st component
// Note: we don't allow components to add new contigs, as currently BAM expects header contigs to be first in the RNAME, and there might already
// be additional contigs after (eg reference contigs added by ctx_populate_zf_ctx_from_contigs and/or "*"). To do: fix this
static void sam_header_verify_contig (STRp (hdr_contig), PosType hdr_contig_LN, void *unused)
{
    ASSINP (sam_hdr_contigs->cp_next < sam_hdr_contigs->cp_num_contigs, 
            "Error: txt header: contigs mismatch between files: first file has %u contigs, but %s has more",
             (unsigned)sam_hdr_contigs->cp_num_contigs, txt_name);

    STR(sam_hdr_contig_name);
    sam_hdr_contig_name = contigs_get_name (sam_hdr_contigs, sam_hdr_contigs->cp_next, &sam_hdr_contig_name_len);
    PosType sam_hdr_contig_LN = contigs_get_LN (sam_hdr_contigs, sam_hdr_contigs->cp_next);
    
    ASSINP (str_issame (hdr_contig, sam_hdr_contig_name),
            "Error: %s header contig=%u: contig name mismatch between files: in first file: \"%s\", in %s: \"%.*s\"",
            dt_name(txt_file->data_type), (unsigned)sam_hdr_contigs->cp_next, sam_hdr_contig_name, txt_name, STRf(hdr_contig));

    ASSINP (hdr_contig_LN == sam_hdr_contig_LN, "Error: %s header in \"%s\": contig length mismatch between files: in first file: LN:%"PRId64", in %s: LN:%"PRId64,
            dt_name(txt_file->data_type), sam_hdr_contig_name, sam_hdr_contig_LN, txt_name, hdr_contig_LN);

    contigs_get_by_index (sam_hdr_contigs, sam_hdr_contigs->cp_next++);
}

// call a callback for each SQ line (contig). Note: callback function is the same as foreach_contig
static void sam_foreach_SQ_line (rom txt_header, // nul-terminated string
                                 ContigsIteratorCallback callback, void *callback_param,
                                 Buffer *new_txt_header) // optional - if provided, all lines need to be copied (modified or not) to this buffer
{
    rom line = txt_header;

    while (1) {
        line = strchr (line, '@');
        if (!line) break;

        char *newline = strchr (line, '\n');
        *newline = 0;

        bool add_line = !!new_txt_header; // whether we need to add the line to new_txt_header

        if (line[1] == 'S' && line[2] == 'Q') { // this test will always be in the buffer - possible in the overflow area

            // line looks something like: @SQ\tSN:chr10\tLN:135534747 (no strict order between SN and LN and there could be more parameters)
            rom chrom_name = strstr (line, "SN:"); 
            rom pos_str    = strstr (line, "LN:"); 
            if (chrom_name && pos_str) {
                unsigned chrom_name_len = strcspn (&chrom_name[3], "\t\n\r");
                rom after_chrom = chrom_name + chrom_name_len + STRLEN("SN:");

                PosType last_pos = (PosType)strtoull (&pos_str[3], NULL, 10);

                ASSINP (last_pos <= MAX_POS_SAM, "Error: @SQ record in header contains LN:%"PRId64" which is beyond the maximum permitted by the SAM spec of %"PRId64,
                        last_pos, MAX_POS_SAM);

                WordIndex ref_index = WORD_INDEX_NONE; // of reference contig with matched name
                callback (&chrom_name[3], chrom_name_len, last_pos, add_line ? &ref_index : callback_param);

                // case --match-chrom-to-reference: we need to update the contig name from the reference, if this contig is in the reference
                if (add_line && ref_index != WORD_INDEX_NONE) {
                    add_line = false; // done
                    bufprintf (evb, new_txt_header, "%.*s%s%.*s\n",
                               (int)(chrom_name - line + STRLEN("SN:")), line,      // line beginning up to the contig name
                               ref_contigs_get_name (gref, ref_index, NULL), // matched contig name from the reference
                               (int)(newline - after_chrom), after_chrom);      // line after contig name
                }
            }
        }

        *newline = '\n'; // restore

        if (add_line)
            buf_add_more (NULL, new_txt_header, line, newline-line+1, NULL);

        line = newline+1;
    }
}

// call a callback for each SQ line (contig). Note: callback function is the same as foreach_contig
static void bam_foreach_SQ_line (rom txt_header, uint32_t txt_header_len, // binary BAM header
                                 ContigsIteratorCallback callback, void *callback_param,
                                 Buffer *new_txt_header) // optional - if provided, all lines need to be copied (modified or not) to this buffer
{
    #define HDRSKIP(n) if (txt_header_len < next + n) goto incomplete_header; next += n
    #define HDR32 (next + 4 <= txt_header_len ? GET_UINT32 (&txt_header[next]) : 0) ; if (txt_header_len < next + 4) goto incomplete_header; next += 4;

    uint32_t next=0;

    // skip magic, l_text, text
    HDRSKIP(4);
    uint32_t l_text = HDR32;
    rom plain_sam_header = &txt_header[next];
    HDRSKIP(l_text);

    uint32_t n_ref = HDR32;

    // case --match-chrom-to-reference: 
    if (new_txt_header) {
        // copy magic and l_text - we will update l_text later
        buf_add (new_txt_header, txt_header, 8); 
        
        // update SAM plain header with new contigs, but DONT add to hdr_contigs - we will do that in the loop below
        SAFE_NUL (&plain_sam_header[l_text]); // per SAM spec, plain SAM header within a BAM header is not nul terminated
        sam_foreach_SQ_line (plain_sam_header, sam_header_get_ref_index, NULL, new_txt_header);
        SAFE_RESTORE;

        // update l_text
        *B32 (*new_txt_header, 1) = LTEN32 ((uint32_t)new_txt_header->len - 8); // this assignment is word-aligned

        // add n_ref
        uint32_t n_ref_lten = LTEN32 (n_ref);
        buf_add_more (NULL, new_txt_header, (char*)&n_ref_lten, sizeof (uint32_t), NULL); // not necessarily word-aligned
    }

    for (uint32_t ref_i=0; ref_i < n_ref; ref_i++) {
        
        uint32_t start_line = next;

        uint32_t l_name = HDR32; // name length including \0
        rom name = &txt_header[next];
        HDRSKIP (l_name);
        uint32_t l_ref = HDR32;

        WordIndex ref_index = WORD_INDEX_NONE; // of reference contig with matched name
        callback (name, l_name-1, l_ref, new_txt_header ? &ref_index : callback_param);

        // case --match-chrom-to-reference: we need to update the contig name from the reference, if this contig is in the reference
        if (new_txt_header) {
            if (ref_index != WORD_INDEX_NONE) {            
                uint32_t new_l_name;
                name = ref_contigs_get_name (gref, ref_index, &new_l_name); 
                new_l_name++; // now its the length including \0

                uint32_t new_l_name_lten = LTEN32 (new_l_name);
                buf_add_more (NULL, new_txt_header, (char*)&new_l_name_lten, sizeof (uint32_t), NULL); // use memcpy as data is not word-aligned
                buf_add_more (NULL, new_txt_header, name, new_l_name_lten, NULL);

                uint32_t l_ref_lten = LTEN32 (l_ref);
                buf_add_more (NULL, new_txt_header, (char*)&l_ref_lten, sizeof (uint32_t), NULL); // use memcpy as data is not word-aligned
            }
            else 
                buf_add_more (NULL, new_txt_header, &txt_header[start_line], l_name + 8, NULL); // use memcpy as data is not word-aligned
        }
    }

    return;

incomplete_header:
    ABORT ("incomplete BAM header (next=%u)", next);

    #undef HDRSKIP
    #undef HDR32
}   

typedef struct { uint32_t num_contigs, dict_len; } ContigsCbParam;
static void sam_header_count_contigs_cb (rom chrom_name, unsigned chrom_name_len, PosType last_pos, void *cb_param)
{
    // note: for simplicity, we consider contigs to be in the range [0, last_pos] even though POS=0 doesn't exist - last_pos+1 loci
    ((ContigsCbParam*)cb_param)->num_contigs++;
    ((ContigsCbParam*)cb_param)->dict_len += chrom_name_len + 1;
}

// ZIP
static void sam_header_alloc_contigs (BufferP txt_header)
{
    ContigsCbParam contigs_count = {}; 
    foreach_SQ_line (sam_header_count_contigs_cb, &contigs_count, NULL);

    // note: In BAM, a contig-less header is considered an unaligned BAM in which all RNAME=*
    // In SAM, a contig-less header just means we don't know the contigs and they are defined in RNAME 
    if (!IS_BAM_ZIP && !contigs_count.num_contigs) return; // "headerless" SAM (not really headerless, just no contigs in the header)

    sam_hdr_contigs = CALLOC (sizeof (ContigPkg));
    sam_hdr_contigs->name = "sam_hdr_contigs";
    
    buf_alloc (evb, &sam_hdr_contigs->contigs, 0, contigs_count.num_contigs, Contig, 1, "JJContigPkg->contigs"); 
    buf_alloc (evb, &sam_hdr_contigs->dict,    0, contigs_count.dict_len, char, 1, "JJContigPkg->dict"); 
}

void sam_header_finalize (void)
{
    contigs_free (sam_hdr_contigs);
    FREE (sam_hdr_contigs);
}

// ZIP
bool sam_header_has_string (ConstBufferP txt_header, rom substr)
{
    uint32_t next=0;
    #define HDRSKIP(n) if (txt_header->len < next + n) return false; next += n
    #define HDR32 (next + 4 <= txt_header->len ? GET_UINT32 (&txt_header->data[next]) : 0) ; if (txt_header->len < next + 4) return false; next += 4;

    if (IS_BAM_ZIP) {
        HDRSKIP(4); // magic
        uint32_t l_text = HDR32;
        rom text = Bc (*txt_header, next);
        SAFE_NUL (&text[l_text]);
        bool found = !!strstr (text, substr);
        SAFE_RESTORE;
        return found;
    }
    else  // SAM - already nul-terminated in sam_header_inspect
        return !!strstr (txt_header->data, substr);

    #undef HDRSKIP
    #undef HDR32
}

// constructs hdr_contigs, and in ZIP, also initialzes refererence ranges and random_access
bool sam_header_inspect (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags)
{    
    if (flag.show_txt_contigs && exe_type == EXE_GENOCAT) exit_ok();

    if (command != ZIP) return true; // nothing to inspect in in PIZ or SA components, all good

    if (!IS_BAM_ZIP) *BAFTc (*txt_header) = 0; // nul-terminate as required by sam_foreach_SQ_line

    segconf.has_bsseeker2 = sam_header_has_string (txt_header, "PN:BS Seeker 2"); 
    segconf.sam_bowtie2   = sam_header_has_string (txt_header, "PN:bowtie2");

    // if there is no external reference provided, then we create our internal one, and store it
    // (if an external reference IS provided, the user can decide whether to store it or not, with --reference vs --REFERENCE)
    if (flag.reference == REF_NONE) flag.reference = REF_INTERNAL;

    // in case of internal reference, we need to initialize. in case of --reference, it was initialized by ref_load_external_reference()
    if (!flag.reference || flag.reference == REF_INTERNAL) 
        ref_initialize_ranges (gref, RT_DENOVO); // it will be REF_INTERNAL if this is the 2nd+ non-conatenated file

    // evb buffers must be alloced by the main thread, since other threads cannot modify evb's buf_list
    random_access_alloc_ra_buf (evb, 0, 0);

    // if its the 2nd+ file while binding in ZIP - we just check that the contigs are the same
    // (or a subset) of previous file's contigs (to do: merge headers, bug 327) 
    if (sam_hdr_contigs && flag.bind && z_file->num_txts_so_far /* 2nd+ file */) {
        sam_hdr_contigs->cp_next = 0; 
        foreach_SQ_line (sam_header_verify_contig, NULL, NULL);
    }
    
    // first bound file, or non-binding ; PIZ - 1st header when concatenating or all headers when unbinding 
    else {
        sam_header_alloc_contigs (txt_header); 
        if (sam_hdr_contigs) {

            if (flag.match_chrom_to_reference) {
                #define new_txt_header txt_header_vb->codec_bufs[0]
                buf_alloc (txt_header_vb, &new_txt_header, 0, txt_header->len + 100, char, 0, "codec_bufs[0]"); // initial allocation (might be a bit bigger due to label changes)

                foreach_SQ_line (sam_header_add_contig, NULL, &new_txt_header);

                // replace txt_header with lifted back one
                buf_free (*txt_header);
                buf_copy (txt_header_vb, txt_header, &new_txt_header, char, 0, 0, "txt_data");
                buf_free (new_txt_header);
                #undef new_txt_header
            }
            else
                foreach_SQ_line (sam_header_add_contig, NULL, NULL);

            contigs_create_index (sam_hdr_contigs, SORT_BY_NAME); // used by sam_sa_add_sa_group
        }
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
            "%s doesn't have a BAM magic - it doesn't seem to be a BAM file: magic=%2.2x%2.2x%2.2x%2.2x ", txt_name, (uint8_t)evb->txt_data.data[0], (uint8_t)evb->txt_data.data[1], (uint8_t)evb->txt_data.data[2], (uint8_t)evb->txt_data.data[3]);

    // sam header text
    uint32_t l_text = HDR32;
    rom text = Bc (evb->txt_data, next);
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

    uint32_t unique_id = adler32 (1, &tb, sizeof (tb));

    // the command line length is unbound, careful not to put it in a bufprintf
    bufprintf (txtheader_buf->vb, txtheader_buf, "@PG\tID:genozip-%u\tPN:genozip\tDS:%s\tVN:%s\tCL:", 
               unique_id, GENOZIP_URL, GENOZIP_CODE_VERSION);
    buf_add_string (txtheader_buf->vb, txtheader_buf, flags_command_line());
    buf_add_string (txtheader_buf->vb, txtheader_buf, "\n");
}

// PIZ main thread: make the txt header either SAM or BAM according to flag.out_dt, and regardless of the source file
TXTHEADER_TRANSLATOR (sam_header_bam2sam)
{
    if (flag.no_header || !txtheader_buf->len /* SA component */) {
        txtheader_buf->len = 0; // remove the BAM header data
        return;
    }

    ASSERT0 (buf_is_alloc (txtheader_buf), "txtheader_buf not allocated");

    uint32_t l_text = GET_UINT32 (Bc (*txtheader_buf, 4));
    memmove (txtheader_buf->data, Bc (*txtheader_buf, 8), l_text);
    txtheader_buf->len = l_text;
    
    sam_header_add_PG (txtheader_buf);
}

static void sam_header_sam2bam_count_sq (rom chrom_name, unsigned chrom_name_len, PosType last_pos, void *callback_param)
{
    (*(uint32_t *)callback_param)++;
}

static void sam_header_sam2bam_ref_info (STRp (ref_contig_name), PosType last_pos, void *callback_param)
{
    Buffer *txtheader_buf = (Buffer *)callback_param;

    buf_alloc (txtheader_buf->vb, txtheader_buf, ref_contig_name_len+9, 0, char, 1, 0);

    // l_name
    ref_contig_name_len++; // inc. nul terminator
    *(uint32_t *)BAFTc (*txtheader_buf) = LTEN32 (ref_contig_name_len); 
    txtheader_buf->len += sizeof (uint32_t);

    ASSINP (ref_contig_name_len <= INT32_MAX, "Error: cannot convert to BAM because l_name=%u exceeds BAM format maximum of %u", ref_contig_name_len, INT32_MAX);

    // name
    buf_add (txtheader_buf, ref_contig_name, ref_contig_name_len-1);                  
    BNXTc (*txtheader_buf) = 0;

    // l_ref
    uint32_t last_pos32 = (uint32_t)last_pos;
    *(uint32_t *)BAFTc (*txtheader_buf) = LTEN32 (last_pos32); 
    txtheader_buf->len += sizeof (uint32_t);

    ASSINP (last_pos <= INT32_MAX, "Error: cannot convert to BAM because contig %.*s has length=%"PRId64" that exceeds BAM format maximum of %u", 
            ref_contig_name_len-1, ref_contig_name, last_pos, INT32_MAX);
}

// prepare BAM header from SAM header, according to https://samtools.github.io/hts-specs/SAMv1.pdf section 4.2
TXTHEADER_TRANSLATOR (sam_header_sam2bam)
{
    // grow buffer to accommodate the BAM header fixed size and text (inc. nul terminator)
    buf_alloc (comp_vb, txtheader_buf, 12 + 1, 0, char, 1, "txt_data");

    // nul-terminate text - required by sam_foreach_SQ_line - but without enlengthening buffer
    *BAFTc (*txtheader_buf) = 0;

    // n_ref = count SQ lines in text
    uint32_t n_ref=0;
    sam_foreach_SQ_line (txtheader_buf->data, sam_header_sam2bam_count_sq, &n_ref, NULL);

    // we can't convert to BAM if its a SAM file with aligned reads but without SQ records, compressed with REF_INTERNAL - as using the REF_INTERNAL
    // contigs would produce lengths that don't match actual reference files - rendering the BAM file useless for downstream
    // analysis. Better give an error here than create confusion downstream. 
    ASSINP (n_ref ||  // has SQ records
            !z_file->z_flags.dts_ref_internal || // has external reference
            ZCTX(SAM_POS)->word_list.len==1, // has only one POS word = "Delta 0" = unaligned SAM that doesn't need contigs 
            "Failed to convert %s from SAM to BAM: the SAM header must have SQ records (see https://samtools.github.io/hts-specs/SAMv1.pdf section 1.3), or alternatively the file should be genozipped with --reference or --REFERENCE", z_name);

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
    memcpy (B1STc (*txtheader_buf), BAM_MAGIC, 4);               // magic
    *B32 (*txtheader_buf, 1) = LTEN32 (l_text);                 // l_text
    *(uint32_t *)Bc (*txtheader_buf, l_text + 8) = LTEN32 (n_ref); // n_ref
    txtheader_buf->len += 12; // BAM header length so far, before adding reference information

    // option 1: copy reference information from SAM SQ lines, if available
    if (from_SQ) 
        sam_foreach_SQ_line (text, sam_header_sam2bam_ref_info, txtheader_buf, NULL);

    // option 2: copy reference information from ref_contigs if available - i.e. if pizzed with external or stored reference
    else if (ref_num_contigs (gref)) 
        foreach_contig (ref_get_ctgs (gref), sam_header_sam2bam_ref_info, txtheader_buf);
}
