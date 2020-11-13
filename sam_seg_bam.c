// ------------------------------------------------------------------
//   sam_bam.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "md5.h"
#include "buffer.h"
#include "vblock.h"
#include "txtfile.h"
#include "file.h"
#include "endianness.h"
#include "sam_private.h"
#include "seg.h"
#include "strings.h"
#include "random_access.h"
#include "dict_id.h"
#include "codec.h"

// copy header ref data during reading header
static uint32_t header_n_ref=0;
static char **header_ref_contigs = NULL;
static uint32_t *header_ref_contigs_len = NULL;

// returns header length if header read is complete + sets lines.len, -1 not complete yet 
// note: usually a BAM header fits into a single 512KB READ BUFFER, so this function is called only twice (without and then with data).
// callback from DataTypeProperties.is_header_done
int32_t bam_is_header_done (void)
{
    uint32_t next=0;
    #define HDRLEN evb->txt_data.len
    #define HDRSKIP(n) if (HDRLEN < next + n) return -1; next += n
    #define HDR32 (next + 4 <= HDRLEN ? GET_UINT32 (&evb->txt_data.data[next]) : 0) ; if (HDRLEN < next + 4) return -1; next += 4;

    HDRSKIP(4); // magic
    ASSERT (!memcmp (evb->txt_data.data, BAM_MAGIC, 4), // magic
            "Error in bam_read_txt_header: %s doesn't have a BAM magic - it doesn't seem to be a BAM file", txt_name);

    // sam header text
    uint32_t l_text = HDR32;
    const char *text = ENT (const char, evb->txt_data, next);
    HDRSKIP(l_text);

    header_n_ref = HDR32;

    if (!header_ref_contigs) {
        header_ref_contigs     = MALLOC (header_n_ref * sizeof (char *));
        header_ref_contigs_len = MALLOC (header_n_ref * sizeof (uint32_t));
    }

    for (uint32_t ref_i=0; ref_i < header_n_ref; ref_i++) {
        header_ref_contigs_len[ref_i] = HDR32;
        header_ref_contigs[ref_i] = ENT (char, evb->txt_data, next);

        HDRSKIP(header_ref_contigs_len[ref_i]);
        header_ref_contigs_len[ref_i]--;

        HDRSKIP(4); // l_ref
    }

    // we have the entire header - count the text lines in the SAM header
    evb->lines.len = 0;
    for (unsigned i=0; i < l_text; i++)
        if (text[i] == '\n') evb->lines.len++;

    return next; // return BAM header length
}   

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
// if first_i > 0, we attempt to heuristically detect the start of a BAM alignment.
int32_t bam_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i)
{
    *i = MIN (*i, vb->txt_data.len - sizeof(BAMAlignmentFixed));

    // find the first alignment in the data (going backwards) that is entirely in the data - 
    // we identify and alignment by l_read_name and read_name
    for (; *i >= first_i; (*i)--) {
        BAMAlignmentFixed *aln = (BAMAlignmentFixed *)ENT (char, vb->txt_data, *i);
        if (&aln->read_name[aln->l_read_name] <= AFTERENT (char, vb->txt_data) && // we have the entire propose read name
            aln->read_name[aln->l_read_name-1] == 0 && // nul-terminated
            *i + LTEN32 (aln->block_size) + 4 <= vb->txt_data.len && // proposed block fits in data (+4 bc block_size value excludes itself)
            str_is_in_range (aln->read_name, aln->l_read_name-1, '!', '~')) // all printable ascii (per SAM spec)
            // TODO - add aln->bin calculation to increase confidence
            
            return vb->txt_data.len - (*i + LTEN32 (aln->block_size) + 4); // everything after this alignment is "unconsumed"
    }

    return -1; // we can't find any alignment - need more data (lower first_i)
}

void bam_seg_initialize (VBlock *vb)
{
    sam_seg_initialize (vb);

    // copy contigs from header in vb=1 for RNAME and RNEXT dictionaries
    if (vb->vblock_i == 1) {
        for (unsigned i=0; i < header_n_ref; i++) {
            ctx_evaluate_snip_seg (vb, &vb->contexts[SAM_RNAME], header_ref_contigs[i], header_ref_contigs_len[i], NULL);
            ctx_evaluate_snip_seg (vb, &vb->contexts[SAM_RNEXT], header_ref_contigs[i], header_ref_contigs_len[i], NULL);
        }

        // keep the contigs in the order as in BAM header as word_index of these also
        // reference to ref_id in the reference list in the BAM header
        vb->contexts[SAM_RNAME].inst |= CTX_INST_NO_VB1_SORT; 
        vb->contexts[SAM_RNEXT].inst |= CTX_INST_NO_VB1_SORT;
    }
}

void bam_zip_finalize (void)
{
    FREE (header_ref_contigs);
    FREE (header_ref_contigs_len);
    header_n_ref = 0;
}

void bam_seg_bin (VBlockSAM *vb, uint16_t bin /* used only in bam */, uint16_t flag, PosType this_pos)
{
    bool segment_unmapped = (flag & 0x4);
    bool is_bam = IS_BAM;

    PosType last_pos = segment_unmapped ? this_pos : (this_pos + vb->ref_consumed - 1);
    uint16_t reg2bin = bam_reg2bin (this_pos-1, last_pos-1);

    if (!is_bam || (last_pos <= (PosType)0x7ffffff /* largest POS is SAM */ && reg2bin == bin))
        seg_by_did_i (vb, ((char []){ SNIP_SPECIAL, SAM_SPECIAL_BIN }), 2, SAM_BAM_BIN, is_bam ? sizeof (uint16_t) : 0)
    
    else {
        WARN ("FYI: bad bin value in vb=%u vb->line_i=%u: this_pos=%"PRId64" ref_consumed=%u flag=%u last_pos=%"PRId64": bin=%u but reg2bin=%u. No harm.",
              vb->vblock_i, vb->line_i, this_pos, vb->ref_consumed, flag, last_pos, bin, reg2bin);
        seg_integer (vb, SAM_BAM_BIN, bin, is_bam);
        vb->contexts[SAM_BAM_BIN].flags = CTX_FL_STORE_INT;
    }
}

static inline void bam_seg_ref_id (VBlockP vb, DidIType did_i, int32_t ref_id, int32_t compare_to_ref_i)
{
    ASSERT (ref_id >= -1 && ref_id < (int32_t)header_n_ref, "Error in bam_seg_ref_id: vb=%u: encountered ref_id=%d but header has only %u contigs",
            vb->vblock_i, ref_id, header_n_ref);

    // get snip and snip_len
    DECLARE_SNIP;
    if (ref_id >= 0) {
        if (ref_id == compare_to_ref_i) {
            snip = "=";
            snip_len = 1;
        }
        else {
            snip = header_ref_contigs[ref_id];
            snip_len = header_ref_contigs_len[ref_id];
        }
    }
    else { 
        snip = "*";
        snip_len = 1;
    }

    WordIndex chrom_index = seg_by_did_i (vb, snip, snip_len, did_i, sizeof (int32_t));
        
    if (did_i==CHROM) 
        random_access_update_chrom (vb, chrom_index, snip, snip_len);
}

static inline void bam_rewrite_cigar (VBlockSAM *vb, uint16_t n_cigar_op, const uint32_t *cigar)
{
    if (!n_cigar_op) {
        buf_alloc (vb, &vb->textual_cigar, 2, 2, "textual_cigar");
        NEXTENT (char, vb->textual_cigar) = '*';
        goto finish;
    }

    // calculate length
    unsigned len=0;
    for (uint16_t i=0; i < n_cigar_op; i++) {
        uint32_t op_len = LTEN32 (cigar[i]) >> 4;
        if      (op_len < 10)      len += 2; // 1 for the op, 1 for the number
        else if (op_len < 100)     len += 3; 
        else if (op_len < 1000)    len += 4; 
        else if (op_len < 10000)   len += 5; 
        else if (op_len < 100000)  len += 6; 
        else if (op_len < 1000000) len += 7;
        else ABORT ("op_len=%u too long in vb=%u: ", vb->vblock_i, op_len); 
    }

    buf_alloc (vb, &vb->textual_cigar, len + 1 /* for \0 */, 2, "textual_cigar");

    for (uint16_t i=0; i < n_cigar_op; i++) {
        uint32_t subcigar = LTEN32 (cigar[i]);
        uint32_t op_len = subcigar >> 4;

        vb->textual_cigar.len += str_int (op_len, AFTERENT (char, vb->textual_cigar));

        static const char op_table[16] = "MIDNSHP=X";
        NEXTENT (char, vb->textual_cigar) = op_table[subcigar & 0xf];
    }

finish:
    *AFTERENT (char, vb->textual_cigar) = 0; // nul terminate
    vb->last_cigar = FIRSTENT (char, vb->textual_cigar);
}

static void bam_seg_cigar_field (VBlockSAM *vb, ZipDataLineSAM *dl, uint32_t l_seq, uint32_t n_cigar_op)
{
    char cigar_snip[vb->textual_cigar.len + 20];
    cigar_snip[0] = SNIP_SPECIAL;
    cigar_snip[1] = SAM_SPECIAL_CIGAR;
    unsigned cigar_snip_len=2;

    // case: we have no sequence - we add a '-' to the CIGAR
    if (!l_seq) cigar_snip[cigar_snip_len++] = '-';

    // case: CIGAR is "*" - we get the dl->seq_len directly from SEQ or QUAL, and add the length to CIGAR eg "151*"
    if (!dl->seq_len) { // CIGAR is not available
        dl->seq_len = l_seq;
        cigar_snip_len += str_int (MAX (dl->seq_len, 1), &cigar_snip[cigar_snip_len]); // if seq_len=0, then we add "1" because we have the * in seq and qual (and consistency with sam_seg_cigar_field)
    }
    
    memcpy (&cigar_snip[cigar_snip_len], vb->textual_cigar.data, vb->textual_cigar.len);
    cigar_snip_len += vb->textual_cigar.len;

    seg_by_did_i (vb, cigar_snip, cigar_snip_len, SAM_CIGAR, 
                  n_cigar_op * sizeof (uint32_t) /* cigar */ + sizeof (uint16_t) /* n_cigar_op */ );
}

static inline void bam_rewrite_seq (VBlockSAM *vb, uint32_t l_seq, const char *next_field)
{
    buf_alloc (vb, &vb->textual_seq, l_seq+1 /* +1 for last half-byte */, 1.5, "textual_seq");

    if (!l_seq) {
        NEXTENT (char, vb->textual_seq) = '*';
        return;        
    }

    char *next = FIRSTENT (char, vb->textual_seq);

    for (uint32_t i=0; i < (l_seq+1) / 2; i++) {
        static const char base_codes[16] = "=ACMGRSVTWYHKDBN";

        *next++ = base_codes[*(uint8_t*)next_field >> 4];
        *next++ = base_codes[*(uint8_t*)next_field & 0xf];
        next_field++;
    }

    vb->textual_seq.len = l_seq;

    ASSERTW (!(l_seq % 2) || (*AFTERENT(char, vb->textual_seq)=='='), 
            "Warning in bam_rewrite_seq vb=%u: expecting the unused lower 4 bits of last seq byte in an odd-length seq_len=%u to be 0, but its not. This will cause an incorrect MD5",
             vb->vblock_i, l_seq);
}

// Rewrite the QUAL field - add +33 to Phred scores to make them ASCII
static inline bool bam_rewrite_qual (uint8_t *qual, uint32_t l_seq)
{
    if (qual[0] == 0xff) return false; // in case SEQ is present but QUAL is omitted, all qual is 0xFF

    for (uint32_t i=0; i < l_seq; i++)
        qual[i] += 33;

    return true;
}

static inline const char *bam_rewrite_one_optional_number (VBlockSAM *vb, const char *next_field, uint8_t type)
{
    static const char special_float[2] = { SNIP_SPECIAL, SAM_SPECIAL_FLOAT };

    switch (type) {
        case 'c': { uint8_t  n = NEXT_UINT8 ; vb->textual_opt.len += str_int ((int8_t )n, AFTERENT (char, vb->textual_opt)); break; }
        case 'C': { uint8_t  n = NEXT_UINT8 ; vb->textual_opt.len += str_int (         n, AFTERENT (char, vb->textual_opt)); break; }
        case 's': { uint16_t n = NEXT_UINT16; vb->textual_opt.len += str_int ((int16_t)n, AFTERENT (char, vb->textual_opt)); break; }
        case 'S': { uint16_t n = NEXT_UINT16; vb->textual_opt.len += str_int (         n, AFTERENT (char, vb->textual_opt)); break; }
        case 'i': { uint32_t n = NEXT_UINT32; vb->textual_opt.len += str_int ((int32_t)n, AFTERENT (char, vb->textual_opt)); break; }
        case 'I': { uint32_t n = NEXT_UINT32; vb->textual_opt.len += str_int (         n, AFTERENT (char, vb->textual_opt)); break; }
        case 'f': { buf_add (&vb->textual_opt, special_float, 2); 
                    uint32_t n = NEXT_UINT32; // n is the 4 bytes of the little endian float, construct and int, and switch to machine endianity 
                    n = LTEN32 (n);           // switch back to Little Endian as it was in the BAM file
                    /* integer as text */     vb->textual_opt.len += str_int (         n, AFTERENT (char, vb->textual_opt)); break; }
        default: ABORT ("Error in bam_rewrite_one_optional_number: enrecognized Optional field type '%c' (ASCII %u) in vb=%u", 
                        type, type, vb->vblock_i);
    }    

    return next_field;
} 

const char *bam_get_one_optional (VBlockSAM *vb, const char *next_field,
                                  const char **tag, char *type, const char **value, unsigned *value_len) // out
{
    *tag  = next_field;
    *type = next_field[2];
    next_field += 3;

    if (*type == 'Z' || *type == 'H') {
        *value = next_field;
        *value_len = strlen (*value);
        return next_field + *value_len + 1; // +1 for \0
    }
    
    else if (*type == 'A') {
        *value = next_field;
        *value_len = 1;
        return next_field + 1;
    } 

    uint32_t max_len = (*type == 'B') ? (12 * GET_UINT32 (next_field) + 10) : // worst case scenario for item: "-1000000000,"
                       30;
    buf_alloc (vb, &vb->textual_opt, max_len, 2, "textual_opt"); // a rather inefficient max in case of arrays, to do: tighten the calculation

    if (*type != 'B')
        next_field = bam_rewrite_one_optional_number (vb, next_field, *type);

    else { // 'B'
        char subtype = *next_field++; // type of elements of array
        uint32_t num_elements = NEXT_UINT32;

        NEXTENT (uint8_t, vb->textual_opt) = subtype;
        
        for (uint32_t i=0; i < num_elements; i++) {
            NEXTENT (char, vb->textual_opt) = ',';
            next_field = bam_rewrite_one_optional_number (vb, next_field, subtype);
        }
    }    

    *value     = vb->textual_opt.data;
    *value_len = vb->textual_opt.len;

    return next_field;
}

const char *bam_seg_txt_line (VBlock *vb_, const char *alignment /* BAM terminology for one line */, bool *has_13_unused)   
{
    VBlockSAM *vb = (VBlockSAM *)vb_;
    ZipDataLineSAM *dl = DATA_LINE (vb->line_i);
    const char *next_field = alignment;

    // *** ingest BAM alignment fixed-length fields ***
    uint32_t block_size = NEXT_UINT32;
    int32_t ref_id      = (int32_t)NEXT_UINT32;     // corresponding to CHROMs in the BAM header
    PosType this_pos    = 1 + (int32_t)NEXT_UINT32; // pos in BAM is 0 based, -1 for unknown 
    uint8_t l_read_name = NEXT_UINT8;               // QNAME length
    uint8_t mapq        = NEXT_UINT8;
    uint16_t bin        = NEXT_UINT16;
    uint16_t n_cigar_op = NEXT_UINT16;
    uint16_t flag       = NEXT_UINT16;
    uint32_t l_seq      = NEXT_UINT32;              // note: we stick with the same logic as SAM for consistency - dl->seq_len is determined by CIGAR 
    int32_t next_ref_id = (int32_t)NEXT_UINT32;     // corresponding to CHROMs in the BAM header
    PosType next_pos    = 1 + (int32_t)NEXT_UINT32; // pos in BAM is 0 based, -1 for unknown
    int32_t tlen        = (int32_t)NEXT_UINT32;

    // *** segment fixed-length fields ***

    bam_seg_ref_id (vb_, SAM_RNAME, ref_id, -1); // ref_id (RNAME)

    // note: pos can have a value even if ref_id=-1 (RNAME="*") - this happens if a SAM with a RNAME that is not in the header is converted to BAM with samtools
    seg_pos_field (vb_, SAM_POS, SAM_POS, false, 0, 0, this_pos, sizeof (uint32_t)); // POS

    seg_integer (vb, SAM_MAPQ, mapq, true); // MAPQ

    seg_integer (vb, SAM_FLAG, flag, true); // FLAG
    
    bam_seg_ref_id (vb_, SAM_RNEXT, next_ref_id, ref_id); // RNEXT

    seg_pos_field (vb_, SAM_PNEXT, SAM_POS, false, 0, 0, next_pos, sizeof (uint32_t)); // PNEXT

    sam_seg_tlen_field (vb, 0, 0, (int64_t)tlen, vb->contexts[SAM_PNEXT].last_delta, dl->seq_len); // TLEN

    seg_compound_field ((VBlockP)vb, &vb->contexts[SAM_QNAME], next_field, // QNAME
                        l_read_name-1, false, 0, 2 /* account for \0 and l_read_name */); 
    next_field += l_read_name; // inc. \0

    // *** ingest & segment variable-length fields ***

    // CIGAR
    bam_rewrite_cigar (vb, n_cigar_op, (uint32_t*)next_field); // re-write BAM format CIGAR as SAM textual format in vb->textual_cigar
    sam_analyze_cigar (vb, vb->textual_cigar.data, vb->textual_cigar.len, &dl->seq_len, &vb->ref_consumed, &vb->ref_and_seq_consumed);
    next_field += n_cigar_op * sizeof (uint32_t);

    // Segment BIN after we've gathered bin, flags, pos and vb->ref_confumed (and before sam_seg_seq_field which ruins vb->ref_consumed)
    bam_seg_bin (vb, bin, flag, this_pos);

    // SEQ - calculate diff vs. reference (denovo or loaded)
    bam_rewrite_seq (vb, l_seq, next_field);
    sam_seg_seq_field (vb, SAM_SQBITMAP, vb->textual_seq.data, vb->textual_seq.len, this_pos, vb->last_cigar, 0, vb->textual_seq.len, 
                       vb->last_cigar, (l_seq+1)/2 + sizeof (uint32_t) /* account for l_seq and seq fields */);
    next_field += (l_seq+1)/2; 

    // QUAL
    bool has_qual=false;
    if (l_seq)
        has_qual = bam_rewrite_qual ((uint8_t *)next_field, l_seq); // add 33 to Phred scores to make them ASCII
    
    if (has_qual) // case we have both SEQ and QUAL
        sam_seg_qual_field (vb, dl, next_field, l_seq, l_seq /* account for qual field */ );

    else { // cases 1. were both SEQ and QUAL are '*' (seq_len=0) and 2. SEQ exists, QUAL not (bam_rewrite_qual returns false)
        *(char *)alignment = '*'; // overwrite as we need it somewhere in txt_data
        sam_seg_qual_field (vb, dl, alignment, 1, l_seq /* account of l_seq 0xff */);
    }
    next_field += l_seq; 

    // finally we can segment the textual CIGAR now (including if n_cigar_op=0)
    bam_seg_cigar_field (vb, dl, l_seq, n_cigar_op);

    // OPTIONAL fields - up to MAX_SUBFIELDS of them
    next_field = sam_seg_optional_all (vb, dl, next_field, 0,0,0, alignment + block_size + sizeof (uint32_t));
    
    buf_free (&vb->textual_cigar);
    buf_free (&vb->textual_seq);

    return next_field;
}
