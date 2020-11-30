// ------------------------------------------------------------------
//   sam_bam.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "digest.h"
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
#include "flags.h"
#include "profiler.h"

// set an estimated number of lines, so that seg_all_data_lines doesn't call seg_estimate_num_lines which
// won't work for BAM as it scans for newlines. Then proceed to sam_seg_initialize
void bam_seg_initialize (VBlock *vb)
{
#   define NUM_LINES_IN_TEST 10 // take an average of this number of lines

    uint32_t next=0, line_i=0;
    while (next < vb->txt_data.len && line_i < NUM_LINES_IN_TEST) {
        next += 4 + GET_UINT32 (ENT (char, vb->txt_data, next)); // block_len excludes the block_len field length itself (4)
        line_i++;
    }

    // estimated number of BAM alignments in the VB, based on the average length of the first few
    vb->lines.len = line_i ? (1.2 * vb->txt_data.len / (next / line_i)) : 0;

    sam_seg_initialize (vb);
}

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
// if first_i > 0, we attempt to heuristically detect the start of a BAM alignment.
int32_t bam_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i)
{
    ASSERT (*i >= 0 && *i < vb->txt_data.len, "Error in def_unconsumed: *i=%d is out of range [0,%"PRIu64"]", *i, vb->txt_data.len);

    *i = MIN (*i, vb->txt_data.len - sizeof(BAMAlignmentFixed));

    // find the first alignment in the data (going backwards) that is entirely in the data - 
    // we identify and alignment by l_read_name and read_name
    for (; *i >= (int32_t)first_i; (*i)--) {
        const BAMAlignmentFixed *aln = (const BAMAlignmentFixed *)ENT (char, vb->txt_data, *i);

        uint32_t block_size = LTEN32 (aln->block_size);
        uint32_t l_seq      = LTEN32 (aln->l_seq);
        uint16_t n_cigar_op = LTEN16 (aln->n_cigar_op);

        // test to see block_size makes sense
        if ((uint64_t)*i + (uint64_t)block_size + 4 > (uint64_t)vb->txt_data.len || // 64 bit arith to catch block_size=-1 that will overflow in 32b
            block_size < sizeof (BAMAlignmentFixed) + 4*n_cigar_op  + aln->l_read_name + l_seq + (l_seq+1)/2)
            continue;

        // test to see l_read_name makes sense
        if (LTEN32 (aln->l_read_name) < 2 ||
            &aln->read_name[aln->l_read_name] > AFTERENT (char, vb->txt_data)) continue;

        // test pos
        int32_t pos = LTEN32 (aln->pos);
        if (pos < -1) continue;

        // test read_name    
        if (aln->read_name[aln->l_read_name-1] != 0 || // nul-terminated
            !str_is_in_range (aln->read_name, aln->l_read_name-1, '!', '~')) continue;  // all printable ascii (per SAM spec)

        // test l_seq vs seq_len implied by cigar
        if (aln->l_seq && n_cigar_op) {
            uint32_t seq_len_by_cigar=0;
            uint32_t *cigar = (uint32_t *)((uint8_t *)(aln+1) + aln->l_read_name);
            for (uint16_t cigar_op_i=0; cigar_op_i < n_cigar_op; cigar_op_i++) {
                uint8_t cigar_op = *(uint8_t *)&cigar[cigar_op_i] & 0xf; // LSB by Little Endian - take 4 LSb
                uint32_t op_len = cigar[cigar_op_i] >> 4;
                if (cigar_lookup_bam[cigar_op] & CIGAR_CONSUMES_QUERY) seq_len_by_cigar += op_len; 
            }
            if (l_seq != seq_len_by_cigar) continue;
        }

        // Note: we don't use add aln->bin calculation because in some files we've seen data that doesn't
        // agree with our formula. see comment in bam_reg2bin

        // all tests passed - this is indeed an alignment
        return vb->txt_data.len - (*i + LTEN32 (aln->block_size) + 4); // everything after this alignment is "unconsumed"
    }

    return -1; // we can't find any alignment - need more data (lower first_i)
}

void bam_seg_bin (VBlockSAM *vb, uint16_t bin /* used only in bam */, uint16_t sam_flag, PosType this_pos)
{
    bool segment_unmapped = (sam_flag & 0x4);
    bool is_bam = IS_BAM;

    PosType last_pos = segment_unmapped ? this_pos : (this_pos + vb->ref_consumed - 1);
    uint16_t reg2bin = bam_reg2bin (this_pos, last_pos); // zero-based, half-closed half-open [start,end)

    if (!is_bam || (last_pos <= MAX_POS_SAM && reg2bin == bin))
        seg_by_did_i (vb, ((char []){ SNIP_SPECIAL, SAM_SPECIAL_BIN }), 2, SAM_BAM_BIN, is_bam ? sizeof (uint16_t) : 0)
    
    else {
#ifdef DEBUG // we show this warning only in DEBUG because I found actual files that have edge cases that don't work with our formula (no harm though)
        static bool warning_shown = false;
        if (!warning_shown) { 
            WARN ("FYI: bad bin value in vb=%u vb->line_i=%u: this_pos=%"PRId64" ref_consumed=%u flag=%u last_pos=%"PRId64": bin=%u but reg2bin=%u. No harm. This warning will not be shown again for this file.",
                  vb->vblock_i, vb->line_i, this_pos, vb->ref_consumed, sam_flag, last_pos, bin, reg2bin);
            warning_shown = true;
        }
#endif
        seg_integer (vb, SAM_BAM_BIN, bin, is_bam);
        vb->contexts[SAM_BAM_BIN].flags.store = STORE_INT;
    }
}

static inline void bam_seg_ref_id (VBlockP vb, DidIType did_i, int32_t ref_id, int32_t compare_to_ref_i)
{
    ASSERT (ref_id >= -1 && ref_id < (int32_t)header_contigs.len, "Error in bam_seg_ref_id: vb=%u line_i=%u: encountered ref_id=%d but header has only %u contigs",
            vb->vblock_i, vb->line_i, ref_id, (uint32_t)header_contigs.len);

    // get snip and snip_len
    DECLARE_SNIP;
    if (ref_id >= 0) {
        if (ref_id == compare_to_ref_i) {
            snip = "=";
            snip_len = 1;
        }
        else {
            RefContig *rc = ENT (RefContig, header_contigs, ref_id);
            snip = ENT (char, header_contigs_dict, rc->char_index);
            snip_len = rc->snip_len;
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

const char *bam_seg_txt_line (VBlock *vb_, const char *alignment /* BAM terminology for one line */,
                              uint32_t remaining_txt_len, bool *has_13_unused)   
{
    VBlockSAM *vb = (VBlockSAM *)vb_;
    ZipDataLineSAM *dl = DATA_LINE (vb->line_i);
    const char *next_field = alignment;

    // *** ingest BAM alignment fixed-length fields ***
    uint32_t block_size = NEXT_UINT32;

    // a non-sensical block_size might indicate an false-positive identification of a BAM alignment in bam_unconsumed
    ASSERT (block_size + 4 >= sizeof (BAMAlignmentFixed) && block_size + 4 <= remaining_txt_len, 
            "Error in bam_seg_txt_line: vb=%u line_i=%u (block_size+4)=%u is out of range - too small, or goes beyond end of txt data: remaining_txt_len=%u",
            vb->vblock_i, vb->line_i, block_size+4, remaining_txt_len);

    int32_t ref_id      = (int32_t)NEXT_UINT32;     // corresponding to CHROMs in the BAM header
    PosType this_pos    = 1 + (int32_t)NEXT_UINT32; // pos in BAM is 0 based, -1 for unknown 
    uint8_t l_read_name = NEXT_UINT8;               // QNAME length
    uint8_t mapq        = NEXT_UINT8;
    uint16_t bin        = NEXT_UINT16;
    uint16_t n_cigar_op = NEXT_UINT16;
    uint16_t sam_flag   = NEXT_UINT16;              // not to be confused with our global var "flag"
    uint32_t l_seq      = NEXT_UINT32;              // note: we stick with the same logic as SAM for consistency - dl->seq_len is determined by CIGAR 
    int32_t next_ref_id = (int32_t)NEXT_UINT32;     // corresponding to CHROMs in the BAM header
    PosType next_pos    = 1 + (int32_t)NEXT_UINT32; // pos in BAM is 0 based, -1 for unknown
    int32_t tlen        = (int32_t)NEXT_UINT32;

    // *** segment fixed-length fields ***

    bam_seg_ref_id (vb_, SAM_RNAME, ref_id, -1); // ref_id (RNAME)

    // note: pos can have a value even if ref_id=-1 (RNAME="*") - this happens if a SAM with a RNAME that is not in the header is converted to BAM with samtools
    seg_pos_field (vb_, SAM_POS, SAM_POS, false, 0, 0, this_pos, sizeof (uint32_t)); // POS
    if (ref_id >= 0) sam_seg_verify_pos (vb_, this_pos);

    seg_integer (vb, SAM_MAPQ, mapq, true); // MAPQ

    seg_integer (vb, SAM_FLAG, sam_flag, true); // FLAG
    
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
    bam_seg_bin (vb, bin, sam_flag, this_pos);

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
