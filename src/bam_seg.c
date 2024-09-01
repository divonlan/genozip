// ------------------------------------------------------------------
//   sam_bam.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "sam_private.h"
#include "txtfile.h"
#include "libdeflate_1.19/libdeflate.h"

void bam_seg_initialize (VBlockP vb)
{
    sam_seg_initialize (vb);

    // we store special quality values at the end of txt_data, as they need to be somewhere in txt_data
    buf_alloc (vb, &vb->txt_data, 1, 0, char, 0, 0); // add 1 character after the end of txt_data
    *BAFTtxt = '*'; // missing qual;

    if (!segconf.running && line_textual_cigars_used)
        buf_alloc (vb, &VB_SAM->line_textual_cigars, 0, segconf.sam_cigar_len * vb->lines.len32 / (segconf.is_long_reads ? 4 : 1),/*divide in case sam_cigar_len is not representative*/
                   char, CTX_GROWTH, "line_textual_cigars");
}

// detect if a generic file is actually a BAM
bool is_bam (STRp(header), bool *need_more)
{
    if (header_len < STRLEN(BAM_MAGIC)) {
        *need_more = true;
        return false;
    }

    return str_isprefix_(STRa(header), _S(BAM_MAGIC));
}

// detect if a generic file is actually a CRAM
bool is_cram (STRp(header), bool *need_more)
{
    if (header_len < STRLEN(CRAM_MAGIC)) {
        *need_more = true;
        return false;
    }

    return str_isprefix_(STRa(header), _S(CRAM_MAGIC));
}

static int32_t bam_unconsumed_scan_forwards (VBlockP vb)
{
    ARRAY (char, txt, vb->txt_data);
    
    if (txt_len < sizeof (BAMAlignmentFixed)) return -1; // this VB doesn't not even contain one single full alignment

    uint32_t aln_size=0, i;
    for (i=0 ; i < txt_len-3; i += aln_size) 
        aln_size = GET_UINT32_((BAMAlignmentFixed *)&txt[i], block_size) + 4;

    if (aln_size > txt_len) 
        return -1; // this VB doesn't not even contain one single full alignment

    else if (i == txt_len)
        return 0;  // we will consume all data - nothing to pass to next VB
    
    else 
        return aln_size - (i - txt_len); // we pass the data of the final, partial, alignment to the next VB
}

static int32_t bam_unconsumed_scan_backwards (VBlockP vb, uint32_t first_i)
{
    int32_t last_i = Ltxt - sizeof(BAMAlignmentFixed);

    // find the first alignment in the data (going backwards) that is entirely in the data - 
    // we identify and alignment by l_read_name and read_name
    for (; last_i >= (int32_t)first_i; (last_i)--) {
        const BAMAlignmentFixed *aln = (const BAMAlignmentFixed *)Btxt (last_i);

        uint32_t block_size = LTEN32 (aln->block_size);
        if (block_size > 100000000) continue; // quick short-circuit - more than 100M for one alignment - clearly wrong

        uint32_t l_seq      = LTEN32 (aln->l_seq);
        uint16_t n_cigar_op = LTEN16 (aln->n_cigar_op);

        // test to see block_size makes sense
        if ((uint64_t)last_i + (uint64_t)block_size + 4 > (uint64_t)vb->txt_data.len || // 64 bit arith to catch block_size=-1 that will overflow in 32b
            block_size + 4 < sizeof (BAMAlignmentFixed) + 4*n_cigar_op  + aln->l_read_name + l_seq + (l_seq+1)/2)
            continue;

        // test to see l_read_name makes sense
        if (LTEN32 (aln->l_read_name) < 2 ||
            &aln->read_name[aln->l_read_name] > BAFTtxt) continue;

        // test pos
        int32_t pos = LTEN32 (aln->pos);
        if (pos < -1) continue;

        // test read_name    
        if (aln->read_name[aln->l_read_name-1] != 0 || // nul-terminated
            !str_is_in_range (aln->read_name, aln->l_read_name-1, '!', '~')) continue;  // all printable ascii (per SAM spec)

        // final test option 1: test l_seq vs seq_len implied by cigar
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
        
        // final test option 2: if QNAME is of a known flavor - test for that flavor
        else if (segconf.qname_flavor[QNAME1]) {
            if (qname_test_flavor (aln->read_name, aln->l_read_name-1, QNAME1, segconf.qname_flavor[QNAME1], true)) // >0 if not qname_flavor
                continue;
        }

        // we don't have a final test - skip
        else continue;

        // Note: we don't use add aln->bin calculation because in some files we've seen data that doesn't
        // agree with our formula. see comment in bam_reg2bin

        // all tests passed - this is indeed an alignment
        return Ltxt - (last_i + LTEN32 (aln->block_size) + 4); // everything after this alignment is "unconsumed"
    }

    return -1; // we can't find any alignment - need more data (lower first_i)    
}  

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
// if first_i > 0, we attempt to heuristically detect the start of a BAM alignment.
int32_t bam_unconsumed (VBlockP vb, uint32_t first_i)
{
    ASSERTNOTZERO (Ltxt);

    int32_t result;

    // if we have all the data - search forward - faster in BAM as we have aln.block_size
    if (first_i == 0)
        result = bam_unconsumed_scan_forwards (vb);

    // stringent -either CIGAR needs to match seq_len, or qname needs to match flavor
    else
        result = bam_unconsumed_scan_backwards (vb, first_i); 

    return result; // if -1 - we will be called again with more data
}

static rom bam_dump_alignment (VBlockSAMP vb, rom alignment, rom after)
{
    buf_free (vb->scratch); // feel free to use scratch bc this is called in an ASSERT before aborting

    buf_set_shared (&vb->txt_data);
    buf_overlay_partial (vb, &vb->scratch, &vb->txt_data, BNUMtxt(alignment), "alignment_buf");
    vb->scratch.len = after - alignment;

    rom fn = "bad_alignment.bam";
    buf_dump_to_file (fn, &vb->scratch, 1, false, false, false, false);
    buf_free (vb->scratch);

    return fn;
}

uint32_t bam_split_aux (VBlockSAMP vb, rom alignment, rom aux, rom after_aux, rom *auxs, uint32_t *aux_lens)
{
    uint32_t n_auxs = 0;
    while (aux < after_aux) {
        auxs[n_auxs] = aux;

        static unsigned const size[256] = { ['A']=1, ['c']=1, ['C']=1, ['s']=2, ['S']=2, ['i']=4, ['I']=4, ['f']=4 };
        
        ASSERT (after_aux - aux >= 4, "%s: Failed to parse BAM AUX fields after n_aux=%u. Only %u bytes remaining in the field, but expecting at least 4 bytes to contain any auxilliary field. Possibly this BAM file has an incorrect value in the block_size field of this alignment. Dumped alignment to %s", 
                LN_NAME, n_auxs, (int)(after_aux-aux), bam_dump_alignment (vb, alignment, after_aux));

        if (aux[2] == 'Z' || aux[2] == 'H') {
            SAFE_NUL (after_aux);
            aux_lens[n_auxs] = strlen (aux+3) + 4; // add tag[2], type and \0
            SAFE_RESTORE;
        }
        else if (aux[2] == 'B') {
            uint32_t array_len = GET_UINT32 (aux+4);
            aux_lens[n_auxs] = 8 + array_len * size[(int)aux[3]]; 
        }
        
        else if (size[(int)aux[2]])
            aux_lens[n_auxs] = 3 + size[(int)aux[2]];
            
        else
            ABORT ("%s: Unrecognized aux type '%c' (ASCII %u). Dumped alignment to %s", LN_NAME, aux[2], (uint8_t)aux[2], bam_dump_alignment (vb, alignment, after_aux));
        
        aux += aux_lens[n_auxs++];
    }

    ASSERT (aux == after_aux, "%s: overflow while parsing auxilliary fields. Dumped alignment to %s", LN_NAME, bam_dump_alignment (vb, alignment, after_aux));

    return n_auxs;
}

void bam_seg_BIN (VBlockSAMP vb, ZipDataLineSAMP dl, uint16_t bin /* used only in bam */, bool is_bam)
{
    PosType32 this_pos = dl->POS;
    PosType32 last_pos = dl->FLAG.unmapped ? this_pos : (this_pos + vb->ref_consumed - 1);
    uint16_t reg2bin = bam_reg2bin (this_pos, last_pos); // zero-based, half-closed half-open [start,end)

    if (!is_bam || (last_pos <= MAX_POS_SAM && reg2bin == bin))
        seg_by_did (VB, ((char []){ SNIP_SPECIAL, SAM_SPECIAL_BIN }), 2, SAM_BAM_BIN, is_bam ? sizeof (uint16_t) : 0);
    
    else {
#ifdef DEBUG // we show this warning only in DEBUG because I found actual files that have edge cases that don't work with our formula (no harm though)
        WARN_ONCE ("FYI: %s: bad bin value in: this_pos=%d ref_consumed=%u flag=%u last_pos=%d: bin=%u but reg2bin=%u. No harm. This warning will not be shown again for this file.",
                    LN_NAME, this_pos, vb->ref_consumed, dl->FLAG.value, last_pos, bin, reg2bin);
#endif
        seg_integer_as_snip (vb, SAM_BAM_BIN, bin, is_bam);
        CTX(SAM_BAM_BIN)->flags.store = STORE_INT;
    }
}

static inline void bam_seg_ref_id (VBlockSAMP vb, ZipDataLineSAMP dl, Did did_i, int32_t ref_id, int32_t compare_to_ref_i)
{
    ASSERT (ref_id == -1 || (sam_hdr_contigs && IN_RANGE (ref_id, 0, sam_hdr_contigs->contigs.len32)), 
            "%s: encountered %s.ref_id=%d but header has only %u contigs%s", 
            LN_NAME, CTX(did_i)->tag_name, ref_id, sam_hdr_contigs ? sam_hdr_contigs->contigs.len32 : 0,
            MP(LONGRANGER) ? ". This is a known longranger bug (samtools won't accept this file either)." : "");

    // get snip and snip_len
    STR0(snip);
    if (ref_id >= 0) {
        if (ref_id == compare_to_ref_i) {
            snip = "=";
            snip_len = 1;            
        }
        else 
            snip = contigs_get_name (sam_hdr_contigs, ref_id, &snip_len);
    }
    else { 
        snip = "*";
        snip_len = 1;
    }

    if (did_i == SAM_RNAME) sam_seg_RNAME (vb, dl, STRa(snip), true, sizeof (int32_t));
    else                    sam_seg_RNEXT (vb, dl, STRa(snip),       sizeof (int32_t));
}

// Rewrite the QUAL field - add +33 to Phred scores to make them ASCII
static inline bool bam_rewrite_qual (uint8_t *qual, uint32_t l_seq)
{
    if (qual[0] == 0xff) return false; // in case SEQ is present but QUAL is omitted, all qual is 0xFF

    for (uint32_t i=0; i < l_seq; i++)
        qual[i] += 33;

    return true;
}

static inline void bam_set_missing_qual_type (VBlockSAMP vb, bytes qual, uint32_t l_seq)
{
    if (l_seq <= 1) return;

    uint8_t filler = qual[1]; // SAM spec - needs to be 0xff. Genozip also supports - as created by pysam - 0x00.
    ASSERT (filler == 0 || filler == 0xff, "%s: Non-compliant QUAL field. l_seq=%u filler=%u (only 0 and 0xff are accepted)",
            LN_NAME, l_seq, filler);

    for (uint32_t i=2; i < l_seq; i++)
        ASSERT (qual[i] == filler, "%s: Non-compliant QUAL field. l_seq=%u filler=%u qual[%u]=%u",
                LN_NAME, l_seq, filler, i, qual[i]);
    
    // note: in older versions of pysam, missing QUAL was 0xff followed by 0s (SAM specficiation requires all bytes to be 0xff)
    if (!filler) {
        if (!segconf.pysam_qual)
            __atomic_store_n (&segconf.pysam_qual, (bool)true, __ATOMIC_RELAXED); 
    }
    else
        // note: this test is not air-tight related to another thread setting pysam_qual after this test, but we really don't expect to see files that are mixed pysam/standard
        ASSINP0 (!segconf.pysam_qual, "Unsupported BAM format: some lines have pysam-style missing QUAL, and some standard missing QUAL");
}

void bam_get_one_aux (VBlockSAMP vb, int16_t idx,
                      rom *tag, char *type, char *array_subtype, // out
                      pSTRp(value), ValueType *numeric) // out - one of these depending on the type
{
    rom aux = vb->auxs[idx];
    *tag  = aux;
    *type = aux[2];     // c, C, s, S, i, I, f, A, Z, H or B
    aux += 3;
    *array_subtype = 0; // c, C, s, S, i, I, f
    *value = NULL;

    switch (*type) {
        // in case of an numeric type, we pass the value as a ValueType
        case 'i': *value_len = 4; numeric->i     = (int32_t)GET_UINT32 (aux); break;
        case 'I': *value_len = 4; numeric->i     = GET_UINT32 (aux);          break; 
        case 'f': *value_len = 4; numeric->f32.f = GET_FLOAT32 (aux);         break; // note: this DOES NOT result in the correct value in last_value.f
        case 's': *value_len = 2; numeric->i     = (int16_t)GET_UINT16 (aux); break;
        case 'S': *value_len = 2; numeric->i     = GET_UINT16 (aux);          break;
        case 'c': *value_len = 1; numeric->i     = (int8_t)*aux;              break;
        case 'C': *value_len = 1; numeric->i     = (uint8_t)*aux;             break;
        case 'Z': 
        case 'H': *value_len = vb->aux_lens[idx] - 4; *value = aux;           break; // value_len excludes the terminating \0
        case 'A': *value_len = 1; *value = aux;                               break;

        // in case of a numerical value we pass the data as is, in machine endianity
        case 'B':
            *array_subtype = *aux++; // type of elements of array
            *value_len = GET_UINT32(aux);  // number of elements, not number of bytes
            *value = aux + sizeof (uint32_t);
            
            // switch from BAM's little endian to machine endianity
            if (!flag.is_lten) 
                switch (*array_subtype) {
                    case 'i': case 'I': case 'f': 
                        for (int i=0; i < *value_len; i++) ((uint32_t*)(*value))[i] = LTEN32(((uint32_t*)(*value))[i]);
                        break;

                    case 's': case 'S': 
                        for (int i=0; i < *value_len; i++) ((uint32_t*)(*value))[i] = LTEN16(((uint32_t*)(*value))[i]);
                        break;

                    default: break;
                }

            ASSERT (aux_width[(uint8_t)*array_subtype], "%s: Invalid array type: '%c' for field \"%c%c\"", 
                    LN_NAME, *array_subtype, (*tag)[0], (*tag)[1]);

            break;

        default:
            ABORT ("%s: Invalid field type: '%c' for field \"%c%c\"", LN_NAME, *type, (*tag)[0], (*tag)[1]);
    }
}

// segconf
static uint32_t bam_segconf_get_transated_sam_line_len (VBlockSAMP vb, ZipDataLineSAMP dl, SamTlenType tlen)
{
    unsigned rname_len, rnext_len;

    if (dl->RNAME == -1) 
        rname_len = 1; // "*" in SAM
    else 
        contigs_get_name (sam_hdr_contigs, dl->RNAME, &rname_len);

    if (dl->RNAME == dl->RNEXT || dl->RNEXT == -1)
        rnext_len = 1; // "=" or "*" in SAM
    else
        contigs_get_name (sam_hdr_contigs, dl->RNEXT, &rnext_len);
    
    uint32_t sam_line_len =  
        dl->QNAME.len + rname_len + rnext_len +
        2 * 9 +             // POS and PNEXT: account for 9 characters - as segconf is biased to small values that are not representative of the file 
        str_int_len (dl->FLAG.value) + str_int_len (dl->MAPQ) + vb->textual_cigar.len32 + 
        str_int_len (tlen) + dl->SEQ.len + (dl->no_qual ? 1 : dl->SEQ.len) + 
        (11 + vb->n_auxs) + // \t or \n seperators
        (5 * vb->n_auxs);   // AUX field prefix, eg AS:i: (we account for the extra 2 characters of B below)

    for (uint32_t idx=0; idx < vb->n_auxs; idx++) {
        STR0(value);
        ValueType numeric;
        rom tag;
        char bam_type = 0, array_subtype = 0;

        bam_get_one_aux (vb, idx, &tag, &bam_type, &array_subtype, pSTRa (value), &numeric);

        #define ASSUMED_FLOAT_LEN 7 // est SAM representation of a BAM floating point
        if (bam_type == 'f')
            sam_line_len += ASSUMED_FLOAT_LEN; 

        else if (!value) // c, C, s, S, i, I
            sam_line_len += str_int_len (numeric.i); 

        else if (bam_type == 'Z' || bam_type == 'H' || bam_type == 'A')
            sam_line_len += value_len;

        else if (bam_type == 'B') {
            sam_line_len += 2 + value_len; // added eg ":c" and commas

            switch (array_subtype) {
                // TO DO: support calculation of Big Endian - numbers need to be flipped before taking the length
                case 'f' : sam_line_len += value_len * ASSUMED_FLOAT_LEN; break; 
                case 'I' : for (int i=0; i < value_len; i++) sam_line_len += str_int_len (GET_UINT32 (value + i * sizeof(uint32_t))); break;
                case 'i' : for (int i=0; i < value_len; i++) sam_line_len += str_int_len ((int32_t)GET_UINT32 (value + i * sizeof(int32_t))); break;
                case 'S' : for (int i=0; i < value_len; i++) sam_line_len += str_int_len (GET_UINT16 (value + i * sizeof(uint16_t))); break;
                case 's' : for (int i=0; i < value_len; i++) sam_line_len += str_int_len ((int16_t)GET_UINT16 (value + i * sizeof(int16_t))); break;
                case 'C' : for (int i=0; i < value_len; i++) sam_line_len += str_int_len (*(uint8_t *)(value + i)); break;
                case 'c' : for (int i=0; i < value_len; i++) sam_line_len += str_int_len (*(int8_t  *)(value + i)); break;
                default  : ABORT ("unrecognized BAM array_subtype='%c'(%u)", array_subtype, array_subtype);
            }
        }

        else
            ABORT ("unrecognized BAM type='%c'(%u)", bam_type, bam_type);
    }

    return sam_line_len;
}

rom bam_seg_txt_line (VBlockP vb_, rom alignment /* BAM terminology for one line */,
                      uint32_t remaining_txt_len, bool *has_13_unused)   
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    sam_reset_line (VB);

    // case: --show - only print BAM lines, no segging
    if (flag.show_bam) return bam_show_line (vb, alignment, remaining_txt_len);

    ZipDataLineSAMP dl = DATA_LINE (vb->line_i);
    rom next_field = alignment;
   
    // BAM alignment fixed-length fields 
    uint32_t block_size = NEXT_UINT32;

    WordIndex prev_line_chrom = vb->chrom_node_index;

    // a non-sensical block_size might indicate an false-positive identification of a BAM alignment in bam_unconsumed
    ASSERT (block_size + 4 >= sizeof (BAMAlignmentFixed) && block_size + 4 <= remaining_txt_len, 
            "%s: (block_size+4)=%u is out of range - too small, or goes beyond end of txt data: txt_data.len=%u remaining_txt_len=%u",
            LN_NAME, block_size+4, vb->txt_data.len32, remaining_txt_len);

    rom after = alignment + block_size + sizeof (uint32_t);

    // note: in BAM, we calculate here instead of in seg_all_data_lines, bc we rewrite qual
    if (flag.debug_lines) 
        vb->debug_line_hash = adler32 (2, alignment, after - alignment); 

    dl->RNAME = vb->chrom_node_index = (int32_t)NEXT_UINT32; // corresponding to CHROMs in the BAM header. -1 in BAM means '*' (no RNAME) - which luckily is WORD_INDEX_NONE.    
    
    // expecting all contigs to be defined in SAM header, and hence in ol_nodes/ol_dict 
    ASSERT (vb->chrom_node_index >= -1 && vb->chrom_node_index < (int32_t)CTX(SAM_RNAME)->ol_nodes.len32, 
            "%s: RNAME=%d ∉ [-1,%d]", LN_NAME, vb->chrom_node_index, (int)CTX(SAM_RNAME)->ol_nodes.len32-1);

    dl->POS              = 1 + (int32_t)NEXT_UINT32; // pos in BAM is 0 based, -1 for unknown 
    uint8_t l_read_name  = NEXT_UINT8;               // QNAME length
    dl->MAPQ             = NEXT_UINT8;
    uint16_t bin         = NEXT_UINT16;
    uint16_t n_cigar_op  = NEXT_UINT16;
    dl->FLAG.value       = NEXT_UINT16;              // not to be confused with our global var "flag"
    uint32_t l_seq       = NEXT_UINT32;              // note: we stick with the same logic as SAM for consistency - dl->SEQ.len is determined by CIGAR 
    dl->RNEXT            = (int32_t)NEXT_UINT32;     // corresponding to CHROMs in the BAM header
    PosType32 next_pos  = 1 + (int32_t)NEXT_UINT32;  // pos in BAM is 0 based, -1 for unknown
    SamTlenType tlen     = (SamTlenType)NEXT_UINT32;
    rom read_name        = next_field;
    dl->QNAME            = (TxtWord){ .index = BNUMtxt (read_name), .len = l_read_name-1 }; // -1 don't count \0
    BamCigarOp *cigar    = (BamCigarOp *)(read_name + l_read_name); // note: the "cigar" pointer might be mis-aligned, but we don't de-reference it
    bytes seq            = (uint8_t *)(cigar + n_cigar_op);
    dl->SEQ.index        = BNUMtxt (seq);
    rom qual             = (rom)seq + (l_seq+1)/2;
    dl->QUAL             = (TxtWord){ .index = BNUMtxt (qual), .len = l_seq };
    rom aux              = qual + l_seq;
    
    vb->RNEXT_is_equal = (dl->RNAME == dl->RNEXT);
    vb->cigar_missing  = n_cigar_op==0; // needed for sam_seg_is_gc_line
 
    // split auxillary fields
    STR_ARRAY (aux,MAX_FIELDS); 
    vb->n_auxs   = bam_split_aux (vb, alignment, aux, after, auxs, aux_lens);
    vb->auxs     = auxs;    // note: pointers to data on the stack
    vb->aux_lens = aux_lens;
                        
    sam_seg_idx_aux (vb);
    
    if (vb->chrom_node_index != WORD_INDEX_NONE) 
        ctx_get_vb_snip_ex (CTX(SAM_RNAME), vb->chrom_node_index, pSTRa(vb->chrom_name));

    else {
        vb->chrom_name = "";
        vb->chrom_name_len = 0;
    }

    // convert BAM seq format to SAM
    buf_alloc (vb, &vb->textual_seq, 0, l_seq+2/* +1 for last half-byte and \0 */, char, 1.5, "textual_seq");    
    {
    START_TIMER; // time here and not in the function, because this is also called from the LongR codec
    bam_seq_to_sam (vb, seq, l_seq, false, true, &vb->textual_seq);
    vb->textual_seq_str = B1STc(vb->textual_seq); 
    COPY_TIMER(bam_seq_to_sam);
    }

    // analyze (but not seg yet) cigar
    buf_append (vb, vb->binary_cigar, BamCigarOp, cigar, n_cigar_op, "binary_cigar");
    uint32_t seq_len;
    bam_seg_cigar_analyze (vb, dl, &seq_len);
    dl->SEQ.len = seq_len; // do it this way to avoid compiler warning

    // if this is a secondary / supplamentary read (aka Dependent) or a read that has an associated sec/sup read 
    // (aka Primary) - move the line to the appropriate component and skip it here (no segging done yet)
    if (vb->check_for_gc && sam_seg_is_gc_line (vb, dl, alignment, after - alignment, true)) {
        vb->debug_line_hash_skip = true;
        memset (dl, 0, sizeof (ZipDataLineSAM));
        goto done;
    }

    // case: biopsy (only arrives here in MAIN VBs if gencomp) biopsy_line: 
    // we just need to pass sam_seg_is_gc_line and we're done
    if ((flag.biopsy && !segconf.running) ||
        (flag.has_biopsy_line && sam_seg_test_biopsy_line (VB, alignment, block_size + 4)) )
        goto done;  

    sam_cigar_binary_to_textual (vb, n_cigar_op, B1ST(BamCigarOp, vb->binary_cigar), // binary_cigar and not "cigar", as the latter is mis-aligned 
                                 &vb->textual_cigar); // re-write BAM format CIGAR as SAM textual format in vb->textual_cigar

    // SEQ - calculate diff vs. reference (denovo or loaded)
    ASSERT (dl->SEQ.len == l_seq || (vb->textual_cigar.len == 1 && *B1STc(vb->textual_cigar) == '*') || !l_seq, 
            "seq_len implied by CIGAR=%s is %u, but actual SEQ length is %u, SEQ=%.*s", 
            B1STc(vb->textual_cigar), dl->SEQ.len, l_seq, l_seq, B1STc(vb->textual_seq));

    // convert BAM QUAL format to SAM
    if (!l_seq || !bam_rewrite_qual ((uint8_t *)qual, l_seq)) { // add 33 to Phred scores to make them ASCII
        bam_set_missing_qual_type (vb, (uint8_t *)qual, l_seq);
        vb->qual_missing = dl->no_qual = true;
    }

    if (has(NM_i)) 
        dl->NM_len = sam_seg_get_aux_int (vb, vb->idx_NM_i, &dl->NM, true, MIN_NM_i, MAX_NM_i, HARD_FAIL);

    // set dl->AS needed by sam_seg_prim_add_sag (in PRIM) and several fields that delta against it
    if (has(AS_i))
        sam_seg_get_aux_int (vb, vb->idx_AS_i, &dl->AS, true, MIN_AS_i, MAX_AS_i, HARD_FAIL);

    if (!IS_MAIN(vb)) 
        sam_seg_sag_stuff (vb, dl, STRb(vb->textual_cigar), B1STc(vb->textual_seq), true);

    // seg QNAME first, as it will find the buddy
    sam_seg_QNAME (vb, dl, read_name, l_read_name-1, 2); // QNAME. account for \0 and l_read_name

    bam_seg_ref_id (vb, dl, SAM_RNAME, vb->chrom_node_index, -1); // ref_id (RNAME)

    // note: pos can have a value even if ref_id=-1 (RNAME="*") - this happens if a SAM with a RNAME that is not in the header is converted to BAM with samtools
    sam_seg_POS (vb, dl, prev_line_chrom, sizeof (uint32_t)); // POS
    
    if (vb->chrom_node_index >= 0) sam_seg_verify_RNAME (vb);

    sam_seg_MAPQ (vb, dl, sizeof (uint8_t));

    sam_seg_FLAG (vb, dl, sizeof (uint16_t));
    
    bam_seg_ref_id (vb, dl, SAM_RNEXT, dl->RNEXT, vb->chrom_node_index); // RNEXT

    sam_seg_PNEXT (vb, dl, 0, 0, next_pos, sizeof (uint32_t));

    // Segment BIN after we've gathered bin, flags, pos and vb->ref_confumed (and before sam_seg_SEQ which ruins vb->ref_consumed)
    bam_seg_BIN (vb, dl, bin, true);

    sam_seg_init_bisulfite (vb, dl); 

    // we search forward for MD:Z now, XG:Z as we will need it for SEQ if it exists
    if (has_MD)
        sam_seg_MD_Z_analyze (vb, dl, STRauxZ(MD_Z,true), dl->POS);
        
    sam_seg_SEQ (vb, dl, STRb(vb->textual_seq), (l_seq+1)/2 + sizeof (uint32_t) /* account for l_seq and seq fields */);

    // finally we can segment the textual CIGAR now (including if n_cigar_op=0)
    sam_seg_CIGAR (vb, dl, vb->textual_cigar.len32, STRb(vb->textual_seq), qual, l_seq, 
                           ((uint32_t)n_cigar_op * sizeof (uint32_t) /* cigar */ + sizeof (uint16_t) /* n_cigar_op */));

    // QUAL. note: can only be called after sam_seg_CIGAR updates SEQ.len
    if (!vb->qual_missing) // case we have both SEQ and QUAL
        sam_seg_QUAL (vb, dl, qual, l_seq, l_seq /* account for qual field */ );

    else { // cases 1. were both SEQ and QUAL are '*' (seq_len=0) and 2. SEQ exists, QUAL not (bam_rewrite_qual returns false)
        dl->QUAL = (TxtWord){ .index=Ltxt, .len = 1 }; // a '*' was placed after txt_data in bam_seg_initialize 

        sam_seg_QUAL (vb, dl, alignment, 1, l_seq /* account of l_seq 0xff */);
    }

    // AUX fields - up to MAX_FIELDS of them
    sam_seg_aux_all (vb, dl);
    
    // TLEN - must be after AUX as we might need data from MC:Z
    sam_seg_TLEN (vb, dl, 0, 0, tlen, vb->chrom_node_index == dl->RNEXT); // TLEN

    if (IS_PRIM(vb)) {
        if      (IS_SAG_NH)   sam_seg_prim_add_sag_NH (vb, dl, dl->NH);
        else if (IS_SAG_CC)   sam_seg_prim_add_sag_CC (vb, dl, dl->NH);
        else if (IS_SAG_FLAG) sam_seg_prim_add_sag (vb, dl, 0, true);
        else if (IS_SAG_SOLO) sam_seg_prim_add_sag_SOLO (vb, dl);
    }

    if (dl->SEQ.len > vb->longest_seq_len) vb->longest_seq_len = dl->SEQ.len;

    if (segconf.running)
        segconf.est_segconf_sam_size += bam_segconf_get_transated_sam_line_len (vb, dl, tlen);

done: {}
    buf_free (vb->textual_cigar);
    buf_free (vb->textual_seq);

    return after;
}
