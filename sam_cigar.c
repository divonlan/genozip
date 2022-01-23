// ------------------------------------------------------------------
//   sam_cigar.c
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

// a module for handling CIGAR and MC:Z

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "reference.h"
#include "segconf.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "codec.h"
#include "profiler.h"

static const bool cigar_valid_op[256] = { ['M']=true, ['I']=true, ['D']=true, ['N']=true, ['S']=true, ['H']=true, ['P']=true, ['=']=true, ['X']=true }; 

//---------
// Shared
//---------

void sam_cigar_binary_to_textual (VBlockSAM *vb, uint16_t n_cigar_op, const uint32_t *cigar, Buffer *textual_cigar /* out */)
{
    if (!n_cigar_op) {
        buf_alloc (vb, textual_cigar, 2, 0, char, 100, textual_cigar->name ? NULL : "textual_cigar");
        NEXTENT (char, *textual_cigar) = '*';
        goto finish;
    }

    // calculate length
    unsigned len=0;
    for (uint16_t i=0; i < n_cigar_op; i++) {
        uint32_t op_len = LTEN32 (cigar[i]) >> 4; // maximum is 268,435,455
        if      (op_len < 10)       len += 2; // 1 for the op, 1 for the number
        else if (op_len < 100)      len += 3; 
        else if (op_len < 1000)     len += 4; 
        else if (op_len < 10000)    len += 5; 
        else if (op_len < 100000)   len += 6; 
        else if (op_len < 1000000)  len += 7;
        else if (op_len < 10000000) len += 8;
        else                        len += 9;
    }

    buf_alloc (vb, textual_cigar, len + 1 /* for \0 */, 100, char, 0, textual_cigar->name ? NULL : "textual_cigar");

    for (uint16_t i=0; i < n_cigar_op; i++) {
        uint32_t subcigar = LTEN32 (cigar[i]);
        uint32_t op_len = subcigar >> 4;

        textual_cigar->len += str_int (op_len, AFTERENT (char, *textual_cigar));

        static const char op_table[16] = "MIDNSHP=X";
        NEXTENT (char, *textual_cigar) = op_table[subcigar & 0xf];
    }

finish:
    *AFTERENT (char, *textual_cigar) = 0; // nul terminate
    
    if (command == ZIP)
        vb->last_cigar = FIRSTENT (char, *textual_cigar);
}

// calculate the expected length of SEQ and QUAL from the CIGAR string
// A CIGAR looks something like: "109S19M23S", See: https://samtools.github.io/hts-specs/SAMv1.pdf 
void sam_cigar_analyze (VBlockSAMP vb, STRp(cigar),  unsigned *seq_consumed)
{
    *seq_consumed            = 0;
    vb->ref_consumed         = 0;
    vb->ref_and_seq_consumed = 0;
    vb->soft_clip            = 0;
    vb->mismatch_bases       = 0;

    ASSERT (cigar[0] != '*' || cigar_len == 1, "Invalid CIGAR: %.*s", cigar_len, cigar); // a CIGAR start with '*' must have 1 character

    if (command == PIZ) 
        vb->textual_cigar.len = 0;

    // ZIP case: if the CIGAR is "*", later sam_cigar_seg_textual uses the length from SEQ and store it as eg "151*". 
    // In PIZ it will be eg "151*" or "1*" if both SEQ and QUAL are "*", so this condition is always false
    if (cigar[0] == '*') return;

    // PIZ case: CIGAR string starts with '-' (indicating missing SEQ) - just skip the '-' for now
    if (*cigar == '-') {
        cigar++;
        cigar_len--;
    }

    // store original textual CIGAR for use of sam_piz_special_MD, as in BAM it will be translated ; also cigar might point to buddy data in ctx->per_line - ctx->per_line might be realloced as we store this line's CIGAR in it 
    if (command == PIZ) {
        buf_add_more (VB, &vb->textual_cigar, STRa(cigar), "textual_cigar");
        *AFTERENT (char, vb->textual_cigar) = 0; // nul-terminate (buf_add_more allocated space for it)
    }

    // if we're reconstructing a BAM, we will create the BAM cigar data in binary_cigar. 
    bool bam_piz = (command == PIZ && flag.out_dt == DT_BAM);
    if (bam_piz) buf_alloc (vb, &vb->binary_cigar, 0, cigar_len/2 /* max possible n_cigar_op */, uint32_t, 2, "binary_cigar");

    unsigned n=0;
    for (unsigned i=0; i < cigar_len; i++) {

        char c = cigar[i];
        char lookup = cigar_lookup_sam[(uint8_t)c];

        ASSINP (lookup, "Invalid CIGAR in %s: invalid operation %c. CIGAR=%.*s", txt_name, cigar[i], cigar_len, cigar);
        lookup &= 0x0f; // remove validity bit

        if (lookup == CIGAR_DIGIT) 
            n = n*10 + (c - '0');
        
        else {
            ASSINP (n, "Invalid CIGAR in %s: operation %c not preceded by a number. CIGAR=%.*s", 
                    txt_name, c, cigar_len, cigar);
            
            if ((lookup & CIGAR_CONSUMES_QUERY))    *seq_consumed += n;
            if ((lookup & CIGAR_CONSUMES_REFERENCE)) vb->ref_consumed += n;
            if ((lookup & CIGAR_CONSUMES_QUERY) && (lookup & CIGAR_CONSUMES_REFERENCE)) vb->ref_and_seq_consumed += n;
            if (c == 'I' || c == 'D') vb->mismatch_bases += n;
            else if (c == 'S') vb->soft_clip += n;

            // note: piz: in case of eg "151*" - *seq_consumed will be updated to the length, but binary_cigar will be empty
            if (bam_piz && c != '*') { 
                // convert character CIGAR op to BAM cigar field op  "MIDNSHP=X" -> 012345678 
                static const uint8_t cigar_char_to_op[256] = { ['M']=0, ['I']=1, ['D']=2, ['N']=3, ['S']=4, 
                                                               ['H']=5, ['P']=6, ['=']=7, ['X']=8           }; 
                uint32_t bam_cigar = (n << 4) | (uint32_t)cigar_char_to_op[(uint8_t)c];
                NEXTENT (uint32_t, vb->binary_cigar) = LTEN32 (bam_cigar);
            }
            n = 0;
        }
    }                          

    if (command == ZIP)
        DATA_LINE (vb->line_i)->ref_consumed = vb->ref_consumed; // consumed by sam_seg_predict_TLEN 

    // PIZ: we store ref_consumed in ctx->piz_ctx_specific_buf because ctx->history is already taken for storing the CIGAR string
    else {
        uint64_t line_i = vb->line_i - vb->first_line;
        if (!line_i)
            buf_alloc_zero (vb, &CTX(SAM_CIGAR)->piz_ctx_specific_buf, 0, vb->lines.len, uint32_t, 0, "piz_ctx_specific_buf"); // initialize to exactly one per line.

        *ENT (uint32_t, CTX(SAM_CIGAR)->piz_ctx_specific_buf, line_i) = vb->ref_consumed;
    }
    
    ASSINP (!n, "Invalid CIGAR in %s: expecting it to end with an operation character. CIGAR=%.*s", 
            txt_name, cigar_len, cigar);

    ASSINP (!seq_consumed || *seq_consumed, "Invalid CIGAR in %s: CIGAR implies 0-length SEQ. CIGAR=%.*s", txt_name, cigar_len, cigar);
}

bool sam_cigar_is_valid (STRp(cigar))
{
    unsigned i=0;
    while (i < cigar_len) {

        unsigned num_digits=0;
        for (; i < cigar_len && IS_DIGIT(cigar[i]) ; i++) num_digits++;

        if (!num_digits) return false;

        if (i == cigar_len || !cigar_valid_op[(int)cigar[i++]])
            return false;
    }
    return true;
}

// reverses a CIGAR, eg "40S111M"->"111M40S". returns false if not a valid CIGAR string.
bool sam_cigar_reverse (char *dst, STRp(cigar))
{
    const char *c = &cigar[cigar_len-1];
    while (c >= cigar) {

        char cigar_op = *(c--);
        if (!cigar_valid_op[(int)cigar_op]) return false;

        unsigned num_digits=0;
        for (; c >= cigar && IS_DIGIT(*c) ; c--) num_digits++;
        
        if (!num_digits) return false;

        memcpy (dst, c+1, num_digits);
        dst[num_digits] = cigar_op;
        
        dst += num_digits + 1;
    }

    return true;
}

//---------
// SEG
//---------

void sam_cigar_seg_textual (VBlockSAM *vb, ZipDataLineSAM *dl, unsigned last_cigar_len, STRp(seq_data), STRp(qual_data))
{
    START_TIMER
    
    bool qual_is_available = (qual_data_len != 1 || *qual_data != '*');
    bool seq_is_available  = (seq_data_len  != 1 || *seq_data  != '*');

    ASSSEG (!(seq_is_available && *seq_data=='*'), seq_data, "seq_data=%.*s (seq_len=%u), but expecting a missing seq to be \"*\" only (1 character)", 
            seq_data_len, seq_data, seq_data_len);

    char cigar_snip[last_cigar_len + 50];
    cigar_snip[0] = SNIP_SPECIAL;
    cigar_snip[1] = SAM_SPECIAL_CIGAR;
    unsigned cigar_snip_len=2;

    // case: SEQ is "*" - we add a '-' to the CIGAR
    if (!seq_is_available) cigar_snip[cigar_snip_len++] = '-';

    // case: CIGAR is "*" - we get the dl->seq_len directly from SEQ or QUAL, and add the length to CIGAR eg "151*"
    if (!dl->seq_len) { // CIGAR is not available
        ASSSEG (!seq_data_len || !qual_is_available || seq_data_len==dl->QUAL.len, seq_data,
                "Bad line: SEQ length is %u, QUAL length is %u, unexpectedly differ. SEQ=%.*s QUAL=%.*s", 
                seq_data_len, dl->QUAL.len, seq_data_len, seq_data, dl->QUAL.len, qual_data);    

        dl->seq_len = MAX_(seq_data_len, dl->QUAL.len); // one or both might be not available and hence =1

        cigar_snip_len += str_int (dl->seq_len, &cigar_snip[cigar_snip_len]);
    } 
    else { // CIGAR is available - just check the seq and qual lengths
        ASSSEG (!seq_is_available || seq_data_len == dl->seq_len, seq_data,
                "Bad line: according to CIGAR, expecting SEQ length to be %u but it is %u. SEQ=%.*s", 
                dl->seq_len, seq_data_len, seq_data_len, seq_data);

        ASSSEG (!qual_is_available || qual_data_len == dl->seq_len, qual_data,
                "Bad line: according to CIGAR, expecting QUAL length to be %u but it is %u. QUAL=%.*s", 
                dl->seq_len, dl->QUAL.len, dl->QUAL.len, qual_data);    
    }

    // case: we buddy non-trival CIGARs with MC:Z. We don't buddy eg "151M" bc this will add rather than reduce entropy.
    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // an invalid pointer if buddy_line_i is -1
    bool seg_done=false;

    if (last_cigar_len > 4 && vb->buddy_line_i != -1 && segconf.has[OPTION_MC_Z] && !segconf.running && 
        cigar_snip_len == 2 && // we don't buddy if CIGAR or SEQ are "*"
        buddy_dl->MC.len == last_cigar_len && 
        !memcmp (vb->last_cigar, ENT (char, vb->txt_data, buddy_dl->MC.index), last_cigar_len)) 

        cigar_snip[cigar_snip_len++] = COPY_BUDDY; // always at cigar_snip[2]

    // case: normal cigar (not buddy copy) 
    else {
        memcpy (&cigar_snip[cigar_snip_len], vb->last_cigar, last_cigar_len);
        cigar_snip_len += last_cigar_len;

        if (last_cigar_len > 7){
            seg_add_to_local_text (VB, CTX(SAM_CIGAR), STRa(cigar_snip), true, last_cigar_len+1);
            seg_done = true;
        }
    }
    
    // store the CIGAR in DataLine for use by a buddy MC:Z and SA:Z
    dl->CIGAR = (TxtWord){ .index = ENTNUM (vb->txt_data, vb->last_cigar), .len = last_cigar_len }; // in SAM (but not BAM) vb->last_cigar points into txt_data

    if (segconf.running) {
        segconf.sam_cigar_len += last_cigar_len;
        segconf.sam_seq_len   += seq_data_len;
    }

    if (!seg_done)
        seg_by_did_i (VB, STRa(cigar_snip), SAM_CIGAR, last_cigar_len+1); // +1 for \t

    COPY_TIMER(sam_cigar_seg);
}


void sam_cigar_seg_binary (VBlockSAM *vb, ZipDataLineSAM *dl, uint32_t l_seq, uint32_t n_cigar_op)
{
    START_TIMER

    char cigar_snip[vb->textual_cigar.len + 20];
    cigar_snip[0] = SNIP_SPECIAL;
    cigar_snip[1] = SAM_SPECIAL_CIGAR;
    unsigned cigar_snip_len=2;

    // case: we have no sequence - we add a '-' to the CIGAR
    if (!l_seq) cigar_snip[cigar_snip_len++] = '-';

    // case: CIGAR is "*" - we get the dl->seq_len directly from SEQ or QUAL, and add the length to CIGAR eg "151*"
    if (!dl->seq_len) { // CIGAR is not available
        dl->seq_len = l_seq;
        cigar_snip_len += str_int (MAX_(dl->seq_len, 1), &cigar_snip[cigar_snip_len]); // if seq_len=0, then we add "1" because we have the * in seq and qual (and consistency with sam_cigar_seg_textual)
    }
    
    // case: we buddy non-trival CIGARs with MC:Z. We don't buddy eg "151M" bc this will add rather than reduce entropy.
    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // an invalid pointer if buddy_line_i is -1
    bool seg_done=false;
    unsigned add_bytes = n_cigar_op * sizeof (uint32_t) /* cigar */ + sizeof (uint16_t) /* n_cigar_op */;

    if (vb->textual_cigar.len > 4 && vb->buddy_line_i != -1 && segconf.has[OPTION_MC_Z] && !segconf.running && 
        buddy_dl->MC.len == vb->textual_cigar.len && 
        cigar_snip_len == 2 && // we don't buddy if CIGAR or SEQ are "*"
        !memcmp (vb->textual_cigar.data, ENT (char, vb->txt_data, buddy_dl->MC.index), vb->textual_cigar.len)) 

        cigar_snip[cigar_snip_len++] = COPY_BUDDY; // always at cigar_snip[2]

    // case: normal cigar (not buddy copy) 
    else {
        memcpy (&cigar_snip[cigar_snip_len], vb->textual_cigar.data, vb->textual_cigar.len);
        cigar_snip_len += vb->textual_cigar.len;

        if (vb->textual_cigar.len > 7) {
            seg_add_to_local_text (VB, CTX(SAM_CIGAR), STRa(cigar_snip), true, add_bytes);
            seg_done = true;
        }
    }

    // store a copy of the CIGAR in buddy_textual_cigars for use by a buddy MC:Z
    if (segconf.has[OPTION_MC_Z] && !segconf.running) {
        dl->CIGAR =(TxtWord){ .index = vb->buddy_textual_cigars.len, .len = vb->textual_cigar.len }; // in BAM dl->CIGAR points into buddy_textual_cigars
        buf_add_buf (VB, &vb->buddy_textual_cigars, &vb->textual_cigar, char, "buddy_textual_cigars");
    }

    if (segconf.running) {
        segconf.sam_cigar_len += cigar_snip_len;
        segconf.sam_seq_len   += dl->seq_len;
    }

    if (!seg_done)
        seg_by_did_i (VB, STRa(cigar_snip), SAM_CIGAR, add_bytes);

    COPY_TIMER(sam_cigar_seg);
}

unsigned sam_cigar_get_MC_ref_consumed (STRp(mc))
{
    // get ref_and_seq_consumed
    unsigned n=0;
    unsigned ref_and_seq_consumed=0;
    for (unsigned i=0; i < mc_len; i++) {

        char c = mc[i];
        char lookup = cigar_lookup_sam[(uint8_t)c];
        if (!lookup) return 0; // invalid CIGAR - unrecognized character

        lookup &= 0x0f; // remove validity bit

        if (lookup == CIGAR_DIGIT) 
            n = n*10 + (c - '0');
        
        else {
            if (!n) return 0; // invalid CIGAR - no number before op

            if ((lookup & CIGAR_CONSUMES_REFERENCE)) 
                ref_and_seq_consumed += n;

            n=0;
        }
    }
    return ref_and_seq_consumed;
}

// MC:Z "CIGAR string for mate/next segment" (https://samtools.github.io/hts-specs/SAMtags.pdf)
void sam_cigar_seg_MC (VBlockSAM *vb, ZipDataLineSAM *dl, STRp(mc), unsigned add_bytes)
{
    segconf_set_has (OPTION_MC_Z);

    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // an invalid pointer if buddy_line_i is -1
    Buffer *buddy_cigar_buf = IS_BAM ? &vb->buddy_textual_cigars : &vb->txt_data; // buddy_dl->MC points into this buffer

    // we buddy non-trival CIGARs. We don't buddy eg "151M" bc this will add rather than reduce entropy.
    if (mc_len > 4 && !segconf.running && vb->buddy_line_i != -1 && 
        mc[0] != '*' &&
        buddy_dl->CIGAR.snip_len == mc_len && !memcmp (mc, ENT (char, *buddy_cigar_buf, buddy_dl->CIGAR.char_index), mc_len)) 

        seg_by_did_i (VB, STRa(MC_buddy_snip), OPTION_MC_Z, add_bytes); // copy MC from earlier-line buddy CIGAR
    
    else if (mc_len > 7) 
        seg_add_to_local_text (VB, CTX(OPTION_MC_Z), STRa(mc), true, add_bytes);
    else
        seg_by_did_i (VB, STRa(mc), OPTION_MC_Z, add_bytes);    

    dl->MC = WORD_IN_TXT_DATA(mc); 

    ctx_set_last_value (VB, CTX(OPTION_MC_Z), (ValueType){ .i = sam_cigar_get_MC_ref_consumed (STRa(mc)) } );
}

//---------
// PIZ
//---------

// CIGAR - calculate vb->seq_len from the CIGAR string, and if original CIGAR was "*" - recover it
SPECIAL_RECONSTRUCTOR (sam_cigar_special_CIGAR)
{
    VBlockSAMP vb_sam = (VBlockSAMP)vb;
    const uint16_t sam_flag = vb->last_int(SAM_FLAG);

    // case: we copy the snip from buddy MC:Z
    if (snip[0] == COPY_BUDDY) 
        reconstruct_from_buddy_get_textual_snip (vb, CTX (OPTION_MC_Z), pSTRa(snip));

    // calculate seq_len (= l_seq, unless l_seq=0), ref_consumed and (if bam) vb->textual_cigar and vb->binary_cigar
    sam_cigar_analyze (vb_sam, STRa(snip), &vb->seq_len); 

    if ((flag.out_dt == DT_SAM || (flag.out_dt == DT_FASTQ && flag.extended_translation)) 
    &&  reconstruct) {
        if (snip[snip_len-1] == '*') // eg "151*" - zip added the "151" to indicate seq_len - we don't reconstruct it, just the '*'
            RECONSTRUCT1 ('*');
        
        else if (snip[0] == '-') // eg "-151M" or "-151*" - zip added the "-" to indicate a '*' SEQ field - we don't reconstruct it
            RECONSTRUCT (snip + 1, snip_len - 1);

        else
            RECONSTRUCT (snip, snip_len);    
    }

    // BAM - output vb->binary_cigar generated in sam_cigar_analyze
    else if (flag.out_dt == DT_BAM) {
        // now we have the info needed to reconstruct bin, l_read_name, n_cigar_op and l_seq
        BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)ENT (char, vb->txt_data, vb->line_start);
        alignment->l_read_name = AFTERENT (char, vb->txt_data) - alignment->read_name;
        alignment->n_cigar_op  = LTEN16 (vb_sam->binary_cigar.len);
        alignment->l_seq       = (snip[0] == '-') ? 0 : LTEN32 (vb->seq_len);

        RECONSTRUCT (vb_sam->binary_cigar.data, vb_sam->binary_cigar.len * sizeof (uint32_t));
        buf_free (&vb_sam->binary_cigar);

        // if BIN is SAM_SPECIAL_BIN, inst.semaphone is set by bam_piz_special_BIN - a signal to us to calculate
        ContextP sam_bam_bin_ctx = CTX(SAM_BAM_BIN);
        if (sam_bam_bin_ctx->semaphore) {
            sam_bam_bin_ctx->semaphore = false;

            PosType pos = CTX(SAM_POS)->last_value.i;
            bool segment_unmapped = (sam_flag & SAM_FLAG_UNMAPPED);
            PosType last_pos = segment_unmapped ? pos : (pos + vb_sam->ref_consumed - 1);
            
            uint16_t bin = bam_reg2bin (pos, last_pos); // zero-based, half-closed half-open [start,end)
            alignment->bin = LTEN16 (bin); // override the -1 previously set by the translator
        }
    }
    
    else if (flag.out_dt == DT_FASTQ) {
        // only analyze, but don't reconstruct CIGAR in FASTQ
    }

    vb_sam->last_cigar = snip;

    return false; // no new value
}   

// copy from buddy CIGAR to MC. If reconstructing to BAM, we convert the binary CIGAR to textual MC:Z.*

SPECIAL_RECONSTRUCTOR (sam_piz_special_COPY_BUDDY_MC)
{
    if (!reconstruct) return false;

    VBlockSAMP sam_vb = (VBlockSAMP)vb;

    // fall back to normal COPY_BUDDY in case of SAM
    if (flag.out_dt != DT_BAM) return reconstruct_from_buddy (vb, ctx, STRa(snip), reconstruct, new_value); // SAM or FASTQ output

    // get CIGAR field value previously reconstructed in BAM binary format
    STR(bam_cigar);
    reconstruct_from_buddy_get_textual_snip (vb, CTX (SAM_CIGAR), pSTRa(bam_cigar));
    
    // convert binary CIGAR to textual MC:Z
    uint32_t n_cigar_op = bam_cigar_len / sizeof (uint32_t);
    sam_cigar_binary_to_textual (sam_vb, n_cigar_op, (uint32_t*)bam_cigar, &vb->txt_data);

    return false; // no new value 
}

// invoked from TOP2FQ (but not TOP2FQEX, bc it reconstructs AUX) to consume MC if it exists in this line, in case this line is
// a buddy line of a future line in which case this MC will be copied to the future line's CIGAR 
SPECIAL_RECONSTRUCTOR (sam_piz_special_CONSUME_MC_Z)
{
    ContextP mc_ctx = CTX(OPTION_MC_Z);
    if (!mc_ctx->flags.store_per_line) goto done; // MC is not buddied
    
    ContextP opt_ctx = CTX(SAM_AUX);
    WordIndex opt_word_index = WORD_INDEX_NONE;

    // get AUX container
    snip_len=0;
    if (opt_ctx->b250.len ||
        (!opt_ctx->b250.len && !opt_ctx->local.len && opt_ctx->dict.len)) {  // all_the_same case - no b250 or local, but have dict      
        opt_word_index = LOAD_SNIP(opt_ctx->did_i); // note: if we have no b250, local but have dict, this will be word_index=0 (see ctx_get_next_snip)

        if (snip_len==1 && *snip == SNIP_LOOKUP)
            snip_len=0;
    }
    
    // case: a singleton (SNIP_LOOKUP) or all data is local
    if (!snip_len) 
        LOAD_SNIP_FROM_LOCAL (opt_ctx);

    if (!snip_len || snip[0] != SNIP_CONTAINER) goto done; // not a container

    ContainerP opt_con = container_retrieve (vb, opt_ctx, opt_word_index, snip+1, snip_len-1, 0, 0);

    // check if this line has an optional tag MC:Z
    bool found = false;
    for (uint32_t item_i=0; item_i < con_nitems (*opt_con); item_i++)
        if (opt_con->items[item_i].dict_id.num == _OPTION_MC_Z) {
            found = true;
            break;
        }
    if (!found) goto done; // AUX has no MC:Z tag

    // store MC:Z in history to be used by future line CIGAR. Note: since MC:Z is not reconstructed, ctx->history will point
    // to either dict or local (depending on where this snip originates) rather than the normal txt_data.
    snip_len = 0;

    // case: MC:Z is in dict (refered to from a b250)
    if (mc_ctx->b250.len ||
        (!mc_ctx->b250.len && !mc_ctx->local.len && mc_ctx->dict.len)) {  // all_the_same case - no b250 or local, but have dict      
        LOAD_SNIP(mc_ctx->did_i); // note: if we have no b250, local but have dict, this will be word_index=0 (see ctx_get_next_snip)

        if (snip_len==1 && *snip == SNIP_LOOKUP)
            snip_len=0;

        else 
            *ENT (HistoryWord, mc_ctx->history, vb->line_i - vb->first_line) = 
                (HistoryWord){ .char_index = ENTNUM (mc_ctx->dict, snip), .snip_len = snip_len, .lookup = LookupDict };
    }
    
    // case: MC:Z is in local
    if (!snip_len) {
        uint32_t char_index = LOAD_SNIP_FROM_LOCAL (mc_ctx);
        *ENT (HistoryWord, mc_ctx->history, vb->line_i - vb->first_line) = 
            (HistoryWord){ .char_index = char_index, .snip_len = snip_len, .lookup = LookupLocal }; 
    }

done:
    return false; // no new value 
}
