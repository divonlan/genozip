// ------------------------------------------------------------------
//   sam_piz.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <math.h>
#include "sam_private.h"
#include "seg.h"
#include "context.h"
#include "piz.h"
#include "reconstruct.h"
#include "strings.h"
#include "dict_id.h"
#include "codec.h" // must be included before reference.h
#include "reference.h"
#include "regions.h"
#include "aligner.h"
#include "file.h"
#include "container.h"

// returns true if section is to be skipped reading / uncompressing
bool sam_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    if (!vb) return false; // we don't skip reading any SEC_DICT sections

    return false;
}

// PIZ: SEQ reconstruction 
void sam_reconstruct_seq (VBlock *vb_, Context *bitmap_ctx, const char *unused, unsigned unused2)
{
#define ROUNDUP_TO_NEAREST_4(x) ((uint32_t)(x) + 3) & ~((uint32_t)0x3)

    VBlockSAMP vb = (VBlockSAMP)vb_;
    ASSERT0 (bitmap_ctx && bitmap_ctx->did_i == SAM_SQBITMAP, "Error in sam_reconstruct_seq: context is not SAM_SQBITMAP");

    if (piz_is_skip_section (vb, SEC_LOCAL, bitmap_ctx->dict_id)) return; // if case we need to skip the SEQ field (for the entire file)

    Context *nonref_ctx      = &vb->contexts[SAM_NONREF];
    const char *nonref       = ENT (const char, nonref_ctx->local, nonref_ctx->next_local); // possibly, this VB has no nonref (i.e. everything is ref), in which case nonref would be an invalid pointer. That's ok, as it will not be accessed.
    const char *nonref_start = nonref;
    unsigned subcigar_len    = 0;
    char cigar_op            = 0;
    const PosType pos        = vb->contexts[SAM_POS].last_value.i;
    const Range *range       = NULL;
    unsigned seq_consumed=0, ref_consumed=0;

    // case: unaligned sequence - pos is 0 
    if (!pos || (vb->chrom_name_len==1 && vb->chrom_name[0]=='*')) {
        // case: compressed with a reference, using our aligner
        if (z_file->z_flags.aligner) {
            aligner_reconstruct_seq ((VBlockP)vb, bitmap_ctx, vb->seq_len, false);
            nonref_ctx->next_local = ROUNDUP_TO_NEAREST_4 (nonref_ctx->next_local);
        }
        // case: no reference was used - in this case, the sequence is not encoded in the bitmap at all. we just copy it from NONREF
        else {
            RECONSTRUCT (nonref, vb->seq_len); 
            nonref_ctx->next_local += ROUNDUP_TO_NEAREST_4 (vb->seq_len);
        }
        return;
    }

    // case: missing sequence - sequence is '*' (which zip marked with a '-' in the cigar) - just reconstruct the '*'
    if (*vb->last_cigar == '-') {
        RECONSTRUCT1 ('*');
        vb->last_cigar++; // skip the '-' so it doesn't affect a subsequent E2 on the same line
        return;
    }

    const char *next_cigar = vb->last_cigar; // don't change vb->last_cigar as we may still need it, eg if we have an E2 optional field
    range = ref_piz_get_range (vb_, pos, vb->ref_consumed);
    
    while (seq_consumed < vb->seq_len || ref_consumed < vb->ref_consumed) {
        
        if (!subcigar_len) {
            subcigar_len = strtod (next_cigar, (char **)&next_cigar); // get number and advance next_cigar
        
            cigar_op = cigar_lookup_sam[(uint8_t)*(next_cigar++)];
            ASSERT (cigar_op, "Error in sam_reconstruct_seq: Invalid CIGAR op while reconstructing line %u: '%c' (ASCII %u)", vb->line_i, *(next_cigar-1), *(next_cigar-1));
            cigar_op &= 0x0f; // remove validity bit
        }

        if (cigar_op & CIGAR_CONSUMES_QUERY) {

            if ((cigar_op & CIGAR_CONSUMES_REFERENCE) && NEXTLOCALBIT (bitmap_ctx)) /* copy from reference */ {

                if (!vb->dont_show_curr_line) { // note: if this line is excluded with --regions, then the reference section covering it might not be loaded
                    uint32_t idx = (pos - range->first_pos) + ref_consumed ;

                    if (!ref_is_nucleotide_set (range, idx)) { 
                        ref_print_is_set (range, pos + ref_consumed, stderr);
                        ABORT ("Error in sam_reconstruct_seq: while reconstructing line %u (vb_i=%u: last_txt_line=%u num_lines=%u): reference is not set: chrom=%u \"%.*s\" pos=%"PRId64" range=[%"PRId64"-%"PRId64"]"
                               " (cigar=%s seq_start_pos=%"PRId64" ref_consumed=%u seq_consumed=%u)",
                               vb->line_i, vb->vblock_i, (uint32_t)(vb->first_line + vb->lines.len - 1), (uint32_t)vb->lines.len, range->chrom, range->chrom_name_len, range->chrom_name, pos + ref_consumed, range->first_pos, range->last_pos, vb->last_cigar, pos, ref_consumed, seq_consumed);
                    }

                    char ref = ref_get_nucleotide (range, idx);
                    RECONSTRUCT1 (ref); 
                }
            }
            else 
                RECONSTRUCT1 (*nonref++);

            seq_consumed++;
        }

        if (cigar_op & CIGAR_CONSUMES_REFERENCE) 
            ref_consumed++;

        subcigar_len--;
    }

    ASSERT (seq_consumed == vb->seq_len,      "Error in sam_reconstruct_seq: expecting seq_consumed(%u) == vb->seq_len(%u)", seq_consumed, vb->seq_len);
    ASSERT (ref_consumed == vb->ref_consumed, "Error in sam_reconstruct_seq: expecting ref_consumed(%u) == vb->ref_consumed(%u)", ref_consumed, vb->ref_consumed);

    bitmap_ctx->last_value.i = bitmap_ctx->next_local; // for SEQ, we use last_value for storing the beginning of the sequence
    
    nonref_ctx->next_local += ROUNDUP_TO_NEAREST_4 (nonref - nonref_start);
}

// CIGAR - calculate vb->seq_len from the CIGAR string, and if original CIGAR was "*" - recover it
SPECIAL_RECONSTRUCTOR (sam_piz_special_CIGAR)
{
    VBlockSAMP vb_sam = (VBlockSAMP)vb;

    // calculate seq_len (= l_seq, unless l_seq=0), ref_consumed and (if bam) vb->textual_cigar
    sam_analyze_cigar (vb_sam, snip, snip_len, &vb->seq_len, &vb_sam->ref_consumed, NULL); 

    if (flag.out_dt == DT_SAM) {
        if (snip[snip_len-1] == '*') // eg "151*" - zip added the "151" to indicate seq_len - we don't reconstruct it, just the '*'
            RECONSTRUCT1 ('*');
        
        else if (snip[0] == '-') // eg "-151M" or "-151*" - zip added the "-" to indicate a '*' SEQ field - we don't reconstruct it
            RECONSTRUCT (snip + 1, snip_len - 1)

        else
            RECONSTRUCT (snip, snip_len);    
    }

    // BAM - output vb->textual_cigar generated in sam_analyze_cigar
    else if (flag.out_dt == DT_BAM) {
        // now we have the info needed to reconstruct bin, l_read_name, n_cigar_op and l_seq
        BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)ENT (char, vb->txt_data, vb->line_start);
        alignment->l_read_name = AFTERENT (char, vb->txt_data) - alignment->read_name;
        alignment->n_cigar_op  = LTEN16 (vb_sam->textual_cigar.len);
        alignment->l_seq       = (snip[0] == '-') ? 0 : LTEN32 (vb->seq_len);

        RECONSTRUCT (vb_sam->textual_cigar.data, vb_sam->textual_cigar.len * sizeof (uint32_t));
        buf_free (&vb_sam->textual_cigar);

        // if BIN is SAM_SPECIAL_BIN, inst.semaphone is set by bam_piz_special_BIN - a signal to us to calculate
        if (vb->contexts[SAM_BAM_BIN].semaphore) {
            vb->contexts[SAM_BAM_BIN].semaphore = false;

            uint16_t flag = vb->contexts[SAM_FLAG].last_value.i;
            PosType pos   = vb->contexts[SAM_POS ].last_value.i;
            bool segment_unmapped = (flag & 0x4);
            PosType last_pos = segment_unmapped ? pos : (pos + vb_sam->ref_consumed - 1);
            
            uint16_t bin = bam_reg2bin (pos, last_pos); // zero-based, half-closed half-open [start,end)
            alignment->bin = LTEN16 (bin); // override the -1 previously set by the translator
        }
    }
    
    else if (flag.out_dt == DT_FASTQ) {
        // only analyze, but don't reconstruct CIGAR in FASTQ
    }
    
    vb_sam->last_cigar = snip;

    if (flag.regions && vb->chrom_node_index != WORD_INDEX_NONE && vb->contexts[SAM_POS].last_value.i && 
        !regions_is_range_included (vb->chrom_node_index, vb->contexts[SAM_POS].last_value.i, vb->contexts[SAM_POS].last_value.i + vb_sam->ref_consumed - 1, true))
        vb->dont_show_curr_line = true;

    return false; // no new value
}   

// Case 1: BIN is set to SPECIAL, we will set new_value here to -1 and wait for CIGAR to calculate it, 
//         as we need vb->ref_consumed - sam_piz_special_CIGAR will update the reconstruced value
// Case 2: BIN is an textual integer snip - its BIN.last_value will be set as normal and transltor will reconstruct it
SPECIAL_RECONSTRUCTOR (bam_piz_special_BIN)
{
    ctx->semaphore = true; // signal to sam_piz_special_CIGAR to calculate
    return false; // no new value
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_TLEN)
{
    ASSERT0 (snip_len, "Error in sam_piz_special_TLEN: snip_len=0");

    int32_t tlen_by_calc = atoi (snip);
    int32_t tlen_val = tlen_by_calc + vb->contexts[SAM_PNEXT].last_delta + vb->seq_len;

    new_value->i = tlen_val;

    if (reconstruct) { RECONSTRUCT_INT (tlen_val); }

    return true; // new value
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_AS)
{
    new_value->i = vb->seq_len - atoi (snip);
    if (reconstruct) { RECONSTRUCT_INT (new_value->i) };
    
    return true; // new value
}

// logic: snip is eg "119C" (possibly also "") - we reconstruct the original, eg "119C31" 
// by concating a number which is (seq_len - partial_seq_len_by_md_field)
SPECIAL_RECONSTRUCTOR (sam_piz_special_MD)
{
    if (!reconstruct) return false;
    
    if (snip_len) RECONSTRUCT (snip, snip_len);

    unsigned partial_seq_len_by_md_field = sam_seg_get_seq_len_by_MD_field (snip, snip_len);
    RECONSTRUCT_INT (vb->seq_len - partial_seq_len_by_md_field);

    return false; // no new value
}

// BD and BI - reconstruct from BD_BI context which contains interlaced BD and BI data. 
SPECIAL_RECONSTRUCTOR (sam_piz_special_BD_BI)
{
    if (!vb->seq_len || !reconstruct) goto done;

    Context *bdbi_ctx = ctx_get_existing_ctx (vb, dict_id_OPTION_BD_BI);

    // note: bd and bi use their own next_local to retrieve data from bdbi_ctx. the actual index
    // in bdbi_ctx.local is calculated given the interlacing
    ASSERT (ctx->next_local + vb->seq_len * 2 <= bdbi_ctx->local.len, "Error reading txt_line=%u: unexpected end of %s data", vb->line_i, dis_dict_id (ctx->dict_id).s);

    char *dst        = AFTERENT (char, vb->txt_data);
    const char *src  = ENT (char, bdbi_ctx->local, ctx->next_local * 2);
    uint32_t seq_len = vb->seq_len; // automatic var for effeciency

    if (ctx->dict_id.num == dict_id_OPTION_BD)
        for (uint32_t i=0; i < seq_len; i++, src+=2, dst++) *dst = *src;
    else
        for (uint32_t i=0; i < seq_len; i++, src+=2, dst++) *dst = *src + *(src+1);
    
    vb->txt_data.len += vb->seq_len;    
    ctx->next_local  += vb->seq_len;

done:
    return false; // no new value
}

// note of float reconstruction:
// When compressing SAM, floats are stored as a textual string, reconstruced natively for SAM and via sam_piz_sam2bam_FLOAT for BAM.
//    Done this way so when reconstructing SAM, the correct number of textual digits is reconstructed.
// When compressing BAM, floats are stored as 32-bit binaries, encoded as uint32, and stringified to a snip. They are reconstructed,
//    either as textual for SAM or binary for BAM via bam_piz_special_FLOAT. Done this way so BAM binary float is reconstructd precisely.
SPECIAL_RECONSTRUCTOR (bam_piz_special_FLOAT)
{
    // get Little Endian n
    int64_t n;
    ASSERT (str_get_int (snip, snip_len, &n), "Error in bam_piz_special_FLOAT: failed to read integer in %s", ctx->name);

    uint32_t lten_n = (uint32_t)n;         // n is now little endian, uint32 
    
    union { // n and f in machine endianity (4 bytes)
        uint32_t n;
        float f;
    } machine_en = { .n = LTEN32 (lten_n) };

    if (!reconstruct) goto finish;

    // binary reconstruction - BAM format
    if (flag.out_dt == DT_BAM)
        RECONSTRUCT (&machine_en.f, sizeof (float))
    
    // textual reconstruction - SAM format 
    else { 
        #define NUM_SIGNIFICANT_DIGITS 6 // 6 significant digits, as samtools does
        
        // calculate digits before and after the decimal point
        double log_f = log10 (machine_en.f >= 0 ? machine_en.f : -machine_en.f);
        unsigned int_digits = (log_f >= 0) + (unsigned)log_f;
        unsigned dec_digits = MAX (0, NUM_SIGNIFICANT_DIGITS - int_digits);
        
        // reconstruct number with exactly NUM_SIGNIFICANT_DIGITS digits
        sprintf (AFTERENT (char, vb->txt_data), "%.*f", dec_digits, machine_en.f); 
        unsigned len = strlen (AFTERENT (char, vb->txt_data)); 
        vb->txt_data.len += len;

        // remove trailing decimal zeros:  "5.500"->"5.5" ; "5.0000"->"5" ; "50"->"50"
        if (dec_digits) {
            unsigned trailing_zeros=0;
            for (int i=vb->txt_data.len-1; i >= vb->txt_data.len-dec_digits; i--)
                if (*ENT (char, vb->txt_data, i) == '0') 
                    trailing_zeros++;
                else
                    break;
            
            vb->txt_data.len -= (dec_digits==trailing_zeros) ? dec_digits+1 : trailing_zeros;
        }
    }

finish:
    new_value->d = (double)machine_en.f;
    return true; // have new value
}

//-----------------------------------------------------------------
// Translator functions for reconstructing SAM data into BAM format
//-----------------------------------------------------------------

// translate SAM ASCII sequence characters to BAM's 4-bit characters:
TRANSLATOR_FUNC (sam_piz_sam2bam_SEQ)
{
    // the characters "=ACMGRSVTWYHKDBN" are mapped to BAM 0->15, in this matrix we add 0x80 as a validity bit. All other characters are 0x00 - invalid
    static const uint8_t sam2bam_seq_map[256] = { ['=']=0x80, ['A']=0x81, ['C']=0x82, ['M']=0x83, ['G']=0x84, ['R']=0x85, ['S']=0x86, ['V']=0x87, 
                                                  ['T']=0x88, ['W']=0x89, ['Y']=0x8a, ['H']=0x8b, ['K']=0x8c, ['D']=0x8d, ['B']=0x8e, ['N']=0x8f };
    
    if (vb->dont_show_curr_line) return 0; // sequence was not reconstructed - nothing to translate
    
    static bool invalid_char_warning_shown = false; // we show this warning up to once per executipn

    BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)ENT (char, vb->txt_data, vb->line_start);
    uint32_t l_seq = LTEN32 (alignment->l_seq);

    // if l_seq=0, just remove the '*'
    if (!l_seq) {
        vb->txt_data.len--;
        return 0;
    }

    // if l_seq is odd, 0 the next byte that will be half of our last result byte
    if (l_seq % 2) *AFTERENT (char, vb->txt_data) = 0; 

    uint8_t *seq_before=(uint8_t *)reconstructed, *seq_after=(uint8_t *)reconstructed; 
    for (uint32_t i=0; i < (l_seq+1)/2; i++, seq_after++, seq_before += 2) {
        uint8_t base[2] = { sam2bam_seq_map[(uint8_t)seq_before[0]], sam2bam_seq_map[(uint8_t)seq_before[1]] };
        
        // check for invalid characters - issue warning (only once per execution), and make then into an 'N'
        for (unsigned b=0; b < 2; b++)
            if (!base[b] && !invalid_char_warning_shown && !(b==1 && (i+1)*2 > l_seq)) {
                WARN ("Warning when converting SAM sequence data to BAM: invalid character encodered, it will be converted as 'N': '%c' (ASCII %u)", base[b], base[b]);
                invalid_char_warning_shown = true;
                base[b] = 0x0f;
            }

        *seq_after = (base[0] << 4) | (base[1] & 0x0f);
    }

    vb->txt_data.len = vb->txt_data.len - l_seq + (l_seq+1)/2;

    return 0;
}

// translate SAM ASCII (33-based) Phread values to BAM's 0-based
TRANSLATOR_FUNC (sam_piz_sam2bam_QUAL)
{
    // if QUAL is "*" there are two options:
    // 1. If l_seq is 0, the QUAL is empty
    // 2. If not (i.e. we have SEQ data but not QUAL) - it is a string of 0xff, length l_seq
    if (reconstructed_len==1 && *reconstructed == '*') {
        BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)ENT (char, vb->txt_data, vb->line_start);
        uint32_t l_seq = LTEN32 (alignment->l_seq);

        if (!l_seq) // option 1
            vb->txt_data.len--;
        else {      // option 2
            memset (LASTENT (uint8_t, vb->txt_data), 0xff, l_seq); // override the '*' and l_seq-1 more
            vb->txt_data.len += l_seq - 1;
        }
    }
    
    else // we have QUAL - update Phred values
        for (uint32_t i=0; i < reconstructed_len; i++)
            reconstructed[i] -= 33; 

    return 0;
}

// output the word_index of RNAME, which is verified in ref_contigs_get_ref_chrom during seg
// to be the same as the reference id 
TRANSLATOR_FUNC (sam_piz_sam2bam_RNAME)
{
    DECLARE_SNIP;
    ctx_get_snip_by_word_index (&ctx->word_list, &ctx->dict, ctx->last_value.i, &snip, &snip_len);

    // if it is '*', reconstruct -1
    if (snip_len == 1 && *snip == '*') 
        RECONSTRUCT_BIN32 (-1)    

    // if its RNEXT and =, emit the last index of RNAME
    else if (ctx->did_i != CHROM && snip_len == 1 && *snip == '=') 
        RECONSTRUCT_BIN32 (vb->contexts[CHROM].last_value.i)

    // otherwise - output the word_index which was stored here because of flags.store=STORE_INDEX set in seg 
    else     
        RECONSTRUCT_BIN32 (ctx->last_value.i); 
    
    return 0;
}

// output, in binary form, POS-1 as BAM uses 0-based POS
TRANSLATOR_FUNC (sam_piz_sam2bam_POS)
{
    RECONSTRUCT_BIN32 (ctx->last_value.i - 1);
    return 0;
}

// place value in correct location in alignment
TRANSLATOR_FUNC (sam_piz_sam2bam_TLEN)
{
    BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)ENT (char, vb->txt_data, vb->line_start);
    alignment->tlen = LTEN32 (ctx->last_value.i);
    return 0;
}

// translate OPTIONAL SAM->BAM - called as translator-only item on within the Optional reconstruction
// fix prefix eg MX:i: -> MXs
TRANSLATOR_FUNC (sam_piz_sam2bam_OPTIONAL_SELF)
{
    Container *con = (Container *)reconstructed;

    if (reconstructed_len == -1) return 0; // no Optional data in this alignment

    char *prefixes_before = &reconstructed[sizeof_container (*con)] + 2; // +2 to skip the empty prefixes of container wide, and item[0]
    char *prefixes_after = prefixes_before;

    for (unsigned i=1; i < con->num_items; i++, prefixes_before+=6, prefixes_after+=4) {
        prefixes_after[0] = prefixes_before[0]; // tag[0] 
        prefixes_after[1] = prefixes_before[1]; // tag[1]
        prefixes_after[2] = prefixes_before[3]; // type
        prefixes_after[3] = SNIP_CONTAINER; // end of prefix

        // a SAM 'i' translate to one of several BAM types using the translator code
        // that may be 0->6 (NONE to SAM2BAM_LTEN_U32)
        if (prefixes_after[2] == 'i') 
            prefixes_after[2] = "\0cCsSiI"[con->items[i].translator];
    }

    return -2 * (con->num_items-1); // change in prefixes_len
}

// translate OPTIONAL SAM->BAM - called after Optional reconstruction is done
// sets block_size
TRANSLATOR_FUNC (sam_piz_sam2bam_OPTIONAL)
{
    BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)ENT (char, vb->txt_data, vb->line_start);
    alignment->block_size = vb->txt_data.len - vb->line_start - sizeof (uint32_t); // block_size doesn't include the block_size field itself
    alignment->block_size = LTEN32 (alignment->block_size);
    return 0;
}

//------------------------------------------------------------------------------------
// Translator and filter functions for reconstructing SAM / BAM data into FASTQ format
//------------------------------------------------------------------------------------

TXTHEADER_TRANSLATOR (txtheader_sam2fq)
{
    txtheader_buf->len = 0; // fastq has no header
}

// filtering during reconstruction: called by container_reconstruct_do for each sam alignment (repeat)
CONTAINER_FILTER_FUNC (sam_piz_sam2fq_filter)
{
    uint16_t flag = (uint16_t)vb->contexts[SAM_FLAG].last_value.i;
    return !(flag & 0x100); // show only if this is a primary alignment (don't show its secondary alignments)
}

// reverse-complement the sequence if needed, and drop if "*"
TRANSLATOR_FUNC (sam_piz_sam2fastq_SEQ)
{
    // full printable ascii with only AGCT->TCGA and agct->tcga
    static const char complement[256] = { [33]='!','"','#','$','%','&','\'','(',')','*','+',',','-','.','/',
                                          '0','1','2','3','4','5','6','7','8','9',':',';','<','=','>','?','@',
                                          'T','B','G','D','E','F','C','H','I','J','K','L','M','N','O','P','Q','R','S','A','U','V','W','X','Y','Z',
                                          '[','\\',']','^','_','`',
                                          't','b','g','d','e','f','c','h','i','j','k','l','m','n','o','p','q','r','s','a','u','v','w','x','y','z'
                                          ,'{','|','}','~' };

    uint16_t flag = (uint16_t)vb->contexts[SAM_FLAG].last_value.i;
    
    // case: SEQ is "*" - don't show this fastq record
    if (reconstructed_len==1 && *reconstructed == '*') 
        vb->dont_show_curr_line = true;

    // case: this sequence is reverse complemented - reverse-complement it
    else if (flag & 0x10) {

        // we move from the outside in, switching the left and right bases while also complementing both
        for (unsigned i=0; i < reconstructed_len / 2; i++) {
            char l_base = reconstructed[i];
            char r_base = reconstructed[reconstructed_len-1-i];

            reconstructed[i]                     = complement[(uint8_t)r_base];
            reconstructed[reconstructed_len-1-i] = complement[(uint8_t)l_base];
        }

        if (reconstructed_len % 2) // we have an odd number of bases - now complement the middle one
            reconstructed[reconstructed_len/2] = complement[(uint8_t)reconstructed[reconstructed_len/2]];
    }

    return 0;
}

// reverse the sequence if needed, and drop if "*"
TRANSLATOR_FUNC (sam_piz_sam2fastq_QUAL)
{
    uint16_t flag = (uint16_t)vb->contexts[SAM_FLAG].last_value.i;
    
    // case: QUAL is "*" - don't show this fastq record
    if (reconstructed_len==1 && *reconstructed == '*') 
        vb->dont_show_curr_line = true;

    // case: this sequence is reverse complemented - reverse the QUAL string
    else if (flag & 0x10) {

        // we move from the outside in, switching the left and right bases 
        for (unsigned i=0; i < reconstructed_len / 2; i++) {
            char tmp = reconstructed[i];
            reconstructed[i] = reconstructed[reconstructed_len-1-i];
            reconstructed[reconstructed_len-1-i] = tmp;
        }
    }

    return 0;
}

// emit 1 if (FLAGS & 0x40) or 2 of (FLAGS & 0x80)
TRANSLATOR_FUNC (sam_piz_sam2fastq_FLAG)
{
    uint16_t flag = (uint16_t)vb->contexts[SAM_FLAG].last_value.i;

    if (flag & 0x40) RECONSTRUCT ("1\n", 2); // usually R1
    if (flag & 0x80) RECONSTRUCT ("2\n", 2); // usually R2

    return 0;
}

// note of float reconstruction:
// When compressing SAM, floats are stored as a textual string, reconstruced natively for SAM and via sam_piz_sam2bam_FLOAT for BAM.
//    Done this way so when reconstructing SAM, the correct number of textual digits is reconstructed.
// When compressing BAM, floats are stored as 32-bit binaries, encoded as uint32, and stringified to a snip. They are reconstructed,
//    either as textual for SAM or binary for BAM via bam_piz_special_FLOAT. Done this way so BAM binary float is reconstructd precisely.
TRANSLATOR_FUNC (sam_piz_sam2bam_FLOAT)
{
    union {
        float f; // 32 bit float
        uint32_t i;
    } value;
    
    ASSERT0 (sizeof (value)==4, "Error in sam_piz_sam2fastq_FLOAT: expecting value to be 32 bits"); // should never happen

    value.f = (float)ctx->last_value.d;
    RECONSTRUCT_BIN32 (value.i);

    return 0;
}

// remove the comma from the prefix that contains the type, eg "i,"->"i"
TRANSLATOR_FUNC (sam_piz_sam2bam_ARRAY_SELF)
{
    Container *con = (Container *)reconstructed;
    char *prefixes = &reconstructed[sizeof_container (*con)];

    // remove the ',' from the prefix, and terminate with CON_PREFIX_SEP_SHOW_REPEATS - this will cause
    // the number of repeats (in LTEN32) to be outputed after the prefix
    prefixes[1] = CON_PREFIX_SEP_SHOW_REPEATS; // prefixes is now { type, CON_PREFIX_SEP_SHOW_REPEATS }
    
    return -1; // change in prefixes length
}
