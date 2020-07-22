// ------------------------------------------------------------------
//   fast.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "fast_private.h"
#include "seg.h"
#include "context.h"
#include "file.h"
#include "strings.h"
#include "piz.h"
#include "optimize.h"
#include "dict_id.h"
#include "refhash.h"
#include "endianness.h"
#include "arch.h"

void fastq_seg_initialize (VBlockFAST *vb)
{
    // thread safety: this will be initialized by vb_i=1, while it holds a mutex in zip_compress_one_vb
    static bool structured_initialized = false;
    if (!structured_initialized) {
        seg_initialize_compound_structured ((VBlockP)vb, "D?ESC", &structured_DESC); 
        structured_initialized = true;
    }

    vb->contexts[FASTQ_SEQ_BITMAP].ltype = CTX_LT_SEQ_BITMAP; 
    vb->contexts[FASTQ_SEQ_STRAND].ltype = CTX_LT_SEQ_BITMAP;
    vb->contexts[FASTQ_SEQ_STRAND].flags = CTX_FL_LOCAL_NONE; // bz2 and lzma only make it biug
    vb->contexts[FASTQ_SEQ_GPOS].ltype   = CTX_LT_UINT32;
    vb->contexts[FASTQ_SEQ_GPOS].flags   = CTX_FL_LOCAL_LZMA;
    vb->contexts[FASTQ_SEQ_NOREF].flags  = CTX_FL_LOCAL_ACGT;
    vb->contexts[FASTQ_SEQ_NOREF].ltype  = CTX_LT_SEQUENCE;
    vb->contexts[FASTQ_QUAL].ltype       = CTX_LT_SEQUENCE;
}

void fastq_seg_finalize (VBlockFAST *vb)
{
/*    printf ("vb_i=%u bases=%u nonref=%u\n", vb->vblock_i, 
            (unsigned)buf_get_bitarray (&vb->contexts[FASTQ_SEQ_BITMAP].local)->num_of_bits, 
            (unsigned)vb->contexts[FASTQ_SEQ_NOREF].local.len);*/
}

static void fastq_seg_seq (VBlockFAST *vb, const char *seq, uint32_t seq_len)
{
    Context *nonref_ctx = &vb->contexts[FASTQ_SEQ_NOREF];
    Context *bitmap_ctx = &vb->contexts[FASTQ_SEQ_BITMAP];
    Context *gpos_ctx =   &vb->contexts[FASTQ_SEQ_GPOS];
    Context *strand_ctx = &vb->contexts[FASTQ_SEQ_STRAND];

    BitArray *bitmap = buf_get_bitarray (&bitmap_ctx->local);

    // case: compressing without a reference - all data goes to "nonref", and we have no bitmap
    if (!flag_reference) {
        buf_alloc (vb, &nonref_ctx->local, MAX (nonref_ctx->local.len + seq_len + 3, vb->lines.len * (seq_len + 5)), CTX_GROWTH, "context->local", nonref_ctx->did_i); 
        buf_add (&nonref_ctx->local, seq, seq_len);
        return;
    }

    // allocate bitmaps - provide name only if buffer is not allocated, to avoid re-writing param which would overwrite num_of_bits that overlays it + param must be 0
    buf_alloc (vb, &bitmap_ctx->local, MAX (bitmap_ctx->local.len + roundup_bits2words64 (seq_len) * sizeof(int64_t), vb->lines.len * (seq_len+5) / 8), CTX_GROWTH, 
               buf_is_allocated (&bitmap_ctx->local) ? NULL : "context->local", 0); 

    buf_alloc (vb, &strand_ctx->local, MAX (nonref_ctx->local.len + sizeof (int64_t), roundup_bits2words64 (vb->lines.len) * sizeof (int64_t)), CTX_GROWTH, 
               buf_is_allocated (&strand_ctx->local) ? NULL : "context->local", 0); 

    buf_alloc (vb, &nonref_ctx->local, MAX (nonref_ctx->local.len + seq_len + 3, vb->lines.len * seq_len / 4), CTX_GROWTH, "context->local", nonref_ctx->did_i); 
    buf_alloc (vb, &gpos_ctx->local,   MAX (nonref_ctx->local.len + sizeof (uint32_t), vb->lines.len * sizeof (uint32_t)), CTX_GROWTH, "context->local", gpos_ctx->did_i); 

    int64_t 
    
    gpos;
    bool is_forward, is_all_ref;
    bool has_match = refhash_best_match ((VBlockP)vb, seq, seq_len, &gpos, &is_forward, &is_all_ref);

    buf_add_bit (&strand_ctx->local, is_forward);
    
    bit_index_t next_bit = buf_extend_bits (&bitmap_ctx->local, seq_len);

    ASSSEG (gpos >= 0 && gpos <= 0xffffffff, seq, "gpos=%"PRId64" is out of uint32_t range", gpos);
    NEXTENT (uint32_t, gpos_ctx->local) = BGEN32 ((uint32_t)gpos);

    // shortcut if there's no reference match
    if (!has_match) {
        bit_array_clear_region (bitmap, next_bit, seq_len); // no bases match the reference
        buf_add (&nonref_ctx->local, seq, seq_len);
        return;
    }

    // we might have a nested spin lock for seq spans over a boundary. we always lock from the lower index to the higer
    if (flag_reference == REF_EXT_STORE) {
        mutex_lock (genome_muteces[GPOS2MUTEX(gpos)]);

        if (GPOS2MUTEX(gpos) != GPOS2MUTEX(gpos + seq_len - 1)) // region might span two spinlock 64K-regions, but not more than that   
            mutex_lock (genome_muteces[GPOS2MUTEX(gpos + seq_len - 1)]);
    }

    // shortcut if we have a full reference match
    if (is_all_ref) {
        bit_array_set_region (bitmap, next_bit, seq_len); // all bases match the reference
        
        if (flag_reference == REF_EXT_STORE) 
            bit_array_set_region (&genome->is_set, gpos, seq_len); // this region of the reference is used (in case we want to store it with REF_EXT_STORE)

        goto done;
    }

    int64_t room_fwd = genome_size - gpos; // how much reference forward might contain a match

    for (uint32_t i=0; i < seq_len; i++) {
                
        bool use_reference = false;

        // case our seq is identical to the reference at this site
        if (i < room_fwd) {
            char seq_base = is_forward ? seq[i] : complement[(uint8_t)seq[i]];
            
            int64_t ref_i = gpos + (is_forward ? i : seq_len-1-i);
            char ref_base = ACGT_DECODE (&genome->ref, ref_i);

            if (seq_base == ref_base) {
                
                if (flag_reference == REF_EXT_STORE) 
                    bit_array_set (&genome->is_set, ref_i); // we will need this ref to reconstruct

                use_reference = true;
            }
        }

        // case: we can't use the reference (different value than base or we have passed the end of the reference)
        if (!use_reference) 
            NEXTENT (char, nonref_ctx->local) = seq[i];

        bit_array_assign (bitmap, next_bit, use_reference);
        next_bit++;
    }

done:
    if (flag_reference == REF_EXT_STORE) {
        // unlock in reverse order
        if (GPOS2MUTEX(gpos) != GPOS2MUTEX(gpos + seq_len - 1)) // region might span to spinlock 64K-regions, but not more than that   
            mutex_unlock (genome_muteces[GPOS2MUTEX(gpos + seq_len - 1)]);

        mutex_unlock (genome_muteces[GPOS2MUTEX(gpos)]);
    }
}

// concept: we treat every 4 lines as a "line". the Description/ID is stored in DESC dictionary and segmented to subfields D?ESC.
// The sequence is stored in SEQ data. In addition, we utilize the TEMPLATE dictionary for metadata on the line, namely
// the length of the sequence and whether each line has a \r.
const char *fastq_seg_txt_line (VBlockFAST *vb, const char *field_start_line, bool *has_13)     // index in vb->txt_data where this line starts
{
    ZipDataLineFAST *dl = DATA_LINE (vb->line_i);

    const char *next_field, *field_start=field_start_line;
    unsigned field_len=0;
    char separator;

    int32_t len = (int32_t)(AFTERENT (char, vb->txt_data) - field_start_line);

    // the leading @ - just verify it (it will be included in D0ESC subfield)
    ASSSEG (*field_start != '\n', field_start, "%s: Invalid FASTQ file format: unexpected newline", global_cmd);

    ASSSEG (*field_start == '@', field_start, "%s: Invalid FASTQ file format: expecting description line to start with @ but it starts with %c",
            global_cmd, *field_start);

    // DESC - the description/id line is vendor-specific. example:
    // @A00910:85:HYGWJDSXX:1:1101:3025:1000 1:N:0:CAACGAGAGC+GAATTGAGTG (<-- this is Illumina format)
    // See here for details of Illumina subfields: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
    next_field = seg_get_next_line (vb, field_start, &len, &field_len, has_13, "DESC");
 
    // we segment it using / | : and " " as separators. 
    seg_compound_field ((VBlockP)vb, &vb->contexts[FASTQ_DESC], field_start, field_len, &vb->desc_mapper, structured_DESC, true, 0);
    SEG_EOL (FASTQ_E1L, true);

    // SEQ - just get the whole line
    const char *seq_start = next_field;
    dl->seq_data_start = next_field - vb->txt_data.data;
    next_field = seg_get_next_item (vb, next_field, &len, true, false, false, &dl->seq_len, &separator, has_13, "SEQ");

    fastq_seg_seq (vb, seq_start, dl->seq_len);

    // Add LOOKUP snip with seq_len
    char snip[10];
    snip[0] = SNIP_LOOKUP;
    unsigned seq_len_str_len = str_int (dl->seq_len, &snip[1]);
    seg_by_did_i (vb, snip, 1 + seq_len_str_len, FASTQ_SEQ_BITMAP, dl->seq_len); 

    SEG_EOL (FASTQ_E1L, true);

    // PLUS - next line is expected to be a "+"
    GET_LAST_ITEM ("+");
    ASSSEG (*field_start=='+' && field_len==1, field_start, "%s: Invalid FASTQ file format: expecting middle line to be a \"+\" (with no spaces) but it is \"%.*s\"",
            global_cmd, field_len, field_start);

    seg_by_did_i (vb, "+", 1, FASTQ_PLUS, 1);
    SEG_EOL (FASTQ_E1L, true);

    // QUAL - just get the whole line and make sure its length is the same as SEQ
    dl->qual_data_start = next_field - vb->txt_data.data;
    GET_LAST_ITEM ("QUAL");
    vb->contexts[FASTQ_QUAL].local.len += dl->seq_len;
    vb->contexts[FASTQ_QUAL].txt_len   += dl->seq_len;

    // End Of Line    
    SEG_EOL (FASTQ_E1L, true);

    ASSSEG (field_len == dl->seq_len, field_start, "%s: Invalid FASTQ file format: sequence_len=%u and quality_len=%u. Expecting them to be the same.\nSEQ=%.*s\nQUAL==%.*s",
            global_cmd, dl->seq_len, field_len, dl->seq_len, seq_start, field_len, field_start);
 
    return next_field;
}

// callback function for compress to get data of one line (called by comp_compress_bzlib)
void fastq_zip_get_start_len_line_i_qual (VBlock *vb, uint32_t vb_line_i, 
                                          char **line_qual_data, uint32_t *line_qual_len, // out
                                          char **unused_data,   uint32_t *unused_len) 
{
    ZipDataLineFAST *dl = DATA_LINE (vb_line_i);
     
    *line_qual_data = ENT (char, vb->txt_data, dl->qual_data_start);
    *line_qual_len  = dl->seq_len;
    *unused_data    = NULL;
    *unused_len     = 0;

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    if (flag_optimize_QUAL) optimize_phred_quality_string (*line_qual_data, *line_qual_len);
}

// returns true if section is to be skipped reading / uncompressing
bool fastq_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    if (!vb) return false; // we don't skip reading any SEC_DICT sections

    // note that piz_read_global_area rewrites --header-only as flag_header_one
    if (flag_header_one && 
        (dict_id.num == dict_id_fields[FASTQ_SEQ_NOREF] || dict_id.num == dict_id_fields[FASTQ_QUAL] || dict_id.num == dict_id_fields[FASTQ_PLUS]))
        return true;
        
    // when grepping by I/O thread - skipping all sections but DESC
    if (vb && flag_grep && (vb->grep_stages == GS_TEST) && 
        dict_id.num != dict_id_fields[FASTQ_DESC] && !dict_id_is_fast_desc_sf (dict_id))
        return true;

    // if grepping, compute thread doesn't need to decompressed DESC again
    if (vb && flag_grep && (vb->grep_stages == GS_UNCOMPRESS) && 
        (dict_id.num == dict_id_fields[FASTQ_DESC] || dict_id_is_fast_desc_sf (dict_id)))
        return true;

    return false;
}

// PIZ: SEQ reconstruction 
void fastq_piz_reconstruct_seq (VBlock *vb_, Context *bitmap_ctx, const char *seq_len_str, unsigned seq_len_str_len)
{
    VBlockFAST *vb = (VBlockFAST *)vb_;
    ASSERT0 (bitmap_ctx && bitmap_ctx->did_i == FASTQ_SEQ_BITMAP, "Error in fastq_piz_reconstruct_seq: context is not FASTQ_SEQ_BITMAP");

    if (piz_is_skip_section (vb, SEC_LOCAL, bitmap_ctx->dict_id)) return; // if case we need to skip the SEQ field (for the entire file)

    Context *nonref_ctx = &vb->contexts[FASTQ_SEQ_NOREF];
    Context *gpos_ctx =   &vb->contexts[FASTQ_SEQ_GPOS];
    Context *strand_ctx = &vb->contexts[FASTQ_SEQ_STRAND];

    // get seq_len, gpos and direction
    vb->seq_len = atoi (seq_len_str); // terminated by SNIP_SEP=0
    
    if (buf_is_allocated (&bitmap_ctx->local)) { // not all non-ref
        uint32_t gpos_be = NEXTLOCAL (uint32_t, gpos_ctx);
        int64_t gpos = BGEN32 (gpos_be); 
        bool forward_strand = NEXTLOCALBIT (strand_ctx);

        if (forward_strand)  // normal (note: this condition test is outside of the tight loop)
            for (uint32_t i=0; i < vb->seq_len; i++)
                if (NEXTLOCALBIT (bitmap_ctx))  // get base from reference
                    RECONSTRUCT1 (ACGT_DECODE (&genome->ref, gpos + i));
                else  // get base from nonref
                    RECONSTRUCT1 (NEXTLOCAL (char, nonref_ctx));

        else // reverse complement
            for (uint32_t i=0; i < vb->seq_len; i++) 
                if (NEXTLOCALBIT (bitmap_ctx))  // case: get base from reference
                    RECONSTRUCT1 (complement [(uint8_t)ACGT_DECODE (&genome->ref, gpos + vb->seq_len-1 - i)]);
                else  // case: get base from nonref
                    RECONSTRUCT1 (NEXTLOCAL (char, nonref_ctx));
    }
    else {
        RECONSTRUCT (ENT (char, nonref_ctx->local, nonref_ctx->next_local), vb->seq_len);
        nonref_ctx->next_local += vb->seq_len;
    }
}

void fastq_piz_reconstruct_vb (VBlockFAST *vb)
{
    ASSERT0 (!flag_reference || genome, "Error in fastq_piz_reconstruct_vb: reference is not loaded correctly");

    if (!flag_grep) piz_map_compound_field ((VBlockP)vb, dict_id_is_fast_desc_sf, &vb->desc_mapper); // it not already done during grep

    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        vb->line_i = 4 * (vb->first_line + vb_line_i); // each vb line is a fastq record which is 4 txt lines
        vb->dont_show_curr_line = false; // might become true due --regions or --grep
        
        uint32_t txt_data_start_line = vb->txt_data.len;

        piz_reconstruct_from_ctx (vb, FASTQ_DESC, 0);
        piz_reconstruct_from_ctx (vb, FASTQ_E1L,  0);

        // minor bug here: since all the EOL fields are aliases, in case header_one, we dont consume the
        // missing line EOLs and they will show wrong ones to the remaining lines. minor bc in virtually all files,
        // the EOL is identical in all lines
        if (!flag_header_one) { // note that piz_read_global_area rewrites --header-only as flag_header_one
            piz_reconstruct_from_ctx (vb, FASTQ_SEQ_BITMAP, 0);
            piz_reconstruct_from_ctx (vb, FASTQ_E2L,  0);
            piz_reconstruct_from_ctx (vb, FASTQ_PLUS, 0);
            piz_reconstruct_from_ctx (vb, FASTQ_E3L,  0);
            piz_reconstruct_from_ctx (vb, FASTQ_QUAL, 0);
            piz_reconstruct_from_ctx (vb, FASTQ_E4L,  0);
        }

        // case: we're grepping, and this line doesn't match
        *AFTERENT (char, vb->txt_data) = 0; // for strstr
        if (vb->dont_show_curr_line || (flag_grep && !strstr (ENT (char, vb->txt_data, txt_data_start_line), flag_grep)))
            vb->txt_data.len = txt_data_start_line; // rollback
    }
}
