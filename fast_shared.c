// ------------------------------------------------------------------
//   fast_shared.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "fast_private.h"
#include "file.h"
#include "endianness.h"
#include "piz.h"
#include "dict_id.h"
#include "reference.h"
#include "refhash.h"
#include "arch.h"
#include "strings.h"
#include "context.h"
#include "seg.h"
#include "compressor.h"

Structured structured_DESC;

unsigned fast_vb_size (void) { return sizeof (VBlockFAST); }
unsigned fast_vb_zip_dl_size (void) { return sizeof (ZipDataLineFAST); }

void fast_vb_release_vb (VBlockFAST *vb)
{
    vb->last_line = 0;
    vb->contig_grepped_out = false;
    vb->pair_num_lines = vb->pair_vb_i = 0;
    memset (&vb->desc_mapper, 0, sizeof (vb->desc_mapper));
}
   
// called by I/O thread in fast_piz_read_one_vb, in case of --grep, to decompress and reconstruct the desc line, to 
// see if this vb is included. 
bool fast_piz_test_grep (VBlockFAST *vb)
{
    ARRAY (const unsigned, section_index, vb->z_section_headers);

    SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)(vb->z_data.data + section_index[0]);
    vb->first_line       = BGEN32 (header->first_line);
    vb->lines.len        = BGEN32 (header->num_lines);
    vb->vb_data_size     = BGEN32 (header->vb_data_size);
    vb->longest_line_len = BGEN32 (header->longest_line_len);

    // in case of --unbind, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every sam component
    if (flag_unbind) vb->vblock_i = BGEN32 (header->h.vblock_i);
    
    // we only need room for one line for now 
    buf_alloc (vb, &vb->txt_data, vb->longest_line_len, 1.1, "txt_data", vb->vblock_i);

    // uncompress & map desc field (filtered by piz_is_skip_section)
    vb->grep_stages = GS_TEST; // tell piz_is_skip_section to skip decompressing sections not needed for determining the grep
    piz_uncompress_all_ctxs ((VBlockP)vb, 0);
    vb->grep_stages = GS_UNCOMPRESS; // during uncompress in the compute thread, uncompress only what was not already uncompressed here

    piz_map_compound_field ((VBlockP)vb, dict_id_is_fast_desc_sf, &vb->desc_mapper);

    // reconstruct each description line and check for string matching with flag_grep
    bool found = false, match = false;

    Context *desc_ctx =  &vb->contexts[vb->data_type == DT_FASTQ ? FASTQ_DESC : FASTA_DESC];
    desc_ctx->iterator.next_b250 = FIRSTENT (uint8_t, desc_ctx->b250); 

    vb->line_i = vb->data_type == DT_FASTQ ? 4 * vb->first_line : vb->first_line;

    while (desc_ctx->iterator.next_b250 < AFTERENT (uint8_t, desc_ctx->b250) ||
           desc_ctx->next_local < desc_ctx->local.len) {
        piz_reconstruct_from_ctx (vb, desc_ctx->did_i, 0);

        *AFTERENT (char, vb->txt_data) = 0; // terminate the desc string

        match = !!strstr (vb->txt_data.data, flag_grep);

        vb->txt_data.len = 0; // reset

        if (match) { // 
            found = true; // we've found a match to the grepped string
            if (vb->data_type == DT_FASTQ) break; // for FASTA, we need to go until the last line, for FASTQ, we can break here
        }

        if (vb->data_type == DT_FASTQ) vb->line_i += 4; // note: for FASTA we have no idea what txt line we're on, because we're only tracking DESC lines
    }

    // last FASTA - carry over whether its grepped to the next VB - in case next VB starts not from the description line
    // similarly, note whether the previous VB ended with a grepped sequence. If previous VB didn't have any description
    // i.e the entire VB was a sequence that started in an earlier VB - the grep status of the easier VB is carried forward
    if (vb->data_type == DT_FASTA) 
        // if the last contig of the previous vb was grepped in - then include this VB anyway
        found = fasta_initialize_contig_grepped_out (vb, desc_ctx->b250.len > 0, match) || found;

    // reset iterators - piz_fast*_reconstruct_vb will use them again 
    mtf_init_iterator (desc_ctx);
    for (unsigned sf_i=0; sf_i < vb->desc_mapper.num_subfields; sf_i++) 
        mtf_init_iterator (&vb->contexts[vb->desc_mapper.did_i[sf_i]]);

    return found; // no match found
}

void fast_seg_seq (VBlockFAST *vb, const char *seq, uint32_t seq_len, int seq_bitmap_field)
{
    // these 4 contexts are consecutive and in the same order for all relevant data_types in data_types.h
    Context *bitmap_ctx = &vb->contexts[seq_bitmap_field + 0]; // FASTx_SEQ_BITMAP
    Context *nonref_ctx = &vb->contexts[seq_bitmap_field + 1]; // FASTx_NONREF
    Context *gpos_ctx   = &vb->contexts[seq_bitmap_field + 2]; // FASTx_GPOS
    Context *strand_ctx = &vb->contexts[seq_bitmap_field + 3]; // FASTx_STRAND

    BitArray *bitmap = buf_get_bitarray (&bitmap_ctx->local);

    // case: compressing without a reference - all data goes to "nonref", and we have no bitmap
    if (flag_reference != REF_EXTERNAL && flag_reference != REF_EXT_STORE) {
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

    int64_t gpos;
    bool is_forward, is_all_ref;
    bool has_match = refhash_best_match ((VBlockP)vb, seq, seq_len, &gpos, &is_forward, &is_all_ref);

    // case: we're the 2nd of the pair - the bit represents whether this strand is equal to the pair's strand (expecting
    // it to be 1 in most cases - making the bitmap highly compressible)
    if (gpos_ctx->inst & CTX_INST_PAIR_LOCAL) {
        const BitArray *pair_strand = buf_get_bitarray (&strand_ctx->pair);
        bool pair_is_forward = bit_array_get (pair_strand, vb->line_i); // same location, in the pair's local
        buf_add_bit (&strand_ctx->local, is_forward == pair_is_forward);
    }
    // case: not 2nd in a pair - just store the strange
    else 
        buf_add_bit (&strand_ctx->local, is_forward);
    
    bit_index_t next_bit = buf_extend_bits (&bitmap_ctx->local, seq_len);

    ASSSEG (gpos >= 0 && gpos <= 0xffffffff, seq, "gpos=%"PRId64" is out of uint32_t range", gpos);
    
    // case: we're the 2nd of the pair - store a delta if its small enough, or a lookup from local if not
    bool store_local = true;
    if (gpos_ctx->inst & CTX_INST_PAIR_LOCAL) {
        int64_t pair_gpos = (int64_t) BGEN32 (*ENT (uint32_t, gpos_ctx->pair, vb->line_i)); // same location, in the pair's local
        int64_t gpos_delta = gpos - pair_gpos;

        if (gpos_delta <= MAX_POS_DELTA && gpos_delta >= -MAX_POS_DELTA) {
            store_local = false;      

            char delta_snip[30] = { SNIP_PAIR_DELTA };
            unsigned delta_str_len = str_int (gpos_delta, &delta_snip[1]);
            seg_by_ctx ((VBlockP)vb, delta_snip, delta_str_len + 1, gpos_ctx, 0, NULL);
        }
        else {
            static const char lookup[1] = { SNIP_LOOKUP }; // lookup from local
            seg_by_ctx ((VBlockP)vb, lookup, 1, gpos_ctx, 0, NULL);
        }
    }
    
    // store the GPOS in local if its not a 2nd pair, or if it is, but the delta is not small enough
    if (store_local)
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
