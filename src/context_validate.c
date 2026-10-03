// ------------------------------------------------------------------
//   context_validate.c
//   Copyright (C) 2026-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying
//   the source code is strictly prohibited and subject to penalties
//   specified in the license.
// ------------------------------------------------------------------

#include <stddef.h>
#include "context_struct.h"

#define ASSERT_CONTEXT_OFFSET(field, expected_word) \
    ASSERT (offsetof(Context, field) == 8 * (expected_word), \
            "Context " #field ": expected word %zu (byte %zu), " \
            "actual word %zu (byte %zu)", \
            (size_t)(expected_word), \
            8 * (size_t)(expected_word), \
            offsetof(Context, field) / 8, \
            offsetof(Context, field))

#define ASSERT_BUFFER_ALIGNED_64(field) \
    ASSERT (offsetof(Context, field) % 64 == 0, \
            "Context Buffer " #field ": byte offset %zu (word %zu) " \
            "is not 64-byte aligned", \
            offsetof(Context, field), \
            offsetof(Context, field) / 8)

void context_validate (void)
{
    ASSERT_CONTEXT_OFFSET(tag_name, 0);
    ASSERT_CONTEXT_OFFSET(dict_id, 8);
    ASSERT_CONTEXT_OFFSET(dyn_int_min, 9);
    ASSERT_CONTEXT_OFFSET(did_i, 11);
    ASSERT_CONTEXT_OFFSET(local_in_z_len, 12);
    ASSERT_CONTEXT_OFFSET(b250_in_z_len, 13);
    ASSERT_CONTEXT_OFFSET(unused_word_14, 14);
    ASSERT_CONTEXT_OFFSET(dict, 16);
    ASSERT_CONTEXT_OFFSET(b250R1, 40);
    ASSERT_CONTEXT_OFFSET(counts, 48);
    ASSERT_CONTEXT_OFFSET(con_cache, 64);
    ASSERT_CONTEXT_OFFSET(ol_chrom2ref_map, 72);
    ASSERT_CONTEXT_OFFSET(last_value, 80);
    ASSERT_CONTEXT_OFFSET(last_delta, 81);
    ASSERT_CONTEXT_OFFSET(last_txt, 82);
    ASSERT_CONTEXT_OFFSET(last_line_i, 83);
    ASSERT_CONTEXT_OFFSET(ctx_specific, 84);
    ASSERT_CONTEXT_OFFSET(iterator, 85);
    ASSERT_CONTEXT_OFFSET(next_local, 86);
    ASSERT_CONTEXT_OFFSET(pair_flags, 87);
    ASSERT_CONTEXT_OFFSET(nodes, 88);
    ASSERT_CONTEXT_OFFSET(global_hash, 96);
    ASSERT_CONTEXT_OFFSET(txt_len, 104);
    ASSERT_CONTEXT_OFFSET(num_new_entries_prev_merged_vb, 106);
    ASSERT_CONTEXT_OFFSET(seg_to_local, 107);
    // §108: STR declaration macro -- offsetof check unavailable
    ASSERT_CONTEXT_OFFSET(ol_dict, 112);
    ASSERT_CONTEXT_OFFSET(unused136, 136);
    // §108: bit-field -- offsetof check unavailable
    ASSERT_CONTEXT_OFFSET(num_failed_singletons, 109);
    ASSERT_CONTEXT_OFFSET(ctx_mutex, 110);
    ASSERT_CONTEXT_OFFSET(word_list, 88);
    ASSERT_CONTEXT_OFFSET(piz_ctx_specific_buf, 112);
    ASSERT_CONTEXT_OFFSET(curr_container, 120);
    ASSERT_CONTEXT_OFFSET(last_wi, 121);
    ASSERT_CONTEXT_OFFSET(other_did_i, 122);

    /* Every Buffer in Context must begin on a 64-byte boundary. */
    ASSERT_BUFFER_ALIGNED_64(dict);
    ASSERT_BUFFER_ALIGNED_64(b250);
    ASSERT_BUFFER_ALIGNED_64(local);
    ASSERT_BUFFER_ALIGNED_64(b250R1);
    ASSERT_BUFFER_ALIGNED_64(alts);
    ASSERT_BUFFER_ALIGNED_64(last_samples);
    ASSERT_BUFFER_ALIGNED_64(sample_copied);
    ASSERT_BUFFER_ALIGNED_64(lookback);
    ASSERT_BUFFER_ALIGNED_64(width_count);
    ASSERT_BUFFER_ALIGNED_64(vep_spec);
    ASSERT_BUFFER_ALIGNED_64(counts);
    ASSERT_BUFFER_ALIGNED_64(con_index);
    ASSERT_BUFFER_ALIGNED_64(con_cache);
    ASSERT_BUFFER_ALIGNED_64(ctx_cache);
    ASSERT_BUFFER_ALIGNED_64(chrom2ref_map);
    ASSERT_BUFFER_ALIGNED_64(snip_cache);
    ASSERT_BUFFER_ALIGNED_64(packed);
    ASSERT_BUFFER_ALIGNED_64(subdicts);
    ASSERT_BUFFER_ALIGNED_64(template);
    ASSERT_BUFFER_ALIGNED_64(value_to_bin);
    ASSERT_BUFFER_ALIGNED_64(longr_state);
    ASSERT_BUFFER_ALIGNED_64(qual_line);
    ASSERT_BUFFER_ALIGNED_64(normalize_buf);
    ASSERT_BUFFER_ALIGNED_64(qname_hash);
    ASSERT_BUFFER_ALIGNED_64(interlaced);
    ASSERT_BUFFER_ALIGNED_64(mi_history);
    ASSERT_BUFFER_ALIGNED_64(XG);
    ASSERT_BUFFER_ALIGNED_64(deep_nonref);
    ASSERT_BUFFER_ALIGNED_64(deep_cigar);
    ASSERT_BUFFER_ALIGNED_64(bamass_cigar);
    ASSERT_BUFFER_ALIGNED_64(format_mapper_buf);
    ASSERT_BUFFER_ALIGNED_64(last_format);
    ASSERT_BUFFER_ALIGNED_64(id_hash);
    ASSERT_BUFFER_ALIGNED_64(info_items);
    ASSERT_BUFFER_ALIGNED_64(deferred_snip);
    ASSERT_BUFFER_ALIGNED_64(ol_chrom2ref_map);
    ASSERT_BUFFER_ALIGNED_64(ref2chrom_map);
    ASSERT_BUFFER_ALIGNED_64(con_len);
    ASSERT_BUFFER_ALIGNED_64(localR1);
    ASSERT_BUFFER_ALIGNED_64(format_contexts);
    ASSERT_BUFFER_ALIGNED_64(sf_i);
    ASSERT_BUFFER_ALIGNED_64(insertion);
    ASSERT_BUFFER_ALIGNED_64(huffman);
    ASSERT_BUFFER_ALIGNED_64(piz_is_set);
    ASSERT_BUFFER_ALIGNED_64(nodes);
    ASSERT_BUFFER_ALIGNED_64(global_hash);
    ASSERT_BUFFER_ALIGNED_64(ol_dict);
    ASSERT_BUFFER_ALIGNED_64(ol_nodes);
    ASSERT_BUFFER_ALIGNED_64(local_hash);
    ASSERT_BUFFER_ALIGNED_64(ston_hash);
    ASSERT_BUFFER_ALIGNED_64(ston_ents);
    ASSERT_BUFFER_ALIGNED_64(word_list);
    ASSERT_BUFFER_ALIGNED_64(history);
    ASSERT_BUFFER_ALIGNED_64(dropped_txt);
    ASSERT_BUFFER_ALIGNED_64(piz_ctx_specific_buf);
    ASSERT_BUFFER_ALIGNED_64(piz_word_list_hash);
    ASSERT_BUFFER_ALIGNED_64(cigar_anal_history);
    ASSERT_BUFFER_ALIGNED_64(line_sqbitmap);
    ASSERT_BUFFER_ALIGNED_64(domq_denorm);
    ASSERT_BUFFER_ALIGNED_64(channel_data);
    ASSERT_BUFFER_ALIGNED_64(homopolymer);

    /* Verify sizes of project types actually used as Context fields. */
    ASSERT (sizeof (Buffer) == 64, "sizeof(Buffer)=%zu, expected 64", sizeof (Buffer));

    ASSERT (sizeof (DictId) == 8, "sizeof(DictId)=%zu, expected 8", sizeof (DictId));
    ASSERT (sizeof (Did) == 2, "sizeof(Did)=%zu, expected 2", sizeof (Did));
    ASSERT (sizeof (LocalType) == 1, "sizeof(LocalType)=%zu, expected 1", sizeof (LocalType));
    ASSERT (sizeof (struct FlagsCtx) == 1, "sizeof(struct FlagsCtx)=%zu, expected 1", sizeof (struct FlagsCtx));
    ASSERT (sizeof (struct FlagsDict) == 1, "sizeof(struct FlagsDict)=%zu, expected 1", sizeof (struct FlagsDict));
    ASSERT (sizeof (B250Size) == 1, "sizeof(B250Size)=%zu, expected 1", sizeof (B250Size));
    ASSERT (sizeof (Codec) == 1, "sizeof(Codec)=%zu, expected 1", sizeof (Codec));

    ASSERT (sizeof (ValueType) == 8, "sizeof(ValueType)=%zu, expected 8", sizeof (ValueType));
    ASSERT (sizeof (WordIndex) == 4, "sizeof(WordIndex)=%zu, expected 4", sizeof (WordIndex));
    ASSERT (sizeof (TxtWord) == 8, "sizeof(TxtWord)=%zu, expected 8", sizeof (TxtWord));
    ASSERT (sizeof (LineIType) == 4, "sizeof(LineIType)=%zu, expected 4", sizeof (LineIType));
    ASSERT (sizeof (IdType) == 1, "sizeof(IdType)=%zu, expected 1", sizeof (IdType));
    ASSERT (sizeof (SamFlags) == 2, "sizeof(SamFlags)=%zu, expected 2", sizeof (SamFlags));
    ASSERT (sizeof (PosType32) == 4, "sizeof(PosType32)=%zu, expected 4", sizeof (PosType32));
    ASSERT (sizeof (thool) == 1, "sizeof(thool)=%zu, expected 1", sizeof (thool));
    ASSERT (sizeof (Ploidy) == 1, "sizeof(Ploidy)=%zu, expected 1", sizeof (Ploidy));
    ASSERT (sizeof (ContextP) == 8, "sizeof(ContextP)=%zu, expected 8", sizeof (ContextP));
    ASSERT (sizeof (SnipIterator) == 8, "sizeof(SnipIterator)=%zu, expected 8", sizeof (SnipIterator));
    ASSERT (sizeof (BamAssTrimCigarTreatment) == 1, "sizeof(BamAssTrimCigarTreatment)=%zu, expected 1", sizeof (BamAssTrimCigarTreatment));
    ASSERT (sizeof (StoreType) == 1, "sizeof(StoreType)=%zu, expected 1", sizeof (StoreType));
    ASSERT (sizeof (LocalDepType) == 1, "sizeof(LocalDepType)=%zu, expected 1", sizeof (LocalDepType));
    ASSERT (sizeof (ContainerP) == 8, "sizeof(ContainerP)=%zu, expected 8", sizeof (ContainerP));
    ASSERT (sizeof (SectionType) == 1, "sizeof(SectionType)=%zu, expected 1", sizeof (SectionType));
    ASSERT (sizeof (SpecialResult) == 1, "sizeof(SpecialResult)=%zu, expected 1", sizeof (SpecialResult));

    /* Types defined inside Context. */
    ASSERT (sizeof (struct ctx_tp) == 8,
            "sizeof(struct ctx_tp)=%zu, expected 8", sizeof (struct ctx_tp));

    ASSERT (sizeof (Context) % 64 == 0,
            "sizeof(Context)=%zu (8B-words=%1.3f) is not a multiple of 64: %u %% 64 = %u",
            sizeof (Context), (double)sizeof (Context) / 8.0, sizeof (Context), sizeof (Context) % 64);            

    ASSERT (RECON_STATE_SIZE <= 61, "Expecting RECON_STATE_SIZE=%u <= 61 (see Ice struct)", RECON_STATE_SIZE);

    if (flag.show_memory)
        iprintf ("Context-related sizes: Context=%u (%u x 64B) Mutex=%u RECON_STATE_SIZE=%u\n", 
                 (unsigned)sizeof(Context), (unsigned)sizeof(Context)/64, (unsigned)sizeof(Mutex), (unsigned)RECON_STATE_SIZE);
}
