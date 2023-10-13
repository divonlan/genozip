// ------------------------------------------------------------------
//   fastq_deep.c
//   Copyright (C) 2023-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "fastq_private.h"
#include "deep.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "bits.h"
#include "huffman.h"
#include "container.h"
#include "codec.h"
#include "htscodecs/arith_dynamic.h"

sSTRl(con_decanonize1_snip,96);
sSTRl(con_decanonize2_snip,96);

//----------------
// ZIP
//----------------

void fastq_deep_zip_initialize (void)
{
    SmallContainer con1 = { .repeats = 1, .nitems_lo = 2, .items = { { .dict_id.num = _FASTQ_Q0NAME }, { .dict_id.num = _FASTQ_Q1NAME } }};
    container_prepare_snip ((ContainerP)&con1, NULL, 0, qSTRa(con_decanonize1_snip));

    SmallContainer con2 = { .repeats = 1, .nitems_lo = 2, .items = { { .dict_id.num = _FASTQ_Q0NAME2 }, { .dict_id.num = _FASTQ_Q1NAME2 } }};
    container_prepare_snip ((ContainerP)&con2, NULL, 0, qSTRa(con_decanonize2_snip));
}

void fastq_deep_zip_after_compute (VBlockFASTQP vb)
{
    for (int i=0; i < NUM_DEEP_STATS; i++)
        z_file->deep_stats[i] += vb->deep_stats[i];
}

void fastq_deep_seg_initialize (VBlockFASTQP vb)
{
    ctx_set_store (VB, STORE_INT, FASTQ_DEEP, DID_EOL); // actually, we store a pointer into one of the Buffers in z_file->deep_ents, but we treat it as an int
    ctx_set_ltype (VB, LT_DYN_INT, FASTQ_DEEP, DID_EOL); 

    if (flag.pair == PAIR_R1 || flag.pair == NOT_PAIRED) 
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_set_deep, '0' }, 3, FASTQ_DEEP, 0);  // all-the-same for FASTQ_DEEP

    else { // pair-2
        ctx_set_ltype (VB, LT_INT8, FASTQ_DEEP_DELTA, DID_EOL); 
        ctx_consolidate_stats (VB, FASTQ_DEEP, FASTQ_DEEP_DELTA, DID_EOL);
    }
}

void fastq_deep_zip_show_entries_stats (void) 
{
    static rom names[] = NO_DEEP_NAMES;

    uint64_t total = z_file->deep_stats[NDP_FQ_READS];
    iprint0 ("\nFASTQ reads breakdown by deepability:\n");
    
    for (int i=0; i < NUM_DEEP_STATS; i++) 
        if (z_file->deep_stats[i])
            iprintf ("%-11.11s: %"PRIu64" (%.1f%%)\n", names[i], z_file->deep_stats[i], 100.0 * (double)z_file->deep_stats[i] / (double)total);

    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    uint64_t count_not_consumed=0, count_dups=0;
    for (uint64_t ent_i=0; ent_i < deep_ents_len; ent_i++) 
        if (deep_ents[ent_i].dup) count_dups++;
        else if (!deep_ents[ent_i].consumed) count_not_consumed++;

    #define NUM_SHOW 20
    iprintf ("\nNumber of %s deepable alignments not consumed by FASTQ reads: %"PRIu64"\n"
             "- Unusable due to multiple alignments with same hash: %"PRIu64"\n"
             "- No matching FASTQ read: %"PRIu64" (showing first %d)\n",
             z_dt_name(), count_not_consumed + count_dups, count_dups, count_not_consumed, (int)MIN_(count_not_consumed, NUM_SHOW));

    // show the first few unconsumed
    if (count_not_consumed)
        for (uint64_t ent_i=0, count=0; ent_i < deep_ents_len && count < NUM_SHOW; ent_i++)
            if (!deep_ents[ent_i].consumed && !deep_ents[ent_i].dup) {
                iprintf ("sam_vb=%u line_i=%u hash=%u,%u,%u\n", deep_ents[ent_i].vb_i, deep_ents[ent_i].line_i, DEEPHASHf(deep_ents[ent_i].hash));
                count++;
            }
}

void fastq_deep_seg_finalize_segconf (uint32_t n_lines)
{
    unsigned n_mch = n_lines - segconf.n_no_mch; // exclude lines that don't match any SAM line - perhaps they were filtered out and excluded from the SAM file 
    unsigned threashold = flag.force_deep ? 1 : MAX_(n_mch/2, 1); // at least half of the lines need this level of matching (or at least 1 if --force-deep)

    // test: at least half of the reads (excluding reads that had no matching SEQ in the SAM) have a matching QNAME
    segconf.deep_qtype = (segconf.n_seq_qname_mch[0] + segconf.n_full_mch[0] >= threashold) ? QNAME1
                       : (segconf.n_seq_qname_mch[1] + segconf.n_full_mch[1] >= threashold) ? QNAME2
                       : QNONE;

    if (segconf.deep_qtype == QNONE) segconf.deep_has_trimmed = false;

    if (!segconf.deep_has_trimmed) buf_destroy (z_file->deep_index_by[BY_QNAME]);

    if (flag.show_deep) 
        iprintf ("segconf: n_lines=%u threashold=%u n_full_mch=(%u,%u) n_seq_qname_mch=(%u,%u) n_seq_qual_mch=%u n_seq_only=%u n_no_match=%u has_trimmed=%s\n",
                 n_lines, threashold, segconf.n_full_mch[0], segconf.n_full_mch[1], segconf.n_seq_qname_mch[0], segconf.n_seq_qname_mch[1], 
                 segconf.n_seq_qual_mch, segconf.n_seq_mch, segconf.n_no_mch, TF(segconf.deep_has_trimmed));

    // test: likewise, matching QUAL
    if ((segconf.n_seq_qual_mch + segconf.n_full_mch[0] < threashold) &&
        (segconf.n_seq_qual_mch + segconf.n_full_mch[1] < threashold))
        segconf.deep_no_qual = true;

    // we need at least 2 of the 3 (QNANE,SEQ,QUAL) to be comparable between FASTQ and SAM to have a total of 64 bit hash (32 from each)    
    ASSINP (segconf.deep_qtype != QNONE || !segconf.deep_no_qual, 
            "Error: cannot use --deep with this file: based on testing the first %u reads of %s, only:\n"
            "%u reads had a matching alignment in the SAM/BAM file matching QNAME, SEQ, QUAL\n"
            "%u reads had a matching alignment in the SAM/BAM file matching QNAME2, SEQ, QUAL\n"
            "%u matching QNAME and SEQ\n"
            "%u matching QNAME2 and SEQ\n"
            "%u matching SEQ and QUAL\n"
            "This is below the threashold needed for meaningful Deep compression. You may compress these files without --deep, or override with --force-deep",
            n_lines, txt_name, segconf.n_full_mch[0], segconf.n_full_mch[1], segconf.n_seq_qname_mch[0], segconf.n_seq_qname_mch[1], segconf.n_seq_qual_mch);
}

static void fastq_deep_seg_segconf (VBlockFASTQP vb, STRp(qname), STRp(qname2), STRp(seq), STRp(qual))
{
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    int num_qnames = 1 + (qname2_len > 0);

    uint32_t qname_hash[2] = { deep_qname_hash (QNAME1, STRa(qname),  NULL),
                               ((num_qnames==2) ? deep_qname_hash (QNAME2, STRa(qname2), NULL) : 0) };
    uint32_t seq_hash  = deep_seq_hash (VB, STRa(seq), false);
    uint32_t qual_hash = deep_qual_hash (VB, STRa(qual), false);

    bool seq_qual_matches=false, seq_qname_matches[2]={}, seq_matches=false;
    uint32_t hash = seq_hash & bitmask32 (num_hash_bits);
    
    for (uint32_t ent_i = *B32(z_file->deep_index_by[BY_SEQ], hash); ent_i != NO_NEXT; ent_i = deep_ents[ent_i].next[BY_SEQ]) {
        ZipZDeep *e = &deep_ents[ent_i];
        
        if (e->hash.seq != seq_hash) continue; // possible, because hash uses less bits that e->hash.seq

        for (int qn=0; qn < num_qnames; qn++) 
            if (e->hash.qname == qname_hash[qn]) {
                if (e->hash.qual == qual_hash) {
                    segconf.n_full_mch[qn]++;
                    return;
                }
                else {
                    seq_qname_matches[qn] = true;
                    goto next_ent;
                }
            }

        if (e->hash.qual == qual_hash)
            seq_qual_matches = true;
        else
            seq_matches = true;

        next_ent: continue;
    }

    // each segconf read increments exactly one category
    if      (seq_qname_matches[0]) segconf.n_seq_qname_mch[0]++; // if one match is seq+qname and another is seq+qual, then count the seq+qname 
    else if (seq_qname_matches[1]) segconf.n_seq_qname_mch[1]++; 
    else if (seq_qual_matches)     segconf.n_seq_qual_mch++;
    else if (seq_matches)          segconf.n_seq_mch++;        // no match for qname and qual, only seq alone
    
    // no match based on SEQ index, perhaps the read appears trimmed (beyond cropped) in the SAM data - search based on QNAME index
    else {
        for (int qn=0; qn < num_qnames; qn++) {
            hash = qname_hash[qn] & bitmask32 (num_hash_bits);
            
            for (uint32_t ent_i = *B32(z_file->deep_index_by[BY_QNAME], hash); ent_i != NO_NEXT; ent_i = deep_ents[ent_i].next[BY_QNAME]) {
                ZipZDeep *e = &deep_ents[ent_i];
                
                if (e->hash.qname != qname_hash[qn] || // possible, because hash uses less bits that e->hash.qname
                    e->seq_len >= seq_len) continue;   // definitely not a match if FASTQ seq_len is shorter than SAM's and also not equal seq_len, as that would have been a no-trimming situation already inspected above
            
                seq_hash = deep_seq_hash (VB, seq, e->seq_len, false);    // trimmed SEQ
                if (seq_hash != e->hash.seq) continue;

                segconf.deep_has_trimmed = true;

                qual_hash = deep_qual_hash (VB, qual, e->seq_len, false); // trimmed QUAL
                if (qual_hash == e->hash.qual) {
                    segconf.n_full_mch[qn]++;
                    return;
                }
                else
                    seq_qname_matches[qn] = true;
            }
        }
        
        if      (seq_qname_matches[0]) segconf.n_seq_qname_mch[0]++;
        else if (seq_qname_matches[1]) segconf.n_seq_qname_mch[1]++; 
        else                           segconf.n_no_mch++; // no match - likely bc this read was filtered out in the SAM file
    }
}

static inline uint64_t fastq_seg_deep_consume_unique_matching_ent (VBlockFASTQP vb, ZipZDeep *ent, DeepHash deep_hash,
                                                                   STRp(qname), STRp(seq), STRp(qual)) 
{
    START_TIMER;

    // copy current place (other threads might consume it in the mean time)
    union ZipZDeepPlace sam_place = { .place = ent->place };
    union ZipZDeepPlace contention_place; // vb and line_i of another thread that matches
    
    if (sam_place.consumed) {
        contention_place = sam_place;
        goto ent_consumed_by_another_thread;
    }

    union ZipZDeepPlace my_place = { .consumed = true,
                                     .vb_i     = vb->vblock_i,
                                     .line_i   = vb->line_i };

    // consume: atomicly set and verify that no other thread beat us to it. note: we will NOT modify if already consumed.
    bool i_consumed = __atomic_compare_exchange_n (&ent->place,      // set value at this pointer
                                                   &sam_place.place, // expected old value
                                                   my_place.place,   // set to this, but only if old value is as expected
                                                   false, __ATOMIC_RELAXED, __ATOMIC_RELAXED);

    if (!i_consumed) { // darn! another thread beat us to consuming
        contention_place = (union ZipZDeepPlace){ .place = ent->place }; // value of consuming FASTQ place as set by the other thread
        goto ent_consumed_by_another_thread;
    }

    ASSERT (sam_place.vb_i <= z_file->vb_start_deep_line.len32, "Expecting vb_i=%u <= z_file->vb_start_deep_line.len=%u", sam_place.vb_i, z_file->vb_start_deep_line.len32);

    // translate vb_i/line_i to a single, 0-based, file-wide number. This is a sequential counter of deepable lines - 
    // i.e. that are not SUPP/SEC, not monochar, and SEQ.len>0

    // MAX is max_uint64-1 bc we +1 this number before storing it. note: uint64_t despite DYN_INT limited to int64: we will interpret the number as uint
    #define MAX_DEEP_LINE (uint64_t)0xfffffffffffffffeULL  
    uint64_t vb_start_deep_line = *B64(z_file->vb_start_deep_line, sam_place.vb_i);
    uint64_t txt_deep_line_i = vb_start_deep_line + sam_place.line_i; // note: in ZIP, vb_start_deep_line is indexed by vb_i

    ASSERT (vb_start_deep_line <= MAX_DEEP_LINE && txt_deep_line_i <= MAX_DEEP_LINE &&
            txt_deep_line_i >= vb_start_deep_line, // txt_deep_line_i didn't overflow beyond the MAX of uint64_t
            "txt_deep_line_i=%"PRIu64" beyond MAX_DEEP_LINE=%"PRIu64": vb_start_deep_line=%"PRIu64" sam_place.line_i=%u", 
            txt_deep_line_i, MAX_DEEP_LINE, vb_start_deep_line, sam_place.line_i);

    COPY_TIMER (fastq_seg_deep_consume_unique_matching_ent);
    return txt_deep_line_i;

ent_consumed_by_another_thread:
    // case: this ent was already consumed by a previous FASTQ read. This can happen for two reasons:
    // 1) two different (QNAME,SEQ,QUAL) triplets have the same 96-bit deep_hash AND only one of the two 
    //    corresponding SAM lines actually exists in the SAM file. Therefore, both reads map to the single SAM
    //    line, and we don't know which is the true match. We have no choice but to abort.
    // 2) two identical (QNAME,SEQ,QUAL), and only one SAM line exists. conceivably, this could happen for example
    //    if SEQ of the two mates sharing the QNAME is trivial but non-monochar (eg an STR), and QUAL is mono-value.
    //    We can't distinguish between this case and case 1, so we need to abort as well.
    // TO DO: recompress normally by executing another instance of genozip
    ABORTINP ("Deep: We hit a rare edge case: two FASTQ reads: current line %s and line FQ\?\?/%u/%u (vb_size=%s)\n"
              "map to the same SAM line (deep_hash=%u,%u,%u qtype=%s no_qual=%s).\n"
              "Please kindly report it to " EMAIL_SUPPORT ". Work around: compress without using --deep.\n"
              "This_line:\nQNAME=\"%.*s\"\nSEQ=\"%.*s\"\nQUAL=\"%.*s\"",
              LN_NAME, contention_place.vb_i, contention_place.line_i, str_size (segconf.vb_size).s, 
              DEEPHASHf(deep_hash), qtype_name(segconf.deep_qtype), TF(segconf.deep_no_qual), 
              STRf(qname), STRf(seq), STRf(qual));
    
    COPY_TIMER (fastq_seg_deep_consume_unique_matching_ent);
    return 0;
}

static int64_t fastq_get_pair_deep_value (VBlockFASTQP vb, ContextP ctx)
{
    ASSERT (ctx->localR1.next < ctx->localR1.len32, "%s: not enough data DEEP.localR1 (len=%u)", LN_NAME, ctx->localR1.len32); 

    switch (ctx->pair_ltype) {
        case LT_INT64:  // note: if UINT64 will appear as INT64 as it was DYN_INT. we treat it as UINT64.
        case LT_UINT64: return *B64(ctx->localR1, ctx->localR1.next++); 
        case LT_UINT32: return *B32(ctx->localR1, ctx->localR1.next++); 
        case LT_UINT16: return *B16(ctx->localR1, ctx->localR1.next++); 
        case LT_UINT8:  return *B8 (ctx->localR1, ctx->localR1.next++); 
        default:   
            ABORT ("Unexpected pair_ltype=%u", ctx->pair_ltype);
    }
}

// called to mark entries as dup in case multiple SAM alignments have the same hash
static void fastq_seg_deep_dup_detected (VBlockFASTQP vb, int by, ZipZDeep *dup1, ZipZDeep *dup2, uint32_t index_hash)
{
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    // first encounter with this dup - traverse the rest of the linked list looking for dups
    // note: we don't bother about atomicity of visibility of "dup" is just used for optimizing the loop and for stats
    if (!dup1->dup) {
        dup1->dup = dup2->dup = true;

        for (uint32_t ent_i = dup2->next[by]; ent_i != NO_NEXT; ent_i = deep_ents[ent_i].next[by]) 
            if (!memcmp (&deep_ents[ent_i].hash, &dup1->hash, sizeof (dup1->hash))) 
                deep_ents[ent_i].dup = true;                
    }

    vb->deep_stats[by==BY_SEQ ? NDP_MULTI_MATCH : NDP_MULTI_TRIMMED]++;

    if (flag.show_deep) 
        iprintf ("%s: Two %s alignments with same hashes (%s), skipping: this=[hash=%u,%u,%u sam_vb_i=%u sam_line_i=%u] prev=[%u,%u,%u sam_vb_i=%u sam_line_i=%u] [deep_no_qual=%s deep_qtype=%s]\n", 
                 VB_NAME, z_dt_name(), by_names[by], DEEPHASHf(dup2->hash), dup2->vb_i, dup2->line_i, 
                 DEEPHASHf(dup1->hash), dup1->vb_i, dup1->line_i,
                 TF(segconf.deep_no_qual), qtype_name (segconf.deep_qtype));            
}

void fastq_seg_deep (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(qname), STRp(qname2), STRp(seq), STRp(qual), 
                     bool *deep_qname, bool *deep_seq, bool *deep_qual, // out - set to true, or left unchanged
                     uint32_t *uncanonical_suffix_len) // out - suffix of deeped QNAME beyond canonical, i.e. not copied from Deep
{
    START_TIMER;

    decl_ctx (FASTQ_DEEP);

    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    if (flag.deep) { // skip tests if only --show-deep without --deep
        vb->deep_stats[NDP_FQ_READS]++;
        
        if (!deep_ents_len) {
            vb->deep_stats[NDP_NO_ENTS]++;
            goto done;
        }

        if (segconf.running) {
            fastq_deep_seg_segconf (vb, STRa(qname), STRa(qname2), STRa(seq), STRa(qual));
            goto no_match;
        }

        // we don't Deep reads monochar SEQ - too much contention
        if (dl->monochar) {
            vb->deep_stats[NDP_MONOSEQ]++;
            goto no_match; 
        }

        // If QNONE (i.e. we're hashing SEQ and QUAL only), we don't Deep reads with monoqual, as it has
        // an elevated chance that two reads will have entirely the same hash (in which case Deep will fail)
        if (segconf.deep_qtype == QNONE && str_is_monochar (STRa(qual))) {
            vb->deep_stats[NDP_MONOQUAL]++;
            goto no_match; 
        }
    }

    DeepHash deep_hash = { .qname = (segconf.deep_qtype == QNAME1) ? deep_qname_hash (QNAME1, STRa(qname),  uncanonical_suffix_len)
                                  : (segconf.deep_qtype == QNAME2) ? deep_qname_hash (QNAME2, STRa(qname2), uncanonical_suffix_len)
                                  :                                  0,
                           .seq   =                             deep_seq_hash (VB, STRa(seq), false),
                           .qual  = segconf.deep_no_qual  ? 0 : deep_qual_hash (VB, STRa(qual), false) };

    if (flag.show_deep == 2) {
        if (deephash_issame (flag.debug_deep_hash, deep_hash)) 
            iprintf ("%s Found deep_hash=%u,%u,%u\nQNAME=\"%.*s\"\nSEQ=\"%.*s\"\nQUAL=\"%.*s\"\n",
                     LN_NAME, DEEPHASHf(deep_hash), STRf(qname), STRf(seq), STRf(qual));

        if (!flag.deep) goto done; // only debug
    }

    // find matching entry - and make sure there is only one matching entry. 
    // case: If there are multiple matching entries - we are not sure which SAM line this FASTQ read relates to, so we don't Deep.
    // case: All matching entries are already used by previous FASTQ lines - so we have multiple FASTQ lines claiming the same
    //       SAM entry - since we can't undo the previous lines at this point - we fail the execution and re-compress without deep
    ZipZDeep *matching_ent = NULL;
    DeepStatsFastq match_type;

    for (uint32_t ent_i = *B32 (z_file->deep_index_by[BY_SEQ], deep_hash.seq & bitmask32 (num_hash_bits)); 
         ent_i != NO_NEXT; ent_i = deep_ents[ent_i].next[BY_SEQ]) {

        ZipZDeep *e = &deep_ents[ent_i];
        
        // case: SEQ matches, and QUAL and QNAME match if they are required to match
        bool qual_matches=false, qname_matches=false;
        
        if (e->hash.seq == deep_hash.seq &&  // SEQ matches
            (segconf.deep_no_qual || !e->hash.qual || (qual_matches = (deep_hash.qual == e->hash.qual))) &&  // QUAL matches (or not relevant)
            (segconf.deep_qtype == QNONE || (qname_matches = (deep_hash.qname == e->hash.qname)))) {  // QNAME matches (or not relevant)
            
            if (e->dup) {
                vb->deep_stats[NDP_MULTI_MATCH]++; // already detected before as a dup entry
                goto no_match;
            }

            // case: two or more SAM lines entries match this read according to the hashes - we don't know which is the true one, so we seg without deep.
            if (matching_ent) {
                fastq_seg_deep_dup_detected (vb, BY_SEQ, matching_ent, e, deep_hash.seq);
                goto no_match;
            }

            matching_ent = e;
            *deep_seq    = true;
            *deep_qual   = qual_matches;
            *deep_qname  = qname_matches;
            match_type   = NDP_DEEPABLE;
        }
    }

    // case: no matching searching by SEQ - perhaps SAM SEQ/QUAL are trimmed - search by QNAME (but skip search if no trimming was found by segconf)
    if (!matching_ent && segconf.deep_has_trimmed)
        for (uint32_t ent_i = *B32 (z_file->deep_index_by[BY_QNAME], deep_hash.qname & bitmask32 (num_hash_bits)); 
            ent_i != NO_NEXT; ent_i = deep_ents[ent_i].next[BY_QNAME]) {

            ZipZDeep *e = &deep_ents[ent_i];
            if (e->hash.qname != deep_hash.qname || // possible, because hash uses less bits that e->hash.qname
                e->seq_len >= seq_len) continue;    // definitely not a match if FASTQ seq_len is shorter than SAM's and also not equal seq_len, as that would have been a no-trimming situation already inspected above

            // case: SEQ matches, and QUAL and QNAME match if they are required to match
            if (deep_seq_hash (VB, seq, e->seq_len, false) != e->hash.seq) // trimmed SEQ
                continue;

            *deep_seq = *deep_qname = true;

            *deep_qual = !segconf.deep_no_qual && e->hash.qual && 
                         (deep_qual_hash (VB, qual, e->seq_len, false) == e->hash.qual); // trimmed QUAL match

            if (e->dup) {
                vb->deep_stats[NDP_MULTI_TRIMMED]++; // already detected before as a dup entry
                goto no_match;
            }

            if (matching_ent) {
                fastq_seg_deep_dup_detected (vb, BY_QNAME, matching_ent, e, deep_hash.qname);
                goto no_match;
            }

            matching_ent = e;
            match_type   = NDP_DEEPABLE_TRIM;
        }

    // single match ent
    uint64_t deep_value;
    if (matching_ent) {
        deep_value = 1 + fastq_seg_deep_consume_unique_matching_ent (vb, matching_ent, deep_hash, STRa(qname), STRa(seq), STRa(qual)); // +1 as 0 means "no match"
        vb->deep_stats[match_type]++;

        dl->sam_seq_len = matching_ent->seq_len;
    }

    // case: no matching ent
    else {
        vb->deep_stats[NDP_NO_MATCH]++;

        if (flag.show_deep) {
            static int count = 0;
            if (__atomic_fetch_add (&count, (int)1, __ATOMIC_RELAXED) < 20) {
                rom once = "";
                DO_ONCE once = "\n\n(Showing first ~20 no-matches)\n";
                iprintf ("%s%s: no matching SAM line (has_match=%s, deep_hash=%u,%u,%u \nQNAME=\"%.*s\"\nSEQ=  \"%.*s\"\nQUAL= \"%.*s\"\n", 
                         once, LN_NAME, TF(!!matching_ent), DEEPHASHf(deep_hash), STRf(qname), STRf(seq), STRf(qual));        
            }
        }

        no_match: 
        deep_value = 0;
        *deep_seq = *deep_qual = *deep_qname = false; // reset
    }

    if (flag.pair == NOT_PAIRED || flag.pair == PAIR_R1)
        seg_add_to_local_resizable (VB, ctx, deep_value, 0);

    else { // PAIR_R2
        uint64_t pair_1_deep_value = fastq_get_pair_deep_value (vb, ctx); // consume whether or not used

        // delta carefully to avoid overflow if values are large uint64's
        uint64_t abs_delta = (pair_1_deep_value > deep_value) ? (pair_1_deep_value - deep_value)
                                                              : (deep_value - pair_1_deep_value);                                       
        if (abs_delta > 127) { 
            seg_add_to_local_resizable (VB, ctx, deep_value, 0);
            seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_set_deep, '0' }, 3, FASTQ_DEEP, 0);
        }   

        else {
            int8_t delta8 = (pair_1_deep_value > deep_value) ? -(int8_t)abs_delta : (int8_t)abs_delta;
            seg_integer_fixed (VB, CTX(FASTQ_DEEP_DELTA), &delta8, false, 0);
            seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_set_deep, '1' }, 3, FASTQ_DEEP, 0);
        }
    }

done:
    COPY_TIMER (fastq_seg_deep);
}

void fastq_deep_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(seq), ContextP bitmap_ctx, ContextP nonref_ctx)
{
    uint32_t trim_len = seq_len - dl->sam_seq_len;

    if (!trim_len)
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_deep_copy_SEQ }, 2, bitmap_ctx, seq_len); 
    
    else {
        SNIPi2 (SNIP_SPECIAL, FASTQ_SPECIAL_deep_copy_SEQ, trim_len);
        seg_by_ctx (VB, STRa(snip), bitmap_ctx, dl->sam_seq_len); 

        seg_add_to_local_fixed (VB, nonref_ctx, &seq[dl->sam_seq_len], trim_len, LOOKUP_NONE, trim_len);
    }
}

void fastq_deep_seg_QUAL (VBlockFASTQP vb, ZipDataLineFASTQ *dl, ContextP qual_ctx, uint32_t qual_len)
{
    uint32_t trim_len = qual_len - dl->sam_seq_len;

    seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_deep_copy_QUAL }, 2, qual_ctx, qual_len);
    
    if (trim_len) 
        qual_ctx->local.len32 += trim_len; // fastq_zip_qual will get the trimmed part of the qual string
    else
        dl->dont_compress_QUAL = true;
}

// used to seg QNAME or QNAME2 as a copy from deep, possibly with a suffix
void fastq_deep_seg_QNAME (VBlockFASTQP vb, Did did_i, STRp(qname), uint32_t uncanonical_suffix_len, unsigned add_bytes)
{
    // case: qname has an uncanonical suffix - seg a container consisting on the canonical and uncanonical substrings
    if (uncanonical_suffix_len) {
        if (did_i == FASTQ_QNAME) seg_by_did (VB, STRa(con_decanonize1_snip), did_i, 0);
        else /*QNAME2*/           seg_by_did (VB, STRa(con_decanonize2_snip), did_i, 0);
        
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_deep_copy_QNAME }, 2, did_i + 1, add_bytes - uncanonical_suffix_len);
        seg_by_did (VB, qname + qname_len - uncanonical_suffix_len, uncanonical_suffix_len, did_i + 2, uncanonical_suffix_len);
    }

    // case: entire qname is canonical - seg as "copy from deep"
    else
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_deep_copy_QNAME }, 2, did_i, add_bytes);
}

//----------------
// PIZ
//----------------

// PIZ: called from compute threads of FASTQ component - after decompressing contexts but before reconstructing
void fastq_deep_piz_wait_for_deep_data (void)
{
    // busy-wait until sam_piz_deep_finalize_ents finalizes the deep data
    while (!z_file->vb_start_deep_line.param) 
        usleep (100000); // 100ms

    // flag.no_writer (eg --test) with gencomp files is forced in-order for SAM bc digest is calculated by writer for gencomp files,
    // but can be out of order in the FQ part as if no_writer=true FASTQ VBs have needs_write=false (see writer_init_vb_info)
    if (flag.no_writer)
        piz_allow_out_of_order(); 
}

// PIZ: binary search for (vb_idx, vb_deepable_line_i) from file-wide txt_deepable_line_i
static void txt_line_to_vb_line (uint64_t txt_deepable_line_i, 
                                 uint32_t first_vb_idx, uint32_t last_vb_idx,   // for binary search
                                 uint32_t *vb_idx, uint32_t *vb_deepable_line_i) // out
{
    if (first_vb_idx > last_vb_idx) {
        *vb_idx = last_vb_idx;
        *vb_deepable_line_i = txt_deepable_line_i - *B64(z_file->vb_start_deep_line, last_vb_idx);
        return;
    }
        
    VBIType mid_vb_idx = (first_vb_idx + last_vb_idx) / 2;

    uint64_t mid_vb_first_line = *B64(z_file->vb_start_deep_line, mid_vb_idx);

    // case: txt_deepable_line_i is in a VB lower than mid_vb_idx
    if (txt_deepable_line_i < mid_vb_first_line) 
        txt_line_to_vb_line (txt_deepable_line_i, first_vb_idx, mid_vb_idx-1, vb_idx, vb_deepable_line_i);

    // case: txt_deepable_line_i is in mid_vb_idx or higher VB. note: if it is in mid_vb_idx, the stop condition will trigger in the recursive call 
    else 
        txt_line_to_vb_line (txt_deepable_line_i, mid_vb_idx+1, last_vb_idx, vb_idx, vb_deepable_line_i);
}

SPECIAL_RECONSTRUCTOR_DT (fastq_special_set_deep)
{
    START_TIMER;

    VBlockFASTQP vb = (VBlockFASTQP)vb_;

    ASSERTNOTEMPTY (z_file->vb_start_deep_line);
    uint64_t txt_deepable_line_i;

    uint64_t pair_1_deep_value = vb->pair_vb_i ? fastq_get_pair_deep_value (vb, ctx) : 0; // consume whether or not used

    if (snip[0] == '0') // no delta
        txt_deepable_line_i = reconstruct_from_local_int (VB, ctx, 0, RECON_OFF);

    else {
        int8_t delta = reconstruct_from_local_int (VB, CTX(FASTQ_DEEP_DELTA), 0, RECON_OFF);

        txt_deepable_line_i = (delta > 0) ? (pair_1_deep_value + delta) : (pair_1_deep_value - (-delta)); // careful with type coverstions and very large uint64 overflow 
    }
    
    // case: this FQ line doesn't match any SAM line - txt_deepable_line_i==0 and we didn't seg against it
    if (!txt_deepable_line_i)
        return NO_NEW_VALUE;

    uint32_t vb_idx;             // SAM VB where this deepable line resides - not counting DEPN VBs, so this is NOT vblock_i!
    uint32_t vb_deepable_line_i; // deepable_line_i within vb_idx - not counting non-deepable lines, so this is NOT line_i!
    txt_line_to_vb_line (txt_deepable_line_i - 1, 0, z_file->deep_index.len32 - 1, &vb_idx, &vb_deepable_line_i);
    
    // get this VB's data
    BufferP deep_ents  = B(Buffer, z_file->deep_ents,  vb_idx);
    BufferP deep_index = B(Buffer, z_file->deep_index, vb_idx);

    ASSPIZ (vb_deepable_line_i < deep_index->len32, "Expecting vb_deepable_line_i=%u < deep_index->len=%u (txt_deepable_line_i=%"PRIu64" vb_idx=%u)", 
            vb_deepable_line_i, deep_index->len32, txt_deepable_line_i, vb_idx);

    uint32_t ents_index = *B32 (*deep_index, vb_deepable_line_i);

    new_value->p = B8 (*deep_ents, ents_index);

    COPY_TIMER (fastq_special_set_deep);
    return HAS_NEW_VALUE; 
}

SPECIAL_RECONSTRUCTOR (fastq_special_deep_copy_QNAME)
{
    START_TIMER;

    uint8_t *deep_ent = CTX(FASTQ_DEEP)->last_value.p;
    ASSERTNOTNULL (deep_ent);

    PizZDeepFlags f = *(PizZDeepFlags *)deep_ent;
    int qname_len = deep_ent[1];
    int prfx_len  = deep_ent[2];
    ASSPIZ (qname_len >= prfx_len, "Expecting qname_len=%d >= prfx_len=%d", qname_len, prfx_len);

    STRli (suffix, qname_len - prfx_len);

    RECONSTRUCT (z_file->master_qname, prfx_len);

    int comp_len = 0;
    if (f.is_qname_comp) {
        comp_len = huffman_decompress (z_file->qname_huf, &deep_ent[3], (uint8_t *)suffix, suffix_len);

        for (int i=0; i < suffix_len; i++)
            suffix[i] ^= z_file->master_qname[prfx_len + i];
    }

    if (reconstruct)
        RECONSTRUCT (f.is_qname_comp ? suffix : (rom)&deep_ent[3], suffix_len);

    // update qual_len to be entire amount of data consumed - inc "prfx_len" and compressed (or not) suffix
    deep_ent[1] = 1/*prfx_len*/ + (f.is_qname_comp ? comp_len : suffix_len);

    COPY_TIMER (fastq_special_deep_copy_QNAME);
    return NO_NEW_VALUE;    
}

SPECIAL_RECONSTRUCTOR (fastq_special_deep_copy_SEQ)
{
    START_TIMER;

    if (!reconstruct) return NO_NEW_VALUE;

    uint32_t trim_len = snip_len ? atoi (snip) : 0;

    uint8_t *deep_ent = CTX(FASTQ_DEEP)->last_value.p;
    ASSERTNOTNULL (deep_ent);

    PizZDeepFlags f = *(PizZDeepFlags *)deep_ent;

    deep_ent += (segconf.deep_qtype == QNONE || flag.seq_only || flag.qual_only) 
        ? 1/*skip flags*/ : (2 + deep_ent[1])/*skip flags, qname_len and qname*/; 

    uint32_t deep_seq_len = f.is_long_seq ? GET_UINT32 (deep_ent) : *deep_ent;
    
    vb->seq_len = deep_seq_len + trim_len; 
    
    ASSERT (vb->seq_len <= vb->longest_seq_len, "%s: unexpectedly seq_len(as read from deep_ent)=%u > vb->longest_seq_len=%u (deep_seq_len=%u trim_len=%u)",
            LN_NAME, vb->seq_len, vb->longest_seq_len, deep_seq_len,trim_len); // sanity check

    deep_ent += f.is_long_seq ? 4 : 1;  

    #define ENCODING(x) (f.seq_encoding == ZDEEP_SEQ_##x)

    // case: SEQ is compressed (bc it contains non-ACGT)
    if (ENCODING(COMPRESSED) || ENCODING(COMPRESSED_LONG_LEN)) {
        uint32_t comp_len = (ENCODING(COMPRESSED) ? *deep_ent : GET_UINT32 (deep_ent));
        deep_ent += (ENCODING(COMPRESSED) ? 1 : 4);

        unsigned out_len = deep_seq_len;        
        ASSPIZ (arith_uncompress_to (vb, deep_ent, comp_len, (uint8_t*)BAFTtxt, &out_len) && out_len == deep_seq_len,
                "Failed arith_uncompress_to SEQ copied from SAM: compressed_len=%u uncompressed_len=%u (expected: %u)", 
                comp_len, out_len, deep_seq_len);
    }

    // case: SEQ is packed
    else if (ENCODING(PACKED)) {
        // case: shortish reads - can fit on stack 
        #define SHORTISH_READ_SIZE 1024
        if ((ROUNDUP32(deep_seq_len) / 4 <= SHORTISH_READ_SIZE)) {
            uint8_t data[SHORTISH_READ_SIZE]; // not too long, so code path gets mileage
            memcpy (data, deep_ent, ROUNDUP4(deep_seq_len) / 4); // copy to word-align
            
            Bits pack = bits_init (deep_seq_len * 2, data, sizeof(data), false);
            bits_2bit_to_ACGT (BAFTtxt, &pack, 0, deep_seq_len);
        }

        // case: too long for stack - allocate on heap
        else {
            ASSERTNOTINUSE (vb->scratch);
            Bits pack = bits_alloc (deep_seq_len*2, true);
            memcpy (pack.words, deep_ent, ROUNDUP4(deep_seq_len) / 4);

            bits_2bit_to_ACGT (BAFTtxt, &pack, 0, deep_seq_len);
            bits_free (&pack);
        }
    }

    else { // ZDEEP_SEQ_PERFECT, ZDEEP_SEQ_PERFECT_REV, ZDEEP_SEQ_MIS_1, ZDEEP_SEQ_MIS_1_REV      
        PizDeepSeq ds = *(PizDeepSeq *)deep_ent; // not nessesarily word-aligned, but I think its ok, since this is just a memcpy, not loading integers
        deep_ent += (ENCODING(PERFECT) || ENCODING(PERFECT_REV)) ? 4 : 6;

        const Bits *genome=NULL;
        ref_get_genome (gref, &genome, NULL, NULL);
        
        BASE_ITER_INIT (genome, ds.gpos, deep_seq_len, (ENCODING(PERFECT) || ENCODING(MIS_1)));

        char *next = BAFTtxt;
        decl_acgt_decode;
        for (uint32_t i=0; i < deep_seq_len; i++) 
            *next++ = acgt_decode(BASE_NEXT);                    

        if (ENCODING(MIS_1)) 
            *(BAFTtxt + ds.mismatch_offset) = ds.mismatch_base; // reconstruct mismatch 

        else if (ENCODING(MIS_1_REV))
            *(BAFTtxt + deep_seq_len-1 - ds.mismatch_offset) = COMPLEM[(int)ds.mismatch_base]; // reconstruct mismatch 
    }

    Ltxt += deep_seq_len;

    if (trim_len)
        reconstruct_from_local_sequence (vb, CTX(SAM_NONREF), trim_len, reconstruct);

    // case we are pair-2: advance pair-1 SQBITMAP iterator, and if pair-1 is aligned - also its GPOS iterator
    if (VB_FASTQ->pair_vb_i/*we are R2*/ && fastq_piz_R1_test_aligned (VB_FASTQ))
        CTX(FASTQ_GPOS)->localR1.next++; // gpos_ctx->localR1.next is an iterator for both gpos and strand

    COPY_TIMER (fastq_special_deep_copy_SEQ);
    return NO_NEW_VALUE;    
}

SPECIAL_RECONSTRUCTOR (fastq_special_deep_copy_QUAL)
{
    START_TIMER;

    uint8_t *deep_ent = CTX(FASTQ_DEEP)->last_value.p;
    ASSERTNOTNULL (deep_ent);

    PizZDeepFlags f = *(PizZDeepFlags *)deep_ent;

    deep_ent += (segconf.deep_qtype == QNONE || flag.seq_only || flag.qual_only) 
        ? 1/*skip flags*/ : (2 + deep_ent[1])/*skip flags, qname_len and qname*/; 

    uint32_t deep_seq_len = f.is_long_seq ? GET_UINT32 (deep_ent) : *deep_ent;

    if (reconstruct) {
        deep_ent += f.is_long_seq ? 4 : 1;
        
        // skip seq
        if (!flag.qual_only) // in qual_only, sam_piz_deep_add_seq didn't store the SEQ
            switch (f.seq_encoding) {
                case ZDEEP_SEQ_COMPRESSED          : deep_ent += 1 + *deep_ent;              break;
                case ZDEEP_SEQ_COMPRESSED_LONG_LEN : deep_ent += 4 + GET_UINT32 (deep_ent);  break;
                case ZDEEP_SEQ_PACKED              : deep_ent += ROUNDUP4(deep_seq_len) / 4; break;
                case ZDEEP_SEQ_PERFECT             : deep_ent += 4;                          break;
                case ZDEEP_SEQ_PERFECT_REV         : deep_ent += 4;                          break;
                case ZDEEP_SEQ_MIS_1               : deep_ent += 6;                          break;  
                case ZDEEP_SEQ_MIS_1_REV           : deep_ent += 6;                          break;  
                default: ABORT0 ("Unknown SEQ encoding");
            }

        uint32_t qual_comp_len = f.is_long_qual_comp ? GET_UINT32 (deep_ent) : *deep_ent;
        deep_ent += (f.is_long_qual_comp ? 4 : 1);

        ASSPIZ0 (qual_comp_len, "qual_comp_len=0 but expecting it to be >0 because QUAL of this line is not monochar");

        unsigned out_len = deep_seq_len;        
        ASSPIZ (arith_uncompress_to (vb, deep_ent, qual_comp_len, (uint8_t *)BAFTtxt, &out_len) && out_len == deep_seq_len,
                "Failed arith_uncompress_to QUAL copied from SAM: compressed_len=%u uncompressed_len=%u (expected: %u)", 
                qual_comp_len, out_len, deep_seq_len);

        Ltxt += out_len;
    }

    if (deep_seq_len < vb->seq_len) // we have trimmed bases
        switch (ctx->ltype) { // the relevant subset of ltypes from reconstruct_from_ctx_do
            case LT_CODEC: 
                codec_args[ctx->lcodec].reconstruct (VB, ctx->lcodec, ctx, vb->seq_len - deep_seq_len); break;

            case LT_SEQUENCE: 
                reconstruct_from_local_sequence (VB, ctx, vb->seq_len - deep_seq_len, reconstruct); break;

            default: ASSPIZ (false, "Invalid ltype=%s for %s", lt_name (ctx->ltype), ctx->tag_name);
        }

    COPY_TIMER (fastq_special_deep_copy_QUAL);
    return NO_NEW_VALUE;
}
