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
#include "htscodecs/arith_dynamic.h"

//----------------
// ZIP
//----------------

void fastq_deep_zip_initialize (void)
{
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

    if (flag.pair == PAIR_R1 || flag.pair == NOT_PAIRED_END) 
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_set_deep, '0' }, 3, FASTQ_DEEP, 0);  // all-the-same for FASTQ_DEEP

    else { // pair-2
        ctx_set_ltype (VB, LT_INT8, FASTQ_DEEP_DELTA, DID_EOL); 
        ctx_consolidate_stats (VB, FASTQ_DEEP, FASTQ_DEEP_DELTA, DID_EOL);
    }
}

void fastq_zip_deep_show_entries_stats (void) 
{
    static rom names[] = NO_DEEP_NAMES;

    uint64_t total = z_file->deep_stats[NDP_FQ_READS];
    iprint0 ("\nFASTQ reads breakdown by deepability:\n");
    
    for (int i=0; i < NUM_DEEP_STATS; i++) 
        if (z_file->deep_stats[i])
            iprintf ("%-11.11s: %"PRIu64" (%.1f%%)\n", names[i], z_file->deep_stats[i], 100.0 * (double)z_file->deep_stats[i] / (double)total);

    ARRAY (uint32_t, hash_table, z_file->deep_hash);
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    iprintf ("\nz_file.deep_hash.len=%"PRIu64" z_file.deep_ents.len=%"PRIu64"\n", z_file->deep_hash.len, z_file->deep_ents.len);

    #define NUM_LENS 16384
    #define PRINT_NUN_LENS 64
    uint64_t used=0;
    uint32_t lens[NUM_LENS] = {}, longer=0; // each entry i contains number linked lists which are of length i; last entry i - linked lists of i or longer

    for (uint64_t hash=0; hash < hash_table_len; hash++) {

        uint32_t this_len=0;
        for (uint32_t ent_i = hash_table[hash]; ent_i != NO_NEXT; ent_i = deep_ents[ent_i].next) 
            this_len++;
        
        if (this_len) used++;

        lens[MIN_(this_len,NUM_LENS-1)]++; // last item reported contains all NUM_LENS-1 *or longer*
    
        if (this_len > PRINT_NUN_LENS) longer++;
    }

    iprintf ("\nDeep hash entries used: %"PRIu64"/%"PRIu64" (%.1f%%)\n", used, hash_table_len, 100.0 * (double)used / (double)hash_table_len);

    iprint0 ("Deep hash linked list length histogram:\n");
    uint64_t total_ents_traversed = 0; // total ents traversed on linked lists when segging FASTQ reads
    for (int this_len=0; this_len < NUM_LENS; this_len++) {
        if (!lens[this_len]) continue;

        if (this_len < PRINT_NUN_LENS)
            iprintf ("%3u: %u\n", this_len, lens[this_len]);

        // our algorithm requires traversal of entire linked list for each entry. calculate the length of the list
        // traversed for the average read (simplifying assumption: there are the same )
        total_ents_traversed += this_len        // number of entries that will access this linked list
                              * this_len        // length of linked list traversed for each entry accessing this list
                              * lens[this_len]; // number of linked lists of this length
    }

    if (longer) iprintf ("longer: %u\n", longer);

    // this is approximate: assuming # of ents hashed is equal to the number of reads
    iprintf ("\nFASTQ SEQ: Linked list length traversed on average: %.2f\n", (double)total_ents_traversed / (double)deep_ents_len);

    uint64_t count_not_consumed=0;
    for (uint64_t ent_i=0; ent_i < deep_ents_len; ent_i++)
        if (!deep_ents[ent_i].consumed) count_not_consumed++;

    iprintf ("\nNumber of %s deepable alignments not consumed by FASTQ reads: %"PRIu64" (showing first %u)\n", 
             dt_name(z_file->data_type), count_not_consumed, (int)MIN_(count_not_consumed, 20));

    if (count_not_consumed)
        for (uint64_t ent_i=0, count=0; ent_i < deep_ents_len && count < 20; ent_i++)
            if (!deep_ents[ent_i].consumed) {
                iprintf ("sam_vb=%u line_i=%u hash=%u,%u,%u\n", deep_ents[ent_i].vb_i, deep_ents[ent_i].line_i, DEEPHASHf(deep_ents[ent_i].hash));
                count++;
            }
}

void fastq_deep_seg_finalize_segconf (uint32_t n_lines)
{
    if (flag.show_deep) 
        iprintf ("segconf: n_lines=%u n_full_mch=%u n_seq_qname_mch=%u n_seq_qual_mch=%u n_seq_only=%u n_no_match=%u\n",
                 n_lines, segconf.n_full_mch, segconf.n_seq_qname_mch, segconf.n_seq_qual_mch, segconf.n_seq_mch,segconf.n_no_mch);

    unsigned n_mch = n_lines - segconf.n_no_mch; // exclude lines that don't match any SAM line - perhaps they were filtered out and excluded from the SAM file 
    unsigned threashold = MAX_(n_mch/2, 1);      // at least half of the lines need this level of matching

    // test: at least half of the reads (excluding reads that had no matching SEQ in the SAM) have a matching QNAME
    if (segconf.n_seq_qname_mch + segconf.n_full_mch < threashold) 
        segconf.deep_no_qname= true;

    // test: likewise, matching QUAL
    if (segconf.n_seq_qual_mch + segconf.n_full_mch < threashold)  
        segconf.deep_no_qual = true;

    // we need at least 2 of the 3 (QNANE,SEQ,QUAL) to be comparable between FASTQ and SAM to have a total of 64 bit hash (32 from each)    
    ASSINP (!segconf.deep_no_qname || !segconf.deep_no_qual, "Error: cannot use --deep with this file: based on testing the first %u reads of %s, only %u reads "
            "had a matching alignment in the SAM/BAM file matching QNAME, SEQ, QUAL, %u matching QNAME and SEQ and %u matching SEQ and QUAL. "
            "This is below the threashold need for meaningful Deep compression. You may compress these files without --deep.",
            n_lines, txt_name, segconf.n_full_mch, segconf.n_seq_qname_mch, segconf.n_seq_qual_mch);
}

static void fastq_deep_seg_segconf (VBlockFASTQP vb, STRp(qname), STRp(seq), STRp(qual))
{
    ARRAY (uint32_t, hash_table, z_file->deep_hash);
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    uint32_t qname_hash = deep_qname_hash (STRa(qname));
    uint32_t seq_hash   = deep_seq_hash (VB, STRa(seq), false);
    uint32_t qual_hash  = deep_qual_hash (VB, STRa(qual), false);

    uint32_t hash = seq_hash & bitmask32 (num_hash_bits);

    bool full_matches=false, seq_qual_matches=false, seq_qname_matches=false, seq_matches=false;
    
    for (uint32_t ent_i = hash_table[hash]; ent_i != NO_NEXT; ent_i = deep_ents[ent_i].next) {
        ZipZDeep *e = &deep_ents[ent_i];
        
        if (e->hash.seq == seq_hash) {
            if (e->hash.qname == qname_hash) {
                if (e->hash.qual == qual_hash)
                    full_matches = true;
                else
                    seq_qname_matches = true;
            }
            else if (e->hash.qual == qual_hash)
                seq_qual_matches = true;
            else
                seq_matches = true;
        }
    }

    // each segconf read increments exactly one category
    if      (full_matches)      segconf.n_full_mch++;      // at least one match was full
    else if (seq_qname_matches) segconf.n_seq_qname_mch++; // if one match is seq+qname and another is seq+qual, then count the seq+qname 
    else if (seq_qual_matches)  segconf.n_seq_qual_mch++;
    else if (seq_matches)       segconf.n_seq_mch++;       // no match for qname and qual, only seq alone
    else                        segconf.n_no_mch++;        // no match - likely bc this read was filtered out in the SAM file
}

static inline uint64_t fastq_seg_deep_consume_unique_matching_ent (VBlockFASTQP vb, ZipZDeep *ent, DeepHash deep_hash,
                                                                   STRp(qname), STRp(seq), STRp(qual)) 
{
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
    ABORTINP ("Deep: We hit an extremely rare edge case: two FASTQ reads: %s and vb_i=%u line_i=%u (vb_size=%s)\n"
              "map to the same SAM line (deep_hash=%u,%u,%u). This is a theorical edge case that has never before\n"
              "been observed, please kindly report it to support@genozip.com. Work around: compress without using --deep.\n"
              "This_line:\nQNAME=\"%.*s\"\nSEQ=\"%.*s\"\nQUAL=\"%.*s\"",
              LN_NAME, contention_place.vb_i, contention_place.line_i, str_size (segconf.vb_size).s, 
              DEEPHASHf(deep_hash), STRf(qname), STRf(seq), STRf(qual));
    return 0;
}

static int64_t fastq_get_pair_deep_value (VBlockFASTQP vb, ContextP ctx)
{
    ASSERT (ctx->localR1.next < ctx->localR1.len32, "%s: not enough data DEEP.localR1 (len=%u)", LN_NAME, ctx->localR1.len32); 

    switch (ctx->pair_ltype) {
        case LT_INT64:  // note: if UINT64 will appear as INT64 is it was DYN_INT. we treat it as UINT64.
        case LT_UINT64: return *B64(ctx->localR1, ctx->localR1.next++); 
        case LT_UINT32: return *B32(ctx->localR1, ctx->localR1.next++); 
        case LT_UINT16: return *B16(ctx->localR1, ctx->localR1.next++); 
        case LT_UINT8:  return *B8 (ctx->localR1, ctx->localR1.next++); 
        default:   
            ABORT_R ("Unexpected pair_ltype=%u", ctx->pair_ltype);
    }
}

void fastq_seg_deep (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(qname), STRp(seq), STRp(qual), 
                     bool *deep_desc, bool *deep_seq, bool *deep_qual) // out - set to true, or left unchanged
{
    decl_ctx (FASTQ_DEEP);

    ARRAY (uint32_t, hash_table, z_file->deep_hash);
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    if (flag.deep) { // skip tests if only --show-deep without --deep
        vb->deep_stats[NDP_FQ_READS]++;
        
        if (!deep_ents_len) {
            vb->deep_stats[NDP_NO_ENTS]++;
            return;
        }

        if (segconf.running) {
            fastq_deep_seg_segconf (vb, STRa(qname), STRa(seq), STRa(qual));
            goto no_match;
        }

        if (dl->monochar) {
            vb->deep_stats[NDP_MONOCHAR]++;
            goto no_match; // we don't deep monochar reads - too much contention
        }
    }

    DeepHash deep_hash = { .qname = segconf.deep_no_qname ? 0 : deep_qname_hash (STRa(qname)),
                           .seq   =                             deep_seq_hash (VB, STRa(seq), false),
                           .qual  = segconf.deep_no_qual  ? 0 : deep_qual_hash (VB, STRa(qual), false) };

    if (flag.show_deep == 2) {
        if (deephash_issame (flag.debug_deep_hash, deep_hash)) 
            iprintf ("%s Found deep_hash=%u,%u,%u\nQNAME=\"%.*s\"\nSEQ=\"%.*s\"\nQUAL=\"%.*s\"\n",
                     LN_NAME, DEEPHASHf(deep_hash), STRf(qname), STRf(seq), STRf(qual));

        if (!flag.deep) return; // only debug
    }

    uint32_t hash = deep_hash.seq & bitmask32 (num_hash_bits);

    // find matching entry - and make sure there is only one matching entry. 
    // case: If there are multiple matching entries - we are not sure which SAM line this FASTQ read relates to, so we don't Deep.
    // case: All matching entries are already used by previous FASTQ lines - so we have multiple FASTQ lines claiming the same
    //       SAM entry - since we can't undo the previous lines at this point - we fail the execution and re-compress without deep
    ZipZDeep *matching_ent = NULL;
    for (uint32_t ent_i = hash_table[hash]; ent_i != NO_NEXT; ent_i = deep_ents[ent_i].next) {
        ZipZDeep *e = &deep_ents[ent_i];
        
        // case: SEQ matches, and QUAL and QNAME match if they are required to match
        bool qual_matches=false, qname_matches=false;
        
        if (e->hash.seq == deep_hash.seq &&  // SEQ matches
            (segconf.deep_no_qual || !e->hash.qual || (qual_matches  = (deep_hash.qual  == e->hash.qual ))) &&  // QUAL matches (or not relevant)
            (segconf.deep_no_qname                 || (qname_matches = (deep_hash.qname == e->hash.qname)))) {  // QNAME matches (or not relevant)
            
            // case: two or more SAM lines entries match this read according to the hashes - we don't know which is the true
            // one, so we seg without deep.
            if (matching_ent) {
                if (flag.show_deep) 
                    iprintf ("Two SAM lines with same hashes, skipping: this=[hash=%u,%u,%u vb_i=%u line_i=%u] other=[%u,%u,%u vb_i=%u line_i=%u] [qual_matches=%s qname_matches=%s deep_no_qual=%s deep_no_qname=%s]\n", 
                             DEEPHASHf(e->hash), e->vb_i, e->line_i, 
                             DEEPHASHf(matching_ent->hash), matching_ent->vb_i, matching_ent->line_i,
                             TF(qual_matches), TF(qname_matches), TF(segconf.deep_no_qual), TF(segconf.deep_no_qname));            

                vb->deep_stats[NDP_MULTI_MATCH]++;
                goto no_match;
            }

            matching_ent = e;
            *deep_seq  = true;
            *deep_qual = qual_matches;
            *deep_desc = qname_matches;
        }
    }

    // single match ent
    uint64_t deep_value;
    if (matching_ent) {
        deep_value = 1 + fastq_seg_deep_consume_unique_matching_ent (vb, matching_ent, deep_hash, STRa(qname), STRa(seq), STRa(qual)); // +1 as 0 means "no match"
        vb->deep_stats[NDP_DEEPABLE]++;
    }

    // case: no matching ent
    else {
        vb->deep_stats[NDP_NO_MATCH]++;

        if (flag.show_deep) {
            static int count = 0;
            if (count < 20) { //  approximate, not exact, due to multi threading
                rom once = "";
                DO_ONCE once = "\n\n(Showing first ~20 no-matches)\n";
                iprintf ("%s%s: no matching SAM line (has_match=%s, deep_hash=%u,%u,%u \nQNAME=\"%.*s\"\nSEQ=\"%.*s\"\nQUAL=\"%.*s\"\n", 
                         once, LN_NAME, TF(!!matching_ent), DEEPHASHf(deep_hash), STRf(qname), STRf(seq), STRf(qual));        
            
                __atomic_add_fetch (&count, (int)1,  __ATOMIC_RELAXED);
            }
        }

        no_match: 
        deep_value = 0;
        *deep_seq = *deep_qual = *deep_desc = false; // reset
    }

    if (flag.pair == NOT_PAIRED_END || flag.pair == PAIR_R1)
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
    BufferP deep_ents  = B(Buffer, z_file->deep_ents, vb_idx);
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

    uint8_t *deep_ent = CTX(FASTQ_DEEP)->last_value.p;
    PizZDeepFlags f = *(PizZDeepFlags *)deep_ent;

    deep_ent += segconf.deep_no_qname ? 1                // skip flags
                                      : 2 + deep_ent[1]; // skip flags, qname_len and qname 

    vb->seq_len = f.is_long_seq ? GET_UINT32 (deep_ent) : *deep_ent;
    deep_ent += f.is_long_seq ? 4 : 1;

    #define ENCODING(x) (f.seq_encoding == ZDEEP_SEQ_##x)

    // case: SEQ is compressed (bc it contains non-ACGT)
    if (ENCODING(COMPRESSED) || ENCODING(COMPRESSED_LONG_LEN)) {
        uint32_t comp_len = (ENCODING(COMPRESSED) ? *deep_ent : GET_UINT32 (deep_ent));
        deep_ent += (ENCODING(COMPRESSED) ? 1 : 4);

        unsigned out_len = vb->seq_len;        
        ASSPIZ (arith_uncompress_to (vb, deep_ent, comp_len, (uint8_t*)BAFTtxt, &out_len) && out_len == vb->seq_len,
                "Failed arith_uncompress_to SEQ copied from SAM: compressed_len=%u uncompressed_len=%u (expected: %u)", 
                comp_len, out_len, vb->seq_len);
    }

    // case: SEQ is packed
    else if (ENCODING(PACKED)) {
        // case: shortish reads - can fit on stack 
        #define SHORTISH_READ_SIZE 1024
        if ((ROUNDUP32(vb->seq_len) / 4 <= SHORTISH_READ_SIZE)) {
            uint8_t data[SHORTISH_READ_SIZE]; // not too long, so code path gets mileage
            memcpy (data, deep_ent, ROUNDUP4(vb->seq_len) / 4); // copy to word-align
            
            Bits pack = bits_init (vb->seq_len * 2, data, sizeof(data), false);
            bits_2bit_to_ACGT (BAFTtxt, &pack, 0, vb->seq_len);
        }

        // case: too long for stack - allocate on heap
        else {
            ASSERTNOTINUSE (vb->scratch);
            Bits pack = bits_alloc (vb->seq_len*2, true);
            memcpy (pack.words, deep_ent, ROUNDUP4(vb->seq_len) / 4);

            bits_2bit_to_ACGT (BAFTtxt, &pack, 0, vb->seq_len);
            bits_free (&pack);
        }
    }

    else { // ZDEEP_SEQ_PERFECT, ZDEEP_SEQ_PERFECT_REV, ZDEEP_SEQ_MIS_1, ZDEEP_SEQ_MIS_1_REV      
        PizDeepSeq ds = *(PizDeepSeq *)deep_ent; // not nessesarily word-aligned, but I think its ok, since this is just a memcpy, not loading integers
        deep_ent += (ENCODING(PERFECT) || ENCODING(PERFECT_REV)) ? 4 : 6;

        const Bits *genome=NULL;
        ref_get_genome (gref, &genome, NULL, NULL);
        
        BASE_ITER_INIT (genome, ds.gpos, vb->seq_len, (ENCODING(PERFECT) || ENCODING(MIS_1)));

        char *next = BAFTtxt;
        decl_acgt_decode;
        for (uint32_t i=0; i < vb->seq_len; i++) 
            *next++ = acgt_decode(BASE_NEXT);                    

        if (ENCODING(MIS_1)) 
            *(BAFTtxt + ds.mismatch_offset) = ds.mismatch_base; // reconstruct mismatch 

        else if (ENCODING(MIS_1_REV))
            *(BAFTtxt + vb->seq_len-1 - ds.mismatch_offset) = COMPLEM[(int)ds.mismatch_base]; // reconstruct mismatch 
    }

    Ltxt += vb->seq_len;

    // case we are pair-2: advance pair-1 SQBITMAP iterator, and if pair-1 is aligned - also its GPOS iterator
    if (VB_FASTQ->pair_vb_i/*we are R2*/ && fastq_piz_R1_test_aligned (VB_FASTQ))
        CTX(FASTQ_GPOS)->localR1.next++; // gpos_ctx->localR1.next is an iterator for both gpos and strand

    COPY_TIMER (fastq_special_deep_copy_SEQ);
    return NO_NEW_VALUE;    
}

SPECIAL_RECONSTRUCTOR (fastq_special_deep_copy_QUAL)
{
    START_TIMER;

    if (!reconstruct) return NO_NEW_VALUE;

    uint8_t *deep_ent = CTX(FASTQ_DEEP)->last_value.p;
    PizZDeepFlags f = *(PizZDeepFlags *)deep_ent;

    deep_ent += segconf.deep_no_qname ? 1                // skip flags
                                      : 2 + deep_ent[1]; // skip flags, qname_len and qname

    // skip seq_len
    deep_ent += f.is_long_seq ? 4 : 1;
    
    // skip seq
    switch (f.seq_encoding) {
        case ZDEEP_SEQ_COMPRESSED          : deep_ent += 1 + *deep_ent;             break;
        case ZDEEP_SEQ_COMPRESSED_LONG_LEN : deep_ent += 4 + GET_UINT32 (deep_ent); break;
        case ZDEEP_SEQ_PACKED              : deep_ent += ROUNDUP4(vb->seq_len) / 4; break;
        case ZDEEP_SEQ_PERFECT             : deep_ent += 4;                         break;
        case ZDEEP_SEQ_PERFECT_REV         : deep_ent += 4;                         break;
        case ZDEEP_SEQ_MIS_1               : deep_ent += 6;                         break;  
        case ZDEEP_SEQ_MIS_1_REV           : deep_ent += 6;                         break;  
        default: ABORT0 ("Unknown SEQ encoding");
    }

    uint32_t qual_comp_len = f.is_long_qual_comp ? GET_UINT32 (deep_ent) : *deep_ent;
    deep_ent += (f.is_long_qual_comp ? 4 : 1);

    unsigned out_len = vb->seq_len;        
    ASSPIZ (arith_uncompress_to (vb, deep_ent, qual_comp_len, (uint8_t *)BAFTtxt, &out_len) && out_len == vb->seq_len,
            "Failed arith_uncompress_to QUAL copied from SAM: compressed_len=%u uncompressed_len=%u (expected: %u)", 
            qual_comp_len, out_len, vb->seq_len);

    Ltxt += out_len;

    COPY_TIMER (fastq_special_deep_copy_QUAL);
    return NO_NEW_VALUE;
}
