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

//-----------------------------------------------------------
// Convert file-wide txt_deep_line_i to/from vb_i+deep_line_i
//-----------------------------------------------------------

// returns 0-based txt_deep_line_i - SAM file-wide 0-based line_i, but counting only lines that are not SUPP/SEC, not monochar, and SEQ.len>0.
// populated in sam_deep_zip_finalize.
static uint64_t vb_line_to_txt_line (VBIType vb_i, uint32_t line_i)
{
    ASSERT (vb_i <= z_file->vb_start_deep_line.len32, "Expecting vb_i=%u <= z_file->vb_start_deep_line.len=%u", vb_i, z_file->vb_start_deep_line.len32);

    return VB_start_deep_line(vb_i) + line_i;
}

static void txt_line_to_vb_line (uint64_t txt_deepable_line_i, 
                                 VBIType first_vb_i, VBIType last_vb_i,       // for binary search
                                 VBIType *vb_i, uint32_t *vb_deepable_line_i) // out
{
    if (first_vb_i > last_vb_i || VB_start_deep_line(first_vb_i) > txt_deepable_line_i) {
        *vb_i = first_vb_i - 1;
        *vb_deepable_line_i = txt_deepable_line_i - VB_start_deep_line (first_vb_i) - 1;
    }
        
    VBIType mid_vb_i = (first_vb_i + last_vb_i) / 2;

    uint64_t mid_vb_first_line = VB_start_deep_line(mid_vb_i);

    // case: txt_deepable_line_i is in a VB lower than mid_vb_i
    if (txt_deepable_line_i < mid_vb_first_line) 
        txt_line_to_vb_line (txt_deepable_line_i, first_vb_i, mid_vb_i-1, vb_i, vb_deepable_line_i);

    // case: txt_deepable_line_i is in mid_vb_i or higher VB. if it is in mid_vb_i, the stop condition will trigger
    // in the recursive call 
    else 
        txt_line_to_vb_line (txt_deepable_line_i, mid_vb_i+1, last_vb_i, vb_i, vb_deepable_line_i);
}

//----------------
// ZIP
//----------------

void fastq_deep_seg_initialize (VBlockP vb)
{
    // all-the-same for FASTQ_DEEP
    seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_set_deep }, 2, FASTQ_DEEP, 0);
    ctx_set_ltype (vb, LT_DYN_INT, FASTQ_DEEP, DID_EOL);
    ctx_set_store (vb, STORE_INT, DID_EOL); // actually, we store a pointer into one of the Buffers in z_file->deep_ents, but we treat it as an int
}

void fastq_deep_show_entries_stats (void) 
{
    ARRAY (uint32_t, hash_table, z_file->deep_hash);
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    printf ("z_file.deep_hash.len=%u z_file.deep_ents.len=%u\n", z_file->deep_hash.len32, z_file->deep_ents.len32);

    #define NUM_LENS 1000
    uint64_t used=0;
    uint32_t lens[NUM_LENS] = {}; // each entry i contains number linked lists which are of length i; last entry i - linked lists of i or longer
    for (uint64_t hash=0; hash < hash_table_len; hash++) {

        uint32_t this_len=0;
        for (uint32_t ent_i = hash_table[hash]; ent_i != NO_NEXT; ent_i = deep_ents[ent_i].next) 
            this_len++;
        
        if (this_len) used++;

        lens[MIN_(this_len,NUM_LENS-1)]++; // one more chain of length this_len
    }

    printf ("Deep hash entries used: %"PRIu64"/%"PRIu64" (%.1f%%)\n", used, hash_table_len, 100.0 * (double)used / (double)hash_table_len);

    printf ("Deep hash linked list length histogram:\n");
    uint64_t total_ents_traversed = 0; // total ents traversed on linked lists when segging FASTQ reads
    for (int this_len=0; this_len < NUM_LENS; this_len++) {
        if (!lens[this_len]) continue;

        printf ("%3u: %u\n", this_len, lens[this_len]);

        // our algorithm requires traversal of entire linked list for each entry. calculate the length of the list
        // traversed for the average read
        total_ents_traversed += this_len        // number of entries that will access this linked list
                              * this_len        // length of linked list traversed for each entry accessing this list
                              * lens[this_len]; // number of linked lists of this length
    }

    uint64_t count_not_consumed = 0;
    for (uint64_t ent_i=0; ent_i < deep_ents_len; ent_i++)
        if (!(deep_ents[ent_i].line_i >> 31)) {
            count_not_consumed++;

        if (count_not_consumed <= 50)
           printf ("Not consumed (first 50): hash=%u,%u,%u sam_vb=%u line_i=%u\n", DEEPHASHf(deep_ents[ent_i].hash), deep_ents[ent_i].vb_i, deep_ents[ent_i].line_i);
    }

    // this is approximate: assuming # of ents hashed is equal to the number of reads
    printf ("FASTQ SEQ: Linked list length traversed on average: %.2f\n", (double)total_ents_traversed / (double)deep_ents_len);
    
    printf ("Number of SAM ents not consumed by FASTQ reads: %"PRIu64"\n", count_not_consumed);
}

void fastq_deep_seg_finalize_segconf (uint32_t n_lines)
{
    if (flag.debug_deep) 
        printf ("segconf: n_lines=%u n_full_mch=%u n_seq_qname_mch=%u n_seq_qual_mch=%u n_seq_only=%u n_no_match=%u\n",
                n_lines, segconf.n_full_mch, segconf.n_seq_qname_mch, segconf.n_seq_qual_mch, segconf.n_seq_mch,segconf.n_no_mch);

    unsigned n_mch = n_lines - segconf.n_no_mch; // exclude lines that don't match any SAM line - perhaps they were filtered out and excluded from the SAM file 
    unsigned threashold = MAX_(n_mch/2, 1);      // at least half of the lines need this level of matching

    // test: at least have of the reads (excluding reads that had no matching SEQ in the SAM) have a matching QNAME
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

void fastq_seg_deep (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(qname), STRp(seq), STRp(qual), 
                     bool *deep_desc, bool *deep_seq, bool *deep_qual) // out - set to true, or left unchanged
{
    ARRAY (uint32_t, hash_table, z_file->deep_hash);
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    if (flag.deep) { // skip tests if only --debug-deep without --deep
        if (!deep_ents_len) return;

        if (segconf.running) {
            fastq_deep_seg_segconf (vb, STRa(qname), STRa(seq), STRa(qual));
            goto no_match;
        }

        if (dl->monobase) goto no_match; // we don't deep monobase reads - too much contention
    }

    DeepHash deep_hash = { .qname = segconf.deep_no_qname ? 0 : deep_qname_hash (STRa(qname)),
                           .seq   =                             deep_seq_hash (VB, STRa(seq), false),
                           .qual  = segconf.deep_no_qual  ? 0 : deep_qual_hash (VB, STRa(qual), false) };

    if (flag.debug_deep == 2) {
        if (deephash_issame (flag.debug_deep_hash, deep_hash)) 
            printf ("%s Found deep_hash=%u,%u,%u\nQNAME=\"%.*s\"\nSEQ=\"%.*s\"\nQUAL=\"%.*s\"\n",
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
                if (flag.debug_deep) 
                    printf ("Two SAM lines with same hashes, skipping: this=[hash=%u,%u,%u vb_i=%u line_i=%u] other=[%u,%u,%u vb_i=%u line_i=%u] [qual_matches=%s qname_matches=%s deep_no_qual=%s deep_no_qname=%s]\n", 
                             DEEPHASHf(e->hash), e->vb_i, e->line_i, 
                             DEEPHASHf(matching_ent->hash), matching_ent->vb_i, matching_ent->line_i,
                             TF(qual_matches), TF(qname_matches), TF(segconf.deep_no_qual), TF(segconf.deep_no_qname));            
                goto no_match;
            }

            matching_ent = e;
            *deep_seq  = true;
            *deep_qual = qual_matches;
            *deep_desc = qname_matches;
        }
    }

    // single match ent
    if (matching_ent) {

        // atomicly set "consumed" (MSb of line_i) to true, and verify that no other thread beat us to it.
        uint32_t expected = matching_ent->line_i & 0x7fffffffUL; // expected - not consumed
        uint32_t set_to   = vb->line_i | 0x80000000UL; // set to - consumed + store consuming FASTQ line
        // xxxuint32_t set_to   = matching_ent->line_i | 0x80000000UL; // set to - consumed
        bool is_first_consumer = __atomic_compare_exchange_n (&matching_ent->line_i, &expected, set_to, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED);

        // case: this ent was already consumed by a previous FASTQ read - we don't know which read this SAM line belongs to, and we can't reverse the previous segged line
        // TO DO: recompress normally by executing another instance of genozip
        ASSINP (is_first_consumer, "Deep %s: Another FASTQ read (vb=%u line_i=%u) maps to the same SAM line - cannot --deep these files - please compress normally:\n"
                                   "deep_hash=%u,%u,%u this_line:\nQNAME=\"%.*s\"\nSEQ=\"%.*s\"\nQUAL=\"%.*s\"",
                LN_NAME, matching_ent->vb_i, matching_ent->line_i, DEEPHASHf(deep_hash), STRf(qname), STRf(seq), STRf(qual));

        int64_t txt_deep_line_i = vb_line_to_txt_line (matching_ent->vb_i, (matching_ent->line_i & 0x7fffffffUL));
        seg_add_to_local_resizable (VB, CTX(FASTQ_DEEP), 1 + txt_deep_line_i, 0); // +1 as 0 means "no match"
// printf ("xxx FQ_vb=%u FQ_line_i=%u SAM_vb=%u SAM_line_i=%u txt_deep_line_i=%u deep_seq=%u deep_qual=%u deep_desc=%u\n", 
// vb->vblock_i, vb->line_i, matching_ent->vb_i, matching_ent->line_i & 0x7fffffffUL, txt_deep_line_i, *deep_seq, *deep_qual, *deep_desc);
        matching_ent->vb_i = vb->vblock_i; // after consumption - store consuming FASTQ vb
    }

    // case: no matching ent
    else {
        if (flag.debug_deep) 
            printf ("%s: no uniquely matching SAM line (has_match=%s, deep_hash=%u,%u,%u) \nQNAME=\"%.*s\"\nSEQ=\"%.*s\"\nQUAL=\"%.*s\"\n", 
                    LN_NAME, TF(!!matching_ent), DEEPHASHf(deep_hash), STRf(qname), STRf(seq), STRf(qual));        
        
        no_match: 
        seg_add_to_local_resizable (VB, CTX(FASTQ_DEEP), 0, 0);
        *deep_seq = *deep_qual = *deep_desc = false; // reset
    }
}

//----------------
// PIZ
//----------------

SPECIAL_RECONSTRUCTOR (fastq_special_set_deep)
{
    uint64_t txt_deepable_line_i = reconstruct_from_local_int (vb, ctx, 0, RECON_OFF);

    VBIType vb_i;
    uint32_t vb_deepable_line_i;
    txt_line_to_vb_line (txt_deepable_line_i, 0, z_file->vb_start_deep_line.len32-1, &vb_i, &vb_deepable_line_i);

    // get this VB's data
    BufferP deep_ents  = B(Buffer, z_file->deep_ents, vb_i-1);
    BufferP deep_index = B(Buffer, z_file->deep_index, vb_i-1);

    ASSPIZ (vb_deepable_line_i < deep_index->len32, "Unexpected vb_deepable_line_i=%u < deep_index->len=%u", vb_deepable_line_i, deep_index->len32);
    uint32_t ents_index = *B32 (*deep_index, vb_deepable_line_i);

    new_value->p = B8 (*deep_ents, ents_index);

    return HAS_NEW_VALUE; 
}

SPECIAL_RECONSTRUCTOR (fastq_special_copy_deep)
{
}
