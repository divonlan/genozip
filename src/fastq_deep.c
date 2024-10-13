// ------------------------------------------------------------------
//   fastq_deep.c
//   Copyright (C) 2023-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "fastq_private.h"
#include "deep.h"
#include "huffman.h"
#include "codec.h"
#include "zip_dyn_int.h"
#include "reconstruct.h"
#include "htscodecs/arith_dynamic.h"

sSTRl(con_decanonize1_snip,96);
sSTRl(con_decanonize2_snip,96);

#define DF vb->piz_deep_flags

//----------------
// ZIP
//----------------

void fastq_deep_zip_initialize (void)
{
    DO_ONCE {
        SmallContainer con1 = { .repeats = 1, .nitems_lo = 2, .items = { { .dict_id.num = _FASTQ_Q0NAME  }, { .dict_id.num = _FASTQ_QmNAME  } }};
        container_prepare_snip ((ContainerP)&con1, NULL, 0, qSTRa(con_decanonize1_snip));

        SmallContainer con2 = { .repeats = 1, .nitems_lo = 2, .items = { { .dict_id.num = _FASTQ_Q0NAME2 }, { .dict_id.num = _FASTQ_QmNAME2 } }};
        container_prepare_snip ((ContainerP)&con2, NULL, 0, qSTRa(con_decanonize2_snip));
    }
}

void fastq_deep_zip_after_compute (VBlockFASTQP vb)
{
    for (int i=0; i < NUM_DEEP_STATS_ZIP; i++)
        z_file->deep_stats[i] += vb->deep_stats[i];
}

void fastq_deep_seg_initialize (VBlockFASTQP vb)
{
    ctx_set_dyn_int (VB, FASTQ_DEEP, DID_EOL); // this also sets STORE_INT. actually, we store a pointer into one of the Buffers in z_file->deep_ents, but we treat it as an int

    if (IS_R1 || flag.pair == NOT_PAIRED) 
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_set_deep, '0' }, 3, FASTQ_DEEP, 0);  // all-the-same for FASTQ_DEEP

    else { // pair-2
        ctx_set_ltype (VB, LT_INT8, FASTQ_DEEP_DELTA, DID_EOL); 
        ctx_consolidate_stats (VB, FASTQ_DEEP, FASTQ_DEEP_DELTA, DID_EOL);
    }

    CTX(FASTQ_SQBITMAP)->no_stons = true;
}

// called by main thread after ALL FASTQ files are done compressing with Deep
void fastq_deep_zip_finalize (void) 
{
    if (flag.show_deep) {
        uint64_t total = z_file->deep_stats[NDP_FQ_READS];
        iprint0 ("\nFASTQ reads breakdown by deepability:\n");
        
        for (int i=0; i < NUM_DEEP_STATS_ZIP; i++) 
            if (z_file->deep_stats[i])
                iprintf ("%-11.11s: %"PRIu64" (%.1f%%)\n", (rom[])DEEP_STATS_NAMES_ZIP[i], z_file->deep_stats[i], 100.0 * (double)z_file->deep_stats[i] / (double)total);
    }

    ASSERTISALLOCED (z_file->deep_ents);
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    // Count deepable alignments in SAM that were did not appear in FASTQ. This would be an user
    // error as we require that FASTQ files covering all SAM alignments (except supplementary, 
    // secondary and consensus) must be provided when using --deep 
    uint64_t count_unconsumed=0, count_dups=0, first_unconsumed=0;
    
    for (uint64_t ent_i=0; ent_i < deep_ents_len; ent_i++) 
        if (deep_ents[ent_i].dup) count_dups++;
        else if (!deep_ents[ent_i].consumed) {
            if (!count_unconsumed) first_unconsumed = ent_i;
            count_unconsumed++;
        }

    if (flag.show_deep) {
        #define NUM_SHOW 20
        iprintf ("\nNumber of %s deepable alignments not consumed by FASTQ reads: %"PRIu64"\n"
                "- Unusable due to multiple alignments with same hash: %"PRIu64"\n"
                "- No matching FASTQ read: %"PRIu64" (showing first %d)\n",
                z_dt_name(), count_unconsumed + count_dups, count_dups, count_unconsumed, (int)MIN_(count_unconsumed, NUM_SHOW));

        // show the first few unconsumed
        if (count_unconsumed)
            for (uint64_t ent_i=0, count=0; ent_i < deep_ents_len && count < NUM_SHOW; ent_i++)
                if (!deep_ents[ent_i].consumed && !deep_ents[ent_i].dup) {
                    iprintf ("sam_vb=%u line_i=%u hash=%016"PRIx64",%08x,%08x\n", deep_ents[ent_i].vb_i, deep_ents[ent_i].line_i, DEEPHASHf(deep_ents[ent_i].hash));
                    count++;
                }
    }

    // technical explanation: if a BAM alignment that has no correspoding FASTQ read it is an indication 
    // that not all BAM alignments are covered by FASTQ reads (contrary to Genozip user instructions).
    // The real problem however, would not be detected here - it is a case that a BAM alignment that is not
    // covered by a FASTQ read, has, by chance, the same hashes as a FASTQ read that is not present in the 
    // BAM data (which is allowed). If this occurs, Deep will incorrectly map the FASTQ read to the redundant 
    // BAM alignment. This BAM alignment will be marked as "consumed" and hence not counted in 
    // "count_unconsumed". However, --test would catch this and error on reconstruction. 
    else if (count_unconsumed) {
        WARN ("WARNING: detected %"PRIu64" %s alignments (other than supplementary, secondary and consensus alignments) which "
              "are absent in the FASTQ file(s). Genozip requires that the FASTQ files included in --deep cover all alignments "
              "in the %s file, or otherwise, in rare cases, the resulting compressed file might be corrupted. If such corruption occurs, "
              "the testing that will follow now will detect it. If the testing completes successfully, then there is no problem. "
              "An example of an alignment present in the %s file but missing in the FASTQ file(s) can be obtained by running:\n"
              "   %s\n",
              count_unconsumed, z_dt_name(), z_dt_name(), z_dt_name(),
              piz_advise_biopsy_line (COMP_NONE, deep_ents[first_unconsumed].vb_i, deep_ents[first_unconsumed].line_i, segconf.sam_deep_filename).s);

        flag.test = true; // force testing in this case
    }
}

static void fastq_deep_show_segconf (uint32_t n_lines, uint32_t proof_of_first, uint32_t proof_of_last, uint32_t proof_of_qname1, uint32_t proof_of_qname2, uint32_t proof_of_qual)
{
    iprintf ("\nFASTQ Deep segconf stats: n_lines=%u proof_of_first=%u proof_of_last=%u proof_of_qname1=%u proof_of_qname2=%u proof_of_qual=%u\n"
             "FIRST:{n_full_mch=(%u,%u) n_seq_qname_mch=(%u,%u)} LAST:{n_full_mch=(%u,%u) n_seq_qname_mch=(%u,%u)} NOT_PAIR={n_full_mch=(%u,%u) n_seq_qname_mch=(%u,%u)} n_no_mch=%u trimming=%s n_full_mch_trimmed=%u\n",
             n_lines, proof_of_first, proof_of_last, proof_of_qname1, proof_of_qname2, proof_of_qual,
             segconf.n_full_mch[2], segconf.n_full_mch[3], segconf.n_seq_qname_mch[2], segconf.n_seq_qname_mch[3], 
             segconf.n_full_mch[4], segconf.n_full_mch[5], segconf.n_seq_qname_mch[4], segconf.n_seq_qname_mch[5], 
             segconf.n_full_mch[0], segconf.n_full_mch[1], segconf.n_seq_qname_mch[0], segconf.n_seq_qname_mch[1],                  
             segconf.n_no_mch, segconf_deep_trimming_name(), segconf.n_full_mch_trimmed);
    
    iprintf ("\nFASTQ Deep segconf settings: paired_qname=%s is_last=%s qtype=%s has_trimmed=%s has_trimmed_left=%s no_qual=%s\n",
             TF(segconf.deep_paired_qname), YN(segconf.deep_is_last), qtype_name (segconf.deep_qtype), TF(segconf.deep_has_trimmed), TF(segconf.deep_has_trimmed_left), TF(segconf.deep_no_qual));
}

void fastq_deep_seg_finalize_segconf (uint32_t n_lines)
{
    if (zip_is_biopsy) return;

    uint32_t proof_of_first=0, proof_of_last=0, proof_of_qname1=0, proof_of_qname2=0, proof_of_qual=0;

    if (segconf.deep_paired_qname) {
        proof_of_first = segconf.n_seq_qname_mch[2] + segconf.n_seq_qname_mch[3] + segconf.n_full_mch[2] + segconf.n_full_mch[3];  
        proof_of_last  = segconf.n_seq_qname_mch[4] + segconf.n_seq_qname_mch[5] + segconf.n_full_mch[4] + segconf.n_full_mch[5];  
        
        if      (proof_of_first && proof_of_first >= 9 * proof_of_last)  segconf.deep_is_last = no;
        
        else if (proof_of_last  && proof_of_last  >= 9 * proof_of_first) segconf.deep_is_last = yes;

        // usually fastq_verify_and_sort_pairs would cause the first (i.e. segconfed) FASTQ file to be
        // the true R1, and usually the user creates the SAM/BAM file with an aligner command line with R1 listed first.
        // If we are wrong, then no FASTQ read will be segged against Deep, which is not great, but not an error.  
        else segconf.deep_is_last = no; 
    }

    unsigned qname_i = segconf.deep_is_last==no?2 : segconf.deep_is_last==yes? 4 : 0/*not paired*/;
    unsigned qname2_i = qname_i + 1;

    // test: at least half of the reads (excluding reads that had no matching SEQ in the SAM) have a matching QNAME
    proof_of_qname1 = segconf.n_seq_qname_mch[qname_i]  + segconf.n_full_mch[qname_i];
    proof_of_qname2 = segconf.n_seq_qname_mch[qname2_i] + segconf.n_full_mch[qname2_i];

    if (proof_of_qname1 || proof_of_qname2) {
        segconf.deep_qtype = (proof_of_qname1 >= proof_of_qname2) ? QNAME1 : QNAME2;

        segconf.deep_has_trimmed = (segconf.n_full_mch_trimmed > 0); // true is some reads might be trimmed (shorter in SAM than FASTQ)
        
        // test: likewise, matching QUAL
        proof_of_qual = segconf.nontrivial_qual &&  // don't deep trivial qual: to avoid bloating deep_ents in piz with qual data without any compression gain 
                        !flag.deep_no_qual      &&  // don't deep qual if user specified --deep=no-qual
                        (segconf.deep_qtype == QNAME1 ? segconf.n_full_mch[qname_i] : segconf.n_full_mch[qname2_i]);

        if (!proof_of_qual) {
            segconf.deep_no_qual = true;
            segconf.deep_N_fq_score = segconf.deep_N_sam_score = 0;
        }
    }

    else { 
        // no luck finding matching lines in segconf - but continue if user insists...
        if (flag.force_deep) {
            segconf.deep_qtype       = QNAME1; // assume match is of QNAME1
            segconf.deep_no_qual     = false;  // assume QUAL matches too
            segconf.deep_has_trimmed = true;   // be liberal, allow trimming (since we require QNAME anyway)
        }
        
        else {
            fastq_deep_show_segconf (n_lines, proof_of_first, proof_of_last, 0, 0, 0);

            STR(master_qname);
            huffman_get_master (SAM_QNAME, pSTRa(master_qname));            
        
            ABORTINP ("Error: cannot use --deep with this file: based on testing the first %u reads of %s.\n"
                      "You may compress these files without --deep, or override with --force-deep\n"
                      "First SAM: %.*s First FASTQ: %s%s",
                      n_lines, txt_name, STRf(master_qname), segconf.deep_1st_desc, report_support_if_unexpected());
        }
    }        

    if (flag.show_deep) 
        fastq_deep_show_segconf (n_lines, proof_of_first, proof_of_last, proof_of_qname1, proof_of_qname2, proof_of_qual);
}

// find a subset of FASTQ's seq that matches the hash of the SAM's seq. 
static bool fastq_deep_seg_find_subseq (VBlockFASTQP vb, STRp (fastq_seq), 
                                        uint32_t sam_seq_len, uint32_t sam_seq_hash, 
                                        bool allow_offset,        // search all possible offsets, not just 0
                                        uint32_t *sam_seq_offset) // out if successful
{
    START_TIMER;

    uint32_t max_offset = allow_offset ? (fastq_seq_len - sam_seq_len) : 0;
    uint32_t found = 0;

    for (uint32_t i=0; i <= max_offset; i++) {

        // hash of the subsequence of FASTQ at offset i and length sam_seq_len
        uint32_t subseq_hash = deep_seq_hash (VB, fastq_seq + i, sam_seq_len, false); 
        
        if (subseq_hash == sam_seq_hash) { 
            // case: first match of a FASTA sub-seq to the SAM seq - this will be our result 
            if (!found) { 
                *sam_seq_offset = i;
                found = true;
            }

            // case: a previous *different* subseq was found, that by chance has the same hash:
            // we don't know which is the correct subseq that matches the SAM seq, so we fail 
            else if (memcmp (fastq_seq + i, fastq_seq + *sam_seq_offset, sam_seq_len)) {
                COPY_TIMER(fastq_deep_seg_find_subseq);
                return false;
            }
            
            // case: this subseq is identical to a previous subseq found - no issues
            else {}
        }
    }

    COPY_TIMER(fastq_deep_seg_find_subseq);

    return found;
}

// if in SAM, all qual score of 'N's is set to a particular value, we adjust qual here, so that the hash may match
static rom fastq_deep_set_N_qual (VBlockFASTQP vb, STRp(seq), STRp(qual))
{
    char *n;
    if (!(n = memchr (seq, 'N', seq_len))) return qual; // quick return in common case: no 'N's, so no change in qual

    // set deep_N_fq_score if not set already. note: if set, this is definitely the value. if 0, another thread might be attempting to set in concurrently.
    if (!segconf.deep_N_fq_score) { 
        char expected = 0;
        __atomic_compare_exchange_n (&segconf.deep_N_fq_score, &expected, qual[n-seq], false, __ATOMIC_ACQ_REL, __ATOMIC_RELAXED);
    }

    // at this point we can trust that segconf.deep_N_fq_score is non zero and immutable

    char deep_N_sam_score = segconf.deep_N_sam_score; // make a copy in case another thread cancels in the middle of the if ↓. so we can avoid atomic.

    // in case deep_N_fq_score is already the requested value - cancel for future lines and threads.
    if (!deep_N_sam_score || // another thread canceled - nothing for us to do
        deep_N_sam_score == segconf.deep_N_fq_score) {
        segconf.deep_N_sam_score = 0; // release whenever, doesn't need to be atomic
        return qual;
    }

    ASSERTNOTINUSE (vb->scratch);
    ARRAY_alloc (char, new_qual, seq_len, false, vb->scratch, vb, "scratch");

    for (uint32_t i=0; i < seq_len; i++)
        if (seq[i] == 'N') {
            if (qual[i] == segconf.deep_N_fq_score)
                new_qual[i] = deep_N_sam_score;
            else 
                return NULL; // unexpected inconsistent quality score of an 'N' base - don't Deep this read since in reconstruction quality of 'N' will be set to deep_N_fq_score          
        }
        else
            new_qual[i] = qual[i];
    
    return new_qual;
}

static void fastq_deep_seg_segconf (VBlockFASTQP vb, STRp(qname), STRp(qname2), STRp(seq), STRp(qual))
{
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    struct { bool cond; QType q; STR(qname); thool is_last; } inst[NUM_INSTS] = 
          // QNAME1                                                       QNAME2
        { { !segconf.deep_paired_qname, QNAME1, STRa(qname), unknown }, { !segconf.deep_paired_qname && qname2_len > 0, QNAME2, STRa(qname2), unknown },   // not paired
          {  segconf.deep_paired_qname, QNAME1, STRa(qname), false   }, { segconf.deep_paired_qname  && qname2_len > 0, QNAME2, STRa(qname2), false   },   // is_first
          {  segconf.deep_paired_qname, QNAME1, STRa(qname), true    }, { segconf.deep_paired_qname  && qname2_len > 0, QNAME2, STRa(qname2), true    } }; // is_last

    uint64_t qname_hash[NUM_INSTS] = {};
    for (int i=0; i < NUM_INSTS; i++)
        if (inst[i].cond) 
            qname_hash[i] = deep_qname_hash (inst[i].q, STRa(inst[i].qname), inst[i].is_last, NULL);
        
    if (segconf.deep_N_sam_score) {
        qual = fastq_deep_set_N_qual (vb, STRa(seq), STRa(qual));
        if (!qual) {
            segconf.n_no_mch++; // some qual scores of 'N' bases are of an inconsistent value - we don't be able to reconstruct them from segconf.deep_N_fq_score
            goto done;
        }
    }

    bool seq_qname_matches[6]={};

    for (int i=0; i < NUM_INSTS; i++) {
        if (!inst[i].cond) continue;

        uint32_t hash = qname_hash[i] & bitmask64 (num_hash_bits);
        
        for (uint32_t ent_i = *B32(z_file->deep_index, hash); ent_i != NO_NEXT; ent_i = deep_ents[ent_i].next) {
            ZipZDeep *e = &deep_ents[ent_i];
            
            if (e->hash.qname != qname_hash[i] || // possible, because hash uses less bits that e->hash.qname
                e->seq_len > seq_len) continue;   // definitely not a match if FASTQ seq_len is shorter than SAM's 
        
            uint32_t sam_seq_offset = 0;
            if (!fastq_deep_seg_find_subseq (vb, STRa(seq), e->seq_len, e->hash.seq, true, &sam_seq_offset))
                continue; 

            if (sam_seq_offset > 0)
                segconf.deep_has_trimmed_left = true; // segconf observed trimming left bases 
            
            uint32_t qual_hash = deep_qual_hash (VB, qual + sam_seq_offset, e->seq_len, false); // possibly trimmed QUAL
            if (qual_hash == e->hash.qual) {
                segconf.n_full_mch[i]++;
                if (seq_len != e->seq_len) segconf.n_full_mch_trimmed++;
                goto done; // QNAME, SEQ and QUAL match (possibly trimmed)
            }
            else
                seq_qname_matches[i] = true; // only QNAME and SEQ match - continue searching, perhaps a full match will be found
        }
    }
    
    for (int i=0; i < NUM_INSTS; i++)
        if (seq_qname_matches[i]) {
            segconf.n_seq_qname_mch[i]++; 
            goto done;
        } 

    segconf.n_no_mch++; // no match - likely bc this read was filtered out in the SAM file

done:
    if (segconf.deep_N_sam_score) 
        buf_free (vb->scratch);
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
    // case: this ent was already consumed by a previous FASTQ read. This can happen if two different FASTQ lines 
    // have the same (QNAME,SEQ,QUAL) or (QNAME,SEQ,0) hashes AND only one of the two corresponding SAM lines 
    // actually exists in the SAM file (if both existed in the SAM they would have been marked as Dup and not deeped). 
    // Therefore, both reads map to the single SAM line, and we don't know which is the true match. We have no choice but to abort.
    ABORTINP ("Deep: We hit a rare edge case: two FASTQ reads: current line %s and line %s/%u/%u (vb_size=%s)\n"
              "map to the same SAM line (deep_hash=%016"PRIx64",%08x,%08x qtype=%s no_qual=%s).\n"
              "%sWork around: compress without using --deep.\n"
              "This_line:\nQNAME=\"%.*s\"\nSEQ  =\"%.*s\"\nQUAL =\"%.*s\"\n"
              "To see other line, run: %s",
              LN_NAME, comp_name (sections_vb_header (contention_place.vb_i)->comp_i), contention_place.vb_i, contention_place.line_i, str_size (segconf.vb_size).s, 
              DEEPHASHf(deep_hash), qtype_name(segconf.deep_qtype), TF(segconf.deep_no_qual), 
              report_support(), STRf(qname), STRf(seq), STRf(qual),
              piz_advise_biopsy_line (COMP_NONE, contention_place.vb_i, contention_place.line_i, NULL).s);
    
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
static bool fastq_seg_deep_is_dup (VBlockFASTQP vb,  
                                   ZipZDeep *e) // first entry on link list that matches criteria
{
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);
    ZipZDeep *e2 = NULL;

    if (!e->dup) // not already marked as dup by another FASTQ read that matched the same criteria
        for (uint32_t ent_i = e->next; ent_i != NO_NEXT; ent_i = deep_ents[ent_i].next) {
            e2 = &deep_ents[ent_i];
            if ( e->hash.qname == e2->hash.qname &&
                 e->hash.seq   == e2->hash.seq   && 
                (e->hash.qual  == e2->hash.qual || segconf.deep_no_qual))
                // first encounter with this dup - traverse the rest of the linked list looking for dups
                // note: no worries about atomicity of visibility of "dup" because all matching threads
                // will discover and set this dup for all dup entries, and thereafter the entries will be immutable 
                // (because they will not be consumed and hence e->consumed will not be set)
                e->dup = e2->dup = true; // possibly setting e->dup multiple times, that's ok
        }        

    if (e->dup && flag.show_deep && (flag.show_deep != SHOW_DEEP_ONE_HASH || deephash_issame (flag.debug_deep_hash, e->hash))) {
        iprintf ("%s: Two or more %s alignments with same hashes (%s,no_qual=%s), skipping: this=[MAIN/%u/%u hash=%016"PRIx64",%08x,%08x]", 
                 LN_NAME, z_dt_name(), qtype_name (segconf.deep_qtype), TF(segconf.deep_no_qual), 
                 e ->vb_i, e ->line_i, DEEPHASHf(e ->hash));
        if (e2) iprintf (" another=[MAIN/%u/%u %016"PRIx64",%08x,%08x]", e2->vb_i, e2->line_i, DEEPHASHf(e2->hash));
        iprint0("\n");
    }

    return e->dup;
}

static DeepStatsZip fastq_seg_find_deep (VBlockFASTQP vb, ZipDataLineFASTQ *dl, DeepHash *deep_hash, STRp(seq), STRp(qual), 
                                           ZipZDeep **matching_ent) // out
{
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    for (uint32_t ent_i = *B32 (z_file->deep_index, deep_hash->qname & bitmask64 (num_hash_bits)); 
         ent_i != NO_NEXT; 
         ent_i = deep_ents[ent_i].next) {

        ZipZDeep *e = &deep_ents[ent_i];

        if (e->hash.qname != deep_hash->qname || // possible, because hash uses less bits that e->hash.qname
            e->seq_len > seq_len             || // definitely not a match if FASTQ seq_len is shorter than SAM's   
            (e->seq_len < seq_len && !segconf.deep_has_trimmed)) // not a match if FASTQ seq_len is too long, and we don't allow trimming
             continue;     

        // case: QNAME matches: see if we can find a subsequence of SEQ that matches, and if needed, also a subsequence of QUAL
        if (!fastq_deep_seg_find_subseq (vb, STRa(seq), e->seq_len, e->hash.seq, segconf.deep_has_trimmed_left, &dl->sam_seq_offset) ||
            (!segconf.deep_no_qual && deep_qual_hash (VB, qual + dl->sam_seq_offset, e->seq_len, false) != e->hash.qual))
            continue; 

        // update to subsequence found (for messages)
        deep_hash->seq  = e->hash.seq;
        deep_hash->qual = e->hash.qual;
        
        // case: two or more SAM lines entries match this read according to the hashes - we don't know which is the true one, so we seg without deep.
        if (fastq_seg_deep_is_dup (vb, e))
            return (e->seq_len == seq_len) ? NDP_SAM_DUP : NDP_SAM_DUP_TRIM;

        *matching_ent = e;
        return (e->seq_len == seq_len) ? NDP_DEEPABLE : NDP_DEEPABLE_TRIM;
    }

    return NDP_NO_MATCH;
}

static void fastq_seg_deep_do (VBlockFASTQP vb, uint64_t deep_value)
{
    decl_ctx (FASTQ_DEEP);

    if (flag.pair == NOT_PAIRED || IS_R1)
        dyn_int_append (VB, ctx, deep_value, 0);

    else { // PAIR_R2
        uint64_t pair_1_deep_value = fastq_get_pair_deep_value (vb, ctx); // consume whether or not used

        // delta carefully to avoid overflow if values are large uint64's
        uint64_t abs_delta = (pair_1_deep_value > deep_value) ? (pair_1_deep_value - deep_value)
                                                              : (deep_value - pair_1_deep_value);                                       
        if (abs_delta > 127) { 
            dyn_int_append (VB, ctx, deep_value, 0);
            seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_set_deep, '0' }, 3, FASTQ_DEEP, 0);
        }   

        else {
            int8_t delta8 = (pair_1_deep_value > deep_value) ? -(int8_t)abs_delta : (int8_t)abs_delta;
            seg_integer_fixed (VB, CTX(FASTQ_DEEP_DELTA), &delta8, false, 0);
            seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_set_deep, '1' }, 3, FASTQ_DEEP, 0);
        }
    }
}

bool fastq_seg_deep (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(qname), STRp(qname2), STRp(seq), STRp(qual), 
                     uint32_t *uncanonical_suffix_len) // out - suffix of deeped QNAME beyond canonical, i.e. not copied from Deep
{
    START_TIMER;

    uint64_t deep_value = 0;
    DeepStatsZip reason;
    #define NO_MATCH(r) ({ reason=(r); goto no_match; })

    if (flag.deep) { // skip tests if only --show-deep without --deep
        vb->deep_stats[NDP_FQ_READS]++;
        
        // SAM file did not produce any Deep entries (e.g. because all alignments seconday/supplementary)
        if (!z_file->deep_ents.len) 
            NO_MATCH (NDP_NO_ENTS);

        if (segconf.running) {
            fastq_deep_seg_segconf (vb, STRa(qname), STRa(qname2), STRa(seq), STRa(qual));
            goto do_seg;
        }

        // we don't Deep reads monochar SEQ - too much contention
        if (dl->monochar) 
            NO_MATCH (NDP_MONOSEQ);
    }

    if (segconf.deep_N_sam_score && !segconf.deep_no_qual) {
        qual = fastq_deep_set_N_qual (vb, STRa(seq), STRa(qual));
        if (!qual) 
            NO_MATCH(NDP_BAD_N_QUAL); // some qual scores of 'N' bases are of an inconsistent value - we don't be able to reconstruct them from segconf.deep_N_fq_score
    }

    bool calc_hash_for_show = (flag.show_deep && !segconf.deep_has_trimmed);
    
    DeepHash deep_hash = { .qname = (segconf.deep_qtype == QNAME1) ? deep_qname_hash (QNAME1, STRa(qname),  segconf.deep_is_last, uncanonical_suffix_len)
                                  :                                  deep_qname_hash (QNAME2, STRa(qname2), segconf.deep_is_last, uncanonical_suffix_len),
                           .seq   = calc_hash_for_show ? deep_seq_hash (VB, STRa(seq), false) : 0,
                           .qual  = (calc_hash_for_show && !segconf.deep_no_qual) ? deep_qual_hash (VB, STRa(qual), false) : 0 };

    if (flag.show_deep == SHOW_DEEP_ONE_HASH || flag.show_deep == SHOW_DEEP_ALL) {
        if (deephash_issame (flag.debug_deep_hash, deep_hash) || flag.show_deep == SHOW_DEEP_ALL) 
            iprintf ("%s Searching for deep_hash=%016"PRIx64",%08x,%08x\nQNAME=\"%.*s\"\nSEQ=\"%.*s\"\nQUAL=\"%.*s\"\n\n",
                     LN_NAME, DEEPHASHf(deep_hash), STRf(qname), STRf(seq), STRf(qual));

        if (!flag.deep) goto done; // only debug
    }

    // find matching entry - and make sure there is only one matching entry. 
    // case: If there are multiple matching entries - we are not sure which SAM line this FASTQ read relates to, so we don't Deep.
    // case: All matching entries are already used by previous FASTQ lines - so we have multiple FASTQ lines claiming the same
    //       SAM entry - since we can't undo the previous lines at this point - we fail the execution and re-compress without deep
    ZipZDeep *matching_ent = NULL;

    // match by QNAME1 or QNAME2 ; SEQ must match ; QUAL must match unless deep_no_qual ; trimming allow if deep_has_trimmed
    reason = fastq_seg_find_deep (vb, dl, &deep_hash, STRa(seq), STRa(qual), &matching_ent);

    // case: single matching ent
    if (matching_ent) {
        deep_value = 1 + fastq_seg_deep_consume_unique_matching_ent (vb, matching_ent, deep_hash, STRa(qname), STRa(seq), STRa(qual)); // +1 as 0 means "no match"
        vb->deep_stats[reason]++;

        dl->sam_seq_len = matching_ent->seq_len;
    }

    // case: no matching ent
    else {
        no_match: 
        vb->deep_stats[reason]++;

        if (reason == NDP_NO_MATCH && flag.show_deep && (flag.show_deep != SHOW_DEEP_ONE_HASH || deephash_issame (flag.debug_deep_hash, deep_hash))) {
            static int count = 0;

            if (flag.show_deep != SHOW_DEEP_SUMMARY || __atomic_fetch_add (&count, (int)1, __ATOMIC_RELAXED) < 20) {
                rom once = "";
                DO_ONCE once = (flag.show_deep == SHOW_DEEP_SUMMARY ? "\n\n(Showing first ~20 no-matches)\n" : "");
                iprintf ("%s%s: no matching SAM line (deep_hash=%016"PRIx64",%08x,%08x is_last=%s trimming=%s)\nQNAME=\"%.*s\"\nSEQ  =\"%.*s\"\nQUAL =\"%.*s\"\n", 
                         once, LN_NAME, DEEPHASHf(deep_hash), YN(segconf.deep_is_last), TF(segconf.deep_has_trimmed), STRf(qname), STRf(seq), STRf(qual));        
            }
        }
    }

do_seg:
    fastq_seg_deep_do (vb, deep_value);

done:
    if (segconf.deep_N_sam_score && !segconf.deep_no_qual) 
        buf_free (vb->scratch);

    COPY_TIMER (fastq_seg_deep);

    return deep_value > 0;
}

void fastq_deep_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(seq), ContextP bitmap_ctx, ContextP nonref_ctx)
{
    uint32_t trim_len = seq_len - dl->sam_seq_len; // left trim + right trim;

    if (!trim_len)
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_deep_copy_SEQ }, 2, bitmap_ctx, seq_len); 
    
    else {
        SNIPi2_2 (SNIP_SPECIAL, FASTQ_SPECIAL_deep_copy_SEQ, trim_len, dl->sam_seq_offset);
        seg_by_ctx (VB, STRa(snip), bitmap_ctx, dl->sam_seq_len); 

        if (dl->sam_seq_offset) // left trim
            seg_add_to_local_fixed (VB, nonref_ctx, seq, dl->sam_seq_offset, LOOKUP_NONE, 0);

        if (dl->sam_seq_offset + dl->sam_seq_len < seq_len) // right trim
            seg_add_to_local_fixed (VB, nonref_ctx, &seq[dl->sam_seq_offset + dl->sam_seq_len], trim_len - dl->sam_seq_offset, LOOKUP_NONE, 0);
    
        nonref_ctx->txt_len += trim_len;
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
        seg_by_did (VB, qname + qname_len - uncanonical_suffix_len, uncanonical_suffix_len, (did_i == FASTQ_QNAME) ? FASTQ_QmNAME : FASTQ_QmNAME2, uncanonical_suffix_len);
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
    // param=1 during SAM VB reconstruction and param=2 after z_file->deep_ents are ready
    ASSERT0 (z_file->vb_start_deep_line.param, "deep_ent not initialized"); // SAM VBs were not reconstructed - perhaps writer_z_initialize issue  

    // busy-wait until sam_piz_deep_finalize_ents finalizes the deep data
    while (z_file->vb_start_deep_line.param != 2) 
        usleep (100000); // 100ms

    // flag.no_writer (eg --test) with gencomp files is forced in-order for SAM bc digest is calculated by writer for gencomp files,
    // but can be out of order in the FQ part as if no_writer=true FASTQ VBs have needs_write=false (see writer_z_initialize)
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

    uint64_t pair_1_deep_value = vb->R1_vb_i ? fastq_get_pair_deep_value (vb, ctx) : 0; // consume whether or not used

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

    new_value->p = B8(*deep_ents, ents_index);

    COPY_TIMER (fastq_special_set_deep);
    return HAS_NEW_VALUE; 
}

SPECIAL_RECONSTRUCTOR (fastq_special_deep_copy_QNAME)
{
    START_TIMER;

    if (flag.qual_only || flag.seq_only) goto done; // if seq/qual_only, SAM did not add qname to deep_ents

    uint8_t *deep_ent = CTX(FASTQ_DEEP)->last_value.p;
    ASSERTNOTNULL (deep_ent);
        
    int qname_uncomp_len = deep_ent[0];
    bytes qname_comp = &deep_ent[1];

    uint32_t qname_comp_len = reconstruct ? RECONSTRUCT_huffman_or_copy (vb, SAM_QNAME, qname_uncomp_len, qname_comp)
                                          : huffman_uncompress_len (SAM_QNAME, qname_comp, qname_uncomp_len);

    CTX(FASTQ_DEEP)->last_value.p += 1 + qname_comp_len;
    
done:
    COPY_TIMER (fastq_special_deep_copy_QNAME);
    return NO_NEW_VALUE;    
}

#define GET_NUMBER_1B_or_5B ({                                  \
    uint32_t value = *next++;                                   \
    if (value == 255) { value = GET_UINT32 (next); next += 4; } \
    value; })

static uint8_t *fastq_special_deep_copy_SEQ_by_arith (VBlockFASTQP vb, uint8_t *next)
{
    uint32_t comp_len = GET_NUMBER_1B_or_5B;

    unsigned out_len = vb->deep_seq_len;        
    ASSPIZ (arith_uncompress_to (VB, next, comp_len, (uint8_t*)BAFTtxt, &out_len) && out_len == vb->deep_seq_len,
            "Failed arith_uncompress_to SEQ copied from SAM: compressed_len=%u uncompressed_len=%u (expected: %u)", 
            comp_len, out_len, vb->deep_seq_len);
    
    return next + comp_len;
}

static uint8_t *fastq_special_deep_copy_SEQ_by_ref (VBlockFASTQP vb, uint8_t *next, uint32_t deep_ref_consumed)
{
    char *recon = BAFTtxt;

    PosType64 gpos = DF.is_long_gpos ? GET_UINT64 (next) : GET_UINT32 (next);
    next += DF.is_long_gpos ? 8 : 4;  
    
    const Bits *genome=NULL;
    ref_get_genome (gref, &genome, NULL, NULL);
    decl_acgt_decode;
    
    BASE_ITER_INIT (genome, gpos, deep_ref_consumed, DF.is_forward);

    // case: no cigar - just copy the reference
    if (!DF.has_cigar) 
        for (uint32_t i=0; i < vb->deep_seq_len; i++) 
            *recon++ = acgt_decode(BASE_NEXT);                    
    
    // case: follow instructions from cigar
    else {
        #define get_buf(buf, name)                                  \
            ASSERTNOTINUSE (vb->buf);                               \
            uint32_t name##_len = GET_NUMBER_1B_or_5B;              \
            buf_alloc_exact (vb, vb->buf, name##_len, char, #buf);  \
            char *name = B1STc(vb->buf);                            \

        get_buf (scratch, cigar);

        next += huffman_uncompress_or_copy (SAM_CIGAR, next, STRa(cigar)); // note: next is advanced to point to nonref data

        get_buf (codec_bufs[0], nonref);

        if (nonref_len) next += str_unpack_bases (nonref, next, nonref_len);

        rom after = cigar + cigar_len;
        rom c = cigar;
        uint32_t D_and_N_bases = 0, S_and_I_bases = 0;

        while (c < after) {
            ASSERT (IS_DIGIT(*c), "%s: bad cigar: %.*s", LN_NAME, STRf(cigar)); // must have at least one digit (may be 0)

            uint32_t n=0; 
            while (IS_DIGIT(*c) && (c < after-1)) { 
                n = 10*n + *c - '0'; 
                c++; 
            }
            
            switch (*c++) {
                case 'M' : case '=' : case 'X' :
                    for (uint32_t i=0; i < n; i++) *recon++ = acgt_decode(BASE_NEXT);                    
                    break;

                case 'S' : case 'I':
                    for (uint32_t i=0; i < n; i++) *recon++ = *nonref++;
                    S_and_I_bases += n;
                    break;

                case 'D' : case 'N':                    
                    for (uint32_t i=0; i < n; i++) BASE_NEXT; // TODO: not very effecient to do this in a loop... we should advance with a direct calculation
                    D_and_N_bases += n;
                    break;

                case 'P' : case 'H': break;

                default: ABORT ("%s: bad cigar: %.*s", LN_NAME, STRf(cigar));
            }
        }

        vb->deep_seq_len = deep_ref_consumed + S_and_I_bases - D_and_N_bases; // finally, we can calculate deep_seq_len
        
        buf_free (vb->scratch);
        buf_free (vb->codec_bufs[0]);
    }

    // overwrite the mismatching bases
    if (DF.has_mismatches) {
        uint32_t num_mismatches = GET_NUMBER_1B_or_5B;
        recon = BAFTtxt;

        for (uint32_t i=0; i < num_mismatches; i++) {
            uint32_t skip = GET_NUMBER_1B_or_5B;
            ASSERT (recon + skip < Btxt(vb->txt_data.size), "%s: skip=%u goes beyond txt_data", LN_NAME, skip);

            recon += skip; // offset - skip this many matching bases
            *recon++ = *next++; // overwrite reference base with mismatching base
        }
    }

    return next;
}

SPECIAL_RECONSTRUCTOR_DT (fastq_special_deep_copy_SEQ)
{
    VBlockFASTQP vb = (VBlockFASTQP)vb_;

    START_TIMER;

    uint8_t *next = CTX(FASTQ_DEEP)->last_value.p; // points to start of SEQ stuff within deep_ents

    if (flag.header_only_fast) goto done; // if --header-only, SAM did not add seq to deep_ents

    uint32_t trim_len = 0;   // how much longer is the FASTQ SEQ than the SAM SEQ

    // 3 cases:
    // 0 items - no trimming
    // 1 item  - trim_len - up to 15.0.22 - trimming with offset 0
    // 2 items - (trim_len, seq_offset), from 15.0.23 - trimming supports offset >= 0 (i.e. left-trimming)
    if (snip_len) {
        str_split_ints (snip, snip_len, 2, ',', item, false);
        trim_len = items[0];
        if (n_items == 2) vb->sam_seq_offset = items[1];
    }

    ASSERTNOTNULL (next);

    vb->piz_deep_flags = *(PizDeepSeqFlags *)next++;
    
    if (flag.qual_only) {
        vb->deep_seq_len = GET_NUMBER_1B_or_5B;  // if qual-only, SAM stores deep_seq_len in deep_ents rather than ref_consumed
        vb->seq_len = vb->deep_seq_len + trim_len; 
        goto done;
    }

    // reconstruct left trim
    if (vb->sam_seq_offset)
        reconstruct_from_local_sequence (VB, CTX(FASTQ_NONREF), vb->sam_seq_offset, reconstruct);

    uint32_t deep_ref_consumed = GET_NUMBER_1B_or_5B; 

    if (!vb->piz_deep_flags.has_cigar)        
        vb->deep_seq_len = deep_ref_consumed; // note: if has_cigar, deep_seq_len is set by fastq_special_deep_copy_SEQ_by_ref based on cigar

    // reconstruct SEQ copied from SAM
    if (reconstruct) {
        #define ENCODING(x) (DF.seq_encoding == ZDEEP_SEQ_##x)
        next = ENCODING(BY_REF) ? fastq_special_deep_copy_SEQ_by_ref (vb, next, deep_ref_consumed)
             : ENCODING(PACKED) ? next + str_unpack_bases (BAFTtxt, next, vb->deep_seq_len)
             : /* ARITH */        fastq_special_deep_copy_SEQ_by_arith (vb, next);
        Ltxt += vb->deep_seq_len;
    }
    
    vb->seq_len = vb->deep_seq_len + trim_len; 

    ASSPIZ (vb->seq_len <= vb->longest_seq_len, "%s: unexpectedly seq_len(as read from deep_ent)=%u > vb->longest_seq_len=%u (deep_seq_len=%u trim_len=%u sam_seq_offset=%u)",
            LN_NAME, vb->seq_len, vb->longest_seq_len, vb->deep_seq_len, trim_len, vb->sam_seq_offset); // sanity check

    // reconstruct right trim
    if (trim_len - vb->sam_seq_offset)
        reconstruct_from_local_sequence (VB, CTX(FASTQ_NONREF), trim_len - vb->sam_seq_offset, reconstruct);

    // case we are pair-2: advance pair-1 SQBITMAP iterator, and if pair-1 is aligned - also its GPOS iterator
    if (vb->R1_vb_i/*we are R2*/ && fastq_piz_R1_test_aligned (vb))
        CTX(FASTQ_GPOS)->localR1.next++; // gpos_ctx->localR1.next is an iterator for both gpos and strand

done:
    CTX(FASTQ_DEEP)->last_value.p = next;
    
    COPY_TIMER (fastq_special_deep_copy_SEQ);
    return NO_NEW_VALUE;    
}

SPECIAL_RECONSTRUCTOR_DT (fastq_special_deep_copy_QUAL)
{
    VBlockFASTQP vb = (VBlockFASTQP)vb_;
    char *qual = BAFTtxt; // save recon start address

    START_TIMER;

    if (flag.header_only_fast || flag.seq_only) goto done; // if --header-only or --seq-only, SAM did not add seq to deep_ents

    // reconstruct left and right trim together (see in zip: fastq_zip_qual)
    ASSERT (vb->seq_len >= vb->deep_seq_len, "Expecting vb->seq_len=%u >= deep_seq_len=%u", vb->seq_len, vb->deep_seq_len);
    uint32_t total_trim = vb->seq_len - vb->deep_seq_len;

    ASSERT (total_trim >= vb->sam_seq_offset, "Expecting total_trim=%u >= vb->sam_seq_offset=%u", total_trim, vb->sam_seq_offset);
    uint32_t right_trim = total_trim - vb->sam_seq_offset; 

    // reconstruct left and right trim, or just consume local data even if !reconstruct
    if (total_trim) { 
        if (ctx->ltype == LT_CODEC)         
            codec_args[ctx->lcodec].reconstruct (VB, ctx->lcodec, ctx, total_trim, reconstruct); 

        else if (ctx->ltype == LT_BLOB) 
            reconstruct_from_local_sequence (VB, ctx, total_trim, reconstruct); 

        else
            ABORT_PIZ ("Invalid ltype=%s for %s", lt_name (ctx->ltype), ctx->tag_name);
    }
    
    if (!reconstruct) goto done;

    // move right trim to its place after the portion to be copied from Deep
    if (right_trim) {
        Ltxt -= right_trim;
        memmove (BAFTtxt + vb->deep_seq_len, BAFTtxt, right_trim); // careful: if there is any QUAL codec that temporarily writes beyond, this won't work
    }

    uint8_t *next = CTX(FASTQ_DEEP)->last_value.p;
    ASSERTNOTNULL (next);

    // case: monochar qual
    if (vb->piz_deep_flags.qual_is_monochar) 
        memset (BAFTtxt, *next, vb->deep_seq_len);

    // case: reconstruct QUAL copied from SAM
    else {
        uint32_t qual_comp_len = GET_NUMBER_1B_or_5B;
        ASSPIZ0 (qual_comp_len, "qual_comp_len=0 but  it to be >0");

        unsigned out_len = vb->deep_seq_len;        
        if (huffman_exists (SAM_QUAL))
            huffman_uncompress (SAM_QUAL, next, BAFTtxt, out_len);
        else
            ASSPIZ (arith_uncompress_to (VB, next, qual_comp_len, (uint8_t *)BAFTtxt, &out_len) && out_len == vb->deep_seq_len,
                    "Failed arith_uncompress_to QUAL copied from SAM: compressed_len=%u uncompressed_len=%u (expected: %u)", 
                    qual_comp_len, out_len, vb->deep_seq_len);
    }

    Ltxt += vb->deep_seq_len + right_trim;

    // case: correct qualities of 'N' bases
    if (segconf.deep_N_fq_score) {
        STRlast (seq, FASTQ_SQBITMAP);

        for (unsigned i=0; i < seq_len; i++)
            if (seq[i] == 'N') qual[i] = segconf.deep_N_fq_score;
    }

done:
    COPY_TIMER (fastq_special_deep_copy_QUAL);
    return NO_NEW_VALUE;
}
