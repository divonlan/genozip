// ------------------------------------------------------------------
//   deep_sam.c
//   Copyright (C) 2023-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#include "threads.h"
#include "huffman.h"
#include "refhash.h"
#include "htscodecs/arith_dynamic.h"

#define PUT_NUMBER_1B_or_5B(n) ({                   \
    if (__builtin_expect ((n) <= 254, true))        \
        *next++ = (n);                              \
    else {                                          \
        *next++ = 255;                              \
        PUT_UINT32 (next, (n));                     \
        next += 4;                                  \
    } })

// --------
// ZIP side
// --------

static void sam_deep_zip_show_index_stats (void) 
{
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);
    ARRAY (uint32_t, index, z_file->deep_index);

    iprintf ("\nz_file.deep_index.len=%"PRIu64" (%s) z_file.deep_ents.len=%"PRIu64" (%s)\n", 
             z_file->deep_index.len, str_size (z_file->deep_index.len * sizeof (uint32_t)).s, 
             z_file->deep_ents.len,  str_size (z_file->deep_ents.len  * sizeof (ZipZDeep)).s);

    #define NUM_LENS 64
    uint64_t used=0;
    uint64_t lens[NUM_LENS] = {}, longer=0; // each entry i contains number linked lists which are of length i; last entry i - linked lists of i or longer
    uint64_t longest_len = 0;
    uint32_t longest_len_hash;

    uint64_t total_ents_traversed = 0; // total ents traversed on linked lists when segging FASTQ reads

    for (uint32_t hash=0; hash < index_len; hash++) {
        uint32_t this_len=0;
        for (uint32_t ent_i = index[hash]; ent_i != NO_NEXT; ent_i = deep_ents[ent_i].next) {
            ASSERT (ent_i < z_file->deep_ents.len, "ent_i=%u >= deep_ents.len=%u in element %u linked list of hash=%u", 
                    ent_i, z_file->deep_ents.len32, this_len, hash);
            this_len++;
        }
        
        if (this_len) used++;

        if (this_len < NUM_LENS) 
            lens[this_len]++;

        else {
            longer++;
            total_ents_traversed += this_len        // number of entries that will access this linked list
                                  * this_len;       // length of linked list traversed for each entry accessing this list

            if (this_len > longest_len) {
                longest_len = this_len;
                longest_len_hash = hash;
            }
        }
    }

    if (longest_len > 16384)
        iprintf ("\nFYI: deep_index entry %u has %"PRIu64 " deep entries on its linked list. This is excessively large and might cause slowness.%s",
                 longest_len_hash, longest_len, report_support_if_unexpected()); // note: error only shows with --show-deep... not very useful

    iprintf ("\ndeep_index entries used: %"PRIu64" / %"PRIu64" (%.1f%%)\n", 
             used, index_len, 100.0 * (double)used / (double)index_len);

    iprint0 ("\ndeep_index linked list length histogram:\n");
    for (int this_len=0; this_len < NUM_LENS; this_len++) {
        if (!lens[this_len]) continue;

        iprintf ("%3u: %"PRIu64"\n", this_len, lens[this_len]);

        // our algorithm requires traversal of entire linked list for each entry. calculate the length of the list
        // traversed for the average read (simplifying assumption: there are the same)
        total_ents_traversed += (uint64_t)this_len        // number of entries that will access this linked list
                              * (uint64_t)this_len        // length of linked list traversed for each entry accessing this list
                              * (uint64_t)lens[this_len]; // number of linked lists of this length
    }

    if (longer) iprintf ("longer: %"PRIu64"\n", longer);

    // this is approximate: assuming # of ents hashed is equal to the number of reads 
    iprintf ("\nFASTQ SEG: deep_index linked list length traversed on average: %.2f (total_ents_traversed=%"PRIu64")\n", 
             (double)total_ents_traversed / (double)deep_ents_len, total_ents_traversed);
}

static void sam_deep_zip_display_reasons (void)
{
    uint64_t total = z_file->num_lines;

    iprintf ("%s Alignments breakdown by deepability:\n", dt_name (txt_file->data_type));
    iprintf ("%-13.13s: %"PRIu64"\n", "Total", total);

    for (int i=0; i < NUM_DEEP_STATS_ZIP; i++) 
        if (z_file->deep_stats[i])
            iprintf ("%-13.13s: %"PRIu64" (%.1f%%)\n", (rom[])DEEP_STATS_NAMES_ZIP[i], z_file->deep_stats[i], 100.0 * (double)z_file->deep_stats[i] / (double)total);
}

// main thread: Called during zip_finalize of the SAM component for a Deep compression
void sam_deep_zip_finalize (void)
{
    threads_log_by_vb (evb, "main_thread", "sam_deep_zip_finalize", 0);
    
    if (zip_is_biopsy) return;

    // return unused memory to libc
    buf_trim (z_file->deep_ents, ZipZDeep);
    
    // Build vb_start_deep_line: first deep_line (0-based SAM-wide line_i, but not counting SUPP/SEC lines, (up 15.0.68: monochar reads), SEQ.len=0) of each SAM VB
    buf_alloc_exact_zero (evb, z_file->vb_start_deep_line, z_file->num_vbs + 1, uint64_t, "z_file->vb_start_deep_line");

    uint64_t count=0;
    
    for (VBIType vb_i=1; vb_i < z_file->vb_num_deep_lines.len32; vb_i++) { // note: z_file->vb_num_deep_lines.len might be less than z_file->num_vbs, as we don't extend the buffer for DEPN VBs
        *B64(z_file->vb_start_deep_line, vb_i) = count; // note: in ZIP, vb_start_deep_line index is vb_i (unlike PIZ, which is vb_idx)        
        count += *B32(z_file->vb_num_deep_lines, vb_i);
    }

    if (flag_show_deep) {
        sam_deep_zip_display_reasons();
        sam_deep_zip_show_index_stats();
    }
        
    memset (z_file->deep_stats, 0, sizeof (z_file->deep_stats)); // we will reuse this field to gather stats while segging the FASTQ files

    // return a bunch of memory to the kernel before moving on to FASTQ
    buf_destroy (z_file->vb_num_deep_lines);
    vb_dehoard_memory (true);
}

// SAM seg: calculate hash values into dl->deep_*

// ZIP compute thread while segging SEQ
void sam_deep_set_QNAME_hash (VBlockSAMP vb, ZipDataLineSAMP dl, QType q, STRp(qname))
{
    // note: we must drop consensus reads, because Deep has an assumption, used to prove uniqueness by hash,
    // that there are no BAM alignments except those in the FASTQ files provided
    if (!dl->FLAG.secondary && !dl->FLAG.supplementary && !segconf.flav_prop[q].is_consensus) {
        dl->deep_hash.qname = deep_qname_hash (VB, q, STRa(qname), segconf.deep_paired_qname ? dl->FLAG.is_last : unknown, NULL);
        dl->is_deepable = true; // note: we haven't tested SEQ yet, so this might still change

        vb->lines.count++;      // counts deepable lines
    }

    else
        vb->deep_stats[segconf.flav_prop[q].is_consensus ? RSN_CONSENSUS 
                     : dl->FLAG.secondary                ? RSN_SECONDARY 
                     :                                     RSN_SUPPLEMENTARY]++; // counting non-deepability reasons
}

// ZIP compute thread while segging SEQ: set hash of forward SEQ (i.e. as it would appear in the FASTQ file) 
void sam_deep_set_SEQ_hash (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(textual_seq))
{
    if (!dl->is_deepable) return; 
    
    // case: sam_deep_set_QNAME_hash happened, but now we discover that SEQ is missing 
    // (up tp 15.0.68 this was: or mono-char seq)
    if (IS_ASTERISK(textual_seq)) {     // note: we don't test SEQ.len, bc SEQ.len=0 for unmapped reads, but they are still deepable   
        dl->is_deepable = false;
        vb->lines.count--;

        vb->deep_stats[RSN_NO_SEQ]++; // counting non-deepability reasons        
    }

    else {
        dl->deep_hash.seq = deep_seq_hash (VB, STRa(textual_seq), dl->FLAG.rev_comp);
        vb->deep_stats[RSN_DEEPABLE]++;
    }
}

// ZIP: compute thread while segging QUAL: hash of forward QUAL (i.e. as it would appear in the FASTQ file) for Deep 
void sam_deep_set_QUAL_hash (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(qual))
{
    if (!dl->is_deepable || 
        !vb->has_qual || segconf.has_bqsr) return; // case: deepable, but without QUAL: dl->deep_qual_hash remains 0

    dl->deep_hash.qual = deep_qual_hash (VB, STRa(qual), dl->FLAG.rev_comp); 

    if (flag.show_deep == SHOW_DEEP_ALL || 
        (flag.show_deep == SHOW_DEEP_ONE_HASH && deephash_issame (dl->deep_hash, flag.debug_deep_hash)))
        iprintf ("\n%s Recording SAM alignment: deep_hash=%016"PRIx64",%08x,%08x\nQNAME=\"%.*s\"\nSEQ=\"%.*s\"\nQUAL=\"%.*s\"\n",
                 LN_NAME, DEEPHASHf(dl->deep_hash), dl->QNAME_len, dl_qname(dl), dl->SEQ.len, vb->textual_seq_str, STRfw(dl->QUAL));
}

// deep and bamass
uint64_t sam_deep_calc_hash_bits (void)
{
    uint64_t est_num_lines = segconf.line_len ? (txt_file->est_seggable_size / segconf.line_len) : 10000;

    // number of bits of hash table size (10-31) (note: minimum cannot be less than 10 due to BamAssEnt.z_qname_hash_hi)
    int clzll = (int)__builtin_clzll (est_num_lines * 2); 
        num_hash_bits = MAX_(10, MIN_(31, 63 - clzll));
    
    if (flag.low_memory)
        num_hash_bits = MAX_(10, num_hash_bits - 3); // 8x smaller index at the expense of longer linked lists

    return est_num_lines;
}

// ZIP compute thread: mutex-protected callback from ctx_merge_in_vb_ctx during merge: add VB's deep_hash to z_file->deep_index/deep_ents
void sam_deep_zip_merge (VBlockP vb_)
{
    if (!flag.deep || zip_is_biopsy) return;

    VBlockSAMP vb = (VBlockSAMP)vb_;
    START_TIMER;

    // initialize - first VB merging
    if (!z_file->deep_index.len) {
        // pointers into deep_ents - initialize to NO_NEXT
        uint64_t est_num_lines = sam_deep_calc_hash_bits();

        // pointers into deep_ents - initialize to NO_NEXT
        z_file->deep_index.can_be_big = true;
        buf_alloc_exact_255 (evb, z_file->deep_index, ((uint64_t)1 << num_hash_bits), uint32_t, "z_file->deep_index"); 

        // note: no initialization of deep_ents - not needed and it takes too long (many seconds to get pages from kernel) while all threads are waiting for vb=1
        z_file->deep_ents.can_be_big = true;
        buf_alloc_exact (evb, z_file->deep_ents, MAX_(1.1 * est_num_lines, vb->lines.count), ZipZDeep, "z_file->deep_ents");

        if (flag_show_deep) 
            printf ("num_hash_bits=%u est_num_lines=%"PRIu64" deep_index.len=%"PRIu64" deep_ents.len=%"PRIu64"\n", 
                    num_hash_bits, est_num_lines, z_file->deep_index.len, z_file->deep_ents.len);

        z_file->deep_ents.len = 0;
    }

    if (vb->comp_i != SAM_COMP_DEPN) {
        buf_alloc_zero (NULL, &z_file->vb_num_deep_lines, 0, vb->vblock_i+1, uint32_t, 0, NULL); // initial allocation is in sam_set_sag_type
        *B32(z_file->vb_num_deep_lines, vb->vblock_i) = vb->lines.count;
        MAXIMIZE (z_file->vb_num_deep_lines.len32, vb->vblock_i+1);
    }

    uint32_t ent_i = z_file->deep_ents.len;

    // we can only handle 4G entries, because hash entry and "next" fields are 32 bit.
    int64_t excess = (uint64_t)ent_i + vb->deep_stats[RSN_DEEPABLE] - MAX_DEEP_ENTS;  
    if (excess > 0) {
        WARN_ONCE ("The number of deepable alignments in %s exceeds %u - the maximum supported for Deep compression. The excess alignments will be compressed normally without Deep. %s", 
                   txt_name, MAX_DEEP_ENTS, report_support());

        vb->deep_stats[RSN_DEEPABLE] -= excess;
        vb->deep_stats[RSN_OVERFLOW] += excess;

        if (!vb->deep_stats[RSN_DEEPABLE]) goto done; // cannot add *any* alignment from this VB
    }

    buf_alloc (evb, &z_file->deep_ents, vb->lines.count/*num deepable lines in VB*/, 0, ZipZDeep, CTX_GROWTH, "z_file->deep_ents"); 

    uint32_t *deep_index = B1ST32 (z_file->deep_index);
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    uint64_t mask = bitmask64 (num_hash_bits); // num_hash_bits is calculated upon first merge
    uint32_t deep_line_i = 0; // line_i within the VB, but counting only deepable lines that have SEQ.len > 0 and are not monochar
    for_buf (ZipDataLineSAM, dl, vb->lines) {
        if (!dl->is_deepable) continue; // case: non-deepable: secondary or supplementary line, no sequence etc

        uint32_t hash = dl->deep_hash.qname & mask;                                 
        
        // create new entry - which is now the head of linked list
        deep_ents[ent_i] = (ZipZDeep){ .next    = deep_index[hash], // previous head of linked list is is now the 2nd element on the list
                                       .seq_len = dl->SEQ.len, 
                                       .hash    = dl->deep_hash,
                                       .vb_i    = vb->vblock_i, 
                                       .line_i  = deep_line_i++ };

        deep_index[hash] = (uint32_t)ent_i; // this is entry is now head of linked list 

        if (++ent_i > MAX_DEEP_ENTS) break;
    }

    z_file->deep_ents.len = ent_i;

done:
    for (int i=0; i < NUM_DEEP_STATS_ZIP; i++)
        z_file->deep_stats[i] += vb->deep_stats[i];

    COPY_TIMER(sam_deep_zip_merge);
} 

//-----------------------------------------------------------------------------
// SAM PIZ: step 1: output QNAME,SEQ,QUAL into vb->deep_ents / vb->deep_index
//-----------------------------------------------------------------------------

// main thread: called for after reading global area
void sam_piz_deep_initialize (void)
{
    store_release (z_file->vb_start_deep_line.param, DEEP_PIZ_ENTS_RECONSTRUCTING);  
}

// called during reconstruction of SEQ to determine if SEQ can be stored vs reference 
void sam_piz_set_deep_seq (VBlockSAMP vb, 
                           bool not_end_of_contig, // doesn't round robin at edge of the chromosome
                           PosType64 gpos)         // true if forward relative to the reference
{
    vb->deep_stats[EXPL_DEEPABLE]++;

    #define NOT_BY_REF(reason) ({ vb->deep_stats[reason]++; \
                                  vb->piz_deep_flags.has_cigar = false; \
                                  return; }) // vb->piz_deep_flags.seq_encoding remains the default ZDEEP_SEQ_PACKED

    if (gpos == NO_GPOS64) 
        NOT_BY_REF(EXPL_SEQ_COPY_VERBATIM);

    if (vb->bisulfite_strand)
        NOT_BY_REF(EXPL_BISULFITE);

    ARRAY32 (char, deep_nonref, CTX(SAM_NONREF)->deep_nonref);

    bool used_aligner = ctx_has_value_in_line_(vb, CTX(SAM_STRAND));

    ASSERT (used_aligner || vb->binary_cigar.len32, "%s: missing cigar", LN_NAME);

    // note: we don't need to store the cigar if it is "perfect" (eg "151M" or aligner used) 
    vb->piz_deep_flags.has_cigar = !used_aligner && !(vb->binary_cigar.len32 == 1 && B1ST(BamCigarOp, vb->binary_cigar)->op == BC_M); 

    if (vb->piz_deep_flags.has_cigar && !str_is_ACGT (STRa(deep_nonref), NULL))
        NOT_BY_REF(EXPL_SEQ_NONREF_NON_ACGT);
    
    if (!not_end_of_contig) 
        NOT_BY_REF(EXPL_SEQ_END_OF_CONTIG);

    if (vb->num_deep_mismatches > MAX_DEEP_SEQ_MISMATCHES) 
        NOT_BY_REF(EXPL_TOO_MANY_MISMATCHES);

    int32_t ref_and_seq_consumed = used_aligner ? vb->seq_len : vb->ref_and_seq_consumed;

    // sanity
    ASSERT (deep_nonref_len == vb->seq_len - ref_and_seq_consumed, "%s: expecting deep_nonref_len=%u == seq_len=%u - ref_and_seq_consumed=%u",
            LN_NAME, deep_nonref_len, vb->seq_len, ref_and_seq_consumed);

    // Yes! we will store SEQ as encoded vs the reference
    
    // calculate is_forward
    bool aligner_rev = used_aligner && !CTX(SAM_STRAND)->last_value.i;

    vb->piz_deep_flags.is_forward     = (last_flags.rev_comp == aligner_rev); // true if FASTQ's sequence is forward relative to the reference
    vb->piz_deep_flags.has_mismatches = (vb->num_deep_mismatches > 0);

    // revcomp the mismatches if FASTQ's and SAM's seq are reversed (regardless of orientation vs reference) 
    if (last_flags.rev_comp) {
        for (int i=0; i < vb->num_deep_mismatches / 2; i++) {
            SWAP (vb->deep_mismatch_offset[i], vb->deep_mismatch_offset[vb->num_deep_mismatches-1-i]);
            SWAP (vb->deep_mismatch_base[i], vb->deep_mismatch_base[vb->num_deep_mismatches-1-i]);
        }

        for (int i=0; i < vb->num_deep_mismatches; i++) {
            vb->deep_mismatch_offset[i] = (vb->seq_len - 1 - vb->deep_mismatch_offset[i]);
            vb->deep_mismatch_base[i] = COMPLEM[(uint8_t)vb->deep_mismatch_base[i]];
        }
    }

    // make offsets relative to each other: offset will now become the number of matching bases (within ref_consumed)
    // since the previous offset or beginning of ref_consumed
    for (int i=vb->num_deep_mismatches-1; i >= 1; i--) 
        vb->deep_mismatch_offset[i] -= vb->deep_mismatch_offset[i-1] + 1;

    vb->piz_deep_flags.seq_encoding = ZDEEP_SEQ_BY_REF; // note: if we didn't reach here, encoding remains the default ZDEEP_SEQ_PACKED
    vb->deep_gpos = gpos;

    // collapse X,= to M ; replace N with D ; I with S ; remove H,P ; reverse if rev_comp 
    if (vb->piz_deep_flags.has_cigar) 
        sam_prepare_deep_cigar (VB, CIG(vb->binary_cigar), last_flags.rev_comp); 

    vb->deep_stats[EXPL_SEQ_AS_REF]++;
}

// strive to alloc exactly twice: 128K initialy, and then a good guess of the total
static void inline sam_piz_alloc_deep_ents (VBlockSAMP vb, uint32_t size)
{
    // case: initial allocation
    if (!vb->deep_ents.len) {
        buf_alloc (vb, &vb->deep_ents, size, 128 KB, uint8_t, 0, "deep_ents"); 
        buf_alloc (vb, &vb->deep_index, 0, vb->lines.len, uint32_t, 0, "deep_index"); 
    }

    // case: have enough memory, no need to allocate
    else if (vb->deep_ents.len + size <= vb->deep_ents.size)
        return;

    // case: we already have stats from a few lines - allocate now for the entire VB
    else {
        uint32_t avg_per_line_so_far = vb->deep_index.len ? (vb->deep_ents.len32 / vb->deep_index.len32) : 0;
        buf_alloc (vb, &vb->deep_ents, size, avg_per_line_so_far * vb->lines.len, uint8_t, 1.1, "deep_ents"); 
    }
}

// copy QNAME (huffman-compressed or uncompressed) to deep_ents
static void sam_piz_deep_add_qname (VBlockSAMP vb)
{
    START_TIMER;

    STRlast (qname, IS_PRIM(vb) ? SAM_QNAMESA : SAM_QNAME);

    BNXT32(vb->deep_index) = vb->deep_ents.len32; // index in deep_ents of data of current deepable line
         
    if (segconf.deep_qtype == QNONE || // back comp - QNONE was possible up to 15.0.66 (deep match was by non-trimmed SEQ and QUAL)
        flag.seq_only || flag.qual_only)
        goto done;         

    ASSPIZ (qname_len <= SAM_MAX_QNAME_LEN, "QNAME.len=%u, but per BAM specfication it is limited to %u characters. QNAME=\"%.*s\"", SAM_MAX_QNAME_LEN, qname_len, STRf(qname));

    qname_canonize (QNAME1, qSTRa(qname), vb->comp_i); // might reduce qname_len

    BNXT8(vb->deep_ents) = qname_len; // uncompressed length

    // since 15.0.65, a SEC_HUFFMAN section is available and we can compress. 
    if (huffman_exists (SAM_QNAME)) {
        uint32_t comp_len = huffman_get_theoretical_max_comp_len (SAM_QNAME, qname_len);
        uint8_t comp[comp_len];
        huffman_compress (VB, SAM_QNAME, STRa(qname), qSTRa(comp)); // using huffman sent via SEC_HUFFMAN

        buf_add (&vb->deep_ents, (rom)comp, comp_len);

        vb->deep_stats[QNAME_BYTES] += 2 + comp_len;
    }
    
    // For files up to 15.0.64. Note that the code used to build the huffman table on-the-fly here.- 
    // For simplicity, we no longer do that, and qname for older files shall remain uncompressed in memory.
    else { 
        buf_add (&vb->deep_ents, (rom)qname, qname_len);

        vb->deep_stats[QNAME_BYTES] += 2 + qname_len;
    }

done:
    COPY_TIMER (sam_piz_deep_add_qname);
}

// compress QUAL or SEQ - adding a 1B or 5B length
static uint32_t sam_piz_deep_compress (VBlockSAMP vb, uint8_t *next, STRp(data), bool is_seq)
{
    START_TIMER;
    uint8_t *start = next;

    // reverse if needed
    if (last_flags.rev_comp) {
        ASSERTNOTINUSE (vb->scratch);
        buf_alloc_exact (vb, vb->scratch, data_len, char, "scratch");

        if (is_seq) str_revcomp (B1STc(vb->scratch), STRa(data));
        else        str_reverse (B1STc(vb->scratch), STRa(data));
    
        data = B1STc(vb->scratch);
    }

    uint32_t comp_len = is_seq ? arith_compress_bound (data_len, X_NOSZ) :
                                 huffman_get_theoretical_max_comp_len (SAM_QUAL, data_len);
    
    bool len_is_1_byte = (data_len < (is_seq ? 1024 : 512)); // just a guess for now
    uint8_t *guess_addr = next + (len_is_1_byte ? 1 : 5);   // we don't know comp_len yet, so just a guess
    
    if (!is_seq)
        huffman_compress (VB, SAM_QUAL, STRa(data), guess_addr, &comp_len);
    else
        ASSPIZ (arith_compress_to (evb, (uint8_t *)STRa(data), guess_addr, &comp_len, X_NOSZ) && comp_len,
                "%s: Failed to arith-compress data: is_seq=%s comp_len=%u data_len=%u data=\"%.*s\"", 
                LN_NAME, TF(is_seq), comp_len, data_len, STRf(data));

    // we guessed wrong how many bytes the comp_len integer will take - move the data
    if (len_is_1_byte != (comp_len <= 254)) {
        memmove (next + (comp_len <= 254 ? 1 : 5), guess_addr, comp_len);
        len_is_1_byte = false;
    }

    PUT_NUMBER_1B_or_5B (comp_len);
    next += comp_len;

    if (last_flags.rev_comp)   
        buf_free (vb->scratch);

    COPY_TIMER (sam_piz_deep_compress);

    return next - start;
}

static uint8_t *sam_piz_deep_add_seq_by_ref (VBlockSAMP vb, uint8_t *next, PizDeepSeqFlags *deep_flags)
{
    // set GPOS and .is_long_gpos (4 or 5 bytes little endian)
    PUT_UINT32 (next, (uint32_t)vb->deep_gpos);
    next += sizeof (PosType32);

    if (vb->deep_gpos > 0xffffffffULL) {
        *next++ = (vb->deep_gpos >> 32);
        deep_flags->is_long_gpos = true;
    }

    // set cigar and nonref
    if (deep_flags->has_cigar) {
        next += nico_compress_cigar (VB, SAM_CIGAR, CIG(CTX(SAM_CIGAR)->deep_cigar), next, vb->deep_ents.size - BNUM64 (vb->deep_ents, next));

        BufferP deep_nonref = &CTX(SAM_NONREF)->deep_nonref;
        PUT_NUMBER_1B_or_5B (deep_nonref->len32); // always output if has_cigar, even if 0

        // note: sam_piz_set_deep_seq guarantees that NONREF is only A,C,G,T
        if (deep_nonref->len32)
            next += str_pack_bases (next, STRb(*deep_nonref), false);
    }

    // set mismatches
    if (vb->num_deep_mismatches) {
        PUT_NUMBER_1B_or_5B (vb->num_deep_mismatches);

        for (int i=0; i < vb->num_deep_mismatches; i++) {
            PUT_NUMBER_1B_or_5B (vb->deep_mismatch_offset[i]); // note: currently always 1B as number of mismatches is limited to 63.
            *next++ = vb->deep_mismatch_base[i];
        }
    }

    return next;
}

// container item callback for SQBITMAP: pack SEQ into 2-bit and into deep_ents or compress if it has non-ACGT (eg N)
static void sam_piz_deep_add_seq (VBlockSAMP vb, STRp(seq))
{
    START_TIMER;

    uint8_t *start = BAFT8(vb->deep_ents); 
    uint8_t *next  = start;

    vb->piz_deep_flags_index = vb->deep_ents.len32;
    PizDeepSeqFlags *deep_flags = (PizDeepSeqFlags *)next++;
    *deep_flags = vb->piz_deep_flags; // .seq_encoding and .is_forward were set already by sam_piz_set_deep_seq 

    // if has_cigar (and not qual-only): we store ref_consumed, otherwise seq_len
    PUT_NUMBER_1B_or_5B ((deep_flags->has_cigar && !flag.qual_only) ? vb->ref_consumed : seq_len); // used for SEQ and QUAL

    if (!flag.qual_only) {
        if (deep_flags->seq_encoding == ZDEEP_SEQ_BY_REF) // as determined by sam_piz_set_deep_seq
            next = sam_piz_deep_add_seq_by_ref (vb, next, deep_flags);

        // if only A,C,G,T - pack in 2-bits: set .seq_encoding + add packed seq
        else if (str_is_ACGT (STRa(seq), NULL)) 
            next += str_pack_bases (next, STRa(seq), last_flags.rev_comp);

        // if has an N or other IUPACs - arith-compress: set .seq_encoding, .is_long_seq_comp + add compressed seq
        else {        
            deep_flags->seq_encoding = ZDEEP_SEQ_ARITH;
            next += sam_piz_deep_compress (vb, next, STRa(seq), true);
        }
    }

    vb->deep_ents.len += (next - start);
    vb->deep_stats[deep_flags->seq_encoding - ZDEEP_SEQ_PACKED + SEQ_PACKED_BYTES] += next - start;

    COPY_TIMER (sam_piz_deep_add_seq);
}

// container item callback for QUAL/QUALSA: compress QUAL into deep_ents
static void sam_piz_deep_add_qual (VBlockSAMP vb, STRp(qual))
{
    START_TIMER;

    // case: monochar qual
    if (str_is_monochar (STRa(qual))) {
        ((PizDeepSeqFlags *)B8(vb->deep_ents, vb->piz_deep_flags_index))->qual_is_monochar = true;
        BNXT8(vb->deep_ents) = qual[0];

        vb->deep_stats[QUAL_MONOCHAR_BYTES]++;
    }

    // case: not monochar qual
    else {
        // back comp: old file with no SEC_HUFFMAN section for SAM_QUAL - create it based on first QUAL reconstructed
        if (!huffman_exists (SAM_QUAL)) // note: non-atomic: if it says it exists, then it exists. If it says it doesn't, it might or might not exist.
            huffman_piz_backcomp_produce_qual (STRa(qual)); // also handles thread-safety

        uint32_t comp_len = sam_piz_deep_compress (vb, BAFT8(vb->deep_ents), STRa(qual), false);
        vb->deep_ents.len32 += comp_len;

        vb->deep_stats[QUAL_BYTES] += comp_len;
    }


    COPY_TIMER (sam_piz_deep_add_qual);
}


// For SAM: called for top-level SAM container if compressed with --deep (defined as con_item_cb in data_types.h)
// for BAM: called from SEQ and QUAL translators (before translation) 
CONTAINER_ITEM_CALLBACK (sam_piz_con_item_cb)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    
    START_TIMER;

    switch (con_item->dict_id.num) {
        case _SAM_SQBITMAP:
            // backcomp note: these must be precisely the same deepability conditions as in ZIP
            if (flag.deep_sam_only                               || // Deep file (otherwise this callback will not be set), but SAM/BAM-only reconstruction
                last_flags.secondary || last_flags.supplementary || // SECONDARY or SUPPLEMENTARY alignment
                IS_ASTERISK(recon)                               || // Missing SEQ
                (z_file->flav_prop[vb->comp_i][QNAME2].is_consensus && ctx_encountered_in_line (VB, SAM_QNAME2)) || // consensus alignment 
                (!VER2(15,69) && str_is_monochar (STRa(recon))))    // mono-char SEQ is excluded in files generated up to 15.0.68
                
                vb->line_not_deepable = true;
            
            else {
                // note: we do QNAME during the SQBITMAP callback, because we need last_flags to be set 
                // (in SAM QNAME reconstructs before flags). On the other hand, SQBITMAP and QUAL must be 
                // handled right after reconstruction, because in BAM, immediately after they get converted to binary format.
                sam_piz_alloc_deep_ents (vb, 4 KB + 8 * recon_len + (1 + vb->binary_cigar.len32) * 6/*6=theo. max. of nico op(1)+N(5)*/); // more than enough for compresssed QNAME, SEQ and QUAL

                sam_piz_deep_add_qname (vb);
                if (!flag.header_only_fast) sam_piz_deep_add_seq (vb, STRa(recon));
            }
            break;
    
        case _SAM_QUAL:
        case _SAM_QUALSA:
            if (!vb->line_not_deepable && !vb->qual_missing && !segconf.deep_no_qual &&
                !flag.seq_only && !flag.header_only_fast) // conditions that may only appear in genocat --fastq of a Deep file
                sam_piz_deep_add_qual (vb, STRa(recon)); 
            break;

        default: break;
    }

    COPY_TIMER (sam_piz_con_item_cb);
}

//-----------------------------------------------------------------------------
// SAM PIZ: step 2: move VB deep* buffers to z_file
//-----------------------------------------------------------------------------

// Called by main thread after completion of pizzing of the SAM/BAM txt file: 
// Finalize z_file->deep_index/deep_ents/vb_start_deep_line to be handed to FASTQ to PIZ against the SAM data
static ASCENDING_SORTER (sort_buffer_array_by_vb_i, Buffer, prm32[0])
static void sam_piz_deep_finalize_ents (void)
{
    START_TIMER;

    if (!z_file->deep_index.len) return;

    qsort (STRb(z_file->deep_ents),  sizeof (Buffer), sort_buffer_array_by_vb_i);
    qsort (STRb(z_file->deep_index), sizeof (Buffer), sort_buffer_array_by_vb_i);

    // create vb_start_deep_line to allow easy binary searching on the FASTQ PIZ side
    uint64_t next_deepable_line = 0;
    uint64_t actual_size = 0;

    for (uint32_t vb_idx=0; vb_idx < z_file->deep_index.len32; vb_idx++) {
        BufferP deep_index = B(Buffer, z_file->deep_index, vb_idx);
        BufferP deep_ents  = B(Buffer, z_file->deep_ents , vb_idx); 

        if (flag.show_deep == SHOW_DEEP_ALL) {
            iprintf ("vb_idx=%u vb=%s/%u start_deepable_line=%"PRIu64" num_deepable_lines=%u ents=(%p, %u)\n", 
                     vb_idx, comp_name(deep_ents->prm32[1]), deep_ents->prm32[0], next_deepable_line, deep_index->len32, deep_ents->data, deep_ents->len32);
        }

        actual_size += deep_ents->len;

        BNXT64(z_file->vb_start_deep_line) = next_deepable_line;
        next_deepable_line += deep_index->len;
    }

    // flush (REL) z_file->vb_start_deep_line before setting DEEP_PIZ_ENTS_READY, and prevent hoisting (ACQ) of setting to vb_start_deep_line.param
    __atomic_thread_fence (__ATOMIC_ACQ_REL); 

    if (flag_show_deep) {
        uint64_t total = 0;
        for (DeepStatsPiz i=QNAME_BYTES; i <= QUAL_BYTES; i++) total += z_file->deep_stats[i];

        iprintf ("\nPIZ deep_ents (actual_size=%s) RAM consumption breakdown by field:\n", str_size (actual_size).s);
        iprintf ("Total: %s\n", str_size (total).s); // expected to be the same as actual_size
        
        for (DeepStatsPiz i=QNAME_BYTES; i <= QUAL_BYTES; i++)
            if (z_file->deep_stats[i])
                iprintf ("%-14.14s: %s (%.1f%%)\n", (rom[])DEEP_STATS_NAMES_PIZ[i], str_size (z_file->deep_stats[i]).s, 100.0 * (double)z_file->deep_stats[i] / (double)total);

        if (z_file->deep_stats[SEQ_PACKED_BYTES] || z_file->deep_stats[SEQ_ARITH_BYTES]) {
            iprint0 ("\nPIZ: Reasons for storing packed SEQ in deep_ents rather than a reference:\n");
            for (DeepStatsPiz i=EXPL_DEEPABLE; i < NUM_DEEP_STATS_PIZ; i++)
                if (z_file->deep_stats[i])
                    iprintf ("%-24.24s: %s (%.1f%%)\n", (rom[])DEEP_STATS_NAMES_PIZ[i], str_int_commas (z_file->deep_stats[i]).s, 100.0 * (double)z_file->deep_stats[i] / (double)z_file->deep_stats[EXPL_DEEPABLE]);
        }
    }

    // We need to set param only after deep_ents/index is finalized.
    // sam_deep_piz_wait_for_deep_data will unlock when vb_start_deep_line.param = DEEP_PIZ_ENTS_READY
    store_release (z_file->vb_start_deep_line.param, DEEP_PIZ_ENTS_READY); 

    COPY_TIMER_EVB (sam_piz_deep_finalize_ents);
}

// PIZ: main thread: might be called out of order in case of --test (because dispatcher is allowed out-of-order)
// even if gencomp (bc while dispatcher is in-order, recon_plan is not in order of VBs)
void sam_piz_deep_grab_deep_ents (VBlockSAMP vb)
{
    START_TIMER;

    // two cases: 1. in FASTQ-only reconstruction (eg --R1), DEPN items have need_recon=false, and hence will not recontructed and
    // this function will not be called for them. When unbinding (SAM+FASTQ), it will be called for DEPN VBs, and return here.
    if (IS_DEPN(vb)) return; // DEPN VBs have only SEC/SUPP alignments, and hence don't contribute to Deep data

    // first call in this z_file - allocate
    if (!z_file->deep_index.len) {
        num_deepable_sam_vbs = sections_get_num_vbs (SAM_COMP_MAIN) + sections_get_num_vbs (SAM_COMP_PRIM); 
        
        buf_alloc (evb, &z_file->vb_start_deep_line, 0, num_deepable_sam_vbs, uint64_t, 0, "z_file->vb_start_deep_line");  
        buf_alloc_zero (evb, &z_file->deep_ents,     0, num_deepable_sam_vbs, Buffer,   0, "z_file->deep_ents"); 
        buf_alloc_zero (evb, &z_file->deep_index,    0, num_deepable_sam_vbs, Buffer,   0, "z_file->deep_index");        
    }

    BufferP vb_deep_ents  = &BNXT(Buffer, z_file->deep_ents);  // the next Buffer in an array of Buffers
    BufferP vb_deep_index = &BNXT(Buffer, z_file->deep_index); 

    // grab VB Buffers as-is, and store in a Buffer array (one ents/index Buffer per MAIN/PRIM VB)
    buf_grab (evb, *vb_deep_ents,  "z_file->deep_ents",  vb->deep_ents);
    buf_grab (evb, *vb_deep_index, "z_file->deep_index", vb->deep_index);    

    // set vb_i to allow sorting and display
    vb_deep_ents->prm32[0] = vb_deep_index->prm32[0] = vb->vblock_i;
    vb_deep_ents->prm32[1] = vb_deep_index->prm32[1] = vb->comp_i;

    for (int i=0; i < NUM_DEEP_STATS_PIZ; i++)
        z_file->deep_stats[i] += vb->deep_stats[i];

    if (z_file->deep_index.len == num_deepable_sam_vbs) // this was the last VB
        sam_piz_deep_finalize_ents();

    COPY_TIMER_EVB (sam_piz_deep_grab_deep_ents);
}
