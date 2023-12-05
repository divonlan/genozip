// ------------------------------------------------------------------
//   deep_sam.c
//   Copyright (C) 2023-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#include "bits.h"
#include "piz.h"
#include "threads.h"
#include "huffman.h"
#include "htscodecs/arith_dynamic.h"

// --------
// ZIP side
// --------

static void sam_deep_zip_show_index_stats (void) 
{
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    for (int i=BY_SEQ; i <= BY_QNAME; i++) {
        ARRAY (uint32_t, index, z_file->deep_index_by[i]);

        iprintf ("\nz_file.deep_index_by[%s].len=%"PRIu64" (%s) z_file.deep_ents.len=%"PRIu64" (%s)\n", 
                 by_names[i], 
                 z_file->deep_index_by[i].len, str_size (z_file->deep_index_by[i].len * sizeof (uint32_t)).s, 
                 z_file->deep_ents.len, str_size (z_file->deep_ents.len * sizeof (ZipZDeep)).s);

        #define NUM_LENS 64
        uint64_t used=0;
        uint64_t lens[NUM_LENS] = {}, longer=0; // each entry i contains number linked lists which are of length i; last entry i - linked lists of i or longer
        uint64_t longest_len = 0;
        uint32_t longest_len_hash;

        uint64_t total_ents_traversed = 0; // total ents traversed on linked lists when segging FASTQ reads

        for (uint64_t hash=0; hash < index_len; hash++) {

            uint64_t this_len=0;
            for (uint32_t ent_i = index[hash]; ent_i != NO_NEXT; ent_i = deep_ents[ent_i].next[i]) 
                this_len++;
            
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
            iprintf ("\nFYI: %s index entry %u has %"PRIu64 " deep entries on its linked list. This is excessively large and might cause slowness.%s",
                     by_names[i], longest_len_hash, longest_len, SUPPORT); // note: error only shows with --show-deep... not very useful

        iprintf ("\nDeep index[%s] entries used: %"PRIu64" / %"PRIu64" (%.1f%%)\n", 
                 by_names[i], used, index_len, 100.0 * (double)used / (double)index_len);

        iprintf ("Deep index[%s] linked list length histogram:\n", by_names[i]);
        for (int this_len=0; this_len < NUM_LENS; this_len++) {
            if (!lens[this_len]) continue;

            iprintf ("%3u: %"PRIu64"\n", this_len, lens[this_len]);

            // our algorithm requires traversal of entire linked list for each entry. calculate the length of the list
            // traversed for the average read (simplifying assumption: there are the same )
            total_ents_traversed += (uint64_t)this_len        // number of entries that will access this linked list
                                  * (uint64_t)this_len        // length of linked list traversed for each entry accessing this list
                                  * (uint64_t)lens[this_len]; // number of linked lists of this length
        }

        if (longer) iprintf ("longer: %"PRIu64"\n", longer);

        // this is approximate: assuming # of ents hashed is equal to the number of reads (only for BY_SEQ - we traverse BY_QNAME only for trimmed reads)
        if (i==BY_SEQ)
            iprintf ("\nFASTQ SEG: BY_SEQ index linked list length traversed on average: %.2f (total_ents_traversed=%"PRIu64")\n", 
                     (double)total_ents_traversed / (double)deep_ents_len, total_ents_traversed);
    }
}

static void sam_deep_zip_display_reasons (void)
{
    static rom rsn_names[] = RSN_NAMES;

    uint64_t total = z_file->num_lines;

    iprintf ("%s Alignments breakdown by deepability:\n", dt_name (txt_file->data_type));
    iprintf ("%-13.13s: %"PRIu64"\n", "Total", total);

    for (int i=0; i < NUM_DEEP_STATS_ZIP; i++) 
        if (z_file->deep_stats[i])
            iprintf ("%-13.13s: %"PRIu64" (%.1f%%)\n", rsn_names[i], z_file->deep_stats[i], 100.0 * (double)z_file->deep_stats[i] / (double)total);
}

// Called during zip_finalize of the SAM component for a Deep compression
void sam_deep_zip_finalize (void)
{
    threads_log_by_vb (evb, "main_thread", "sam_deep_zip_finalize", 0);
    
    if (flag.biopsy || flag.biopsy_line.line_i != NO_LINE) return;

    // return unused memory to libc
    buf_trim (z_file->deep_ents, ZipZDeep);
    
    // Build vb_start_deep_line: first deep_line (0-based SAM-wide line_i, but not counting SUPP/SEC lines, monochar reads, SEQ.len=0) of each SAM VB
    buf_alloc_exact_zero (evb, z_file->vb_start_deep_line, z_file->num_vbs + 1, uint64_t, "z_file->vb_start_deep_line");

    uint64_t count=0;
    
    for (VBIType vb_i=1; vb_i <= z_file->num_vbs; vb_i++) {
        Section sec = sections_vb_header (vb_i);
        if (IS_DEPN(sec)) continue;

        *B64(z_file->vb_start_deep_line, vb_i) = count; // note: in ZIP, vb_start_deep_line index is vb_i (unlike PIZ, which is vb_idx)
        count += sec->num_deep_lines;
    }

    if (flag.show_deep) {
        sam_deep_zip_display_reasons();
        sam_deep_zip_show_index_stats();
    }

    memset (z_file->deep_stats, 0, sizeof (z_file->deep_stats)); // we will reuse this field to gather stats while segging the FASTQ files

    // return a bunch of memory to the kernel before moving on to FASTQ
    vb_dehoard_memory (true);
}

// SAM seg: calculate hash values into dl->deep_*

#define MAX_ENTRIES 0xfffffffe // our current implementation is limited to 4G reads

// ZIP compute thread while segging SEQ
void sam_deep_set_QNAME_hash (VBlockSAMP vb, ZipDataLineSAM *dl, QType q, STRp(qname))
{
    // note: we must drop consensus reads, because Deep has an assumption, used to prove uniqueness by hash,
    // that there are no BAM alignments except those in the FASTQ files provided
    if (!dl->FLAG.secondary && !dl->FLAG.supplementary && !segconf.flav_prop[q].is_consensus) {
        dl->deep_hash.qname = deep_qname_hash (q, STRa(qname), NULL);
        dl->is_deepable = true; // note: we haven't tested SEQ yet, so this might still change

        vb->lines.count++;      // counts deepable lines
    }

    else
        vb->deep_stats[segconf.flav_prop[q].is_consensus ? RSN_CONSENSUS 
                     : dl->FLAG.secondary                ? RSN_SECONDARY 
                     :                                     RSN_SUPPLEMENTARY]++; // counting non-deepability reasons
}

// ZIP compute thread while segging SEQ: set hash of forward SEQ (i.e. as it would appear in the FASTQ file) 
void sam_deep_set_SEQ_hash (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(textual_seq))
{
    if (!dl->is_deepable) return; 
    
    // case: sam_deep_set_QNAME_hash happened, but now we discover that SEQ is missing or mono-char
    if (IS_ASTERISK(textual_seq) ||     // note: we don't test SEQ.len, bc SEQ.len=0 for unmapped reads, but they are still deepable
        str_is_monochar (STRa(textual_seq))) { // note: not deepable bc mono-char an elevated chance of hash contention
     
        dl->is_deepable = false;
        vb->lines.count--;

        vb->deep_stats[IS_ASTERISK(textual_seq) ? RSN_NO_SEQ : RSN_MONOCHAR]++; // counting non-deepability reasons        
    }

    else {
        dl->deep_hash.seq = deep_seq_hash (VB, STRa(textual_seq), dl->FLAG.rev_comp);
        vb->deep_stats[RSN_DEEPABLE]++;
    }
}

// ZIP: compute thread while segging QUAL: hash of forward QUAL (i.e. as it would appear in the FASTQ file) for Deep 
void sam_deep_set_QUAL_hash (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qual))
{
    if (!dl->is_deepable || 
        !vb->has_qual || segconf.has_bqsr) return; // case: deepable, but without QUAL: dl->deep_qual_hash remains 0

    dl->deep_hash.qual = deep_qual_hash (VB, STRa(qual), dl->FLAG.rev_comp); 

    if (flag.show_deep == 2 && deephash_issame (dl->deep_hash, flag.debug_deep_hash)) 
        iprintf ("%s Found deep_hash=%u,%u,%u\nQNAME=\"%.*s\"\nSEQ=\"%.*s\"\nQUAL=\"%.*s\"\n",
                 LN_NAME, DEEPHASHf(dl->deep_hash), STRfw(dl->QNAME), dl->SEQ.len, vb->textual_seq_str, STRfw(dl->QUAL));
}

// ZIP compute thread: mutex-protected callback from ctx_merge_in_vb_ctx during merge: add VB's deep_hash to z_file->deep_index_by_*/deep_ents
void sam_deep_merge (VBlockP vb_)
{
    if (!flag.deep || flag.biopsy || flag.biopsy_line.line_i != NO_LINE) return;

    VBlockSAMP vb = (VBlockSAMP)vb_;
    START_TIMER;

    // initialize - first VB merging
    if (!z_file->deep_index_by[BY_SEQ].len) {
        double est_num_vbs = MAX_(1, (double)txtfile_get_seggable_size() / (double)Ltxt * 1.05/* additonal PRIM VBs */);
        uint64_t est_num_lines = MIN_(MAX_ENTRIES, (uint64_t)(est_num_vbs * (double)vb->lines.len * 1.15)); // inflate the estimate by 15% - to reduce hash contention, and to reduce realloc chances for z_file->deep_ents

        // number of bits of hash table size (maximum 31)
        int clzll = (int)__builtin_clzll (est_num_lines * 2); 
        num_hash_bits = MIN_(31, 64 - clzll); // num hash bits

        for (int by=BY_SEQ; by <= BY_QNAME; by++) {
            z_file->deep_index_by[by].can_be_big = true;

            // pointer into deep_ents - initialize to NO_NEXT
            buf_alloc_exact_255 (evb, z_file->deep_index_by[by], ((uint64_t)1 << num_hash_bits), uint32_t, "z_file->deep_index_by"); 
        }

        // note: no initialization of deep_ents - not needed and it takes too long (many seconds to get pages from kernel) while all threads are waiting for vb=1
        z_file->deep_ents.can_be_big = true;
        buf_alloc_exact (evb, z_file->deep_ents, MAX_(est_num_lines, vb->lines.count), ZipZDeep, "z_file->deep_ents");

        if (flag.show_deep) 
            printf ("num_hash_bits=%u est_num_lines*1.15=%"PRIu64" clzll=%d deep_index_by[BY_SEQ].len=%"PRIu64" deep_index_by[BY_QNAME].len=%"PRIu64" deep_ents.len=%"PRIu64"\n", 
                    num_hash_bits, est_num_lines, clzll, z_file->deep_index_by[BY_SEQ].len, z_file->deep_index_by[BY_QNAME].len, z_file->deep_ents.len);

        z_file->deep_ents.len = 0;
    }

    buf_alloc (evb, &z_file->deep_ents, vb->lines.count/*num deepable lines in VB*/, 0, ZipZDeep, CTX_GROWTH, "z_file->deep_ents"); 

    uint32_t ent_i = z_file->deep_ents.len32;

    uint32_t *deep_index_by[2] = { B1ST32 (z_file->deep_index_by[BY_SEQ]), B1ST32 (z_file->deep_index_by[BY_QNAME]) };
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    uint32_t mask = bitmask32 (num_hash_bits); // num_hash_bits is calculated upon first merge
    uint32_t deep_line_i = 0; // line_i within the VB, but counting only deepable lines that have SEQ.len > 0 and are not monochar
    for_buf2 (ZipDataLineSAM, dl, line_i, vb->lines) {
        if (!dl->is_deepable) continue; // case: non-deepable: secondary or supplementary line, no sequence or mono-base

        // we can only handle 4G entries, because hash entry and "next" fields are 32 bit.
        if (ent_i == 0xffffffff) {
            WARN_ONCE ("The number of non-secondary non-supplementary alignments in %s exceeds %u - the maximum supported for Deep compression. The excess alignments will be compressed normally without Deep. Please report this to " EMAIL_SUPPORT, 
                       txt_name, 0xffffffff);
            break;
        }

        // add to indices (by SEQ and by QNAME)
        uint32_t prev_head[2];
        for (int by=BY_SEQ; by <= BY_QNAME; by++) {
            uint32_t hash = (by==BY_SEQ ? dl->deep_hash.seq : dl->deep_hash.qname) & mask;                                 
            prev_head[by] = deep_index_by[by][hash];
            deep_index_by[by][hash] = ent_i; // new head of linked list 
        }
        
        // create new entry - which is now the head of both linked links (BY_SEQ, BY_QNAME)
        deep_ents[ent_i++] = (ZipZDeep){ .next    = { prev_head[0], prev_head[1] }, // previous head, if there is one, will now be the 2nd element
                                         .seq_len = dl->SEQ.len, 
                                         .hash    = dl->deep_hash,
                                         .vb_i    = vb->vblock_i, 
                                         .line_i  = deep_line_i++ };
    }

    z_file->deep_ents.len32 = ent_i;

    for (int i=0; i < NUM_DEEP_STATS_ZIP; i++)
        z_file->deep_stats[i] += vb->deep_stats[i];

    COPY_TIMER(sam_deep_merge);
} 

//-----------------------------------------------------------------------------
// SAM PIZ: step 1: output QNAME,SEQ,QUAL into vb->deep_ents / vb->deep_index
//-----------------------------------------------------------------------------

// main thread: called for after reading global area
void sam_piz_deep_initialize (void)
{
    mutex_initialize (z_file->qname_huf_mutex);

    // In Deep, if we have gencomp, we digest the SAM VBs in writer (during which the serializer is not used), and the 
    // FASTQ VBs in their compute threads (with the serializer). Therefore, the serializer needs to be
    // initialized to skip the SAM VBs 
    if (!z_file->z_flags.adler && // serializer is only used for md5 (considering Deep files only exist since v15)
         z_has_gencomp)           // if no gencomp, we digest all VBs in the compute thread
        z_file->digest_serializer.vb_i_last = sections_get_first_vb_i (SAM_COMP_FQ00) - 1;
}

void sam_piz_deep_init_vb (VBlockSAMP vb, ConstSectionHeaderVbHeaderP header)
{
    vb->arith_compress_bound_longest_seq_len = arith_compress_bound (vb->longest_seq_len, X_NOSZ);
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

#define NUM_QNAMES_TO_SAMPLE 100 // number of QNAMEs to sample for producing the huffman comressor

static void sam_piz_deep_sample_qname (STRp(suffix), // if prfx_len>0 then this is just the suffix the suffix 
                                       int prfx_len)
{
    // note: the mutex is not necessary locked in the order of z_file->qnames_sampled
    mutex_lock (z_file->qname_huf_mutex);

    // first thread that locked the mutex - initialize
    if (!z_file->qname_huf) {
        // this lucky qname is now our master_qname (note: at this point suffix is the whole qname)
        memcpy (z_file->master_qname, suffix, suffix_len); // note: we verified already that qname_len <= SAM_MAX_QNAME_LEN

        // We need to set qname_huf only after master_qname is is set. 
        __atomic_store_n (&z_file->qname_huf, huffman_initialize(), __ATOMIC_RELEASE); 
    }

    uint8_t sample[suffix_len];
    for (int i=0; i < suffix_len; i++)
        sample[i] = suffix[i] ^ z_file->master_qname[prfx_len + i];

    int sampled_so_far = huffman_chew_one_sample (z_file->qname_huf, sample, suffix_len);

    if (sampled_so_far == NUM_QNAMES_TO_SAMPLE) 
        huffman_produce_compressor (z_file->qname_huf);
        
    mutex_unlock (z_file->qname_huf_mutex);
}

// copy QNAME as is to deep_ents
static void sam_piz_deep_add_qname (VBlockSAMP vb)
{
    START_TIMER;

    STRlast (qname, IS_PRIM(vb) ? SAM_QNAMESA : SAM_QNAME);

    sam_piz_alloc_deep_ents (vb, 3 + qname_len); // +3: PizZDeepFlags, qname_len, prfx_len

    BNXT32 (vb->deep_index) = vb->deep_ents.len32; // index in deep_ents of data of current deepable line
    
    PizZDeepFlags *deep_flags = &BNXT(PizZDeepFlags, vb->deep_ents); // also, skip past flags
    bool suffix_is_compressed = false;
    
    if (segconf.deep_qtype == QNONE || flag.seq_only || flag.qual_only)
        goto done;         

    ASSPIZ (qname_len <= SAM_MAX_QNAME_LEN, "QNAME.len=%u, but per BAM specfication it is limited to %u characters. QNAME=\"%.*s\"", SAM_MAX_QNAME_LEN, qname_len, STRf(qname));

    qname_canonize (QNAME1, qSTRa(qname)); // might reduce qname_len

    // get common prefix length of qname and segconf.qame1
    int prfx_len=0; 
    if (z_file->qname_huf) // case: master_qname is initialized
        for (; prfx_len < qname_len; prfx_len++)
            if (qname[prfx_len] != z_file->master_qname[prfx_len]) break;
    
    rom suffix = qname + prfx_len;
    int suffix_len = qname_len - prfx_len;

    if (z_file->qnames_sampled < NUM_QNAMES_TO_SAMPLE && // first test without atomic - quick fail without the atomicity overhead
        __atomic_fetch_add (&z_file->qnames_sampled, (int)1, __ATOMIC_RELAXED) < NUM_QNAMES_TO_SAMPLE)
        sam_piz_deep_sample_qname (STRa(suffix), prfx_len);

    BNXT8(vb->deep_ents) = qname_len;
    BNXT8(vb->deep_ents) = prfx_len;

    // if the huffman compressor is ready, we compress the suffix
    if (huffman_is_produced (z_file->qname_huf)) { 
        // compress suffix XOR master_qname
        uint8_t uncomp[suffix_len];
        for (int i=0; i < suffix_len; i++)
            uncomp[i] = suffix[i] ^ z_file->master_qname[prfx_len + i];

        uint32_t comp_len = huffman_comp_len_required_allocation(suffix_len);
        uint8_t comp[comp_len];
        memset (comp, 0, comp_len);
        
        huffman_compress (z_file->qname_huf, uncomp, suffix_len, qSTRa(comp));

        if (comp_len >= suffix_len) goto uncompressed_suffix;

        suffix_is_compressed = true;
        buf_add (&vb->deep_ents, comp, comp_len);

        vb->deep_stats[QNAME_BYTES] += 2 + comp_len;
    }

    else uncompressed_suffix: {
        buf_add (&vb->deep_ents, suffix, suffix_len);

        vb->deep_stats[QNAME_BYTES] += 2 + suffix_len;
    }

done:
    // update flags
    *deep_flags = (PizZDeepFlags){ .is_qname_comp = suffix_is_compressed }; 

    COPY_TIMER (sam_piz_deep_add_qname);
}

// compress QUAL or SEQ - adding a 1B or 4B length
static bool sam_piz_deep_compress (VBlockSAMP vb, STRp(data), bool is_seq)
{
    START_TIMER;

    // reverse if needed
    if (last_flags.rev_comp) {
        ASSERTNOTINUSE (vb->scratch);
        buf_alloc_exact (vb, vb->scratch, data_len, char, "scratch");

        if (is_seq) str_revcomp (B1STc(vb->scratch), STRa(data));
        else        str_reverse (B1STc(vb->scratch), STRa(data));
    
        data = B1STc(vb->scratch);
    }

    sam_piz_alloc_deep_ents (vb, 4 + vb->arith_compress_bound_longest_seq_len);

    uint8_t *next = BAFT8 (vb->deep_ents);

    uint32_t comp_len = vb->arith_compress_bound_longest_seq_len;
    
    bool len_is_1_byte = (data_len < (is_seq ? 1024 : 512)); // just a guess for now
    uint8_t *guess_addr = next + (len_is_1_byte ?  1 : 4);   // we don't know comp_len yet, so just a guess
    ASSPIZ (arith_compress_to (evb, (uint8_t *)STRa(data), guess_addr, &comp_len, X_NOSZ) && comp_len,
            "%s: Failed to compress data: data_len=%u data=\"%.*s\"", LN_NAME, data_len, STRf(data));

    // we guessed wrong that we need 1 byte for the length - actually need 4 
    if (len_is_1_byte && comp_len > 255) {
        memmove (guess_addr+3, guess_addr, comp_len); // make room for another 3 bytes
        len_is_1_byte = false;
    }

    if (len_is_1_byte)  // 1B length
        *next = comp_len;

    else  // 4B length
        memcpy (next, &comp_len, 4); // memcpy bc non-aligned

    uint32_t n_bytes = comp_len + (len_is_1_byte ? 1 : 4);
    vb->deep_ents.len32 += n_bytes;
    vb->deep_stats[is_seq ? SEQ_BYTES : QUAL_BYTES] += n_bytes;

    if (last_flags.rev_comp)   
        buf_free (vb->scratch);

    COPY_TIMER (sam_piz_deep_compress);

    return !len_is_1_byte; // true long length
}

// container item callback for SQBITMAP: pack SEQ into 2-bit and into deep_ents or compress if it has non-ACGT (eg N)
static void sam_piz_deep_add_seq (VBlockSAMP vb, STRp(seq))
{
    START_TIMER;

    if (str_is_monochar (STRa(seq))) 
        vb->seq_is_monochar = true; 

    // case: not deepable bc no SEQ or monochar - rollback
    if (vb->seq_missing || vb->seq_is_monochar) {
        vb->deep_ents.len32 = *BLST32 (vb->deep_index);
        vb->deep_index.len32--;
        return;
    }

    PizZDeepFlags *deep_flags = B(PizZDeepFlags, vb->deep_ents, *BLST32(vb->deep_index));

    if (flag.header_only_fast) return;

    if (seq_len <= 255) 
        BNXT8 (vb->deep_ents) = seq_len;
    
    else {
        deep_flags->is_long_seq = true;
        memcpy (BAFT8 (vb->deep_ents), &seq_len, sizeof (uint32_t)); // not word-aligned
        vb->deep_ents.len32 += sizeof (uint32_t);
    }

    if (flag.qual_only) return;

    // case: perfect match with perfect CIGAR (verbatim copy of the reference) - we need only GPOS
    if (vb->deep_seq.gpos) { // note possible edge case: gpos exists and is actually zero. we will fail this condition and fallback to the normal treatment of this alignment.
        bool perfect = !vb->deep_seq.mismatch_base;
        memcpy (BAFT8 (vb->deep_ents), &vb->deep_seq, (perfect ? 4 : 6));
        vb->deep_ents.len32 += (perfect ? 4 : 6);

        bool used_aligner = ctx_has_value_in_line_(vb, CTX(SAM_STRAND));
        bool aligner_rev = used_aligner && !CTX(SAM_STRAND)->last_value.i;
        bool rev_comp = last_flags.rev_comp ^ aligner_rev;
        
        deep_flags->seq_encoding = perfect ? (rev_comp ? ZDEEP_SEQ_PERFECT_REV : ZDEEP_SEQ_PERFECT)
                                           : (rev_comp ? ZDEEP_SEQ_MIS_1_REV   : ZDEEP_SEQ_MIS_1);
    }

    // if only A,C,G,T - pack in 2-bits
    else if (str_is_ACGT (STRa(seq), NULL)) {
        sam_piz_alloc_deep_ents (vb, 4 + seq_len + 2*sizeof(uint64_t)/*padding for Bits*/); 

        uint32_t pack_index = ROUNDUP8 (vb->deep_ents.len32); // 64b word-align as required by Bits
        Bits pack = bits_init (seq_len * 2, B8(vb->deep_ents, pack_index), vb->deep_ents.size - pack_index, false);
        sam_seq_pack (vb, &pack, 0, STRa(seq), false, last_flags.rev_comp, HARD_FAIL);

        uint32_t n_bytes = roundup_bits2bytes (pack.nbits);

        if (BAFT8 (vb->deep_ents) != (uint8_t*)pack.words)
            memmove (BAFT8 (vb->deep_ents), pack.words, n_bytes); // move back to (unaligned) place

        vb->deep_ents.len32 += n_bytes;
        vb->deep_stats[SEQ_BYTES] += n_bytes;
        deep_flags->seq_encoding = ZDEEP_SEQ_PACKED;
    }

    // if has an N or other IUPACs - compress
    else {
        bool is_long = sam_piz_deep_compress (vb, STRa(seq), true);
        
        deep_flags = B(PizZDeepFlags, vb->deep_ents, *BLST32(vb->deep_index)); // cafeful! sam_piz_deep_compress reallocs vb->deep_ents
        deep_flags->seq_encoding = (is_long ? ZDEEP_SEQ_COMPRESSED_LONG_LEN : ZDEEP_SEQ_COMPRESSED);
    }

    COPY_TIMER (sam_piz_deep_add_seq);
}

// container item callback for QUAL/QUALSA: compress QUAL into deep_ents
static void sam_piz_deep_add_qual (VBlockSAMP vb, STRp(qual))
{
    START_TIMER;

    if (vb->seq_missing || vb->seq_is_monochar || VB_SAM->qual_missing || // not a deepable alignment or no QUAL
        flag.seq_only || flag.header_only_fast) return; // FASTQ is not going to reconstruct QUAL

    uint32_t index = *BLST32(vb->deep_index);
    ASSERT (index < vb->deep_ents.len32, "index=%u out of range: vb->deep_ents.len=%u", index, vb->deep_ents.len32);

    bool is_long = sam_piz_deep_compress (vb, STRa(qual), false);

    PizZDeepFlags *deep_flags = B(PizZDeepFlags, vb->deep_ents, *BLST32(vb->deep_index)); // cafeful! sam_piz_deep_compress reallocs vb->deep_ents
    deep_flags->is_long_qual_comp = is_long;

    COPY_TIMER (sam_piz_deep_add_qual);
}


// For SAM: called for top-level SAM container if compressed with --deep (defined as con_item_cb in data_types.h)
// for BAM: called from SEQ and QUAL translators (before translation) 
CONTAINER_ITEM_CALLBACK (sam_piz_con_item_cb)
{
    START_TIMER;

    if (con_item->dict_id.num == _SAM_SQBITMAP) {
        if (!flag.deep                                       || // Deep file (otherwise this callback will not be set), but SAM/BAM-only reconstruction
            last_flags.secondary || last_flags.supplementary || // SECONDARY or SUPPLEMENTARY alignment
            (segconf.flav_prop[QNAME2].is_consensus && ctx_encountered_in_line (vb, SAM_QNAME2)) || // consensus alignment 
            IS_ASTERISK(recon)                               || // Missing SEQ
            str_is_monochar (STRa(recon))) {                    // mono-char SEQ

            VB_SAM->line_not_deepable = true;
        }
        else {
            // note: we do QNAME during the SQBITMAP callback, because we need last_flags to be set 
            // (in SAM QNAME reconstructs before flags). On the other hand, SQBITMAP and QUAL must be 
            // handled right after reconstruction, because in BAM, immediately after they get converted to binary format.
            sam_piz_deep_add_qname (VB_SAM);
            sam_piz_deep_add_seq (VB_SAM, STRa(recon));
        }
    }

    else if (!VB_SAM->line_not_deepable && 
             (con_item->dict_id.num == _SAM_QUAL || con_item->dict_id.num == _SAM_QUALSA))
        sam_piz_deep_add_qual (VB_SAM, STRa(recon)); 

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
   
    for (uint32_t vb_idx=0; vb_idx < z_file->deep_index.len32; vb_idx++) {
        BufferP deep_index = B(Buffer, z_file->deep_index, vb_idx);

        // too much output, uncomment when needed
        // BufferP deep_ents = B(Buffer, z_file->deep_ents, vb_idx); 
        // if (flag.show_deep)
        //     iprintf ("vb_idx=%u vb=%s/%u start_deepable_line=%"PRIu64" num_deepable_lines=%u ents=(%p, %u)\n", 
        //              vb_idx, comp_name(deep_ents->prm32[1]), deep_ents->prm32[0], next_deepable_line, deep_index->len32, deep_ents->data, deep_ents->len32);

        BNXT64(z_file->vb_start_deep_line) = next_deepable_line;
        next_deepable_line += deep_index->len;
    }

    if (flag.show_deep) {
        static rom names[] = DEEP_ENT_NAMES;

        uint64_t total = z_file->deep_stats[QNAME_BYTES] + z_file->deep_stats[SEQ_BYTES] + z_file->deep_stats[QUAL_BYTES];
        iprint0 ("\ndeep_ents RAM consumption breakdown by field:\n");
        iprintf ("Total: %s\n", str_size (total).s);
        
        for (int i=0; i < NUM_DEEP_STATS_PIZ; i++)
            iprintf ("%-5.5s: %s (%.1f%%)\n", names[i], str_size (z_file->deep_stats[i]).s, 100.0 * (double)z_file->deep_stats[i] / (double)total);
    }

    // We need to set param only after deep_ents/index is finalized.
    // sam_deep_piz_wait_for_deep_data will unlock when vb_start_deep_line.param = true
    __atomic_store_n (&z_file->vb_start_deep_line.param, (uint64_t)true, __ATOMIC_RELEASE); 

    COPY_TIMER_EVB (sam_piz_deep_finalize_ents);
}

// PIZ: main thread: might be called out of order in case of --test (because dispatcher is allowed out-of-order)
// even if gencomp (bc while dispatcher is in-order, recon_plan is not in order of VBs)
void sam_piz_deep_grab_deep_ents (VBlockSAMP vb)
{
    START_TIMER;

    #define num_deepable_sam_vbs z_file->deep_index.param

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
