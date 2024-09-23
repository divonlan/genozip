// ------------------------------------------------------------------
//   sam_ingest_grps.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#include "writer.h"
#include "qname.h"
#include "compressor.h"
#include "huffman.h"
#include "htscodecs/rANS_static4x16.h"

//------------------------------------
// Ingesting SA Group from PRIM VBs
//------------------------------------

static Mutex seq_mutex={}, qual_mutex={}, qname_mutex={}, aln_mutex={}, grp_mutex={};

typedef struct { uint32_t qname_hash, grp_i; } SAGroupIndexEntry; 

void sam_sa_prim_initialize_ingest (void)
{
    mutex_initialize (seq_mutex);
    mutex_initialize (qual_mutex);
    mutex_initialize (qname_mutex);
    mutex_initialize (aln_mutex);
    mutex_initialize (grp_mutex);

    // suppress buf_alloc warning of big allocations
    z_file->sag_grps.can_be_big = z_file->sag_seq.can_be_big = z_file->sag_qual.can_be_big = z_file->sag_qnames.can_be_big = true;
}

static ASCENDING_SORTER (group_index_sorter, SAGroupIndexEntry, qname_hash)

// ZIP: called from main thread by sam_zip_after_compute after final PRIM vb
void sam_sa_prim_finalize_ingest (void)
{
    START_TIMER;
    if (!z_file->sag_grps.len) return; 

    ASSERT (z_file->sag_grps.len <= 0xffffffffULL, "sag_grps.len=%"PRIu64" exceeds 32 bits - need to widen SAGroupIndexEntry.grp_i", z_file->sag_grps.len);

    // build index by hash(qname)
    qsort (z_file->sag_grps_index.data, z_file->sag_grps_index.len, sizeof(SAGroupIndexEntry), group_index_sorter);

    if (flag.show_sag) sam_show_sag();

    buf_trim (z_file->sag_qnames, char);
    buf_trim (z_file->solo_data, char);
    
    if      (IS_SAG_SA)   buf_trim (z_file->sag_alns, SAAln);
    else if (IS_SAG_CC)   buf_trim (z_file->sag_alns, CCAln);
    else if (IS_SAG_SOLO) buf_trim (z_file->sag_alns, SoloAln);

    // TODO: trim z_file->sag_seq (bits buffer)
    
    COPY_TIMER_EVB (sam_sa_prim_finalize_ingest);
}

// pack seq (into 2bit ACGT format) of each PRIM line separately 
static void sam_zip_prim_ingest_vb_pack_seq (VBlockSAMP vb, Sag *vb_grps, uint32_t vb_grps_len,
                                             BufferP packed_seq_buf, bool is_bam_format)
{
    START_TIMER;

    uint32_t total_seq_len=0;
    for (uint32_t grp_i=0; grp_i < vb_grps_len; grp_i++) 
        total_seq_len += vb_grps[grp_i].seq_len;

    // allocate memory but don't extend the bitmap yet
    ASSERTNOTINUSE (*packed_seq_buf);
    BitsP sag_seq = buf_alloc_bits (vb, packed_seq_buf, total_seq_len * 2, NOINIT, 0, "packed_seq_buf");
    uint64_t next_bit = 0;

    for (uint32_t vb_grp_i=0; vb_grp_i < vb_grps_len; vb_grp_i++) {    
        Sag *vb_grp = &vb_grps[vb_grp_i];

        sam_seq_pack (vb, sag_seq, next_bit, Btxt(vb_grp->seq), vb_grp->seq_len, is_bam_format, false, HARD_FAIL);
        vb_grp->seq = next_bit / 2; // update from an index into txt_data to an index (bases not bits) into sag_seq
    
        next_bit += vb_grp->seq_len * 2;
    }

    bits_clear_excess_bits_in_top_word ((BitsP)packed_seq_buf, true);

    COPY_TIMER (sam_zip_prim_ingest_vb_pack_seq);
}

// separately compress qual of each PRIM line 
static void sam_zip_prim_ingest_vb_compress_qual (VBlockSAMP vb, Sag *vb_grps, uint32_t vb_grps_len, BufferP comp_qual_buf)
{
    START_TIMER;

    // count total qual scores in this VB
    uint32_t total_qual_len=0;
    for (uint32_t grp_i=0; grp_i < vb_grps_len; grp_i++) 
        total_qual_len += (vb_grps[grp_i].no_qual ? 0 : vb_grps[grp_i].seq_len);

    if (!total_qual_len) return;

    // note: in an unlikely case, the compressed size might be beyond this - in which case we will abandon the compression.
    ASSERTNOTINUSE (*comp_qual_buf);
    uint32_t est_total_comp_qual_len = MIN_(total_qual_len / 1.3 + rans_compress_bound_4x16 (segconf.longest_seq_len, X_NOSZ), MAX_SA_SEQ_LEN); // heuristic (we don't want this to be unnecessarily too big)
    buf_alloc (vb, comp_qual_buf, 0, est_total_comp_qual_len, char, 0, "comp_qual_buf"); // likely already allocated

    for (uint32_t vb_grp_i=0; vb_grp_i < vb_grps_len; vb_grp_i++) {
        Sag *vb_grp = &vb_grps[vb_grp_i];
        if (vb_grp->no_qual) continue;

        uint32_t comp_len = rans_compress_bound_4x16 (vb_grp->seq_len, X_NOSZ);
        buf_alloc (vb, comp_qual_buf, comp_len, 0, char, CTX_GROWTH, "comp_qual_buf"); // likely already allocated
        
        uint8_t *qual = B8 (vb->txt_data, vb_grp->qual);  // always in SAM format
        ASSERT (rans_compress_to_4x16 (VB, qual, vb_grp->seq_len, BAFT8(*comp_qual_buf), &comp_len, X_NOSZ) && comp_len,
                "%s: Failed to compress PRIM qual: qual_len=%u", LN_NAME, vb_grp->seq_len);

        // TO DO: this is not expected to occur as the limit is quite high. If this does occur, we can avoid failing by
        // marking this vb_grp_i for removal and removing it prior to appending to z_file.
        ASSERT (comp_len <= MAX_SA_QUAL_COMP_LEN, "%s: while ingesting vb_grp_i=%u: qual_comp_len=%u > MAX_SA_QUAL_COMP_LEN=%u. seq_len=%u. %s",
                VB_NAME, vb_grp_i, comp_len, MAX_SA_QUAL_COMP_LEN, vb_grp->seq_len, report_support());

        vb_grp->qual          = comp_qual_buf->len; // if compressed - index in comp_qual_buf, if not - remains index in txt_data
        vb_grp->qual_comp_len = comp_len; 
        comp_qual_buf->len   += comp_len;    
    }

    vb->comp_qual_len = comp_qual_buf->len; // goes into SectionHeaderVbHeader.sam_prim_comp_qual_len 

    COPY_TIMER (sam_zip_prim_ingest_vb_compress_qual);
}

static void sam_zip_prim_ingest_vb_compress_qnames (VBlockSAMP vb, Sag *vb_grps, uint32_t vb_grps_len, BufferP comp_qname_buf, BufferP qname_hash_buf)
{
    START_TIMER;

    ASSERTNOTINUSE (*comp_qname_buf);
    ASSERTNOTINUSE (*qname_hash_buf);

    ARRAY_alloc (uint32_t, qname_hashes, vb_grps_len, false, *qname_hash_buf, vb, "qname_hash_buf");
    
    // count total qname length in this VB
    uint32_t total_qnames_len=0;
    for (uint32_t grp_i=0; grp_i < vb_grps_len; grp_i++) 
        total_qnames_len += vb_grps[grp_i].qname_len;
    
    buf_alloc (vb, comp_qname_buf, 0, total_qnames_len / 2, uint8_t, 0, "comp_qname_buf");

    for (uint32_t vb_grp_i=0; vb_grp_i < vb_grps_len; vb_grp_i++) {
        Sag *g = &vb_grps[vb_grp_i];
        ZipDataLineSAMP dl = DATA_LINE(vb_grp_i);

        rom qname = Btxt (g->qname);
        uint32_t qname_len = g->qname_len;
        
        // case: segged as COPY_BUDDY (for mate) from mate, with flavor that has the same exactly qname for both mates (i.e. no 1/2)
        if (dl->qname_mate_copied_exactly) {
            g->qname = vb_grps[dl->mate_line_i].qname;
         
            qname_hashes[vb_grp_i] = qname_hash_change_last (qname_hashes[dl->mate_line_i], g->is_last);
        }    
        
        else {
            g->qname = comp_qname_buf->len; // update from index into vb->txt_data to index into z_file->sa_qname

            uint32_t comp_len = huffman_comp_len_required_allocation (qname_len);
            buf_alloc (vb, comp_qname_buf, comp_len, 0, uint8_t, 0, NULL);
            
            huffman_compress (SAM_QNAME, STRa(qname), BAFT8(*comp_qname_buf), &comp_len);
            comp_qname_buf->len += comp_len;

            qname_hashes[vb_grp_i] = qname_calc_hash (QNAME1, COMP_NONE, STRa(qname), g->is_last, false, CRC32, NULL);
        }
    }

    CTX(SAM_QNAME)->huffman.comp_len = comp_qname_buf->len32; // goes into VbHeader.sam_prim_comp_qname_len

    COPY_TIMER (sam_zip_prim_ingest_vb_compress_qnames);
}

static void sam_zip_prim_ingest_vb_create_index (VBlockSAMP vb, BufferP index_buf, ConstBufferP qname_hash_buf)
{
    START_TIMER;

    buf_alloc_exact (vb, *index_buf, qname_hash_buf->len, SAGroupIndexEntry, "sag_index");
    SAGroupIndexEntry *first = B1ST (SAGroupIndexEntry, *index_buf);

    for_buf_tandem (SAGroupIndexEntry, index, *index_buf, uint32_t, hash, *qname_hash_buf)
        *index = (SAGroupIndexEntry){ .grp_i = index - first, .qname_hash = *hash };

    qsort (first, index_buf->len32, sizeof(SAGroupIndexEntry), group_index_sorter);

    COPY_TIMER (sam_zip_prim_ingest_vb_create_index);
}

static inline ZWord copy_solo (VBlockSAMP vb, char **next, TxtWord txtw, bool check_copy)
{
    uint64_t start = BNUM64(z_file->solo_data, *next);

    // if identical, UR and UB, CR and CB, point to the same data in solo_data
    if (check_copy && start >= txtw.len && !memcmp (*next - txtw.len, Btxt (txtw.index), txtw.len)) 
        return (ZWord){ .index = start - txtw.len, .len = txtw.len };
    
    else {
        *next = mempcpy (*next, Btxt (txtw.index), txtw.len); 
        return (ZWord){ .index = start, .len = txtw.len };
    }
}

// runs in compute thread, serialized by mutex
void sam_zip_prim_ingest_solo_data (VBlockSAMP vb)
{
    START_TIMER;

    // allocate memory for solo_data
    ASSERT (vb->sag_grps.len32 == vb->lines.len32, "%s: Expecting sag_grps.len=%u == lines.len=%u", 
            VB_NAME, vb->sag_grps.len32, vb->lines.len32);

    uint32_t total_solo_len=0;    
    for_buf (ZipDataLineSAM, dl, vb->lines) 
        for (int tag_i = 0; tag_i < NUM_SOLO_TAGS; tag_i++)
            total_solo_len += dl->solo_z_fields[tag_i].len;

    z_file->solo_data.can_be_big = z_file->sag_alns.can_be_big = true; // suppress warnings
    buf_alloc (evb, &z_file->solo_data, total_solo_len, 0, char, CTX_GROWTH, "solo_data");
    buf_alloc (evb, &z_file->sag_alns, vb->sag_grps.len32, 0, SoloAln, CTX_GROWTH, "solo_aln");

    char *next = BAFTc (z_file->solo_data);
    for_buf (ZipDataLineSAM, dl, vb->lines) {
        SoloAln solo_aln;
        for (int tag_i = 0; tag_i < NUM_SOLO_TAGS; tag_i++)
            solo_aln.word[tag_i] = copy_solo (vb, &next, dl->solo_z_fields[tag_i], solo_props[tag_i].maybe_same_as_prev);
        
        BNXT(SoloAln, z_file->sag_alns) = solo_aln;
    }

    vb->solo_data_len = BNUM64 (z_file->solo_data, next) - z_file->solo_data.len;
    z_file->solo_data.len = BNUM64 (z_file->solo_data, next);

    COPY_TIMER (sam_zip_prim_ingest_solo_data);
}

// ZIP PRIM VB: called by compute thread, from sam_seg_finalize - VBs are out of order
// Called after segging a PRIM VB, to ingest its data into vb->sa_*, for later use by DEPN VBs. 
void sam_zip_prim_ingest_vb (VBlockSAMP vb)
{
    START_TIMER;

    ARRAY (Sag, vb_grps, vb->sag_grps);
    if (!vb_grps_len) return; // no SA groups

    ASSERTNOTINUSE (vb->codec_bufs[1]); // note: not [0], bc it is used by codec_hts_compress which we use to compress reread prescriptions
    ASSERTNOTINUSE (vb->codec_bufs[2]);
    ASSERTNOTINUSE (vb->codec_bufs[3]);
    #define comp_qname_buf vb->codec_bufs[1]
    #define qname_hash_buf vb->codec_bufs[2]
    #define comp_qual_buf  vb->codec_bufs[3]
    buf_free (vb->z_data);
    #define packed_seq_buf vb->z_data // recycle z_data memory instead of allocating more
    #define sag_grps_index_buf comp_qual_buf
    
    // compression happens outside of any mutex
    sam_zip_prim_ingest_vb_compress_qnames (vb, STRa(vb_grps), &comp_qname_buf, &qname_hash_buf);
    
    sam_zip_prim_ingest_vb_pack_seq (vb, STRa(vb_grps), &packed_seq_buf, IS_BAM_ZIP);    
    
    sam_zip_prim_ingest_vb_compress_qual (vb, STRa(vb_grps), &comp_qual_buf);

    // Note: Buffers are pre-INITed in file_initialize_bufs, so we can buf_alloc them from serialized compute threads
    // Note: We must serialize the data populating too, bc otherwise other VBs can re-alloc the buffers under our feet
    // VBs can enter in arbitrary order, and also handle SEQ, QUAL, QUAL, ALN in arbitrary order by availability of the mutex
    
    bool seq_done=false, qual_done=false, qname_done=false, aln_done=false;

    while (!seq_done || !qual_done || !qname_done || !aln_done) {

        bool achieved_something = false;
        
        // concatenate this VB's sequence to z_file->sag_seq (in ACGT format)
        if (!seq_done && mutex_trylock (seq_mutex)) {
            uint64_t start_seq = z_file->sag_seq.nbits / 2; // next_seq in bases, not bits
            
            ASSERT ((z_file->sag_seq.nbits + packed_seq_buf.nbits) / 2 <= MAX_SA_SEQ_INDEX, "%s: while ingesting: Total seq length of all prim sequences exceeds maximum of %"PRIu64" bases. Workaround: re-run with --no-gencomp. %s",
                    VB_NAME, MAX_SA_SEQ_INDEX, report_support());

            buf_alloc_bits (evb, &z_file->sag_seq, z_file->sag_seq.nbits + packed_seq_buf.nbits, NOINIT, CTX_GROWTH, "z_file->sag_seq");
            bits_copy ((BitsP)&z_file->sag_seq, start_seq * 2, (BitsP)&packed_seq_buf, 0, packed_seq_buf.nbits);

            seq_done = achieved_something = true;                
            mutex_unlock (seq_mutex);

            for (uint32_t vb_grp_i=0; vb_grp_i < vb_grps_len; vb_grp_i++)
                vb_grps[vb_grp_i].seq += start_seq; 
        } 

        // concatenate this VB's compressed qual to z_file->sag_qual
        if (!qual_done && mutex_trylock (qual_mutex)) {
            uint64_t start_qual = z_file->sag_qual.len;

            ASSERT (z_file->sag_qual.len + comp_qual_buf.len <= MAX_SA_QUAL_INDEX, "%s: while ingesting: Total qual length of all prim sequences exceeds maximum of %"PRIu64". Workaround: re-run with --no-gencomp. %s",
                    VB_NAME, MAX_SA_QUAL_INDEX, report_support());

            buf_append (evb, z_file->sag_qual, uint8_t, comp_qual_buf.data, comp_qual_buf.len, NULL);
            qual_done = achieved_something = true;
            mutex_unlock (qual_mutex);

            for (uint32_t vb_grp_i=0; vb_grp_i < vb_grps_len; vb_grp_i++) 
                vb_grps[vb_grp_i].qual += start_qual;
        }

        // concatenate this VB's compressed qnames to z_file->sag_qnames and extend the (to-be-sorted) index in z_file->sag_grps_index
        if (!qname_done && mutex_trylock (qname_mutex)) {
            uint64_t start_qname = z_file->sag_qnames.len;

            ASSERT (z_file->sag_qnames.len + comp_qname_buf.len <= MAX_SA_QNAME_INDEX, "%s: while ingesting: Total qname length of all prim sequences exceeds maximum of %"PRIu64". Workaround: re-run with --no-gencomp. %s",
                    VB_NAME, MAX_SA_QNAME_INDEX, report_support());

            buf_append (evb, z_file->sag_qnames, uint8_t, comp_qname_buf.data, comp_qname_buf.len, NULL);
            qname_done = achieved_something = true;
            mutex_unlock (qname_mutex);

            for (uint32_t vb_grp_i=0; vb_grp_i < vb_grps_len; vb_grp_i++) 
                vb_grps[vb_grp_i].qname += start_qname;
        }

        // concatenate this VB's alignment data (depends on SAG type) to z_file->sag_alns
        if (!aln_done && mutex_trylock (aln_mutex)) {
            uint64_t start_alns = 0;
            if (IS_SAG_SA) {
                start_alns = z_file->sag_alns.len;
                buf_append_buf (evb, &z_file->sag_alns, &vb->sag_alns, SAAln, NULL); // copy alignments (also allocs memory)
            
                z_file->sag_alns.count = z_file->sag_alns.len;
            }

            else if (IS_SAG_CC) 
                buf_append_buf (evb, &z_file->sag_alns, &vb->sag_alns, CCAln, NULL); 

            else if (IS_SAG_SOLO) 
                sam_zip_prim_ingest_solo_data (vb);

            if (!IS_SAG_SA)
                for (uint32_t grp_i=0; grp_i < vb_grps_len; grp_i++)
                    z_file->sag_alns.count += vb_grps[grp_i].num_alns; // in Sag types other than SA, we don't store alignments (just the prim) - just count them (for show_sag)
             
            aln_done = achieved_something = true;
            mutex_unlock (aln_mutex);

            // update index from index into vb->sag_alns to z_file->sag_alns
            if (IS_SAG_SA) {
                ASSERT (vb_grps[vb_grps_len-1].first_aln_i + start_alns <= MAX_SA_GRP_ALNS, 
                        "%s: while ingesting: Total number of prim+depn alignments in this file exceeds maximum of %"PRIu64". (vb_grps_len=%"PRIu64" first_aln_i=%"PRIu64" start_alns=%"PRIu64") Workaround: re-run with --no-gencomp. %s",
                        VB_NAME, MAX_SA_GRP_ALNS, vb_grps_len, (uint64_t)vb_grps[vb_grps_len-1].first_aln_i, start_alns, report_support());

                for (uint32_t grp_i=0; grp_i < vb_grps_len; grp_i++)
                    vb_grps[grp_i].first_aln_i += start_alns; 
            }
        }

        if (!achieved_something) {
            START_TIMER;
            usleep (1000); // 1 ms
            COPY_TIMER_EVB (sam_zip_prim_ingest_idle);
            if (!seq_done)   COPY_TIMER_EVB (sam_zip_prim_ingest_wait_for_seq_mutex);
            if (!qual_done)  COPY_TIMER_EVB (sam_zip_prim_ingest_wait_for_qual_mutex);
            if (!qname_done) COPY_TIMER_EVB (sam_zip_prim_ingest_wait_for_qname_mutex);
            if (!aln_done)   COPY_TIMER_EVB (sam_zip_prim_ingest_wait_for_aln_mutex);
        }
    }

    vb_grps[0].first_grp_in_vb = true; // mark the first group in this PRIM VB

    // calculate index outside of mutex (time consuming due to qsort)
    sam_zip_prim_ingest_vb_create_index (vb, &sag_grps_index_buf, &qname_hash_buf);

    // Groups and index (by qname_hash) into groups - must be last, after vb_grps is updated ^. Order of VBs is arbitrary
    mutex_lock (grp_mutex);
    vb->first_grp_i = z_file->sag_grps.len; // the index of first group of this PRIM VB, in z_file->sag_grps. note: might be out-of-order of VBs.
    buf_append_buf (evb, &z_file->sag_grps, &vb->sag_grps, Sag, NULL); 

    for_buf (SAGroupIndexEntry, index, sag_grps_index_buf)
        index->grp_i += z_file->sag_grps_index.len32;

    buf_append_buf (evb, &z_file->sag_grps_index, &sag_grps_index_buf, SAGroupIndexEntry, NULL); 

    mutex_unlock (grp_mutex);    

    buf_destroy (comp_qname_buf);
    buf_destroy (qname_hash_buf);
    buf_destroy (comp_qual_buf);
    buf_free (vb->z_data);

    COPY_TIMER_EVB (sam_zip_prim_ingest_vb); 
}
