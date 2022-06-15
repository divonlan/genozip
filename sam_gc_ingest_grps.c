// ------------------------------------------------------------------
//   sam_ingest_grps.c
//   Copyright (C) 2022-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "chrom.h"
#include "context.h"
#include "profiler.h"
#include "codec.h"
#include "bits.h"
#include "writer.h"
#include "sections.h"
#include "codec.h"
#include "compressor.h"
#include "libdeflate/libdeflate.h"
#include "htscodecs/rANS_static4x16.h"

//------------------------------------
// Ingesting SA Group from PRIM VBs
//------------------------------------

static Mutex seq_mutex={}, qual_mutex={}, qname_mutex={}, aln_mutex={}, grp_mutex={};

typedef struct __attribute__ ((__packed__)) { uint32_t qname_hash, grp_i; } SAGroupIndexEntry; 

static ASCENDING_SORTER (group_index_sorter, SAGroupIndexEntry, qname_hash)

static void sam_sa_create_group_index (void)
{
    ARRAY (SAGroupType, grp, z_file->sa_groups);

    if (!grp_len) return;

    // Note: if SAGroup is 32 bit, verify grp_len 
    ASSERT (sizeof(SAGroup)==8 || grp_len <= 0xffffffffULL, "sa_groups.len=%"PRIu64" exceeds 32 bits - need to extend index to be 64 bits", z_file->sa_groups.len);
    
    buf_alloc (evb, &z_file->sa_groups_index, grp_len, 0, sizeof (SAGroupIndexEntry), 1, "z_file->sa_groups_index");
    z_file->sa_groups_index.len = grp_len;

    ARRAY (SAGroupIndexEntry, index, z_file->sa_groups_index);
    for (uint64_t i=0; i < grp_len; i++) {
        rom grp_qname = GRP_QNAME(&grp[i]);
        index[i] = (SAGroupIndexEntry){ .grp_i = i, .qname_hash = QNAME_HASH (grp_qname, grp[i].qname_len, grp[i].is_last) };
    }

    qsort (index, index_len, sizeof(SAGroupIndexEntry), group_index_sorter);
}

void sam_sa_prim_initialize_ingest (void)
{
    mutex_initialize (seq_mutex);
    mutex_initialize (qual_mutex);
    mutex_initialize (qname_mutex);
    mutex_initialize (aln_mutex);
    mutex_initialize (grp_mutex);

    // suppress buf_alloc warning of big allocations
    z_file->sa_seq.can_be_big = z_file->sa_qual.can_be_big = true;
}

// ZIP: called from main thread by sam_zip_after_compute after final PRIM vb
void sam_sa_prim_finalize_ingest (void)
{
    START_TIMER;

    // build index by crc32(qname). 
    sam_sa_create_group_index();

    if (flag.show_sa) sam_show_sa();

    COPY_TIMER_VB (evb, sam_sa_prim_finalize_ingest);
}

// separately compress qual of each PRIM line - into overlaid buffer. returns total qual length.
static uint32_t sam_zip_prim_ingest_vb_compress_qual (VBlockSAMP vb, SAGroupType *vb_grps, uint32_t vb_grps_len,
                                                      BufferP underlying_buf, BufferP comp_qual_buf)
{
    // initialize total_qual_len to total uncompressed size of QUAL 
    uint32_t total_qual_len=0;
    for (uint32_t grp_i=0; grp_i < vb_grps_len; grp_i++) 
        total_qual_len += (vb_grps[grp_i].no_qual ? 0 : vb_grps[grp_i].seq_len);

    if (!total_qual_len) return 0;

    // note: in an unlikely case, the compressed size might be beyond this - in which case we will abandon the compression.
    uint32_t max_comp_len = total_qual_len / 1.2 + rans_compress_bound_4x16(vb_grps[0].seq_len, X_NOSZ); // heuristic (we don't want this to be unnecessarily too big - bigger than allocated to z_data)
    buf_alloc (vb, underlying_buf, max_comp_len, 0, char, 0, "z_data"); // likely already allocated
    buf_overlay_partial (vb, comp_qual_buf, underlying_buf, underlying_buf->len/*after packed_seq_buf*/, "comp_qual_buf");

    for (uint32_t vb_grp_i=0; vb_grp_i < vb_grps_len; vb_grp_i++) {
        
        SAGroupType *vb_grp = &vb_grps[vb_grp_i];
        if (vb_grp->no_qual) continue;

        uint32_t comp_len = rans_compress_bound_4x16 (vb_grp->seq_len, X_NOSZ);
        if (comp_qual_buf->len + comp_len <= max_comp_len)  { // still room to compress
            uint8_t *qual = B8 (vb->txt_data, vb_grp->qual); // always in SAM format
            ASSERT (rans_compress_to_4x16 (evb, qual, vb_grp->seq_len, BAFT(uint8_t, *comp_qual_buf), &comp_len, X_NOSZ) && comp_len,
                    "%s: Failed to compress PRIM qual: qual_len=%u", LN_NAME, vb_grp->seq_len);
        }

        // add to buffer only if compression actually compresses (otherwise it remains pointing at txt_data)
        if (comp_len < vb_grp->seq_len) {
            vb_grp->qual          = comp_qual_buf->len; // if compressed - index in comp_qual_buf, if not - remains index in txt_data
            vb_grp->qual_comp_len = comp_len; 
            comp_qual_buf->len   += comp_len;    
            total_qual_len        = total_qual_len + comp_len - vb_grp->seq_len;
        } 
    }

    return total_qual_len;
}

static void sam_zip_prim_ingest_vb_copy_qual_vb_to_z (VBlockSAMP vb, SAGroupType *vb_grps, uint32_t vb_grps_len, ConstBufferP comp_qual_buf, uint32_t total_qual_len)
{
    // concatenate this VB's compressed qual to z_file->sa_qual
    buf_alloc (evb, &z_file->sa_qual, total_qual_len, 0, uint8_t, 1, "z_file->sa_qual");
    
    // short cut if all groups have their QUAL compressed
    if (total_qual_len == comp_qual_buf->len) {
        // update index to be into z_file->sa_qual instead of comp_qual_buf
        for (uint32_t vb_grp_i=0; vb_grp_i < vb_grps_len; vb_grp_i++) 
            vb_grps[vb_grp_i].qual += z_file->sa_qual.len;

        buf_add_buf (evb, &z_file->sa_qual, comp_qual_buf, char, NULL);
    }
    
    else 
        for (uint32_t vb_grp_i=0; vb_grp_i < vb_grps_len; vb_grp_i++) {
            SAGroupType *vb_grp = &vb_grps[vb_grp_i];
            uint64_t z_qual = z_file->sa_qual.len;

            // case: qual is compressed - copy from comp_qual_buf
            if (vb_grp->qual_comp_len) 
                buf_add (&z_file->sa_qual, Bc(*comp_qual_buf, vb_grp->qual), vb_grp->qual_comp_len);

            // case: qual is not compressed - copy from txt_data
            else
                buf_add (&z_file->sa_qual, Btxt(vb_grp->qual), vb_grp->seq_len); // note: qual in txt_data is always textual (in BAM, qual was re-written to textual in bam_rewrite_qual)

            vb_grp->qual = z_qual; // now its an index into z_file->sa_qual
        }
}

static void sam_zip_prim_ingest_vb_copy_qname_vb_to_z (VBlockSAMP vb, SAGroupType *vb_grps, uint32_t vb_grps_len, uint32_t total_qname_len)
{
    buf_alloc (evb, &z_file->sa_qnames, total_qname_len, 200000, char, CTX_GROWTH, "z_file->sa_qnames");

    for (uint32_t grp_i=0; grp_i < vb_grps_len; grp_i++) {
        
        SAGroupType *g = &vb_grps[grp_i];

        // qname. case: collated files - mates are consecutive - in case both of primary, we can just use the existing qname
        if (segconf.sam_is_collated && grp_i && str_issameR_(STRtxt(g->qname), STRtxt((g-1)->qname)))
            g->qname = (g-1)->qname; // same index into z_file->sa_qnames like previous qname
        
        else {
            uint64_t z_qname = z_file->sa_qnames.len;;
            memcpy (BAFTc (z_file->sa_qnames), Btxt (g->qname), g->qname_len);
            g->qname = z_qname; // update from index into vb->txt_data to index into z_file->sa_qname
            z_file->sa_qnames.len += g->qname_len;
        }
    }
}

// ZIP PRIM VB: called by compute thread, from sam_seg_finalize - VBs are out of order
// Called after segging a PRIM VB, to ingest its data into vb->sa_*, for later use by DEPN VBs. 
void sam_zip_prim_ingest_vb (VBlockSAMP vb)
{
    START_TIMER;

    ARRAY (SAGroupType, vb_grps, vb->sa_groups);

    if (!vb_grps_len) return; // no SA groups

    // overlay buffers - we can use automatic vars as overlaid buffers are not added to the buffer list
    Buffer comp_qual_buf={}, packed_seq_buf={};

    // using z_data buf that is already allocated but still empty - we're likely not going to allocate additional memory
    ASSERT (!vb->z_data.len, "expecting vb->z_data.len=0 but it is %"PRIu64, vb->z_data.len);
    
    // compress seq and qual of all lines into packed_seq_buf & comp_qual_buf (overlaid on z_data, to save allocating memory)
    sam_zip_prim_ingest_vb_pack_seq (vb, STRa(vb_grps), &vb->z_data, &packed_seq_buf, IS_BAM_ZIP);
    
    uint32_t total_qual_len = sam_zip_prim_ingest_vb_compress_qual (vb, STRa(vb_grps), &vb->z_data, &comp_qual_buf);
    vb->comp_qual_len = total_qual_len; // goes into SectionHeaderVbHeader.sam_prim_comp_qual_len 

    // no mutex locked: calculate qname memory requirements
    uint32_t total_qname_len=0;
    for (uint32_t grp_i=0; grp_i < vb_grps_len; grp_i++)
        total_qname_len += vb_grps[grp_i].qname_len;

    // Note: Buffers are pre-INITed in file_initialize_bufs, so we can buf_alloc them from serialized compute threads
    // Note: We must serialize the data populating too, bc otherwise other VBs can re-alloc the buffers under our feet
    // VBs can enter in arbitrary order, and also handle SEQ, QUAL, QUAL, ALN in arbitrary order by availability of the mutex
    
    bool seq_done=false, qual_done=false, qname_done=false, aln_done=false;

    while (!seq_done || !qual_done || !qname_done || !aln_done) {

        bool achieved_something = false;
        
        if (!seq_done && mutex_trylock (seq_mutex)) {
            // concatenate this VB's sequence to z_file->sa_seq (in ACGT format)
            uint64_t z_seq = z_file->sa_seq.nbits / 2; // next_seq in bases, not bits
            bits_concat (buf_get_bitarray(&z_file->sa_seq), buf_get_bitarray (&packed_seq_buf), 0); // also allocate memory
            seq_done = achieved_something = true;                
            mutex_unlock (seq_mutex);

            for (uint32_t grp_i=0; grp_i < vb_grps_len; grp_i++)
                vb_grps[grp_i].seq += z_seq; // update index from index into vb->sa_seq to z_file->sa_seq
        } 

        if (!qual_done && mutex_trylock (qual_mutex)) {
            sam_zip_prim_ingest_vb_copy_qual_vb_to_z (vb, vb_grps, vb_grps_len, &comp_qual_buf, total_qual_len);
            qual_done = achieved_something = true;
            mutex_unlock (qual_mutex);
        }

        if (!qname_done && mutex_trylock (qname_mutex)) {
            sam_zip_prim_ingest_vb_copy_qname_vb_to_z (vb, vb_grps, vb_grps_len, total_qname_len);
            qname_done = achieved_something = true;
            mutex_unlock (qname_mutex);
        }

        if (!aln_done && mutex_trylock (aln_mutex)) {
            uint64_t z_first_aln = z_file->sa_alns.len;
            buf_add_buf (evb, &z_file->sa_alns, &vb->sa_alns, SAAlnType, NULL); // copy alignments (also allocs memory)
            aln_done = achieved_something = true;
            mutex_unlock (aln_mutex);

            for (uint32_t grp_i=0; grp_i < vb_grps_len; grp_i++)
                vb_grps[grp_i].first_aln_i += z_first_aln; // update index from index into vb->sa_alns to z_file->sa_alns
        }

        if (!achieved_something) usleep (10000); // 10 ms
    }

    vb_grps[0].first_grp_in_vb = true; // mark the first group in this PRIM VB

    // Groups - must be last, after vb_grps is updated ^. Order of VBs is arbitrary
    mutex_lock (grp_mutex);
    vb->first_grp_i = z_file->sa_groups.len; // the index of first group of this PRIM VB, in z_file->sa_groups. note: might be out-of-order of VBs.
    buf_add_buf (evb, &z_file->sa_groups, &vb->sa_groups, SAGroupType, NULL); 
    mutex_unlock (grp_mutex);    

    buf_free (packed_seq_buf);
    buf_free (comp_qual_buf);
    buf_free (vb->z_data);

    COPY_TIMER_VB (evb, sam_zip_prim_ingest_vb); 
}
