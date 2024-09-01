// ------------------------------------------------------------------
//   vblock.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

// vb stands for VBlock - it started its life as VBlockVCF when genozip could only compress VCFs, but now
// it means a block of lines from the text file. 

#include "vblock.h"
#include "file.h"
#include "digest.h"
#include "mgzip.h"
#include "threads.h"
#include "writer.h"

// pool of VBs allocated based on number of threads
static VBlockPoolP pools[NUM_POOL_TYPES] = {};

VBlockP evb = NULL; // outside a pool

VBlockPool *vb_get_pool (VBlockPoolType type, FailType soft_fail)
{
    ASSERT (pools[type] || soft_fail, "VB Pool type=%u is not allocated", type);
    return pools[type];
}

VBlockP vb_get_from_pool (VBlockPoolP pool, VBID vb_id) 
{
    ASSERTNOTNULL (pool);

    if (!(vb_id == VB_ID_EVB || (vb_id >= 0 && vb_id < pool->num_vbs)))
        return NULL; // soft fail on invalid vb_id

    return (vb_id == VB_ID_EVB) ? evb : pool->vb[vb_id];
}

static VBlockP nonpool_vbs[NUM_NONPOOL_VBs] = {}; 

static inline bool is_in_use (VBlockP vb)
{
    return __atomic_load_n (&vb->in_use, __ATOMIC_ACQUIRE);
}

static void set_in_use (VBlockP vb, bool in_use)
{
    __atomic_store_n (&vb->in_use, in_use, __ATOMIC_RELEASE);   
}

void vb_release_vb_do (VBlockP *vb_p, rom task_name, rom func)
{
    START_TIMER;

    VBlockP vb = *vb_p;
    if (!vb) return; // nothing to release
    
    ASSERT (is_in_use (vb) || vb==evb, "Cannot release VB because it is not in_use (called from %s): vb->id=%d vb->vblock_id=%u", 
            func, vb->id, vb->vblock_i);

    threads_log_by_vb (vb, vb->compute_task ? vb->compute_task : func, "RELEASING VB", 0);

    if (flag.show_time_comp_i == vb->comp_i || flag.show_time_comp_i == COMP_ALL)
        profiler_add (vb);

    if (vb->id != VB_ID_EVB) // cannot test evb, see comment in buflist_test_overflows_do
        buflist_test_overflows(vb, func); 

    // verify that gzip_compressor was released after use
    ASSERT (!vb->gzip_compressor, "vb=%s: expecting gzip_compressor=NULL", VB_NAME);

    // release all buffers in vb->buffer_list, and zero the space between these buffers
    buflist_free_vb (vb); 

    // this release can be run by either the main or writer thread. we make sure to update in_use as the very
    // last change, and do so atomically

    // case: this VB is from the pool (i.e. not evb)
    int32_t num_in_use = -1;
    if (vb != evb) {
        // Logic: num_in_use is always AT LEAST sum(vb)->in_use. i.e. pessimistic. (it can be mometarily less between these two updates)
        set_in_use (vb, false);   // released the VB back into the pool - it may now be reused 
        if (vb->id >= 0) 
            num_in_use = __atomic_sub_fetch (&pools[vb->pool]->num_in_use, 1, __ATOMIC_ACQ_REL); // atomic to prevent concurrent update by writer thread and main thread (must be after update of in_use)
        *vb_p = NULL;
    }

    if (flag_is_show_vblocks (task_name)) 
        iprintf (vb->id >= 0 ? "VB_RELEASE(task=%s id=%d) vb=%s caller=%s in_use=%d/%d\n"
                             : "VB_RELEASE(task=%s id=%d) vb=%s caller=%s\n", 
                 task_name, vb->id, VB_NAME, func, num_in_use, (vb->id >= 0 ? pools[vb->pool]->num_vbs : -1));

    buflist_compact (vb);

    if (vb->id >= 0) COPY_TIMER_EVB (vb_release_vb_do);
}


void vb_destroy_vb_do (VBlockP *vb_p, rom func)
{
    ASSERTMAINTHREAD;
    START_TIMER;

    VBlockP vb = *vb_p;
    if (!vb) return;

    if (flag_is_show_vblocks (vb->compute_task)) 
        iprintf ("VB_DESTROY(id=%d) vb_i=%d caller=%s\n", vb->id, vb->vblock_i, func);

    bool is_evb = vb->id == VB_ID_EVB;

    buflist_destroy_vb_bufs (vb, false);

    // case: this is a nonpool VB
    for (int i=0; i < NUM_NONPOOL_VBs; i++)
        if (nonpool_vbs[i] == vb) 
            nonpool_vbs[i] = NULL;

    FREE (*vb_p);

    if (!is_evb) COPY_TIMER_EVB (vb_destroy_vb);
}

// return all VBlocks memory and unused evb memory to libc and optionally to the kernel
void vb_dehoard_memory (bool release_to_kernel)
{
    vb_destroy_pool_vbs (POOL_MAIN, false); 
    buflist_destroy_vb_bufs (evb, true); // destroys all unused buffers

    if (release_to_kernel)
        buf_low_level_release_memory_back_to_kernel();
}

void vb_create_pool (VBlockPoolType type, rom name)
{
    // only main-thread dispatcher can create a pool. other dispatcher (eg writer's bgzf compression) can must existing pool
    uint32_t num_vbs = (type == POOL_MAIN) ? MAX_(1, global_max_threads) +               // compute thread VBs
                                             (IS_PIZ ? 2 : 0)    +                       // txt header VB and wvb (for PIZ)
                                             (IS_PIZ ? z_file->max_conc_writing_vbs : 0) // SAM: max number of thread-less VBs handed over to the writer thread which the writer must load concurrently 
                                           : writer_get_max_bgzf_threads();
    if (flag_is_show_vblocks (NULL)) 
        iprintf ("CREATING_VB_POOL: type=%s global_max_threads=%u max_conc_writing_vbs=%u num_vbs=%u\n", 
                 name, global_max_threads, z_file->max_conc_writing_vbs, num_vbs); 

    uint32_t size = sizeof (VBlockPool) + num_vbs * sizeof (VBlockP);

    if (!pools[type])  
        // allocation includes array of pointers (initialized to NULL)
        pools[type] = (VBlockPool *)CALLOC (size); // note we can't use Buffer yet, because we don't have VBs yet...

    // case: old pool is too small - realloc it (the pool contains only pointers to VBs, so the VBs themselves are not realloced)
    else if (pools[type]->num_vbs < num_vbs) {
        REALLOC (&pools[type], size, "VBlockPool"); 
        memset (&pools[type]->vb[pools[type]->num_vbs], 0, (num_vbs - pools[type]->num_vbs) * sizeof (VBlockP)); // initialize new entries
    }

    pools[type]->name    = name;
    pools[type]->size    = size; 
    pools[type]->num_vbs = MAX_(num_vbs, pools[type]->num_vbs); 
}

VBlockP vb_initialize_nonpool_vb (VBID vb_id, DataType dt, rom task)
{
    VBlockP vb            = CALLOC (get_vb_size (dt));
    vb->data_type         = DT_NONE;
    vb->id                = vb_id;
    vb->compute_task      = task;
    vb->data_type         = dt;
    vb->data_type_alloced = dt;
    vb->comp_i            = COMP_NONE;
    vb->pool              = NO_POOL;
    init_dict_id_to_did_map (vb->d2d_map); 
    
    if (!vb->buffer_list.vb) {
        vb->buffer_list.name = "buffer_list";
        buf_init_lock (&vb->buffer_list);
        buflist_add_buf (vb, &vb->buffer_list, __FUNCLINE);
        vb->buffer_list.vb = vb; // indication buffer was added to buffer list
    }

    nonpool_vbs[NUM_NONPOOL_VBs + vb_id] = vb; // vb_id is a negative integer
    set_in_use (vb, true);

    return vb;
}

VBlockP vb_get_nonpool_vb (VBID vb_id)
{
    return nonpool_vbs[NUM_NONPOOL_VBs + vb_id]; // may be NULL
}

// used to change segconf VB data_type is seg_initiatlize (FASTA->FASTQ)
void vb_change_datatype_nonpool_vb (VBlockP *vb_p, DataType new_dt)
{
    REALLOC (vb_p, get_vb_size (new_dt), "VBlock");
    VBlockP vb = *vb_p;

    vb->data_type = vb->data_type_alloced = new_dt;
    init_dict_id_to_did_map (vb->d2d_map); 
    
    buflist_update_vb_addr_change (vb, vb->buffer_list.vb);

    nonpool_vbs[NUM_NONPOOL_VBs + vb->id] = vb;
}

static VBlockP vb_update_data_type (VBlockP vb, DataType dt, DataType alloc_dt, uint64_t sizeof_vb)
{
    if (z_file && vb->data_type == dt) return vb;

    // the new data type has a private section in its VB, that is different that the one of alloc_dt - realloc private section
    if (vb->data_type_alloced != alloc_dt && alloc_dt != DT_NONE) {

        // destroy private part of previous data_type. we also destroy all contexts as new data type is going
        // to allocate different contexts with a different memory usage profile (eg b250 vs local) for each
        if (vb->data_type_alloced != DT_NONE) 
            buflist_destroy_private_and_context_vb_bufs (vb); 
        
        buflist_compact (vb); // remove buffer_list entries marked for removal

        // calloc private part of new data_type
        VBlockP old_vb = vb;
        REALLOC (&vb, sizeof_vb, "VBlock");

        // initialize private part of new data_type, keeping common part intact
        memset ((char*)vb + sizeof (VBlock), 0, sizeof_vb - sizeof (VBlock));

        // update buf->vb in all buffers of this VB to new VB address
        if (vb != old_vb) buflist_update_vb_addr_change (vb, old_vb);
        vb->data_type_alloced = alloc_dt;
    }

    vb->data_type = dt;
    return vb;
}

// allocate an unused vb from the pool. seperate pools for zip and unzip
VBlockP vb_get_vb (VBlockPoolType type, rom task_name, VBIType vblock_i, CompIType comp_i)
{
    START_TIMER;

    VBlockPool *pool = vb_get_pool (type, HARD_FAIL);

#ifdef DEBUG
    // if GFF VB becauses larger than FASTA, then we need to change the dt assignment conditions below
    ASSERT0 (get_vb_size (DT_FASTA) > get_vb_size (DT_GFF), "Failed assumption that FASTA has larger VB than GFF");
#endif

    DataType dt = (type == POOL_BGZF)                             ? DT_NONE
                : (flag.deep && flag.zip_comp_i >= SAM_COMP_FQ00) ? DT_FASTQ
                : (IS_ZIP && segconf.has_embedded_fasta)          ? DT_FASTA // GFF3 with embedded FASTA (allocate a FASTA VB as it larger than GFF and can accommodate both)
                : (IS_ZIP && txt_file)                            ? txt_file->data_type
                : (IS_PIZ && z_file && Z_DT(GFF))                 ? DT_FASTA // vb is allocated before reading the VB_HEADER, so we don't yet know if it is a embdedded fasta VB. To be safe, we allocate enough memory for a VBlockFASTA (which is larger), so we can change the dt in gff_piz_init_vb() without needing to realloc
                : (IS_PIZ && z_file && flag.deep && comp_i != COMP_NONE && comp_i >= SAM_COMP_FQ00) ? DT_FASTQ
                : (IS_PIZ && z_file)                              ? z_file->data_type  
                :                                                   DT_NONE;
    
    uint64_t sizeof_vb = get_vb_size (dt);

    DataType alloc_dt = sizeof_vb == sizeof (VBlock) ? DT_NONE
                      : (dt == DT_REF && IS_PIZ)     ? DT_NONE
                      : dt == DT_BAM                 ? DT_SAM
                      : dt == DT_BCF                 ? DT_VCF 
                      :                                dt;

    if (type == POOL_MAIN && IS_PIZ && z_file && Z_DT(GFF)) 
        dt = DT_GFF; // return GFF dt to its true dt after getting the size, otherwise it won't work   

    // circle around until a VB becomes available (busy wait)
    VBID vb_id; for (vb_id=0; ; vb_id = (vb_id+1) % pool->num_vbs) {
        if (!pool->vb[vb_id]) { // VB is not allocated - allocate it
            pool->vb[vb_id] = CALLOC (sizeof_vb); 
            pool->num_allocated_vbs++;
            pool->vb[vb_id]->data_type_alloced = alloc_dt;
            break;
        }

        else if (!is_in_use (pool->vb[vb_id])) {
            pool->vb[vb_id] = vb_update_data_type (pool->vb[vb_id], dt, alloc_dt, sizeof_vb); // possibly realloc
            break;
        }

        // case: we've checked all the VBs and none is available - wait a bit and check again
        // in PIZ, this happens when a lot VBs are handed over to the writer thread which has not processed them yet.
        // for example, if writer is blocking on write(), waiting for a pipe counterpart to read.
        // it will be released when the writer thread completes one VB.
        if (vb_id == pool->num_vbs-1) usleep (50000); // 50 ms
    }

    VBlockP vb = pool->vb[vb_id]; 

    // Logic: num_in_use is always AT LEAST sum(vb)->in_use. i.e. pessimistic. (it can be mometarily more between these two updates)
    uint32_t num_in_use = __atomic_add_fetch (&pool->num_in_use, 1, __ATOMIC_ACQ_REL); // atomic to prevent concurrent update by writer thread and main thread (must be before update of in_use)
    set_in_use (vb, true);

    // initialize VB fields that need to be a value other than 0
    vb->id                = vb_id;
    vb->data_type         = dt;
    vb->vblock_i          = vblock_i;
    vb->comp_i            = comp_i;
    vb->compute_thread_id = THREAD_ID_NONE;
    vb->compute_task      = task_name;
    vb->pool              = type;
    init_dict_id_to_did_map (vb->d2d_map);

    if (!vb->buffer_list.vb) {
        vb->buffer_list.name = "buffer_list";
        buf_init_lock (&vb->buffer_list);
        buflist_add_buf (vb, &vb->buffer_list, __FUNCLINE);
        vb->buffer_list.vb = vb; // indication buffer was added to buffer list
    }
    
    if (flag_is_show_vblocks (task_name)) 
        iprintf ("VB_GET_VB(task=%s id=%u) vb_i=%s/%d num_in_use=%u/%u%s\n", 
                  task_name, vb->id, comp_name (vb->comp_i), vb->vblock_i, num_in_use, pool->num_vbs,
                  flag.preprocessing ? " preprocessing" : "");

    threads_log_by_vb (vb, task_name, "GET VB", 0);

    if (flag.debug_memory)
        iprintf ("vb_get_vb: got vb_i=%d id=%d task=%s dt=%s address=[%p - %p]\n", 
                 vblock_i, vb_id, task_name, dt_name(alloc_dt), vb, (char*)vb + sizeof_vb);

    COPY_TIMER_EVB (vb_get_vb);
    return vb;
}

uint32_t vb_pool_get_num_in_use (VBlockPoolType type, VBID *id/*optional out*/)
{
    VBlockPool *pool = vb_get_pool (type, HARD_FAIL);
    int num_in_use = __atomic_load_n (&pool->num_in_use, __ATOMIC_ACQUIRE); // atomic, bc for POOL_MAIN, writer thread might update concurrently.

    if (id) {
        *id = VB_ID_NONE;
        if (num_in_use) {
            for (VBID vb_id=0; vb_id < pool->num_vbs; vb_id++)
                if (pool->vb[vb_id] && pool->vb[vb_id]->in_use) 
                    return *id = vb_id;

            if (*id == VB_ID_NONE) // all lost use while we were checking
                num_in_use = 0;
        }
    }

    return num_in_use;
}

// Note: num_in_use is always AT LEAST sum(vb)->in_use (it can be mometarily less than sum(vb)->in_use as they are getting updated)
// therefore, this function may return true when pool is actually no longer full.
bool vb_pool_is_full (VBlockPoolType type)
{
    return vb_pool_get_num_in_use (type, NULL) == vb_get_pool(type, HARD_FAIL)->num_vbs;
}
 
// Note: As in vb_pool_is_full, if the function returns false (not empty) the pool might in fact already be empty
bool vb_pool_is_empty (VBlockPoolType type)
{
    return vb_pool_get_num_in_use (type, NULL) == 0;
}

bool vb_is_valid (VBlockP vb)
{
    if (vb->pool != NO_POOL)
        for (VBID vb_id=0; vb_id < pools[vb->pool]->num_vbs; vb_id++)
            if (vb == pools[vb->pool]->vb[vb_id]) return true;

    for (int i=0; i < NUM_NONPOOL_VBs; i++)
        if (nonpool_vbs[i] == vb) return true;

    return false;
}

// frees memory of all VBs, except for non-pool VBs (evb, segconf, writer,...)
void vb_destroy_pool_vbs (VBlockPoolType type, bool destroy_pool)
{
    if (!pools[type]) return;

    for (VBID vb_id=0; vb_id < pools[type]->num_vbs; vb_id++) 
        vb_destroy_vb (&pools[type]->vb[vb_id]);

    if (destroy_pool)
        FREE (pools[type]);
}

StrText err_vb_pos (void *vb)
{
    StrText s;

    snprintf (s.s, sizeof (s), "vb i=%u position in %s file=%"PRIu64, 
             (VB)->vblock_i, dt_name (txt_file->data_type), (VB)->vb_position_txt_file);
    return s;
}

unsigned def_vb_size (DataType dt) { return sizeof (VBlock); }

void vb_set_is_processed (VBlockP vb)
{
    __atomic_store_n (&vb->is_processed, (bool)true, __ATOMIC_RELEASE); 
}

bool vb_is_processed (VBlockP vb)
{
    return __atomic_load_n (&vb->is_processed, __ATOMIC_ACQUIRE);
}

bool vb_buf_locate (VBlockP vb, ConstBufferP buf)
{
    if (!vb) return false;

    unsigned sizeof_vb = get_vb_size (vb->data_type_alloced);

    return vb && is_p_in_range (buf, vb, sizeof_vb);
}

rom textual_assseg_line (VBlockP vb)
{
    if (vb->line_start >= Ltxt) return "Invalid line_start";

    char *nl = memchr (Btxt(vb->line_start), '\n', Ltxt - vb->line_start);
    if (!nl) nl = BAFTtxt; // possibly overwriting txt_data's overflow fence

    *nl = 0; // terminate string
    return Btxt(vb->line_start);
}

static void vb_deferred_q_reorder (VBlockP vb, Did did_i, int q_index, int depth)
{
    ASSERT (depth < 10, "Cyclic seg order requirements, ctx=%s", CTX(did_i)->tag_name);

    for (int i=0; i < q_index; i++)
        if (vb->deferred_q[i].seg_after_did_i == did_i) {
            // move element i to after element q_index
            // example start: [A, B, C, D] (B.seg_after_did_i=D). result: [A, C, D, B]
            DeferredField df = vb->deferred_q[i];
            memmove (&vb->deferred_q[i], &vb->deferred_q[i+1], (q_index - i) * sizeof (DeferredField));
            vb->deferred_q[q_index] = df;

            // now, move any element that needs to be after i (B in the example)
            vb_deferred_q_reorder (vb, df.did_i, q_index, depth+1);
        }
}

void vb_add_to_deferred_q (VBlockP vb, ContextP ctx, DeferredSeg seg, int16_t idx,
                           Did seg_after_did_i) // optional (DID_NONE if not) - ctx cannot be segged before seg_after_did_i (= another context that might be in the deferred queue)
{
    ASSERT (vb->deferred_q_len+1 < DEFERRED_Q_SZ, "%s: deferred queue is full (deferred_q_len=%u) when adding %s", VB_NAME, vb->deferred_q_len, ctx->tag_name);
    ASSERT (idx >= 0, "Invalid idx=%d when adding %s", idx, ctx->tag_name);

    vb->deferred_q[vb->deferred_q_len++] = (DeferredField){ .did_i=ctx->did_i, .seg=seg, .idx=idx, .seg_after_did_i = seg_after_did_i };

    // change order of segging if needed
    vb_deferred_q_reorder (vb, ctx->did_i, vb->deferred_q_len-1, 1);
}

void vb_display_deferred_q (VBlockP vb, rom func)
{
    if (!vb->deferred_q_len) return;

    iprintf ("%s: %s Deferred seg queue: ", func, LN_NAME);

    for (int i=0; i < vb->deferred_q_len; i++)
        iprintf ("%s ", CTX(vb->deferred_q[i].did_i)->tag_name);

    iprint0 ("\n");
}
