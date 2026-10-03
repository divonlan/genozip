// ------------------------------------------------------------------
//   vblock.c
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
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
#include "dispatcher.h"

// pool of VBs allocated based on number of threads
static VBlockPoolP pools[NUM_POOLS] = {};

VBlockP evb = NULL; // a copy of pools[POOL_NONPOOL].vb[VB_ID_EVB]

rom pool_name (VBlockPoolType type)
{
    if (IN_RANGE (type, 0, NUM_POOLS)) 
        return (rom[])POOL_NAMES[type];

    else
        return "(Invalid pool name)";
}

StrText dis_vb_id (VBID vb_id)
{
    StrText s;
    
    switch (vb_id.pool) {
        case NO_POOL:
            strcpy (s.s, "NONE");
            break;

        case POOL_NONPOOL:
            if(IN_RANGE(vb_id.index, 0, NUM_NONPOOL_VBs))
                snprintf (s.s, sizeof(s), "%s", (rom[])NONPOOL_NAMES[vb_id.index]);
            else
                snprintf (s.s, sizeof(s), "(Invalid NONPOOL index=%u)", vb_id.pool);
            break;

        case POOL_MAIN: case POOL_BGZF: case POOL_MISC:
            snprintf (s.s, sizeof(s), "%s:%u", pool_name (vb_id.pool), vb_id.index);
            break;

        default:
            snprintf (s.s, sizeof(s), "(Invalid pool=%u)", vb_id.pool);
    }

    return s;
}

VBlockPool *vb_get_pool (VBlockPoolType type, FailType soft_fail)
{
    ASSERTINRANGE(type, 1, NUM_POOLS);

    ASSERT (pools[type] || soft_fail, "VB Pool %s is not allocated", pool_name (type));

    return pools[type];
}

VBlockP vb_get_from_pool (VBlockPoolP pool, VBIDIndex index) 
{
    ASSERTNOTNULL (pool);

    if (index >= pool->num_vbs)
        return NULL; // soft fail on invalid index

    return pool->vb[index];
}

static inline bool is_in_use (VBlockP vb)
{
    return load_acquire (vb->in_use);
}

static void set_in_use (VBlockP vb, bool in_use)
{
    store_release (vb->in_use, in_use);   
}

void vb_release_vb_do (VBlockP *vb_p, rom func)
{
    START_TIMER;

    VBlockP vb = *vb_p;
    if (!vb) return; // nothing to release

    ASSERT (is_in_use (vb) || vb==evb, "Cannot release VB because it is not in_use (called from %s): vb->id.index=%d vb->vblock_id=%u", 
            func, vb->id.index, vb->vblock_i);

    Task task = vb->compute_task;
    threads_log_by_vb (vb, task ? task_name (task) : func, "RELEASING VB", 0);

    if (flag.show_time_comp_i == vb->comp_i || flag.show_time_comp_i == COMP_ALL)
        𝓅𝓇ℴ𝒻𝒾𝓁ℯ (profiler_add (vb));

    if (!(vb->id.pool == POOL_NONPOOL && vb->id.index == VB_ID_EVB)) // cannot test evb, see comment in buflist_test_overflows_do
        buflist_test_overflows(vb, func); 

    // verify that bgzf memory was released after use
    ASSERTISNULL (vb->libdef_decomp_mem);
    ASSERTISNULL (vb->gz_deflate_mem);

    // release all buffers in vb->buffer_list, and zero the space between these buffers
    buflist_free_vb (vb); 

    buflist_compact (vb);

    // IMPORTANT: this release can be run by either the main or writer thread. 
    // we make sure to update in_use as the very last change, and do so atomically

    VBID vb_id = vb->id; // copy before releasing

    // case: this VB is from the pool (i.e. not evb)
    int32_t num_in_use = -1;
    if (vb_id.pool != POOL_NONPOOL) {
        set_in_use (vb, false);   // released the VB back into the pool - it may now be reused 

        __atomic_thread_fence (__ATOMIC_RELEASE); 

        // **** NO ACCESS TO vb AFTER THIS POINT (in PIZ, it might be realloced by the next user) ****

        // Logic: num_in_use is always AT LEAST sum(vb)->in_use. i.e. pessimistic. (it can be mometarily less between these two updates)
        num_in_use = decrement_relaxed (pools[vb_id.pool]->num_in_use); // atomic to prevent concurrent update by writer thread and main thread (must be after update of in_use)
        *vb_p = NULL;
    }

    if (flag_is_show_vblocks (task)) 
        iprintf ("VB_RELEASE(task=%s id=%d) vb=%s caller=%s%s%s\n",
                 task_name (task), vb_id.index, VB_NAME, func, 
                 cond_int (vb_id.pool != POOL_NONPOOL, " in_use=", num_in_use), 
                 cond_int (vb_id.pool != POOL_NONPOOL, "/", pools[vb_id.pool]->num_vbs));

    if (vb_id.pool != POOL_NONPOOL) 
        COPY_TIMER_EVB (vb_release_vb_do);
}


void vb_destroy_vb_do (VBlockP *vb_p, rom func)
{
    ASSERTMAINTHREAD;
    START_TIMER;

    VBlockP vb = *vb_p;
    if (!vb) return;

    if (flag_is_show_vblocks (vb->compute_task)) 
        iprintf ("VB_DESTROY(id=%d) vb_i=%d caller=%s\n", vb->id.index, vb->vblock_i, func);

    pools[vb->id.pool]->vb[vb->id.index] = NULL; // remove from pool

    bool is_evb = (vb->id.pool == POOL_NONPOOL && vb->id.index == VB_ID_EVB);

    buflist_destroy_vb_bufs (vb, false);

    FREE_ALIGN_64B (*vb_p);

    if (!is_evb) COPY_TIMER_EVB (vb_destroy_vb); //can't store profiling in evb after it is destroyed...
}

// return all VBlocks memory and unused evb memory to libc and optionally to the kernel
void vb_dehoard_memory (bool release_to_kernel)
{
    vb_destroy_pool (POOL_MAIN, false); 
    buflist_destroy_vb_bufs (evb, true); // destroys all unused buffers

    if (release_to_kernel)
        return_freed_memory_to_kernel();
}

void vb_create_pool (VBlockPoolType type)
{
    // only main-thread dispatcher can create a pool. other dispatcher (eg writer's bgzf compression) can must existing pool
    uint32_t num_vbs = 
        (type == POOL_BGZF)    ? writer_get_max_bgzf_threads()
      : (type == POOL_MISC)    ? MAX_(1, global_max_threads)
      : (type == POOL_NONPOOL) ? NUM_NONPOOL_VBs
      : /*POOL_MAIN*/            MAX_(1, global_max_threads) + // compute thread VBs
                                 (IS_PIZ ? 2 : 0) + 
                                 (IS_PIZ && !flag.no_writer_thread ? z_file->max_conc_writing_vbs : 0); // SAM: max number of thread-less VBs handed over to the writer thread which the writer must load concurrently 
    
    num_vbs = MIN_(num_vbs, MAX_POOL_VBS);

    if (flag_is_show_vblocks (TASK_NONE)) 
        iprintf ("CREATING_VB_POOL: type=%s global_max_threads=%u max_conc_writing_vbs=%u num_vbs=%u\n", 
                 pool_name (type), global_max_threads, z_file->max_conc_writing_vbs, num_vbs); 

    uint32_t size = sizeof (VBlockPool) + num_vbs * sizeof (VBlockP);

    if (!pools[type])  
        // allocation includes array of pointers (initialized to NULL)
        pools[type] = (VBlockPool *)CALLOC(size); // note we can't use Buffer yet, because we don't have VBs yet...

    // case: old pool is too small - realloc it (the pool contains only pointers to VBs, so the VBs themselves are not realloced)
    else if (pools[type]->num_vbs < num_vbs) {
        REALLOC (&pools[type], size, "VBlockPool"); 
        memset ((char *)&pools[type]->vb[pools[type]->num_vbs], 0, (num_vbs - pools[type]->num_vbs) * sizeof (VBlockP)); // initialize new entries
    }

    pools[type]->name    = pool_name (type);
    pools[type]->size    = size; 
    pools[type]->num_vbs = MAX_(num_vbs, pools[type]->num_vbs); 
}

extern void context_validate (void); // function autogenerated by context_validate.sh

VBlockP vb_initialize_nonpool_vb (VBIDIndex nonpool, DataType dt, Task task)
{
    // verify that VBlock buffers and contexts are 64B-aligned
    DO_ONCE {
        context_validate();
        
        ASSERT (sizeof(Context) % 64 == 0, "Expecting sizeof(Context)=%u (%1.3f 64B-blocks) to be a multiple of 64",
                (int)sizeof(Context), sizeof(Context) / 64.0);

        ASSERT (sizeof(ContextArray) % 64 == 0, "Expecting sizeof(ContextArray)=%u (%1.3f 64B-blocks) to be a multiple of 64",
                (int)sizeof(ContextArray), sizeof(ContextArray) / 64.0);

        // verify 64B alignment for CA and buffers (assuming that VBlock itself is aligned as it is allocated with CALLOC_ALIGN_64B)
        #define ASSERT_ALIGNED_64(field) \
            ASSERT (offsetof(VBlock, field) % 64 == 0, "VBlock field " #field ": byte offset %zu (64B block: %1.3f) is not 64-byte aligned", \
                    offsetof(VBlock, field), offsetof(VBlock, field) / 64.0)
        ASSERT_ALIGNED_64(ca);
        ASSERT_ALIGNED_64(txt_data); // first field after ca
        ASSERT_ALIGNED_64(codec_bufs[1]);
    }

    ASSERT (IN_RANGE(nonpool, 0, NUM_NONPOOL_VBs), "nonpool ∉ [0,%u]", NUM_NONPOOL_VBs-1);

    if (!pools[POOL_NONPOOL])
        vb_create_pool (POOL_NONPOOL);

    VBlockP vb = pools[POOL_NONPOOL]->vb[nonpool] = CALLOC_ALIGN_64B (get_vb_size (dt));

    vb->data_type         = DT_NONE;
    vb->id                = (VBID){ .pool = POOL_NONPOOL, .index = nonpool };
    vb->compute_task      = task;
    vb->data_type         = dt;
    vb->data_type_alloced = dt;
    vb->comp_i            = COMP_NONE;
    ca_init_d2d_map (&vb->ca); 
    
    if (!vb->buffer_list.vb) {
        vb->buffer_list.nameר = ר("buffer_list");
        buf_init_lock (&vb->buffer_list);
        buflist_add_buf (vb, &vb->buffer_list, THIS_CODE_LINE);
        vb->buffer_list.vb = vb; // indication buffer was added to buffer list
    }

    set_in_use (vb, true);

    return vb;
}

VBlockP vb_get_nonpool_vb (VBIDIndex nonpool)
{
    ASSERT (IN_RANGE(nonpool, 0, NUM_NONPOOL_VBs), "nonpool ∉ [0,%u]", NUM_NONPOOL_VBs-1);

    return pools[POOL_NONPOOL]->vb[nonpool]; // may be NULL
}

static VBlockP vb_update_data_type (VBlockP vb, DataType dt, DataType alloc_dt, uint64_t sizeof_vb)
{
    if (z_file && vb->data_type == dt) return vb;

    // the new data type has a private section in its VB, that is different that the one of alloc_dt - realloc private section
    if (vb->data_type_alloced != alloc_dt) {

        // destroy private part of previous data_type. we also destroy all contexts as new data type is going
        // to allocate different contexts with a different memory usage profile (eg b250 vs local) for each
        if (vb->data_type_alloced != DT_NONE) 
            buflist_destroy_private_and_context_vb_bufs (vb); 
        
        buflist_compact (vb); // remove buffer_list entries marked for removal

        // unfortunately there is no aligned_realloc: we allocated a new block and copy
        VBlockP old_vb = vb;
        vb = CALLOC_ALIGN_64B (sizeof_vb);

        // copy common part, leave new private part zeroed
        memcpy (vb, old_vb, sizeof (VBlock));

        // update buf->vb in all buffers of this VB to new VB address
        buflist_update_vb_addr_change (vb, old_vb);
        vb->data_type_alloced = alloc_dt;

        FREE_ALIGN_64B (old_vb);
    }

    vb->data_type = dt;
    return vb;
}

// used to change segconf VB data_type is seg_initiatlize (FASTA->FASTQ)
void vb_change_datatype_nonpool_vb (VBlockP *vb_p, DataType new_dt)
{
    *vb_p = vb_update_data_type (*vb_p, new_dt, new_dt, get_vb_size (new_dt));

    pools[POOL_NONPOOL]->vb[(*vb_p)->id.index] = *vb_p;
}

// allocate an unused vb from the pool. separate pools for zip and unzip
VBlockP vb_get_vb (VBlockPoolType type, Task task, VBIType vblock_i, CompIType comp_i)
{
    START_TIMER;

    VBlockPoolP pool = vb_get_pool (type, HARD_FAIL);

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
    VBlockP vb;
    VBIDIndex index; for (index=0; ; index = (index+1) % pool->num_vbs) {
        if (!pool->vb[index]) { // VB is not allocated - allocate it
            vb = pool->vb[index] = CALLOC_ALIGN_64B (sizeof_vb);
            pool->num_allocated_vbs++;
            vb->id = (VBID){ .index = index, .pool = type };
            vb->data_type_alloced = alloc_dt;

            vb->buffer_list.nameר = ר("buffer_list");
            buf_init_lock (&vb->buffer_list);
            buflist_add_buf (vb, &vb->buffer_list, THIS_CODE_LINE);
            vb->buffer_list.vb = vb; // indication buffer was added to buffer list
            break;
        }

        else if (!is_in_use (pool->vb[index])) {
            vb = pool->vb[index] = vb_update_data_type (pool->vb[index], dt, alloc_dt, sizeof_vb); // possibly realloc
            break;
        }

        // case: we've checked all the VBs and none is available - wait a bit and check again
        // in PIZ, this happens when a lot VBs are handed over to the writer thread which has not processed them yet.
        // for example, if writer is blocking on write(), waiting for a pipe counterpart to read.
        // it will be released when the writer thread completes one VB.
        if (index == pool->num_vbs-1) usleep (50000); // 50 ms
    }

    // Logic: num_in_use is always AT LEAST sum(vb)->in_use. i.e. pessimistic. (it can be mometarily more between these two updates)
    uint32_t num_in_use = __atomic_add_fetch (&pool->num_in_use, 1, __ATOMIC_ACQ_REL); // atomic to prevent concurrent update by writer thread and main thread (must be before update of in_use)
    set_in_use (vb, true);

    // initialize VB fields that need to be a value other than 0
    vb->data_type         = dt;
    vb->vblock_i          = vblock_i;
    vb->comp_i            = comp_i;
    vb->compute_thread_id = THREAD_ID_NONE;
    vb->compute_task      = task;
    ca_init_d2d_map (&vb->ca);
    
    if (flag_is_show_vblocks (task)) 
        iprintf ("VB_GET_VB(task=%s id=%u) vb_i=%s/%d num_in_use=%u/%u%s\n", 
                  task_name (task), vb->id.index, comp_name (vb->comp_i), vb->vblock_i, num_in_use, pool->num_vbs,
                  flag.preprocessing ? " preprocessing" : "");

    threads_log_by_vb (vb, task_name (task), "GET VB", 0);

    if (flag.debug_memory)
        iprintf ("vb_get_vb: got vb_i=%d id=%d task=%s dt=%s address=[%p - %p]\n", 
                 vblock_i, index, task_name (task), dt_name(alloc_dt), vb, (char*)vb + sizeof_vb);

    COPY_TIMER_EVB (vb_get_vb);
    return vb;
}

uint32_t vb_pool_get_num_in_use (VBlockPoolType type, VBID *id/*optional out*/)
{
    VBlockPool *pool = vb_get_pool (type, HARD_FAIL);
    int num_in_use = load_acquire (pool->num_in_use); // atomic, bc for POOL_MAIN, writer thread might update concurrently.

    if (id) {
        *id = (VBID){ .pool = NO_POOL };
        if (num_in_use) {
            for (VBIDIndex index=0; index < pool->num_vbs; index++)
                if (pool->vb[index] && pool->vb[index]->in_use) { // not thread safe!
                    *id = (VBID){ .pool = type, .index = index };
                    goto done;
                }

            // all lost use while we were checking
            num_in_use = 0;
        }
    }
done:
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
    ASSERTNOTNULL(vb);

    return IN_RANGE(vb->id.pool, 0, NUM_POOLS) &&
           IN_RANGE(vb->id.index, 0, pools[vb->id.pool]->num_vbs) &&
           vb == pools[vb->id.pool]->vb[vb->id.index];
}

// frees memory of all VBs in a pool and the pool itself
void vb_destroy_pool (VBlockPoolType type, bool destroy_pool)
{
    if (!pools[type]) return;

    for (VBIDIndex index=0; index < pools[type]->num_vbs; index++) 
        vb_destroy_vb (&pools[type]->vb[index]);

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
    store_release (vb->is_processed, (bool)true); 
}

bool vb_is_processed (VBlockP vb)
{
    return load_acquire (vb->is_processed);
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

    iprint_newline();
}

