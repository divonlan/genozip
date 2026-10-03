// ------------------------------------------------------------------
//   buffer.c
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#ifndef _WIN32
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#endif
#ifdef _WIN32
#include <windows.h>
#elif defined __linux__
#include <malloc.h>
#elif defined __APPLE__
#include <malloc/malloc.h>
#endif
#include <math.h>
#include <fcntl.h> 
#include "profiler.h"
#include "buf_struct.h"
#include "buf_list.h"
#include "bits.h"
#include "file.h"
#include "threads.h"
#include "arch.h"

#define DISPLAY_ALLOCS_AFTER 0 // display allocations, except the first X allocations. reallocs are always displayed

// We mark a buffer list entry as removed, by setting its LSb. This keeps the buffer list sorted, in case it is sorted.
// Normally the LSb is 0 as buffers are word-aligned (verified in buflist_add_buf)
#define BL_SET_REMOVED(bl_ent) bl_ent = ((BufferP)((uint64_t)(bl_ent) | 1))
#define BL_IS_REMOVED(bl_ent)  ((uint64_t)(bl_ent) & 1) 

𝓌𝒾𝓃 (static HANDLE heap;)

void buf_increment_user_count (BufferP buf)
{
    ASSERT (buf_user_count (buf) < 0xffff, "user_count at max for buf=%s", buf_desc(buf).s);
    BOLCOUNTER(buf)++;
}

uint16_t buf_decrement_user_count (BufferP buf)
{
    ASSERT (buf_user_count (buf) > 0, "user_count is 0 for buf=%s", buf_desc(buf).s);
    return --BOLCOUNTER(buf);
}

#define buf_unlock_and_decrement_lock_links                                             \
    if (spinlock) {                                                                     \
        ASSERT ((buf)->spinlock[0]->link_count, "spinlock->link_count is 0 for %s", unר(buf->nameר)); \
        if (decrement_relaxed (spinlock->link_count) == 0) {                            \
            FREE ((buf)->spinlock[0]);                                                     \
            (buf)->shared = false;                                                      \
        }                                                                               \
        else                                                                            \
            __atomic_clear (&spinlock->lock, __ATOMIC_RELEASE);                         \
        spinlock = NULL;                                                                \
    }

#define reset_memory_pointer(buf)                                                       \
    /* careful to make sure that instructions are not re-ordered in either direction */ \
    char *old_memory UNUSED = (buf)->memory;                           \
    store_release ((buf)->memory, (void*)0);                                            \
    __atomic_thread_fence (__ATOMIC_ACQ_REL); 

void buf_initialize()
{
#ifdef _WIN32
    heap = GetProcessHeap();
    ASSERT (heap, "GetProcessHeap failed: %s", str_win_error());
#endif
}

rom buf_type_name (ConstBufferP buf)
{
    if (IN_RANGE (buf->type, 0, BUF_NUM_TYPES)) 
        return (rom[])BUFTYPE_NAMES[buf->type];

    else {
        char *s = malloc (32); // used for error printing
        snprintf (s, 32, "invalid_buf_type=%u", buf->type);
        return s;
    }
}

const StrText1K buf_desc (ConstBufferP buf)
{
    #define IS_CONTEXT(obj) ((rom)buf >= (rom)obj->ca.contexts && (rom)buf < (rom)&obj->ca.contexts[MAX_DICTS])
    #define TAG_NAME(obj) obj->ca.contexts[((rom)buf - (rom)obj->ca.contexts) / (sizeof (obj->ca.contexts) / MAX_DICTS /*note: may be different than sizeof(Context) due to word alignment*/)].tag_name

    if (!buf) return (StrText1K){ .s = "NULL" };

    // case: buffer is one of the Buffers within a Context - show tag_name
    rom tag_name = NULL;
    VBlockP vb = buf->vb; 
    if (vb && vb != evb && IS_CONTEXT(vb))
        tag_name = TAG_NAME (vb);
    else if (vb && vb == evb && IS_CONTEXT(z_file))
        tag_name = TAG_NAME (z_file);

    StrText1K desc; // use static memory instead of malloc since we could be in the midst of a memory issue when this is called
    snprintf (desc.s, sizeof (desc.s), "{ \"%s\"%.20s memory=%p data=%p param=%"PRId64"(0x%016"PRIx64") len=%"PRIu64" size=%"PRId64" type=%s shared=%s%.16s promiscuous=%s spinlock=%p%.20s%.20s allocated in %s:%u%.20s }", 
             unר(buf->nameר), cond_str (tag_name, " ctx=", tag_name), 
             buf->memory, buf->data, buf->param, buf->param, buf->len, (uint64_t)buf->size, buf_type_name (buf), 
             TF(buf->shared), cond_int (buf->memory && vb, " users=", BOLCOUNTER(buf)), TF(buf->promiscuous),
             buf->spinlock[0], // note: access spinlock fields only if vb is set. If vb=NULL, buf was destroyed, and spinlock should be NULL - if its not, its likely a memory corruption
             cond_int (buf->spinlock[0] && vb, " locked=", buf->spinlock[0]->lock), cond_int (buf->spinlock[0] && vb, " link_count=", buf->spinlock[0]->link_count),
             unר(buf->funcר), buf->code_line, cond_int (vb, " by vb=", vb->vblock_i));
    return desc;
}

// show underflow, overflow, user count and padding size
const StrText1K buf_desc_CTL (ConstBufferP buf)
{
    if (!buf->memory) return (StrText1K){ "Unallocated" };

    StrText1K s;
    uint64_t overflow  = BOVERFLOW(buf);
    uint64_t underflow = BUNDERFLOW(buf);
    uint8_t users      = BOLCOUNTER(buf);
    int padding        = buf->data ? (int)(buf->memory - buf->data - sizeof(uint64_t)) : 0;

    snprintf (s.s, sizeof(s), "UNDERFLOW=\"%.8s\" (%016"PRIx64") OVERFLOW=\"%.8s\" (%016"PRIx64") users=%u padding=%d",
              (rom)&underflow, underflow, (rom)&overflow, overflow, users, padding);

    return s;
}

// quick inline for internal buf_struct.c use check overflow and underflow in an allocated buffer
static inline void no_integrity (ConstBufferP buf, Caller caller, rom buf_func)
{
    flag.quiet = false;

    ASSERTW (BUNDERFLOW(buf) == UNDERFLOW_TRAP, _ERR "called from %s:%u to %s %s%s: Error in %s: buffer has corrupt underflow trap: %s",
             CALLERf, buf_func, version_str().s, license_get_number().s, buf_desc(buf).s, str_to_printable_(buf->memory, 8).s);

    ASSERTW (BOVERFLOW(buf) == OVERFLOW_TRAP, _ERR "called from %s:%u to %s %s%s: Error in %s: buffer has corrupt overflow trap: %s",
            CALLERf, buf_func, version_str().s, license_get_number().s, buf_desc(buf).s, str_to_printable_(buf->memory + buf->size + sizeof(uint64_t), 8).s);
    
    bool corruption_detected = buflist_test_overflows (buf->vb, buf_func);
    if (corruption_detected) buflist_test_overflows_all_other_vb (buf->vb, buf_func, true, false); // corruption not from this VB - test the others
    exit_on_error (true);
}

// quick inline wrapper
static inline void buf_verify_integrity (ConstBufferP buf, Caller caller, rom buf_func)
{
    if (buf->memory && buf->type == BUF_REGULAR && 
        (BUNDERFLOW(buf) != UNDERFLOW_TRAP || BOVERFLOW(buf) != OVERFLOW_TRAP))
        
        no_integrity (buf, caller, buf_func);
}

// called from other modules for debugging memory issues. 
void buf_verify_do (ConstBufferP buf, rom msg, Caller caller)
{
    if (!buf || !buf->memory) return;

    BufferSpinlockP spinlock = buf->promiscuous ? buf_lock_promiscuous (buf, caller) : NULL; // prevent frees or reallocs while we're testing
    if (buf->promiscuous && !spinlock) return; // by the time we acquired the lock, buf was already freed

    buf_verify_integrity (buf, THIS_CODE_LINE, msg);

    buf_unlock;
}

static void buf_reset (BufferP buf)
{
    // first set "memory" to 0 - so buf_test_overflow doesn't test this buffer 
    reset_memory_pointer (buf);

    VBlockP save_vb    = buf->vb;        // preserve vb because still in vb->buffer_list
    int32_t save_funcר = buf->funcר;     // preserve func and code_line for buf_list* error reporting
    int save_line      = buf->code_line;

    memset (buf, 0, sizeof (Buffer)); // make this buffer UNALLOCATED

    buf->vb        = save_vb;
    buf->funcר     = save_funcר;
    buf->code_line = save_line;
}

static void buf_init (BufferP buf, char *memory, uint64_t size, bool aligned_64B, Caller caller, rom name)
{
    // set some parameters before allocation so they can go into the error message in case of failure
    buf->funcר     = caller.funcר;
    buf->code_line = caller.code_line;

    if (name) 
        buf->nameר = ר(name);
    else
        ASSERT (buf->nameר, "buffer has no name. func=%s:%u", unר(buf->funcר), buf->code_line);

    if (!memory) { // malloc or realloc failed
        buflist_show_memory (true, 0, 0);

        ABORT ("%s: Out of memory%s. Details: %s:%u failed to allocate %s bytes. Buffer: %s", 
               global_cmd, 
               cond_int (IS_ZIP, ". Try running with a lower vblock size using --vblock. Current vblock size is: ", segconf.vb_size >> 20),
               CALLERf, str_int_commas (size + CTL_SIZE).s, buf_desc(buf).s);
    }

    buf->data = memory + sizeof (uint64_t);
    buf->size = size;

    // data needs to be 64B aligned
    if (aligned_64B) 
        buf_aligned_64B (buf); // adjusts data and size

    BUNDERFLOW_(memory) = UNDERFLOW_TRAP; // underflow protection
    BOVERFLOW_(memory, buf->data, buf->size) = OVERFLOW_TRAP;  // overflow prortection (underflow protection was copied with realloc)
    BOLCOUNTER_(memory, buf->data, buf->size) = 1;  // 1 when memory is first allocated

    // only when we're done initializing - we update memory - that buf_test_overflow running concurrently doesn't test
    // half-initialized buffers
    store_release (buf->memory, memory); 
}

// allocates or enlarges buffer
// if it needs to enlarge a buffer fully overlaid by an overlay buffer - it abandons its memory (leaving it to
// the overlaid buffer) and allocates new memory
void buf_alloc_do (VBlockP vb, BufferP buf, uint64_t requested_size,
                   float grow_at_least_factor, // IF we need to allocate or reallocate physical memory, we get this much more than requested
                   bool aligned_64B, // if true, buf->data will be aligned to 64B (cache_line).  
                   rom name, Caller caller)      
{
    START_TIMER; // don't account time for these ^ calls - we're interested in actual allocations 

    if (!vb) vb = buf->vb;
    ASSERT (vb, "called from %s:%u: null vb", CALLERf);

    // **** sanity checks ****
    ASSERT ((int64_t)requested_size > 0, "called from %s:%u: negative requested_size=%"PRId64" for name=%s", CALLERf, requested_size, name);

#define REQUEST_TOO_BIG_THREADSHOLD (3 GB)
    if (requested_size > REQUEST_TOO_BIG_THREADSHOLD && !buf->can_be_big) // use WARN instead of ASSERTW to have a place for breakpoint
        WARN (_WRN "buf_alloc called from %s:%u %s for \"%s\" requested %s. This is suspiciously high and might indicate a bug %s. z_dt=%s vb->vblock_i=%u buf=%s line_i=%d vb_size=%s RAM=%3.1lf GB txt_file->disk_size=%s",
              CALLERf, version_str().s, name ? name : unר(buf->nameר), str_size (requested_size).s, report_support(), z_dt_name(), vb->vblock_i, buf_desc (buf).s, vb->line_i, str_size (segconf.vb_size).s, arch_get_physical_mem_size(), txt_file ? str_size (txt_file->disk_size).s : "N/A");

    ASSERT (buf->type == BUF_REGULAR || buf->type == BUF_UNALLOCATED, "called from %s:%u: cannot buf_alloc a buffer of type %s. details: %s", 
            CALLERf, buf_type_name (buf), buf_desc (buf).s);

    // if this happens: either 1. the wrong VB was given now, or when initially allocating this buffer OR
    // 2. VB was REALLOCed in vb_get_vb, but for some reason this buf->vb was not updated because it was not on the buffer list
    ASSERT (!buf->vb || vb == buf->vb, "called from %s:%u: buffer=%p has wrong VB: vb=%p (id=%s vblock_i=%u) but buf->vb=%p", 
            CALLERf, buf, vb, dis_vb_id(vb->id).s, vb->vblock_i, buf->vb);

    ASSERT (vb != evb || buf->promiscuous || threads_am_i_main_thread(), "called from %s:%u: A non-main thread is attempting to allocate an evb buffer \"%s\" with promiscuous=false", 
            CALLERf, name ? name : unר(buf->nameר));

    // **** initial memory allocation ****

    // CASE 2: initial allocation: exactly at requested size
    if (!buf->memory) {
        // round up to 64 bit boundary to avoid aliasing errors with the overflow indicator
        // for data to be 64B aligned (if requested) the worst-case scenario is if "memory" is 64B-aligned and we need to add 56 bytes after UNDERFLOW to get to the next 64B
        // (given Linux, Mac and Windows allocate on 16B-boundary, verified in arch_initialize)
        uint64_t new_size = MIN_(MAX_BUFFER_SIZE, ROUNDUP8(requested_size + (aligned_64B ? 56 : 0)));

        ASSERT (new_size >= requested_size, "called from %s:%u: requested too much memory=%s for buf=%s. vb->vblock_i=%u", 
                CALLERf, buf_desc(buf).s, str_int_commas (requested_size).s, vb->vblock_i); 

        char *memory = (char *)buf_low_level_malloc (new_size + CTL_SIZE, false, false, caller);
        buf->type = BUF_REGULAR;

        buf_init (buf, memory, new_size, aligned_64B, caller, name);
        
        if (buf != &vb->buffer_list) { // buffer_list buffer is added in vb_get_vb / vb_initialize_nonpool_vb
            if (!buf->promiscuous) // if promiscuous or buffer_list, already added
                buflist_add_buf (vb, buf, caller);
            else
                ASSERT (buf->vb, "called from %s:%u: Expecting promiscuous buffer to be on the buffer_list: %s", CALLERf, buf_desc (buf).s);
        }
        
        goto done;
    }

    // **** realloc: calculate size include "growth" ****

    // add an epsilon to avoid floating point multiplication ending up slightly less than the integer
    grow_at_least_factor = MAX_(1.0001, grow_at_least_factor); 

    // grow us requested - rounding up to 64 bit boundary to avoid aliasing errors with the overflow indicator
    uint64_t new_size = MIN_(MAX_BUFFER_SIZE, ROUNDUP8((uint64_t)(requested_size * grow_at_least_factor) + (aligned_64B ? 56 : 0)));

    ASSERT (new_size >= requested_size, "called from %s:%u: allocated too little memory for buffer %s: requested=%s, allocated=%s. vb->vblock_i=%u", 
            CALLERf, buf_desc (buf).s, str_int_commas (requested_size).s, str_int_commas (new_size).s, vb->vblock_i); // floating point paranoia

    // CASE 3: realloc of a non-shared - use realloc that will extend instead of malloc & copy if possible 
    if (!buf->shared) {    
        buf_lock_if (buf, buf->spinlock[0]); // promiscous or buf_list or caller lock...

        buf_verify_integrity (buf, caller, "buf_alloc_do");

        reset_memory_pointer (buf);
         
        char *new_memory = (char *)buf_low_level_realloc (old_memory, new_size + CTL_SIZE, name, caller);
        buf_init (buf, new_memory, new_size, aligned_64B, caller, name);
    
        buf_unlock;
    }

    else { // shared cases
        ASSERT (!aligned_64B, "aligned_64B not supported for overlaying: %s", name);

        buf_lock (buf);

        buf_verify_integrity (buf, caller, "buf_alloc_do(shared)");

        // CASE 4: shared: currently no overlayers and standard "data" value - we can realloc
        if (buf_user_count (buf) == 1 && 
            (!buf->data || buf->data == buf->memory + sizeof(uint64_t))) {
            reset_memory_pointer (buf);

            char *new_memory = (char *)buf_low_level_realloc (old_memory, new_size + CTL_SIZE, name, caller);
            buf_init (buf, new_memory, new_size, false, caller, name);
        }

        // CASE 5: shared: have overlayers, or non-standard data value - malloc & copy
        else {
            char *new_memory = (char *)buf_low_level_malloc (new_size + CTL_SIZE, false, false, caller);

            memcpy (new_memory + sizeof (uint64_t), buf->data, buf->size); // copy old data
            uint16_t user_count = buf_decrement_user_count (buf); // we are no longer using the old memory

            reset_memory_pointer (buf);
            
            if (!user_count)
                buf_low_level_free (old_memory, false, caller);

            buf_init (buf, new_memory, new_size, false, caller, name);
        }

        buf_unlock; // note: spinlock stays the same even if memory is realloced
    }

done:
    if (flag.debug_memory && !buflist_locate (buf, NULL)) 
        // not in any VB, file, reference or gencomp structure
        iprintf ("buf_alloc_do: allocated independent buf %p: %s\n", buf, buf_desc(buf).s);

    if (vb == evb) COPY_TIMER_EVB (buf_alloc_main); // works even for promiscuous bc uses atomic 
    else           COPY_TIMER (buf_alloc_compute); 
}

// shrink buffer down to size, returning memory to libc (but not to kernel)
// main thread only, but just because of the COPY_TIMER_EVB
void buf_trim_do (BufferP buf, uint64_t size, Caller caller)
{
    START_TIMER;

    size = ROUNDUP8 (size); // buffer size is required to be a multiple of 8 

    if (size >= buf->size || !buf->memory) return; // nothing to do - size if already smaller

    ASSERT (!buf->shared, "%s:%u trimming is not currently supported on shared buffers: buf=%s", CALLERf, buf_desc(buf).s);

    buf_lock_if (buf, buf->spinlock[0]); // promiscous or buf_list or caller lock...

    buf_verify_integrity (buf, caller, "buf_alloc_do");

    reset_memory_pointer (buf);
        
    char *new_memory = (char *)buf_low_level_realloc (old_memory, size + CTL_SIZE, unר(buf->nameר), caller);
    buf_init (buf, new_memory, size, false, caller, unר(buf->nameר));

    buf_unlock;

    COPY_TIMER_EVB (buf_trim_do);
}

void buf_set_shared (BufferP buf)
{
    if (!buf->shared) {
        buf_init_lock (buf);
        buf->shared = true;
    }
}

void buf_remove_spinlock (BufferP buf)
{
    if (!buf->spinlock[0]) return;

    ASSERT (!buf->memory, "cannot remove spinlock: buf has memory: %s", buf_desc(buf).s);

    if (!decrement_relaxed (buf->spinlock[0]->link_count))
        FREE (buf->spinlock[0]); // we were the only uses of the spinlock

    buf->shared = buf->promiscuous = false;
}

// an overlay buffer is a buffer using some of the memory of another buffer - it doesn't have its own memory
void buf_overlay_do (VBlockP vb, 
                     BufferP top_buf, // dst 
                     BufferP bottom_buf, 
                     Caller caller, rom name)
{   
    START_TIMER;

    // if this buffer was used by a previous VB as a bottom buffer - we need to "destroy" it first
    if (top_buf->type == BUF_REGULAR && top_buf->data == NULL && 
        (top_buf->memory || (top_buf->spinlock[0] && top_buf->spinlock[0] != bottom_buf->spinlock[0]))) 
        buf_destroy (*top_buf);

    ASSERT (top_buf->type == BUF_UNALLOCATED, "%s: Call from %s:%u: cannot buf_overlay to a buffer %s already in use", VB_NAME, CALLERf, buf_desc (top_buf).s);

    // overlaying a SHM buffer, just creates another SHM buffer
    if (bottom_buf->type == BUF_SHM) {
        buf_attach_to_shm_do (vb, top_buf, 
                              bottom_buf->memory, bottom_buf->size, 0, 
                              caller, name);

        top_buf->len = bottom_buf->len;
        return;
    }

    ASSERT (bottom_buf->shared && bottom_buf->spinlock[0], 
            "%s: Call from %s:%u: expecting bottom_buf %s to have a spinlock and shared=true", VB_NAME, CALLERf, buf_desc (bottom_buf).s);

    ASSERT (bottom_buf->type == BUF_REGULAR,
            "%s: Call from %s:%u: bottom_buf %s in buf_overlay must be a bottom or shm buffer", VB_NAME, CALLERf, buf_desc (bottom_buf).s);

    top_buf->type      = BUF_REGULAR;
    top_buf->nameר     = ר(name);
    top_buf->len       = bottom_buf->len;
    top_buf->funcר     = caller.funcר;
    top_buf->code_line = caller.code_line;
    top_buf->shared    = true;
    top_buf->spinlock[0]  = bottom_buf->spinlock[0];

    // note: while we're waiting for the lock, bottom_buf may be reallocing - but this doesn't change spin_lock
    buf_lock (bottom_buf); // locking spinlock shared between top, bottom buffers and all other overlayers

    buf_verify_integrity (bottom_buf, caller, "buf_overlay_do");

    // note: data+size MUST be at the control region, as we have the overlay counter there
    top_buf->size = bottom_buf->size;
    top_buf->data = bottom_buf->data;

    // increment spinlock users and memory users (note: spinlock link cou t >= memory_users bc a spinlock can be used by multiple memories in case of a realloc-induced split)
    increment_relaxed (bottom_buf->spinlock[0]->link_count);    

    buf_increment_user_count (bottom_buf);
    
    // final step
    store_release (top_buf->memory, bottom_buf->memory); 

    buf_unlock;

    buflist_add_buf (vb, top_buf, caller); 

    COPY_TIMER (buf_overlay_do);
}

// note: top_buf struct is not assumed to be initialized - it is fully overwritten
// note: caller guarantees that top_buf is NOT in buffer_list: must buf_destroy before calling if it might be
void buf_superimpose_do (VBlockP vb, BufferP top_buf, BufferP bottom_buf, uint64_t start_in_bottom, Caller caller, rom name)
{
    ASSERT (bottom_buf->type == BUF_REGULAR || bottom_buf->type == BUF_SHM,
            "%s: Call from %s:%u: bottom_buf %s in buf_superimpose must be a bottom or shm buffer", VB_NAME, CALLERf, buf_desc (bottom_buf).s);

    ASSERT (start_in_bottom < bottom_buf->size, 
            "called from %s:%u: not enough room in bottom buffer for superimposed buf: start_in_bottom=%"PRIu64" but bottom_buf.size=%"PRIu64,
            CALLERf, start_in_bottom, (uint64_t)bottom_buf->size);
            
    *top_buf = (Buffer){
        .type      = BUF_SUPERIMPOSED,
        .nameר     = ר(name),
        .funcר     = caller.funcר,
        .code_line = caller.code_line,
        .size      = bottom_buf->size - start_in_bottom,
        .data      = bottom_buf->data + start_in_bottom,
        .len       = start_in_bottom ? 0 : bottom_buf->len, // copy len if superimposing the entire buffer
        .memory    = NULL
    };
}

void buf_attach_to_shm_do (VBlockP vb, BufferP buf, void *memory, uint64_t size, uint64_t start, Caller caller, rom name)
{
    // preserve len and param (i.e. nbits and nwords)
    uint64_t save_len   = buf->len;
    uint64_t save_param = buf->param;

    // if this buffer was used by a previous VB as a regular buffer - we need to "destroy" it first
    if (buf->vb) 
        buf_destroy (*buf);

    *buf = (Buffer){
        .type      = BUF_SHM,
        .nameר     = ר(name),
        .funcר     = caller.funcר,
        .code_line = caller.code_line,
        .vb        = vb,
        .memory    = memory,
        .data      = memory + start, // note: no control area in shm buffers
        .size      = size,
        .len       = save_len,
        .param     = save_param
    };
}

void buf_free_do (BufferP buf, Caller caller) 
{
    switch (buf->type) {

        case BUF_REGULAR: {
            START_TIMER;
        
            buf_lock_if (buf, buf->spinlock[0]);

            ASSERT (!buf->spinlock[0] || buf->spinlock[0]->link_count == 1, "Cannot buf_free an overlaid buffer: use buf_destroy: %s", buf_desc (buf).s);

            buf_verify_integrity (buf, caller, "buf_free_do");

            uint16_t user_count = buf_user_count (buf);
            ASSERT (user_count, "%s:%u: user_count=0 in buffer %s", CALLERf, buf_desc(buf).s);

            // case: if user_count >= 2: reset buffer and leave memory to other user 
            // (often: compute thread resets vb buffer and leaves memory to main thread z_file buffer)
            if (user_count >= 2) {
                user_count = buf_decrement_user_count (buf);

                reset_memory_pointer (buf);
                buf->size = 0;
                buf->nameר = buf->funcר = 0;
                buf->code_line = 0;
                buf->type = BUF_UNALLOCATED;
                // preserved: vb, promiscuous, shared, spinlock + still on its VB's buffer_list
            }

            // if last remaining user is a partial overlay or 64B alignment - increase size to be based on entire memory
            else if (buf->data && (buf->data - buf->memory != sizeof (uint64_t)))
                buf->size += (buf->data - buf->memory - sizeof (uint64_t));

            buf->data        = NULL; 
            buf->can_be_big  = false;
            buf->len         = 0;
            buf->param       = 0;

            buf_unlock;

#ifdef PROFILE
            if (buf->vb == evb) COPY_TIMER_EVB (buf_free_main); // works even for promiscuous bc uses atomic 
            else { VBlockP vb = buf->vb; COPY_TIMER (buf_free_compute); };
#endif
            break;
        }
        case BUF_UNALLOCATED: // reset len and param that may be used even without allocating the buffer
            buf->len         = 0;
            buf->param       = 0;
            break;

        case BUF_SHM:
        case BUF_SUPERIMPOSED:
            *buf = (Buffer){}; // note: SHM and SUPERIMPOSED buffers are not in buffer_list, so we must zero buf->vb
            break;

        default:
            ABORT ("invalid buf->type=%s", buf_type_name (buf));
    }
} 

void buf_destroy_do_do (BufListEnt *ent, Caller caller)
{
    if (!ent) return;
    START_TIMER;

    BufferP buf = ent->buf;
    VBlockP vb  = buf->vb;

    if (flag.debug_memory==1) 
        iprintf ("Destroy %s: buf_addr=%p vb->id=%s buf_i=%u\n", buf_desc (buf).s, buf, dis_vb_id(buf->vb->id).s, BNUM (buf->vb->buffer_list, ent));

    // remove from buffer list
    BL_SET_REMOVED (ent->buf);
    buf->vb = NULL;

    switch (buf->type) {
        case BUF_REGULAR : { 
            buf_lock_if (buf, buf->spinlock[0]);

            if (buf->memory) { // NULL if buffer was disowned
                buf_verify_integrity (buf, caller, "buf_destroy_do");
                uint16_t remaining_user_count = buf_decrement_user_count (buf);

                // first set "memory" to 0 - so buf_test_overflow doesn't test this buffer after destroyed
                reset_memory_pointer (buf);

                if (!remaining_user_count) 
                    buf_low_level_free (old_memory, false, caller); 
            }

            buf_unlock_and_decrement_lock_links; // also frees spinlock if we're the last user
            break;
        }

        case BUF_SHM : 
        case BUF_SUPERIMPOSED:
            buf_free (*buf); 
            break;
        
        case BUF_UNALLOCATED : {
            buf_lock_if (buf, buf->spinlock[0]); // possibly set as promiscuous but never allocated
            buf_unlock_and_decrement_lock_links;
            break;
        }

        default : ABORT ("called from %s:%u: Error in buf_destroy_do: invalid buffer type %s", CALLERf, buf_type_name (buf));
    }

    buf_reset (buf);

    if (vb==evb) COPY_TIMER_EVB (buf_destroy_do_do_main);
    else         COPY_TIMER (buf_destroy_do_do_compute);
}

void buf_destroy_do (BufferP buf, Caller caller)
{
    if (!buf || 
        (!buf->vb && buf->type == BUF_UNALLOCATED) || // never allocated 
        flag.let_OS_cleanup_on_exit) return; // nothing to do (we don't destroy on exit, as the exiting thread may not be able to remove from buf_list)

    BufListEnt *ent;

    if (buf->type == BUF_DISOWNED) {
        buf_low_level_free (buf->memory, false, caller);
        *buf = (Buffer){};  
    }

    else if ((ent = buflist_find_buf (buf->vb, buf, SOFT_FAIL))) 
        buf_destroy_do_do (ent, caller);
    
    else
        *buf = (Buffer){};  
}

// similar to buf_move, but also moves buffer between buf_lists. can be run by the main thread only.
// IMPORTANT: only works when called from main thread, when BOTH src and dst VB are in full control of main thread, so that there
// no chance another thread is concurrently modifying the buf_list of the src or dst VBs 
void buf_grab_do (VBlockP dst_vb, BufferP dst_buf, rom dst_name/*optional*/, BufferP src_buf, Caller caller)
{
    ASSERTMAINTHREAD;
    ASSERT (src_buf, "called from %s:%u: buf is NULL", CALLERf);
    if (src_buf->type == BUF_UNALLOCATED) return; // nothing to grab

    ASSERT (src_buf->type == BUF_REGULAR && !src_buf->shared, "called from %s:%u: this function can only be called for a non-shared REGULAR buf", CALLERf);
    ASSERT (dst_buf->type == BUF_UNALLOCATED, "called from %s:%u: expecting dst_buf to be UNALLOCATED", CALLERf);

    reset_memory_pointer (src_buf);

    dst_buf->type     = BUF_REGULAR;
    dst_buf->len      = src_buf->len;
    dst_buf->param    = src_buf->param;
    dst_buf->spinlock[0] = src_buf->spinlock[0];
    buf_init (dst_buf, old_memory, src_buf->size, false, caller, dst_name ? dst_name : unר(src_buf->nameר));

    buflist_add_buf (dst_vb, dst_buf, caller);

    // remove src_buf from buffer list
    buf_lock (&src_buf->vb->buffer_list);
    BL_SET_REMOVED (buflist_find_buf (src_buf->vb, src_buf, HARD_FAIL)->buf);
    buf_unlock;
    
    src_buf->vb = NULL;

    buf_reset (src_buf);
}

// moves all the data between buffers in the same VB, keeping the same buffer_list entry. 
void buf_move_do (VBlockP vb, BufferP dst_buf, rom dst_name/*optional*/, BufferP src_buf, Caller caller)
{
    ASSERT (!dst_buf->data, "%s:%u: dst_buf is not empty: %s", CALLERf, buf_desc(dst_buf).s);
    ASSERT (src_buf->type == BUF_REGULAR, "%s:%u: src_buf is %s", CALLERf, buf_type_name (src_buf));
    ASSERT (src_buf->vb == vb, "%s:%u: src_buf has wrong vb", CALLERf);
    ASSERT (buf_user_count (src_buf) == 1, "%s:%u: expecting src_buf to have 1 user but it has %u", CALLERf, buf_user_count (src_buf));
    buf_verify_integrity (src_buf, caller, "buf_move");

    if (dst_buf->vb) buf_destroy (*dst_buf); // also remove from buffer_list

    if (!dst_name) dst_name = unר(src_buf->nameר);

    *dst_buf = (Buffer){ .funcר      = caller.funcר,
                         .code_line  = caller.code_line,
                         .nameר      = ר(dst_name),
                         .data       = src_buf->data,
                         .size       = src_buf->size,
                         .param      = src_buf->param,
                         .len        = src_buf->len,
                         .type       = BUF_REGULAR,
                         .can_be_big = src_buf->can_be_big,
                         .vb         = vb,
                         .memory     = src_buf->memory }; // no need for atomic_store, bc buflist_move_buf unlocks with ATOMIC_RELEASE

    // make the buffer_list entry of src_buf now point to dst_buf     
    buflist_move_buf (vb, dst_buf, dst_name, src_buf, caller);

    // promiscuous, shared and spinlock are NOT moved. 
    buf_remove_spinlock (src_buf);

    reset_memory_pointer (src_buf);
    *src_buf = (Buffer){};
}

// move buffer struct to a new location, without adding it the buffer list 
void buf_disown_do (VBlockP vb, BufferP src_buf, BufferP dst_buf, bool make_a_copy, Caller caller)
{
    ASSERT (vb == src_buf->vb, "called from %s:%u: buf->vb mismatches vb. vb->vblock_i=%u", CALLERf, vb->vblock_i);
    ASSERT (src_buf->type == BUF_REGULAR, "called from %s:%u: expecting src_buf to be BUF_REGULAR, buf_type=%s", CALLERf, buf_type_name(src_buf));
    ASSERT (!src_buf->promiscuous, "called from %s:%u: src_buf cannot be promiscuous", CALLERf);
    ASSERT (!src_buf->shared, "called from %s:%u: src_buf cannot be shared", CALLERf);
    ASSERT (dst_buf->type == BUF_UNALLOCATED, "called from %s:%u: expecting dst_buf to be BUF_UNALLOCATED, buf_type=%s", CALLERf, buf_type_name(dst_buf));

    buf_verify_integrity (src_buf, caller, "buf_disown_do");

    *dst_buf = (Buffer) {
        .memory     = src_buf->memory,
        .data       = src_buf->data,
        .len        = src_buf->len,
        .param      = src_buf->param,
        .can_be_big = src_buf->can_be_big,
        .nameר      = src_buf->nameר,
        .size       = src_buf->size,
        .type       = BUF_DISOWNED,
        .funcר      = caller.funcר,
        .code_line  = caller.code_line };
    
    // dst is a disowned copy of src
    if (make_a_copy) {
        dst_buf->memory = MALLOC (dst_buf->size);
        dst_buf->data   = src_buf->data ? (dst_buf->memory + (src_buf->data - src_buf->memory)) : 0;
        memcpy (dst_buf->memory, src_buf->memory, src_buf->size);
    }

    // src is moved to dst and disowned
    else {
        dst_buf->memory = src_buf->memory;
        dst_buf->data   = src_buf->data;

        // first set "memory" to 0 - so buf_test_overflow doesn't test this buffer after destroyed
        reset_memory_pointer (src_buf);
        buf_destroy (*src_buf);
    }
}

// removes data and len from buffer, returning them to caller, and replenishes the buffer memory
void buf_extract_data_do (BufferP buf, char **data_p, uint64_t *len_p, uint32_t *len32_p, char **memory_p, Caller caller)
{
    ASSERT (buf->type == BUF_REGULAR, "called from %s:%u: expecting buf to be BUF_REGULAR, buf_type=%s", CALLERf, buf_type_name(buf));

    if (data_p)  *data_p  = buf->data;
    if (len_p)   *len_p   = buf->len;
    if (len32_p) *len32_p = buf->len32;

    *memory_p = buf->memory; // caller must free this
    
    char *new_memory = (char *)buf_low_level_malloc (buf->size + CTL_SIZE, false, false, caller);
    buf_init (buf, new_memory, buf->size, false, caller, NULL);

    buf->len = 0;
}

//-------------------------
// low-level functions
//-------------------------

void buf_low_level_free (void *p, bool aligned_64B, Caller caller)
{
    if (!p) return; // nothing to do

    START_TIMER;
    bool p_is_evb = (p == evb);

    if (flag.debug_memory==1) 
        iprintf ("Memory freed by free(): %p %s:%u\n", p, CALLERf);

#ifndef _WIN32
    free (p);
#else
    if (aligned_64B)
        _aligned_free (p);
    else
        ASSERT (HeapFree (heap, 0, p), "HeapFree failed: %s", str_win_error());
#endif

    if (!p_is_evb) // unless we just freed evb...
        COPY_TIMER_EVB (buf_low_level_free);
}

static StrText1K oom_tip (void)
{
    StrText1K s;
    snprintf (s.s, sizeof(s), "\n" _TIP "use --low-memory or alternatively limit the number of concurrent threads with --threads "
             "(currently %d - affects speed) and/or reduce the amount of data processed by each thread with --vblock "
             "(currently %d - affects compression ratio)", global_max_threads, (int)(segconf.vb_size / (1 MB)));

    return s;
} 

void *buf_low_level_realloc (void *p, size_t size, rom name, Caller caller)
{
    void *new = X𝓌𝒾𝓃 (realloc (p, size))
                𝓌𝒾𝓃  (HeapReAlloc (heap, 0, p, size));

    ASSERTW (new, _ERR "Out of memory in %s:%u: realloc failed (name=%s size=%"PRIu64" bytes). %s", 
             CALLERf, name, (uint64_t)size, IS_ZIP ? oom_tip().s : "");

    if (flag.debug_memory && size >= flag.debug_memory) {
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wpragmas"         // avoid warning if "-Wuse-after-free" is not defined in this version of gcc
#pragma GCC diagnostic ignored "-Wunknown-warning-option" // same
#pragma GCC diagnostic ignored "-Wuse-after-free"  // avoid compiler warning of using p after it is freed
        iprintf ("realloc(): old=%p new=%p name=%s size=%"PRIu64" %s:%u\n", p, new, name, (uint64_t)size, CALLERf);
#pragma GCC diagnostic pop
    }

    return new;
}

void *buf_low_level_malloc (size_t size, bool zero, bool aligned_64B, Caller caller)
{
    void *new;

    if (aligned_64B) {
        size = ROUNDUP64 (size);
        new = X𝓌𝒾𝓃(aligned_alloc (64, size))
              𝓌𝒾𝓃(_aligned_malloc (size, 64));
        if (zero) memset (new, 0, size);
    }

    else
        new = X𝓌𝒾𝓃 (zero ? calloc (size, 1) : malloc (size)) 
              𝓌𝒾𝓃 (HeapAlloc (heap, zero ? HEAP_ZERO_MEMORY : 0, size));

    ASSERT (new, "Out of memory in %s:%u: malloc failed (size=%"PRIu64" bytes). %s", 
            CALLERf, (uint64_t)size, IS_ZIP ? oom_tip().s : "");

    if (flag.debug_memory && size >= flag.debug_memory) 
        iprintf ("malloc(): %p size=%"PRIu64" %s:%u\n", new, (uint64_t)size, CALLERf);
    
    return new;
}

void return_freed_memory_to_kernel (void)
{
    ℓ𝒾𝓃𝓊𝓍 (malloc_trim (0);)                      // return whole free pages to the kernel
    𝓂𝒶𝒸 (malloc_zone_pressure_relief (NULL, 0);)  // tell OS that this process is interested in participating in "pressure relief" - freeing memory. OS will decide if and when to actually release the memory.
    𝓌𝒾𝓃 (HeapCompact (heap, 0);)                  // return blocks marked for "deferred free" to the kernel
}

uint64_t buf_mem_size (ConstBufferP buf) 
{ 
    // note: this calculation does not discount for memory overlaid in multiple buffers (with shared=true)
    return buf->type != BUF_REGULAR ? 0
         : !buf->memory             ? 0
         : buf->data                ? ((buf->data - buf->memory) + buf->size + sizeof (uint64_t) + sizeof(uint16_t)) // might be "partial overlay"
         :                            (buf->size + CTL_SIZE); 
}

//---------------------
// Bits stuff
//---------------------

// adds "more_bits" to the bitmap, and optionally allocates memory beyond the end of the bitmap 
// to avoid future allocations. optionally initializes (only) the new bits added.
BitsP buf_alloc_bits_do (VBlockP vb, BufferP buf, uint64_t more_bits, uint64_t preallocate_at_least_bits, BitsInitType init_to, float grow_at_least_factor, rom name, Caller caller)
{
    ASSERT0 (buf->type == BUF_UNALLOCATED || buf->type == BUF_REGULAR, "buf needs to be BUF_UNALLOCATED or BUF_REGULAR");

    uint64_t old_nbits  = buf->nbits;
    uint64_t old_nwords = buf->nwords;
    buf->nbits += more_bits;   
    buf->nwords = roundup_bits2words64 (buf->nbits);
    
    // case: we're adding more words (if not, top bits in unchanged top word are expected to be zero already, per Bits invariant)
    if (buf->nwords != old_nwords) {
        buf_alloc_(vb, buf, buf->nwords - old_nwords, preallocate_at_least_bits / 64, 
                   sizeof(uint64_t), grow_at_least_factor, name, caller);

        // we need to clear the unused bits in the high word, however we can't use bits_clear_excess_bits_in_top_word because
        // it is read-modify-write, in which read, reads uninitialized memory. instead, we just the zero the entire new word
        buf->words[buf->nwords-1] = 0;
    }
    
    // clear / set only added bits
    if (init_to == CLEAR && buf->nbits > old_nbits) 
        bits_clear_region (buf, old_nbits, buf->nbits - old_nbits);           
    
    else if (init_to == SET && buf->nbits > old_nbits) 
        bits_set_region (buf, old_nbits, buf->nbits - old_nbits);           

    DEBUG_VALIDATE_BITS (buf);
    
    return buf;
}

// creates and optionally initializes a bitmap of requested "exact_bits" 
BitsP buf_alloc_bits_exact_do (VBlockP vb, BufferP buf, uint64_t exact_bits, BitsInitType init_to, float grow_at_least_factor, rom name, Caller caller)
{
    ASSERT0 (buf->type == BUF_UNALLOCATED || buf->type == BUF_REGULAR, "buf needs to be BUF_UNALLOCATED or BUF_REGULAR");

    buf->nwords = roundup_bits2words64 (exact_bits);
    buf->nbits  = exact_bits;   
    
    buf_alloc_(vb, buf, 0, buf->nwords, sizeof(uint64_t), grow_at_least_factor, name, caller);

    // clear / set all bits
    if (init_to == CLEAR) 
        memset (buf->data, 0, buf->nwords * sizeof(uint64_t));

    else if (init_to == SET) {
        memset (buf->data, 0xff, buf->nwords * sizeof(uint64_t));
        bits_clear_excess_bits_in_top_word (buf); // zero unused bits
    }
    else
        buf->words[buf->nwords-1] = 0; // zero top word

    return buf;
}

// convert a Buffer from a z_file section whose len is in char to a bits
Bits *buf_zfile_buf_to_bits (BufferP buf, uint64_t nbits)
{
    ASSERT (roundup_bits2bytes (nbits) <= buf->len, "nbits=%"PRId64" indicating a length of at least %"PRId64", but buf->len=%"PRId64,
            nbits, roundup_bits2bytes (nbits), buf->len);

    Bits *bits = buf;
    bits->nbits  = nbits;
    bits->nwords = roundup_bits2words64 (bits->nbits);

    ASSERT (roundup_bits2bytes64 (nbits) <= buf->size, "buffer to small: buf->size=%"PRId64" but bits has %"PRId64" words and hence requires %"PRId64" bytes",
            (uint64_t)buf->size, bits->nwords, bits->nwords * sizeof(uint64_t));

    LTEN_bits (bits);

    bits_clear_excess_bits_in_top_word (bits);

    return bits;
}

