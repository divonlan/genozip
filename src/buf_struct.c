// ------------------------------------------------------------------
//   buffer.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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
#include <fcntl.h> 
#include "genozip.h"
#include "profiler.h"
#include "buf_struct.h"
#include "buf_list.h"
#include "bits.h"
#include "file.h"
#include "threads.h"
#include "version.h"

#define DISPLAY_ALLOCS_AFTER 0 // display allocations, except the first X allocations. reallocs are always displayed

// We mark a buffer list entry as removed, by setting its LSb. This keeps the buffer list sorted, in case it is sorted.
// Normally the LSb is 0 as buffers are word-aligned (verified in buflist_add_buf)
#define BL_SET_REMOVED(bl_ent) bl_ent = ((BufferP)((uint64_t)(bl_ent) | 1))
#define BL_IS_REMOVED(bl_ent)  ((uint64_t)(bl_ent) & 1) 

#ifdef _WIN32
static HANDLE heap;
#endif

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
        ASSERT ((buf)->spinlock->link_count, "spinlock->link_count is 0 for %s", buf->name); \
        if (__atomic_sub_fetch (&spinlock->link_count, 1, __ATOMIC_RELAXED) == 0)       \
            FREE ((buf)->spinlock);                                                     \
        else                                                                            \
            __atomic_clear (&spinlock->lock, __ATOMIC_RELEASE);                         \
        spinlock = NULL;                                                                \
    }

#define reset_memory_pointer(buf)                                                       \
    /* careful to make sure that instructions are not re-ordered in either direction */ \
    char *old_memory __attribute__((unused)) = (buf)->memory;                           \
    __atomic_store_n (&(buf)->memory, (void*)0, __ATOMIC_RELEASE);                      \
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
    static rom names[] = BUFTYPE_NAMES;
    if (buf->type >= 0 && buf->type < BUF_NUM_TYPES) 
        return names[buf->type];
    else {
        char *s = malloc (32); // used for error printing
        sprintf (s, "invalid_buf_type=%u", buf->type);
        return s;
    }
}

const BufDescType buf_desc (ConstBufferP buf)
{
    #define IS_CONTEXT(ctxs) ((rom)buf >= (rom)ctxs && (rom)buf < (rom)&ctxs[MAX_DICTS])
    #define TAG_NAME(ctxs) ctxs[((rom)buf - (rom)ctxs) / (sizeof (ctxs) / MAX_DICTS /*note: may be different than sizeof(Context) due to word alignment*/)].tag_name

    if (!buf) return (BufDescType){ .s = "NULL" };

    // case: buffer is one of the Buffers within a Context - show tag_name
    rom tag_name = NULL;
    if (buf->vb && buf->vb != evb && IS_CONTEXT(buf->vb->contexts))
        tag_name = TAG_NAME (buf->vb->contexts);
    else if (buf->vb && buf->vb == evb && IS_CONTEXT(z_file->contexts))
        tag_name = TAG_NAME (z_file->contexts);

    BufDescType desc; // use static memory instead of malloc since we could be in the midst of a memory issue when this is called
    sprintf (desc.s, "\"%s\"%.20s memory=%p data=%p param=%"PRId64"(0x%016"PRIx64") len=%"PRIu64" size=%"PRId64" type=%s shared=%s%.16s promiscuous=%s spinlock=%p%.20s%.20s allocated in %s:%u%.20s", 
             buf->name ? buf->name : "(no name)", cond_str (tag_name, " ctx=", tag_name), 
             buf->memory, buf->data, buf->param, buf->param, buf->len, (uint64_t)buf->size, buf_type_name (buf), 
             TF(buf->shared), cond_int (buf->memory && buf->vb, " users=", BOLCOUNTER(buf)), TF(buf->promiscuous),
             buf->spinlock, // note: access spinlock fields only if buf->vb is set. If vb=NULL, buf was destroyed, and spinlock should be NULL - if its not, its likely a memory corruption
             cond_int (buf->spinlock && buf->vb, " locked=", buf->spinlock->lock), cond_int (buf->spinlock && buf->vb, " link_count=", buf->spinlock->link_count),
             buf->func ? buf->func : "(no func)", buf->code_line, cond_int (buf->vb, " by vb=", buf->vb->vblock_i));
    return desc;
}

// quick inline for internal buf_struct.c use check overflow and underflow in an allocated buffer
static void no_integrity (ConstBufferP buf, FUNCLINE, rom buf_func)
{
    flag.quiet = false;

    ASSERTW (BUNDERFLOW(buf) == UNDERFLOW_TRAP, "called from %s:%u to %s: Error in %s: buffer has corrupt underflow trap",
             func, code_line, buf_func, buf_desc(buf).s);

    ASSERTW (BOVERFLOW(buf) == OVERFLOW_TRAP, "called from %s:%u to %s: Error in %s: buffer has corrupt overflow trap",
            func, code_line, buf_func, buf_desc(buf).s);
    
    bool corruption_detected = buflist_test_overflows (buf->vb, buf_func);
    if (corruption_detected) buflist_test_overflows_all_other_vb (buf->vb, buf_func, true); // corruption not from this VB - test the others
    exit_on_error (false);
}

// quick inline wrapper
static inline void buf_verify_integrity (ConstBufferP buf, FUNCLINE, rom buf_func)
{
    if (buf->memory && buf->type == BUF_REGULAR && 
        (BUNDERFLOW(buf) != UNDERFLOW_TRAP || BOVERFLOW(buf) != OVERFLOW_TRAP))
        
        no_integrity (buf, func, code_line, buf_func);
}

// called from other modules for debugging memory issues. 
void buf_verify_do (ConstBufferP buf, rom msg, FUNCLINE)
{
    if (!buf || !buf->memory) return;

    BufferSpinlock *spinlock = buf->promiscuous ? buf_lock_promiscuous (buf, func, code_line) : NULL; // prevent frees or reallocs while we're testing
    if (buf->promiscuous && !spinlock) return; // by the time we acquired the lock, buf was already freed

    buf_verify_integrity (buf, __FUNCLINE, msg);

    buf_unlock;
}

static void buf_reset (BufferP buf)
{
    // first set "memory" to 0 - so buf_test_overflow doesn't test this buffer 
    reset_memory_pointer (buf);

    VBlockP save_vb = buf->vb;        // preserve vb because still in vb->buffer_list
    rom save_func   = buf->func;      // preserve func and code_line for buf_list* error reporting
    int save_line   = buf->code_line;

    memset (buf, 0, sizeof (Buffer)); // make this buffer UNALLOCATED

    buf->vb        = save_vb;
    buf->func      = save_func;
    buf->code_line = save_line;
}

static void buf_init (BufferP buf, char *memory, uint64_t size, FUNCLINE, rom name)
{
    // set some parameters before allocation so they can go into the error message in case of failure
    buf->func      = func;
    buf->code_line = code_line;

    if (name) 
        buf->name = name;
    else
        ASSERT (buf->name, "buffer has no name. func=%s:%u", buf->func, buf->code_line);

    if (!memory) { // malloc or realloc failed
        buflist_show_memory (true, 0, 0);

        ABORT ("%s: Out of memory%s. Details: %s:%u failed to allocate %s bytes. Buffer: %s", 
               global_cmd, 
               cond_int (IS_ZIP, ". Try running with a lower vblock size using --vblock. Current vblock size is: ", segconf.vb_size >> 20),
               func, code_line, str_int_commas (size + CTL_SIZE).s, buf_desc(buf).s);
    }

    buf->data = memory + sizeof (uint64_t);
    buf->size = size;

    BUNDERFLOW_(memory)      = UNDERFLOW_TRAP; // underflow protection
    BOVERFLOW_(memory, size) = OVERFLOW_TRAP;  // overflow prortection (underflow protection was copied with realloc)
    BOLCOUNTER_(memory, buf->data, size) = 1;  // 1 when memory is first allocated

    // only when we're done initializing - we update memory - that buf_test_overflow running concurrently doesn't test
    // half-initialized buffers
    __atomic_store_n (&buf->memory, memory, __ATOMIC_RELEASE); 
}

// allocates or enlarges buffer
// if it needs to enlarge a buffer fully overlaid by an overlay buffer - it abandons its memory (leaving it to
// the overlaid buffer) and allocates new memory
void buf_alloc_do (VBlockP vb, BufferP buf, uint64_t requested_size,
                   float grow_at_least_factor, // IF we need to allocate or reallocate physical memory, we get this much more than requested
                   rom name, FUNCLINE)      
{
    START_TIMER; // don't account time for these ^ calls - we're interested in actual allocations 

    // **** sanity checks ****
    ASSERT ((int64_t)requested_size > 0, "called from %s:%u: negative requested_size=%"PRId64" for name=%s", func, code_line, requested_size, name);

#define REQUEST_TOO_BIG_THREADSHOLD (3 GB)
    if (requested_size > REQUEST_TOO_BIG_THREADSHOLD && !buf->can_be_big) // use WARN instead of ASSERTW to have a place for breakpoint
        WARN ("Warning: buf_alloc called from %s:%u %s for \"%s\" requested %s. This is suspiciously high and might indicate a bug - please report to " EMAIL_SUPPORT ". vb->vblock_i=%u buf=%s line_i=%d",
              func, code_line, GENOZIP_CODE_VERSION, name, str_size (requested_size).s, vb->vblock_i, buf_desc (buf).s, vb->line_i);

    ASSERT (buf->type == BUF_REGULAR || buf->type == BUF_UNALLOCATED, "called from %s:%u: cannot buf_alloc a buffer of type %s. details: %s", 
            func, code_line, buf_type_name (buf), buf_desc (buf).s);

    if (!vb) vb = buf->vb;
    ASSERT (vb, "called from %s:%u: null vb", func, code_line);

    // if this happens: either 1. the wrong VB was given now, or when initially allocating this buffer OR
    // 2. VB was REALLOCed in vb_get_vb, but for some reason this buf->vb was not updated because it was not on the buffer list
    ASSERT (!buf->vb || vb == buf->vb, "called from %s:%u: buffer=%p has wrong VB: vb=%p (id=%u vblock_i=%u) but buf->vb=%p", 
            func, code_line, buf, vb, vb->id, vb->vblock_i, buf->vb);

    ASSERT (vb != evb || buf->promiscuous || threads_am_i_main_thread(), "called from %s:%u: A non-main thread is attempting to allocate an evb buffer \"%s\" with promiscuous=false", 
            func, code_line, name ? name : buf->name);

    // **** initial memory allocation ****

    // CASE 2: initial allocation: exactly at requested size
    if (!buf->memory) {
        ASSERT (requested_size <= MAX_BUFFER_SIZE, "called from %s:%u: Requested %s bytes which is beyond the Buffer maximum of %s",
                func, code_line, str_int_commas (requested_size).s, str_int_commas (MAX_BUFFER_SIZE).s);

        char *memory = (char *)buf_low_level_malloc (requested_size + CTL_SIZE, false, func, code_line);

        buf->type = BUF_REGULAR;

        buf_init (buf, memory, requested_size, func, code_line, name);
        
        if (buf != &vb->buffer_list) { // buffer_list buffer is added in vb_get_vb / vb_initialize_nonpool_vb
            if (!buf->promiscuous) // if promiscuous or buffer_list, already added
                buflist_add_buf (vb, buf, func, code_line);
            else
                ASSERT (buf->vb, "called from %s:%u: Expecting promiscuous buffer to be on the buffer_list: %s", func, code_line, buf_desc (buf).s);
        }
        
        goto done;
    }

    // **** realloc: calculate size include "growth" ****

    // add an epsilon to avoid floating point multiplication ending up slightly less than the integer
    grow_at_least_factor = MAX_(1.0001, grow_at_least_factor); 

    // grow us requested - rounding up to 64 bit boundary to avoid aliasing errors with the overflow indicator
    uint64_t new_size = MIN_(MAX_BUFFER_SIZE, ROUNDUP8((uint64_t)(requested_size * grow_at_least_factor)));

    ASSERT (new_size >= requested_size, "called from %s:%u: allocated too little memory for buffer %s: requested=%"PRIu64", allocated=%"PRIu64". vb->vblock_i=%u", 
            func, code_line, buf_desc (buf).s, requested_size, new_size, vb->vblock_i); // floating point paranoia

    // CASE 3: realloc of a non-shared - use realloc that will extend instead of malloc & copy if possible 
    if (!buf->shared) {    
        buf_lock_if (buf, buf->spinlock); // promiscous or buf_list or caller lock...

        buf_verify_integrity (buf, func, code_line, "buf_alloc_do");

        reset_memory_pointer (buf);
         
        char *new_memory = (char *)buf_low_level_realloc (old_memory, new_size + CTL_SIZE, name, func, code_line);
        buf_init (buf, new_memory, new_size, func, code_line, name);
    
        buf_unlock;
    }

    else { // shared cases
        buf_lock (buf);

        buf_verify_integrity (buf, func, code_line, "buf_alloc_do(shared)");

        // CASE 4: shared: currently no overlayers and standard "data" value - we can realloc
        if (buf_user_count (buf) == 1 && 
            (!buf->data || buf->data == buf->memory + sizeof(uint64_t))) {
            reset_memory_pointer (buf);

            char *new_memory = (char *)buf_low_level_realloc (old_memory, new_size + CTL_SIZE, name, func, code_line);
            buf_init (buf, new_memory, new_size, func, code_line, name);
        }

        // CASE 5: shared: have overlayers, or non-standard data value - malloc & copy
        else {
            char *new_memory = (char *)buf_low_level_malloc (new_size + CTL_SIZE, false, func, code_line);

            memcpy (new_memory + sizeof (uint64_t), buf->data, buf->size); // copy old data
            uint16_t user_count = buf_decrement_user_count (buf); // we are no longer using the old memory

            reset_memory_pointer (buf);
            
            if (!user_count)
                buf_low_level_free (old_memory, func, code_line);

            buf_init (buf, new_memory, new_size, func, code_line, name);
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
void buf_trim_do (BufferP buf, uint64_t size, FUNCLINE)
{
    if (size >= buf->size || !buf->memory) return; // nothing to do - size if already smaller

    ASSERT (!buf->shared, "%s:%u trimming is not currently supported on shared buffers: buf=%s", func, code_line, buf_desc(buf).s);

    buf_lock_if (buf, buf->spinlock); // promiscous or buf_list or caller lock...

    buf_verify_integrity (buf, func, code_line, "buf_alloc_do");

    reset_memory_pointer (buf);
        
    char *new_memory = (char *)buf_low_level_realloc (old_memory, size + CTL_SIZE, buf->name, func, code_line);
    buf_init (buf, new_memory, size, func, code_line, buf->name);

    buf_unlock;
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
    if (!buf->spinlock) return;

    ASSERT (!buf->memory, "cannot remove spinlock: buf has memory: %s", buf_desc(buf).s);

    if (!__atomic_sub_fetch (&buf->spinlock->link_count, 1, __ATOMIC_RELAXED))
        FREE (buf->spinlock); // we were the only uses of the spinlock

    buf->shared = buf->promiscuous = false;
}

// an overlay buffer is a buffer using some of the memory of another buffer - it doesn't have its own memory
void buf_overlay_do (VBlockP vb, 
                     BufferP top_buf, // dst 
                     BufferP bottom_buf, 
                     uint64_t start_in_bottom, // 0 means full overlay, and copy len 
                     FUNCLINE, rom name)
{   
    START_TIMER;

    // if this buffer was used by a previous VB as a bottom buffer - we need to "destroy" it first
    if (top_buf->type == BUF_REGULAR && top_buf->data == NULL && 
        (top_buf->memory || (top_buf->spinlock && top_buf->spinlock != bottom_buf->spinlock))) 
        buf_destroy (*top_buf);

    ASSERT (top_buf->type == BUF_UNALLOCATED, "%s: Call from %s:%u: cannot buf_overlay to a buffer %s already in use", VB_NAME, func, code_line, buf_desc (top_buf).s);

    // overlaying a SHM buffer, just creates another SHM buffer
    if (bottom_buf->type == BUF_SHM) {
        ASSERT (start_in_bottom < bottom_buf->size, "%s: Call from %s:%u: expecting start_in_bottom=%"PRIu64" < bottom_buf->size=%"PRIu64, 
                VB_NAME, func, code_line, start_in_bottom, (uint64_t)bottom_buf->size);

        buf_attach_to_shm_do (vb, top_buf, 
                              bottom_buf->memory, bottom_buf->size, start_in_bottom, 
                              func, code_line, name);

        if (!start_in_bottom) top_buf->len = bottom_buf->len;
        return;
    }

    ASSERT (bottom_buf->shared && bottom_buf->spinlock, 
            "%s: Call from %s:%u: expecting bottom_buf %s to have a spinlock and shared=true", VB_NAME, func, code_line, buf_desc (bottom_buf).s);

    ASSERT (bottom_buf->type == BUF_REGULAR,
            "%s: Call from %s:%u: bottom_buf %s in buf_overlay must be a bottom or shm buffer", VB_NAME, func, code_line, buf_desc (bottom_buf).s);

    top_buf->type      = BUF_REGULAR;
    top_buf->name      = name;
    top_buf->len       = start_in_bottom ? 0ULL : bottom_buf->len;
    top_buf->func      = func;
    top_buf->code_line = code_line;
    top_buf->shared    = true;
    top_buf->spinlock  = bottom_buf->spinlock;

    // note: while we're waiting for the lock, bottom_buf may be reallocing - but this doesn't change spin_lock
    buf_lock (bottom_buf); // locking spinlock shared between top, bottom buffers and all other overlayers

    buf_verify_integrity (bottom_buf, func, code_line, "buf_overlay_do");

    ASSERT (start_in_bottom < bottom_buf->size, 
            "called from %s:%u: not enough room in bottom buffer for overlaid buf: start_in_bottom=%"PRIu64" but bottom_buf.size=%"PRIu64,
            func, code_line, start_in_bottom, (uint64_t)bottom_buf->size);

    // note: data+size MUST be at the control region, as we have the overlay counter there
    top_buf->size = bottom_buf->size - start_in_bottom;
    top_buf->data = bottom_buf->data + start_in_bottom;

    // increment spinlock users and memory users (note: spinlock link cou t >= memory_users bc a spinlock can be used by multiple memories in case of a realloc-induced split)
    __atomic_fetch_add (&bottom_buf->spinlock->link_count, 1, __ATOMIC_RELAXED);    

    buf_increment_user_count (bottom_buf);
    
    // final step
    __atomic_store_n (&top_buf->memory, bottom_buf->memory, __ATOMIC_RELEASE); 

    buf_unlock;

    buflist_add_buf (vb, top_buf, func, code_line); 

    COPY_TIMER (buf_overlay_do);
}

void buf_attach_to_shm_do (VBlockP vb, BufferP buf, void *memory, uint64_t size, uint64_t start, FUNCLINE, rom name)
{
    // preserve len and param (i.e. nbits and nwords)
    uint64_t save_len   = buf->len;
    uint64_t save_param = buf->param;

    // if this buffer was used by a previous VB as a regular buffer - we need to "destroy" it first
    if (buf->vb) 
        buf_destroy (*buf);

    *buf = (Buffer){
        .type      = BUF_SHM,
        .name      = name,
        .func      = func,
        .code_line = code_line,
        .vb        = vb,
        .memory    = memory,
        .data      = memory + start, // note: no control area in shm buffers
        .size      = size,
        .len       = save_len,
        .param     = save_param
    };
}

void buf_free_do (BufferP buf, FUNCLINE) 
{
    switch (buf->type) {

        case BUF_REGULAR: {
            START_TIMER;
        
            buf_lock_if (buf, buf->spinlock);

            buf_verify_integrity (buf, func, code_line, "buf_free_do");

            uint16_t user_count = buf_user_count (buf);
            ASSERT (user_count, "%s:%u: user_count=0 in buffer %s", func, code_line, buf_desc(buf).s);

            // case: if user_count >= 2: reset buffer and leave memory to other user 
            // (often: compute thread resets vb buffer and leaves memory to main thread z_file buffer)
            if (user_count >= 2) {
                user_count = buf_decrement_user_count (buf);

                reset_memory_pointer (buf);
                buf->size = 0;
                buf->name = buf->func = 0;
                buf->code_line = 0;
                buf->type = BUF_UNALLOCATED;
                // preserved: vb, promiscuous, shared, spinlock + still on its VB's buffer_list
            }

            // if last remaining user is a partial overlay - increase size to be based on entire memory
            else if (buf->data && (buf->data - buf->memory != sizeof (uint64_t)))
                buf->size += (buf->data - buf->memory - sizeof (uint64_t));

            buf->data        = NULL; 
            buf->can_be_big  = false;
            buf->len         = 0;
            buf->param       = 0;

            buf_unlock;

            if (buf->vb == evb) COPY_TIMER_EVB (buf_free_main); // works even for promiscuous bc uses atomic 
            else { VBlockP vb = buf->vb; COPY_TIMER (buf_free_compute); };

            break;
        }
        case BUF_UNALLOCATED: // reset len and param that may be used even without allocating the buffer
            buf->len         = 0;
            buf->param       = 0;
            break;

        case BUF_SHM:
            buf_reset (buf);
            break;

        default:
            ABORT ("Error: invalid buf->type=%s", buf_type_name (buf));
    }
} 

void buf_destroy_do_do (BufListEnt *ent, FUNCLINE)
{
    if (!ent) return;
    START_TIMER;

    BufferP buf = ent->buf;
    VBlockP vb = buf->vb;

    if (flag.debug_memory==1) 
        iprintf ("Destroy %s: buf_addr=%p vb->id=%d buf_i=%u\n", buf_desc (buf).s, buf, buf->vb->id, BNUM (buf->vb->buffer_list, ent));

    // remove from buffer list
    BL_SET_REMOVED (ent->buf);
    buf->vb = NULL;

    switch (buf->type) {
        case BUF_REGULAR : { 
            buf_lock_if (buf, buf->spinlock);

            if (buf->memory) { // NULL if buffer was disowned
                buf_verify_integrity (buf, func, code_line, "buf_destroy_do");
                uint16_t remaining_user_count = buf_decrement_user_count (buf);

                // first set "memory" to 0 - so buf_test_overflow doesn't test this buffer after destroyed
                reset_memory_pointer (buf);

                if (!remaining_user_count) 
                    buf_low_level_free (old_memory, func, code_line); 
            }

            buf_unlock_and_decrement_lock_links; // also frees spinlock if we're the last user
            break;
        }

        case BUF_SHM : 
            buf_free (*buf); 
            break;
        
        case BUF_UNALLOCATED : {
            buf_lock_if (buf, buf->spinlock); // possibly set as promiscuous but never allocated
            buf_unlock_and_decrement_lock_links;
            break;
        }

        default : ABORT ("called from %s:%u: Error in buf_destroy_do: invalid buffer type %s", func, code_line, buf_type_name (buf));
    }

    buf_reset (buf);

    if (vb==evb) COPY_TIMER_EVB (buf_destroy_do_do_main);
    else         COPY_TIMER (buf_destroy_do_do_compute);
}

void buf_destroy_do (BufferP buf, FUNCLINE)
{
    if (!buf || 
        (!buf->vb && buf->type == BUF_UNALLOCATED) || // never allocated 
        flag.let_OS_cleanup_on_exit) return; // nothing to do (we don't destroy on exit, as the exiting thread may not be able to remove from buf_list)

    BufListEnt *ent;

    if (buf->type == BUF_DISOWNED) {
        buf_low_level_free (buf->memory, func, code_line);
        *buf = (Buffer){};  
    }

    else if ((ent = buflist_find_buf (buf->vb, buf, SOFT_FAIL))) 
        buf_destroy_do_do (ent, func, code_line);
    
    else
        *buf = (Buffer){};  
}

// similar to buf_move, but also moves buffer between buf_lists. can be run by the main thread only.
// IMPORTANT: only works when called from main thread, when BOTH src and dst VB are in full control of main thread, so that there
// no chance another thread is concurrently modifying the buf_list of the src or dst VBs 
void buf_grab_do (VBlockP dst_vb, BufferP dst_buf, rom dst_name/*optional*/, BufferP src_buf, FUNCLINE)
{
    ASSERTMAINTHREAD;
    ASSERT (src_buf, "called from %s:%u: buf is NULL", func, code_line);
    if (src_buf->type == BUF_UNALLOCATED) return; // nothing to grab

    ASSERT (src_buf->type == BUF_REGULAR && !src_buf->shared, "called from %s:%u: this function can only be called for a non-shared REGULAR buf", func, code_line);
    ASSERT (dst_buf->type == BUF_UNALLOCATED, "called from %s:%u: expecting dst_buf to be UNALLOCATED", func, code_line);

    reset_memory_pointer (src_buf);

    dst_buf->type     = BUF_REGULAR;
    dst_buf->len      = src_buf->len;
    dst_buf->param    = src_buf->param;
    dst_buf->spinlock = src_buf->spinlock;
    buf_init (dst_buf, old_memory, src_buf->size, func, code_line, dst_name ? dst_name : src_buf->name);

    buflist_add_buf (dst_vb, dst_buf, func, code_line);

    // remove src_buf from buffer list
    buf_lock (&src_buf->vb->buffer_list);
    BL_SET_REMOVED (buflist_find_buf (src_buf->vb, src_buf, HARD_FAIL)->buf);
    buf_unlock;
    
    src_buf->vb = NULL;

    buf_reset (src_buf);
}

// moves all the data between buffers in the same VB, keeping the same buffer_list entry. 
void buf_move_do (VBlockP vb, BufferP dst_buf, rom dst_name/*optional*/, BufferP src_buf, FUNCLINE)
{
    ASSERT (!dst_buf->data, "%s:%u: dst_buf is not empty: %s", func, code_line, buf_desc(dst_buf).s);
    ASSERT (src_buf->type == BUF_REGULAR, "%s:%u: src_buf is %s", func, code_line, buf_type_name (src_buf));
    ASSERT (src_buf->vb == vb, "%s:%u: src_buf has wrong vb", func, code_line);
    ASSERT (buf_user_count (src_buf) == 1, "%s:%u: expecting src_buf to have 1 user but it has %u", func, code_line, buf_user_count (src_buf));
    buf_verify_integrity (src_buf, func, code_line, "buf_move");

    if (dst_buf->vb) buf_destroy (*dst_buf); // also remove from buffer_list

    if (!dst_name) dst_name = src_buf->name;

    *dst_buf = (Buffer){ .func       = func,
                         .code_line  = code_line,
                         .name       = dst_name,
                         .data       = src_buf->data,
                         .size       = src_buf->size,
                         .param      = src_buf->param,
                         .len        = src_buf->len,
                         .type       = BUF_REGULAR,
                         .can_be_big = src_buf->can_be_big,
                         .vb         = vb,
                         .memory     = src_buf->memory }; // no need for atomic_store, bc buflist_move_buf unlocks with ATOMIC_RELEASE

    // make the buffer_list entry of src_buf now point to dst_buf     
    buflist_move_buf (vb, dst_buf, dst_name, src_buf, func, code_line);

    // promiscuous, shared and spinlock are NOT moved. 
    buf_remove_spinlock (src_buf);

    reset_memory_pointer (src_buf);
    *src_buf = (Buffer){};
}

// move buffer struct to a new location, without adding it the buffer list 
void buf_disown_do (VBlockP vb, BufferP src_buf, BufferP dst_buf, bool make_a_copy, FUNCLINE)
{
    ASSERT (vb == src_buf->vb, "called from %s:%u: buf->vb mismatches vb. vb->vblock_i=%u", func, code_line, vb->vblock_i);
    ASSERT (src_buf->type == BUF_REGULAR, "called from %s:%u: expecting src_buf to be BUF_REGULAR, buf_type=%s", func, code_line, buf_type_name(src_buf));
    ASSERT (!src_buf->promiscuous, "called from %s:%u: src_buf cannot be promiscuous", func, code_line);
    ASSERT (!src_buf->shared, "called from %s:%u: src_buf cannot be shared", func, code_line);
    ASSERT (dst_buf->type == BUF_UNALLOCATED, "called from %s:%u: expecting dst_buf to be BUF_UNALLOCATED, buf_type=%s", func, code_line, buf_type_name(dst_buf));

    buf_verify_integrity (src_buf, func, code_line, "buf_disown_do");

    *dst_buf = (Buffer) {
        .memory     = src_buf->memory,
        .data       = src_buf->data,
        .len        = src_buf->len,
        .param      = src_buf->param,
        .can_be_big = src_buf->can_be_big,
        .name       = src_buf->name,
        .size       = src_buf->size,
        .type       = BUF_DISOWNED,
        .func       = func,
        .code_line  = code_line };
    
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

//-------------------------
// low-level functions
//-------------------------

void buf_low_level_free (void *p, FUNCLINE)
{
    if (!p) return; // nothing to do

    START_TIMER;
    bool p_is_evb = (p == evb);

    if (flag.debug_memory==1) 
        iprintf ("Memory freed by free(): %p %s:%u\n", p, func, code_line);

#ifndef _WIN32
    free (p);
#else
    ASSERT (HeapFree (heap, 0, p), "HeapFree failed: %s", str_win_error());
#endif

    if (!p_is_evb) // unless we just freed evb...
        COPY_TIMER_EVB (buf_low_level_free);
}

void *buf_low_level_realloc (void *p, size_t size, rom name, FUNCLINE)
{
#ifndef _WIN32
    void *new = realloc (p, size);
#else
    void *new = HeapReAlloc (heap, 0, p, size);
#endif

    ASSERTW (new, "Out of memory in %s:%u: realloc failed (name=%s size=%"PRIu64" bytes). %s", func, code_line, name, (uint64_t)size, 
             IS_ZIP ? "Try limiting the number of concurrent threads with --threads (affects speed) or reducing the amount of data processed by each thread with --vblock (affects compression ratio)" : "");

    if (flag.debug_memory && size >= flag.debug_memory) {
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wpragmas"         // avoid warning if "-Wuse-after-free" is not defined in this version of gcc
#pragma GCC diagnostic ignored "-Wunknown-warning-option" // same
#pragma GCC diagnostic ignored "-Wuse-after-free"  // avoid compiler warning of using p after it is freed
        iprintf ("realloc(): old=%p new=%p name=%s size=%"PRIu64" %s:%u\n", p, new, name, (uint64_t)size, func, code_line);
#pragma GCC diagnostic pop
    }

    return new;
}

void *buf_low_level_malloc (size_t size, bool zero, FUNCLINE)
{
#ifndef _WIN32
    void *new = malloc (size);
#else
    void *new = HeapAlloc (heap, zero ? HEAP_ZERO_MEMORY : 0, size);
#endif
    ASSERT (new, "Out of memory in %s:%u: malloc failed (size=%"PRIu64" bytes). %s", func, code_line, (uint64_t)size,
            IS_ZIP ? "Try limiting the number of concurrent threads with --threads (affects speed) or reducing the amount of data processed by each thread with --vblock (affects compression ratio)" : "");

    if (flag.debug_memory && size >= flag.debug_memory) 
        iprintf ("malloc(): %p size=%"PRIu64" %s:%u\n", new, (uint64_t)size, func, code_line);

#ifndef _WIN32
    if (zero) memset (new, 0, size);
#endif
    
    return new;
}

void buf_low_level_release_memory_back_to_kernel (void)
{
#ifdef __linux__
    malloc_trim (0);                        // return whole free pages to the kernel
#elif defined __APPLE__
    malloc_zone_pressure_relief (NULL, 0);  // tell OS that this process is interested in participating in "pressure relief" - freeing memory. OS will decide if and when to actually release the memory.
#elif defined _WIN32
    HeapCompact (heap, 0);                  // return blocks marked for "deferred free" to the kernel
#endif
}

uint64_t buf_mem_size (ConstBufferP buf) 
{ 
    // note: this calculation does not discount for memory overlaid in multiple buffers (with shared=true)
    return buf->type != BUF_REGULAR ? 0
         : !buf->memory             ? 0
         : buf->data                ? ((buf->data - buf->memory) + buf->size + sizeof (uint64_t) + sizeof(uint16_t)) // might be "partial overlay"
         :                            (buf->size + CTL_SIZE); 
}
