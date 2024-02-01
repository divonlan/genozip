// ------------------------------------------------------------------
//   buf_list.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifdef _WIN32
#include <windows.h>
#endif
#include "genozip.h"
#include "vblock.h"
#include "buf_struct.h"
#include "buf_list.h"
#include "context.h"
#include "file.h"
#include "threads.h"
#include "reference.h"
#include "gencomp.h"
#include "arch.h"

// We mark a buffer list entry as removed, by setting its LSb. This keeps the buffer list sorted, in case it is sorted.
// Normally the LSb is 0 as buffers are word-aligned (verified in buflist_add_buf)
#define BL_SET_REMOVED(bl_ent) bl_ent = ((BufferP)((uint64_t)bl_ent | 1))
#define BL_IS_REMOVED(bl_ent)  ((uint64_t)bl_ent & 1) 

#define is_sorted prm8[0] 

//-----------------------------
// Iterators
//-----------------------------

// sorts buffer list ahead of multiple calls to buflist_find_buf to speed it up
static ASCENDING_SORTER (buflist_sorter, BufListEnt, buf)
void buflist_sort (VBlockP vb, bool already_locked)
{
    START_TIMER;

    if (flag.let_OS_cleanup_on_exit ||      // buf_destroy is disabled, so no need to sort either
        vb->buffer_list.is_sorted) return;   // already sorted

    buf_lock_if (&vb->buffer_list, !already_locked);

    qsort (B1ST(BufListEnt, vb->buffer_list), vb->buffer_list.len32, sizeof(BufListEnt), buflist_sorter);
    vb->buffer_list.is_sorted = true; // now it's sorted

    buf_unlock;

    COPY_TIMER (buflist_sort);
}

// remove from buffer_list of this vb
// thread safety: only the thread owning the VB of the buffer (main thread of evb) can remove a buffer
// from the buf list, OR the main thread may remove for all VBs, IF no compute thread is running
static BINARY_SEARCHER (buflist_finder, BufListEnt, ConstBufferP, buf, true)
BufListEnt *buflist_find_buf (VBlockP vb, ConstBufferP buf, FailType soft_fail)
{
    START_TIMER;

    ASSERTNOTNULL (buf);
    BufListEnt *ent = NULL;

    if (buf->vb) { // this buffer is on the buf_list 
        if (!vb->buffer_list.is_sorted) { // linear search if not sorted
            for_buf_back (BufListEnt, my_ent, vb->buffer_list) // back - we often destroy recently added
                if (my_ent->buf == buf) {
                    ent = my_ent;
                    break;
                }
        }
        else // binary search if sorted
            ent = binary_search (buflist_finder, BufListEnt, vb->buffer_list, buf);
    }
    
    // note: it is possible that the buffer is not found in the list if it is never allocated or destroyed more than once. that's fine.
    ASSERT (ent || soft_fail == SOFT_FAIL, "Cannot find buf=%p on buffer_list: %s", buf, buf_desc(buf).s);

    COPY_TIMER (buflist_find_buf);
    return ent;
}

static void buflist_foreach_buffer (VBlockP vb, bool (*callback)(ConstBufferP, VBlockPoolType pool_type, int vb_id, void *arg), VBlockPoolType pool_type, void *arg)
{
    buf_lock (&vb->buffer_list);

    for_buf2 (BufListEnt, bl, buf_i, vb->buffer_list) {
        #ifdef __linux__
            // if (!BL_IS_REMOVED (*bl) && access ((rom)*bl, F_OK) == -1 && errno == EFAULT)
            //     ABORT ("buffer structure inaccessible (invalid pointer=%p) in buf_i=%u of vb_id=%d", *bl, buf_i, vb_id); 
        #elif defined WIN32
            if (!BL_IS_REMOVED (bl->buf) && IsBadReadPtr (bl->buf, sizeof (Buffer))) 
                ABORT ("buffer structure inaccessible (invalid pointer=%p) in buf_i=%u of vb_id=%d", bl->buf, buf_i, vb->id); 
        #endif

        if (BL_IS_REMOVED (bl->buf) || !bl->buf->memory) continue; // exclude destroyed, not-yet-allocated, overlay buffers and buffers that were src in buf_move

        if (!callback (bl->buf, pool_type, vb->id, arg))
            break;
    }

    buf_unlock;
}

static void buflist_foreach_buffer_in_each_vbs (bool (*callback)(ConstBufferP, VBlockPoolType pool_type, int vb_id, void *arg), void *arg)
{
    for (VBlockPoolType pool_type=0; pool_type < NUM_POOL_TYPES; pool_type++) {

        VBlockPool *vb_pool = vb_get_pool (pool_type, SOFT_FAIL);
        if (!vb_pool) continue;

        // note: we don't cover EVB (only in DEBUG and --show-mem) as it segfaults for an unknown reason (likely the Buffer structure itself resides with a structure
        // that is freed (eg File)). TODO: debug this.
        for (int vb_id=0; vb_id < (int)vb_pool->num_allocated_vbs; vb_id++) {

            VBlockP vb = vb_get_from_pool (vb_pool, vb_id);
            if (!vb) continue;

            buflist_foreach_buffer (vb, callback, pool_type, arg);
        }
    }

    // non-pool VBs 
    if (threads_am_i_main_thread() || flag.debug || flag.show_memory) 
        for (int vb_id=-1; vb_id >= -NUM_NONPOOL_VBs; vb_id--) {
            VBlockP vb = vb_get_nonpool_vb (vb_id);
            if (vb) buflist_foreach_buffer (vb, callback, -1, arg);
        }
}

// reports locked buflists - used in case of a deadlock
void buflist_who_is_locked (void)
{
    for (VBlockPoolType pool_type=0; pool_type < NUM_POOL_TYPES; pool_type++) {

        VBlockPool *vb_pool = vb_get_pool (pool_type, SOFT_FAIL);
        if (!vb_pool) continue;

        // note: we don't cover EVB (only in DEBUG and --show-mem) as it segfaults for an unknown reason (likely the Buffer structure itself resides with a structure
        // that is freed (eg File)). TODO: debug this.
        for (int vb_id=0; vb_id < (int)vb_pool->num_allocated_vbs; vb_id++) {

            VBlockP vb = vb_get_from_pool (vb_pool, vb_id);
            if (!vb) continue;

            if (vb->buffer_list.spinlock && vb->buffer_list.spinlock->lock)
                fprintf (stderr, "vb=%d vb_id=%d buffer_list is locked\n", vb->vblock_i, vb_id);
        }
    }

    // non-pool VBs 
    if (threads_am_i_main_thread() || flag.debug || flag.show_memory) 
        for (int vb_id=-1; vb_id >= -NUM_NONPOOL_VBs; vb_id--) {
            VBlockP vb = vb_get_nonpool_vb (vb_id);
            if (vb && vb->buffer_list.spinlock && vb->buffer_list.spinlock->lock)
                fprintf (stderr, "non-pool vb_id=%d buffer_list is locked\n", vb_id);
        }    
}

void buflist_display (VBlockP vb) // for debugging
{
    uint32_t vb_size = get_vb_size (vb->data_type);

    iprintf ("Buffers of vblock_i=%d id=%u range=[%p - %p]:\n", vb->vblock_i, vb->id, vb, (rom)vb + vb_size - 1);
    for_buf2 (BufListEnt, ent, ent_i, vb->buffer_list)
        iprintf ("%3u\tin_vb=%s\t%-20s\t%s:%u\t%p%s\n", 
                 ent_i, TF(is_p_in_range (ent->buf, vb, vb_size)), ent->name, ent->func, ent->code_line, ent->buf,
                 cond_str (BL_IS_REMOVED (ent->buf), "\t", "REMOVED"));
}

//---------------------------------------------
// add, change and remove items in buflist
//---------------------------------------------

// thread safety: only the thread owning the VB of the buffer (main thread of evb) can add a buffer
// to the buf list OR it may be added by the main thread IF the compute thread of this VB is not running yet
void buflist_add_buf (VBlockP vb, BufferP buf, FUNCLINE)
{
    if (buf->vb == vb) return; // already in buf_list - nothing to do
    START_TIMER;

    ASSERT (buf->name, "buf has no name: %s", buf_desc(buf).s);

    ASSERT (vb != evb || threads_am_i_main_thread(), "called from %s:%u: A thread other than main thread is attempting to modify evb buffer list: buf=%s", 
            func, code_line, buf_desc(buf).s);

    ASSERT (!buf->vb, "called from %s:%u: cannot add buffer %s to buf_list of vb->vblock_i=%u because it is already in buf_list of vb_i=%u.",
             func, code_line, buf_desc(buf).s, vb->vblock_i, buf->vb->vblock_i);    

    ASSERT (((uint64_t)buf & 7)==0, "called from %s:%u: Expecting buffer %s of vb->vblock_i=%u to be word-aligned, but its address is %p",
             func, code_line, buf_desc(buf).s, vb->vblock_i, buf); // BL_SET_REMOVED expects this

#define INITIAL_MAX_MEM_NUM_BUFFERS 10000 /* for files that have ht,gt,phase,variant,and line - the factor would be about 5.5 so there will be 1 realloc per vb, but most files don't */
    BufferP bl = &vb->buffer_list;

    buf_alloc (vb, bl, 1, INITIAL_MAX_MEM_NUM_BUFFERS, BufListEnt, 2, "buffer_list");

    buf_lock (bl);

    BNXT (BufListEnt, *bl) = (BufListEnt){ .buf = buf, .func = func, .code_line = code_line, .name = buf->name };
    if (bl->len32 >= 2 && buf < (BLST(BufListEnt, *bl)-1)->buf)
        bl->is_sorted = false; // buffer_list is not sorted anymore

    buf->vb = vb; // successfully added to buf list

    buf_unlock;
    COPY_TIMER (buflist_add_buf);

    #define DISPLAY_ALLOCS_AFTER 0 // display allocations, except the first X allocations. reallocs are always displayed
    if (flag.debug_memory==1 && vb->buffer_list.len > DISPLAY_ALLOCS_AFTER)
        iprintf ("buflist_add_buf: &buf=%p (%s): %s: size=%"PRIu64" buffer=%p func=%s:%u vb->id=%d buf_i=%u\n", 
                 buf, func, buf_desc(buf).s, (uint64_t)buf->size, buf, func, code_line, vb->id, vb->buffer_list.len32-1);    
}

void buflist_remove_buf (BufferP buf, FUNCLINE)
{
    START_TIMER;

    VBlockP vb = buf->vb;

    ASSERTISALLOCED (vb->buffer_list);

    if (!vb_is_valid (vb)) { // this is bug 576
        WARN ("Unexpected, but unharmful: cannot remove buf=%p from buffer list, because buf->vb=%p refers to a VB that no longer exists", buf, vb);
        WARN ("%s allocated in %s:%u", buf->name, buf->func, buf->code_line); // separate line in case it segfaults if buf is corrupt
        return;
    }

    ASSERT (vb != evb || threads_am_i_main_thread(), "A thread other than main thread is attempting to modify evb buffer list: buf=%s", buf_desc (buf).s);

    BufListEnt *ent = buflist_find_buf (vb, buf, HARD_FAIL); 

    buf_lock (&vb->buffer_list);

    BL_SET_REMOVED (ent->buf);
    buf->vb = NULL;

    buf_unlock;

    COPY_TIMER (buflist_remove_buf);
}

// change entry from pointing to one buffer to pointing to another
void buflist_move_buf (VBlockP vb, BufferP new_buf, rom new_name, ConstBufferP old_buf, FUNCLINE)
{
    BufListEnt *ent = buflist_find_buf (vb, old_buf, HARD_FAIL);
    
    buf_lock (&vb->buffer_list);
   
    *ent = (BufListEnt){ .buf       = new_buf,
                         .func      = func,
                         .code_line = code_line,
                         .name      = new_name  };

    if (vb->buffer_list.is_sorted && 
          ( (ent != B1ST(BufListEnt, vb->buffer_list) && new_buf < (ent-1)->buf) ||
            (ent != BLST(BufListEnt, vb->buffer_list) && new_buf > (ent+1)->buf)))
        vb->buffer_list.is_sorted = false;

    buf_unlock;
}

// called from main thread when preparing the VB, so it is safe to change
// the individual buffers' "vb" field. however, buflist is still locked in
// case live VBs verify concurrently.
void buflist_update_vb_addr_change (VBlockP new_vb, ConstVBlockP old_vb)
{
    DataType old_dt = new_vb->data_type_alloced; // still not updated
    uint32_t old_vb_size = get_vb_size (old_dt);

    buf_lock (&new_vb->buffer_list);

    for_buf (BufListEnt, ent, new_vb->buffer_list) {

        // buf_list must be compacted before calling this function.
        ASSERT (!BL_IS_REMOVED (ent->buf), "buf_list not compacted: vb_i=%u name=%s %s:%u", 
                new_vb->vblock_i, ent->name, ent->func, ent->code_line);

        // if this buffer is within common area of the moved VB - update its address to the same offset relative to the new vb address
        if (is_p_in_range (ent->buf, old_vb, sizeof (VBlock))) {
            ent->buf = (BufferP)((char*)new_vb + ((char*)ent->buf - (char*)old_vb)); 

            // update VB address within the buffer
            (ent->buf)->vb = new_vb; 
        }
        else 
            ASSERT (!is_p_in_range (ent->buf, old_vb, old_vb_size), "Expecting private buffer \"%s\" allocated=%s:%u of old allocated_data_type=%s to have been destroyed and removed from buffer_list", 
                    ent->name, ent->func, ent->code_line, dt_name (old_dt));
    }

    buf_unlock;
}

// after removing marking buffers as removed, actually remove them from the list
void buflist_compact (VBlockP vb)
{
    if (!vb) return;

    buf_lock_(&vb->buffer_list);
    buf_remove_items_except_(BufListEnt, vb->buffer_list, !BL_IS_REMOVED (ent->buf));
    buf_unlock;
}

bool buflist_locate (ConstBufferP buf, rom prefix /*NULL if silent*/)
{
    for (VBlockPoolType pool_type=0; pool_type < NUM_POOL_TYPES; pool_type++) {

        VBlockPoolP vb_pool = vb_get_pool (pool_type, SOFT_FAIL);
        if (!vb_pool) continue;

        // note: we don't cover EVB (only in DEBUG and --show-mem) as it segfaults for an unknown reason (likely the Buffer structure itself resides with a structure
        // that is freed (eg File)). TODO: debug this.
        for (int vb_id=0; vb_id < (int)vb_pool->num_allocated_vbs; vb_id++) {

            VBlockP vb = vb_get_from_pool (vb_pool, vb_id);
            if (!vb) continue;

            if (vb_buf_locate (vb, buf)) {
                if (prefix) 
                    fprintf (stderr, "%s %p located in pool=%s vb_id=%d: %s\n", 
                             prefix, buf, vb_pool->name, vb_id, buf_desc (buf).s);
                return true;
            }
        }
    }

    #define OBJECT_TEST(locator, obj)                                                                   \
        if (locator (obj, buf)) {                                                                \
            if (prefix)                                                                                \
                fprintf (stderr, "%s %p located in " #obj ": %s\n", prefix, buf, buf_desc (buf).s); \
            return true;                                                                                \
        }

    for (int i=1; i <= NUM_NONPOOL_VBs; i++) 
        OBJECT_TEST(vb_buf_locate, vb_get_nonpool_vb (-i));

    OBJECT_TEST(file_buf_locate, z_file);
    OBJECT_TEST(file_buf_locate, txt_file);
    OBJECT_TEST(ref_buf_locate, gref);
    OBJECT_TEST(ref_buf_locate, prim_ref);

    void *depn=0, *componentsP=0, *queueP=0, *preabsorb_queue=0; // dummies
    OBJECT_TEST(gencomp_buf_locate_depn, depn);
    OBJECT_TEST(gencomp_buf_locate_componentsP, componentsP);
    OBJECT_TEST(gencomp_buf_locate_queueP, queueP);
    OBJECT_TEST(gencomp_buf_locate_preabsorb_queue, preabsorb_queue);

    if (prefix) 
        fprintf (stderr, "%s %p is not located in any of the objects tested.", prefix, buf);

    return false;
}

// returns "" - just to make it easy to include in an ASSERT
static rom buflist_locate_s (ConstBufferP buf)
{
    buflist_locate (buf, "Duplicate Buffer");
    return "";
}

//-----------------------------
// Freers and Destroyers
//-----------------------------

// frees all buffers on vb->buffer_list, and also zeros the spaces between Buffers in the VB.
// does not affect "fields that survive buflist_free_vb" as defined in vblock.h
void buflist_free_vb (VBlockP vb) 
{
    if (flag.let_OS_cleanup_on_exit) return;

    buf_lock_(&vb->buffer_list);
    buflist_sort (vb, true); // need buffers to be sorted to erase the space between them
    
    uint32_t sizeof_vb = get_vb_size (vb->data_type); // note: buffers are expected in current data_type part of buffer, not in the unused portion due to an historical alloced_data_type
    
    ASSERT (sizeof_vb <= get_vb_size (vb->data_type_alloced), "Expecting vb_i=%u vb->id=%u data_type=%s to be of size >= %u, but it is allocated as %s of size %u",
            vb->vblock_i, vb->id, dt_name (vb->data_type), sizeof_vb, dt_name (vb->data_type_alloced), get_vb_size (vb->data_type_alloced));

    VBIType vb_i = vb->vblock_i; // save before it is erased

    char *start_erase = (char *)&vb->in_use + 1; // in_use is last in "fields that survive buflist_free_vb"

    for_buf2 (BufListEnt, ent, buf_i, vb->buffer_list) 
        if (!BL_IS_REMOVED(ent->buf) &&
            ent->buf != &vb->buffer_list && // don't free buffer_list
            // if evb, free only buffers contained in the VBlock structure
            (vb != evb || is_p_in_range (ent->buf, vb, sizeof_vb))) { 

            // test with next, not with freed previous
            ASSERT (buf_i+1 == vb->buffer_list.len32 || ent->buf != (ent+1)->buf, 
                    "Duplicate buffer_list entries in buf_i=%u vb_i=%u evb=%s ent1:{%s:%u %s} ent2:{%s:%u %s} for buf: %s%s", 
                    buf_i, vb_i, TF(vb==evb), 
                    ent->func, ent->code_line, ent->name, (ent+1)->func, (ent+1)->code_line, (ent+1)->name,
                    buf_desc(ent->buf).s, buflist_locate_s (ent->buf));
                
            char *save_buf = (char *)ent->buf; // before buf_destroy changes the ent->buf pointer

            if (ent->buf->shared && vb != evb)
                buf_destroy_do_do (ent, __FUNCLINE); 
            else
                buf_free (*ent->buf);

            // erase the space between the end of the previous buffer and the beginning of this one
            if (is_p_in_range (save_buf, vb, sizeof_vb)) {
                ASSERT (save_buf >= start_erase, "expecting save_buf=%p > start_erase=%p (vb_i=%d buf_i=%u/%u)", save_buf, start_erase, vb_i, buf_i, vb->buffer_list.len32);
                memset (start_erase, 0, (rom)save_buf - start_erase);
                start_erase = save_buf + sizeof(Buffer);
            }
        }

    // erase space after the last Buffer in the VB
    memset (start_erase, 0, (rom)vb + sizeof_vb - start_erase);

    buf_unlock;
}

// frees all buffers in a Context, and also erases the spaces between Buffers in the Context.
void buflist_free_ctx (VBlockP vb, ContextP ctx) 
{
    buf_lock_(&vb->buffer_list);

    buflist_sort (vb, true); // need buffers to be sorted to erase the space between them

    // find first buffer than belongs to ctx in vb->buffer_list
    char *start_erase = (char *)ctx; 

    BufListEnt *first_ent = MAX_(buflist_find_buf (vb, &ctx->FIRST_BUFFER_IN_Context, SOFT_FAIL), // binary search
                                 B1ST(BufListEnt, vb->buffer_list)); // if binary search returns NULL
    if (!first_ent) goto done;

    for (BufListEnt *ent = first_ent; (rom)ent->buf < (rom)(ctx+1); ent++)
        if (!BL_IS_REMOVED(ent->buf) && is_p_in_range (ent->buf, ctx, sizeof(Context))) { 

            char *save_buf = (char *)ent->buf; // before buf_destroy changes the ent->buf pointer

            if (ent->buf->shared && vb != evb)
                buf_destroy_do_do (ent, __FUNCLINE); 
            else
                buf_free (*ent->buf);

            // erase the space between the end of the previous buffer and the beginning of this one
            ASSERT (save_buf >= start_erase, "expecting ent->buf=%p > start_erase=%p", save_buf, start_erase);
            memset (start_erase, 0, save_buf - start_erase);
            start_erase = save_buf + sizeof (Buffer);
        }

    // erase space after the last Buffer in the VB
    memset (start_erase, 0, (rom)(ctx+1) - start_erase);

done:
    buf_unlock;
}

// destroys all buffers in the VB's buffer list - even if not contained in the VBlock structure
void buflist_destroy_vb_bufs (VBlockP vb, bool only_if_unused)
{
    if (flag.let_OS_cleanup_on_exit) return;
    
    buf_lock_(&vb->buffer_list);

    BufListEnt *buf_list_ent = NULL;
    for_buf (BufListEnt, ent, vb->buffer_list) {
        if (BL_IS_REMOVED(ent->buf) || 
            (only_if_unused && (ent->buf->data || ent->buf->len || ent->buf->param || ent->buf->promiscuous))) 
            continue;
    
        if (ent->buf != &vb->buffer_list) 
            buf_destroy_do_do (ent, __FUNCLINE);

        else 
            buf_list_ent = ent;
    }

    buf_unlock;

    buf_destroy_do_do (buf_list_ent, __FUNCLINE); // unlock before, bc this also destroys the buf_list lock
}

// destroy all private bufs of this VB, i.e. excluding buffers in the common area
// note: destroys REGULAR buffers, but leaves UNALLOCATED, SHM, STANDALONE_BITS buffers untouched
void buflist_destroy_private_and_context_vb_bufs (VBlockP vb)
{
    uint32_t sizeof_alloced_vb  = get_vb_size (vb->data_type_alloced); // note: data_type_alloced, not data_type! 
    uint32_t sizeof_common_area = (rom)(&vb->final_member + 1) - (rom)vb;
    buf_lock_(&vb->buffer_list);

    for_buf (BufListEnt, ent, vb->buffer_list) 
        if (!BL_IS_REMOVED(ent->buf) &&                             // not marked for removal
            is_p_in_range (ent->buf, vb, sizeof_alloced_vb) &&      // is in VB struct, but:
              (!is_p_in_range (ent->buf, vb, sizeof_common_area) || // either: not in common area
               is_p_in_range (ent->buf, vb->contexts, sizeof (ContextArray)) || // or: context data
               ent->buf == &vb->lines))                             // or: vb->lines
            buf_destroy_do_do (ent, __FUNCLINE);

    buf_unlock;
}

// main thread only: destroy all evb buffers with a specific name
// note: destroys REGULAR  buffers, but leaves UNALLOCATED, SHM, STANDALONE_BITS buffers untouched
void buflist_destroy_file_bufs (FileP file)
{
    if (flag.let_OS_cleanup_on_exit) return;
    
    ASSERTMAINTHREAD;

    buf_lock (&evb->buffer_list);

    buflist_sort (evb, true); // for effeciency (likely already sorted)

    for_buf (BufListEnt, ent, evb->buffer_list) 
        if (!BL_IS_REMOVED(ent->buf) && is_p_in_range (ent->buf, file, sizeof (File))) {
            // user_count=0 if promiscuous buffer was added to buffer list but never allocated, or 1 if buffer is allocated
            ASSERT (buf_user_count (ent->buf) <= 1, "%s has count=%u > 1", ent->buf->name, buf_user_count (ent->buf));
            buf_destroy_do_do (ent, __FUNCLINE);
        }

        else if ((rom)ent->buf > (rom)(file+1))
            break; 

    buf_unlock;

    buflist_compact (evb); // this locks for itself
}

// main thread only: destroy all evb buffers with a specific name
// note: destroys REGULAR buffers, but leaves UNALLOCATED, SHM, STANDALONE_BITS buffers untouched
void buflist_destroy_bufs_by_name (rom name, bool compact_after_destroying)
{
    if (flag.let_OS_cleanup_on_exit) return;

    ASSERTMAINTHREAD;

    buf_lock_(&evb->buffer_list);

    for_buf (BufListEnt, ent, evb->buffer_list) 
        if (!BL_IS_REMOVED(ent->buf) && ent->buf->name && !strcmp (ent->buf->name, name))
            buf_destroy_do_do (ent, __FUNCLINE);

    buf_unlock;

    if (compact_after_destroying) buflist_compact (evb);
}

//-----------------------------------------------
// testing integrity of buffers on buflist
//-----------------------------------------------

static void buflist_find_underflow_culprit (ConstVBlockP calling_vb, rom memory, rom msg)
{
    Buffer highest = {};
    rom highest_after = 0;
    uint64_t highest_of = 0;
    VBIType highest_vb_i;
    int highest_vb_id;
    uint32_t highest_buf_i;

    for (VBlockPoolType type=POOL_MAIN; type <= POOL_BGZF; type++) {
        VBlockPool *vb_pool = vb_get_pool (type, SOFT_FAIL);
        if (!vb_pool) continue;
    
        for (int vb_id=-1; vb_id < (int)vb_pool->num_vbs; vb_id++) {
            VBlockP vb = vb_get_from_pool (vb_pool, vb_id);
            if (!vb) continue;
            
            char func[50];
            sprintf (func, "%s(vb_i=%u)\n", __FUNCTION__, vb->vblock_i);

            buf_lock_if (&vb->buffer_list, vb != calling_vb); // calling_vb->buf_list is already locked by caller

            for_buf2 (BufListEnt, ent, i, vb->buffer_list) {
                if (!ent->buf || BL_IS_REMOVED(ent->buf)) continue;

                BufferP buf = ent->buf;         
                BufferSpinlock *spinlock = buf->promiscuous ? buf_lock_promiscuous (buf, __FUNCLINE) : NULL; // prevent frees or reallocs while we're testing
                if (buf->promiscuous && !spinlock) continue; // by the time we acquired the lock, buf was already freed
                
                // memory=0 could be a buffer that has been buf_moved, or a promiscous evb buffer current being realloced by a compute thread
                if (buf->memory) {
                    // find higher overflower that is lower than memory
                    rom after_buf = buf->memory + buf_mem_size (buf);

                    if (after_buf <= memory &&               // lower than us
                        highest_after > after_buf  &&        // highest so far
                        BOVERFLOW(buf) != OVERFLOW_TRAP &&   // overflowing
                        BUNDERFLOW(buf) == UNDERFLOW_TRAP) { // but not underflowing (overflowing AND underflowing is an indication that an lower buf completely overwrote this buffer)
                        
                        highest = *buf; // copy as it might change after unlocking
                        highest_of = BOVERFLOW(&highest);
                        highest_after = after_buf;
                        highest_vb_i = buf->vb->vblock_i;
                        highest_vb_id = buf->vb->id;
                        highest_buf_i = i;
                    }
                }

                buf_unlock; // promiscuous lock
            }

            buf_unlock; // buf_list lock
        }
    }

    if (highest_after) { // found
        rom of = (rom)&highest_of;
        fprintf (stderr,
                "Candidate culprit: vb->id=%d vblock_i=%d type=%s: buffer: %s memory: %p-%p buf_i=%u Overflow fence=%c%c%c%c%c%c%c%c\n",
                highest_vb_id, highest_vb_i, buf_type_name (&highest), buf_desc(&highest).s,
                highest.memory, highest.memory + buf_mem_size (&highest) - 1, highest_buf_i, 
                of[0], of[1], of[2], of[3], of[4], of[5], of[6], of[7]);
    }
    else 
        fprintf (stderr, "Cannot find a Buffer which has overflowed, and located just before the underflowed buffer\n\n");
}

// this function cannot contain ASSERT or ABORT as exit_on_error calls it
static bool buflist_test_overflows_do (VBlockP vb, bool primary, rom msg)
{
    if (!vb) return false;
    START_TIMER; 

    rom corruption = NULL;

    buf_lock_(&vb->buffer_list);

    BufferP buf; // declare outside the loop, so it is observable in the debugger in case of a crash
    uint32_t code_line;

    for_buf2 (BufListEnt, ent, buf_i, vb->buffer_list) {
        static rom nl[2] = {"", "\n\n"};
        #define fcn ent->func, ent->code_line, ent->name

        buf = ent->buf;
        func = ent->func;
        code_line = ent->code_line;

        // case: buf was 'buf_destroy'd
        if (BL_IS_REMOVED (ent->buf)) continue;  

        if (!buf) {
            // known issue: bug 912
            buflist_display (vb);
            fprintf (stderr, "%s: buf=NULL for buf_i=%u in vb->buffer_list ^^^. buffer_list=%s", VB_NAME, buf_i, buf_desc(&vb->buffer_list).s);
            corruption = "buffer_list has entry with buf=NULL";
            goto done;
        }

#ifdef WIN32
        if (IsBadReadPtr (ent->buf, sizeof (Buffer))) {
            fprintf (stderr, "%s%s: Memory corruption in vb->id=%d (vb->vblock_i=%d) buffer=%p func=%s:%u \"%s\" (buf_i=%u): buffer structure inaccessible (invalid pointer)\n", 
                     nl[primary], msg, vb->id, vb->vblock_i, ent->buf, fcn, buf_i);
            corruption = "buf has an invalid address";
        }
#endif

        BufferSpinlock *spinlock = buf->promiscuous ? buf_lock_promiscuous (buf, __FUNCLINE) : NULL; // prevent frees or reallocs while we're testing
        if (buf->promiscuous && !spinlock) return false; // by the time we acquired the lock, buf was already freed

        // memory=0 could be a buffer that has been buf_moved, or a promiscous evb buffer current being 
        // realloced by a compute thread (for example, Context.global_hash is realloced by all threads (under mutex protection)
        if (!buf->memory) goto done_one_buf;

        if (vb && buf->vb != vb) {
            fprintf (stderr, "%s%s: Memory corruption in vb->id=%d (vb->vblock_i=%d) buffer=%p func=%s:%u \"%s\" (buf_i=%u): Corrupt Buffer structure OR invalid buffer pointer - buf->vb=%p != vb=%p\n", 
                     nl[primary], msg, vb->id, vb->vblock_i, buf, fcn, buf_i, buf->vb, vb);
            buflist_locate (buf, "Corrupt Buffer");
            corruption = "VB mismatch";
        }
        else if (buf->data && buf->vb->vblock_i != vb->vblock_i) { // buffers might still be here from the previous incarnation of this vb - its ok if they're not allocated yet
                    fprintf (stderr, "%s%s: Memory corruption in vb_id=%d: buf_vb_i=%d differs from thread_vb_i=%d: buffer: %s %p func: %s:%u \"%s\" memory: %p-%p name: %s vb_i=%u buf_i=%u\n",
                             nl[primary], msg, vb ? vb->id : -999, buf->vb->vblock_i, vb->vblock_i, buf_type_name(buf), 
                             buf, fcn, buf->memory, buf->memory + buf_mem_size (buf)-1,
                             buf_desc (buf).s, buf->vb->vblock_i, buf_i);
            corruption = "vblock_i mismatch";
        }
        else if (buf->type < 0 || buf->type > BUF_NUM_TYPES) {
            fprintf (stderr, "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d) buffer=%p func=%s:%u \"%s\" (buf_i=%u): Corrupt Buffer structure OR invalid buffer pointer - invalid buf->type", 
                        nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, buf, fcn, buf_i);
            fprintf (stderr, " Buffer=%s\n", buf_desc(buf).s);  // separate fprintf in case it seg faults
            corruption = "invalid type";
        }
        else if (!buf->name) {
            fprintf (stderr, "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d): buffer=%p func=%s:%u \"%s\" (buf_i=%u): Corrupt Buffer structure - null name", 
                     nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, buf, fcn, buf_i);
            fprintf (stderr, " Buffer=%s\n", buf_desc(buf).s);  // separate fprintf in case it seg faults
            corruption = "missing name";
        }
        else if (BUNDERFLOW(buf) != UNDERFLOW_TRAP) {
            fprintf (stderr, 
                     "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d): Underflow: buffer: %s %p func: %s:%u memory: %p-%p name: %s vb_i=%u buf_i=%u. Fence=%c%c%c%c%c%c%c%c\n",
                     nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, buf_type_name(buf), buf, func, code_line, buf->memory, buf->memory+buf_mem_size(buf)-1, 
                     buf_desc (buf).s, buf->vb->vblock_i, buf_i, 
                     buf->memory[0], buf->memory[1], buf->memory[2], buf->memory[3], buf->memory[4], buf->memory[5], buf->memory[6], buf->memory[7]);

            buf_unlock;
            buflist_find_underflow_culprit (vb, buf->memory, msg);
            buf_lock (&vb->buffer_list); // relock

            if (primary) buflist_test_overflows_all_other_vb (vb, msg, false);
            primary = false;

            corruption = "underflow";
        }
        else if (BOVERFLOW(buf) != OVERFLOW_TRAP) {
            char *of = &buf->memory[buf->size + sizeof(uint64_t)];
            fprintf (stderr,
                     "%s%s: Memory corruption in vb_id=%d (vb_i=%d): Overflow: buffer: %s %p func: %s:%u memory: %p-%p name: %s vb_i=%u buf_i=%u Fence=%c%c%c%c%c%c%c%c\n",
                     nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, buf_type_name(buf), buf, func, code_line, buf->memory, buf->memory+buf_mem_size(buf)-1, 
                     buf_desc (buf).s, buf->vb->vblock_i, buf_i, of[0], of[1], of[2], of[3], of[4], of[5], of[6], of[7]);
            
            if (primary) buflist_test_overflows_all_other_vb (vb, msg, false);
            primary = false;

            corruption = "overflow";
        }

    done_one_buf:
        if (buf) buf_unlock; // promiscuous lock
    }
    
    done: 
    buf_unlock; // buf_list lock

    ASSERT (!primary || !corruption, "Aborting due to memory corruption \"%s\"", corruption); // primary will exit on corruption

    COPY_TIMER (buflist_test_overflows_do);
    return corruption > 0;
}

void buflist_test_overflows_all_other_vb (VBlockP caller_vb, rom msg, bool force)
{
    // IMPORTANT: this function is not thread safe - while buflist_test_overflows_do
    // locks the buffer_list of each VB, the individual buffers are not locked
    // unless shared or promiscuous - and may be destroyed or realloced concurrently
    // with the test causing a segfault. That's why the "return" is normally here. (bug 828)
    if (!force && !flag.debug_memory && !flag.xthreads) 
        return; 

    fprintf (stderr, "\nTesting all other VBs (WARNING: NOT thread safe - might segfault; activated by certain flags (see code)):\n");
    bool corruption_detected = false;
    for (VBlockPoolType type=POOL_MAIN; type <= POOL_BGZF; type++) {
        VBlockPool *vb_pool = vb_get_pool (type, SOFT_FAIL);
        if (!vb_pool) continue;

        for (int vb_id=-1; vb_id < (int)vb_pool->num_vbs; vb_id++) {
            VBlockP vb = vb_get_from_pool (vb_pool, vb_id);
            if (!vb || vb == caller_vb) continue; // skip caller's VB
            corruption_detected |= buflist_test_overflows_do (vb, false, msg);
        }
    }

    if (!corruption_detected) fprintf (stderr, "No issues found in other VBs\n");
}

// this function cannot contain ASSERT or ABORT as exit_on_error calls it
void buflist_test_overflows_all_vbs (rom msg)
{
    buflist_test_overflows_all_other_vb (NULL, msg, false);
}

// this function cannot contain ASSERT or ABORT as exit_on_error calls it
bool buflist_test_overflows (VBlockP vb, rom msg)
{
    return buflist_test_overflows_do (vb, true, msg); // true if corruption detected
}

//-----------------------------------------------
// reporting memory usage
//-----------------------------------------------

static DESCENDING_SORTER (buf_stats_sort_by_bytes, MemStats, bytes)

static bool buf_count_mem_usage (ConstBufferP buf, VBlockPoolType pool_type, int vb_id, void *mem_usage)
{
    *((uint64_t *)mem_usage) += buf_mem_size (buf);

    return true; // continue
}

uint64_t buflist_get_memory_usage (void)
{
    uint64_t mem_usage = 0;
    buflist_foreach_buffer_in_each_vbs (buf_count_mem_usage, &mem_usage);
    return mem_usage;
}

#define MAX_MEMORY_STATS 150
static MemStats stats[MAX_MEMORY_STATS]; // must be pre-allocated, because buflist_show_memory is called when malloc fails, so it cannot malloc
static unsigned num_stats=0, num_buffers=0;

static bool buflist_show_memory_add_buf (ConstBufferP buf, VBlockPoolType pool_type, int vb_id, void *unused)
{
    ASSERTW (buf->name && strlen (buf->name) > 0, "FYI: buffer allocated in %s:%u has no name", buf->func, buf->code_line);

    bool found = false;
    for (unsigned st_i=0; st_i < num_stats && !found; st_i++) 
        if (!strcmp (stats[st_i].name, buf->name)) {
            if (!buf->shared || buf->vb == evb) { // for a shared buf, increment only in evb so we don't double-count
                stats[st_i].buffers++;
                stats[st_i].bytes += buf_mem_size (buf);
            }
            found = true;
        }

    if (!found) {
        stats[num_stats++] = (MemStats){ .buffers = 1, 
                                         .name    = buf->name, 
                                         .bytes   = buf->size + CTL_SIZE };
        ASSERT (num_stats < MAX_MEMORY_STATS, "# memory stats exceeded %u, consider increasing MAX_MEMORY_STATS", MAX_MEMORY_STATS);
    }

    num_buffers++;

    return true; // continue
}

void buflist_show_memory (bool memory_full, unsigned max_threads, unsigned used_threads)
{
    // if memory is full - only one thread needs to show this - other threads stall
    static bool once=0;
    if (memory_full && __atomic_test_and_set (&once, __ATOMIC_RELAXED)) 
        while (1) usleep (1000000000);
    
    memset (stats, 0, sizeof(stats));
    num_stats = num_buffers = 0;
    buflist_foreach_buffer_in_each_vbs (buflist_show_memory_add_buf, 0);
    
    if (memory_full)
        fprintf (stderr, "\n\nError memory is full:\n");
    else
        iprint0 ("\n-------------------------------------------------------------------------------------\n");

    // sort stats by bytes
    qsort (stats, num_stats, sizeof (MemStats), buf_stats_sort_by_bytes);

    uint64_t total_bytes=0;
    for (unsigned i=0; i< num_stats; i++) total_bytes += stats[i].bytes;

    uint32_t num_allocated_vbs = (vb_get_pool (POOL_MAIN, SOFT_FAIL) ? vb_get_pool (POOL_MAIN, HARD_FAIL)->num_allocated_vbs : 0)
                               + (vb_get_pool (POOL_BGZF, SOFT_FAIL) ? vb_get_pool (POOL_BGZF, HARD_FAIL)->num_allocated_vbs : 0);

    fprintf (memory_full ? stderr : info_stream, "Total bytes: %s in %u buffers in %u buffer lists. global_max_threads=%u max_resident_size=%s\n", 
             str_size (total_bytes).s, num_buffers, num_allocated_vbs, global_max_threads, str_size (arch_get_max_resident_set()).s);
    if (IS_ZIP) 
        fprintf (memory_full ? stderr : info_stream, "vb_size = %u MB\n", (unsigned)(segconf.vb_size >> 20));
    
    if (max_threads)
        fprintf (memory_full ? stderr : info_stream, "Compute threads: max_permitted=%u actually_used=%u\n", max_threads, used_threads);

    for (unsigned i=0; i < num_stats; i++)
        fprintf (memory_full ? stderr : info_stream, "%-30s: %-8s (%4.1f%%) in %u buffers\n", stats[i].name, str_size (stats[i].bytes).s, 100.0 * (float)stats[i].bytes / (float)total_bytes, stats[i].buffers);

    fprintf (memory_full ? stderr : info_stream, "\n");
}

