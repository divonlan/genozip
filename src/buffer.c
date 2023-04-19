// ------------------------------------------------------------------
//   buffer.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// memory management - when running the same code by the same thread for another VB - we reuse
// the previous variant's block memory. this way we save repetitive malloc/free cycles which might
// be very time consuming.

#ifndef _WIN32
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#else
#include <windows.h>
#endif
#include <fcntl.h> 

#include "genozip.h"
#include "profiler.h"
#include "buffer.h"
#include "vblock.h"
#include "strings.h"
#include "reference.h"
#include "bits.h"
#include "mutex.h"
#include "file.h"
#include "endianness.h"
#include "segconf.h"
#include "threads.h"
#include "gencomp.h"

#define DISPLAY_ALLOCS_AFTER 0 // display allocations, except the first X allocations. reallocs are always displayed

#define UNDERFLOW_TRAP 0x574F4C4652444E55ULL // "UNDRFLOW" - inserted at the begining of each memory block to detected underflows
#define OVERFLOW_TRAP  0x776F6C667265766FULL // "overflow" - inserted at the end of each memory block to detected overflows

#define BUFFER_BEING_MODIFIED ((char*)0x777)

// We mark a buffer list entry as removed, by setting its LSb. This keeps the buffer list sorted, in case it is sorted.
// Normally the LSb is 0 as buffers are word-aligned (verified in buf_add_to_buffer_list_do)
#define BL_SET_REMOVED(bl_ent) bl_ent = ((BufferP)((uint64_t)bl_ent | 1))
#define BL_IS_REMOVED(bl_ent)  ((uint64_t)bl_ent & 1) 

static const unsigned control_size = 2*sizeof (uint64_t) + sizeof(uint16_t); // underflow, overflow and user counter

static Mutex overlay_mutex   = {}; // used to thread-protect overlay counters (note: not initializing here - different in different OSes)
static Mutex only_once_mutex = {};
static uint64_t abandoned_mem_current = 0;
static uint64_t abandoned_mem_high_watermark = 0;

static bool cleanup_on_exit = false; // if true, the memory release is left for the OS process exit() and not done explicitly

#ifdef _WIN32
static HANDLE heap;
#endif

typedef struct {
    BufferP buf;  // address of Buffer structure (not the data!). Buffer is removed if removed bit is set with BL_SET_REMOVED 
    rom func;     // function in which this Buffer was first allocated (at which time it was added to the buffer_list)
    rom name;     // buffer name
    uint32_t code_line;
} BufListEnt;

static inline bool is_promiscuous (ConstBufferP buf)
{
    return buf->promiscuous && (global_max_threads > 1); // if --threads=1, even promiscuous buffers cannot behave promiscuously 
}

void buf_initialize()
{
    mutex_initialize (overlay_mutex);
    mutex_initialize (only_once_mutex);

#ifdef _WIN32
    heap = GetProcessHeap();
    ASSERT (heap, "GetProcessHeap failed: %s", str_win_error());
#endif
}

// note that the system is cleanup up before an exit(). this is irreversible
void buf_set_cleanup_on_exit (void)
{
    __atomic_store_n (&cleanup_on_exit, (bool)true, __ATOMIC_SEQ_CST);
}

static rom buf_display_type (ConstBufferP buf)
{
    static rom names[] = BUFTYPE_NAMES;
    if (buf->type >= 0 && buf->type < BUF_NUM_TYPES) 
        return names[buf->type];
    else
        return "INVALID";
}

// get string with buffer's metadata for debug message. this function is NOT thread-safe
char *buf_display (ConstBufferP buf)
{
    static char str[256]; // NOT thread-safe

    sprintf (str, "Buffer %s (%"PRId64"): size=%"PRIu64" len=%"PRIu64" data=%p memory=%p",
             buf->name, buf->param, (uint64_t)buf->size, buf->len, buf->data, buf->memory);
    return str;    
}

const BufDescType buf_desc (ConstBufferP buf)
{
    #define IS_CONTEXT(ctxs) ((rom)buf >= (rom)ctxs && (rom)buf < (rom)&ctxs[MAX_DICTS])
    #define TAG_NAME(ctxs) ctxs[((rom)buf - (rom)ctxs) / (sizeof (ctxs) / MAX_DICTS /*note: may be different than sizeof(Context) due to word alignment*/)].tag_name

    // case: buffer is one of the Buffers within a Context - show tag_name
    rom tag_name = NULL;
    if (buf->vb && buf->vb != evb && IS_CONTEXT(buf->vb->contexts))
        tag_name = TAG_NAME (buf->vb->contexts);
    else if (buf->vb && buf->vb == evb && IS_CONTEXT(z_file->contexts))
        tag_name = TAG_NAME (z_file->contexts);

    BufDescType desc; // use static memory instead of malloc since we could be in the midst of a memory issue when this is called
    sprintf (desc.s, "\"%s\"%s%s param=%"PRId64" len=%"PRIu64" size=%"PRId64" type=%s allocated in %s:%u by vb=%d", 
             buf->name ? buf->name : "(no name)", 
             tag_name ? " ctx=" : "", tag_name ? tag_name : "",
             buf->param, buf->len, (uint64_t)buf->size, buf_display_type (buf), buf->func ? buf->func : "(no func)", buf->code_line, 
             (buf->vb ? buf->vb->vblock_i : -999));
    return desc;
}

static inline void buf_reset (BufferP buf)
{
    VBlockP save_vb = buf->vb;
    memset (buf, 0, sizeof (Buffer)); // make this buffer UNALLOCATED

    buf->vb = save_vb;
}

static void buf_foreach_buffer_in_vb (VBlockP vb, bool (*callback)(ConstBufferP, VBlockPoolType pool_type, int vb_id, void *arg), VBlockPoolType pool_type, int vb_id, void *arg)
{
    for_buf2 (BufListEnt, bl, buf_i, vb->buffer_list) {
        #ifdef __linux__
            // if (!BL_IS_REMOVED (*bl) && access ((rom)*bl, F_OK) == -1 && errno == EFAULT)
            //     ABORT ("buffer structure inaccessible (invalid pointer=%p) in buf_i=%u of vb_id=%d", *bl, buf_i, vb_id); 
        #elif defined WIN32
            if (!BL_IS_REMOVED (bl->buf) && IsBadReadPtr (bl->buf, sizeof (Buffer))) 
                ABORT ("buffer structure inaccessible (invalid pointer=%p) in buf_i=%u of vb_id=%d", bl->buf, buf_i, vb->id); 
        #endif

        if (BL_IS_REMOVED (bl->buf) || !bl->buf->memory) continue; // exclude destroyed, not-yet-allocated, overlay buffers and buffers that were src in buf_move

        if (!callback (bl->buf, pool_type, vb_id, arg))
            break;
    }
}

static void buf_foreach_buffer (bool (*callback)(ConstBufferP, VBlockPoolType pool_type, int vb_id, void *arg), void *arg)
{
    for (VBlockPoolType pool_type=0; pool_type < NUM_POOL_TYPES; pool_type++) {

        VBlockPool *vb_pool = vb_get_pool (pool_type, SOFT_FAIL);
        if (!vb_pool) continue;

        // note: we don't cover EVB (only in DEBUG and --show-mem) as it segfaults for an unknown reason (likely the Buffer structure itself resides with a structure
        // that is freed (eg File)). TODO: debug this.
        for (int vb_id=0; vb_id < (int)vb_pool->num_allocated_vbs; vb_id++) {

            VBlockP vb = vb_get_from_pool (vb_pool, vb_id);
            if (!vb) continue;

            buf_foreach_buffer_in_vb (vb, callback, pool_type, vb_id, arg);
        }
    }

    // non-pool VBs: 
    // TO DO: add the other non-pool VBs (txt_header_vb...)
    // note: we don't cover EVB (only in DEBUG and --show-mem) as it segfaults for an unknown reason (likely the Buffer structure itself resides with a structure
    // that is freed (eg File)). TODO: debug this.
    if (flag.debug || flag.show_memory)          
        buf_foreach_buffer_in_vb (evb, callback, -1, VB_ID_EVB, arg);
}

static inline bool buf_has_overflowed (ConstBufferP buf, rom msg)
{
    // note on promiscuous evb buffers: if another thread reallocs the memory concurrently, the test (memory != OVERFLOW_TRAP) might seg-fault
    if (is_promiscuous (buf) || buf->type != BUF_REGULAR) return false;

    // memory==BUFFER_BEING_MODIFIED if an evb buffer is currently being allocated by another thread, and hence memory is not set yet.
    // for example, global_hash_prime is realloced by all threads (under mutex protection)
    // note: memory can become 0x777 at this point, even though we've tested it before at it was not yet
    char *memory = buf->memory;
    if (buf->vb && buf->vb->id==-1 && memory == BUFFER_BEING_MODIFIED) return false;

    ASSERT (memory != BUFFER_BEING_MODIFIED, "%s: buf->memory=BUFFER_BEING_MODIFIED. buffer %s size=%"PRIu64,
            msg, buf_desc(buf).s, (uint64_t)buf->size);

    // note on promiscuous evb buffers: if another thread reallocs the memory concurrently, this might seg-fault
    return *((uint64_t*)(memory + buf->size + sizeof(uint64_t))) != OVERFLOW_TRAP; 
}

static inline bool buf_has_underflowed (ConstBufferP buf, rom msg)
{
    // note: overlayed buffers might be partial and start at a different address 
    // note on promiscuous evb buffers: if another thread reallocs the memory concurrently, the test (memory != UNDERFLOW_TRAP) might seg-fault
    if (is_promiscuous (buf) || buf->type != BUF_REGULAR) return false; 

    // see comment in buf_has_overflowed
    char *memory = buf->memory;
    if (buf->vb && buf->vb->id==-1 && memory == BUFFER_BEING_MODIFIED) return false;

    ASSERT (memory != BUFFER_BEING_MODIFIED, "%s: buf->memory=BUFFER_BEING_MODIFIED. buffer %s size=%"PRIu64,
            msg, buf_desc(buf).s, (uint64_t)buf->size);

    return *(uint64_t*)memory != UNDERFLOW_TRAP; 
}

// not thread-safe, used in emergency 
static void buf_find_underflow_culprit (rom memory, rom msg)
{
    bool found=false;
    for (VBlockPoolType type=POOL_MAIN; type <= POOL_BGZF; type++) {
        VBlockPool *vb_pool = vb_get_pool (type, SOFT_FAIL);
        if (!vb_pool) continue;
    
        for (int vb_id=-1; vb_id < (int)vb_pool->num_vbs; vb_id++) {
            VBlockP vb = vb_get_from_pool (vb_pool, vb_id);
            if (!vb) continue;
            
            for_buf2 (BufferP , buf_p, i, vb->buffer_list) {
                BufferP buf = *buf_p;
                if (buf) {
                    char *after_buf = buf->memory + buf->size + control_size;
                    if (after_buf <= memory && (after_buf + 100 > memory) && buf_has_overflowed (buf, msg)) {
                        char *of = &buf->memory[buf->size + sizeof(uint64_t)];
                        fprintf (stderr,
                                "Candidate culprit: vb->id=%d (vb->vblock_i=%d): buffer: %s %p memory: %p-%p name: %s vb_id=%u buf_i=%u Overflow fence=%c%c%c%c%c%c%c%c\n",
                                vb ? vb->id : -999, vb->vblock_i, buf_display_type(buf), 
                                buf, buf->memory, after_buf-1, 
                                buf_desc (buf).s, buf->vb->vblock_i, i, 
                                of[0], of[1], of[2], of[3], of[4], of[5], of[6], of[7]);
                        found = true;
                    }
                }
            }
        }
    }

    if (!found) fprintf (stderr, "Cannot find a Buffer which has overflowed, and located just before the underflowed buffer\n\n");
}

// this function cannot contain ASSERT or ABORT as exit_on_error calls it
static bool buf_test_overflows_do (ConstVBlockP vb, bool primary, rom msg);
static void buf_test_overflows_all_other_vb(ConstVBlockP caller_vb, rom msg, bool force)
{
    // IMPORTANT: this function is not thread safe as it checks OTHER thread's memories which
    // may be dealloced as it checking them causing weird errors.
    // That's why the "return" is normally here.
    if (!force) return; 

    fprintf (stderr, "Testing all other VBs (WARNING: NOT thread safe - might cause all kinds of unrelated errors, activated by --debug or GENOZIP_TEST):\n");
    bool corruption_detected = false;
    for (VBlockPoolType type=POOL_MAIN; type <= POOL_BGZF; type++) {
        VBlockPool *vb_pool = vb_get_pool (type, SOFT_FAIL);
        if (!vb_pool) continue;

        for (int vb_id=-1; vb_id < (int)vb_pool->num_vbs; vb_id++) {
            VBlockP vb = vb_get_from_pool (vb_pool, vb_id);
            if (!vb || vb == caller_vb) continue; // skip caller's VB
            corruption_detected |= buf_test_overflows_do (vb, false, msg);
        }
    }

    if (!corruption_detected) fprintf (stderr, "No issues found in other VBs\n");
}

static bool buf_locate (ConstBufferP buf, bool verbose)
{
    for (VBlockPoolType pool_type=0; pool_type < NUM_POOL_TYPES; pool_type++) {

        VBlockPool *vb_pool = vb_get_pool (pool_type, SOFT_FAIL);
        if (!vb_pool) continue;

        // note: we don't cover EVB (only in DEBUG and --show-mem) as it segfaults for an unknown reason (likely the Buffer structure itself resides with a structure
        // that is freed (eg File)). TODO: debug this.
        for (int vb_id=0; vb_id < (int)vb_pool->num_allocated_vbs; vb_id++) {

            VBlockP vb = vb_get_from_pool (vb_pool, vb_id);
            if (!vb) continue;

            if (vb_buf_locate (vb, buf)) {
                if (verbose) 
                    fprintf (stderr, "Corrupt Buffer %p located in pool=%s vb_id=%d: %s\n", 
                             buf, pool_name(pool_type), vb_id, buf_desc (buf).s);
                return true;
            }
        }
    }

    #define OBJECT_TEST(locator, obj)                                                                   \
        if (locator (obj, buf)) {                                                                \
            if (verbose)                                                                                \
                fprintf (stderr, "Corrupt Buffer %p located in " #obj ": %s\n", buf, buf_desc (buf).s); \
            return true;                                                                                \
        }
    
    extern VBlockP compress_depn_vb, scan_vb, segconf_vb;
    void *depn=0, *componentsP=0, *queueP=0; // dummies
    OBJECT_TEST(vb_buf_locate, evb);
    OBJECT_TEST(vb_buf_locate, compress_depn_vb);
    OBJECT_TEST(vb_buf_locate, scan_vb);
    OBJECT_TEST(vb_buf_locate, segconf_vb);
    OBJECT_TEST(file_buf_locate, z_file);
    OBJECT_TEST(file_buf_locate, txt_file);
    OBJECT_TEST(ref_buf_locate, gref);
    OBJECT_TEST(ref_buf_locate, prim_ref);
    OBJECT_TEST(gencomp_buf_locate_depn, depn);
    OBJECT_TEST(gencomp_buf_locate_componentsP, componentsP);
    OBJECT_TEST(gencomp_buf_locate_queueP, queueP);

    if (verbose) 
        fprintf (stderr, "Corrupt Buffer %p is not located in any of the objects tested.", buf);

    return false;
}

// this function cannot contain ASSERT or ABORT as exit_on_error calls it
static bool buf_test_overflows_do (ConstVBlockP vb, bool primary, rom msg)
{
    if (!vb) return false;

    int corruption = 0;

    ConstBufferP buf; // declare outside the loop, so it is observable in the debugger in case of a crash
    rom func, name;
    uint32_t code_line;

    for_buf2 (BufListEnt, ent, buf_i, vb->buffer_list) {
        static rom nl[2] = {"", "\n\n"};

        buf       = ent->buf;
        func      = ent->func;
        name      = ent->name;
        code_line = ent->code_line;

        if (BL_IS_REMOVED (buf) ||  // buf was 'buf_destroy'd
            buf->type == BUF_OVERLAY) continue; 

#ifdef WIN32
        if (IsBadReadPtr (buf, sizeof (Buffer))) {
            fprintf (stderr, "%s%s: Memory corruption in vb->id=%d (vb->vblock_i=%d) buffer=%p func=%s:%u \"%s\" (buf_i=%u): buffer structure inaccessible (invalid pointer)\n", 
                     nl[primary], msg, vb->id, vb->vblock_i, buf, func, code_line, name, buf_i);
            corruption = 100;
            break;
        }
#endif

        // IMPORTANT NOTE regarding promiscuous evb buffers: testing promiscuous buffers (only evb might have them) might FAIL and should 
        // not be done in production! this if another thread can modify its buffers (under mutex protection) concurrently with this test, 
        // causing inconsistent state between eg data and memory we attempt to prevent many of the cases by not checking buffers that are 
        // BUFFER_BEING_MODIFIED at the onset, but they still may become BUFFER_BEING_MODIFIED mid way through the test
        if (is_promiscuous (buf)) continue;

        if (buf->memory && buf->memory != BUFFER_BEING_MODIFIED) {

            if (vb && buf->vb != vb) {
                fprintf (stderr, "%s%s: Memory corruption in vb->id=%d (vb->vblock_i=%d) buffer=%p func=%s:%u \"%s\" (buf_i=%u): Corrupt Buffer structure OR invalid buffer pointer - buf->vb=%p != vb=%p\n", 
                         nl[primary], msg, vb->id, vb->vblock_i, buf, func, code_line, name, buf_i, buf->vb, vb);
                buf_locate (buf, true);
                corruption = 1;
            }
            else if (buf->data && buf->vb->vblock_i != vb->vblock_i) { // buffers might still be here from the previous incarnation of this vb - its ok if they're not allocated yet
                        fprintf (stderr, "%s%s: Memory corruption in vb_id=%d: buf_vb_i=%d differs from thread_vb_i=%d: buffer: %s %p func: %s:%u \"%s\" memory: %p-%p name: %s vb_i=%u buf_i=%u\n",
                        nl[primary], msg, vb ? vb->id : -999, buf->vb->vblock_i, vb->vblock_i, buf_display_type(buf), 
                        buf, func, code_line, name, buf->memory, buf->memory+buf->size+control_size-1,
                        buf_desc (buf).s, buf->vb->vblock_i, buf_i);
                corruption = 2;
            }
            else if (buf->type < 0 || buf->type > BUF_NUM_TYPES) {
                fprintf (stderr, "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d) buffer=%p func=%s:%u \"%s\" (buf_i=%u): Corrupt Buffer structure OR invalid buffer pointer - invalid buf->type", 
                         nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, buf, func, code_line, name, buf_i);
                fprintf (stderr, " Buffer=%s\n", buf_desc(buf).s);  // separate fprintf in case it seg faults
                corruption = 3;
            }
            else if (!buf->name) {
                fprintf (stderr, "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d): buffer=%p func=%s:%u \"%s\" (buf_i=%u): Corrupt Buffer structure - null name", 
                         nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, buf, func, code_line, name, buf_i);
                fprintf (stderr, " Buffer=%s\n", buf_desc(buf).s);  // separate fprintf in case it seg faults
                corruption = 4;
            }
            else if (buf->data && buf->type == BUF_REGULAR && (buf->data != buf->memory + sizeof(uint64_t))) {
                fprintf (stderr, 
                         "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d): data!=memory+8: allocating_vb_i=%u buf_i=%u buffer=%p func=%s:%u \"%s\" memory=%p name=%s : Corrupt Buffer structure - expecting data+8 == memory. buf->data=%p\n", 
                         nl[primary], msg, vb ? vb->id : -999, vb->vblock_i,  buf->vb->vblock_i, buf_i, buf, func, code_line, name, buf->memory, buf_desc(buf).s, buf->data);
                corruption = 5;
            }
            else if (buf_has_underflowed(buf, msg)) {
                fprintf (stderr, 
                        "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d): Underflow: buffer: %s %p func: %s:%u memory: %p-%p name: %s vb_i=%u buf_i=%u. Fence=%c%c%c%c%c%c%c%c\n",
                        nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, buf_display_type(buf), buf, func, code_line, buf->memory, buf->memory+buf->size+control_size-1, 
                        buf_desc (buf).s, buf->vb->vblock_i, buf_i, 
                        buf->memory[0], buf->memory[1], buf->memory[2], buf->memory[3], buf->memory[4], buf->memory[5], buf->memory[6], buf->memory[7]);

                buf_find_underflow_culprit (buf->memory, msg);

                if (primary) buf_test_overflows_all_other_vb (vb, msg, false);
                primary = false;

                corruption = 6;
            }
            else if (buf_has_overflowed(buf, msg)) {
                char *of = &buf->memory[buf->size + sizeof(uint64_t)];
                fprintf (stderr,
                        "%s%s: Memory corruption in vb_id=%d (vb_i=%d): Overflow: buffer: %s %p func: %s:%u memory: %p-%p name: %s vb_i=%u buf_i=%u Fence=%c%c%c%c%c%c%c%c\n",
                        nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, buf_display_type(buf), buf, func, code_line, buf->memory, buf->memory+buf->size+control_size-1, 
                        buf_desc (buf).s, buf->vb->vblock_i, buf_i, of[0], of[1], of[2], of[3], of[4], of[5], of[6], of[7]);
                
                if (primary) buf_test_overflows_all_other_vb (vb, msg, false);
                primary = false;

                corruption = 7;
            }
        }
    }
    
    ASSERT (!primary || !corruption, "Aborting due to memory corruption #%u", corruption); // primary will exit on corruption

    return corruption > 0;
}

// this function cannot contain ASSERT or ABORT as exit_on_error calls it
bool buf_test_overflows (void *vb, rom msg)
{
    return buf_test_overflows_do ((ConstVBlockP)vb, true, msg); // true if corruption detected
}

// this function cannot contain ASSERT or ABORT as exit_on_error calls it
void buf_test_overflows_all_vbs (rom msg)
{
    buf_test_overflows_all_other_vb (NULL, msg, false);
}

// quick inline for internal buffer.c use check overflow and underflow in an allocated buffer
static inline void buf_verify_integrity (ConstBufferP buf, FUNCLINE, rom buf_func)
{
    // note on promiscuous evb buffers: if another thread reallocs the memory concurrently, the test (memory == OVER/UNDERFLOW_TRAP) might seg-fault
    if (is_promiscuous (buf)) return; 

    // note: // overlayed buffers might be partial and start at a different address 
    if (buf->type == BUF_REGULAR) {
        ASSERTGOTO (*(uint64_t *)(buf->memory) == UNDERFLOW_TRAP, "called from %s:%u to %s: Error in %s: buffer has corrupt underflow trap",
                    func, code_line, buf_func, buf_desc(buf).s);

        ASSERTGOTO (buf->memory + 8 == buf->data || !buf->data, "called from %s:%u to %s: expecting memory=%p + 8 == data=%p", 
                    func, code_line, buf_func, buf->memory, buf->data);
    }

    ASSERTGOTO (*(uint64_t *)(buf->memory + sizeof(uint64_t) + buf->size) == OVERFLOW_TRAP, "called from %s:%u to %s: Error in %s: buffer has corrupt overflow trap",
                func, code_line, buf_func, buf_desc(buf).s);
    
    return; // all good

error: {
    bool corruption_detected = buf_test_overflows (buf->vb, buf_func);
    if (corruption_detected) buf_test_overflows_all_other_vb (buf->vb, buf_func, true); // corruption not from this VB - test the others
    exit_on_error (false);
}}

// called from other modules for debugging memory issues. 
void buf_verify (ConstBufferP buf, rom msg)
{
    if (!buf || !buf->memory) return;

    DO_ONCE if (is_promiscuous (buf) && !threads_am_i_main_thread())
        fprintf (stderr, "FYI: buf_verify on promiscuous buffer \"%s\" might segfault. This is not an indication of an error.\n", buf->name ? buf->name : "(none)");

    bool save_promiscuous = ((BufferP)buf)->promiscuous;
    ((BufferP)buf)->promiscuous = false;

    buf_verify_integrity (buf, __FUNCLINE, msg);

    ((BufferP)buf)->promiscuous = save_promiscuous;
}

static DESCENDING_SORTER (buf_stats_sort_by_bytes, MemStats, bytes)

void buf_update_buf_list_vb_addr_change (VBlockP new_vb, VBlockP old_vb)
{
    if (old_vb == new_vb) return; // no change in address
        
    for_buf (BufListEnt, bl, new_vb->buffer_list) {

        if (BL_IS_REMOVED (bl->buf)) continue;

        // if this buffer is within the moved VB - update its address to the same offset relative to the new vb address
        if ((char*)bl->buf >= (char*)old_vb && (char*)bl->buf < (char *)(old_vb+1)) 
            bl->buf = (BufferP)((char*)new_vb + ((char*)bl->buf - (char*)old_vb)); 

        // update VB address within the buffer
        (bl->buf)->vb = new_vb; 
    }
}

static bool buf_count_mem_usage (ConstBufferP buf, VBlockPoolType pool_type, int vb_id, void *mem_usage)
{
    *((uint64_t *)mem_usage) += buf->size + control_size;

    return true; // continue
}

uint64_t buf_get_memory_usage (void)
{
    uint64_t mem_usage = 0;
    buf_foreach_buffer (buf_count_mem_usage, &mem_usage);
    return mem_usage;
}

#define MAX_MEMORY_STATS 100
static MemStats stats[MAX_MEMORY_STATS]; // must be pre-allocated, because buf_show_memory is called when malloc fails, so it cannot malloc
static unsigned num_stats=0, num_buffers=0;

static bool buf_add_mem_usage_to_stats (ConstBufferP buf, VBlockPoolType pool_type, int vb_id, void *unused)
{
    ASSERTW (buf->name && strlen (buf->name) > 0, "FYI: buffer allocated in %s:%u has no name", buf->func, buf->code_line);

    bool found = false;
    for (unsigned st_i=0; st_i < num_stats && !found; st_i++) 
        if (!strcmp (stats[st_i].name, buf->name)) {
            stats[st_i].buffers++;
            stats[st_i].bytes += buf->size + control_size;
            found = true;
        }

    if (!found) {
        stats[num_stats++] = (MemStats){ .buffers = 1, 
                                         .name    = buf->name, 
                                         .bytes   = buf->size + control_size };
        ASSERT (num_stats < MAX_MEMORY_STATS, "# memory stats exceeded %u, consider increasing MAX_MEMORY_STATS", MAX_MEMORY_STATS);
    }

    num_buffers++;

    return true; // continue
}

void buf_show_memory (bool memory_full, unsigned max_threads, unsigned used_threads)
{
    // if memory is full - only one thread needs to show this - others threads hang until we exit the process
    if (memory_full) mutex_lock (only_once_mutex);
    
    memset (stats, 0, sizeof(stats));
    num_stats = num_buffers = 0;
    buf_foreach_buffer (buf_add_mem_usage_to_stats, 0);
    
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

    fprintf (memory_full ? stderr : info_stream, "Total bytes: %s in %u buffers in %u buffer lists. global_max_threads=%u\n", 
             str_size (total_bytes).s, num_buffers, num_allocated_vbs, global_max_threads);
    if (IS_ZIP) 
        fprintf (memory_full ? stderr : info_stream, "vb_size = %u MB\n", (unsigned)(segconf.vb_size >> 20));
    
    if (max_threads)
        fprintf (memory_full ? stderr : info_stream, "Compute threads: max_permitted=%u actually_used=%u\n", max_threads, used_threads);

    for (unsigned i=0; i < num_stats; i++)
        fprintf (memory_full ? stderr : info_stream, "%-30s: %-8s (%4.1f%%) in %u buffers\n", stats[i].name, str_size (stats[i].bytes).s, 100.0 * (float)stats[i].bytes / (float)total_bytes, stats[i].buffers);

    fprintf (memory_full ? stderr : info_stream, "\n");
}

#ifndef _WIN32
// signal handler of SIGUSR1 
void buf_show_memory_handler (void) 
{
    buf_show_memory (false, 0, 0);
}
#endif

// thread safety: only the thread owning the VB of the buffer (main thread of evb) can add a buffer
// to the buf list OR it may be added by the main thread IF the compute thread of this VB is not running yet
void buf_add_to_buffer_list_do (VBlockP vb, BufferP buf, FUNCLINE)
{
    if (buf->vb == vb) return; // already in buf_list - nothing to do

    ASSERT (vb != evb || threads_am_i_main_thread(), "called from %s:%u: A thread other than main thread is attempting to modify evb buffer list: buf=%s", 
            func, code_line, buf_display (buf));

    ASSERT (!buf->vb, "called from %s:%u: cannot add buffer %s to buf_list of vb->vblock_i=%u because it is already in buf_list of vb_i=%u.",
             func, code_line, buf_desc(buf).s, vb->vblock_i, buf->vb->vblock_i);    

    ASSERT (((uint64_t)buf & 7)==0, "called from %s:%u: Expecting buffer %s of vb->vblock_i=%u to be word-aligned, but its address is %p",
             func, code_line, buf_desc(buf).s, vb->vblock_i, buf); // BL_SET_REMOVED expects this

#define INITIAL_MAX_MEM_NUM_BUFFERS 10000 /* for files that have ht,gt,phase,variant,and line - the factor would be about 5.5 so there will be 1 realloc per vb, but most files don't */
    BufferP bl = &vb->buffer_list;

    buf_alloc (vb, bl, 1, INITIAL_MAX_MEM_NUM_BUFFERS, BufListEnt, 2, "buffer_list");
    BNXT (BufListEnt, *bl) = (BufListEnt){ .buf = buf, .func = func, .code_line = code_line, .name = buf->name };
    bl->param = false; // buffer list is not sorted anymore

    buf->vb = vb; // successfully added to buf list

    if (flag.debug_memory==1 && vb->buffer_list.len > DISPLAY_ALLOCS_AFTER)
        iprintf ("buf_init_promiscuous &buf=%p (%s): %s: size=%"PRIu64" buffer=%p func=%s:%u vb->id=%d buf_i=%u\n", 
                 buf, func, buf_desc(buf).s, (uint64_t)buf->size, buf, func, code_line, vb->id, vb->buffer_list.len32-1);    
}

static void buf_init (BufferP buf, char *memory, uint64_t size, uint64_t old_size, 
                      FUNCLINE, rom name)
{
    // set some parameters before allocation so they can go into the error message in case of failure
    buf->func        = func;
    buf->code_line   = code_line;

    if (name) buf->name  = name;
    ASSERT (buf->name, "buffer has no name. func=%s:%u", buf->func, buf->code_line);

    if (!memory) {
        buf_show_memory (true, 0, 0);

        ABORT ("%s: Out of memory. %sDetails: %s:%u failed to allocate %s bytes. Buffer: %s", 
               global_cmd, 
               (command==ZIP ? "Try running with a lower vblock size using --vblock. " : ""), 
               func, code_line, str_int_commas (size + control_size).s, buf_desc(buf).s);
    }

    buf->data        = memory + sizeof (uint64_t);
    buf->size        = size;
    buf->overlayable = false;

    *(uint64_t *)memory = UNDERFLOW_TRAP;                    // underflow protection
    *(uint64_t *)(buf->data + size) = OVERFLOW_TRAP;         // overflow prortection (underflow protection was copied with realloc)
    *(uint16_t *)(buf->data + size + sizeof (uint64_t)) = 1; // counter of buffers that use of this memory (0 or 1 main buffer + any number of overlays)

    // only when we're done initializing - we update memory - that buf_test_overflow running concurrently doesn't test
    // half-initialized buffers
    __atomic_store_n (&buf->memory, memory, __ATOMIC_RELAXED); 
}

// allocates or enlarges buffer
// if it needs to enlarge a buffer fully overlaid by an overlay buffer - it abandons its memory (leaving it to
// the overlaid buffer) and allocates new memory
uint64_t buf_alloc_do (VBlockP vb, BufferP buf, uint64_t requested_size,
                       float grow_at_least_factor, // IF we need to allocate or reallocate physical memory, we get this much more than requested
                       FUNCLINE, rom name)      
{
    START_TIMER;

    if (!requested_size) return 0; // nothing to do

    ASSERT ((int64_t)requested_size > 0, "called from %s:%u: negative requested_size=%"PRId64" for name=%s", func, code_line, requested_size, name);

    ASSERT (requested_size <= MAX_BUFFER_SIZE, "called from %s:%u: Requested %s bytes which is beyond the Buffer maximum of %s",
            func, code_line, str_int_commas (requested_size).s, str_int_commas (MAX_BUFFER_SIZE).s);

#define REQUEST_TOO_BIG_THREADSHOLD (3 GB)
    if (requested_size > REQUEST_TOO_BIG_THREADSHOLD && !buf->can_be_big) // use WARN instead of ASSERTW to have a place for breakpoint
        WARN ("Warning: buf_alloc called from %s:%u requested %s. This is suspiciously high and might indicate a bug - please report to " EMAIL_SUPPORT ". vb->vblock_i=%u buf=%s line_i=%d",
              func, code_line, str_size (requested_size).s, vb->vblock_i, buf_desc (buf).s, vb->line_i);

    // sanity checks
    ASSERT (buf->type == BUF_REGULAR || buf->type == BUF_UNALLOCATED, "called from %s:%u: cannot buf_alloc a buffer of type %s. details: %s", 
            func, code_line, buf_display_type (buf), buf_desc (buf).s);

    ASSERT (vb, "called from %s:%u: null vb", func, code_line);

    // if this happens: either 1. the wrong VB was given now, or when initially allocating this buffer OR
    // 2. VB was REALLOCed in vb_get_vb, but for some reason this buf->vb was not updated because it was not on the buffer list
    ASSERT (!buf->vb || vb == buf->vb, "called from %s:%u: buffer=%p has wrong VB: vb=%p (id=%u vblock_i=%u) but buf->vb=%p", 
            func, code_line, buf, vb, vb->id, vb->vblock_i, buf->vb);

    ASSERT (vb != evb || buf->promiscuous || threads_am_i_main_thread(), "called from %s:%u: A non-main thread is attempting to allocate an evb buffer \"%s\" with promiscuous=false", 
            func, code_line, name ? name : buf->name);

    // case 1: we have enough memory already
    if (requested_size <= buf->size) {
        if (!buf->data) buf_init (buf, buf->memory, buf->size, buf->size, func, code_line, name);
        goto finish;
    }

    // add an epsilon to avoid floating point multiplication ending up slightly less than the integer
    grow_at_least_factor = MAX_(1.0001, grow_at_least_factor); 

    // grow us requested - rounding up to 64 bit boundary to avoid aliasing errors with the overflow indicator
    uint64_t new_size = MIN_(MAX_BUFFER_SIZE, ROUNDUP8((uint64_t)(requested_size * grow_at_least_factor)));

    ASSERT (new_size >= requested_size, "called from %s:%u: allocated too little memory for buffer %s: requested=%"PRIu64", allocated=%"PRIu64". vb->vblock_i=%u", 
            func, code_line, buf_desc (buf).s, requested_size, new_size, vb->vblock_i); // floating point paranoia

    // case 2: buffer was allocated already in the past - allocate new memory and copy over the data
    if (buf->memory) {

        ASSERT (buf->type != BUF_SHM, "called from %s:%u: SHM buffers cannot be realloced", func, code_line);

        uint64_t old_size = buf->size;

        // special handling if we have an overlaying buffer
        if (buf->overlayable) {
            mutex_lock (overlay_mutex);
            uint16_t *overlay_count = (uint16_t*)(buf->data + buf->size + sizeof(uint64_t));

            rom old_data     = buf->data;
            uint64_t old_len = buf->len;

            // if there is currently an overlay buffer on top of our buffer - abandon the memory
            // (leave it to the overlay buffer(s) that will eventually free() it), and allocate fresh memory
            if (*overlay_count > 1) {

                abandoned_mem_current += buf->size;
                abandoned_mem_high_watermark = MAX_(abandoned_mem_high_watermark, abandoned_mem_current);

                (*overlay_count)--; // overlaying buffers are now on their own - no regular buffer
                buf->memory = buf->data = NULL;
                buf->size = buf->len = 0;
                buf_alloc_do (vb, buf, new_size, 1, func, code_line, name); // recursive call - simple alloc
          
                // copy old data
                memcpy (buf->data, old_data, old_size);
                buf->len = old_len;
            }
            else {
                buf_verify_integrity (buf, func, code_line, "buf_alloc_do(overlayable)");

                // buffer is overlayable - but no current overlayers - regular realloc - however,
                // still within mutex to prevent another thread from overlaying while we're at it
                char *old_memory = buf->memory;
                __atomic_store_n (&buf->memory, BUFFER_BEING_MODIFIED, __ATOMIC_RELAXED);
                char *new_memory = (char *)buf_low_level_realloc (old_memory, new_size + control_size, name, func, code_line);
                buf_init (buf, new_memory, new_size, old_size, func, code_line, name);
            }
            buf->overlayable = true; // renew this, as it was reset by buf_init
            mutex_unlock (overlay_mutex);
        }

        else { // non-overlayable buffer - regular realloc without mutex
            buf_verify_integrity (buf, func, code_line, "buf_alloc_do");

            char *old_memory = buf->memory;
            __atomic_store_n (&buf->memory, BUFFER_BEING_MODIFIED, __ATOMIC_RELAXED);
            char *new_memory = (char *)buf_low_level_realloc (old_memory, new_size + control_size, name, func, code_line);
            buf_init (buf, new_memory, new_size, old_size, func, code_line, name);
        }
    }

    // case 3: we need to allocate memory - buffer is not yet allocated, so no need to copy data
    else {
        __atomic_store_n (&buf->memory, BUFFER_BEING_MODIFIED, __ATOMIC_RELAXED);
        char *memory = (char *)buf_low_level_malloc (new_size + control_size, false, func, code_line);
        ASSERT (memory != BUFFER_BEING_MODIFIED, "called from %s:%u: malloc didn't assign, very weird! buffer %s new_size=%"PRIu64,
                func, code_line, buf_desc(buf).s, new_size);

        buf->type = BUF_REGULAR;

        buf_init (buf, memory, new_size, 0, func, code_line, name);
        buf_add_to_buffer_list_do (vb, buf, func, code_line);
    }

    if (flag.debug_memory && !buf_locate (buf, false)) 
        // not in any VB, file, reference or gencomp structure
        iprintf ("buf_alloc_do: allocated independent buf %p: %s\n", buf, buf_desc(buf).s);

finish:
    if (vb != evb) COPY_TIMER (buf_alloc); // this is not thread-safe for evb as evb buffers might be allocated by any thread (?? is this still the case?)
    return buf->size;
}

BitsP buf_alloc_bits_do (VBlockP vb, BufferP buf, uint64_t nbits, FUNCLINE, rom name)
{
    uint64_t nwords = roundup_bits2words64 (nbits);
    
    if (!buf->data || buf->size < nwords * sizeof (uint64_t))
        buf_alloc_do (vb, buf, nwords * sizeof (uint64_t), 1, func, code_line, name);

    buf->nbits  = nbits;   
    buf->nwords = nwords; 

    bits_clear_excess_bits_in_top_word ((BitsP)buf);

    return (BitsP)buf;
}

void buf_alloc_bits_buffer_do (VBlockP vb, BufferP buf, uint64_t nbits, FUNCLINE, rom name)
{
    uint64_t nwords = roundup_bits2words64 (nbits);
    buf_alloc_do (vb, buf, nwords * sizeof (uint64_t), 1, func, code_line, name);
}

BitsP buf_overlay_bitarr_do (VBlockP vb,
                             BufferP overlaid_buf, BufferP regular_buf,  
                             uint64_t start_byte_in_regular_buf,
                             uint64_t nbits,
                             FUNCLINE, rom name)
{
    uint64_t nwords = roundup_bits2words64 (nbits);

    buf_overlay_do (evb, overlaid_buf, regular_buf, start_byte_in_regular_buf, func, code_line, name);

    overlaid_buf->nbits  = nbits;
    overlaid_buf->nwords = nwords;
    return (BitsP)overlaid_buf;
}

// an overlay buffer is a buffer using some of the memory of another buffer - it doesn't have its own memory
void buf_overlay_do (VBlockP vb, 
                     BufferP overlaid_buf, // dst 
                     BufferP regular_buf, 
                     uint64_t start_in_regular, // 0 means full overlay, and copy len 
                     FUNCLINE, rom name)
{   
     // if this buffer was used by a previous VB as a regular buffer - we need to "destroy" it first
    if (overlaid_buf->type == BUF_REGULAR && overlaid_buf->data == NULL && overlaid_buf->memory) {
        buf_low_level_free (overlaid_buf->memory, func, code_line);
        overlaid_buf->type = BUF_UNALLOCATED;
    }

    ASSERT (overlaid_buf->type == BUF_UNALLOCATED, "%s: Call from %s:%u: cannot buf_overlay to a buffer %s already in use", VB_NAME, func, code_line, buf_desc (overlaid_buf).s);

    // overlaying a SHM buffer, just creates another SHM buffer
    if (regular_buf->type == BUF_SHM) {
        ASSERT (start_in_regular < regular_buf->size, "%s: Call from %s:%u: expecting start_in_regular=%"PRIu64" < regular_buf->size=%"PRIu64, 
                VB_NAME, func, code_line, start_in_regular, (uint64_t)regular_buf->size);

        buf_attach_to_shm_do (vb, overlaid_buf, 
                              regular_buf->memory, regular_buf->size, start_in_regular, 
                              func, code_line, name);

        if (!start_in_regular) overlaid_buf->len = regular_buf->len;
        return;
    }

    ASSERT (regular_buf->type == BUF_REGULAR,
            "%s: Call from %s:%u: regular_buf %s in buf_overlay must be a regular or shm buffer", VB_NAME, func, code_line, buf_desc (regular_buf).s);
    //xxx ASSERT (regular_buf->overlayable, "%s: Call from %s:%u: buf_overlay: buffer %s is not overlayble", VB_NAME, func, code_line, buf_desc (regular_buf).s);

    regular_buf->overlayable  = true;

    overlaid_buf->type        = BUF_OVERLAY;
    overlaid_buf->memory      = 0;
    overlaid_buf->overlayable = false;
    overlaid_buf->name        = name;
    overlaid_buf->len         = start_in_regular ? 0ULL : regular_buf->len;
    overlaid_buf->func        = func;
    overlaid_buf->code_line   = code_line;

    // full or partial buffer overlay - if size=0, copy len too and update overlay counter
    mutex_lock (overlay_mutex);

    ASSERT (start_in_regular < regular_buf->size, 
            "called from %s:%u: not enough room in regular buffer for overlaid buf: start_in_regular=%"PRIu64" but regular_buf.size=%"PRIu64,
            func, code_line, start_in_regular, (uint64_t)regular_buf->size);

    // note: data+size MUST be at the control region, as we have the overlay counter there
    overlaid_buf->size = regular_buf->size - start_in_regular;
    overlaid_buf->data = regular_buf->data + start_in_regular;

    uint16_t *overlay_count = (uint16_t*)(regular_buf->data + regular_buf->size + sizeof(uint64_t));
    (*overlay_count)++; // counter of users of this memory
    
    mutex_unlock (overlay_mutex);

    buf_add_to_buffer_list_do (vb, overlaid_buf, func, code_line);
}

void buf_attach_to_shm_do (VBlockP vb, BufferP buf, void *memory, uint64_t size, uint64_t start, FUNCLINE, rom name)
{
    // if this buffer was used by a previous VB as a regular buffer - we need to "destroy" it first
    if (buf->type == BUF_REGULAR && buf->data == NULL && buf->memory) 
        buf_low_level_free (buf->memory, func, code_line);

    *buf = (Buffer){
        .type      = BUF_SHM,
        .name      = name,
        .func      = func,
        .code_line = code_line,
        .vb        = vb,
        .memory    = memory,
        .data      = memory + start, // note: no control area in shm buffers
        .size      = size
    };
}

static void buf_abandon_overlay (BufferP buf)
{
    mutex_lock (overlay_mutex);
    uint16_t *overlay_count = (uint16_t*)(buf->data + buf->size + sizeof(uint64_t)); // number of buffers (overlay and regular) sharing buf->memory

    if (*overlay_count > 1) { // current overlays exist - abandon memory - leave it to the overlaid buffer(s) which will free() this memory when they're done with it
        (*overlay_count)--;
    
        abandoned_mem_current += buf->size;
        abandoned_mem_high_watermark = MAX_(abandoned_mem_high_watermark, abandoned_mem_current);

        buf_reset (buf);
    }
    else  // case: no overlay exists
        buf->overlayable = false;
    
    mutex_unlock (overlay_mutex);            
}

// free buffer. A future buf_alloc of this buffer will reuse the memory if possible.
void buf_free_do (BufferP buf, FUNCLINE) 
{
    switch (buf->type) {

        case BUF_REGULAR: 

            buf_verify_integrity (buf, func, code_line, "buf_free_do");

            if (buf->overlayable)   
                buf_abandon_overlay (buf);

            // this causes crashes in rare cases, see bug 308 (solved?)

            // In Windows and Mac, we observe that free() operations are expensive and significantly slow down execution - so we
            // just recycle the same memory
            if (flag.is_linux && !buf->overlayable &&     // note: possibly lost overlayability in buf_abandon_overlay
                flag.show_memory != SHOW_MEM_PEAK) { // note: we don't free if --show-memory==PEAK, otherwise we won't report VB memory, only evb
                buf_low_level_free (buf->memory, func, code_line);
                buf->memory = NULL;
                buf->size   = 0;
                buf->type   = BUF_UNALLOCATED; // note: buf->vb is still set, so buf won't be added to buffer_list again in the next alloc
            }

            buf->data        = NULL; 
            buf->overlayable = false;
            buf->can_be_big  = false;

            // fall through (name (and in Windows/Mac also memory and size) are not changed, and buffer is still in buffer list)

        case BUF_UNALLOCATED: // reset len and param that may be used even without allocating the buffer
            buf->len         = 0;
            buf->param       = 0;
            break;

        case BUF_OVERLAY: {
            ASSERTNOTNULL (buf->data); // cannot be NULL if BUF_OVERLAY

            mutex_lock (overlay_mutex);
            uint16_t *overlay_count = (uint16_t*)(buf->data + buf->size + sizeof(uint64_t));
            (*overlay_count)--;
            mutex_unlock (overlay_mutex);            
            
            // we are the last user - we can free the memory now.
            // do this outside of the mutex - free is a system call and can take some time.
            // this is safe because if we ever observe *overlay_count==0, it means that no buffer has this memory,
            // therefore there is no possibility it would be subsequently overlayed between the test and the free().
            if (! (*overlay_count)) {
                buf_low_level_free (buf->data - sizeof(uint64_t), func, code_line); // the original buf->memory
                abandoned_mem_current -= buf->size;
            }
    
            VBlockP save_vb = buf->vb;
            buf_reset (buf);
            buf->vb   = save_vb;     // it is still in the buffer_list
            break;
        }

        case BUF_SHM:
            buf_reset (buf);
            break;

        default:
            ABORT0 ("Error: invalid buf->type");
    }
} 

static ASCENDING_SORTER (buf_list_sorter, BufListEnt, buf)
static BINARY_SEARCHER (buf_list_find_buffer, BufListEnt, BufferP, buf, true)

// after removing marking buffers as removed, actually remove them from the list
void buf_compact_buf_list (VBlockP vb)
{
    if (!vb) return;
    
    ARRAY32 (BufListEnt, bl, vb->buffer_list);

    uint32_t new_i = 0;
    for (uint32_t old_i=0; old_i < bl_len; old_i++)
        if (!BL_IS_REMOVED (bl[old_i].buf)) {
            if (new_i != old_i) bl[new_i] = bl[old_i];
            new_i++;
        }
            
    vb->buffer_list.len32 = new_i;
}

// sorts buffer list ahead of multiple calls to buf_find_in_buffer_list to speed it up
void buf_sort_buffer_list (VBlockP vb)
{
    if (cleanup_on_exit ||              // buf_destroy is disabled, so no need to sort either
       vb->buffer_list.param) return;   // already sorted

    qsort (B1ST(BufListEnt, vb->buffer_list), vb->buffer_list.len32, sizeof(BufListEnt), buf_list_sorter);
    vb->buffer_list.param = true; // now it's sorted
}

// remove from buffer_list of this vb
// thread safety: only the thread owning the VB of the buffer (main thread of evb) can remove a buffer
// from the buf list, OR the main thread may remove for all VBs, IF no compute thread is running
static BufListEnt *buf_find_in_buffer_list (BufferP buf)
{
    if (!buf->vb) return NULL; // this buffer is not on the buf_list - nothing to do

    START_TIMER;
    VBlockP vb = buf->vb;

    ASSERTISALLOCED (vb->buffer_list);

    if (!vb_is_valid (vb)) { // this is bug 576
        if (flag.debug_or_test) {
            WARN ("Unexpected, but unharmful: cannot remove buf=%p from buffer list, because buf->vb=%p refers to a VB that no longer exists", buf, vb);
            WARN ("%s allocated in %s:%u", buf->name, buf->func, buf->code_line); // separate line in case it segfaults if buf is corrupt
        }
        return NULL;
    }

    ASSERT (vb != evb || threads_am_i_main_thread(), "A thread other than main thread is attempting to modify evb buffer list: buf=%s", buf_display (buf));
    
    BufListEnt *ent = NULL;
    if (!vb->buffer_list.param) { // linear search if not sorted
        for_buf_back (BufListEnt, my_ent, vb->buffer_list) // back - we often destroy recently added
            if (my_ent->buf == buf) {
                ent = my_ent;
                break;
            }
    }
    else // binary search if sorted
        ent = binary_search (buf_list_find_buffer, BufListEnt, vb->buffer_list, buf);

    // note: it is possible that the buffer is not found in the list if it is never allocated or destroyed more than once. that's fine.
    
    COPY_TIMER (buf_find_in_buffer_list);
    return ent;
}

void buf_destroy_do_do (BufListEnt *ent, FUNCLINE)
{
    if (!ent) return;
    BufferP buf = ent->buf;
    
    if (flag.debug_memory==1) 
        iprintf ("Destroy %s: buf_addr=%p vb->id=%d buf_i=%u\n", buf_desc (buf).s, buf, buf->vb->id, BNUM (buf->vb->buffer_list, ent));

    BL_SET_REMOVED (ent->buf);
    buf->vb = NULL;

    // make sure that all overlayers have freed (applicable to BUF_REGULAR)
    uint16_t overlay_count = 1;
    if (buf->overlayable) {
        mutex_lock (overlay_mutex);
        overlay_count = (*(uint16_t*)(buf->data + buf->size + sizeof(uint64_t)));
        mutex_unlock (overlay_mutex);            
    }
    ASSERT (overlay_count==1, "called from %s:%u: cannot destroy buffer %s because it is currently overlaid", func, code_line, buf->name);

    switch (buf->type) {
        case BUF_REGULAR     : 
            buf_verify_integrity (buf, func, code_line, "buf_destroy_do");
            buf_low_level_free (buf->memory, func, code_line); 
            break;

        case BUF_OVERLAY     : buf_free (*buf); break;
        case BUF_SHM         : buf_free (*buf); break;
        case BUF_UNALLOCATED : break;
        default              : ABORT ("called from %s:%u: Error in buf_destroy_do: invalid buffer type %u", func, code_line, buf->type);
    }

    *buf = EMPTY_BUFFER; // reset to factory defaults
}

void buf_destroy_do (BufferP buf, FUNCLINE)
{
    if (!buf || cleanup_on_exit) return; // nothing to do (we don't destroy on exit, as the exiting thread may not be able to remove from buf_list)

    BufListEnt *ent = buf_find_in_buffer_list (buf); 

    if (ent) buf_destroy_do_do (ent, func, code_line);
}

void buf_destroy_vb_bufs (VBlockP vb)
{
    if (cleanup_on_exit) return;
    
    uint32_t sizeof_vb = DT_FUNC(vb, sizeof_vb)(vb->data_type_alloced);
    
    BufListEnt *buf_list_ent = NULL;
    for_buf (BufListEnt, ent, vb->buffer_list) {
        if (BL_IS_REMOVED(ent->buf)) continue;

        // all buffers that blong to the VB but are not in the VBlock structure, should be destroyed before calling this function
        ASSERT ((char*)ent->buf >= (char*)vb && (char*)ent->buf < ((char*)vb + sizeof_vb), "Found a buffer belonging to vb=%u, but not in VBlock: %s",
                vb->vblock_i, buf_desc (ent->buf).s);
    
        if (ent->buf != &vb->buffer_list) 
            buf_destroy_do_do (ent, __FUNCLINE);

        else 
            buf_list_ent = ent;
    }

    buf_destroy_do_do (buf_list_ent, __FUNCLINE);
}

void buf_destroy_file_bufs (FileP file)
{
    if (cleanup_on_exit) return;
    
    for_buf (BufListEnt, ent, evb->buffer_list) 
        // all buffers that blong to the VB but are not in the VBlock structure, should be destroyed before calling this function
        if (!BL_IS_REMOVED(ent->buf) &&
            (char*)ent->buf >= (char*)file && (char*)ent->buf < ((char*)file + sizeof (File))) 

            buf_destroy_do_do (ent, __FUNCLINE);

    buf_compact_buf_list (evb);
}

// main thread only: destroy all evb buffers with a specific name
void buf_destroy_by_name (rom name, bool compact_after_destroying)
{
    if (cleanup_on_exit) return;
    
    for_buf (BufListEnt, ent, evb->buffer_list) 
        if (!BL_IS_REMOVED(ent->buf) && ent->buf->name && !strcmp (ent->buf->name, name))
            buf_destroy_do_do (ent, __FUNCLINE);

    if (compact_after_destroying) buf_compact_buf_list (evb);
}

// similar to buf_move, but also moves buffer between buf_lists. can be run by the main thread only.
// IMPORTANT: only works when called from main thread, when BOTH src and dst VB are in full control of main thread, so that there
// no chance another thread is concurrently modifying the buf_list of the src or dst VBs 
void buf_grab_do (VBlockP dst_vb, BufferP dst_buf, rom dst_name, BufferP src_buf, FUNCLINE)
{
    ASSERT (src_buf, "called from %s:%u: buf is NULL", func, code_line);
    if (src_buf->type == BUF_UNALLOCATED) return; // nothing to grab

    ASSERT (src_buf->type == BUF_REGULAR && !src_buf->overlayable, "called from %s:%u: this function can only be called for a non-overlayable REGULAR buf", func, code_line);
    ASSERT (dst_buf->type == BUF_UNALLOCATED, "called from %s:%u: expecting dst_buf to be UNALLOCATED", func, code_line);

    dst_buf->type  = BUF_REGULAR;
    dst_buf->len   = src_buf->len;
    dst_buf->param = src_buf->param;
    buf_init (dst_buf, src_buf->memory, src_buf->size, 0, func, code_line, dst_name);

    buf_add_to_buffer_list_do (dst_vb, dst_buf, func, code_line);

    // reset src buffer, but keep vb - as it remains on the buffer list
    VBlockP save_vb = src_buf->vb;
    memset (src_buf, 0, sizeof (Buffer));
    src_buf->vb = save_vb;
}

// copy data - possibly within the same buffer
void buf_copy_do (VBlockP dst_vb, BufferP dst, ConstBufferP src, 
                  uint64_t bytes_per_entry, // how many bytes are counted by a unit of .len
                  uint64_t src_start_entry, uint64_t max_entries,  // if 0 copies the entire buffer 
                  FUNCLINE,
                  rom dst_name) // dst buffer settings, or take from src if 0
{
    ASSERTNOTNULL (src);
    ASSERTNOTNULL (dst);

    ASSERT (src->data, "called from %s:%u: src->data is NULL", func, code_line);
    
    ASSERT (!max_entries || src_start_entry < src->len, 
            "buf_copy of %s called from %s:%u: src_start_entry=%"PRIu64" is larger than src->len=%"PRIu64, buf_desc(src).s, func, code_line, src_start_entry, src->len);

    uint64_t num_entries = max_entries ? MIN_(max_entries, src->len - src_start_entry) : src->len - src_start_entry;
    if (!bytes_per_entry) bytes_per_entry=1;
    
    if (num_entries) {
        buf_alloc (dst_vb, dst, 0, num_entries * bytes_per_entry, char, 1, dst_name ? dst_name : src->name); 

        if (dst != src || src_start_entry >= num_entries)
            memcpy (dst->data, &src->data[src_start_entry * bytes_per_entry], num_entries * bytes_per_entry);
        else 
            memmove (dst->data, &src->data[src_start_entry * bytes_per_entry], num_entries * bytes_per_entry); // need memmove for overlapping areas
    }

    dst->len = num_entries;  
}   

// moves all the data from one buffer to another, leaving the source buffer unallocated
// both the src and dst buffer remain (or are added) to the buf_lists of their respective VBs
void buf_move (VBlockP dst_vb, BufferP dst, VBlockP src_vb, BufferP src)
{
    if (dst->type == BUF_REGULAR && !dst->data) buf_destroy (*dst);
    
    ASSERT (dst->type == BUF_UNALLOCATED, "attempt to move into an already-allocated buffer. src: %s dst: %s", buf_desc (src).s, buf_desc (dst).s);

    ASSERT (src_vb==dst_vb || dst->vb==dst_vb, "to move a buffer between VBs, the dst buffer needs to be added"
                                               " to the dst_vb buffer_list in advance. If dst_vb=evb the dst buffer must be added to"
                                               " the buffer_list by the main thread only. src: %s dst: %s src_vb->vb_i=%d dst_vb->vb_i=%d",
            buf_desc (src).s, buf_desc (dst).s, (src_vb ? src_vb->vblock_i : -999), (dst_vb ? dst_vb->vblock_i : -999));

    if (!dst->vb && (src->type == BUF_REGULAR || src->type == BUF_OVERLAY)) 
        buf_add_to_buffer_list_do (dst_vb, dst, __FUNCLINE); // this can only happen if src_vb==dst_vb 

    memcpy (dst, src, sizeof(Buffer));    
    dst->vb = dst_vb;

    buf_reset (src); // zero buffer except vb
}

// removes a section from the buffer
void buf_remove_do (BufferP buf, unsigned sizeof_item, uint64_t remove_start, uint64_t remove_len)
{
    if (!remove_len) return;

    ASSERT (remove_start + remove_len <= buf->len, "Out of range: remove_start=%"PRIu64" + remove_len=%"PRIu64" > buf->len=%"PRIu64,
            remove_start, remove_len, buf->len);

    if (remove_len != buf->len) { // skip in common case of deleting entire buffer 
        uint64_t remove_start_byte = remove_start * sizeof_item;
        uint64_t remove_bytes      = remove_len   * sizeof_item;
        memmove (buf->data + remove_start_byte, 
                 buf->data + remove_start_byte + remove_bytes, 
                 (buf->len * sizeof_item) - (remove_start_byte + remove_bytes));
    }

    buf->len -= remove_len;
}

void buf_insert_do (VBlockP vb, BufferP buf, unsigned width, uint64_t insert_at, const void *new_data, uint64_t new_data_len, rom name, FUNCLINE) 
{ 
    if (!new_data_len) return;

    buf_alloc_do (vb ? vb : buf->vb, buf, (buf->len + new_data_len + 1) * width/*room for \0 or separator*/, CTX_GROWTH, func, code_line, name); 

    if (insert_at != buf->len) {
        ASSERT (insert_at < buf->len, "called from %s:%u: expecting insert_at=%"PRIu64" <= buf->len=%"PRIu64" in buf=%s", 
                func, code_line, insert_at, buf->len, buf_desc(buf).s);

        memmove (&buf->data[(insert_at + new_data_len) * width], &buf->data[insert_at * width], (buf->len - insert_at) * width);
    }

    memcpy (&buf->data[insert_at * width], new_data, new_data_len * width);   
    buf->len += new_data_len; 

    buf_verify_integrity (buf, func, code_line, "buf_insert_do");
} 

void buf_append_string (VBlockP vb, BufferP buf, rom str) 
{ 
    uint64_t len = strlen (str); 
    ASSERT (len < 10000000, "len=%"PRIu64" too long, looks like a bug", len);

    buf_add_more (vb, buf, str, len, buf->name ? buf->name : "string_buf"); // allocates one char extra
    *BAFTc (*buf) = '\0'; // string terminator without increasing buf->len
}

void buf_print (BufferP buf, bool add_newline)
{
    for (uint64_t i=0; i < buf->len; i++) 
        fputc (buf->data[i], info_stream);  // safer than printf %.*s ?

    iprint0 (add_newline ? "\n" : "");
}

void buf_low_level_free (void *p, FUNCLINE)
{
    if (!p) return; // nothing to do

    START_TIMER;
    bool p_is_evb = (p == evb);

    if (flag.debug_memory==1) 
        iprintf ("Memory freed by free(): %p %s:%u\n", p, func, code_line);

    if (p == BUFFER_BEING_MODIFIED) {
        fprintf (stderr, "Warning in buf_low_level_free: corrupt pointer = 0x777 while attempting free()\n");
        return; // this can happen if there are memory overflow issues
    }

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

    if (flag.debug_memory && size >= flag.debug_memory) 
        iprintf ("realloc(): old=%p new=%p name=%s size=%"PRIu64" %s:%u\n", p, new, name, (uint64_t)size, func, code_line);

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

// convert a Buffer from a z_file section whose len is in char to a bitarray
Bits *buf_zfile_buf_to_bitarray (BufferP buf, uint64_t nbits)
{
    ASSERT (roundup_bits2bytes (nbits) <= buf->len, "nbits=%"PRId64" indicating a length of at least %"PRId64", but buf->len=%"PRId64,
            nbits, roundup_bits2bytes (nbits), buf->len);

    Bits *bits = (BitsP)buf;
    bits->nbits  = nbits;
    bits->nwords = roundup_bits2words64 (bits->nbits);

    ASSERT (roundup_bits2bytes64 (nbits) <= buf->size, "buffer to small: buf->size=%"PRId64" but bitarray has %"PRId64" words and hence requires %"PRId64" bytes",
            (uint64_t)buf->size, bits->nwords, bits->nwords * sizeof(uint64_t));

    LTEN_bits (bits);

    bits_clear_excess_bits_in_top_word (bits);

    return bits;
}

void buf_add_bit (BufferP buf, int64_t new_bit) 
{
    Bits *bar = (BitsP)buf;

    ASSERT (bar->nbits < buf->size * 8, "no room in Buffer %s to extend the bitmap", buf->name);
    bar->nbits++;     
    if (bar->nbits % 64 == 1) { // starting a new word                
        bar->nwords++;
        bar->words[bar->nwords-1] = new_bit; // LSb is as requested, other 63 bits are 0
    } 
    else
        bits_assign (bar, bar->nbits-1, new_bit);  
}
 
uint64_t buf_extend_bits (BufferP buf, int64_t num_new_bits) 
{
    Bits *bar = (BitsP)buf;

    ASSERT (bar->nbits + num_new_bits <= buf->size * 8, "Error in %s:%u: no room in Buffer %s to extend the bitmap: nbits=%"PRIu64", num_new_bits=%"PRId64", buf->size=%"PRIu64, 
            __FUNCLINE, buf->name, bar->nbits, num_new_bits, (uint64_t)buf->size);
    
    uint64_t next_bit = bar->nbits;

    bar->nbits += num_new_bits;     
    bar->nwords = roundup_bits2words64 (bar->nbits);
    bits_clear_excess_bits_in_top_word (bar);

    return next_bit;
}

void buf_erase_do (BufferP buf, FUNCLINE)
{
    // in Windows, buf_free keeps the buffer as BUF_REGULAR, so we no need to test if its UNALLOCATED
    ASSERT (flag.is_windows || buf->type == BUF_UNALLOCATED, "called from %s:%u: Cannot erase buffer because it is allocated: %s",
            func, code_line, buf_desc(buf).s);
    
    memset (buf, 0, sizeof(Buffer)); 
}


// writes a buffer to a file, return true if successful
// note: this is designed to run in any there, so it cannot create any buffers in evb
bool buf_dump_to_file (rom filename, ConstBufferP buf, unsigned buf_word_width, bool including_control_region, 
                       bool no_dirs, bool verbose, bool do_gzip)
{
    RETURNW (buf->type == BUF_REGULAR || buf->type == BUF_OVERLAY, false, 
             "FYI: failed to dump buffer.type=%s name=%s while putting %s", 
             buf_display_type (buf), buf->name ? buf->name : "(null)", filename);

    int fn_len = strlen(filename);
    char update_filename[fn_len + 10];
    strcpy (update_filename, filename);

    if (no_dirs) {
        for (unsigned i=0; i < fn_len; i++)
            if (filename[i] == '/' || (flag.is_windows && (filename[i] == '\\' || (filename[i] == ':' && i!=1))))
                update_filename[i] = '-';
        filename = update_filename;
    }

    bool success;
    if (including_control_region) {
        ASSERT (*(uint64_t *)(buf->memory)           == UNDERFLOW_TRAP, "dumping to %s: buffer has underflowed", filename);
        ASSERT (*(uint64_t *)(buf->data + buf->size) == OVERFLOW_TRAP,  "dumping to %s: buffer has overflowed",  filename);

        success = file_put_data (update_filename, buf->memory, buf->size + control_size, 0);
    }
    else
        success = file_put_data (update_filename, buf->data, buf->len * buf_word_width, 0);

    if (success && do_gzip) {
        char command[fn_len + 50];
        sprintf (command, "gzip -f \"%s\"", update_filename);
        int ret = system (command);
        ASSERTW (!ret, "FYI: \"%s\" returned %d. No harm.", command, ret); 

        if (!ret) {
            // special case: rename .bam.gz -> .bam
            if (fn_len >= 4 && !memcmp (&update_filename[fn_len-4], ".bam", 4)) {
                char gz_filename[fn_len + 10];
                sprintf (gz_filename, "%s.gz", update_filename);
                file_rename (gz_filename, update_filename, false);
            }
            else 
                strcpy (&update_filename[fn_len], ".gz");
        }
    }        

    if (success && verbose) iprintf ("\nDumped file %s\n", update_filename);

    return success;
}
 
void interlace_d8_buf       (BufferP buf, LocalType *lt) { for_buf (int8_t,  num, *buf) *num =        (INTERLACE(int8_t,  *num)); }
void BGEN_interlace_d16_buf (BufferP buf, LocalType *lt) { for_buf (int16_t, num, *buf) *num = BGEN16 (INTERLACE(int16_t, *num)); }
void BGEN_interlace_d32_buf (BufferP buf, LocalType *lt) { for_buf (int32_t, num, *buf) *num = BGEN32 (INTERLACE(int32_t, *num)); }
void BGEN_interlace_d64_buf (BufferP buf, LocalType *lt) { for_buf (int64_t, num, *buf) *num = BGEN64 (INTERLACE(int64_t, *num)); }
void LTEN_interlace_d16_buf (BufferP buf, LocalType *lt) { for_buf (int16_t, num, *buf) *num = LTEN16 (INTERLACE(int16_t, *num)); }
void LTEN_interlace_d32_buf (BufferP buf, LocalType *lt) { for_buf (int32_t, num, *buf) *num = LTEN32 (INTERLACE(int32_t, *num)); }
void LTEN_interlace_d64_buf (BufferP buf, LocalType *lt) { for_buf (int64_t, num, *buf) *num = LTEN64 (INTERLACE(int64_t, *num)); }

void BGEN_u8_buf  (BufferP buf, LocalType *lt) {}
void BGEN_u16_buf (BufferP buf, LocalType *lt) { if ( flag.is_lten) for_buf (uint16_t, num, *buf) *num = BGEN16 (*num); }
void BGEN_u32_buf (BufferP buf, LocalType *lt) { if ( flag.is_lten) for_buf (uint32_t, num, *buf) *num = BGEN32 (*num); }
void BGEN_u64_buf (BufferP buf, LocalType *lt) { if ( flag.is_lten) for_buf (uint64_t, num, *buf) *num = BGEN64 (*num); }
void LTEN_u16_buf (BufferP buf, LocalType *lt) { if (!flag.is_lten) for_buf (uint16_t, num, *buf) *num = LTEN16 (*num); }
void LTEN_u32_buf (BufferP buf, LocalType *lt) { if (!flag.is_lten) for_buf (uint32_t, num, *buf) *num = LTEN32 (*num); }
void LTEN_u64_buf (BufferP buf, LocalType *lt) { if (!flag.is_lten) for_buf (uint64_t, num, *buf) *num = LTEN64 (*num); }

// number of columns is trasmitted in the count, except if this is a matrix of VCF samples, in which case param=0 and we take 
// the number of columns to be the number of samples in the VCF header
static inline uint32_t BGEN_transpose_num_cols (ConstBufferP buf)
{
    uint32_t cols = buf->count; // cols and rows in terms of the target non-transposed matrix (0 if VCF)

    if (!cols) cols = vcf_header_get_num_samples(); 
    ASSERT0 (cols, "vcf_header_get_num_samples=0");
    
    return cols;
}

void BGEN_transpose_u8_buf (BufferP buf, LocalType *lt)
{
    if (!buf->len) return;

    uint32_t cols = BGEN_transpose_num_cols (buf);
    uint32_t rows = buf->len / cols;

    buf_alloc (buf->vb, &buf->vb->scratch, 0, buf->len, uint8_t, 1, "scratch");
    ARRAY (uint8_t, target, buf->vb->scratch);
    ARRAY (uint8_t, transposed, *buf);

    for (uint32_t c=0; c < cols; c++) 
        for (uint32_t r=0; r < rows; r++) 
            target[r * cols + c] = transposed[c * rows + r];

    buf->vb->scratch.len = buf->len;
    buf_copy (buf->vb, buf, &buf->vb->scratch, uint8_t, 0, 0, "contexts->local"); // copy and not move, so we can keep local's memory for next vb

    buf_free (buf->vb->scratch);

    if (lt) *lt = LT_UINT8; // no longer transposed
}

void BGEN_transpose_u16_buf (BufferP buf, LocalType *lt)
{
    if (!buf->len) return;

    uint32_t cols = BGEN_transpose_num_cols (buf);
    uint32_t rows = buf->len / cols;

    buf_alloc (buf->vb, &buf->vb->scratch, 0, buf->len, uint16_t, 1, "scratch");
    ARRAY (uint16_t, target, buf->vb->scratch);
    ARRAY (uint16_t, transposed, *buf);

    for (uint32_t c=0; c < cols; c++) 
        for (uint32_t r=0; r < rows; r++) 
            target[r * cols + c] = BGEN16 (transposed[c * rows + r]);

    buf->vb->scratch.len = buf->len;
    buf_copy (buf->vb, buf, &buf->vb->scratch, uint16_t, 0, 0, "contexts->local"); // copy and not move, so we can keep local's memory for next vb

    buf_free (buf->vb->scratch);

    *lt = LT_UINT16; // no longer transposed
}

void BGEN_transpose_u32_buf (BufferP buf, LocalType *lt)
{
    if (!buf->len) return;

    uint32_t cols = BGEN_transpose_num_cols (buf);
    uint32_t rows = buf->len / cols;

    buf_alloc (buf->vb, &buf->vb->scratch, 0, buf->len, uint32_t, 1, "scratch");
    ARRAY (uint32_t, target, buf->vb->scratch);
    ARRAY (uint32_t, transposed, *buf);

    for (uint32_t c=0; c < cols; c++) 
        for (uint32_t r=0; r < rows; r++) 
            target[r * cols + c] = BGEN32 (transposed[c * rows + r]);

    buf->vb->scratch.len = buf->len;
    buf_copy (buf->vb, buf, &buf->vb->scratch, uint32_t, 0, 0, "contexts->local"); // copy and not move, so we can keep local's memory for next vb

    buf_free (buf->vb->scratch);

    *lt = LT_UINT32; // no longer transposed
}

void BGEN_transpose_u64_buf (BufferP buf, LocalType *lt)
{
    if (!buf->len) return;

    uint32_t cols = BGEN_transpose_num_cols (buf);
    uint32_t rows = buf->len / cols;

    buf_alloc (buf->vb, &buf->vb->scratch, 0, buf->len, uint64_t, 1, "scratch");
    ARRAY (uint64_t, target, buf->vb->scratch);
    ARRAY (uint64_t, transposed, *buf);

    for (uint32_t c=0; c < cols; c++) 
        for (uint32_t r=0; r < rows; r++) 
            target[r * cols + c] = BGEN64 (transposed[c * rows + r]);

    buf->vb->scratch.len = buf->len;
    buf_copy (buf->vb, buf, &buf->vb->scratch, uint64_t, 0, 0, "contexts->local"); // copy and not move, so we can keep local's memory for next vb

    buf_free (buf->vb->scratch);

    *lt = LT_UINT64; // no longer transposed
}

void BGEN_deinterlace_d8_buf (BufferP buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint8_t unum = *B8 (*buf, i);
        *B(int8_t, *buf, i) = DEINTERLACE(int8_t,unum); 
    }
}

void BGEN_deinterlace_d16_buf (BufferP buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint16_t num_big_en = *B16 (*buf, i);
        uint16_t unum = BGEN16 (num_big_en);
        *B(int16_t, *buf, i) = DEINTERLACE(int16_t,unum); 
    }
}

void BGEN_deinterlace_d32_buf (BufferP buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint32_t num_big_en = *B32 ( *buf, i);
        uint32_t unum = BGEN32 (num_big_en);
        *B(int32_t, *buf, i) = DEINTERLACE(int32_t,unum); 
    }
}

void BGEN_deinterlace_d64_buf (BufferP buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint64_t num_big_en = *B64 (*buf, i);
        uint64_t unum = BGEN64 (num_big_en);
        *B(int64_t, *buf, i) = DEINTERLACE(int64_t,unum); 
    }
}
