// ------------------------------------------------------------------
//   buffer.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// memory management - when running the same code by the same thread for another variant block - we reuse
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
#include "bit_array.h"
#include "mutex.h"
#include "file.h"
#include "endianness.h"
#include "segconf.h"

#define DISPLAY_ALLOCS_AFTER 0 // display allocations, except the first X allocations. reallocs are always displayed

#define UNDERFLOW_TRAP 0x574F4C4652444E55ULL // "UNDRFLOW" - inserted at the begining of each memory block to detected underflows
#define OVERFLOW_TRAP  0x776F6C667265766FULL // "OVERFLOW" - inserted at the end of each memory block to detected overflows

#define BUFFER_BEING_MODIFIED ((char*)0x777)

static const unsigned control_size = 2*sizeof (uint64_t) + sizeof(uint16_t); // underflow, overflow and user counter

static Mutex overlay_mutex = {}; // used to thread-protect overlay counters (note: not initializing here - different in different OSes)
static Mutex only_once_mutex = {};
static uint64_t abandoned_mem_current = 0;
static uint64_t abandoned_mem_high_watermark = 0;

void buf_initialize()
{
    mutex_initialize (overlay_mutex);
    mutex_initialize (only_once_mutex);
}

static const char *buf_display_type (const Buffer *buf)
{
    static const char *names[] = BUFTYPE_NAMES;
    if (buf->type >= 0 && buf->type < BUF_NUM_TYPES) 
        return names[buf->type];
    else
        return "INVALID";
}

// get string with buffer's metadata for debug message. this function is NOT thread-safe
char *buf_display (const Buffer *buf)
{
    static char str[200]; // NOT thread-safe

    sprintf (str, "Buffer %s (%"PRId64"): size=%"PRIu64" len=%"PRIu64" data=%s memory=%s",
             buf->name, buf->param, buf->size, buf->len, str_pointer(buf->data).s, str_pointer(buf->memory).s);
    return str;    
}

const BufDescType buf_desc (const Buffer *buf)
{
    BufDescType desc; // use static memory instead of malloc since we could be in the midst of a memory issue when this is called
    sprintf (desc.s, "\"%s\" param=%"PRId64" len=%"PRIu64" size=%"PRId64" type=%s allocated in %s:%u by vb_i=%d", 
             buf->name ? buf->name : "(no name)", buf->param, buf->len, buf->size, buf_display_type (buf), buf->func ? buf->func : "(no func)", buf->code_line, 
             (buf->vb ? buf->vb->vblock_i : -999));
    return desc;
}

static inline void buf_reset (Buffer *buf)
{
    VBlockP save_vb = buf->vb;
    memset (buf, 0, sizeof (Buffer)); // make this buffer UNALLOCATED

    buf->vb = save_vb;
}

static inline bool buf_has_overflowed (const Buffer *buf, const char *msg)
{
    // memory==BUFFER_BEING_MODIFIED if an evb buffer is currently being allocated by another thread, and hence memory is not set yet.
    // for example, global_hash_prime is realloced by all threads (under mutex protection)
    // note: memory can become 0x777 at this point, even though we've tested it before at it was not yet
    char *memory = buf->memory;
    if (buf->vb && buf->vb->id==-1 && memory == BUFFER_BEING_MODIFIED) return false;

    ASSERT (memory != BUFFER_BEING_MODIFIED, "%s: buf->memory=BUFFER_BEING_MODIFIED. buffer %s size=%"PRIu64,
            msg, buf_desc(buf).s, buf->size);

    return *((uint64_t*)(memory + buf->size + sizeof(uint64_t))) != OVERFLOW_TRAP; // note on evb: if another thread reallocs the memory concurrently, this might seg-fault
}

static inline bool buf_has_underflowed (const Buffer *buf, const char *msg)
{
    if (buf->type == BUF_OVERLAY) return false; // overlayed buffers might be partial and start at a different address 

    // see comment in buf_has_overflowed
    char *memory = buf->memory;
    if (buf->vb && buf->vb->id==-1 && memory == BUFFER_BEING_MODIFIED) return false;

    ASSERT (memory != BUFFER_BEING_MODIFIED, "%s: buf->memory=BUFFER_BEING_MODIFIED. buffer %s size=%"PRIu64,
            msg, buf_desc(buf).s, buf->size);

    return *(uint64_t*)memory != UNDERFLOW_TRAP; // note on evb: if another thread reallocs the memory concurrently, this might seg-fault
}

// not thread-safe, used in emergency 
static void buf_find_underflow_culprit (const char *memory, const char *msg)
{
    VBlockPool *vb_pool = vb_get_pool();
    
    bool found=false;
    for (int vb_i=-1; vb_i < (int)vb_pool->num_vbs; vb_i++) {
        VBlock *vb = vb_get_from_pool (vb_i);

        if (!vb) continue;
        
        ARRAY (Buffer *, buf_list, vb->buffer_list);

        for (unsigned buf_i=0; buf_i < buf_list_len; buf_i++) {
            const Buffer *buf = buf_list[buf_i];

            if (buf) {
                char *after_buf = buf->memory + buf->size + control_size;
                if (after_buf <= memory && (after_buf + 100 > memory) && buf_has_overflowed (buf, msg)) {
                    char *of = &buf->memory[buf->size + sizeof(uint64_t)];
                    fprintf (stderr,
                            "Candidate culprit: vb_id=%d (vb_i=%d): buffer: %s %s memory: %s-%s name: %s vb_i=%u buf_i=%u Overflow fence=%c%c%c%c%c%c%c%c\n",
                            vb ? vb->id : -999, vb->vblock_i, buf_display_type(buf), 
                            str_pointer(buf).s, str_pointer(buf->memory).s, str_pointer(after_buf-1).s, 
                            buf_desc (buf).s, buf->vb->vblock_i, buf_i, 
                            of[0], of[1], of[2], of[3], of[4], of[5], of[6], of[7]);
                    found = true;
                }
            }
        }
    }
    if (!found) fprintf (stderr, "Cannot find a Buffer which has overflowed, and located just before the underflowed buffer\n\n");
}

// this function cannot contain ASSERT or ABORT as exit_on_error calls it
static bool buf_test_overflows_do (const VBlock *vb, bool primary, const char *msg);
static void buf_test_overflows_all_other_vb(const VBlock *caller_vb, const char *msg)
{
    // IMPORTANT: this function is not thread safe as it checks OTHER thread's memories which
    // may be dealloced as it checking them causing weird errors.
    // That's why the "return" is here. It can be removed when debugging specific cases.
    return; 

    VBlockPool *vb_pool = vb_get_pool();

    fprintf (stderr, "Testing all other VBs:\n");
    for (int vb_i=-1; vb_i < (int)vb_pool->num_vbs; vb_i++) {
        VBlock *vb = vb_get_from_pool (vb_i);
        if (!vb || vb == caller_vb) continue; // skip caller's VB
        buf_test_overflows_do (vb, false, msg);
    }
}

// this function cannot contain ASSERT or ABORT as exit_on_error calls it
static bool buf_test_overflows_do (const VBlock *vb, bool primary, const char *msg)
{
    if (!vb) return false;

    int corruption = 0;

    const Buffer *buf; // declare outside, so it is observable in the debugger in case of a crash
    for (unsigned buf_i=0; buf_i < vb->buffer_list.len; buf_i++) {
        static const char *nl[2] = {"", "\n\n"};

        // IMPORTANT NOTE regarding evb: testing evb might FAIL and should not be done in production! this if another thread 
        // can modify its buffers (under mutex protection) concurrently with this test, causing in consistent state between eg data and memory
        // we attempt to prevent many of the cases by not checking buffers that are BUFFER_BEING_MODIFIED at the onset, but they still may
        // become BUFFER_BEING_MODIFIED mid way through the test

        if (!(buf = *ENT (Buffer *, vb->buffer_list, buf_i)))
             continue; // buf was 'buf_destroy'd

#ifdef WIN32
        if (IsBadReadPtr (buf, sizeof (Buffer))) {
            fprintf (stderr, "%s%s: Memory corruption in vb_id=%d (vb_i=%d) buffer=%s (buf_i=%u): buffer structure inaccessible (invalid pointer)\n", 
                     nl[primary], msg, vb->id, vb->vblock_i, str_pointer(buf).s, buf_i);
            corruption = 100;
            break;
        }
#endif

        if (buf->memory && buf->memory != BUFFER_BEING_MODIFIED) {

            if (vb && buf->vb != vb) {
                fprintf (stderr, "%s%s: Memory corruption in vb_id=%d (vb_i=%d) buffer=%s (buf_i=%u): Corrupt Buffer structure OR invalid buffer pointer - buf->vb=%s != vb=%s\n", 
                         nl[primary], msg, vb->id, vb->vblock_i, str_pointer(buf).s, buf_i, str_pointer(buf->vb).s, str_pointer(vb).s);
                corruption = 0;
            }
            else if (buf->data && buf->vb->vblock_i != vb->vblock_i) { // buffers might still be here from the previous incarnation of this vb - its ok if they're not allocated yet
                        fprintf (stderr, "%s%s: Memory corruption in vb_id=%d: buf_vb_i=%d differs from thread_vb_i=%d: buffer: %s %s memory: %s-%s name: %s vb_i=%u buf_i=%u\n",
                        nl[primary], msg, vb ? vb->id : -999, buf->vb->vblock_i, vb->vblock_i, buf_display_type(buf), 
                        str_pointer(buf).s, str_pointer(buf->memory).s, str_pointer(buf->memory+buf->size+control_size-1).s,
                        buf_desc (buf).s, buf->vb->vblock_i, buf_i);
                corruption = 1;
            }
            else if (buf->type < 0 || buf->type > BUF_NUM_TYPES) {
                fprintf (stderr, "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d) buffer=%s (buf_i=%u): Corrupt Buffer structure OR invalid buffer pointer - invalid buf->type", 
                         nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, str_pointer (buf).s, buf_i);
                fprintf (stderr, " Buffer=%s\n", buf_desc(buf).s);  // separate fprintf in case it seg faults
                corruption = 2;
            }
            else if (!buf->name) {
                fprintf (stderr, "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d): buffer=%s (buf_i=%u): Corrupt Buffer structure - null name", 
                         nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, str_pointer (buf).s, buf_i);
                fprintf (stderr, " Buffer=%s\n", buf_desc(buf).s);  // separate fprintf in case it seg faults
                corruption = 3;
            }
            else if (buf->data && buf->type != BUF_OVERLAY && (buf->data != buf->memory + sizeof(uint64_t))) {
                fprintf (stderr, 
                         "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d): data!=memory+8: allocating_vb_i=%u buf_i=%u buffer=%s memory=%s name=%s : Corrupt Buffer structure - expecting data+8 == memory. buf->data=%s\n", 
                         nl[primary], msg, vb ? vb->id : -999, vb->vblock_i,  buf->vb->vblock_i, buf_i, str_pointer(buf).s, str_pointer(buf->memory).s, buf_desc(buf).s, str_pointer(buf->data).s);
                corruption = 4;
            }
            else if (buf_has_underflowed(buf, msg)) {
                fprintf (stderr, 
                        "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d): Underflow: buffer: %s %s memory: %s-%s name: %s vb_i=%u buf_i=%u. Fence=%c%c%c%c%c%c%c%c\n",
                        nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, buf_display_type(buf), str_pointer(buf).s, str_pointer(buf->memory).s, str_pointer(buf->memory+buf->size+control_size-1).s, 
                        buf_desc (buf).s, buf->vb->vblock_i, buf_i, 
                        buf->memory[0], buf->memory[1], buf->memory[2], buf->memory[3], buf->memory[4], buf->memory[5], buf->memory[6], buf->memory[7]);

                buf_find_underflow_culprit (buf->memory, msg);

                if (primary) buf_test_overflows_all_other_vb (vb, msg);
                primary = false;

                corruption = 5;
            }
            else if (buf_has_overflowed(buf, msg)) {
                char *of = &buf->memory[buf->size + sizeof(uint64_t)];
                fprintf (stderr,
                        "%s%s: Memory corruption in vb_id=%d (vb_i=%d): Overflow: buffer: %s %s memory: %s-%s name: %s vb_i=%u buf_i=%u Fence=%c%c%c%c%c%c%c%c\n",
                        nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, buf_display_type(buf), str_pointer(buf).s, str_pointer(buf->memory).s, str_pointer(buf->memory+buf->size+control_size-1).s, 
                        buf_desc (buf).s, buf->vb->vblock_i, buf_i, of[0], of[1], of[2], of[3], of[4], of[5], of[6], of[7]);
                
                if (primary) buf_test_overflows_all_other_vb (vb, msg);
                primary = false;

                corruption = 6;
            }
        }
    }
    
    ASSERT (!primary || !corruption, "Aborting due to memory corruption #%u", corruption); // primary will exit on corruption

    return corruption;
}

// this function cannot contain ASSERT or ABORT as exit_on_error calls it
void buf_test_overflows (void *vb, const char *msg)
{
    buf_test_overflows_do ((ConstVBlockP)vb, true, msg);
}

// this function cannot contain ASSERT or ABORT as exit_on_error calls it
void buf_test_overflows_all_vbs (const char *msg)
{
    buf_test_overflows_all_other_vb (NULL, msg);
}

static int buf_stats_sort_by_bytes(const void *a, const void *b)  
{ 
    return ((MemStats*)a)->bytes < ((MemStats*)b)->bytes ? 1 : -1;
}

void buf_update_buf_list_vb_addr_change (VBlockP new_vb, VBlockP old_vb)
{
    if (old_vb == new_vb) return; // no change in address
    
    ARRAY (Buffer *, bl, new_vb->buffer_list);
    
    for (uint32_t buf_i=0; buf_i < bl_len; buf_i++) {

        if (!bl[buf_i]) continue;

        // if this buffer is within the moved VB - update its address to the same offset relative to the new vb address
        if ((char*)bl[buf_i] >= (char*)old_vb && (char*)bl[buf_i] < (char *)(old_vb+1)) 
            bl[buf_i] = (Buffer *)((char*)new_vb + ((char*)bl[buf_i] - (char*)old_vb)); 

        // update VB address within the buffer
        bl[buf_i]->vb = new_vb; 
    }
}

static void buf_foreach_buffer (void (*callback)(const Buffer *, void *arg), void *arg)
{
    VBlockPool *vb_pool = vb_get_pool ();
    if (!vb_pool) return;

    // note: we don't cover EVB as it segfaults for an unknown reason (likely the Buffer structure itself resides with a structure
    // that is freed (eg File)). TODO: debug this.
    for (int vb_i=0; vb_i < (int)vb_pool->num_allocated_vbs; vb_i++) {

        VBlockP vb = vb_i >= 0 ? vb_get_from_pool (vb_i) : evb;
        if (!vb) continue;

        ARRAY (const Buffer *, bl, vb->buffer_list);
        
        for (uint32_t buf_i=0; buf_i < bl_len; buf_i++) {
            #ifdef WIN32
            if (IsBadReadPtr (&bl[buf_i], sizeof (Buffer))) 
                ABORT ("buffer structure inaccessible (invalid pointer) in buf_i=%u of vb_i=%d", buf_i, vb_i); 
            #endif

            if (bl[buf_i] && bl[buf_i]->memory) callback (bl[buf_i], arg); // exclude destroyed, not-yet-allocated, overlay buffers and buffers that were src in buf_move
        }
    }
}

static void buf_count_mem_usage (const Buffer *buf, void *mem_usage)
{
    *((uint64_t *)mem_usage) += buf->size + control_size;
}

uint64_t buf_get_memory_usage (void)
{
    uint64_t mem_usage=0;
    buf_foreach_buffer (buf_count_mem_usage, &mem_usage);

    return mem_usage;
}

#define MAX_MEMORY_STATS 100
static MemStats stats[MAX_MEMORY_STATS]; // must be pre-allocated, because buf_show_memory is called when malloc fails, so it cannot malloc
static unsigned num_stats=0, num_buffers=0;

static void buf_add_mem_usage_to_stats (const Buffer *buf, void *unused)
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

    // add non-Buffer reference memory
    if (flag.reference != REF_NONE) 
        stats[num_stats++] = ref_memory_consumption (gref); // only reports for DENOVO references

    // sort stats by bytes
    qsort (stats, num_stats, sizeof (MemStats), buf_stats_sort_by_bytes);

    uint64_t total_bytes=0;
    for (unsigned i=0; i< num_stats; i++) total_bytes += stats[i].bytes;

    fprintf (memory_full ? stderr : info_stream, "Total bytes: %s in %u buffers in %u buffer lists. global_max_threads=%u\n", 
             str_size (total_bytes).s, num_buffers, vb_get_pool() ? vb_get_pool()->num_allocated_vbs : 0, global_max_threads);
    if (command == ZIP) 
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
// to the buf list OR it may be added by the main thread IF the compute thread of this VB is not 
// running yet
void buf_add_to_buffer_list_do (VBlock *vb, Buffer *buf, const char *func)
{
    if (buf->vb == vb) return; // already in buf_list - nothing to do

    ASSERT (!buf->vb, "cannot add buffer %s to buf_list of vb_i=%u because it is already in buf_list of vb_i=%u.",
             buf_desc(buf).s, vb->vblock_i, buf->vb->vblock_i);    

#define INITIAL_MAX_MEM_NUM_BUFFERS 10000 /* for files that have ht,gt,phase,variant,and line - the factor would be about 5.5 so there will be 1 realloc per vb, but most files don't */
    Buffer *bl = &vb->buffer_list;

    buf_alloc (vb, bl, 1, INITIAL_MAX_MEM_NUM_BUFFERS, Buffer *, 2, "buffer_list");

    ((Buffer **)bl->data)[bl->len++] = buf;

    if (flag.debug_memory==1 && vb->buffer_list.len > DISPLAY_ALLOCS_AFTER)
        iprintf ("buf_add_to_buffer_list (%s): %s: size=%"PRIu64" buffer=%s vb->id=%d buf_i=%u\n", 
                 func, buf_desc(buf).s, buf->size, str_pointer(buf).s, vb->id, (uint32_t)vb->buffer_list.len-1);
    
    buf->vb = vb; // successfully added to buf list
}

static void buf_init (Buffer *buf, char *memory, uint64_t size, uint64_t old_size, 
                      const char *func, uint32_t code_line, const char *name)
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
               func, code_line, str_uint_commas (size + control_size).s, buf_desc(buf).s);
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
uint64_t buf_alloc_do (VBlock *vb,
                       Buffer *buf, 
                       uint64_t requested_size,
                       float grow_at_least_factor, // IF we need to allocate or reallocate physical memory, we get this much more than requested
                       const char *func, uint32_t code_line,
                       const char *name)      
{
    START_TIMER;

    if (!requested_size) return 0; // nothing to do

    ASSERT ((int64_t)requested_size > 0, "called from %s:%u: negative requested_size=%"PRId64" for name=%s", func, code_line, requested_size, name);

#define REQUEST_TOO_BIG_THREADSHOLD (3ULL << 30) // 3 GB
    if (requested_size > REQUEST_TOO_BIG_THREADSHOLD && !buf->can_be_big) // use WARN instead of ASSERTW to have a place for breakpoint
        WARN ("Warning: buf_alloc called from %s:%u requested %s. This is suspiciously high and might indicate a bug. buf=%s",
              func, code_line, str_size (requested_size).s, buf_desc (buf).s);

    // sanity checks
    ASSERT (buf->type == BUF_REGULAR || buf->type == BUF_UNALLOCATED, "called from %s:%u: cannot buf_alloc a buffer of type %s. details: %s", 
            func, code_line, buf_display_type (buf), buf_desc (buf).s);

    ASSERT (vb, "called from %s:%u: null vb", func, code_line);

    // if this happens: either 1. the wrong VB was given now, or when initially allocating this buffer OR
    // 2. VB was REALLOCed in vb_get_vb, but for some reason this buf->vb was not updated because it was not on the buffer list
    ASSERT (!buf->vb || vb == buf->vb, "called from %s:%u: buffer=%p has wrong VB: vb=%p (id=%u vblock_i=%u) but buf->vb=%p", 
            func, code_line, buf, vb, vb->id, vb->vblock_i, buf->vb);

    // case 1: we have enough memory already
    if (requested_size <= buf->size) {
        if (!buf->data) buf_init (buf, buf->memory, buf->size, buf->size, func, code_line, name);
        goto finish;
    }

    // add an epsilon to avoid floating point multiplication ending up slightly less than the integer
    grow_at_least_factor = MAX_(1.0001, grow_at_least_factor); 

    // grow us requested - rounding up to 64 bit boundary to avoid aliasing errors with the overflow indicator
    uint64_t new_size = (uint64_t)(requested_size * grow_at_least_factor + 7) & 0xfffffffffffffff8ULL; // aligned to 8 bytes

    ASSERT (new_size >= requested_size, "called from %s:%u: allocated too little memory for buffer %s: requested=%"PRIu64", allocated=%"PRIu64". vb_i=%u", 
            func, code_line, buf_desc (buf).s, requested_size, new_size, vb->vblock_i); // floating point paranoia

    // case 2: buffer was allocated already in the past - allocate new memory and copy over the data
    if (buf->memory) {

        uint64_t old_size = buf->size;

        // special handling if we have an overlaying buffer
        if (buf->overlayable) {
            mutex_lock (overlay_mutex);
            uint16_t *overlay_count = (uint16_t*)(buf->data + buf->size + sizeof(uint64_t));

            char *old_data = buf->data;
            uint32_t old_len = buf->len;

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
            char *old_memory = buf->memory;
            __atomic_store_n (&buf->memory, BUFFER_BEING_MODIFIED, __ATOMIC_RELAXED);
            char *new_memory = (char *)buf_low_level_realloc (old_memory, new_size + control_size, name, func, code_line);
            buf_init (buf, new_memory, new_size, old_size, func, code_line, name);
        }
    }

    // case 3: we need to allocate memory - buffer is not yet allocated, so no need to copy data
    else {
        __atomic_store_n (&buf->memory, BUFFER_BEING_MODIFIED, __ATOMIC_RELAXED);
        char *memory = (char *)malloc (new_size + control_size);
        ASSERT (memory != BUFFER_BEING_MODIFIED, "called from %s:%u: malloc didn't assign, very weird! buffer %s new_size=%"PRIu64,
                func, code_line, buf_desc(buf).s, new_size);

        buf->type = BUF_REGULAR;

        buf_init (buf, memory, new_size, 0, func, code_line, name);
        buf_add_to_buffer_list(vb, buf);
    }

finish:
    if (vb != evb) COPY_TIMER (buf_alloc); // this is not thread-safe for evb as evb buffers might be allocated by any thread (?? is this still the case?)
    return buf->size;
}

BitArray *buf_alloc_bitarr_do (VBlock *vb,
                               Buffer *buf, 
                               uint64_t nbits,
                               const char *func, uint32_t code_line,
                               const char *name)
{
    uint64_t nwords = roundup_bits2words64 (nbits);
    buf_alloc_do (vb, buf, nwords * sizeof (uint64_t), 1, func, code_line, name);
    BitArray *bitarr = buf_get_bitarray (buf);
    bitarr->nbits  = nbits;
    bitarr->nwords = nwords;

    return bitarr;
}

BitArray *buf_overlay_bitarr_do (VBlock *vb,
                                 Buffer *overlaid_buf, Buffer *regular_buf,  
                                 uint64_t start_byte_in_regular_buf,
                                 uint64_t nbits,
                                 const char *func, uint32_t code_line,
                                 const char *name)
{
    uint64_t nwords = roundup_bits2words64 (nbits);

    buf_overlay_do (evb, overlaid_buf, regular_buf, start_byte_in_regular_buf, func, code_line, name);

    BitArray *bitarr = buf_get_bitarray (overlaid_buf);
    bitarr->nbits  = nbits;
    bitarr->nwords = nwords;
    return bitarr;
}

// an overlay buffer is a buffer using some of the memory of another buffer - it doesn't have its own memory
void buf_overlay_do (VBlock *vb, 
                     Buffer *overlaid_buf, // dst 
                     Buffer *regular_buf, 
                     uint64_t start_in_regular, // 0 means full overlay, and copy len 
                     const char *func, uint32_t code_line, const char *name)
{
    // if this buffer was used by a previous VB as a regular buffer - we need to "destroy" it first
    if (overlaid_buf->type == BUF_REGULAR && overlaid_buf->data == NULL && overlaid_buf->memory) {
        buf_low_level_free (overlaid_buf->memory, func, code_line);
        overlaid_buf->type = BUF_UNALLOCATED;
    }
    
    ASSERT (overlaid_buf->type == BUF_UNALLOCATED, "Error in %s:%u: cannot buf_overlay to a buffer %s already in use", func, code_line, buf_desc (overlaid_buf).s);
    ASSERT (regular_buf->type == BUF_REGULAR || regular_buf->type == BUF_MMAP, "Error in %s:%u: regular_buf %s in buf_overlay must be a regular or mmap buffer", func, code_line, buf_desc (regular_buf).s);
    ASSERT (regular_buf->overlayable, "Error in %s:%u: buf_overlay: buffer %s is not overlayble", func, code_line, buf_desc (regular_buf).s);

    overlaid_buf->type        = BUF_OVERLAY;
    overlaid_buf->memory      = 0;
    overlaid_buf->overlayable = false;
    overlaid_buf->name        = name;
    overlaid_buf->len         = start_in_regular ? 0 : regular_buf->len;
    // note: we don't add overlaid_buf to buf_list, and because of that we also don't set its overlaid_buf->vb because it won't get updated in buf_update_buf_list_vb_addr_change

    // full or partial buffer overlay - if size=0, copy len too and update overlay counter
    mutex_lock (overlay_mutex);

    ASSERT (start_in_regular < regular_buf->size, 
            "called from %s:%u: not enough room in regular buffer for overlaid buf: start_in_regular=%"PRIu64" but regular_buf.size=%"PRIu64,
            func, code_line, start_in_regular, regular_buf->size);

    // note: data+size MUST be at the control region, as we have the overlay counter there
    overlaid_buf->size = regular_buf->size - start_in_regular;
    overlaid_buf->data = regular_buf->data + start_in_regular;

    uint16_t *overlay_count = (uint16_t*)(regular_buf->data + regular_buf->size + sizeof(uint64_t));
    (*overlay_count)++; // counter of users of this memory

    mutex_unlock (overlay_mutex);
}

// creates a file mapping: data is mapping from a read-only file, any modifications are private
// to the process and not written back to the file. buf->param is used for mmapping.
bool buf_mmap_do (VBlock *vb, Buffer *buf, const char *filename, 
                  bool read_only_buffer, // if false, make a copy-on-write memory mapping, creating private pages upon write
                  const char *func, uint32_t code_line, const char *name)
{
    START_TIMER;

    int fd = -1;

    if (!file_exists (filename)) return false; 

    // if this buffer was used by a previous VB as a regular buffer - we need to "destroy" it first
    if (buf->type == BUF_REGULAR && buf->data == NULL && buf->memory) 
        buf_low_level_free (buf->memory, func, code_line);

    uint64_t file_size = file_get_size (filename);

    *buf = (Buffer){
        .type      = BUF_MMAP,
        .name      = name,
        .func      = func,
        .code_line = code_line,
        .vb        = vb,
        .size      = file_size - control_size
    };

    fd = open (filename, O_RDONLY);    

#ifdef _WIN32
    HANDLE file = (HANDLE)_get_osfhandle (fd);
    buf->param = (int64_t)CreateFileMapping (file, NULL, PAGE_WRITECOPY, file_size >> 32, file_size & 0xffffffff, NULL);
    ASSERTGOTO (buf->param, "Error in buf_mmap of %s file_size=%"PRIu64": CreateFileMapping failed: %s", filename, file_size, str_win_error());
    
    // note that mmap'ed buffers include the Buffer 
    buf->memory = MapViewOfFile ((HANDLE)buf->param, read_only_buffer ? FILE_MAP_READ : FILE_MAP_COPY, 0, 0, file_size);
    ASSERTGOTO (buf->memory, "Error in buf_mmap of %s file_size=%"PRIu64": MapViewOfFile failed: %s", filename, file_size, str_win_error());

#else
    buf->memory = mmap (NULL, file_size, PROT_READ | (read_only_buffer ? 0 : PROT_WRITE), MAP_PRIVATE, fd, 0);
    ASSERT (buf->memory != MAP_FAILED, "Error in buf_mmap of %s file_size=%"PRIu64": mmap failed: %s", filename, file_size, strerror (errno));
#endif
    close (fd);
    fd=-1;

    // verify buffer integrity
    buf->data = buf->memory + sizeof (uint64_t);
    ASSERTGOTO (*(uint64_t *)(buf->memory) == UNDERFLOW_TRAP, "Error in buf_mmap of %s: mmap'ed buffer has corrupt underflow trap", filename);
    ASSERTGOTO (*(uint64_t *)(buf->data + buf->size) == OVERFLOW_TRAP, "Error in buf_mmap of %s: mmap'ed buffer has corrupt overflow trap - possibly file is trucated", filename);

    // reset overlay counter    
    if (!read_only_buffer)
        *(uint16_t *)(buf->data + buf->size + sizeof (uint64_t)) = 1;
    
    COPY_TIMER (buf_mmap_do);
    
    return true;

error:
    // close mapping and delete file - ignore errors
    if (fd >= 0) close (fd);
#ifdef _WIN32
    UnmapViewOfFile (buf->data);
    CloseHandle ((HANDLE)buf->param);
#else
    munmap (buf->data, buf->size);
#endif

    memset (buf, 0, sizeof(Buffer));
    return false;
}

// free buffer - without freeing memory. A future buf_alloc of this buffer will reuse the memory if possible.
void buf_free_do (Buffer *buf, const char *func, uint32_t code_line) 
{
    uint16_t *overlay_count; // number of buffers (overlay and regular) sharing buf->memory

    switch (buf->type) {

        case BUF_REGULAR: 

            if (buf->overlayable) {

                mutex_lock (overlay_mutex);
                overlay_count = (uint16_t*)(buf->data + buf->size + sizeof(uint64_t));

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

// this causes crashes in rare cases, see bug 308
/*#ifdef __linux__
            // In Windows and Mac, we observe that free() operations are expensive and significantly slow down execution - so we
            // just recycle the same memory
            if (!buf->overlayable) {
                buf_low_level_free (buf->memory, func, code_line);
                buf->memory = NULL;
                buf->size   = 0;
            }
#endif*/
            buf->data        = NULL; 
            buf->overlayable = false;

            // fall through (name (and in Windows also memory and size) are not changed, and buffer is still in buffer list)

        case BUF_UNALLOCATED: // reset len and param that may be used even without allocating the buffer
            buf->len         = 0;
            buf->param       = 0;
            break;

        case BUF_OVERLAY:
            mutex_lock (overlay_mutex);
            overlay_count = (uint16_t*)(buf->data + buf->size + sizeof(uint64_t));
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
    
            buf_reset (buf);
            break;

        case BUF_MMAP:
#ifdef _WIN32
            ASSERT (UnmapViewOfFile (buf->memory), "called from %s:%u: UnmapViewOfFile failed: %s", buf->func, buf->code_line, str_win_error());
            ASSERT (CloseHandle ((HANDLE)buf->param), "called from %s:%u: CloseHandle failed: %s", buf->func, buf->code_line, str_win_error());            
#else
            ASSERT (!munmap (buf->memory, buf->size), "called from %s:%u: munmap failed: %s", buf->func, buf->code_line, strerror (errno));
#endif
            buf_reset (buf);
            break;

        default:
            ABORT0 ("Error: invalid buf->type");
    }
} 

// remove from buffer_list of this vb
// thread safety: only the thread owning the VB of the buffer (main thread of evb) can remove a buffer
// from the buf list, OR the main thread may remove for all VBs, IF no compute thread is running
static void buf_remove_from_buffer_list (Buffer *buf)
{
    if (!buf->vb) return; // this buffer is not on the buf_list - nothing to do

    ARRAY (Buffer *, buf_list, buf->vb->buffer_list);

     for (unsigned i=0; i < buf->vb->buffer_list.len; i++) 

        if (buf_list[i] == buf) {
            
            if (flag.debug_memory==1) 
                iprintf ("Destroy %s: buf_addr=%s buf->vb->id=%d buf_i=%u\n", buf_desc (buf).s, str_pointer(buf).s, buf->vb->id, i);

            buf_list[i] = NULL;
            buf->vb = NULL;

            break;
        }

    // note: it is possible that the buffer is not found in the list if it is never allocated or destroyed more than once. that's fine.
}

void buf_destroy_do (Buffer *buf, const char *func, uint32_t code_line)
{
    if (!buf) return; // nothing to do

    buf_remove_from_buffer_list (buf); 

    // make sure that all overlayers have freed (applicable to BUF_REGULAR and BUF_MMAP)
    uint16_t overlay_count = 1;
    if (buf->overlayable) {
        mutex_lock (overlay_mutex);
        overlay_count = (*(uint16_t*)(buf->data + buf->size + sizeof(uint64_t)));
        mutex_unlock (overlay_mutex);            
    }
    ASSERT (overlay_count==1, "called from %s:%u: cannot destroy buffer %s because it is currently overlaid", func, code_line, buf->name);

    switch (buf->type) {
        case BUF_REGULAR     : buf_low_level_free (buf->memory, func, code_line); break;
        case BUF_OVERLAY     : buf_free (buf);   /* stop overlaying */            break;
        case BUF_MMAP        : buf_free (buf);   /* stop mmap'ing   */            break;
        case BUF_UNALLOCATED :                                                    break;
        default              : ABORT ("called from %s:%u: Error in buf_destroy_do: invalid buffer type %u", func, code_line, buf->type);
    }

    memset (buf, 0, sizeof (Buffer)); // reset to factory defaults
}

// similar to buf_move, but for evb only
// IMPORTANT: only works when called from main thread, when BOTH src and dst VB are in full control of main thread, so that there
// no chance another thread is concurrently modifying the buf_list of the src or dst VBs 
void buf_grab_do (VBlock *dst_vb, Buffer *dst_buf, const char *dst_name, Buffer *src_buf, const char *func, uint32_t code_line)
{
    ASSERT (src_buf, "called from %s:%u: buf is NULL", func, code_line);
    if (src_buf->type == BUF_UNALLOCATED) return; // nothing to grab

    ASSERT (src_buf->type == BUF_REGULAR && !src_buf->overlayable, "called from %s:%u: this function can only be called for a non-overlayable REGULAR buf", func, code_line);
    ASSERT (dst_buf->type == BUF_UNALLOCATED, "called from %s:%u: expecting dst_buf to be UNALLOCATED", func, code_line);

    buf_remove_from_buffer_list (src_buf); 

    dst_buf->type = BUF_REGULAR;
    dst_buf->len = src_buf->len;
    buf_init (dst_buf, src_buf->memory, src_buf->size, 0, func, code_line, dst_name);

    buf_add_to_buffer_list (dst_vb, dst_buf);

    memset (src_buf, 0, sizeof (Buffer)); // reset to factory defaults
}

// copy data - possibly within the same buffer
void buf_copy_do (VBlock *dst_vb, Buffer *dst, const Buffer *src, 
                  uint64_t bytes_per_entry, // how many bytes are counted by a unit of .len
                  uint64_t src_start_entry, uint64_t max_entries,  // if 0 copies the entire buffer 
                  const char *func, unsigned code_line,
                  const char *dst_name) // dst buffer settings, or take from src if 0
{
    ASSERT (src->data, "called from %s:%u: src->data is NULL", func, code_line);
    
    ASSERT (!max_entries || src_start_entry < src->len, 
            "buf_copy of %s called from %s:%u: src_start_entry=%"PRIu64" is larger than src->len=%"PRIu64, buf_desc(src).s, func, code_line, src_start_entry, src->len);

    uint64_t num_entries = max_entries ? MIN_(max_entries, src->len - src_start_entry) : src->len - src_start_entry;
    if (!bytes_per_entry) bytes_per_entry=1;
    
    if (num_entries) {
        buf_alloc (dst_vb, dst, 0, num_entries * bytes_per_entry, char, 1, dst_name ? dst_name : src->name); // use realloc rather than malloc to allocate exact size

        if (dst != src || src_start_entry >= num_entries)
            memcpy (dst->data, &src->data[src_start_entry * bytes_per_entry], num_entries * bytes_per_entry);
        else 
            memmove (dst->data, &src->data[src_start_entry * bytes_per_entry], num_entries * bytes_per_entry); // need memmove for overlapping areas
    }

    dst->len = num_entries;  
}   

// moves all the data from one buffer to another, leaving the source buffer unallocated
// both the src and dst buffer remain (or are added) to the buf_lists of their respective VBs
void buf_move (VBlock *dst_vb, Buffer *dst, VBlock *src_vb, Buffer *src)
{
    if (dst->type == BUF_REGULAR && !dst->data) buf_destroy (dst);
    
    ASSERT (dst->type == BUF_UNALLOCATED, "attempt to move to an already-allocated src: %s dst: %s", buf_desc (src).s, buf_desc (dst).s);

    ASSERT (src_vb==dst_vb || dst->vb==dst_vb, "to move a buffer between VBs, the dst buffer needs to be added"
                                               " to the dst_vb buffer_list in advance. If dst_vb=evb the dst buffer must be added to"
                                               " the buffer_list by the main thread only. src: %s dst: %s src_vb->vb_i=%d dst_vb->vb_i=%d",
            buf_desc (src).s, buf_desc (dst).s, (src_vb ? src_vb->vblock_i : -999), (dst_vb ? dst_vb->vblock_i : -999));

    if (!dst->vb) buf_add_to_buffer_list (dst_vb, dst); // this can only happen if src_vb==dst_vb 

    memcpy (dst, src, sizeof(Buffer));    
    dst->vb = dst_vb;

    buf_reset (src); // zero buffer except vb
}

void buf_add_string (VBlockP vb, Buffer *buf, const char *str) 
{ 
    uint64_t len = strlen (str); 
    ASSERT (len < 10000000, "len=%"PRIu64" too long, looks like a bug", len);

    buf_add_more (vb, buf, str, len, buf->name ? buf->name : "string_buf"); // allocates one char extra
    buf->data[buf->len] = '\0'; // string terminator without increasing buf->len
}

void buf_add_int_as_text (VBlockP vb, Buffer *buf, int64_t value)
{
    char s[20];
    unsigned len = str_int (value, s);
    buf_add (buf, s, len);
}

void buf_print (Buffer *buf, bool add_newline)
{
    for (uint64_t i=0; i < buf->len; i++) 
        fputc (buf->data[i], info_stream);  // safer than printf %.*s ?

    iprint0 (add_newline ? "\n" : "");
}

void buf_low_level_free (void *p, const char *func, uint32_t code_line)
{
    if (!p) return; // nothing to do
    
    if (flag.debug_memory==1) 
        iprintf ("Memory freed by free(): %s %s:%u\n", str_pointer (p).s, func, code_line);

    if (p == BUFFER_BEING_MODIFIED) {
        fprintf (stderr, "Warning in buf_low_level_free: corrupt pointer = 0x777 while attempting free()\n");
        return; // this can happen if there are memory overflow issues
    }

    free (p);
}

void *buf_low_level_realloc (void *p, size_t size, const char *name, const char *func, uint32_t code_line)
{
    void *new = realloc (p, size);
    ASSERTW (new, "Out of memory in %s:%u: realloc failed (name=%s size=%"PRIu64" bytes). %s", func, code_line, name, (uint64_t)size, 
             command == ZIP ? "Try limiting the number of concurrent threads with --threads (affects speed) or reducing the amount of data processed by each thread with --vblock (affects compression ratio)" : "");

    if (flag.debug_memory && size >= flag.debug_memory) 
        iprintf ("realloc(): old=%p new=%p name=%s size=%"PRIu64" %s:%u\n", p, new, name, (uint64_t)size, func, code_line);

    return new;
}

void *buf_low_level_malloc (size_t size, bool zero, const char *func, uint32_t code_line)
{
    void *new = malloc (size);
    ASSERT (new, "Out of memory in %s:%u: malloc failed (size=%"PRIu64" bytes). %s", func, code_line, (uint64_t)size,
            command == ZIP ? "Try limiting the number of concurrent threads with --threads (affects speed) or reducing the amount of data processed by each thread with --vblock (affects compression ratio)" : "");

    if (flag.debug_memory && size >= flag.debug_memory) 
        iprintf ("malloc(): %s size=%"PRIu64" %s:%u\n", str_pointer (new).s, (uint64_t)size, func, code_line);

    if (zero) memset (new, 0, size);
    
    return new;
}

// convert a Buffer from a z_file section whose len is in char to a bitarray
BitArray *buf_zfile_buf_to_bitarray (Buffer *buf, uint64_t nbits)
{
    ASSERT (roundup_bits2bytes (nbits) <= buf->len, "nbits=%"PRId64" indicating a length of at least %"PRId64", but buf->len=%"PRId64,
            nbits, roundup_bits2bytes (nbits), buf->len);

    BitArray *bitarr = buf_get_bitarray (buf);
    bitarr->nbits  = nbits;
    bitarr->nwords = roundup_bits2words64 (bitarr->nbits);

    ASSERT (roundup_bits2bytes64 (nbits) <= buf->size, "buffer to small: buf->size=%"PRId64" but bitarray has %"PRId64" words and hence requires %"PRId64" bytes",
            buf->size, bitarr->nwords, bitarr->nwords * sizeof(uint64_t));

    LTEN_bit_array (bitarr);

    return bitarr;
}

void buf_add_bit (Buffer *buf, int64_t new_bit) 
{
    BitArray *bar = buf_get_bitarray (buf);

    ASSERT (bar->nbits < buf->size * 8, "Error in %s:%u: no room in Buffer %s to extend the bitmap", __FUNCTION__, __LINE__, buf->name);
    bar->nbits++;     
    if (bar->nbits % 64 == 1) { // starting a new word                
        bar->nwords++;
        bar->words[bar->nwords-1] = new_bit; // LSb is as requested, other 63 bits are 0
    } 
    else
        bit_array_assign (bar, bar->nbits-1, new_bit);  
}
 
bit_index_t buf_extend_bits (Buffer *buf, int64_t num_new_bits) 
{
    BitArray *bar = buf_get_bitarray (buf);

    ASSERT (bar->nbits + num_new_bits <= buf->size * 8, "Error in %s:%u: no room in Buffer %s to extend the bitmap: nbits=%"PRIu64", num_new_bits=%"PRId64", buf->size=%"PRIu64, 
            __FUNCTION__, __LINE__, buf->name, bar->nbits, num_new_bits, buf->size);
    
    bit_index_t next_bit = bar->nbits;

    bar->nbits += num_new_bits;     
    bar->nwords = roundup_bits2words64 (bar->nbits);
    bit_array_clear_excess_bits_in_top_word (bar);

    return next_bit;
}

// writes a buffer to a file, return true if successful
// note: this is run as separate thread from ref_create_cache_in_background and refhash_create_cache_in_background
// so it cannot allocate any buffers
bool buf_dump_to_file (const char *filename, const Buffer *buf, unsigned buf_word_width, bool including_control_region, 
                       bool no_dirs, bool verbose, bool do_gzip)
{
    RETURNW (buf->type == BUF_REGULAR, false, "FYI: failed to dump buffer.type=%s name=%s while putting %s", 
             buf_display_type (buf), buf->name ? buf->name : "(null)", filename);

    int fn_len = strlen(filename);
    char update_filename[fn_len + 10];
    strcpy (update_filename, filename);

    if (no_dirs) {
        for (unsigned i=0; i < fn_len; i++)
            if (filename[i] == '/' || (flag.is_windows && (filename[i] == '\\' || filename[i] == ':')))
                update_filename[i] = '-';
        filename = update_filename;
    }

    bool success;
    if (including_control_region) {
        ASSERT (*(uint64_t *)(buf->memory) == UNDERFLOW_TRAP, "dumping to %s: buffer has underflowed", filename);
        ASSERT (*(uint64_t *)(buf->data + buf->size) == OVERFLOW_TRAP, "dumping to %s: buffer has underflowed", filename);

        success = file_put_data (update_filename, buf->memory, buf->size + control_size, 0);
    }
    else
        success = file_put_data (update_filename, buf->data, buf->len * buf_word_width, 0);

    if (success && do_gzip) {
        char command[fn_len + 50];
        sprintf (command, "gzip %s", update_filename);
        system (command); // ignore failure

        // special case: rename .bam.gz -> .bam
        if (fn_len >= 4 && !memcmp (&update_filename[fn_len-4], ".bam", 4)) {
            char gz_filename[fn_len + 10];
            sprintf (gz_filename, "%s.gz", update_filename);
            rename (gz_filename, update_filename);
        }
        else 
            strcpy (&update_filename[fn_len], ".gz");
    }        

    if (success && verbose) iprintf ("\nDumped file %s\n", update_filename);

    return success;
}
 
void BGEN_u8_buf (Buffer *buf, LocalType *lt)
{
}

void BGEN_u16_buf (Buffer *buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint16_t num_big_en = *ENT (uint16_t, *buf, i);
        *ENT (uint16_t, *buf, i) = BGEN16 (num_big_en);            
    }
}

void BGEN_u32_buf (Buffer *buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint32_t num_big_en = *ENT (uint32_t, *buf, i);
        *ENT (uint32_t, *buf, i) = BGEN32 (num_big_en);            
    }
}

void BGEN_u64_buf (Buffer *buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint64_t num_big_en = *ENT (uint64_t, *buf, i);
        *ENT (uint64_t, *buf, i) = BGEN64 (num_big_en);            
    }
}

// number of columns is trasmitted in the param, except if this is a matrix of VCF samples, in which case param=0 and we take 
// the number of columns to be the number of samples in the VCF header
static inline uint32_t BGEN_transpose_num_cols (const Buffer *buf)
{
    uint32_t cols = buf->param; // cols and rows in terms of the target non-transposed matrix (0 if VCF)

    if (!cols) cols = vcf_header_get_num_samples(); 
    ASSERT0 (cols, "vcf_header_get_num_samples=0");
    
    return cols;
}

void BGEN_transpose_u8_buf (Buffer *buf, LocalType *lt)
{
    if (!buf->len) return;

    uint32_t cols = BGEN_transpose_num_cols (buf);
    uint32_t rows = buf->len / cols;

    buf_alloc (buf->vb, &buf->vb->compressed, 0, buf->len, uint8_t, 1, "compressed");
    ARRAY (uint8_t, target, buf->vb->compressed);
    ARRAY (uint8_t, transposed, *buf);

    for (uint32_t c=0; c < cols; c++) 
        for (uint32_t r=0; r < rows; r++) 
            target[r * cols + c] = transposed[c * rows + r];

    buf->vb->compressed.len = buf->len;
    buf_copy (buf->vb, buf, &buf->vb->compressed, uint8_t, 0, 0, "contexts->local"); // copy and not move, so we can keep local's memory for next vb

    buf_free (&buf->vb->compressed);

    *lt = LT_UINT8; // no longer transposed
}

void BGEN_transpose_u16_buf (Buffer *buf, LocalType *lt)
{
    if (!buf->len) return;

    uint32_t cols = BGEN_transpose_num_cols (buf);
    uint32_t rows = buf->len / cols;

    buf_alloc (buf->vb, &buf->vb->compressed, 0, buf->len, uint16_t, 1, "compressed");
    ARRAY (uint16_t, target, buf->vb->compressed);
    ARRAY (uint16_t, transposed, *buf);

    for (uint32_t c=0; c < cols; c++) 
        for (uint32_t r=0; r < rows; r++) 
            target[r * cols + c] = BGEN16 (transposed[c * rows + r]);

    buf->vb->compressed.len = buf->len;
    buf_copy (buf->vb, buf, &buf->vb->compressed, uint16_t, 0, 0, "contexts->local"); // copy and not move, so we can keep local's memory for next vb

    buf_free (&buf->vb->compressed);

    *lt = LT_UINT16; // no longer transposed
}

void BGEN_transpose_u32_buf (Buffer *buf, LocalType *lt)
{
    if (!buf->len) return;

    uint32_t cols = BGEN_transpose_num_cols (buf);
    uint32_t rows = buf->len / cols;

    buf_alloc (buf->vb, &buf->vb->compressed, 0, buf->len, uint32_t, 1, "compressed");
    ARRAY (uint32_t, target, buf->vb->compressed);
    ARRAY (uint32_t, transposed, *buf);

    for (uint32_t c=0; c < cols; c++) 
        for (uint32_t r=0; r < rows; r++) 
            target[r * cols + c] = BGEN32 (transposed[c * rows + r]);

    buf->vb->compressed.len = buf->len;
    buf_copy (buf->vb, buf, &buf->vb->compressed, uint32_t, 0, 0, "contexts->local"); // copy and not move, so we can keep local's memory for next vb

    buf_free (&buf->vb->compressed);

    *lt = LT_UINT32; // no longer transposed
}

void BGEN_transpose_u64_buf (Buffer *buf, LocalType *lt)
{
    if (!buf->len) return;

    uint32_t cols = BGEN_transpose_num_cols (buf);
    uint32_t rows = buf->len / cols;

    buf_alloc (buf->vb, &buf->vb->compressed, 0, buf->len, uint64_t, 1, "compressed");
    ARRAY (uint64_t, target, buf->vb->compressed);
    ARRAY (uint64_t, transposed, *buf);

    for (uint32_t c=0; c < cols; c++) 
        for (uint32_t r=0; r < rows; r++) 
            target[r * cols + c] = BGEN64 (transposed[c * rows + r]);

    buf->vb->compressed.len = buf->len;
    buf_copy (buf->vb, buf, &buf->vb->compressed, uint64_t, 0, 0, "contexts->local"); // copy and not move, so we can keep local's memory for next vb

    buf_free (&buf->vb->compressed);

    *lt = LT_UINT64; // no longer transposed
}

void BGEN_deinterlace_d8_buf (Buffer *buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint8_t unum = *ENT (uint8_t, *buf, i);
        *ENT (int8_t, *buf, i) = DEINTERLACE(int8_t,unum); 
    }
}

void BGEN_deinterlace_d16_buf (Buffer *buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint16_t num_big_en = *ENT (uint16_t, *buf, i);
        uint16_t unum = BGEN16 (num_big_en);
        *ENT (int16_t, *buf, i) = DEINTERLACE(int16_t,unum); 
    }
}

void BGEN_deinterlace_d32_buf (Buffer *buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint32_t num_big_en = *ENT (uint32_t, *buf, i);
        uint32_t unum = BGEN32 (num_big_en);
        *ENT (int32_t, *buf, i) = DEINTERLACE(int32_t,unum); 
    }
}

void BGEN_deinterlace_d64_buf (Buffer *buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint64_t num_big_en = *ENT (uint64_t, *buf, i);
        uint64_t unum = BGEN64 (num_big_en);
        *ENT (int64_t, *buf, i) = DEINTERLACE(int64_t,unum); 
    }
}
