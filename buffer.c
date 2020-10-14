// ------------------------------------------------------------------
//   buffer.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// memory management - when running the same code by the same thread for another variant block - we reuse
// the previous variant's block memory. this way we save repetitive malloc/free cycles which might
// be very time consuming.

#include "genozip.h"
#include "profiler.h"
#include "buffer.h"
#include "vblock.h"
#include "strings.h"
#include "reference.h"
#include "bit_array.h"
#include "mutex.h"

#define DISPLAY_ALLOCS_AFTER 0 // display allocations, except the first X allocations. reallocs are always displayed

#define UNDERFLOW_TRAP 0x574F4C4652444E55ULL // "UNDRFLOW" - inserted at the begining of each memory block to detected underflows
#define OVERFLOW_TRAP  0x776F6C667265766FULL // "OVERFLOW" - inserted at the end of each memory block to detected overflows

#define BUFFER_BEING_MODIFIED ((char*)0x777)

static const unsigned overhead_size = 2*sizeof (uint64_t) + sizeof(uint16_t); // underflow, overflow and user counter

MUTEX (overlay_mutex); // used to thread-protect overlay counters (note: not initializing here - different in different OSes)
static uint64_t abandoned_mem_current = 0;
static uint64_t abandoned_mem_high_watermark = 0;

void buf_initialize()
{
    mutex_initialize (overlay_mutex);
}

static const char *bt_str (const Buffer *buf)
{
    static const char *names[] = BUFTYPE_NAMES;
    if (buf->type >= BUF_UNALLOCATED && buf->type <= BUF_OVERLAY) 
        return names[buf->type];
    else
        return "INVALID";
}

// get string with buffer's metadata for debug message. this function is NOT thread-safe
char *buf_display (const Buffer *buf)
{
    static char str[200]; // NOT thread-safe

    char s1[POINTER_STR_LEN], s2[POINTER_STR_LEN];
    sprintf (str, "Buffer %s (%"PRId64"): size=%"PRIu64" len=%"PRIu64" data=%s memory=%s",
             buf->name, buf->param, buf->size, buf->len, str_pointer(buf->data, s1), str_pointer(buf->memory,s2));
    return str;    
}

// NOT THREAD SAFE - used in error messages terminating executation
const char *buf_desc (const Buffer *buf)
{
    static char desc[300]; // use static memory instead of malloc since we could be in the midst of a memory issue when this is called
    sprintf (desc, "%s:%"PRId64" len=%"PRIu64" allocated in %s:%u by vb_i=%d", 
             buf->name ? buf->name : "(no name)", buf->param, buf->len, buf->func, buf->code_line, (buf->vb ? buf->vb->vblock_i : -999));
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

    ASSERT (memory != BUFFER_BEING_MODIFIED, "%s: Error in buf_has_overflowed: buf->memory=BUFFER_BEING_MODIFIED. buffer %s size=%"PRIu64,
            msg, buf_desc(buf), buf->size);

    return *((uint64_t*)(memory + buf->size + sizeof(uint64_t))) != OVERFLOW_TRAP; // note on evb: if another thread reallocs the memory concurrently, this might seg-fault
}

static inline bool buf_has_underflowed (const Buffer *buf, const char *msg)
{
    // see comment in buf_has_overflowed
    char *memory = buf->memory;
    if (buf->vb && buf->vb->id==-1 && memory == BUFFER_BEING_MODIFIED) return false;

    ASSERT (memory != BUFFER_BEING_MODIFIED, "%s: Error in buf_has_underflowed: buf->memory=BUFFER_BEING_MODIFIED. buffer %s size=%"PRIu64,
            msg, buf_desc(buf), buf->size);

    return *(uint64_t*)memory != UNDERFLOW_TRAP; // note on evb: if another thread reallocs the memory concurrently, this might seg-fault
}

// not thread-safe, used in emergency 
static void buf_find_underflow_culprit (const char *memory, const char *msg)
{
    VBlockPool *vb_pool = vb_get_pool();
    char s1[POINTER_STR_LEN], s2[POINTER_STR_LEN], s3[POINTER_STR_LEN];
    
    bool found=false;
    for (int vb_i=-1; vb_i < (int)vb_pool->num_vbs; vb_i++) {
        VBlock *vb = (vb_i == -1) ? evb : vb_pool->vb[vb_i]; 

        if (!vb) continue;
        
        ARRAY (Buffer *, buf_list, vb->buffer_list);

        for (unsigned buf_i=0; buf_i < vb->buffer_list.len; buf_i++) {
            const Buffer *buf = buf_list[buf_i];

            if (buf) {
                char *after_buf = buf->memory + buf->size + overhead_size;
                if (after_buf <= memory && (after_buf + 100 > memory) && buf_has_overflowed (buf, msg)) {
                    char *of = &buf->memory[buf->size + sizeof(uint64_t)];
                    fprintf (stderr,
                            "Candidate culprit: vb_id=%d (vb_i=%d): buffer: %s %s memory: %s-%s name: %s vb_i=%u buf_i=%u Overflow fence=%c%c%c%c%c%c%c%c\n",
                            vb ? vb->id : -999, vb->vblock_i, bt_str(buf), str_pointer(buf,s1), str_pointer(buf->memory,s2), str_pointer(after_buf-1,s3), 
                            buf_desc (buf), buf->vb->vblock_i, buf_i, 
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
        VBlock *vb = (vb_i == -1) ? evb : vb_pool->vb[vb_i]; 
        if (vb == caller_vb) continue; // skip caller's VB
        buf_test_overflows_do (vb, false, msg);
    }
}

// this function cannot contain ASSERT or ABORT as exit_on_error calls it
static bool buf_test_overflows_do (const VBlock *vb, bool primary, const char *msg)
{
    if (!vb) return false;

    const Buffer *buf_list = &vb->buffer_list;

    int corruption = 0;
    const Buffer *buf; // declare outside, so it is observable in the debugger in case of a crash
    for (unsigned buf_i=0; buf_i < buf_list->len; buf_i++) {

// IMPORTANT NOTE regarding evb: testing evb might FAIL and should not be done in production! this is another thread thread
// can modify its buffers (under mutex protection) concurrently with this test, causing in consistent state between eg data and memory
// we attempt to prevent many of the cases by not checking buffers that are BUFFER_BEING_MODIFIED at the onset, but they still may
// become BUFFER_BEING_MODIFIED mid way through the test

        buf = ((Buffer **)buf_list->data)[buf_i];

        if (!buf) continue; // buf was 'buf_destroy'd

        static const char *nl[2] = {"", "\n\n"};
        char s1[POINTER_STR_LEN], s2[POINTER_STR_LEN], s3[POINTER_STR_LEN];
        if (buf->memory && buf->memory != BUFFER_BEING_MODIFIED) {

            if (buf->data && buf->vb->vblock_i != vb->vblock_i) { // buffers might still be here from the previous incarnation of this vb - its ok if they're not allocated yet
                        fprintf (stderr, "%s%s: Memory corruption in vb_id=%d: buf_vb_i=%d differs from thread_vb_i=%d: buffer: %s %s memory: %s-%s name: %s vb_i=%u buf_i=%u\n",
                        nl[primary], msg, vb ? vb->id : -999, buf->vb->vblock_i, vb->vblock_i, bt_str(buf), str_pointer(buf,s1), str_pointer(buf->memory,s2), str_pointer(buf->memory+buf->size+overhead_size-1,s3),
                        buf_desc (buf), buf->vb->vblock_i, buf_i);
                corruption = 1;
            }
            if (buf->type < BUF_UNALLOCATED || buf->type > BUF_OVERLAY) {
                fprintf (stderr, "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d) buffer=%s (buf_i=%u): Corrupt Buffer structure OR invalid buffer pointer - invalid buf->type\n", 
                         nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, str_pointer (buf, s1), buf_i);
                corruption = 2;
            }
            else if (!buf->name) {
                fprintf (stderr, "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d): buffer=%s (buf_i=%u): Corrupt Buffer structure - null name\n", 
                         nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, str_pointer (buf, s1), buf_i);
                corruption = 3;
            }
            else if (buf->data && (buf->data != buf->memory + sizeof(uint64_t))) {
                fprintf (stderr, 
                         "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d): data!=memory+8: allocating_vb_i=%u buf_i=%u buffer=%s memory=%s name=%s : Corrupt Buffer structure - expecting data+8 == memory. buf->data=%s\n", 
                         nl[primary], msg, vb ? vb->id : -999, vb->vblock_i,  buf->vb->vblock_i, buf_i, str_pointer(buf,s1), str_pointer(buf->memory,s2), buf_desc(buf), str_pointer(buf->data,s3));
                corruption = 4;
            }
            else if (buf_has_underflowed(buf, msg)) {
                fprintf (stderr, 
                        "%s%s: Memory corruption in vb_id=%d (thread vb_i=%d): Underflow: buffer: %s %s memory: %s-%s name: %s vb_i=%u buf_i=%u. Fence=%c%c%c%c%c%c%c%c\n",
                        nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, bt_str(buf), str_pointer(buf,s1), str_pointer(buf->memory,s2), str_pointer(buf->memory+buf->size+overhead_size-1,s3), 
                        buf_desc (buf), buf->vb->vblock_i, buf_i, 
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
                        nl[primary], msg, vb ? vb->id : -999, vb->vblock_i, bt_str(buf), str_pointer(buf,s1), str_pointer(buf->memory,s2), str_pointer(buf->memory+buf->size+overhead_size-1,s3), 
                        buf_desc (buf), buf->vb->vblock_i, buf_i, of[0], of[1], of[2], of[3], of[4], of[5], of[6], of[7]);
                
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

void buf_display_memory_usage (bool memory_full, unsigned max_threads, unsigned used_threads)
{
    #define MAX_MEMORY_STATS 100
    static MemStats stats[MAX_MEMORY_STATS]; // must be pre-allocated, because buf_display_memory_usage is called when malloc fails, so it cannot malloc
    unsigned num_stats = 0, num_buffers = 0;

    if (memory_full)
        fprintf (stderr, "\n\nError memory is full:\n");
    else
        fprintf (stderr, "\n-------------------------------------------------------------------------------------\n");

    VBlockPool *vb_pool = vb_get_pool ();

    for (int vb_i=-1; vb_i < (int)vb_pool->num_allocated_vbs; vb_i++) {

        Buffer *buf_list = (vb_i == -1) ? &evb->buffer_list
                                        : &vb_pool->vb[vb_i]->buffer_list; // a pointer to a buffer, which contains an array of pointers to buffers of a single vb/non-vb

        if (!buf_list->len) continue; // no buffers allocated yet for this VB

        for (unsigned buf_i=0; buf_i < buf_list->len; buf_i++) {
    
            ASSERT (buf_list->memory, "Error: memory of buffer_list of vb_i=%u is not allocated", vb_i); // this should never happen

            Buffer *buf = ((Buffer **)buf_list->data)[buf_i];
            
            if (!buf || !buf->memory) continue; // exclude destroyed, not-yet-allocated, overlay buffers and buffers that were src in buf_move

            bool found = false;
            for (unsigned st_i=0; st_i < num_stats && !found; st_i++) {
                MemStats *st = &stats[st_i];

                if (!strcmp (st->name, buf->name)) {
                    st->buffers++;
                    st->bytes += buf->size + overhead_size;
                    found = true;
                }
            }

            if (!found) {
                stats[num_stats].name    = buf->name;
                stats[num_stats].bytes   = buf->size + overhead_size;
                stats[num_stats].buffers = 1;
                num_stats++;
                ASSERT (num_stats < MAX_MEMORY_STATS, "# memory stats exceeded %u, consider increasing MAX_MEMORY_STATS", MAX_MEMORY_STATS);
            }

            num_buffers++;
        }
    }

    // add non-Buffer reference memory
    if (flag_reference != REF_NONE) 
        stats[num_stats++] = ref_memory_consumption();

    // sort stats by bytes
    qsort (stats, num_stats, sizeof (MemStats), buf_stats_sort_by_bytes);

    uint64_t total_bytes=0;
    for (unsigned i=0; i< num_stats; i++) total_bytes += stats[i].bytes;

    char str[30];
    str_size (total_bytes, str);
    fprintf (stderr, "Total bytes: %s in %u buffers in %u buffer lists:\n", str, num_buffers, vb_pool->num_allocated_vbs);
    fprintf (stderr, "Compute threads: max_permitted=%u actually_used=%u\n", max_threads, used_threads);

    for (unsigned i=0; i < num_stats; i++) {
        str_size (stats[i].bytes, str);
        fprintf (stderr, "%-30s: %-8s (%4.1f%%) in %u buffers\n", stats[i].name, str, 100.0 * (float)stats[i].bytes / (float)total_bytes, stats[i].buffers);
    }
}

// thread safety: only the thread owning the VB of the buffer (I/O thread of evb) can add a buffer
// to the buf list OR it may be added by the I/O thread IF the compute thread of this VB is not 
// running yet
void buf_add_to_buffer_list (VBlock *vb, Buffer *buf)
{
    if (buf->vb == vb) return; // already in buf_list - nothing to do

    ASSERT (!buf->vb, "Error in buf_add_to_buffer_list: cannot add buffer %s to buf_list of vb_i=%u because it is already in buf_list of vb_i=%u.",
            buf_desc(buf), vb->vblock_i, buf->vb->vblock_i);    

#define INITIAL_MAX_MEM_NUM_BUFFERS 10000 /* for files that have ht,gt,phase,variant,and line - the factor would be about 5.5 so there will be 1 realloc per vb, but most files don't */
    Buffer *bl = &vb->buffer_list;

    buf_alloc (vb, bl, MAX (INITIAL_MAX_MEM_NUM_BUFFERS, bl->len+1) * sizeof(Buffer *), 2, "buffer_list", vb->id);

    ((Buffer **)bl->data)[bl->len++] = buf;

    if (flag_debug_memory && vb->buffer_list.len > DISPLAY_ALLOCS_AFTER) {
        char s[POINTER_STR_LEN];
        fprintf (stderr, "buf_add_to_buffer_list: %s: size=%"PRIu64" buffer=%s vb->id=%d buf_i=%u\n", 
                 buf_desc(buf), buf->size, str_pointer(buf,s), vb->id, (uint32_t)vb->buffer_list.len-1);
    }
    
    buf->vb = vb; // successfully added to buf list
}

static void buf_init (Buffer *buf, char *memory, uint64_t size, uint64_t old_size, 
                      const char *func, uint32_t code_line, const char *name, int64_t param)
{
    // set some parameters before allocation so they can go into the error message in case of failure
    buf->func        = func;
    buf->code_line   = code_line;

    if (name) {
        buf->name  = name;
        buf->param = param;
    } 
    ASSERT (buf->name, "Error: buffer has no name. func=%s:%u", buf->func, buf->code_line);

    if (!memory) {
        if (flag_show_memory)
            buf_display_memory_usage (true, 0, 0);

        char s[30];
        ABORT ("%s: Out of memroy. %s%sDetails: %s:%u failed to allocate %s bytes. Buffer: %s", 
               global_cmd, 
               (command==ZIP ? "Try running with a lower vblock size using --vblock. " : ""), 
               (!flag_show_memory ? "To see memory details - run again with --show-memory. " : ""), 
               func, code_line, str_uint_commas (size + overhead_size, s), buf_desc(buf));
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
                       double grow_at_least_factor, // IF we need to allocate or reallocate physical memory, we get this much more than requested
                       const char *func, uint32_t code_line,
                       const char *name, int64_t param)      
{
    START_TIMER;

    if (!requested_size) return 0; // nothing to do

#define REQUEST_TOO_BIG_THREADSHOLD (4ULL*1024*1024*1024) // 4 GB
    char s[30];
    ASSERTW (requested_size < REQUEST_TOO_BIG_THREADSHOLD, "Warning: buf_alloc called from %s:%u requested %s. This is suspeciously high and might indicate a bug",
             func, code_line, str_size (requested_size, s));

    // sanity checks
    ASSERT (buf->type == BUF_REGULAR || buf->type == BUF_UNALLOCATED, "Error in buf_alloc_do called from %s:%u: cannot buf_alloc an overlayed buffer. details: %s", 
            func, code_line, buf_desc (buf));

    ASSERT0 (vb, "Error in buf_alloc_do: null vb");

    // case 1: we have enough memory already
    if (requested_size <= buf->size) {
        if (!buf->data) buf_init (buf, buf->memory, buf->size, buf->size, func, code_line, name, param);
        goto finish;
    }

    // add an epsilon to avoid floating point multiplication ending up slightly less that the integer
    grow_at_least_factor = MAX (1.00000001, grow_at_least_factor); 

    // grow us requested - rounding up to 64 bit boundary to avoid aliasing errors with the overflow indicator
    uint64_t new_size = (uint64_t)(requested_size * grow_at_least_factor + 7) & 0xfffffffffffffff8ULL; // aligned to 8 bytes

    ASSERT (new_size >= requested_size, "Error in buf_alloc_do called from %s:%u: allocated too little memory for buffer %s: requested=%"PRIu64", allocated=%"PRIu64". vb_i=%u", 
            func, code_line, buf_desc (buf), requested_size, new_size, vb->vblock_i); // floating point paranoia

    // case 2: buffer was allocated already in the past - allocate new memory and copy over the data
    if (buf->memory) {

        uint64_t old_size = buf->size;

        // special handling if we have an overlaying buffer
        if (buf->overlayable) {
            pthread_mutex_lock (&overlay_mutex);
            uint16_t *overlay_count = (uint16_t*)(buf->data + buf->size + sizeof(uint64_t));

            char *old_data = buf->data;
            uint32_t old_len = buf->len;

            // if there is currently an overlay buffer on top of our buffer - abandon the memory
            // (leave it to the overlay buffer(s) that will eventually free() it), and allocate fresh memory
            if (*overlay_count > 1) {

                abandoned_mem_current += buf->size;
                abandoned_mem_high_watermark = MAX (abandoned_mem_high_watermark, abandoned_mem_current);

                (*overlay_count)--; // overlaying buffers are now on their own - no regular buffer
                buf->memory = buf->data = NULL;
                buf->size = buf->len = 0;
                buf_alloc_do (vb, buf, new_size, 1, func, code_line, name, param); // recursive call - simple alloc
          
                // copy old data
                memcpy (buf->data, old_data, old_size);
                buf->len = old_len;
            }
            else {
                // buffer is overlayable - but no current overlayers - regular realloc - however,
                // still within mutex to prevent another thread from overlaying while we're at it
                char *old_memory = buf->memory;
                __atomic_store_n (&buf->memory, BUFFER_BEING_MODIFIED, __ATOMIC_RELAXED);
                char *new_memory = (char *)buf_low_level_realloc (old_memory, new_size + overhead_size, func, code_line);
                buf_init (buf, new_memory, new_size, old_size, func, code_line, name, param);
            }
            buf->overlayable = true; // renew this, as it was reset by buf_init
            pthread_mutex_unlock (&overlay_mutex);
        }

        else { // non-overlayable buffer - regular realloc without mutex
            char *old_memory = buf->memory;
            __atomic_store_n (&buf->memory, BUFFER_BEING_MODIFIED, __ATOMIC_RELAXED);
            char *new_memory = (char *)buf_low_level_realloc (old_memory, new_size + overhead_size, func, code_line);
            buf_init (buf, new_memory, new_size, old_size, func, code_line, name, param);
        }
    }

    // case 3: we need to allocate memory - buffer is not yet allocated, so no need to copy data
    else {
        __atomic_store_n (&buf->memory, BUFFER_BEING_MODIFIED, __ATOMIC_RELAXED);
        char *memory = (char *)malloc (new_size + overhead_size);
        ASSERT (memory != BUFFER_BEING_MODIFIED, "Error in buf_alloc_do called from %s:%u: malloc didn't assign, very weird! buffer %s new_size=%"PRIu64,
                func, code_line, buf_desc(buf), new_size);

        buf->type   = BUF_REGULAR;

        buf_init (buf, memory, new_size, 0, func, code_line, name, param);
        buf_add_to_buffer_list(vb, buf);
    }

//    char size_str[20];
//    ASSERTW (new_size < 0x80000000ULL, "FYI: allocated > 2GB of memory buffer for %s: allocated=%s vb_i=%u", 
//             buf_desc (buf), str_size (new_size, size_str), vb->vblock_i);

finish:
    if (vb != evb) COPY_TIMER (buf_alloc); // this is not thread-safe for evb as evb buffers might be allocated by any thread
    return buf->size;
}

// an overlay buffer is a buffer using some of the memory of another buffer - it doesn't have its own memory
void buf_overlay_do (VBlock *vb, Buffer *overlaid_buf, Buffer *regular_buf, const char *func, uint32_t code_line,
                     const char *name, int64_t param)
{
    // if this buffer was used by a previous VB as a regular buffer - we need to "destroy" it first
    if (overlaid_buf->type == BUF_REGULAR && overlaid_buf->data == NULL && overlaid_buf->memory) {
        buf_low_level_free (overlaid_buf->memory, func, code_line);
        overlaid_buf->type = BUF_UNALLOCATED;
    }
    
    ASSERT (overlaid_buf->type == BUF_UNALLOCATED, "Error in %s:%u: cannot buf_overlay to a buffer %s already in use", func, code_line, buf_desc (overlaid_buf));
    ASSERT (regular_buf->type == BUF_REGULAR, "Error in %s:%u: regular_buf %s in buf_overlay must be a regular buffer", func, code_line, buf_desc (regular_buf));
    ASSERT (regular_buf->overlayable, "Error in %s:%u: buf_overlay: buffer %s is not overlayble", func, code_line, buf_desc (regular_buf));

    overlaid_buf->size        = 0;
    overlaid_buf->len         = 0;
    overlaid_buf->type        = BUF_OVERLAY;
    overlaid_buf->memory      = 0;
    overlaid_buf->overlayable = false;
    overlaid_buf->vb          = vb;
    overlaid_buf->name        = name;
    overlaid_buf->param       = param;

    // full buffer overlay - copy len too and update overlay counter
    pthread_mutex_lock (&overlay_mutex);

    overlaid_buf->size = regular_buf->size;
    overlaid_buf->len  = regular_buf->len;
    overlaid_buf->data = regular_buf->data;
    uint16_t *overlay_count = (uint16_t*)(regular_buf->data + regular_buf->size + sizeof(uint64_t));
    (*overlay_count)++; // counter of users of this memory

    pthread_mutex_unlock (&overlay_mutex);
}

// free buffer - without freeing memory. A future buf_alloc of this buffer will reuse the memory if possible.
void buf_free_do (Buffer *buf, const char *func, uint32_t code_line) 
{
    uint16_t *overlay_count; // number of buffers (overlay and regular) sharing buf->memory

    switch (buf->type) {

        case BUF_UNALLOCATED:
            return; // nothing to do

        case BUF_REGULAR: 

            if (buf->overlayable) {

                pthread_mutex_lock (&overlay_mutex);
                overlay_count = (uint16_t*)(buf->data + buf->size + sizeof(uint64_t));

                if (*overlay_count > 1) { // current overlays exist - abandon memory - leave it to the overlaid buffer(s) which will free() this memory when they're done with it
                    (*overlay_count)--;
             
                    abandoned_mem_current += buf->size;
                    abandoned_mem_high_watermark = MAX (abandoned_mem_high_watermark, abandoned_mem_current);

                    buf_reset (buf);
                }
                // if no overlay exists then we just keep .memory and reuse it in future allocations

                pthread_mutex_unlock (&overlay_mutex);            
            }
            
            buf->data        = NULL; 
            buf->len         = 0;
            buf->param       = 0;
            buf->overlayable = false;
            
            // name, param, memory and size are not changed

            break;

        case BUF_OVERLAY:
            pthread_mutex_lock (&overlay_mutex);
            overlay_count = (uint16_t*)(buf->data + buf->size + sizeof(uint64_t));
            (*overlay_count)--;
            pthread_mutex_unlock (&overlay_mutex);            

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

        default:
            ABORT0 ("Error: invalid buf->type");
    }
} 

// remove from buffer_list of this vb
// thread safety: only the thread owning the VB of the buffer (I/O thread of evb) can remove a buffer
// from the buf list, OR the I/O thread may remove for all VBs, IF no compute thread is running
void buf_remove_from_buffer_list (Buffer *buf)
{
    if (!buf->vb) return; // this buffer is not on the buf_list - nothing to do

    ARRAY (Buffer *, buf_list, buf->vb->buffer_list);

     for (unsigned i=0; i < buf->vb->buffer_list.len; i++) 

        if (buf_list[i] == buf) {
            
            if (flag_debug_memory) {
                char s[POINTER_STR_LEN];
                fprintf (stderr, "Destroy %s: buf_addr=%s buf->vb->id=%d buf_i=%u\n", buf_desc (buf), str_pointer(buf,s), buf->vb->id, i);
            }

            buf_list[i] = NULL;
            buf->vb = NULL;

            break;
        }

    // note: it is possible that the buffer is not found in the list if it is never allocated. that's fine.
}

void buf_destroy_do (Buffer *buf, const char *func, uint32_t code_line)
{
    if (!buf) return; // nothing to do

    buf_remove_from_buffer_list (buf); 

    if (buf->memory) {
    
        uint16_t overlay_count = 1;
        if (buf->overlayable) {
            pthread_mutex_lock (&overlay_mutex);
            overlay_count = (*(uint16_t*)(buf->data + buf->size + sizeof(uint64_t)));
            pthread_mutex_unlock (&overlay_mutex);            
        }

        ASSERT (overlay_count==1, "Error: cannot destroy buffer %s because it is currently overlaid", buf->name);

        buf_low_level_free (buf->memory, func, code_line);
    }

    memset (buf, 0, sizeof (Buffer)); // reset to factory defaults
}

void buf_copy_do (VBlock *dst_vb, Buffer *dst, const Buffer *src, 
                  uint64_t bytes_per_entry, // how many bytes are counted by a unit of .len
                  uint64_t src_start_entry, uint64_t max_entries,  // if 0 copies the entire buffer 
                  const char *func, unsigned code_line,
                  const char *name, int64_t param)
{
    ASSERT (src->data, "Error in buf_copy call from %s:%u: src->data is NULL", func, code_line);
    
    ASSERT (!max_entries || src_start_entry < src->len, 
            "Error buf_copy of %s called from %s:%u: src_start_entry=%"PRIu64" is larger than src->len=%"PRIu64, buf_desc(src), func, code_line, src_start_entry, src->len);

    uint64_t num_entries = max_entries ? MIN (max_entries, src->len - src_start_entry) : src->len - src_start_entry;
    if (!bytes_per_entry) bytes_per_entry=1;
    
    buf_alloc_do (dst_vb, dst, num_entries * bytes_per_entry, 1, func, code_line,
                  name ? name : src->name, name ? param : src->param); // use realloc rather than malloc to allocate exact size

    memcpy (dst->data, &src->data[src_start_entry * bytes_per_entry], num_entries * bytes_per_entry);

    dst->len = num_entries;  
}   

// moves all the data from one buffer to another, leaving the source buffer unallocated
// both the src and dst buffer remain (or are added) to the buf_lists of their respective VBs
void buf_move (VBlock *dst_vb, Buffer *dst, VBlock *src_vb, Buffer *src)
{
    if (dst->type == BUF_REGULAR && !dst->data) buf_destroy (dst);
    
    ASSERT (dst->type == BUF_UNALLOCATED, "Error: attempt to move to an already-allocated src: %s dst: %s", buf_desc (src), buf_desc (dst));

    ASSERT (src_vb==dst_vb || dst->vb==dst_vb, "Error in buf_move: to move a buffer between VBs, the dst buffer needs to be added"
                                               " to the dst_vb buffer_list in advance. If dst_vb=evb the dst buffer must be added to"
                                               " the buffer_list by the I/O thread only. src: %s dst: %s src_vb->vb_i=%d dst_vb->vb_i=%d",
            buf_desc (src), buf_desc (dst), (src_vb ? src_vb->vblock_i : -999), (dst_vb ? dst_vb->vblock_i : -999));

    if (!dst->vb) buf_add_to_buffer_list (dst_vb, dst); // this can only happen if src_vb==dst_vb 

    memcpy (dst, src, sizeof(Buffer));    
    dst->vb = dst_vb;

    buf_reset (src); // zero buffer except vb
}

void buf_add_string (VBlock *vb, Buffer *buf, const char *str) 
{ 
    unsigned len = strlen (str); 
    buf_alloc (vb, buf, MAX (1000, buf->len + len + 1), 2, "string_buf", 0);
    buf_add (buf, str, len);
    buf->data[buf->len] = '\0'; // string terminator without increasing buf->len
}

void buf_print (Buffer *buf, bool add_newline)
{
    for (uint64_t i=0; i < buf->len; i++) 
        putchar (buf->data[i]);  // safer than printf %.*s ?

    if (add_newline) putchar ('\n');
}

void buf_low_level_free (void *p, const char *func, uint32_t code_line)
{
    if (flag_debug_memory) {
        char s[POINTER_STR_LEN];
        fprintf (stderr, "Memory freed by free(): %s %s:%u\n", str_pointer (p, s), func, code_line);
    }

    free (p);
}

void *buf_low_level_realloc (void *p, size_t size, const char *func, uint32_t code_line)
{
    void *new = realloc (p, size);

    if (flag_debug_memory && new != p) {
        char s[POINTER_STR_LEN];
        fprintf (stderr, "Memory freed by realloc(): %s %s:%u\n", str_pointer (p, s), func, code_line);
    }

    return new;
}

// convert a Buffer from a z_file section whose len is in char to a bitarray
BitArray *buf_zfile_buf_to_bitarray (Buffer *buf, uint64_t num_of_bits)
{
    ASSERT (roundup_bits2bytes (num_of_bits) <= buf->len, "Error in buf_zfile_buf_to_bitarray: num_of_bits=%"PRId64" indicating a length of at least %"PRId64", but buf->len=%"PRId64,
            num_of_bits, roundup_bits2bytes (num_of_bits), buf->len);

    BitArray *bitarr = buf_get_bitarray (buf);
    bitarr->num_of_bits  = num_of_bits;
    bitarr->num_of_words = roundup_bits2words64 (bitarr->num_of_bits);

    ASSERT (roundup_bits2bytes64 (num_of_bits) <= buf->size, "Error in buf_zfile_buf_to_bitarray: buffer to small: buf->size=%"PRId64" but bitarray has %"PRId64" words and hence requires %"PRId64" bytes",
            buf->size, bitarr->num_of_words, bitarr->num_of_words * sizeof(uint64_t));

    LTEN_bit_array (bitarr);

    return bitarr;
}

void buf_add_bit (Buffer *buf, int64_t new_bit) 
{
    BitArray *bar = buf_get_bitarray (buf);

    ASSERT (bar->num_of_bits < buf->size * 8, "Error in %s:%u: no room in Buffer %s to extend the bitmap", __FUNCTION__, __LINE__, buf->name);
    bar->num_of_bits++;     
    if (bar->num_of_bits % 64 == 1) { // starting a new word                
        bar->num_of_words++;
        bar->words[bar->num_of_words-1] = new_bit; // LSb is as requested, other 63 bits are 0
    } 
    else
        bit_array_assign (bar, bar->num_of_bits-1, new_bit);  
}
 
bit_index_t buf_extend_bits (Buffer *buf, int64_t num_new_bits) 
{
    BitArray *bar = buf_get_bitarray (buf);

    ASSERT (bar->num_of_bits + num_new_bits<= buf->size * 8, "Error in %s:%u: no room in Buffer %s to extend the bitmap", __FUNCTION__, __LINE__, buf->name);
    
    bit_index_t next_bit = bar->num_of_bits;

    bar->num_of_bits += num_new_bits;     
    bar->num_of_words = roundup_bits2words64 (bar->num_of_bits);
    bit_array_clear_excess_bits_in_top_word (bar);

    return next_bit;
}
 
void BGEN_u8_buf (Buffer *buf)
{
}

void BGEN_u16_buf (Buffer *buf)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint16_t num_big_en = *ENT (uint16_t, *buf, i);
        *ENT (uint16_t, *buf, i) = BGEN16 (num_big_en);            
    }
}

void BGEN_u32_buf (Buffer *buf)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint32_t num_big_en = *ENT (uint32_t, *buf, i);
        *ENT (uint32_t, *buf, i) = BGEN32 (num_big_en);            
    }
}

void BGEN_u64_buf (Buffer *buf)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint64_t num_big_en = *ENT (uint64_t, *buf, i);
        *ENT (uint64_t, *buf, i) = BGEN64 (num_big_en);            
    }
}

void BGEN_deinterlace_d8_buf (Buffer *buf)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint8_t unum = *ENT (uint8_t, *buf, i);
        *ENT (int8_t, *buf, i) = DEINTERLACE(int8_t,unum); 
   }
}

void BGEN_deinterlace_d16_buf (Buffer *buf)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint16_t num_big_en = *ENT (uint16_t, *buf, i);
        uint16_t unum = BGEN16 (num_big_en);
        *ENT (int16_t, *buf, i) = DEINTERLACE(int16_t,unum); 
   }
}

void BGEN_deinterlace_d32_buf (Buffer *buf)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint32_t num_big_en = *ENT (uint32_t, *buf, i);
        uint32_t unum = BGEN32 (num_big_en);
        *ENT (int32_t, *buf, i) = DEINTERLACE(int32_t,unum); 
   }
}

void BGEN_deinterlace_d64_buf (Buffer *buf)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint64_t num_big_en = *ENT (uint64_t, *buf, i);
        uint64_t unum = BGEN64 (num_big_en);
        *ENT (int64_t, *buf, i) = DEINTERLACE(int64_t,unum); 
   }
}
