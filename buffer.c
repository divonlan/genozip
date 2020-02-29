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
#include "vb.h"
#ifndef _MSC_VER // Microsoft compiler
#include <pthread.h>
#else
#include "compatability/visual_c_pthread.h"
#endif

//#define DISPLAY_ALLOCS_AFTER 4090 // display allocations, except the first X allocations. reallocs are always displayed

#define UNDERFLOW_TRAP 0x574F4C4652444E55ULL // "UNDRFLOW" - inserted at the begining of each memory block to detected underflows
#define OVERFLOW_TRAP  0x776F6C667265766FULL // "OVERFLOW" - inserted at the end of each memory block to detected overflows

static const unsigned overhead_size = 2*sizeof (uint64_t) + sizeof(uint16_t); // underflow, overflow and user counter

static pthread_mutex_t overlay_mutex; // used to thread-protect overlay counters
static uint64_t abandoned_mem_current = 0;
static uint64_t abandoned_mem_high_watermark = 0;

void buf_initialize()
{
    pthread_mutex_init (&overlay_mutex, NULL);
}

char *buf_human_readable_size (int64_t size, char *str /* out */)
{
    if      (size > (1LL << 40)) sprintf (str, "%3.1lf TB", ((double)size) / (double)(1LL << 40));
    else if (size > (1LL << 30)) sprintf (str, "%3.1lf GB", ((double)size) / (double)(1LL << 30));
    else if (size > (1LL << 20)) sprintf (str, "%3.1lf MB", ((double)size) / (double)(1LL << 20));
    else if (size > (1LL << 10)) sprintf (str, "%3.1lf KB", ((double)size) / (double)(1LL << 10));
    else                         sprintf (str, "%3d B"    ,     (int)size)                       ;

    return str; // for convenience so caller can use in printf directly
}

// get string with buffer's metadata for debug message. this function is NOT thread-safe
char *buf_display (const Buffer *buf)
{
    static char str[200]; // NOT thread-safe

    sprintf (str, 
#ifdef _MSC_VER
             "Buffer %s (%u): size=%u len=%u data=0x%I64x memory=0x%I64x",
#else
             "Buffer %s (%u): size=%u len=%u data=0x%"PRIx64" memory=0x%"PRIx64,
#endif
             buf->name, buf->param, buf->size, buf->len, (uint64_t)buf->data, (uint64_t)buf->memory);
    return str;    
}

static inline bool buf_has_overflowed (const Buffer *buf)
{
    return *((uint64_t*)(buf->memory + buf->size + sizeof(uint64_t))) != OVERFLOW_TRAP;
}

static inline bool buf_has_underflowed (const Buffer *buf)
{
    return *(uint64_t*)buf->memory != UNDERFLOW_TRAP;
}

void buf_test_overflows(const VariantBlock *vb)
{
    const Buffer *buf_list = &vb->buffer_list;

    bool corruption = false;
    for (unsigned buf_i=0; buf_i < buf_list->len; buf_i++) {
        const Buffer *buf = ((Buffer **)buf_list->data)[buf_i];

        if (!buf) continue; // buf was 'buf_destroy'd

        if (buf->memory) {
            if (!buf->name) {
                fprintf (stderr, 
#ifdef _MSC_VER
                         "buffer=0x%I64x : Corrupt Buffer structure - null name\n", 
#else
                         "buffer=0x%"PRIx64" : Corrupt Buffer structure - null name\n", 
#endif
                         (uint64_t)(uintptr_t)buf);
                corruption = true;
            }
            else if (buf->data && buf->data != buf->memory + sizeof(uint64_t)) {
                fprintf (stderr, 
#ifdef _MSC_VER
                         "vb_id=%u buf_i=%u buffer=0x%I64x memory=0x%I64x : Corrupt Buffer structure - expecting data+8 == memory. name=%.30s param=%u buf->data=0x%I64x\n", 
#else
                         "vb_id=%u buf_i=%u buffer=0x%"PRIx64" memory=0x%"PRIx64" : Corrupt Buffer structure - expecting data+8 == memory. name=%.30s param=%u buf->data=0x%"PRIx64"\n", 
#endif
                         vb ? vb->id : 0, buf_i, (uint64_t)(uintptr_t)buf, (uint64_t)(uintptr_t)buf->memory, buf->name, buf->param, (uint64_t)(uintptr_t)buf->data);
                corruption = true;
            }
            else if (buf_has_underflowed(buf)) {
                fprintf (stderr, 
#ifdef _MSC_VER
                        "vb_id=%u buf_i=%u buffer=0x%I64x memory=0x%I64x : Underflow in buffer %.30s param=%u fence=\"%c%c%c%c%c%c%c%c\"\n", 
#else
                        "vb_id=%u buf_i=%u buffer=0x%"PRIx64" memory=0x%"PRIx64" : Underflow in buffer %.30s param=%u fence=\"%c%c%c%c%c%c%c%c\"\n", 
#endif
                         vb ? vb->id : 0, buf_i, (uint64_t)(uintptr_t)buf, (uint64_t)(uintptr_t)buf->memory, buf->name, buf->param, 
                         buf->memory[0], buf->memory[1], buf->memory[2], buf->memory[3], buf->memory[4], buf->memory[5], buf->memory[6], buf->memory[7]);
                corruption = true;
            }
            else if (buf_has_overflowed(buf)) {
                char *of = &buf->memory[buf->size + sizeof(uint64_t)];
                fprintf (stderr,
#ifdef _MSC_VER
                         "vb_id=%u buf_i=%u buffer=0x%I64x memory=0x%I64x size=%u : Overflow in buffer %.30s param=%u fence=\"%c%c%c%c%c%c%c%c\"\n", 
#else
                         "vb_id=%u buf_i=%u buffer=0x%"PRIx64" memory=0x%"PRIx64" size=%u : Overflow in buffer %.30s param=%u fence=\"%c%c%c%c%c%c%c%c\"\n", 
#endif
                         vb ? vb->id : 0, buf_i, (uint64_t)(uintptr_t)buf, (uint64_t)(uintptr_t)buf->memory, buf->size, buf->name, buf->param, 
                         of[0], of[1], of[2], of[3], of[4], of[5], of[6], of[7]);
                
                corruption = true;
            }
        }
    }
    ASSERT0 (!corruption, "Aborting due to memory corruption");
}

typedef struct {
    const char *name;
    uint64_t bytes; 
    unsigned buffers;
} MemStats;

static int buf_stats_sort_by_bytes(const void *a, const void *b)  
{ 
    return ((MemStats*)a)->bytes < ((MemStats*)b)->bytes ? 1 : -1;
}

void buf_display_memory_usage (PoolId pool_id, bool memory_full)
{
    #define MAX_MEMORY_STATS 100
    static MemStats stats[MAX_MEMORY_STATS]; // must be pre-allocated, because buf_display_memory_usage is called when malloc fails, so it cannot malloc
    unsigned num_stats = 0, num_buffers = 0;

    if (memory_full)
        fprintf (stderr, "\n\nError memory is full:\n");
    else
        fprintf (stderr, "\n-------------------------------------------------------------------------------------\n");

    VariantBlockPool *vb_pool = vb_get_pool (pool_id);

    for (unsigned vb_i=0; vb_i < vb_pool->num_vbs; vb_i++) {

        Buffer *buf_list = &vb_pool->vb[vb_i].buffer_list; // a pointer to a buffer, which contains an array of pointers to buffers of a single vb/non-vb

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

    stats[num_stats].name    = "abandoned_mem_high_watermark";
    stats[num_stats].bytes   = abandoned_mem_high_watermark;
    stats[num_stats].buffers = 0;
    num_stats++;

    // sort stats by bytes
    qsort (stats, num_stats, sizeof (MemStats), buf_stats_sort_by_bytes);

    uint64_t total_bytes=0;
    for (unsigned i=0; i< num_stats; i++) total_bytes += stats[i].bytes;

    char str[30];
    buf_human_readable_size (total_bytes, str);
    fprintf (stderr, "Total bytes: %s in %u buffers in %u buffer lists:\n", str, num_buffers, vb_pool->num_vbs);

    for (unsigned i=0; i < num_stats; i++) {
        buf_human_readable_size (stats[i].bytes, str);
        fprintf (stderr, "%-30s: %-8s (%4.1f%%) in %u buffers\n", stats[i].name, str, 100.0 * (float)stats[i].bytes / (float)total_bytes, stats[i].buffers);
    }
}

long long buf_vb_memory_consumption (const VariantBlock *vb)
{
    const Buffer *buf_list   = &vb->buffer_list;

    // memory of the structure itself
    long long vb_memory = sizeof (*vb);

    if (vb->data_lines) vb_memory += global_max_lines_per_vb * sizeof (DataLine);

    // memory allocated outside of Buffer (direct calloc)
    if (vb->haplotype_sections_data) vb_memory += vb->num_sample_blocks * sizeof (Buffer);
    if (vb->genotype_sections_data)  vb_memory += vb->num_sample_blocks * sizeof (Buffer);
    if (vb->phase_sections_data)     vb_memory += vb->num_sample_blocks * sizeof (Buffer);

    // memory allocated via Buffer (the vast majority, usually)
    for (unsigned buf_i=0; buf_i < buf_list->len; buf_i++) {
        Buffer *buf = ((Buffer **)buf_list->data)[buf_i]; 
        
        if (!buf || !buf->memory) continue; // exclude destroyed or not-yet-allocated buffers
        
        vb_memory += buf->size;
    }

    return vb_memory;
}

static inline void buf_add_to_buffer_list (VariantBlock *vb, Buffer *buf)
{
#define INITIAL_MAX_MEM_NUM_BUFFERS 10000 /* for files that have ht,gt,phase,variant,and line - the factor would be about 5.5 so there will be 1 realloc per vb, but most files don't */
    Buffer *bl = &vb->buffer_list;

    buf_alloc (vb, bl, MAX (INITIAL_MAX_MEM_NUM_BUFFERS, bl->len+1) * sizeof(Buffer *), 2, "buffer_list", vb->id);

    ((Buffer **)bl->data)[bl->len++] = buf;
}

static void buf_init (VariantBlock *vb, Buffer *buf, unsigned size, unsigned old_size, const char *name, unsigned param)
{
    if (!buf->memory) {
        buf_test_overflows(vb);
#ifdef DEBUG
        buf_display_memory_usage (vb->pool_id, true);
#endif
        ABORT ("Error: Failed to allocate %u bytes name=%s param=%u", size + overhead_size, name, param);
    }

    buf->data        = buf->memory + sizeof (uint64_t);
    buf->size        = size;
    buf->overlayable = false;

    if (name) {
        buf->name = name;
        buf->param = param;
    } 
    ASSERT0 (buf->name, "Error: buffer has no name");

    *(uint64_t *)buf->memory        = UNDERFLOW_TRAP;        // underflow protection
    *(uint64_t *)(buf->data + size) = OVERFLOW_TRAP;         // overflow prortection (underflow protection was copied with realloc)
    *(uint16_t *)(buf->data + size + sizeof (uint64_t)) = 1; // counter of buffers that use of this memory (0 or 1 main buffer + any number of overlays)

#ifdef DISPLAY_ALLOCS_AFTER
    if (vb->buffer_list.len > DISPLAY_ALLOCS_AFTER)
        fprintf (stderr, "%s (%u): old_size=%u new_size=%u\n", buf->name, buf->param, old_size, size);
#endif
}

// allocates or enlarges buffer
// if it needs to enlarge a buffer fully overlaid by an overlay buffer - it abandons its memory (leaving it to
// the overlaid buffer) and allocates new memory
unsigned buf_alloc_do (VariantBlock *vb,
                       Buffer *buf, 
                       unsigned requested_size,
                       float grow_at_least_factor, // IF we need to allocate or reallocate physical memory, we get this much more than requested
                       const char *name, unsigned param)      
{
    START_TIMER;

    if (!requested_size) return 0; // nothing to do

    // sanity checks
    ASSERT (buf->type == BUF_REGULAR || buf->type == BUF_UNALLOCATED, "Error: cannot buf_alloc an overlayed buffer. name=%s", buf->name ? buf->name : "");
    ASSERT0 (vb, "Error: null vb");

    // case 1: we have enough memory already
    if (requested_size <= buf->size) {
        buf->data = buf->memory + sizeof (uint64_t); // allocate if not already allocated
        goto finish;
    }

    // grow us requested - rounding up to 64 bit boundary to avoid aliasing errors with the overflow indicator
    uint32_t new_size = (uint32_t)(requested_size * MAX (grow_at_least_factor, 1) + 7) & 0xfffffff8; 

    // case 2: we need to allocate memory - buffer is already allocated so copy over the data
    if (buf->memory) {

        unsigned old_size = buf->size;

        // special handing if we have an overlaying buffer
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
                buf_alloc (vb, buf, new_size, 1, name, param); // recursive call - simple alloc
                
                // copy old data
                memcpy (buf->data, old_data, old_size);
                buf->len = old_len;
            }
            else {
                // buffer is overlayable - but no current overlayers - regular realloc - however,
                // still within mutex to prevent another thread from overlaying while we're at it
                buf->memory = (char *)realloc (buf->memory, new_size + overhead_size);
                buf_init (vb, buf, new_size, old_size, name, param);
            }
            buf->overlayable = true;
            pthread_mutex_unlock (&overlay_mutex);
        }

        else { // non-overlayable buffer - regular realloc without mutex
            buf->memory = (char *)realloc (buf->memory, new_size + overhead_size);
            buf_init (vb, buf, new_size, old_size, name, param);
        }
    }

    // case 3: we need to allocate memory - buffer is not yet allocated, so no need to copy data
    else {
        buf->memory = (char *)malloc (new_size + overhead_size);
        buf->type  = BUF_REGULAR;

        buf_init (vb, buf, new_size, 0, name, param);
        buf_add_to_buffer_list(vb, buf);
    }

finish:
    if (vb) COPY_TIMER (vb->profile.buf_alloc);
    return buf->size;
}

// an overlay buffer is a buffer using some of the memory of another buffer - it doesn't have its own memory
void buf_overlay (Buffer *overlaid_buf, Buffer *regular_buf, const Buffer *copy_from /* optional */, 
                  unsigned *regular_buf_offset, const char *name, unsigned param)
{
    bool full_overlay = !regular_buf_offset && !copy_from;

//printf ("Overlaying onto buffer old_name=%s old_param=%u new_name=%s new_param=%u\n", 
//        overlaid_buf->name, overlaid_buf->param, regular_buf->name, regular_buf->param);      

    ASSERT (overlaid_buf->type == BUF_UNALLOCATED, "Error: cannot buf_overlay to a buffer already in use. overlaid_buf->name=%s", overlaid_buf->name ? overlaid_buf->name : "");
    ASSERT (regular_buf->type == BUF_REGULAR, "Error: regular_buf in buf_overlay must be a regular buffer. regular_buf->name=%s", regular_buf->name ? regular_buf->name : "");
    ASSERT (!full_overlay || regular_buf->overlayable, "Error: buf_overlay: only overlayble buffers can be fully overlaid. regular_buf->name=%s", regular_buf->name ? regular_buf->name : "");

    overlaid_buf->size   = 0;
    overlaid_buf->len    = copy_from ? copy_from->len : 0;
    overlaid_buf->type   = full_overlay ? BUF_FULL_OVERLAY : BUF_PARTIAL_OVERLAY;
    overlaid_buf->memory = 0;
    overlaid_buf->overlayable = false;

    if (name) {
        overlaid_buf->name   = name;
        overlaid_buf->param  = param;
    }
    else {
        overlaid_buf->name  = regular_buf->name;
        overlaid_buf->param = regular_buf->param;
    }

    if (!full_overlay) {
        overlaid_buf->data = regular_buf->data + (regular_buf_offset ? *regular_buf_offset : 0);

        if (copy_from && copy_from->len) {

            ASSERT ((regular_buf_offset ? *regular_buf_offset : 0) + copy_from->len <= regular_buf->size, 
                    "Error: buf_overlay exceeds the size of the regular buffer: offset=%u size=%u regular_buf.size=%u", 
                    *regular_buf_offset, copy_from->len, regular_buf->size);
        
            memcpy (overlaid_buf->data, copy_from->data, copy_from->len);

            // new data was copied to overlaid buffer - move the offset forward and update the len of the regular buffer
            // to enable to the next buffer to be overlaid subsequently
            if (regular_buf_offset) {
                *regular_buf_offset += overlaid_buf->len;
                regular_buf->len = *regular_buf_offset;
            }
        }
    }

    // full buffer overlay - copy len too and update overlay counter
    else {
        pthread_mutex_lock (&overlay_mutex);

        overlaid_buf->size = regular_buf->size;
        overlaid_buf->len  = regular_buf->len;
        overlaid_buf->data = regular_buf->data;
        uint16_t *overlay_count = (uint16_t*)(regular_buf->data + regular_buf->size + sizeof(uint64_t));
        (*overlay_count)++; // counter of users of this memory

        pthread_mutex_unlock (&overlay_mutex);
    }
}

// free buffer - without freeing memory. A future buf_alloc of this buffer will reuse the memory if possible.
void buf_free (Buffer *buf) 
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

                    memset (buf, 0, sizeof (Buffer)); // make this buffer UNALLOCATED
                }
                // if no overlay exists then we just keep .memory and reuse it in future allocations

                pthread_mutex_unlock (&overlay_mutex);            
            }
            
            buf->data = NULL; 
            buf->len = 0;
            buf->overlayable = false;
            
            // name, param, memory and size are not changed

            break;

        case BUF_FULL_OVERLAY:
//printf ("Freeing overlay buffer name=%s param=%u\n", buf->name, buf->param);      
            pthread_mutex_lock (&overlay_mutex);
            overlay_count = (uint16_t*)(buf->data + buf->size + sizeof(uint64_t));
            (*overlay_count)--;
            pthread_mutex_unlock (&overlay_mutex);            

            // we are the last user - we can free the memory now.
            // do this outside of the mutex - free is a system call and can take some time.
            // this is safe because if we ever observe *overlay_count==0, it means that no buffer has this memory,
            // therefore there is no possibility it would be subsequently overlayed between the test and the free().
            if (! (*overlay_count)) {
                free (buf->data - sizeof(uint64_t)); // the original buf->memory
                abandoned_mem_current -= buf->size;
            }
            
            // fall through

        case BUF_PARTIAL_OVERLAY:
            memset (buf, 0, sizeof (Buffer));
            break;

        default:
            ABORT0 ("Error: invalid buf->type");
    }
} 

void buf_destroy (VariantBlock *vb, Buffer *buf)
{
    if (!buf) return; // nothing to do

    if (buf->memory) {

        uint16_t overlay_count = 1;
        if (buf->overlayable) {
            pthread_mutex_lock (&overlay_mutex);
            overlay_count = (*(uint16_t*)(buf->data + buf->size + sizeof(uint64_t)));
            pthread_mutex_unlock (&overlay_mutex);            
        }

        ASSERT (overlay_count==1, "Error: cannot destroy buffer %s because it is currently overlaid", buf->name);

        free (buf->memory);
    }

    // remove from buffer_list of this vb
    unsigned i=0; for (; i < vb->buffer_list.len; i++) 
        if (((Buffer **)vb->buffer_list.data)[i] == buf) {
            ((Buffer **)vb->buffer_list.data)[i] = NULL;
            break;
        }

    // sanity (note: it is possible that a buffer is not in the list if it was never allocated)
    ASSERT0 (!buf->memory || i < vb->buffer_list.len, "Error: cannot find buffer in buffer_list");
}

void buf_copy (VariantBlock *vb, Buffer *dst, const Buffer *src, 
               unsigned bytes_per_entry, // how many bytes are counted by a unit of .len
               unsigned start_entry, unsigned max_entries,  // if 0 copies the entire buffer 
               const char *name, unsigned param)
{
    ASSERT0 (src->data, "Error in buf_copy: src->data is NULL");
    
    unsigned num_entries = max_entries ? MIN (max_entries, src->len - start_entry) : src->len - start_entry;

    buf_alloc(vb, dst, num_entries * bytes_per_entry, 1, 
              name ? name : src->name, name ? param : src->param); // use realloc rather than malloc to allocate exact size

    memcpy (dst->data, &src->data[start_entry * bytes_per_entry], num_entries * bytes_per_entry);

    dst->len = num_entries;  
}   

// moves all the data from one buffer to another, leaving the source buffer unallocated
void buf_move (VariantBlock *vb, Buffer *dst, Buffer *src)
{
    memcpy (dst, src, sizeof(Buffer));
    memset (src, 0, sizeof(Buffer));

    buf_add_to_buffer_list (vb, dst);
}

void buf_set_overlayable (Buffer *buf) 
{ 
    buf->overlayable = true;
}

void buf_add_string (VariantBlock *vb, Buffer *buf, const char *str) 
{ 
    unsigned len = strlen (str); 
    buf_alloc (vb, buf, MAX (1000, buf->len + len + 1), 2, "string_buf", 0);
    buf_add (buf, str, len);
    buf->data[buf->len] = '\0'; // string terminator without increasing buf->len
}
