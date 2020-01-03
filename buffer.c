// ------------------------------------------------------------------
//   buffer.c
//   Copyright (C) 2019 Divon Lan <genozip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// memory management - when running the same code by the same thread for another variant block - we reuse
// the previous variant's block memory. this way we save repetitive malloc/free cycles which might
// be very time consuming.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <inttypes.h>

#include "genozip.h"

//#define DISPLAY_ALLOCS_AFTER 4100 // display allocations, except the first X allocations. reallocs are always displayed

#define UNDERFLOW_TRAP 0x574F4C4652444E55ULL // "UNDRFLOW" - inserted at the begining of each memory block to detected underflows
#define OVERFLOW_TRAP  0x776F6C667265766FULL // "OVERFLOW" - inserted at the end of each memory block to detected overflows

static Buffer **buffer_lists = NULL; 
static unsigned num_buffer_lists = 0;

char *buf_human_readable_size (uint64_t size, char *str /* out */)
{
    if      (size > (1ULL << 40)) sprintf (str, "%3.1lf TB", ((double)size) / (double)(1ULL << 40));
    else if (size > (1ULL << 30)) sprintf (str, "%3.1lf GB", ((double)size) / (double)(1ULL << 30));
    else if (size > (1ULL << 20)) sprintf (str, "%3.1lf MB", ((double)size) / (double)(1ULL << 20));
    else if (size > (1ULL << 10)) sprintf (str, "%3.1lf KB", ((double)size) / (double)(1ULL << 10));
    else                          sprintf (str, "%3u B"    ,(unsigned)size)                        ;

    return str; // for convenience so caller can use in printf directly
}

bool buf_has_overflowed(const Buffer *buf)
{
    return *((long long*)(buf->memory + buf->size + sizeof(long long))) != OVERFLOW_TRAP;
}

bool buf_has_underflowed(const Buffer *buf)
{
    return *(long long*)buf->memory != UNDERFLOW_TRAP;
}

void buf_test_overflows(const VariantBlock *vb)
{
    const Buffer *buf_list = &vb->buffer_list;

    bool corruption = false;
    for (unsigned buf_i=0; buf_i < buf_list->len; buf_i++) {
        const Buffer *buf = ((Buffer **)buf_list->data)[buf_i];

        if (buf->memory) {
            if (!buf->name) {
                fprintf (stderr, "buffer=0x%"PRIx64" : Corrupt Buffer structure - null name\n", (uint64_t)(uintptr_t)buf);
                corruption = true;
            }
            else if (buf->data && buf->data != buf->memory + sizeof(long long)) {
                fprintf (stderr, "vb_id=%u buf_i=%u buffer=0x%"PRIx64" memory=0x%"PRIx64" : Corrupt Buffer structure - expecting data+8 == memory. name=%.30s param=%u buf->data=0x%"PRIx64, 
                         vb ? vb->id : 0, buf_i, (uint64_t)(uintptr_t)buf, (uint64_t)(uintptr_t)buf->memory, buf->name, buf->param, (uint64_t)(uintptr_t)buf->data);
                corruption = true;
            }
            else if (buf_has_underflowed(buf)) {
                fprintf (stderr, "vb_id=%u buf_i=%u buffer=0x%"PRIx64" memory=0x%"PRIx64" : Underflow in buffer %.30s param=%u \"%.8s\"", 
                         vb ? vb->id : 0, buf_i, (uint64_t)(uintptr_t)buf, (uint64_t)(uintptr_t)buf->memory, buf->name, buf->param, buf->memory);
                corruption = true;
            }
            else if (buf_has_overflowed(buf)) {
                fprintf (stderr,"vb_id=%u buf_i=%u buffer=0x%"PRIx64" memory=0x%"PRIx64" size=%u : Overflow in buffer %.30s param=%u \"%.8s\"", 
                         vb ? vb->id : 0, buf_i, (uint64_t)(uintptr_t)buf, (uint64_t)(uintptr_t)buf->memory, buf->size, buf->name, buf->param, &buf->memory[buf->size + sizeof(long long)]);
                
                corruption = true;
            }
        }
    }
    ASSERT0 (!corruption, "Aborting due to memory corruption");
}

typedef struct {
    const char *name;
    unsigned bytes, buffers;
} MemStats;

static int buf_stats_sort_by_bytes(const void *a, const void *b)  
{ 
    return ((MemStats*)a)->bytes < ((MemStats*)b)->bytes ? 1 : -1;
}

void buf_display_memory_usage(bool memory_full)
{
    #define MAX_MEMORY_STATS 1000
    static MemStats stats[MAX_MEMORY_STATS]; // must be pre-allocated, because buf_display_memory_usage is called when malloc fails, so it cannot malloc
    unsigned num_stats = 0, num_buffers = 0;

    if (memory_full)
        printf ("\n\nError memory is full:\n");
    else
        printf ("\n-------------------------------------------------------------------------------------\n");

    for (unsigned buf_list_i=0; buf_list_i < num_buffer_lists; buf_list_i++ ) {

        Buffer *buf_list = buffer_lists[buf_list_i]; // a pointer to a buffer, which contains an array of pointers to buffers of a single vb/non-vb

        if (!buf_list->len) continue; // no buffers allocated yet for this VB

        for (unsigned buf_i=0; buf_i < buf_list->len; buf_i++) {
    
            ASSERT (buf_list->memory, "Error: buf_list[%u] memory is not allocated", buf_list_i); // this should never happen

            Buffer *buf = ((Buffer **)buf_list->data)[buf_i];
            
            if (!buf->memory) continue; // exclude destroyed or not-yet-allocated buffers

            bool found = false;
            for (unsigned st_i=0; st_i < num_stats && !found; st_i++) {
                MemStats *st = &stats[st_i];

                if (!strcmp (st->name, buf->name)) {
                    st->buffers++;
                    st->bytes += buf->size + 2*(sizeof(long long));
                    found = true;
                }
            }

            if (!found) {
                stats[num_stats].name    = buf->name;
                stats[num_stats].bytes   = buf->size + 2*(sizeof(long long));
                stats[num_stats].buffers = 1;
                num_stats++;
                ASSERT (num_stats < MAX_MEMORY_STATS, "# memory stats exceeded %u, consider increasing MAX_MEMORY_STATS", MAX_MEMORY_STATS);
            }

            num_buffers++;
        }
    }

    // sort stats by bytes
    qsort (stats, num_stats, sizeof (MemStats), buf_stats_sort_by_bytes);

    unsigned total_bytes=0;
    for (unsigned i=0; i< num_stats; i++) total_bytes += stats[i].bytes;

    char str[30];
    buf_human_readable_size (total_bytes, str);
    printf ("Total bytes: %s in %u buffers in %u buffer lists:\n", str, num_buffers, num_buffer_lists);

    for (unsigned i=0; i < num_stats; i++) {
        buf_human_readable_size (stats[i].bytes, str);
        printf ("%-30s: %-8s (%4.1f%%) in %u buffers\n", stats[i].name, str, 100.0 * (float)stats[i].bytes / (float)total_bytes, stats[i].buffers);
    }
}

long long buf_vb_memory_consumption (const VariantBlock *vb)
{
    const Buffer *buf_list   = &vb->buffer_list;

    // memory of the structure itself
    long long vb_memory = sizeof (*vb);

    // memory allocated outside of Buffer (direct calloc)
    if (vb->haplotype_sections_data) vb_memory += vb->num_sample_blocks * sizeof (Buffer);
    if (vb->genotype_sections_data)  vb_memory += vb->num_sample_blocks * sizeof (Buffer);
    if (vb->phase_sections_data)     vb_memory += vb->num_sample_blocks * sizeof (Buffer);

    // memory allocated via Buffer (the vast majority, usually)
    for (unsigned buf_i=0; buf_i < buf_list->len; buf_i++) 
        vb_memory += ((Buffer **)buf_list->data)[buf_i]->size;

    return vb_memory;
}

static inline void buf_add(VariantBlock *vb, Buffer *buf)
{
    Buffer *buffer_list   = &vb->buffer_list;
    unsigned max_num_buffers = buffer_list->size / sizeof (Buffer *);

#define INITIAL_MAX_MEM_NUM_BUFFERS (VARIANTS_PER_BLOCK*3) /* for files that have ht,gt,phase,variant,and line - the factor would be about 5.5 so there will be 1 realloc per vb, but most files don't */
    if (!buf_is_allocated (buffer_list)) {

        if (!buffer_lists) buffer_lists = malloc (sizeof (Buffer *) * (global_max_threads + 1)); // +1 for psuedo vbs

        ASSERT (num_buffer_lists < MAX(2,global_max_threads) + 1, "Error: buffer_lists maxed out num_buffer_lists=%u", num_buffer_lists);

        // malloc - this will call this function recursively - that's ok bc that point buffer_list is already allocated
        buf_alloc (vb, buffer_list, INITIAL_MAX_MEM_NUM_BUFFERS * sizeof(Buffer *), false, "buffer_list", vb ? vb->id : 0);
        
        buffer_lists[num_buffer_lists++] = buffer_list;
    }
    else if (buffer_list->len && buffer_list->len == max_num_buffers) 
        buf_alloc (vb, buffer_list, buffer_list->size * 2, 1, "buffer_list realloc", vb ? vb->id : 0);

    ((Buffer **)buffer_list->data)[buffer_list->len++] = buf;
}

// allocates or enlarges buffer. if this buffer is not already allocated, then it behaves like malloc.
unsigned buf_alloc (VariantBlock *vb,
                    Buffer *buf, 
                    unsigned requested_size,
                    float grow_at_least_factor, // IF we need to allocate or reallocate physical memory, we get this much more than requested
                    const char *name, unsigned param)      // for debugging
{
    START_TIMER;

    if (!requested_size) return 0; // nothing to do

    // sanity checks
    ASSERT (buf->type != BUF_OVERLAYED, "Error: cannot buf_alloc an overlayed buffer. name=%s", buf->name ? buf->name : "");
    ASSERT0 (vb, "Error: null vb");

    // case 1: we have enough memory already
    if (requested_size <= buf->size) 
        buf->data = buf->memory + sizeof (long long); // allocate if not already allocated

    else { // not enough memory
        unsigned new_size = requested_size * MAX (grow_at_least_factor, 1);

        // case 2: we need to allocate memory - buffer is already allocated so copy over the data
        if (buf->memory) {
#ifdef DISPLAY_ALLOCS_AFTER
            unsigned  old_size = buf->size;
#endif
            buf->memory = realloc (buf->memory, new_size + 2*sizeof (long long));
            ASSERT (buf->memory, "Error: buf_alloc failed to realloc %u bytes. name=%s param=%u", new_size + 2*(unsigned)sizeof (long long), name, param);

            buf->data = buf->memory + sizeof (long long);
            buf->size = new_size;

            if (!buf->memory) {
                buf_test_overflows(vb);
#ifdef DEBUG
                    buf_display_memory_usage(true);
#endif
            }

            *(long long *)(buf->data + new_size) = OVERFLOW_TRAP; // overflow prortection (underflow protection was copied with realloc)

#ifdef DISPLAY_ALLOCS_AFTER
                if (vb->buffer_list.len > DISPLAY_ALLOCS_AFTER)
                    printf ("%s (%u): old_size=%u requested_size=%u growth_factor=%f new_size=%u\n", name, param, old_size, requested_size, grow_at_least_factor, new_size);
#endif
        }

        // case 3: we need to allocate memory - buffer is not yet allocated, so no need to copy data
        else {
            buf->memory = malloc (new_size + 2*sizeof (long long));

            if (!buf->memory) {
                buf_test_overflows(vb);
    
#               ifdef DEBUG
                    buf_display_memory_usage(true);
#               endif
            }

            ASSERT (buf->memory, "Error: buf_alloc failed to malloc %u bytes. name=%s param=%u", new_size + 2*(unsigned)sizeof (long long), name, param);

            buf->data = buf->memory + sizeof (long long);
            buf->size = new_size;

            *(long long *)buf->memory            = UNDERFLOW_TRAP; // underflow protection
            *(long long *)(buf->data + new_size) = OVERFLOW_TRAP;  // overflow protection

            buf->name  = name;
            buf->param = param;
            buf->type  = BUF_REGULAR;
            
            buf_add(vb, buf);

#ifdef DISPLAY_ALLOCS_AFTER
            if (vb->buffer_list.len > DISPLAY_ALLOCS_AFTER)
                printf ("%s (%u) (malloc): %u\n", name, param, new_size);
#endif
        }
    }

    if (vb) COPY_TIMER (vb->profile.buf_alloc);

    return buf->size;
}

// an overlay buffer is a buffer using some of the memory of another buffer - it doesn't have its own memory
void buf_overlay (Buffer *overlaid_buf, Buffer *regular_buf, const Buffer *copy_from /* optional */, 
                  unsigned *regular_buf_offset, const char *name, unsigned param)
{
    ASSERT (overlaid_buf->type == BUF_UNALLOCATED, "Error: cannot buf_overlay to a buffer already in use. overlaid_buf->name=%s", overlaid_buf->name ? overlaid_buf->name : "");
    ASSERT (regular_buf->type == BUF_REGULAR, "Error: regular_buf in buf_overlay must be a regular buffer. regular_buf->name=%s", regular_buf->name ? regular_buf->name : "");
    ASSERT (!copy_from || (regular_buf_offset ? *regular_buf_offset : 0) + copy_from->len <= regular_buf->size, "Error: buf_overlay exceeds the size of the regular buffer: offset=%u size=%u regular_buf.size=%u", *regular_buf_offset, copy_from->len, regular_buf->size);

    overlaid_buf->data   = regular_buf->data + (regular_buf_offset ? *regular_buf_offset : 0);

    if (copy_from && copy_from->len)
        memcpy (overlaid_buf->data, copy_from->data, copy_from->len);

    overlaid_buf->size   = 0;
    overlaid_buf->len    = copy_from ? copy_from->len : 0;
    overlaid_buf->type   = BUF_OVERLAYED;
    overlaid_buf->memory = 0;
    overlaid_buf->name   = name;
    overlaid_buf->param  = param;

    // new data was copied to overlaid buffer - move the offset forward and update the len of the regular buffer
    // to enable to the next buffer to be overlaid subsequently
    if (regular_buf_offset) {
        *regular_buf_offset += overlaid_buf->len;
        regular_buf->len = *regular_buf_offset;
    }
    
    // full buffer overlay - copy len too
    if (!regular_buf_offset && !copy_from)
        overlaid_buf->len = regular_buf->len;
}

static Buffer abandoned_memories = EMPTY_BUFFER; // no thread safety issues - only used by dispatcher thread

// if needed, we extend the buf, but without reallocting memory which could disturb other threads which have their
// own buffers overlaid on buf. Instead, we abandon the memory used by buf, but without freeing it.
// We add the abandoned memory to an abandoned memory list to be freed at the end of this operation (piz). Then,
// we allocate new memory (at least double) and copy the old data and len to it.
void buf_extend (VariantBlock *vb, Buffer *buf, 
                 unsigned num_new_units, unsigned sizeof_unit, /* len is measured in these units */
                 const char *name, unsigned param)
{
    if (!buf_is_allocated (buf)) 
        buf_alloc (vb, buf, num_new_units * sizeof_unit, 1, name, param);

    // if we need more space, we abandon the memory used by buf, but without freeing it as VBs are overlaying on it.
    // instead, we add it to an abandoned memory list to be freed at the end of this piz
    if ((buf->len + num_new_units) * sizeof_unit > buf->size) { 

        if (!buf_is_allocated (&abandoned_memories) || abandoned_memories.len * sizeof(char *) == abandoned_memories.size)
            buf_alloc (vb, &abandoned_memories, 1000 * sizeof(char*), 2, "abandoned_memories", 0);
        
        unsigned old_len       = buf->len;
        const char *old_memory = buf->memory;
        const char *old_data   = buf->data;

        ((const char**)abandoned_memories.data)[abandoned_memories.len++] = old_memory;

        memset (buf, 0, sizeof (Buffer));
        buf_alloc (vb, buf, sizeof_unit * MAX (old_len * 2, old_len + num_new_units), 1, name, param);
        memcpy (buf->data, old_data, old_len * sizeof_unit);
        buf->len = old_len;
    }
}

// if there is space, the fragment is just appended to the base and len is updated. 
// If not, we extend the memory without disturbing other VBs that are overlaying on it, and add
// the fragment
void buf_append (VariantBlock *vb, Buffer *base, 
                 Buffer *fragment, unsigned sizeof_unit, // len is measured in these units 
                 const char *name, unsigned param)
{
    buf_extend (vb, base, fragment->len, sizeof_unit, name, param);

    memcpy (&base->data[base->len * sizeof_unit], fragment->data, fragment->len * sizeof_unit);
    base->len += fragment->len;
}

void buf_free_abandoned_memories()
{
    for (unsigned i=0; i < abandoned_memories.len; i++)
        free (((char**)abandoned_memories.data)[i]);
}

// free buffer - without freeing memory. A future buf_malloc of this buffer will reuse the memory if possible.
void buf_free(Buffer *buf) 
{
    if (buf->type == BUF_REGULAR) {
        buf->data = NULL ; 
        buf->len = 0;
        // name, memory and size are not changed
    }
    else
        memset (buf, 0, sizeof (Buffer));
} 

void buf_copy (VariantBlock *vb, Buffer *dst, Buffer *src, 
               unsigned bytes_per_entry, // how many bytes are counted by a unit of .len
               unsigned start_entry, unsigned max_entries,  // if 0 copies the entire buffer 
               const char *name, unsigned param)
{
    ASSERT0 (src->data, "Error in buf_copy: src->data is NULL");
    
    unsigned num_entries = max_entries ? MIN (max_entries, src->len - start_entry) : src->len - start_entry;

    buf_alloc(vb, dst, num_entries * bytes_per_entry, 1, name, param); // use realloc rather than malloc to allocate exact size

    memcpy (dst->data, &src->data[start_entry * bytes_per_entry], num_entries * bytes_per_entry);

    dst->len = num_entries;  
}   
