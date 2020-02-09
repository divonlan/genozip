// ------------------------------------------------------------------
//   buffer.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef BUFFER_INCLUDED
#define BUFFER_INCLUDED

#include "genozip.h"

struct variant_block_; 

#define NUM_POOLS 2 
typedef enum { POOL_ID_UNIT_TEST=-1, POOL_ID_ZIP=0, POOL_ID_UNZIP=1 } PoolId;

typedef enum {BUF_UNALLOCATED=0, BUF_REGULAR, BUF_FULL_OVERLAY, BUF_PARTIAL_OVERLAY} BufferType; // BUF_UNALLOCATED must be 0

typedef struct {
    BufferType type;
    bool overlayable; // this buffer may be fully overlaid by one or more overlay buffers
    const char *name; // name of allocator - used for memory debugging & statistics
    unsigned param;   // parameter provided by allocator - used for memory debugging & statistics
    unsigned size;    // number of bytes allocated to memory
    unsigned len;     // used by the buffer user according to its internal logic. not modified by malloc/realloc, zeroed by buf_free
    char *data;       // =memory+2*sizeof(long long) if buffer is allocated or NULL if not
    char *memory;     // memory allocated to this buffer - amount is: size + 2*sizeof(longlong) to allow for OVERFLOW and UNDERFLOW)
} Buffer;
#define EMPTY_BUFFER {BUF_UNALLOCATED,false,NULL,0,0,0,NULL,NULL}

extern void buf_initialize();
extern unsigned buf_alloc (VariantBlockP vb,
                           Buffer *buf, 
                           unsigned requested_size, // whether contents of memory should be zeroed
                           float grow_at_least_factor, // grow more than new_size    
                           const char *name, unsigned param); // for debugging
static inline void buf_set_overlayable (Buffer *buf) { buf->overlayable = true;}
extern void buf_overlay (Buffer *overlaid_buf, Buffer *regular_buf, const Buffer *copy_from, unsigned *regular_buf_offset, const char *name, unsigned param);
extern void buf_free (Buffer *buf); // free buffer - without freeing memory. A future buf_alloc of this buffer will reuse the memory if possible.
extern void buf_destroy (VariantBlockP vb, Buffer *buf);

static inline bool buf_is_allocated (const Buffer *buf) {return buf->data != NULL && (buf->type == BUF_REGULAR || buf->type == BUF_FULL_OVERLAY || buf->type == BUF_PARTIAL_OVERLAY);}

extern void buf_copy (VariantBlockP vb, Buffer *dst, const Buffer *src, unsigned bytes_per_entry,
                      unsigned start_entry, unsigned max_entries, // if 0 copies the entire buffer
                      const char *name, unsigned param);

extern void buf_move (VariantBlockP vb, Buffer *dst, Buffer *src);

static inline void buf_add (Buffer *buf, const void *data, unsigned len) { memcpy (&buf->data[buf->len], data, len);  buf->len += len; }
#define buf_add_string(buf,str) buf_add (buf, str, strlen (str));

extern void buf_test_overflows(ConstVariantBlockP vb);

extern long long buf_vb_memory_consumption (ConstVariantBlockP vb);
extern void buf_display_memory_usage (PoolId pool_id, bool memory_full);

extern char *buf_human_readable_size (uint64_t size, char *str /* out */);

#endif