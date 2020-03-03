// ------------------------------------------------------------------
//   buffer.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef BUFFER_INCLUDED
#define BUFFER_INCLUDED

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "genozip.h"

struct variant_block_; 

typedef enum {BUF_UNALLOCATED=0, BUF_REGULAR, BUF_FULL_OVERLAY, BUF_PARTIAL_OVERLAY} BufferType; // BUF_UNALLOCATED must be 0

typedef struct buffer_ {
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

#define buf_is_allocated(buf_p) ((buf_p)->data != NULL && ((buf_p)->type == BUF_REGULAR || (buf_p)->type == BUF_FULL_OVERLAY || (buf_p)->type == BUF_PARTIAL_OVERLAY))

extern unsigned buf_alloc_do (VariantBlockP vb,
                              Buffer *buf, 
                              unsigned requested_size, 
                              float grow_at_least_factor, // grow more than new_size    
                              const char *name, unsigned param); // for debugging

// efficient wrapper
#define buf_alloc(vb, buf, requested_size, grow_at_least_factor, name, param) \
  ((!(buf)->data || (buf)->size < (requested_size)) ? buf_alloc_do ((vb), (buf), (requested_size), (grow_at_least_factor), (name), (param)) \
                                                    : (buf)->size) 

#define buf_set_overlayable(buf) (buf)->overlayable = true

extern void buf_overlay (Buffer *overlaid_buf, Buffer *regular_buf, const Buffer *copy_from, unsigned *regular_buf_offset, const char *name, unsigned param);
extern void buf_free (Buffer *buf); // free buffer - without freeing memory. A future buf_alloc of this buffer will reuse the memory if possible.
extern void buf_destroy (VariantBlockP vb, Buffer *buf);

#define buf_is_large_enough(buf_p, requested_size) (buf_is_allocated ((buf_p)) && (buf_p)->size >= requested_size)

extern void buf_copy (VariantBlockP vb, Buffer *dst, const Buffer *src, unsigned bytes_per_entry,
                      unsigned start_entry, unsigned max_entries, // if 0 copies the entire buffer
                      const char *name, unsigned param);

extern void buf_move (VariantBlockP vb, Buffer *dst, Buffer *src);

#define buf_add(buf, new_data, new_data_len) { memcpy (&(buf)->data[(buf)->len], (new_data), (new_data_len));  (buf)->len += (new_data_len); }
extern void buf_add_string (VariantBlockP vb, Buffer *buf, const char *str);
#define bufprintf(vb, buf, format, ...)  { char s[5000]; sprintf (s, (format), __VA_ARGS__); buf_add_string ((vb), (buf), s); }

extern void buf_test_overflows(ConstVariantBlockP vb);

extern int64_t buf_vb_memory_consumption (ConstVariantBlockP vb);
extern void buf_display_memory_usage (bool memory_full);

extern char *buf_human_readable_size (int64_t size, char *str /* out */);
extern char *buf_human_readable_uint (int64_t n, char *str /* out */);

#define buf_zero(buf_p) { memset ((buf_p)->data, 0, (buf_p)->size); }

#endif