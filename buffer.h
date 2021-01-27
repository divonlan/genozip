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
#include "endianness.h"

struct variant_block_; 

typedef enum { BUF_UNALLOCATED=0, BUF_REGULAR, BUF_OVERLAY, BUF_MMAP, BUF_NUM_TYPES } BufferType; // BUF_UNALLOCATED must be 0, must be identical to BitArrayType
#define BUFTYPE_NAMES { "UNALLOCATED", "REGULAR", "OVERLAY", "MMAP" }

typedef struct Buffer {
    bool overlayable; // this buffer may be fully overlaid by one or more overlay buffers
    
    const char *name; // name of allocator - used for memory debugging & statistics
    uint64_t size;    // number of bytes available to the user (i.e. not including the allocated overhead)

    //------------------------------------------------------------------------------------------------------
    // these 4 fields can be overlayed with a BitArray - their order and size is identical to BitArray
    BufferType type;
    char *data;       // ==memory+8 if buffer is allocated or NULL if not
    int64_t param;    // parameter provided by allocator 
    uint64_t len;     // used by the buffer user according to its internal logic. not modified by malloc/realloc, zeroed by buf_free
    //------------------------------------------------------------------------------------------------------

    char *memory;     // memory allocated to this buffer - amount is: size + 2*sizeof(longlong) to allow for OVERFLOW and UNDERFLOW)

    // info on the allocator of this buffer
    VBlockP vb;       // vb that owns this buffer, and which this buffer is in its buf_list
    const char *func; // the allocating function
    uint32_t code_line;
} Buffer;

#define EMPTY_BUFFER {}

#define ARRAY(element_type, name, buf) element_type *name = ((element_type *)((buf).data)) 

#define ENT(type, buf, index) ((type *)(&(buf).data[(index) * sizeof(type)]))
#define FIRSTENT(type, buf)   ((type *)( (buf).data))
#define LASTENT(type, buf)    ((type *)(&(buf).data[((buf).len-1) * sizeof(type)]))
#define AFTERENT(type, buf)   ((type *)(&(buf).data[((buf).len  ) * sizeof(type)]))

typedef struct { char s[300]; } BufDescType;
extern const BufDescType buf_desc (const Buffer *buf);

static inline uint64_t NEXTENT_get_index (Buffer *buf, size_t size, const char *func, uint32_t code_line) 
{ 
    uint64_t index = (buf->len++) * size;
    ASSERTE (index + size <= buf->size, "called from %s:%u: NEXTENT went beyond end of buffer: size=%u index=%"PRIu64": %s", 
             func, code_line, (unsigned)size, index, buf_desc (buf).s);
    return index;
}
#define NEXTENT(type, buf)    (*(type *)(&(buf).data[NEXTENT_get_index (&(buf), sizeof(type), __FUNCTION__, __LINE__)]))
#define ENTNUM(buf, ent)      ((uint32_t)((((char*)(ent)) - ((buf).data)) / sizeof (*ent)))
#define ISLASTENT(buf,ent)    (ENTNUM((buf),(ent)) == (buf).len - 1)

extern void buf_initialize(void);

#define buf_is_allocated(buf_p) ((buf_p)->data != NULL && (buf_p)->type != BUF_UNALLOCATED)

extern uint64_t buf_alloc_do (VBlockP vb,
                              Buffer *buf, 
                              uint64_t requested_size, 
                              double grow_at_least_factor, // grow more than new_size   
                              const char *func, uint32_t code_line,
                              const char *name);

// efficient wrapper
#define buf_alloc(vb, buf, requested_size, grow_at_least_factor, name) { \
    uint64_t new_req_size = (requested_size); /* make copy to allow ++ */ \
    ((!(buf)->data || (buf)->size < (new_req_size)) ? buf_alloc_do ((VBlockP)(vb), (buf), (new_req_size), (grow_at_least_factor), __FUNCTION__, __LINE__, (name)) \
                                                    : (buf)->size); \
}

#define buf_alloc_more(vb, buf, more, at_least, type, grow_at_least_factor,name) \
    buf_alloc ((vb), (buf), MAX((at_least), ((buf)->len+(more)))*sizeof(type), (grow_at_least_factor), (name))

#define buf_alloc_more_zero(vb, buf, more, at_least, type, grow_at_least_factor,name) { \
    uint64_t size_before = (buf)->size; \
    buf_alloc_more((vb), (buf), (more), (at_least), type, (grow_at_least_factor), (name)); \
    if ((buf)->size > size_before) memset (&(buf)->data[size_before], 0, (buf)->size - size_before); \
}

extern bool buf_mmap_do (VBlockP vb, Buffer *buf, const char *filename, const char *func, uint32_t code_line, const char *name);
#define buf_mmap(vb, buf, filename, name) \
    buf_mmap_do((vb), (buf), (filename), __FUNCTION__, __LINE__, (name))

extern BitArrayP buf_alloc_bitarr_do (VBlockP vb, Buffer *buf, uint64_t nbits, const char *func, uint32_t code_line, const char *name);
#define buf_alloc_bitarr(vb, buf, nbits, name) buf_alloc_bitarr_do((vb), (buf), (nbits), __FUNCTION__, __LINE__, (name))

extern BitArrayP buf_overlay_bitarr_do (VBlockP vb, Buffer *overlaid_buf, Buffer *regular_buf, uint64_t start_byte_in_regular_buf, uint64_t nbits, const char *func, uint32_t code_line, const char *name);
#define buf_overlay_bitarr(vb, overlaid_buf, regular_buf, start_byte_in_regular_buf, nbits, name) \
    buf_overlay_bitarr_do((vb), (overlaid_buf), (regular_buf), (start_byte_in_regular_buf), (nbits), __FUNCTION__, __LINE__, (name))

#define buf_set_overlayable(buf) (buf)->overlayable = true

extern void buf_overlay_do (VBlockP vb, Buffer *overlaid_buf, Buffer *regular_buf,  uint64_t start_in_regular,  
                            const char *func, uint32_t code_line, const char *name);
#define buf_overlay(vb, overlaid_buf, regular_buf, name) \
    buf_overlay_do((vb), (overlaid_buf), (regular_buf), 0, __FUNCTION__, __LINE__, (name)) 

#define buf_overlay_partial(vb, overlaid_buf, regular_buf, start, name) \
    buf_overlay_do((vb), (overlaid_buf), (regular_buf), (start), __FUNCTION__, __LINE__, (name)) 

extern void buf_free_do (Buffer *buf, const char *func, uint32_t code_line);
#define buf_free(buf) buf_free_do ((buf), __FUNCTION__, __LINE__);

extern void buf_destroy_do (Buffer *buf, const char *func, uint32_t code_line);
#define buf_destroy(buf) buf_destroy_do ((buf), __FUNCTION__, __LINE__)

#define buf_is_large_enough(buf_p, requested_size) (buf_is_allocated ((buf_p)) && (buf_p)->size >= requested_size)

extern void buf_copy_do (VBlockP dst_vb, Buffer *dst, const Buffer *src, uint64_t bytes_per_entry,
                         uint64_t src_start_entry, uint64_t max_entries, // if 0 copies the entire buffer
                         const char *func, uint32_t code_line,
                         const char *name);
#define buf_copy(dst_vb,dst,src,bytes_per_entry,src_start_entry,max_entries,name) \
  buf_copy_do ((VBlockP)(dst_vb),(dst),(src),(bytes_per_entry),(src_start_entry),(max_entries),__FUNCTION__,__LINE__,(name))

extern void buf_move (VBlockP dst_vb, Buffer *dst, VBlockP src_vb, Buffer *src);

#define buf_has_space(buf, new_len) ((buf)->len + (new_len) <= (buf)->size)

#define buf_add(buf, new_data, new_data_len) { uint32_t new_len = (new_data_len); /* copy in case caller uses ++ */ \
                                               ASSERTE (buf_has_space(buf, new_len), \
                                                        "Error in buf_add: buffer %s is out of space: len=%u size=%u new_data_len=%u", \
                                                        buf_desc (buf).s, (uint32_t)(buf)->len, (uint32_t)(buf)->size, new_len);\
                                               memcpy (&(buf)->data[(buf)->len], (new_data), new_len);   \
                                               (buf)->len += new_len; }

extern void buf_add_string (VBlockP vb, Buffer *buf, const char *str);
extern void buf_add_int (VBlockP vb, Buffer *buf, int64_t value);

#define BUFPRINTF_MAX_LEN 5000
#define bufprintf(vb, buf, format, ...)  { char __s[BUFPRINTF_MAX_LEN]; sprintf (__s, (format), __VA_ARGS__); buf_add_string ((vb), (buf), __s); }

extern void buf_print (Buffer *buf, bool add_newline);

extern void buf_test_overflows (void *vb, const char *msg);
extern void buf_test_overflows_all_vbs (const char *msg);

typedef struct {
    const char *name;
    uint64_t bytes; 
    unsigned buffers;
} MemStats;

extern void buf_display_memory_usage (bool memory_full, unsigned max_threads, unsigned used_threads);

#define buf_set(buf_p,value) { if ((buf_p)->data) memset ((buf_p)->data, value, (buf_p)->size); }
#define buf_zero(buf_p) buf_set(buf_p, 0)

extern void buf_add_to_buffer_list (VBlockP vb, Buffer *buf);
extern void buf_remove_from_buffer_list (Buffer *buf);

extern void buf_low_level_free (void *p, const char *func, uint32_t code_line);
#define FREE(p) { if (p) { buf_low_level_free (((void*)p), __FUNCTION__, __LINE__); p=NULL; } }

extern void *buf_low_level_realloc (void *p, size_t size, const char *func, uint32_t code_line);
#define REALLOC(p,size) buf_low_level_realloc (p, size, __FUNCTION__, __LINE__)

extern void *buf_low_level_malloc (size_t size, bool zero, const char *func, uint32_t code_line);
#define MALLOC(size) buf_low_level_malloc (size, false, __FUNCTION__, __LINE__)
#define CALLOC(size) buf_low_level_malloc (size, true,  __FUNCTION__, __LINE__)

extern bool buf_dump_to_file (const char *filename, const Buffer *buf, unsigned buf_word_width, bool including_control_region, bool no_dirs);

// bitmap stuff
extern uint64_t buf_extend_bits (Buffer *buf, int64_t num_new_bits);
extern void buf_add_bit (Buffer *buf, int64_t new_bit);
extern BitArrayP buf_zfile_buf_to_bitarray (Buffer *buf, uint64_t nbits);
#define buf_get_bitarray(buf) ((BitArrayP)(&(buf)->type))
#define buf_add_set_bit(buf)   buf_add_bit (buf, 1)
#define buf_add_clear_bit(buf) buf_add_bit (buf, 0)

// endianness
extern BgEnBufFunc BGEN_u8_buf, BGEN_u16_buf, BGEN_u32_buf, BGEN_u64_buf, 
                   BGEN_transpose_u8_buf, BGEN_transpose_u16_buf, BGEN_transpose_u32_buf, BGEN_transpose_u64_buf,
                   BGEN_deinterlace_d8_buf, BGEN_deinterlace_d16_buf, BGEN_deinterlace_d32_buf, BGEN_deinterlace_d64_buf;

#endif