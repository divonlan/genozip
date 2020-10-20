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

typedef enum {BUF_UNALLOCATED=0, BUF_REGULAR, BUF_OVERLAY} BufferType; // BUF_UNALLOCATED must be 0
#define BUFTYPE_NAMES { "UNALLOCATED", "REGULAR", "OVERLAY" }

typedef struct Buffer {
    BufferType type;
    bool overlayable; // this buffer may be fully overlaid by one or more overlay buffers
    
    const char *name; // name of allocator - used for memory debugging & statistics
    uint64_t size;    // number of bytes available to the user (i.e. not including the allocated overhead)

    //------------------------------------------------------------------------------------------------------
    // these 3 fields can be overlayed with a BitArray - their order and size is identical to BitArray
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

#define EMPTY_BUFFER { .type        = BUF_UNALLOCATED,\
                       .overlayable = false,\
                       .name        = NULL,\
                       .param       = 0,\
                       .size        = 0,\
                       .len         = 0,\
                       .data        = NULL,\
                       .memory      = NULL,\
                       .vb          = NULL,\
                       .func        = NULL,\
                       .code_line   =0 }

#define ARRAY(element_type, name, buf) element_type *name = ((element_type *)((buf).data)) 

#define ENT(type, buf, index) ((type *)(&(buf).data[(index) * sizeof(type)]))
#define FIRSTENT(type, buf)   ((type *)( (buf).data))
#define LASTENT(type, buf)    ((type *)(&(buf).data[((buf).len-1) * sizeof(type)]))
#define AFTERENT(type, buf)   ((type *)(&(buf).data[((buf).len  ) * sizeof(type)]))
#define NEXTENT(type, buf)    (*(type *)(&(buf).data[((buf).len++) * sizeof(type)]))
#define ENTNUM(buf, ent)      ((((char*)(ent)) - ((buf).data)) / sizeof (*ent))
#define ISLASTENT(buf,ent)    (ENTNUM((buf),(ent)) == (buf).len - 1)

extern void buf_initialize(void);

#define buf_is_allocated(buf_p) ((buf_p)->data != NULL && ((buf_p)->type == BUF_REGULAR || (buf_p)->type == BUF_OVERLAY))

extern uint64_t buf_alloc_do (VBlockP vb,
                              Buffer *buf, 
                              uint64_t requested_size, 
                              double grow_at_least_factor, // grow more than new_size   
                              const char *func, uint32_t code_line,
                              const char *name, int64_t param);

// efficient wrapper
#define buf_alloc(vb, buf, requested_size, grow_at_least_factor, name, param) { \
    uint64_t new_req_size = (requested_size); /* make copy to allow ++ */ \
    ((!(buf)->data || (buf)->size < (new_req_size)) ? buf_alloc_do ((VBlockP)(vb), (buf), (new_req_size), (grow_at_least_factor), __FUNCTION__, __LINE__, (name), (param)) \
                                                    : (buf)->size); \
}

#define buf_alloc_more(vb, buf, more, at_least, type, grow_at_least_factor) \
    buf_alloc ((vb), (buf), MAX(at_least, ((buf)->len+(more)))*sizeof(type), (grow_at_least_factor), (buf)->name, (buf)->param)

#define buf_alloc_more_zero(vb, buf, more, at_least, type, grow_at_least_factor) { \
    uint64_t size_before = (buf)->size; \
    buf_alloc_more((vb), (buf), (more), (at_least), type, (grow_at_least_factor)); \
    if ((buf)->size > size_before) memset (&(buf)->data[size_before], 0, (buf)->size - size_before); \
}

#define buf_set_overlayable(buf) (buf)->overlayable = true

extern void buf_overlay_do (VBlockP vb, Buffer *overlaid_buf, Buffer *regular_buf, const char *func, uint32_t code_line, const char *name, int64_t param);
#define buf_overlay(vb, overlaid_buf, regular_buf, name, param) \
    buf_overlay_do(vb, overlaid_buf, regular_buf, __FUNCTION__, __LINE__, name, param) 

extern void buf_free_do (Buffer *buf, const char *func, uint32_t code_line);
#define buf_free(buf) buf_free_do (buf, __FUNCTION__, __LINE__);

extern void buf_destroy_do (Buffer *buf, const char *func, uint32_t code_line);
#define buf_destroy(buf) buf_destroy_do (buf, __FUNCTION__, __LINE__)

#define buf_is_large_enough(buf_p, requested_size) (buf_is_allocated ((buf_p)) && (buf_p)->size >= requested_size)

extern void buf_copy_do (VBlockP dst_vb, Buffer *dst, const Buffer *src, uint64_t bytes_per_entry,
                         uint64_t src_start_entry, uint64_t max_entries, // if 0 copies the entire buffer
                         const char *func, uint32_t code_line,
                         const char *name, int64_t param);
#define buf_copy(dst_vb,dst,src,bytes_per_entry,src_start_entry,max_entries,name,param) \
  buf_copy_do ((VBlockP)(dst_vb),(dst),(src),(bytes_per_entry),(src_start_entry),(max_entries),__FUNCTION__,__LINE__,(name),(param))

extern void buf_move (VBlockP dst_vb, Buffer *dst, VBlockP src_vb, Buffer *src);

#define buf_has_space(buf, new_len) ((buf)->len + (new_len) <= (buf)->size)

#define buf_add(buf, new_data, new_data_len) { uint32_t new_len = (new_data_len); /* copy in case caller uses ++ */ \
                                               ASSERT (buf_has_space(buf, new_len), \
                                                       "Error in buf_add called from %s:%u: buffer %s is out of space: len=%u size=%u new_data_len=%u", \
                                                       __FUNCTION__, __LINE__, buf_desc (buf), (uint32_t)(buf)->len, (uint32_t)(buf)->size, new_len);\
                                               memcpy (&(buf)->data[(buf)->len], (new_data), new_len);   \
                                               (buf)->len += new_len; }

extern void buf_add_string (VBlockP vb, Buffer *buf, const char *str);
#define bufprintf(vb, buf, format, ...)  { char __s[5000]; sprintf (__s, (format), __VA_ARGS__); buf_add_string ((vb), (buf), __s); }

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

extern void *buf_low_level_malloc (size_t size, const char *func, uint32_t code_line);
#define MALLOC(size) buf_low_level_malloc (size, __FUNCTION__, __LINE__)

extern const char *buf_desc (const Buffer *buf);

// bitmap stuff
extern uint64_t buf_extend_bits (Buffer *buf, int64_t num_new_bits);
extern void buf_add_bit (Buffer *buf, int64_t new_bit);
extern BitArrayP buf_zfile_buf_to_bitarray (Buffer *buf, uint64_t num_of_bits);
#define buf_get_bitarray(buf) ((BitArray*)(&(buf)->data))
#define buf_add_set_bit(buf)   buf_add_bit (buf, 1)
#define buf_add_clear_bit(buf) buf_add_bit (buf, 0)

// endianness
extern BgEnBufFunc BGEN_u8_buf, BGEN_u16_buf, BGEN_u32_buf, BGEN_u64_buf, 
                   BGEN_deinterlace_d8_buf, BGEN_deinterlace_d16_buf, BGEN_deinterlace_d32_buf, BGEN_deinterlace_d64_buf;

#endif