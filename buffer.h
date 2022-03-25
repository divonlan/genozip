// ------------------------------------------------------------------
//   buffer.h
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "genozip.h"

struct variant_block_; 

typedef enum { BUF_UNALLOCATED=0, BUF_REGULAR, BUF_OVERLAY, BUF_MMAP, BUF_NUM_TYPES } BufferType; // BUF_UNALLOCATED must be 0, must be identical to BitArrayType
#define BUFTYPE_NAMES { "UNALLOCATED", "REGULAR", "OVERLAY", "MMAP" }

typedef struct Buffer {
    bool overlayable; // this buffer may be fully overlaid by one or more overlay buffers
    bool can_be_big;  // do not display warning if buffer grows very big
    
    const char *name; // name of allocator - used for memory debugging & statistics
    uint64_t size;    // number of bytes available to the user (i.e. not including the allocated overhead)

    //------------------------------------------------------------------------------------------------------
    // these 4 fields can be overlayed with a BitArray - their order and size is identical to BitArray
    BufferType type;
    char *data;       // ==memory+8 if buffer is allocated or NULL if not
    int64_t param;    // parameter provided by allocator (in BitArray - nbits)
    uint64_t len;     // used by the buffer user according to its internal logic. not modified by malloc/realloc, zeroed by buf_free (in BitArray - nwords)
    //------------------------------------------------------------------------------------------------------

    char *memory;     // memory allocated to this buffer - amount is: size + 2*sizeof(longlong) to allow for OVERFLOW and UNDERFLOW)

    // info on the allocator of this buffer
    VBlockP vb;       // vb that owns this buffer, and which this buffer is in its buf_list
    const char *func; // the allocating function
    uint32_t code_line;
} Buffer;

#define EMPTY_BUFFER {}

#define ARRAY(element_type, name, buf) \
    element_type *name = ((element_type *)((buf).data)); \
    const uint64_t name##_len __attribute__((unused)) = (buf).len; // read-only copy of len 
    
#define ENT(type, buf, index) ((type *)(&(buf).data[(index) * sizeof(type)]))
#define FIRSTENT(type, buf)   ((type *)( (buf).data))
#define LASTENT(type, buf)    ((type *)(&(buf).data[((buf).len-1) * sizeof(type)]))
#define AFTERENT(type, buf)   ((type *)(&(buf).data[((buf).len  ) * sizeof(type)]))

typedef struct { char s[300]; } BufDescType;
extern const BufDescType buf_desc (const Buffer *buf);

static inline uint64_t NEXTENT_get_index (Buffer *buf, size_t size, const char *func, uint32_t code_line) 
{ 
    uint64_t index = (buf->len++) * size;
    ASSERT (index + size <= buf->size, "called from %s:%u: NEXTENT went beyond end of buffer: size=%u index=%"PRIu64": %.*s", 
            func, code_line, (unsigned)size, index, (int)sizeof(BufDescType)-1, &buf_desc (buf).s[0]);
    return index;
}
#define NEXTENT(type, buf)    (*(type *)(&(buf).data[NEXTENT_get_index (&(buf), sizeof(type), __FUNCTION__, __LINE__)]))
#define ENTNUM(buf, ent)      ((int32_t)((((char*)(ent)) - ((buf).data)) / (int32_t)sizeof (*(ent)))) // signed integer
#define REMAINING(buf, ent)   ((buf).len - ENTNUM ((buf),(ent)))
#define ISLASTENT(buf,ent)    (ENTNUM((buf),(ent)) == (buf).len - 1)
#define ISAFTERENT(buf,ent)   (!(buf).data || ENTNUM((buf),(ent)) >= (buf).len)
#define ISFIRSTENT(buf,ent)   (ENTNUM((buf),(ent)) == 0)
#define ISBEFOREENT(buf,ent)  (!(buf).data || ENTNUM((buf),(ent)) <= -1)

extern void buf_initialize(void);

#define buf_is_alloc(buf_p) ((buf_p)->data != NULL && (buf_p)->type != BUF_UNALLOCATED)
#define ASSERTNOTINUSE(buf) ASSERT (!buf_is_alloc (&(buf)) && !(buf).len, "expecting "#buf" to be free, but it's not: %s", buf_desc (&(buf)).s)
#define ASSERTISALLOCED(buf) ASSERT0 (buf_is_alloc (&(buf)), #buf" is not allocated")
#define ASSERTNOTEMPTY(buf) ASSERT (buf_is_alloc (&(buf)) && (buf).len, "expecting "#buf" to be contain some data, but it doesn't: %s", buf_desc (&(buf)).s)

extern uint64_t buf_alloc_do (VBlockP vb,
                              Buffer *buf, 
                              uint64_t requested_size, 
                              float grow_at_least_factor, // grow more than new_size   
                              const char *func, uint32_t code_line,
                              const char *name);

#define buf_alloc(alloc_vb, buf, more, at_least, type, grow_at_least_factor, name) do {\
    uint64_t new_req_size = MAX_((at_least), ((buf)->len+(more)))*sizeof(type); /* make copy to allow ++ */ \
    ((!(buf)->data || (buf)->size < (new_req_size)) ? buf_alloc_do (((alloc_vb) ? ((VBlockP)alloc_vb) : (buf)->vb), (buf), (new_req_size), (grow_at_least_factor), __FUNCTION__, __LINE__, (name)) \
                                                    : (buf)->size); \
} while(0)

#define buf_alloc_zero(vb, buf, more, at_least, element_type, grow_at_least_factor,name) do { \
    uint64_t size_before = (buf)->data ? (buf)->size : 0; /* always zero the whole buffer in an initial allocation */ \
    buf_alloc((vb), (buf), (more), (at_least), element_type, (grow_at_least_factor), (name)); \
    if ((buf)->data && (buf)->size > size_before) memset (&(buf)->data[size_before], 0, (buf)->size - size_before); \
} while(0)

#define buf_alloc_255(vb, buf, more, at_least, element_type, grow_at_least_factor,name) do { \
    uint64_t size_before = (buf)->data ? (buf)->size : 0; /* always zero the whole buffer in an initial allocation */ \
    buf_alloc((vb), (buf), (more), (at_least), element_type, (grow_at_least_factor), (name)); \
    if ((buf)->data && (buf)->size > size_before) memset (&(buf)->data[size_before], 255, (buf)->size - size_before); \
} while(0)

extern bool buf_mmap_do (VBlockP vb, Buffer *buf, const char *filename, bool read_only_buffer, const char *func, uint32_t code_line, const char *name);
#define buf_mmap(vb, buf, filename, read_only_buffer, name) \
    buf_mmap_do((VBlockP)(vb), (buf), (filename), (read_only_buffer), __FUNCTION__, __LINE__, (name))

extern void buf_free_do (Buffer *buf, const char *func, uint32_t code_line);
#define buf_free(buf) buf_free_do ((buf), __FUNCTION__, __LINE__)

extern void buf_destroy_do (Buffer *buf, const char *func, uint32_t code_line);
#define buf_destroy(buf) buf_destroy_do ((buf), __FUNCTION__, __LINE__)

#define buf_is_large_enough(buf_p, requested_size) (buf_is_alloc ((buf_p)) && (buf_p)->size >= requested_size)

extern void buf_copy_do (VBlockP dst_vb, Buffer *dst, const Buffer *src, uint64_t bytes_per_entry,
                         uint64_t src_start_entry, uint64_t max_entries, // if 0 copies the entire buffer
                         const char *func, uint32_t code_line,
                         const char *name);
#define buf_copy(dst_vb,dst,src,type,src_start_entry,max_entries,dst_name) \
  buf_copy_do ((VBlockP)(dst_vb),(dst),(src),sizeof(type),(src_start_entry),(max_entries),__FUNCTION__,__LINE__,(dst_name))

extern void buf_move (VBlockP dst_vb, Buffer *dst, VBlockP src_vb, Buffer *src);
extern void buf_grab_do (VBlockP dst_vb, Buffer *dst_buf, const char *dst_name, Buffer *src_buf, const char *func, uint32_t code_line);
#define buf_grab(dst_vb, dst_buf, dst_name, src_buf) buf_grab_do ((dst_vb), (dst_buf), (dst_name), (src_buf), __FUNCTION__, __LINE__)

extern void buf_cut_out_do (BufferP buf, unsigned sizeof_item, uint64_t remove_start, uint64_t remove_len);
#define buf_remove(buf, type, remove_start, remove_len) buf_cut_out_do ((buf), sizeof(type), (remove_start), (remove_len))

#define buf_has_space(buf, new_len) ((buf)->len + (new_len) <= (buf)->size)

#define buf_add(buf, new_data, new_data_len) do { \
    uint32_t new_len = (uint32_t)(new_data_len); /* copy in case caller uses ++ */ \
    ASSERT (buf_has_space(buf, new_len), \
            "Error in buf_add: buffer %s is out of space: len=%u size=%u new_data_len=%u", \
            buf_desc (buf).s, (uint32_t)(buf)->len, (uint32_t)(buf)->size, new_len);\
    memcpy (&(buf)->data[(buf)->len], (new_data), new_len);   \
    (buf)->len += new_len; \
} while(0)

static inline void buf_add_more_ex (VBlockP vb, BufferP buf, unsigned width, const void *new_data, unsigned new_data_len, const char *name) 
{ 
    if (!new_data_len) return;

    buf_alloc (vb ? vb : buf->vb, buf, new_data_len, (buf->len + new_data_len) * width +1 /* +1 - room for \0 or separator */, char, 1.5, name); 
    memcpy (&buf->data[buf->len * width], new_data, new_data_len * width);   
    buf->len += new_data_len; 
} 
static inline void buf_add_more (VBlockP vb, BufferP buf, STRp(new_data), const char *name) 
{   buf_add_more_ex (vb, buf, 1, STRa(new_data), name); }

#define buf_add_more_(vb, buf, type, new_data, new_data_len, name) \
    buf_add_more_ex ((VBlockP)(vb), (buf), sizeof(type), (new_data), (new_data_len), (name))

#define buf_add_moreC(vb_, buf, literal_str, name) buf_add_more ((VBlockP)(vb_), (buf), literal_str, sizeof literal_str-1, (name))
#define buf_add_moreS(vb_, buf, str, name) buf_add_more ((VBlockP)(vb_), (buf), str, str##_len, (name))
#define buf_add_buf(vb_,dst_buf,src_buf,type,name) do { \
    buf_alloc ((vb_) ? (vb_) : (dst_buf)->vb, (dst_buf), (src_buf)->len, 0, type, CTX_GROWTH, (name)); \
    memcpy (AFTERENT(type, *(dst_buf)), (src_buf)->data, (src_buf)->len * sizeof (type));   \
    (dst_buf)->len += (src_buf)->len; \
} while (0)

extern void buf_add_string (VBlockP vb, Buffer *buf, const char *str);
extern void buf_add_int_as_text (VBlockP vb, Buffer *buf, int64_t value);

#define buf_add_int(vb, buf, n) ({ buf_alloc ((VBlockP)(vb), &(buf), 1, 0, typeof(n), 0, NULL); NEXTENT(typeof(n), (buf)) = n; })

#define buf_add_vector(vb, type, dst, src, name) do { \
    if ((src).len) { \
        if (!(dst).len) \
            buf_copy ((vb), &(dst), &(src), sizeof (type), 0, 0, name); \
        else \
            for (uint64_t i=0; i < (src).len; i++) \
                *ENT (type, (dst), i) += *ENT (type, (src), i); \
    } \
} while(0)

#define BUFPRINTF_MAX_LEN 5000
#define bufprintf(vb, buf, format, ...)  do { char __s[BUFPRINTF_MAX_LEN]; sprintf (__s, (format), __VA_ARGS__); buf_add_string ((VBlockP)(vb), (buf), __s); } while (0)
#define bufprint0 buf_add_string 

extern void buf_print (Buffer *buf, bool add_newline);

extern void buf_test_overflows (void *vb, const char *msg);
extern void buf_test_overflows_all_vbs (const char *msg);

typedef struct {
    const char *name;
    uint64_t bytes; 
    unsigned buffers;
} MemStats;

uint64_t buf_get_memory_usage (void);
extern void buf_show_memory (bool memory_full, unsigned max_threads, unsigned used_threads);
extern void buf_show_memory_handler (void);

#define buf_set(buf_p,value) do { if ((buf_p)->data) memset ((buf_p)->data, value, (buf_p)->size); } while(0)
#define buf_zero(buf_p) buf_set(buf_p, 0)

extern void buf_add_to_buffer_list_do (VBlockP vb, Buffer *buf, const char *func);
#define buf_add_to_buffer_list(vb,buf) buf_add_to_buffer_list_do ((vb), (buf), __FUNCTION__)

extern void buf_update_buf_list_vb_addr_change (VBlockP new_vb, VBlockP old_vb);

extern void buf_compact_buf_list (VBlockP vb);

extern void buf_low_level_free (void *p, const char *func, uint32_t code_line);
#define FREE(p) do { if (p) { buf_low_level_free (((void*)(p)), __FUNCTION__, __LINE__); p=NULL; } } while(0)

extern void *buf_low_level_malloc (size_t size, bool zero, const char *func, uint32_t code_line);
#define MALLOC(size) buf_low_level_malloc (size, false, __FUNCTION__, __LINE__)
#define CALLOC(size) buf_low_level_malloc (size, true,  __FUNCTION__, __LINE__)

extern void *buf_low_level_realloc (void *p, size_t size, const char *name, const char *func, uint32_t code_line);
#define REALLOC(p,size,name) if (!(*(p) = buf_low_level_realloc (*(p), (size), (name), __FUNCTION__, __LINE__))) ABORT0 ("REALLOC failed")

extern bool buf_dump_to_file (const char *filename, const Buffer *buf, unsigned buf_word_width, bool including_control_region, bool no_dirs, bool verbose, bool do_gzip);

// ------------
// bitmap stuff
// ------------
extern BitArrayP buf_alloc_bitarr_do (VBlockP vb, Buffer *buf, uint64_t nbits, const char *func, uint32_t code_line, const char *name);
#define buf_alloc_bitarr(vb, buf, nbits, name) \
    buf_alloc_bitarr_do((VBlockP)(vb), (buf), (nbits), __FUNCTION__, __LINE__, (name))

extern BitArrayP buf_overlay_bitarr_do (VBlockP vb, Buffer *overlaid_buf, Buffer *regular_buf, uint64_t start_byte_in_regular_buf, uint64_t nbits, const char *func, uint32_t code_line, const char *name);
#define buf_overlay_bitarr(vb, overlaid_buf, regular_buf, start_byte_in_regular_buf, nbits, name) \
    buf_overlay_bitarr_do((VBlockP)(vb), (overlaid_buf), (regular_buf), (start_byte_in_regular_buf), (nbits), __FUNCTION__, __LINE__, (name))

#define buf_set_overlayable(buf) (buf)->overlayable = true

extern void buf_overlay_do (VBlockP vb, Buffer *overlaid_buf, Buffer *regular_buf,  uint64_t start_in_regular,  
                            const char *func, uint32_t code_line, const char *name);
#define buf_overlay(vb, overlaid_buf, regular_buf, name) \
    buf_overlay_do((VBlockP)(vb), (overlaid_buf), (regular_buf), 0, __FUNCTION__, __LINE__, (name)) 

#define buf_overlay_partial(vb, overlaid_buf, regular_buf, start, name) \
    buf_overlay_do((VBlockP)(vb), (overlaid_buf), (regular_buf), (start), __FUNCTION__, __LINE__, (name)) 

extern uint64_t buf_extend_bits (Buffer *buf, int64_t num_new_bits);
extern void buf_add_bit (Buffer *buf, int64_t new_bit);
extern BitArrayP buf_zfile_buf_to_bitarray (Buffer *buf, uint64_t nbits);
#define buf_get_bitarray(buf) ((BitArrayP)(&(buf)->type))
#define buf_add_set_bit(buf)   buf_add_bit (buf, 1)
#define buf_add_clear_bit(buf) buf_add_bit (buf, 0)

// endianness
extern BgEnBufFunc BGEN_u8_buf, BGEN_u16_buf, BGEN_u32_buf, BGEN_u64_buf, 
                   BGEN_transpose_u8_buf, BGEN_transpose_u16_buf, BGEN_transpose_u32_buf, BGEN_transpose_u64_buf,
                   BGEN_deinterlace_d8_buf, BGEN_deinterlace_d16_buf, BGEN_deinterlace_d32_buf, BGEN_deinterlace_d64_buf,
                   interlace_d8_buf, BGEN_interlace_d16_buf, BGEN_interlace_d32_buf, BGEN_interlace_d64_buf,
                   LTEN_u16_buf, LTEN_u32_buf, LTEN_u64_buf, LTEN_interlace_d16_buf, LTEN_interlace_d32_buf, LTEN_interlace_d64_buf;
                   
