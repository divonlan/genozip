// ------------------------------------------------------------------
//   buffer.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include "genozip.h"
#include "strings.h"
#include "endianness.h"

typedef enum __attribute__ ((__packed__)) { // 1 byte 
    BUF_UNALLOCATED=0, BUF_REGULAR, BUF_OVERLAY, BUF_OVERLAY_RO, BUF_MMAP, BUF_MMAP_RO, BUF_STANDALONE_BITARRAY, BUF_NUM_TYPES 
} BufferType; // BUF_UNALLOCATED must be 0, must be identical to BitsType

#define BUFTYPE_NAMES { "UNALLOCATED", "REGULAR", "OVERLAY", "MMAP" }

typedef struct Buffer { // 64 bytes
    //------------------------------------------------------------------------------------------------------
    // these 4 fields, must be first, and can be overlayed with a Bits - their order and size is identical to Bits
    union {
        char *data;            // ==memory+8 if buffer is allocated or NULL if not
        uint64_t *words;       // for Bits
    };
    union {                    // the "parameter" field is for discretionary use by the caller. some options to use the parameter are provided.
        int64_t param;    
        int64_t next;     
        int64_t count;    
        uint64_t nbits;        // for Bits
        uint32_t prm32[2];
        uint16_t prm16[4];
        uint8_t  prm8 [8];
#ifdef _WIN32
        void *mmap_handle;     // handle to memory mmaped object, used by BUF_MMAP and BUF_MMAP_RO (Windows)
#endif        
    };
    union {
        uint64_t len;          // used by the buffer user according to its internal logic. not modified by malloc/realloc, zeroed by buf_free (in Bits - nwords)
        uint64_t nwords;       // for Bits
        struct {
#ifdef __LITTLE_ENDIAN__
            uint32_t len32, len32hi; // we use len32 instead of len if easier, in cases where we are certain len is smaller than 4B
#else                
            uint32_t len32hi, len32;
#endif  
        };
    };
    uint64_t type        : 3;  // BufferType
    //------------------------------------------------------------------------------------------------------
    uint64_t overlayable : 1;  // this buffer may be fully overlaid by one or more overlay buffers
    uint64_t can_be_big  : 1;  // do not display warning if buffer grows very big
    #define BG_LOADING_BIT 5   // bit # of bg_loading within the uint64_t 
    uint64_t bg_loading  : 1;  // mmap thread is loading in the background
    uint64_t code_line   : 12; // the allocating line number in source code file (up to 4095)
    #define MAX_BUFFER_SIZE ((1ULL<<46)-1) // according to bits of "size" (64 TB)
    uint64_t size        : 46; // number of bytes available to the user (i.e. not including the allocated overhead). 

    VBlockP vb;                // vb that owns this buffer, and which this buffer is in its buf_list

    rom name;                  // name of allocator - used for memory debugging & statistics.
    rom func;                  // the allocating function
    pthread_t bg_loader;       // thread loading mmap data (used by BUF_MMAP and BUF_MMAP_RO)

    char *memory;              // memory allocated to this buffer - amount is: size + 2*sizeof(longlong) to allow for OVERFLOW and UNDERFLOW)

} Buffer; 

#define EMPTY_BUFFER ((Buffer){})

extern void buf_initialize(void);

#define buf_is_alloc(buf_p) ((buf_p)->data != NULL && (buf_p)->type != BUF_UNALLOCATED)
#define ASSERTNOTINUSE(buf) ASSERT (!buf_is_alloc (&(buf)) && !(buf).len && !(buf).param, "expecting "#buf" to be free, but it's not: %s", buf_desc (&(buf)).s)
#define ASSERTISALLOCED(buf) ASSERT0 (buf_is_alloc (&(buf)), #buf" is not allocated")
#define ASSERTISEMPTY(buf) ASSERT (buf_is_alloc (&(buf)) && !(buf).len, "expecting "#buf" to be be allocated an empty, but it isn't: %s", buf_desc (&(buf)).s)
#define ASSERTNOTEMPTY(buf) ASSERT (buf_is_alloc (&(buf)) && (buf).len, "expecting "#buf" to be contain some data, but it doesn't: %s", buf_desc (&(buf)).s)

extern uint64_t buf_alloc_do (VBlockP vb,
                              BufferP buf, 
                              uint64_t requested_size, 
                              float grow_at_least_factor, // grow more than new_size   
                              FUNCLINE,
                              rom name);

#define buf_alloc(alloc_vb, buf, more, at_least, type, grow_at_least_factor, name) ({\
    uint64_t new_req_size = MAX_((at_least), ((buf)->len+(more)))*sizeof(type); /* make copy to allow ++ */ \
    ((!(buf)->data || (buf)->size < (new_req_size)) ? buf_alloc_do (((alloc_vb) ? ((VBlockP)alloc_vb) : (buf)->vb), (buf), (new_req_size), (grow_at_least_factor), __FUNCLINE, (name)) \
                                                    : (buf)->size); })

// alloc a set amount of bytes, and set buf.len
#define buf_alloc_exact(alloc_vb, buf, exact_len, type, name) ({  \
    buf_alloc((alloc_vb), &(buf), 0, (exact_len), type, 1, name); \
    (buf).len = (exact_len); })

#define buf_alloc_exact_zero(alloc_vb, buf, exact_len, type, name) ({   \
    buf_alloc_exact (alloc_vb, buf, exact_len, type, name); \
    memset ((buf).data, 0, (exact_len) * sizeof(type)); })

#define buf_alloc_exact_255(alloc_vb, buf, exact_len, type, name) ({   \
    buf_alloc_exact (alloc_vb, buf, exact_len, type, name); \
    memset ((buf).data, 255, (exact_len) * sizeof(type)); })

// allocates exactly the requested amount and sets let, and declares ARRAY
#define ARRAY_alloc(element_type, array_name, array_len, init_zero, buf, alloc_vb, buf_name) \
    (buf).len = (array_len); \
    if (!(buf).data || (buf).size < ((buf).len) * sizeof(element_type)) \
        buf_alloc_do (((alloc_vb) ? ((VBlockP)alloc_vb) : (buf).vb), &(buf), ((buf).len * sizeof(element_type)), 1, __FUNCLINE, (buf_name)); \
    if (init_zero) memset ((buf).data, 0, (buf).len * sizeof(element_type)); /* resets the entire buffer, not just newly allocated memory */ \
    element_type *array_name = ((element_type *)((buf).data)); \
    const uint64_t array_name##_len __attribute__((unused)) = (buf).len; // read-only copy of len 

#define buf_alloc_zero(vb, buf, more, at_least, element_type, grow_at_least_factor,name) ({ \
    uint64_t size_before = (buf)->data ? (buf)->size : 0; /* always zero the whole buffer in an initial allocation */ \
    buf_alloc((vb), (buf), (more), (at_least), element_type, (grow_at_least_factor), (name)); \
    if ((buf)->data && (buf)->size > size_before) memset (&(buf)->data[size_before], 0, (buf)->size - size_before); })

#define buf_alloc_255(vb, buf, more, at_least, element_type, grow_at_least_factor,name) ({ \
    uint64_t size_before = (buf)->data ? (buf)->size : 0; /* always zero the whole buffer in an initial allocation */ \
    buf_alloc((vb), (buf), (more), (at_least), element_type, (grow_at_least_factor), (name)); \
    if ((buf)->data && (buf)->size > size_before) memset (&(buf)->data[size_before], 255, (buf)->size - size_before); })

#define ARRAY(element_type, name, buf) \
    element_type *name = ((element_type *)((buf).data)); \
    const uint64_t name##_len __attribute__((unused)) = (buf).len; // read-only copy of len 

#define ARRAY32(element_type, name, buf) \
    element_type *name = ((element_type *)((buf).data)); \
    const uint32_t name##_len __attribute__((unused)) = (buf).len32; // read-only copy of len 

// entry #index in Buffer
#define B(type, buf, index) ((type *)(&(buf).data[(index) * sizeof(type)]))
#define Bc(buf, index)      B(char,     (buf), (index))
#define B8(buf, index)      B(uint8_t,  (buf), (index))
#define B16(buf, index)     B(uint16_t, (buf), (index))
#define B32(buf, index)     B(uint32_t, (buf), (index))
#define B64(buf, index)     B(uint64_t, (buf), (index))
#define Btxt(index)         Bc (vb->txt_data, index)

// first entry in Buffer
#define B1ST(type, buf)     ((type     *)(buf).data)
#define B1STc(buf)          ((char     *)(buf).data)
#define B1ST8(buf)          ((uint8_t  *)(buf).data)
#define B1ST16(buf)         ((uint16_t *)(buf).data)
#define B1ST32(buf)         ((uint32_t *)(buf).data)
#define B1ST64(buf)         ((uint64_t *)(buf).data)
#define B1STtxt             (vb->txt_data.data)

// last entry in Buffer
#define BLST(type, buf)     ((type *)(&(buf).data[((buf).len-1) * sizeof(type)]))
#define BLSTc(buf)          BLST(char,    (buf))
#define BLST8(buf)          BLST(uint8_t, (buf))
#define BLST16(buf)         BLST(uint16_t,(buf))
#define BLST32(buf)         BLST(uint32_t,(buf))
#define BLST64(buf)         BLST(uint64_t,(buf))
#define BLSTtxt             (&vb->txt_data.data[(vb->txt_data.len-1)])

// entry after the end of the Buffer 
#define BAFT(type, buf)     ((type *)(&(buf).data[((buf).len) * sizeof(type)]))
#define BAFTc(buf)          BAFT(char,    (buf))
#define BAFT8(buf)          BAFT(uint8_t, (buf))
#define BAFT16(buf)         BAFT(uint16_t,(buf))
#define BAFT32(buf)         BAFT(uint32_t,(buf))
#define BAFT64(buf)         BAFT(uint64_t,(buf))
#define BAFTtxt             (&vb->txt_data.data[vb->txt_data.len])

#define for_buf(element_type, iterator, buf)  \
    for (element_type *iterator=B1ST(element_type, (buf)); iterator < BAFT(element_type, (buf)); iterator++)

#define for_buf_back(element_type, iterator, buf)  \
    for (element_type *iterator=BLST(element_type, (buf)); iterator >= B1ST(element_type, (buf)); iterator--)

// loop with two concurrent iterators "iter_p" (pointer to element_type) and "iter_i" (32bit) 
#define for_buf2(element_type, iter_p, iter_i, buf) \
    for (uint32_t iter_i=0; iter_i < (buf).len32;)  \
        for (element_type *iter_p=B1ST(element_type, (buf)); iter_p < BAFT(element_type, (buf)); iter_p++, iter_i++)

#define for_buf2_back(element_type, iter_p, iter_i, buf) \
    for (int32_t iter_i=(buf).len32-1; iter_i >= 0;)  \
        for (element_type *iter_p=BLST(element_type, (buf)); iter_p >= B1ST(element_type, (buf)); iter_p--, iter_i--)

typedef struct { char s[300]; } BufDescType;
extern const BufDescType buf_desc (ConstBufferP buf);

static inline uint64_t BNXT_get_index (BufferP buf, size_t size, FUNCLINE) 
{ 
    uint64_t index = (buf->len++) * size;
    ASSERT (index + size <= buf->size, "called from %s:%u: BNXT went beyond end of buffer: size=%u index=%"PRIu64": %.*s", 
            func, code_line, (unsigned)size, index, (int)sizeof(BufDescType)-1, &buf_desc (buf).s[0]);
    return index;
}

// next entry in Buffer
#define BNXT(type, buf)     (*(type *)(&(buf).data[BNXT_get_index (&(buf), sizeof(type), __FUNCLINE)]))
#define BNXTc(buf)          (buf).data[BNXT_get_index (&(buf), 1, __FUNCLINE)]
#define BNXT8(buf)          BNXT(uint8_t, (buf))
#define BNXT16(buf)         BNXT(uint16_t, (buf))
#define BNXT32(buf)         BNXT(uint32_t, (buf))
#define BNXT64(buf)         BNXT(uint64_t, (buf))
#define BNUM(buf, ent)      ((int32_t)((((char*)(ent)) - ((buf).data)) / (int32_t)sizeof (*(ent)))) // signed integer
#define BNUM64(buf, ent)    ((int64_t)((((char*)(ent)) - ((buf).data)) / (int64_t)sizeof (*(ent)))) // signed integer
#define BNUMtxt(ent)        BNUM(vb->txt_data, (ent))
#define BREMAINS(buf, ent)  ((buf).len - BNUM ((buf),(ent)))
#define BIS1ST(buf,ent)     (BNUM((buf),(ent)) == 0)
#define BISLST(buf,ent)     (BNUM((buf),(ent)) == (buf).len - 1)
#define BISBEFORE(buf,ent)  (BNUM((buf),(ent)) <= -1)
#define BISAFT(buf,ent)     (BNUM((buf),(ent)) >= (buf).len)
#define BISVALID(buf, ent)  (!BISBEFORE((buf),(ent)) && !BISAFT((buf),(ent)))
#define ASSBISVALID(buf,ent) ASSERT (BISVALID((buf), (ent)), #ent "=%p is not in its Buffer " #buf, (ent))

#define INSERBtxtAFTER(type, buf, index) ({             \
    buf_alloc ((buf).vb, &(buf), 1, 0, type, 1, NULL);  \
    memmove (B(type, (buf), (index)+2), B(type, (buf), (index)+1), ((buf).len - ((index)+1)) * sizeof(type)); \
    (buf).len++;                                        \
    B(type, (buf), (index)+1);                          \
})

extern bool buf_mmap_do (VBlockP vb, BufferP buf, rom filename, bool read_only_buffer, FUNCLINE, rom name);
#define buf_mmap(vb, buf, filename, read_only_buffer, name) \
    buf_mmap_do((VBlockP)(vb), (buf), (filename), (read_only_buffer), __FUNCLINE, (name))

extern void buf_free_do (BufferP buf, FUNCLINE);
#define buf_free(buf) buf_free_do (&(buf), __FUNCLINE)

extern void buf_destroy_do (BufferP buf, FUNCLINE);
#define buf_destroy(buf) buf_destroy_do (&(buf), __FUNCLINE)

#define buf_is_large_enough(buf_p, requested_size) (buf_is_alloc ((buf_p)) && (buf_p)->size >= requested_size)

extern void buf_copy_do (VBlockP dst_vb, BufferP dst, ConstBufferP src, uint64_t bytes_per_entry,
                         uint64_t src_start_entry, uint64_t max_entries, // if 0 copies the entire buffer
                         FUNCLINE,
                         rom name);
#define buf_copy(dst_vb,dst,src,type,src_start_entry,max_entries,dst_name) \
  buf_copy_do ((VBlockP)(dst_vb),(dst),(src),sizeof(type),(src_start_entry),(max_entries),__FUNCTION__,__LINE__,(dst_name))

extern void buf_move (VBlockP dst_vb, BufferP dst, VBlockP src_vb, BufferP src);
extern void buf_grab_do (VBlockP dst_vb, BufferP dst_buf, rom dst_name, BufferP src_buf, FUNCLINE);
#define buf_grab(dst_vb, dst_buf, dst_name, src_buf) buf_grab_do ((dst_vb), &(dst_buf), (dst_name), &(src_buf), __FUNCLINE)

extern void buf_remove_do (BufferP buf, unsigned sizeof_item, uint64_t remove_start, uint64_t remove_len);
#define buf_remove(buf, type, remove_start, remove_len) buf_remove_do (&(buf), sizeof(type), (remove_start), (remove_len))

extern void buf_insert_do (VBlockP vb, BufferP buf, unsigned width, uint64_t insert_at, const void *new_data, uint64_t new_data_len, rom name);

#define buf_has_space(buf, new_len) ((buf)->len + (new_len) <= (buf)->size)

#define buf_add(buf, new_data, new_data_len)  \
    ({ uint32_t new_len = (uint32_t)(new_data_len); /* copy in case caller uses ++ */ \
       ASSERT (buf_has_space(buf, new_len), \
            "buf_add: buffer %s is out of space: len=%u size=%u new_data_len=%u", \
            buf_desc (buf).s, (uint32_t)(buf)->len, (uint32_t)(buf)->size, new_len);\
       memcpy (&(buf)->data[(buf)->len], (new_data), new_len);   \
       (buf)->len += new_len; })

static inline void buf_add_more (VBlockP vb, BufferP buf, STRp(new_data), rom name) 
{   buf_insert_do (vb, buf, 1, buf->len, STRa(new_data), name); }

#define buf_append(vb, buf, type, new_data, new_data_len, name) \
    buf_insert_do ((VBlockP)(vb), &(buf), sizeof(type), (buf).len, (new_data), (new_data_len), (name))

#define buf_append_one(buf, item) ({ \
    typeof(item) item_ = (item); /* evaluate once */\
    buf_insert_do (NULL, &(buf), sizeof(typeof(item)), (buf).len, &item_, 1, NULL);\
})

#define buf_insert(vb, buf, type, insert_at, new_data, new_data_len, name) \
    buf_insert_do ((VBlockP)(vb), &(buf), sizeof(type), (insert_at), (new_data), (new_data_len), (name))

#define buf_add_moreC(vb_, buf, literal_str, name) buf_add_more ((VBlockP)(vb_), (buf), literal_str, sizeof literal_str-1, (name))
#define buf_add_moreS(vb_, buf, str, name) buf_add_more ((VBlockP)(vb_), (buf), str, str##_len, (name))
#define buf_add_buf(vb_,dst_buf,src_buf,type,name) ({ \
    buf_alloc ((vb_) ? (VBlockP)(vb_) : (dst_buf)->vb, (dst_buf), (src_buf)->len, 0, type, CTX_GROWTH, (name)); \
    memcpy (BAFT(type, *(dst_buf)), (src_buf)->data, (src_buf)->len * sizeof (type));   \
    (dst_buf)->len += (src_buf)->len; })

extern void buf_append_string (VBlockP vb, BufferP buf, rom str);

// adds a textual int at the end of a buffer. caller needs to allocate enough space in the buffer.
static inline unsigned buf_add_int_as_text (BufferP buf, int64_t n)
{
    unsigned n_len = str_int_ex (n, BAFTc (*buf), false); 
    buf->len += n_len; 
    return n_len;
}

// adds a textual int at the end of a buffer. caller needs to allocate enough space in the buffer.
static inline unsigned buf_add_hex_as_text (BufferP buf, int64_t n, bool uppercase)
{
    unsigned n_len = str_hex_ex (n, BAFTc (*buf), uppercase, false); 
    buf->len += n_len; 
    return n_len;
}

#define buf_add_int(vb, buf, n) ({ buf_alloc ((VBlockP)(vb), &(buf), 1, 0, typeof(n), 0, NULL); BNXT(typeof(n), (buf)) = n; })

#define buf_add_vector(vb, type, dst, src, name)  \
({  if ((src).len) { \
        if (!(dst).len) \
            buf_copy ((vb), &(dst), &(src), sizeof (type), 0, 0, name); \
        else \
            for (uint64_t i=0; i < (src).len; i++) \
                *B(type, (dst), i) += *B(type, (src), i); \
    } })

#define BUFPRINTF_MAX_LEN 5000
#define bufprintf(vb, buf, format, ...)  ({ char __s[BUFPRINTF_MAX_LEN+1]; \
                                            __s[BUFPRINTF_MAX_LEN] = '#';  /* overflow fence */ \
                                            sprintf (__s, (format), __VA_ARGS__); \
                                            ASSERT (__s[BUFPRINTF_MAX_LEN] == '#', "bufprintf: String too long - stack overflow: len=%u", (int)strlen (__s)); \
                                            buf_append_string ((VBlockP)(vb), (buf), __s); })
#define bufprint0 buf_append_string  // note: the string isn't a printf format, so escaping works differently (eg % doesn't need to be %%)

extern void buf_print (BufferP buf, bool add_newline);

extern void buf_test_overflows (void *vb, rom msg);
extern void buf_test_overflows_all_vbs (rom msg);

typedef struct {
    rom name;
    uint64_t bytes; 
    unsigned buffers;
} MemStats;

uint64_t buf_get_memory_usage (void);
extern void buf_show_memory (bool memory_full, unsigned max_threads, unsigned used_threads);
extern void buf_show_memory_handler (void);

#define buf_set(buf_p,value) ({ if ((buf_p)->data) memset ((buf_p)->data, value, (buf_p)->size); })
#define buf_zero(buf_p) buf_set(buf_p, 0)

extern void buf_add_to_buffer_list_do (VBlockP vb, BufferP buf, FUNCLINE);
#define buf_add_to_buffer_list(vb,buf) buf_add_to_buffer_list_do ((VBlockP)(vb), (buf), __FUNCLINE)
#define buf_add_to_buffer_list_(vb,buf,buf_name) ({ buf_add_to_buffer_list_do ((VBlockP)(vb), (buf), __FUNCLINE); (buf)->name = (buf_name); })

extern void buf_update_buf_list_vb_addr_change (VBlockP new_vb, VBlockP old_vb);

extern void buf_compact_buf_list (VBlockP vb);

extern void buf_low_level_free (void *p, FUNCLINE);
#define FREE(p) ({ if (p) { buf_low_level_free (((void*)(p)), __FUNCLINE); p=NULL; } })

extern void *buf_low_level_malloc (size_t size, bool zero, FUNCLINE);
#define MALLOC(size) buf_low_level_malloc (size, false, __FUNCLINE)
#define CALLOC(size) buf_low_level_malloc (size, true,  __FUNCLINE)

extern void *buf_low_level_realloc (void *p, size_t size, rom name, FUNCLINE);
#define REALLOC(p,size,name) if (!(*(p) = buf_low_level_realloc (*(p), (size), (name), __FUNCLINE))) ABORT0 ("REALLOC failed")

extern bool buf_dump_to_file (rom filename, ConstBufferP buf, unsigned buf_word_width, bool including_control_region, bool no_dirs, bool verbose, bool do_gzip);

// ------------
// bitmap stuff
// ------------
// allocate bit array and set nbits
extern BitsP buf_alloc_bits_do (VBlockP vb, BufferP buf, uint64_t nbits, FUNCLINE, rom name);
#define buf_alloc_bits(vb, buf, nbits, name) \
    buf_alloc_bits_do((VBlockP)(vb), (buf), (nbits), __FUNCLINE, (name))

// allocate bit array without setting nbits
extern void buf_alloc_bits_buffer_do (VBlockP vb, BufferP buf, uint64_t nbits, FUNCLINE, rom name);
#define buf_alloc_bits_buffer(vb, buf, nbits, name) \
    buf_alloc_bits_buffer_do((VBlockP)(vb), (buf), (nbits), __FUNCLINE, (name))

extern BitsP buf_overlay_bitarr_do (VBlockP vb, BufferP overlaid_buf, BufferP regular_buf, uint64_t start_byte_in_regular_buf, uint64_t nbits, FUNCLINE, rom name);
#define buf_overlay_bitarr(vb, overlaid_buf, regular_buf, start_byte_in_regular_buf, nbits, name) \
    buf_overlay_bitarr_do((VBlockP)(vb), (overlaid_buf), (regular_buf), (start_byte_in_regular_buf), (nbits), __FUNCLINE, (name))

#define buf_set_overlayable(buf) (buf)->overlayable = true

extern void buf_overlay_do (VBlockP vb, BufferP overlaid_buf, BufferP regular_buf,  uint64_t start_in_regular,  
                            FUNCLINE, rom name);
#define buf_overlay(vb, overlaid_buf, regular_buf, name) \
    buf_overlay_do((VBlockP)(vb), (overlaid_buf), (regular_buf), 0, __FUNCLINE, (name)) 

#define buf_overlay_partial(vb, overlaid_buf, regular_buf, start_in_regular, name) \
    buf_overlay_do((VBlockP)(vb), (overlaid_buf), (regular_buf), (start_in_regular), __FUNCLINE, (name)) 

extern uint64_t buf_extend_bits (BufferP buf, int64_t num_new_bits);
extern void buf_add_bit (BufferP buf, int64_t new_bit);
extern BitsP buf_zfile_buf_to_bitarray (BufferP buf, uint64_t nbits);
#define buf_add_set_bit(buf)   buf_add_bit (buf, 1)
#define buf_add_clear_bit(buf) buf_add_bit (buf, 0)

#define buf_get_buffer_from_bits(bitarr) ((BufferP)(bitarr))

extern char *buf_display (ConstBufferP buf);

// endianness
extern BgEnBufFunc BGEN_u8_buf, BGEN_u16_buf, BGEN_u32_buf, BGEN_u64_buf, 
                   BGEN_transpose_u8_buf, BGEN_transpose_u16_buf, BGEN_transpose_u32_buf, BGEN_transpose_u64_buf,
                   BGEN_deinterlace_d8_buf, BGEN_deinterlace_d16_buf, BGEN_deinterlace_d32_buf, BGEN_deinterlace_d64_buf,
                   interlace_d8_buf, BGEN_interlace_d16_buf, BGEN_interlace_d32_buf, BGEN_interlace_d64_buf,
                   LTEN_u16_buf, LTEN_u32_buf, LTEN_u64_buf, LTEN_interlace_d16_buf, LTEN_interlace_d32_buf, LTEN_interlace_d64_buf;
                   
    
