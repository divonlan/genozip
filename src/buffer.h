// ------------------------------------------------------------------
//   buffer.h
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "endianness.h"
#include "buf_struct.h"
#include "local_type.h"

#define ARRAY(element_type, name, buf) \
    element_type *name = ((element_type *)((buf).data)); \
    const uint64_t name##_len __attribute__((unused)) = (buf).len   // read-only copy of len 

#define ARRAY32(element_type, name, buf) \
    element_type *name = ((element_type *)((buf).data)); \
    const uint32_t name##_len __attribute__((unused)) = (buf).len32 // read-only copy of len 

#pragma GCC diagnostic ignored "-Wnonnull"  // avoid warning when using B macros in memcpy, memmove and atoi - warning bc they are sometimes NULL

// entry #index in Buffer
#define B(type, buf, index) ((buf).data ? ((type *)(&(buf).data[(index) * sizeof(type)])) : NULL)
#define Bc(buf, index)      B(char,     (buf), (index))
#define B8(buf, index)      B(uint8_t,  (buf), (index))
#define B16(buf, index)     B(uint16_t, (buf), (index))
#define B32(buf, index)     B(uint32_t, (buf), (index))
#define B64(buf, index)     B(uint64_t, (buf), (index))
#define Btxt(index)         Bc (vb->txt_data, index)
#define Ltxt                (vb->txt_data.len32)
#define Stxt                (vb->txt_data.size)
#define Rtxt                (vb->txt_data.len < Stxt ? (Stxt - vb->txt_data.len) : 0) // remaining (non-negative)

// first entry in Buffer
#define B1ST(type, buf)     ((type     *)(buf).data)
#define B1STc(buf)          ((char     *)(buf).data)
#define B1ST8(buf)          ((uint8_t  *)(buf).data)
#define B1ST16(buf)         ((uint16_t *)(buf).data)
#define B1ST32(buf)         ((uint32_t *)(buf).data)
#define B1ST40(buf)         ((uint40_t *)(buf).data)
#define B1ST64(buf)         ((uint64_t *)(buf).data)
#define B1STtxt             (vb->txt_data.data)

// last entry in Buffer
#define BLST(type, buf)     ((buf).data ? ((type *)(&(buf).data[((buf).len-1) * sizeof(type)])) : NULL)
#define BLSTc(buf)          BLST(char,    (buf))
#define BLST8(buf)          BLST(uint8_t, (buf))
#define BLST16(buf)         BLST(uint16_t,(buf))
#define BLST32(buf)         BLST(uint32_t,(buf))
#define BLST64(buf)         BLST(uint64_t,(buf))
#define BLSTtxt             (vb->txt_data.data ? (&vb->txt_data.data[(Ltxt-1)]) : NULL)

// entry after the end of the Buffer 
#define BAFT(type, buf)     ((buf).data ? ((type *)(&(buf).data[((buf).len) * sizeof(type)])) : NULL)
#define BAFTc(buf)          BAFT(char,    (buf))
#define BAFT8(buf)          BAFT(uint8_t, (buf))
#define BAFT16(buf)         BAFT(uint16_t,(buf))
#define BAFT32(buf)         BAFT(uint32_t,(buf))
#define BAFT64(buf)         BAFT(uint64_t,(buf))
#define BAFTtxt             (&vb->txt_data.data[Ltxt])

#define for_buf(element_type, iterator, buf)  \
    for (element_type *iterator=B1ST(element_type, (buf)), *after_##iterator=BAFT(element_type, (buf)); iterator && iterator < after_##iterator; iterator++)

#define for_buf_back(element_type, iterator, buf)  \
    for (element_type *iterator=BLST(element_type, (buf)), *first_##iterator=B1ST(element_type, (buf)); iterator && iterator >= first_##iterator; iterator--)

// loop with two concurrent iterators "iter_p" (pointer to element_type) and "iter_i" (32bit) 
#define for_buf2(element_type, iter_p, iter_i, buf) \
    for (uint32_t iter_i=0, after_##iter_i=(buf).len32; iter_i < after_##iter_i;)  \
        for (element_type *iter_p=B1ST(element_type, (buf)), *after_##iter_p=BAFT(element_type, (buf)); iter_p && iter_p < after_##iter_p; iter_p++, iter_i++)

#define for_buf2_back(element_type, iter_p, iter_i, buf) \
    for (int32_t iter_i=(buf).len32-1; iter_i >= 0;)  \
        for (element_type *iter_p=BLST(element_type, (buf)), *first_##iterator=B1ST(element_type, (buf)); iter_p && iter_p >= first_##iterator; iter_p--, iter_i--)

#define for_buf_tandem(element_type1, iterator1, buf1, element_type2, iterator2, buf2)  \
    ASSERT ((buf1).len32 == (buf2).len32, "expecting %s.len=%u == %s.len=%u", (buf1).name, (buf1).len32, (buf2).name, (buf2).len32);\
    element_type1 *iterator1=B1ST(element_type1, (buf1)), *after_##iterator1=BAFT(element_type1, (buf1)); \
    element_type2 *iterator2=B1ST(element_type2, (buf2)); \
    for (; iterator1 && iterator1 < after_##iterator1; iterator1++, iterator2++)

#define for_buf_tandem_back(element_type1, iterator1, buf1, element_type2, iterator2, buf2)  \
    ASSERT ((buf1).len32 == (buf2).len32, "expecting %s.len=%u == %s.len=%u", (buf1).name, (buf1).len32, (buf2).name, (buf2).len32);\
    element_type1 *iterator1=BLST(element_type1, (buf1)), *first_##iterator1=B1ST(element_type1, (buf1)); \
    element_type2 *iterator2=BLST(element_type2, (buf2)); \
    for (; iterator1 && iterator1 >= first_##iterator1; iterator1--, iterator2--)

// remove entries from buffer that fail to meet the condition
#define buf_remove_items_except_(type, buf, keep_predicate) \
    type *new_e = B1ST(type, (buf));                        \
    for_buf (type, ent, (buf)) {                            \
        if (keep_predicate) {                               \
            if (new_e != ent) *new_e = *ent;                \
            new_e++;                                        \
        }                                                   \
    }                                                       \
    (buf).len = BNUM((buf), new_e);

#define buf_remove_items_except(type, buf, field, must_be)  \
    buf_remove_items_except_(type, (buf), (ent->field == (must_be)))

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
#define BNXTf(buf)          BNXT(float, (buf))
#define BNUM(buf, ent)      ((int32_t)((((char*)(ent)) - ((buf).data)) / (int32_t)sizeof (*(ent)))) // signed integer
#define BNUM64(buf, ent)    ((int64_t)((((char*)(ent)) - ((buf).data)) / (int64_t)sizeof (*(ent)))) // signed integer
#define BNUMtxt(ent)        BNUM(vb->txt_data, (ent))
#define BREMAINS(buf, ent)  ((buf).len - BNUM ((buf),(ent)))
#define BFREE8(buf)         ((buf).size - (buf).len) // free space, for char or uint8_t buffers
#define BIS1ST(buf,ent)     (BNUM((buf),(ent)) == 0)
#define BISLST(buf,ent)     (BNUM((buf),(ent)) == (buf).len - 1)
#define BISBEFORE(buf,ent)  (BNUM((buf),(ent)) <= -1)
#define BISAFT(buf,ent)     (BNUM((buf),(ent)) >= (buf).len)
#define BISVALID(buf, ent)  (!BISBEFORE((buf),(ent)) && !BISAFT((buf),(ent)))
#define ASSBISVALID(buf,ent) ASSERT (BISVALID((buf), (ent)), #ent "=%p is not in its Buffer " #buf, (ent))

#define IN_BUF(ent,buf) (BNUM64((buf),(ent)) >= 0 && BNUM64((buf),(ent)) < (buf).len) // is item contained in buffer

#define INSERTAFTER(type, buf, index) ({                \
    buf_alloc ((buf).vb, &(buf), 1, 0, type, 1, NULL);  \
    memmove (B(type, (buf), (index)+2), B(type, (buf), (index)+1), ((buf).len - ((index)+1)) * sizeof(type)); \
    (buf).len++;                                        \
    B(type, (buf), (index)+1);                          \
})

extern void buf_remove_do (BufferP buf, unsigned sizeof_item, uint64_t remove_start, uint64_t remove_len, FUNCLINE);
#define buf_remove(buf, type, remove_start, remove_len) buf_remove_do (&(buf), sizeof(type), (remove_start), (remove_len), __FUNCLINE)

extern void buf_insert_do (VBlockP vb, BufferP buf, unsigned width, uint64_t insert_at, const void *new_data, uint64_t new_data_len, rom name, FUNCLINE);

#define buf_has_space(buf, new_len) ((buf)->data && (buf)->len + (new_len) <= (buf)->size)

static inline void buf_add_do (BufferP buf, STRp(data))
{
    memcpy (&buf->data[buf->len], data, data_len);
    buf->len += data_len;
}

extern void buf_add (BufferP buf, STRp(data));

#define buf_add_more(vb, buf, new_data, new_data_len, name)                             \
    buf_insert_do ((VBlockP)(vb), (buf), 1, (buf)->len, (new_data), (new_data_len), (name), __FUNCLINE)

#define buf_append(vb, buf, type, new_data, new_data_len, name)                         \
    buf_insert_do ((VBlockP)(vb), &(buf), sizeof(type), (buf).len, (new_data), (new_data_len), (name), __FUNCLINE)

#define buf_append_one(buf, item) ({                                                    \
    typeof(item) item_ = (item); /* evaluate once */                                    \
    buf_insert_do (NULL, &(buf), sizeof(typeof(item)), (buf).len, &item_, 1, NULL, __FUNCLINE);\
    BLST(typeof(item), (buf));                                                          \
})

#define buf_append_buf(dst_vb,dst_buf,src_buf,type,name) ({ \
    buf_alloc ((dst_vb) ? (VBlockP)(dst_vb) : (dst_buf)->vb, (dst_buf), (src_buf)->len, 0, type, CTX_GROWTH, (name)); \
    memcpy ((char *)BAFT(type, *(dst_buf)), (src_buf)->data, (src_buf)->len * sizeof (type));   \
    (dst_buf)->len += (src_buf)->len; })

#define buf_insert(vb, buf, type, insert_at, new_data, new_data_len, name) \
    buf_insert_do ((VBlockP)(vb), &(buf), sizeof(type), (insert_at), (new_data), (new_data_len), (name), __FUNCLINE)

#define buf_add_moreC(vb_, buf, literal_str, name) buf_add_more ((VBlockP)(vb_), (buf), literal_str, sizeof literal_str-1, (name))
#define buf_add_moreS(vb_, buf, str, name) buf_add_more ((VBlockP)(vb_), (buf), str, str##_len, (name))

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

#define buf_add_int(vb, buf, n) ({ buf_alloc ((VBlockP)(vb), &(buf), 1, CTX_GROWTH, typeof(n), 0, NULL); BNXT(typeof(n), (buf)) = n; })

#define buf_add_vector(vb, type, dst, src, name)  \
({  if ((src).len) { \
        if (!(dst).len) \
            buf_copy ((vb), &(dst), &(src), sizeof (type), 0, 0, name); \
        else \
            for (uint64_t i=0; i < (src).len; i++) \
                *B(type, (dst), i) += *B(type, (src), i); \
    } })

extern void buf_swap (BufferP buf1, BufferP buf2);

#define BUFPRINTF_MAX_LEN 5000
#define bufprintf(vb, buf, format, ...)  ({ char __s[BUFPRINTF_MAX_LEN+2]; \
                                            int __s_len = snprintf (__s, sizeof(__s), (format), __VA_ARGS__); \
                                            ASSERT (__s_len <= BUFPRINTF_MAX_LEN , "bufprintf: String too long - stack overflow: len=%u", __s_len); \
                                            buf_append_string ((VBlockP)(vb), (buf), __s); })
#define bufprint0 buf_append_string  // note: the string isn't a printf format, so escaping works differently (eg % doesn't need to be %%)

static inline bool buf_issame (BufferP a, BufferP b, unsigned width)
{
    return a->len == b->len && !memcmp (a->data, b->data, a->len * width);
}

#define buf_set(buf_p,value) ({ if ((buf_p)->data) memset ((buf_p)->data, value, (buf_p)->size); })
#define buf_zero(buf_p) buf_set(buf_p, 0)

extern void buf_copy_do (VBlockP dst_vb, BufferP dst, ConstBufferP src, uint64_t bytes_per_entry,
                         uint64_t src_start_entry, uint64_t max_entries, // if 0 copies the entire buffer
                         FUNCLINE,
                         rom name);
#define buf_copy(dst_vb,dst,src,type,src_start_entry,max_entries,dst_name) \
  buf_copy_do ((VBlockP)(dst_vb),(dst),(src),sizeof(type),(src_start_entry),(max_entries),__FUNCLINE,(dst_name))

typedef bool (*TxtIteratorCallback)(rom line, unsigned line_len, void *cb_param1, void *cb_param2, unsigned cb_param3);
extern char *buf_foreach_line (BufferP buf, bool reverse, TxtIteratorCallback callback, void *cb_param1, void *cb_param2, unsigned cb_param3, int64_t *line_len);

extern void buf_print (BufferP buf, bool add_newline);

extern bool buf_dump_to_file (rom filename, ConstBufferP buf, unsigned buf_word_width, bool including_control_region, bool no_dirs, bool verbose, bool do_gzip);

//-----------------
// bits
//-----------------

// allocate bit array and set nbits
typedef enum { CLEAR=0, SET=1, NOINIT=2 } BitsInitType;
extern BitsP buf_alloc_bits_exact_do (VBlockP vb, BufferP buf, uint64_t exact_bits, BitsInitType init_to, float grow_at_least_factor, rom name, FUNCLINE);
#define buf_alloc_bits_exact(vb, buf, exact_nbits, init_to, grow_at_least_factor, name) \
    buf_alloc_bits_exact_do ((VBlockP)(vb), (buf), (exact_nbits), (init_to), (grow_at_least_factor), (name), __FUNCLINE)

extern BitsP buf_alloc_bits_do (VBlockP vb, BufferP buf, uint64_t nbits, uint64_t at_east_bits, BitsInitType init_to, float grow_at_least_factor, rom name, FUNCLINE);
#define buf_alloc_bits(vb, buf, more_bits, at_least_bits, init_to, grow_at_least_factor, name) \
    buf_alloc_bits_do ((VBlockP)(vb), (buf), (more_bits), (at_least_bits), (init_to), (grow_at_least_factor), (name), __FUNCLINE)

extern BitsP buf_overlay_bits_do (VBlockP vb, BufferP overlaid_buf, BufferP regular_buf, uint64_t start_byte_in_regular_buf, uint64_t nbits, FUNCLINE, rom name);
#define buf_overlay_bits(vb, overlaid_buf, regular_buf, start_byte_in_regular_buf, nbits, name) \
    buf_overlay_bits_do((VBlockP)(vb), (overlaid_buf), (regular_buf), (start_byte_in_regular_buf), (nbits), __FUNCLINE, (name))

extern uint64_t buf_extend_bits (BufferP buf, int64_t num_new_bits);
extern void buf_add_bit (BufferP buf, int64_t new_bit);
extern BitsP buf_zfile_buf_to_bits (BufferP buf, uint64_t nbits);
#define buf_add_set_bit(buf)   buf_add_bit (buf, 1)
#define buf_add_clear_bit(buf) buf_add_bit (buf, 0)

//-----------------
// Endianity stuff
//-----------------

typedef void BgEnBufFunc (BufferP buf, LocalType *lt); 

typedef BgEnBufFunc (*BgEnBuf);

extern BgEnBufFunc BGEN_u8_buf, BGEN_u16_buf, BGEN_u32_buf, BGEN_u64_buf, 
                   BGEN_transpose_u8_buf, BGEN_transpose_u16_buf, BGEN_transpose_u32_buf,
                   BGEN_ptranspose_u8_buf, BGEN_ptranspose_u16_buf, BGEN_ptranspose_u32_buf,
                   BGEN_deinterlace_d8_buf, BGEN_deinterlace_d16_buf, BGEN_deinterlace_d32_buf, BGEN_deinterlace_d64_buf,
                   interlace_d8_buf, BGEN_interlace_d16_buf, BGEN_interlace_d32_buf, BGEN_interlace_d64_buf,
                   LTEN_u16_buf, LTEN_u32_buf, LTEN_u64_buf, LTEN_interlace_d16_buf, LTEN_interlace_d32_buf, LTEN_interlace_d64_buf;
                   
