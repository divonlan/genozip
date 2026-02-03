// ------------------------------------------------------------------
//   buffer.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include <stdint.h>
#include <pthread.h>
#include "strings.h"
#include "buf_list.h"
#include "endianness.h"

// Notes: 1 byte. BUF_UNALLOCATED must be 0. all values must be identical to 
#define BUF_TYPE_BITS 3
typedef enum          { BUF_UNALLOCATED=0, BUF_REGULAR,   BITS_STANDALONE, BUF_SHM, BUF_DISOWNED, BUF_NUM_TYPES } BufferType; 
#define BUFTYPE_NAMES {    "UNALLOCATED",     "REGULAR", "BITS_STANDALONE",   "SHM",   "DISOWNED"                                                   }

typedef struct { 
    bool lock;
    uint16_t link_count;       // # of buffers pointing to this spinlock. >= Buffer user count, bc if buffer splits, spinlock doesn't.
} BufferSpinlock, *BufferSpinlockP; // see internal-docs/overlay-logic.txt

typedef struct Buffer {        // 72 bytes
    //------------------------------------------------------------------------------------------------------
    // these 4 fields, must be first, and can be overlayed with a Bits - their order and size is identical to Bits
    union {
        char *data;            // memory+8 when initially allocated or NULL if not, but allowed to have a bigger offset vs memory (eg partial overlay)
        uint64_t *words;       // for Bits
    };
    union {                    // the "parameter" field is for discretionary use by the caller. some options to use the parameter are provided.
        int64_t param;    
        int64_t next;     
        int64_t count;    
        uint64_t n_cols;       // for matrices
        uint64_t nbits;        // for Bits
        int32_t newest_index;  // for lookback
        uint32_t prm32[2];
        uint16_t prm16[4];
        uint8_t  prm8 [8];
        void *pointer;
        struct { int32_t uncomp_len, comp_len;   }; // (signed) used for compressed buffer: uncomp_len is the length of the uncompressing comp_len data from the buffer
        struct { uint32_t consumed_by_prev_vb;      // vb->gz_blocks: bytes of the first BGZF block consumed by the prev VB or txt_header
                 uint32_t current_bb_i;          }; // index into vb->gz_blocks of first bgzf block of current line
        struct { uint32_t next_index, next_line; }; // used by z_file->gencomp_vb_lines for building recon plan
        CompIType prev_comp_i;                      // used by vb->vb_plan
    };
    union {
        uint64_t len;          // used by the buffer user according to its internal logic. not modified by malloc/realloc, zeroed by buf_free (in Bits - nwords)
        uint64_t nwords;       // for Bits
        uint32_t gap_index;    // for lookback
        struct {
#ifdef __LITTLE_ENDIAN__
            uint32_t len32, len32hi; // we use len32 instead of len if easier, in cases where we are certain len is smaller than 4B
#else                
            uint32_t len32hi, len32;
#endif  
        };
    };

    BufferType type      : BUF_TYPE_BITS;  // BufferType
    //------------------------------------------------------------------------------------------------------
    uint64_t can_be_big  : 1;  // do not display warning if buffer grows very big
    uint64_t promiscuous : 1;  // used only in evb buffers: if true, a compute thread may allocate the buffer (not just the main thread as with usual evb buffers), 
                               // but it needs to be first buf_set_promiscuous by the main thread
    uint64_t shared      : 1;  // "memory" of this Buffer may be shared with other Buffers
    uint64_t code_line   : 12; // the allocating line number in source code file (up to 4096)
    #define BUF_SIZE_BITS (64 - BUF_TYPE_BITS - 15)
    #define MAX_BUFFER_SIZE ((1ULL<<BUF_SIZE_BITS)-1) // according to bits of "size" 
    uint64_t size        : BUF_SIZE_BITS;             // number of bytes available to the user (i.e. not including the allocated overhead). 
    VBlockP vb;                // vb that owns this buffer, and which this buffer is in its buf_list

    rom name;                  // name of allocator - used for memory debugging & statistics.
    rom func;                  // the allocating function

    char *memory;              // memory allocated to this buffer - amount is: size + CTL_SIZE to allow for UNDERFLOW, OVERFLOW and user_count)

    BufferSpinlockP spinlock;  // see internal-docs/overlay-logic.txt
} Buffer; 

#define UNDERFLOW_TRAP 0x574F4C4652444E55ULL // "UNDRFLOW" - inserted at the begining of each memory block to detected underflows
#define OVERFLOW_TRAP  0x776F6C667265766FULL // "overflow" - inserted at the end of each memory block to detected overflows

#define BUNDERFLOW_(memory) (*(uint64_t *)(memory))
#define BUNDERFLOW(buf) BUNDERFLOW_((buf)->memory)
#define BOVERFLOW_(memory,size)  (*(uint64_t *)((memory) + (size) + sizeof(uint64_t)))
#define BOVERFLOW(buf) BOVERFLOW_((buf)->memory, (buf)->size + ((buf)->data/*might be partial overlay*/ ? ((buf)->data - (buf)->memory - sizeof(uint64_t)) : 0))
#define BOLCOUNTER_(memory,data,size) (*(uint16_t *)(&BOVERFLOW_(memory, size + (data ? (data - memory - sizeof(uint64_t)) : 0)) + 1))
#define BOLCOUNTER(buf) BOLCOUNTER_((buf)->memory, (buf)->data, (buf)->size)

#define CTL_SIZE (2*sizeof (uint64_t) + sizeof(uint16_t)) // underflow, overflow and user counter

extern void buf_initialize(void);

#define buf_is_alloc(buf_p) ((buf_p)->data != NULL && (buf_p)->type != BUF_UNALLOCATED)
#define ASSERTNOTINUSE(buf)  ASSERT (!buf_is_alloc (&(buf)) && !(buf).len && !(buf).param, "expecting %s to be free, but it's not: %s", #buf, buf_desc (&(buf)).s)
#define ASSERTISALLOCED(buf) ASSERT (buf_is_alloc (&(buf)), "%s is not allocated", #buf)
#define ASSERTISEMPTY(buf)   ASSERT (buf_is_alloc (&(buf)) && !(buf).len, "expecting %s to be be allocated and empty, but it isn't: %s", #buf, buf_desc (&(buf)).s)
#define ASSERTNOTEMPTY(buf)  ASSERT ((buf).len && (buf).data, "expecting %s to be contain some data, but it doesn't: %s", #buf, buf_desc (&(buf)).s)

extern void buf_alloc_do (VBlockP vb, BufferP buf, uint64_t requested_size, float grow_at_least_factor, rom name, FUNCLINE);

static inline void buf_alloc_quick (BufferP buf, uint64_t req_size, rom name, FUNCLINE)
{
    if (__builtin_expect (!buf->data && req_size, false)) {
        buf->func      = func; 
        buf->code_line = code_line; 
        buf->data      = buf->memory + sizeof (uint64_t); 
        if (name) buf->name = name; 
        else ASSERT (buf->name, "%s:%u: no name", func, code_line); 
    }
}

#define buf_alloc_(alloc_vb, buf, more, at_least, width, grow_at_least_factor, name, func, code_line) ({\
    uint64_t new_more = (more); /* avoid evaluating twice */                                \
    uint64_t if_more = new_more ? ((buf)->len + new_more) : 0; /* in units of type */       \
    uint64_t new_req_size = MAX_((uint64_t)(at_least), if_more) * width; /* make copy to allow ++ */  \
    if (__builtin_expect(new_req_size <= (buf)->size, true))                                \
        buf_alloc_quick ((buf), new_req_size, (name), func, code_line);                     \
    else                                                                                    \
        buf_alloc_do (((alloc_vb) ? ((VBlockP)alloc_vb) : (buf)->vb), (buf), new_req_size, (grow_at_least_factor), (name), func, code_line); \
})

#define buf_alloc(alloc_vb, buf, more, at_least, type, grow_at_least_factor, name) \
    buf_alloc_((alloc_vb), (buf), (more), (at_least), sizeof(type), (grow_at_least_factor), (name), __FUNCTION__, __LINE__)

#define buf_alloc_zero(vb, buf, more, at_least, element_type, grow_at_least_factor,name) ({ \
    uint64_t size_before = (buf)->data ? (buf)->size : 0; /* always zero the whole buffer in an initial allocation */ \
    buf_alloc((vb), (buf), (more), (at_least), element_type, (grow_at_least_factor), (name)); \
    if ((buf)->data && (buf)->size > size_before) memset (&(buf)->data[size_before], 0, (buf)->size - size_before); })

#define buf_alloc_255(vb, buf, more, at_least, element_type, grow_at_least_factor,name) ({ \
    uint64_t size_before = (buf)->data ? (buf)->size : 0; /* always zero the whole buffer in an initial allocation */ \
    buf_alloc((vb), (buf), (more), (at_least), element_type, (grow_at_least_factor), (name)); \
    if ((buf)->data && (buf)->size > size_before) memset (&(buf)->data[size_before], 255, (buf)->size - size_before); })

// alloc a set amount of bytes, and set buf.len
#define buf_alloc_exact(alloc_vb, buf, exact_len, type, name) ({  \
    buf_alloc((alloc_vb), &(buf), 0, (exact_len), type, 1, name); \
    (buf).len = (exact_len); })

// note: the entire buffer is zeroed, not just the added bytes
#define buf_alloc_exact_zero(alloc_vb, buf, exact_len, type, name) ({   \
    buf_alloc_exact (alloc_vb, buf, exact_len, type, name); \
    memset ((buf).data, 0, (exact_len) * sizeof(type)); })

// note: the entire buffer is set to 255, not just the added bytes
#define buf_alloc_exact_255(alloc_vb, buf, exact_len, type, name) ({   \
    buf_alloc_exact (alloc_vb, buf, exact_len, type, name); \
    memset ((buf).data, 255, (exact_len) * sizeof(type)); })

// allocates exactly the requested amount and sets let, and declares ARRAY
#define ARRAY_alloc(element_type, array_name, array_len, init_zero, buf, alloc_vb, buf_name) \
    buf_alloc_exact (((alloc_vb) ? ((VBlockP)alloc_vb) : (buf).vb), (buf), (array_len), element_type, (buf_name)); \
    if (init_zero) memset ((buf).data, 0, (buf).len * sizeof(element_type)); /* resets the entire buffer, not just newly allocated memory */ \
    element_type *array_name = ((element_type *)((buf).data)); \
    const uint64_t array_name##_len __attribute__((unused)) = (buf).len; // read-only copy of len 

extern void buf_attach_to_shm_do (VBlockP vb, BufferP buf, void *data, uint64_t size, uint64_t start, FUNCLINE, rom name);

#define buf_attach_to_shm(vb, buf, data, size, name)                                     \
    buf_attach_to_shm_do ((VBlockP)(vb), (buf), (data), (size), 0, __FUNCLINE, (name))

#define buf_attach_bits_to_shm(vb, buf, data, n_bits, name)                              \
    ({ (buf)->nwords = roundup_bits2words64 (n_bits);                                    \
       (buf)->nbits  = (n_bits);                                                         \
       buf_attach_to_shm_do ((VBlockP)(vb), (buf), (data), (buf)->nwords * sizeof(uint64_t), 0, __FUNCLINE, (name));  })                                                                              \
       
extern void buf_free_do (BufferP buf, FUNCLINE);
#define buf_free(buf) buf_free_do (&(buf), __FUNCLINE)

extern void buf_destroy_do_do (BufListEnt *ent, FUNCLINE);
extern void buf_destroy_do (BufferP buf, FUNCLINE);
#define buf_destroy(buf) buf_destroy_do (&(buf), __FUNCLINE)

#define buf_is_large_enough(buf_p, requested_size) (buf_is_alloc ((buf_p)) && (buf_p)->size >= requested_size)

extern void buf_move_do (VBlockP vb, BufferP dst_buf, rom dst_name, BufferP src_buf, FUNCLINE);
#define buf_move(vb, dst_buf, dst_name, src_buf) buf_move_do ((VBlockP)(vb), &(dst_buf), (dst_name), &(src_buf), __FUNCLINE)

extern void buf_grab_do (VBlockP dst_vb, BufferP dst_buf, rom dst_name, BufferP src_buf, FUNCLINE);
#define buf_grab(dst_vb, dst_buf, dst_name, src_buf) buf_grab_do ((VBlockP)(dst_vb), &(dst_buf), (dst_name), &(src_buf), __FUNCLINE)

extern void buf_disown_do (VBlockP vb, BufferP src_buf, BufferP dst_buf, bool make_a_copy, FUNCLINE);
#define buf_disown(vb, src_buf, dst_buf, make_a_copy) buf_disown_do ((vb), &(src_buf), &(dst_buf), (make_a_copy), __FUNCLINE)

extern void buf_extract_data_do (BufferP buf, char **data_p, uint64_t *len_p, uint32_t *len32_p, char **memory_p, FUNCLINE);
#define buf_extract_data(buf, data_p, len_p, len32_p, memory_p) buf_extract_data_do (&(buf), (data_p), (len_p), (len32_p), (memory_p), __FUNCLINE)

extern void buf_verify_do (ConstBufferP buf, rom msg, FUNCLINE);
#define buf_verify(buf, msg) buf_verify_do (&(buf), (msg), __FUNCLINE)

extern void buf_trim_do (BufferP buf, uint64_t size, FUNCLINE);
#define buf_trim(buf, type) buf_trim_do (&(buf), (buf).len * sizeof(type), __FUNCLINE)

typedef struct {
    rom name;
    uint64_t bytes; 
    unsigned buffers;
} MemStats;

extern void buf_set_promiscuous_do (VBlockP vb, BufferP buf, rom buf_name, FUNCLINE);
#define buf_set_promiscuous(buf, buf_name) buf_set_promiscuous_do (evb, (buf), (buf_name), __FUNCLINE)

extern void buf_low_level_free (void *p, FUNCLINE);
#define FREE(p) ({ if (p) { buf_low_level_free (((void*)(p)), __FUNCLINE); (p)=NULL; } })

extern void *buf_low_level_malloc (size_t size, bool zero, FUNCLINE);
#define MALLOC(size) buf_low_level_malloc (size, false, __FUNCLINE)
#define CALLOC(size) buf_low_level_malloc (size, true,  __FUNCLINE)

extern void *buf_low_level_realloc (void *p, size_t size, rom name, FUNCLINE);
#define REALLOC(p,size,name) if (!(*(p) = buf_low_level_realloc (*(p), (size), (name), __FUNCLINE))) ABORT0 ("REALLOC failed")

extern void buf_low_level_release_memory_back_to_kernel (void);

extern void buf_set_shared (BufferP buf);
extern void buf_remove_spinlock (BufferP buf);

extern void buf_overlay_do (VBlockP vb, BufferP top_buf, BufferP bottom_buf, uint64_t start_in_bottom, bool copy_len, FUNCLINE, rom name);
#define buf_overlay(vb, top_buf, bottom_buf, name) \
    buf_overlay_do((VBlockP)(vb), (top_buf), (bottom_buf), 0, true, __FUNCLINE, (name)) 

// note: partially OVERLAY buffers MUST be freed before their bottom_buf 
#define buf_overlay_partial(vb, top_buf, bottom_buf, start_in_bottom, name) \
    buf_overlay_do((VBlockP)(vb), (top_buf), (bottom_buf), (start_in_bottom), false, __FUNCLINE, (name)) 

extern uint64_t buf_mem_size (ConstBufferP buf);

//--------------------------
// thread synchronization
//--------------------------

extern void buf_init_lock (BufferP buf);

#define buf_lock_if(buf, cond) \
    BufferSpinlockP spinlock = (cond) ? (buf)->spinlock : NULL; \
    ASSERT (!(cond) || spinlock, "spinlock not initialized for %s", buf_desc(buf).s); \
    if (spinlock) while (({ bool expected = (bool)false; !__atomic_compare_exchange_n (&spinlock->lock, &expected, (bool)true, false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE); })); /* spinlock */ 

#define buf_lock(buf) buf_lock_if ((buf), true)
#define buf_lock_(buf) rom func __attribute__((unused)) = __FUNCTION__; buf_lock_if ((buf), true)

#define buf_unlock ({ if (spinlock) { __atomic_clear (&spinlock->lock, __ATOMIC_RELEASE); \
                                      spinlock = NULL; \
                                      /* printf ("unlocked %s\n", func); */}; })

extern BufferSpinlockP buf_lock_promiscuous (ConstBufferP buf, FUNCLINE);

static inline uint16_t buf_user_count (ConstBufferP buf) 
{
    return buf->memory ? BOLCOUNTER(buf) : 0;
}

typedef struct { char s[300]; } BufDescType;
extern const BufDescType buf_desc (ConstBufferP buf);

extern rom buf_type_name (ConstBufferP buf);
