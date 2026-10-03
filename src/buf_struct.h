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

// Notes: 1 byte. BUF_UNALLOCATED must be 0. 
typedef enum          { BUF_UNALLOCATED=0, BUF_REGULAR,   BITS_OVERLAY,   BITS_SELF_ALLOC, BUF_SHM, BUF_DISOWNED, BUF_SUPERIMPOSED, BUF_NUM_TYPES } BufferType; 
#define BUFTYPE_NAMES {    "UNALLOCATED",     "REGULAR", "BITS_OVERLAY", "BITS_SELF_ALLOC",   "SHM",   "DISOWNED"    "SUPERIMPOSED"                                                   }

typedef struct { 
    bool lock;
    uint16_t link_count;       // # of buffers pointing to this spinlock. >= Buffer user count, bc if buffer splits, spinlock doesn't.
} BufferSpinlock, *BufferSpinlockP; // see internal-docs/overlay-logic.txt

// note: we would like to alignas(32) Buffer, but a gcc bug causes Buffer to be sometimes 
// incorrectly unaligned when on the stack, but nevertheless SMID instructions are used when copying 
// a local Buffer (e.g buf = (Buffer){...}), attempting aligned memory load, causing segfault
typedef struct Buffer {        // 64 bytes: the first 32B are frequently accessed, 2nd 32B much less frequently
    // the two most commonly accessed field (data, len) are adjacent
    union { // 8 bytes
        char *data;            // memory+8 when initially allocated or NULL if not, but allowed to have a bigger offset vs memory (eg partial overlay)
        uint64_t *words;       // for Bits
    };
    union { // 8 bytes
        uint64_t len;          // used by the buffer user according to its internal logic. not modified by malloc/realloc, zeroed by buf_free (in Bits - nwords)
        uint64_t nwords;       // for Bits
        uint32_t gap_index;    // for lookback
        struct {
            ℒ𝒾𝓉ℰ (uint32_t len32;) // we use len32 instead of len if easier, in cases where we are certain len is smaller than 4B
            uint32_t len32hi;
            ℬ𝒾ℊℰ (uint32_t len32;) 
        };
    };

    union { /*8 bytes */       // the "parameter" field is for discretionary use by the caller. some options to use the parameter are provided.
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
        
        struct QBits { // used by gencomp: queue[gct].gc_txts[buf_i]
            uint32_t num_lines : 27; // VB size is limited to 1GB and BAM has a theorical minium of 38 bytes per line, and SAM 22.
            uint32_t comp_i    : 2;
            uint32_t unused    : 3;
            uint16_t next;           // queue: towards head ; stack: down stack
            uint16_t prev;           // queue: towards tail ; stack: always END_OF_LIST
        } qbits;
    };

    // --- 8 bytes starts here ---
    #define MAX_BUFFER_SIZE ((1ULL<<46)-1) // according to bits of "size" 
    uint64_t size        : 46; // number of bytes available to the user (i.e. not including the allocated overhead). 
    BufferType type      : 3;  
    uint64_t can_be_big  : 1;  // do not display warning if buffer grows very big
    uint64_t promiscuous : 1;  // used only in evb buffers: if true, a compute thread may allocate the buffer (not just the main thread as with usual evb buffers), 
                               // but it needs to be first buf_set_promiscuous by the main thread
    uint64_t shared      : 1;  // "memory" of this Buffer may be shared with other Buffers
    uint64_t code_line   : 12; // the allocating line number in source code file (up to 4096)
    // --- 8 bytes ends here ---

    char *memory;              // memory allocated to this buffer - amount is: size + CTL_SIZE to allow for UNDERFLOW, OVERFLOW and user_count)

    VBlockP vb;                // vb that owns this buffer, and which this buffer is in its buf_list

    // store these two static strings as relative (ר) to string_anchor (int32_t instead of a 64bit pointer)
    Pointeר nameר;             // buffer name - used for memory debugging & statistics.
    Pointeר funcר;             // the allocating function 

    BufferSpinlockP spinlock[1];// Used for overlay top/bottom buffers, promiscuous buffers and buffer_list. see internal-docs/overlay-logic.txt
} Buffer; 

#define UNDERFLOW_TRAP 0x574F4C4652444E55ULL // "UNDRFLOW" - inserted at the begining of each memory block to detected underflows
#define OVERFLOW_TRAP  0x776F6C667265766FULL // "overflow" - inserted at the end of each memory block to detected overflows

// extra padding between UNDERFLOW and data: partial overlay, or 64B-aligned "data" pointer.
#define BUNDERFLOW_(memory) (*(uint64_t *)(memory)) 
#define BUNDERFLOW(buf) BUNDERFLOW_((buf)->memory)
#define BOVERFLOW_(memory,data,size)  (*(uint64_t *)((memory) + sizeof(uint64_t)/*underflow*/ + (data ? (data - memory - sizeof(uint64_t)) : 0)/*padding: partial overlay, or 64B alignment*/ + (size)))
#define BOVERFLOW(buf) BOVERFLOW_((buf)->memory, (buf)->data, (buf)->size)
#define BOLCOUNTER_(memory,data,size) (*(uint16_t *)(&BOVERFLOW_(memory, data, size) + 1))
#define BOLCOUNTER(buf) BOLCOUNTER_((buf)->memory, (buf)->data, (buf)->size)

#define CTL_SIZE (2*sizeof (uint64_t) + sizeof(uint16_t)) // underflow, overflow and user counter

extern void buf_initialize(void);

#define buf_is_alloc(buf_p) ((buf_p)->data != NULL && (buf_p)->type != BUF_UNALLOCATED)
#define ASSERTNOTINUSE(buf)  ASSERT (!buf_is_alloc (&(buf)) && !(buf).len && !(buf).param, "expecting %s to be free, but it's not: %s", #buf, buf_desc (&(buf)).s)
#define ASSERTISALLOCED(buf) ASSERT (buf_is_alloc (&(buf)), "%s is not allocated", #buf)
#define ASSERTISEMPTY(buf)   ASSERT (buf_is_alloc (&(buf)) && !(buf).len, "expecting %s to be be allocated and empty, but it isn't: %s", #buf, buf_desc (&(buf)).s)
#define ASSERTNOTEMPTY(buf)  ASSERT ((buf).len && (buf).data, "expecting %s to be contain some data, but it doesn't: %s", #buf, buf_desc (&(buf)).s)

extern void buf_alloc_do (VBlockP vb, BufferP buf, uint64_t requested_size, float grow_at_least_factor, bool aligned_64B, rom name, Caller caller);

static inline void buf_alloc_quick (BufferP buf, uint64_t req_size, rom name, Caller caller)
{
    if (__builtin_expect (!buf->data && req_size, false)) {
        buf->funcר     = caller.funcר; 
        buf->code_line = caller.code_line; 
        buf->data      = buf->memory + sizeof (uint64_t); 
        if (name) buf->nameר = ר(name); 
        else ASSERT (buf->nameר, "%s:%u: no name", CALLERf); 
    }
}

#define buf_alloc_(alloc_vb, buf, more, at_least, width, grow_at_least_factor, name, caller) ({\
    uint64_t new_more = (more); /* avoid evaluating twice */                                \
    uint64_t if_more = new_more ? ((buf)->len + new_more) : 0; /* in units of type */       \
    uint64_t new_req_size = MAX_((uint64_t)(at_least), if_more) * width; /* make copy to allow ++ */  \
    if (__builtin_expect(new_req_size <= (buf)->size, true))                                \
        buf_alloc_quick ((buf), new_req_size, (name), caller);                              \
    else                                                                                    \
        buf_alloc_do ((VBlockP)alloc_vb, (buf), new_req_size, (grow_at_least_factor), false, (name), caller); \
})

#define buf_alloc(alloc_vb, buf, more, at_least, type, grow_at_least_factor, name) \
    buf_alloc_((alloc_vb), (buf), (more), (at_least), sizeof(type), (grow_at_least_factor), (name), THIS_CODE_LINE)

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
    if ((buf).data) memset ((buf).data, 0, (exact_len) * sizeof(type)); }) // if protects against undefined behaviour (per C spec) if data=NULL and len=0

// note: the entire buffer is set to 255, not just the added bytes
#define buf_alloc_exact_255(alloc_vb, buf, exact_len, type, name) ({   \
    buf_alloc_exact (alloc_vb, buf, exact_len, type, name); \
    if ((buf).data) memset ((buf).data, 255, (exact_len) * sizeof(type)); })

static inline void buf_aligned_64B (BufferP buf)
{
    uint64_t padding = ROUNDUP64 ((uintptr_t)buf->data) - (uintptr_t)buf->data;
    buf->data += padding;
    buf->size -= padding;
}

static inline void buf_alloc_quick_aligned_64B (BufferP buf, uint64_t req_size, rom name, Caller caller)
{
    if (__builtin_expect (!buf->data && req_size, false)) {
        buf->funcר     = caller.funcר; 
        buf->code_line = caller.code_line; 
        buf->data      = buf->memory + sizeof (uint64_t); 

        // this will be the same padding as when the buffer was originally alloced,
        // so "size" refers to memory available with this padding
        buf_aligned_64B (buf);

        if (name) buf->nameר = ר(name); 
        else ASSERT (buf->nameר, "%s:%u: no name", CALLERf); 
    }
}

// buf->data is guaranteed to be 64B-aligned. buf->size is NOT guaranteed to be a multiple of 64.
#define buf_alloc_aligned_64B_(alloc_vb, buf, more, at_least, width, grow_at_least_factor, name, caller) ({ \
    uint64_t new_more = (more); /* avoid evaluating twice */                                \
    uint64_t if_more = new_more ? ((buf)->len + new_more) : 0; /* in units of type */       \
    uint64_t new_req_size = MAX_((uint64_t)(at_least), if_more) * width; /* make copy to allow ++ */  \
    if (__builtin_expect(new_req_size <= (buf)->size, true))                                \
        buf_alloc_quick_aligned_64B ((buf), new_req_size, (name), caller);                  \
    else                                                                                    \
        buf_alloc_do ((VBlockP)alloc_vb, (buf), new_req_size, (grow_at_least_factor), true, (name), caller); \
})

// note: cannot realloc an unaligned buffer into aligned. must buf_free first.
#define buf_alloc_aligned_64B(alloc_vb, buf, more, at_least, type, grow_at_least_factor, name) \
    buf_alloc_aligned_64B_((alloc_vb), (buf), (more), (at_least), sizeof(type), (grow_at_least_factor), (name), THIS_CODE_LINE)

// allocates exactly the requested amount and sets let, and declares ARRAY
#define ARRAY_alloc(element_type, array_name, array_len, init_zero, buf, alloc_vb, buf_name) \
    buf_alloc_exact (((alloc_vb) ? ((VBlockP)alloc_vb) : (buf).vb), (buf), (array_len), element_type, (buf_name)); \
    if (init_zero && (buf).data) memset ((buf).data, 0, (buf).len * sizeof(element_type)); /* resets the entire buffer, not just newly allocated memory */ \
    element_type *array_name = ((element_type *)((buf).data)); \
    const uint64_t array_name##_len UNUSED = (buf).len; // read-only copy of len 

#define ARRAY_alloc𐤐(element_type, array_name, array_len, init_zero, buf, alloc_vb, buf_name) \
    buf_alloc_exact (((alloc_vb) ? ((VBlockP)alloc_vb) : (buf).vb), (buf), (array_len), element_type, (buf_name)); \
    if (init_zero && (buf).data) memset ((buf).data, 0, (buf).len * sizeof(element_type)); /* resets the entire buffer, not just newly allocated memory */ \
    element_type *restrict array_name = ((element_type *)((buf).data)); \
    const uint64_t array_name##_len UNUSED = (buf).len; // read-only copy of len 

extern void buf_attach_to_shm_do (VBlockP vb, BufferP buf, void *data, uint64_t size, uint64_t start, Caller caller, rom name);

#define buf_attach_to_shm(vb, buf, data, size, name)                                     \
    buf_attach_to_shm_do ((VBlockP)(vb), (buf), (data), (size), 0, THIS_CODE_LINE, (name))

#define buf_attach_bits_to_shm(vb, buf, data, n_bits, name)                              \
    ({ (buf)->nwords = roundup_bits2words64 (n_bits);                                    \
       (buf)->nbits  = (n_bits);                                                         \
       buf_attach_to_shm_do ((VBlockP)(vb), (buf), (data), (buf)->nwords * sizeof(uint64_t), 0, THIS_CODE_LINE, (name));  })                                                                              \
       
extern void buf_free_do (BufferP buf, Caller caller);
#define buf_free(buf) buf_free_do (&(buf), THIS_CODE_LINE)

extern void buf_destroy_do_do (BufListEnt *ent, Caller caller);
extern void buf_destroy_do (BufferP buf, Caller caller);
#define buf_destroy(buf) buf_destroy_do (&(buf), THIS_CODE_LINE)

#define buf_is_large_enough(buf_p, requested_size) (buf_is_alloc ((buf_p)) && (buf_p)->size >= requested_size)

extern void buf_move_do (VBlockP vb, BufferP dst_buf, rom dst_name, BufferP src_buf, Caller caller);
#define buf_move(vb, dst_buf, dst_name, src_buf) buf_move_do ((VBlockP)(vb), &(dst_buf), (dst_name), &(src_buf), THIS_CODE_LINE)

extern void buf_grab_do (VBlockP dst_vb, BufferP dst_buf, rom dst_name, BufferP src_buf, Caller caller);
#define buf_grab(dst_vb, dst_buf, dst_name, src_buf) buf_grab_do ((VBlockP)(dst_vb), &(dst_buf), (dst_name), &(src_buf), THIS_CODE_LINE)

extern void buf_disown_do (VBlockP vb, BufferP src_buf, BufferP dst_buf, bool make_a_copy, Caller caller);
#define buf_disown(vb, src_buf, dst_buf, make_a_copy) buf_disown_do ((vb), &(src_buf), &(dst_buf), (make_a_copy), THIS_CODE_LINE)

extern void buf_extract_data_do (BufferP buf, char **data_p, uint64_t *len_p, uint32_t *len32_p, char **memory_p, Caller caller);
#define buf_extract_data(buf, data_p, len_p, len32_p, memory_p) buf_extract_data_do (&(buf), (data_p), (len_p), (len32_p), (memory_p), THIS_CODE_LINE)

extern void buf_verify_do (ConstBufferP buf, rom msg, Caller caller);
#define buf_verify(buf, msg) buf_verify_do (&(buf), (msg), THIS_CODE_LINE)

extern void buf_trim_do (BufferP buf, uint64_t size, Caller caller);
#define buf_trim(buf, type) buf_trim_do (&(buf), (buf).len * sizeof(type), THIS_CODE_LINE)

typedef struct {
    int32_t nameר;
    unsigned buffers;
    uint64_t bytes; 
} MemStats;

extern void buf_set_promiscuous_do (VBlockP vb, BufferP buf, rom buf_name, Caller caller);
#define buf_set_promiscuous(buf, buf_name) buf_set_promiscuous_do (evb, (buf), (buf_name), THIS_CODE_LINE)

extern void buf_low_level_free (void *p, bool aligned_64B, Caller caller);
#define FREE(p) ({ if (p) { buf_low_level_free (((void*)(p)), false, THIS_CODE_LINE); (p)=NULL; } })
#define FREE_ALIGN_64B(p) ({ if (p) { buf_low_level_free (((void*)(p)), true, THIS_CODE_LINE); (p)=NULL; } })

extern void *buf_low_level_malloc (size_t size, bool zero, bool aligned_64B, Caller caller);
#define MALLOC(size) buf_low_level_malloc (size, false, false, THIS_CODE_LINE)
#define CALLOC(size) buf_low_level_malloc (size, true,  false, THIS_CODE_LINE)
#define CALLOC_ALIGN_64B(size) buf_low_level_malloc (size, true, true, THIS_CODE_LINE)

extern void *buf_low_level_realloc (void *p, size_t size, rom name, Caller caller);
#define REALLOC(p,size,name) if (!(*(p) = buf_low_level_realloc (*(p), (size), (name), THIS_CODE_LINE))) ABORT0 ("REALLOC failed")

extern void return_freed_memory_to_kernel (void);

// overlaying: bottom buffer can be realloced, and overlayed buffer still using old data
// Both top and buttom buffers must be ShareableBuffer
extern void buf_set_shared (BufferP buf);
extern void buf_remove_spinlock (BufferP buf);

extern void buf_overlay_do (VBlockP vb, BufferP top_buf, BufferP bottom_buf, Caller caller, rom name);
#define buf_overlay(vb, top_buf, bottom_buf, name) \
    buf_overlay_do((VBlockP)(vb), (top_buf), (bottom_buf), THIS_CODE_LINE, (name)) 

// superimpose: 
// 1. caller guarantees that bottom buffer is not realloced or freed while having superimposed buffers
// 2. caller guarantees that top_buf is NOT in buf_list: must buf_destroy first if it might have been
// 3. top_buf is not added to the buf_list, and therefore can be an automatic variable
// 4. If Buffer is an automatic variable, it is not necessary to free or destroy it
extern void buf_superimpose_do (VBlockP vb, BufferP top_buf, BufferP bottom_buf, uint64_t start_in_bottom, Caller caller, rom name);
#define buf_superimpose(vb, top_buf, bottom_buf, start_in_bottom, name) \
    buf_superimpose_do((VBlockP)(vb), (top_buf), (bottom_buf), (start_in_bottom), THIS_CODE_LINE, (name)) 

extern uint64_t buf_mem_size (ConstBufferP buf);

//-----------------
// bits
//-----------------

// allocate bit array and set nbits
typedef enum { CLEAR=0, SET=1, NOINIT=2 } BitsInitType;
extern BitsP buf_alloc_bits_exact_do (VBlockP vb, BufferP buf, uint64_t exact_bits, BitsInitType init_to, float grow_at_least_factor, rom name, Caller caller);
#define buf_alloc_bits_exact(vb, buf, exact_nbits, init_to, grow_at_least_factor, name) \
    buf_alloc_bits_exact_do ((VBlockP)(vb), (buf), (exact_nbits), (init_to), (grow_at_least_factor), (name), THIS_CODE_LINE)

extern BitsP buf_alloc_bits_do (VBlockP vb, BufferP buf, uint64_t nbits, uint64_t at_east_bits, BitsInitType init_to, float grow_at_least_factor, rom name, Caller caller);
#define buf_alloc_bits(vb, buf, more_bits, at_least_bits, init_to, grow_at_least_factor, name) \
    buf_alloc_bits_do ((VBlockP)(vb), (buf), (more_bits), (at_least_bits), (init_to), (grow_at_least_factor), (name), THIS_CODE_LINE)

extern BitsP buf_zfile_buf_to_bits (BufferP buf, uint64_t nbits);

//--------------------------
// thread synchronization
//--------------------------

extern void buf_init_lock (BufferP buf);

#define buf_lock_if(buf, cond) \
    BufferSpinlockP spinlock = (cond) ? (buf)->spinlock[0] : NULL; \
    ASSERT (!(cond) || spinlock, "spinlock not initialized for %s", buf_desc(buf).s); \
    if (spinlock) while (({ bool expected = (bool)false; !cas_strong_rel_acq (spinlock->lock, expected, (bool)true); })); /* spinlock */ 

#define buf_lock(buf) buf_lock_if ((buf), true)
#define buf_lock_(buf) rom func UNUSED = __FUNCTION__; buf_lock_if ((buf), true)

#define buf_unlock ({ if (spinlock) { __atomic_clear (&spinlock->lock, __ATOMIC_RELEASE); \
                                      spinlock = NULL; \
                                      /* printf ("unlocked %s\n", func); */}; })

extern BufferSpinlockP buf_lock_promiscuous (ConstBufferP buf, Caller caller);

static inline uint16_t buf_user_count (ConstBufferP buf) 
{
    return buf->memory ? BOLCOUNTER(buf) : 0;
}

extern const StrText1K buf_desc (ConstBufferP buf);
extern const StrText1K buf_desc_CTL (ConstBufferP buf);

extern rom buf_type_name (ConstBufferP buf);
