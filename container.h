// ------------------------------------------------------------------
//   container.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"
#include "sections.h"

#pragma pack(1)
#define CONTAINER_MAX_PREFIXES_LEN (32 * MAX_FIELDS)    // max len of just the names string, without the data eg "INFO1=INFO2=INFO3="
#define CON_PX_SEP              '\x4'  // starts the prefix string and terminates every prefix within it
#define CON_PX_SEP_             "\4"   // string version (careful not \x4 as it can combine with the next character to eg \x4F)
#define CON_PX_SEP_SHOW_REPEATS '\x5'  // an alternative terminator - outputs the number of repeats in LTEN32 after the prefix (used for BAM 'B' array count field)

#define CONTAINER_MAX_REPEATS 0xfffffe     // 3 byte unsigned int (16M) (minus 1 for easier detection of overflows)
#define CONTAINER_MAX_SELF_TRANS_CHANGE 50

#define CONTAINER_MAX_DICTS 2047           // 11 bits -matches CONTAINER_FIELDS.nitems_lo+nitems_hi

typedef struct ContainerItem {
    DictId dict_id;                        // note: the code counts on this field being first (assigning "item = { dict_id }")
    uint8_t did_i_small;                   // PIZ only: can store dids 0->254, 255 means did_i too large to store

    // special values of seperator[0]
    #define CI0_NONE        ((uint8_t)0x00) // no seperator 
    #define CI0_INVISIBLE   ((uint8_t)0x01) // this item does not appear in the original or reconstructed text. it should be consumed with reconstruct=false
    #define CI0_FIXED_0_PAD ((uint8_t)0x02) // fixed width, zero-left-padded, width in sep[2] (introduced v13)

    // separator[0] values with bit 7 set (0x80) are interpreted as flags rather than a separator, in 
    // which case separator[1] is a parameter of the flags
    #define CI0_ITEM_HAS_FLAG(item) ((uint8_t)(item)->separator[0] & 0x80)
    #define CI0_TRANS_NUL   ((uint8_t)0x81) // In translated mode: '\0' separator 
    #define CI0_TRANS_NOR   ((uint8_t)0x82) // In translated mode: reconstruct prefix, consume but don't reconstruct value. If item is a sub-container - NOT propagated to this sub-container.
    #define CI0_TRANS_MOVE  ((uint8_t)0x84) // (ORed) in addition: in translated: txt_data.len is moved separator[1] bytes (0-255), after all recontruction and/or translation
    #define CI0_NATIVE_NEXT ((uint8_t)0x88) // in non-translated mode: separator is in separator[1]

    // special values of separator[1]
    #define CI1_NONE        ((uint8_t)0x00) // no seperator 
    #define CI1_ITEM_CB     ((uint8_t)0x01) // item callback

    uint8_t separator[2];                   // 2 byte separator reconstructed after the item (or flags)
    
    TranslatorId translator;                // instructions how to translate this item, if this Container is reconstructed translating from one data type to another
} ContainerItem;

// container snip: it starts with SNIP_CONTAINER, following by a base64 of a big endian Container, with the number of
// items actually used and optoinally followed by a prefix array.
// The prefix array, if it exists, starts with CON_PX_SEP followed by the container-wide prefix, followed
// by a prefix for each item. Every prefix, if provided, is terminated by CON_PX_SEP.
// Only the container-wide prefix may alternatively be terminated by CON_PX_SEP_SHOW_REPEATS.

#define CONTAINER_FIELDS(nitems)        \
    uint32_t nitems_hi            : 3;  /* nitems_lo+hi=11 bits, matches CONTAINER_MAX_DICTS. MSB of num_items (until 9.0.22 it was MSB of repeats (after BGEN)), until 12.0.27 nitems_hi was 8 bits */ \
    /* container flags set during reconstruction */ \
    uint32_t unused               : 4;  /* can be used to enlarge nitems_hi or to add flags */ \
    uint32_t no_translation       : 1;  /* Cancel translation for this container and all of its items */\
    \
    uint32_t repeats              : 24; /* number of "repeats" (array elements) */ \
    uint8_t nitems_lo;                  /* LSB of num_items */  \
    /* container flags set during Seg */               \
    uint8_t drop_final_item_sep_of_final_repeat : 1; /* Deprecated - should not be used in new code. drop separator of final item of FINAL repeat */  \
    uint8_t drop_final_repeat_sep : 1;  \
    uint8_t filter_repeats        : 1; /* filter called before reconstruction of each repeat to determine if it should be reconstructed */ \
    uint8_t filter_items          : 1; /* filter called before reconstruction of each item to determine if it should be reconstructed */ \
    uint8_t is_toplevel           : 1;  \
    uint8_t keep_empty_item_sep   : 1; /* normally, we delete the separator preceding an empty item. this flag supprnor its repeat separator is reconstructed */ \
    uint8_t callback              : 1; /* callback called after reconstruction of each repeat (introduced 10.0.6) */ \
    uint8_t drop_final_item_sep   : 1; /* drop separator of final item of each repeat (introduced v12) */ \
    char repsep[2];                    /* repeat separator - two bytes that appear at the end of each repeat (ignored if 0) */ \
    ContainerItem items[nitems];

typedef struct Container       { CONTAINER_FIELDS(MAX_FIELDS)        } Container;
typedef struct MiniContainer   { CONTAINER_FIELDS(1)                 } MiniContainer;

#define SMALL_CON_NITEMS 16
typedef struct SmallContainer  { CONTAINER_FIELDS(SMALL_CON_NITEMS)  } SmallContainer;

#define MEDIUM_CON_NITEMS 80
typedef struct MediumContainer { CONTAINER_FIELDS(MEDIUM_CON_NITEMS) } MediumContainer;

#pragma pack()
#define con_nitems(con) ((con).nitems_hi * 256 + (con).nitems_lo)
#define con_set_nitems(con, n) { ASSERT (n <= CONTAINER_MAX_DICTS, "Container nitems=%u too large, max is %u", (unsigned)(n), CONTAINER_MAX_DICTS);\
                                 (con).nitems_hi = (n) >> 8; (con).nitems_lo = ((n) & 0xff); }
#define con_inc_nitems(con) con_set_nitems ((con), con_nitems (con) + 1)
#define con_dec_nitems(con) con_set_nitems ((con), con_nitems (con) - 1)
#define con_sizeof(con) (sizeof(con) - sizeof((con).items) + con_nitems (con) * sizeof((con).items[0]))

// prefixes may be NULL, or:
// [0] CON_PX_SEP - start
//     container-wide-prefix (may be empty)  + CON_PX_SEP
//     a prefix for each item (may be empty) + CON_PX_SEP
// empty prefixes of trailing items may be omitted
extern void container_prepare_snip (ConstContainerP con, STRp(prefixes), char *snip, unsigned *snip_len);
extern WordIndex container_seg_do (VBlockP vb, ContextP ctx, ConstContainerP con, STRp(prefixes), STRp(ren_prefixes), unsigned add_bytes, bool *is_new);
#define container_seg(vb, ctx, con, prefixes, prefixes_len, add_bytes) container_seg_do ((VBlockP)(vb), (ctx), (con), (prefixes), (prefixes_len), 0, 0, (add_bytes), NULL)
#define container_seg_with_rename(vb, ctx, con, prefixes, prefixes_len, ren_prefixes, ren_prefixes_len, add_bytes, is_new) container_seg_do ((VBlockP)(vb), (ctx), (con), (prefixes), (prefixes_len), (ren_prefixes), (ren_prefixes_len), (add_bytes), (is_new))
#define container_seg_by_dict_id(vb,dict_id,con,add_bytes) container_seg (vb, ctx_get_ctx (vb, dict_id), con, NULL, 0, add_bytes)

extern ValueType container_reconstruct (VBlockP vb, ContextP ctx, ConstContainerP con, STRp(prefixes));
extern ContainerP container_retrieve (VBlockP vb, ContextP ctx, WordIndex word_index, STRp(snip), pSTRp(out_prefixes));

extern void container_display (ConstContainerP con);

CONTAINER_FILTER_FUNC (container_no_filter);

// Translators reconstructing last_value as a little endian binary
TRANSLATOR_FUNC (container_translate_I8);   
TRANSLATOR_FUNC (container_translate_U8);   
TRANSLATOR_FUNC (container_translate_LTEN_I16);   
TRANSLATOR_FUNC (container_translate_LTEN_U16);   
TRANSLATOR_FUNC (container_translate_LTEN_I32);   
TRANSLATOR_FUNC (container_translate_LTEN_U32);   
