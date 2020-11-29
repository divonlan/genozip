// ------------------------------------------------------------------
//   container.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef CONTAINER_INCLUDED
#define CONTAINER_INCLUDED

#include "genozip.h"
#include "sections.h"

#pragma pack(1)
#define CONTAINER_MAX_PREFIXES_LEN 1000 // max len of just the names string, without the data eg "INFO1=INFO2=INFO3="
#define CON_PREFIX_SEP              '\x4'   // starts the prefix string and terminates every prefix within it
#define CON_PREFIX_SEP_SHOW_REPEATS '\x5'   // an alternative terminator - outputs the number of repeats in LTEN32 after the prefix (used for BAM 'B' array count field)

#define CONTAINER_MAX_REPEATS (UINT32_MAX-1) // one less than maxuint32 to make it easier to loop with con.repeats without overflow 
#define CONTAINER_MAX_SELF_TRANS_CHANGE 50

typedef struct ContainerItem {
    DictId dict_id;  
    DidIType did_i;                       // Used only in PIZ, must remain DID_I_NONE in ZIP

    // seperator[0] values with bit 7 set (0x80) are interpreted as flags rather than a seperator, in 
    // which case seperator[1] is a parameter of the flags
    #define CI_ITEM_HAS_FLAG(item) ((uint8_t)(item)->seperator[0] & 0x80)
    #define CI_TRANS_NUL   ((uint8_t)0x81) // In translated mode: '\0' seperator 
    #define CI_TRANS_NOR   ((uint8_t)0x82) // In translated mode: NO RECONSTRUCT: don't reconstruct number, just store it in last_value (not implemented for LT_SEQUENCE, LT_BITMAP, Containers, Sequences)
    #define CI_TRANS_MOVE  ((uint8_t)0x84) // (ORed) in addition: in translated: txt_data.len is moved seperator[1] bytes (0-255), after all recontruction and/or translation
    #define CI_NATIVE_NEXT ((uint8_t)0x88) // in non-translated mode: seperator is in seperator[1]
    uint8_t seperator[2];                 // 2 byte seperator reconstructed after the item (or flags)
    
    TranslatorId translator;              // instructions how to translate this item, if this Container is reconstructed translating from one data type to another
} ContainerItem;

// container snip: it starts with SNIP_CONTAINER, following by a base64 of a big endian Container, with the number of
// items actually used and optoinally followed by a prefix array.
// The prefix array, if it exists, starts with CON_PREFIX_SEP followed by the container-wide prefix, followed
// by a prefix for each item. Every prefix, if provided, is terminated by CON_PREFIX_SEP.
// Only the container-wide prefix may alternatively be terminated by CON_PREFIX_SEP_SHOW_REPEATS.

typedef struct Container {
    uint32_t repeats;     // number of "repeats" (array elements)
    uint8_t num_items;    // 1 to MAX_CONTAINER_ITEMS

    // container flags
    uint8_t drop_final_item_sep   : 1;
    uint8_t drop_final_repeat_sep : 1;
    uint8_t filter_repeats        : 1;
    uint8_t filter_items          : 1;
    uint8_t is_toplevel           : 1;
    
    char repsep[2];       // repeat seperator - two bytes that appear at the end of each repeat (ignored if 0)
    ContainerItem items[MAX_SUBFIELDS];
} Container;

// identical to Container, but with only one item
typedef struct MiniContainer {
    uint32_t repeats;     // number of "repeats" (array elements)
    uint8_t num_items;    // must be 1

    // container flags
    uint8_t drop_final_item_sep   : 1;
    uint8_t drop_final_repeat_sep : 1;
    uint8_t filter_repeats        : 1;
    uint8_t filter_items          : 1;
    uint8_t is_toplevel           : 1;

    char repsep[2];       // repeat seperator - two bytes that appear at the end of each repeat (ignored if 0)
    ContainerItem items[1];
} MiniContainer;

#pragma pack()
#define sizeof_container(con) (sizeof(con) - sizeof((con).items) + (con).num_items * sizeof((con).items[0]))

extern WordIndex container_seg_by_ctx (VBlockP vb, ContextP ctx, ContainerP con, const char *prefixes, unsigned prefixes_len, unsigned add_bytes);
#define container_seg_by_dict_id(vb,dict_id,con,add_bytes) container_seg_by_ctx ((VBlockP)vb, ctx_get_ctx (vb, dict_id), con, NULL, 0, add_bytes)

extern void container_reconstruct (VBlockP vb, ContextP ctx, WordIndex word_index, const char *snip, unsigned snip_len);

// Translators reconstructing last_value as a little endian binary
TRANSLATOR_FUNC (container_translate_I8);   
TRANSLATOR_FUNC (container_translate_U8);   
TRANSLATOR_FUNC (container_translate_LTEN_I16);   
TRANSLATOR_FUNC (container_translate_LTEN_U16);   
TRANSLATOR_FUNC (container_translate_LTEN_I32);   
TRANSLATOR_FUNC (container_translate_LTEN_U32);   

#endif
