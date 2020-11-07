// ------------------------------------------------------------------
//   container.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef CONTAINER_INCLUDED
#define CONTAINER_INCLUDED

#include "genozip.h"

// container snip: it starts with SNIP_CONTAINER, following by a base64 of a big endian Container
#pragma pack(1)

#define CONTAINER_DROP_FINAL_ITEM_SEP   0x01
#define CONTAINER_DROP_FINAL_REPEAT_SEP 0x02
#define CONTAINER_FILTER_REPEATS        0x04
#define CONTAINER_FILTER_ITEMS          0x08
#define CONTAINER_TOPLEVEL              0x10

#define CONTAINER_MAX_REPEATS 4294967294UL // one less than maxuint32 to make it easier to loop with con.repeats without overflow 
#define CONTAINER_MAX_PREFIXES_LEN 1000 // max len of just the names string, without the data eg "INFO1=INFO2=INFO3="

typedef struct ContainerItem {
    DictId dict_id;  
    DidIType did_i;                       // Used only in PIZ, must remain DID_I_NONE in ZIP

    // seperator[0] values with bit 7 set (0x80) are interpreted as flags rather than a seperator, in 
    // which case seperator[1] is a parameter of the flags
    #define CI_TRANS_NUL   ((uint8_t)0x81) // In translated mode: '\0' seperator 
    #define CI_TRANS_NOR   ((uint8_t)0x82) // In translated mode: NO RECONSTRUCT: don't reconstruct number, just store it in last_value (not implemented for LT_SEQUENCE, LT_BITMAP, Containers, Sequences)
    #define CI_TRANS_MOVE  ((uint8_t)0x84) // (ORed) in addition: in translated: txt_data.len is moved seperator[1] bytes (0-255), after all recontruction and/or translation
    #define CI_NATIVE_NEXT ((uint8_t)0x88) // in non-translated mode: seperator is in seperator[1]

    uint8_t seperator[2];                 // 2 byte seperator reconstructed after the item (or flags)
    TranslatorId translator;              // instructions how to translate this item, if this Container is reconstructed translating from one data type to another
} ContainerItem;

typedef struct Container {
    uint32_t repeats;     // number of "repeats" (array elements)
    uint8_t num_items;    // 1 to MAX_CONTAINER_ITEMS
    uint8_t flags;
    char repsep[2];       // repeat seperator - two bytes that appear at the end of each repeat (ignored if 0)
    ContainerItem items[MAX_SUBFIELDS];
} Container;

// identical to Container, but with only one item
typedef struct MiniContainer {
    uint32_t repeats;     // number of "repeats" (array elements)
    uint8_t num_items;    // must be 1
    uint8_t flags;
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
