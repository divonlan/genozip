// ------------------------------------------------------------------
//   gencomp.h - "generated component"
//   Copyright (C) 2022-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"
#include "file.h"

// vb->gencomp is an array of these:
#define GENCOMP_LINE_LEN_BITS 30
#define GENCOMP_MAX_LINE_LEN ((1<<GENCOMP_LINE_LEN_BITS)-1)
typedef struct __attribute__ ((__packed__)) {
    uint32_t line_i;
    uint32_t line_index; // index into txt_data
    uint32_t line_len : GENCOMP_LINE_LEN_BITS;
    uint32_t comp_i   : 32-GENCOMP_LINE_LEN_BITS;
} GencompLineIEntry;

typedef enum { 
    GCT_NONE,
    // VBs of generated component are accumulated in memory, and compressed (out-of-band) when have enough data
    // in parallel with compression of MAIN data (DVCF PRIM and LUFT components, SAM PRIM component)
    GCT_OOB,  // out-of-band

    // VBs are accumulated in memory, and when the in-memory queue is full, they are offloaded to disk. They
    // are compressed only after all the main and out-of-band VBs are compressed. (SAM DEPN component)
    GCT_DEPN, // Dependent. Compressing requires some kind of ingestion of data of the previous VBs (eg SAM_COMP_DEPN)
    NUM_GC_TYPES
} GencompType;

extern void gencomp_seg_add_line (VBlockP vb, CompIType comp_i, STRp(line));
extern void gencomp_initialize (CompIType comp_i, GencompType gc_type);
extern void gencomp_destroy(void);
extern void gencomp_absorb_vb_gencomp_lines (VBlockP vb);
extern bool gencomp_get_txt_data (VBlockP vb);
extern void gencomp_a_main_vb_has_been_dispatched (void);
extern void gencomp_sam_prim_vb_has_been_ingested (VBlockP vb);

extern bool gencomp_comp_eligible_for_digest (VBlockP vb);
extern bool gencomp_am_i_expecting_more_txt_data (void);

extern int64_t gencomp_get_num_lines (CompIType comp_i);
