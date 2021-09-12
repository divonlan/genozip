// ------------------------------------------------------------------
//   hash.c
//   Copyright (C) 2021-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef SEGCONF_INCLUDED
#define SEGCONF_INCLUDED

#include "genozip.h"

typedef struct {
    bool running;               // currently in segconf_calculate()

    // Seg parameters - general
    uint64_t vb_size;
    
    // Seg parameters - FASTA
    bool fasta_has_contigs;     // true if the sequences in this FASTA represent contigs (as opposed to reads) - in which case we have a FASTA_CONTIG dictionary and RANDOM_ACCESS
} SegConf;

extern SegConf segconf;

extern void segconf_calculate (void);

#endif
