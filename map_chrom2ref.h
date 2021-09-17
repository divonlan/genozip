// ------------------------------------------------------------------
//   map_chrom2ref.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef MAP_CHROM2REF_INCLUDED
#define MAP_CHROM2REF_INCLUDED

#include "reference.h"

extern void map_chrom2ref_load (Reference ref);
extern void map_chrom2ref_compress (Reference ref);

// gets the index of the ref_contig index against which this chrom was compressed with --reference 
#define map_chrom2ref(chrom_index) \
    (buf_is_alloc (&z_file->chrom2ref_map) ? *ENT (WordIndex, z_file->chrom2ref_map, (chrom_index)) : (chrom_index))

#endif