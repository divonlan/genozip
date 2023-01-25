// ------------------------------------------------------------------
//   deep.h
//   Copyright (C) 2023-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "sam_private.h"

// Stuff that happens during SAM seg
extern void deep_sam_seg_initialize (VBlockSAMP vb);
extern void deep_sam_set_QNAME_hash (VBlockSAMP vb, ZipDataLineSAM *dl, uint32_t qname_hash);
extern void deep_sam_set_SEQ_hash (VBlockSAMP vb,ZipDataLineSAM *dl, STRp(textual_seq));
extern void deep_sam_set_QUAL_hash (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qual));

#define IS_DEEP_ALN (flag.deep && !segconf.running && FLAG_IS_PRIMARY(dl->FLAG))

