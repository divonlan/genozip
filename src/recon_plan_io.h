// ------------------------------------------------------------------
//   recon_plan_io.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "genozip.h"

// ZIP
extern void recon_plan_compress (uint32_t my_conc_writing_vbs);

// PIZ
extern void recon_plan_uncompress (Section sec, uint32_t *out_conc_writing_vbs, uint32_t *out_vblock_mb);

// Misc
extern void recon_plan_show (FileP file, uint32_t conc_writing_vbs, uint32_t vblock_mb);

extern rom recon_plan_flavors[8];
