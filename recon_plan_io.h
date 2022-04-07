// ------------------------------------------------------------------
//   recon_plan_io.h
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

// ZIP
extern void recon_plan_compress (uint32_t my_conc_writing_vbs, bool my_is_luft);

// PIZ
extern void recon_plan_uncompress (Section sec, uint32_t *out_conc_writing_vbs, uint32_t *out_vblock_mb);

// Misc
extern void recon_plan_show (FileP file, bool is_luft, uint32_t conc_writing_vbs, uint32_t vblock_mb);

extern rom recon_plan_flavors[8];
