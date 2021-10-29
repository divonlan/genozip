// ------------------------------------------------------------------
//   coords.h
//   Copyright (C) 2021-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#define NUM_COORDS 4
typedef enum __attribute__ ((__packed__)) { DC_NONE, DC_PRIMARY, DC_LUFT, DC_BOTH } Coords; // 2 bits, part of the file format, goes into FlagsTxtHeader, FlagsVbHeader

#define OTHER_COORDS(c) ((c)==DC_PRIMARY ? DC_LUFT : DC_PRIMARY)
#define SEL(prim,luft) ((vb->line_coords == DC_PRIMARY) ? (prim) : (luft))

extern const char *coords_name (Coords coord);
