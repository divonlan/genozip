// ------------------------------------------------------------------
//   coords.c
//   Copyright (C) 2021-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "coords.h"

const char *coords_name (Coords coord)
{
    static const char *coords_names[4] = { "NONE", "PRIM", "LUFT", "BOTH" };
    
    return (coord < 0 || coord >= NUM_COORDS) ? "(invalid coord)" : coords_names[coord];
}
