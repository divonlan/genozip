// ------------------------------------------------------------------
//   platform.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "../genozip.h"

#define ALPHABET_SIZE 256

extern void *(* bsc_malloc)(void *vb, size_t size);
extern void  (* bsc_free)(void *vb, void *address);
extern void *bsc_zero_malloc (void *vb, size_t size);

