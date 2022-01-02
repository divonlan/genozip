// ------------------------------------------------------------------
//   platform.h
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "../genozip.h"

#define ALPHABET_SIZE 256

extern void *(* bsc_malloc_do)(void *vb, size_t size, const char *func, uint32_t code_line);
#define bsc_malloc(vb, size) bsc_malloc_do (vb, size, __FUNCTION__, __LINE__)

extern void (* bsc_free_do)(void *vb, void *address, const char *func, uint32_t code_line);
#define bsc_free(vb, address) bsc_free_do ((vb), (address), __FUNCTION__, __LINE__)

extern void *bsc_zero_malloc (void *vb, size_t size);

