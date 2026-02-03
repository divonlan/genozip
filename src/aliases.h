// ------------------------------------------------------------------
//   aliases.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "sections.h"

typedef struct __attribute__ ((packed)) { 
    AliasType alias_type; // added v15 - up to v14 DictIdAlias contained only alias, dst
    DictId alias, dst; 
} DictIdAlias;

extern void aliases_compress (void);
extern ConstBufferP aliases_get (void);
extern void show_aliases (void);
