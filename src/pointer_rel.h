// ------------------------------------------------------------------
//   pointeר.h
//   Copyright (C) 2026-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include <stdint.h>
// don't #include ANY genozip files here!

// separate header file it can be included in the 3rd party library, where inclusion
// of the whole genozip.h would cause conflicts.

// Relative strings - marked with ר (Hebrew R): static strings (e.g. in .rodata)

typedef int32_t Pointeר;          // a pointer into the static area (instructions, static data) extressed as delta vs string_anchor - 32b instead of 64b
extern const char *string_anchor; // base address, relative to which other strings are expressed

// ר is only safe with static strings (string literals or pointers to static memory), bc Pointeר is only int32_t. 
// for string literals, it will evaluate to a constant by the compiler+linker
#define ר(static_string) ((Pointeר)((uintptr_t)(void *)static_string - (uintptr_t)(void *)string_anchor)) // static strings only!

#define unר(strר) ((rom)((int64_t)(strר) + (int64_t)string_anchor))

typedef struct Caller {
    Pointeר funcר; // function name is always a static += 1MB from anchor with .rodata (Linux) / .rdata (Windows) / __TEXT,__cstring (Mac)
    uint16_t code_line;
} Caller;

#define CALLERff(caller) unר(caller.funcר), caller.code_line /* for printf arguments */
#define CALLERf CALLERff(caller)
#define THIS_CODE_LINE ((Caller){ .funcר = ר(__FUNCTION__), .code_line = __LINE__})
