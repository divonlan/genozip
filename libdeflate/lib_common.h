/*
 * lib_common.h - internal header included by all library code
 */

#ifndef LIB_LIB_COMMON_H
#define LIB_LIB_COMMON_H

#ifdef LIBDEFLATE_H
#  error "lib_common.h must always be included before libdeflate.h"
   /* because BUILDING_LIBDEFLATE must be set first */
#endif

#define BUILDING_LIBDEFLATE

#include "common_defs.h"

/*
 * Prefix with "_libdeflate_" all global symbols which are not part of the API
 * and don't already have a "libdeflate" prefix.  This avoids exposing overly
 * generic names when libdeflate is built as a static library.
 *
 * Note that the chosen prefix is not really important and can be changed
 * without breaking library users.  It was just chosen so that the resulting
 * symbol names are unlikely to conflict with those from any other software.
 * Also note that this fixup has no useful effect when libdeflate is built as a
 * shared library, since these symbols are not exported.
 */
#define SYM_FIXUP(sym)			_libdeflate_##sym
#define deflate_get_compression_level	SYM_FIXUP(deflate_get_compression_level)
#define _cpu_features			SYM_FIXUP(_cpu_features)
#define setup_cpu_features		SYM_FIXUP(setup_cpu_features)

void *libdeflate_malloc_do(size_t size, void *opaque, const char *func, uint32_t code_line);
#define libdeflate_malloc(size, opaque) libdeflate_malloc_do ((size), (opaque), __FUNCTION__, __LINE__)

void libdeflate_free_do(void *ptr, void *opaque, const char *func, uint32_t code_line);
#define libdeflate_free(ptr, opaque) libdeflate_free_do ((ptr), (opaque), __FUNCTION__, __LINE__)

void *libdeflate_aligned_malloc(size_t alignment, size_t size, void *opaque);
void libdeflate_aligned_free(void *ptr, void *opaque);

#ifdef FREESTANDING
/*
 * With -ffreestanding, <string.h> may be missing, and we must provide
 * implementations of memset(), memcpy(), memmove(), and memcmp().
 * See https://gcc.gnu.org/onlinedocs/gcc/Standards.html
 *
 * Also, -ffreestanding disables interpreting calls to these functions as
 * built-ins.  E.g., calling memcpy(&v, p, WORDBYTES) will make a function call,
 * not be optimized to a single load instruction.  For performance reasons we
 * don't want that.  So, declare these functions as macros that expand to the
 * corresponding built-ins.  This approach is recommended in the gcc man page.
 * We still need the actual function definitions in case gcc calls them.
 */
void *memset(void *s, int c, size_t n);
#define memset(s, c, n)		__builtin_memset((s), (c), (n))

void *memcpy(void *dest, const void *src, size_t n);
#define memcpy(dest, src, n)	__builtin_memcpy((dest), (src), (n))

void *memmove(void *dest, const void *src, size_t n);
#define memmove(dest, src, n)	__builtin_memmove((dest), (src), (n))

int memcmp(const void *s1, const void *s2, size_t n);
#define memcmp(s1, s2, n)	__builtin_memcmp((s1), (s2), (n))
#else
#include <string.h>
#endif

#endif /* LIB_LIB_COMMON_H */
