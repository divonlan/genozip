/*
 * libdeflate.h - public header for libdeflate
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// #define LIBDEFLATE_VERSION_MAJOR	1
// #define LIBDEFLATE_VERSION_MINOR	7
// #define LIBDEFLATE_VERSION_STRING	"1.7"

#include <stddef.h>
#include <stdint.h>

/*
 * On Windows, if you want to link to the DLL version of libdeflate, then
 * #define LIBDEFLATE_DLL.  Note that the calling convention is "stdcall".
 */
#ifdef LIBDEFLATE_DLL
#  ifdef BUILDING_LIBDEFLATE
#    define LIBDEFLATEEXPORT	LIBEXPORT
#  elif defined(_WIN32) || defined(__CYGWIN__)
#    define LIBDEFLATEEXPORT	__declspec(dllimport)
#  endif
#endif
#ifndef LIBDEFLATEEXPORT
#  define LIBDEFLATEEXPORT
#endif

#if defined(_WIN32) && !defined(_WIN64)
#  define LIBDEFLATEAPI_ABI	__stdcall
#else
#  define LIBDEFLATEAPI_ABI
#endif

#if defined(BUILDING_LIBDEFLATE) && defined(__GNUC__) && \
	defined(_WIN32) && !defined(_WIN64)
    /*
     * On 32-bit Windows, gcc assumes 16-byte stack alignment but MSVC only 4.
     * Realign the stack when entering libdeflate to avoid crashing in SSE/AVX
     * code when called from an MSVC-compiled application.
     */
#  define LIBDEFLATEAPI_STACKALIGN	__attribute__((force_align_arg_pointer))
#else
#  define LIBDEFLATEAPI_STACKALIGN
#endif

#define LIBDEFLATEAPI_1_7	LIBDEFLATEAPI_ABI LIBDEFLATEAPI_STACKALIGN

/* ========================================================================== */
/*                             Compression                                    */
/* ========================================================================== */

struct libdeflate_compressor;

/*
 * libdeflate_alloc_compressor() allocates a new compressor that supports
 * DEFLATE, zlib, and gzip compression.  'compression_level' is the compression
 * level on a zlib-like scale but with a higher maximum value (1 = fastest, 6 =
 * medium/default, 9 = slow, 12 = slowest).  Level 0 is also supported and means
 * "no compression", specifically "create a valid stream, but only emit
 * uncompressed blocks" (this will expand the data slightly).
 *
 * The return value is a pointer to the new compressor, or NULL if out of memory
 * or if the compression level is invalid (i.e. outside the range [0, 12]).
 *
 * Note: for compression, the sliding window size is defined at compilation time
 * to 32768, the largest size permissible in the DEFLATE format.  It cannot be
 * changed at runtime.
 *
 * A single compressor is not safe to use by multiple threads concurrently.
 * However, different threads may use different compressors concurrently.
 */
LIBDEFLATEEXPORT struct libdeflate_compressor * LIBDEFLATEAPI_1_7
libdeflate_alloc_compressor_1_7(int compression_level, void *opaque);

/*
 * libdeflate_deflate_compress() performs raw DEFLATE compression on a buffer of
 * data.  The function attempts to compress 'in_nbytes' bytes of data located at
 * 'in' and write the results to 'out', which has space for 'out_nbytes_avail'
 * bytes.  The return value is the compressed size in bytes, or 0 if the data
 * could not be compressed to 'out_nbytes_avail' bytes or fewer.
 */
LIBDEFLATEEXPORT size_t LIBDEFLATEAPI_1_7
libdeflate_deflate_compress_1_7(struct libdeflate_compressor *compressor,
			    const void *in, size_t in_nbytes,
			    void *out, size_t out_nbytes_avail);

/*
 * libdeflate_deflate_compress_bound() returns a worst-case upper bound on the
 * number of bytes of compressed data that may be produced by compressing any
 * buffer of length less than or equal to 'in_nbytes' using
 * libdeflate_deflate_compress() with the specified compressor.  Mathematically,
 * this bound will necessarily be a number greater than or equal to 'in_nbytes'.
 * It may be an overestimate of the true upper bound.  The return value is
 * guaranteed to be the same for all invocations with the same compressor and
 * same 'in_nbytes'.
 *
 * As a special case, 'compressor' may be NULL.  This causes the bound to be
 * taken across *any* libdeflate_compressor that could ever be allocated with
 * this build of the library, with any options.
 *
 * Note that this function is not necessary in many applications.  With
 * block-based compression, it is usually preferable to separately store the
 * uncompressed size of each block and to store any blocks that did not compress
 * to less than their original size uncompressed.  In that scenario, there is no
 * need to know the worst-case compressed size, since the maximum number of
 * bytes of compressed data that may be used would always be one less than the
 * input length.  You can just pass a buffer of that size to
 * libdeflate_deflate_compress() and store the data uncompressed if
 * libdeflate_deflate_compress() returns 0, indicating that the compressed data
 * did not fit into the provided output buffer.
 */
LIBDEFLATEEXPORT size_t LIBDEFLATEAPI_1_7
libdeflate_deflate_compress_bound_1_7(struct libdeflate_compressor *compressor,
				  size_t in_nbytes);

/*
 * Like libdeflate_deflate_compress(), but stores the data in the zlib wrapper
 * format.
 */
LIBDEFLATEEXPORT size_t LIBDEFLATEAPI_1_7
libdeflate_zlib_compress_1_7(struct libdeflate_compressor *compressor,
			 const void *in, size_t in_nbytes,
			 void *out, size_t out_nbytes_avail);

/*
 * Like libdeflate_deflate_compress_bound(), but assumes the data will be
 * compressed with libdeflate_zlib_compress() rather than with
 * libdeflate_deflate_compress().
 */
LIBDEFLATEEXPORT size_t LIBDEFLATEAPI_1_7
libdeflate_zlib_compress_bound_1_7(struct libdeflate_compressor *compressor,
			       size_t in_nbytes);

/*
 * Like libdeflate_deflate_compress(), but stores the data in the gzip wrapper
 * format.
 */
LIBDEFLATEEXPORT size_t LIBDEFLATEAPI_1_7
libdeflate_gzip_compress_1_7(struct libdeflate_compressor *compressor,
			 const void *in, size_t in_nbytes,
			 void *out, size_t out_nbytes_avail);

/*
 * Like libdeflate_deflate_compress_bound(), but assumes the data will be
 * compressed with libdeflate_gzip_compress() rather than with
 * libdeflate_deflate_compress().
 */
LIBDEFLATEEXPORT size_t LIBDEFLATEAPI_1_7
libdeflate_gzip_compress_bound_1_7(struct libdeflate_compressor *compressor,
			       size_t in_nbytes);

/*
 * libdeflate_free_compressor() frees a compressor that was allocated with
 * libdeflate_alloc_compressor().  If a NULL pointer is passed in, no action is
 * taken.
 */
LIBDEFLATEEXPORT void LIBDEFLATEAPI_1_7
libdeflate_free_compressor_1_7(struct libdeflate_compressor *compressor);

/* ========================================================================== */
/*                           Custom memory allocator                          */
/* ========================================================================== */

/*
 * Install a custom memory allocator which libdeflate will use for all memory
 * allocations.  'malloc_func' is a function that must behave like malloc(), and
 * 'free_func' is a function that must behave like free().
 *
 * There must not be any libdeflate_compressor or libdeflate_decompressor
 * structures in existence when calling this function.
 */
LIBDEFLATEEXPORT void LIBDEFLATEAPI_1_7
libdeflate_set_memory_allocator_1_7 (void *(*malloc_func)(void *, unsigned, unsigned, const char*, uint32_t),
				void (*free_func)(void *, void *, const char*, uint32_t));

#ifdef __cplusplus
}
#endif

