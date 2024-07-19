// ------------------------------------------------------------------
//   endianness.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.
//   Copyright claimed on additions and modifications.
//
// This file is derived from https://gist.github.com/jtbr/7a43e6281e6cca353b33ee501421860c, see unmodified license statement below

/** 
* @file   endianness.h
* @brief  Convert Endianness of shorts, longs, long longs, regardless of architecture/OS
*
* Defines (without pulling in platform-specific network include headers):
* bswap16, bswap32, bswap64, ntoh16, hton16, ntoh32 hton32, ntoh64, hton64
*
* Should support linux / macos / solaris / windows.
* Supports GCC (on any platform, including embedded), MSVC2015, and clang,
* and should support intel, solaris, and ibm compilers as well.
*
* Released under MIT license
*/

#pragma once

#include <stdlib.h>
#include <stdint.h>

/* Detect platform endianness at compile time */

// If boost were available on all platforms, could use this instead to detect endianness
// #include <boost/predef/endian.h>

// When available, these headers can improve platform endianness detection
#ifdef __has_include // C++17, supported as extension to C++11 in clang, GCC 5+, vs2015
#  if __has_include(<endian.h>)
#    include <endian.h> // gnu libc normally provides, linux
#  elif __has_include(<machine/endian.h>)
#    include <machine/endian.h> //open bsd, macos
#  elif __has_include(<sys/param.h>)
#    include <sys/param.h> // mingw, some bsd (not open/macos)
#  elif __has_include(<sys/isadefs.h>)
#    include <sys/isadefs.h> // solaris
#  endif
#endif

#if !defined(__LITTLE_ENDIAN__) && !defined(__BIG_ENDIAN__)
#  if (defined(__BYTE_ORDER__)  && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__) || \
     (defined(__BYTE_ORDER) && __BYTE_ORDER == __BIG_ENDIAN) || \
	 (defined(_BYTE_ORDER) && _BYTE_ORDER == _BIG_ENDIAN) || \
	 (defined(BYTE_ORDER) && BYTE_ORDER == BIG_ENDIAN) || \
     (defined(__sun) && defined(__SVR4) && defined(_BIG_ENDIAN)) || \
     defined(__ARMEB__) || defined(__THUMBEB__) || defined(__AARCH64EB__) || \
     defined(_MIBSEB) || defined(__MIBSEB) || defined(__MIBSEB__) || \
     defined(_M_PPC)
#    define __BIG_ENDIAN__
#  elif (defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__) || /* gcc */\
     (defined(__BYTE_ORDER) && __BYTE_ORDER == __LITTLE_ENDIAN) /* linux header */ || \
	 (defined(_BYTE_ORDER) && _BYTE_ORDER == _LITTLE_ENDIAN) || \
	 (defined(BYTE_ORDER) && BYTE_ORDER == LITTLE_ENDIAN) /* mingw header */ ||  \
     (defined(__sun) && defined(__SVR4) && defined(_LITTLE_ENDIAN)) || /* solaris */ \
     defined(__ARMEL__) || defined(__THUMBEL__) || defined(__AARCH64EL__) || \
     defined(_MIPSEL) || defined(__MIPSEL) || defined(__MIPSEL__) || \
     defined(_M_IX86) || defined(_M_X64) || defined(_M_IA64) || /* msvc for intel processors */ \
     defined(_M_ARM) /* msvc code on arm executes in little endian mode */
#    define __LITTLE_ENDIAN__
#  endif
#endif

#if defined(bswap16) || defined(bswap32) || defined(bswap64)
#  error "unexpected define!" // freebsd may define these; probably just need to undefine them
#endif

/* Define byte-swap functions, using fast processor-native built-ins where possible */
#if defined(_MSC_VER) // needs to be first because msvc doesn't short-circuit after failing defined(__has_builtin)
#  define bswap16(x)     _byteswap_ushort((x))
#  define bswap32(x)     _byteswap_ulong((x))
#  define bswap64(x)     _byteswap_uint64((x))
#elif (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 8)
#  define bswap16(x)     __builtin_bswap16((x))
#  define bswap32(x)     __builtin_bswap32((x))
#  define bswap64(x)     __builtin_bswap64((x))
#elif defined(__has_builtin) 
#  if __has_builtin(__builtin_bswap64)  /* for clang; gcc 5 fails on this and && shortcircuit fails; must be after GCC check */
#     define bswap16(x)     __builtin_bswap16((x))
#     define bswap32(x)     __builtin_bswap32((x))
#     define bswap64(x)     __builtin_bswap64((x))
#  endif
#endif

#define CONST_BSWAP64(x) \
               ((( x & 0xff00000000000000ull ) >> 56 ) | \
                (( x & 0x00ff000000000000ull ) >> 40 ) | \
                (( x & 0x0000ff0000000000ull ) >> 24 ) | \
                (( x & 0x000000ff00000000ull ) >> 8  ) | \
                (( x & 0x00000000ff000000ull ) << 8  ) | \
                (( x & 0x0000000000ff0000ull ) << 24 ) | \
                (( x & 0x000000000000ff00ull ) << 40 ) | \
                (( x & 0x00000000000000ffull ) << 56 ))

#ifndef bswap64 
    /* even in this case, compilers often optimize by using native instructions */
    static inline uint16_t bswap16(uint16_t x) {
		return ((( x  >> 8 ) & 0xffu ) | (( x  & 0xffu ) << 8 ));
	}
    static inline uint32_t bswap32(uint32_t x) {
        return ((( x & 0xff000000u ) >> 24 ) |
                (( x & 0x00ff0000u ) >> 8  ) |
                (( x & 0x0000ff00u ) << 8  ) |
                (( x & 0x000000ffu ) << 24 ));
    }
    static inline uint64_t bswap64(uint64_t x) {
        return CONST_BSWAP64(x);
    }
#endif

#ifdef __LITTLE_ENDIAN__
#define BGEN16(x) bswap16(x)
#define BGEN24(x) (bswap32(x) >> 8)
#define BGEN32(x) bswap32(x)
#define BGEN64(x) bswap64(x)
#define CONST_BGEN64(x) CONST_BSWAP64(x) // when applied to a constant, it is computed in compilation time
#define LTEN16(x) (x)
#define LTEN32(x) (x)
#define LTEN64(x) (x)
#define CONST_LTEN64(x) (x)
#elif defined __BIG_ENDIAN__
#define LTEN16(x) bswap16(x)
#define LTEN32(x) bswap32(x)
#define LTEN64(x) bswap64(x)
#define CONST_LTEN64(x) CONST_BSWAP64(x) // when applied to a constant, it is computed in compilation time
#define BGEN16(x) (x)
#define BGEN24(x) (x)
#define BGEN32(x) (x)
#define BGEN64(x) (x)
#define CONST_BGEN64(x) (x)
#else
#error "Neither __LITTLE_ENDIAN__ nor __BIG_ENDIAN__ is defined"
#endif

#define LTEN32F(x) ({ union { float f;  uint32_t i; } u = {.f = (x)}; u.i = LTEN32(u.i); u.f; })
#define BGEN32F(x) ({ union { float f;  uint32_t i; } u = {.f = (x)}; u.i = BGEN32(u.i); u.f; })
#define LTEN64F(x) ({ union { double f; uint64_t i; } u = {.f = (x)}; u.i = LTEN64(u.i); u.f; })
#define BGEN64F(x) ({ union { double f; uint64_t i; } u = {.f = (x)}; u.i = BGEN64(u.i); u.f; })

// getting and putting unaligned Little Endian words
#ifdef GENOZIP_ALLOW_UNALIGNED_ACCESS
    #define GET_UINT16(p)  LTEN16 (*((uint16_t *)(p)))
    #define GET_UINT32(p)  LTEN32 (*((uint32_t *)(p)))
    #define GET_UINT64(p)  LTEN64 (*((uint64_t *)(p)))
    #define GET_FLOAT32(p) LTEN32F(*((float    *)(p)))

    #define GET_UINT32_(st_p, member) ((st_p)->member)

    #define PUT_UINT16(p,n) *((uint16_t *)(p)) = LTEN16(n)
    #define PUT_UINT32(p,n) *((uint32_t *)(p)) = LTEN32(n)

    #define PUT_UINT16_(st_p, member, n) (st_p)->member = LTEN16(n)
    #define PUT_UINT32_(st_p, member, n) (st_p)->member = LTEN32(n)
#else
    // loading a Little Endian uint32_t from an unaligned memory location
    #define GET_UINT16(p)  ((uint16_t)((uint8_t*)(p))[0] | ((uint16_t)((uint8_t*)(p))[1] << 8))
    #define GET_UINT32(p)  ((uint32_t)((uint8_t*)(p))[0] | ((uint32_t)((uint8_t*)(p))[1] << 8) | ((uint32_t)((uint8_t*)(p))[2] << 16) | ((uint32_t)((uint8_t*)(p))[3] << 24))
    #define GET_UINT64(p)  ((uint64_t)((uint8_t*)(p))[0] | ((uint64_t)((uint8_t*)(p))[1] << 8) | ((uint64_t)((uint8_t*)(p))[2] << 16) | ((uint64_t)((uint8_t*)(p))[3] << 24) | ((uint64_t)((uint8_t*)(p))[4] << 32) | ((uint64_t)((uint8_t*)(p))[5] << 40) | ((uint64_t)((uint8_t*)(p))[6] << 48) | ((uint64_t)((uint8_t*)(p))[7] << 56)))
    #define GET_FLOAT32(p) ({ union { uint32_t i; float f; } n= {.i = GET_UINT32(p)}; n.f; })

    #define GET_UINT32_(st_p, member) ({ typeof(*st_p) dummy; bytes _P=(bytes)(st_p) + ((bytes)&dummy.member - (bytes)&dummy); ((uint32_t)_P[0] | ((uint32_t)_P[1] << 8) | ((uint32_t)_P[2] << 16) | ((uint32_t)_P[3] << 24)); })

    // storing a Little Endian integer in an unaligned memory location
    #define PUT_UINT16(p,n) ({ uint16_t _N=(n); uint8_t *_P=(uint8_t *)(p); _P[0]=_N; _P[1]=_N>>8; }) 
    #define PUT_UINT32(p,n) ({ uint32_t _N=(n); uint8_t *_P=(uint8_t *)(p); _P[0]=_N; _P[1]=_N>>8; _P[2]=_N>>16; _P[3]=_N>>24; })

    // storing as a struct member
    #define PUT_UINT16_(st_p, member, n) ({ typeof(*st_p) dummy; PUT_UINT16 ((rom)(st_p) + ((rom)&dummy.member - (rom)&dummy), LTEN16(n)); })
    #define PUT_UINT32_(st_p, member, n) ({ typeof(*st_p) dummy; PUT_UINT32 ((rom)(st_p) + ((rom)&dummy.member - (rom)&dummy), LTEN32(n)); })
#endif

#define GET_UINT8(p)   ((uint8_t)(((uint8_t*)(p))[0]))
#define GET_UINT24(p)  ((uint32_t)(((uint8_t*)(p))[0] | (((uint8_t*)(p))[1] << 8))| (((uint8_t*)(p))[2] << 16))

#define PUT_UINT8(p,n)  ({ ((uint8_t*)(p))[0] = (n); })
#define PUT_UINT24(p,n) ({ uint32_t _N=(n); uint8_t *_P=(uint8_t *)(p); _P[0]=_N; _P[1]=_N>>8; _P[2]=_N>>16; })
