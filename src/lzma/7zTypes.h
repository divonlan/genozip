/* 7zTypes.h -- Basic types
2018-08-04 : Igor Pavlov : Public domain */

#pragma once

#include "stddef.h" // ptrdiff_t
#include "../genozip.h"

#define SZ_OK 0

#define SZ_ERROR_DATA 1
#define SZ_ERROR_MEM 2
#define SZ_ERROR_CRC 3
#define SZ_ERROR_UNSUPPORTED 4
#define SZ_ERROR_PARAM 5
#define SZ_ERROR_INPUT_EOF 6
#define SZ_ERROR_OUTPUT_EOF 7
#define SZ_ERROR_READ 8
#define SZ_ERROR_WRITE 9
#define SZ_ERROR_PROGRESS 10
#define SZ_ERROR_FAIL 11
#define SZ_ERROR_THREAD 12

#define SZ_ERROR_ARCHIVE 16
#define SZ_ERROR_NO_ARCHIVE 17

typedef int SRes;


#ifdef _WIN32

/* typedef DWORD WRes; */
typedef unsigned WRes;
#define MY_SRes_HRESULT_FROM_WRes(x) HRESULT_FROM_WIN32(x)

#else

typedef int WRes;
#define MY__FACILITY_WIN32 7
#define MY__FACILITY__WRes MY__FACILITY_WIN32
#define MY_SRes_HRESULT_FROM_WRes(x) ((HRESULT)(x) <= 0 ? ((HRESULT)(x)) : ((HRESULT) (((x) & 0x0000FFFF) | (MY__FACILITY__WRes << 16) | 0x80000000)))

#endif


#ifndef RINOK
#define RINOK(x) { int __result__ = (x); if (__result__ != 0) return __result__; }
#endif

#ifndef BYTE_DEFINED // also defined (identically) in zlib/zconf.h --divon
#define BYTE_DEFINED
typedef unsigned char Byte;
#endif
typedef short Int16;
typedef unsigned short UInt16;

typedef int32_t Int32;
typedef uint32_t UInt32;

typedef int64_t Int64;
typedef uint64_t UInt64;

#if defined(_MSC_VER) || defined(__BORLANDC__)
#define UINT64_CONST(n) n
#else
#define UINT64_CONST(n) n ## ULL
#endif

typedef uint64_t SizeT; // Divon change - always 64bit

typedef int BoolInt;

#define True 1
#define False 0

#ifdef _WIN32
#define MY_STD_CALL __stdcall
#else
#define MY_STD_CALL
#endif

#ifdef _MSC_VER

#if _MSC_VER >= 1300
#define MY_NO_INLINE __declspec(noinline)
#else
#define MY_NO_INLINE
#endif

#define MY_FORCE_INLINE __forceinline

#define MY_CDECL __cdecl
#define MY_FAST_CALL __fastcall

#else

#define MY_NO_INLINE
#define MY_FORCE_INLINE
#define MY_CDECL
#define MY_FAST_CALL

/* inline keyword : for C++ / C99 */

/* GCC, clang: */
/*
#if defined (__GNUC__) && (__GNUC__ >= 4)
#define MY_FORCE_INLINE __attribute__((always_inline))
#define MY_NO_INLINE __attribute__((noinline))
#endif
*/

#endif


/* The following interfaces use first parameter as pointer to structure */

typedef struct IByteIn IByteIn;
struct IByteIn
{
  Byte (*Read)(const IByteIn *p); /* reads one byte, returns 0 in case of EOF or error */
};
#define IByteIn_Read(p) (p)->Read(p)


typedef struct IByteOut IByteOut;
struct IByteOut
{
  void (*Write)(const IByteOut *p, Byte b);
};
#define IByteOut_Write(p, b) (p)->Write(p, b)


typedef struct ISeqInStream ISeqInStream;
struct ISeqInStream
{
  SRes (*Read)(const ISeqInStream *p, void *buf, size_t *size);
    /* if (input(*size) != 0 && output(*size) == 0) means end_of_stream.
       (output(*size) < input(*size)) is allowed */

  // all the remaining fields (except for Read) were added by Divon
  void *vb;        
  void *ctx;            
  unsigned line_i;        // current partially-consumed line, or next line if no partially consumed line
  unsigned avail_in;      // total bytes remaining
  
  // we refer to two line segments (eg SEQ and E2 or QUAL and U2)
  char *next_in_1;        // if there is a line partially-consumed, this is the next byte
  char *next_in_2;        // if there is a line partially-consumed, and next_in==NULL, this is the next byte
  unsigned avail_in_1;    // if there is a line partially-consumed (0 if not)
  unsigned avail_in_2;    // if there is a line partially-consumed (0 if not)
  void (*callback)();     // callback to fetch more data
};
#define ISeqInStream_Read(p, buf, size) (p)->Read(p, buf, size)

typedef struct ISeqOutStream ISeqOutStream;
struct ISeqOutStream
{
  size_t (*Write)(const ISeqOutStream *p, const void *buf, size_t size);
    /* Returns: result - the number of actually written bytes.
       (result < size) means error */
  char *next_out;                     // added by Divon: next compressed data goes here
  unsigned avail_out;                 // added by Divon: remaining space
};

#define ISeqOutStream_Write(p, buf, size) (p)->Write(p, buf, size)

#ifndef MY_offsetof
  #ifdef offsetof
    #define MY_offsetof(type, m) offsetof(type, m)
    /*
    #define MY_offsetof(type, m) FIELD_OFFSET(type, m)
    */
  #else
    #define MY_offsetof(type, m) ((size_t)&(((type *)0)->m))
  #endif
#endif

#ifndef MY_container_of

/*
#define MY_container_of(ptr, type, m) container_of(ptr, type, m)
#define MY_container_of(ptr, type, m) CONTAINING_RECORD(ptr, type, m)
#define MY_container_of(ptr, type, m) ((type *)((char *)(ptr) - offsetof(type, m)))
#define MY_container_of(ptr, type, m) (&((type *)0)->m == (ptr), ((type *)(((char *)(ptr)) - MY_offsetof(type, m))))
*/

/*
  GCC shows warning: "perhaps the 'offsetof' macro was used incorrectly"
    GCC 3.4.4 : classes with constructor
    GCC 4.8.1 : classes with non-public variable members"
*/

#define MY_container_of(ptr, type, m) ((type *)((char *)(1 ? (ptr) : &((type *)0)->m) - MY_offsetof(type, m)))


#endif

#define CONTAINER_FROM_VTBL_SIMPLE(ptr, type, m) ((type *)(ptr))

/*
#define CONTAINER_FROM_VTBL(ptr, type, m) CONTAINER_FROM_VTBL_SIMPLE(ptr, type, m)
*/
#define CONTAINER_FROM_VTBL(ptr, type, m) MY_container_of(ptr, type, m)

#define CONTAINER_FROM_VTBL_CLS(ptr, type, m) CONTAINER_FROM_VTBL_SIMPLE(ptr, type, m)
/*
#define CONTAINER_FROM_VTBL_CLS(ptr, type, m) CONTAINER_FROM_VTBL(ptr, type, m)
*/



#ifdef _WIN32

#define CHAR_PATH_SEPARATOR '\\'
#define WCHAR_PATH_SEPARATOR L'\\'
#define STRING_PATH_SEPARATOR "\\"
#define WSTRING_PATH_SEPARATOR L"\\"

#else

#define CHAR_PATH_SEPARATOR '/'
#define WCHAR_PATH_SEPARATOR L'/'
#define STRING_PATH_SEPARATOR "/"
#define WSTRING_PATH_SEPARATOR L"/"

#endif

