/* LzFind.h -- Match finder for LZ algorithms
2017-06-10 : Igor Pavlov : Public domain */

#pragma once

#include "7zTypes.h"

typedef UInt32 CLzRef;

typedef struct _CMatchFinder
{
  Byte *buffer;
  UInt32 pos;
  UInt32 posLimit;
  UInt32 streamPos;
  UInt32 lenLimit;

  UInt32 cyclicBufferPos;
  UInt32 cyclicBufferSize; /* it must be = (historySize + 1) */

  Byte streamEndWasReached;
  Byte btMode;
  Byte bigHash;
  Byte directInput;

  UInt32 matchMaxLen;
  CLzRef *hash;
  CLzRef *son;
  UInt32 hashMask;
  UInt32 cutValue;

  Byte *bufferBase;
  ISeqInStream *stream;
  
  UInt32 blockSize;
  UInt32 keepSizeBefore;
  UInt32 keepSizeAfter;

  UInt32 numHashBytes;
  size_t directInputRem;
  UInt32 historySize;
  UInt32 fixedHashSize;
  UInt32 hashSizeSum;
  SRes result;
  UInt32 crc[256];
  size_t numRefs;

  UInt64 expectedDataSize;
} CMatchFinder;

#define Inline_MatchFinder_GetPointerToCurrentPos(p) ((p)->buffer)

#define Inline_MatchFinder_GetNumAvailableBytes(p) ((p)->streamPos - (p)->pos)

#define Inline_MatchFinder_IsFinishedOK(p) \
    ((p)->streamEndWasReached \
        && (p)->streamPos == (p)->pos \
        && (!(p)->directInput || (p)->directInputRem == 0))

void MatchFinder_Construct(CMatchFinder *p);

/* Conditions:
     historySize <= 3 GB
     keepAddBufferBefore + matchMaxLen + keepAddBufferAfter < 511MB
*/
int MatchFinder_Create(VBlockP vb, CMatchFinder *p, UInt32 historySize,
    UInt32 keepAddBufferBefore, UInt32 matchMaxLen, UInt32 keepAddBufferAfter);
void MatchFinder_Free(VBlockP vb, CMatchFinder *p);

/*
Conditions:
  Mf_GetNumAvailableBytes_Func must be called before each Mf_GetMatchLen_Func.
  Mf_GetPointerToCurrentPos_Func's result must be used only before any other function
*/

typedef void (*Mf_Init_Func)(void *object);
typedef UInt32 (*Mf_GetNumAvailableBytes_Func)(void *object);
typedef const Byte * (*Mf_GetPointerToCurrentPos_Func)(void *object);
typedef UInt32 (*Mf_GetMatches_Func)(void *object, UInt32 *distances);
typedef void (*Mf_Skip_Func)(void *object, UInt32);

typedef struct _IMatchFinder
{
  Mf_Init_Func Init;
  Mf_GetNumAvailableBytes_Func GetNumAvailableBytes;
  Mf_GetPointerToCurrentPos_Func GetPointerToCurrentPos;
  Mf_GetMatches_Func GetMatches;
  Mf_Skip_Func Skip;
} IMatchFinder;

void MatchFinder_CreateVTable(CMatchFinder *p, IMatchFinder *vTable);

