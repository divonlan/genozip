// ------------------------------------------------------------------
//   digest.h
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

typedef packed_enum { DIGEST_MD5=0, DIGEST_ADLER=1, DIGEST_XXH3=2 } DigestAlg; // part of the file format: GenozipHeader.ref_genome_digest_alg (compatable with "bool is_adler" of previous versions: 1=Adler32 0=MD5)
#define DIGEST_UNKNOWN 127 // not part of the enum, and never in DigestStateHeader.alg

typedef union { 
    struct { uint32_t adler32, unused32[3]; };
    struct { uint64_t xxh3, unused64; }; 
    uint8_t md5[16]; 
    uint32_t words[4];
    uint128_t w128;  
} Digest;

#define DIGEST_NONE ((Digest){})

typedef struct {
    uint64_t bytes_digested : 60;
    DigestAlg alg           : 2; 
    bool initialized        : 1;
    bool log                : 1;
} DigestStateHeader;

typedef struct { // 96 bytes
    DigestStateHeader;
    uint64_t     n_bytes;
    uint32_t     a, b, c, d;
    union {
        uint8_t  bytes[64];
        uint32_t words[16];
    };
} Md5State, *Md5StateP;

typedef struct {
    DigestStateHeader;
    uint32_t adler; // native endianity
} AdlerState;

typedef union {
    DigestStateHeader;
    Md5State md5_state;
    AdlerState adler_state;
} DigestState;

#define DIGEST_CONTEXT_NONE (DigestState){}

extern Digest digest_snapshot (const DigestState *ctx, rom msg);
extern Digest digest_txt_header (BufferP data, Digest piz_expected_digest, CompIType comp_i);
extern bool digest_one_vb (VBlockP vb, bool is_compute_thread, BufferP data);
extern void digest_piz_verify_one_txt_file (unsigned txt_file_i);
extern bool digest_piz_has_it_failed (void);

extern Digest digest_do (rom data, uint64_t data_len, DigestAlg alg, rom show_digest_msg);

typedef struct { char s[64]; } DigestDisplay;
extern DigestDisplay digest_display (Digest digest);
extern DigestDisplay digest_display_(Digest digest, DigestAlg alg);

typedef enum { DD_NORMAL, DD_MD5, DD_MD5_IF_MD5, DD_SHORT } DigestDisplayMode;
extern DigestDisplay digest_display_ex_(const Digest digest, DigestDisplayMode mode, DigestAlg alg, bool show_alg);
extern DigestDisplay digest_display_ex (const Digest digest, DigestDisplayMode mode);
extern rom digest_name (void);
extern rom digest_alg_name(DigestAlg alg);

#define digest_is_equal(digest1,digest2) ((digest1).w128 == (digest2).w128)
extern void digest_verify_correct_external_reference_is_loaded (void);

#define digest_is_zero(digest) (!(digest).w128)

static inline DigestAlg digest_guess_alg (Digest digest)
{
    return digest.unused64    ? DIGEST_MD5    
         : digest.unused32[0] ? DIGEST_XXH3   // 2^-64 chance of guessing XXH3 when it is actually MD5
         :                      DIGEST_ADLER; // 2^-32 chance of guessing ADLER when it is actually XXH3 
}

// backward compatability note: in v8 files compressed without --md5 or --test, we had no digest. starting v9, we have Adler32
#define piz_need_digest (z_file->z_flags.has_digest && !flag.piz_txt_modified && !flag_loading_auxiliary) 

#define zip_need_digest (!flag.make_reference && !zip_is_biopsy)
