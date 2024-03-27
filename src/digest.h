// ------------------------------------------------------------------
//   digest.h
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

typedef union { 
    // used if the digest is MD5
    uint8_t bytes[16]; 
    uint32_t words[4];
    uint128_t w128;

    // used if the digest is Adler32
#define adler_bgen words[0] // Big Endian
} Digest;

#define DIGEST_NONE ((Digest){})

typedef struct {
    bool is_adler; // always false
    bool initialized;
    bool log;
    uint64_t bytes_digested;
    uint32_t     lo, hi;
    uint32_t     a, b, c, d;
    union {
        uint8_t  bytes[64];
        uint32_t words[16];
    } buffer;
} Md5Context;

typedef struct {
    bool is_adler; // always true
    bool initialized;
    bool log;
    uint64_t bytes_digested;
    uint32_t adler; // native endianity
} AdlerContext;

typedef union {
    struct {
        bool is_adler;
        bool initialized;
        bool log;
        uint64_t bytes_digested; 
    };
    Md5Context md5_ctx;
    AdlerContext adler_ctx;
} DigestContext;

#define DIGEST_CONTEXT_NONE (DigestContext){}

extern Digest digest_snapshot (const DigestContext *ctx, rom msg);
extern Digest digest_txt_header (BufferP data, Digest piz_expected_digest, CompIType comp_i);
extern bool digest_one_vb (VBlockP vb, bool is_compute_thread, BufferP data);
extern void digest_piz_verify_one_txt_file (unsigned txt_file_i);
extern bool digest_piz_has_it_failed (void);

extern Digest digest_do (STRp(data), bool is_adler, rom show_digest_msg);

typedef struct { char s[64]; } DigestDisplay;
extern DigestDisplay digest_display (Digest digest);
extern DigestDisplay digest_display_(Digest digest, bool is_adler);

typedef enum { DD_NORMAL, DD_MD5, DD_MD5_IF_MD5, DD_SHORT } DigestDisplayMode;
extern DigestDisplay digest_display_ex_(const Digest digest, DigestDisplayMode mode, thool is_adler, bool show_alg);
extern DigestDisplay digest_display_ex (const Digest digest, DigestDisplayMode mode);
extern rom digest_name (void);
extern rom digest_name_(bool is_adler);

#define digest_is_equal(digest1,digest2) ((digest1).w128 == (digest2).w128)
extern void digest_verify_ref_is_equal (const Reference ref, rom header_ref_filename, const Digest header_md5);

#define digest_is_zero(digest) (!(digest).w128)

// backward compatability note: in v8 files compressed without --md5 or --test, we had no digest. starting v9, we have Adler32
#define piz_need_digest (z_file->z_flags.has_digest && !flag.piz_txt_modified && !flag.no_writer && !flag_loading_auxiliary) 

#define zip_need_digest (!flag.make_reference && !flag.zip_txt_modified && !flag.biopsy && flag.no_biopsy_line)
