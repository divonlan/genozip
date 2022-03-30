// ------------------------------------------------------------------
//   digest.c
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "libdeflate/libdeflate.h"
#include "genozip.h"
#include "digest.h"
#include "buffer.h"
#include "md5.h"
#include "vblock.h"
#include "mutex.h"
#include "file.h"
#include "codec.h"
#include "txtfile.h"
#include "gencomp.h"
#include "strings.h"
#include "endianness.h"
#include "profiler.h"

static Mutex vb_digest_mutex = {}; // ZIP: used for serializing MD5ing of VBs

#define IS_ADLER ((command == ZIP) ? !flag.md5 : z_file->z_flags.adler)
#define DIGEST_NAME (IS_ADLER ? "Adler32" : "MD5")

#define DIGEST_LOG_FILENAME (command==ZIP ? "digest.zip.log" : "digest.piz.log")

void digest_initialize (void)
{
    // reset if we're starting a new z_file
    if (!z_file->num_txts_so_far) 
        mutex_initialize (vb_digest_mutex); 
}

// get digest of data so far
Digest digest_snapshot (const DigestContext *ctx, rom msg)
{
    // make a copy of ctx, and finalize it, keeping the original copy unfinalized
    DigestContext ctx_copy = *ctx;

    Digest digest = IS_ADLER ? (Digest) { .words = { BGEN32 (ctx_copy.adler_ctx.adler) } }
                             : md5_finalize (&ctx_copy.md5_ctx);

    if (msg && flag.show_digest) 
        iprintf ("%s finalize %s: %s\n", DIGEST_NAME, msg, digest_display_ex (digest, DD_NORMAL).s);

    return digest;
}

void digest_start_log (DigestContext *ctx)
{
    ctx->common.log = true;
    file_remove (DIGEST_LOG_FILENAME, true);
}

Digest digest_do (const void *data, uint32_t len)
{
    if (IS_ADLER)
        return (Digest) { .words = { BGEN32 (adler32 (1, data, len)) } };
    else
        return md5_do (data, len);
}

void digest_update_do (VBlockP vb, DigestContext *ctx, rom data, uint64_t data_len, rom msg)
{
    if (!data_len) return;
    
    START_TIMER;
    DigestContext before;

    if (IS_ADLER) {
        if (!ctx->adler_ctx.initialized) {
            ctx->adler_ctx.adler = 1;
            ctx->adler_ctx.initialized = true;
        }
        if (flag.show_digest) before = *ctx;
        ctx->adler_ctx.adler = adler32 (ctx->adler_ctx.adler, STRa(data));
    }

    else {
        if (!ctx->md5_ctx.initialized) {
            bool save_log = ctx->md5_ctx.log;
            ctx->md5_ctx.log = 0;
            md5_initialize (&ctx->md5_ctx);
            ctx->md5_ctx.log = save_log;
            ctx->md5_ctx.initialized = true;
        }
        if (flag.show_digest) before = *ctx;
        md5_update (&ctx->md5_ctx, STRa(data));
    }

    ctx->common.bytes_digested += data_len;

    if (flag.show_digest) {
        char str[65];
        iprintf ("vb=%05d %s update %s (len=%"PRIu64" so_far=%"PRIu64") 32chars=\"%s\": before=%s after=%s\n", 
                 vb->vblock_i, DIGEST_NAME, msg, data_len, ctx->common.bytes_digested, 
                 str_to_printable (data, MIN_(32, data_len), str), 
                 digest_display_ex (digest_snapshot (&before, NULL), DD_NORMAL).s, 
                 digest_display_ex (digest_snapshot (ctx, NULL), DD_NORMAL).s);
    }

    if (ctx->common.log) {
        FILE *f=fopen (DIGEST_LOG_FILENAME, "ab");
        fwrite (data, data_len, 1, f);
        fclose (f);
    }

    COPY_TIMER_VB (vb, digest);
}

void digest_piz_verify_one_vb (VBlockP vb)
{
    static bool failed = false; 

    // if testing, compare digest up to this VB to that calculated on the original file and transferred through SectionHeaderVbHeader
    // note: we cannot test this unbind mode, because the digests are commulative since the beginning of the bound file
    if (!failed && !flag.unbind && !v8_digest_is_zero (vb->digest_so_far)) {
        Digest piz_digest_so_far = digest_snapshot (&z_file->digest_ctx, NULL);

        // warn if VB is bad, but don't exit, so file reconstruction is complete and we can debug it
        if (!digest_recon_is_equal (piz_digest_so_far, vb->digest_so_far) && !digest_is_equal (vb->digest_so_far, DIGEST_NONE)) {

            TEMP_FLAG (quiet, flag.quiet && !flag.show_digest);

            // dump bad vb to disk
            WARN ("reconstructed vblock=%u, component=%u (%s=%s) differs from original file (%s=%s).\n"
                  "Note: genounzip is unable to check the %s subsequent vblocks once a vblock is bad\n"
                  "Bad reconstructed vblock has been dumped to: %s\n"
                  "To see the same data in the original file:\n"
                  "genozip --biopsy %u %s (+any parameters used to compress this file)\n"
                  "If this is unexpected, please contact support@genozip.com.\n", 
                  vb->vblock_i, vb->comp_i,  
                  DIGEST_NAME, digest_display (piz_digest_so_far).s, 
                  DIGEST_NAME, digest_display (vb->digest_so_far).s, 
                  DIGEST_NAME, txtfile_dump_vb (vb, z_name), vb->vblock_i, file_guess_original_filename (txt_file));

            RESTORE_FLAG (quiet);

            failed = true; // no point in test the rest of the vblocks as they will all fail - MD5 is commulative
        }
    }
}

// ZIP and PIZ: called by compute thread to calculate MD5 for one VB - need to serialize VBs using a mutex
void digest_one_vb (VBlockP vb)
{
    #define WAIT_TIME_USEC 1000
    #define DIGEST_TIMEOUT (30*60) // 30 min

    // serialize VBs in order
    mutex_lock_by_vb_order (vb->vblock_i, vb_digest_mutex);

    // note: we don't digest generated components, but we need to enter the mutex anyway as it serializes VBs in order
    if (gencomp_comp_eligible_for_digest(vb)) {
    
        if (command == ZIP) {
            digest_update (&z_file->digest_ctx, &vb->txt_data, "vb");

            // take a snapshot of digest as per the end of this VB - this will be used to test for errors in piz after each VB  
            vb->digest_so_far = digest_snapshot (&z_file->digest_ctx, NULL);
        }
        
        else { // PIZ
            digest_update (&z_file->digest_ctx, &vb->txt_data, "vb"); // labels consistent with ZIP so we can easily diff PIZ vs ZIP

            digest_piz_verify_one_vb (vb);        
        }
    }

    mutex_unlock (vb_digest_mutex);
}

// in v6-v11 we had a bug were the 2nd 32b of an MD5 digest was a copy of the first 32b
bool digest_recon_is_equal (const Digest recon_digest, const Digest expected_digest) 
{
    if (z_file->genozip_version >= 12) return digest_is_equal (recon_digest, expected_digest);

    return recon_digest.words[0] == expected_digest.words[0] &&
           (recon_digest.words[0] == expected_digest.words[1] /* buggy md5 */ || !recon_digest.words[1] /* adler */) &&
           recon_digest.words[2] == expected_digest.words[2] &&
           recon_digest.words[3] == expected_digest.words[3]; 
}

void digest_verify_ref_is_equal (const Reference ref, rom header_ref_filename, const Digest header_md5)
{
    Digest ref_md5 = ref_get_file_md5 (ref);
    uint8_t version = ref_get_genozip_version (ref);

    // in v6-v11 we had a bug where the 2nd 32b of an MD5 digest was a copy of the first 32b.
    // we don't know the genozip version of the reference file used for compressing - therefore we accept both MD5 options.
    // this way, if a reference is re-generated from FASTA using a new genozip, it can be used to decompress an old file
    ASSINP (//digest_is_equal (ref_md5, header_md5) ||
            (ref_md5.words[0] == header_md5.words[0] && ref_md5.words[2] == header_md5.words[2] && ref_md5.words[3] == header_md5.words[3]),
            "%s: Bad reference file:\n%s (MD5=%s) was used for compressing\n%s (MD5=%s version=%u) has a different MD5",
            z_name, header_ref_filename, digest_display_ex (header_md5, DD_MD5).s, 
            ref_get_filename (ref), digest_display_ex (ref_md5, DD_MD5).s, version);
}

DigestDisplay digest_display_ex (const Digest digest, DigestDisplayMode mode)
{
    DigestDisplay dis = {};

    bytes b = digest.bytes; 
    
    if ((mode == DD_NORMAL && !IS_ADLER && md5_is_zero (digest)) ||
        (mode == DD_MD5                 && md5_is_zero (digest)) || 
        (mode == DD_MD5_IF_MD5 && (IS_ADLER || md5_is_zero (digest)))) 
        sprintf (dis.s, "N/A                             ");

    else if ((mode == DD_NORMAL || mode == DD_SHORT) && IS_ADLER)
        sprintf (dis.s, "%-10u", BGEN32 (digest.adler_bgen));

    else if ((mode == DD_NORMAL && !IS_ADLER && !md5_is_zero (digest)) ||
        (mode == DD_MD5                 && !md5_is_zero (digest)) || 
        (mode == DD_MD5_IF_MD5 && !IS_ADLER && !md5_is_zero (digest)))
        sprintf (dis.s, "%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x", 
                 b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12], b[13], b[14], b[15]);    
    
    else if (mode == DD_SHORT && !IS_ADLER)
        sprintf (dis.s, "%2.2x%2.2x%2.2x%2.2x", b[0], b[1], b[2], b[3]);    

    return dis;
}
DigestDisplay digest_display (const Digest digest) { return digest_display_ex (digest, DD_NORMAL); }

rom digest_name (void) { return DIGEST_NAME; }
