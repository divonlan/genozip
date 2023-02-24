// ------------------------------------------------------------------
//   digest.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "libdeflate/libdeflate.h"
#include "genozip.h"
#include "digest.h"
#include "buffer.h"
#include "md5.h"
#include "vblock.h"
#include "mutex.h"
#include "filename.h"
#include "file.h"
#include "codec.h"
#include "txtfile.h"
#include "gencomp.h"
#include "strings.h"
#include "endianness.h"
#include "profiler.h"
#include "progress.h"

#define IS_ADLER (IS_ZIP ? !flag.md5 : z_file->z_flags.adler)
#define DIGEST_NAME (IS_ADLER ? "Adler32" : "MD5")

#define DIGEST_LOG_FILENAME (command==ZIP ? "digest.zip.log" : "digest.piz.log")

static bool digest_recon_is_equal (const Digest recon_digest, const Digest expected_digest) 
{
    if (VER(12)) return digest_is_equal (recon_digest, expected_digest);

    // in v6-v11 we had a bug were the 2nd 32b of an MD5 digest was a copy of the first 32b
    return recon_digest.words[0] == expected_digest.words[0] &&
           (recon_digest.words[0] == expected_digest.words[1] /* buggy md5 */ || !recon_digest.words[1] /* adler */) &&
           recon_digest.words[2] == expected_digest.words[2] &&
           recon_digest.words[3] == expected_digest.words[3]; 
}

static bool piz_digest_failed = false;
bool digest_piz_has_it_failed (void) { return piz_digest_failed; }

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

static Digest digest_do (const void *data, uint32_t len, rom show_digest_msg)
{
    Digest digest = IS_ADLER ? (Digest){ .words = { BGEN32 (adler32 (1, data, len)) } }
                             : md5_do (data, len);

    if (flag.show_digest)
        iprintf ("%s digest_one_vb %s: %s\n", DIGEST_NAME, show_digest_msg, digest_display_ex (digest, DD_NORMAL).s);

    return digest;
}

#define digest_update(ctx, buf, msg) digest_update_do ((buf)->vb, (ctx), STRb(*(buf)), (msg))
static void digest_update_do (VBlockP vb, DigestContext *ctx, rom data, uint64_t data_len, rom msg)
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
        iprintf ("vb=%10s %s update %s (len=%"PRIu64" so_far=%"PRIu64") 32chars=\"%s\": before=%s after=%s\n", 
                 VB_NAME, DIGEST_NAME, msg, data_len, ctx->common.bytes_digested, 
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

void digest_piz_verify_one_txt_file (unsigned txt_file_i/* 0-based */)
{
    char s[200];  

    // since v14, if Alder32, we verify each TxtHeader and VB, but we don't create a cumulative digest for the entire file. 
    // now, we just confirm that all VBs were verified as expected.
    if (VER(14) && z_file->z_flags.adler) {
        
        CompIType comp_i = (flag.deep && txt_file_i==1)             ? SAM_COMP_FQ00
                         : (flag.deep && txt_file_i==2)             ? SAM_COMP_FQ01
                         : (fastq_piz_is_paired() && txt_file_i==1) ? FQ_COMP_R2
                         :                                            COMP_MAIN;

        uint32_t expected_vbs_verified = sections_get_num_vbs (comp_i);

        ASSERT (z_file->num_vbs_verified == expected_vbs_verified ||  // success
                txt_file->vb_digest_failed,                           // failure already announced
                "Expected to have verified (adler32) all %u VBlocks, but verified only %u",
                expected_vbs_verified, z_file->num_vbs_verified);

        if (flag.show_digest)
            iprintf ("Txt file #%u: %u VBs verified\n", txt_file_i, z_file->num_vbs_verified);

        if (flag.test) { 
            sprintf (s, "verified as identical to the original %s", dt_name (txt_file->data_type));
            progress_finalize_component (s); 
        }

        z_file->num_vbs_verified = 0; // reset for next component
    }

    // txt file commulative digest: MD5, and up to v13, also Adler32
    else {
        Digest decompressed_file_digest = digest_snapshot (&z_file->digest_ctx, "file"); 

        if (digest_recon_is_equal (decompressed_file_digest, z_file->digest)) {
            if (flag.test) { 
                sprintf (s, "verified as identical to the original %s (%s=%s)", 
                        dt_name (txt_file->data_type), digest_name(), digest_display (decompressed_file_digest).s);
                progress_finalize_component (s); 
            }
        }
        
        else if (flag.test) {
            progress_finalize_component ("FAILED!");
            ABORT ("Error: %s of original file=%s is different than decompressed file=%s\n",
                digest_name(), digest_display (z_file->digest).s, digest_display (decompressed_file_digest).s);
        }

        // if decompressed incorrectly - warn, but still give user access to the decompressed file
        else { 
            piz_digest_failed = true; // inspected by main_genounzip
            WARN ("File integrity error: %s of decompressed file %s is %s, but %s of the original %s file was %s", 
                digest_name(), txt_file->name, digest_display (decompressed_file_digest).s, digest_name(), 
                dt_name (txt_file->data_type), digest_display (z_file->digest).s);
        }
    }
}

static void digest_piz_verify_one_vb (VBlockP vb)
{
    // Compare digest up to this VB transmitted through SectionHeaderVbHeader. If Adler32, it is a stand-alone
    // digest of the VB, and if MD5, it is a commulative digest up to this VB.
    if ((!txt_file->vb_digest_failed || IS_ADLER) && // note: for MD5, we report only the first failed VB, bc the digest is commulative, so all subsequent VBs will fail for sure
        (VER(14) || !flag.unbind)) {                 // note: for files <= v13, we cannot test per-VB digest in unbind mode, because the digests (MD5 and Adler32) are commulative since the beginning of the bound file. However, we still test component-wide digest in piz_verify_digest_one_txt_file.

        Digest piz_digest = (VER(14) && IS_ADLER) ? vb->digest  // stand-alone digest of this VB
                                                  : digest_snapshot (&z_file->digest_ctx, NULL); // commulative digest so far

        // warn if VB is bad, but don't exit, so file reconstruction is complete and we can debug it
        if (!digest_recon_is_equal (piz_digest, vb->expected_digest)) { 

            TEMP_FLAG (quiet, flag.quiet && !flag.show_digest);

            char recon_size_warn[100] = "";
            if (vb->recon_size != vb->txt_data.len)
                sprintf (recon_size_warn, "Expecting: VB_HEADER.recon_size=%u == txt_data.len=%"PRIu64"\n", vb->recon_size, vb->txt_data.len);

            WARN ("reconstructed vblock=%s/%u, (%s=%s) differs from original file (%s=%s).\n%s",
                  comp_name (vb->comp_i), vb->vblock_i, 
                  DIGEST_NAME, digest_display (piz_digest).s, 
                  DIGEST_NAME, digest_display (vb->expected_digest).s, 
                  recon_size_warn);

            // case: first bad VB: dump bad VB to disk
            if (!txt_file->vb_digest_failed)
                WARN ("Bad reconstructed vblock has been dumped to: %s.gz\n"
                      "To see the same data in the original file:\n"
                      "genozip --biopsy %u %s (+any parameters used to compress this file)\n"
                      "If this is unexpected, please contact support@genozip.com.\n", 
                    txtfile_dump_vb (vb, z_name), vb->vblock_i, filename_guess_original (txt_file));

            if (flag.test) exit_on_error (false);

            RESTORE_FLAG (quiet);

            txt_file->vb_digest_failed = true;
        }

        else
            __atomic_fetch_add (&z_file->num_vbs_verified, (uint32_t)1, __ATOMIC_RELAXED); // count VBs verified so far
    }
}


// ZIP and PIZ: called by compute thread to calculate MD5 or Adler32 of the txt header for one VB - possibly serializing VBs using a mutex
bool digest_one_vb (VBlockP vb, bool is_compute_thread, 
                    BufferP data) // if NULL, data digested is vb->txt_data
{
    #define WAIT_TIME_USEC 1000
    #define DIGEST_TIMEOUT (30*60) // 30 min

    if (!data) data = &vb->txt_data;

    bool digestable = gencomp_comp_eligible_for_digest(vb);

    // starting V14, if adler32, we digest each VB stand-alone.
    if (IS_ADLER && VER(14)) {
        if (digestable) {
            vb->digest = digest_do (STRb(*data), VB_NAME);
            
            if (IS_PIZ) digest_piz_verify_one_vb (vb);  
        }
    }

    else {
        // serialize VBs in order
        if (is_compute_thread)
            serializer_lock (z_file->digest_serializer, vb->vblock_i);

        // note: we don't digest generated components, but we need to enter the mutex anyway as otherwise the serialization will go out of sync
        if (digestable && IS_ZIP) {
            digest_update (&z_file->digest_ctx, data, "vb");

            // take a snapshot of the commulative digest as per the end of this VB 
            vb->digest = digest_snapshot (&z_file->digest_ctx, NULL);
        }
            
        else if (digestable && IS_PIZ) {
            digest_update (&z_file->digest_ctx, data, "vb"); // labels consistent with ZIP so we can easily diff PIZ vs ZIP

            digest_piz_verify_one_vb (vb);        
        }

        if (is_compute_thread)
            serializer_unlock (z_file->digest_serializer);
    }

    return digestable;
}

// ZIP and PIZ: called by main thread to calculate MD5 or Adler32 of the txt header
Digest digest_txt_header (BufferP data, Digest piz_expected_digest)
{
    START_TIMER;

    // start dogest log if need
    if (flag.log_digest) {
        z_file->digest_ctx.common.log = true;
        file_remove (DIGEST_LOG_FILENAME, true);    
    }

    if (!data->len) return DIGEST_NONE;

    Digest digest;
    
    // starting V14, if adler32, we digest each TXT_HEADER stand-alone.
    if (IS_ADLER && VER(14))
        digest = digest_do (STRb(*data), "TXT_HEADER");

    // if MD5 or v13 (or earlier) - digest is for commulative for the whole file, and we take a snapshot
    else {
        digest_update (&z_file->digest_ctx, data, "txt_header");
        digest = digest_snapshot (&z_file->digest_ctx, NULL);
    }

    if (IS_PIZ) {  
        TEMP_FLAG (quiet, flag.quiet && !flag.show_digest);

        if (!digest_recon_is_equal (digest, piz_expected_digest)) {
            
            WARN ("%s of reconstructed %s header (%s) differs from original file (%s)\n"
                  "Bad reconstructed header has been dumped to: %s\n"
                  "To see the same data in the original file:\n"
                  "genozip --biopsy 0 %s\n"
                  "If this is unexpected, please contact support@genozip.com.\n", 
                  digest_name(),
                  dt_name (z_file->data_type), digest_display (digest).s, digest_display (piz_expected_digest).s,
                  txtfile_dump_vb (data->vb, z_name), filename_guess_original (txt_file));

            if (flag.test) exit_on_error(false);
        }

        RESTORE_FLAG (quiet);
    }

    COPY_TIMER_VB (evb, digest_txt_header);

    return digest;
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
    
    if ((mode == DD_NORMAL && !IS_ADLER && digest_is_zero (digest)) ||
        (mode == DD_MD5                 && digest_is_zero (digest)) || 
        (mode == DD_MD5_IF_MD5 && (IS_ADLER || digest_is_zero (digest)))) 
        sprintf (dis.s, "N/A                             ");

    else if ((mode == DD_NORMAL || mode == DD_SHORT) && IS_ADLER)
        sprintf (dis.s, "%-10u", BGEN32 (digest.adler_bgen));

    else if ((mode == DD_NORMAL && !IS_ADLER && !digest_is_zero (digest)) ||
        (mode == DD_MD5                 && !digest_is_zero (digest)) || 
        (mode == DD_MD5_IF_MD5 && !IS_ADLER && !digest_is_zero (digest)))
        sprintf (dis.s, "%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x", 
                 b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12], b[13], b[14], b[15]);    
    
    else if (mode == DD_SHORT && !IS_ADLER)
        sprintf (dis.s, "%2.2x%2.2x%2.2x%2.2x", b[0], b[1], b[2], b[3]);    

    return dis;
}
DigestDisplay digest_display (const Digest digest) { return digest_display_ex (digest, DD_NORMAL); }

rom digest_name (void) { return DIGEST_NAME; }

