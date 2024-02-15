// ------------------------------------------------------------------
//   digest.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "libdeflate_1.19/libdeflate.h"
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
#include "writer.h"
#include "txtheader.h"

#define IS_ADLER (IS_ZIP ? !flag.md5 : z_file->z_flags.adler)
#define IS_MD5 (!IS_ADLER)

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

    Digest digest = ctx->is_adler ? (Digest){ .words = { BGEN32 (ctx_copy.adler_ctx.adler) } }
                                  : md5_finalize (&ctx_copy.md5_ctx);

    if (msg && flag.show_digest) 
        iprintf ("%s finalize %s: %s\n", DIGEST_NAME, msg, digest_display_ex (digest, DD_NORMAL).s);

    return digest;
}

Digest digest_do (STRp(data), bool is_adler, rom show_digest_msg)
{
    Digest digest = is_adler ? (Digest){ .words = { BGEN32 (adler32 (1, STRa(data))) } }
                             : md5_do (STRa(data));

    if (flag.show_digest)
        iprintf ("%s digest_do %s: %s\n", DIGEST_NAME, show_digest_msg, digest_display_ex_(digest, DD_NORMAL, is_adler, false).s);

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
            ctx->is_adler = true;
            ctx->adler_ctx.adler = 1;
            ctx->adler_ctx.initialized = true;
        }
        if (flag.show_digest) before = *ctx;
        ctx->adler_ctx.adler = adler32 (ctx->adler_ctx.adler, STRa(data));
    }

    else {
        if (!ctx->md5_ctx.initialized) {
            ctx->is_adler = false;
            bool save_log = ctx->md5_ctx.log;
            ctx->md5_ctx.log = 0;
            md5_initialize (&ctx->md5_ctx);
            ctx->md5_ctx.log = save_log;
            ctx->md5_ctx.initialized = true;
        }
        if (flag.show_digest) before = *ctx;
        md5_update (&ctx->md5_ctx, STRa(data));
    }

    ctx->bytes_digested += data_len;

    if (flag.show_digest) {
        iprintf ("vb=%10s %s update %s (len=%"PRIu64" so_far=%"PRIu64") 32chars=\"%s\": before=%s after=%s\n", 
                 VB_NAME, DIGEST_NAME, msg, data_len, ctx->bytes_digested, 
                 str_to_printable_(data, MIN_(32, data_len)).s, 
                 digest_display_ex (digest_snapshot (&before, NULL), DD_NORMAL).s, 
                 digest_display_ex (digest_snapshot (ctx, NULL),     DD_NORMAL).s);
    }

    if (ctx->log) {
        FILE *f=fopen (DIGEST_LOG_FILENAME, "ab");
        fwrite (data, data_len, 1, f);
        fclose (f);
    }

    COPY_TIMER (digest);
}

void digest_piz_verify_one_txt_file (unsigned txt_file_i/* 0-based */)
{
    char s[200];  

    // note: in case of --test, we abort in case of a digest error (here or in digest_piz_verify_one_vb)
    //       but in case of decompression, we print an error but allow the decompression to proceed, so
    //       that we can detect where the incorrect reconstruction occurred + reduce damage to minimum

    // since v14, if Alder32, we verify each TxtHeader and VB, but we don't create a cumulative digest for the entire file. 
    // now, we just confirm that all VBs were verified as expected.
    if (VER(14) && z_file->z_flags.adler) {
        
        CompIType comp_i = (flag.deep && txt_file_i >= 1) ? (txt_file_i - 1 + SAM_COMP_FQ00)
                         : (flag.pair && txt_file_i ==1 ) ? FQ_COMP_R2
                         :                                COMP_MAIN;

        uint32_t expected_vbs_verified = sections_get_num_vbs (comp_i);

        ASSERT (z_file->num_vbs_verified == expected_vbs_verified ||  // success
                txt_file->vb_digest_failed,                           // failure already announced
                "Expected to have verified (adler32) all %u VBlocks, but verified %u (txt_file=%s txt_file_i=%u)",
                expected_vbs_verified, z_file->num_vbs_verified, txt_name, txt_file_i);

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
            ABORT ("Error: %s of original file=%s is different than decompressed file=%s (txt_file=%s txt_file_i=%u)\n",
                digest_name(), digest_display (z_file->digest).s, digest_display (decompressed_file_digest).s, txt_name, txt_file_i);
        }

        // if decompressed incorrectly - warn, but still give user access to the decompressed file
        else { 
            piz_digest_failed = true; // inspected by main_genounzip
            WARN ("File integrity error: %s of decompressed file %s is %s, but %s of the original %s file was %s (txt_file=%s txt_file_i=%u)", 
                digest_name(), txt_file->name, digest_display (decompressed_file_digest).s, digest_name(), 
                dt_name (txt_file->data_type), digest_display (z_file->digest).s, txt_name, txt_file_i);
        }
    }
}

static void digest_piz_verify_one_vb (VBlockP vb)
{
    // Compare digest up to this VB transmitted through SectionHeaderVbHeader. If Adler32, it is a stand-alone
    // digest of the VB, and if MD5, it is a commulative digest up to this VB.
    if ((!txt_file->vb_digest_failed || IS_ADLER) && // note: for MD5, we report only the first failed VB, bc the digest is commulative, so all subsequent VBs will fail for sure
        (VER(14) || !flag.unbind)) {                 // note: for files <= v13, we cannot test per-VB digest in unbind mode, because the digests (MD5 and Adler32) are commulative since the beginning of the bound file. However, we still test component-wide digest in piz_verify_digest_one_txt_file.

        // add VB to commulative digests as needed
        Digest single_comp_commulative_digest = (!VER(14) || IS_MD5)      ? digest_snapshot (&z_file->digest_ctx, NULL) : DIGEST_NONE; // up to v13 Adler was commulative too
        Digest multi_comp_commulative_digest  = (!VER(14) && flag.unbind) ? digest_snapshot (&z_file->v13_commulative_digest_ctx, NULL) : DIGEST_NONE;

        Digest piz_digest = (VER(14) && IS_ADLER)     ? vb->digest  // stand-alone digest of this VB
                          : (!VER(14) && flag.unbind) ? multi_comp_commulative_digest
                          :                             single_comp_commulative_digest;

        TEMP_FLAG (quiet, false);

        // warn if VB is bad, but don't exit, so file reconstruction is complete and we can debug it
        if (!digest_recon_is_equal (piz_digest, vb->expected_digest)) { 

            char recon_size_warn[100] = "";
            if (vb->recon_size != vb->txt_data.len) // note: leave vb->txt_data.len 64bit to detect bugs
                sprintf (recon_size_warn, "Expecting: VB_HEADER.recon_size=%u == txt_data.len=%"PRIu64"\n", vb->recon_size, vb->txt_data.len);

            WARN ("reconstructed vblock=%s/%u (vb_line_i=0 -> txt_line_i(1-based)=%"PRIu64" num_lines=%u), (%s=%s) differs from original file (%s=%s).\n%s",
                  comp_name (vb->comp_i), vb->vblock_i, writer_get_txt_line_i (vb, 0), vb->lines.len32,
                  DIGEST_NAME, digest_display (piz_digest).s, 
                  DIGEST_NAME, digest_display (vb->expected_digest).s, 
                  recon_size_warn);

            // case: first bad VB: dump bad VB to disk and maybe exit
            if (!__atomic_test_and_set (&txt_file->vb_digest_failed, __ATOMIC_RELAXED)) { // not WARN_ONCE because we might be genounzipping multiple files - we want to show this for every failed file (see also note in digest_piz_verify_one_txt_file)
                WARN ("Bad reconstructed vblock has been dumped to: %s.gz\n"
                      "To see the same data in the original file:\n"
                      "genozip --biopsy %u -B%u %s%s%s",
                      txtfile_dump_vb (vb, z_name), vb->vblock_i, (unsigned)(segconf.vb_size >> 20), 
                      (Z_DT(SAM) && !z_file->z_flags.has_gencomp) ? "--no-gencomp " : "",
                      (txt_file && txt_file->name) ? filename_guess_original (txt_file) : IS_PIZ ? txtheader_get_txt_filename_from_section().s : "(uncalculable)",
                      SUPPORT);

                if (flag.test) exit_on_error (false); // must be inside the atomic test, otherwise another thread will exit before we completed dumping
            }
        }

        else {
            if (flag.show_digest)
                iprintf ("%s: verfied digest\n", VB_NAME);

            __atomic_fetch_add (&z_file->num_vbs_verified, (uint32_t)1, __ATOMIC_RELAXED); // count VBs verified so far
        }

        RESTORE_FLAG (quiet);
    }
}


// ZIP and PIZ: called by compute thread to calculate MD5 or Adler32 of one VB - possibly serializing VBs using a mutex
bool digest_one_vb (VBlockP vb, bool is_compute_thread, 
                    BufferP data) // if NULL, data digested is vb->txt_data
{
    #define WAIT_TIME_USEC 1000
    #define DIGEST_TIMEOUT (30*60) // 30 min

    if (!data) data = &vb->txt_data;

    bool digestable = gencomp_comp_eligible_for_digest(vb);

    // starting V14, if adler32, we digest each VB stand-alone.
    if (IS_ADLER && (IS_ZIP || VER(14))) {
        if (digestable) {
            vb->digest = digest_do (STRb(*data), IS_ADLER, VB_NAME);
            
            if (IS_PIZ) digest_piz_verify_one_vb (vb);  
        }
    }

    else {
        // serialize VBs in order. note: we don't serialize when called from writer, as writer already serializes
        // the output in the correct order, but digestable=false VBs (e.g. PRIM and DEPN in SAM) are called in 
        // the order dictated by PLAN_END_OF_VB in the recon_plan which is not neccesarily their sequential order.
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

            if (!VER(14) && flag.unbind)
                digest_update (&z_file->v13_commulative_digest_ctx, data, "vb"); // up to v13, multi-component, VB digest is commulative across all components (in MD5 and Adler)

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
        z_file->digest_ctx.log = true;
        file_remove (DIGEST_LOG_FILENAME, true);    
    }

    if (!data->len) return DIGEST_NONE;

    Digest digest;
    
    // starting V14, if adler32, we digest each TXT_HEADER stand-alone.
    if (IS_ADLER && VER(14))
        digest = digest_do (STRb(*data), IS_ADLER, "TXT_HEADER");

    // if MD5 or v13 (or earlier) - digest is for commulative for the whole file, and we take a snapshot
    else {
        digest_update (&z_file->digest_ctx, data, "txt_header");
        digest = digest_snapshot (&z_file->digest_ctx, NULL);

        // up to v13, in unbind, VB digest relied on commulative digest of all componets
        if (!VER(14) && flag.unbind)
            digest_snapshot (&z_file->v13_commulative_digest_ctx, NULL);
    }

    if (IS_PIZ) {  
        TEMP_FLAG (quiet, flag.quiet && !flag.show_digest);

        if (!digest_recon_is_equal (digest, piz_expected_digest)) {
            
            WARN ("%s of reconstructed %s header (%s) differs from original file (%s)\n"
                  "Bad reconstructed header has been dumped to: %s\n"
                  "To see the same data in the original file:\n"
                  "genozip --biopsy 0 %s%s",
                  digest_name(),
                  dt_name (z_file->data_type), digest_display (digest).s, digest_display (piz_expected_digest).s,
                  txtfile_dump_vb (data->vb, z_name), filename_guess_original (txt_file), SUPPORT);

            if (flag.test) exit_on_error(false);
        }

        RESTORE_FLAG (quiet);
    }

    COPY_TIMER_EVB (digest_txt_header);

    return digest;
}

void digest_verify_ref_is_equal (const Reference ref, rom header_ref_filename, const Digest header_ref_genome_digest)
{
    Digest ref_digest   = ref_get_genome_digest (ref);
    uint8_t ref_version = ref_get_genozip_version (ref);
    bool is_adler       = ref_is_digest_adler (ref);

    // in v6-v11 we had a bug where the 2nd 32b of an MD5 digest was a copy of the first 32b.
    // we don't know the genozip version of the reference file used for compressing - therefore we accept both MD5 options.
    // this way, if a reference is re-generated from FASTA using a new genozip, it can be used to decompress an old file
    if (ref_version < 15)
        ASSINP ((ref_digest.words[0] == header_ref_genome_digest.words[0] && 
                 ref_digest.words[2] == header_ref_genome_digest.words[2] && 
                 ref_digest.words[3] == header_ref_genome_digest.words[3]),
                "%s: Incorrect reference file:\n%s (%s) was used for compressing\nNot %s (%s version=%u)",
                z_name, header_ref_filename, digest_display_ex_(header_ref_genome_digest, DD_NORMAL, unknown, true).s, 
                ref_get_filename (ref), digest_display_ex_(ref_digest, DD_NORMAL, no, true).s, ref_version);

    // reference files produced since v15 can be either Adler or MD5, and digest is of in-memory genome, not FASTA file
    else {
        if (!digest_is_equal (header_ref_genome_digest, ref_digest)) { 
            // cannot decompress file compressed with Genozip <= v14 with a reference generated with >= v15, because
            // in the former uses a reference file with a the digest is of the FASTA, and in references >= v15 the digest is of the in-memory genome
            ASSINP (VER(15), "Error: The reference file %s was generated with Genozip version %u, but the file %s was compressed with a reference file generated by an earlier version of Genozip. To decompress this file, use a reference file generated with Genozip <= 14.%s",
                    ref_get_filename (ref), ref_version, z_file->basename, 
                    cond_str (!VER(14), " ", "Note: if the reference file used to compress was generated with Genozip version 13.0.20 or earlier, then the original reference file must be used and it is not possible to re-generate it."));

            ASSINP (!VER(15), 
                    "%s: Incorrect reference file:\n%s with %s was used for compressing\nBut %s has %s (ref_genozip_version=%u)",
                    z_name, header_ref_filename, digest_display_ex_(header_ref_genome_digest, DD_NORMAL, unknown, true).s, 
                    ref_get_filename (ref), digest_display_ex_(ref_digest, DD_NORMAL, is_adler, true).s, ref_version);
        }
    }
}

DigestDisplay digest_display_ex_(const Digest digest, DigestDisplayMode mode, thool is_adler, bool show_alg)
{
    DigestDisplay dis = {};

    bytes b = digest.bytes; 
    
    if (is_adler == unknown)
        is_adler = (!digest.words[1] && !digest.words[2] && !digest.words[3]); // if not known - guess

    rom alg = show_alg ? digest_name_(is_adler) : "";
    rom alg_eq = show_alg ? "=" : "";

    if ((mode == DD_NORMAL && !is_adler && digest_is_zero (digest)) ||
        (mode == DD_MD5                 && digest_is_zero (digest)) || 
        (mode == DD_MD5_IF_MD5 && (is_adler || digest_is_zero (digest)))) 
        sprintf (dis.s, "N/A                                    ");

    else if ((mode == DD_NORMAL || mode == DD_SHORT) && is_adler)
        sprintf (dis.s, "%s%s%-10u", alg, alg_eq, BGEN32 (digest.adler_bgen));

    else if ((mode == DD_NORMAL && !is_adler     && !digest_is_zero (digest)) ||
             (mode == DD_MD5                     && !digest_is_zero (digest)) || 
             (mode == DD_MD5_IF_MD5 && !is_adler && !digest_is_zero (digest)))
        sprintf (dis.s, "%s%s%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x", 
                 alg, alg_eq,
                 b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12], b[13], b[14], b[15]);    
    
    else if (mode == DD_SHORT && !is_adler)
        sprintf (dis.s, "%s%s%2.2x%2.2x%2.2x%2.2x", alg, alg_eq, b[0], b[1], b[2], b[3]);    

    return dis;
}

DigestDisplay digest_display_ex (const Digest digest, DigestDisplayMode mode)
{
    return digest_display_ex_(digest, mode, IS_ADLER, false);
}

DigestDisplay digest_display (const Digest digest) 
{ 
    return digest_display_ex (digest, DD_NORMAL); 
}

DigestDisplay digest_display_(const Digest digest, bool is_adler) 
{ 
    return digest_display_ex_(digest, DD_NORMAL, is_adler, false); 
}

rom digest_name (void) { return DIGEST_NAME; }
rom digest_name_(bool is_adler) { return is_adler ? "Adler32" : "MD5"; }

