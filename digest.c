// ------------------------------------------------------------------
//   digest.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
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
#include "strings.h"
#include "endianness.h"

static Mutex vb_digest_mutex = {};   // ZIP: used for serializing MD5ing of VBs
static uint32_t vb_digest_last=0; // last vb to be MD5ed 

#define IS_ADLER ((command == ZIP) ? !flag.md5 : z_file->z_flags.adler)
#define DIGEST_NAME (IS_ADLER ? "Adler32" : "MD5")

void digest_initialize (void)
{
    mutex_initialize (vb_digest_mutex);
    if (!z_file->num_txt_components_so_far) vb_digest_last = 0; // reset if we're starting a new z_file
}

Digest digest_finalize (DigestContext *ctx, const char *msg)
{
    Digest digest = IS_ADLER ? (Digest) { .words = { BGEN32 (ctx->adler_ctx.adler) } }
                             : md5_finalize (&ctx->md5_ctx);

    if (flag.show_digest) 
        iprintf ("%s finalize %s: %s\n", DIGEST_NAME, msg, digest_display_ex (digest, DD_NORMAL).s);

    return digest;
}

// get digest of data so far, without finalizing
Digest digest_snapshot (const DigestContext *ctx)
{
    // make a copy of ctx, and finalize it, keeping the original copy unfinalized
    DigestContext ctx_copy = *ctx;

    return IS_ADLER ? (Digest) { .words = { BGEN32 (ctx_copy.adler_ctx.adler) } }
                    : md5_finalize (&ctx_copy.md5_ctx);
}

Digest digest_do (const void *data, uint32_t len)
{
    if (IS_ADLER)
        return (Digest) { .words = { BGEN32 (libdeflate_adler32 (1, data, len)) } };
    else
        return md5_do (data, len);
}

void digest_update (DigestContext *ctx, const Buffer *buf, const char *msg)
{
    DigestContext before;

    if (IS_ADLER) {
        if (!ctx->adler_ctx.initialized) {
            ctx->adler_ctx.adler = 1;
            ctx->adler_ctx.initialized = true;
        }
        if (flag.show_digest) before = *ctx;
        ctx->adler_ctx.adler = libdeflate_adler32 (ctx->adler_ctx.adler, buf->data, buf->len);
    }

    else {
        if (!ctx->md5_ctx.initialized) {
            md5_initialize (&ctx->md5_ctx);
            ctx->md5_ctx.initialized = true;
        }
        if (flag.show_digest) before = *ctx;
        md5_update (&ctx->md5_ctx, buf->data, buf->len);
    }

    ctx->common.bytes_digested += buf->len;

    if (flag.show_digest) {
        char str[65];
        iprintf ("vb=%05d %s update %s (len=%"PRIu64" so_far=%"PRIu64") 32chars=\"%s\": before=%s after=%s\n", 
                 buf->vb->vblock_i, DIGEST_NAME, msg, buf->len, ctx->common.bytes_digested, 
                 str_to_single_line_printable (buf->data, MIN_(32, (int)buf->len), str), 
                 digest_display_ex (digest_snapshot (&before), DD_NORMAL).s, 
                 digest_display_ex (digest_snapshot (ctx), DD_NORMAL).s);
    }
}

// ZIP and PIZ: called by compute thread to calculate MD5 for one VB - need to serialize VBs using a mutex
void digest_one_vb (VBlock *vb)
{
    #define WAIT_TIME_USEC 1000
    #define DIGEST_TIMEOUT (30*60) // 30 min

    // wait for our turn
    for (unsigned i=0; ; i++) {
        mutex_lock (vb_digest_mutex);

        ASSERT (vb_digest_last < vb->vblock_i, "Expecting vb_digest_last=%u < vb->vblock_i=%u", vb_digest_last, vb->vblock_i);
        
        if (vb_digest_last == vb->vblock_i - 1) break; // its our turn now

        // not our turn, wait 1ms and try again
        mutex_unlock (vb_digest_mutex);
        usleep (WAIT_TIME_USEC);

        // timeout after approx 30 seconds
        ASSERT (i < DIGEST_TIMEOUT*(1000000/WAIT_TIME_USEC), "Timeout (%u sec) while waiting for vb_digest_mutex in vb=%u. vb_digest_last=%u", 
                DIGEST_TIMEOUT, vb->vblock_i, vb_digest_last);
    }

    if (command == ZIP) {
        // digest of all data up to and including this VB is just the total digest of the file so far (as there is no unconsumed data)
        if (flag.bind) digest_update (&z_file->digest_ctx_bound, &vb->txt_data, "vb:digest_ctx_bound");
        digest_update (&z_file->digest_ctx_single, &vb->txt_data, "vb:digest_ctx_single");

        // take a snapshot of digest as per the end of this VB - this will be used to test for errors in piz after each VB  
        vb->digest_so_far = digest_snapshot (flag.bind ? &z_file->digest_ctx_bound : &z_file->digest_ctx_single);
    }
    
    else { // PIZ
        static bool failed = false; // note: when testing multiple files, if a file fails the test, we don't test subsequent files, so no need to reset this variable

        digest_update (&txt_file->digest_ctx_bound, &vb->txt_data, flag.unbind ? "vb:digest_ctx_single" : "vb:digest_ctx_bound"); // labels consistent with ZIP so we can easily diff PIZ vs ZIP

        // if testing, compare digest up to this VB to that calculated on the original file and transferred through SectionHeaderVbHeader
        // note: we cannot test this unbind mode, because the digests are commulative since the beginning of the bound file
        if (!failed && !flag.unbind && !v8_digest_is_zero (vb->digest_so_far)) {
            Digest piz_digest_so_far = digest_snapshot (&txt_file->digest_ctx_bound);

            // warn if VB is bad, but don't exit, so file reconstruction is complete and we can debug it
            if (!digest_recon_is_equal (piz_digest_so_far, vb->digest_so_far) && !digest_is_equal (vb->digest_so_far, DIGEST_NONE)) {

                TEMP_FLAG (quiet, flag.quiet && !flag.show_digest);

                // dump bad vb to disk
                WARN ("%s of reconstructed vblock=%u,component=%u (%s) differs from original file (%s).\n"
                      "Note: genounzip is unable to check the %s subsequent vblocks once a vblock is bad\n"
                      "Bad reconstructed vblock has been dumped to: %s\n"
                      "To see the same data in the original file:\n"
                      "   %s %s | head -c %"PRIu64" | tail -c %u > %s", DIGEST_NAME,
                      vb->vblock_i, z_file->num_txt_components_so_far, digest_display (piz_digest_so_far).s, digest_display (vb->digest_so_far).s, 
                      DIGEST_NAME,
                      txtfile_dump_vb (vb, z_name),
                      codec_args[txt_file->codec].viewer, file_guess_original_filename (txt_file),
                      vb->vb_position_txt_file + vb->txt_data.len, // head argument
                      (uint32_t)vb->txt_data.len,                  // tail argument
                      txtfile_dump_filename (vb, z_name, "good")); // redirection filename

                RESTORE_FLAG (quiet);

                failed = true; // no point in test the rest of the vblocks as they will all fail - MD5 is commulative
            }
        }
    }

    vb_digest_last++; // next please
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

// in v6-v11 we had a bug were the 2nd 32b of an MD5 digest was a copy of the first 32b
void digest_verify_ref_is_equal (const Reference ref, const char *header_ref_filename, const Digest header_md5)
{
    Digest ref_md5 = ref_get_file_md5 (ref);
    uint8_t version = ref_get_genozip_version (ref);

    ASSINP ((version >= 12 && digest_is_equal (ref_md5, header_md5)) ||
            (version <= 11 && ref_md5.words[0] == header_md5.words[0] && ref_md5.words[2] == header_md5.words[2] && ref_md5.words[3] == header_md5.words[3]),
            "%s: Bad reference file:\n%s (MD5=%s) was used for compressing\n%s (MD5=%s) has a different MD5",
            z_name, header_ref_filename, digest_display_ex (header_md5, DD_MD5).s, 
            ref_get_filename (ref), digest_display_ex (ref_md5, DD_MD5).s);
}

DigestDisplay digest_display_ex (const Digest digest, DigestDisplayMode mode)
{
    DigestDisplay dis = {};

    const uint8_t *b = digest.bytes; 
    
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

const char *digest_name (void) { return DIGEST_NAME; }
