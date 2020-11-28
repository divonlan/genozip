// ------------------------------------------------------------------
//   digest.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <libdeflate.h>
#include "genozip.h"
#include "digest.h"
#include "buffer.h"
#include "md5.h"
#include "vblock.h"
#include "mutex.h"
#include "file.h"
#include "codec.h"
#include "txtfile.h"

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
        fprintf (stderr, "%s finalize %s: %s\n", DIGEST_NAME, msg, digest_display (digest).s);

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
    if (IS_ADLER) {
        if (!ctx->initialized) {
            ctx->adler_ctx.adler = 1;
            ctx->initialized = true;
        }
        ctx->adler_ctx.adler = libdeflate_adler32 (ctx->adler_ctx.adler, buf->data, buf->len);
    }

    else {
        if (!ctx->initialized) {
            md5_initialize (&ctx->md5_ctx);
            ctx->initialized = true;
        }
        md5_update (&ctx->md5_ctx, buf->data, buf->len);
    }

    if (flag.show_digest) 
        fprintf (stderr, "%s update %s (len=%"PRIu64"): %s\n", DIGEST_NAME, msg, buf->len, digest_display (digest_snapshot (ctx)).s);
}

// ZIP: called by compute thread to calculate MD5 for one VB - need to serialize VBs using a mutex
void digest_one_vb (VBlock *vb)
{
    // wait for our turn
    while (1) {
        mutex_lock (vb_digest_mutex);
        if (vb_digest_last == vb->vblock_i - 1) break; // its our turn now

        // not our turn, wait 10ms and try again
        mutex_unlock (vb_digest_mutex);
        usleep (10000);
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

        digest_update (&txt_file->digest_ctx_bound, &vb->txt_data, "vb:digest_ctx_bound");

        // if testing, compare digest up to this VB to that calculated on the original file and transferred through SectionHeaderVbHeader
        // note: we cannot test this unbind mode, because the digests are commulative since the beginning of the bound file
        if (!failed && !flag.unbind && !v8_digest_is_zero (vb->digest_so_far)) {
            Digest piz_digest_so_far = digest_snapshot (&txt_file->digest_ctx_bound);

            // warn if VB is bad, but don't exit, so file reconstruction is complete and we can debug it
            if (!digest_is_equal (vb->digest_so_far, piz_digest_so_far)) {
uint8_t *b=vb->digest_so_far.bytes;
printf ("XXXX  vb->digest_so_far %2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x\n", 
b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12], b[13], b[14], b[15]);
b=piz_digest_so_far.bytes;
printf ("XXXX  vb->digest_so_far %2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x\n", 
b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12], b[13], b[14], b[15]);

                // dump bad vb to disk
                WARN ("%s of reconstructed vblock=%u (%s) differs from original file (%s).\n"
                      "Bad reconstructed vblock has been dumped to: %s\n"
                      "To see the same data in the original file:\n"
                      "   %s %s | head -c %"PRIu64" | tail -c %u > %s", DIGEST_NAME,
                      vb->vblock_i, digest_display (piz_digest_so_far).s, digest_display (vb->digest_so_far).s, txtfile_dump_vb (vb, z_name),
                      codec_args[txt_file->codec].viewer, file_guess_original_filename (txt_file),
                      vb->vb_position_txt_file + vb->txt_data.len, (uint32_t)vb->txt_data.len, txtfile_dump_filename (vb, z_name, "good"));

                failed = true; // no point in test the rest of the vblocks as they will all fail - MD5 is commulative
            }
        }
    }

    vb_digest_last++; // next please
    mutex_unlock (vb_digest_mutex);
}

DigestDisplay digest_display_ex (const Digest digest, DigestDisplayMode mode)
{
    DigestDisplay dis;

    const uint8_t *b = digest.bytes; 
    
    if ((mode == DD_NORMAL && !IS_ADLER && md5_is_zero (digest)) ||
        (mode == DD_MD5                 && md5_is_zero (digest)) || 
        (mode == DD_MD5_IF_MD5 && (IS_ADLER || md5_is_zero (digest)))) 
        sprintf (dis.s, "N/A                             ");

    if (mode == DD_NORMAL && IS_ADLER)
        sprintf (dis.s, "%-10u", BGEN32 (digest.adler_bgen));

    if ((mode == DD_NORMAL && !IS_ADLER && !md5_is_zero (digest)) ||
        (mode == DD_MD5                 && !md5_is_zero (digest)) || 
        (mode == DD_MD5_IF_MD5 && !IS_ADLER && !md5_is_zero (digest)))
        sprintf (dis.s, "%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x", 
                 b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12], b[13], b[14], b[15]);    
    
    return dis;
}
DigestDisplay digest_display (const Digest digest) { return digest_display_ex (digest, DD_NORMAL); }

const char *digest_name (void) { return DIGEST_NAME; }
