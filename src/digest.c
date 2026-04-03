// ------------------------------------------------------------------
//   digest.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "libdeflate_1.19/libdeflate.h"
#include "xxhash/xxhash.h"
#include "digest.h"
#include "md5.h"
#include "filename.h"
#include "codec.h"
#include "txtfile.h"
#include "gencomp.h"
#include "progress.h"
#include "writer.h"
#include "txtheader.h"
#include "piz.h"
#include "reference.h"

#define IS_MD5 (IS_ZIP ? flag.md5 : !z_file->z_flags.not_md5)

#define DIGEST_NAME digest_alg_name (IS_MD5?DIGEST_MD5 : (IS_ZIP || VER2(15,81))?DIGEST_XXH3 : DIGEST_ADLER)

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

// get digest of data so far in a commulative digest: MD5, and up to v13 - also Adler32 (never XXH3)
Digest digest_snapshot (const DigestState *state, rom msg)
{
    // make a copy of state, and finalize it, keeping the original copy unfinalized
    DigestState ctx_copy = *state;

    Digest digest = (state->alg == DIGEST_ADLER) ? (Digest){ .words = { BGEN32 (ctx_copy.adler_state.adler) } }
                                                    : md5_finalize (&ctx_copy.md5_state);

    if (msg && flag.show_digest) 
        iprintf ("%s finalize %s: %s\n", DIGEST_NAME, msg, digest_display_ex (digest, DD_NORMAL).s);

    return digest;
}

Digest digest_do (rom data, uint64_t data_len, DigestAlg alg, rom show_digest_msg)
{
    Digest digest;
    switch (alg) {
        case DIGEST_ADLER : digest = (Digest){ .adler32 = BGEN32 (adler32 (1, STRa(data))) };  break;
        case DIGEST_XXH3  : digest = (Digest){ .xxh3   = BGEN64 (XXH3_64bits (STRa(data))) }; break;
        case DIGEST_MD5   : digest = md5_do (STRa(data)); break;
        default           : ABORT ("Invalid digest alg=%d", alg);
    }

    if (flag.show_digest)
        iprintf ("%s digest_do(%"PRIu64" bytes) %s: %s\n", DIGEST_NAME, data_len, show_digest_msg,  digest_display_ex_(digest, DD_NORMAL, alg, false).s);

    return digest;
}

#define digest_update(state, buf, msg) digest_update_do ((buf)->vb, (state), STRb(*(buf)), (msg))
static void digest_update_do (VBlockP vb, DigestState *state, rom data, uint64_t data_len, rom msg)
{
    if (!data_len) return;
    
    START_TIMER;
    DigestState before;

    if (!IS_MD5) { // used for files up to v13 where adler digest was commulative
        if (!state->initialized) {
            state->alg = DIGEST_ADLER;
            state->adler_state.adler = 1;
            state->adler_state.initialized = true;
        }
        if (flag.show_digest) before = *state;
        state->adler_state.adler = adler32 (state->adler_state.adler, STRa(data));
    }

    else {
        if (!state->initialized) 
            md5_initialize (&state->md5_state, state->md5_state.log); // preserve log

        if (flag.show_digest) before = *state;
        md5_update (&state->md5_state, STRa(data));
    }

    state->bytes_digested += data_len;

    if (flag.show_digest) {
        iprintf ("vb=%10s %s update %s (len=%"PRIu64" so_far=%"PRIu64") 32chars=\"%s\": before=%s after=%s\n", 
                 VB_NAME, DIGEST_NAME, msg, data_len, state->bytes_digested, 
                 str_to_printable_(data, MIN_(32, data_len)).s, 
                 digest_display_ex (digest_snapshot (&before, NULL),  DD_NORMAL).s, 
                 digest_display_ex (digest_snapshot (state, NULL), DD_NORMAL).s);
    }

    if (state->log) {
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
    if (VER(14) && z_file->z_flags.not_md5) {
        
        CompIType comp_i = ((flag.deep || flag.pair) && flag.one_component) ? flag.one_component-1
                         : (flag.deep && txt_file_i >= 1) ? (txt_file_i - 1 + SAM_COMP_FQ00)
                         : (flag.pair && txt_file_i ==1 ) ? FQ_COMP_R2
                         :                                  COMP_MAIN;

        uint32_t expected_vbs_verified = sections_get_num_vbs (comp_i);

        ASSERT (z_file->num_vbs_verified == expected_vbs_verified ||  // success
                txt_file->vb_digest_failed,                           // failure already announced
                "Expected to have verified (not_md5) all %u VBlocks, but verified %u (txt_file=%s txt_file_i=%u)",
                expected_vbs_verified, z_file->num_vbs_verified, txt_name, txt_file_i);

        if (flag.show_digest)
            iprintf ("Txt file #%u: %u VBs verified\n", txt_file_i, z_file->num_vbs_verified);

        if (flag.test) { 
            snprintf (s, sizeof (s), "verified as identical to the %s %s", 
                      segconf.zip_txt_modified ? "modified" : "original", // in case ZIP modified, e.g. with --optimize
                      dt_name_faf (txt_file->data_type));
            progress_finalize_component (s); 
        }

        z_file->num_vbs_verified = 0; // reset for next component
    }

    // txt file commulative digest: MD5, and up to v13, also Adler32
    else {
        Digest decompressed_file_digest = digest_snapshot (&z_file->digest_state, "file"); 

        if (digest_recon_is_equal (decompressed_file_digest, z_file->digest)) {
            if (flag.test || IS_MD5) { 
                snprintf (s, sizeof (s), "verified as identical to the original %s (%s=%s)", 
                          dt_name_faf (txt_file->data_type), digest_name(), digest_display (decompressed_file_digest).s);
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
                  dt_name_faf (txt_file->data_type), digest_display (z_file->digest).s, txt_name, txt_file_i);
        }
    }
}

static void digest_piz_verify_one_vb (VBlockP vb, 
                                      BufferP txt_data) // may differ from vb->txt_data if integrated PRIM/DEPN lines
{
    // Compare digest up to this VB transmitted through SectionHeaderVbHeader. If Adler32, it is a stand-alone
    // digest of the VB, and if MD5, it is a commulative digest up to this VB.
    if ((!txt_file->vb_digest_failed || !IS_MD5) && // note: for MD5, we report only the first failed VB, bc the digest is commulative, so all subsequent VBs will fail for sure
        (VER(14) || !flag.unbind)) {                // note: for files <= v13, we cannot test per-VB digest in unbind mode, because the digests (MD5 and Adler32) are commulative since the beginning of the bound file. However, we still test component-wide digest in piz_verify_digest_one_txt_file.

        // add VB to commulative digests as needed
        Digest single_comp_commulative_digest = (!VER(14) || IS_MD5)      ? digest_snapshot (&z_file->digest_state, NULL) : DIGEST_NONE; // up to v13 Adler was commulative too
        Digest multi_comp_commulative_digest  = (!VER(14) && flag.unbind) ? digest_snapshot (&z_file->v13_commulative_digest_state, NULL) : DIGEST_NONE;

        Digest piz_digest = (VER(14) && !IS_MD5)      ? vb->digest  // stand-alone digest of this VB
                          : (!VER(14) && flag.unbind) ? multi_comp_commulative_digest
                          :                             single_comp_commulative_digest;

        // warn if VB is bad, but don't exit, so file reconstruction is complete and we can debug it
        if (!digest_recon_is_equal (piz_digest, vb->expected_digest)) { 

            // compare recon_size (note: this excludes integrated PRIM/DEPN lines)
            char recon_size_warn[128] = "";
            if (vb->recon_size != vb->txt_data.len) // note: leave vb->txt_data.len 64bit to detect bugs
                snprintf (recon_size_warn, sizeof (recon_size_warn), "Expecting: VB_HEADER.recon_size=%u == txt_data.len=%"PRIu64"%s\n", 
                          vb->recon_size, vb->txt_data.len, (z_sam_gencomp ? " (note these sizes exclude gencomp data)" : ""));  

            NOISYWARN ("reconstructed vblock=%s/%u (vb_line_i=0 -> txt_line_i(1-based)=%"PRId64" num_lines=%u), (%s=%s) differs from the %s file (%s=%s). File compressed with version=%s uncompressed with=%s\n%s",
                       comp_name (vb->comp_i), vb->vblock_i, writer_get_txt_line_i (vb, 0), vb->lines.len32,
                       DIGEST_NAME, digest_display (piz_digest).s, 
                       segconf.zip_txt_modified ? "modified" : "original",
                       DIGEST_NAME, digest_display (vb->expected_digest).s, 
                       STRver(file_version()).s , STRver(code_version()).s, 
                       recon_size_warn);

            // case: first bad VB: dump bad VB to disk and maybe exit
            if (!test_and_set_relaxed (txt_file->vb_digest_failed)) { // not WARN_ONCE because we might be genounzipping multiple files - we want to show this for every failed file (see also note in digest_piz_verify_one_txt_file)
                NOISYWARN ("Bad reconstructed %s has been dumped to: %s. vb_position_txt_file=%"PRIu64", txt_data.len=%u recon_size(expected)=%u\n%s%s",
                           VB_NAME, txtfile_dump_vb (vb, z_name, txt_data).s, // note: txt_data dumped INCLUDES integrated gencomp lines (if there are any)
                           vb->vb_position_txt_file, vb->txt_data.len32, vb->recon_size, 
                           piz_advise_biopsy (vb).s, report_support_if_unexpected());

                if (flag.test) exit_on_error (false); // must be inside the atomic test, otherwise another thread will exit before we completed dumping
            }
        }

        else {
            if (flag.show_digest)
                iprintf ("%s: verfied digest\n", VB_NAME);

            increment_relaxed (z_file->num_vbs_verified); // count VBs verified so far
        }
    }
}


// ZIP and PIZ: called by compute thread to calculate MD5 or Adler32 of one VB - possibly serializing VBs using a mutex
bool digest_one_vb (VBlockP vb, bool is_compute_thread, 
                    BufferP txt_data) // if NULL, txt_data digested is vb->txt_data (if not NULL: this might be PIZ of a SAM MAIN vb with integrated PRIM and DEPN lines)
{
    #define WAIT_TIME_USEC 1000
    #define DIGEST_TIMEOUT (30*60) // 30 min

    if (!txt_data) txt_data = &vb->txt_data;

    bool digestable = gencomp_comp_eligible_for_digest(vb);

    // starting V14, unless md5, we digest each VB stand-alone.
    if (!IS_MD5 && (IS_ZIP || VER(14))) {
        if (digestable) {
            vb->digest = digest_do (STRb(*txt_data), (IS_ZIP || VER2(15,81)) ? DIGEST_XXH3 : DIGEST_ADLER, VB_NAME);

            if (IS_PIZ) digest_piz_verify_one_vb (vb, txt_data);  
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
            digest_update (&z_file->digest_state, txt_data, "vb");

            // take a snapshot of the commulative digest as per the end of this VB 
            vb->digest = digest_snapshot (&z_file->digest_state, NULL);
        }
            
        else if (digestable && IS_PIZ) {
            digest_update (&z_file->digest_state, txt_data, "vb"); // labels consistent with ZIP so we can easily diff PIZ vs ZIP

            if (!VER(14) && flag.unbind)
                digest_update (&z_file->v13_commulative_digest_state, txt_data, "vb"); // up to v13, multi-component, VB digest is commulative across all components (in MD5 and Adler)

            digest_piz_verify_one_vb (vb, txt_data);        
        }

        if (is_compute_thread)
            serializer_unlock (z_file->digest_serializer);
    }

    return digestable;
}

// ZIP and PIZ: called by main thread to calculate MD5 or Adler32 of the txt header
Digest digest_txt_header (BufferP txt_data, Digest piz_expected_digest, CompIType comp_i)
{
    START_TIMER;

    // start digest log if need
    if (flag.log_digest) {
        z_file->digest_state.log = true;
        file_remove (DIGEST_LOG_FILENAME, true);    
    }

    // serialize VBs of this txt file (if MD5, or if v13 or earlier)
    if (IS_PIZ && (IS_MD5 || !VER(14)) && sections_get_num_vbs (comp_i)) 
        z_file->digest_serializer.vb_i_last = sections_get_first_vb_i (comp_i) - 1; 

    if (!txt_data->len) return DIGEST_NONE;

    Digest digest;
    
    // in 15.0.81 we switched from Adler32 to XXH3
    if (!IS_MD5 && VER2(15,81))
        digest = digest_do (STRb(*txt_data), DIGEST_XXH3, "TXT_HEADER");

    // starting V14, if adler32, we digest each TXT_HEADER stand-alone.
    else if (!IS_MD5 && VER(14))
        digest = digest_do (STRb(*txt_data), DIGEST_ADLER, "TXT_HEADER");

    // if MD5 or v13 (or earlier) - digest is for commulative for the whole file, and we take a snapshot
    else {
        digest_update (&z_file->digest_state, txt_data, "txt_header");
        digest = digest_snapshot (&z_file->digest_state, NULL);

        // up to v13, in unbind, VB digest relied on commulative digest of all componets
        if (!VER(14) && flag.unbind)
            digest_snapshot (&z_file->v13_commulative_digest_state, NULL);
    }

    if (IS_PIZ) {  
        TEMP_FLAG (quiet, flag.quiet && !flag.show_digest);

        if (!digest_recon_is_equal (digest, piz_expected_digest)) {
            
            WARN ("%s of reconstructed %s header (%s) differs from %s file (%s)\n"
                  "Bad reconstructed header has been dumped to: %s\n"
                  "To see the same data in the original file:\n"
                  "genozip --biopsy 0 %s%s",
                  digest_name(), dt_name (z_file->data_type), digest_display (digest).s, 
                  segconf.zip_txt_modified ? "modified" : "original", // in case ZIP modified, e.g. with --optimize
                  digest_display (piz_expected_digest).s,
                  txtfile_dump_vb (txt_data->vb, z_name, txt_data).s, filename_guess_original (txt_file), report_support_if_unexpected());

            if (flag.test) exit_on_error(false);
        }

        RESTORE_FLAG (quiet);
    }

    COPY_TIMER_EVB (digest_txt_header);

    return digest;
}

// test for matching digest between loaded external reference and reference specified in SectionHeaderGenozipHeader
void digest_verify_correct_external_reference_is_loaded (void)
{
    // case: we loaded an external ref, even though the file was not comprssed with one. E.g. when translating 23andMe to VCF.    
    if (digest_is_zero (z_file->ref_genome_digest)) return;

    Digest loaded_ref_digest = ref_get_genome_digest();
    Version loaded_ref_ver   = ref_get_genozip_ver();
    DigestAlg loaded_ref_alg = ref_get_genome_digest_alg();

    // in v6-v11 we had a bug where the 2nd 32b of an MD5 digest was a copy of the first 32b.
    // we don't know the genozip version of the reference file used for compressing - therefore we accept both MD5 options.
    // this way, if a reference is re-generated from FASTA using a new genozip, it can be used to decompress an old file
    
    // loaded reference is v14 or older: alg is always MD5, but digest is of the FASTA, no the in-memory genome
    if (loaded_ref_ver.major < 15) 
        ASSINP ((loaded_ref_digest.words[0] == z_file->ref_genome_digest.words[0] && 
                 loaded_ref_digest.words[2] == z_file->ref_genome_digest.words[2] && 
                 loaded_ref_digest.words[3] == z_file->ref_genome_digest.words[3]),
                 "%s: Incorrect reference file:\n%s (%s) was used for compressing\nNot %s (%s version=%u)",
                 z_name, z_file->ref_filename_used_in_zip, digest_display_ex_(z_file->ref_genome_digest, DD_NORMAL, DIGEST_UNKNOWN, true).s, 
                 ref_get_filename(), digest_display_ex_(loaded_ref_digest, DD_NORMAL, DIGEST_MD5, true).s, loaded_ref_ver.major);

    // reference files produced since v15 can be either MD5 or (XXH3 since 15.0.81 and Adler32 prior), and digest is of in-memory genome Buffer
    else {
        if (VER(15)) { // file and loaded reference are both v15
            // case: file was compressed with a different external reference digested with a different alg
            // but it might be still ok if the genome is the same. re-digest the genome with the the alg recorded in the compressed file
            if (z_file->ref_genome_digest_alg != loaded_ref_alg) 
                loaded_ref_digest = reference_re_digest_genome (z_file->ref_genome_digest_alg);

            ASSINP (digest_is_equal (z_file->ref_genome_digest, loaded_ref_digest), 
                    "%s: Incorrect reference file:\n%s with %s was used for compressing (ref_genozip_version=%s)\nBut %s has %s (ref_genozip_version=%s)",
                    z_name, z_file->ref_filename_used_in_zip, digest_display_ex_(z_file->ref_genome_digest, DD_NORMAL, z_file->ref_genome_digest_alg, true).s, 
                    STRver(z_file->ref_genozip_ver).s, // available since 15.0.81 
                    ref_get_filename(), digest_display_ex_(loaded_ref_digest, DD_NORMAL, z_file->ref_genome_digest_alg, true).s, STRver(loaded_ref_ver).s);
        }

        else { // loaded reference is v15, file is older
            // cannot uncompress file compressed with Genozip <= v14 with a reference generated with >= v15, because
            // in the former uses a reference file with a the digest is of the FASTA, and in references >= v15 the digest is of the in-memory genome
            ABORTINP ("Error: The reference file %s was generated with Genozip version %u, but the file %s was compressed with a reference file generated by an earlier version of Genozip. To uncompress this file, use a reference file generated with Genozip <= 14.%s",
                      ref_get_filename(), loaded_ref_ver.major, z_file->basename, 
                      cond_str (!VER(14), " ", "Note: if the reference file used to compress was generated with Genozip version 13.0.20 or earlier, then the original reference file must be used and it is not possible to re-generate it.")); // defect 2022-08-22
        }
    }
}

DigestDisplay digest_display_ex_(const Digest digest, DigestDisplayMode mode, DigestAlg alg, bool show_alg)
{
    DigestDisplay dis = {};

    bytes b = digest.md5; 
    
    if (alg == DIGEST_UNKNOWN) 
        alg = digest_guess_alg (digest);

    rom alg_name = show_alg ? digest_alg_name(alg) : "";
    rom alg_eq   = show_alg ? "="                  : "";

    bool md5_output = (mode == DD_MD5) || ((mode == DD_NORMAL || mode == DD_MD5_IF_MD5) && alg == DIGEST_MD5);

    if (md5_output && digest_is_zero (digest))
        snprintf (dis.s, sizeof (dis.s), "N/A                                    ");

    else if (mode == DD_MD5_IF_MD5 && alg != DIGEST_MD5)
        snprintf (dis.s, sizeof (dis.s), "N/A                                    ");

    else if ((mode == DD_NORMAL || mode == DD_SHORT) && alg == DIGEST_ADLER)
        snprintf (dis.s, sizeof (dis.s), "%s%s%08x", alg_name, alg_eq, BGEN32 (digest.adler32));

    else if ((mode == DD_NORMAL || mode == DD_SHORT) && alg == DIGEST_XXH3)
        snprintf (dis.s, sizeof (dis.s), "%s%s%016"PRIx64, alg_name, alg_eq, BGEN64 (digest.xxh3));

    else if (md5_output && !digest_is_zero (digest))
        snprintf (dis.s, sizeof (dis.s), "%s%s%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x", 
                  alg_name, alg_eq,
                  b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12], b[13], b[14], b[15]);    
    
    else if (mode == DD_SHORT && alg==DIGEST_MD5)
        snprintf (dis.s, sizeof (dis.s), "%s%s%2.2x%2.2x%2.2x%2.2x", alg_name, alg_eq, b[0], b[1], b[2], b[3]);    

    return dis;
}

DigestDisplay digest_display_ex (const Digest digest, DigestDisplayMode mode)
{
    return digest_display_ex_(digest, mode, IS_MD5?DIGEST_MD5 : VER2(15,81)?DIGEST_XXH3 : DIGEST_ADLER, false);
}

DigestDisplay digest_display (const Digest digest) 
{ 
    return digest_display_ex (digest, DD_NORMAL); 
}

DigestDisplay digest_display_(const Digest digest, DigestAlg alg) 
{ 
    return digest_display_ex_(digest, DD_NORMAL, alg, false); 
}

rom digest_name (void) { return DIGEST_NAME; }

rom digest_alg_name (DigestAlg alg) 
{ 
    switch (alg) {
        case DIGEST_ADLER   : return "Adler32";
        case DIGEST_MD5     : return "MD5";
        case DIGEST_XXH3    : return "XXH3";
        case DIGEST_UNKNOWN : return "Unknown";
        default             : return "Invalid";
    }
}

