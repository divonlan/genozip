// ------------------------------------------------------------------
//   mgzip_sections.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "mgzip.h"
#include "vblock.h"
#include "file.h"
#include "codec.h"
#include "bai.h"
#include "compressor.h"
#include "zfile.h"
#include "filename.h"

const FlagsMgzip bgzf_recompression_levels[1+MAX_FLAG_BGZF] = {
    { .library = BGZF_LIBDEFLATE19, .level = 0 },  // --bgzf=0 : BGZF blocks with no compression     
    { .library = BGZF_IGZIP,        .level = 1 },  // --bgzf=1 : note: this is IGZIP LVL0 
    { .library = BGZF_IGZIP,        .level = 2 },  // --bgzf=2 : note: this is IGZIP LVL1
    { .library = BGZF_LIBDEFLATE19, .level = 1 },  // --bgzf=3 
    { .library = BGZF_LIBDEFLATE19, .level = 7 },  // --bgzf=4 
    { .library = BGZF_LIBDEFLATE19, .level = 9 },  // --bgzf=5
};
#define BGZF_DEFAULT_LEVEL 2 // PIZ: used if --bgzf is not specified (it is actually faster than 1 if also writing to disk)

#define bgzf_no_recompression (FlagsMgzip){ .library = BGZF_NO_LIBRARY, .level = BGZF_NO_BGZF }

// PIZ: data read from SEC_GZ_* sections, to be moved to txt_file after it opens 
static struct {
    Buffer isizes, gz_digests; 
    uint64_t gz_file_size;
    Codec effective_codec;
} orig = {}; 

//-----------------------------------------------------
// Share ZIP/PIZ SIDE
//-----------------------------------------------------

// ZIP/PIZ output of --show-isizes
static void mgzip_show_isizes (CompIType comp_i)
{
    if (!flag.explicit_quiet) 
        iprint0 ("blk_i\tisize\tdigest\n");

    BufferP isizes_buf  = IS_ZIP ? &txt_file->mgzip_isizes  : &orig.isizes;
    BufferP digests_buf = IS_ZIP ? &txt_file->mgzip_digests : &orig.gz_digests;

    for_buf_tandem (uint32_t, isize_p, *isizes_buf, uint32_t, digest_p, *digests_buf)
        iprintf ("%i\t%u\t%08x\n", BNUM(*isizes_buf, isize_p), *isize_p, digests_buf->len ? *digest_p : 0);

    if (!flag.explicit_quiet) {
        if (IS_ZIP || VER2(15,80))
            iprintf ("comp_i=%u gz_file_size=%"PRIu64"effective_codec=%s\n", comp_i, 
                     IS_ZIP ? txt_file->disk_so_far : orig.gz_file_size, codec_name (orig.effective_codec));
        else
            iprintf ("comp_i=%u\n", comp_i);
    }

    if (is_genocat) exit_ok;
}

//-----------------------------------------------------
// ZIP SIDE
//-----------------------------------------------------

// ZIP: BGZF section per txt_file component
void mgzip_compress_SEC_GZ_sections (void)
{
    // cases where we don't write the BGZF blocks section
    if (!txt_file->is_exactable)  // --bgzf=exact will not work (not a MGZIP codec we can reconstruct, GZ library not idenfied, truncated, data-modifying command line options...)
        return;
    
    z_file->comp_num_EOF_blocks[flag.zip_comp_i] = txt_file->num_EOF_blocks;
    
    evb->comp_i = flag.zip_comp_i; // this goes into SectionEntFileFormat.comp_i via sections_add_to_list

    // sanity check
    uint64_t total_isize=0;
    for_buf (uint32_t, isize_p, txt_file->mgzip_isizes)
        total_isize += *isize_p;
    
    ASSERT (txt_file->mgzip_isizes.len == txt_file->mgzip_digests.len, 
            "expecting mgzip_isizes.len=%"PRIu64" == mgzip_digests.len=%"PRIu64,  
            txt_file->mgzip_isizes.len, txt_file->mgzip_digests.len);

    ASSERT (total_isize == txt_file->txt_data_so_far_single + txt_file->last_truncated_line_len, 
            "Expecting total_isize=%"PRId64" == txt_data_so_far_single=%"PRId64,
            total_isize, txt_file->txt_data_so_far_single);

    if (flag.show_isizes) mgzip_show_isizes (flag.zip_comp_i);

    BGEN_u32_buf (&txt_file->mgzip_isizes,  NULL);
    BGEN_u32_buf (&txt_file->mgzip_digests, NULL);

    txt_file->mgzip_isizes.len  *= sizeof (uint32_t); // len is now a count of bytes
    txt_file->mgzip_digests.len *= sizeof (uint32_t); 
    
    // compress SEC_GZ_ISIZES
    Codec codec = codec_assign_best_codec (evb, NULL, &txt_file->mgzip_isizes, SEC_GZ_ISIZES);

    SectionHeader isizes_header = (SectionHeader){ 
        .magic                 = BGEN32 (GENOZIP_MAGIC),
        .flags                 = (SectionFlags)txt_file->mgzip_flags,
        .section_type          = SEC_GZ_ISIZES,
        .data_uncompressed_len = BGEN32 (txt_file->mgzip_isizes.len32),
        .codec                 = codec,
    };

    comp_compress (evb, NULL, &evb->z_data, &isizes_header, txt_file->mgzip_isizes.data, NO_CALLBACK, "SEC_GZ_ISIZES");

    // compress SEC_GZ_DIGESTS
    SectionHeaderGzDigests digests_header = (SectionHeaderGzDigests){ 
        .magic                 = BGEN32 (GENOZIP_MAGIC),
        .section_type          = SEC_GZ_DIGESTS,
        .data_uncompressed_len = BGEN32 (txt_file->mgzip_digests.len32),
        .codec                 = CODEC_NONE, // digests are not compressible
        .effective_codec       = txt_file->effective_codec, // this is passed to piz, not related to compression of this section
        .gz_file_size          = BGEN64 (txt_file->disk_so_far), // works with stdin/url too
    };
    
    comp_compress (evb, NULL, &evb->z_data, &digests_header, txt_file->mgzip_digests.data, NO_CALLBACK, "SEC_GZ_DIGESTS");

    txt_file->mgzip_isizes.len  /= sizeof (uint32_t); // restore
    txt_file->mgzip_digests.len /= sizeof (uint32_t);
}

// ZIP: the is-exactable decision is made at four points, of them only #4 and #5 may report "true":
// 1. file_open_txt_read:    "false" if file is not .gz, .bam or CODEC_NONE (which is sometimes gz anyway)
// 2. file_open_txt_read_gz: "false" if CODEC_NONE file is confirmed as indeed not .gz, or for non-FASTQ/A, the gz file is confirmed as non-BGZF
// 3. bgzf_initialize_discovery : "false" in case fastq.gz, and not BGZF
// 4. bgzf_finalize_discovery: in case of a short file which ended before completing the BGZF test: may report true or false 
// 5. mgzip_uncompress_one_block: after discovering the BGZF compression libraries, we report whether the discovery was successful 
void mgzip_set_is_exactable (FileP file, bool is_exactable, rom reason_why_not)
{
    if (file->is_exactable == unknown || 
        (file->is_exactable == yes && is_exactable == no)) { // discovery said yes, but down the line we found out that no

        if (flag.show_is_exactable) {
            printf (flag.explicit_quiet ? "%s" : "Can 'genounzip --bgzf=exact' be used to re-compress this file with the same exact gz: %s", YN(is_exactable)); // stdout
            if (!is_exactable && !flag.explicit_quiet) printf (" (%s)", reason_why_not);
            printf ("\n");
            exit_ok;
        }
    
        file->is_exactable = is_exactable;
    }
}

//-----------------------------------------------------
// PIZ SIDE - setting up BGZF for a particular txt file
//-----------------------------------------------------

rom NON_EXACT_ERROR = "file cannot be gz-recompressed with the exact same gz-compression as the original file.";
static rom NON_EXACT_NOT_GZ = "file cannot be gz-recompressed with the exact same gz-compression because original file was not gz-compressed";

static FlagsMgzip bgzf_load_isizes (CompIType comp_i, bool is_CODEC_NONE, bool show_only) 
{
    Section sec = sections_get_comp_GZ_ISIZES_sec (comp_i);
    if (!sec) ignore: {
        ASSINP (flag.bgzf != BGZF_EXACT_STRICT, "%s %s", is_CODEC_NONE ? NON_EXACT_NOT_GZ : NON_EXACT_ERROR, WEBSITE_GZ);
        
        ASSERTW (show_only, "%s Using default gz-recompression.", is_CODEC_NONE ? NON_EXACT_NOT_GZ : NON_EXACT_ERROR);

        flag.bgzf = BGZF_EXACT_FAILED; // so bai_initialize does write an index file after all
        goto fallback; // this component doesn't contain a BGZF section
    }

    int32_t offset = zfile_read_section (z_file, evb, 0, &evb->z_data, "z_data", SEC_GZ_ISIZES, sec);
    SectionHeaderP header = (SectionHeaderP)Bc(evb->z_data, offset);

    // if we don't know the compression level (in older Genozip versions we wrote the SEC_GZ_ISIZES even 
    // if level discovery failed)
    if (header->flags.mgzip.level == BGZF_COMP_LEVEL_UNKNOWN) 
        goto ignore;

    zfile_uncompress_section (evb, header, &orig.isizes, "txt_file->mgzip_isizes", 0, SEC_GZ_ISIZES);

    if (show_only) goto fallback;
    
    // up to 15.0.63 section contains an array of 32bit isizes to eventually support MGZIP compression schemes beyond BGZF (still, only BGZF supported)
    if (VER2(15,63)) { 
        if (!VER2(15,80) && (Z_DT(FASTQ) || Z_DT(FASTA) || (flag.deep && comp_i >= SAM_COMP_FQ00)) && 
            header->flags.mgzip.level == 0 && header->flags.mgzip.level == 0) goto ignore; // defect 2026-03-04

        orig.isizes.len /= sizeof (uint32_t);
        BGEN_u32_buf (&orig.isizes, NULL);
    }

    // up to 15.0.62 buffer was 16 bit, values were (isize-1), and EOF was indicated by header.has_eof_block
    else {
        orig.isizes.len /= sizeof (uint16_t);
        buf_alloc (evb, &orig.isizes, 0, orig.isizes.len + 1, uint32_t, 0, NULL);

        for_buf_tandem_back (uint16_t, isize16, orig.isizes, uint32_t, isize32, orig.isizes)
            *isize32 = (uint32_t)BGEN16 (*isize16) + 1;

        if (header->flags.mgzip.OLD_has_eof_block)
            BNXT32(orig.isizes) = 0; // append EOF block
    }

    return header->flags.mgzip; // mgzip_isizes successfully loaded

fallback:
    buf_destroy (orig.isizes);
    return bgzf_recompression_levels[BGZF_DEFAULT_LEVEL];
}          

static void bgzf_load_digests (CompIType comp_i, bool show_only) 
{
    if (!orig.isizes.len && !show_only) return; // since no isizes section, there will be no digests either
    if (!VER2(15,80)) goto done; // SECTION_GZ_DIGESTS only available in exactable files since 15.0.80

    Section sec = sections_get_comp_GZ_ISIZES_sec (comp_i) + 1; // SECTION_GZ_DIGESTS always follows SECTION_GZ_ISIZES

    int32_t offset = zfile_read_section (z_file, evb, 0, &evb->z_data, "z_data", SEC_GZ_DIGESTS, sec);
    SectionHeaderGzDigestsP header = (SectionHeaderGzDigestsP)Bc(evb->z_data, offset);

    zfile_uncompress_section (evb, header, &orig.gz_digests, "txt_file->mgzip_digests", 0, SEC_GZ_DIGESTS);

    if (show_only) {
        buf_destroy (orig.gz_digests);
        return;
    }

    // these will all be grabbed/copied to txt_file after it opens
    orig.gz_file_size    = BGEN64(header->gz_file_size);  
    orig.effective_codec = header->effective_codec;
    orig.gz_digests.len /= sizeof (uint32_t);
    BGEN_u32_buf (&orig.gz_digests, NULL);
    
done:
    if (flag.show_isizes && !show_only) mgzip_show_isizes (comp_i); // show isizes and digests
}          

// PIZ: called from main thread after reading txt_header's header
FlagsMgzip mgzip_piz_calculate_mgzip_flags (CompIType comp_i, Codec src_codec)
{
    #define C(cdc) (src_codec == CODEC_##cdc)
    FlagsMgzip mgzip_flags;

    #define HAS_EXT(x) filename_has_ext (flag.out_filename, #x)
    bool bgzf_implied_by_out_filename = flag.out_filename && (HAS_EXT(.gz) || HAS_EXT(.bgz) || HAS_EXT(.bam));
    bool no_bgzf_implied_by_out_filename = file_piz_get_dt_of_out_filename() == flag.out_dt && !(HAS_EXT(.gz) || HAS_EXT(.bgz) || HAS_EXT(.bam));
    bool mgzip_section_was_read = false;

    // cases where there is no BGZF re-compression 
    if (flag.test    || 
        OUT_DT(CRAM) ||
        (flag.bgzf == 0 && !OUT_DT(BAM) && !OUT_DT(BCF)) || // note: in BCF and BAM --bgzf=0 means BGZF blocks with no compression (as opposed to no BGZF at all)
        (flag.bgzf == BGZF_EXACT && C(NONE) && flag.reconstruct_as_src))  // case: --bgzf=exact and source codec was CODEC_NONE
        
        mgzip_flags = bgzf_no_recompression; 
    
    // case: reconstructing BCF: piz sends VCF to bcftools in CODEC_NONE, and bcftools compressed by the level given in mgzip_flags
    else if (OUT_DT(BCF))
        mgzip_flags = (FlagsMgzip){ .library = BGZF_EXTERNAL_LIB, 
                                    .level = (flag.bgzf < 0) ? 4 : (int[]){0, 2, 4, 6, 8, 9 }[flag.bgzf] }; // convert Genozip level 0-5 to bcftools level 0-9
    
    // case: --bgzf=exact[-strict]
    else if (IS_PIZ_EXACT || flag.show_isizes) {
        mgzip_flags = bgzf_load_isizes (comp_i, C(NONE), false);
        bgzf_load_digests (comp_i, false); 
        mgzip_section_was_read = true; // section read and header displayed if requested, regardless if successfully loaded or not
    }
    
    // case: --bgzf=0 to 5
    else if (flag.bgzf >= 0) {
        mgzip_flags = bgzf_recompression_levels[flag.bgzf]; // set to --bgzf command line value

        // if user specified --bgzf and --output - make sure output filename is .gz, .bam or .bcf
        ASSINP (flag.force || !flag.out_filename || bgzf_implied_by_out_filename || HAS_EXT(.bcf) || mgzip_flags.level==0, 
                "using %s in combination with %s for outputting a %s file, requires the output filename to end with %s (override with --force)", 
                OT("output", "o"), OT("bgzf", "z"), dt_name(flag.out_dt), OUT_DT(BAM)?".bam" : OUT_DT(BCF)?".bcf" : ".gz");
    
        ASSINP (!OUT_DT(BCF) || !IS_PIZ_EXACT, // because we have no control over bcftools' BGZF block generation
                "cannot use --bgzf=exact%s when outputing a BCF file", flag.bgzf == BGZF_EXACT_STRICT ? "-strict" : ""); 
    }

    // case: genocat to stdout without --bgzf: - no re-compression. 
    else if (is_genocat && !flag.out_filename)
        mgzip_flags = OUT_DT(BAM) ? bgzf_recompression_levels[0] : bgzf_no_recompression; // file_open_txt_write interprets level=0 as CODEC_BGZF without compression for BAM, and CODEC_NONE for other types 

    // case: genocat or genounzip out_filename and no --bgzf: - determine by file name (except BAM - bgzf regardless of filename)
    else if (flag.out_filename)
        mgzip_flags = (bgzf_implied_by_out_filename || OUT_DT(BAM) || (!no_bgzf_implied_by_out_filename && !C(NONE))) ? bgzf_recompression_levels[BGZF_DEFAULT_LEVEL] : bgzf_no_recompression;
    
    // case: genounzip without explicit filename, and no --bgzf: default compression or no compression
    else
        // note: for bz2, xz, and zip - we reconstruct as gz too. better choice than plain.
        mgzip_flags = (IS_GZIP(src_codec) || C(BAM) || C(BZ2) || C(XZ) || C(ZIP)) ? bgzf_recompression_levels[BGZF_DEFAULT_LEVEL] : bgzf_no_recompression; // note: similar logic to txtheader_piz_get_filename
    
    // case: user wants to see this section header, despite not needing BGZF data
    if (!mgzip_section_was_read && (flag.only_headers == SEC_GZ_ISIZES+1 || (flag.only_headers && flag.show_sec_headers[SEC_GZ_ISIZES])))
        bgzf_load_isizes (comp_i, C(NONE), true); 

    if (!mgzip_section_was_read && (flag.only_headers == SEC_GZ_DIGESTS+1 || (flag.only_headers && flag.show_sec_headers[SEC_GZ_DIGESTS])))
        bgzf_load_digests (comp_i, true); 
        
    if (flag_show_bgzf)
        iprintf ("comp_i=%u with src_codec=%s out_dt=%s: calculated mgzip_flags={%s, %d}\n",
                 comp_i, codec_name (src_codec), dt_name (flag.out_dt), bgzf_library_name (mgzip_flags.library, true), mgzip_flags.level);

    return mgzip_flags;
    #undef C
}

// PIZ main thread: update txt_file with BGZF info calculated earlier
void mgzip_piz_set_txt_file_info (FlagsMgzip mgzip_flags, uint32_t OLD_gz_size_3LSB)
{
    if (orig.isizes.len) {
        buf_grab (evb, txt_file->mgzip_isizes, "txt_file->mgzip_isizes", orig.isizes);

        if (VER2(15,80)) {
            buf_grab (evb, txt_file->mgzip_digests, "txt_file->mgzip_digests", orig.gz_digests);
            txt_file->effective_codec   = orig.effective_codec;
            txt_file->orig_gz_file_size = orig.gz_file_size;
        }

        else {
            txt_file->effective_codec   = CODEC_BGZF; 
            txt_file->orig_gz_file_size = OLD_gz_size_3LSB; // for files up to 15.0.79: only 3 LSB of file size
        }

        // reinitialize for next file
        buf_destroy (orig.isizes);
        buf_destroy (orig.gz_digests);
        memset (&orig, 0, sizeof (orig)); 
    }

    txt_file->mgzip_flags  = mgzip_flags;

    // sanity
    ASSERT (txt_file->mgzip_flags.level >= 0 && txt_file->mgzip_flags.level <= BGZF_MAX_LEVEL, "txt_file->mgzip_flags.level=%u ∉ [0,%u]", 
            txt_file->mgzip_flags.level, BGZF_MAX_LEVEL);

    ASSERT (txt_file->mgzip_flags.library >= 0 && txt_file->mgzip_flags.library < NUM_BGZF_LIBRARIES, "txt_file->mgzip_flags.library=%u ∉ [0,%u]", 
            txt_file->mgzip_flags.level, NUM_BGZF_LIBRARIES-1);
}

rom bgzf_library_name (MgzipLibraryType library, bool long_name)
{
    return (library < 0 || library >= NUM_ALL_BGZF_LIBRARIES) ? "INVALID_BGZF_LIBRARY"
         : long_name ? (rom[])BGZF_LIB_NAMES_LONG[library]
         :             (rom[])BGZF_LIB_NAMES_SHRT[library];
}

StrText bgzf_lib_name_level (FlagsMgzip mgzip_flags)
{
    StrText s;
    snprintf (s.s, sizeof(s), "%s⁀%u", bgzf_library_name (mgzip_flags.library, false), mgzip_flags.level);

    return s;
}
