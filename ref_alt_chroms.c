// ------------------------------------------------------------------
//   ref_alt_chroms.c
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "ref_private.h"
#include "sections.h"
#include "data_types.h"
#include "file.h"
#include "zfile.h"
#include "endianness.h"
#include "strings.h"
#include "vblock.h"
#include "mutex.h"
#include "seg.h"
#include "txtheader.h"

// ZIP of a file with an external reference: 
// When a chrom index in the txt file matches an a different chrom index in the reference file, 
// we create a mapping here pass it to Piz as a SEC_REF_ALT_CHROMS section. It contains at most as many entries 
// as the number of contigs in the reference file.

// this is needed in two cases:
// 1. In a file without a header: in case the chrom_name in the txt is an alternate name to the one in the header (eg "22"->"chr22")
// 2. In a file with a header: the reference index will be different from the txt if the header and reference chroms are not
//    in the same order

void ref_alt_chroms_compress (Reference ref)
{
    // case: when compressing SAM, BAM or VCF with a header (including BAM header with no SQs - unaligned BAM), 
    // alt chroms were already prepared in ref_contigs_ref_chrom_from_header_chrom
    if (txtheader_get_contigs()) goto just_compress;

    Context *ctx = &z_file->contexts[CHROM];
    uint32_t num_chroms = ctx->nodes.len;
    uint32_t num_contigs = ref->loaded_contigs.len;   // chroms that are in the reference file

    ASSERT (num_chroms >= num_contigs, "expecting num_chroms=%u >= num_contigs=%u", num_chroms, num_contigs);

    uint32_t num_alt_chroms = num_chroms - num_contigs; // chroms that are only in the txt file, not in the reference

    if (!num_alt_chroms) return; // no need for an alt chroms sections as we have none

    buf_alloc (evb, &z_file->alt_chrom_map, 0, num_alt_chroms, AltChrom, 1, "z_file->alt_chrom_map");

    if (flag.show_ref_alts) 
        iprint0 ("\nAlternative chrom indices (output of --show-ref-alts): chroms that are in the file and are mapped to a different name in the reference\n");

    for (uint32_t i=0; i < num_alt_chroms; i++) {
        WordIndex chrom_index = num_contigs + i;
        const CtxNode *chrom_node = ENT (CtxNode, ctx->nodes, chrom_index);
        const char *chrom_name = ENT (const char, ctx->dict, chrom_node->char_index);
        
        WordIndex alt_index = ref_alt_chroms_get_alt_index (ref, chrom_name, chrom_node->snip_len, 0, WORD_INDEX_NONE);

        // an alt_index might be missing for chrom snips like '=' or '*' or sequence-less chroms that don't appear in the reference
        if (alt_index != WORD_INDEX_NONE)           
            NEXTENT (AltChrom, z_file->alt_chrom_map) = (AltChrom){ .txt_chrom = BGEN32 (chrom_index), 
                                                                    .ref_chrom = BGEN32 (alt_index) };
    
        if (flag.show_ref_alts) {
            const CtxNode *alt_node = ENT (CtxNode, ctx->nodes, alt_index);
            const char *alt_name = ENT (const char, ctx->dict, alt_node->char_index);

            iprintf ("In file: '%s' (%d) In reference: '%s' (%d)\n", chrom_name, chrom_index, alt_name, alt_index);
        }
    }

just_compress:
    if (z_file->alt_chrom_map.len) {
        z_file->alt_chrom_map.len *= sizeof (AltChrom);
        zfile_compress_section_data_ex (evb, SEC_REF_ALT_CHROMS, &z_file->alt_chrom_map, 0,0, CODEC_LZMA, SECTION_FLAGS_NONE); // compresses better with LZMA than BZLIB
    }
    
    buf_free (&z_file->alt_chrom_map);
}

void ref_alt_chroms_load (Reference ref)
{
    Section sl = sections_last_sec (SEC_REF_ALT_CHROMS, true);
    if (!sl) return; // we don't have alternate chroms

    zfile_get_global_section (SectionHeader, SEC_REF_ALT_CHROMS, sl, &evb->compressed, "compressed");

    if (flag.show_ref_alts) 
        iprint0 ("\nAlternative chrom indices (output of --show-ref-alts): chroms that are in the txt file and are mapped to a different index in the reference\n");

    evb->compressed.len /= sizeof (AltChrom);
    Context *ctx = &z_file->contexts[CHROM];

    // create mapping user index -> reference index
    buf_alloc (evb, &z_file->alt_chrom_map, 0, ctx->word_list.len, WordIndex, 1, "z_file->alt_chrom_map");
    z_file->alt_chrom_map.len = ctx->word_list.len;

    // initialize with unity mapping
    ARRAY (WordIndex, map, z_file->alt_chrom_map);
    for (uint32_t i=0; i < ctx->word_list.len; i++)
        map[i] = i;

    // the indices of chroms that are NOT in the reference (they are only in the user file), will be mapped to ref chroms
    uint32_t num_ref_contigs = ref->loaded_contigs.len;
    for (uint32_t i=0; i < evb->compressed.len; i++) {
        AltChrom *ent = ENT (AltChrom, evb->compressed, i);
        WordIndex txt_chrom_index = BGEN32 (ent->txt_chrom);
        WordIndex ref_chrom_index = BGEN32 (ent->ref_chrom);

        ASSERT (txt_chrom_index >= 0 && txt_chrom_index < ctx->word_list.len, "txt_chrom_index=%d out of range [0,%d]", txt_chrom_index, (int32_t)ctx->word_list.len);
        ASSERT (!num_ref_contigs /* ref not loaded */ || (ref_chrom_index >= 0 && ref_chrom_index < num_ref_contigs), 
                "ref_chrom_index=%d out of range [0,%u]", ref_chrom_index, (unsigned)ref->loaded_contigs.len);

        map[txt_chrom_index] = ref_chrom_index;

        if (flag.show_ref_alts) {
            const char *chrom_name = ctx_get_words_snip (ctx, txt_chrom_index);
            const char *alt_name   = ctx_get_words_snip (ctx, ref_chrom_index);
            iprintf ("In file: '%s' (%d) In reference: '%s' (%d)\n", chrom_name, txt_chrom_index, alt_name, ref_chrom_index);
        }
    }

    if (flag.show_ref_alts && exe_type == EXE_GENOCAT) exit (EXIT_OK); // in genocat this, not the data

    buf_free (&evb->compressed);
}

// Used in ZIP and in PIZ - only in regions_make_chregs. In PIZ, alternative chroms found in the file data are already mapped in SEC_REF_ALT_CHROMS
// looks in Reference contigs if ref is provided, or in CHROM if ref=NULL
WordIndex ref_alt_chroms_get_alt_index (Reference ref, const char *chrom, unsigned chrom_len, PosType chrom_LN, WordIndex fallback_index)
{
    if (!ref_is_external_loaded (ref)) return WORD_INDEX_NONE;

    WordIndex alt_chrom_index = WORD_INDEX_NONE;
    
    // if LN is known, look for a chromosome with the same length
    if (chrom_LN)
        alt_chrom_index = ref_contigs_get_by_uniq_len (ref, chrom_LN);

    // 22 -> chr22 (1->22, X, Y, M, MT chromosomes)
    else if ((chrom_len == 1 && (IS_DIGIT (chrom[0]) || chrom[0]=='X' || chrom[0]=='Y')) ||
        (chrom_len == 2 && ((IS_DIGIT (chrom[0]) && IS_DIGIT (chrom[1]))))) {

        char chr_chrom[5] = "chr";
        chr_chrom[3] = chrom[0];
        chr_chrom[4] = (chrom_len == 2 ? chrom[1] : 0);

        alt_chrom_index = ref_contigs_get_by_name (ref, chr_chrom, chrom_len+3, true); 
    }

    // M, MT, chrM, chrMT 
    else if ((chrom_len==4 && chrom[0]=='c' && chrom[1]=='h' && chrom[2]=='r' && chrom[3]=='M') || 
             (chrom_len==5 && chrom[0]=='c' && chrom[1]=='h' && chrom[2]=='r' && chrom[3]=='M' && chrom[4]=='T') || 
             (chrom_len==1 && chrom[0]=='M') ||
             (chrom_len==2 && chrom[0]=='M' && chrom[1]=='T')) {

        alt_chrom_index = ref_contigs_get_by_name (ref, "chrMT", 5, true); 
        if (alt_chrom_index == WORD_INDEX_NONE) alt_chrom_index = ref_contigs_get_by_name (ref, "chrM", 4, true); 
        if (alt_chrom_index == WORD_INDEX_NONE) alt_chrom_index = ref_contigs_get_by_name (ref, "MT", 2, true); 
        if (alt_chrom_index == WORD_INDEX_NONE) alt_chrom_index = ref_contigs_get_by_name (ref, "M", 1, true); 
    }
    
    // Chr? or Chr?? -> ? or ??
    else if ((chrom_len == 4 || chrom_len == 5) && !memcmp (chrom, "chr", 3))
        alt_chrom_index = ref_contigs_get_by_name (ref, &chrom[3], chrom_len-3, true); 

    // AC subfield in DESC in reference FASTAs, eg GRCh37/38, eg "GL000207.1" -> "chr18_gl000207_random"
    // https://www.ncbi.nlm.nih.gov/Sequin/acc.html
    else if (ref && chrom_len >= 6 && IS_CLETTER (chrom[0]) && chrom[chrom_len-2]=='.' && IS_DIGIT(chrom[chrom_len-1])) 
        alt_chrom_index = ref_contigs_get_by_accession_number (ref, chrom, chrom_len);

    // final check - see that the LN is as expected (if caller provided a LN)
    if (alt_chrom_index != WORD_INDEX_NONE && chrom_LN && 
        chrom_LN != ENT (RefContig, ref->loaded_contigs, alt_chrom_index)->max_pos)
        alt_chrom_index = WORD_INDEX_NONE;

    return (alt_chrom_index != WORD_INDEX_NONE) ? alt_chrom_index : fallback_index;
}

