// ------------------------------------------------------------------
//   map_chrom2ref.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "context.h"
#include "file.h"
#include "zfile.h"
#include "endianness.h"
#include "contigs.h"
#include "reference.h"

// ZIP of a file with an external reference: 
// The CHROM dictionary includes: first - the txtheader contig names, then reference file contig names that are not in txtheader, 
// and finally names encountered in the file that are not in either the header or the reference (see zip_prepopulate_contig_ctxs()). 
// Since the CHROM context is prepopulated from the txtheader and the reference, often not all the entries are used by the file data.
//
// Here, we create a mapping between those chrom entries that are used (count>0) and the index of the contig in the reference file
// against which this chrom data was compressed. We rely on contigs_get_matching() returning the same result here as it did in Seg.

void map_chrom2ref_compress (Reference ref)
{
    Context *ctx = ZCTX(CHROM);
    ARRAY (CtxNode, nodes, ctx->nodes);

    buf_alloc (evb, &z_file->chrom2ref_map, 0, nodes_len, ChromRefMap, 1, "z_file->chrom2ref_map");

    if (flag.show_ref_alts) 
        iprint0 ("\nAlternative chrom indices (output of --show-ref-alts): chroms that are in the file and are mapped to a different index in the reference\n");

    for (uint32_t chrom_index=0; chrom_index < nodes_len; chrom_index++) {

        if (! *ENT (int64_t, ctx->counts, chrom_index)) continue; // case: this chrom_name doesn't appear in the file data - no need to map it

        const char *chrom_name = ENT (const char, ctx->dict, nodes[chrom_index].char_index);

        bool is_alt; // using is_alt = search for exact match too
        ConstContigPkgP ctgs = ref_get_ctgs (ref); 
        WordIndex ref_index = contigs_get_matching (ctgs, chrom_name, nodes[chrom_index].snip_len, 0, &is_alt);

        // an alt_index might be missing for chrom snips like '=' or '*' or SAM sequence-less chroms that don't appear in the reference
        if (ref_index != WORD_INDEX_NONE) { 
            
            if (chrom_index != ref_index) // chrom_index==ref_index for reference contigs in the absence of header contigs or if header contigs are in the same order of ref contigs - see zip_prepopulate_contig_ctxs
                NEXTENT (ChromRefMap, z_file->chrom2ref_map) = 
                    (ChromRefMap){ .txt_chrom = BGEN32 (chrom_index), .ref_chrom = BGEN32 (ref_index) };
    
            if (flag.show_ref_alts) {
                const char *ref_name = ref_contigs_get_name (ref, ref_index, NULL);
                iprintf ("In file: '%s' (%d) In reference: '%s' (%d) %s\n", chrom_name, chrom_index, ref_name, ref_index,
                         chrom_index != ref_index ? "INDEX_CHANGE" : "");
            }
        }
    }

    if (z_file->chrom2ref_map.len) {
        z_file->chrom2ref_map.len *= sizeof (ChromRefMap);
        zfile_compress_section_data_ex (evb, SEC_CHROM2REF_MAP, &z_file->chrom2ref_map, 0,0, CODEC_LZMA, SECTION_FLAGS_NONE); // compresses better with LZMA than BZLIB
    }
    
    buf_free (&z_file->chrom2ref_map);
}

void map_chrom2ref_load (Reference ref)
{
    Section sl = sections_last_sec (SEC_CHROM2REF_MAP, true);
    if (!sl) return; // we don't have alternate chroms

    zfile_get_global_section (SectionHeader, SEC_CHROM2REF_MAP, sl, &evb->compressed, "compressed");

    if (flag.show_ref_alts) 
        iprint0 ("\nAlternative chrom indices (output of --show-ref-alts): chroms that are in the txt file and are mapped to a different index in the reference\n");

    evb->compressed.len /= sizeof (ChromRefMap);
    Context *ctx = &z_file->contexts[CHROM];

    // create mapping user index -> reference index
    buf_alloc (evb, &z_file->chrom2ref_map, 0, ctx->word_list.len, WordIndex, 1, "z_file->chrom2ref_map");
    z_file->chrom2ref_map.len = ctx->word_list.len;

    // initialize with unity mapping
    ARRAY (WordIndex, map, z_file->chrom2ref_map);
    for (uint32_t i=0; i < ctx->word_list.len; i++)
        map[i] = i;

    // the indices of chroms that are NOT in the reference (they are only in the user file), will be mapped to ref chroms
    ConstContigPkgP ctgs = ref_get_ctgs (ref); 
    uint32_t num_ref_contigs = ctgs->contigs.len;
    for (uint32_t i=0; i < evb->compressed.len; i++) {
        ChromRefMap *ent = ENT (ChromRefMap, evb->compressed, i);
        WordIndex txt_chrom_index = BGEN32 (ent->txt_chrom);
        WordIndex ref_chrom_index = BGEN32 (ent->ref_chrom);

        ASSERT (txt_chrom_index >= 0 && txt_chrom_index < ctx->word_list.len, "txt_chrom_index=%d out of range [0,%d]", txt_chrom_index, (int32_t)ctx->word_list.len);
        ASSERT (!num_ref_contigs /* ref not loaded */ || (ref_chrom_index >= 0 && ref_chrom_index < num_ref_contigs), 
                "ref_chrom_index=%d out of range [0,%u]", ref_chrom_index, num_ref_contigs-1);

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
