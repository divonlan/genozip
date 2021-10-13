// ------------------------------------------------------------------
//   ref_contigs.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

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
#include "refhash.h"
#include "profiler.h"
#include "fasta.h"
#include "contigs.h"
#include "chrom.h"
#include "buffer.h"

#ifdef __linux__ 
extern int strncasecmp (const char *s1, const char *s2, size_t n); // defined in <strings.h>, but file name conflicts with "strings.h" (to do: sort this out in the Makefile)
#endif

static void BGEN_ref_contigs (Buffer *contigs_buf)
{
    for (uint32_t i=0; i < contigs_buf->len; i++) {
        Contig *cn   = ENT (Contig, *contigs_buf, i);
        cn->gpos        = BGEN64 (cn->gpos);
        cn->min_pos     = BGEN64 (cn->min_pos);
        cn->max_pos     = BGEN64 (cn->max_pos);
        cn->char_index  = BGEN64 (cn->char_index);
        cn->snip_len    = BGEN32 (cn->snip_len);
        cn->ref_index   = BGEN32 (cn->ref_index);
    };
}

// in the case of data compressed with the aligner (FASTQ or unaligned SAM) it might not have all the contigs in CHROM
// (we might have some, eg in case of a SAM with mixed aligned and unaligned reads). If this is REF_EXT_STORE, the a 
// SEC_REF_CONTIGS section is created, which refers to the CHROM dict for names, hence we add these chroms. 
void ref_contigs_populate_aligned_chroms (void)
{
    // create sorted index into CHROM
    chrom_index_by_name (CHROM);

    for (uint32_t range_i=0; range_i < gref->ranges.len; range_i++) {
        Range *r = ENT (Range, gref->ranges, range_i);

        bool already_exists = (chrom_get_by_name (STRa(r->chrom_name)) != WORD_INDEX_NONE); // it might exist, eg in a case of a SAM with mixed aligned and unaligned reads
        if (already_exists) continue;

        // add if chrom doesn't already exist ()
        if ((r->num_set = bit_array_num_bits_set (&r->is_set))) {
            
            WordIndex chrom_index = ctx_populate_zf_ctx (CHROM, STRa(r->chrom_name), r->range_i); // add to dictionary

            buf_alloc_255 (evb, &z_file->chrom2ref_map, 0, MAX_(1000, chrom_index+1), WordIndex, CTX_GROWTH, "z_file->chrom2ref_map");

            *ENT (WordIndex, z_file->chrom2ref_map, chrom_index) = r->chrom;
        }
    }
}

static void ref_contigs_show (const Buffer *contigs_buf, bool created)
{
    ARRAY (const Contig, cn, *contigs_buf);

    iprintf ("\nContigs as they appear in the reference%s:\n", created ? " created (note: contig names are as they appear in the txt data, not the reference)" : "");
    for (uint32_t i=0; i < cn_len; i++) {

        const char *chrom_name = ENT (const char, ZCTX(CHROM)->dict, cn[i].char_index);
        bool ext_ref = (flag.reference & REF_ZIP_LOADED) || Z_DT(DT_REF);

        if (ext_ref && created)
            iprintf ("\"%s\" length=%"PRId64" ref_index=%d gpos=%"PRId64" metadata=\"%s\"\n",
                     chrom_name,  cn[i].max_pos, cn[i].ref_index, cn[i].gpos, cn[i].metadata.str);

        else if (ext_ref) 
            iprintf ("\"%s\" length=%"PRId64" ref_index=%d gpos=%"PRId64" %s\n",
                     chrom_name,  cn[i].max_pos, cn[i].ref_index, cn[i].gpos, display_acc_num (&cn[i].metadata.parsed.ac).s);

        else if (cn[i].snip_len)
            iprintf ("i=%d '%s' gpos=%"PRId64" min_pos=%"PRId64" max_pos=%"PRId64" ref_index=%d char_index=%"PRIu64" snip_len=%u\n",
                     i, chrom_name, cn[i].gpos, cn[i].min_pos, cn[i].max_pos, cn[i].ref_index, cn[i].char_index, cn[i].snip_len);
        else
            iprintf ("i=%d (unused - not present in txt data)\n", cn[i].ref_index);
    }
}

static void ref_contigs_compress_do (BufferP created_contigs)
{
    if (flag.show_ref_contigs) ref_contigs_show (created_contigs, true);

    BGEN_ref_contigs (created_contigs);

    created_contigs->len *= sizeof (Contig);
    zfile_compress_section_data_ex (evb, SEC_REF_CONTIGS, created_contigs, 0,0, CODEC_BSC, SECTION_FLAGS_NONE); 

    buf_free (created_contigs);
}

// REF_INTERNAL (SAM/BAM): create and compress contigs - sorted by chrom_index
void ref_contigs_compress_internal (Reference ref)
{
    START_TIMER;

    static Buffer created_contigs = EMPTY_BUFFER;  

    if (!buf_is_alloc (&ZCTX(CHROM)->nodes)) return; // no contigs (note: non-empty SAM/BAM files always have contigs, even if only a "*")

    // the number of contigs is at most the number of chroms - but could be less if some chroms have no sequence
    // in SAM with header and -E, we don't copy the contigs to CHROM, so we take the length from CHROM.
    buf_alloc (evb, &created_contigs, 0, ZCTX(CHROM)->nodes.len, Contig, 1, "created_contigs");

    // create sorted index into CHROM
    chrom_index_by_name (CHROM);

    Contig *last = NULL;

    ConstBufferP contig_metadata = ref_make_get_contig_metadata(); // empty unless this is --make-reference

    for (uint32_t range_i=0; range_i < ref->ranges.len; range_i++) {
        Range *r = ENT (Range, ref->ranges, range_i);

        // chrom_word_index might still be WORD_INDEX_NONE. We get it now from the z_file data
        if (r->chrom == WORD_INDEX_NONE)
            r->chrom = chrom_get_by_name (STRa(r->chrom_name));

        // first range of a contig
        if (!last || r->chrom != last->ref_index) {

            // we assign 64-aligned gpos
            r->gpos = range_i ? ROUNDUP64 ((r-1)->gpos + (r-1)->last_pos - (r-1)->first_pos + 1) : 0;
                
            CtxNode *chrom_node = (r->chrom != WORD_INDEX_NONE) ? ENT (CtxNode, ZCTX(CHROM)->nodes, r->chrom) : NULL;

            NEXTENT (Contig, created_contigs) = (Contig){
                .gpos        = r->gpos, 
                .min_pos     = r->first_pos,
                .max_pos     = r->last_pos,
                .ref_index   = r->chrom,  // index into ZCTX(CHROM)
                .char_index  = chrom_node ? chrom_node->char_index : 0, 
                .snip_len    = chrom_node ? chrom_node->snip_len   : 0
            };

            last = LASTENT (Contig, created_contigs);

            if (flag.make_reference)
                memcpy (last->metadata.str, ENT (ContigMetadata, *contig_metadata, r->chrom), REFCONTIG_MD_LEN);
        }
        else {
            // set gpos of the next range relative to this chrom (note: ranges might not be contiguous)
            r->gpos = (r-1)->gpos + (r->first_pos - (r-1)->first_pos); // note: only the first range in the chrom needs to be 64-aligned            
            last->max_pos = r->last_pos; // update
        }
    }
    
    COPY_TIMER_VB (evb, ref_contigs_compress); // we don't count the compression itself and the disk writing

    ref_contigs_compress_do (&created_contigs);
}

// create and compress contigs (when --make-reference or when compressing with REF_EXT_STORE or REF_INTERNAL)
void ref_contigs_compress_ext_store (Reference ref)
{
    START_TIMER;
    ContextP ctx = ZCTX(CHROM);
    uint64_t num_chroms = ctx->nodes.len;

    static Buffer created_contigs = EMPTY_BUFFER;          

    if (!buf_is_alloc (&ctx->nodes)) return; // no contigs
    ARRAY (CtxNode, nodes, ctx->nodes);

    // NOTE: ref_get_range_by_chrom access ranges by chrom, but ref_initialize_loaded_ranges (incorrectly) allocates ranges.len by ref->ctgs.contigs.len 
    // (for non REF_INTERNAL) therefore, we must have ref_contigs aligned with CHROM
    
    buf_alloc_zero (evb, &created_contigs, 0, num_chroms, Contig, 1, "created_contigs");
    created_contigs.len = num_chroms;
    ARRAY (Contig, cn, created_contigs);

    for (uint64_t i=0; i < num_chroms; i++) {
        cn[i].ref_index  = WORD_INDEX_NONE; // initialize (not all CHROMs have contigs, eg "*" in SAM)
        cn[i].char_index = nodes[i].char_index;
        cn[i].snip_len   = nodes[i].snip_len;
        cn[i].ref_index  = i; // ref_index always refers to the CHROM index of the file in which it is stored (eg REF_CONTIG in References files, RNAME in SAM etc)
        cn[i].gpos       = 0; // ref_initialize_ranges overlays even empty contigs on the genome, so GPOS must be valid number
    }

    // create sorted index into CHROM
    chrom_index_by_name (CHROM);

    for (uint32_t range_i=0; range_i < ref->ranges.len; range_i++) {
        Range *r = ENT (Range, ref->ranges, range_i);

        // if gpos of the first range in this contig has been increased due to removing flanking regions and is no longer 64-aligned,
        // we start the contig a little earlier to be 64-aligned (note: this is guaranteed to be within the original contig before 
        // compacting, because the original contig had a 64-aligned gpos)
        PosType delta = r->gpos % 64;

        WordIndex chrom = *ENT (WordIndex, z_file->ref2chrom_map, r->chrom); // the CHROM corresponding to this ref_index, even if a different version of the chrom name 
        if (chrom == WORD_INDEX_NONE) continue; // this contig was not used in the data 
        
        cn[chrom].gpos    = r->gpos - delta;
        cn[chrom].min_pos = r->first_pos - delta;
        cn[chrom].max_pos = r->last_pos;
    }

    COPY_TIMER_VB (evb, ref_contigs_compress); // we don't count the compression itself and the disk writing

    ref_contigs_compress_do (&created_contigs);
}

// read and uncompress a contigs section (when pizzing the reference file or pizzing a data file with a stored reference)
void ref_contigs_load_contigs (Reference ref)
{
    Section sl = sections_last_sec (SEC_REF_CONTIGS, true);
    if (!sl) return; // section doesn't exist

    zfile_get_global_section (SectionHeader, SEC_REF_CONTIGS, sl, &ref->ctgs.contigs, "ContigPkg->contigs");

    ref->ctgs.contigs.len /= sizeof (Contig);
    BGEN_ref_contigs (&ref->ctgs.contigs);

    ASSERT0 (ZCTX(CHROM)->dict.len, "CHROM dictionary is empty");

    buf_copy (evb, &ref->ctgs.dict, &ZCTX(CHROM)->dict, char, 0, 0, "ContigPkg->dict");

    contigs_create_index (&ref->ctgs, SORT_BY_NAME | SORT_BY_AC);

    if (flag.show_ref_contigs) { 
        ref_contigs_show (&ref->ctgs.contigs, false);
        if (exe_type == EXE_GENOCAT) exit_ok();  // in genocat this, not the data
    }

    // check chromosome naming - assume the first contig is a chromosome
    const char *c = ref->ctgs.dict.data;
    ref->chrom_style = (IS_DIGIT(c[0]) && (!c[1] || (IS_DIGIT(c[1]) && !c[2]))) ? CHROM_STYLE_22
                     : (c[0]=='c' && c[1]=='h' && c[2]=='r')                    ? CHROM_STYLE_chr22
                     :                                                            CHROM_STYLE_UNKNOWN;
}

// binary search for this chrom in ref->ctgs. we count on gcc tail recursion optimization to keep this fast.
// looks in Reference contigs if ref is provided, or in CHROM if ref=NULL
WordIndex ref_contigs_get_by_name (const Reference ref, const char *chrom_name, uint32_t chrom_name_len, bool alt_ok, bool soft_fail)
{
    return alt_ok ? contigs_get_matching (&ref->ctgs, STRa(chrom_name), 0, false, NULL)
                  : contigs_get_by_name (&ref->ctgs, STRa(chrom_name));
}

// ZIP: gets matching contig in the reference (same or alt). Returns ref_index if a match was found, or WORD_INDEX_NONE
WordIndex ref_contigs_get_matching (const Reference ref, PosType LN, STRp(txt_chrom), // in
                                    STRp(*ref_contig),          // out - contig name
                                    bool strictly_alt,          // in - only tests for differnet names
                                    bool *is_alt,               // out - true if ref_contig is an alternate chrom different than txt_chrom
                                    int32_t *chrom_name_growth) // optional out
{
    WordIndex ref_index = contigs_get_matching (&ref->ctgs, STRa(txt_chrom), 0, strictly_alt, is_alt);

    if (ref_index != WORD_INDEX_NONE) {
        *ref_contig = contigs_get_name (&ref->ctgs, ref_index, ref_contig_len);
        if (chrom_name_growth) *chrom_name_growth = (int32_t)*ref_contig_len - (int32_t)txt_chrom_len;
    }

    return ref_index;
}

const char *ref_contigs_get_name (Reference ref, WordIndex ref_index, unsigned *contig_name_len /* optional */)
{
    return contigs_get_name (&ref->ctgs, ref_index, contig_name_len);
}

const char *ref_contigs_get_name_by_ref_index (Reference ref, WordIndex ref_index, const char **snip, uint32_t *snip_len)
{
    const Contig *contig = ref_contigs_get_contig_by_ref_index (ref, ref_index, false);
    const char *my_snip = ENT (const char, ref->ctgs.dict, contig->char_index);
    
    if (snip) *snip = my_snip;
    if (snip_len) *snip_len = contig->snip_len;

    return my_snip; 
}

void ref_contigs_generate_data_if_denovo (Reference ref)
{
        // copy data from the reference FASTA's CONTIG context, so it survives after we finish reading the reference and close z_file
    Context *chrom_ctx = ZCTX(CHROM);

    ASSINP (flag.reference == REF_INTERNAL || (buf_is_alloc (&chrom_ctx->dict) && buf_is_alloc (&chrom_ctx->word_list)),
            "Error: cannot use %s as a reference as it is missing a CONTIG dictionary", z_name);

    // we copy from the z_file context after we completed segging a file
    // note: we can't rely on chrom_ctx->nodes for the correct order as vb_i=1 resorted. rather, by mimicking the
    // word_list generation as done in PIZ, we guarantee that we will get the same ref_index
    // in case of multiple bound files, we re-do this in every file in case of additional chroms (not super effecient, but good enough because the context is small)
    contigs_free (&ref->ctgs);
    contigs_build_contig_pkg_from_ctx (&ref->ctgs, chrom_ctx, SORT_BY_NAME);
}

// ZIP SAM/BAM and VCF: verify that we have the specified chrom (name & last_pos) loaded from the reference at the same index. 
// called by sam_header_add_contig, vcf_header_consume_contig. returns true if contig is in reference
WordIndex ref_contigs_ref_chrom_from_header_chrom (Reference ref, const char *chrom_name, unsigned chrom_name_len, 
                                                   PosType *hdr_LN) // if 0, set from reference, otherwise verify
{               
    WordIndex ref_contig_index = ref_contigs_get_by_name (ref, chrom_name, chrom_name_len, true, true); // including alts
        
    // if its not found, we ignore it. sequences that have this chromosome will just be non-ref
    if (ref_contig_index == WORD_INDEX_NONE) {
        if (command == ZIP)
            WARN_ONCE ("FYI: header of %s has contig '%.*s' (and maybe others, too), missing in %s. No harm.",
                        txt_file->basename, chrom_name_len, chrom_name, ref->filename);
        return WORD_INDEX_NONE;
    }

    if (buf_is_alloc (&ref->ranges)) { // it is not allocated in --sex/coverage
        PosType ref_LN = ENT (Range, ref->ranges, ref_contig_index)->last_pos; // get from ranges because Contig.LN=0 - we don't populate it at reference creation
        
        // case: file header is missing length, update from reference
        if (! *hdr_LN) *hdr_LN = ref_LN;

        // case: file header specifies length - it must be the same as the reference
        else if (*hdr_LN != ref_LN) {
            const char *ref_contig_name = ref_contigs_get_name (ref, ref_contig_index, NULL); // might be different that chrom_name if it matches an alt_name
        
            char fmt[200]; // no unbound names in fmt
            sprintf (fmt, "Error: wrong reference file: %%s has a \"%s\", but in %%s '%%s' has LN=%%"PRId64, DTPT (hdr_contigs));
            ASSINP (false, fmt, txt_name, chrom_name_len, chrom_name, *hdr_LN, ref->filename, ref_contig_name, ref_LN);
        }
    }

    return ref_contig_index;
}

// get length of contig according to ref_contigs (loaded or stored reference) - used in piz when generating
// a file header of a translated data type eg SAM->BAM or ME23->VCF
PosType ref_contigs_get_contig_length (const Reference ref,
                                       WordIndex ref_index, // option 1 
                                       const char *chrom_name, unsigned chrom_name_len, // option 2
                                       bool enforce)
{
    if (ref_index == WORD_INDEX_NONE && chrom_name)
        ref_index = ref_contigs_get_by_name (ref, STRa(chrom_name), true, true);

    ASSERT (!enforce || (ref_index >= 0 && ref_index < ref->ranges.len), "contig=\"%.*s\" ref_index=%d not found in reference contigs (range=[0,%d])",
            chrom_name_len, chrom_name ? chrom_name : "", ref_index, (int)ref->ranges.len-1);

    if (ref_index == WORD_INDEX_NONE) return -1; // chrom_name not found in ref_contigs
    
    // get info as it appears in reference
    return ENT (Range, ref->ranges, ref_index)->last_pos;
}

// get contig by ref_index, by binary searching for it
static const Contig *ref_contigs_get_contig_do (const Reference ref, WordIndex ref_index, int32_t start_i, int32_t end_i)
{
    int32_t mid_i = (start_i + end_i) / 2;
    Contig *mid_rc = ENT (Contig, ref->ctgs.contigs, mid_i);
    if (mid_rc->ref_index == ref_index) return mid_rc; // found
    if (start_i >= end_i) return NULL; // not found

    // continue searching
    if (ref_index < mid_rc->ref_index)
        return ref_contigs_get_contig_do (ref, ref_index, start_i, mid_i-1);
    else    
        return ref_contigs_get_contig_do (ref, ref_index, mid_i+1, end_i);
}

const Contig *ref_contigs_get_contig_by_ref_index (const Reference ref, WordIndex ref_index, bool soft_fail)
{
    ASSERTISALLOCED (ref->ctgs.contigs);

    const Contig *rc = ref_contigs_get_contig_do (ref, ref_index, 0, ref->ctgs.contigs.len-1);

    ASSERT (rc || soft_fail, "cannot find contig for ref_index=%d", ref_index);

    return rc;
}

PosType ref_contigs_get_genome_nbases (const Reference ref)
{
    if (!ref->ctgs.contigs.len) return 0;

    // note: contigs are sorted by chrom and then pos within chrom, NOT by gpos! (chroms might have been added out of order during reference creation)
    ARRAY (const Contig, rc, ref->ctgs.contigs);
    const Contig *rc_with_largest_gpos = NULL;

    for (uint64_t i=0; i < rc_len; i++)  {

        if (!rc[i].min_pos && !rc[i].max_pos) continue; // not real contig, eg "*" in SAM

        if (!rc_with_largest_gpos || rc[i].gpos > rc_with_largest_gpos->gpos) rc_with_largest_gpos = &rc[i];
    }
    
    ASSERT (rc_with_largest_gpos, "There are %"PRIu64" elements in the contigs array, but none are actual contigs", ref->ctgs.contigs.len);

    // note: contigs with gpos beyond 4Gb will not be used by the aligner

    return rc_with_largest_gpos->gpos + (rc_with_largest_gpos->max_pos - rc_with_largest_gpos->min_pos + 1);
}

WordIndex ref_contig_get_by_gpos (const Reference ref, PosType gpos,
                                  PosType *pos) // optional out, POS within the CHROM matching gpos
{
    // note: contigs are sorted by chrom and then pos within chrom, NOT by gpos! (chroms might have been added out of order during reference creation)
    for (WordIndex ref_index=0 ; ref_index < ref->ctgs.contigs.len; ref_index++) {
        Contig *rc = ENT (Contig, ref->ctgs.contigs, ref_index);
        if (gpos >= rc->gpos && gpos <= rc->gpos + (rc->max_pos - rc->min_pos)) {
            if (pos) {
                ASSERT (ref_index < ref->ranges.len, "Unexpected ref_index=%d >= ref->ranges.len=%"PRIu64, ref_index, ref->ranges.len);  
                const Range *r = ENT (Range, ref->ranges, ref_index);
                
                *pos = r->first_pos + gpos - r->gpos;
            }
            return ref_index;
        }
    }

    return WORD_INDEX_NONE; // not found
}
