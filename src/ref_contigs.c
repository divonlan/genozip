// ------------------------------------------------------------------
//   ref_contigs.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "ref_private.h"
#include "zfile.h"
#include "seg.h"
#include "refhash.h"
#include "chrom.h"

#ifdef __linux__ 
extern int strncasecmp (rom s1, rom s2, size_t n); // defined in <strings.h>, but file name conflicts with "strings.h" (to do: sort this out in the Makefile)
#endif

static void BGEN_ref_contigs_not_compacted (BufferP contigs_buf)
{
    for (uint32_t i=0; i < contigs_buf->len; i++) {
        Contig *cn   = B(Contig, *contigs_buf, i);
        cn->gpos        = BGEN64 (cn->gpos);
        cn->min_pos     = BGEN64 (cn->min_pos);
        cn->max_pos     = BGEN64 (cn->max_pos);
        cn->char_index  = BGEN64 (cn->char_index);
        cn->snip_len    = BGEN32 (cn->snip_len);
        cn->ref_index   = BGEN32 (cn->ref_index);
    };
}

// main thread from zip_write_global_area: in the case of data compressed with the aligner (FASTQ or unaligned reads in SAM) 
// it might not have all the contigs in CHROM (we might have some, eg in case of a SAM with mixed aligned and unaligned reads). 
// If this is REF_EXT_STORE, a SEC_REF_CONTIGS is created, which refers to the CHROM dict for names, hence we add these chroms. 
void ref_contigs_populate_aligned_chroms (void)
{
    // create sorted index into CHROM
    chrom_index_by_name (CHROM);

    for_buf (Range, r, gref->ranges) {
        r->num_set = bits_num_set_bits (&r->is_set);
        if (!r->num_set) continue;

        WordIndex chrom_index = chrom_get_by_name (STRa(r->chrom_name)); 
        decl_zctx(CHROM);

        // contig r->chrom_name was used by aligner, but it is not in CHROM.dict (it could have been in CHROM in a mixed mapped/mapped SAM file)
        if (chrom_index == WORD_INDEX_NONE) 
            // add it to the the CHROM.dict
            chrom_index = ctx_populate_zf_ctx (CHROM, STRa(r->chrom_name), r->range_id); 

        // make sure count is at least 1 for ref_contigs_compress_stored to store 
        buf_alloc_zero (evb, &zctx->counts, 0, MAX_(chrom_index + 1, gref->ranges.len32), uint64_t, CTX_GROWTH, "zctx->counts");
        *B64(zctx->counts, chrom_index) = MAX_(*B64(zctx->counts, chrom_index), 1);
    }
}

static void ref_contigs_show (ConstBufferP contigs_buf, bool created)
{
    iprintf ("\nContigs as they appear in the reference%s:\n", created ? " created (note: contig names are as they appear in the txt data, not the reference)" : "");
    for_buf2 (const Contig, cn, i, *contigs_buf) {

        if (!cn->max_pos) continue; // unused contig

        rom chrom_name = B(const char, ZCTX(CHROM)->dict, cn->char_index);
        bool ext_ref = (IS_ZIP && (flag.reference & REF_ZIP_LOADED)) || Z_DT(REF);

        if (ext_ref && created)
            iprintf ("\"%s\" length=%"PRId64" ref_index=%d gpos=%"PRId64" metadata=\"%s\"\n",
                     chrom_name, cn->max_pos, cn->ref_index, cn->gpos, cn->metadata.str);

        else if (ext_ref) 
            iprintf ("\"%s\" length=%"PRId64" ref_index=%d gpos=%"PRId64" %s\n",
                     chrom_name, cn->max_pos, cn->ref_index, cn->gpos, display_acc_num (&cn->metadata.parsed.ac).s);

        else if (cn->snip_len)
            iprintf ("[%d] \"%s\" gpos=%"PRId64" min_pos=%"PRId64" max_pos=%"PRId64" ref_index=%d char_index=%"PRIu64" snip_len=%u\n",
                     i, chrom_name, cn->gpos, cn->min_pos, cn->max_pos, cn->ref_index, cn->char_index, cn->snip_len);
        else
            iprintf ("[%d] (unused - not present in txt data)\n", cn->ref_index);
    }
}

static void ref_contigs_compress_do (BufferP created_contigs, bool sequential_ref_index, bool compacted)
{
    if (!compacted) // note: compacted contigs are already big endian
        BGEN_ref_contigs_not_compacted (created_contigs);

    created_contigs->len *= compacted ? sizeof (CompactContig) : sizeof (Contig);
    zfile_compress_section_data_ex (evb, NULL, SEC_REF_CONTIGS, created_contigs, 0, 0, CODEC_BSC, 
                                    (struct FlagsRefContigs){ .sequential_ref_index = sequential_ref_index,
                                                              .compacted            = compacted             }, 
                                    NULL); 

    buf_free (*created_contigs);
}

// Used by: 1. REF_INTERNAL (SAM/BAM): create and compress contigs - sorted by chrom_index 2. --make-reference
void ref_contigs_compress_ref_make (Reference ref)
{
    START_TIMER;

    static Buffer created_contigs = {};  

    if (!buf_is_alloc (&ZCTX(CHROM)->nodes)) return; // no contigs (note: non-empty SAM/BAM files always have contigs, even if only a "*")

    // the number of contigs is at most the number of chroms - but could be less if some chroms have no sequence
    // in SAM with header and -E, we don't copy the contigs to CHROM, so we take the length from CHROM.
    buf_alloc_zero (evb, &created_contigs, 0, ZCTX(CHROM)->nodes.len, Contig, 1, "created_contigs");

    // create sorted index into CHROM
    chrom_index_by_name (CHROM);

    Contig *last = NULL;

    ConstBufferP contig_metadata = ref_make_get_contig_metadata(); // empty unless this is --make-reference

    Range *prev_r = NULL;
    for_buf (Range, r, ref->ranges) {

        if (r->first_pos > r->last_pos) continue; // range not actually used

        // chrom_word_index might still be WORD_INDEX_NONE. We get it now from the z_file data
        if (r->chrom == WORD_INDEX_NONE) {
            r->chrom = chrom_get_by_name (STRa(r->chrom_name));
            
            if (r->chrom == WORD_INDEX_NONE) {
                buf_dump_to_file ("RNAME.dict", &ZCTX(CHROM)->dict, 1, false, true, true, false);
                ABORT ("Unable to find chrom index for contig \"%.*s\"", STRf(r->chrom_name));
            }
        }

        // first range of a contig
        if (!last || r->chrom != last->ref_index) {

            // we assign 64-aligned gpos
            r->gpos = prev_r ? ROUNDUP64 (prev_r->gpos + prev_r->last_pos - prev_r->first_pos + 1) : 0;
                
            BNXT (Contig, created_contigs) = (Contig){
                .gpos        = r->gpos, 
                .min_pos     = r->first_pos,
                .max_pos     = r->last_pos,
                .ref_index   = r->chrom,  // word_index in ZCTX(CHROM)
            };

            last = BLST (Contig, created_contigs);

            if (flag.make_reference)
                memcpy (last->metadata.str, B(ContigMetadata, *contig_metadata, r->chrom), REFCONTIG_MD_LEN);
        }
        else {
            // set gpos of the next range relative to this chrom (note: ranges might not be contiguous)
            r->gpos = prev_r->gpos + (r->first_pos - prev_r->first_pos); // note: only the first range in the chrom needs to be 64-aligned            
            last->max_pos = r->last_pos; // update
        }

        prev_r = r;
    }
    
    COPY_TIMER_EVB (ref_contigs_compress); // we don't count the compression itself and the disk writing

    if (flag.show_ref_contigs) ref_contigs_show (&created_contigs, true);

    ref_contigs_compress_do (&created_contigs, false, false);

    buf_destroy (created_contigs);
}

// convert Contig to CompactContig if possible, ahead of writing to disk
static bool ref_contigs_compress_compact (BufferP created_contigs)
{
    // verify that contigs can be compacted
    for_buf2 (Contig, cn, i, *created_contigs) 
        if (cn->min_pos > 0xffffffffULL || cn->max_pos > 0xffffffffULL) 
            return false;

    // do the compacting
    CompactContig *next = B1ST (CompactContig, *created_contigs); // replace in-place
    for_buf2 (Contig, cn, ref_index, *created_contigs) {
        if (!cn->max_pos) continue; // skip unused (cannot happen if sequential_ref_index)

        CompactContig cc = { .ref_index = BGEN32 (ref_index),
                             .min_pos   = BGEN32 ((uint32_t)cn->min_pos),
                             .max_pos   = BGEN32 ((uint32_t)cn->max_pos),
                             .gpos      = BGEN64 (cn->gpos)              };                                
        *next++ = cc;
    }

    created_contigs->len = BNUM (*created_contigs, next);

    return true;
}

// create and compress contigs (when compressing with REF_EXT_STORE or REF_INTERNAL)
void ref_contigs_compress_stored (Reference ref)
{
    START_TIMER;
    decl_zctx (CHROM);

    ARRAY (uint64_t, chrom_counts, ZCTX(CHROM)->counts);

    // NOTE: ref_get_range_by_chrom access ranges by chrom, but ref_initialize_loaded_ranges (incorrectly) allocates ranges.len by ref->ctgs.contigs.len 
    // (for non REF_INTERNAL) therefore, we must have ref_contigs aligned with CHROM
    
    static Buffer created_contigs = {};          
    ARRAY_alloc (Contig, cn, chrom_counts_len, true, created_contigs, evb, "created_contigs");

    // create sorted index into CHROM
    chrom_index_by_name (CHROM);
    
    bool has_some_contigs = false;

    for_buf2 (Range, r, range_i, ref->ranges) {
        // if gpos of the first range in this contig has been increased due to removing flanking regions and is no longer 64-aligned,
        // we start the contig a little earlier to be 64-aligned (note: this is guaranteed to be within the original contig before 
        // compacting, because the original contig had a 64-aligned gpos)
        PosType64 delta = r->gpos % 64;

        WordIndex chrom = IS_REF_EXT_STORE ? *B(WordIndex, zctx->ref2chrom_map, r->chrom) : r->chrom; // the CHROM corresponding to this ref_index, even if a different version of the chrom name 

        // don't store min_pos/max_pos/gpos (better compression) is this contig was not used explicitly (note: aligner doesn't use contigs, but GPOS)
        if (chrom == WORD_INDEX_NONE || !chrom_counts[chrom] || 
            ((IS_REF_INTERNAL || IS_REF_EXT_STORE) && !r->is_set.nbits)) continue; 

        cn[chrom].gpos      = r->gpos - delta;
        cn[chrom].min_pos   = r->first_pos - delta;
        cn[chrom].max_pos   = r->last_pos;

        if (flag.show_ref_contigs) {
            cn[chrom].char_index = ctx_get_char_index_of_snip (zctx, STRa(r->chrom_name), false);
            cn[chrom].snip_len   = r->chrom_name_len;
        }

        has_some_contigs = true;
    }

    if (flag.show_ref_contigs) {
        ref_contigs_show (&created_contigs, true); // note: cannot call this after compacting

        for_buf (Contig, cn, created_contigs) {
            cn->char_index = 0; // these don't need to be written to the file
            cn->snip_len   = 0;        
        }
    }

    // compact contigs if possible - 20 bytes per contigs vs 132
    bool is_compacted = ref_contigs_compress_compact (&created_contigs);

    if (has_some_contigs) // Could be false, e.g. in FASTQ if no read was successfully aligned to the reference
        ref_contigs_compress_do (&created_contigs, true, is_compacted); // note: sequential_ref_index introduced v14.0.0, is_compacted introduce v14.0.10.

    buf_destroy (created_contigs);
    
    COPY_TIMER_EVB (ref_contigs_compress); 
}

static void ref_contigs_load_set_contig_names (Reference ref)
{
    decl_zctx (CHROM);
    ASSERT0 (zctx->dict.len, "CHROM dictionary is empty");

    buf_copy (evb, &ref->ctgs.dict, &zctx->dict, char, 0, 0, "ContigPkg->dict");

    ARRAY (CtxWord, chrom, zctx->word_list);
    ARRAY (Contig, contig, ref->ctgs.contigs);

    for (uint32_t i=0; i < contig_len; i++) {
        if (!contig[i].max_pos) continue;

        WordIndex chrom_index = contig[i].ref_index;
        ASSERTINRANGE (chrom_index, 0, chrom_len);
        contig[i].char_index = chrom[chrom_index].char_index;
        contig[i].snip_len   = chrom[chrom_index].snip_len;
    }
}

static void ref_contigs_load_uncompact (BufferP compact, BufferP contigs)
{
    compact->len /= sizeof (CompactContig);

    ARRAY_alloc (Contig, cn, 
                 BGEN32 (BLST(CompactContig, *compact)->ref_index) + 1, 
                 true, *contigs, evb, "ContigPkg->contigs");

    for_buf (CompactContig, compact_cn, *compact) {
        uint32_t ref_index = BGEN32 (compact_cn->ref_index); 
        cn[ref_index] = (Contig){
            .ref_index = ref_index, // chrom
            .min_pos   = BGEN32 (compact_cn->min_pos),
            .max_pos   = BGEN32 (compact_cn->max_pos),
            .gpos      = BGEN64 (compact_cn->gpos)
        };
    }
}

// read and uncompress a contigs section (when pizzing the reference file or pizzing a data file with a stored reference)
void ref_contigs_load_contigs (Reference ref)
{
    Section sl = sections_last_sec (SEC_REF_CONTIGS, SOFT_FAIL);
    if (!sl) return; // section doesn't exist

    if (sl->flags.ref_contigs.compacted) { // non-reference-file generated in v14.0.10 or later
        zfile_get_global_section (SectionHeader, sl, &evb->scratch, "scratch");

        ref_contigs_load_uncompact (&evb->scratch, &ref->ctgs.contigs);

        buf_free (evb->scratch);
    }
    else {
        zfile_get_global_section (SectionHeader, sl, &ref->ctgs.contigs, "ContigPkg->contigs");

        ref->ctgs.contigs.len /= sizeof (Contig);
        BGEN_ref_contigs_not_compacted (&ref->ctgs.contigs);
    }

    // case: ref_index is set to 0 on disk to improve compression, we re-created it here (sequentially)
    // note: we do this even if compacted, since we also need ref_index set for the skipped contigs
    if (sl->flags.ref_contigs.sequential_ref_index)
        for_buf2 (Contig, cn, ref_index, ref->ctgs.contigs) 
            cn->ref_index = ref_index; 

    // get contig names from CHROM 
    ref_contigs_load_set_contig_names (ref);

    contigs_create_index (&ref->ctgs, SORT_BY_NAME | SORT_BY_AC);

    if (flag.show_ref_contigs) { 
        ref_contigs_show (&ref->ctgs.contigs, false);
        if (is_genocat) exit_ok;  // in genocat this, not the data
    }

    // check chromosome naming - assume the first contig is a chromosome
    rom c = B1STc(ref->ctgs.dict);
    ref->chrom_style = (IS_DIGIT(c[0]) && (!c[1] || (IS_DIGIT(c[1]) && !c[2]))) ? CHROM_STYLE_22
                     : (c[0]=='c' && c[1]=='h' && c[2]=='r')                    ? CHROM_STYLE_chr22
                     :                                                            CHROM_STYLE_UNKNOWN;
}

// binary search for this chrom in ref->ctgs. we count on gcc tail recursion optimization to keep this fast.
// looks in Reference contigs if ref is provided, or in CHROM if ref=NULL
WordIndex ref_contigs_get_by_name (const Reference ref, rom chrom_name, uint32_t chrom_name_len, bool alt_ok, FailType soft_fail)
{
    return alt_ok ? contigs_get_matching (&ref->ctgs, STRa(chrom_name), 0, false, NULL)
                  : contigs_get_by_name (&ref->ctgs, STRa(chrom_name));
}

// ZIP: gets matching contig in the reference (same or alt). Returns ref_index if a match was found, or WORD_INDEX_NONE
WordIndex ref_contigs_get_matching (const Reference ref, PosType64 LN, STRp(txt_chrom), // in
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

rom ref_contigs_get_name (Reference ref, WordIndex ref_index, unsigned *contig_name_len /* optional */)
{
    return contigs_get_name (&ref->ctgs, ref_index, contig_name_len);
}

rom ref_contigs_get_name_by_ref_index (Reference ref, WordIndex ref_index, pSTRp(snip), PosType64 *gpos)
{
    const Contig *contig = ref_contigs_get_contig_by_ref_index (ref, ref_index, HARD_FAIL);
    rom my_snip = B(const char, ref->ctgs.dict, contig->char_index);
    
    if (snip) *snip = my_snip;
    if (snip_len) *snip_len = contig->snip_len;
    if (gpos) *gpos = contig->gpos;

    return my_snip; 
}

// ZIP SAM/BAM and VCF: verify that we have the specified chrom (name & last_pos) loaded from the reference at the same index. 
// called by sam_header_add_contig, vcf_header_consume_contig. returns true if contig is in reference
WordIndex ref_contigs_ref_chrom_from_header_chrom (Reference ref, STRp(chrom_name), 
                                                   PosType64 *hdr_LN) // if 0, set from reference, otherwise verify
{               
    WordIndex ref_contig_index = ref_contigs_get_by_name (ref, STRa(chrom_name), true, true); // including alts
        
    // if its not found, we ignore it. sequences that have this chromosome will just be non-ref
    if (ref_contig_index == WORD_INDEX_NONE) {
        if (IS_ZIP)
            WARN_ONCE ("FYI: header of %s has contig \"%.*s\" (and maybe others, too), missing in %s. If the file contains many %ss with this contig, it might compress a bit less than when using a reference file that contains all contigs.",
                        txt_file->basename, STRf(chrom_name), ref->filename, DTPT(line_name));
        return WORD_INDEX_NONE;
    }

    if (buf_is_alloc (&ref->ranges)) { // it is not allocated in --coverage
        PosType64 ref_LN = B(Range, ref->ranges, ref_contig_index)->last_pos; // get from ranges because Contig.LN=0 - we don't populate it at reference creation
        
        // case: file header is missing length, update from reference
        if (! *hdr_LN) *hdr_LN = ref_LN;

        // case: file header specifies length - it must be the same as the reference
        else if (*hdr_LN != ref_LN && IS_ZIP) { // note: in PIZ, the REF_INTERNAL contigs might differ in length from the header contigs
            rom ref_contig_name = ref_contigs_get_name (ref, ref_contig_index, NULL); // might be different that chrom_name if it matches an alt_name
        
            ASSINP (false, "Error: wrong reference file - different contig length: in %s \"%.*s\" has LN=%"PRIu64", but in %s \"%s\" has LN=%"PRId64, 
                    txt_name, STRf(chrom_name), *hdr_LN, ref->filename, ref_contig_name, ref_LN);
        }
    }

    return ref_contig_index;
}

// ZIP CRAM: verify that a contig name from the SAM header (in CRAM) has identical name and length to a contig in the reference
void ref_contigs_verify_same_contig_as_ref (Reference ref, rom cram_filename, STRp(chrom_name), PosType64 hdr_LN)
{            
    ASSINP0 (ref->filename, "when compressing a CRAM file that has @SQ fields, --reference or --REFERENCE must be specified");

    // get contig with exact name in reference (no alts allowed)
    WordIndex ref_contig_index = ref_contigs_get_by_name (ref, STRa(chrom_name), false, true); 
        
    // if its not found, error, as samtools might hang
    ASSINP (ref_contig_index != WORD_INDEX_NONE, "header of %s has contig \"%.*s\" (and maybe others, too), missing in %s. ",
            cram_filename, STRf(chrom_name), ref->filename);

    if (hdr_LN) { // 0 if header doesn't specify LN
        PosType64 ref_LN = B(Range, ref->ranges, ref_contig_index)->last_pos; // get from ranges because Contig.LN=0 - we don't populate it at reference creation
    
        // case: file header specifies length - it must be the same as the reference
        ASSINP (hdr_LN == ref_LN, "Error: wrong reference file - different contig length: in %s \"%.*s\" has LN=%"PRIu64", but in %s it has LN=%"PRId64, 
                cram_filename, STRf(chrom_name), hdr_LN, ref->filename, ref_LN);
    }
}

// get length of contig according to ref_contigs (loaded or stored reference) - used in piz when generating
// a file header of a translated data type eg SAM->BAM or ME23->VCF
PosType64 ref_contigs_get_contig_length (const Reference ref,
                                       WordIndex ref_index, // option 1 
                                       STRp(chrom_name),    // option 2
                                       bool enforce)
{
    if (ref_index == WORD_INDEX_NONE && chrom_name)
        ref_index = ref_contigs_get_by_name (ref, STRa(chrom_name), true, true);

    ASSERT (!enforce || (ref_index >= 0 && ref_index < ref->ranges.len), "contig=\"%.*s\" ref_index=%d not found in reference contigs (range=[0,%d])",
            chrom_name_len, chrom_name ? chrom_name : "", ref_index, (int)ref->ranges.len-1);

    if (ref_index == WORD_INDEX_NONE) return -1; // chrom_name not found in ref_contigs
    
    // get info as it appears in reference
    return B(Range, ref->ranges, ref_index)->last_pos;
}

// get contig by ref_index, by binary searching for it
static const Contig *ref_contigs_get_contig_do (const Reference ref, WordIndex ref_index, int32_t start_i, int32_t end_i)
{
    int32_t mid_i = (start_i + end_i) / 2;
    Contig *mid_rc = B(Contig, ref->ctgs.contigs, mid_i);
    if (mid_rc->ref_index == ref_index) return mid_rc; // found
    if (start_i >= end_i) return NULL; // not found

    // continue searching
    if (ref_index < mid_rc->ref_index)
        return ref_contigs_get_contig_do (ref, ref_index, start_i, mid_i-1);
    else    
        return ref_contigs_get_contig_do (ref, ref_index, mid_i+1, end_i);
}

const Contig *ref_contigs_get_contig_by_ref_index (const Reference ref, WordIndex ref_index, FailType soft_fail)
{
    ASSERTISALLOCED (ref->ctgs.contigs);

    const Contig *rc = ref_contigs_get_contig_do (ref, ref_index, 0, ref->ctgs.contigs.len-1);

    ASSERT (rc || soft_fail, "cannot find contig for ref_index=%d", ref_index);

    return rc;
}

uint32_t ref_contigs_get_num_contigs (Reference ref)
{
    return flag.make_reference ? ZCTX(CHROM)->nodes.len32 : ref->ctgs.contigs.len32;
}

PosType64 ref_contigs_get_genome_nbases (const Reference ref)
{
    if (flag.make_reference)
        return B(Range, ref->ranges, ref->ranges.len-1)->gpos +
               B(Range, ref->ranges, ref->ranges.len-1)->last_pos;

    if (!ref->ctgs.contigs.len) return 0;

    // note: contigs are sorted by chrom and then pos within chrom, NOT by gpos! 
    // (up to 13.0.20, chroms might have been added out of order during reference creation)
    const Contig *rc_with_largest_gpos = NULL;
    for_buf (const Contig, rc, ref->ctgs.contigs)  {
        if (!rc->min_pos && !rc->max_pos) continue; // not a real contig, eg "*" in SAM
        if (!rc_with_largest_gpos || rc->gpos > rc_with_largest_gpos->gpos) rc_with_largest_gpos = rc;
    }
    
    ASSERT (rc_with_largest_gpos, "There are %u elements in the contigs array, but none are actual contigs", ref->ctgs.contigs.len32);

    // note: contigs with gpos beyond 4Gb will not be used by the aligner

    return rc_with_largest_gpos->gpos + (rc_with_largest_gpos->max_pos - rc_with_largest_gpos->min_pos + 1);
}

// up to 13.0.20, chroms might have been added out of order during reference creation, resulting in 
// contigs, which are sorted by chrom, being NOT necessarily sorted by gpos 
static ContigP ref_contig_search_by_gpos_v13 (const Reference ref, PosType64 gpos)
{
    // note: contigs are sorted by chrom and then pos within chrom, NOT by gpos! (chroms might have been added out of order during reference creation)
    for_buf (Contig, rc, ref->ctgs.contigs) 
        if (gpos >= rc->gpos && gpos <= rc->gpos + (rc->max_pos - rc->min_pos)) 
            return rc;

    return NULL; // not found
}

static ContigP ref_contig_search_by_gpos (const Reference ref, PosType64 gpos, 
                                          WordIndex first_ctg_i, WordIndex last_ctg_i,
                                          bool next_contig_if_in_gap)
{
    if (first_ctg_i > last_ctg_i) 
        return next_contig_if_in_gap ? B(Contig, ref->ctgs.contigs, first_ctg_i)
                                     : NULL; // gpos is after all contigs, or in the gaps between contigs

    WordIndex mid_ctg_i = (first_ctg_i + last_ctg_i) / 2;
    ContigP rc = B(Contig, ref->ctgs.contigs, mid_ctg_i);

    if (gpos >= rc->gpos && gpos <= rc->gpos + (rc->max_pos - rc->min_pos))     
        return rc;

    else if (gpos < rc->gpos)
        return ref_contig_search_by_gpos (ref, gpos, first_ctg_i, mid_ctg_i - 1, next_contig_if_in_gap);

    else
        return ref_contig_search_by_gpos (ref, gpos, mid_ctg_i + 1, last_ctg_i, next_contig_if_in_gap);
}

WordIndex ref_contig_get_by_gpos (const Reference ref, PosType64 gpos,
                                  int32_t seq_len, // if non-0 succeed only if range is entirely with the contig (may be positive or negative number)
                                  PosType32 *pos,  // optional out, POS within the CHROM matching gpos
                                  bool next_contig_if_in_gap)
{
    ContigP rc = VER(14) ? ref_contig_search_by_gpos (ref, gpos, 0, ref->ctgs.contigs.len32 - 1, next_contig_if_in_gap)
                         : ref_contig_search_by_gpos_v13 (ref, gpos);
    
    if (!rc) 
        return WORD_INDEX_NONE; // gpos is after all contigs, or in the gaps between contigs
        
    WordIndex ref_index = BNUM (ref->ctgs.contigs, rc);
    ASSERT (ref_index < ref->ranges.len, "Unexpected ref_index=%d >= ref->ranges.len=%"PRIu64, ref_index, ref->ranges.len);  

    if (pos) 
        *pos = rc->min_pos + gpos - rc->gpos;

    // case: gpos itself is not in gap, but the sequence starting a gpos, of length seq_len, is not entirely on the contig
    if ((seq_len > 0 && *pos + seq_len - 1 > rc->max_pos) || // beyond end of contig
        (seq_len < 0 && *pos + seq_len + 1 < rc->max_pos))   // before the start of the contig 
        return WORD_INDEX_NONE;

    return ref_index;
}
