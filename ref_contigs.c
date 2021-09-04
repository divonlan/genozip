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
#include "txtheader.h"
#include "fasta.h"

#ifdef __linux__ 
extern int strncasecmp (const char *s1, const char *s2, size_t n); // defined in <strings.h>, but file name conflicts with "strings.h" (to do: sort this out in the Makefile)
#endif

static void BGEN_ref_contigs (Buffer *contigs_buf)
{
    for (uint32_t i=0; i < contigs_buf->len; i++) {
        RefContig *cn   = ENT (RefContig, *contigs_buf, i);
        cn->gpos        = BGEN64 (cn->gpos);
        cn->min_pos     = BGEN64 (cn->min_pos);
        cn->max_pos     = BGEN64 (cn->max_pos);
        cn->char_index  = BGEN64 (cn->char_index);
        cn->snip_len    = BGEN32 (cn->snip_len);
        cn->chrom_index = BGEN32 (cn->chrom_index);
    };
}

static void ref_contigs_show (const Buffer *contigs_buf, bool created)
{
    ARRAY (const RefContig, cn, *contigs_buf);

    iprintf ("\nContigs as they appear in the reference%s:\n", created ? " created (note: contig names are as they appear in the txt data, not the reference)" : "");
    for (uint32_t i=0; i < cn_len; i++) {

        const char *chrom_name = ENT (const char, ZCTX(CHROM)->dict, cn[i].char_index);

        if (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE || flag.reference == REF_LIFTOVER || 
            z_file->data_type == DT_REF) 
            iprintf ("\"%s\" length=%"PRId64" i=%u gpos=%"PRId64" \"%s\"\n",
                     chrom_name,  cn[i].max_pos, cn[i].chrom_index, cn[i].gpos, cn[i].metadata);

        else if (cn[i].snip_len)
            iprintf ("i=%u '%s' gpos=%"PRId64" min_pos=%"PRId64" max_pos=%"PRId64" chrom_index=%d char_index=%"PRIu64" snip_len=%u\n",
                     i, chrom_name, cn[i].gpos, cn[i].min_pos, cn[i].max_pos, cn[i].chrom_index, cn[i].char_index, cn[i].snip_len);
        else
            iprintf ("i=%u chrom_index=%d (unused - not present in txt data)\n", i, cn[i].chrom_index);
    }
}

// create and compress contigs - sorted by chrom_index
void ref_contigs_compress (Reference ref)
{
    START_TIMER;

    static Buffer created_contigs = EMPTY_BUFFER;  

    if (!buf_is_alloc (&ZCTX(CHROM)->nodes)) return; // no contigs

    // the number of contigs is at most the number of chroms - but could be less if some chroms have no sequence
    // in SAM with header and -E, we don't copy the contigs to CHROM, so we take the length from loaded_contigs.
    buf_alloc (evb, &created_contigs, 0, MAX_(ZCTX(CHROM)->nodes.len, ref->loaded_contigs.len), RefContig, 1, "created_contigs");

    RefContig *last = NULL;

    ConstBufferP contig_metadata = ref_make_get_contig_metadata(); // empty unless this is --make-reference

    for (uint32_t range_i=0; range_i < ref->ranges.len; range_i++) {
        Range *r = ENT (Range, ref->ranges, range_i);

        // in case of RT_DENOVO, chrom_word_index might still be WORD_INDEX_NONE. We get it now from the z_file data
        if (r->chrom == WORD_INDEX_NONE)
            r->chrom = ref_contigs_get_by_name (ref, r->chrom_name, r->chrom_name_len, false, false);

        // first range of a contig
        if (!last || r->chrom != last->chrom_index) {

            // if gpos of the first range in this contig has been increased due to removing flanking regions
            // and is no longer 64-aligned, we start the contig a little earlier to be 64-aligned
            // (note: this is guaranteed to be within the original contig before compacting, because the original
            // contig had a 64-aligned gpos)
            PosType delta = (ranges_type(ref) == RT_LOADED || ranges_type(ref) == RT_CACHED) ? r->gpos % 64 : 0;

            // in case of RT_DENOVO, we assign 64-aligned gpos (in case of RT_LOADED - gposes are loaded)
            if (ranges_type(ref) == RT_DENOVO)
                r->gpos = range_i ? ROUNDUP64 ((r-1)->gpos + (r-1)->last_pos - (r-1)->first_pos + 1) : 0;

            WordIndex txt_chrom = WORD_INDEX_NONE;

            // case: we have header contigs. In this case, the contigs in reference (including in ranges) have a different index
            // as their index is relative to the reference
            const Buffer *header_contigs = txtheader_get_contigs();
            if (header_contigs) {
                // search for the reference chrom index r->chrom in header_contigs to find the equivalent chrom (same or alt name) in txt
                for (uint32_t i=0; i < header_contigs->len; i++) 
                    if (ENT (RefContig, *header_contigs, i)->chrom_index == r->chrom) {
                        txt_chrom = i;
                        break; // found
                    }
                // note: this is O(num_chrom^2). but tests on my PC with 3000 contigs ref x 300 contig file result in 3 millisec. we could build a sorted index and binary search, but not worth it
            }
            
            // case: absent header contigs, reference chroms are copied to contexts[CHROM] and hence the indices are the same
            else
                txt_chrom = r->chrom;
                
            CtxNode *chrom_node = (txt_chrom != WORD_INDEX_NONE) ? ENT (CtxNode, ZCTX(CHROM)->nodes, txt_chrom) : NULL;

            NEXTENT (RefContig, created_contigs) = (RefContig){
                .gpos        = r->gpos - delta, 
                .min_pos     = r->first_pos - delta,
                .max_pos     = r->last_pos,
                .chrom_index = r->chrom,  // reference index (except if REF_INTERNAL - index into ZCTX(CHROM))
                .char_index  = chrom_node ? chrom_node->char_index : 0, 
                .snip_len    = chrom_node ? chrom_node->snip_len   : 0
            };

            last = LASTENT (RefContig, created_contigs);

            if (flag.make_reference)
                memcpy (last->metadata, ENT (ContigMetadata, *contig_metadata, r->chrom), REFCONTIG_MD_LEN);
        }
        else {
            // set gpos of the next range relative to this chrom (note: ranges might not be contiguous)
            if (ranges_type(ref) == RT_DENOVO) 
                r->gpos = (r-1)->gpos + (r->first_pos - (r-1)->first_pos); // note: only the first range in the chrom needs to be 64-aligned
            
            last->max_pos = r->last_pos; // update
        }
    }
    
    if (flag.show_ref_contigs) ref_contigs_show (&created_contigs, true);

    BGEN_ref_contigs (&created_contigs);

    COPY_TIMER_VB (evb, ref_contigs_compress); // we don't count the compression itself and the disk writing

    created_contigs.len *= sizeof (RefContig);
    zfile_compress_section_data_ex (evb, SEC_REF_CONTIGS, &created_contigs, 0,0, CODEC_BSC, SECTION_FLAGS_NONE); 

    buf_free (&created_contigs);
}

static Reference sorted_ref = 0; // ref_contigs_create_sorted_index is always called from the main thread, so no thread safety issues
static int ref_contigs_sort_contigs_alphabetically (const void *a, const void *b)
{
    uint32_t index_a = *(uint32_t *)a;
    uint32_t index_b = *(uint32_t *)b;

    RefContig *contig_a = ENT (RefContig, sorted_ref->loaded_contigs, index_a);
    RefContig *contig_b = ENT (RefContig, sorted_ref->loaded_contigs, index_b);
    
    return strcmp (ENT (char, sorted_ref->loaded_contigs_dict, contig_a->char_index),
                   ENT (char, sorted_ref->loaded_contigs_dict, contig_b->char_index));
}

static int ref_contigs_sort_contigs_by_length (const void *a, const void *b)
{
    uint32_t index_a = *(uint32_t *)a;
    uint32_t index_b = *(uint32_t *)b;

    RefContig *contig_a = ENT (RefContig, sorted_ref->loaded_contigs, index_a);
    RefContig *contig_b = ENT (RefContig, sorted_ref->loaded_contigs, index_b);
    
    // note: explicit comparison and not substraction, as LN is 64bit and return value is 32.
    return contig_a->max_pos > contig_b->max_pos ? 1
        :  contig_a->max_pos < contig_b->max_pos ? -1
        :                                          0;
}

// create to sorted indices for reference contigs - one by name, one by LN
static void ref_contigs_create_sorted_index (Reference ref)
{
    if (buf_is_alloc (&ref->loaded_contigs_by_name)) return; // already done

    buf_alloc (evb, &ref->loaded_contigs_by_name, 0, ref->loaded_contigs.len, uint32_t, 1, "loaded_contigs_by_name");
    buf_alloc (evb, &ref->loaded_contigs_by_LN,   0, ref->loaded_contigs.len, uint32_t, 1, "loaded_contigs_by_LN");

    for (uint32_t i=0; i < ref->loaded_contigs.len; i++)
        NEXTENT (uint32_t, ref->loaded_contigs_by_name) = NEXTENT (uint32_t, ref->loaded_contigs_by_LN) = i;

    sorted_ref = ref;

    qsort (ref->loaded_contigs_by_name.data, ref->loaded_contigs_by_name.len, sizeof(uint32_t), ref_contigs_sort_contigs_alphabetically);
    qsort (ref->loaded_contigs_by_LN.data, ref->loaded_contigs_by_LN.len, sizeof(uint32_t), ref_contigs_sort_contigs_by_length);
}

static int ref_contigs_sort_chroms_alphabetically (const void *a, const void *b)
{
    Context *ctx = ZCTX (CHROM);

    uint32_t index_a = *(uint32_t *)a;
    uint32_t index_b = *(uint32_t *)b;

    CtxWord *word_a = ENT (CtxWord, ctx->word_list, index_a);
    CtxWord *word_b = ENT (CtxWord, ctx->word_list, index_b);
    
    return strcmp (ENT (char, ctx->dict, word_a->char_index),
                   ENT (char, ctx->dict, word_b->char_index));
}

// called by main thread from piz_read_global_area of user file (not reference file)
void ref_contigs_sort_chroms (void)
{
    uint32_t num_chroms = ZCTX(CHROM)->word_list.len;

    // z_file->chroms_sorted_index - an array of uint32 of indexes into ZCTX(CHROM)->word_list - sorted by alphabetical order of the snip in ZCTX(CHROM)->dict
    buf_alloc (evb, &z_file->chroms_sorted_index, 0, num_chroms, uint32_t, 1, "z_file->chroms_sorted_index");
    for (uint32_t i=0; i < num_chroms; i++)
        NEXTENT (uint32_t, z_file->chroms_sorted_index) = i;

    qsort (z_file->chroms_sorted_index.data, num_chroms, sizeof(uint32_t), ref_contigs_sort_chroms_alphabetically);
}

// read and uncompress a contigs section (when pizzing the reference file or pizzing a data file with a stored reference)
void ref_contigs_load_contigs (Reference ref)
{
    Section sl = sections_last_sec (SEC_REF_CONTIGS, true);
    if (!sl) return; // section doesn't exist

    zfile_get_global_section (SectionHeader, SEC_REF_CONTIGS, sl, &ref->loaded_contigs, "loaded_contigs");

    ref->loaded_contigs.len /= sizeof (RefContig);
    BGEN_ref_contigs (&ref->loaded_contigs);

    ASSERT0 (ZCTX(CHROM)->dict.len, "CHROM dictionary is empty");

    buf_copy (evb, &ref->loaded_contigs_dict, &ZCTX(CHROM)->dict, char, 0, 0, "loaded_contigs_dict");

    ref_contigs_create_sorted_index(ref);

    if (flag.show_ref_contigs) {
        ref_contigs_show (&ref->loaded_contigs, false);
        if (exe_type == EXE_GENOCAT) exit_ok();  // in genocat this, not the data
    }
}

// called by ctx_build_zf_ctx_from_contigs when initializing ZIP for a new file using a pre-loaded external reference
void ref_contigs_get (Reference ref, ConstBufferP *out_contig_dict, ConstBufferP *out_contigs)
{
    if (out_contig_dict) *out_contig_dict = &ref->loaded_contigs_dict;
    if (out_contigs) *out_contigs = &ref->loaded_contigs;
}

uint32_t ref_num_loaded_contigs (Reference ref)
{
    return (uint32_t)ref->loaded_contigs.len;
}

// binary search for this chrom in loaded_contigs_by_name. we count on gcc tail recursion optimization to keep this fast.
static WordIndex ref_contigs_get_by_name_do (Reference ref, const char *chrom_name, unsigned chrom_name_len, 
                                             WordIndex first_sorted_index, WordIndex last_sorted_index)
{
    if (first_sorted_index > last_sorted_index) return WORD_INDEX_NONE; // not found

    WordIndex mid_sorted_index = (first_sorted_index + last_sorted_index) / 2;
    
    const char *snip;
    uint32_t snip_len;
    WordIndex word_index;
    if (ref) {
        word_index = *ENT (WordIndex, ref->loaded_contigs_by_name, mid_sorted_index);
        RefContig *mid_word = ENT (RefContig, ref->loaded_contigs, word_index);
        snip = ENT (const char, ref->loaded_contigs_dict, mid_word->char_index);
        snip_len = mid_word->snip_len;
    }
    else {
        word_index = *ENT (WordIndex, z_file->chroms_sorted_index, mid_sorted_index);
        CtxWord *mid_word = ENT (CtxWord, ZCTX(CHROM)->word_list, word_index);
        snip = ENT (const char, ZCTX(CHROM)->dict, mid_word->char_index);
        snip_len = mid_word->snip_len;
    }
 
    int cmp = strncmp (snip, chrom_name, chrom_name_len);
    if (!cmp && snip_len != chrom_name_len) // identical prefix but different length
        cmp = snip_len - chrom_name_len;

    if (cmp < 0) return ref_contigs_get_by_name_do (ref, chrom_name, chrom_name_len, mid_sorted_index+1, last_sorted_index);
    if (cmp > 0) return ref_contigs_get_by_name_do (ref, chrom_name, chrom_name_len, first_sorted_index, mid_sorted_index-1);

    return word_index;
}             

// binary search for this chrom in loaded_contigs_by_name. we count on gcc tail recursion optimization to keep this fast.
// looks in Reference contigs if ref is provided, or in CHROM if ref=NULL
WordIndex ref_contigs_get_by_name (Reference ref, const char *chrom_name, unsigned chrom_name_len, bool alt_ok, bool soft_fail)
{
    WordIndex ref_index = ref_contigs_get_by_name_do (ref, chrom_name, chrom_name_len, 0, 
                                                      ref ? ref->loaded_contigs_by_name.len-1 : z_file->chroms_sorted_index.len-1);

    // case alt_ok: if not found, try common alternative names
    if (alt_ok && ref_index == WORD_INDEX_NONE) 
        ref_index = ref_alt_chroms_get_alt_index (ref, chrom_name, chrom_name_len, 0, WORD_INDEX_NONE);

    ASSINP (soft_fail || ref_index != WORD_INDEX_NONE, "Error: contig \"%.*s\" is observed in %s but is not found in the reference %s",
            chrom_name_len, chrom_name, txt_name, ref->filename);

    return ref_index;
}

static WordIndex ref_contigs_get_by_uniq_len_do (Reference ref, PosType LN, WordIndex first_sorted_index, WordIndex last_sorted_index)
{
    #define contig_LN(sorted_index) ENT (RefContig, ref->loaded_contigs, *ENT (WordIndex, ref->loaded_contigs_by_LN, (sorted_index)))->max_pos

    if (first_sorted_index > last_sorted_index) return WORD_INDEX_NONE; // not found

    WordIndex mid_sorted_index = (first_sorted_index + last_sorted_index) / 2;
    PosType mid_LN = contig_LN (mid_sorted_index);
 
    if (mid_LN < LN) return ref_contigs_get_by_uniq_len_do (ref, LN, mid_sorted_index+1, last_sorted_index);
    if (mid_LN > LN) return ref_contigs_get_by_uniq_len_do (ref, LN, first_sorted_index, mid_sorted_index-1);

    // verify that this is the only contig with this length
    if ((mid_sorted_index > 0 && contig_LN (mid_sorted_index-1) == LN) ||
        (mid_sorted_index+1 < ref->loaded_contigs_by_LN.len && contig_LN (mid_sorted_index+1) == LN))
        return WORD_INDEX_NONE; // not unique

    return *ENT (WordIndex, ref->loaded_contigs_by_LN, mid_sorted_index);
    #undef contig_LN
}             

// binary search for this chrom in loaded_contigs_by_LN. 
WordIndex ref_contigs_get_by_uniq_len (Reference ref, PosType LN)
{
    return ref_contigs_get_by_uniq_len_do (ref, LN, 0, ref->loaded_contigs_by_LN.len-1);
}

const char *ref_contigs_get_chrom_snip (Reference ref, WordIndex chrom_index, const char **snip, uint32_t *snip_len)
{
    const RefContig *contig = ref_contigs_get_contig (ref, chrom_index, false);
    const char *my_snip = ENT (const char, ref->loaded_contigs_dict, contig->char_index);
    
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

    buf_copy (evb, &ref->loaded_contigs_dict, &chrom_ctx->dict, char, 0, 0, "contig_dict");
    
    // we copy from the z_file context after we completed segging a file
    // note: we can't rely on chrom_ctx->nodes for the correct order as vb_i=1 resorted. rather, by mimicking the
    // word_list generation as done in PIZ, we guarantee that we will get the same chrom_index
    // in case of multiple bound files, we re-do this in every file in case of additional chroms (not super effecient, but good enough because the context is small)
    ref->loaded_contigs.len = chrom_ctx->nodes.len;
    buf_alloc (evb, &ref->loaded_contigs, 0, ref->loaded_contigs.len, RefContig, 1, "loaded_contigs");

    // similar logic to ctx_dict_build_word_lists
    char *start = ref->loaded_contigs_dict.data;
    for (uint32_t snip_i=0; snip_i < ref->loaded_contigs.len; snip_i++) {

        RefContig *contig = ENT (RefContig, ref->loaded_contigs, snip_i);

        char *c=start; while (*c) c++;
        contig->snip_len   = c - start;
        contig->char_index = start - ref->loaded_contigs_dict.data;

        start = c+1; // skip over the \0
    }

    // contig_words_sorted_index - an array of uint32 of indexes into contig_words - sorted by alphabetical order of the snip in contig_dict
    ref_contigs_create_sorted_index(ref);
}

// ZIP SAM/BAM and VCF: verify that we have the specified chrom (name & last_pos) loaded from the reference at the same index. 
// called by txtheader_add_contig. returns true if contig is in reference
WordIndex ref_contigs_ref_chrom_from_header_chrom (Reference ref, const char *chrom_name, unsigned chrom_name_len, 
                                                   PosType *last_pos, // if 0, set from reference, otherwise verify
                                                   WordIndex header_chrom)
{               
    WordIndex ref_chrom = ref_contigs_get_by_name (ref, chrom_name, chrom_name_len, true, true);
        
    // case: found (original or alternate name), but in a different index than header - create a alt_chrom entry to map indeces
    // note: in case of an alt name, but in the correct index - no need to create an alt chrom entry as no index mapping is needed
    if (ref_chrom != WORD_INDEX_NONE && ref_chrom != header_chrom) {
        buf_alloc (evb, &z_file->alt_chrom_map, 1, 100, AltChrom, 2, "z_file->alt_chrom_map");
        NEXTENT (AltChrom, z_file->alt_chrom_map) = (AltChrom){ .txt_chrom = BGEN32 (header_chrom), 
                                                                .ref_chrom = BGEN32 (ref_chrom) };

        if (flag.show_ref_alts) 
            iprintf ("In %s header: '%.*s' (index=%d) In reference: '%s' (index=%d)\n", 
                     dt_name (txt_file->data_type), chrom_name_len, chrom_name, header_chrom, 
                     ref_contigs_get_chrom_snip (ref, ref_chrom, 0, 0), ref_chrom);
    }

    // if its not found, we ignore it. sequences that have this chromosome will just be non-ref
    if (ref_chrom == WORD_INDEX_NONE) {
        static bool shown = false;
        ASSERTW (command != ZIP || shown, "FYI: header of %s has contig '%.*s' (and maybe others, too), missing in %s. No harm.",
                 txt_file->basename, chrom_name_len, chrom_name, ref->filename);
        shown = true; // show only once 
        return WORD_INDEX_NONE;
    }

    // get info as it appears in reference
    RefContig *contig = ENT (RefContig, ref->loaded_contigs, ref_chrom);

    const char *ref_chrom_name = ENT (const char, ref->loaded_contigs_dict, contig->char_index); // might be different that chrom_name if we used an alt_name

    if (buf_is_alloc (&ref->ranges)) { // it is not allocated in --sex/coverage
        PosType ref_last_pos = ENT (Range, ref->ranges, ref_chrom)->last_pos; // get from ranges because RefContig.LN=0 - we don't populate it at reference creation
        
        // case: file header is missing length, update from reference
        if (! *last_pos) *last_pos = ref_last_pos;

        // case: file header specifies length - it must be the same as the reference
        else if (*last_pos != ref_last_pos) {
            char fmt[200]; // no unbound names in fmt
            sprintf (fmt, "Error: wrong reference file: %%s has a \"%s\", but in %%s '%%s' has LN=%%"PRId64, DTPT (header_contigs));
            ASSINP (false, fmt, txt_name, chrom_name_len, chrom_name, *last_pos, ref->filename, ref_chrom_name, ref_last_pos);
        }
    }

    return ref_chrom;
}

// get length of contig according to ref_contigs (loaded or stored reference) - used in piz when generating
// a file header of a translated data type eg SAM->BAM or ME23->VCF
PosType ref_contigs_get_contig_length (const Reference ref,
                                       WordIndex chrom_index, // option 1 
                                       const char *chrom_name, unsigned chrom_name_len, // option 2
                                       bool enforce)
{
    if (chrom_index == WORD_INDEX_NONE && chrom_name)
        chrom_index = ref_contigs_get_by_name (ref, chrom_name, chrom_name_len, false, true);

    // if not found, try common alternative names
    if ((chrom_index == WORD_INDEX_NONE || chrom_index >= ref->ranges.len) && chrom_name)
        chrom_index = ref_alt_chroms_get_alt_index (ref, chrom_name, chrom_name_len, 0, WORD_INDEX_NONE);

    ASSERT (!enforce || (chrom_index >= 0 && chrom_index < ref->ranges.len), "contig=\"%.*s\" chrom_index=%d not found in reference contigs (range=[0,%d])",
            chrom_name_len, chrom_name ? chrom_name : "", chrom_index, (int)ref->ranges.len-1);

    if (chrom_index == WORD_INDEX_NONE) return -1; // chrom_name not found in ref_contigs
    
    // get info as it appears in reference
    return ENT (Range, ref->ranges, chrom_index)->last_pos;
}

// get contig by chrom_index, by binary searching for it
static const RefContig *ref_contigs_get_contig_do (const Reference ref, WordIndex chrom_index, int32_t start_i, int32_t end_i)
{
    int32_t mid_i = (start_i + end_i) / 2;
    RefContig *mid_rc = ENT (RefContig, ref->loaded_contigs, mid_i);
    if (mid_rc->chrom_index == chrom_index) return mid_rc; // found
    if (start_i >= end_i) return NULL; // not found

    // continue searching
    if (chrom_index < mid_rc->chrom_index)
        return ref_contigs_get_contig_do (ref, chrom_index, start_i, mid_i-1);
    else    
        return ref_contigs_get_contig_do (ref, chrom_index, mid_i+1, end_i);
}

const RefContig *ref_contigs_get_contig (const Reference ref, WordIndex chrom_index, bool soft_fail)
{
    const RefContig *rc = ref_contigs_get_contig_do (ref, chrom_index, 0, ref->loaded_contigs.len-1);

    ASSERT (rc || soft_fail, "cannot find contig for chrom_index=%d", chrom_index);

    return rc;
}

// call a callback for each accessed contig. Note: callback function is the same as sam_foreach_SQ_line
void ref_contigs_iterate (const Reference ref, RefContigsIteratorCallback callback, void *callback_param)
{
    ARRAY (const RefContig, rc, ref->loaded_contigs);

    for (uint64_t i=0; i < ref->loaded_contigs.len; i++) {
        const char *ref_chrom_name = ENT (const char, ref->loaded_contigs_dict, rc[i].char_index); 
        callback (ref_chrom_name, strlen (ref_chrom_name), rc[i].max_pos, callback_param);
    }
}

PosType ref_contigs_get_genome_nbases (const Reference ref)
{
    if (!ref->loaded_contigs.len) return 0;

    // note: contigs are sorted by chrom and then pos within chrom, NOT by gpos! (chroms might have been added out of order during reference creation)
    ARRAY (const RefContig, rc, ref->loaded_contigs);
    const RefContig *rc_with_largest_gpos = &rc[0];
    for (uint64_t i=1; i < ref->loaded_contigs.len; i++) 
        if (rc[i].gpos > rc_with_largest_gpos->gpos) rc_with_largest_gpos = &rc[i];

    // note: gpos can exceed MAX_GPOS if compressed with REF_INTERNAL (and will, if there are a lot of tiny contigs
    // to which we grant 1M gpos space) - this is ok because denovo doesn't use gpos, rather the POS from the SAM alighment

    ASSERT ((rc_with_largest_gpos->gpos >= 0 && rc_with_largest_gpos->gpos <= MAX_GPOS) || IS_REF_INTERNAL (z_file),
            "gpos=%"PRId64" out of range 0-%"PRId64, rc_with_largest_gpos->gpos, MAX_GPOS);

    return rc_with_largest_gpos->gpos + (rc_with_largest_gpos->max_pos - rc_with_largest_gpos->min_pos + 1);
}

WordIndex ref_contigs_get_by_accession_number (const Reference ref, const char *ac, unsigned ac_len)
{
    // "GL000207.1" -> "000207"
    char numeric[ac_len+1];
    unsigned numeric_len=0, letter_len=0;
    for (unsigned i=0; i < ac_len; i++) {
        if (ac[i] == '.') break;
        if (!IS_DIGIT (ac[i])) {
            letter_len++;
            continue;
        }
        numeric[numeric_len++] = ac[i];
    }
    numeric[numeric_len] = 0;

    // search by numeric only - easier as we don't have to worry about case. when candidate is found, compare letters
    for (WordIndex alt_chrom=0 ; alt_chrom < ref->loaded_contigs.len; alt_chrom++) {
        const char *alt_chrom_name = ENT (const char, ref->loaded_contigs_dict, ENT (RefContig, ref->loaded_contigs, alt_chrom)->char_index);
        const char *substr = strstr (alt_chrom_name, numeric);
        if (substr && (substr - alt_chrom_name >= letter_len) && !strncasecmp (substr-letter_len, ac, (size_t)(numeric_len + letter_len))) 
            return alt_chrom;
    }

    return NODE_INDEX_NONE; // not found
}

WordIndex ref_contig_get_by_gpos (const Reference ref, PosType gpos)
{
    // note: contigs are sorted by chrom and then pos within chrom, NOT by gpos! (chroms might have been added out of order during reference creation)
    for (WordIndex chrom_index=0 ; chrom_index < ref->loaded_contigs.len; chrom_index++) {
        RefContig *rc = ENT (RefContig, ref->loaded_contigs, chrom_index);
        if (gpos >= rc->gpos && gpos <= rc->gpos + (rc->max_pos - rc->min_pos)) 
            return chrom_index;
    }

    return WORD_INDEX_NONE; // not found
}

// returns txt_chrom_index if its in the reference OR an alt chrom index (eg 22->chr22) in the reference OR WORD_INDEX_NONE
WordIndex ref_contig_get_by_chrom (ConstVBlockP vb, const Reference ref, WordIndex txt_chrom_index, const char *chrom_name, unsigned chrom_name_len,
                                   PosType *max_pos) // optional out
{
    WordIndex contig_index = WORD_INDEX_NONE;

    // case: chrom is part of the reference (same index)
    const RefContig *ctg;
    if (txt_chrom_index < ref->loaded_contigs.len &&
           (!chrom_name_len ||
              ((ctg = ENT (RefContig, ref->loaded_contigs, txt_chrom_index))->snip_len == chrom_name_len &&
                !memcmp (ENT (char, ref->loaded_contigs_dict, ctg->char_index), chrom_name, chrom_name_len))))  {
        
        // note: we compare the name too, because txt_chrom isn't always the same as ref_chrom. example: VCF file in --chain, where VCF_CHROM is populated
        // from chain file, but if VCF variant is an alt name, and VCF has no header contigs, txt_chrom_index will be a new snip after the chain contigs - not
        // matching the reference 
        contig_index = txt_chrom_index; 
    }

    // case: chrom is not in the reference as is, test if it is in the reference using an alternative name (eg "22"->"chr22")
    else if (chrom_name) 
        contig_index = ref_contigs_get_by_name (ref, chrom_name, chrom_name_len, true, true);

    if (max_pos && contig_index != WORD_INDEX_NONE) 
        *max_pos = ENT (RefContig, ref->loaded_contigs, contig_index)->max_pos;

    return contig_index;
}

