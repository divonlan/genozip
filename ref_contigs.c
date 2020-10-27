// ------------------------------------------------------------------
//   ref_contigs.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
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

static Buffer loaded_contigs              = EMPTY_BUFFER; // array of RefContig
static Buffer loaded_contig_dict          = EMPTY_BUFFER;
static Buffer loaded_contigs_sorted_index = EMPTY_BUFFER; // an array of uint32 of indexes into loaded_contigs - sorted by alphabetical order of the snip in contig_dict

#ifdef __linux__ 
extern int strncasecmp (const char *s1, const char *s2, size_t n); // defined in <strings.h>, but file name conflicts with "strings.h" (to do: sort this out in the Makefile)
#endif

void ref_contigs_free (void)
{
    buf_free (&loaded_contigs);          
    buf_free (&loaded_contig_dict);
    buf_free (&loaded_contigs_sorted_index);
}

uint32_t ref_contigs_num_contigs(void) { return (uint32_t)loaded_contigs.len; }

static void BGEN_ref_contigs (Buffer *contigs_buf)
{
    ARRAY (RefContig, cn, *contigs_buf);

    for (uint32_t i=0; i < contigs_buf->len; i++) 
        cn[i] = (RefContig){ .gpos        = BGEN64 (cn[i].gpos),
                             .min_pos     = BGEN64 (cn[i].min_pos),
                             .max_pos     = BGEN64 (cn[i].max_pos),
                             .char_index  = BGEN64 (cn[i].char_index),
                             .snip_len    = BGEN32 (cn[i].snip_len),
                             .chrom_index = BGEN32 (cn[i].chrom_index)     };
}

static void ref_contigs_show (const Buffer *contigs_buf, bool created)
{
    ARRAY (const RefContig, cn, *contigs_buf);
    Context *chrom_ctx = &z_file->contexts[CHROM];

    fprintf (stderr, "\nContigs as they appear in the reference%s:\n", created ? " created" : "");
    for (uint32_t i=0; i < contigs_buf->len; i++) {

        const char *chrom_name = ENT (const char, chrom_ctx->dict, cn[i].char_index);

        fprintf (stderr, "i=%u %s gpos=%"PRId64" min_pos=%"PRId64" max_pos=%"PRId64" chrom_index=%d char_index=%"PRIu64" snip_len=%u\n",
                 i, chrom_name, cn[i].gpos, cn[i].min_pos, cn[i].max_pos, cn[i].chrom_index, cn[i].char_index, cn[i].snip_len);
    }
}

// create and compress contigs - sorted by chrom_index
void ref_contigs_compress (void)
{
    static Buffer created_contigs = EMPTY_BUFFER;  

    if (!buf_is_allocated (&z_file->contexts[CHROM].mtf)) return; // no contigs

    // the number of contigs is at most the number of chroms - but could be less if some chroms have no sequence
    buf_alloc (evb, &created_contigs, sizeof (RefContig) * z_file->contexts[CHROM].mtf.len, 1, "created_contigs", 0);

    RefContig *last = NULL;

    for (uint32_t range_i=0; range_i < ranges.len; range_i++) {
        Range *r = ENT (Range, ranges, range_i);

        // in case of REF_INTERNAL, chrom_word_index might still be WORD_INDEX_NONE. We get it now from the z_file data
        if (r->chrom == WORD_INDEX_NONE)
            r->chrom = ref_contigs_get_word_index (r->chrom_name, r->chrom_name_len, WI_REF_CONTIG, false);

        // first range of a contig
        if (!last || r->chrom != last->chrom_index) {

            // if gpos of the first range in this contig has been increased due to removing flaking regions
            // and is no longer 64-aligned, we start the contig a little earlier to be 64-aligned
            // (note: this is guaranteed to be within the original contig before compacting, because the original
            // contig had a 64-aligned gpos)
            PosType delta = (ranges.param == RT_LOADED) ? r->gpos % 64 : 0;

            // in case of RT_DENOVO, we assign 64-aligned gpos (in case of RT_LOADED - gposes are loaded)
            if (ranges.param == RT_DENOVO)
                r->gpos = r->chrom ? ROUNDUP64 ((r-1)->gpos + (r-1)->last_pos - (r-1)->first_pos + 1) : 0;

            MtfNode *chrom_node = ENT (MtfNode, z_file->contexts[CHROM].mtf, r->chrom);

            NEXTENT (RefContig, created_contigs) = (RefContig){
                .gpos        = r->gpos - delta, 
                .min_pos     = r->first_pos - delta,
                .max_pos     = r->last_pos,
                .chrom_index = r->chrom, 
                .char_index  = chrom_node->char_index, 
                .snip_len    = chrom_node->snip_len 
            };

            last = LASTENT (RefContig, created_contigs);
        }
        else {
            // set gpos of the next range relative to this chrom (note: ranges might not be contiguous)
            if (ranges.param == RT_DENOVO) 
                r->gpos = (r-1)->gpos + (r->first_pos - (r-1)->first_pos); // note: only the first range in the chrom needs to be 64-aligned
            
            last->max_pos = r->last_pos; // update
        }
    }
    
    if (flag_show_ref_contigs) ref_contigs_show (&created_contigs, true);

    BGEN_ref_contigs (&created_contigs);

    created_contigs.len *= sizeof (RefContig);
    zfile_compress_section_data_codec (evb, SEC_REF_CONTIGS, &created_contigs, 0,0, CODEC_LZMA); // compresses better with LZMA than BZLIB
    
    buf_free (&created_contigs);
}

static int ref_contigs_sort_contigs_alphabetically (const void *a, const void *b)
{
    uint32_t index_a = *(uint32_t *)a;
    uint32_t index_b = *(uint32_t *)b;

    RefContig *contig_a = ENT (RefContig, loaded_contigs, index_a);
    RefContig *contig_b = ENT (RefContig, loaded_contigs, index_b);
    
    return strcmp (ENT (char, loaded_contig_dict, contig_a->char_index),
                   ENT (char, loaded_contig_dict, contig_b->char_index));
}

static void ref_contigs_create_sorted_index (void)
{
    if (buf_is_allocated (&loaded_contigs_sorted_index)) return; // already done

    // contig_words_sorted_index - an array of uint32 of indexes into contig_words - sorted by alphabetical order of the snip in contig_dict
    buf_alloc (evb, &loaded_contigs_sorted_index, sizeof(uint32_t) * loaded_contigs.len, 1, "loaded_contigs_sorted_index", 0);
    for (uint32_t i=0; i < loaded_contigs.len; i++)
        NEXTENT (uint32_t, loaded_contigs_sorted_index) = i;

    qsort (loaded_contigs_sorted_index.data, loaded_contigs_sorted_index.len, sizeof(uint32_t), ref_contigs_sort_contigs_alphabetically);
}

static int ref_contigs_sort_chroms_alphabetically (const void *a, const void *b)
{
    Context *ctx = &z_file->contexts[CHROM];

    uint32_t index_a = *(uint32_t *)a;
    uint32_t index_b = *(uint32_t *)b;

    MtfWord *word_a = ENT (MtfWord, ctx->word_list, index_a);
    MtfWord *word_b = ENT (MtfWord, ctx->word_list, index_b);
    
    return strcmp (ENT (char, ctx->dict, word_a->char_index),
                   ENT (char, ctx->dict, word_b->char_index));
}

// called by I/O thread from piz_read_global_area of user file (not reference file)
void ref_contigs_sort_chroms (void)
{
    uint32_t num_chroms = z_file->contexts[CHROM].word_list.len;

    // z_file->chroms_sorted_index - an array of uint32 of indexes into z_file->contexts[CHROM].word_list - sorted by alphabetical order of the snip in z_file->contexts[CHROM].dict
    buf_alloc (evb, &z_file->chroms_sorted_index, sizeof(uint32_t) * num_chroms, 1, "z_file->chroms_sorted_index", 0);
    for (uint32_t i=0; i < num_chroms; i++)
        NEXTENT (uint32_t, z_file->chroms_sorted_index) = i;

    qsort (z_file->chroms_sorted_index.data, num_chroms, sizeof(uint32_t), ref_contigs_sort_chroms_alphabetically);
}

// read and uncompress a contigs section (when pizzing the reference file or pizzing a data file with a stored reference)
void ref_contigs_load_contigs (void)
{
    SectionListEntry *sl = sections_get_first_section_of_type (SEC_REF_CONTIGS, false);
    zfile_read_section (z_file, evb, 0, &evb->z_data, "z_data", sizeof (SectionHeader), SEC_REF_CONTIGS, sl);

    if (flag_show_headers && exe_type == EXE_GENOCAT) goto done;

    zfile_uncompress_section (evb, evb->z_data.data, &loaded_contigs, "loaded_contigs", 0, SEC_REF_CONTIGS);

    loaded_contigs.len /= sizeof (RefContig);
    BGEN_ref_contigs (&loaded_contigs);

    buf_copy (evb, &loaded_contig_dict, &z_file->contexts[CHROM].dict, 1, 0, 0, "contig_dict", 0);

    ref_contigs_create_sorted_index();

    if (flag_show_ref_contigs) {
        ref_contigs_show (&loaded_contigs, false);
        if (exe_type == EXE_GENOCAT) exit_ok;  // in genocat this, not the data
    }

done:
    buf_free (&evb->z_data);
}

// called by mtf_copy_reference_contig_to_chrom_ctx when initializing ZIP for a new file using a pre-loaded external reference
void ref_contigs_get (ConstBufferP *out_contig_dict, ConstBufferP *out_contigs)
{
    if (out_contig_dict) *out_contig_dict = &loaded_contig_dict;
    if (out_contigs) *out_contigs = &loaded_contigs;
}

uint32_t ref_num_loaded_contigs (void)
{
    return (uint32_t)loaded_contigs.len;
}

// binary search for this chrom in loaded_contigs_sorted_index. we count on gcc tail recursion optimization to keep this fast.
static WordIndex ref_contigs_get_word_index_do (const char *chrom_name, unsigned chrom_name_len, GetWordIndexType wi_type, 
                                                WordIndex first_sorted_index, WordIndex last_sorted_index)
{
    if (first_sorted_index > last_sorted_index) return WORD_INDEX_NONE; // not found

    WordIndex mid_sorted_index = (first_sorted_index + last_sorted_index) / 2;
    
    const char *snip;
    uint32_t snip_len;
    WordIndex word_index;
    if (wi_type == WI_REF_CONTIG) {
        word_index = *ENT (WordIndex, loaded_contigs_sorted_index, mid_sorted_index);
        RefContig *mid_word = ENT (RefContig, loaded_contigs, word_index);
        snip = ENT (const char, loaded_contig_dict, mid_word->char_index);
        snip_len = mid_word->snip_len;
    }
    else {
        word_index = *ENT (WordIndex, z_file->chroms_sorted_index, mid_sorted_index);
        MtfWord *mid_word = ENT (MtfWord, z_file->contexts[CHROM].word_list, word_index);
        snip = ENT (const char, z_file->contexts[CHROM].dict, mid_word->char_index);
        snip_len = mid_word->snip_len;
    }
 
    int cmp = strncmp (snip, chrom_name, chrom_name_len);
    if (!cmp && snip_len != chrom_name_len) // identical prefix but different length
        cmp = snip_len - chrom_name_len;

    if (cmp < 0) return ref_contigs_get_word_index_do (chrom_name, chrom_name_len, wi_type, mid_sorted_index+1, last_sorted_index);
    if (cmp > 0) return ref_contigs_get_word_index_do (chrom_name, chrom_name_len, wi_type, first_sorted_index, mid_sorted_index-1);

    return word_index;
}             

// binary search for this chrom in loaded_contigs_sorted_index. we count on gcc tail recursion optimization to keep this fast.
WordIndex ref_contigs_get_word_index (const char *chrom_name, unsigned chrom_name_len, GetWordIndexType wi_type, bool soft_fail)
{
    WordIndex word_index = ref_contigs_get_word_index_do (chrom_name, chrom_name_len, wi_type, 0, 
                                                          wi_type == WI_REF_CONTIG ? loaded_contigs_sorted_index.len-1 : z_file->chroms_sorted_index.len-1);

    ASSERT (soft_fail || word_index != WORD_INDEX_NONE, "Error: contig \"%.*s\" is observed in %s but is not found in the reference %s",
            chrom_name_len, chrom_name, txt_name, ref_filename);

    return word_index;
}

const char *ref_contigs_get_chrom_snip (WordIndex chrom_index, const char **snip, uint32_t *snip_len)
{
    const RefContig *contig = ref_contigs_get_contig (chrom_index, false);
    const char *my_snip = ENT (const char, loaded_contig_dict, contig->char_index);
    
    if (snip) *snip = my_snip;
    if (snip_len) *snip_len = contig->snip_len;

    return my_snip; 
}

void ref_contigs_generate_data_if_denovo (void)
{
    // copy data from the reference FASTA's CONTIG context, so it survives after we finish reading the reference and close z_file
    Context *chrom_ctx = &z_file->contexts[CHROM];

    ASSERT (flag_reference == REF_INTERNAL || (buf_is_allocated (&chrom_ctx->dict) && buf_is_allocated (&chrom_ctx->word_list)),
            "Error: cannot use %s as a reference as it is missing a CONTIG dictionary", z_name);

    buf_copy (evb, &loaded_contig_dict, &chrom_ctx->dict, 1, 0, 0, "contig_dict", 0);
    
    // we copy from the z_file context after we completed segging a file
    // note: we can't rely on chrom_ctx->mtf for the correct order as vb_i=1 resorted. rather, by mimicking the
    // word_list generation as done in PIZ, we guarantee that we will get the same chrom_index
    // in case of multiple bound files, we re-do this in every file in case of additional chroms (not super effecient, but good enough because the context is small)
    loaded_contigs.len = chrom_ctx->mtf.len;
    buf_alloc (evb, &loaded_contigs, loaded_contigs.len * sizeof (RefContig), 1, "loaded_contigs", 0);

    // similar logic to mtf_integrate_dictionary_fragment
    char *start = loaded_contig_dict.data;
    for (uint32_t snip_i=0; snip_i < loaded_contigs.len; snip_i++) {

        RefContig *contig = ENT (RefContig, loaded_contigs, snip_i);

        char *c=start; while (*c != SNIP_SEP) c++;
        contig->snip_len   = c - start;
        contig->char_index = start - loaded_contig_dict.data;

        start = c+1; // skip over the SNIP_SEP
    }

    // contig_words_sorted_index - an array of uint32 of indexes into contig_words - sorted by alphabetical order of the snip in contig_dict
    ref_contigs_create_sorted_index();
}

// ZIP SAM: verify that we have the specified chrom (name & last_pos) loaded from the reference. called by sam_inspect_txt_header
// returns the chrom_index of the chrom
WordIndex ref_contigs_verify_identical_chrom (const char *chrom_name, unsigned chrom_name_len, PosType last_pos,
                                              WordIndex must_be_chrom_index) // provided if chrom must appear with this index or WORD_INDEX_NONE                                    
{
    WordIndex chrom_index = ref_contigs_get_word_index (chrom_name, chrom_name_len, WI_REF_CONTIG, true);

    // if not found, try common alternative names
    if (chrom_index == WORD_INDEX_NONE)
        chrom_index = ref_alt_chroms_zip_get_alt_index (chrom_name, chrom_name_len, WI_REF_CONTIG, WORD_INDEX_NONE);

    ASSERT (chrom_index != WORD_INDEX_NONE || must_be_chrom_index == WORD_INDEX_NONE, "Error: contig %s appears in file header, but is absent from the reference file %s",
            chrom_name, ref_filename);

    ASSERT (chrom_index == must_be_chrom_index || must_be_chrom_index == WORD_INDEX_NONE, "Error: contig %s appears in the file header at index %d, but appears at index %d in the reference file %s",
            chrom_name, must_be_chrom_index, chrom_index, ref_filename);

    // if its not found, we ignore it. sequences that have this chromosome will just be non-ref
    RETURNW (chrom_index != WORD_INDEX_NONE, chrom_index, "Warning: header of %s lists contig '%.*s', which is absent from %s. This might result in sub-optimal compression (but no harm otherwise).",
             txt_file->basename, chrom_name_len, chrom_name, ref_filename);

    // get info as it appears in reference
    PosType ref_last_pos = ENT (Range, ranges, chrom_index)->last_pos;
    RefContig *contig = ENT (RefContig, loaded_contigs, chrom_index);
    const char *ref_chrom_name = ENT (const char, loaded_contig_dict, contig->char_index); // might be different that chrom_name if we used an alt_name

    ASSERT (last_pos == ref_last_pos, "Error: wrong reference file: file %s has a @SQ header field specifying contig '%.*s' with LN=%"PRId64", however in reference %s the last position of '%s' is %"PRId64,
            txt_name, chrom_name_len, chrom_name, last_pos, ref_filename, ref_chrom_name, ref_last_pos);

    return chrom_index;
}

// get contig by chrom_index, by binary searching for it
static const RefContig *ref_contigs_get_contig_do (WordIndex chrom_index, int32_t start_i, int32_t end_i)
{
    int32_t mid_i = (start_i + end_i) / 2;
    RefContig *mid_rc = ENT (RefContig, loaded_contigs, mid_i);
    if (mid_rc->chrom_index == chrom_index) return mid_rc; // found
    if (start_i >= end_i) return NULL; // not found

    // continue searching
    if (chrom_index < mid_rc->chrom_index)
        return ref_contigs_get_contig_do (chrom_index, start_i, mid_i-1);
    else    
        return ref_contigs_get_contig_do (chrom_index, mid_i+1, end_i);
}

const RefContig *ref_contigs_get_contig (WordIndex chrom_index, bool soft_fail)
{
    const RefContig *rc = ref_contigs_get_contig_do (chrom_index, 0, loaded_contigs.len-1);

    ASSERT (rc || soft_fail, "Error in ref_contigs_get_contig: cannot find contig for chrom_index=%d", chrom_index);

    return rc;
}

PosType ref_contigs_get_genome_size (void)
{
    if (!loaded_contigs.len) return 0;

    // note: contigs are sorted by chrom and then pos within chrom, NOT by gpos! (chroms might have been added out of order during reference creation)
    ARRAY (RefContig, rc, loaded_contigs);
    RefContig *rc_with_largest_gpos = &rc[0];
    for (uint64_t i=1; i < loaded_contigs.len; i++) 
        if (rc[i].gpos > rc_with_largest_gpos->gpos) rc_with_largest_gpos = &rc[i];

    return rc_with_largest_gpos->gpos + (rc_with_largest_gpos->max_pos - rc_with_largest_gpos->min_pos + 1);
}

WordIndex ref_contigs_get_by_accession_number (const char *ac, unsigned ac_len)
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
    for (WordIndex alt_chrom=0 ; alt_chrom < loaded_contigs.len; alt_chrom++) {
        const char *alt_chrom_name = ENT (const char, loaded_contig_dict, ENT (RefContig, loaded_contigs, alt_chrom)->char_index);
        const char *substr = strstr (alt_chrom_name, numeric);
        if (substr && (substr - alt_chrom_name >= letter_len) && !strncasecmp (substr-letter_len, ac, (size_t)(numeric_len + letter_len))) 
            return alt_chrom;
    }

    return NODE_INDEX_NONE; // not found
}
