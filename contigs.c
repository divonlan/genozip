// ------------------------------------------------------------------
//   contigs.c
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "contigs.h"
#include "buffer.h"
#include "strings.h"
#include "context.h"

//------------------------------------------------------
// Calculating Accession Numbers 
// see: https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/
//------------------------------------------------------

// turn contig name into a canonical AC: Supported formats:
// hg19:   "chr4_gl383528_alt"                  -> { "GL383528",     '1' }
// GRCh38: "chrUn_JTFH01001867v2_decoy"         -> { "JTFH01001867", '2' }
// hs37d5: "GL000192.1"                         -> { "GL000192",     '1' }
// GRCh37_latest_genomic.fna.gz: "NC_000002.11" -> { "NC_000002",    '11' }
static bool contig_name_to_acc_num (STRp(contig), AccessionNumber *ac) // out - canonical AC or zeroed
{
    *ac = (AccessionNumber){};

    // test for hs37d5 style
    if (contig_len >= 6 && contig_len <= ACCESSION_LEN + 2 && 
        contig[contig_len-2] == '.' && IS_DIGIT (contig[contig_len-1])) {

        *ac = (AccessionNumber){};
    
        for (unsigned i=0; i < contig_len-2; i++)
            ac->AC[i] = UPPER_CASE (contig[i]);

        ac->version = contig[contig_len-1];
        return true;
    }

    // hs37d5 style, 2 digit version
    if (contig_len >= 6 && contig_len <= ACCESSION_LEN + 3 && 
        contig[contig_len-3] == '.' && IS_DIGIT (contig[contig_len-2]) && IS_DIGIT (contig[contig_len-1])) {

        *ac = (AccessionNumber){};
    
        for (unsigned i=0; i < contig_len-3; i++)
            ac->AC[i] = UPPER_CASE (contig[i]);

        ac->version  = contig[contig_len-2];
        ac->version2 = contig[contig_len-1];

        return true;
    }

    // test for hg19 / GRCh38 style
    if (contig_len < 11 || contig[0]!='c' || contig[1]!='h' || contig[2]!='r') return false;

    str_split (contig, contig_len, 3, '_', item, true);
    if (!n_items || item_lens[1] < 6 || item_lens[1] > ACCESSION_LEN) return false; // the shortest an AC can be is 6

    if (items[1][item_lens[1]-2] == 'v') {
        ac->version = items[1][item_lens[1]-1];
        item_lens[1] -= 2;
    }
    else
        ac->version = '1';

    for (unsigned i=0; i < item_lens[1]; i++)
        ac->AC[i] = UPPER_CASE (items[1][i]);

    return true;
}

static void contigs_calculate_accession_numbers (BufferP contigs, ConstBufferP contigs_dict)
{
    ARRAY (Contig, ctg, *contigs);

    for (uint64_t i=0; i < ctg_len; i++) {

        // first check for eg "AC:GL949752.1" in the metadata
        char *s = strstr (ctg[i].metadata.str, "AC:");
        if (s) {
            s += 3;
            char ac[ACCESSION_LEN] = "";
            
            for (unsigned j=0; j < ACCESSION_LEN && (IS_DIGIT(*s) || IS_LETTER(*s)); j++, s++)
                ac[j] = UPPER_CASE (*s);

            char version  = (s - ctg[i].metadata.str <= REFCONTIG_MD_LEN-2 && *s == '.')      ? s[1] : '1';
            char version2 = (s - ctg[i].metadata.str <= REFCONTIG_MD_LEN-3 && IS_DIGIT(s[2])) ? s[2] : 0;
                 
            // overwrite string metadata with parsed meta data
            memcpy (ctg[i].metadata.parsed.ac.AC, ac, ACCESSION_LEN);
            ctg[i].metadata.parsed.ac.version = version;
            ctg[i].metadata.parsed.ac.version = version2;
        }

        // second, try to extract the AC from the contig name. 
        else {
            const char *contig_name = ENT (const char, *contigs_dict, ctg[i].char_index);
            contig_name_to_acc_num (contig_name, ctg[i].snip_len, &ctg[i].metadata.parsed.ac);
        }
    }
}

AccNumText display_acc_num (const AccessionNumber *ac)
{
    AccNumText s = {};
    if (ac->AC[0]) sprintf (s.s, "AC=%s.%c%c", ac->AC, ac->version, ac->version2);
    return s;
}

//-------------------
// Sorting
//-------------------

static ConstBufferP sorter_contigs=0, sorter_contigs_dicts=0;
#define DEF_SORTER_CONTIGS                                      \
    uint32_t index_a = *(uint32_t *)a;                          \
    uint32_t index_b = *(uint32_t *)b;                          \
    Contig *contig_a = ENT (Contig, *sorter_contigs, index_a);  \
    Contig *contig_b = ENT (Contig, *sorter_contigs, index_b);

static int contigs_alphabetical_sorter (const void *a, const void *b)
{
    DEF_SORTER_CONTIGS;
    return strcmp (ENT (char, *sorter_contigs_dicts, contig_a->char_index),
                   ENT (char, *sorter_contigs_dicts, contig_b->char_index));
}

static int contigs_accession_number_sorter (const void *a, const void *b)
{
    DEF_SORTER_CONTIGS;
    int cmp = memcmp (&contig_a->metadata.parsed.ac.AC, &contig_b->metadata.parsed.ac.AC, sizeof (AccessionNumber));
    if (!cmp) cmp = contig_a->metadata.parsed.ac.version - contig_b->metadata.parsed.ac.version;
    if (!cmp) cmp = contig_a->metadata.parsed.ac.version2 - contig_b->metadata.parsed.ac.version2;
    return cmp;
}

static int contigs_ref_index_sorter (const void *a, const void *b)
{
    DEF_SORTER_CONTIGS;
    return contig_a->ref_index - contig_b->ref_index;
}

// not thread safe: can only run from the main thread
static void contigs_sort (ConstBufferP contigs, ConstBufferP contigs_dict, BufferP index_buf, SortBy sort_by)
{
    if (!contigs->len) return; // nothing to do
    
    sorter_contigs = contigs;
    sorter_contigs_dicts = contigs_dict;    

    buf_alloc (evb, index_buf, 0, contigs->len, uint32_t, 1, contigs->name);
    index_buf->len = contigs->len;
    ARRAY (uint32_t, index, *index_buf);

    for (uint32_t i=0; i < index_len; i++)
        index[i] = i;

    switch (sort_by) {
        case SORT_BY_NAME      : qsort (index, index_len, sizeof(uint32_t), contigs_alphabetical_sorter);     break;
        case SORT_BY_REF_INDEX : qsort (index, index_len, sizeof(uint32_t), contigs_ref_index_sorter);        break;
        case SORT_BY_AC        : qsort (index, index_len, sizeof(uint32_t), contigs_accession_number_sorter); break;
        default: ABORT ("Invalid sort_by=%u", sort_by);
    }
}

//-------------------
// Finding
//-------------------

// binary search for this chrom in contigs. we count on gcc tail recursion optimization to keep this fast.
WordIndex contigs_get_by_name_do (ConstBufferP contigs, ConstBufferP contigs_dict, ConstBufferP index_buf, STRp(contig_name),
                                  WordIndex first_sorted_index, WordIndex last_sorted_index)
{
    if (first_sorted_index > last_sorted_index) return WORD_INDEX_NONE; // not found

    WordIndex mid_sorted_index = (first_sorted_index + last_sorted_index) / 2;
    
    WordIndex word_index = *ENT (WordIndex, *index_buf, mid_sorted_index);
    Contig *mid_word = ENT (Contig, *contigs, word_index);
    const char *snip = ENT (const char, *contigs_dict, mid_word->char_index);
    uint32_t snip_len = mid_word->snip_len;
 
    int cmp = strncmp (snip, contig_name, contig_name_len);
    if (!cmp && snip_len != contig_name_len) // identical prefix but different length
        cmp = snip_len - contig_name_len;

    if (cmp < 0) return contigs_get_by_name_do (contigs, contigs_dict, index_buf, STRa(contig_name), mid_sorted_index+1, last_sorted_index);
    if (cmp > 0) return contigs_get_by_name_do (contigs, contigs_dict, index_buf, STRa(contig_name), first_sorted_index, mid_sorted_index-1);

    return word_index;
}             

// binary search for this chrom in contigs. we count on gcc tail recursion optimization to keep this fast.
// looks in Reference contigs if ref is provided, or in CHROM if ref=NULL
WordIndex contigs_get_by_name (ConstContigPkgP ctgs, STRp(contig_name))
{
    ASSERT (ctgs->contigs.len == ctgs->by_name.len, "%s: expecting contigs.len=%"PRIu64" == by_name.len=%"PRIu64, ctgs->name, ctgs->contigs.len, ctgs->by_name.len);

    return contigs_get_by_name_do (&ctgs->contigs, &ctgs->dict, &ctgs->by_name, STRa(contig_name), 0, ctgs->contigs.len-1);
}

// binary search for this chrom in contigs. we count on gcc tail recursion optimization to keep this fast.
static WordIndex contigs_get_by_accession_number_do (ConstBufferP contigs, ConstBufferP index_buf, const AccessionNumber *ac, 
                                                     WordIndex first_sorted_index, WordIndex last_sorted_index)
{
    if (first_sorted_index > last_sorted_index) return WORD_INDEX_NONE; // not found

    WordIndex mid_sorted_index = (first_sorted_index + last_sorted_index) / 2;
    
    WordIndex word_index = *ENT (WordIndex, *index_buf, mid_sorted_index);
    Contig *mid_word = ENT (Contig, *contigs, word_index);
 
    int cmp = memcmp (&mid_word->metadata.parsed.ac.AC, ac->AC, sizeof (ac->AC));
    
    // also compare by version, version2, unless ac->version is 0 to indicate no searching by version is needed
    if (!cmp && ac->version) cmp = mid_word->metadata.parsed.ac.version  - ac->version;
    if (!cmp && ac->version) cmp = mid_word->metadata.parsed.ac.version2 - ac->version2;
    
    if (cmp < 0) return contigs_get_by_accession_number_do (contigs, index_buf, ac, mid_sorted_index+1, last_sorted_index);
    if (cmp > 0) return contigs_get_by_accession_number_do (contigs, index_buf, ac, first_sorted_index, mid_sorted_index-1);

    return word_index; // found
}             

static WordIndex contigs_get_by_accession_number (ConstContigPkgP ctgs, STRp (contig))
{
    ASSERT (ctgs->contigs.len == ctgs->by_AC.len, "%s: expecting contigs.len=%"PRIu64" == by_AC.len=%"PRIu64, ctgs->name, ctgs->contigs.len, ctgs->by_AC.len);

    AccessionNumber ac;
    if (!contig_name_to_acc_num (STRa(contig), &ac)) return WORD_INDEX_NONE;

    return contigs_get_by_accession_number_do (&ctgs->contigs, &ctgs->by_AC, &ac, 0, ctgs->contigs.len-1);
}

// search by ac.AC only, ac.version must be set to 0. 
static WordIndex contigs_get_by_accession_number_no_version (ConstContigPkgP ctgs, const AccessionNumber *ac)
{
    ASSERT (ctgs->contigs.len == ctgs->by_AC.len, "%s: expecting contigs.len=%"PRIu64" == by_AC.len=%"PRIu64, ctgs->name, ctgs->contigs.len, ctgs->by_AC.len);

    return contigs_get_by_accession_number_do (&ctgs->contigs, &ctgs->by_AC, ac, 0, ctgs->contigs.len-1);
}

// binary search for this chrom in contigs. we count on gcc tail recursion optimization to keep this fast.
static WordIndex contigs_get_by_ref_index_do (ConstBufferP contigs, ConstBufferP index_buf, WordIndex ref_index, 
                                              WordIndex first_sorted_index, WordIndex last_sorted_index)
{
    if (first_sorted_index > last_sorted_index) return WORD_INDEX_NONE; // not found

    WordIndex mid_sorted_index = (first_sorted_index + last_sorted_index) / 2;
    
    WordIndex word_index = *ENT (WordIndex, *index_buf, mid_sorted_index);
    Contig *mid_word = ENT (Contig, *contigs, word_index);
 
    int cmp = mid_word->ref_index - ref_index;

    if (cmp < 0) return contigs_get_by_ref_index_do (contigs, index_buf, ref_index, mid_sorted_index+1, last_sorted_index);
    if (cmp > 0) return contigs_get_by_ref_index_do (contigs, index_buf, ref_index, first_sorted_index, mid_sorted_index-1);

    return word_index; // found
}             

WordIndex contigs_get_by_ref_index (ConstContigPkgP ctgs,WordIndex ref_index)
{
    ASSERT (ctgs->contigs.len == ctgs->by_ref_index.len, "%s: expecting contigs.len=%"PRIu64" == by_ref_index.len=%"PRIu64, ctgs->name, ctgs->contigs.len, ctgs->by_ref_index.len);
    return contigs_get_by_ref_index_do (&ctgs->contigs, &ctgs->by_ref_index, ref_index, 0, ctgs->contigs.len-1);
}

// get index of a contig matching the given name: check for exact match, then naming style, then accession number.
// if length is given, it is verified.
// note that the function must always return the same result with the same name (so txtheader, seg and chrom_2ref_compress 
// are all consistent), whether or not LN is provided, so we cannot use LN in the decision making.
WordIndex contigs_get_matching (ConstContigPkgP ctgs, STRp(name), PosType LN, /* optional */
                                bool strictly_alt, // only tests for differnet names
                                bool *is_alt) // if not NULL, also search for as-is, and return whether its alt or not
{
    // human chromosomes 1-22,X,Y as they appear, eg, on:  ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz
    static const AccessionNumber NC[25] = { {},
        {.AC="NC_000001"}, {.AC="NC_000002"}, {.AC="NC_000003"}, {.AC="NC_000004"}, {.AC="NC_000005"}, {.AC="NC_000006"},
        {.AC="NC_000007"}, {.AC="NC_000008"}, {.AC="NC_000009"}, {.AC="NC_000010"}, {.AC="NC_000011"}, {.AC="NC_000012"},
        {.AC="NC_000013"}, {.AC="NC_000014"}, {.AC="NC_000015"}, {.AC="NC_000016"}, {.AC="NC_000017"}, {.AC="NC_000018"},
        {.AC="NC_000019"}, {.AC="NC_000020"}, {.AC="NC_000021"}, {.AC="NC_000022"}, {.AC="NC_000023"}, {.AC="NC_000024"}
        };

    #define CHECK_IF_DONE if (ctg_i != WORD_INDEX_NONE) goto finalize;

    #define IS_CHROMOSOME(s,len) (((len) == 1 && (IS_DIGIT ((s)[0]) || ((s)[0]>='W' && (s)[0]<='Z'))) || \
                                  ((len) == 2 && ((IS_DIGIT ((s)[0]) && IS_DIGIT ((s)[1])))))

    #define NC_NUM(s,len) ((len)==2)                   ? (((s)[0]-'0') * 10 + ((s)[1]-'0'))   \
                        : ((s)[0]=='X')                ? 23                 \
                        : ((s)[0]=='Y')                ? 24                 \
                        : ((s)[0]=='W' || (s)[0]=='Z') ? 1000 /* invalid */ \
                        :                                ((s)[0]-'0');

    // if is_alt provided, search for exact name match
    WordIndex ctg_i = !strictly_alt ? contigs_get_by_name (ctgs, STRa(name)) : WORD_INDEX_NONE;

    if (is_alt) *is_alt = (ctg_i == WORD_INDEX_NONE);
    CHECK_IF_DONE;

    // 22 -> chr22 (1->22, X, Y, W, Z chromosomes)
        if (IS_CHROMOSOME (name, name_len)) {

        char chr_chrom[5] = "chr";
        chr_chrom[3] = name[0];
        chr_chrom[4] = (name_len == 2 ? name[1] : 0);

        ctg_i = contigs_get_by_name (ctgs, chr_chrom, name_len+3); 
        CHECK_IF_DONE;

        if (name[0] != 'W' && name[0] != 'Z') {
            int nc_num = NC_NUM (name, name_len);
            if (nc_num <= 24) {
                ctg_i = contigs_get_by_accession_number_no_version (ctgs, &NC[nc_num]); 
                CHECK_IF_DONE;
            }
        }
    }    

    // chr? or chr?? -> ? or ??
    if ((name_len == 4 || name_len == 5) && !memcmp (name, "chr", 3) && IS_CHROMOSOME(&name[3], name_len-3)) {
        ctg_i = contigs_get_by_name (ctgs, &name[3], name_len-3);
        CHECK_IF_DONE;

        int nc_num = NC_NUM (&name[3], name_len-3);
        if (nc_num <= 24) {
            ctg_i = contigs_get_by_accession_number_no_version (ctgs, &NC[nc_num]); 
            CHECK_IF_DONE;
        }
    }

    // NC_.. -> ? or chr?
    if ((name_len == 11 || name_len == 12) && name[0]=='N' && name[1]=='C' && name[2]=='_' && name[9]=='.') {
        uint64_t nc_num; 
        str_get_int_dec (&name[3], 6, &nc_num);

        static const char nc[25][3] = { "", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
                                        "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y" };

        if (nc_num > 0 && nc_num <= 24) {

            // NC_... -> ?
            ctg_i = contigs_get_by_name (ctgs, nc[nc_num], 1 + (nc_num>=10 && nc_num<=22));
            CHECK_IF_DONE;

            char name[6] = { 'c','h','r',nc[nc_num][0],nc[nc_num][1],0 };
            ctg_i = contigs_get_by_name (ctgs, name, 4 + (nc_num>=10 && nc_num<=22));
            CHECK_IF_DONE;
        }
    }

    // M, MT, chrM, chrMT 
    if ((name_len==4 && name[0]=='c' && name[1]=='h' && name[2]=='r' && name[3]=='M') || 
        (name_len==5 && name[0]=='c' && name[1]=='h' && name[2]=='r' && name[3]=='M' && name[4]=='T') || 
        (name_len==1 && name[0]=='M') ||
        (name_len==2 && name[0]=='M' && name[1]=='T')) {

        ctg_i = contigs_get_by_name (ctgs, "chrMT", 5); CHECK_IF_DONE; 
        ctg_i = contigs_get_by_name (ctgs, "chrM",  4); CHECK_IF_DONE; 
        ctg_i = contigs_get_by_name (ctgs, "MT",    2); CHECK_IF_DONE; 
        ctg_i = contigs_get_by_name (ctgs, "M",     1); CHECK_IF_DONE; 
    }

    // for contigs eg "GL000192.1", "chrUn_JTFH01001867v2_decoy", "chr4_gl383528_alt". 
    if (ctg_i == WORD_INDEX_NONE && name_len >= 6) { // 6 is an accession number minimum length 
        ctg_i = contigs_get_by_accession_number (ctgs, name, name_len); 
        CHECK_IF_DONE; 
    }

finalize:
    // if we found a match, check the the length matches
    if (ctg_i != WORD_INDEX_NONE && LN && CONTIG(*ctgs, ctg_i)->max_pos != LN)
        ctg_i = WORD_INDEX_NONE;

    return ctg_i;
}

const char *contigs_get_name (ConstContigPkgP ctgs, WordIndex index, unsigned *contig_name_len /* optional */)
{
    ASSERT (index >= 0 && index < ctgs->contigs.len, "expecting 0 <= index=%d < %s.len=%"PRIu64, index, ctgs->name, ctgs->contigs.len);

    Contig *ctg = ENT (Contig, ctgs->contigs, index);

    if (contig_name_len) *contig_name_len = ctg->snip_len;
    
    return ENT (const char, ctgs->dict, ctg->char_index);
}

// -----------------------------
// iterator
// -----------------------------

// call a callback for each accessed contig. Note: callback function is the same as sam_foreach_SQ_line
void foreach_contig (ConstContigPkgP ctgs, ContigsIteratorCallback callback, void *callback_param)
{
    ASSERTNOTNULL (ctgs);

    ARRAY (const Contig, ctg, ctgs->contigs);

    for (uint64_t i=0; i < ctg_len; i++) {
        const char *ref_chrom_name = ENT (const char, ctgs->dict, ctg[i].char_index); 
        callback (ref_chrom_name, strlen (ref_chrom_name), ctg[i].max_pos, callback_param);
    }
}

// -----------------------------
// initialization & finalization
// -----------------------------

void contigs_create_index (ContigPkg *ctgs, SortBy sort_by)
{
    ASSERTNOTNULL (ctgs);

    if (sort_by & SORT_BY_NAME) 
        contigs_sort (&ctgs->contigs, &ctgs->dict, &ctgs->by_name, SORT_BY_NAME);

    if (sort_by & SORT_BY_REF_INDEX) 
        contigs_sort (&ctgs->contigs, &ctgs->dict, &ctgs->by_ref_index, SORT_BY_REF_INDEX);

    if (sort_by & SORT_BY_AC) {
        contigs_calculate_accession_numbers (&ctgs->contigs, &ctgs->dict);
        contigs_sort (&ctgs->contigs, &ctgs->dict, &ctgs->by_AC, SORT_BY_AC);
    }  

    ctgs->sorted_by |= sort_by;
}

void contigs_build_contig_pkg_from_ctx (ContigPkg *ctgs, ConstContextP ctx, SortBy sort_by)
{
    ASSERTNOTNULL (ctgs);

    buf_copy (evb, &ctgs->dict, &ctx->dict, char, 0, 0, "ContigPkg->dict");
    
    // note: we take num_contigs from nodes (ZIP) or word_list (PIZ), but in the case of nodes, we can't rely
    // on chrom_ctx->nodes for the correct order as vb_i=1 resorted. so we generate directly from dict
    uint32_t num_contigs = MAX_(ctx->nodes.len, ctx->word_list.len);
    uint32_t num_counts  = ctx->counts.len;

    ASSERT (!num_counts || ctx->counts.len == num_contigs, "expecting ctx=%s to have num_contigs=%u == num_counts=%u", ctx->tag_name, num_contigs, num_counts);

    ctgs->contigs.len = ctgs->dict.param = num_contigs;
    ctgs->has_counts = (num_counts > 0);

    buf_alloc_zero (evb, &ctgs->contigs, 0, num_contigs, Contig, 1, "ContigPkg->contigs");

    // similar logic to ctx_dict_build_word_lists
    char *start = ctgs->dict.data;
    Contig *contig = FIRSTENT (Contig, ctgs->contigs);
    for (uint32_t i=0; i < num_contigs; i++, contig++) {

        if (num_counts)
            contig->metadata.parsed.count = *ENT (int64_t, ctx->counts, i);

        char *c=start; while (*c) c++;
        contig->snip_len   = c - start;
        contig->char_index = start - ctgs->dict.data;

        const char *contig_name = ENT (char, ctgs->dict, contig->char_index);
        contig_name_to_acc_num (contig_name, contig->snip_len, &contig->metadata.parsed.ac);

        start = c+1; // skip over the \0
    }

    contigs_create_index (ctgs, sort_by);    
} 

void contigs_free (ContigPkg *ctgs)
{
    if (!ctgs) return;

    buf_free (&ctgs->contigs);
    buf_free (&ctgs->dict);
    buf_free (&ctgs->by_name);
    buf_free (&ctgs->by_LN);
    buf_free (&ctgs->by_AC);
    buf_free (&ctgs->by_ref_index);
    ctgs->sorted_by = 0;
}

void contigs_destroy (ContigPkg *ctgs)
{
    if (!ctgs) return;

    buf_destroy (&ctgs->contigs);
    buf_destroy (&ctgs->dict);
    buf_destroy (&ctgs->by_name);
    buf_destroy (&ctgs->by_LN);
    buf_destroy (&ctgs->by_AC);
    buf_destroy (&ctgs->by_ref_index);
}

