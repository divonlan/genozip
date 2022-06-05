// ------------------------------------------------------------------
//   kraken.c
//   Copyright (C) 2021-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
// handle kraken output format: https://github.com/DerrickWood/kraken/wiki/Manual#output-formats

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "context.h"
#include "file.h"
#include "endianness.h"
#include "strings.h"
#include "piz.h"
#include "dict_id.h"
#include "stats.h"
#include "codec.h"
#include "version.h"
#include "kraken.h"
#include "reconstruct.h"
#include "compressor.h"
#include "zfile.h"
#include "hash.h"
#include "progress.h"
#include "website.h"
#include "qname.h"

// Search algorithm:
// A. During kraken loading:
//    1. qname_dict is reconstructed directly from TOP2TAXID - concatenation of txt_data from all VBs
//    2. qname_nodes is created per VB during loading of kraken data - during reconstruction  
//       -- only items that that pass the --taxid filter are entered into qname_dict and qname_nodes
//    3. during kraken loading, after reconstruction, A. qsort qname_nodes by hash
//       B. generate qname_hashtab - an array indexed by hash, containing an index into qname_nodes.
// B. At toplevel callback of SAM/BAM, FASTQ, FASTA and Kraken:
//    Call kraken_is_included_loaded - this will check if the qname is in qname_dict by using the hash table and sorted list

// Used by Load

#define QNAME_NODE_NONE 0xffffffff // entry in qname_hashtab indicating this qname is not in dict
#define MAX_QNAME_NODES 0xfffffffe

#define BITS_CHAR_INDEX 40
#define BITS_SNIP_LEN   7
#define BITS_TAXID      24
#define BITS_TAXID_LO   16

#define MAX_SNIP_LEN   ((1 << BITS_SNIP_LEN) - 1)
#define MAX_TAXID      ((1 << BITS_TAXID) - 1)
#define MAX_CHAR_INDEX ((1ULL << BITS_CHAR_INDEX) -1)

// 13 bytes
typedef struct __attribute__((packed, aligned(1))) {
    uint32_t hash;                          // hash of the QNAME - qname_nodes is sorted by hash, and it also the index into qname_hashtab
    uint64_t char_index : BITS_CHAR_INDEX;  // index into qname_dict of the qname
    uint64_t snip_len   : BITS_SNIP_LEN;
    uint64_t is_paired  : 1;
    uint64_t taxid_lo   : BITS_TAXID_LO;
    uint8_t  taxid_hi;
} QnameNode;

char *kraken_filename = NULL;               // global - kraken filename. Non-NULL indicates kraken is loaded.
static bool taxid_negative = false;         // user specified ^
static bool taxid_also_0   = false;         // user specified +0

static Buffer qname_dict    = EMPTY_BUFFER;
static Buffer qname_nodes   = EMPTY_BUFFER; // array of QnameNode
static Buffer qname_hashtab = EMPTY_BUFFER; // array of uint32_t - index into qname_nodes
static unsigned one_qname_length = 0;       // average length of one QNAME
static TaxonomyId dom_taxid = TAXID_NONE;   // most common taxid in kraken file - eliminated from qname_nodes

// Used by Seg
static char copy_taxid_snip[30];
static unsigned copy_taxid_snip_len;

void kraken_set_show_kraken (rom optarg)
{
    if (!optarg) 
        flag.show_kraken = KRK_ALL;
    else if (str_case_compare (optarg, "included", NULL)) 
        flag.show_kraken = KRK_INCLUDED;
    else if (str_case_compare (optarg, "excluded", NULL)) 
        flag.show_kraken = KRK_EXCLUDED;
    else
        ASSINP0 (false, "--show-kraken argument error: see "WEBSITE_KRAKEN" for valid options");
}

//-----------------------
// Segmentation functions
//-----------------------

void kraken_zip_initialize (void)
{
    seg_prepare_snip_other (SNIP_COPY, _KRAKEN_TAXID, 0, 0, copy_taxid_snip);

    qname_zip_initialize (KRAKEN_QNAME);
}

void kraken_seg_initialize (VBlockP vb)
{
    CTX(KRAKEN_TAXID)->flags.store    = STORE_INT;
    CTX(KRAKEN_TAXID)->no_stons       = true; // must be no_stons the SEC_COUNTS data needs to mirror the dictionary words
    CTX(KRAKEN_TAXID)->counts_section = true; 

    stats_set_consolidation (vb, KRAKEN_KMERS, 2, KRAKEN_KMERTAX, KRAKEN_KMERLEN);

    qname_seg_initialize (VB, KRAKEN_QNAME);
}

void kraken_seg_finalize (VBlockP vb)
{
    // top level snip - IMPORTNAT - kraken_piz_filter needs to changed and support backward compatability
    SmallContainer top_level = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .callback     = true,
        .nitems_lo    = 6,
        .items        = { { .dict_id = { _KRAKEN_CU },     .separator = {'\t'} },
                          { .dict_id = { _KRAKEN_QNAME },  .separator = {'\t'} },
                          { .dict_id = { _KRAKEN_TAXID },  .separator = {'\t'} },
                          { .dict_id = { _KRAKEN_SEQLEN }, .separator = {'\t'} },
                          { .dict_id = { _KRAKEN_KMERS },                      },
                          { .dict_id = { _KRAKEN_EOL },                        } }
    };

    container_seg (vb, CTX(KRAKEN_TOPLEVEL), (ContainerP)&top_level, 0, 0, 0);

    // top level container when loading a kraken file with --kraken
    SmallContainer top2taxid = { 
        .repeats        = vb->lines.len,
        .is_toplevel    = true,
        .callback       = true,
        .nitems_lo      = 2,
        .items          = { { .dict_id = { _KRAKEN_QNAME },  .separator = { CI0_TRANS_NUL /* '\0' */} },
                            { .dict_id = { _KRAKEN_TAXID },  .separator = { CI0_TRANS_NOR /* no reconstruct */ } } },
    };

    container_seg (vb, CTX(KRAKEN_TOP2TAXID), (ContainerP)&top2taxid, 0, 0, 0);
}

// ZIP: called by main thread after compute has finished. 
void kraken_zip_after_compute (VBlockP vb)
{
    // add up the total length of all QNAMEs in the file - will be transferred to PIZ via SectionHeaderCounts.nodes_param of TAXID
    ZCTX(KRAKEN_TAXID)->nodes.count += vb->last_int(KRAKEN_QNAME); 
}

bool kraken_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return 
        // typically small 
        dict_id.num == _KRAKEN_CU       ||
        dict_id.num == _KRAKEN_QNAME    || // container
        dict_id.num == _KRAKEN_TAXID    || // usually a few thousand species
        dict_id.num == _KRAKEN_SEQLEN   ||
        dict_id.num == _KRAKEN_KMERS    || // container
        dict_id.num == _KRAKEN_EOL      ||
        dict_id.num == _KRAKEN_TOPLEVEL ||
        dict_id.num == _KRAKEN_TOP2TAXID;
}

static void kraken_seg_kmers (VBlockP vb, rom value, int32_t value_len, rom taxid, unsigned taxid_len) // must be signed
{
    bool final_space = (value_len && value[value_len-1] == ' '); // sometimes there is, sometimes there isn't...
    value_len -= final_space;

    unsigned num_kmers=1;
    for (unsigned i=0; i < value_len; i++)
        if (value[i] == ' ') {
            if (i && value[i-1] == ' ') { // not valid khmers
                seg_by_did_i (VB, value, value_len, KRAKEN_KMERS, value_len);
                return;
            }
            num_kmers++;
        }

    SmallContainer kmers_con = (SmallContainer){ 
        .repeats   = num_kmers, 
        .nitems_lo = 2,
        .repsep    = {' '},
        .items     = { { .dict_id = { _KRAKEN_KMERTAX }, .separator = {':'} },
                       { .dict_id = { _KRAKEN_KMERLEN }                     } },
        .drop_final_repeat_sep = !final_space
    };


    str_split_enforce (value, value_len, num_kmers, ' ', kmer, true, "kmers");

    for (unsigned k=0; k < num_kmers; k++) {
        str_split_enforce (kmers[k], kmer_lens[k], 2, ':', item, true, "kmer items");

        // KMER taxid - either copy from TAXID or seg a normal snip
        if (taxid_len == item_lens[0] && !memcmp (taxid, items[0], taxid_len)) 
            seg_by_did_i (VB, copy_taxid_snip, copy_taxid_snip_len, KRAKEN_KMERTAX, item_lens[0]);
        else 
            seg_by_did_i (VB, items[0], item_lens[0], KRAKEN_KMERTAX, item_lens[0]);

        // KMER length
        seg_by_did_i (VB, items[1], item_lens[1], KRAKEN_KMERLEN, item_lens[1]);
    }

    container_seg (vb, CTX(KRAKEN_KMERS), (ContainerP)&kmers_con, 0, 0, num_kmers*2-1 + final_space); // account for ':' within kmers and ' ' betweem them
}

// example: "C       ST-E00180:535:HCNW2CCX2:8:1101:16183:1221       570     150|150 A:1 570:21 543:5 91347:8 0:81 |:| A:1 543:4 570:25 0:32 28384:5 0:49"
rom kraken_seg_txt_line (VBlockP vb, rom field_start_line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    rom next_field=field_start_line, field_start;
    unsigned field_len=0;
    char separator;

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    SEG_NEXT_ITEM (KRAKEN_CU);

    // QNAME - We break down the QNAME into subfields separated by / and/or : - these are vendor-defined strings. See examples in sam_seg_txt_line()
    GET_NEXT_ITEM (KRAKEN_QNAME);
    qname_seg (VB, CTX(KRAKEN_QNAME), field_start, field_len, 1 /* \t */);

    vb->last_int(KRAKEN_QNAME) += field_len+1; // count total QNAME lengths in this VB (+1 for separator)

    SEG_NEXT_ITEM (KRAKEN_TAXID);
    GET_NEXT_ITEM (KRAKEN_SEQLEN); // for paired files we will have eg "150|149", for non-paired eg "150"
    if (memchr (field_start, '|', field_len)) { 
        // struct and not array, bc SEQLEN_2 is often a lot more random than SEQLEN_1 so better use different contexts
        static const MediumContainer con_SEQLEN = { .nitems_lo = 2,      
                                                    .items     = { { .dict_id = { _KRAKEN_SEQLEN_1 }, .separator = {'|'} },  
                                                                   { .dict_id = { _KRAKEN_SEQLEN_2 },                    } } };

        seg_array_of_struct (VB, CTX(KRAKEN_SEQLEN), con_SEQLEN, field_start, field_len, (SegCallback[]){seg_pos_field_cb, 0}); // first element is good to delta, second is not
    }
    else 
        seg_pos_field_cb (VB, CTX(KRAKEN_SEQLEN), field_start, field_len, 0);

    CTX(KRAKEN_SEQLEN)->txt_len++; // \t

    GET_LAST_ITEM (KRAKEN_KMERS);
    kraken_seg_kmers (vb, field_start, field_len, KRAKEN_TAXID_str, KRAKEN_TAXID_len);
    
    SEG_EOL (KRAKEN_EOL, true);

    return next_field;
}

//----------------------------------------------------------------------------------------------------
// PIZ functions when loading as a result of --kraken
// Note: when just pizzing a kraken normally, it is a plain vanilla piz without any special functions
//----------------------------------------------------------------------------------------------------

bool kraken_is_translation (VBlockP vb)
{
    return flag.reading_kraken;
    //return flag.kraken_taxid > 0;
}

// returns true if section is to be skipped reading / uncompressing
IS_SKIP (kraken_piz_is_skip_section)
{
    if (!dict_id.num) return false; // only attempt to skip B250/LOCAL/COUNT sections

    if (flag.reading_kraken && // note: we only need some of the fields when ingesting loading kraken data
        dict_id.num != _KRAKEN_QNAME     &&
        !dict_id_is_qname_sf(dict_id)    && 
        dict_id.num != _KRAKEN_TAXID     &&
        dict_id.num != _KRAKEN_EOL       &&
        dict_id.num != _KRAKEN_TOP2TAXID &&
        dict_id_typeless (dict_id).id[0] != 'Q') // array item and compound items of QNAME (in files compressed with 12.0.42 or earlier)  
        return true;

    return false;
}

void kraken_piz_handover_data (VBlockP vb)
{
    uint64_t dict_start = qname_dict.len;
    buf_add_buf (evb, &qname_dict, &vb->txt_data, char, "qname_dict");

    ARRAY (QnameNode, nodes, CTX(KRAKEN_QNAME)->qname_nodes);
    for (uint64_t i=0; i < nodes_len; i++)
        nodes[i].char_index += dict_start;

buf_add_buf (evb, &qname_nodes, &CTX(KRAKEN_QNAME)->qname_nodes, QnameNode, "qname_nodes");
}

// callback called after every repeat of TOP2HASH, i.e. when run with --taxid
// inspect QNAME and decide if to include in the qname_nodes or not
CONTAINER_CALLBACK (kraken_piz_container_cb)
{
    // when loading kraken to filter another file, remove lines according to the requested filter
    if (dict_id.num == _KRAKEN_TOP2TAXID) {

        uint32_t this_taxid = vb->last_int (KRAKEN_TAXID);

        // we store everything EXCEPT the lines with the most dominant taxid. Missing qnode for a line would indicate
        // it has the dominant taxid (under the assumption that all lines in the filtered file are covered by lines in the kraken file)
        if (this_taxid != dom_taxid) {

            // if this is a paired read - ending with /1 or /2 - remove the suffix
            bool is_paired = false;
            unsigned snip_len = recon_len - 1; // excluding the \0 separator
            
            if (recon_len > 3 && recon[snip_len-2] == '/' &&
                (recon[recon_len-2] == '1' || recon[snip_len-1] == '2')) {
                is_paired = true;
                snip_len -= 2;
            }

            buf_alloc (vb, &CTX(KRAKEN_QNAME)->qname_nodes, 1, vb->lines.len, QnameNode, 1.5, "contexts->qname_nodes");
            
            ASSERT (this_taxid <= MAX_TAXID, "taxid=%u exceeds maximum of %u", this_taxid, MAX_TAXID);
            ASSERT (snip_len <= MAX_SNIP_LEN, "Length of QNAME \"%.*s\" exceeds maximum of %u", snip_len, recon, MAX_SNIP_LEN);

            BNXT (QnameNode, CTX(KRAKEN_QNAME)->qname_nodes) = (QnameNode){ 
                .char_index = recon - vb->txt_data.data, // will be updated in kraken_piz_handover_data
                .snip_len   = snip_len,
                .hash       = hash_do (qname_hashtab.len, recon, snip_len),
                .is_paired  = z_file->num_components == 2 || is_paired,
                .taxid_lo   = this_taxid & ((1<<BITS_TAXID_LO)-1),
                .taxid_hi   = this_taxid >> BITS_TAXID_LO
            };
        }
        else
            vb->drop_curr_line = "dom_taxid"; // txt_data will turn into qname_dict - so we eliminate unused qnames
    } 

    // when pizzing a kraken to be filtered by itself (loaded or stored) (useful only for testing) - apply :kraken_is_included_loaded" filter 
    else if (flag.kraken_taxid != TAXID_NONE && 
             dict_id.num == _KRAKEN_TOPLEVEL && 
             (   ( kraken_is_loaded && !kraken_is_included_loaded (vb, last_txt (vb, KRAKEN_QNAME), vb->last_txt_len (KRAKEN_QNAME)))
              || (!kraken_is_loaded && !kraken_is_included_stored (vb, KRAKEN_TAXID, true)))) // TAXID was recon in the TOPLEVEL container
        
        vb->drop_curr_line = "taxid";
}

// Loading kraken: called after global area is read (including SEC_COUNTS), before reconstruction of VBs
bool kraken_piz_initialize (void)
{
    if (!flag.reading_kraken) return true; // proceed with simple PIZ of the kraken file

    ASSERTNOTINUSE (qname_hashtab);
    ASSERTNOTINUSE (qname_nodes);
    ASSERTNOTINUSE (qname_dict);
    
    ASSERT (z_file->num_components <= 2, "--kraken requires a .kraken.genozip file with up to 2 components, but %s has %u components",
            z_name, z_file->num_components);

    Context *zctx = ZCTX(KRAKEN_TAXID);
    ARRAY (uint64_t, counts, zctx->counts);

    // verify that the user selected taxonomy ID is in the kraken data
    if (exe_type == EXE_GENOCAT) {
        ASSINP0 (flag.kraken_taxid != TAXID_NONE, "--taxid must be provided if --kraken is used");

        WordIndex taxid_word_i = ctx_get_word_index_by_snip (evb, zctx, str_int_s (flag.kraken_taxid).s, 0);
        if (taxid_word_i == WORD_INDEX_NONE) {
            progress_finalize_component ("Skipped");
            WARN ("FYI: %s has no sequences with a Taxanomic ID of \"%d\"", z_name, flag.kraken_taxid);
            return false;
        }
    }

    // calculate the total number of sequences in the kraken file (and the FASTQ file(s) from which it was generated)
    uint64_t total_sequences_in_file=0; 
    uint64_t max_count = 0;
    WordIndex max_count_word_index = WORD_INDEX_NONE;

    for (uint32_t i=0; i < counts_len; i++) {
        total_sequences_in_file += counts[i];
        if (counts[i] > max_count) {
            max_count = counts[i];
            max_count_word_index = i;
        }
    }
    uint64_t num_non_dom_seqs = total_sequences_in_file - max_count;

    dom_taxid = atoi (ctx_get_words_snip (zctx, max_count_word_index));

    ASSERT0 (total_sequences_in_file, "unexpectedly, total_sequences_in_file=0");

    // calculate the average length of a QNAME string in the file
    int64_t total_qname_length = zctx->nodes.count; // sent ZIP->PIZ via SectionHeaderCounts.node_param
    one_qname_length = (unsigned)(1 + total_qname_length / total_sequences_in_file); // average (rounded up) QNAME length in the file
  
    // allocate
    qname_nodes.can_be_big = qname_dict.can_be_big = true; // we're might be allocating a huge amount of memory - suppress buf_alloc warning
    buf_alloc (evb, &qname_nodes, 0, num_non_dom_seqs, QnameNode, 0, "qname_nodes"); 
    buf_alloc (evb, &qname_dict,  0, num_non_dom_seqs, char[one_qname_length+2],  0, "qname_dict"); // approximate (+2), might grow if needed in kraken_piz_handover_data

    qname_hashtab.len = hash_next_size_up (3 * num_non_dom_seqs, true);

    buf_alloc (evb, &qname_hashtab, 0, qname_hashtab.len, uint32_t, 1, "qname_hashtab");
    buf_set (&qname_hashtab, 0xff); // empty

//printf ("total_sequences_in_file=%u num_non_dom_seqs=%u max_count=%u qname_hashtab.size=%"PRIu64" one_qname_length=%u qname_nodes.size=%"PRIu64" qname_dict.size=%"PRIu64"\n", (int)total_sequences_in_file, (int)num_non_dom_seqs, (int)max_count, qname_hashtab.size, (int)one_qname_length, qname_nodes.size, qname_dict.size);    

    return true; // proceed with PIZ of kraken file
}

static ASCENDING_SORTER (kraken_qname_nodes_cmp, QnameNode, hash)

// genocat: load kraken file as a result of genocat --kraken
void kraken_load (void)
{
    ASSERTNOTNULL (flag.reading_kraken);
    SAVE_FLAGS_AUX;

    flag.maybe_vb_modified_by_reconstructor = true;    // we drop the lines not matching --taxid
    flag.data_modified    = true;    // the reconstructed kraken file is not the same as the original...
    flag.no_writer = flag.no_writer_thread = true;
    flag.genocat_no_reconstruct = false;
    flag.genocat_global_area_only = false;
    flag.out_dt           = DT_NONE; // needed for dt_get_translation to find the translation defined in TRANSLATIONS
    flag.reference        = REF_NONE;
    TEMP_VALUE (command, PIZ);

    z_file = file_open (flag.reading_kraken, READ, Z_FILE, DT_KRAKEN);    
    zfile_read_genozip_header (0);
    
    ASSINP (Z_DT(DT_KRAKEN), "expected %s to be a genozip'ed kraken output file, but its a %s file. Tip: compress the output file generated by kraken2 with \"genozip --input kraken\"", 
            z_name, dt_name (z_file->data_type));

    z_file->basename = file_basename (flag.reading_kraken, false, "(kraken-file)", NULL, 0);

    Dispatcher dispachter = piz_z_file_initialize();
    bool kraken_loaded = piz_one_txt_file (dispachter, false, false, COMP_NONE);

    kraken_filename = file_make_unix_filename (z_name); // full-path unix-style filename, allocates memory

    file_close (&z_file, false, false);
    file_close (&txt_file, false, false); // close the txt_file object we created (even though we didn't open the physical file). it was created in file_open called from txtheader_piz_read_and_reconstruct.
    
    // note: kraken_loaded is false if flag.kraken_taxid is absent from kraken data (then kraken_piz_initialize aborts pizzing the data) 
    if (kraken_loaded) {
        qsort (STRb(qname_nodes), sizeof (QnameNode), kraken_qname_nodes_cmp);

        ASSERT (qname_nodes.len <= MAX_QNAME_NODES, "qname_nodes.len=%"PRIu64" exceeds the maximum of %u", qname_nodes.len, MAX_QNAME_NODES);

        // build the hash table
        for (uint32_t i=0; i < qname_nodes.len32; i++) {
            uint32_t hash = B(QnameNode, qname_nodes, i)->hash;
            *B32 (qname_hashtab, hash) = i; // overwriting an empty or populated hash entry 
        }
    }

    // recover globals
    RESTORE_VALUE (command);
    RESTORE_FLAGS;

    flag.reading_kraken = NULL;
}

// reset to factory defaults
void kraken_destroy (void)
{
    buf_destroy (qname_dict);
    buf_destroy (qname_nodes);
    buf_destroy (qname_hashtab);
    kraken_filename = NULL;
    taxid_negative = taxid_also_0 = false;
    one_qname_length = 0;
}

// -------------------------------------------
// using the kraken data in genocat --kraken
// -------------------------------------------

void kraken_set_taxid (rom optarg)
{
    unsigned len = strlen (optarg);

    taxid_negative = optarg[0] == '^';
    taxid_also_0 = (len > 3 && optarg[len-2] == '+' && optarg[len-1] == '0');

    str_get_int_range32 (optarg + taxid_negative, len - taxid_negative - 2 * taxid_also_0, 
                         0, 9999999, &flag.kraken_taxid); // NCBI Taxonomic ID is a 7 digit positive integer (and "0" means unclassified in kraken)
}

static inline TaxonomyId get_taxid (const QnameNode *l)
{
    return (l->taxid_hi << BITS_TAXID_LO) | l->taxid_lo;
}

static inline const QnameNode *kraken_search_up_and_down (const QnameNode *listent, uint32_t hash, STRp(qname))
{
    // search down
    for (const QnameNode *l=listent-1; l >= B1ST (QnameNode, qname_nodes) && l->hash == hash; l--)
        if (l->snip_len == qname_len && !memcmp (qname, Bc (qname_dict, l->char_index), qname_len))
            return l;

    // search up
    for (const QnameNode *l=listent+1; l <= BLST (QnameNode, qname_nodes) && l->hash == hash; l++)
        if (l->snip_len == qname_len && !memcmp (qname, Bc (qname_dict, l->char_index), qname_len))
            return l;

    return NULL;
}

// search for an additional entry for the same qname - this could either originate from a file with
// 2 components, or reads with /1 /2
static TaxonomyId kraken_get_taxid_with_pair (const QnameNode *l1, uint32_t hash, STRp(qname))
{
    TaxonomyId l1_taxid = get_taxid (l1);

    if (exe_type == EXE_GENOCAT &&  l1_taxid == flag.kraken_taxid) 
        return l1_taxid; // we found the taxid in l1, no need to check l2

    const QnameNode *l2 = kraken_search_up_and_down (l1, hash, STRa(qname));
    TaxonomyId l2_taxid = l2 ? get_taxid (l2) : dom_taxid; // if paired qnode is missing - it is the dominant tax id;

    return l2_taxid != TAXID_UNCLASSIFIED ? l2_taxid  // if t2 is classified (in genocat, possibly kraken_taxid), we return it
         : l1_taxid != TAXID_UNCLASSIFIED ? l1_taxid  // or, if t1 is classified (definitely not kraken_taxid), we return it
         :                                  TAXID_UNCLASSIFIED; // both are unclassified - result is unclassified
}

// genozip --kraken: returns the taxid of the read. in case of paired - returns the read which is classified.
//                   if both are classified, returns r2.
// genocat --kraken: same, except that if both reads are classified, returns flag.kraken_taxid if one of the reads equals it
static TaxonomyId kraken_get_taxid (STRp(qname))
{
    uint32_t hash = hash_do (qname_hashtab.len32, STRa(qname));

    uint32_t list_index = *B32 (qname_hashtab, hash);
    if (list_index == QNAME_NODE_NONE) 
        return dom_taxid; // taxid not found at all - it is dom_taxid

    ASSERT (list_index < qname_nodes.len, "expecting list_index=%u < qname_nodes.len=%"PRIu64, list_index, qname_nodes.len);

    // hash found, it could be the requested qname, or a different one with the same hash. let's check
    const QnameNode *listent = B(QnameNode, qname_nodes, list_index);
    rom snip = Bc (qname_dict, listent->char_index);

    if (listent->snip_len == qname_len && !memcmp (qname, snip, qname_len))
        return listent->is_paired ? kraken_get_taxid_with_pair (listent, hash, STRa(qname)) : get_taxid (listent);

    // we have a hash table contention - the hash file points to a different QNAME with the same hash. 
    // we look for all other entries with the same hash. these are adajent to listent as qname_nodes is sorted by hash

    const QnameNode *l = kraken_search_up_and_down (listent, hash, STRa(qname));
    if (!l) return dom_taxid; // taxid not found at all - it is dom_taxid

    return l->is_paired ? kraken_get_taxid_with_pair (l, hash, STRa(qname)) : get_taxid (l);
}

// Find whether QNAME is included in the taxid filter in O(1) (using the loaded kraken file)
bool kraken_is_included_loaded (VBlockP vb, STRp(qname))
{
    ASSINP0 (flag.kraken_taxid != TAXID_NONE, "When using --kraken, you must specify --taxid <number> (positive filter) or --taxid ^<number> (negative filter)");
    ASSINP0 (!z_file->z_flags.has_taxid || VB_DT(DT_KRAKEN), "You cannot use --kraken, because the file already contains taxonomic information (it was compressed with --kraken) - just use --taxid");

    // case: kraken data is not loaded, because kraken_piz_initialize determined that kraken_taxid is absent
    // from the kraken data
    if (!qname_hashtab.len) 
        return taxid_negative;

    TaxonomyId qname_taxid = kraken_get_taxid (STRa(qname));
    bool is_requested_taxid = (qname_taxid == flag.kraken_taxid)
                           || (taxid_also_0 && !qname_taxid);

    bool res = (taxid_negative == !is_requested_taxid);

    if (flag.show_kraken) {
        if (res && flag.show_kraken != KRK_EXCLUDED)
            iprintf ("%.*s\tINCLUDED\n", STRf(qname));
        else if (!res && flag.show_kraken != KRK_INCLUDED)
            iprintf ("%.*s\tEXCLUDED\n", STRf(qname));
    }

    return res;
}

// Find whether QNAME is included in the taxid filter in O(1) (using the loaded kraken file)
bool kraken_is_included_stored (VBlockP vb, DidIType did_i_taxid, bool already_reconstructed)
{
    if (flag.kraken_taxid == TAXID_NONE) return true; // user didn't specify --taxid - everything's included

    ASSINP0 (z_file->z_flags.has_taxid, "To use --taxid, either use \"genocat --kraken <file> --taxid <taxid>\" to load a kraken file, or use \"genozip --kraken <file>\" to include taxonomic information in the genozip file");

    if (!already_reconstructed) 
        reconstruct_from_ctx (vb, did_i_taxid, 0, false); // update last_int without reconstructing
    
    uint32_t this_taxid = vb->last_int (did_i_taxid);

    return  (( taxid_negative && flag.kraken_taxid != this_taxid) ||
             (!taxid_negative && flag.kraken_taxid == this_taxid));
}

//----------------------------------------------------------------------------------------------------
// SEG helper functions for other data types segging with genozip --kraken
//----------------------------------------------------------------------------------------------------

// returns snip_len if successful, or 0 if failed (only happens if fail_if_missing)
unsigned kraken_seg_taxid_do (VBlockP vb, DidIType did_i_taxid, STRp(qname), 
                              char *snip, // caller-allocated out
                              bool fail_if_missing)
{
    TaxonomyId taxid = kraken_get_taxid (STRa(qname));

    ctx_set_last_value (vb, CTX(did_i_taxid), (ValueType){.i = taxid });

    ASSINP (!fail_if_missing || taxid != TAXID_NONE, "Cannot find taxonomy id in kraken file for QNAME \"%.*s\"", STRf(qname));

    if (taxid == TAXID_NONE) 
        return 0; // failed

    else {
        unsigned snip_len = str_int (taxid, snip);

        if (flag.show_kraken)
            iprintf ("%.*s\t%d\n", STRf(qname), taxid);

        seg_by_did_i (VB, STRa(snip), did_i_taxid, 0);

        return snip_len;
    }
}

// always segs if even_if_missing ; or if even_if_missing=false, returns true if taxid found and segged
// returns snip_len if successful, or 0 if failed (only happens if fail_if_missing)
unsigned kraken_seg_taxid (VBlockP vb, DidIType did_i_taxid, STRp(qname), bool fail_if_missing)
{
    char snip[20];
    return kraken_seg_taxid_do (vb, did_i_taxid, STRa(qname), snip, fail_if_missing);
}
