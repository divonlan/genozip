// ------------------------------------------------------------------
//   kraken.c
//   Copyright (C) 2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
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

// Data structures used during Kraken PIZ 

char *kraken_filename = NULL; // global - kraken filename

// Search algorithm:
// A. During kraken loading:
//    1. qname_dict is reconstructed directly from TOP2TAXID - concatenation of txt_data from all VBs
//    2. qname_nodes is created per VB during loading of kraken data - during reconstruction  
//       -- only items that that pass the --taxid filter are entered into qname_dict and qname_nodes
//    3. during kraken loading, after reconstruction, A. qsort qname_nodes by hash
//       B. generate qname_hashtab - an array indexed by hash, containing an index into qname_nodes.
// B. At toplevel callback of SAM/BAM, FASTQ, FASTA and Kraken:
//    Call kraken_is_included_loaded - this will check if the qname is in qname_dict by using the hash table and sorted list

#define MAX_QNAME_NODES 0xfffffffe
static Buffer qname_dict    = EMPTY_BUFFER;
static Buffer qname_nodes   = EMPTY_BUFFER; // array of QnameNode
static enum { QNAME_INCLUDE, QNAME_EXCLUDE } qname_type;

static Buffer qname_hashtab = EMPTY_BUFFER; // array of uint32_t - index into qname_nodes
static unsigned one_qname_length = 0; // average length of one QNAME

// 16 bytes
typedef struct __attribute__((packed, aligned(4))) {
    uint32_t hash;            // hash of the QNAME - qname_nodes is sorted by hash, and it also the index into qname_hashtab
    uint32_t taxid;
    uint64_t char_index : 48; // index into qname_dict of the qname
    uint64_t snip_len   : 16;
} QnameNode;

#define QNAME_NODE_NONE 0xffffffff // entry in qname_hashtab indicating this qname is not in dict

static char copy_taxid_snip[30];
static unsigned copy_taxid_snip_len;

//-----------------------
// Segmentation functions
//-----------------------

void kraken_zip_initialize (void)
{
    copy_taxid_snip_len = sizeof (copy_taxid_snip);
    seg_prepare_snip_other (SNIP_OTHER_COPY, (DictId)dict_id_fields[KRAKEN_TAXID], 0, 0, copy_taxid_snip, &copy_taxid_snip_len);
}

void kraken_seg_initialize (VBlock *vb)
{
    vb->contexts[KRAKEN_TAXID].flags.store    = STORE_INT;
    vb->contexts[KRAKEN_TAXID].no_stons       = true; // must be no_stons the SEC_COUNTS data needs to mirror the dictionary words
    vb->contexts[KRAKEN_TAXID].counts_section = true; 

    vb->contexts[KRAKEN_KMERTAX].st_did_i   = KRAKEN_KMERS;
    vb->contexts[KRAKEN_KMERLEN].st_did_i   = KRAKEN_KMERS;

    vb->contexts[KRAKEN_KMERS].no_stons     =  // container context
    vb->contexts[KRAKEN_QNAME].no_stons     =  // container context
    vb->contexts[KRAKEN_TOPLEVEL].no_stons  = 
    vb->contexts[KRAKEN_TOP2TAXID].no_stons = 
    vb->contexts[KRAKEN_EOL].no_stons       = 
    vb->contexts[KRAKEN_CU].no_stons        = true; // no singletons, so b250 can be optimized away 
}

void kraken_seg_finalize (VBlockP vb)
{
    // top level snip - IMPORTNAT - kraken_piz_filter needs to changed and support backward compatability
    SmallContainer top_level = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .callback     = true,
        .nitems_lo    = 6,
        .items        = { { .dict_id = (DictId)dict_id_fields[KRAKEN_CU],     .seperator = {'\t'} },
                          { .dict_id = (DictId)dict_id_fields[KRAKEN_QNAME],  .seperator = {'\t'} },
                          { .dict_id = (DictId)dict_id_fields[KRAKEN_TAXID],  .seperator = {'\t'} },
                          { .dict_id = (DictId)dict_id_fields[KRAKEN_SEQLEN], .seperator = {'\t'} },
                          { .dict_id = (DictId)dict_id_fields[KRAKEN_KMERS],                      },
                          { .dict_id = (DictId)dict_id_fields[KRAKEN_EOL],                        } }
    };

    container_seg_by_ctx (vb, &vb->contexts[KRAKEN_TOPLEVEL], (ContainerP)&top_level, 0, 0, 0);

    // top level container when loading a kraken file with --kraken
    SmallContainer top2taxid = { 
        .repeats        = vb->lines.len,
        .is_toplevel    = true,
        .callback       = true,
        .nitems_lo      = 2,
        .items          = { { .dict_id = (DictId)dict_id_fields[KRAKEN_QNAME],  .seperator = { CI_TRANS_NUL /* '\0' */} },
                            { .dict_id = (DictId)dict_id_fields[KRAKEN_TAXID],  .seperator = { CI_TRANS_NOR /* no reconstruct */ } } },
    };

    container_seg_by_ctx (vb, &vb->contexts[KRAKEN_TOP2TAXID], (ContainerP)&top2taxid, 0, 0, 0);
}

// ZIP: called by main thread after compute has finished. 
void kraken_zip_after_compute (VBlockP vb)
{
    // add up the total length of all QNAMEs in the file - will be transferred to PIZ via SectionHeaderCounts.nodes_param of TAXID
    z_file->contexts[KRAKEN_TAXID].nodes.param += vb->last_int(KRAKEN_QNAME); 
}

bool kraken_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return 
        // typically small 
        dict_id.num == dict_id_fields[KRAKEN_CU]       ||
        dict_id.num == dict_id_fields[KRAKEN_QNAME]    || // container
        dict_id.num == dict_id_fields[KRAKEN_TAXID]    || // usually a few thousand species
        dict_id.num == dict_id_fields[KRAKEN_SEQLEN]   ||
        dict_id.num == dict_id_fields[KRAKEN_KMERS]    || // container
        dict_id.num == dict_id_fields[KRAKEN_EOL]      ||
        dict_id.num == dict_id_fields[KRAKEN_TOPLEVEL] ||
        dict_id.num == dict_id_fields[KRAKEN_TOP2TAXID];
}

static void kraken_seg_kmers (VBlock *vb, const char *value, int32_t value_len, const char *taxid, unsigned taxid_len) // must be signed
{
    bool final_space = (value_len && value[value_len-1] == ' '); // sometimes there is, sometimes there isn't...
    value_len -= final_space;

    unsigned num_kmers=1;
    for (unsigned i=0; i < value_len; i++)
        if (value[i] == ' ') {
            if (i && value[i-1] == ' ') { // not valid khmers
                seg_by_did_i (vb, value, value_len, KRAKEN_KMERS, value_len);
                return;
            }
            num_kmers++;
        }

    SmallContainer kmers_con = (SmallContainer){ 
        .repeats   = num_kmers, 
        .nitems_lo = 2,
        .repsep    = {' '},
        .items     = { { .dict_id = (DictId)dict_id_fields[KRAKEN_KMERTAX], .seperator = {':'} },
                       { .dict_id = (DictId)dict_id_fields[KRAKEN_KMERLEN]                     } },
        .drop_final_repeat_sep = !final_space
    };


    const char *kmers[num_kmers+1];
    unsigned kmer_lens[num_kmers+1];
    str_split (value, value_len, num_kmers, ' ', kmers, kmer_lens, "kmers");

    for (unsigned k=0; k < num_kmers; k++) {
        const char *items[3];
        unsigned item_lens[3];

        str_split (kmers[k], kmer_lens[k], 2, ':', items, item_lens, "kmer items");

        // KMER taxid - either copy from TAXID or seg a normal snip
        if (taxid_len == item_lens[0] && !memcmp (taxid, items[0], taxid_len)) 
            seg_by_did_i (vb, copy_taxid_snip, copy_taxid_snip_len, KRAKEN_KMERTAX, item_lens[0]);
        else 
            seg_by_did_i (vb, items[0], item_lens[0], KRAKEN_KMERTAX, item_lens[0]);

        // KMER length
        seg_by_did_i (vb, items[1], item_lens[1], KRAKEN_KMERLEN, item_lens[1]);
    }

    container_seg_by_ctx (vb, &vb->contexts[KRAKEN_KMERS], (ContainerP)&kmers_con, 0, 0, num_kmers*2-1 + final_space); // account for ':' within kmers and ' ' betweem them
}

// example: "C       ST-E00180:535:HCNW2CCX2:8:1101:16183:1221       570     150|150 A:1 570:21 543:5 91347:8 0:81 |:| A:1 543:4 570:25 0:32 28384:5 0:49"
const char *kraken_seg_txt_line (VBlock *vb, const char *field_start_line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    const char *next_field=field_start_line, *field_start;
    unsigned field_len=0;
    char separator;

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    SEG_NEXT_ITEM (KRAKEN_CU);

    // QNAME - We break down the QNAME into subfields separated by / and/or : - these are vendor-defined strings. Examples:
    // Illumina: <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> for example "A00488:61:HMLGNDSXX:4:1101:15374:1031" see here: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
    // PacBio BAM: {movieName}/{holeNumber}/{qStart}_{qEnd} see here: https://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
    GET_NEXT_ITEM ("QNAME");
    SegCompoundArg arg = { .slash = true, .pipe = true, .dot = true, .colon = true };
    seg_compound_field ((VBlockP)vb, &vb->contexts[KRAKEN_QNAME], field_start, field_len, arg, 0, 1 /* \t */);

    vb->last_int(KRAKEN_QNAME) += field_len+1; // count total QNAME lengths in this VB (+1 for separator)

    SEG_NEXT_ITEM (KRAKEN_TAXID);
    const char *taxid = field_start;
    unsigned taxid_len = field_len;

    GET_NEXT_ITEM ("SEQLEN"); // for paired files we will have eg "150|149", for non-paired eg "150"
    if (memchr (field_start, '|', field_len)) { 
        SegCompoundArg arg = { .pipe = true };
        // compound and not array, bc pair_2's seq_len is often a lot more random than pair_1's so use different contexts
        seg_compound_field ((VBlockP)vb, &vb->contexts[KRAKEN_SEQLEN], field_start, field_len, arg, 0, 1 /* \t */); 
    }
    else 
        seg_by_did_i (vb, field_start, field_len, KRAKEN_SEQLEN, field_len+1);
    
    GET_LAST_ITEM ("KMERS");
    kraken_seg_kmers (vb, field_start, field_len, taxid, taxid_len);
    
    SEG_EOL (KRAKEN_EOL, true);

    return next_field;
}

//----------------------------------------------------------------------------------------------------
// PIZ functions when loading as a result of --kraken
// Note: when just pizzing a kraken normally, it is a plain vanilla piz without any special functions
//----------------------------------------------------------------------------------------------------

// returns true if section is to be skipped reading / uncompressing
bool kraken_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    if (flag.kraken_taxid && // note: we only need some of the fields when ingesting with --taxid
        dict_id.num != dict_id_fields[KRAKEN_QNAME]     &&
        dict_id.num != dict_id_fields[KRAKEN_TAXID]     &&
        dict_id.num != dict_id_fields[KRAKEN_EOL]       &&
        dict_id.num != dict_id_fields[KRAKEN_TOP2TAXID] &&
        dict_id_typeless (dict_id).id[0] != 'Q') // array item and compound items of QNAME  
        return true;

    return false;
}

void kraken_piz_handover_data (VBlockP vb)
{
    uint64_t dict_start = qname_dict.len;
    buf_add_buf (evb, &qname_dict, &vb->txt_data, char, "qname_dict");

    ARRAY (QnameNode, nodes, vb->contexts[KRAKEN_QNAME].nodes);
    for (uint64_t i=0; i < nodes_len; i++)
        nodes[i].char_index += dict_start;

    buf_add_buf (evb, &qname_nodes, &vb->contexts[KRAKEN_QNAME].nodes, QnameNode, "qname_nodes");
}

// callback called after every repeat of TOP2HASH, i.e. when run with --taxid
// inspect QNAME and decide if to include in the qname_nodes or not
CONTAINER_CALLBACK (kraken_piz_container_cb)
{
    // when loading kraken to filter another file, remove lines according to the requested filter
    if (dict_id.num == dict_id_fields[KRAKEN_TOP2TAXID]) {

        uint32_t this_taxid = vb->last_int (KRAKEN_TAXID);

        // case: line is filtered out in genocat (in genozip - all lines are loaded)
        if (flag.kraken_taxid >= 0 && // it will be -1 when "genozip --kraken" as it is --taxid is unavailable for genozip
            (    (qname_type == (QNAME_INCLUDE != flag.kraken_taxid_negative) && flag.kraken_taxid != this_taxid)
              || (qname_type == (QNAME_EXCLUDE != flag.kraken_taxid_negative) && flag.kraken_taxid == this_taxid)))
            vb->drop_curr_line = true;

        // case line is included
        else {
            buf_alloc (vb, &vb->contexts[KRAKEN_QNAME].nodes, 1, vb->lines.len, QnameNode, 1.5, "contexts->nodes");
            
            NEXTENT (QnameNode, vb->contexts[KRAKEN_QNAME].nodes) = (QnameNode){ 
                .char_index = reconstructed - vb->txt_data.data, // will be updated at merge
                .snip_len   = reconstructed_len - 1,             // excluding the \0 seperator
                .hash       = hash_do (qname_hashtab.len, reconstructed, reconstructed_len-1),
                .taxid      = this_taxid
            };
        } 
    }

    // when pizzing a kraken to be filtered by itself (loaded or stored) (useful only for testing) - apply :kraken_is_included_loaded" filter 
    else if (flag.kraken_taxid && 
             dict_id.num == dict_id_fields[KRAKEN_TOPLEVEL] && 
             (   ( kraken_is_loaded && !kraken_is_included_loaded (vb, last_txt (vb, KRAKEN_QNAME), vb->last_txt_len (KRAKEN_QNAME)))
              || (!kraken_is_loaded && !kraken_is_included_stored (vb, KRAKEN_TAXID, true)))) // TAXID was reconstructed in the TOPLEVEL container
        
        vb->drop_curr_line = true;
}

// Loading kraken: called after global area is read (including SEC_COUNTS), before reconstruction of VBs
void kraken_piz_initialize (void)
{
    if (!flag.reading_kraken) return;

    ASSERTNOTINUSE (qname_hashtab);
    ASSERTNOTINUSE (qname_nodes);
    ASSERTNOTINUSE (qname_dict);
    
    Context *ctx = &z_file->contexts[KRAKEN_TAXID];
    ARRAY (int64_t, counts, ctx->counts);

    // calculate the total number of sequences in the kraken file (and the FASTQ file(s) from which it was generated)
    int64_t total_sequences_in_file=0; 
    for (uint64_t i=0; i < counts_len; i++)
        total_sequences_in_file += counts[i];

    ASSERT0 (total_sequences_in_file, "unexpectedly, total_sequences_in_file=0");

    // calculate the average length of a QNAME string in the file
    int64_t sequences_in_list, 
            total_qname_length = ctx->nodes.param; // sent ZIP->PIZ via SectionHeaderCounts.node_param
    one_qname_length = (unsigned)(1 + total_qname_length / total_sequences_in_file); // average (rounded up) QNAME length in the file
    
    // get info of the user selected taxonomy ID 
    if (exe_type == EXE_GENOCAT) {

        ASSINP0 (flag.kraken_taxid != -1, "--taxid must be provided if --kraken is used");

        WordIndex taxid_word_i = ctx_get_word_index_by_snip (ctx, str_int_s (flag.kraken_taxid).s);
        RETURNW (taxid_word_i != WORD_INDEX_NONE,, "FYI: %s has no sequences with a Taxanomic ID of \"%d\"", z_name, flag.kraken_taxid);

        int64_t taxid_count = counts[taxid_word_i]; // number of QNAMEs in the file with this taxid (note: in SAM, usually there are multiple lines with the same QNAME)
    
        // case: we have a single component (the common case): decide whether qname_nodes will be the 
        // included or excluded qnames based on the smaller of the two - to save memory
        if (z_file->num_components == 1)
            qname_type = taxid_count < total_sequences_in_file/2 ? QNAME_INCLUDE : QNAME_EXCLUDE;     

        // case: we have multiple components (eg components generated by krakening R1 and R2 FASTQs).
        // The components might have a different classification for the same QNAME.
        // In a positive --taxid - we always use QNAME_INCLUDE it is sufficient that one of the components classifies 
        // as the requested taxid.
        // In a negative --taxid (eg ^9606) - we always uses QNAME_EXCLUDE as it is sufficient that one of the components
        // classifies as the undesired taxid (9606 in this example), that it will be excluded
        else
            qname_type = flag.kraken_taxid_negative ? QNAME_EXCLUDE : QNAME_INCLUDE;
            
        sequences_in_list = (qname_type == QNAME_INCLUDE) ? taxid_count : (total_sequences_in_file - taxid_count);
    }
    else { // genozip - we load everything
        qname_type = QNAME_INCLUDE;
        sequences_in_list = total_sequences_in_file;
    }

    // allocate
    TEMP_FLAG (quiet, true); // we're allocating a huge amount of memory (~30GB for a 40x BAM) - temporarily suppress buf_alloc warning
    buf_alloc (evb, &qname_nodes, 0, sequences_in_list, QnameNode, 0, "qname_nodes"); 
    buf_alloc (evb, &qname_dict,  0, sequences_in_list, char[one_qname_length+2],  0, "qname_dict"); // approximate (+2), might grow if needed in kraken_piz_handover_data
    RESTORE_FLAG (quiet);

    qname_hashtab.len = hash_next_size_up (3 * sequences_in_list, true);

    //printf ("total_sequences_in_file=%u sequences_in_list=%u qname_hashtab.len=%u one_qname_length=%u\n", (int)total_sequences_in_file, (int)sequences_in_list, (int)qname_hashtab.len, (int)one_qname_length);    

    buf_alloc (evb, &qname_hashtab, 0, qname_hashtab.len, uint32_t, 1, "qname_hashtab");
    buf_set (&qname_hashtab, 0xff); // empty
}

static int kraken_qname_nodes_cmp (const void *a, const void *b)
{
    return ((QnameNode *)a)->hash - ((QnameNode *)b)->hash;
}

// genocat: load kraken file as a result of genocat --kraken
void kraken_load (void)
{
    ASSERTNOTNULL (flag.reading_kraken);
    SAVE_FLAGS_AUX;

    flag.may_drop_lines   = true;    // we drop the lines not matching --taxid
    flag.trans_containers = true;    // execute TOP2TAXID in translate mode
    flag.data_modified    = true;    // the reconstructed kraken file is not the same as the original...
    flag.out_dt           = DT_NONE; // needed for dt_get_translation to find the translation defined in TRANSLATIONS
    flag.reference        = REF_NONE;
    TEMP_VALUE (command, PIZ);

    z_file = file_open (flag.reading_kraken, READ, Z_FILE, DT_KRAKEN);    
    zfile_read_genozip_header (0,0,0);
    
    ASSINP (z_file->data_type == DT_KRAKEN, "expected %s to be a genozip'ed kraken output file, but its a %s file. Tip: compress the output file generated by kraken2 with \"genozip --input kraken\"", 
            z_name, dt_name (z_file->data_type));

    z_file->basename = file_basename (flag.reading_kraken, false, "(kraken-file)", NULL, 0);

    Dispatcher dispachter = piz_z_file_initialize (false);
    piz_one_txt_file (dispachter, false);

    kraken_filename = file_make_unix_filename (z_name); // full-path unix-style filename, allocates memory

    file_close (&z_file, false, false);
    file_close (&txt_file, false, false); // close the txt_file object we created (even though we didn't open the physical file). it was created in file_open called from txtheader_piz_read_and_reconstruct.
    
    qsort (qname_nodes.data, qname_nodes.len, sizeof (QnameNode), kraken_qname_nodes_cmp);

    ASSERT (qname_nodes.len <= MAX_QNAME_NODES, "qname_nodes.len=%"PRIu64" exceeds the maximum of %u", qname_nodes.len, MAX_QNAME_NODES);

    // build the hash table
    for (uint32_t i=0; i < (uint32_t)qname_nodes.len; i++) {
        uint32_t hash = ENT (QnameNode, qname_nodes, i)->hash;
        *ENT (uint32_t, qname_hashtab, hash) = i; // overwriting an empty or populated hash entry 
    }

    // recover globals
    RESTORE_VALUE (command);
    RESTORE_FLAGS;

    flag.reading_kraken = NULL;
}

// -------------------------------------------
// using the kraken data in genocat --kraken
// -------------------------------------------

void kraken_set_taxid (const char *optarg)
{
    flag.kraken_taxid_negative = optarg[0] == '^';

    str_get_int_range32 (optarg + flag.kraken_taxid_negative, 0, 0, 9999999, &flag.kraken_taxid); // NCBI Taxonomic ID is a 7 digit positive integer (and "0" means unclassified in kraken)
}

static TaxonomyId kraken_get_taxid (const char *qname, unsigned qname_len)
{
    uint32_t hash = hash_do (qname_hashtab.len, qname, qname_len);

    uint32_t list_index = *ENT (uint32_t, qname_hashtab, hash);
    if (list_index == QNAME_NODE_NONE) 
        return TAXID_NONE;

    ASSERT (list_index < qname_nodes.len, "expecting list_index=%u < qname_nodes.len=%"PRIu64, list_index, qname_nodes.len);

    // hash found, it could be the requested qname, or a different one with the same hash. let's check
    QnameNode *listent = ENT (QnameNode, qname_nodes, list_index);
    const char *snip = ENT (char, qname_dict, listent->char_index);

    if (listent->snip_len == qname_len && !memcmp (qname, snip, qname_len))
        return listent->taxid;

    // we have a hash table contention - the hash file points to a different QNAME with the same hash. 
    // we look for all other entries with the same hash. these are adajent to listent as qname_nodes is sorted by hash

    // search down
    for (QnameNode *l=listent-1; l >= FIRSTENT (QnameNode, qname_nodes) && l->hash == hash; l--)
        if (l->snip_len == qname_len && !memcmp (qname, ENT (char, qname_dict, l->char_index), qname_len))
            return l->taxid;

    // search up
    for (QnameNode *l=listent+1; l <= LASTENT (QnameNode, qname_nodes) && l->hash == hash; l++)
        if (l->snip_len == qname_len && !memcmp (qname, ENT (char, qname_dict, l->char_index), qname_len))
            return l->taxid;

    return TAXID_NONE;
}

// Find whether QNAME is included in the taxid filter in O(1) (using the loaded kraken file)
bool kraken_is_included_loaded (VBlockP vb, const char *qname, unsigned qname_len)
{
    ASSINP0 (flag.kraken_taxid != -1, "When using --kraken, you must specify --taxid <number> (positive filter) or --taxid ^<number> (negative filter)");
    ASSINP0 (!z_file->z_flags.has_taxid || vb->data_type == DT_KRAKEN, "You cannot use --kraken, because the file already contains taxonomic information (it was compressed with --kraken) - just use --taxid");

    if (!qname_hashtab.len) {
        WARN_ONCE ("FYI: %s contains no sequence with taxid=%d", kraken_filename, flag.kraken_taxid);
        return false; // this happens with the kraken file contains no sequences with the requested taxid
    }

    bool res;
    #define INCLUDE_LINE_IF(x) res = (qname_type == (x))

    if (kraken_get_taxid (qname, qname_len) == TAXID_NONE) 
        INCLUDE_LINE_IF (QNAME_EXCLUDE);
    else 
        INCLUDE_LINE_IF (QNAME_INCLUDE);  // qname found in dict

    if (flag.show_kraken)
        iprintf ("%.*s\t%s\n", qname_len, qname, res ? "INCLUDED" : "EXCLUDED");

    return res;
}

// Find whether QNAME is included in the taxid filter in O(1) (using the loaded kraken file)
bool kraken_is_included_stored (VBlockP vb, DidIType did_i_taxid, bool already_reconstructed)
{
    if (flag.kraken_taxid == -1) return true; // use didn't specify --taxid - everything's included

    ASSINP0 (z_file->z_flags.has_taxid, "To use --taxid, either use \"genocat --kraken <file> --taxid <taxid>\" to load a kraken file, or use \"genozip --kraken <file>\" to include taxonomic information in the genozip file");

    if (!already_reconstructed) 
        reconstruct_from_ctx (vb, did_i_taxid, 0, false); // update last_int without reconstructing
    
    uint32_t this_taxid = vb->last_int (did_i_taxid);

    return  (( flag.kraken_taxid_negative && flag.kraken_taxid != this_taxid) ||
             (!flag.kraken_taxid_negative && flag.kraken_taxid == this_taxid));
}

//----------------------------------------------------------------------------------------------------
// SEG helper functions for other data types segging with genozip --kraken
//----------------------------------------------------------------------------------------------------

// returns snip_len if successful, or 0 if failed (only happens if fail_if_missing)
unsigned kraken_seg_taxid_do (VBlockP vb, DidIType did_i_taxid, const char *qname, unsigned qname_len, 
                              char *snip, // caller-allocated out
                              bool fail_if_missing)
{
    TaxonomyId taxid = kraken_get_taxid (qname, qname_len);

    ASSINP (!fail_if_missing || taxid != TAXID_NONE, 
            "Cannot find taxonomy id in kraken file for QNAME \"%.*s\"", qname_len, qname);

    if (taxid == TAXID_NONE) 
        return 0; // failed

    else {
        unsigned snip_len = str_int (taxid, snip);

        if (flag.show_kraken)
            iprintf ("%.*s\t%d\n", qname_len, qname, taxid);

        seg_by_did_i (vb, snip, snip_len, did_i_taxid, 0);

        return snip_len;
    }
}

// always segs if even_if_missing ; or if even_if_missing=false, returns true if taxid found and segged
// returns snip_len if successful, or 0 if failed (only happens if fail_if_missing)
unsigned kraken_seg_taxid (VBlockP vb, DidIType did_i_taxid, 
                           const char *qname, unsigned qname_len, 
                           bool fail_if_missing)
{
    char snip[20];
    return kraken_seg_taxid_do (vb, did_i_taxid, qname, qname_len, snip, fail_if_missing);
}
