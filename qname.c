// ------------------------------------------------------------------
//   qname.c
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "qname.h"
#include "tokenizer.h"
#include "buffer.h"
#include "vblock.h"
#include "segconf.h"
#include "strings.h"
#include "seg.h"
#include "sam.h"
#include "fastq.h"
#include "stats.h"
#include "qname_cons.h"
#include "file.h"

#define MAX_QNAME_ITEMS 11 // mate + 10 others (matching Q?NAME defined in sam.h, fastq.h, kraken.h) 

static STRl(copy_qname, 50);

//----------------------------------------------------------------

typedef struct { STR(s); int32_t fixed_i; } FixedStr;

typedef struct QnameFlavorStruct {
    char name[16], example[256];
    SeqTech tech;                              // The sequencing technology used to generate this data
    int fq_only;                               // this QF is only for FASTQ, not SAM/BAM/KRAKEN. FQ: 1 if /1 goes after the first string, 2 if after the second space-separated string
    SmallContainer *con_template;              // container template - copied to con in qname_zip_initialize
    unsigned num_seps;                         // number of printable separator and prefix characters
    int integer_items[MAX_QNAME_ITEMS+1];      // array of item_i (0-based) that are expected to be integer - no leading zeros (+1 for termiantor; terminated by -1)
    int numeric_items[MAX_QNAME_ITEMS+1];      // array of item_i (0-based) that are expected to be numeric - leading zeros ok (+1 for termiantor; terminated by -1)
    int in_local[MAX_QNAME_ITEMS+1];           // int or numeric items that should be stored in local. if not set, item is segged to dict.
    int hex_items[MAX_QNAME_ITEMS+1];          // array of item_i (0-based) that are expected to be hexademical (+1 for termiantor; terminated by -1)
    int ordered_item1, ordered_item2;          // an item that may be delta'd in non-sorted files. this should be the fast-moving item which consecutive values are close
    int range_end_item;                        // item containing end of range - delta vs preceding item in container
    unsigned fixed_len;                        // fixed length (0 if it is not fixed length)
    FixedStr px_strs[MAX_QNAME_ITEMS+1];       // fixed substrings in template (=prefixes of items in fixed locations). (+1 for termiantor; terminated by empty entry) 
    int qname2;                                // embedded qname2 item (generated in qname_zip_initialize)
    STRl (con_snip,  200);                     // container snips (generated in qname_zip_initialize)
    STRl (con_snip2, 200);                     // same as con_snip, just with QNAME2 dict_ids
    #define MAX_PREFIX_LEN 30
    STRl (con_prefix, MAX_PREFIX_LEN);         // prefix of container
    SmallContainer con;                        // container
} QnameFlavorStruct;

static QnameFlavorStruct qf[] = { 
/*  mate   name             example                                   tech     fq_only   con_template     #sp int_items       numeric_items   in-local        hex_items       ord1,2 rng len px_strs{str,str_len,fixed_i}                               */
         { "Illumina-fastq","A00488:61:HMLGNDSXX:4:1101:4345:1000 2:N:0:CTGAAGCT+ATAGAGGC",
                                                                      TECH_ILLUM_7, 1, &con_illumina_7_fq, 7, {1,3,4,5,6,-1}, {-1},           {-1},           {-1},           5,6,   -1,                                                                 },
    {},  { "Illumina",      "A00488:61:HMLGNDSXX:4:1101:4345:1000",   TECH_ILLUM_7, 0, &con_illumina_7,    6, {1,3,4,5,6,-1}, {-1},           {-1},           {-1},           5,6,   -1,                                                                 },
    {},  { "BGI-E",         "E100020409L1C001R0030000801",            TECH_BGI,     0, &con_bgi_E,         4, {-1},           {0,1,2,3,4,-1}, {-1},           {-1},           4,-1,  -1, 27, { {"E",1,0},  {"L",1,10},{"C",1,12},{"R",1,16},{"",0,20}, {(char[]){CI0_SKIP},0,27} } },
    {},  { "BGI-CL",        "CL100025298L1C002R050_244547",           TECH_BGI,     0, &con_bgi_CL,        5, {-1},           {0,1,2,3,4,-1}, {-1},           {-1},           4,-1,  -1, 0,  { {"CL",2,0}, {"L",1,11},{"C",1,13},{"R",1,17},{"_",1,21} } },
    {},  { "IonTorrent",    "ZEWTM:10130:07001",                      TECH_IONTORR, 0, &con_ion_torrent_3, 2, {-1},           {1,2,-1},       {-1},           {-1},           -1,-1, -1, 17                                                              },
    {},  { "Illumina-old#", "HWI-ST550_0201:3:1101:1626:2216#ACAGTG", TECH_ILLUM_5, 0, &con_illumina_5i,   5, {1,2,3,4,-1},   {-1},           {-1},           {-1},           -1,-1, -1,                                                                 },
    {},  { "Illumina-old",  "SOLEXA-1GA-1_4_FC20ENL:7:258:737:870",   TECH_ILLUM_5, 0, &con_illumina_5,    4, {1,2,3,4,-1},   {-1},           {-1},           {-1},           -1,-1, -1,                                                                 },
    {},  { "Roche-454",     "000050_1712_0767",                       TECH_454,     0, &con_roche_454,     2, {-1},           {0,1,2,-1},     {-1},           {-1},           -1,-1, -1, 16, { {"",0,0},{"_",1,6},{"_",1,11},{(char[]){CI0_SKIP},0,16} } },
    {},  { "Helicos",       "VHE-242383071011-15-1-0-2",              TECH_HELICOS, 0, &con_helicos,       5, {2,3,4,5,-1},   {1,-1},         {-1},           {-1},           -1,-1, -1,                                                                 },
    {},  { "PacBio-3",      "56cdb76f_70722_4787",                    TECH_PACBIO,  0, &con_pacbio_3,      2, {1,2,-1},       {-1},           {-1},           {0,-1},         -1,-1, -1,                                                                 },
    {},  { "PacBio-Range",  "m130802_221257_00127_c100560082550000001823094812221334_s1_p0/128361/872_4288",
                                                                      TECH_PACBIO,  0, &con_pacbio_range,  4, {1,2,3,-1},     {-1},           {-1},           {-1},           1,-1,   3, 0,  { { "m",1,0} }                                              },
    {},  { "PacBio-Label",  "m64136_200621_234916/18/ccs",            TECH_PACBIO,  0, &con_pacbio_label,  3, {1,-1},         {-1},           {-1},           {-1},           1,-1,  -1, 0,  { { "m",1,0} }                                              },
    {},  { "PacBio-Plain",  "m64136_200621_234916/18",                TECH_PACBIO,  0, &con_pacbio_plain,  2, {1,-1},         {-1},           {-1},           {-1},           1,-1,  -1, 0,  { { "m",1,0} }                                              },
    {},  { "Nanopore",      "af84b0c1-6945-4323-9193-d9f6f2c38f9a",   TECH_ONP,     0, &con_nanopore,      4, {-1},           {-1},           {0,1,2,3,4,-1}, {0,1,2,3,4,-1}, -1,-1, -1, 36,                                                             },
    {},  { "Nanopore-ext",  "2a228edf-d8bc-45d4-9c96-3d613b8530dc_Basecall_2D_000_template",                               
                                                                      TECH_ONP,     0, &con_nanopore_ext,  5, {-1},           {-1},           {0,1,2,3,4,-1}, {0,1,2,3,4,-1}, -1,-1, -1,                                                                 },
    {},  { "NCBI-SRA2+-FQ", "ERR2708427.1.1 51e7525d-fa50-4b1a-ad6d-4f4ae25c1df7 someextradata length=1128",                         
                                                                      TECH_UNKNOWN, 2, &con_ncbi_sra2P_fq, 6, {2,3,8,-1},     {1, -1},        {-1},           {-1},           3,-1,  -1,                                                                 },
    {},  { "NCBI-SRA+-FQ",  "ERR1111170.1 07dc4948-eb0c-45f2-9b40-a933a9bd5cf7_Basecall_2D_000_template BOWDEN04_20151016_MN15199_FAA67113_BOWDEN04_MdC_MARC_Phase2a_4833_1_ch19_file1_strand length=52",
                                                                      TECH_UNKNOWN, 2, &con_ncbi_sraP_fq,  5, {2,7,-1},       {1, -1},        {-1},           {-1},           2,-1,  -1,                                                                 },
    {},  { "NCBI-SRA2-FQ",  "ERR2708427.1.1 51e7525d-fa50-4b1a-ad6d-4f4ae25c1df7 length=1128",                         
                                                                      TECH_UNKNOWN, 2, &con_ncbi_sra2_fq,  5, {2,3,7,-1},     {1, -1},        {-1},           {-1},           3,-1,  -1,                                                                 },
    {},  { "NCBI-SRA-FQ",   "ERR1111170.1 07dc4948-eb0c-45f2-9b40-a933a9bd5cf7_Basecall_2D_000_template length=52",
                                                                      TECH_UNKNOWN, 2, &con_ncbi_sra_fq,   4, {2,6,-1},       {1, -1},        {-1},           {-1},           2,-1,  -1,                                                                 },
    {},  { "NCBI-SRA2",     "ERR2708427.1.1",                         TECH_UNKNOWN, 0, &con_ncbi_sra2,     2, {2,3,-1},       {1, -1},        {-1},           {-1},           3,-1,  -1,                                                                 },
    {},  { "NCBI-SRA",      "SRR001666.1",                            TECH_UNKNOWN, 0, &con_ncbi_sra,      1, {2,-1},         {1, -1},        {-1},           {-1},           2,-1,  -1,                                                                 },
    {},  { "Genozip-opt",   "basic.1",  /* must be last */            TECH_UNKNOWN, 0, &con_genozip_opt,   1, {1,-1},         {-1},           {-1},           {-1},           1,-1,  -1,                                                                 },
};

#define NUM_QFs (sizeof(qf)/sizeof(qf[0]))

static inline unsigned qname_get_px_str_len (QnameFlavorStruct *qfs)
{
    unsigned i=0; 
    while (qfs->px_strs[i].s) i++;
    return i;
}

static void qname_remove_skip (SmallContainerP con_no_skip, uint32_t *prefix_no_skip_len, 
                               ConstSmallContainerP con, STRp (con_prefix)) // out
{
    // all container data, except for items
    memcpy (con_no_skip, con, sizeof (MiniContainer) - sizeof (ContainerItem));
    con_no_skip->nitems_lo = 0;

    // copy items, except SKIP items    
    for (int item_i=0; item_i < con->nitems_lo; item_i++)
        if (con->items[item_i].separator[0] != CI0_SKIP) 
            con_no_skip->items[con_no_skip->nitems_lo++] = con->items[item_i];

    // copy prefix, except last item if it is a SKIP (we only support SKIP for final prefix items)
    const char *px_skip = con_prefix ? memchr (con_prefix, CI0_SKIP, con_prefix_len) : 0;
    *prefix_no_skip_len = px_skip ? (px_skip - con_prefix) : *prefix_no_skip_len;
}

static void qname_genarate_qfs_with_mate (QnameFlavorStruct *qfs)
{
    // find mate item (can't be item 0)
    int mate_item_i;
    for (mate_item_i=1; mate_item_i < qfs->con.nitems_lo; mate_item_i++)
        if (qfs->con.items[mate_item_i].dict_id.num == _SAM_QmNAME) break;

    ASSERT (mate_item_i < qfs->con.nitems_lo, "Can't find mate item in QFS=%s", qfs->name);

    // add new item as numeric item (at beginning of array - easier, and the order doesn't matter)
    memmove (&qfs->numeric_items[1], qfs->numeric_items, sizeof(qfs->numeric_items) - sizeof(qfs->numeric_items[0])); // make room
    qfs->numeric_items[0] = mate_item_i;

    // name
    strcpy (&qfs->name[strlen(qfs->name)], "/");
    
    // example
    unsigned example_len = strlen(qfs->example);
    char *after = &qfs->example[example_len];
    if (!qfs->fq_only)
        strcpy (after, "/1");
    else {
        str_split (qfs->example, example_len, 0, ' ', str, false);
        char *slash = (char *)strs[qfs->fq_only-1] + str_lens[qfs->fq_only-1];
        memmove (slash+2, slash, after - slash);
        slash[0] = '/';
        slash[1] = '1';
    }

    // case 1: previous item is not fixed - move its separator (possibly 0) to the mate item and make it '/'
    if (qfs->con.items[mate_item_i-1].separator[0] != CI0_FIXED_0_PAD) {
        qfs->con.items[mate_item_i].separator[0] = qfs->con.items[mate_item_i-1].separator[0]; // 0 or ' '
        qfs->con.items[mate_item_i-1].separator[0] = '/';
        qfs->num_seps++;
    }

    // case 2: previous item is fixed - add separator as a prefix. 
    else { // eg con_bgi_E, con_roche_454
        // qfs must already have prefix item for the mate item - initially empty string
        ASSERT (qname_get_px_str_len (qfs) > mate_item_i, "QnameFlavor=%s: expecting prefix to exist for mate_item_i=%u",
                qfs->name, mate_item_i);

        qfs->px_strs[mate_item_i].s = "/";
        qfs->px_strs[mate_item_i].s_len = 1;
    }

    if (qfs->fixed_len) // note: qfs can have fixed_len even if items aren't fixed (eg IonTorrent)
        qfs->fixed_len += 2;
}

// note: we need to re-initialize in each file, lest the data type has changed
void qname_zip_initialize (DidIType qname_did_i)
{
    ContextP qname_zctx = ZCTX(qname_did_i);

    // we need to prepare the containers only once and they can serve all files, as the containers
    // are data-type independent (since all data types use the dict_id Q?NAME for qnames)
    DO_ONCE {
        for (QnameFlavorStruct *qfs=&qf[0]; qfs < &qf[NUM_QFs]; qfs++) {

            // case: this qfs os version of the next qfs, just with room for /1 /2 - generate it now
            if (!qfs->name[0]) { 
                *qfs = *(qfs+1);
                qfs->con = *qfs->con_template;
                qname_genarate_qfs_with_mate (qfs);
            }
            else
                qfs->con = *qfs->con_template; // create container

            qfs->qname2 = -1; // initialize pessimistically
            ContextP next_zctx = qname_zctx + 1;
            for (int item_i=0; item_i < qfs->con.nitems_lo; item_i++) {
                uint64_t dnum = qfs->con.items[item_i].dict_id.num;
                // set qname2
                if (dnum == _FASTQ_QNAME2) {
                    ASSERT0 (qfs->fq_only && qfs->tech==TECH_UNKNOWN, "Bad embedded qname2 definition"); // a con with FASTQ_QNAME2 must have these properties
                    qfs->qname2 = item_i;
                }
                // verify dict_id
                else if (dnum != _SAM_QmNAME) {
                    ASSERT (next_zctx->dict_id.num == dnum, "Expecting item #%u of %s to have to be %s", item_i, qfs->name, next_zctx->tag_name);
                    next_zctx++;
                }
            }

            // case: this QF needs a container prefix - generate it
            if (qfs->px_strs[0].s) {
                qfs->con_prefix[qfs->con_prefix_len++] = CON_PX_SEP; 
                qfs->con_prefix[qfs->con_prefix_len++] = CON_PX_SEP; // no container-wide prefix

                for (const FixedStr *px = qfs->px_strs; px->s ; px++) {
                    memcpy (&qfs->con_prefix[qfs->con_prefix_len], px->s, px->s_len);
                    qfs->con_prefix_len += px->s_len;
                    qfs->con_prefix[qfs->con_prefix_len++] = CON_PX_SEP; 
                }
            } 

            // remove SKIP item, if one exists, from container and prefix
            SmallContainer con_no_skip = qfs->con; // copy
            uint32_t prefix_no_skip_len=qfs->con_prefix_len;
            qname_remove_skip (&con_no_skip, &prefix_no_skip_len, &qfs->con, STRa (qfs->con_prefix));
            
            // prepare container snip for this QF
            qfs->con_snip_len = sizeof (qfs->con_snip);            
            container_prepare_snip ((Container*)&con_no_skip, qfs->con_prefix, prefix_no_skip_len, qfs->con_snip, &qfs->con_snip_len);

            // prepare snip of identical container, except dict_ids are of QNAME2
            DidIType next_qname2_did_i = FASTQ_Q0NAME2;
            for (unsigned item_i=0; item_i < con_no_skip.nitems_lo; item_i++) 
                if (con_no_skip.items[item_i].dict_id.num != _SAM_QmNAME && con_no_skip.items[item_i].dict_id.num != _FASTQ_QNAME2) {
                    con_no_skip.items[item_i].dict_id = dt_fields[DT_FASTQ].predefined[next_qname2_did_i].dict_id; // only FASTQ has QNAME2 dict_ids, this is not nessarilly the current data type
                    next_qname2_did_i++;
                }

            qfs->con_snip2_len = sizeof (qfs->con_snip2);            
            container_prepare_snip ((Container*)&con_no_skip, qfs->con_prefix, prefix_no_skip_len, qfs->con_snip2, &qfs->con_snip2_len);
        }
    }

    // prepare copy_qname snip - every time, as container_did_i changes between data types
    seg_prepare_snip_other (SNIP_COPY, qname_zctx->dict_id, false, 0, copy_qname); // QNAME dict_id is the same for SAM, FASTQ, KRAKEN
}

static void qname_seg_initialize_do (VBlockP vb, QnameFlavor qfs, DidIType qname_did_i, DidIType st_did_i, unsigned num_qname_items)
{
    if (!qfs) return;

    // consolidate all items under QNAME in --stats
    ContextP qname_ctxs[num_qname_items+1];
    for (int i=0; i < num_qname_items+1; i++) {
        qname_ctxs[i] = CTX(qname_did_i + i);

        if (flag.pair) {
            qname_ctxs[i]->no_stons = true; // prevent singletons, so pair_1 and pair_2 are comparable based on b250 only
            
            if (flag.pair == PAIR_READ_2)
                ctx_create_node (vb, qname_did_i + i, (char[]){ SNIP_MATE_LOOKUP }, 1); // required by ctx_convert_generated_b250_to_mate_lookup
        }
    }

    stats_set_consolidation_(vb, st_did_i, num_qname_items+1, qname_ctxs); 

    // set STORE_INT as appropriate
    if (qfs->ordered_item1  != -1) CTX(qname_did_i + 1 + qfs->ordered_item1)->flags.store = STORE_INT; 
    if (qfs->ordered_item2  != -1) CTX(qname_did_i + 1 + qfs->ordered_item2)->flags.store = STORE_INT; 
    if (qfs->range_end_item != -1) CTX(qname_did_i + 1 + qfs->range_end_item - 1)->flags.store = STORE_INT; 

    if (qfs->qname2 != -1)
        qname_seg_initialize_do (vb, segconf.qname_flavor2, FASTQ_QNAME2, st_did_i, num_qname_items-1/*-1 bc no mate*/);
}

void qname_seg_initialize (VBlockP vb, DidIType qname_did_i) 
{   
    qname_seg_initialize_do (vb, segconf.qname_flavor, qname_did_i, qname_did_i, MAX_QNAME_ITEMS); 
}

static inline bool qname_has_px_strs (STRp(qname), QnameFlavor qfs)
{
    for (const FixedStr *px = qfs->px_strs; px->s ; px++) {
        int i = (px->fixed_i >= 0) ? px->fixed_i : ((int)qname_len + px->fixed_i);
        if (i + px->s_len > qname_len ||          // case: fixed str goes beyond the end of qname
            memcmp (px->s, &qname[i], px->s_len)) // case: fixed str does not appear in qname
            return false;
    }

    return true; // all fixed strings appear as expected
}

// note: we run this function only in discovery, not in segging, because it is quite expesive - checking all numerics.
bool qname_is_flavor (STRp(qname), QnameFlavor qfs,
                      pSTRp (qname2)) // out
{
    if (!qname_len) return false;

    // test fixed length, if applicable
    if (qfs->fixed_len && qname_len != qfs->fixed_len) return false;

    // check that the fixed strings expected by the flavor appear in qname
    if (!qname_has_px_strs (STRa(qname), qfs)) return false;

    str_split_by_container (qname, qname_len, &qfs->con, qfs->con_prefix, qfs->con_prefix_len, item);
    if (!n_items) return false; // failed to split according to container - qname is not of this flavor

    // items cannot include whitespace or non-printable chars
    for (uint32_t item_i=0; item_i < n_items; item_i++)
        if (!str_is_no_ws (STRi(item, item_i))) return false;

    // check that all the items expected to be numeric (leading zeros ok) are indeed so
    for (const int *item_i = qfs->numeric_items; *item_i != -1; item_i++)
        if (!str_is_numeric (STRi(item, *item_i))) return false;

    // check that all the items expected to be hexadecimal are indeed so
    for (const int *item_i = qfs->hex_items; *item_i != -1; item_i++)
        if (!str_is_hexlo (STRi(item, *item_i))) return false;

    // check that all the items expected to be integer (no leading zeros) are indeed so
    for (const int *item_i = qfs->integer_items; *item_i != -1; item_i++)
        if (!str_is_int (STRi(item, *item_i))) return false;

    // return embedded qname2, eg "82a6ce63-eb2d-4812-ad19-136092a95f3d" in "@ERR3278978.1 82a6ce63-eb2d-4812-ad19-136092a95f3d/1"
    if (qname2 && qfs->qname2 != -1) {
        *qname2 = items[qfs->qname2];
        *qname2_len = item_lens[qfs->qname2];
    }

    return true; // yes, qname is of this flavor
}

// called for the first line in segconf.running
void qname_segconf_discover_flavor (VBlockP vb, DidIType qname_did_i, STRp(qname))
{
    segconf.qname_flavor = 0; // unknown
    STR0(qname2);

    for (QnameFlavor qfs=&qf[0]; qfs < &qf[NUM_QFs]; qfs++) 
        if ((VB_DT(DT_FASTQ) || !qfs->fq_only) && qname_is_flavor (STRa(qname), qfs, pSTRa(qname2))) {
            segconf.qname_flavor = qfs;
            segconf.tech = qfs->tech;
            qname_seg_initialize (vb, qname_did_i); // so the rest of segconf.running can seg fast using the discovered container

            // if it has an embedded qname2 - find the true tech from qname2 (only possible for FASTQ, in SAM we can't find the tech from NCBI qname)
            if (qname2) 
                for (QnameFlavor qf2=&qf[0]; qf2 < &qf[NUM_QFs]; qf2++) 
                    if (!qf2->fq_only && qname_is_flavor (STRa (qname2), qf2, 0, 0)) { // qname2 cannot be "fq_only"
                        segconf.qname_flavor2 = qf2;
                        segconf.tech = qf2->tech;
                    }
            break;
        }

    // when optimizing qname with --optimize_DESC - capture the correct TECH ^ but set flavor to Genozip-opt
    if (flag.optimize_DESC) {
        segconf.qname_flavor = &qf[ARRAY_LEN(qf)-1]; // Genozip-opt is last
        segconf.qname_flavor2 = NULL;
    }

    if (flag.debug) // unit test
        for (QnameFlavor qfs=&qf[0]; qfs < &qf[NUM_QFs]; qfs++) 
            ASSERT (qname_is_flavor (qfs->example, strlen (qfs->example), qfs, 0, 0), 
                    "Failed to identify qname \"%s\" as %s", qfs->example, qfs->name);
}

// attempt to seg according to the qf - return true if successful
static inline bool qname_seg_qf (VBlockP vb, ContextP qname_ctx, QnameFlavor qfs, STRp(qname), unsigned add_additional_bytes)
{
    // check that the fixed strings expected by the flavor appear in qname. note: we don't verify numeric fields as
    // it would be expensive. instead, we fallback on segging them as text if they turn out to be not numeric
    if (!qname_has_px_strs (STRa(qname), qfs)) return false;

    str_split_by_container (qname, qname_len, &qfs->con, qfs->con_prefix, qfs->con_prefix_len, item);
    if (!n_items) return false; // failed to split according to container - qname is not of this flavor

    // seg
    int64_t value, prev_value;

    bool is_int[n_items], is_hex[n_items], in_local[n_items], is_numeric[n_items];

    memset (is_int, 0, n_items * sizeof(bool));
    for (unsigned i=0; qfs->integer_items[i] != -1; i++)
        is_int[qfs->integer_items[i]] = true;

    memset (is_numeric, 0, n_items * sizeof(bool));
    for (unsigned i=0; qfs->numeric_items[i] != -1; i++)
        is_numeric[qfs->numeric_items[i]] = true;

    memset (is_hex, 0, n_items * sizeof(bool));
    for (unsigned i=0; qfs->hex_items[i] != -1; i++)
        is_hex[qfs->hex_items[i]] = true;
    
    memset (in_local, 0, n_items * sizeof(bool));
    for (unsigned i=0; qfs->in_local[i] != -1; i++)
        in_local[qfs->in_local[i]] = true;

    // seg container
    if (qfs == segconf.qname_flavor)
        seg_by_ctx (vb, STRa(qfs->con_snip), qname_ctx, qfs->num_seps + add_additional_bytes); // account for container separators, prefixes and caller-requested add_additional_bytes 
    else
        seg_by_ctx (vb, STRa(qfs->con_snip2), qname_ctx, qfs->num_seps + add_additional_bytes); 

    int encountered_mName=0, encountered_QNAME2=0;
    for (unsigned item_i=0; item_i < qfs->con.nitems_lo; item_i++) {

        // get item ctx - Q?NAME did_i are immediately following QNAME, and QmNAME si the last.
        const ContainerItem *item = &qfs->con.items[item_i];
        
        // calculate context - a bit messy, but faster than looking up by dict_id
        ContextP item_ctx = (item->dict_id.num == _SAM_QmNAME)   ? (qname_ctx + MAX_QNAME_ITEMS) 
                          : (item->dict_id.num == _FASTQ_QNAME2) ? CTX(FASTQ_QNAME2)
                          :                                        (qname_ctx + 1+item_i - encountered_mName - encountered_QNAME2);
                          
        encountered_mName  |= item->dict_id.num == _SAM_QmNAME;
        encountered_QNAME2 |= item->dict_id.num == _FASTQ_QNAME2;

        // case: this is the file is sorted by qname - delta against previous
        if ((item_i == qfs->ordered_item1 || item_i == qfs->ordered_item2) && 
            !segconf.sam_is_sorted && 
            str_get_int_dec (STRi(item, item_i), (uint64_t*)&value))

            seg_self_delta (vb, item_ctx, value, item_lens[item_i]);
        
        // case: end-of-range item, seg as delta vs previous item which is start-of-range
        else if (qfs->range_end_item == item_i &&
                 !flag.pair && // we don't use delta against other if we're pairing - too complicated
                 str_get_int (STRi(item, item_i), &value) && 
                 str_get_int (STRi(item, item_i-1), &prev_value)) {
            SNIP(32);
            seg_prepare_snip_other (SNIP_OTHER_DELTA, qfs->con.items[qfs->range_end_item-1].dict_id, true, value - prev_value, snip);
            seg_by_ctx (vb, STRa(snip), item_ctx, item_lens[item_i]);      
        }

        else if (in_local[item_i] && !flag.pair) { // note: we can't store in local if pairing        
            if (is_int[item_i] || is_numeric[item_i])
                seg_integer_or_not (vb, item_ctx, STRi(item, item_i), item_lens[item_i]);

            else 
                seg_add_to_local_text (vb, item_ctx, STRi(item, item_i), true, item_lens[item_i]);
        }

        // TO DO - add field xor_diff allowing specifying xor_diff method (if we find a case where its useful)
        // else if (xor_diff[item_i] && !flag.pair) {
        //    seg_xor_diff (vb, item_ctx, STRi(item, item_i), item_ctx->flags.all_the_same, item_lens[item_i]);
        //     seg_set_last_txt (vb, item_ctx, STRi(item, item_i));
        // }
        
        else if (item_i == qfs->qname2 && segconf.qname_flavor2) 
            qname_seg_qf (vb, CTX (FASTQ_QNAME2), segconf.qname_flavor2, STRi(item, item_i), 0);

        else if (item->separator[0] == CI0_SKIP)
            {} // no segging a skipped item

        // case: textual item
        else
            seg_by_ctx (vb, STRi(item, item_i), item_ctx, item_lens[item_i]);   
    }

    return true;
}

void qname_seg (VBlock *vb, Context *qname_ctx, STRp (qname), unsigned add_additional_bytes)  // account for characters in addition to the field
{
    START_TIMER;
    
    // copy if identical to previous (> 50% of lines in collated) - small improvement in compression and compression time
    // no need in sam_is_sorted, as already handled in sam_seg_QNAME with buddy 
    if (!segconf.sam_is_sorted && vb->line_i && !flag.optimize_DESC && is_same_last_txt (vb, qname_ctx, STRa(qname))) {
        seg_by_ctx (vb, STRa(copy_qname), qname_ctx, qname_len + add_additional_bytes);
        goto done;
    }

    bool success = false;
    if (segconf.qname_flavor)
        success = qname_seg_qf (VB, qname_ctx, segconf.qname_flavor, STRa(qname), add_additional_bytes);

    if (!success) 
        tokenizer_seg (VB, qname_ctx, STRa(qname), 
                       VB_DT(DT_FASTQ) ? sep_with_space : sep_without_space, 
                       add_additional_bytes);

done:
    seg_set_last_txt (vb, qname_ctx, STRa(qname)); // needed also for kraken in sam
    COPY_TIMER (qname_seg);
}

const char *qf_name (QnameFlavor qf)
{
    return qf ? qf->name : "";
}
