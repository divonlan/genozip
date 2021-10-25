// ------------------------------------------------------------------
//   qname.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "qname.h"
#include "tokenizer.h"
#include "buffer.h"
#include "vblock.h"
#include "segconf.h"
#include "strings.h"
#include "seg.h"
#include "sam.h"
#include "stats.h"

#define MAX_QNAME_ITEMS 9 // mate + 8 others (matching Q?NAME defined in sam.h, fastq.h etc) 
#define item_0_did_i (qname_did_i+1)
#define item_1_did_i (qname_did_i+2)
#define item_2_did_i (qname_did_i+3)
#define item_3_did_i (qname_did_i+4)
#define item_4_did_i (qname_did_i+5)
#define item_5_did_i (qname_did_i+6)
#define item_6_did_i (qname_did_i+7)
#define item_7_did_i (qname_did_i+8)
#define mate_did_i   (qname_did_i+9) // the mate if the qname ends with eg /1 or /2

static STRl(copy_qname, 50);

//-----------------------------------------------------------------------------------------
// Illumina-7 format: <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>  
// Example: HWI-D00360:5:H814YADXX:1:1106:10370:52569
// Example FASTQ: A00488:61:HMLGNDSXX:4:1101:1940:1000 2:N:0:CTGAAGCT+ATAGAGGC
// See here: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
// interpretation of the first field here: https://github.com/10XGenomics/supernova/blob/master/tenkit/lib/python/tenkit/illumina_instrument.py#L12-L45
//-----------------------------------------------------------------------------------------
static SmallContainer con_illumina_7 = {
    .repeats   = 1,
    .nitems_lo = 7,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":" },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q6NAME },                  } } 
};

static SmallContainer con_illumina_7_fq = {
    .repeats   = 1,
    .nitems_lo = 8,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":" },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q6NAME }, .separator = " " },
                   { .dict_id = { _SAM_Q7NAME },                  } } // entire part after space is segged together
};

//--------------------------------------------------------------------------------------------------------------
// BGI_E format: 27 characters: E, FlowCellSerialNumber[9], L, Lane[1], C, Column[3], R, Row[3] Tile[7]
// Example: E100020409L1C001R0030000234
// See: https://github.com/IMB-Computational-Genomics-Lab/BGIvsIllumina_scRNASeq
//--------------------------------------------------------------------------------------------------------------
static SmallContainer con_bgi_E = {
    .repeats             = 1,
    .nitems_lo           = 5,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI_FIXED_0_PAD, 9 } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI_FIXED_0_PAD, 1 } }, // note: CI_FIXED_0_PAD is only useful if the field is segged as a number (which may be shorter than the field length). if it is segged as a text, the separator will have no effect.
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI_FIXED_0_PAD, 3 } },
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI_FIXED_0_PAD, 3 } },
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI_FIXED_0_PAD, 7 } } } 
};

//--------------------------------------------------------------------------------------------------------------------
// BGI_CL format: CL, FlowCellSerialNumber[9], L, Lane[1], C, Column[3], R, Row[3], _, variable-length number
//  CL100025298L1C002R050_244547 reported by https://en.wikipedia.org/wiki/File:BGI_seq_platform_read_name_description.png
// @CL100072652L2C001R001_12/1 for FASTQ reported by https://www.biostars.org/p/395715/
//--------------------------------------------------------------------------------------------------------------------
static SmallContainer con_bgi_CL = {
    .repeats             = 1,
    .nitems_lo           = 5,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI_FIXED_0_PAD, 9 } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI_FIXED_0_PAD, 1 } },
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI_FIXED_0_PAD, 3 } },
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI_FIXED_0_PAD, 3 } },
                             { .dict_id = { _SAM_Q4NAME },                                  } } 
};

//--------------------------------------------------------------------------------------------------------------
// ION_TORRENT_3 format: 17 characters: RunID[5] : Row[5] : Column[5]
// Example: ZEWTM:10130:07001
// See: https://www.biostars.org/p/330644/
//--------------------------------------------------------------------------------------------------------------
static SmallContainer con_ion_torrent_3 = {
    .repeats   = 1,
    .nitems_lo = 3,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":" },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q2NAME },                  } } 
};

//--------------------------------------------------------------------------------------------------------------
// Illumina-5 (old, including Solexa) format: <machine_id>:<lane>:<tile>:<x_coord>:<y_coord> 
// Example: SOLEXA-1GA-1_4_FC20ENL:7:258:737:870 (all fields are variable-length)
// Example: HWI-ST550_0201:3:1101:1626:2216#ACAGTG
// See: https://www.biostars.org/p/330644/
//--------------------------------------------------------------------------------------------------------------
static SmallContainer con_illumina_5 = {
    .repeats   = 1,
    .nitems_lo = 5,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":" },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q4NAME },                  } } 
};
static SmallContainer con_illumina_5i = {
    .repeats   = 1,
    .nitems_lo = 6,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":" },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q4NAME }, .separator = "#" },
                   { .dict_id = { _SAM_Q5NAME },                  } } 
};

//--------------------------------------------------------------------------------------------------------------------
// ROCHE_454 format: 
// Example: 000050_1712_0767
// See: https://www.sequencing.uio.no/illumina-services/4.%20Data%20delivery/results_454.html
//--------------------------------------------------------------------------------------------------------------------
static SmallContainer con_roche_454 = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI_FIXED_0_PAD, 6 } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI_FIXED_0_PAD, 4 } },
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI_FIXED_0_PAD, 4 } } } 
};

// See: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats
static SmallContainer con_helicos = {
    .repeats             = 1,
    .nitems_lo           = 6,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "-" }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-" },
                             { .dict_id = { _SAM_Q2NAME }, .separator = "-" },
                             { .dict_id = { _SAM_Q3NAME }, .separator = "-" },
                             { .dict_id = { _SAM_Q4NAME }, .separator = "-" },
                             { .dict_id = { _SAM_Q5NAME }                   } } 
};

static SmallContainer con_pacbio_3 = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "_" }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "_" },
                             { .dict_id = { _SAM_Q2NAME }                   } } 
};

// <MovieName>/<ZMW_number>/<subread-start>_<subread-end>
// example: m130802_221257_00127_c100560082550000001823094812221334_s1_p0/128361/872_4288
static SmallContainer con_pacbio_range = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "/" }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "/" },
                             { .dict_id = { _SAM_Q2NAME }, .separator = "_" },
                             { .dict_id = { _SAM_Q3NAME }                   } } 
};

// <MovieName>/<ZMW_number>/ccs
// example: m64136_200621_234916/18/ccs
static SmallContainer con_pacbio_label = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "/" }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "/" },
                             { .dict_id = { _SAM_Q2NAME }                   } } 
};

// @<MovieName>/<ZMW_number>
static SmallContainer con_pacbio_plain = {
    .repeats             = 1,
    .nitems_lo           = 2,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "/" }, 
                             { .dict_id = { _SAM_Q1NAME }                   } } 
};

// example: af84b0c1-6945-4323-9193-d9f6f2c38f9a
static SmallContainer con_nanopore = {
    .repeats             = 1,
    .nitems_lo           = 5,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "-" }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-" }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = "-" }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = "-" }, 
                             { .dict_id = { _SAM_Q4NAME }                   } } 
};

static SmallContainer con_ncbi_sra = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI_FIXED_0_PAD, 3 } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "."    }, 
                             { .dict_id = { _SAM_Q2NAME }                      } } 
};

static SmallContainer con_genozip_opt = {
    .repeats             = 1,
    .nitems_lo           = 2,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "."    }, 
                             { .dict_id = { _SAM_Q1NAME }                      } } 
};
//----------------------------------------------------------------

typedef struct { STR(s); int32_t fixed_i; } FixedStr;

typedef struct QnameFlavorStruct {
    char name[16], example[112];
    SeqTech tech;                              // The sequencing technology used to generate this data
    bool fq_only;                              // this QF is only for FASTQ, not SAM/BAM/KRAKEN
    SmallContainer *con_template;              // container template - copied to con in qname_zip_initialize
    unsigned num_seps;                         // number of printable separator and prefix characters
    int integer_items[MAX_QNAME_ITEMS+1];      // array of item_i (0-based) that are expected to be integer - no leading zeros (+1 for termiantor; terminated by -1)
    int numeric_items[MAX_QNAME_ITEMS+1];      // array of item_i (0-based) that are expected to be numeric - leading zeros ok (+1 for termiantor; terminated by -1)
    int hex_items[MAX_QNAME_ITEMS+1];          // array of item_i (0-based) that are expected to be hexademical (+1 for termiantor; terminated by -1)
    int ordered_item1, ordered_item2;          // an item that may be delta'd in non-sorted files. this should be the fast-moving item which consecutive values are close
    int range_end_item;                        // item containing end of range - delta vs preceding item in container
    unsigned fixed_len;                        // fixed length (0 if it is not fixed length)
    FixedStr px_strs[MAX_QNAME_ITEMS+1];       // fixed substrings in template (=prefixes of items in fixed locations). (+1 for termiantor; terminated by empty entry) 
    STRl (con_snip, 200);                      // container snips (generated in qname_zip_initialize)
    STRl (con_prefix, 30);                     // prefix of container
    SmallContainer con;                        // container
} QnameFlavorStruct;

static QnameFlavorStruct qf[] = { 
/*  mate   name             example                                   tech     fq_only   con_template     #sp int_items      numeric_items    hex_items       ord1,2 rng len px_strs{str,str_len,fixed_i}                               */
         { "Unknown"                                                                                                                                                                                               }, 
         { "Illumina-fastq","A00488:61:HMLGNDSXX:4:1101:4345:1000"
           " 2:N:0:CTGAAGCT+ATAGAGGC",                                TECH_ILLUM_7, 1, &con_illumina_7_fq, 7, {1,3,4,5,6,-1}, {-1},           {-1},           5,6,   -1,                                                                 },
    {},  { "Illumina",      "A00488:61:HMLGNDSXX:4:1101:4345:1000",   TECH_ILLUM_7, 0, &con_illumina_7,    6, {1,3,4,5,6,-1}, {-1},           {-1},           5,6,   -1,                                                                 },
    {},  { "BGI-E",         "E100020409L1C001R0030000801",            TECH_BGI,     0, &con_bgi_E,         4, {-1},           {0,1,2,3,4,-1}, {-1},           4,-1,  -1, 27, { {"E",1,0},  {"L",1,10},{"C",1,12},{"R",1,16},{"",0,20} }  },
    {},  { "BGI-CL",        "CL100025298L1C002R050_244547",           TECH_BGI,     0, &con_bgi_CL,        5, {-1},           {0,1,2,3,4,-1}, {-1},           4,-1,  -1, 0,  { {"CL",2,0}, {"L",1,11},{"C",1,13},{"R",1,17},{"_",1,21} } },
    {},  { "IonTorrent",    "ZEWTM:10130:07001",                      TECH_IONTORR, 0, &con_ion_torrent_3, 2, {-1},           {1,2,-1},       {-1},           -1,-1, -1, 17                                                              },
    {},  { "Illumina-old#", "HWI-ST550_0201:3:1101:1626:2216#ACAGTG", TECH_ILLUM_5, 0, &con_illumina_5i,   5, {1,2,3,4,-1},   {-1},           {-1},           -1,-1, -1,                                                                 },
    {},  { "Illumina-old",  "SOLEXA-1GA-1_4_FC20ENL:7:258:737:870",   TECH_ILLUM_5, 0, &con_illumina_5,    4, {1,2,3,4,-1},   {-1},           {-1},           -1,-1, -1,                                                                 },
    {},  { "Roche-454",     "000050_1712_0767",                       TECH_454,     0, &con_roche_454,     2, {-1},           {0,1,2,-1},     {-1},           -1,-1, -1, 16, { {"",0,0}, {"_",1,6},{"_",1,11} }                          },
    {},  { "Helicos",       "VHE-242383071011-15-1-0-2",              TECH_HELICOS, 0, &con_helicos,       5, {2,3,4,5,-1},   {1,-1},         {-1},           -1,-1, -1                                                                  },
    {},  { "PacBio-3",      "56cdb76f_70722_4787",                    TECH_PACBIO,  0, &con_pacbio_3,      2, {1,2,-1},       {-1},           {0,-1},         -1,-1, -1                                                                  },
    {},  { "PacBio-Range",  "m130802_221257_00127_"
           "c100560082550000001823094812221334_s1_p0/128361/872_4288",TECH_PACBIO,  0, &con_pacbio_range,  4, {1,2,3,-1},     {-1},           {-1},           -1,-1, 3,  0,  { { "m",1,0} }                                              },
    {},  { "PacBio-Label",  "m64136_200621_234916/18/ccs",            TECH_PACBIO,  0, &con_pacbio_label,  3, {1,-1},         {-1},           {-1},           -1,-1, -1, 0,  { { "m",1,0} }                                              },
    {},  { "PacBio-Plain",  "m64136_200621_234916/18",                TECH_PACBIO,  0, &con_pacbio_plain,  2, {1,-1},         {-1},           {-1},           -1,-1, -1, 0,  { { "m",1,0} }                                              },
    {},  { "Nanopore",      "af84b0c1-6945-4323-9193-d9f6f2c38f9a",   TECH_ONP,     0, &con_nanopore,      4, {-1},           {-1},           {0,1,2,3,4,-1}, -1,-1, -1, 36,                                                             },
    {},  { "NCBI-SRA",      "SRR001666.1",                            TECH_UNKNOWN, 0, &con_ncbi_sra,      1, {2,-1},         {1,-1},         {-1},           2,-1,  -1,                                                                 },
    {},  { "Genozip-opt",   "basic.1",                                TECH_UNKNOWN, 1, &con_genozip_opt,   1, {1,-1},         {-1},           {-1},           1,-1,  -1,                                                                 },
};

#define NUM_QFs (sizeof(qf)/sizeof(qf[0]))

static inline unsigned qname_get_px_str_len (QnameFlavorStruct *qfs)
{
    unsigned i=0; 
    while (qfs->px_strs[i].s) i++;
    return i;
}

// note: we need to re-initialize in each file, lest the data type has changed
void qname_zip_initialize (DictId qname_dict_id)
{
    // we need to prepare the containers only once and they can serve all files, as the containers
    // are data-type independent (since all data types use the dict_id Q?NAME for qnames)
    DO_ONCE {
        for (unsigned qf_i=1; qf_i < NUM_QFs; qf_i++) {
            QnameFlavorStruct *qfs = &qf[qf_i];
            
            // case: version of the next entry with room for /1 /2 - populate now
            bool is_paired = !qfs->name[0]; // a mate {} means we should prepare a version that can end with /1 /2
            if (is_paired) *qfs = *(qfs+1);

            // create container
            qfs->con = *qfs->con_template;

            if (is_paired) {
                // add new item as numeric item (at beginning of array - easier, and the order doesn't matter)
                memmove (&qfs->numeric_items[1], qfs->numeric_items, sizeof(qfs->numeric_items) - sizeof(qfs->numeric_items[0])); // make room
                qfs->numeric_items[0] = qfs->con.nitems_lo;
                qfs->num_seps++; // one more separator - the '/'

                if (qfs->fixed_len) qfs->fixed_len += 2;

                strcpy (&qfs->name[strlen(qfs->name)], "/");
                strcpy (&qfs->example[strlen(qfs->example)], "/1");
                
                // if last item has no separator: add a '/' separator for the last item
                if (!qfs->con.items[qfs->con.nitems_lo-1].separator[0])
                    qfs->con.items[qfs->con.nitems_lo-1].separator[0] = '/';
                
                // otherwise: add a '/' as a prefix for the next item
                else {
                    // items with a set last separator, MUST have px_strs set for ALL items
                    ASSERT (qname_get_px_str_len (qfs) == qfs->con.nitems_lo, "QnameFlavor=%s: Since separator[0] is set for the last item, px_str MUST contain exactly the same number of elements as the container",
                            qfs->name);

                    qfs->px_strs[qfs->con.nitems_lo] = (FixedStr){ .s="/", .s_len=1, .fixed_i=-2 /* qname_len-2 */ };
                }
                
                // add the new item
                qfs->con.items[qfs->con.nitems_lo++] = (ContainerItem){ .dict_id = { _SAM_QmatNAME } };
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

            // prepare container snip for this QF
            qfs->con_snip_len = sizeof (qfs->con_snip);
            container_prepare_snip ((Container*)&qfs->con, STRa(qfs->con_prefix), qfs->con_snip, &qfs->con_snip_len);
        }
    }

    // prepare copy_qname snip - every time, as container_did_i changes between data types
    seg_prepare_snip_other (SNIP_COPY, qname_dict_id, false, 0, copy_qname); // QNAME dict_id is the same for SAM, FASTQ, KRAKEN

    segconf.tech = TECH_UNKNOWN; // initialize
}

void qname_seg_initialize (VBlockP vb, DidIType qname_did_i)
{
    if (!segconf.qname_flavor) return;
    QnameFlavorStruct *qfs = &qf[segconf.qname_flavor];

    // consolidate all items under QNAME in --stats
    stats_set_consolidation (vb, qname_did_i, MAX_QNAME_ITEMS, 
                             item_0_did_i, item_1_did_i, item_2_did_i, item_3_did_i, 
                             item_4_did_i, item_5_did_i, item_6_did_i, item_7_did_i, mate_did_i); 

    // set STORE_INT as appropriate
    if (qfs->ordered_item1  != -1) CTX(item_0_did_i + qfs->ordered_item1)->flags.store = STORE_INT; 
    if (qfs->ordered_item2  != -1) CTX(item_0_did_i + qfs->ordered_item2)->flags.store = STORE_INT; 
    if (qfs->range_end_item != -1) CTX(item_0_did_i + qfs->range_end_item - 1)->flags.store = STORE_INT; 
}

static inline bool qname_has_px_strs (STRp(qname), const QnameFlavorStruct *qfs)
{
    if (!qfs->px_strs) return true; // flavor has no px_strs requirements - all good

    for (const FixedStr *px = qfs->px_strs; px->s ; px++) {
        int i = (px->fixed_i >= 0) ? px->fixed_i : ((int)qname_len + px->fixed_i);
        if (i + px->s_len > qname_len ||          // case: fixed str goes beyond the end of qname
            memcmp (px->s, &qname[i], px->s_len)) // case: fixed str does not appear in qname
            return false;
    }

    return true; // all fixed strings appear as expected
}

// note: we run this function only in discovery, not in segging, because it is quite expesive - checking all numerics.
static bool qname_is_flavor (STRp(qname), const QnameFlavorStruct *qfs)
{
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

    return true; // yes, qname is of this flavor
}

// called for the first line in segconf.running
void qname_segconf_discover_flavor (VBlockP vb, DidIType qname_did_i, STRp(qname))
{
    segconf.qname_flavor = 0; // unknown

    for (unsigned qf_i=1; qf_i < NUM_QFs; qf_i++) 
        if ((VB_DT(DT_FASTQ) || !qf[qf_i].fq_only) && qname_is_flavor (qname, qname_len, &qf[qf_i])) {
            segconf.qname_flavor = qf_i;
            segconf.tech = qf[qf_i].tech;
            qname_seg_initialize (vb, qname_did_i); // so the rest of segconf.running can seg fast using the discovered container
            break;
        }

    if (flag.debug) // unit test
        for (unsigned qf_i=1; qf_i < NUM_QFs; qf_i++)
            ASSERT (qname_is_flavor (qf[qf_i].example, strlen (qf[qf_i].example), &qf[qf_i]), 
                    "Failed to identify qname \"%s\" as %s (qf_i=%u)", qf[qf_i].example, qf[qf_i].name, qf_i);
}

// attempt to seg according to the qf - return true if successful
static inline bool qname_seg_qf (VBlockP vb, ContextP qname_ctx, const QnameFlavorStruct *qfs, STRp(qname), unsigned add_additional_bytes)
{
    // check that the fixed strings expected by the flavor appear in qname. note: we don't verify numeric fields as
    // it would be expensive. instead, we fallback on segging them as text if they turn out to be not numeric
    if (!qname_has_px_strs (STRa(qname), qfs)) return false;

    str_split_by_container (qname, qname_len, &qfs->con, qfs->con_prefix, qfs->con_prefix_len, item);
    if (!n_items) return false; // failed to split according to container - qname is not of this flavor

    // seg
    int64_t value, prev_value;

    for (unsigned item_i=0; item_i < qfs->con.nitems_lo; item_i++) {

        // get item ctx - Q?NAME did_i are immediately following QNAME, and QmatNAME si the last.
        ContextP item_ctx = qname_ctx + ((qfs->con.items[item_i].dict_id.num == _SAM_QmatNAME) ? MAX_QNAME_ITEMS : (1+item_i)); 
        
        // case: this is the file is sorted by qname - delta against previous
        if ((item_i == qfs->ordered_item1 || item_i == qfs->ordered_item2) && !segconf.sam_is_sorted && str_get_int_dec (STRi(item, item_i), (uint64_t*)&value)) 
            seg_self_delta (vb, item_ctx, value, item_lens[item_i]);
        
        // case: end-of-range item, seg as delta vs previous item which is start-of-range
        else if (qfs->range_end_item == item_i &&
                 str_get_int (STRi(item, item_i), &value) && 
                 str_get_int (STRi(item, item_i-1), &prev_value)) {
            SNIP(50);
            seg_prepare_snip_other (SNIP_OTHER_DELTA, qfs->con.items[qfs->range_end_item-1].dict_id, true, value - prev_value, snip);
            seg_by_ctx (vb, STRa(snip), item_ctx, item_lens[item_i]);      
        }
        // case: textual item
        else
            seg_by_ctx (vb, STRi(item, item_i), item_ctx, item_lens[item_i]);      
    }

    // seg container
    seg_by_ctx (vb, STRa(qfs->con_snip), qname_ctx, qfs->num_seps + add_additional_bytes); // account for container separators, prefixes and caller-requested add_additional_bytes 

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
        success = qname_seg_qf (VB, qname_ctx, &qf[segconf.qname_flavor], STRa(qname), add_additional_bytes);

    if (!success) 
        tokenizer_seg (VB, qname_ctx, STRa(qname), 
                       VB_DT(DT_FASTQ) ? sep_with_space : sep_without_space, 
                       add_additional_bytes);

done:
    seg_set_last_txt (vb, qname_ctx, STRa(qname), STORE_NONE); // needed also for kraken in sam
    COPY_TIMER (qname_seg);
}

const char *qf_name (unsigned qf_i)
{
    return qf[qf_i].name;
}
