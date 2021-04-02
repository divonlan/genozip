// ------------------------------------------------------------------
//   dict_id.c
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "dict_id.h"
#include "data_types.h"
#include "file.h"
#include "zfile.h"
#include "sections.h"
#include "vblock.h"
#include "libdeflate/libdeflate.h"
#include "vcf.h"

// globals externed in dict_id.h and initialized in dict_id_initialize
static Buffer dict_id_aliases_buf  = EMPTY_BUFFER;
const DictIdAlias *dict_id_aliases = NULL;
uint32_t dict_id_num_aliases = 0;

uint64_t dict_id_fields[MAX_NUM_FIELDS_PER_DATA_TYPE];

// VCF stuff
uint64_t dict_id_FORMAT_PL=0, dict_id_FORMAT_GL=0, dict_id_FORMAT_GP=0, dict_id_FORMAT_DP=0, dict_id_FORMAT_MIN_DP=0, 
         dict_id_FORMAT_PS=0, dict_id_FORMAT_GT=0, dict_id_FORMAT_GT_HT=0, dict_id_FORMAT_GT_HT_INDEX=0,
         dict_id_PBWT_RUNS=0, dict_id_PBWT_FGRC=0,
         dict_id_FORMAT_GT_SHARK_DB=0, dict_id_FORMAT_GT_SHARK_GT=0, dict_id_FORMAT_GT_SHARK_EX=0,
         dict_id_FORMAT_AD=0, dict_id_FORMAT_ADF=0, dict_id_FORMAT_ADR=0, dict_id_FORMAT_ADALL=0, 
         dict_id_FORMAT_GQ=0, dict_id_FORMAT_DS=0,
         dict_id_INFO_AC=0, dict_id_INFO_AF=0, dict_id_INFO_AN=0, dict_id_INFO_DP=0, dict_id_INFO_VQSLOD=0,
         dict_id_INFO_END=0, dict_id_INFO_SVLEN=0, dict_id_INFO_DP4=0, dict_id_INFO_SF=0,
         dict_id_INFO_LIFTOVER=0, dict_id_INFO_LIFTBACK=0, dict_id_INFO_LIFTREJT=0, 
         dict_id_INFO_BaseCounts=0,

         // tags from VEP (Varient Effect Predictor) and similar tools
         dict_id_INFO_CSQ=0, dict_id_INFO_vep=0, dict_id_INFO_DP_HIST=0, dict_id_INFO_GQ_HIST=0, 
         dict_id_INFO_AGE_HISTOGRAM_HET=0, dict_id_INFO_AGE_HISTOGRAM_HOM=0; 
   
// SAM stuff
uint64_t dict_id_OPTION_AM=0, dict_id_OPTION_AS=0, dict_id_OPTION_CM=0, dict_id_OPTION_LB=0, dict_id_OPTION_FI=0, dict_id_OPTION_H0=0,
         dict_id_OPTION_H1=0, dict_id_OPTION_H2=0, dict_id_OPTION_MD=0, dict_id_OPTION_MQ=0, dict_id_OPTION_NH=0, dict_id_OPTION_NM=0, 
         dict_id_OPTION_OA=0, dict_id_OPTION_OC=0, dict_id_OPTION_PG=0, dict_id_OPTION_E2=0, dict_id_OPTION_U2=0,
         dict_id_OPTION_PQ=0, dict_id_OPTION_PU=0, dict_id_OPTION_RG=0, dict_id_OPTION_SA=0, dict_id_OPTION_SM=0, dict_id_OPTION_TC=0, 
         dict_id_OPTION_UQ=0, dict_id_OPTION_CC=0, dict_id_OPTION_MC=0,
         dict_id_OPTION_X0=0, dict_id_OPTION_X1=0, dict_id_OPTION_XA=0, dict_id_OPTION_XA_RNAME=0, dict_id_OPTION_XN=0, dict_id_OPTION_XM=0, dict_id_OPTION_XO=0,
         dict_id_OPTION_XG=0, dict_id_OPTION_XS=0, dict_id_OPTION_XE=0,
         dict_id_OPTION_mc=0, dict_id_OPTION_ms=0,
         dict_id_OPTION_BD=0, dict_id_OPTION_BI=0, dict_id_OPTION_BD_BI=0,
         dict_id_OPTION_ZM=0,
  
         // private genozip dict
         dict_id_OPTION_STRAND=0, dict_id_OPTION_RNAME=0, dict_id_OPTION_POS=0, dict_id_OPTION_CIGAR=0, dict_id_OPTION_MAPQ=0; 

// GVF stuff
uint64_t dict_id_ATTR_ID=0, dict_id_ATTR_Variant_seq=0, dict_id_ATTR_Reference_seq=0, dict_id_ATTR_Variant_freq=0,
         dict_id_ATTR_Dbxref=0, // from from GRCh37/38 - example: "dbSNP_151:rs1282280967"
         dict_id_ATTR_ancestral_allele=0, // from from GRCh37/38 - example ancestral_allele=GTTA
         dict_id_ATTR_Variant_effect=0, // example: "Variant_effect=non_coding_transcript_variant 0 ncRNA ENST00000431238,intron_variant 0 primary_transcript ENST00000431238"
         dict_id_ATTR_sift_prediction=0, dict_id_ATTR_polyphen_prediction=0, dict_id_ATTR_variant_peptide=0,

         dict_id_ENSTid=0; // private genozip dict

// our stuff used in multiple data types
uint64_t dict_id_WindowsEOL=0;         

DictId dict_id_make (const char *str, unsigned str_len, DictIdType dict_id_type) 
{ /*
    DictId dict_id = DICT_ID_NONE; 

    if (!str_len) str_len = strlen (str);

    if (str_len <= DICT_ID_LEN) 
        memcpy (dict_id.id, str, str_len);
    
    // case: name is too long - compress it
    else { 
        union { uint32_t i; uint8_t s[4]; } adler = {
            .i = libdeflate_adler32 (1, str, str_len)
        };
        
        // compressed dict_id uses 2 characters from the Adler32 to disambiguate, including one in id[1] that
        // goes into the map, hopefully avoiding a map contention
        dict_id.id[0] = str[0];
        dict_id.id[1] = adler.s[0];
        dict_id.id[2] = str[2];
        dict_id.id[3] = str[3];
        dict_id.id[4] = adler.s[1];
        dict_id.id[5] = str[str_len-3];
        dict_id.id[6] = str[str_len-2];
        dict_id.id[7] = str[str_len-1];
    }

    return dict_id;
*/
    DictId dict_id = DICT_ID_NONE; 

    if (!str_len) str_len = strlen (str);

    if (str_len <= DICT_ID_LEN) 
        memcpy (dict_id.id, str, str_len);
    
    else { 
        #define half1_len (DICT_ID_LEN/2)
        #define half2_len (DICT_ID_LEN - DICT_ID_LEN/2)

        memcpy (dict_id.id, str, half1_len); // take 1/2 from the start and 1/2 from then end
        memcpy (dict_id.id + half1_len, str+str_len-half2_len, half2_len);
    }

    switch (dict_id_type) {
        case DTYPE_FIELD  : dict_id.id[0] = dict_id.id[0] & 0x3f; break;
        case DTYPE_1      : dict_id.id[0] = dict_id.id[0] | 0xc0; break;
        case DTYPE_2      : break;
        default: ABORT ("Error in dict_id_make: invalid type %d", dict_id_type);
    }
    return dict_id;
}

void dict_id_initialize (DataType data_type) 
{   
    ASSERT0 (data_type != DT_NONE, "data_type is DT_NONE");

    for (int f=0; f < dt_fields[data_type].num_fields; f++) {
        const char *field_name = dt_fields[data_type].names[f];

        ASSERT (field_name, "Data type %s is missing a field name in DATA_TYPE_FIELDS for field %u", dt_name (data_type), f);
        dict_id_fields[f] = dict_id_make (field_name, strlen (field_name), DTYPE_FIELD).num; 
    }

    dict_id_WindowsEOL = dict_id_make ("#", 1, DTYPE_1).num; 

    switch (data_type) { 
    case DT_VCF:
    case DT_BCF:
        dict_id_FORMAT_DP     = dict_id_make ("DP", 2, DTYPE_VCF_FORMAT).num;
        dict_id_FORMAT_DS     = dict_id_make ("DS", 2, DTYPE_VCF_FORMAT).num;
        dict_id_FORMAT_GT     = dict_id_make ("GT", 2, DTYPE_VCF_FORMAT).num;
        dict_id_FORMAT_GT_HT  = dict_id_make ("@HT", 3, DTYPE_VCF_FORMAT).num; // different first 2 letters than GT, for lookup table
        dict_id_FORMAT_GT_HT_INDEX = dict_id_make ("@INDEXHT", 8, DTYPE_VCF_FORMAT).num; // different first 2 letters
        dict_id_FORMAT_GT_SHARK_DB = dict_id_make ("@1SHRKDB", 8, DTYPE_VCF_FORMAT).num; // different first 2 letters, 
        dict_id_FORMAT_GT_SHARK_GT = dict_id_make ("@2SHRKGT", 8, DTYPE_VCF_FORMAT).num;
        dict_id_FORMAT_GT_SHARK_EX = dict_id_make ("@3SHRKEX", 8, DTYPE_VCF_FORMAT).num;  
        dict_id_PBWT_RUNS     = dict_id_make ("@1BWTRUN", 8, DTYPE_VCF_FORMAT).num; // PBWT runs
        dict_id_PBWT_FGRC     = dict_id_make ("@2BWTFGR", 8, DTYPE_VCF_FORMAT).num; // PBWT foreground run count
        dict_id_FORMAT_PL     = dict_id_make ("PL", 2, DTYPE_VCF_FORMAT).num;
        dict_id_FORMAT_GP     = dict_id_make ("GP", 2, DTYPE_VCF_FORMAT).num;
        dict_id_FORMAT_GL     = dict_id_make ("GL", 2, DTYPE_VCF_FORMAT).num;
        dict_id_FORMAT_GQ     = dict_id_make ("GQ", 2, DTYPE_VCF_FORMAT).num;
        dict_id_FORMAT_AD     = dict_id_make ("AD", 2, DTYPE_VCF_FORMAT).num;
        dict_id_FORMAT_ADF    = dict_id_make ("ADF", 3, DTYPE_VCF_FORMAT).num;
        dict_id_FORMAT_ADR    = dict_id_make ("ADR", 3, DTYPE_VCF_FORMAT).num;
        dict_id_FORMAT_ADALL  = dict_id_make ("^ADALL", 6, DTYPE_VCF_FORMAT).num; // different 2 letters than AD
        dict_id_INFO_AC       = dict_id_make ("AC", 2, DTYPE_VCF_INFO).num;
        dict_id_INFO_AF       = dict_id_make ("AF", 2, DTYPE_VCF_INFO).num;
        dict_id_INFO_AGE_HISTOGRAM_HET = dict_id_make ("AGE_HISTOGRAM_HET", 17, DTYPE_VCF_INFO).num; 
        dict_id_INFO_AGE_HISTOGRAM_HOM = dict_id_make ("AGE_HISTOGRAM_HOM", 17, DTYPE_VCF_INFO).num;
        dict_id_INFO_AN       = dict_id_make ("AN", 2, DTYPE_VCF_INFO).num;
        dict_id_INFO_BaseCounts = dict_id_make ("BaseCounts", 10, DTYPE_VCF_INFO).num;
        dict_id_INFO_CSQ      = dict_id_make ("CSQ", 3, DTYPE_VCF_INFO).num;
        dict_id_INFO_DP       = dict_id_make ("DP", 2, DTYPE_VCF_INFO).num;
        dict_id_INFO_DP4      = dict_id_make ("DP4", 3, DTYPE_VCF_INFO).num;
        dict_id_INFO_DP_HIST  = dict_id_make ("DP_HIST", 7, DTYPE_VCF_INFO).num; // unfortunately there's a 2-letter conflict with DP, but we can't change names of INFO fields 
        dict_id_INFO_GQ_HIST  = dict_id_make ("GQ_HIST", 7, DTYPE_VCF_INFO).num;
        dict_id_INFO_SF       = dict_id_make ("SF", 2, DTYPE_VCF_INFO).num;
        dict_id_INFO_vep      = dict_id_make ("vep", 3, DTYPE_VCF_INFO).num;
        dict_id_INFO_VQSLOD   = dict_id_make ("VQSLOD", 6, DTYPE_VCF_INFO).num;
        dict_id_INFO_LIFTOVER = dict_id_make (INFO_LIFTOVER, INFO_LIFTOVER_LEN, DTYPE_VCF_INFO).num;
        dict_id_INFO_LIFTBACK = dict_id_make (INFO_LIFTBACK, INFO_LIFTBACK_LEN, DTYPE_VCF_INFO).num;
        dict_id_INFO_LIFTREJT = dict_id_make (INFO_LIFTREJT, INFO_LIFTREJT_LEN, DTYPE_VCF_INFO).num;

        // Added by GATK HaplotypeCaller in a gVCF: https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format
        dict_id_INFO_END      = dict_id_make ("END", 3, DTYPE_VCF_INFO).num;
        dict_id_INFO_SVLEN    = dict_id_make ("SVLEN", 5, DTYPE_VCF_INFO).num;
        dict_id_FORMAT_MIN_DP = dict_id_make ("MIN_DP", 6, DTYPE_VCF_FORMAT).num;

        // This appears if the VCF line has a Windows-style \r\n line ending
        break;

    case DT_SAM:
    case DT_BAM:
        dict_id_OPTION_AM = dict_id_make ("AM:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_AS = dict_id_make ("AS:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_CC = dict_id_make ("CC:Z", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_BD = dict_id_make ("BD:Z", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_BI = dict_id_make ("BI:Z", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_BD_BI = dict_id_make ("BD_BI", 5, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_CM = dict_id_make ("CM:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_E2 = dict_id_make ("E2:Z", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_FI = dict_id_make ("FI:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_H0 = dict_id_make ("H0:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_H1 = dict_id_make ("H1:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_H2 = dict_id_make ("H2:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_LB = dict_id_make ("LB:Z", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_MC = dict_id_make ("MC:Z", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_MD = dict_id_make ("MD:Z", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_MQ = dict_id_make ("MQ:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_NH = dict_id_make ("NH:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_NM = dict_id_make ("NM:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_OA = dict_id_make ("OA:Z", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_OC = dict_id_make ("OC:Z", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_PG = dict_id_make ("PG:Z", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_PQ = dict_id_make ("PQ:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_PU = dict_id_make ("PU:Z", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_RG = dict_id_make ("RG:Z", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_SA = dict_id_make ("SA:Z", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_SM = dict_id_make ("SM:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_TC = dict_id_make ("TC:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_UQ = dict_id_make ("UQ:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_U2 = dict_id_make ("U2:Z", 4, DTYPE_SAM_OPTIONAL).num;
                
        // Ion Torrent flow signal array
        dict_id_OPTION_ZM = dict_id_make ("ZM:B", 4, DTYPE_SAM_OPTIONAL).num;

        // bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
        dict_id_OPTION_X0 = dict_id_make ("X0:i", 4, DTYPE_SAM_OPTIONAL).num; 
        dict_id_OPTION_X1 = dict_id_make ("X1:i", 4, DTYPE_SAM_OPTIONAL).num; 
        dict_id_OPTION_XA = dict_id_make ("XA:Z", 4, DTYPE_SAM_OPTIONAL).num; 
        dict_id_OPTION_XA_RNAME = dict_id_make ("X0ARNAME", 8, DTYPE_SAM_OPTIONAL).num; 
        dict_id_OPTION_XN = dict_id_make ("XN:i", 4, DTYPE_SAM_OPTIONAL).num; 
        dict_id_OPTION_XM = dict_id_make ("XM:i", 4, DTYPE_SAM_OPTIONAL).num; 
        dict_id_OPTION_XO = dict_id_make ("XO:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_XG = dict_id_make ("XG:i", 4, DTYPE_SAM_OPTIONAL).num; 
        dict_id_OPTION_XS = dict_id_make ("XS:i", 4, DTYPE_SAM_OPTIONAL).num; 
        dict_id_OPTION_XE = dict_id_make ("XE:i", 4, DTYPE_SAM_OPTIONAL).num;

        // biobambam tags
        dict_id_OPTION_mc = dict_id_make ("mc:i", 4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_ms = dict_id_make ("ms:i", 4, DTYPE_SAM_OPTIONAL).num;

        // added by GATK's BQSR (Base Quality Score Recalibration)
        dict_id_OPTION_BD = dict_id_make ("BD:Z", 4, DTYPE_SAM_OPTIONAL).num; // not used in newer versions of GATK
        dict_id_OPTION_BI = dict_id_make ("BI:Z", 4, DTYPE_SAM_OPTIONAL).num; // not used in newer versions of GATK

        // our private dictionary for + or 0 strands
        dict_id_OPTION_STRAND = dict_id_make ("@STRAND", 7, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_RNAME  = dict_id_make ("@RNAME",  6, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_POS    = dict_id_make ("@POS",    4, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_CIGAR  = dict_id_make ("@CIGAR",  6, DTYPE_SAM_OPTIONAL).num;
        dict_id_OPTION_MAPQ   = dict_id_make ("@MAPQ",   5, DTYPE_SAM_OPTIONAL).num;

        break;

    case DT_FASTA:
        break;

    case DT_GFF3:
        // standard GVF fields (ID is also a standard GFF3 field)
        dict_id_ATTR_ID               = dict_id_make ("ID", 2, DTYPE_GFF3_ATTR).num;
        dict_id_ATTR_Variant_seq      = dict_id_make ("Variant_seq", 0, DTYPE_GFF3_ATTR).num;
        dict_id_ATTR_Reference_seq    = dict_id_make ("Reference_seq", 0, DTYPE_GFF3_ATTR).num;
        dict_id_ATTR_Variant_freq     = dict_id_make ("Variant_freq", 0, DTYPE_GFF3_ATTR).num;

        // fields added in the GVFs of GRCh37/38
        dict_id_ATTR_Dbxref           = dict_id_make ("Dbxref", 6, DTYPE_GFF3_ATTR).num;
        dict_id_ATTR_ancestral_allele = dict_id_make ("ancestral_allele", 0, DTYPE_GFF3_ATTR).num;
        dict_id_ATTR_Variant_effect   = dict_id_make ("Variant_effect", 0, DTYPE_GFF3_ATTR).num;
        dict_id_ATTR_sift_prediction  = dict_id_make ("sift_prediction", 0, DTYPE_GFF3_ATTR).num;
        dict_id_ATTR_polyphen_prediction = dict_id_make ("polyphen_prediction", 0, DTYPE_GFF3_ATTR).num;
        dict_id_ATTR_variant_peptide  = dict_id_make ("variant_peptide", 0, DTYPE_GFF3_ATTR).num;

        dict_id_ENSTid                = dict_id_make ("ENSTid", 0, DTYPE_GFF3_ENST).num; 
        break;

    default:
        break; // no special fields for the other data types
    }
}

static void dict_id_show_aliases (void)
{
    iprint0 ("Contents of SEC_DICT_ID_ALIASES section:\n");
    for (unsigned i=0; i < dict_id_num_aliases; i++) 
        iprintf ("alias=%s dst=%s\n", 
                 dis_dict_id (dict_id_aliases[i].alias).s, dis_dict_id (dict_id_aliases[i].dst).s);

    if (exe_type == EXE_GENOCAT) exit_ok;
}

// called by ZIP main thread for writing to global section
Buffer *dict_id_create_aliases_buf (void)
{
    static struct { DataType dt; uint64_t *dict_id_alias; uint64_t *dict_id_dst; } aliases_def[] = DICT_ID_ALIASES;

    // count
    dict_id_num_aliases = 0;
    for (unsigned i=0; i < sizeof(aliases_def)/sizeof(aliases_def[0]); i++)
        if (aliases_def[i].dt == z_file->data_type)
            dict_id_num_aliases++;

    // build global alias reference, which will be immutable until the end of this z_file
    dict_id_aliases_buf.len = dict_id_num_aliases * sizeof (DictIdAlias);
    buf_alloc_old (evb, &dict_id_aliases_buf, dict_id_aliases_buf.len, 1, "dict_id_aliases_buf");

    DictIdAlias *next = FIRSTENT (DictIdAlias, dict_id_aliases_buf);
    for (unsigned i=0; i < sizeof(aliases_def)/sizeof(aliases_def[0]); i++)
        if (aliases_def[i].dt == z_file->data_type) {
            next->alias = (DictId)*aliases_def[i].dict_id_alias;
            next->dst   = (DictId)*aliases_def[i].dict_id_dst;
            next++;
        }

    if (flag.show_aliases) dict_id_show_aliases();

    return &dict_id_aliases_buf;
}

// PIZ main thread: read all dict_id aliaeses, if there are any
void dict_id_read_aliases (void) 
{ 
    if (!sections_next_sec1 (NULL, SEC_DICT_ID_ALIASES, false, true)) return; // no aliases section

    buf_free (&dict_id_aliases_buf); // needed in case this is the 2nd+ file being pizzed

    zfile_get_global_section (SectionHeader, SEC_DICT_ID_ALIASES, NULL, &dict_id_aliases_buf, "dict_id_aliases_buf");

    dict_id_aliases = FIRSTENT (DictIdAlias, dict_id_aliases_buf);
    dict_id_num_aliases = dict_id_aliases_buf.len / sizeof (DictIdAlias);

    for (unsigned i=0; i < dict_id_num_aliases; i++) 
        ASSERT0 (dict_id_aliases[i].dst.id[0] && dict_id_aliases[i].alias.id[0], "corrupted aliases buffer");
    
    if (flag.show_aliases) dict_id_show_aliases();
}

// template can be 0 - anything OR a type - must 2 MSb of id[0] are used OR a specific dict_id
// candidate is a specific dict_id that we test for its matching of the template
bool dict_id_is_match (DictId template, DictId candidate)
{
    if (!template.num) return true;

    if (template.num == candidate.num) return true;

    // template if it is 0 except for the 2 MSb of byte 0 - i.e. we're searching for a type of dict_id rather than a specific one
    DictId copy = template;
    copy.id[0] &= 0x3f; // remove the two MSb */

    if (!copy.num && ((template.id[0] & 0xc0) == (candidate.id[0] & 0xc0)))
        return true;

    return false;
}

const char *dict_id_display_type (DataType dt, DictId dict_id)
{
    if (dict_id_is_field  (dict_id)) return dt_props[dt].stat_dict_types[0];
    if (dict_id_is_type_1 (dict_id)) return dt_props[dt].stat_dict_types[1]; 
    if (dict_id_is_type_2 (dict_id)) return dt_props[dt].stat_dict_types[2]; 

    ABORT0_R ("Error in dict_id_display_type");
}

// print the dict_id - NOT thread safe, for use in execution-termination messages
DisplayPrintId dis_dict_id_ex (DictId dict_id, bool with_type_if_vcf)
{
    DisplayPrintId s = {};

    if (!dict_id.num) return (DisplayPrintId){ .s = "<none>" };

    unsigned start=0;
    if (with_type_if_vcf && z_file && (z_file->data_type == DT_VCF || z_file->data_type == DT_BCF)) {
        if      (dict_id_is_vcf_format_sf (dict_id)) { memcpy (s.s, "FORMAT/", 7); start = 7; }
        else if (dict_id_is_vcf_info_sf   (dict_id)) { memcpy (s.s, "INFO/",   5); start = 5; }
    }

    s.s[start] = (dict_id.id[0] & 0x7f) | 0x40;  // set 2 Msb to 01

    for (unsigned i=1; i < DICT_ID_LEN; i++) 
        s.s[start+i] = dict_id.id[i] == 0  ? ' '
                     : dict_id.id[i] < 32  ? '?' // non printable char is converted to '?'
                     : dict_id.id[i] > 126 ? '?'
                     :                       dict_id.id[i];

    // trim final ' '
    for (int i=DICT_ID_LEN-1; i >= 0; i--) {
        if (s.s[start+i] != ' ') break;
        s.s[start+i] = 0;
    }

    return s;
}
