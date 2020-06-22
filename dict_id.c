// ------------------------------------------------------------------
//   dict_id.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "dict_id.h"
#include "data_types.h"
#include "file.h"
#include "zfile.h"
#include "sections.h"

// globals externed in dict_id.h and initialized in dict_id_initialize
static Buffer dict_id_aliases_buf  = EMPTY_BUFFER;
const DictIdAlias *dict_id_aliases = NULL;
uint32_t dict_id_num_aliases = 0;

uint64_t dict_id_fields[MAX_NUM_FIELDS_PER_DATA_TYPE];

// VCF stuff
uint64_t dict_id_FORMAT_PL=0, dict_id_FORMAT_GL=0, dict_id_FORMAT_GP=0, dict_id_FORMAT_DP=0, dict_id_FORMAT_MIN_DP=0, 
         dict_id_INFO_AC=0, dict_id_INFO_AF=0, dict_id_INFO_AN=0, dict_id_INFO_DP=0, dict_id_INFO_VQSLOD=0,
         dict_id_INFO_END=0;

// SAM stuff
uint64_t dict_id_OPTION_AM=0, dict_id_OPTION_AS=0, dict_id_OPTION_CM=0, dict_id_OPTION_LB=0, dict_id_OPTION_FI=0, dict_id_OPTION_H0=0,
         dict_id_OPTION_H1=0, dict_id_OPTION_H2=0, dict_id_OPTION_MD=0, dict_id_OPTION_MQ=0, dict_id_OPTION_NH=0, dict_id_OPTION_NM=0, 
         dict_id_OPTION_OA=0, dict_id_OPTION_OC=0, dict_id_OPTION_PG=0, dict_id_OPTION_E2=0, dict_id_OPTION_U2=0,
         dict_id_OPTION_PQ=0, dict_id_OPTION_PU=0, dict_id_OPTION_RG=0, dict_id_OPTION_SA=0, dict_id_OPTION_SM=0, dict_id_OPTION_TC=0, 
         dict_id_OPTION_UQ=0, dict_id_OPTION_CC=0, dict_id_OPTION_MC=0,
         dict_id_OPTION_X0=0, dict_id_OPTION_X1=0, dict_id_OPTION_XA=0, dict_id_OPTION_XN=0, dict_id_OPTION_XM=0, dict_id_OPTION_XO=0,
         dict_id_OPTION_XG=0, dict_id_OPTION_XS=0, dict_id_OPTION_XE=0,
         dict_id_OPTION_mc=0, dict_id_OPTION_ms=0,
         dict_id_OPTION_BD=0, dict_id_OPTION_BI=0,
         dict_id_OPTION_ZM=0,
  
         // private genozip dict
         dict_id_OPTION_STRAND=0, dict_id_OPTION_RNAME=0, dict_id_OPTION_POS=0, dict_id_OPTION_CIGAR=0, dict_id_OPTION_MAPQ=0,
         dict_id_SAM_SQnonref=0; 

// FASTA stuff
uint64_t dict_id_FASTA_DESC=0, dict_id_FASTA_SEQ=0, dict_id_FASTA_COMMENT=0;

// GVF stuff
uint64_t dict_id_ATTR_ID=0, dict_id_ATTR_Variant_seq=0, dict_id_ATTR_Reference_seq=0, dict_id_ATTR_Variant_freq=0,
         dict_id_ATTR_Dbxref=0, // from from GRCh37/38 - example: "dbSNP_151:rs1282280967"
         dict_id_ATTR_ancestral_allele=0, // from from GRCh37/38 - example ancestral_allele=GTTA
         dict_id_ATTR_Variant_effect=0, // example: "Variant_effect=non_coding_transcript_variant 0 ncRNA ENST00000431238,intron_variant 0 primary_transcript ENST00000431238"
         dict_id_ATTR_sift_prediction=0, dict_id_ATTR_polyphen_prediction=0, dict_id_ATTR_variant_peptide=0,

         dict_id_ENSTid=0; // private genozip dict

// our stuff used in multiple data types
uint64_t dict_id_WindowsEOL=0;         

DictId dict_id_make(const char *str, unsigned str_len) 
{ 
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

    return dict_id;
}

void dict_id_initialize (DataType data_type) 
{   
    ASSERT0 (data_type != DT_NONE, "Error in dict_id_initialize: data_type is DT_NONE");

    for (int f=0; f < dt_fields[data_type].num_fields; f++) {
        const char *field_name = dt_fields[data_type].names[f];
        dict_id_fields[f] = dict_id_field (dict_id_make (field_name, strlen (field_name))).num; 
    }

    dict_id_WindowsEOL = dict_id_type_1 (dict_id_make ("#", 1)).num; 

    switch (data_type) { 
    case DT_VCF:
        dict_id_FORMAT_PL     = dict_id_vcf_format_sf (dict_id_make ("PL", 2)).num;
        dict_id_FORMAT_GP     = dict_id_vcf_format_sf (dict_id_make ("GP", 2)).num;
        dict_id_FORMAT_GL     = dict_id_vcf_format_sf (dict_id_make ("GL", 2)).num;
        dict_id_FORMAT_DP     = dict_id_vcf_format_sf (dict_id_make ("DP", 2)).num;
        
        dict_id_INFO_AC       = dict_id_vcf_info_sf   (dict_id_make ("AC", 2)).num;
        dict_id_INFO_AF       = dict_id_vcf_info_sf   (dict_id_make ("AF", 2)).num;
        dict_id_INFO_AN       = dict_id_vcf_info_sf   (dict_id_make ("AN", 2)).num;
        dict_id_INFO_DP       = dict_id_vcf_info_sf   (dict_id_make ("DP", 2)).num;
        dict_id_INFO_VQSLOD   = dict_id_vcf_info_sf   (dict_id_make ("VQSLOD", 6)).num;

        // Added by GATK HaplotypeCaller in a gVCF: https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format
        dict_id_INFO_END      = dict_id_vcf_info_sf   (dict_id_make ("END", 3)).num;
        dict_id_FORMAT_MIN_DP = dict_id_vcf_format_sf (dict_id_make ("MIN_DP", 6)).num;

        // This appears if the VCF line has a Windows-style \r\n line ending
        break;

    case DT_SAM:
        dict_id_OPTION_AM = sam_dict_id_optnl_sf (dict_id_make ("AM:i", 4)).num;
        dict_id_OPTION_AS = sam_dict_id_optnl_sf (dict_id_make ("AS:i", 4)).num;
        dict_id_OPTION_CC = sam_dict_id_optnl_sf (dict_id_make ("CC:Z", 4)).num;
        dict_id_OPTION_BD = sam_dict_id_optnl_sf (dict_id_make ("BD:Z", 4)).num;
        dict_id_OPTION_BI = sam_dict_id_optnl_sf (dict_id_make ("BI:Z", 4)).num;
        dict_id_OPTION_CM = sam_dict_id_optnl_sf (dict_id_make ("CM:i", 4)).num;
        dict_id_OPTION_E2 = sam_dict_id_optnl_sf (dict_id_make ("E2:Z", 4)).num;
        dict_id_OPTION_FI = sam_dict_id_optnl_sf (dict_id_make ("FI:i", 4)).num;
        dict_id_OPTION_H0 = sam_dict_id_optnl_sf (dict_id_make ("H0:i", 4)).num;
        dict_id_OPTION_H1 = sam_dict_id_optnl_sf (dict_id_make ("H1:i", 4)).num;
        dict_id_OPTION_H2 = sam_dict_id_optnl_sf (dict_id_make ("H2:i", 4)).num;
        dict_id_OPTION_LB = sam_dict_id_optnl_sf (dict_id_make ("LB:Z", 4)).num;
        dict_id_OPTION_MC = sam_dict_id_optnl_sf (dict_id_make ("MC:Z", 4)).num;
        dict_id_OPTION_MD = sam_dict_id_optnl_sf (dict_id_make ("MD:Z", 4)).num;
        dict_id_OPTION_MQ = sam_dict_id_optnl_sf (dict_id_make ("MQ:i", 4)).num;
        dict_id_OPTION_NH = sam_dict_id_optnl_sf (dict_id_make ("NH:i", 4)).num;
        dict_id_OPTION_NM = sam_dict_id_optnl_sf (dict_id_make ("NM:i", 4)).num;
        dict_id_OPTION_OA = sam_dict_id_optnl_sf (dict_id_make ("OA:Z", 4)).num;
        dict_id_OPTION_OC = sam_dict_id_optnl_sf (dict_id_make ("OC:Z", 4)).num;
        dict_id_OPTION_PG = sam_dict_id_optnl_sf (dict_id_make ("PG:Z", 4)).num;
        dict_id_OPTION_PQ = sam_dict_id_optnl_sf (dict_id_make ("PQ:i", 4)).num;
        dict_id_OPTION_PU = sam_dict_id_optnl_sf (dict_id_make ("PU:Z", 4)).num;
        dict_id_OPTION_RG = sam_dict_id_optnl_sf (dict_id_make ("RG:Z", 4)).num;
        dict_id_OPTION_SA = sam_dict_id_optnl_sf (dict_id_make ("SA:Z", 4)).num;
        dict_id_OPTION_SM = sam_dict_id_optnl_sf (dict_id_make ("SM:i", 4)).num;
        dict_id_OPTION_TC = sam_dict_id_optnl_sf (dict_id_make ("TC:i", 4)).num;
        dict_id_OPTION_UQ = sam_dict_id_optnl_sf (dict_id_make ("UQ:i", 4)).num;
        dict_id_OPTION_U2 = sam_dict_id_optnl_sf (dict_id_make ("U2:Z", 4)).num;
                
        // Ion Torrent flow signal array
        dict_id_OPTION_ZM = sam_dict_id_optnl_sf (dict_id_make ("ZM:B", 4)).num;

        // bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
        dict_id_OPTION_X0 = sam_dict_id_optnl_sf (dict_id_make ("X0:i", 4)).num; 
        dict_id_OPTION_X1 = sam_dict_id_optnl_sf (dict_id_make ("X1:i", 4)).num; 
        dict_id_OPTION_XA = sam_dict_id_optnl_sf (dict_id_make ("XA:Z", 4)).num; 
        dict_id_OPTION_XN = sam_dict_id_optnl_sf (dict_id_make ("XN:i", 4)).num; 
        dict_id_OPTION_XM = sam_dict_id_optnl_sf (dict_id_make ("XM:i", 4)).num; 
        dict_id_OPTION_XO = sam_dict_id_optnl_sf (dict_id_make ("XO:i", 4)).num;
        dict_id_OPTION_XG = sam_dict_id_optnl_sf (dict_id_make ("XG:i", 4)).num; 
        dict_id_OPTION_XS = sam_dict_id_optnl_sf (dict_id_make ("XS:i", 4)).num; 
        dict_id_OPTION_XE = sam_dict_id_optnl_sf (dict_id_make ("XE:i", 4)).num;

        // biobambam tags
        dict_id_OPTION_mc = sam_dict_id_optnl_sf (dict_id_make ("mc:i", 4)).num;
        dict_id_OPTION_ms = sam_dict_id_optnl_sf (dict_id_make ("ms:i", 4)).num;

        // added by GATK's BQSR (Base Quality Score Recalibration)
        dict_id_OPTION_BD = sam_dict_id_optnl_sf (dict_id_make ("BD:Z", 4)).num; // not used in newer versions of GATK
        dict_id_OPTION_BI = sam_dict_id_optnl_sf (dict_id_make ("BI:Z", 4)).num; // not used in newer versions of GATK

        // our private dictionary for + or 0 strands
        dict_id_SAM_SQnonref  = dict_id_field        (dict_id_make ("SQnonref",8)).num;
        dict_id_OPTION_STRAND = sam_dict_id_optnl_sf (dict_id_make ("@STRAND", 7)).num;
        dict_id_OPTION_RNAME  = sam_dict_id_optnl_sf (dict_id_make ("@RNAME",  6)).num;
        dict_id_OPTION_POS    = sam_dict_id_optnl_sf (dict_id_make ("@POS",    4)).num;
        dict_id_OPTION_CIGAR  = sam_dict_id_optnl_sf (dict_id_make ("@CIGAR",  6)).num;
        dict_id_OPTION_MAPQ   = sam_dict_id_optnl_sf (dict_id_make ("@MAPQ",   5)).num;

        break;

    case DT_FASTA:
        dict_id_FASTA_DESC    = dict_id_field (dict_id_make ("DESC",    4)).num;
        dict_id_FASTA_SEQ     = dict_id_field (dict_id_make ("SEQ",     3)).num;
        dict_id_FASTA_COMMENT = dict_id_field (dict_id_make ("COMMENT", 7)).num;
        break;

    case DT_GFF3:
        // standard GVF fields (ID is also a standard GFF3 field)
        dict_id_ATTR_ID               = dict_id_gff3_attr_sf (dict_id_make ("ID", 2)).num;
        dict_id_ATTR_Variant_seq      = dict_id_gff3_attr_sf (dict_id_make ("Variant_seq", 0)).num;
        dict_id_ATTR_Reference_seq    = dict_id_gff3_attr_sf (dict_id_make ("Reference_seq", 0)).num;
        dict_id_ATTR_Variant_freq     = dict_id_gff3_attr_sf (dict_id_make ("Variant_freq", 0)).num;

        // fields added in the GVFs of GRCh37/38
        dict_id_ATTR_Dbxref           = dict_id_gff3_attr_sf (dict_id_make ("Dbxref", 6)).num;
        dict_id_ATTR_ancestral_allele = dict_id_gff3_attr_sf (dict_id_make ("ancestral_allele", 0)).num;
        dict_id_ATTR_Variant_effect   = dict_id_gff3_attr_sf (dict_id_make ("Variant_effect", 0)).num;
        dict_id_ATTR_sift_prediction  = dict_id_gff3_attr_sf (dict_id_make ("sift_prediction", 0)).num;
        dict_id_ATTR_polyphen_prediction = dict_id_gff3_attr_sf (dict_id_make ("polyphen_prediction", 0)).num;
        dict_id_ATTR_variant_peptide  = dict_id_gff3_attr_sf (dict_id_make ("variant_peptide", 0)).num;

        dict_id_ENSTid                = dict_id_type_2 (dict_id_make ("ENSTid", 0)).num; // type 2 is ENST subsubfields
        break;

    default:
        break; // no special fields for the other data types
    }
}

// called by ZIP I/O thread for writing to global section
Buffer *dict_id_create_aliases_buf (void)
{
    static struct { DataType dt; uint64_t *dict_id_alias; uint64_t *dict_id_dst; } aliases_def[] = DICT_ID_ALIASES;

    // count
    dict_id_aliases_buf.len = 0;
    for (unsigned i=0; i < sizeof(aliases_def)/sizeof(aliases_def[0]); i++)
        if (aliases_def[i].dt == z_file->data_type)
            dict_id_aliases_buf.len += sizeof (DictIdAlias);

    // build global alias reference, which will be immutable until the end of this z_file
    buf_alloc (evb, &dict_id_aliases_buf, dict_id_aliases_buf.len, 1, "dict_id_aliases_buf", 0);

    DictIdAlias *next = FIRSTENT (DictIdAlias, dict_id_aliases_buf);
    for (unsigned i=0; i < sizeof(aliases_def)/sizeof(aliases_def[0]); i++)
        if (aliases_def[i].dt == z_file->data_type) {
            next->alias = (DictId)*aliases_def[i].dict_id_alias;
            next->dst   = (DictId)*aliases_def[i].dict_id_dst;
            next++;
        }

    return &dict_id_aliases_buf;
}

// PIZ I/O thread: read all dict_id aliaeses, if there are any
void dict_id_read_aliases (void) 
{ 
    if (!sections_seek_to (SEC_DICT_ID_ALIASES)) return; // no aliases section

    static Buffer compressed_aliases = EMPTY_BUFFER;

    buf_free (&dict_id_aliases_buf); // needed in case this is the 2nd+ file being pizzed
    
    zfile_read_section (evb, 0, NO_SB_I, &dict_id_aliases_buf, "dict_id_aliases_buf", 
                        sizeof(SectionHeader), SEC_DICT_ID_ALIASES, NULL);    

    SectionHeader *header = (SectionHeader *)dict_id_aliases_buf.data;
    zfile_uncompress_section (evb, header, &dict_id_aliases_buf, "dict_id_aliases_buf", SEC_DICT_ID_ALIASES);

    buf_destroy (&compressed_aliases);

    dict_id_aliases = FIRSTENT (DictIdAlias, dict_id_aliases_buf);
    dict_id_num_aliases = dict_id_aliases_buf.len / sizeof (DictIdAlias);
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

    ABORT0 ("Error in dict_id_display_type");
    return 0;    
}

// print the dict_id - NOT thread safe, for use in execution-termination messages
const char *err_dict_id (DictId dict_id)
{
    static char s[DICT_ID_LEN+1];
    sprintf (s, "%.*s", DICT_ID_LEN, dict_id.num ? (char*)dict_id_printable(dict_id).id : "<null>");
    return s;
}
