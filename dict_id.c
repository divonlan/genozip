// ------------------------------------------------------------------
//   dict_id.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "dict_id.h"
#include "header.h"
#include "file.h"

// globals externed in dict_id.h and initialized in dict_id_initialize

uint64_t dict_id_fields[MAX_NUM_FIELDS_PER_DATA_TYPE];

// VCF stuff
uint64_t dict_id_FORMAT_PL=0, dict_id_FORMAT_GL=0, dict_id_FORMAT_GP=0, 
         dict_id_INFO_AC=0, dict_id_INFO_AF=0, dict_id_INFO_AN=0, dict_id_INFO_DP=0, dict_id_INFO_VQSLOD=0,
         dict_id_INFO_13=0;

// SAM stuff
uint64_t dict_id_OPTION_AM=0, dict_id_OPTION_AS=0, dict_id_OPTION_CM=0, dict_id_OPTION_LB=0, dict_id_OPTION_FI=0, dict_id_OPTION_H0=0,
         dict_id_OPTION_H1=0, dict_id_OPTION_H2=0, dict_id_OPTION_MD=0, dict_id_OPTION_MQ=0, dict_id_OPTION_NH=0, dict_id_OPTION_NM=0, 
         dict_id_OPTION_OA=0, dict_id_OPTION_OC=0, dict_id_OPTION_PG=0, dict_id_OPTION_E2=0, dict_id_OPTION_U2=0,
         dict_id_OPTION_PQ=0, dict_id_OPTION_PU=0, dict_id_OPTION_RG=0, dict_id_OPTION_SA=0, dict_id_OPTION_SM=0, dict_id_OPTION_TC=0, 
         dict_id_OPTION_UQ=0, dict_id_OPTION_CC=0, dict_id_OPTION_MC=0,
         dict_id_OPTION_X0=0, dict_id_OPTION_X1=0, dict_id_OPTION_XA=0, dict_id_OPTION_XN=0, dict_id_OPTION_XM=0, dict_id_OPTION_XO=0,
         dict_id_OPTION_XG=0, dict_id_OPTION_XS=0, dict_id_OPTION_XE=0,
         dict_id_OPTION_mc=0, dict_id_OPTION_ms=0,
         dict_id_OPTION_STRAND=0;

DictIdType DICT_ID_NONE = {0};

void dict_id_initialize (void) 
{   // note: this uint64_t values will be different in big and little endian machines 
    // (it's ok, they never get stored in the file)

    for (int f=0; f <= datatype_last_field[z_file->data_type]; f++) {
        const char *field_name = field_names[z_file->data_type][f];
        dict_id_fields[f] = dict_id_field (dict_id_make (field_name, strlen (field_name))).num; 
    }

    switch (z_file->data_type) { 
    case DT_VCF:
        dict_id_FORMAT_PL   = dict_id_vcf_format_sf (dict_id_make ("PL", 2)).num;
        dict_id_FORMAT_GP   = dict_id_vcf_format_sf (dict_id_make ("GP", 2)).num;
        dict_id_FORMAT_GL   = dict_id_vcf_format_sf (dict_id_make ("GL", 2)).num;
        
        dict_id_INFO_AC     = dict_id_vcf_info_sf   (dict_id_make ("AC", 2)).num;
        dict_id_INFO_AF     = dict_id_vcf_info_sf   (dict_id_make ("AF", 2)).num;
        dict_id_INFO_AN     = dict_id_vcf_info_sf   (dict_id_make ("AN", 2)).num;
        dict_id_INFO_DP     = dict_id_vcf_info_sf   (dict_id_make ("DP", 2)).num;
        dict_id_INFO_VQSLOD = dict_id_vcf_info_sf   (dict_id_make ("VQSLOD", 6)).num;

        dict_id_INFO_13     = dict_id_vcf_info_sf   (dict_id_make ("#", 1)).num; // This appears if the VCF line has a Windows-style \r\n line ending
        break;

    case DT_SAM:
        dict_id_OPTION_AM = dict_id_sam_optnl_sf (dict_id_make ("AM:i", 4)).num;
        dict_id_OPTION_AS = dict_id_sam_optnl_sf (dict_id_make ("AS:i", 4)).num;
        dict_id_OPTION_CC = dict_id_sam_optnl_sf (dict_id_make ("CC:Z", 4)).num;
        dict_id_OPTION_CM = dict_id_sam_optnl_sf (dict_id_make ("CM:i", 4)).num;
        dict_id_OPTION_E2 = dict_id_sam_optnl_sf (dict_id_make ("E2:Z", 4)).num;
        dict_id_OPTION_FI = dict_id_sam_optnl_sf (dict_id_make ("FI:i", 4)).num;
        dict_id_OPTION_H0 = dict_id_sam_optnl_sf (dict_id_make ("H0:i", 4)).num;
        dict_id_OPTION_H1 = dict_id_sam_optnl_sf (dict_id_make ("H1:i", 4)).num;
        dict_id_OPTION_H2 = dict_id_sam_optnl_sf (dict_id_make ("H2:i", 4)).num;
        dict_id_OPTION_LB = dict_id_sam_optnl_sf (dict_id_make ("LB:Z", 4)).num;
        dict_id_OPTION_MC = dict_id_sam_optnl_sf (dict_id_make ("MC:Z", 4)).num;
        dict_id_OPTION_MD = dict_id_sam_optnl_sf (dict_id_make ("MD:Z", 4)).num;
        dict_id_OPTION_MQ = dict_id_sam_optnl_sf (dict_id_make ("MQ:i", 4)).num;
        dict_id_OPTION_NH = dict_id_sam_optnl_sf (dict_id_make ("NH:i", 4)).num;
        dict_id_OPTION_NM = dict_id_sam_optnl_sf (dict_id_make ("NM:i", 4)).num;
        dict_id_OPTION_OA = dict_id_sam_optnl_sf (dict_id_make ("OA:Z", 4)).num;
        dict_id_OPTION_OC = dict_id_sam_optnl_sf (dict_id_make ("OC:Z", 4)).num;
        dict_id_OPTION_PG = dict_id_sam_optnl_sf (dict_id_make ("PG:Z", 4)).num;
        dict_id_OPTION_PQ = dict_id_sam_optnl_sf (dict_id_make ("PQ:i", 4)).num;
        dict_id_OPTION_PU = dict_id_sam_optnl_sf (dict_id_make ("PU:Z", 4)).num;
        dict_id_OPTION_RG = dict_id_sam_optnl_sf (dict_id_make ("RG:Z", 4)).num;
        dict_id_OPTION_SA = dict_id_sam_optnl_sf (dict_id_make ("SA:Z", 4)).num;
        dict_id_OPTION_SM = dict_id_sam_optnl_sf (dict_id_make ("SM:i", 4)).num;
        dict_id_OPTION_TC = dict_id_sam_optnl_sf (dict_id_make ("TC:i", 4)).num;
        dict_id_OPTION_UQ = dict_id_sam_optnl_sf (dict_id_make ("UQ:i", 4)).num;
        dict_id_OPTION_U2 = dict_id_sam_optnl_sf (dict_id_make ("U2:Z", 4)).num;
                
        // bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
        dict_id_OPTION_X0 = dict_id_sam_optnl_sf (dict_id_make ("X0:i", 4)).num; 
        dict_id_OPTION_X1 = dict_id_sam_optnl_sf (dict_id_make ("X1:i", 4)).num; 
        dict_id_OPTION_XA = dict_id_sam_optnl_sf (dict_id_make ("XA:Z", 4)).num; 
        dict_id_OPTION_XN = dict_id_sam_optnl_sf (dict_id_make ("XN:i", 4)).num; 
        dict_id_OPTION_XM = dict_id_sam_optnl_sf (dict_id_make ("XM:i", 4)).num; 
        dict_id_OPTION_XO = dict_id_sam_optnl_sf (dict_id_make ("XO:i", 4)).num;
        dict_id_OPTION_XG = dict_id_sam_optnl_sf (dict_id_make ("XG:i", 4)).num; 
        dict_id_OPTION_XS = dict_id_sam_optnl_sf (dict_id_make ("XS:i", 4)).num; 
        dict_id_OPTION_XE = dict_id_sam_optnl_sf (dict_id_make ("XE:i", 4)).num;

        // biobambam tags
        dict_id_OPTION_mc = dict_id_sam_optnl_sf (dict_id_make ("mc:i", 4)).num;
        dict_id_OPTION_ms = dict_id_sam_optnl_sf (dict_id_make ("ms:i", 4)).num;

        dict_id_OPTION_STRAND = dict_id_sam_optnl_sf (dict_id_make ("STRAND", 6)).num;

        break;

    default:
        break; // no special fields for the other data types
    }
}

const char *dict_id_display_type (DictIdType dict_id)
{
    static const char *dict_type_by_data_type[NUM_DATATYPES][3] = STAT_DICT_TYPES;

    if (dict_id_is_field (dict_id))  return dict_type_by_data_type[z_file->data_type][0]; 
    if (dict_id_is_type_1 (dict_id)) return dict_type_by_data_type[z_file->data_type][1]; 
    if (dict_id_is_type_2 (dict_id)) return dict_type_by_data_type[z_file->data_type][2]; 

    return "Bug!";
}

// print the dict_id - NOT thread safe, for use in execution-termination messages
const char *err_dict_id (DictIdType dict_id)
{
    static char s[DICT_ID_LEN+1];
    sprintf (s, "%.*s", DICT_ID_LEN, dict_id_printable(dict_id).id);
    return s;
}
