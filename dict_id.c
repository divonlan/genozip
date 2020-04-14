// ------------------------------------------------------------------
//   dict_id.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "dict_id.h"
#include "header.h"

// globals externed in dict_id.h and initialized in dict_id_initialize

// VCF stuff
uint64_t dict_id_vcf_fields[NUM_VCF_FIELDS] = {0,0,0,0,0,0,0,0},
         dict_id_FORMAT_PL=0, dict_id_FORMAT_GL=0, dict_id_FORMAT_GP=0, 
         dict_id_INFO_AC=0, dict_id_INFO_AF=0, dict_id_INFO_AN=0, dict_id_INFO_DP=0, dict_id_INFO_VQSLOD=0,
         dict_id_INFO_13=0;

// SAM stuff
uint64_t dict_id_sam_fields[NUM_SAM_FIELDS] = {0,0,0,0,0,0,0,0,0,0},
         dict_id_OPTION_AM=0, dict_id_OPTION_AS=0, dict_id_OPTION_CM=0, dict_id_OPTION_LB=0, dict_id_OPTION_FI=0, dict_id_OPTION_H0=0,
         dict_id_OPTION_H1=0, dict_id_OPTION_H2=0, dict_id_OPTION_MQ=0, dict_id_OPTION_NH=0, dict_id_OPTION_NM=0, dict_id_OPTION_PG=0, 
         dict_id_OPTION_PQ=0, dict_id_OPTION_PU=0, dict_id_OPTION_RG=0, dict_id_OPTION_SM=0, dict_id_OPTION_TC=0, dict_id_OPTION_UQ=0,
         dict_id_OPTION_CC=0, dict_id_OPTION_CG=0, dict_id_OPTION_MC=0;

static DataType last_data_type = DATA_TYPE_NONE;

void dict_id_initialize (DataType data_type) 
{   // note: this uint64_t values will be different in big and little endian machines 
    // (it's ok, they never get stored in the file)

    switch (data_type) { 
    case DATA_TYPE_VCF:
        for (VcfFields f=VCF_CHROM; f <= VCF_FORMAT; f++)
            dict_id_vcf_fields[f] = dict_id_field (dict_id_make (vcf_field_names[f], strlen (vcf_field_names[f]))).num; 
        
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

    case DATA_TYPE_SAM:
        for (SamFields f=SAM_QNAME; f <= SAM_OPTIONAL; f++)
            dict_id_sam_fields[f] = dict_id_field (dict_id_make (sam_field_names[f], strlen (sam_field_names[f]))).num; 

        dict_id_OPTION_AM    = dict_id_vcf_info_sf   (dict_id_make ("AM", 2)).num;
        dict_id_OPTION_AS    = dict_id_vcf_info_sf   (dict_id_make ("AS", 2)).num;
        dict_id_OPTION_CM    = dict_id_vcf_info_sf   (dict_id_make ("CM", 2)).num;
        dict_id_OPTION_LB    = dict_id_vcf_info_sf   (dict_id_make ("LB", 2)).num;
        dict_id_OPTION_FI    = dict_id_vcf_info_sf   (dict_id_make ("FI", 2)).num;
        dict_id_OPTION_H0    = dict_id_vcf_info_sf   (dict_id_make ("H0", 2)).num;
        dict_id_OPTION_H1    = dict_id_vcf_info_sf   (dict_id_make ("H1", 2)).num;
        dict_id_OPTION_H2    = dict_id_vcf_info_sf   (dict_id_make ("H2", 2)).num;
        dict_id_OPTION_MQ    = dict_id_vcf_info_sf   (dict_id_make ("MQ", 2)).num;
        dict_id_OPTION_NH    = dict_id_vcf_info_sf   (dict_id_make ("NH", 2)).num;
        dict_id_OPTION_NM    = dict_id_vcf_info_sf   (dict_id_make ("NM", 2)).num;
        dict_id_OPTION_PG    = dict_id_vcf_info_sf   (dict_id_make ("PG", 2)).num;
        dict_id_OPTION_PQ    = dict_id_vcf_info_sf   (dict_id_make ("PQ", 2)).num;
        dict_id_OPTION_PU    = dict_id_vcf_info_sf   (dict_id_make ("PU", 2)).num;
        dict_id_OPTION_RG    = dict_id_vcf_info_sf   (dict_id_make ("RG", 2)).num;
        dict_id_OPTION_SM    = dict_id_vcf_info_sf   (dict_id_make ("SM", 2)).num;
        dict_id_OPTION_TC    = dict_id_vcf_info_sf   (dict_id_make ("TC", 2)).num;
        dict_id_OPTION_UQ    = dict_id_vcf_info_sf   (dict_id_make ("UQ", 2)).num;
        dict_id_OPTION_CC    = dict_id_vcf_info_sf   (dict_id_make ("CC", 2)).num;
        dict_id_OPTION_CG    = dict_id_vcf_info_sf   (dict_id_make ("CG", 2)).num;
        dict_id_OPTION_MC    = dict_id_vcf_info_sf   (dict_id_make ("MC", 2)).num;

        break;

    default:
        ABORT ("Error in dict_id_initialize: unknown data_type: %d", data_type);
    }

    last_data_type = data_type;
}

const char *dict_id_display_type (DictIdType dict_id)
{
    switch (last_data_type) { 
    case DATA_TYPE_VCF:
        if (dict_id_is_field (dict_id))     return "FIELD";
        if (dict_id_is_vcf_info_sf (dict_id))   return "INFO";
        if (dict_id_is_vcf_format_sf (dict_id)) return "FORMAT";
        break;

    case DATA_TYPE_SAM:
        if (dict_id_is_field (dict_id))     return "FIELD";
        if (dict_id_is_sam_qname_sf (dict_id))  return "QNAME";
        if (dict_id_is_sam_optnl_sf (dict_id))  return "OPTION";
        break;
    
    default: break; 
    }
    return "ERROR!";
}

// returns field to which this dict_id belongs, if its a main field dict_id, or minus the next field if not
int dict_id_get_field (DictIdType dict_id)
{
    static const uint64_t *dict_id_datatype_fields[NUM_DATATYPES] = { dict_id_vcf_fields, dict_id_sam_fields };

    if (!dict_id_is_field (dict_id)) 
        return -(datatype_last_field[last_data_type] + 1); // not a main field - returning minus the next field

    for (int f=0; f <= datatype_last_field[last_data_type]; f++)
        if (dict_id.num == dict_id_datatype_fields[last_data_type][f]) return f;

    ABORT ("Error in dict_id_get_field: dict_id=%.*s is not a field dictionary despite appeared to be so", DICT_ID_LEN, dict_id.id); 
    return 0; // quieten compiler warning - never reaches here
}
