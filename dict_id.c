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
uint64_t dict_id_sam_fields[NUM_SAM_FIELDS] = {0,0,0,0,0,0,0,0,0},
         dict_id_OPTION_AM=0, dict_id_OPTION_AS=0, dict_id_OPTION_CM=0, dict_id_OPTION_LB=0, dict_id_OPTION_FI=0, dict_id_OPTION_H0=0,
         dict_id_OPTION_H1=0, dict_id_OPTION_H2=0, dict_id_OPTION_MQ=0, dict_id_OPTION_NH=0, dict_id_OPTION_NM=0, dict_id_OPTION_OA=0,
         dict_id_OPTION_OC=0, dict_id_OPTION_PG=0, dict_id_OPTION_E2=0, dict_id_OPTION_U2=0,
         dict_id_OPTION_PQ=0, dict_id_OPTION_PU=0, dict_id_OPTION_RG=0, dict_id_OPTION_SA=0, dict_id_OPTION_SM=0, dict_id_OPTION_TC=0, 
         dict_id_OPTION_UQ=0, dict_id_OPTION_CC=0, dict_id_OPTION_CG=0, dict_id_OPTION_MC=0,
         dict_id_OPTION_X0=0, dict_id_OPTION_X1=0, dict_id_OPTION_XA=0, dict_id_OPTION_XN=0, dict_id_OPTION_XM=0, dict_id_OPTION_XO=0,
         dict_id_OPTION_XG=0, dict_id_OPTION_XS=0, dict_id_OPTION_XE=0,
         dict_id_OPTION_mc=0,
         dict_id_OPTION_STRAND=0;
          
static DataType last_data_type = DATA_TYPE_NONE;

DictIdType DICT_ID_NONE = {0};

void dict_id_initialize (DataType data_type) 
{   // note: this uint64_t values will be different in big and little endian machines 
    // (it's ok, they never get stored in the file)

    switch (data_type) { 
    case DATA_TYPE_VCF:
        for (VcfFields f=VCF_CHROM; f <= VCF_FORMAT; f++)
            dict_id_vcf_fields[f] = dict_id_field (dict_id_make (field_names[DATA_TYPE_VCF][f], strlen (field_names[DATA_TYPE_VCF][f]))).num; 
        
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
            dict_id_sam_fields[f] = dict_id_field (dict_id_make (field_names[DATA_TYPE_SAM][f], strlen (field_names[DATA_TYPE_SAM][f]))).num; 

        dict_id_OPTION_AM = dict_id_sam_optnl_sf (dict_id_make ("AM", 2)).num;
        dict_id_OPTION_AS = dict_id_sam_optnl_sf (dict_id_make ("AS", 2)).num;
        dict_id_OPTION_CC = dict_id_sam_optnl_sf (dict_id_make ("CC", 2)).num;
        dict_id_OPTION_CG = dict_id_sam_optnl_sf (dict_id_make ("CG", 2)).num;
        dict_id_OPTION_CM = dict_id_sam_optnl_sf (dict_id_make ("CM", 2)).num;
        dict_id_OPTION_E2 = dict_id_sam_optnl_sf (dict_id_make ("E2", 2)).num;
        dict_id_OPTION_FI = dict_id_sam_optnl_sf (dict_id_make ("FI", 2)).num;
        dict_id_OPTION_H0 = dict_id_sam_optnl_sf (dict_id_make ("H0", 2)).num;
        dict_id_OPTION_H1 = dict_id_sam_optnl_sf (dict_id_make ("H1", 2)).num;
        dict_id_OPTION_H2 = dict_id_sam_optnl_sf (dict_id_make ("H2", 2)).num;
        dict_id_OPTION_LB = dict_id_sam_optnl_sf (dict_id_make ("LB", 2)).num;
        dict_id_OPTION_MC = dict_id_sam_optnl_sf (dict_id_make ("MC", 2)).num;
        dict_id_OPTION_MQ = dict_id_sam_optnl_sf (dict_id_make ("MQ", 2)).num;
        dict_id_OPTION_NH = dict_id_sam_optnl_sf (dict_id_make ("NH", 2)).num;
        dict_id_OPTION_NM = dict_id_sam_optnl_sf (dict_id_make ("NM", 2)).num;
        dict_id_OPTION_OA = dict_id_sam_optnl_sf (dict_id_make ("OA", 2)).num;
        dict_id_OPTION_OC = dict_id_sam_optnl_sf (dict_id_make ("OC", 2)).num;
        dict_id_OPTION_PG = dict_id_sam_optnl_sf (dict_id_make ("PG", 2)).num;
        dict_id_OPTION_PQ = dict_id_sam_optnl_sf (dict_id_make ("PQ", 2)).num;
        dict_id_OPTION_PU = dict_id_sam_optnl_sf (dict_id_make ("PU", 2)).num;
        dict_id_OPTION_RG = dict_id_sam_optnl_sf (dict_id_make ("RG", 2)).num;
        dict_id_OPTION_SA = dict_id_sam_optnl_sf (dict_id_make ("SA", 2)).num;
        dict_id_OPTION_SM = dict_id_sam_optnl_sf (dict_id_make ("SM", 2)).num;
        dict_id_OPTION_TC = dict_id_sam_optnl_sf (dict_id_make ("TC", 2)).num;
        dict_id_OPTION_UQ = dict_id_sam_optnl_sf (dict_id_make ("UQ", 2)).num;
        dict_id_OPTION_U2 = dict_id_sam_optnl_sf (dict_id_make ("U2", 2)).num;
                
        dict_id_OPTION_X0 = dict_id_sam_optnl_sf (dict_id_make ("X0", 2)).num; 
        dict_id_OPTION_X1 = dict_id_sam_optnl_sf (dict_id_make ("X1", 2)).num; 
        dict_id_OPTION_XA = dict_id_sam_optnl_sf (dict_id_make ("XA", 2)).num; 
        dict_id_OPTION_XN = dict_id_sam_optnl_sf (dict_id_make ("XN", 2)).num; 
        dict_id_OPTION_XM = dict_id_sam_optnl_sf (dict_id_make ("XM", 2)).num; 
        dict_id_OPTION_XO = dict_id_sam_optnl_sf (dict_id_make ("XO", 2)).num;
        dict_id_OPTION_XG = dict_id_sam_optnl_sf (dict_id_make ("XG", 2)).num; 
        dict_id_OPTION_XS = dict_id_sam_optnl_sf (dict_id_make ("XS", 2)).num; 
        dict_id_OPTION_XE = dict_id_sam_optnl_sf (dict_id_make ("XE", 2)).num;

        dict_id_OPTION_mc = dict_id_sam_optnl_sf (dict_id_make ("mc", 2)).num;

        dict_id_OPTION_STRAND = dict_id_sam_optnl_sf (dict_id_make ("STRAND", 6)).num;

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

