// ------------------------------------------------------------------
//   dict_id.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "dict_id.h"
#include "header.h"

// globals externed in dict_id.h and initialized in dict_id_initialize
uint64_t dict_id_vardata_fields[8] = {0,0,0,0,0,0,0,0},
         dict_id_FORMAT_PL=0, dict_id_FORMAT_GL=0, dict_id_FORMAT_GP=0, 
         dict_id_INFO_AC=0, dict_id_INFO_AF=0, dict_id_INFO_AN=0, dict_id_INFO_DP=0, dict_id_INFO_VQSLOD=0,
         dict_id_INFO_13=0;

void dict_id_initialize(void) 
{   // note: this uint64_t values will be different in big and little endian machines 
    // (it's ok, they never get stored in the file)
    for (VcfFields f=VCF_CHROM; f <= VCF_FORMAT; f++)
        dict_id_vardata_fields[f] = dict_id_vardata_field (dict_id_make (vcf_field_names[f], strlen (vcf_field_names[f]))).num; 
    
    dict_id_FORMAT_PL   = dict_id_format_subfield (dict_id_make ("PL", 2)).num;
    dict_id_FORMAT_GP   = dict_id_format_subfield (dict_id_make ("GP", 2)).num;
    dict_id_FORMAT_GL   = dict_id_format_subfield (dict_id_make ("GL", 2)).num;
    
    dict_id_INFO_AC     = dict_id_info_subfield   (dict_id_make ("AC", 2)).num;
    dict_id_INFO_AF     = dict_id_info_subfield   (dict_id_make ("AF", 2)).num;
    dict_id_INFO_AN     = dict_id_info_subfield   (dict_id_make ("AN", 2)).num;
    dict_id_INFO_DP     = dict_id_info_subfield   (dict_id_make ("DP", 2)).num;
    dict_id_INFO_VQSLOD = dict_id_info_subfield   (dict_id_make ("VQSLOD", 6)).num;

    dict_id_INFO_13     = dict_id_info_subfield   (dict_id_make ("#", 1)).num; // This appears if the VCF line has a Windows-style \r\n line ending
}

const char *dict_id_display_type (DictIdType dict_id)
{
    if (dict_id_is_format_subfield (dict_id)) return "FORMAT";
    if (dict_id_is_info_subfield (dict_id)) return "INFO";
    if (dict_id_is_vardata_field (dict_id)) return "FIELD";
    return "ERROR!n";
}

// returns field to which this dict_id belongs, if its a vardata dict_id, or -1 if not
int dict_id_get_field (DictIdType dict_id)
{
    for (VcfFields f=VCF_CHROM; f <= VCF_FORMAT; f++)
        if (dict_id.num == dict_id_vardata_fields[f]) return f;

    return -1;
}