// ------------------------------------------------------------------
//   regions.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "regions.h"
#include "buffer.h"

#define REGIONS_NO_POS 0xffffffff

typedef struct {
    const char *chrom;           // NULL if no chrom
    uint32_t start_pos, end_pos; // start/end either both have a value or are both REGIONS_NO_POS (0xffffffff) if no start/end pos. 
                                 // the region searched will include both the start and the end
} Region;

static Buffer regions_buf = EMPTY_BUFFER;

typedef enum {RPT_SINGLTON, RPT_START_ONLY, RPT_END_ONLY, RPT_BOTH_START_END} RegionPosType;

static bool regions_parse_pos (const char *str, 
                               RegionPosType *type, uint32_t *start_pos, uint32_t *end_pos) // optional outs - only if case of true
{
    unsigned len = strlen (str);

    // pos needs to be one of 4 formats: N -N N- or N-N. Where N is a non-negative integer less than 0xffffffff
    // if there are two numbers - the first must be not larger than the second.
    int hyphen_found = false;
    bool digit_found  = false;
    for (unsigned i=0; i < len; i++) {
        if (str[i] == '-') {
            if (hyphen_found) goto fail; // only one hyphen is permitted
            hyphen_found = true;
        }
        else if (str[i] >= '0' && str[i] <= '9') 
            digit_found = true;
        else
            goto fail; // only digits and hyphen are allowed
    }
    if (!digit_found) goto fail; // at least one digit is required

    // this doesn't have the same power as regex in validating the format exactly, for example,
    // if the user gives a number longer than 10 digits, we will take the first 10. but that's fine -
    // i don't want to deal with finding a good cross-platform regex library just for this

    uint32_t num1, num2;
    bool hyphen_first = str[0]=='-';
    bool hyphen_last  = str[len-1]=='-';

    int count = sscanf (&str[hyphen_first], "%10u-%10u", &num1, &num2); // if first char is -, start from 2nd char
    
    if (count == 2) {
        if (start_pos > end_pos) goto fail; // a range such 100-1 is not permitted
        if (type) *type = RPT_BOTH_START_END;
        if (start_pos) { *start_pos = num1; *end_pos = num2; }
    }
    else if (hyphen_first) {
        if (type) *type = RPT_END_ONLY;
        if (start_pos) { *start_pos = 0; *end_pos = num1; }
    }
    else if (hyphen_last) {
        if (type) *type = RPT_START_ONLY;
        if (start_pos) { *start_pos = num1; *end_pos = 0xfffffffe; }
    }
    else {
        if (type) *type = RPT_SINGLTON;
        if (start_pos) { *start_pos = *end_pos = num1; }
    }

    return true;

fail:
    if (start_pos) { *start_pos = *end_pos = REGIONS_NO_POS; }
    return false;
}

static bool regions_is_valid_chrom (const char *str)
{
    // if it looks like a non-singleton range, we take it as not being a pos
    RegionPosType pos_type;
    if (regions_parse_pos (str, &pos_type, NULL, NULL) && pos_type != RPT_SINGLTON) return false;

    // per VCF 4.2 specification, the limitations on CHROM are it is not allowed to contain whitespace or a colon
    unsigned len = strlen (str);

    for (unsigned i=0; i < len; i++)
        if (str[i] < 33 || str[i] > 127 || str[i] == ':') return false;

    return true;
}


void regions_add (char *reg_str)
{
    char *next_region_token = reg_str;

    while (1) {
        char *one_reg_str = strtok_s (next_region_token, ",", &next_region_token);
        if (!one_reg_str) break;

        buf_alloc (external_vb, &regions_buf, MAX (regions_buf.len + 1, 100), 2, "regions_buf", 0);

        char *after_colon;
        char *before_colon = strtok_s (one_reg_str, ":", &after_colon);

        ASSERT (before_colon, "Error: invalid region string: %s", reg_str);

        Region *reg = &((Region *)regions_buf.data)[regions_buf.len++]; // update after possible realloc
        reg->start_pos = reg->end_pos = REGIONS_NO_POS;
        reg->chrom = NULL;

        // case: we have both chrom and pos - easy!
        if (after_colon[0]) {
            ASSERT (regions_is_valid_chrom (before_colon), "Error: Invalid CHROM in region string: %s", one_reg_str);
            reg->chrom = before_colon;
            
            ASSERT (regions_parse_pos (after_colon, NULL, &reg->start_pos, &reg->end_pos), 
                    "Error: Invalid position range in region string: %s", one_reg_str);
        }

        // case: only one substring. we need to determine if the single substring is a pos or a chrom. if it
        // is both a valid chrom and pos (i.e. it is a single number e.g. 22), then we create
        // two regions - one for the chrom and one for the pos
        else {
            if (regions_is_valid_chrom (before_colon)) 
                reg->chrom = before_colon;

            regions_parse_pos (before_colon, NULL, &reg->start_pos, &reg->end_pos);

            // make sure at least one of them is valid
            ASSERT (reg->chrom || reg->start_pos != REGIONS_NO_POS, "Error: Invalid region string: %s", one_reg_str);
        }

        printf ("region: chrom=%s start=%u end=%u\n", reg->chrom ? reg->chrom : "NONE" , reg->start_pos, reg->end_pos);
    }
}
