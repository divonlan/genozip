// ------------------------------------------------------------------
//   regions.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "regions.h"
#include "buffer.h"
#include "move_to_front.h"
#include "vcf_header.h"
#include "vb.h"
#include "file.h"

typedef struct {
    const char *chrom;   // NULL if no chrom
    int32_t chrom_node_index; // NIL if chrom is NULL
    uint32_t start_pos;  // if the user did specify pos then start_pos=0 and end_pos=0xffffffff
    uint32_t end_pos;    // the region searched will include both the start and the end
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
        if (type) *type = RPT_BOTH_START_END;
        if (num1 > num2) goto fail;   // a range such 100-1 is not permitted
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
    if (start_pos) { *start_pos = 0; *end_pos = 0xffffffff; }
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


void regions_add (const char *region_str)
{
    ASSERT0 (region_str, "Error: region_str is NULL");

    // make a copy of the string and leave the original one for error message. 
    // we don't free this memory as chrom fields in regions will be pointing to it
    char *rs = malloc (strlen (region_str)+1); 
    strcpy (rs, region_str);

    char *next_region_token = rs;

    while (1) {
        char *one_rs = strtok_s (next_region_token, ",", &next_region_token);
        if (!one_rs) break;

        buf_alloc (evb, &regions_buf, MAX (regions_buf.len + 1, 100), 2, "regions_buf", 0);

        char *after_colon;
        char *before_colon = strtok_s (one_rs, ":", &after_colon);

        ASSERT (before_colon, "Error: invalid region string: %s", region_str);

        Region *reg = &((Region *)regions_buf.data)[regions_buf.len++]; // update after possible realloc
        Region *reg2 = NULL;

        reg->start_pos = 0;
        reg->end_pos = 0xffffffff;
        reg->chrom = NULL;
        reg->chrom_node_index = NIL;

        // case: we have both chrom and pos - easy!
        if (after_colon[0]) {
            ASSERT (regions_is_valid_chrom (before_colon), "Error: Invalid CHROM in region string: %s", region_str);
            reg->chrom = before_colon;
            
            ASSERT (regions_parse_pos (after_colon, NULL, &reg->start_pos, &reg->end_pos), 
                    "Error: Invalid position range in region string: %s", region_str);
        }

        // case: only one substring. we need to determine if the single substring is a pos or a chrom. if it
        // is both a valid chrom and pos (i.e. it is a single number e.g. 22), then we create
        // two regions - one for the chrom and one for the pos
        else {
            if (regions_is_valid_chrom (before_colon)) 
                reg->chrom = before_colon;

            bool has_pos = regions_parse_pos (before_colon, NULL, &reg->start_pos, &reg->end_pos);

            // make sure at least one of them is valid
            ASSERT (reg->chrom || has_pos, "Error: Invalid region string: %s", region_str);

            // if both are valid, but the number is <= MAX_NUM_THAT_WE_ASSUME_IS_A_CHROM_AND_NOT_POS, we assume it is a chromosome.
            // Otherwise, we create two regions. Note: the user can always force a region with 10-10 or 1:10
#define MAX_NUM_THAT_WE_ASSUME_IS_A_CHROM_AND_NOT_POS 50
            if (reg->chrom && has_pos) {

                if (reg->start_pos > MAX_NUM_THAT_WE_ASSUME_IS_A_CHROM_AND_NOT_POS) {
                    reg2 = &((Region *)regions_buf.data)[regions_buf.len++];
                    reg2->chrom = reg->chrom;
                    reg ->chrom = NULL;
                    reg2->chrom_node_index = NIL;
                    reg2->start_pos = 0;
                    reg2->end_pos = 0xffffffff;
                }
                else {
                    reg->start_pos = 0;
                    reg->end_pos = 0xffffffff;
                }
            }
        }

        printf ("region: chrom=%s start=%u end=%u\n", reg->chrom ? reg->chrom : "NONE" , reg->start_pos, reg->end_pos);
        if (reg2) printf ("region: chrom=%s start=%u end=%u\n", reg2->chrom ? reg2->chrom : "NONE" , reg2->start_pos, reg2->end_pos);
    }
}

void regions_get_chrom_index()
{
    if (!flag_regions) return; // nothing to do

    Region *regions = (Region *)regions_buf.data;
    MtfContext *chrom_ctx = &evb->z_file->mtf_ctx[CHROM];

    for (int i=0; i < regions_buf.len; i++) {

        if (regions[i].chrom) { // this region is for a specific chrom
            regions[i].chrom_node_index =  mtf_search_for_node_index (chrom_ctx, regions[i].chrom, strlen (regions[i].chrom));

            // if the requested chrom does not exist in the VCF, we remove this region
            if (regions[i].chrom_node_index == NIL) {
                memcpy (&regions[i], &regions[i+1], (regions_buf.len-i-1) * sizeof (Region));
                regions_buf.len--;
                i--;
            }
        }
        else // this region is for all chroms
            regions[i].chrom_node_index = NIL;
    }
}

// we check if if a range within a vb (represented by the parameters of this funciton) intersects with any of the -r/-R regions
bool regions_is_region_included_in_requested_regions (uint32_t chrom_node_index, uint32_t min_pos, uint32_t max_pos)
{
    if (!flag_regions) return true; // nothing to do

    Region *regions = (Region *)regions_buf.data;

    // if any -r region intersects with this VB random_access region - this this VB should be considered
    for (unsigned i=0; i < regions_buf.len; i++)
        if ((regions[i].chrom_node_index == chrom_node_index || regions[i].chrom_node_index == NIL) && // chrom matches, or -r didn't specify which chrom=we match against all
            !(regions[i].start_pos > max_pos || regions[i].end_pos < min_pos)) // regions are intersecting
            return true;

    return false;
}
