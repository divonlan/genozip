// ------------------------------------------------------------------
//   regions.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "regions_vcf.h"
#include "buffer.h"
#include "move_to_front.h"
#include "header.h"
#include "vblock.h"
#include "file.h"

// region as parsed from the --regions option
typedef struct {
    const char *chrom;        // NULL means all chromosomes (i.e. not a specific chromosome)
    uint32_t start_pos;       // if the user did specify pos then start_pos=0 and end_pos=0xffffffff
    uint32_t end_pos;         // the region searched will include both the start and the end
} Region;

// region of a specific chromosome
typedef struct {
    uint32_t start_pos;       // if the user did specify pos then start_pos=0 and end_pos=0xffffffff
    uint32_t end_pos;         // the region searched will include both the start and the end
} Chreg; // = Chromosome Region

static Buffer regions_buf = EMPTY_BUFFER; // all regions together
static Buffer *chregs = NULL;     // one entry per chrom

static uint32_t num_chroms;

static bool is_negative_regions = false; // true if the user used ^ to negate the regions

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

// called from main when parsing the command line to add the argument of --regions
void regions_add (const char *region_str)
{
    ASSERT0 (region_str, "Error: region_str is NULL");

    bool is_negated = region_str[0] == '^';

    bool is_conflicting_negation = (regions_buf.len && (is_negative_regions != is_negated));
    ASSERT0 (!is_conflicting_negation, "Error: inconsistent negation - all regions listed must either be negated or not");

    is_negative_regions = is_negated;

    // make a copy of the string and leave the original one for error message. 
    // we don't free this memory as chrom fields in regions will be pointing to it
    char *rs = malloc (strlen (region_str)+1); 
    strcpy (rs, region_str + is_negated); // drop the ^ if there is one

    char *next_region_token = rs;

    while (1) {
        char *one_rs = strtok_r (next_region_token, ",", &next_region_token);
        if (!one_rs) break;

        buf_alloc (evb, &regions_buf, MAX (regions_buf.len + 1, 100), 2, "regions_buf", 0);

        char *after_colon;
        char *before_colon = strtok_r (one_rs, ":", &after_colon);

        ASSERT (before_colon, "Error: invalid region string: %s", region_str);

        Region *reg = &((Region *)regions_buf.data)[regions_buf.len++]; // update after possible realloc
        Region *reg2 = NULL;

        reg->start_pos = 0;
        reg->end_pos = 0xffffffff;
        reg->chrom = NULL;

        // case: we have both chrom and pos - easy!
        if (after_colon && after_colon[0]) {
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
                    reg2->start_pos = 0;
                    reg2->end_pos = 0xffffffff;
                }
                else {
                    reg->start_pos = 0;
                    reg->end_pos = 0xffffffff;
                }
            }
        }

        //regions_display ("After regions_add"); 
    }
}

// convert the list of regions as parsed from --regions, to an array of chregs - one for each chromosome.
// 1. convert the chrom string to a chrom word index
// 2. for "all chrom" regions - include them in all chregs
void regions_make_chregs(void)
{
    if (!flag_regions) return; // nothing to do

    ARRAY (Region, regions, regions_buf);
    MtfContext *chrom_ctx = &z_file->mtf_ctx[VCF_CHROM];

    num_chroms = chrom_ctx->word_list.len;
    chregs = calloc (num_chroms, sizeof (Buffer)); // a module global variable - array of buffers, one for each chrom
    
    for (int i=0; i < regions_buf.len; i++) {

        Region *reg = &regions[i];
        
        int32_t chrom_word_index = NIL; // All chromosomes, unless reg->chrom is defined
        if (reg->chrom) {
            chrom_word_index = mtf_search_for_node_index (chrom_ctx, regions[i].chrom, strlen (regions[i].chrom));

            // if the requested chrom does not exist in the VCF, we remove this region
            if (chrom_word_index == NIL) continue;
        }

        // loop through either a single chromosome or all chromosomes 
        for (unsigned chr_i = (chrom_word_index == NIL ? 0 : chrom_word_index);
             chr_i         <= (chrom_word_index == NIL ? num_chroms-1 : chrom_word_index);
             chr_i++) {
            
            buf_alloc (evb, &chregs[chr_i], (++chregs[chr_i].len) * sizeof (Chreg), 2, "chregs", chr_i);
            
            Chreg *chreg = LASTENT (Chreg, chregs[chr_i]);
            chreg->start_pos = reg->start_pos;
            chreg->end_pos   = reg->end_pos;
        }
    }

    buf_destroy (&regions_buf); // free the memory, we won't need it again

    //regions_display("After regions_make_chregs");
}

// a user can specific a negative region eg ^13:100-200. We convert the set of negative regions
// to the complementary set of positive regions
void regions_transform_negative_to_positive_complement()
{
    if (!is_negative_regions) return; // nothing to do

    Buffer *neg_chregs = chregs;
    chregs = calloc (num_chroms, sizeof (Buffer));

    // initialize regions for each chr - to be the whole chr
    for (unsigned chr_i=0; chr_i < num_chroms; chr_i++) {
        buf_alloc (evb, &chregs[chr_i], sizeof (Chreg), 1, "chregs", chr_i);
        Chreg *chreg = ENT (Chreg, chregs[chr_i], 0);
        chreg->start_pos   = 0;
        chreg->end_pos     = 0xffffffff;
        chregs[chr_i].len  = 1;

        // process each negative regions - substract from positive chreg 
        for (unsigned negreg_i=0; negreg_i < neg_chregs[chr_i].len; negreg_i++) {

            Chreg *neg_chreg = ENT (Chreg, neg_chregs[chr_i], negreg_i);
            
            // for each positive region that overlaps the negative region - fix it to remove the negative region
            for (unsigned posreg_i=0; posreg_i < chregs[chr_i].len; posreg_i++) {

                Chreg *pos_chreg = ENT (Chreg, chregs[chr_i], posreg_i);

                // case: nagative completely eliminates the positive region
                if (neg_chreg->start_pos <= pos_chreg->start_pos && neg_chreg->end_pos >= pos_chreg->end_pos) {
                    memcpy (pos_chreg, ENT (Chreg, chregs[chr_i], posreg_i+1), (chregs[chr_i].len - posreg_i -1) * sizeof (Chreg));
                    chregs[chr_i].len--;
                }

                // case: negative cuts out the lower side of positive
                else if (neg_chreg->start_pos <= pos_chreg->start_pos && neg_chreg->end_pos >= pos_chreg->start_pos)
                    pos_chreg->start_pos = neg_chreg->end_pos + 1;

                // case: negative cuts out the higher side of positive
                else if (neg_chreg->start_pos <= pos_chreg->end_pos && neg_chreg->end_pos >= pos_chreg->end_pos)
                    pos_chreg->end_pos = neg_chreg->start_pos-1;

                // case: negative is strictly within positive - split positive to the two flanking regions
                else if (neg_chreg->start_pos > pos_chreg->start_pos && neg_chreg->end_pos < pos_chreg->end_pos) {
                    chregs[chr_i].len++;
                    buf_alloc (evb, &chregs[chr_i], chregs[chr_i].len * sizeof (Chreg), 2, "chregs", chr_i);

                    Chreg *new_pos_chreg = LASTENT (Chreg, chregs[chr_i]);
                    Chreg *pos_chreg = ENT (Chreg, chregs[chr_i], posreg_i); // update after realloc

                    new_pos_chreg->start_pos = neg_chreg->end_pos + 1;
                    new_pos_chreg->end_pos   = pos_chreg->end_pos;
                    pos_chreg->end_pos       = neg_chreg->start_pos - 1;
                }
            }
        }
    }

    // copy the destroy the negative buffers
    for (unsigned chr_i=0; chr_i < num_chroms; chr_i++) 
        buf_destroy (&neg_chregs[chr_i]);
    
    is_negative_regions = false; // yay!

    FREE (neg_chregs);

    //regions_display("After regions_transform_negative_to_positive_complement"); 
}

// PIZ: we calculate which regions (specified in the command line -r/-R) intersect with 
// the ra (=a range of a single chrome within a vb) (represented by the parameters of this funciton) - 
// filling in a bytemap of the intersection, and returning true if there is any intersection at all
bool regions_get_ra_intersection (uint32_t chrom_word_index, uint32_t min_pos, uint32_t max_pos,
                                  char *intersection_array) // optional out
{
    if (!flag_regions) return true; // nothing to do

    Buffer *chregs_buf = &chregs[chrom_word_index];

    // if any -r region intersects with this VB random_access region - this this VB should be considered
    bool intersection_found = false;
    for (unsigned chreg_i=0; chreg_i < chregs_buf->len; chreg_i++) {

        Chreg *chreg = ENT (Chreg, *chregs_buf, chreg_i);

        if (chreg->start_pos <= max_pos && chreg->end_pos >= min_pos) { // regions are intersecting
            if (intersection_array) intersection_array[chreg_i] = true;
            intersection_found = true;
        }
    }
    return intersection_found;
}

// PIZ: check if a (chrom,pos) that comes from a specific line, is included in any positive region of
// a specific ra (i.e. chromosome)
bool regions_is_site_included (uint32_t chrom_word_index, uint32_t pos)
{
    // it sufficient that the site is included in one (positive) region
    Buffer *chregs_buf = &chregs[chrom_word_index];
    for (unsigned chreg_i=0; chreg_i < chregs_buf->len; chreg_i++) {
        Chreg *chreg = ENT (Chreg, *chregs_buf, chreg_i);
        if (pos >= chreg->start_pos && pos <= chreg->end_pos) return true;
    }
    return false;
}

unsigned regions_max_num_chregs(void) 
{ 
    static int result = 0; // initialize

    if (!result)
        for (unsigned chr_i=0; chr_i < num_chroms; chr_i++)
            if (chregs[chr_i].len > result) result = chregs[chr_i].len;

    return result; 
}

void regions_display(const char *title)
{
    fprintf (stderr, "regions_display: %s\n", title);

    if (buf_is_allocated (&regions_buf)) {

        fprintf (stderr, "Showing %u %s regions:\n", (uint32_t)regions_buf.len, is_negative_regions ? "NEGATIVE" : "POSITIVE");

        for (unsigned reg_i = 0; reg_i < regions_buf.len; reg_i++) {
            Region *reg = &((Region *)regions_buf.data)[reg_i];
            fprintf (stderr, "chrom=%s start=%u end=%u\n", 
                    reg->chrom ? reg->chrom : "ALL", reg->start_pos, reg->end_pos); 
        }
    }

    if (chregs) {
        fprintf (stderr, "Showing %s chromosomeXregions (\"chregs\") across %u chromosomes:\n", is_negative_regions ? "NEGATIVE" : "POSITIVE", num_chroms);

        for (unsigned chr_i = 0; chr_i < num_chroms; chr_i++)
            for (unsigned chreg_i=0; chreg_i < chregs[chr_i].len; chreg_i++) {
                Chreg *chreg = ENT (Chreg, chregs[chr_i], chreg_i);
                fprintf (stderr, "chrom_word_index=%d start=%u end=%u\n", chr_i, chreg->start_pos, chreg->end_pos); 
            }
    }
}