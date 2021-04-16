// ------------------------------------------------------------------
//   regions.c
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "regions.h"
#include "buffer.h"
#include "context.h"
#include "vblock.h"
#include "file.h"
#include "strings.h"

// region as parsed from the --regions option
typedef struct {
    const char *chrom;        // NULL means all chromosomes (i.e. not a specific chromosome)
    PosType start_pos;       // if the user did specify pos then start_pos=0 and end_pos=MAX_POS
    PosType end_pos;         // the region searched will include both the start and the end
} Region;

// region of a specific chromosome
typedef struct {
    PosType start_pos;       // if the user did specify pos then start_pos=0 and end_pos=MAX_POS
    PosType end_pos;         // the region searched will include both the start and the end
} Chreg; // = Chromosome Region

static Buffer regions_buf = EMPTY_BUFFER; // all regions together
static Buffer *chregs = NULL;     // one entry per chrom

static uint32_t num_chroms;

static bool is_negative_regions = false; // true if the user used ^ to negate the regions

// returns true if this is valid pos range string 
static bool regions_parse_pos (const char *str, Region *reg) 
{
    unsigned len = strlen (str);

    reg->start_pos = 0;
    reg->end_pos   = MAX_POS;

    // case: "-1000"
    if (str[0] == '-') return str_get_int (&str[1], len-1, &reg->end_pos);

    // case: "1000-"
    if (str[len-1] == '-') return str_get_int (str, len-1, &reg->start_pos);

    // case: "1000-1500" (start 1000, end 1500)
    const char *sep = strchr (str, '-');
    if (sep) {
        if (!str_get_int (str, sep - str, &reg->start_pos)) return false;
        if (!str_get_int (sep+1, str+len-(sep+1), &reg->end_pos)) return false;

        if (reg->start_pos > reg->end_pos) SWAP (reg->start_pos, reg->end_pos);

        return true;
    }

    // case "1000+500" (start 1000, length 500)
    sep = strchr (str, '+');
    if (sep) {
        if (!str_get_int (str, sep - str, &reg->start_pos)) return false;
        PosType region_len;
        if (!str_get_int (sep+1, str+len-(sep+1), &region_len)) return false;
        reg->end_pos = reg->start_pos + region_len - 1;
        return true;
    }

    // case: "1000"
    if (!str_get_int (str, len, &reg->start_pos)) return false;
    reg->end_pos = reg->start_pos;
    return true;
}

static bool regions_is_valid_chrom (const char *str)
{
    // if it looks like a non-singleton range, we take it as not being a pos
    Region reg;
    if (regions_parse_pos (str, &reg) && reg.start_pos != reg.end_pos) return false;

    // per VCF 4.2 specification, the limitations on CHROM are it is not allowed to contain whitespace or a colon
    unsigned len = strlen (str);

    for (unsigned i=0; i < len; i++)
        if (str[i] < 33 || str[i] > 127 || str[i] == ':') return false;

    return true;
}

// called from main when parsing the command line to add the argument of --regions
void regions_add (const char *region_str)
{
    ASSERTNOTNULL (region_str);

    bool is_negated = region_str[0] == '^';

    bool is_conflicting_negation = (regions_buf.len && (is_negative_regions != is_negated));
    ASSINP0 (!is_conflicting_negation, "Error: inconsistent negation - all regions listed must either be negated or not");

    is_negative_regions = is_negated;

    // make a copy of the string and leave the original one for error message. 
    // we don't free this memory as chrom fields in regions will be pointing to it
    char *rs = MALLOC (strlen (region_str)+1); // heap memoery, as reg->chrom points into this
    strcpy (rs, region_str + is_negated); // drop the ^ if there is one

    char *next_region_token = rs;

    while (1) {
        char *one_rs = strtok_r (next_region_token, ",", &next_region_token);
        if (!one_rs) break;

        buf_alloc (evb, &regions_buf, 1, 100, Region, 2, "regions_buf");

        char *after_colon;
        char *before_colon = strtok_r (one_rs, ":", &after_colon);

        ASSINP (before_colon, "Error: invalid region string: %s", region_str);

        Region *reg = &NEXTENT (Region, regions_buf);
        *reg = (Region){ .chrom = NULL, .start_pos = 0, .end_pos = MAX_POS };

        // case: we have both chrom and pos - easy!
        if (after_colon && after_colon[0]) {
            ASSINP (regions_is_valid_chrom (before_colon), "Error: Invalid CHROM in region string: %s", region_str);
            reg->chrom = before_colon;
            
            ASSINP (regions_parse_pos (after_colon, reg), "Error: Invalid position range in region string: %s", region_str);
        }

        // case: only one substring. we need to determine if the single substring is a pos or a chrom. if it
        // is both a valid chrom and pos (i.e. it is a single number e.g. 22), then we create
        // two regions - one for the chrom and one for the pos
        else {
            if (regions_is_valid_chrom (before_colon)) 
                reg->chrom = before_colon;

            bool has_pos = regions_parse_pos (before_colon, reg);

            // make sure at least one of them is valid
            ASSINP (reg->chrom || has_pos, "Error: Invalid region string: %s", region_str);

            // if both are valid, but the number is <= MAX_NUM_THAT_WE_ASSUME_IS_A_CHROM_AND_NOT_POS, we assume it is a chromosome.
            // Otherwise, we create two regions. Note: the user can always force a region with 10-10 or 1:10
#define MAX_NUM_THAT_WE_ASSUME_IS_A_CHROM_AND_NOT_POS 50
            if (reg->chrom && has_pos) {

                // if large number - have two regions: entire chrom of this number, and all chroms at this pos
                if (reg->start_pos > MAX_NUM_THAT_WE_ASSUME_IS_A_CHROM_AND_NOT_POS) {
                    NEXTENT (Region, regions_buf) = (Region){ .chrom = reg->chrom, .start_pos = 0, .end_pos = MAX_POS };
                    reg->chrom = NULL;
                }

                // if small number - assume it is a single entire chrom
                else { 
                    reg->start_pos = 0;
                    reg->end_pos = MAX_POS;
                }
            }
        }

        //regions_display ("After regions_add"); 
    }
}

// convert the list of regions as parsed from --regions, to an array of chregs - one for each chromosome.
// 1. convert the chrom string to a chrom word index
// 2. for "all chrom" regions - include them in all chregs
void regions_make_chregs (void)
{
    ARRAY (Region, regions, regions_buf);
    Context *chrom_ctx = &z_file->contexts[CHROM];

    num_chroms = chrom_ctx->word_list.len;
    chregs = CALLOC (num_chroms * sizeof (Buffer)); // a module global variable - array of buffers, one for each chrom
    
    for (int i=0; i < regions_len; i++) {

        Region *reg = &regions[i];
        
        int32_t chrom_word_index = WORD_INDEX_NONE; // All chromosomes, unless reg->chrom is defined
        if (reg->chrom) {
            chrom_word_index = ctx_search_for_word_index (chrom_ctx, regions[i].chrom, strlen (regions[i].chrom));

            // if the requested chrom does not exist in the file, we remove this region
            if (chrom_word_index == WORD_INDEX_NONE) continue;
        }

        // loop through either a single chromosome or all chromosomes 
        for (unsigned chr_i = (chrom_word_index == WORD_INDEX_NONE ? 0 : chrom_word_index);
             chr_i         <= (chrom_word_index == WORD_INDEX_NONE ? num_chroms-1 : chrom_word_index);
             chr_i++) {
            
            chregs[chr_i].len++;
            buf_alloc (evb, &chregs[chr_i], 0, chregs[chr_i].len, Chreg, 2, "chregs");
            
            Chreg *chreg = LASTENT (Chreg, chregs[chr_i]);
            chreg->start_pos = reg->start_pos;
            chreg->end_pos   = reg->end_pos;
        }
    }

    //regions_display("After regions_make_chregs");
}

// a user can specific a negative region eg ^13:100-200. We convert the set of negative regions
// to the complementary set of positive regions
void regions_transform_negative_to_positive_complement()
{
    if (!is_negative_regions) return; // nothing to do

    Buffer *neg_chregs = chregs;
    chregs = CALLOC (num_chroms * sizeof (Buffer));

    // initialize regions for each chr - to be the whole chr
    for (unsigned chr_i=0; chr_i < num_chroms; chr_i++) {
        buf_alloc (evb, &chregs[chr_i], 0, 1, Chreg, 1, "chregs");
        Chreg *chreg = ENT (Chreg, chregs[chr_i], 0);
        chreg->start_pos   = 0;
        chreg->end_pos     = MAX_POS;
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
                    buf_alloc (evb, &chregs[chr_i], 0, chregs[chr_i].len, Chreg, 2, "chregs");
                    
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
bool regions_get_ra_intersection (WordIndex chrom_word_index, PosType min_pos, PosType max_pos)
{
    if (!flag.regions) return true; // nothing to do

    ASSERT (chrom_word_index >= 0 && chrom_word_index < num_chroms, "chrom_word_index=%d out of range", chrom_word_index);

    // if any -r region intersects with this VB random_access region - this this VB should be considered
    ARRAY (Chreg, ch, chregs[chrom_word_index]);
    for (uint64_t chreg_i=0; chreg_i < ch_len; chreg_i++) 
        if (ch[chreg_i].start_pos <= max_pos && ch[chreg_i].end_pos >= min_pos)  // regions are intersecting
            return true; //intersection found
    
    return false; // no intersection found
}

// used by ref_display_ref. if an intersection was found - returns the min,max pos and true, otherwise returns false
bool regions_get_range_intersection (WordIndex chrom_word_index, PosType min_pos, PosType max_pos,
                                     PosType *intersect_min_pos, PosType *intersect_max_pos) // out
{
    if (!flag.regions) { // if no regions are specified, the entire range "intersects"
        *intersect_min_pos = min_pos;
        *intersect_max_pos = max_pos;
        return true;
    }

    Buffer *chregs_buf = &chregs[chrom_word_index];
    if (!chregs_buf->len) return false; // no intersection with this chromosome

    ASSINP0 (chregs_buf->len==1, "Error: when using --regions to display a reference, you can specify at most on region per chromosome");

    Chreg *chreg = FIRSTENT (Chreg, *chregs_buf);

    if (chreg->start_pos > max_pos || chreg->end_pos < min_pos) 
        return false; // no intersection with this chromosome

    *intersect_min_pos = MAX (min_pos, chreg->start_pos);
    *intersect_max_pos = MIN (max_pos, chreg->end_pos);

    return true;
}

// PIZ: check if a (chrom,pos) that comes from a specific line, is included in any positive region of
// a specific ra (i.e. chromosome)
bool regions_is_site_included (WordIndex chrom_word_index, PosType pos)
{
    ASSERT (chrom_word_index >= 0 && chrom_word_index < num_chroms, "chrom_word_index=%d out of range", chrom_word_index);

    // it sufficient that the site is included in one (positive) region
    Buffer *chregs_buf = &chregs[chrom_word_index];
    for (unsigned chreg_i=0; chreg_i < chregs_buf->len; chreg_i++) {
        Chreg *chreg = ENT (Chreg, *chregs_buf, chreg_i);
        if (pos >= chreg->start_pos && pos <= chreg->end_pos) return true;
    }
    return false;
}

// PIZ: check if a range (chrom,start_pos,end_pos) overlaps with an included region. used when loading reference ranges.
bool regions_is_range_included (WordIndex chrom_word_index, PosType start_pos, PosType end_pos, bool completely_included)
{
    ASSERT (chrom_word_index >= 0 && chrom_word_index < num_chroms, "chrom_word_index=%d out of range", chrom_word_index);

    // it sufficient that the site is included in one (positive) region
    Buffer *chregs_buf = &chregs[chrom_word_index];
    for (unsigned chreg_i=0; chreg_i < chregs_buf->len; chreg_i++) {
        Chreg *chreg = ENT (Chreg, *chregs_buf, chreg_i);

        // check for complete inclusion
        if (start_pos >= chreg->start_pos && end_pos <= chreg->end_pos) 
            return true;
        
        // check for region overlap
        if (end_pos >= chreg->start_pos && start_pos <= chreg->end_pos) { // requested range overlaps chreg

            // if we only need overlap, we're done
            if (!completely_included) return true;

            // part of the region is covered, decided based on the remaining part (which can be either before or after this chreg or both)
            bool left_flanking_covered=true, right_flanking_covered=true;

            if (start_pos < chreg->start_pos) 
                left_flanking_covered = regions_is_range_included (chrom_word_index, start_pos, chreg->start_pos-1, true);

            if (end_pos > chreg->end_pos) 
                right_flanking_covered = regions_is_range_included (chrom_word_index, chreg->end_pos+1, end_pos, true);

            return left_flanking_covered && right_flanking_covered;
        }
    }
    return false;
}


unsigned regions_max_num_chregs (void) 
{ 
    static int result = 0; // initialize

    if (!result)
        for (unsigned chr_i=0; chr_i < num_chroms; chr_i++)
            if (chregs[chr_i].len > result) result = chregs[chr_i].len;

    return result; 
}

void regions_display(const char *title)
{
    iprintf ("regions_display: %s\n", title);

    if (buf_is_allocated (&regions_buf)) {

        iprintf ("Showing %u %s regions:\n", (uint32_t)regions_buf.len, is_negative_regions ? "NEGATIVE" : "POSITIVE");

        for (unsigned reg_i = 0; reg_i < regions_buf.len; reg_i++) {
            Region *reg = &((Region *)regions_buf.data)[reg_i];
            iprintf ("chrom=%s start=%"PRId64" end=%"PRId64"\n", 
                     reg->chrom ? reg->chrom : "ALL", reg->start_pos, reg->end_pos); 
        }
    }

    if (chregs) {
        iprintf ("Showing %s chromosomeXregions (\"chregs\") across %u chromosomes:\n", is_negative_regions ? "NEGATIVE" : "POSITIVE", num_chroms);

        for (unsigned chr_i = 0; chr_i < num_chroms; chr_i++)
            for (unsigned chreg_i=0; chreg_i < chregs[chr_i].len; chreg_i++) {
                Chreg *chreg = ENT (Chreg, chregs[chr_i], chreg_i);
                iprintf ("chrom_word_index=%d start=%"PRId64" end=%"PRId64"\n", chr_i, chreg->start_pos, chreg->end_pos); 
            }
    }
}
