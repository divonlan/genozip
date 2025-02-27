// ------------------------------------------------------------------
//   regions.c
//   Copyright (C) 2020-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "regions.h"
#include "file.h"
#include "contigs.h"
#include "piz.h"

// region as parsed from the --regions option
typedef struct {
    rom chrom;           // Pointer into regions_data. NULL means all chromosomes (i.e. not a specific chromosome)
    PosType64 start_pos; // if the user did specify pos then start_pos=0 and end_pos=MAX_POS
    PosType64 end_pos;   // the region searched will include both the start and the end
} Region, *RegionP;

// region of a specific chromosome
typedef struct {
    PosType64 start_pos; // if the user did specify pos then start_pos=0 and end_pos=MAX_POS
    PosType64 end_pos;   // the region searched will include both the start and the end
    bool revcomp;        // display the region in reverse complement
} Chreg, *ChregP; // = Chromosome Region

static Buffer regions_data = {};// if using --regions
static BufferP regions_data_from_file = NULL; // if using --regions-file
static Buffer regions_buf = {}; // all regions together
static BufferP chregs = NULL;   // an array of Buffer - one entry per chrom

static WordIndex num_chroms;    // signed value as its compared to chrom

static bool is_negative_regions = false; // true if the user used ^ to negate the regions

// returns true if this is valid pos range string 
static bool regions_parse_pos (rom str, Region *reg) 
{
    unsigned len = strlen (str);

    reg->start_pos = 0;
    reg->end_pos   = MAX_POS;

    // case: "-1000"
    if (str[0] == '-') return str_get_int (&str[1], len-1, &reg->end_pos);

    // case: "1000-"
    if (str[len-1] == '-') return str_get_int (str, len-1, &reg->start_pos);

    // case: "1000-1500" (start 1000, end 1500)
    rom sep = strchr (str, '-');
    if (sep) {
        if (!str_get_int (str, sep - str, &reg->start_pos)) return false;
        if (!str_get_int (sep+1, str+len-(sep+1), &reg->end_pos)) return false;

        // case: "1000000-5" (start 1000000, end 999996)
        if (reg->start_pos > reg->end_pos && reg->end_pos < 1000) 
            reg->end_pos = MAX_(1, reg->start_pos - reg->end_pos + 1);

        return true;
    }

    // case "1000+500" (start 1000, length 500)
    sep = strchr (str, '+');
    if (sep) {
        if (!str_get_int (str, sep - str, &reg->start_pos)) return false;
        PosType64 region_len;
        if (!str_get_int (sep+1, str+len-(sep+1), &region_len)) return false;
        reg->end_pos = reg->start_pos + region_len - 1;
        return true;
    }

    // case: "1000"
    if (!str_get_int (str, len, &reg->start_pos)) return false;
    reg->end_pos = reg->start_pos;
    return true;
}

static bool regions_is_valid_chrom (rom str)
{
    // if it looks like a non-singleton range, we take it as not being a pos
    Region reg;
    if (regions_parse_pos (str, &reg) && reg.start_pos != reg.end_pos) return false;

    // per VCF 4.2 specification, the limitations on CHROM are it is not allowed to contain whitespace or a colon
    unsigned len = strlen (str);

    for (unsigned i=0; i < len; i++)
        if (str[i] < 33 || str[i] > 126 || str[i] == ':') return false;

    return true;
}

// called from main when parsing the command line to add the argument of --regions
void regions_add (rom regions_arg) 
{
    ASSERTNOTNULL (regions_arg);

    buf_append (evb, regions_data, char, regions_arg, strlen (regions_arg)+1, "regions_data");

    bool is_negated = (*B1STc(regions_data) == '^');

    bool is_conflicting_negation = (regions_buf.len && (is_negative_regions != is_negated));
    ASSINP0 (!is_conflicting_negation, "Error: inconsistent negation - all regions listed must either be negated or not");

    is_negative_regions = is_negated;

    char *next_region_token = Bc(regions_data, is_negated); // skip the ^ if there is one

    while (1) {
        char *one_rs = strtok_r (next_region_token, ",", &next_region_token);
        if (!one_rs) break;

        buf_alloc (evb, &regions_buf, 1, 100, Region, 2, "regions_buf");

        char *after_colon;
        char *before_colon = strtok_r (one_rs, ":", &after_colon);

        ASSINP (before_colon, "Error: invalid region string: %s", regions_arg);

        Region *reg = &BNXT (Region, regions_buf);
        *reg = (Region){ .chrom = NULL, .start_pos = 0, .end_pos = MAX_POS };

        // case: we have both chrom and pos - easy!
        if (after_colon && after_colon[0]) {
            ASSINP (regions_is_valid_chrom (before_colon), "Error: Invalid CHROM in region string: %s", regions_arg);
            reg->chrom = before_colon;
            
            ASSINP (regions_parse_pos (after_colon, reg), "Error: Invalid position range in region string: \"%s\"", regions_arg);
        }

        // case: only one substring. we need to determine if the single substring is a pos or a chrom. if it
        // is both a valid chrom and pos (i.e. it is a single number e.g. 22), then we create
        // two regions - one for the chrom and one for the pos
        else {
            if (regions_is_valid_chrom (before_colon)) 
                reg->chrom = before_colon;

            bool has_pos = regions_parse_pos (before_colon, reg);

            // make sure at least one of them is valid
            ASSINP (reg->chrom || has_pos, "Error: Invalid region string: %s", regions_arg);

            // if both are valid, but the number is <= MAX_NUM_THAT_WE_ASSUME_IS_A_CHROM_AND_NOT_POS, we assume it is a chromosome.
            // Otherwise, we create two regions. Note: the user can always force a region with 10-10 or 1:10
            #define MAX_NUM_THAT_WE_ASSUME_IS_A_CHROM_AND_NOT_POS 50
            if (reg->chrom && has_pos) {

                // if large number - have two regions: entire chrom of this number, and all chroms at this pos
                if (reg->start_pos > MAX_NUM_THAT_WE_ASSUME_IS_A_CHROM_AND_NOT_POS) {
                    BNXT (Region, regions_buf) = (Region){ .chrom = reg->chrom, .start_pos = 0, .end_pos = MAX_POS };
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

// called from main when parsing the command line to add the argument of --regions
// file format: tab separated file, each line contains 2 or 3 columns: CHROM POS or CHROM POS END (1-based POS like VCF, inclusive range)
void regions_add_by_file (rom regions_filename)
{
    // common user error -R1 instead of --R1
    ASSINP ((regions_filename[0]!='1' && regions_filename[0]!='2') || regions_filename[1],
            "-R is the short form of --regions-file, and it expects a filename. Did you mean --R%c (double hyphen)?", regions_filename[0]);
    
    if (regions_filename[0] == '^') {
        is_negative_regions = true;
        regions_filename++;
        ASSINP0 (regions_filename[0], "bad --regions-file argument");
    }

    file_split_lines (regions_filename, "regions", VERIFY_ASCII);
    if (!n_lines) return; 

    buf_alloc (evb, &regions_buf, 0, n_lines, Region, 0, "regions_buf");

    for (unsigned i=0; i < n_lines; i++) { // -1 as last line (following terminating newline) is empty

        if (lines[i][0] == '#') continue; // skip comment line
        
        unsigned num_fields = str_count_char (lines[i], line_lens[i], '\t') + 1;
        ASSINP (num_fields == 2 || num_fields == 3, "Expecting either 2 or 3 tab-separated columns in line %u of %s, but found: \"%.*s\"",
                i+1, regions_filename, line_lens[i], lines[i]);

        str_split_enforce (lines[i], line_lens[i], num_fields, '\t', field, true, "fields");

        ((char**)fields)[0][field_lens[0]] = 0; // nul-terminate chrom
        ASSINP (regions_is_valid_chrom (fields[0]), "Invalid chrom (column 1 in a tab-separated file) in %s line %i: \"%s\"", regions_filename, i+1, fields[0]);
        
        bool has_len = num_fields == 3 && fields[2][0] == '+';
        
        PosType64 start_pos, end_pos, len;
        ASSINP (str_get_int (fields[1], field_lens[1], &start_pos), "Invalid start pos (column 2 in a tab-separated file) in %s line %u: \"%.*s\"", regions_filename, i+1, field_lens[1], fields[1]);
        
        ASSINP (has_len || num_fields == 2 || str_get_int_range64 (fields[2], field_lens[2], start_pos, MAX_POS, &end_pos), 
                "Invalid end pos (column 3 in a tab-separated file) in %s line %u: \"%.*s\"", regions_filename, i+1, field_lens[2], fields[2]);

        ASSINP (!has_len || str_get_int_range64 (fields[2]+1, field_lens[2]-1, 0, MAX_POS, &len), 
                "Invalid len (column 3 in a tab-separated file) in %s line %u: \"%.*s\"", regions_filename, i+1, field_lens[2], fields[2]);

        BNXT (Region, regions_buf) = (Region){ 
            .chrom     = fields[0], // points into "data" allocated by file_split_lines
            .start_pos = start_pos, 
            .end_pos   = (num_fields==2) ? start_pos 
                       : has_len         ? (start_pos + len - 1) // user specified length, eg +10, rather than end
                       :                   end_pos 
        };
    }

    regions_data_from_file = &data; // so we can destroy it if needed
}

// convert the list of regions as parsed from --regions, to an array of chregs - one for each chromosome.
// 1. convert the chrom string to a chrom word index
// 2. for "all chrom" regions - include them in all chregs
void regions_make_chregs (ContextP chrom_ctx)
{
    ARRAY (Region, regions, regions_buf);

    num_chroms = chrom_ctx->word_list.len;
    // case: FASTA compressed without contigs, or FASTQ...
    ASSINP (num_chroms, "--regions is not supported for this file because it was not indexed during compression%s",
            (Z_DT(FASTA) || FAF) ? " (to index, compress with --index)" : ""); 

    chregs = CALLOC (num_chroms * sizeof (Buffer)); // a module global variable - array of buffers, one for each chrom
    
    for (int i=0; i < regions_len; i++) {

        Region *reg = &regions[i];
        
        int32_t chrom_word_index = WORD_INDEX_NONE; // All chromosomes, unless reg->chrom is defined
        if (reg->chrom) {
            chrom_word_index = ctx_search_for_word_index (chrom_ctx, regions[i].chrom, strlen (regions[i].chrom));

            // if we have a reference file loaded, try getting an alternative name
            if (chrom_word_index == WORD_INDEX_NONE)
                chrom_word_index = contigs_get_matching (ref_get_ctgs(), regions[i].chrom, strlen (regions[i].chrom), 0, true, NULL);

            // if the requested chrom does not exist in the file, we remove this region
            if (chrom_word_index == WORD_INDEX_NONE) continue;
        }

        // loop through either a single chromosome or all chromosomes 
        for (unsigned chr_i = (chrom_word_index == WORD_INDEX_NONE ? 0 : chrom_word_index);
             chr_i         <= (chrom_word_index == WORD_INDEX_NONE ? num_chroms-1 : chrom_word_index);
             chr_i++) {

            RangeP r = ref_get_range_by_ref_index (evb, chr_i);

            PosType64 reg_start_pos = reg->start_pos - ((r && flag.gpos) ? (r->gpos - 1) : 0);
            PosType64 reg_end_pos   = reg->end_pos   - ((r && flag.gpos) ? (r->gpos - 1) : 0);
            
            chregs[chr_i].len++;
            buf_alloc (evb, &chregs[chr_i], 0, chregs[chr_i].len, Chreg, 2, "chregs");
            
            ChregP chreg = BLST (Chreg, chregs[chr_i]);
            chreg->revcomp   = reg_end_pos < reg_start_pos; 
            chreg->start_pos = chreg->revcomp ? reg_end_pos : reg_start_pos;
            chreg->end_pos   = chreg->revcomp ? reg_start_pos : reg_end_pos;

            // remove chreg if it does not overlap chromosome at all
            if (r && (chreg->start_pos > r->last_pos || chreg->end_pos < 1)) 
                chregs[chr_i].len--;
        }
    }

    //regions_display("After regions_make_chregs");
}

// a user can specify negative region eg ^13:100-200. We convert the set of negative regions
// to the complementary set of positive regions
void regions_transform_negative_to_positive_complement()
{
    if (!is_negative_regions) return; // nothing to do

    BufferP neg_chregs = chregs;
    chregs = CALLOC (num_chroms * sizeof (Buffer));

    // initialize regions for each chr - to be the whole chr
    for (unsigned chr_i=0; chr_i < num_chroms; chr_i++) {
        buf_alloc (evb, &chregs[chr_i], 0, 1, Chreg, 1, "chregs");
        ChregP chreg = B(Chreg, chregs[chr_i], 0);
        chreg->start_pos   = 0;
        chreg->end_pos     = MAX_POS;
        chregs[chr_i].len  = 1;

        // process each negative regions - substract from positive chreg 
        for (unsigned negreg_i=0; negreg_i < neg_chregs[chr_i].len32; negreg_i++) {

            ChregP neg_chreg = B(Chreg, neg_chregs[chr_i], negreg_i);
            
            // for each positive region that overlaps the negative region - fix it to remove the negative region
            for (unsigned posreg_i=0; posreg_i < chregs[chr_i].len32; posreg_i++) {

                ChregP pos_chreg = B(Chreg, chregs[chr_i], posreg_i);

                // case: nagative completely eliminates the positive region
                if (neg_chreg->start_pos <= pos_chreg->start_pos && neg_chreg->end_pos >= pos_chreg->end_pos) {
                    memcpy (pos_chreg, B(Chreg, chregs[chr_i], posreg_i+1), (chregs[chr_i].len - posreg_i -1) * sizeof (Chreg));
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
                    
                    ChregP new_pos_chreg = BLST (Chreg, chregs[chr_i]);
                    ChregP pos_chreg = B(Chreg, chregs[chr_i], posreg_i); // update after realloc

                    new_pos_chreg->start_pos = neg_chreg->end_pos + 1;
                    new_pos_chreg->end_pos   = pos_chreg->end_pos;
                    pos_chreg->end_pos       = neg_chreg->start_pos - 1;
                }
            }
        }
    }

    // copy the destroy the negative buffers
    for (unsigned chr_i=0; chr_i < num_chroms; chr_i++) 
        buf_destroy (neg_chregs[chr_i]);
    
    is_negative_regions = false; // yay!

    FREE (neg_chregs);

    //regions_display("After regions_transform_negative_to_positive_complement"); 
}

// PIZ: we calculate which regions (specified in the command line -r/-R) intersect with 
// the ra (=a range of a single chrome within a vb) (represented by the parameters of this funciton) - 
// filling in a bytemap of the intersection, and returning true if there is any intersection at all
bool regions_get_ra_intersection (WordIndex chrom_word_index, PosType64 min_pos, PosType64 max_pos)
{
    if (!flag.regions) return true; // nothing to do

    ASSERTINRANGE (chrom_word_index, 0, num_chroms);

    // if any -r region intersects with this VB random_access region - this this VB should be considered
    ARRAY (Chreg, ch, chregs[chrom_word_index]);
    for (uint64_t chreg_i=0; chreg_i < ch_len; chreg_i++) 
        if (ch[chreg_i].start_pos <= max_pos && ch[chreg_i].end_pos >= min_pos)  // regions are intersecting
            return true; //intersection found
    
    return false; // no intersection found
}

unsigned regions_get_num_range_intersections (WordIndex chrom_word_index)
{
    return flag.regions ? chregs[chrom_word_index].len : 1; // no --regions = whole chromosome
}

// used by ref_display_ref. if an intersection was found - returns the min,max pos and true, otherwise returns false
bool regions_get_range_intersection (WordIndex chrom_word_index, PosType64 min_pos, PosType64 max_pos, unsigned intersect_i,
                                     PosType64 *intersect_min_pos, PosType64 *intersect_max_pos, bool *revcomp) // out
{
    if (!flag.regions) { // if no regions are specified, the entire range "intersects"
        *intersect_min_pos = min_pos;
        *intersect_max_pos = max_pos;
        return true;
    }

    BufferP chregs_buf = &chregs[chrom_word_index];
    if (intersect_i >= chregs_buf->len) return false; // intersect_i is out of range with this chromosome (possibly because len=0 - no intersections with this chromosome)

    ChregP chreg = B(Chreg, *chregs_buf, intersect_i);

    if (chreg->start_pos > max_pos || chreg->end_pos < min_pos) 
        return false; // no intersection with this chromosome

    *intersect_min_pos = MAX_(min_pos, chreg->start_pos);
    *intersect_max_pos = MIN_(max_pos, chreg->end_pos);
    *revcomp = chreg->revcomp;

    return true;
}

// PIZ: check if a (chrom,pos) that comes from a specific line, is included in any positive region of
// a specific ra (i.e. chromosome)
bool regions_is_site_included (VBlockP vb)
{
    Did pos_did_i   = DTF(pos);

    WordIndex chrom = vb->last_index (CHROM);
    PosType64 pos = (pos_did_i == DID_NONE)          ? 1 
                : CTX(pos_did_i)->pos_last_value ? CTX(pos_did_i)->pos_last_value // use saved value if one exists (used in VCF, bc VCF_POS.last_value might be modified by a INFO/END)
                :                                  CTX(pos_did_i)->last_value.i;
    
    ASSPIZ (IN_RANGE(chrom, 0, num_chroms), "chrom=%d ∉ [0,%u)", chrom, num_chroms);

    // it sufficient that the site is included in one (positive) region
    BufferP chregs_buf = &chregs[chrom];
    for (unsigned chreg_i=0; chreg_i < chregs_buf->len; chreg_i++) {
        ChregP chreg = B(Chreg, *chregs_buf, chreg_i);
        if (pos >= chreg->start_pos && pos <= chreg->end_pos) return true;
    }
    return false;
}

// PIZ: check if a range (chrom,start_pos,end_pos) overlaps with an included region. used when loading reference ranges.
bool regions_is_range_included (WordIndex chrom_word_index, PosType64 start_pos, PosType64 end_pos, bool completely_included)
{
    ASSERTINRANGE (chrom_word_index, 0, num_chroms);

    // it sufficient that the site is included in one (positive) region
    BufferP chregs_buf = &chregs[chrom_word_index];
    for (unsigned chreg_i=0; chreg_i < chregs_buf->len; chreg_i++) {
        ChregP chreg = B(Chreg, *chregs_buf, chreg_i);

        // check for complete inclusion
        if (start_pos >= chreg->start_pos && end_pos <= chreg->end_pos) 
            return true;
        
        // check for region overlap
        if (end_pos >= chreg->start_pos && start_pos <= chreg->end_pos) { // requested range overlaps chreg

            // if we only need overlap, we're done
            if (!completely_included) 
                return true;

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

void regions_display(rom title)
{
    iprintf ("regions_display: %s\n", title);

    if (buf_is_alloc (&regions_buf)) {

        iprintf ("Showing %u %s regions:\n", regions_buf.len32, is_negative_regions ? "NEGATIVE" : "POSITIVE");

        for (unsigned reg_i = 0; reg_i < regions_buf.len32; reg_i++) {
            Region *reg = &((Region *)regions_buf.data)[reg_i];
            iprintf ("chrom=%s start=%"PRId64" end=%"PRId64"\n", 
                     reg->chrom ? reg->chrom : "ALL", reg->start_pos, reg->end_pos); 
        }
    }

    if (chregs) {
        iprintf ("Showing %s chromosomeXregions (\"chregs\") across %u chromosomes:\n", is_negative_regions ? "NEGATIVE" : "POSITIVE", num_chroms);

        for (unsigned chr_i = 0; chr_i < num_chroms; chr_i++)
            for (unsigned chreg_i=0; chreg_i < chregs[chr_i].len32; chreg_i++) {
                ChregP chreg = B(Chreg, chregs[chr_i], chreg_i);
                iprintf ("chrom_word_index=%d start=%"PRId64" end=%"PRId64"\n", chr_i, chreg->start_pos, chreg->end_pos); 
            }
    }
}

void regions_destroy (void)
{
    if (!is_genocat) return; // --regions and --regions-file are available only in genocat

    buf_destroy (regions_data);
    buf_destroy (regions_buf);

    if (chregs) {
        for (unsigned chr_i = 0; chr_i < num_chroms; chr_i++)
            buf_destroy (chregs[chr_i]);

        FREE (chregs);
    }

    if (regions_data_from_file)
        buf_destroy (*regions_data_from_file);
}
