// ------------------------------------------------------------------
//   segregate.c
//   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <unistd.h>
#include <inttypes.h>

#include "vczip.h"

static void segregate_pos_field (VariantBlock *vb, const char *str, 
                                 long long *pos_delta, const char **pos_start, unsigned *pos_len,
                                 unsigned line_i)
{
    *pos_start = str;

    // scan by ourselves - hugely faster the sscanf
    long long this_pos=0;
    const char *s; for (s=str; *s != '\t'; s++)
        this_pos = this_pos * 10 + (*s - '0');
    
    *pos_len     = s - str;
    *pos_delta   = this_pos - vb->last_pos;
    vb->last_pos = this_pos;
}

static void segregate_format_field(DataLine *dl, const char *str, unsigned line_i)
{
    if (!global_num_samples) return; // if we're not expecting any samples, we ignore the FORMAT field

    ASSERT(str[0] && str[1] && str[2], "Error: missing samples in line %u, expected by VCF file sample name header line", line_i);

    if (str[0] == 'G' && str[1] == 'T') // GT field in FORMAT columns - must always appear first per VCF spec (if at appears)
        dl->has_haplotype_data = true; 

    if (str[2] == ':' || !dl->has_haplotype_data) {
        dl->has_genotype_data = true; // no GT subfield, or subfields in addition to GT

        // check which subfield (if any) is GL
        dl->gl_subfield = gl_optimize_get_gl_subfield_index(str);
    }
}

static void segregate_variant_area(VariantBlock *vb, DataLine *dl, const char *str, unsigned len,
                                   long long pos_delta, const char *pos_start, unsigned pos_len)
{
    char pos_delta_str[22];
    sprintf (pos_delta_str, "%"PRId64, pos_delta);
    unsigned pos_delta_len = strlen (pos_delta_str);

    vb->line_variant_data.len = len + (pos_delta_len - pos_len) + 1; // +1 for the \n, but not including the \0

    buf_alloc (vb, &vb->line_variant_data, vb->line_variant_data.len + 1, 1, "line_variant_data", dl->line_i); // +1 for \0

    unsigned len_substr_before_pos = pos_start - str;

    sprintf (vb->line_variant_data.data, "%.*s%s%.*s\n",
             len_substr_before_pos, str, // substring prior to POS
             pos_delta_str,
             len - len_substr_before_pos - pos_len, pos_start + pos_len);

    vb->add_bytes[SEC_VARIANT_DATA] -= pos_len - pos_delta_len; // we saved this number of bytes vs. the original file
}

static bool segregate_haplotype_area(VariantBlock *vb, DataLine *dl, const char *str, unsigned len, 
                                     unsigned line_i, unsigned sample_i)
{
    // check ploidy
    unsigned ploidy=1;
    for (unsigned i=1; i<len-1; i++)
        if (str[i] == '|' || str[i] == '/') ploidy++;

    ASSERT (ploidy <= 65535, "Error: ploidy=%u exceeds the maximum of 65535 in line_i=%u", ploidy, line_i);
    
    // if the ploidy of this sample is bigger than the ploidy of the data in this VB so far, then
    // we have to re-do the VB with a larger ploidy. This can happen for example in the X chromosome
    // if initial samples are male with ploidy=1 and then a female sample with ploidy=2
    if (vb->ploidy && ploidy > vb->ploidy) {
        vb->ploidy = ploidy;
        return true; // ploidy overflow
    }

    if (!vb->ploidy) vb->ploidy = ploidy; // very first sample in the vb

    // note - ploidy of this sample might be smaller than vb->ploidy (eg a male sample in an X chromosesome that was preceded by a female sample)

    // initial allocation, or might be enlarged in case of a re-do following a ploidy overflow
    if (!vb->line_haplotype_data.data || vb->line_haplotype_data.size < vb->ploidy * global_num_samples) {
        buf_alloc (vb, &vb->line_haplotype_data, vb->ploidy * global_num_samples, 1, "line_haplotype_data", line_i);
        dl->phase_type = (vb->ploidy==1 ? PHASE_HAPLO : PHASE_UNKNOWN);
    }

    char *ht_data = &vb->line_haplotype_data.data[vb->ploidy * sample_i];

    vb->add_bytes[SEC_PHASE_DATA] -= (ploidy-1); // we "saved" this number of phase characters vs the vcf file (... we will add them back later if we didn't)

    for (unsigned ht_i=0; ht_i < ploidy; ht_i++) {

        char ht = *(str++); 
        len--;

        ASSERT ((ht >= '0' && ht <= '9') || ht == '.',
                "Error: invalid VCF file - line %u - expecting an allele in a sample to be a number 0-9 or . , but seeing %c", line_i, *str);

        // single-digit allele numbers
        ht_data[ht_i] = ht;

        if (!len) break;

        // handle 2-digit allele numbers
        if (ht != '.' && *str >= '0' && *str <= '9') {
            unsigned allele = 10 * (ht-'0') + (*(str++) - '0');
            len--;

            // make sure there isn't a 3rd digit
            ASSERT (!len || *str < '0' || * str > '9', "Error: VCF file - line %u sample %u - vczip currently supports only alleles up to 99", line_i, sample_i+1);

            ht_data[ht_i] = '0' + allele; // use ascii 48->147
            vb->add_bytes[SEC_HAPLOTYPE_DATA]--; // we saved one byte vs the original file
        }

        // get phase
        if (ploidy > 1 && ht_i < ploidy-1) {
            
            PhaseType cell_phase_type = *(str++);
            PhaseType ht0_phase_type;
            len--;

            ASSERT (cell_phase_type == '|' || cell_phase_type == '/', "Error: invalid VCF file - line %u - unable to parse sample %u: expecting a | or / but seeing %c", line_i, sample_i+1, cell_phase_type);

            // deal with phase - only at the first separator eg 1|1|0|1
            if (ht_i==0) { 
                if (cell_phase_type == dl->phase_type) {} // do nothing

                else if ((cell_phase_type == '|' || cell_phase_type == '/') && 
                         (dl->phase_type == PHASE_UNKNOWN || dl->phase_type == PHASE_HAPLO))
                    dl->phase_type = cell_phase_type;

                else if ((dl->phase_type == '|' && cell_phase_type == '/') || 
                         (dl->phase_type == '/' && cell_phase_type == '|')) {
                    dl->phase_type = PHASE_MIXED_PHASED;
                
                    // realloc and not malloc, bc it can be already allocated in case of a re-do following a ploidy overflow
                    buf_alloc (vb, &vb->line_phase_data, global_num_samples, 1, "line_phase_data", line_i);

                    // fill in the so-far uniform phase (which is different than the current one)
                    memset (vb->line_phase_data.data, (cell_phase_type == '|' ? '/' : '|'), sample_i);
                    vb->line_phase_data.data[sample_i] = (char)cell_phase_type;
                }
                else if (dl->phase_type == PHASE_MIXED_PHASED)
                    vb->line_phase_data.data[sample_i] = (char)cell_phase_type;

                ht0_phase_type = cell_phase_type;
            }
            // subsequent seperators - only make sure they're consistent
            else
                ASSERT (cell_phase_type==ht0_phase_type, "Error: invalid VCF file - line %u - unable to parse sample %u: inconsistent phasing symbol '|' '/'", line_i, sample_i+1);
        }
    } // for characters in a sample

    // if the ploidy of the sample is lower than vb->ploidy, set missing ht as '*' ("ploidy padding")
    if (ploidy != vb->ploidy) {
        
        for (unsigned ht_i=ploidy; ht_i < vb->ploidy; ht_i++) 
            ht_data[ht_i] = '*';
        
        vb->add_bytes[SEC_HAPLOTYPE_DATA] += vb->ploidy - ploidy; // we added some characters that weren't in the original vcf
    }

    if (ploidy==1 && vb->ploidy > 1 && dl->phase_type == PHASE_MIXED_PHASED)
        vb->line_phase_data.data[sample_i] = (char)PHASE_HAPLO;

    return false;
}


static void segregate_genotype_area (VariantBlock *vb, DataLine *dl, 
                                     const char *cell_gt_data, 
                                     unsigned cell_gt_data_len,  // not include the \t or \n 
                                     unsigned line_i)
{
    // add memory if needed - but most chances are that we already have enough and realloc will return immediately
    vb->longest_line_genotype_data = MAX (vb->longest_line_genotype_data, vb->line_genotype_data.len + cell_gt_data_len + 1 /* \t */);

    // allocate buffer - this will allocate existing memory, or realloc it if needed
    buf_alloc (vb, &vb->line_genotype_data, vb->longest_line_genotype_data, 2, "line_genotype_data", line_i); 

    // case - empty genotype - just add \t
    if (!cell_gt_data_len) { 
        vb->line_genotype_data.data[vb->line_genotype_data.len++] = '\t';
        return;
    }

    // copy genotype data
    char *next = vb->line_genotype_data.data + vb->line_genotype_data.len;
    memcpy (next, cell_gt_data, cell_gt_data_len);
    next[cell_gt_data_len] = '\t';
 
    // if we are asked to optimize the GL subfield - zero-out the largest number
    if (dl->gl_subfield) 
        gl_optimize_do (vb, next, cell_gt_data_len, dl->gl_subfield - dl->has_haplotype_data);
    
    vb->line_genotype_data.len += cell_gt_data_len + 1;  // +1 for the \t
}

/* split variant line to: 
   1. variant data (all fields up to FORMAT field) + newline per sample block (if any samples exist):
   2. genotype data (except the 1|1) - space after every piece of data including last one
   3. haplotype data - a string of eg 1 or 0 (no whitespace) - ordered by the permutation order
   4. phase data - only if MIXED phase - a string of | and / - one per sample
*/
static bool segregate_data_line (VariantBlock *vb, /* may be NULL if testing */
                                 DataLine *dl, 
                                 unsigned line_i) // line in original VCF file
{
    // Caller must make sure all the fields of dl are 0 except for line and line_len
    
    dl->phase_type = PHASE_UNKNOWN;

    unsigned num_tabs = 0;
    unsigned sample_i = 0;
    
    enum {DONE, VARIANT, GENOTYPE, HAPLOTYPE} area = VARIANT;
    char *area_start = dl->line.data;
    long long pos_delta; // delta between POS in this row and the previous row
    const char *pos_start; unsigned pos_len; // start and length of pos field
    bool ploidy_overflow = false;

    char *c; for (c = dl->line.data; *c && area != DONE; c++) {

        // skip if we didn't arrive at a meaningful separator - tab, newline or : after a haplotype
        if (*c != '\t' && *c != '\n' && (*c != ':' || area != HAPLOTYPE)) continue;

        if (*c == '\t') num_tabs++;

        // handle the data in the area that we just passed
        switch (area) {
            case VARIANT:
                // we are at at the POS field - re-encode it to be the delta from the previous line
                if (num_tabs == 1) {
                    segregate_pos_field (vb, c+1, &pos_delta, &pos_start, &pos_len, line_i);
                    continue;
                }

                // we arrived at the FORMAT column - figure out if we have haplotype and genotype data for this line
                if (num_tabs == 8) { 
                    segregate_format_field(dl, c+1, line_i);                 
                    continue;
                }

                if (*c != '\n' && num_tabs < 9) continue; // if not yet the end of the variant data

                // we are at the end of the variant area
                segregate_variant_area(vb, dl, area_start, c - area_start, pos_delta, pos_start, pos_len);
                break;

            case HAPLOTYPE:
                ASSERT (dl->has_genotype_data || *c != ':', "Error: unexpected ':' after haplotype info, line %u sample %u", line_i, sample_i+1);

                ploidy_overflow = segregate_haplotype_area (vb, dl, area_start, c-area_start, line_i, sample_i);
                if (ploidy_overflow) goto cleanup; // we need to re-do this VB with the updated, higher, ploidy

                break;

            case GENOTYPE:
                segregate_genotype_area (vb, dl, area_start, c-area_start, line_i);
                break;

            case DONE: 
            default:
                ABORT0 ("in segregate_data_line(), should never reach here");
        }

        if (*c != ':') {
            // missing genotype data despite being declared in FORMAT - this is permitted by VCF spec
            if (area == HAPLOTYPE && dl->has_genotype_data) {
                segregate_genotype_area (vb, dl, NULL, 0, line_i);
                vb->add_bytes[SEC_GENOTYPE_DATA]++; // we added a \t that was not in the original file
                area = GENOTYPE; // we just processed an (empty) genotype area
            }

            if (area == HAPLOTYPE || area == GENOTYPE)
                sample_i++;
        }

        // calculate the next area
        area_start = c+1;

        if (*c == '\n') 
            area = DONE;
        else if (dl->has_haplotype_data && (area != HAPLOTYPE || !dl->has_genotype_data)) 
            area = HAPLOTYPE;
        else 
            area = GENOTYPE;
    }

    // some sanity checks
    ASSERT (!area || *c, "Error: last line %u is missing a newline", line_i);

    ASSERT (!area, "Error: line %u ends pre-maturely", line_i);

    ASSERT(sample_i == global_num_samples, 
           "Error: Invalid VCF file: the number of samples in line %u is %u, different than the VCF column header line which has %u samples",
           line_i, sample_i, global_num_samples);

    // update lengths
    if (dl->has_haplotype_data) {
        unsigned a = global_num_samples * vb->ploidy;
        vb->line_haplotype_data.len = global_num_samples * vb->ploidy;

        if (dl->phase_type == PHASE_MIXED_PHASED) 
            vb->line_phase_data.len = global_num_samples;
    }

    // now, overlay the data over the line memory so we can re-use our working buffers vb->line_* for the next line
    unsigned total_len = vb->line_variant_data.len + vb->line_haplotype_data.len + vb->line_genotype_data.len + vb->line_phase_data.len;
    buf_alloc (vb, &dl->line, total_len, 1, "dl->line", line_i);

    unsigned offset_in_line = 0;

    buf_overlay (&dl->variant_data, &dl->line, &vb->line_variant_data, &offset_in_line, "dl->variant_data", line_i);
    
    if (dl->has_haplotype_data)
        buf_overlay (&dl->haplotype_data, &dl->line, &vb->line_haplotype_data, &offset_in_line, "dl->haplotype_data", line_i);

    if (dl->has_haplotype_data && dl->phase_type == PHASE_MIXED_PHASED)
        buf_overlay (&dl->phase_data, &dl->line, &vb->line_phase_data, &offset_in_line, "dl->phase_data", line_i);    

    if (dl->has_genotype_data)
        buf_overlay (&dl->genotype_data, &dl->line, &vb->line_genotype_data, &offset_in_line, "dl->genotype_data", line_i);

cleanup:
    buf_free (&vb->line_variant_data);
    buf_free (&vb->line_haplotype_data);
    buf_free (&vb->line_phase_data);
    buf_free (&vb->line_genotype_data);

    return ploidy_overflow; 
}

void segregate_data_line_unit_test()
{
    unsigned *n = &global_num_samples;

    // *n=2 ; char *test_data = "20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\n"; // normal
    // *n=0 ; char *test_data = "20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\n";  // no samples (=ok)
    // *n=2 ; char *test_data = "20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1|0:3:5:65,3\n"; // inconsistent ploidy
    // *n=2 ; char *test_data = "20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t1:49:3:58,50\t0:3:5:65,3\n"; // genotype data only
    // *n=2 ; char *test_data = "20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0\t0|1\n"; // missing genotype data
    // *n=2 ; char *test_data = "20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT\t0|0:1\t0|1:0\n"; // unexpected genotype data
    // *n=2 ; char *test_data = "20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT\t0|0\t0/1\n"; // mixed phase
    // *n=2 ; char *test_data = "20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT\n"; // missing haplotype
    // *n=2 ; char *test_data = "20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tHT\n"; // missing genotype data
    // *n=1 ; char *test_data = "20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT\t0|0\t0/1\n"; // too many samples
    *n=3 ; char *test_data = "20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT\t0|0\t0/1\n"; // not enough samples
    
    DataLine *data_line = calloc (1, sizeof(DataLine));

    VariantBlock *vb = vb_get_vb(NULL, NULL, NULL, 0);

    buf_alloc (vb, &data_line->line, strlen(test_data) + 1, 0, "segregate_data_line_unit_test0", 0);
    strcpy (data_line->line.data, test_data);
    
    unsigned line_len = strlen(test_data) + 2; // maximum length of any part of the line (leave room for terminator eg \n or \t and a \0)

    unsigned line_i = 555;

    segregate_data_line (NULL, data_line, line_i);
    unsigned stophere=1;
}

static void segregate_update_vb_from_dl (VariantBlock *vb, DataLine *dl)
{
    // update block data
    vb->has_genotype_data = vb->has_genotype_data || dl->has_genotype_data;

    vb->has_haplotype_data = vb->has_haplotype_data || dl->has_haplotype_data;
    
    if (vb->phase_type == PHASE_UNKNOWN) 
        vb->phase_type = dl->phase_type;
    
    else if ((vb->phase_type == PHASE_PHASED     && dl->phase_type == PHASE_NOT_PHASED) ||
             (vb->phase_type == PHASE_NOT_PHASED && dl->phase_type == PHASE_PHASED) ||
              dl->phase_type  == PHASE_MIXED_PHASED)
        vb->phase_type = PHASE_MIXED_PHASED;    
}

// complete haplotype and genotypes of lines that don't have them, if we found out that some
// other line has them, and hence all lines in the vb must have them
void segregate_complete_missing_lines (VariantBlock *vb)
{
    vb->num_haplotypes_per_line = vb->ploidy * global_num_samples;
    
    for (unsigned vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {

        DataLine *dl = &vb->data_lines[vb_line_i];                                

        if (vb->has_haplotype_data && dl->has_haplotype_data)
            vb->add_bytes[SEC_STATS_HT_SEPERATOR] -= global_num_samples; // we saved this amount by not including the separator (tab or :) that appears after the haplotype data in the original VCF

        if (vb->has_haplotype_data && !dl->has_haplotype_data) {
            // realloc line which is the buffer on which haplotype data is overlaid
            buf_alloc (vb, &dl->line, dl->line.len + vb->num_haplotypes_per_line, 1, dl->line.name, dl->line.param);

            // overlay the haplotype buffer overlaid at the end of the line buffer
            buf_overlay (&dl->haplotype_data, &dl->line, NULL, &dl->line.len, "dl->haplotype_data", vb->first_line + vb_line_i);

            memset (dl->haplotype_data.data, '-', vb->num_haplotypes_per_line); // '-' means missing haplotype - note: ascii 45 (haplotype values start at ascii 48)
            dl->haplotype_data.len = vb->num_haplotypes_per_line;

            vb->add_bytes[SEC_HAPLOTYPE_DATA] += vb->num_haplotypes_per_line; // we added these bytes that are not in the orignial file

            dl->has_haplotype_data = true;
        }

        if (vb->has_genotype_data && !dl->has_genotype_data) {
            // realloc line which is the buffer on which genotype data is overlaid
            buf_alloc (vb, &dl->line, dl->line.len + global_num_samples, 1, dl->line.name, dl->line.param);

            // overlay the genotype buffer overlaid at the end of the line buffer
            buf_overlay (&dl->genotype_data, &dl->line, NULL, &dl->line.len, "dl->genotype_data", vb->first_line + vb_line_i);

            memset (dl->genotype_data.data, '\t', global_num_samples); // empty genotype data
            dl->genotype_data.len = global_num_samples;

            vb->add_bytes[SEC_GENOTYPE_DATA] += global_num_samples; // we added these bytes that are not in the orignial file

            dl->has_genotype_data = true;
        }
    }
}

// split each lines in this variant block to its components
void segregate_all_data_lines (VariantBlock *vb, Buffer *lines_orig /* for testing */)
{
    START_TIMER;

    unsigned num_ploidy_overlows = 0;

    for (int vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {

        DataLine *dl = &vb->data_lines[vb_line_i];

        if (lines_orig) buf_copy (vb, &lines_orig[vb_line_i], &dl->line, 0, dl->line.len+1); // if testing

        long long saved_last_pos = vb->last_pos;
        
        bool ploidy_overflow = segregate_data_line (vb, dl, vb->first_line + vb_line_i);
#       define MAX_PLOIDY_OVERFLOWS 1000 /* an arbitrary large number to avoid an ininifinate loop in case of a bug */

        // ploidy overflow is a situation where earlier samples in the VB have a low ploidy and then
        // a higher ploidy sample is encountered. we re-do the VB with the higher ploidy. Examples where this can happen:
        // 1. a VCF file contains eg a Y chromosome followed by an autosomal chromosome and a VB transcents both
        // 2. a X chromosome has male samples followed by female
        if (ploidy_overflow) {
            // just to be safe, we will limit the number of ploidy increases to 5
            ASSERT (num_ploidy_overlows < MAX_PLOIDY_OVERFLOWS, "Error: too many ploidy overflows, line_i=%u", vb->first_line + vb_line_i);
            vb_line_i = -1; // next line is 0
            num_ploidy_overlows++;
            vb->last_pos = saved_last_pos;

            // reset adders
            memset (vb->add_bytes, 0, sizeof(vb->add_bytes));

            continue;
        }

        segregate_update_vb_from_dl (vb, dl);
    }

    if (vb->has_genotype_data || vb->has_haplotype_data)
        segregate_complete_missing_lines(vb);

    COPY_TIMER(vb->profile.segregate_all_data_lines);
}