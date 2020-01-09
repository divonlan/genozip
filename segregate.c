// ------------------------------------------------------------------
//   segregate.c
//   Copyright (C) 2019-2020 Divon Lan <genozip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"

// returns true if this line has the same chrom as this VB, or if it is the first line
static bool seg_chrom_field (VariantBlock *vb, const char *str)
{
    unsigned i=0;
    if (! *vb->chrom) { // first line in this VB
        do {
            vb->chrom[i] = str[i];
            i++;
        } while (str[i-1] != '\t' && i < MAX_CHROM_LEN-1);
        vb->chrom[i-1] = 0;
       return true;
    }
 
    for (i=0; i < MAX_CHROM_LEN-1 && vb->chrom[i] && str[i] != '\t'; i++)
        if (vb->chrom[i] != str[i]) return false;

    return str[i] == '\t' && vb->chrom[i] == '\0'; // true if they are the same length
}

static void seg_pos_field (VariantBlock *vb, const char *str, 
                                 int64_t *pos_delta, const char **pos_start, unsigned *pos_len)
{
    *pos_start = str;

    // scan by ourselves - hugely faster the sscanf
    long long this_pos=0;
    const char *s; for (s=str; *s != '\t'; s++)
        this_pos = this_pos * 10 + (*s - '0');

    if (this_pos < vb->last_pos) 
        vb->is_sorted_by_pos = false;

    if ((int64_t)this_pos < vb->min_pos || vb->min_pos < 0)
        vb->min_pos = (int64_t)this_pos;

    if ((int64_t)this_pos > vb->max_pos)
        vb->max_pos = (int64_t)this_pos;

    *pos_len     = s - str;
    *pos_delta   = this_pos - vb->last_pos;
    vb->last_pos = this_pos;
}

// traverses the FORMAT field, gets ID of subfield, and moves to the next subfield
SubfieldIdType seg_get_subfield (const char **str, unsigned len, // remaining length of line
                                 unsigned line_i) // line in original vcf file
{
    SubfieldIdType subfield = EMPTY_SUBFIELD_ID;

    for (unsigned i=0; i < len; i++) {
        if ((*str)[i] != ':' && (*str)[i] != '\t' && (*str)[i] != '\n') { // note: in vcf files, we will see \t at end of FORMAT, but in variant data sections of piz files, we will see a \n
            
            if (i < SUBFIELD_ID_LEN)
                subfield.id[i] = (*str)[i];

        } else {
            *str += i+1;
            return subfield;
        }
    }
    // we reached the end of the line without encountering a : or \t OR this is piz and we reached \n
    ASSERT ((*str)[-1]=='\n', "Error: invalid FORMAT field in line %u", line_i); 
    return subfield; // never reach here, just to avoid a compiler warning
}

static void seg_format_field(VariantBlock *vb, DataLine *dl, 
                             const char *str, unsigned len, // remaining length of line
                             unsigned line_i) // line in original vcf file
{
    dl->num_subfields = 0; // in case of re-segmenting after ploidy overflow it might already be set and we need to reset it

    if (!global_num_samples) return; // if we're not expecting any samples, we ignore the FORMAT field

    ASSERT(str[0] && str[1] && str[2], "Error: missing samples in line %u, expected by VCF file sample name header line", line_i);

    if (str[0] == 'G' && str[1] == 'T') // GT field in FORMAT columns - must always appear first per VCF spec (if at appears)
        dl->has_haplotype_data = true; 

    if (str[2] == ':' || !dl->has_haplotype_data) {
        dl->has_genotype_data = true; // no GT subfield, or subfields in addition to GT

        if (dl->has_haplotype_data) str +=3; // skip over GT

        do {
            SubfieldIdType subfield = seg_get_subfield (&str, len, line_i);

            unsigned sf_i = mtf_get_sf_i_by_subfield (vb->mtf_ctx, &vb->num_subfields, subfield);

            dl->sf_i[dl->num_subfields++] = sf_i;
        } 
        while (str[-1] != '\t' && str[-1] != '\n');
    }

    if (dl->has_genotype_data)
        buf_alloc (vb, &vb->line_gt_data, dl->num_subfields * global_num_samples * sizeof(uint32_t), 1, "line_gt_data", line_i); 
}

static void seg_variant_area(VariantBlock *vb, DataLine *dl, const char *str, unsigned len,
                                   int64_t pos_delta, const char *pos_start, unsigned pos_len)
{
    char pos_delta_str[30];
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

static bool seg_haplotype_area(VariantBlock *vb, DataLine *dl, const char *str, unsigned len, unsigned line_i, unsigned sample_i)
{
    // check ploidy
    unsigned ploidy=1;
    for (unsigned i=1; i<len-1; i++)
        if (str[i] == '|' || str[i] == '/') ploidy++;

    ASSERT (ploidy <= MAX_PLOIDY, "Error: ploidy=%u exceeds the maximum of %u in line_i=%u", ploidy, MAX_PLOIDY, line_i);
    
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
    if (!vb->line_ht_data.data || vb->line_ht_data.size < vb->ploidy * global_num_samples) {
        buf_alloc (vb, &vb->line_ht_data, vb->ploidy * global_num_samples, 1, "line_ht_data", line_i);
        dl->phase_type = (vb->ploidy==1 ? PHASE_HAPLO : PHASE_UNKNOWN);
    }

    char *ht_data = &vb->line_ht_data.data[vb->ploidy * sample_i];

    vb->add_bytes[SEC_PHASE_DATA] -= (ploidy-1); // we "saved" this number of phase characters vs the vcf file (... we will add them back later if we didn't)

    bool padding_only = (*str == '*');

    PhaseType ht0_phase_type = PHASE_UNKNOWN;
    for (unsigned ht_i=0; ht_i < ploidy; ht_i++) {

        char ht = *(str++); 
        len--;

        ASSERT ((ht >= '0' && ht <= '9') || ht == '.' || ht == '*',
                "Error: invalid VCF file - line %u - expecting an allele in a sample to be a number 0-9 or . , but seeing %c", line_i, *str);

        // single-digit allele numbers
        ht_data[ht_i] = ht;

        if (!len) break;

        // handle 2-digit allele numbers
        if (ht != '.' && *str >= '0' && *str <= '9') {
            unsigned allele = 10 * (ht-'0') + (*(str++) - '0');
            len--;

            // make sure there isn't a 3rd digit
            ASSERT (!len || *str < '0' || * str > '9', "Error: VCF file - line %u sample %u - genozip currently supports only alleles up to 99", line_i, sample_i+1);

            ht_data[ht_i] = '0' + allele; // use ascii 48->147
            vb->add_bytes[SEC_HAPLOTYPE_DATA]--; // we saved one byte vs the original file
        }

        // get phase
        if (ploidy > 1 && ht_i < ploidy-1) {
            
            PhaseType cell_phase_type = *(str++);
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

    if (padding_only)
        vb->add_bytes[SEC_HAPLOTYPE_DATA] += 1; // we adding a '*' that was not in the original file

    if (ploidy==1 && vb->ploidy > 1 && dl->phase_type == PHASE_MIXED_PHASED)
        vb->line_phase_data.data[sample_i] = (char)PHASE_HAPLO;

    return false;
}

// returns length of the snip, not including the separator
static inline unsigned seg_snip_len_tnc (const char *snip)
{
    const char *s; for (s=snip ; *s != '\t' && *s != ':' && *s != '\n'; s++);
    return s - snip;
}

static void seg_genotype_area (VariantBlock *vb, DataLine *dl, 
                               const char *cell_gt_data, 
                               unsigned cell_gt_data_len,  // not including the \t or \n 
                               unsigned line_i)
{
    uint32_t *next = (uint32_t *)(&vb->line_gt_data.data[vb->line_gt_data.len]);
 
    bool end_of_cell = !cell_gt_data_len;
    for (unsigned sf=0; sf < dl->num_subfields; sf++) { // iterate on the order as in the line

        // move next to the beginning of the subfield data, if there is any
        unsigned len = end_of_cell ? 0 : seg_snip_len_tnc (cell_gt_data);

        MtfContext *ctx = &vb->mtf_ctx[dl->sf_i[sf]];
        uint32_t node_index = mtf_evaluate_snip (vb, ctx, cell_gt_data, len);
        *(next++) = node_index;

        if (node_index != SEG_MISSING_SF) // don't skip the \t if we might have more missing subfields
            cell_gt_data += len + 1; // skip separator too

        end_of_cell = end_of_cell || cell_gt_data[-1] != ':'; // a \t or \n encountered

        if (vb->variant_block_i == 1 && node_index <= SEG_MAX_INDEX) 
            ((MtfNode *)ctx->mtf.data)[node_index].count++;
    }
    ASSERT0 (end_of_cell, "Error: invalid reading of genotype data");

    vb->line_gt_data.len += dl->num_subfields * sizeof(uint32_t);  // len is number of bytes

    vb->add_bytes[SEC_GENOTYPE_DATA] += dl->num_subfields * sizeof(uint32_t) - (cell_gt_data_len ? cell_gt_data_len+1 : 0);
}

/* split variant line to: 
   1. variant data (all fields up to FORMAT field) + newline per sample block (if any samples exist):
   2. genotype data (except the 1|1) - space after every piece of data including last one
   3. haplotype data - a string of eg 1 or 0 (no whitespace) - ordered by the permutation order
   4. phase data - only if MIXED phase - a string of | and / - one per sample
*/
static bool seg_data_line (VariantBlock *vb, /* may be NULL if testing */
                           DataLine *dl, 
                           unsigned line_i) // line in original VCF file
{
    // Caller must make sure all the fields of dl are 0 except for line and line_len
    
    dl->phase_type = PHASE_UNKNOWN;

    unsigned num_tabs = 0;
    unsigned sample_i = 0;
    
    enum {DONE, VARIANT, GENOTYPE, HAPLOTYPE} area = VARIANT;
    char *area_start = dl->line.data;
    int64_t pos_delta=0; // delta between POS in this row and the previous row
    const char *pos_start=NULL; unsigned pos_len=0; // start and length of pos field
    bool ploidy_overflow = false;
    unsigned gt_line_len=0;

    char *c; for (c = dl->line.data; *c && area != DONE; c++) {

        // skip if we didn't arrive at a meaningful separator - tab, newline or : after a haplotype
        if (c != dl->line.data && *c != '\t' && *c != '\n' && (*c != ':' || area != HAPLOTYPE)) continue;

        if (*c == '\t' || *c == '\n') num_tabs++;

        // handle the data in the area that we just passed
        switch (area) {
            case VARIANT:
                // we are at the CHROM field - make sure the VB has consisent CHROM
                if (num_tabs == 0) {
                    //bool same_chrom = TO DO: finish this vb and start a new one when chrom changes
                    seg_chrom_field (vb, c);
                    continue;
                }

                // we are at at the POS field - re-encode it to be the delta from the previous line
                if (num_tabs == 1) {
                    seg_pos_field (vb, c+1, &pos_delta, &pos_start, &pos_len);
                    continue;
                }

                // we arrived at the FORMAT column - figure out if we have haplotype and genotype data for this line
                if (num_tabs == 8) { 
                    seg_format_field (vb, dl, c+1, dl->line.len - ((c+1)-dl->line.data),line_i);                 
                    continue;
                }

                if (*c != '\n' && num_tabs < 9) continue; // if not yet the end of the variant data

                // we are at the end of the variant area
                seg_variant_area(vb, dl, area_start, c - area_start, pos_delta, pos_start, pos_len);
                break;

            case HAPLOTYPE:
                ASSERT (dl->has_genotype_data || *c != ':', "Error: unexpected ':' after haplotype info, line %u sample %u", line_i, sample_i+1);

                ploidy_overflow = seg_haplotype_area (vb, dl, area_start, c-area_start, line_i, sample_i);
                if (ploidy_overflow) goto cleanup; // we need to re-do this VB with the updated, higher, ploidy

                break;

            case GENOTYPE:
                seg_genotype_area (vb, dl, area_start, c-area_start, line_i);
                gt_line_len += c-area_start + 1; // including the \t or \n
                break;

            case DONE: 
            default:
                ABORT0 ("in seg_data_line(), should never reach here");
        }

        if (*c != ':') {
            // missing genotype data despite being declared in FORMAT - this is permitted by VCF spec
            if (area == HAPLOTYPE && dl->has_genotype_data) {
                seg_genotype_area (vb, dl, NULL, 0, line_i);
                gt_line_len++; // adding the SEG_MISSING_SF
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

    // some real-world files I encountered have too-short lines due to human errors. we pad them
    if (sample_i < global_num_samples) {
        
        ASSERTW(false, "Warning: the number of samples in line %u is %u, different than the VCF column header line which has %u samples",
                line_i, sample_i, global_num_samples);

        for (; sample_i < global_num_samples; sample_i++) {

            if (dl->has_haplotype_data) {
                // '*' (haplotype padding) with ploidy 1
                seg_haplotype_area (vb, dl, "*", 1, line_i, sample_i);
                vb->add_bytes[SEC_HAPLOTYPE_DATA]++;
            }

            if (dl->has_genotype_data) {
                seg_genotype_area (vb, dl, NULL, 0, line_i);
                gt_line_len++; // adding the SEG_MISSING_SF
            }

        }
    }

    // update lengths
    if (dl->has_haplotype_data) {
        vb->line_ht_data.len = global_num_samples * vb->ploidy;

        if (dl->phase_type == PHASE_MIXED_PHASED) 
            vb->line_phase_data.len = global_num_samples;
    } else 
        vb->line_ht_data.len = 0;

    vb->max_gt_line_len = MAX (vb->max_gt_line_len, gt_line_len);

    // now, overlay the data over the line memory so we can re-use our working buffers vb->line_* for the next line
    unsigned total_len = vb->line_variant_data.len + vb->line_ht_data.len + vb->line_gt_data.len + vb->line_phase_data.len;
    buf_alloc (vb, &dl->line, total_len, 1, "dl->line", line_i);

    unsigned offset_in_line = 0;

    buf_overlay (&dl->variant_data, &dl->line, &vb->line_variant_data, &offset_in_line, "dl->variant_data", line_i);

    if (dl->has_haplotype_data) {
        buf_overlay (&dl->haplotype_data, &dl->line, &vb->line_ht_data, &offset_in_line, "dl->haplotype_data", line_i);
        
        if (flag_show_alleles)
            printf ("%.*s\n", dl->haplotype_data.len, dl->haplotype_data.data);
    }

    if (dl->has_haplotype_data && dl->phase_type == PHASE_MIXED_PHASED)
        buf_overlay (&dl->phase_data, &dl->line, &vb->line_phase_data, &offset_in_line, "dl->phase_data", line_i);    

    if (dl->has_genotype_data) 
        buf_overlay (&dl->genotype_data, &dl->line, &vb->line_gt_data, &offset_in_line, "dl->genotype_data", line_i);

cleanup:
    buf_free (&vb->line_variant_data);
    buf_free (&vb->line_ht_data);
    buf_free (&vb->line_phase_data);
    buf_free (&vb->line_gt_data);

    // roll back in case of ploidy overflow
    if (ploidy_overflow) {
        vb->last_pos = 0;
        memset (vb->add_bytes, 0, sizeof(vb->add_bytes));
    }

    return ploidy_overflow; 
}

static void seg_update_vb_from_dl (VariantBlock *vb, DataLine *dl)
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
void seg_complete_missing_lines (VariantBlock *vb)
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

            // NOTE: we DONT set dl->has_haplotype_data to true bc downstream we still
            // count this row as having no GT field when analyzing gt data
        }

        if (vb->has_genotype_data && !dl->has_genotype_data) {
            // realloc line which is the buffer on which haplotype data is overlaid
            buf_alloc (vb, &dl->line, dl->line.len + global_num_samples * sizeof(uint32_t), 1, dl->line.name, dl->line.param);

            // overlay the haplotype buffer overlaid at the end of the line buffer
            buf_overlay (&dl->genotype_data, &dl->line, NULL, &dl->line.len, "dl->genotype_data", vb->first_line + vb_line_i);
            for (unsigned i=0; i < global_num_samples; i++) 
                ((uint32_t*)dl->genotype_data.data)[i] = SEG_MISSING_SF;

            dl->genotype_data.len += global_num_samples * sizeof(uint32_t);
            dl->num_subfields = 1;

            vb->add_bytes[SEC_GENOTYPE_DATA] += global_num_samples * sizeof(uint32_t); // we added these bytes that are not in the orignial file
        }
    }
}

// split each lines in this variant block to its components
void seg_all_data_lines (VariantBlock *vb, Buffer *lines_orig /* for testing */)
{
    START_TIMER;

    unsigned num_ploidy_overlows = 0;

    for (unsigned vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {
        //printf ("vb_line_i=%u\n", vb_line_i);
        DataLine *dl = &vb->data_lines[vb_line_i];

        if (lines_orig) buf_copy (vb, &lines_orig[vb_line_i], &dl->line, 1, 0, dl->line.len+1, "lines_orig", vb->variant_block_i); // if testing
        bool ploidy_overflow = seg_data_line (vb, dl, vb->first_line + vb_line_i);
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

            continue;
        }

        seg_update_vb_from_dl (vb, dl);
    }

    if (/*vb->has_genotype_data || */vb->has_haplotype_data)
        seg_complete_missing_lines(vb);
        
    COPY_TIMER(vb->profile.seg_all_data_lines);
}