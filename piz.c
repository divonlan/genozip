// ------------------------------------------------------------------
//   piz.c
//   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
//   Please see terms and conditions in the file LICENSE.txt

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <inttypes.h>

#include "vczip.h"

// decode the delta-encouded value of the POS field
static inline void piz_decode_pos (VariantBlock *vb,
                                   const char *str,
                                   char *pos_str, const char **delta_pos_start, unsigned *delta_pos_len, int *add_len /* out */)
{
    START_TIMER;

    *delta_pos_start = str;

    long long delta=0;

    // we parse the string ourselves - this is hugely faster than sscanf.
    unsigned negative = *str == '-'; // 1 or 0

    const char *s; for (s=(str + negative); *s != '\t' ; s++)
        delta = delta*10 + (*s - '0');

    if (negative) delta = -delta;

    vb->last_pos += delta;
    sprintf (pos_str, "%"PRIu64, vb->last_pos);
    
    *delta_pos_len = s - str;
    *add_len = strlen (pos_str) - *delta_pos_len;

    COPY_TIMER(vb->profile.piz_decode_pos);
}

static void piz_get_variant_data_line (VariantBlock *vb, unsigned line_i,
                                       unsigned *length_remaining, // for safety
                                       const char **line_start,
                                       unsigned *gl_subfield)
{
    START_TIMER;

    // decoding the delta-encouded value of the POS field
    const char *delta_pos_start; // start of delta-POS field
    unsigned delta_pos_len; // length of encoded delta-pos
    int add_len; // how much more (or less) characters do we need in the string after decoding POS?
    char pos_str[22]; // decoded POS value

    const char *next = *line_start;
    const char *after = *line_start + *length_remaining;
    unsigned column = 1;

    while (next < after) {

        // starting a new field
        if (*next == '\t') {
            column++;

            // if we're at the POS field, decode the delta encoding
            if (column == 2)  
                piz_decode_pos (vb, next+1, pos_str, &delta_pos_start, &delta_pos_len, &add_len);

            // if this is the FORMAT column, and needed, check which subfield (if any) is GL
            else if (column == 9) 
                *gl_subfield = gl_optimize_get_gl_subfield_index (next+1);
        }

        else if (*next == '\n') {
            next++; // past the newline
            unsigned line_strlen = next - *line_start; // inc. the newline, not including \0

            buf_alloc (vb, &vb->line_variant_data, line_strlen + add_len + 1 /* +1 sprintf adds a \0 */, 
                         1, "line_variant_data", 0);

            const char *after_delta = delta_pos_start + delta_pos_len;

            sprintf (vb->line_variant_data.data, "%.*s%s%.*s",
                     delta_pos_start - *line_start, *line_start,    // substring until \t preceding delta
                     pos_str,                           // decoded pos string
                     next - after_delta, after_delta);  // substring starting \t after delta
            
            vb->line_variant_data.len = line_strlen + add_len;

            *line_start = next;
            return;
        }

        next++;
    }
    
    *length_remaining = after - next;

    ASSERT (false, "Error: corrupt dv file - at end of variant_data buffer, and no newline was found", "");

    COPY_TIMER(vb->profile.piz_get_variant_data_line);
}

unsigned piz_get_genotype_sample_starts (VariantBlock *vb)
{
    buf_alloc (vb, &vb->next_gt_in_sample, sizeof(char*) * global_num_samples, 1, "next_in_sample", 0);
    const char **first_line_in_sample = (const char **)vb->next_gt_in_sample.data; // an array of char * - each element pointing to the gt of the first line of a sample

    buf_alloc (vb, &vb->gt_line_lengths, sizeof(unsigned) * vb->num_lines, 1, "gt_line_lengths", 0);
    unsigned *gt_line_lengths = (unsigned *)vb->gt_line_lengths.data; // an array of unsigned  - each element is the length of gts in a line
    memset (gt_line_lengths, 0, sizeof(unsigned) * vb->num_lines);
    
    unsigned sample_i=0;
    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned num_samples_in_sample_block = vb_num_samples_in_sb (vb, sb_i);
        
        const char *next_gt_in_sb = vb->genotype_sections_data[sb_i].data;
        const char *beyond_last_gt_in_sb = vb->genotype_sections_data[sb_i].data + vb->genotype_sections_data[sb_i].len;

        for (unsigned sample_in_sb_i=0; sample_in_sb_i < num_samples_in_sample_block; sample_in_sb_i++) {
            
            first_line_in_sample[sample_i++] = next_gt_in_sb;
            
            // now skip all remaining genotypes of lines in the variant block (each terminated by a \t)
            // until arriving at the genotype of first line of the next sample
            unsigned line_i; for (line_i=0; line_i < vb->num_lines && next_gt_in_sb < beyond_last_gt_in_sb; next_gt_in_sb++) {

                gt_line_lengths[line_i]++;

                if (*next_gt_in_sb == '\t')
                    line_i++;
            }

            // verify that we found a genotype of each line
            ASSERT (line_i==vb->num_lines, "Error: expected to find %u genotypes, but found only %u. sb_i=%u, first_line=%u",
                    vb->num_lines * num_samples_in_sample_block, vb->num_lines * (sb_i * vb->num_samples_per_block + sample_in_sb_i) + line_i,
                    sb_i, vb->first_line); 
        }

        // verify that there is no data left over after getting all the genotypes of this sample block
        ASSERT (next_gt_in_sb == beyond_last_gt_in_sb, "Error: expected to find %genotypes, but found more. sb_i=%u, first_line=%u",
                vb->num_lines * num_samples_in_sample_block, sb_i, vb->first_line);
    }

    // find maximum line length
    unsigned max_line_len = 0;
    for (unsigned line_i=0; line_i < vb->num_lines; line_i++) 
        if (gt_line_lengths[line_i] > max_line_len)
            max_line_len = gt_line_lengths[line_i];

    return max_line_len;
}

// convert genotype data from sample block format to line format
static void piz_get_genotype_data_line (VariantBlock *vb, unsigned line_i, unsigned gl_subfield)
{
    START_TIMER;

    for (unsigned sample_i=0; sample_i < global_num_samples; sample_i++) {

        unsigned len = strcpy_tab (&vb->line_genotype_data.data[vb->line_genotype_data.len],
                                   ((char **)vb->next_gt_in_sample.data)[sample_i]);
        
        if (gl_subfield >= 1) 
            // if gl-optimzed, 90% of the time of this routine, adn 15% of the overall decompress time are in this routine
            gl_optimize_undo (vb, &vb->line_genotype_data.data[vb->line_genotype_data.len], len, 
                              gl_subfield - vb->has_haplotype_data);

        ((char **)vb->next_gt_in_sample.data)[sample_i] += len;
        vb->line_genotype_data.len += len;
    }

    vb->data_lines[line_i].has_genotype_data = vb->line_genotype_data.len > global_num_samples; // not all just \t

    COPY_TIMER(vb->profile.piz_get_genotype_data_line);
}

static void piz_get_phase_data_line (VariantBlock *vb, unsigned line_i)
{
    START_TIMER;

    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned num_samples_in_sample_block = vb_num_samples_in_sb (vb, sb_i);

        memcpy (&vb->line_phase_data.data[sb_i * vb->num_samples_per_block],
                &vb->phase_sections_data[sb_i].data[line_i * num_samples_in_sample_block], 
                num_samples_in_sample_block);
    }

    COPY_TIMER(vb->profile.piz_get_phase_data_line);
}

// for each haplotype column, retrieve its it address in the haplotype sections. Note that since the haplotype sections are
// transposed, each column will be a row, or a contiguous array, in the section data. This function returns an array
// of pointers, each pointer being a beginning of column data within the section array
static const char **piz_get_ht_columns_data (VariantBlock *vb)
{
    START_TIMER;

    buf_alloc (vb, &vb->ht_columns_data, sizeof (char *) * (vb->num_haplotypes_per_line + 7), 1, "ht_columns_data", 0); // realloc for exact size

    const char **ht_columns_data = (const char **)vb->ht_columns_data.data;

    const unsigned *permutatation_index = (const unsigned *)vb->haplotype_permutation_index.data;
    unsigned max_ht_per_block = vb->num_samples_per_block * vb->ploidy; // last sample block may have less, but that's ok for our div/mod calculations below

    for (unsigned ht_i=0; ht_i < vb->num_haplotypes_per_line; ht_i++) {
        unsigned permuted_ht_i = permutatation_index[ht_i];
        unsigned sb_i    = permuted_ht_i / max_ht_per_block; // get haplotype sample block per this ht is
        unsigned row     = permuted_ht_i % max_ht_per_block; // get row transposed haplotype sample block. column=line_i
        unsigned sb_ht_i = row * vb->num_lines;     // index within haplotype block 

        ht_columns_data[ht_i] = &vb->haplotype_sections_data[sb_i].data[sb_ht_i];
    }

    // provide 7 extra zero-columns for the convenience of the permuting loop (suppoting 32bit and 64bit assignments)
    static char *column_of_zeros = NULL; // this static is allocated once and never changed, so no thread safety issues here
    if (!column_of_zeros) column_of_zeros = calloc (VARIANTS_PER_BLOCK, 1);

    for (unsigned ht_i=vb->num_haplotypes_per_line; ht_i < vb->num_haplotypes_per_line + 7; ht_i++)
        ht_columns_data[ht_i] = column_of_zeros;

    COPY_TIMER(vb->profile.piz_get_ht_permutation_lookups);

    return ht_columns_data;
}

// build haplotype for a line - reversing the permutation and the transposal.
static void piz_get_haplotype_data_line (VariantBlock *vb, unsigned line_i, const char **ht_columns_data)
{
    START_TIMER;

    // this loop can consume up to 50% of the entire decompress compute time (tested with 1KGP data)
    // note: we do memory assignment 32 bit at time (its about 10% faster than byte-by-byte)
    // TO DO: this is LITTLE ENDIAN, need an #ifdef to support BIG ENDIAN
    // TO DO: we can also #ifdef between 32bit and 64bit compilers and do 8 at a time for the latter
    unsigned *next = (unsigned *)vb->line_haplotype_data.data;
    for (unsigned ht_i=0; ht_i < vb->num_haplotypes_per_line; ht_i += 4) 
        *(next++) = ((unsigned)(unsigned char)ht_columns_data[ht_i    ][line_i]      ) |  // this is LITTLE ENDIAN order
                    ((unsigned)(unsigned char)ht_columns_data[ht_i + 1][line_i] << 8 ) |
                    ((unsigned)(unsigned char)ht_columns_data[ht_i + 2][line_i] << 16) |
                    ((unsigned)(unsigned char)ht_columns_data[ht_i + 3][line_i] << 24) ;  // no worries if num_haplotypes_per_line is not a multiple of 4 - we have extra columns of zero

    // check if this row has now haplotype data (no GT field) despite some other rows in the VB having data
    DataLine *dl = &vb->data_lines[line_i];
    dl->has_haplotype_data = vb->line_haplotype_data.data[0] != '-'; // either the entire line is '-' or there is no '-' in the line

    COPY_TIMER(vb->profile.piz_get_haplotype_data_line);
}

// merge line components (variant, haplotype, genotype, phase) back into a line
static void piz_merge_line(VariantBlock *vb, unsigned line_i)
{
    DataLine *dl = &vb->data_lines[line_i]; 

    // calculate the line length & allocate it
    unsigned ht_digits_len  = dl->has_haplotype_data ? vb->num_haplotypes_per_line : 0; 
    unsigned var_data_len   = vb->line_variant_data.len;        // includes a \n separator
    unsigned phase_sepr_len = dl->has_haplotype_data ? global_num_samples * (vb->ploidy-1) : 0; // the phase separators (/ or |)
    unsigned gt_colon_len   = ((dl->has_genotype_data && dl->has_haplotype_data) ? global_num_samples : 0); // the colon separating haplotype data from genotype data 
    unsigned gt_data_len    = (dl->has_genotype_data ? vb->line_genotype_data.len // includes accounting for separator (\t or \n) after each sample
                                                     : (global_num_samples ? global_num_samples // separators after haplotype data with no genotype data
                                                                           : 0));  //  only variant data, no samples
    // adjustments to lengths
    if (dl->has_haplotype_data) {
        for (unsigned i=0; i < vb->num_haplotypes_per_line; i++) {
            // adjust for 2-digit alleles represented by 'A'..'Z' in the haplotype data
            if ((unsigned char)vb->line_haplotype_data.data[i] >= '0'+10) 
                ht_digits_len++; // 2-digit haplotype (10 to 99)

            // adjust for samples with ploidy less than vb->ploidy 
            if ((unsigned char)vb->line_haplotype_data.data[i] == '*') { // '*' means "ploidy padding"
                ht_digits_len--;
                phase_sepr_len--;
            }
        }
    }

    dl->line.len = var_data_len + ht_digits_len + phase_sepr_len + gt_colon_len + gt_data_len; // assigning line->line.len before allocating line->line should work
    buf_alloc (vb, &dl->line, dl->line.len + 2, 1, "dl->line", vb->first_line + line_i); // +1 for string terminator, +1 for temporary additonal phase in case of '*'

    char *next    = dl->line.data;
    char *next_gt = vb->line_genotype_data.data;

    // add variant data - change the \n separator to \t if needed
    memcpy (next, vb->line_variant_data.data, vb->line_variant_data.len);

    if (dl->has_genotype_data || dl->has_haplotype_data) 
        next[vb->line_variant_data.len-1] = '\t';

    next += vb->line_variant_data.len;

    // add samples
    for (unsigned sample_i=0; sample_i < global_num_samples; sample_i++) {
        // add haplotype data - ploidy haplotypes per sample 
        if (dl->has_haplotype_data) {

            PhaseType phase = (vb->phase_type == PHASE_MIXED_PHASED ? vb->line_phase_data.data[sample_i]
                                                                    : vb->phase_type);
            ASSERT (phase=='/' || phase=='|' || phase=='*', "Error: invalid phase %c line_i=%u sample_i=%u", phase, line_i, sample_i+1);

            for (unsigned p=0; p < vb->ploidy ; p++) {
                
                unsigned char ht = vb->line_haplotype_data.data[sample_i * vb->ploidy + p];
                if (ht == '.') // unknown 
                    *(next++) = ht;

                else if (ht == '*')  // missing haplotype - delete previous phase
                    *(--next) = 0;

                else { // allele 0 to 99
                    unsigned allele = ht - '0'; // allele 0->99 represented by ascii 48->147
                    ASSERT (allele <= 99, "Error: allele out of range: %u line_i=%u sample_i=%u", allele, line_i, sample_i+1);
                    
                    if (allele >= 10) *(next++) = '0' + allele / 10;
                    *(next++) = '0' + allele % 10;
                }

                // add the phase character between the haplotypes
                if (vb->ploidy >= 2 && p < vb->ploidy-1) 
                    *(next++) = phase;
            }
        }

        // get length of genotype data - we need it now, to see if we need a :
        unsigned gt_len;
        if (dl->has_genotype_data) {
            char *tab = strchr (next_gt, '\t');
            ASSERT (tab, "Error: has_genotype_data=true, but cannot find a tab separator between genotype elements in vb->line_genotype_data", "")

            gt_len = tab - next_gt;
        }

        // add colon separating haplotype from genotype data
        if (dl->has_haplotype_data && dl->has_genotype_data) {
            if (gt_len)
                *(next++) = ':';
            else
                // this sample is in a line with genotype data, but has not genotype data. This is permitted
                // by VCF spce. We adjust the line length to account for this - we couldn't have known in the beginning without 
                // expensive scanning of the line
                dl->line.len--;            
        }

        // add genotype data - dropping the \t
        if (dl->has_genotype_data) {
            memcpy (next, next_gt, gt_len);
            next_gt += gt_len + 1;
            next += gt_len;
        }

        // add tab or newline separator after sample data
        if (dl->has_haplotype_data || dl->has_genotype_data)
            *(next++) = (sample_i == global_num_samples - 1) ? '\n' : '\t';
    }
    *next = '\0'; // end of string;

    // sanity check
    ASSERT (next - dl->line.data == dl->line.len, "Error: unexpected line size: calculated=%u, actual=%u", dl->line.len, next - dl->line.data);
}

// combine all the sections of a variant block to regenerate the variant_data, haplotype_data,
// genotype_data and phase_data for each row of the variant block
void piz_reconstruct_line_components (VariantBlock *vb)
{
    START_TIMER;

    // initialize variant_data stuff
    const char *variant_data_next_line = vb->variant_data_section_data.data;
    unsigned variant_data_length_remaining = vb->variant_data_section_data.len;
    
    // initialize phase data if needed
    if (vb->phase_type == PHASE_MIXED_PHASED) 
        buf_alloc (vb, &vb->line_phase_data, global_num_samples, 1, "line_phase_data", 0);

    // initialize haplotype stuff
    const char **ht_columns_data;
    if (vb->has_haplotype_data) {

        //  memory - realloc for exact size, add 7 because depermuting_loop works on a word (32/64 bit) boundary
        buf_alloc (vb, &vb->line_haplotype_data, vb->num_haplotypes_per_line + 7, 1, "line_haplotype_data", 0);

        ht_columns_data = piz_get_ht_columns_data (vb);
    }

    // initialize genotype stuff
    if (vb->has_genotype_data) {
        unsigned max_line_len = piz_get_genotype_sample_starts(vb);

        buf_alloc (vb, &vb->line_genotype_data, max_line_len, 1, "line_genotype_data", 0);
    }

    for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {

        // reset len for next line - no need to realloc as we have realloced what is needed already
        vb->line_haplotype_data.len = vb->line_genotype_data.len = vb->line_phase_data.len = 0;

        // de-permute variant data into vb->line_variant_data
        unsigned gl_subfield; 
        piz_get_variant_data_line (vb, line_i, &variant_data_length_remaining, &variant_data_next_line, &gl_subfield);

        // transform sample blocks (each block: n_lines x s_samples) into line components (each line: 1 line x ALL_samples)
        if (vb->has_genotype_data) 
            piz_get_genotype_data_line (vb, line_i, gl_subfield);

        if (vb->phase_type == PHASE_MIXED_PHASED) 
            piz_get_phase_data_line (vb, line_i);

        if (vb->has_haplotype_data) 
            piz_get_haplotype_data_line (vb, line_i, ht_columns_data);

        piz_merge_line (vb, line_i);
    }

    COPY_TIMER(vb->profile.piz_reconstruct_line_components);
}

static void piz_uncompress_all_sections (VariantBlock *vb)
{
    unsigned *section_index = (unsigned *)vb->z_section_headers.data;

    // get the variant data - newline-seperated lines, each containing the first 8 (if no FORMAT field) or 9 fields (if FORMAT exists)
    zfile_uncompress_section (vb, vb->z_data.data + section_index[0], &vb->variant_data_section_data, SEC_VARIANT_DATA);
    
    SectionHeaderVariantData *vardata_header = (SectionHeaderVariantData *)(vb->z_data.data + section_index[0]);
    vb->num_lines               = ENDN16 (vardata_header->num_lines);
    vb->first_line              = ENDN32 (vardata_header->first_line);
    vb->has_genotype_data       = vardata_header->has_genotype_data;
    vb->has_haplotype_data      = vardata_header->num_haplotypes_per_line > 0;
    vb->phase_type              = vardata_header->phase_type;
    vb->num_haplotypes_per_line = ENDN32 (vardata_header->num_haplotypes_per_line);
    vb->num_samples_per_block   = ENDN32 (vardata_header->num_samples_per_block);
    vb->ploidy                  = ENDN16 (vardata_header->ploidy);
    vb->num_sample_blocks       = ENDN32 (vardata_header->num_sample_blocks);
    vb->vcf_data_size           = ENDN32 (vardata_header->vcf_data_size);
        
    ASSERT (global_num_samples == ENDN32 (vardata_header->num_samples), "Error: Expecting variant block to have %u samples, but it has %u", global_num_samples, vardata_header->num_samples);

    // we allocate memory for the Buffer arrays only once the first time this VariantBlock
    // is used. Subsequent blocks reusing the memory will have the same number of samples (by VCF spec)
    if (vardata_header->has_genotype_data && !vb->genotype_sections_data) 
        vb->genotype_sections_data  = calloc (vb->num_sample_blocks, sizeof (Buffer));

    if (vardata_header->phase_type == PHASE_MIXED_PHASED && !vb->phase_sections_data) 
        vb->phase_sections_data     = calloc (vb->num_sample_blocks, sizeof (Buffer));
    
    if (vardata_header->num_haplotypes_per_line && !vb->haplotype_sections_data) 
        vb->haplotype_sections_data = calloc (vb->num_sample_blocks, sizeof (Buffer));

    // unsqueeze permutation index - if this VCF has samples
    if (global_num_samples) {
        buf_alloc (vb, &vb->haplotype_permutation_index, vb->num_haplotypes_per_line * sizeof(unsigned), 0, 
                    "haplotype_permutation_index", vb->first_line);

        unsqueeze (vb,
                   (unsigned *)vb->haplotype_permutation_index.data, 
                   vardata_header->haplotype_index, 
                   ENDN16 (vardata_header->haplotype_index_checksum),
                   vb->num_haplotypes_per_line);
    }

  // get data for sample blocks - each block *may* have up to 3 file sections - genotype, phase and haplotype

    unsigned section_i=1;

    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned num_samples_in_sb = (sb_i == vb->num_sample_blocks-1 ? global_num_samples % vb->num_samples_per_block : vb->num_samples_per_block);

        // if genotype data exists, it appears first
        if (vb->has_genotype_data) 
            zfile_uncompress_section (vb, vb->z_data.data + section_index[section_i++], &vb->genotype_sections_data[sb_i], SEC_GENOTYPE_DATA);

        // next, comes phase data
        if (vb->phase_type == PHASE_MIXED_PHASED) {
            
            zfile_uncompress_section (vb, vb->z_data.data + section_index[section_i++], &vb->phase_sections_data[sb_i], SEC_PHASE_DATA);
            
            unsigned expected_size = vb->num_lines * num_samples_in_sb;
            ASSERT (vb->phase_sections_data[sb_i].len==expected_size, 
                    "Error: unexpected size of phase_sections_data[%u]: expecting %u but got %u", sb_i, expected_size, vb->phase_sections_data[sb_i].len)
        }

        // finally, comes haplotype data
        if (vb->has_haplotype_data) {
            
            zfile_uncompress_section (vb, vb->z_data.data + section_index[section_i++], &vb->haplotype_sections_data[sb_i], SEC_HAPLOTYPE_DATA);
            
            unsigned expected_size = vb->num_lines * num_samples_in_sb * vb->ploidy;
            ASSERT (vb->haplotype_sections_data[sb_i].len == expected_size, 
                    "Error: unexpected size of haplotype_sections_data[%u]: expecting %u but got %u", sb_i, expected_size, vb->haplotype_sections_data[sb_i].len)
        }
    }
}

// this is the compute thread entry point. It receives all data of a variant block and processes it
// in memory to the uncompressed format. This thread then terminates the I/O thread writes the output.
static void *piz_uncompress_variant_block (VariantBlock *vb)
{
    START_TIMER;

    piz_uncompress_all_sections (vb);

    // combine all the sections of a variant block to regenerate the variant_data, haplotype_data,
    // genotype_data and phase_data for each row of the variant block
    piz_reconstruct_line_components (vb);

    // merge line components (variant, haplotype, genotype, phase) back into a line
    //piz_merge_all_lines (vb);

    COPY_TIMER (vb->profile.piz_uncompress_variant_block);

#ifdef DEBUG
     buf_test_overflows(vb);
#endif
}

void *piz_uncompress_variant_block_thread_entry (void *vb_)
{
    piz_uncompress_variant_block ((VariantBlock*)vb_);
}

void piz_dispatcher (char *z_basename, File *z_file, File *vcf_file, bool concat_mode, bool test_mode, unsigned max_threads)
{
    START_WALLCLOCK;

    VariantBlockPool *vb_pool = vb_construct_pool (max_threads+1 /* +1 for pseudo-vb */, 100 /* for zip */);

    VariantBlock *pseudo_vb = vb_get_vb (vb_pool, vcf_file, z_file, 0);

    // read and write VCF header
    bool success = vcf_header_vcz_to_vcf ((VariantBlock*)pseudo_vb, concat_mode);
    if (!success) goto finish; // empty file - not an error

    typedef struct {
        pthread_t thread_id;
        VariantBlock *vb;
    } Thread;
    
    static Buffer compute_threads_buf = EMPTY_BUFFER;
    buf_alloc (pseudo_vb, &compute_threads_buf, sizeof(Thread) * (max_threads-1), 1, "compute_threads_buf", 0);
    Thread *compute_threads = (Thread *)compute_threads_buf.data;

    VariantBlock *writer_vb = NULL; // the variant block with completed data ready to be written to disk

    unsigned next_thread_to_dispatched   = 0;
    unsigned next_thread_to_be_joined    = 0;

    unsigned num_running_compute_threads = 0;
    unsigned input_exhausted = false;

    // this is the dispatcher loop. In each iteration, it can do one of 3 things, in this order of priority:
    // 1. In input is not exhausted, and a compute thread is available - read a variant block and compute it
    // 2. If data is available for writing, write it to disk
    // 3. Wait for the first thread (by sequential order) to complete the compute make data available for writing

    unsigned variant_block_i = 0;
    bool show_progress = z_basename && !!isatty(2);

#ifndef PROFILER
    if (show_progress) fprintf (stderr, "dvunzip %s: 0%%", z_basename);
#endif
    double last_pc_complete = 0.0;

    do {
        // PRIORITY 1: In input is not exhausted, and a compute thread is available - read a variant block and compute it
       if (!input_exhausted && num_running_compute_threads < max_threads - 1) {

            variant_block_i++;
            VariantBlock *vb = vb_get_vb (vb_pool, vcf_file, z_file, variant_block_i);

            success = zfile_read_one_vb (vb);
            if (!success) {
                input_exhausted = true;
                vb_release_vb (&vb);
                continue;
            }
            
            compute_threads[next_thread_to_dispatched].vb = vb;

            // note: the first variant block is always done on the main thread, so we can measure
            // memory consumption and reduce the number of compute threads and/or num_lines in subsequent vbs
            if (max_threads > 1 
#if __WIN32__
                && variant_block_i > 1
#endif
               ) {
                unsigned err = pthread_create(&compute_threads[next_thread_to_dispatched].thread_id, NULL, 
                                              piz_uncompress_variant_block_thread_entry, vb);
                ASSERT (!err, "Error: failed to create thread for processing variant block #%u in dv file, err=%u", vb->variant_block_i, err);
            }
            else { // single thread
                piz_uncompress_variant_block(vb);

                if (max_threads > 1)
                    vb_adjust_max_threads_by_vb_size(vb, &max_threads, test_mode);
            }

            num_running_compute_threads++;
            next_thread_to_dispatched = (next_thread_to_dispatched + 1) % (max_threads - 1);

            continue;
        }

        // PRIORITY 2: If data is available for writing, write it to disk
        if (writer_vb) {
            vcffile_write_one_variant_block (vcf_file, writer_vb);

            z_file->vcf_data_so_far += writer_vb->vcf_data_size; 

#ifdef DEBUG
            buf_test_overflows(writer_vb);
#endif
#ifdef PROFILER
            profiler_add (&pseudo_vb->profile, &writer_vb->profile);

            fprintf (stderr, "VUNBLOCK COMPLETED block #%u %s\n", 
                     writer_vb->variant_block_i, profiler_print_short (&writer_vb->profile));
#else
            if (show_progress) {
                int wallclock_ms=0;
                COPY_WALLCLOCK (wallclock_ms);

                vb_show_progress (&last_pc_complete, z_file, vcf_file->vcf_data_so_far,
                                  0, wallclock_ms/1000, 
                                  input_exhausted && !num_running_compute_threads, false);
            }
#endif
            vb_release_vb (&writer_vb); // cleanup vb and get it ready for another usage (without freeing memory)

            continue;
        }
        // PRIORITY 3: Wait for the first thread (by sequential order) to complete the compute make data available 
        // for writing
        if (num_running_compute_threads && !writer_vb) {

            // wait for thread to complete (possibly it completed already)
            Thread *th = &compute_threads[next_thread_to_be_joined];
            
            if (max_threads > 1) pthread_join(th->thread_id, NULL);

            writer_vb = th->vb;
            
            th->vb = NULL;
            th->thread_id = 0;
            num_running_compute_threads--;

            next_thread_to_be_joined = (next_thread_to_be_joined + 1) % (max_threads - 1);

            continue;
        }

    } while (!input_exhausted || num_running_compute_threads || writer_vb);

#ifdef PROFILER
    int wallclock_ms=0;
    COPY_WALLCLOCK (wallclock_ms);

    fprintf (stderr, "\nPROFILER:\nwallclock: %u\n%s\n", wallclock_ms, profiler_print_report (&pseudo_vb->profile));
#endif

finish:
    buf_free (&compute_threads_buf);
    vb_release_vb (&pseudo_vb);
}
