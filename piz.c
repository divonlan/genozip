// ------------------------------------------------------------------
//   piz.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"

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
    sprintf (pos_str, 
#ifdef _MSC_VER
    "%I64u", 
#else
    "%"PRIu64, 
#endif
    vb->last_pos);
    
    *delta_pos_len = s - str;
    *add_len = strlen (pos_str) - *delta_pos_len;

    COPY_TIMER(vb->profile.piz_decode_pos);
}

static void piz_get_line_get_num_subfields (VariantBlock *vb, unsigned line_i, // line in vcf file
                                            const char **line, unsigned *remaining_len,
                                            const char **subfields_start, unsigned *subfields_len, int *num_subfields)
{    
    START_TIMER;

    const char *after = *line + *remaining_len;

    unsigned column=1, i=0; for (; i < *remaining_len && column < 9; i++)
        if      ((*line)[i] == '\t') column++;
        else if ((*line)[i] == '\n') break;

    ASSERT (column==9, "Error: line %u is missing a FORMAT field", line_i); 

    *line += i; // past the tab
    
    if ((*line)[0] == '\n') { // this row has no genotype data, but maybe elsewhere in the VB there is
        *num_subfields = (unsigned)vb->has_genotype_data; // if we have genotype, then zip inserted BASE250_MISSING_SF 
        (*line)++; // past the newline
        goto cleanup;
    }
    if ((*line)[0] == 'G' && (*line)[1] == 'T') { // GT field in FORMAT columns - must always appear first per VCF spec (if at appears)
        *line += 3; // past the GT and : or \n
        if ((*line)[-1] == '\n') { // no subfields in this line
            *num_subfields = (unsigned)vb->has_genotype_data; // if we have genotype, then zip inserted BASE250_MISSING_SF 
            goto cleanup;
        }
    }

    *subfields_start = *line;

    for (*num_subfields = 1; *num_subfields <= MAX_SUBFIELDS; (*num_subfields)++) {
        seg_get_subfield (line, after-*line, line_i);
        if ((*line)[-1] == '\n') break;
    } 

    *subfields_len = (unsigned)(*line - *subfields_start); // length including separator

    ASSERT ((*line)[-1] == '\n', "Error: number of subfields declared in FORMAT exceeds maximum allowed of %u (excluding GT), line=%u", MAX_SUBFIELDS, vb->first_line+line_i);

cleanup:
    *remaining_len = after - *line;

    COPY_TIMER (vb->profile.piz_get_line_get_num_subfields)
}

static void piz_get_line_subfields (VariantBlock *vb, unsigned line_i, // line in vcf file
                                    const char *subfields_start, unsigned subfields_len,
                                    int *line_subfields) // out
{
    START_TIMER;

    for (unsigned i=0; i < MAX_SUBFIELDS; i++) line_subfields[i] = NIL;

    // case: this line has no subfields, despite other lines in the VB having
    if (!subfields_len) {
        line_subfields[0] = -2;
        return;
    }

    const char *after = subfields_start + subfields_len;
    for (unsigned i=0; i < MAX_SUBFIELDS; i++) {
        SubfieldIdType subfield = seg_get_subfield (&subfields_start, after-subfields_start, line_i);

        // the dictionaries were already read, so all subfields are expected to have a ctx
        unsigned sf_i=0 ; for (; sf_i < vb->z_file->num_subfields; sf_i++) 
            if (!memcmp (subfield.id, vb->z_file->mtf_ctx[sf_i].subfield.id, SUBFIELD_ID_LEN)) {
                // entry i corresponds to subfield i in FORMAT (excluding GT), and contains the index in mtf_ctx of this subfield
                line_subfields[i] = sf_i;
                break;
            }
#ifdef DEBUG
        ASSERTW (sf_i < vb->z_file->num_subfields, 
                 "Warning: subfield %.*s not found in dictionaries, line=%u. This can happen legitimately if the subfield is declared in FORMAT, but a value is never provided in any sample", SUBFIELD_ID_LEN, subfield.id, line_i);
#endif
        if (subfields_start[-1] == '\t' || subfields_start[-1] == '\n') break;
    } 

    COPY_TIMER (vb->profile.piz_get_line_subfields)
}

static void piz_get_variant_data_line (VariantBlock *vb, 
                                       unsigned line_i, // line in vcf file
                                       unsigned *length_remaining, // for safety
                                       const char **line_start) // out
{
    START_TIMER;

    // decoding the delta-encouded value of the POS field
    const char *delta_pos_start=0; // start of delta-POS field
    unsigned delta_pos_len=0; // length of encoded delta-pos
    int add_len=0; // how much more (or less) characters do we need in the string after decoding POS?
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
        }

        else if (*next == '\n') {
            next++; // past the newline
            unsigned line_strlen = next - *line_start; // inc. the newline, not including \0

            buf_alloc (vb, &vb->line_variant_data, line_strlen + add_len + 1 /* +1 sprintf adds a \0 */, 
                       1.2, "line_variant_data", 0);

            const char *after_delta = delta_pos_start + delta_pos_len;

            sprintf (vb->line_variant_data.data, "%.*s%s%.*s",
                     (int)(delta_pos_start - *line_start), *line_start,    // substring until \t preceding delta
                     pos_str,                           // decoded pos string
                     (int)(next - after_delta), after_delta);  // substring starting \t after delta
            
            vb->line_variant_data.len = line_strlen + add_len;

            *line_start = next;
            *length_remaining = after - next;

            goto cleanup;
        }

        next++;
    }
    
    ABORT0 ("Error: corrupt genozip file - at end of variant_data buffer, and no newline was found");

cleanup:
    COPY_TIMER(vb->profile.piz_get_variant_data_line);
    return;
}

// number of bytes this base250 number consumes in the gt data
static inline unsigned base250_len (const uint8_t *data) { 
    if (*data < BASE250_2_NUMERALS) return 1; // 1 byte
    else                            return *data - BASE250_2_NUMERALS + 3;
}

void piz_get_genotype_sample_starts (VariantBlock *vb, int *num_subfields)
{
    START_TIMER;
    
    buf_alloc (vb, &vb->next_gt_in_sample, sizeof(uint8_t*) * global_num_samples, 1, "next_gt_in_sample", 0);
    const uint8_t **next_gt_in_sample = (const uint8_t **)vb->next_gt_in_sample.data; // an array of uint8_t * - each element pointing to the gt of the first line of a sample
    
    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned num_samples_in_sb = vb_num_samples_in_sb (vb, sb_i);
        
        const uint8_t *next = (const uint8_t *)vb->genotype_sections_data[sb_i].data;
        const uint8_t *after = next + vb->genotype_sections_data[sb_i].len;

        unsigned sample_after = sb_i * SAMPLES_PER_BLOCK + num_samples_in_sb;
        
        unsigned sample_i = sb_i * SAMPLES_PER_BLOCK; 
        for (;sample_i < sample_after && next < after; sample_i++) {
            
            next_gt_in_sample[sample_i] = next; // line=0 of each sample_i (column)
            
            // now skip all remaining genotypes in this column, arriving at the beginning of the next column
            // (gt data is stored transposed - i.e. column by column)
            for (unsigned line_i=0; line_i < vb->num_lines; line_i++)
                for (unsigned sf=0; sf < num_subfields[line_i]; sf++) 
                    next += base250_len (next);
        }

        // sanity checks to see we read the correct amount of genotypes
        ASSERT (sample_i == sample_after, "Error: expected to find %u genotypes in sb_i=%u of variant_block_i=%u, but found only %u",
                vb->num_lines * num_samples_in_sb, sb_i, vb->variant_block_i, vb->num_lines * (sample_i - sb_i * SAMPLES_PER_BLOCK));

        ASSERT (next == after, "Error: expected to find %u genotypes in sb_i=%u of variant_block_i=%u, but found more. ",
                vb->num_lines * num_samples_in_sb, sb_i, vb->variant_block_i);
    }

    COPY_TIMER (vb->profile.piz_get_genotype_sample_starts)
}

// convert genotype data from sample block format of indices in base-250 to line format
// of tab-separated genotypes
static void piz_get_genotype_data_line (VariantBlock *vb, unsigned line_i, int *line_subfields)
{
    START_TIMER;

    uint8_t **next_gt_in_sample = (uint8_t **)vb->next_gt_in_sample.data; // for convenience

    char *next = vb->line_gt_data.data;
    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned first_sample = sb_i*SAMPLES_PER_BLOCK;
        unsigned num_samples_in_sb = vb_num_samples_in_sb (vb, sb_i);

        for (unsigned sample_i=first_sample; 
             sample_i < first_sample + num_samples_in_sb; 
             sample_i++) {

            char *snip = NULL; // will be set to a pointer into a dictionary
            
            for (unsigned sf_i=0; sf_i < vb->num_subfields; sf_i++) {

                if (line_subfields[sf_i] != NIL) {  // this line has this subfield (according to its FORMAT field)

                    // add a colon before, if needed
                    if (snip) *(next++) = ':'; // this works for empty "" snip too

                    uint8_t *word_index_base250 = next_gt_in_sample[sample_i];
                    unsigned snip_len;
                    mtf_get_snip_by_word_index (vb, &vb->mtf_ctx[line_subfields[sf_i]], // note: line_subfields[sf_i] maybe -2 (set in piz_get_line_subfields()), and this is an invalid value. this is ok, bc in this case word_index_base250 will be a control character
                                                word_index_base250, &snip, &snip_len);

                    if (snip && snip_len) { // it can be a valid empty subfield if snip="" and snip_len=0
                        memcpy (next, snip, snip_len);
                        next += snip_len;
                    }

                    next_gt_in_sample[sample_i] += base250_len (word_index_base250);
                }
            }

            // if we ended with a : - remove it
            next -= (next[-1] == ':');

            // add sample terminator - \t
            *(next++) = '\t';

            // safety
            ASSERT (next <= vb->line_gt_data.data + vb->line_gt_data.size, 
                    "Error: line_gt_data buffer overflow. variant_block_i=%u line_i=%u sb_i=%u sample_i=%u",
                    vb->variant_block_i, line_i + vb->first_line, sb_i, sample_i);
        }
    }

    // change last terminator to a \n
    next[-1] = '\n';

    vb->line_gt_data.len = next - vb->line_gt_data.data;

    vb->data_lines[line_i].has_genotype_data = vb->line_gt_data.len > global_num_samples; // not all just \t

    COPY_TIMER(vb->profile.piz_get_genotype_data_line);
}

static void piz_get_phase_data_line (VariantBlock *vb, unsigned line_i)
{
    START_TIMER;

    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned num_samples_in_sb = vb_num_samples_in_sb (vb, sb_i);

        memcpy (&vb->line_phase_data.data[sb_i * vb->num_samples_per_block],
                &vb->phase_sections_data[sb_i].data[line_i * num_samples_in_sb], 
                num_samples_in_sb);
    }

    COPY_TIMER(vb->profile.piz_get_phase_data_line);
}

// for each haplotype column, retrieve its it address in the haplotype sections. Note that since the haplotype sections are
// transposed, each column will be a row, or a contiguous array, in the section data. This function returns an array
// of pointers, each pointer being a beginning of column data within the section array
static const char **piz_get_ht_columns_data (VariantBlock *vb)
{
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
    if (!column_of_zeros) column_of_zeros = (char *)calloc (VARIANTS_PER_BLOCK, 1);

    for (unsigned ht_i=vb->num_haplotypes_per_line; ht_i < vb->num_haplotypes_per_line + 7; ht_i++)
        ht_columns_data[ht_i] = column_of_zeros;

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
    unsigned *next = (unsigned *)vb->line_ht_data.data;
    for (unsigned ht_i=0; ht_i < vb->num_haplotypes_per_line; ht_i += 4) 
        *(next++) = ((unsigned)(unsigned char)ht_columns_data[ht_i    ][line_i]      ) |  // this is LITTLE ENDIAN order
                    ((unsigned)(unsigned char)ht_columns_data[ht_i + 1][line_i] << 8 ) |
                    ((unsigned)(unsigned char)ht_columns_data[ht_i + 2][line_i] << 16) |
                    ((unsigned)(unsigned char)ht_columns_data[ht_i + 3][line_i] << 24) ;  // no worries if num_haplotypes_per_line is not a multiple of 4 - we have extra columns of zero

    // check if this row has now haplotype data (no GT field) despite some other rows in the VB having data
    DataLine *dl = &vb->data_lines[line_i];
    dl->has_haplotype_data = vb->line_ht_data.data[0] != '-'; // either the entire line is '-' or there is no '-' in the line

    COPY_TIMER(vb->profile.piz_get_haplotype_data_line);
}

// merge line components (variant, haplotype, genotype, phase) back into a line
static void piz_merge_line(VariantBlock *vb, unsigned line_i)
{
    START_TIMER;

    DataLine *dl = &vb->data_lines[line_i]; 

    // calculate the line length & allocate it
    unsigned ht_digits_len  = dl->has_haplotype_data ? vb->num_haplotypes_per_line : 0; 
    unsigned var_data_len   = vb->line_variant_data.len;        // includes a \n separator
    unsigned phase_sepr_len = dl->has_haplotype_data ? global_num_samples * (vb->ploidy-1) : 0; // the phase separators (/ or |)
    unsigned gt_colon_len   = ((dl->has_genotype_data && dl->has_haplotype_data) ? global_num_samples : 0); // the colon separating haplotype data from genotype data 
    unsigned gt_data_len    = (dl->has_genotype_data ? vb->line_gt_data.len // includes accounting for separator (\t or \n) after each sample
                                                     : (global_num_samples ? global_num_samples // separators after haplotype data with no genotype data
                                                                           : 0));  //  only variant data, no samples
    // adjustments to lengths
    if (dl->has_haplotype_data) {
        for (unsigned i=0; i < vb->num_haplotypes_per_line; i++) {
            // adjust for 2-digit alleles represented by 'A'..'Z' in the haplotype data
            if ((unsigned char)vb->line_ht_data.data[i] >= '0'+10) 
                ht_digits_len++; // 2-digit haplotype (10 to 99)

            // adjust for samples with ploidy less than vb->ploidy 
            if ((unsigned char)vb->line_ht_data.data[i] == '*') { // '*' means "ploidy padding"
                ht_digits_len--;
                phase_sepr_len--;
            }
        }
    }

    // this buf_alloc, when just naively called with the size actually needed, is responsible for 63% of the memory
    // consumed, and 34% of the execution time on one of our large test files. the time is mostly due to reallocs by subsequent
    // VBs. By having a 1.1 growth factor, we avoid most of the reallocs and significantly bring down the overall
    // execution time of genounzip
    
    dl->line.len = var_data_len + ht_digits_len + phase_sepr_len + gt_colon_len + gt_data_len; 
    buf_alloc (vb, &dl->line, (dl->line.len + 2), 1.1, "dl->line", vb->first_line + line_i); // +1 for string terminator, +1 for temporary additonal phase in case of '*'

    char *next    = dl->line.data;
    char *next_gt = vb->line_gt_data.data;

    // add variant data - change the \n separator to \t if needed
    memcpy (next, vb->line_variant_data.data, vb->line_variant_data.len);

    if (dl->has_genotype_data || dl->has_haplotype_data) 
        next[vb->line_variant_data.len-1] = '\t';

    next += vb->line_variant_data.len;

    // add samples
    for (unsigned sample_i=0; sample_i < global_num_samples; sample_i++) {

        // add haplotype data - ploidy haplotypes per sample 
        if (dl->has_haplotype_data) {

            PhaseType phase = (vb->phase_type == PHASE_MIXED_PHASED ? (PhaseType)vb->line_phase_data.data[sample_i]
                                                                    : vb->phase_type);
            ASSERT (phase=='/' || phase=='|' || phase=='1' || phase=='*', "Error: invalid phase character '%c' line_i=%u sample_i=%u", phase, line_i, sample_i+1);

            for (unsigned p=0; p < vb->ploidy ; p++) {
                
                unsigned char ht = vb->line_ht_data.data[sample_i * vb->ploidy + p];
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
        if (dl->has_genotype_data) 
            for (gt_len=0; next_gt[gt_len] != '\t' && next_gt[gt_len] != '\n'; gt_len++);

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

    // sanity check (the actual can be smaller in a line with missing samples)
    ASSERT (next - dl->line.data <= dl->line.len, "Error: unexpected line size in line_i=%u: calculated=%u, actual=%u", 
            vb->first_line + line_i, dl->line.len, (unsigned)(next - dl->line.data));

    dl->line.len = next - dl->line.data; // update line len to actual, which will be smaller in case of missing samples

    COPY_TIMER (vb->profile.piz_merge_line);
}

// combine all the sections of a variant block to regenerate the variant_data, haplotype_data,
// genotype_data and phase_data for each row of the variant block
void piz_reconstruct_line_components (VariantBlock *vb)
{
    START_TIMER;

    // initialize phase data if needed
    if (vb->phase_type == PHASE_MIXED_PHASED) 
        buf_alloc (vb, &vb->line_phase_data, global_num_samples, 1, "line_phase_data", 0);

    // initialize haplotype stuff
    const char **ht_columns_data=NULL;
    if (vb->has_haplotype_data) {

        //  memory - realloc for exact size, add 7 because depermuting_loop works on a word (32/64 bit) boundary
        buf_alloc (vb, &vb->line_ht_data, vb->num_haplotypes_per_line + 7, 1, "line_ht_data", 0);

        ht_columns_data = piz_get_ht_columns_data (vb);
    }

    // traverse the variant data first, only processing the FORMAT field - populate
    // the subfield data needed by piz_get_genotype_sample_starts

    const char *variant_data_next_line = vb->variant_data_section_data.data;
    unsigned variant_data_length_remaining = vb->variant_data_section_data.len;
    const char *subfields_start[VARIANTS_PER_BLOCK]; // pointer within the FORMAT field
    unsigned subfields_len[VARIANTS_PER_BLOCK];      // length of the FORMAT field, excluding GT, including the separator
    int num_subfields[VARIANTS_PER_BLOCK];           // number of subfields excluding GT
        
    // initialize genotype stuff
    if (vb->has_genotype_data) {
        for (unsigned line_i=0; line_i < vb->num_lines; line_i++) 
            // get subfields info from the FORMAT field
            piz_get_line_get_num_subfields (vb, vb->first_line + line_i, 
                                            &variant_data_next_line, &variant_data_length_remaining,
                                            &subfields_start[line_i], &subfields_len[line_i], &num_subfields[line_i]);

        piz_get_genotype_sample_starts(vb, num_subfields);

        buf_alloc (vb, &vb->line_gt_data, vb->max_gt_line_len, 1, "line_gt_data", 0);
    }

    // initialize again - for piz_get_variant_data_line
    variant_data_next_line = vb->variant_data_section_data.data;
    variant_data_length_remaining = vb->variant_data_section_data.len;

    for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {

        // de-permute variant data into vb->line_variant_data
        piz_get_variant_data_line (vb, vb->first_line + line_i, &variant_data_length_remaining, &variant_data_next_line);

        // reset len for next line - no need to realloc as we have realloced what is needed already
        vb->line_ht_data.len = vb->line_gt_data.len = vb->line_phase_data.len = 0;

        // transform sample blocks (each block: n_lines x s_samples) into line components (each line: 1 line x ALL_samples)
        if (vb->has_genotype_data)  {
            int line_subfields[MAX_SUBFIELDS]; // entry i corresponds to subfield i in FORMAT (excluding GT), and contains the index in mtf_ctx of this subfield
            piz_get_line_subfields (vb, vb->first_line + line_i,
                                    subfields_start[line_i], subfields_len[line_i], line_subfields);

            piz_get_genotype_data_line (vb, line_i, line_subfields);
        }
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
    vb->first_line              = ENDN32 (vardata_header->first_line);
    vb->num_lines               = ENDN16 (vardata_header->num_lines);
    vb->phase_type              = (PhaseType)vardata_header->phase_type;
    vb->has_genotype_data       = vardata_header->has_genotype_data;
    vb->is_sorted_by_pos        = vardata_header->is_sorted_by_pos;
    vb->num_haplotypes_per_line = ENDN32 (vardata_header->num_haplotypes_per_line);
    vb->has_haplotype_data      = vb->num_haplotypes_per_line > 0;
    vb->num_sample_blocks       = ENDN32 (vardata_header->num_sample_blocks);
    vb->num_samples_per_block   = ENDN32 (vardata_header->num_samples_per_block);
    vb->ploidy                  = ENDN16 (vardata_header->ploidy);
    vb->num_subfields           = ENDN16 (vardata_header->num_subfields);
    // num_dictionary_sections is read in zfile_read_one_vb()
    vb->max_gt_line_len         = ENDN32 (vardata_header->max_gt_line_len);
    memcpy(vb->chrom, vardata_header->chrom, MAX_CHROM_LEN);
    vb->min_pos                 = ENDN64 (vardata_header->min_pos);
    vb->max_pos                 = ENDN64 (vardata_header->max_pos);
    vb->vcf_data_size           = ENDN32 (vardata_header->vcf_data_size);
    
    ASSERT (global_num_samples == ENDN32 (vardata_header->num_samples), "Error: Expecting variant block to have %u samples, but it has %u", global_num_samples, vardata_header->num_samples);

    // we allocate memory for the Buffer arrays only once the first time this VariantBlock
    // is used. Subsequent blocks reusing the memory will have the same number of samples (by VCF spec)
    // BUG: this won't work if we're doing mutiple unrelated VCF on the command line
    if (vb->has_genotype_data && !vb->genotype_sections_data) 
        vb->genotype_sections_data  = (Buffer *)calloc (vb->num_sample_blocks, sizeof (Buffer));

    if (vb->phase_type == PHASE_MIXED_PHASED && !vb->phase_sections_data) 
        vb->phase_sections_data     = (Buffer *)calloc (vb->num_sample_blocks, sizeof (Buffer));
    
    if (vb->num_haplotypes_per_line && !vb->haplotype_sections_data) 
        vb->haplotype_sections_data = (Buffer *)calloc (vb->num_sample_blocks, sizeof (Buffer));

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
static void piz_uncompress_variant_block (VariantBlock *vb)
{
    START_TIMER;

    piz_uncompress_all_sections (vb);

    // combine all the sections of a variant block to regenerate the variant_data, haplotype_data,
    // genotype_data and phase_data for each row of the variant block
    piz_reconstruct_line_components (vb);

    // merge line components (variant, haplotype, genotype, phase) back into a line
    //piz_merge_all_lines (vb);

    COPY_TIMER (vb->profile.compute);

#ifdef DEBUG
    buf_test_overflows(vb);
#endif
}

void piz_dispatcher (const char *z_basename, File *z_file, File *vcf_file, bool test_mode, unsigned max_threads)
{
    Dispatcher dispatcher = dispatcher_init (max_threads, POOL_ID_UNZIP, vcf_file, z_file, test_mode, true, z_basename);

    // read and write VCF header
    bool success = vcf_header_genozip_to_vcf (dispatcher_get_pseudo_vb (dispatcher));
    if (!success) goto finish; // empty file - not an error

    // this is the dispatcher loop. In each iteration, it can do one of 3 things, in this order of priority:
    // 1. In input is not exhausted, and a compute thread is available - read a variant block and compute it
    // 2. Wait for the first thread (by sequential order) to complete and write data

    do {
        // PRIORITY 1: In input is not exhausted, and a compute thread is available - read a variant block and compute it
        if (!dispatcher_is_input_exhausted (dispatcher) && dispatcher_has_free_thread (dispatcher)) {

            bool success = zfile_read_one_vb (dispatcher_generate_next_vb (dispatcher));
            if (success)
                dispatcher_compute (dispatcher, piz_uncompress_variant_block);
            else
                dispatcher_input_exhausted (dispatcher);
        }

        // PRIORITY 2: Wait for the first thread (by sequential order) to complete and write data
        else {
            VariantBlock *processed_vb = dispatcher_get_next_processed_vb (dispatcher); // NULL if there is no computing thread
    
            vcffile_write_one_variant_block (vcf_file, processed_vb);

            z_file->vcf_data_so_far += processed_vb->vcf_data_size; 

            dispatcher_finalize_one_vb (dispatcher, z_file, vcf_file->vcf_data_so_far, 0);
        }

    } while (!dispatcher_is_done (dispatcher));

finish:
    dispatcher_finish (dispatcher);
}
