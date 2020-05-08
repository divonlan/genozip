// ------------------------------------------------------------------
//   piz_vcf.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "zfile.h"
#include "txtfile.h"
#include "header.h"
#include "seg.h"
#include "vblock.h"
#include "base250.h"
#include "dispatcher.h"
#include "move_to_front.h"
#include "file.h"
#include "endianness.h"
#include "squeeze_vcf.h"
#include "piz.h"
#include "sections.h"
#include "random_access.h"
#include "regions.h"
#include "samples.h"
#include "gtshark_vcf.h"
#include "dict_id.h"
#include "strings.h"

#define DATA_LINE(i) ENT (PizDataLineVCF, vb->lines, i)

typedef struct {
    unsigned num_subfields;         // number of subfields this FORMAT has
    MtfContext *ctx[MAX_SUBFIELDS]; // pointer to the ctx of each format subfield
} FormatInfo;


// 1. populates vb->format_info_buf with info about each unique format type in this vb (FormatInfo structure)
// 2. for each line, populates vb->lines[line_i].format_mtf_i - an index into vb->format_info_buf
static void piz_vcf_get_format_info (VBlockVCF *vb)
{    
    START_TIMER;

    MtfContext *format_ctx = &vb->mtf_ctx[VCF_FORMAT];

    // get number of subfields for each FORMAT item in dictionary, by traversing the FORMAT dectionary mtf array

    buf_alloc (vb, &vb->format_info_buf, sizeof(FormatInfo) * format_ctx->word_list.len, 1.5, "format_info_buf", 0);
    ARRAY (FormatInfo, formats, vb->format_info_buf);

    ARRAY (const MtfWord, snip_list, format_ctx->word_list);

    for (unsigned format_i=0; format_i < format_ctx->word_list.len; format_i++)
    {
        const char *format_snip = ENT (const char, format_ctx->dict, snip_list[format_i].char_index);
        uint32_t format_snip_len = snip_list[format_i].snip_len;
        
        // count colons in FORMAT snip
        unsigned num_colons = 0;
        int colons[MAX_SUBFIELDS+1];
        colons[num_colons++] = -1;
        for (unsigned i=0; i < format_snip_len; i++)
            if (format_snip[i] == ':') colons[num_colons++] = i;
        colons[num_colons++] = format_snip_len;

        bool format_has_gt_subfield = format_snip_len >= 2 && format_snip[0] == 'G' && format_snip[1] == 'T' && 
                                      (format_snip_len == 2 || format_snip[2] == ':');

        formats[format_i].num_subfields = num_colons - 1 - format_has_gt_subfield; // if FORMAT has a GT subfield - don't count it
        
        for (unsigned sf_i=0; sf_i < formats[format_i].num_subfields; sf_i++) {
            
            // construct dict_id for this format subfield
            
            const char *start = &format_snip[colons[sf_i + format_has_gt_subfield] + 1];
            unsigned len = colons[sf_i + format_has_gt_subfield + 1] - colons[sf_i + format_has_gt_subfield] - 1;

            DictIdType dict_id = dict_id_make (start, len);

            // get the did_i of this subfield. note: did_i can be NIL if the subfield appeared in a FORMAT field
            // in this VB, but never had any value in any sample on any line in this VB
            uint8_t did_i = mtf_get_existing_did_i_by_dict_id (dict_id);
            formats[format_i].ctx[sf_i] = (did_i != DID_I_NONE) ? &vb->mtf_ctx[did_i] : NULL;
        }
    }

    // now, get the FORMAT type (format_mtf_i) in each line of the VB, by traversing the FORMAT b250 data
    for (unsigned line_i=0; line_i < vb->lines.len; line_i++) 
        DATA_LINE (line_i)->format_mtf_i = mtf_get_next_snip ((VBlockP)vb, format_ctx, NULL, NULL, NULL, vb->first_line + line_i);

    // reset format_ctx reader iterator fields, as we are going to traverse FORMAT again when reconstructing the lines
    mtf_init_iterator (format_ctx);

    COPY_TIMER (vb->profile.piz_vcf_get_format_info)
}

static bool piz_vcf_reconstruct_special_info_subfields (VBlock *vb, uint8_t did_i, DictIdType dict_id, uint32_t txt_line_i)
{
    // starting v5, END is a delta vs POS
    if (is_v5_or_above && dict_id.num == dict_id_INFO_END) {
        DECLARE_SNIP;
        RECONSTRUCT_FROM_DICT_POS (did_i, vb->last_pos, true, NULL, false);
        return false;
    }

    return true; // proceed with normal reconstruction
}
        
static bool piz_vcf_reconstruct_fields (VBlockVCF *vb, unsigned vb_line_i, 
                                        bool *has_13) // out: original vcf line ended with Windows-style \r\n
{
    uint32_t txt_line_i = vb->first_line + vb_line_i;
    bool line_included = true;
    uint64_t txt_data_start = vb->txt_data.len;

    DECLARE_SNIP;
    uint32_t chrom_word_index = RECONSTRUCT_FROM_DICT (VCF_CHROM, true);
    RECONSTRUCT_FROM_DICT_POS (VCF_POS, vb->last_pos, true, NULL, true); // reconstruct from delta

    // in case of --regions - check if this line is needed at all (based on CHROM and POS)
    if (flag_regions && !regions_is_site_included (chrom_word_index, vb->last_pos)) // test here bc vb->last_pos might change if INFO has END
        line_included = false;

    IFNOTSTRIP(".",1) { RECONSTRUCT_ID (VCF_ID, &vb->id_numeric_data, &vb->next_id_numeric, NULL, true); }
    RECONSTRUCT_FROM_DICT (VCF_REFALT, true);
    IFNOTSTRIP(".",1) { RECONSTRUCT_FROM_DICT (VCF_QUAL, true);}
    IFNOTSTRIP(".",1) { RECONSTRUCT_FROM_DICT (VCF_FILTER, true);}

    IFNOTSTRIP(".",1) { 
        uint32_t iname_word_index = LOAD_SNIP (VCF_INFO);        
        piz_reconstruct_info ((VBlockP)vb, iname_word_index, snip, snip_len, piz_vcf_reconstruct_special_info_subfields, 
                              txt_line_i, has_13);
        RECONSTRUCT1 ("\t");
    }

    if (flag_gt_only || flag_strip) {
        if (global_vcf_num_samples) RECONSTRUCT ("GT\t", 3);
    }
    else
        RECONSTRUCT_FROM_DICT (VCF_FORMAT, true);

    // after consuming sections' data, if this line is not to be outputed - shorten txt_data back to start of line
    if (!line_included) vb->txt_data.len = txt_data_start; 

    return line_included;
}

// initialize vb->sample_iterator to the first line in the gt data for each sample (column) 
static void piz_vcf_initialize_sample_iterators (VBlockVCF *vb)
{
    START_TIMER;
    
    buf_alloc (vb, &vb->sample_iterator, sizeof(SnipIterator) * global_vcf_num_samples, 1, "sample_iterator", 0);
    
    ARRAY (SnipIterator, sample_iterator, vb->sample_iterator);
    ARRAY (FormatInfo, formats, vb->format_info_buf);

    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        if (! *ENT(bool, vb->is_sb_included, sb_i)) continue; // skip if this sample block is excluded by --samples

        unsigned num_samples_in_sb = vb_vcf_num_samples_in_sb (vb, sb_i);
        
        const uint8_t *next  = FIRSTENT (const uint8_t, vb->genotype_sections_data[sb_i]);
        const uint8_t *after = AFTERENT (const uint8_t, vb->genotype_sections_data[sb_i]);

        unsigned sample_after = sb_i * vb->num_samples_per_block + num_samples_in_sb;
        
        unsigned sample_i = sb_i * vb->num_samples_per_block; 
        for (;sample_i < sample_after && next < after; sample_i++) {
            
            sample_iterator[sample_i].next_b250 = next; // line=0 of each sample_i (column)
            sample_iterator[sample_i].prev_word_index = 1;

            // now skip all remaining genotypes in this column, arriving at the beginning of the next column
            // (gt data is stored transposed - i.e. column by column)
            for (uint32_t line_i=0; line_i < (uint32_t)vb->lines.len; line_i++) {
                
                FormatInfo *line_format_info = &formats[DATA_LINE (line_i)->format_mtf_i];
                uint32_t num_subfields = line_format_info->num_subfields;
                
                for (unsigned sf=0; sf < num_subfields; sf++) 
                    next += base250_len (next); // if this format has no non-GT subfields, it will not have a ctx 
            }
        }

        // sanity checks to see we read the correct amount of genotypes
        ASSERT (sample_i == sample_after, "Error: expected to find %u genotypes in sb_i=%u of vblock_i=%u, but found only %u",
                (uint32_t)vb->lines.len * num_samples_in_sb, sb_i, vb->vblock_i, (uint32_t)vb->lines.len * (sample_i - sb_i * vb->num_samples_per_block));

        ASSERT (next == after, "Error: unused data remains in buffer after processing genotype data for sb_i=%u of vblock_i=%u (%u lines x %u samples)",
                sb_i, vb->vblock_i, (uint32_t)vb->lines.len, num_samples_in_sb);
    }

    COPY_TIMER (vb->profile.piz_vcf_initialize_sample_iterators)
}

// convert genotype data from sample block format of indices in base-250 to line format
// of tab-separated genotype data string, each string being a colon-seperated list of subfields, 
// the subfields being defined in the FORMAT of this line
static void piz_vcf_reconstruct_genotype_data_line (VBlockVCF *vb, unsigned vb_line_i, bool is_line_included)
{
    START_TIMER;

    PizDataLineVCF *dl = DATA_LINE (vb_line_i);

    ARRAY (SnipIterator, sample_iterator, vb->sample_iterator);
    ARRAY (const FormatInfo, formats, vb->format_info_buf);

    const FormatInfo *line_format_info = &formats[dl->format_mtf_i];

    char *next = vb->line_gt_data.data;
    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned first_sample = sb_i * vb->num_samples_per_block;
        unsigned num_samples_in_sb = vb_vcf_num_samples_in_sb (vb, sb_i);

        for (unsigned sample_i=first_sample; 
             sample_i < first_sample + num_samples_in_sb; 
             sample_i++) {

            if (!samples_am_i_included (sample_i)) continue;

            const char *snip = NULL; // will be set to a pointer into a dictionary
            
            uint32_t dp_value=0;  // used for calculating MIN_DP = dp_value + delta.

            for (unsigned sf_i=0; sf_i < line_format_info->num_subfields; sf_i++) {

                MtfContext *sf_ctx = line_format_info->ctx[sf_i];

                ASSERT (sf_ctx || *sample_iterator[sample_i].next_b250 == BASE250_MISSING_SF, 
                        "Error: line_format_info->ctx[sf_i=%u] for line %u sample %u (both counting from 1) is NULL, indicating that this subfield has no value in the vb in any sample or any line. And yet, it does...", 
                        sf_i, vb_line_i + vb->first_line, sample_i+1);

                // add a colon before, if needed
                if (snip && is_line_included) *(next++) = ':'; // this works for empty "" snip too

                unsigned snip_len;
                mtf_get_next_snip ((VBlockP)vb, sf_ctx, &sample_iterator[sample_i], &snip, &snip_len, vb->first_line + vb_line_i);

                // handle MIN_DP : if its a DP, store it...
                char min_dp[30];
                if (sf_ctx && snip_len && sf_ctx->dict_id.num == dict_id_FORMAT_DP) 
                    dp_value = atoi (snip);
                // ...and if its an MIN_DP - calculate it 
                else if (sf_ctx && sf_ctx->dict_id.num == dict_id_FORMAT_MIN_DP && snip_len && is_v5_or_above) {
                    int32_t delta = atoi (snip);
                    str_int (dp_value - delta, min_dp, &snip_len); // note: we dp_value==0 if no DP subfield preceeds DP_MIN, that's fine
                    snip = min_dp;
                }

                if (snip && snip_len && is_line_included) { // it can be a valid empty subfield if snip="" and snip_len=0
                    memcpy (next, snip, snip_len); 
                    next += snip_len;
                }
            }

            if (is_line_included) {
                // if we ended with a : - remove it
                next -= (next[-1] == ':');

                // add sample terminator - \t
                *(next++) = '\t';

                // safety
                ASSERT (next <= vb->line_gt_data.data + vb->line_gt_data.size, 
                        "Error: line_gt_data buffer overflow. vblock_i=%u vcf_line_i=%u sb_i=%u sample_i=%u",
                        vb->vblock_i, vb_line_i + vb->first_line, sb_i, sample_i);
            }
        } // for sample
    } // for sample block
    
    // change last terminator to a \n
    if (is_line_included) next[-1] = '\n';

    vb->line_gt_data.len = next - vb->line_gt_data.data;

    dl->has_genotype_data = (vb->line_gt_data.len > global_vcf_num_displayed_samples); // not all just \t

    COPY_TIMER(vb->profile.piz_vcf_reconstruct_genotype_data_line);
}

static void piz_vcf_get_phase_data_line (VBlockVCF *vb, unsigned vb_line_i)
{
    START_TIMER;

    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        if (! *ENT(bool, vb->is_sb_included, sb_i)) continue; // skip if this sample block is excluded by --samples

        unsigned num_samples_in_sb = vb_vcf_num_samples_in_sb (vb, sb_i);

        memcpy (ENT (char, vb->line_phase_data, sb_i * vb->num_samples_per_block),
                ENT (char, vb->phase_sections_data[sb_i], vb_line_i * num_samples_in_sb), 
                num_samples_in_sb);
    }

    COPY_TIMER(vb->profile.piz_vcf_get_phase_data_line);
}

// for each haplotype column, retrieve its it address in the haplotype sections. Note that since the haplotype sections are
// transposed, each column will be a row, or a contiguous array, in the section data. This function returns an array
// of pointers, each pointer being a beginning of column data within the section array
static const char **piz_vcf_get_ht_columns_data (VBlockVCF *vb)
{
    buf_alloc (vb, &vb->ht_columns_data, sizeof (char *) * (vb->num_haplotypes_per_line + 7), 1, "ht_columns_data", 0); // realloc for exact size (+15 is padding for 64b operations)

    // each entry is a pointer to the beginning of haplotype column located in vb->haplotype_sections_data
    // note: in genozip v1 haplotype columns were permuted across the entire samples, as of v2 they are permuted
    // only within their own sample block
    ARRAY (const char *, ht_columns_data, vb->ht_columns_data); 
    ARRAY (const unsigned, permutatation_index, vb->haplotype_permutation_index);
    
    const unsigned max_ht_per_block = vb->num_samples_per_block * vb->ploidy; // last sample block may have less, but that's ok for our div/mod calculations below

    // provide 7 extra zero-columns for the convenience of the permuting loop (supporting 64bit assignments)
    // note: txt_file->max_lines_per_vb will be zero if genozip file was created by redirecting output
    buf_alloc (vb, &vb->column_of_zeros, MAX (txt_file->max_lines_per_vb, vb->lines.len), 1, "column_of_zeros", 0);
    buf_zero (&vb->column_of_zeros);

    for (uint32_t sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        bool is_sb_included = *ENT(bool, vb->is_sb_included, sb_i);

        for (unsigned ht_i = sb_i * max_ht_per_block; 
             ht_i < MIN ((sb_i+1) * max_ht_per_block, vb->num_haplotypes_per_line); 
             ht_i++) {

            if (!is_sb_included) { // for sample blocks excluded by --samples, just have columns of 0 for the convenience of the permuting loop
                ht_columns_data[ht_i] = vb->column_of_zeros.data;
                continue;
            }

            unsigned permuted_ht_i = permutatation_index[ht_i];
            //unsigned sb_i    = permuted_ht_i / max_ht_per_block; // get haplotype sample block per this ht is


            unsigned row     = permuted_ht_i % max_ht_per_block; // get row transposed haplotype sample block. column=line_i
            unsigned sb_ht_i = row * (uint32_t)vb->lines.len;              // index within haplotype block 

            ht_columns_data[ht_i] = &vb->haplotype_sections_data[sb_i].data[sb_ht_i];

            ASSERT (ht_columns_data[ht_i], 
                    "Error in piz_vcf_get_ht_columns_data: haplotype column is NULL for vb_i=%u sb_i=%u sb_ht_i=%u", vb->vblock_i, sb_i, sb_ht_i);
        }
    }

    for (unsigned ht_i=vb->num_haplotypes_per_line; ht_i < vb->num_haplotypes_per_line + 7; ht_i++)
        ht_columns_data[ht_i] = vb->column_of_zeros.data;

    return ht_columns_data;
}

// build haplotype for a line - reversing the permutation and the transposal.
static void piz_vcf_get_haplotype_data_line (VBlockVCF *vb, unsigned vb_line_i, const char **ht_columns_data)
{
    START_TIMER;

    const uint32_t max_ht_per_block = vb->num_samples_per_block * vb->ploidy; // last sample block may have less, but that's ok for our div/mod calculations below

    if (flag_samples) memset (vb->line_ht_data.data, 0, vb->num_haplotypes_per_line); // if we're not filling in all samples, initialize to 0;

    uint32_t ht_i = 0;
    for (uint32_t sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        if (! *ENT(bool, vb->is_sb_included, sb_i)) continue;

        // start from the nearest block 8 columns that includes our start column (might include some previous columns too)
        // but if already done (because it overlaps the previous SB that that SB was included) 
        // start from the next one that is not done. last block of 8 columns might overlap the next vb
        ht_i = MAX (ht_i, (sb_i * max_ht_per_block) & 0xfffffff8);

        uint32_t ht_i_after = MIN ((sb_i+1) * max_ht_per_block, vb->num_haplotypes_per_line);

        uint64_t *next = (uint64_t *)&vb->line_ht_data.data[ht_i];

        // this loop can consume up to 25-50% of the entire decompress compute time (tested with 1KGP data)
        // note: we do memory assignment 64 bit at time (its about 10% faster than byte-by-byte)
        
        for (; ht_i < ht_i_after; ht_i += 8) {

#ifdef __LITTLE_ENDIAN__
            *(next++) = ((uint64_t)(uint8_t)ht_columns_data[ht_i    ][vb_line_i]      ) |  // this is LITTLE ENDIAN order
                        ((uint64_t)(uint8_t)ht_columns_data[ht_i + 1][vb_line_i] << 8 ) |
                        ((uint64_t)(uint8_t)ht_columns_data[ht_i + 2][vb_line_i] << 16) |
                        ((uint64_t)(uint8_t)ht_columns_data[ht_i + 3][vb_line_i] << 24) |
                        ((uint64_t)(uint8_t)ht_columns_data[ht_i + 4][vb_line_i] << 32) |
                        ((uint64_t)(uint8_t)ht_columns_data[ht_i + 5][vb_line_i] << 40) |
                        ((uint64_t)(uint8_t)ht_columns_data[ht_i + 6][vb_line_i] << 48) |
                        ((uint64_t)(uint8_t)ht_columns_data[ht_i + 7][vb_line_i] << 56) ;  // no worries if num_haplotypes_per_line is not a multiple of 4 - we have extra columns of zero
#else
            *(next++) = ((uint64_t)(uint8_t)ht_columns_data[ht_i    ][vb_line_i] << 56) |  // this is BIG ENDIAN order
                        ((uint64_t)(uint8_t)ht_columns_data[ht_i + 1][vb_line_i] << 48) |
                        ((uint64_t)(uint8_t)ht_columns_data[ht_i + 2][vb_line_i] << 40) |
                        ((uint64_t)(uint8_t)ht_columns_data[ht_i + 3][vb_line_i] << 32) |
                        ((uint64_t)(uint8_t)ht_columns_data[ht_i + 4][vb_line_i] << 24) |
                        ((uint64_t)(uint8_t)ht_columns_data[ht_i + 5][vb_line_i] << 16) |
                        ((uint64_t)(uint8_t)ht_columns_data[ht_i + 6][vb_line_i] << 8 ) |
                        ((uint64_t)(uint8_t)ht_columns_data[ht_i + 7][vb_line_i]      ) ;  
#endif
        }
    }

    // check if this row has no haplotype data (no GT field) despite some other rows in the VB having data
    PizDataLineVCF *dl = DATA_LINE (vb_line_i);
    
    // check to see if this line has any sample with haplotype info (when not using --samples, this loop
    // usually terminates in the first iteration - so not adding a lot of overhead)
    dl->has_haplotype_data = false;
    for (uint32_t ht_i=0; ht_i < vb->num_haplotypes_per_line; ht_i++) {
        // it can be '-' in three scenarios. 
        // 1. '-' - line has no GT despite other lines in the VB having - seg_vcf_complete_missing_lines sets it all to '-'
        // 2. 0   - initialized above ^ and not filled in due to --samples 
        // 3. 0   - it comes from a column of zeros set in piz_vcf_get_ht_columns_data - due to --samples
        if (vb->line_ht_data.data[ht_i] != '-' && vb->line_ht_data.data[ht_i] != 0) { 
            dl->has_haplotype_data = true; // found one sample that has haplotype
            break;
        }
    }
    COPY_TIMER(vb->profile.piz_vcf_get_haplotype_data_line);
}

// add the samples (haplotype, genotype, phase) to txt_dataa
static void piz_vcf_reconstruct_samples (VBlockVCF *vb, unsigned vb_line_i, bool has_13)
{
    START_TIMER;

    PizDataLineVCF *dl = DATA_LINE (vb_line_i);
    uint32_t line_start = vb->txt_data.len;
    char *next_gt = vb->line_gt_data.data;

    if (dl->has_genotype_data || dl->has_haplotype_data) 
        *LASTENT (char, vb->txt_data) = '\t';

    // add samples
    for (unsigned sample_i=0; sample_i < global_vcf_num_samples; sample_i++) {

        if (!samples_am_i_included (sample_i)) continue;

        // add haplotype data - ploidy haplotypes per sample 
        if (dl->has_haplotype_data) {

            PhaseType phase = (vb->phase_type == PHASE_MIXED_PHASED ? (PhaseType)vb->line_phase_data.data[sample_i]
                                                                    : vb->phase_type);
            ASSERT (phase=='/' || phase=='|' || phase=='1' || phase=='*', 
                    "Error: invalid phase character '%c' line_i=%u sample_i=%u", phase, vb_line_i + vb->first_line, sample_i+1);

            for (unsigned p=0; p < vb->ploidy ; p++) {
                
                uint8_t ht = vb->line_ht_data.data[sample_i * vb->ploidy + p];
                if (ht == '.') // unknown
                    NEXTENT (char, vb->txt_data) = ht;

                else if (ht == '*')  // missing haplotype - delete previous phase
                    vb->txt_data.len--;

                else { // allele 0 to 99
                    unsigned allele = ht - '0'; // allele 0->99 represented by ascii 48->147
                    ASSERT (allele <= VCF_MAX_ALLELE_VALUE, "Error: allele out of range: %u (ht=ascii(%u)) line_i=%u sample_i=%u", allele, (unsigned)ht, vb->first_line + vb_line_i, sample_i+1);
                    
                    if (allele >= 10) NEXTENT (char, vb->txt_data) = '0' + allele / 10;
                    NEXTENT (char, vb->txt_data) = '0' + allele % 10;
                }

                // add the phase character between the haplotypes
                if (vb->ploidy >= 2 && p < vb->ploidy-1) 
                    NEXTENT (char, vb->txt_data) = phase;
            }
        }

        // get length of genotype data - we need it now, to see if we need a :
        unsigned gt_len;
        if (dl->has_genotype_data) 
            for (gt_len=0; next_gt[gt_len] != '\t' && next_gt[gt_len] != '\n'; gt_len++);

        // add colon separating haplotype from genotype data
        if (dl->has_haplotype_data && dl->has_genotype_data && gt_len) 
            NEXTENT (char, vb->txt_data) = ':';

        // add genotype data - dropping the \t
        if (dl->has_genotype_data) {
            RECONSTRUCT (next_gt, gt_len);
            next_gt += gt_len + 1;
        }

        // add tab separator after sample data
        if (dl->has_haplotype_data || dl->has_genotype_data) 
            NEXTENT (char, vb->txt_data) = '\t';
    }

    // trim trailing tabs due to missing data
    for (; vb->txt_data.len >= line_start+2 && *ENT(char, vb->txt_data, vb->txt_data.len - 2) == '\t'; vb->txt_data.len--); // after this loop, next points to the first tab after the last non-tab character

    // now, we need to replace the last tab with \n or with \r\n (depending on has_13)
    vb->txt_data.len--;
    RECONSTRUCT (has_13 ? "\r\n" : "\n", has_13 ? 2 : 1);

    COPY_TIMER (vb->profile.piz_vcf_reconstruct_samples);
}

// combine all the sections of a variant block to regenerate the variant_data, haplotype_data,
// genotype_data and phase_data for each row of the variant block
static void piz_vcf_reconstruct_vb (VBlockVCF *vb)
{
    START_TIMER;

    buf_alloc (vb, &vb->lines, vb->lines.len * sizeof (PizDataLineVCF), 1.2, "lines", vb->vblock_i);
    buf_zero (&vb->lines);

    // initialize phase data if needed
    if (vb->phase_type == PHASE_MIXED_PHASED && !flag_drop_genotypes) 
        buf_alloc (vb, &vb->line_phase_data, global_vcf_num_samples, 1, "line_phase_data", vb->vblock_i);

    // initialize haplotype stuff
    const char **ht_columns_data=NULL;
    if (vb->has_haplotype_data && !flag_drop_genotypes) {

        //  memory - realloc for exact size, add 7 because depermuting_loop works on a word (32/64 bit) boundary
        buf_alloc (vb, &vb->line_ht_data, vb->num_haplotypes_per_line + 7, 1, "line_ht_data", vb->vblock_i);

        ht_columns_data = piz_vcf_get_ht_columns_data (vb);
    }
    
    // initialize genotype stuff
    if (vb->has_genotype_data && !flag_drop_genotypes && !flag_strip && !flag_gt_only) {
        
        // get info about the different types of FORMAT in this vb (vb->format_info_buf)
        // as well as which is used for each line (dl->format_mtf_i)
        piz_vcf_get_format_info (vb);

        // initialize vb->sample_iterator to the first line in the gt data for each sample (column) 
        piz_vcf_initialize_sample_iterators(vb);

        buf_alloc (vb, &vb->line_gt_data, vb->max_gt_line_len, 1, "line_gt_data", vb->vblock_i);
    }

    // these arrays (for fields) and iname_mapper->next (for info subfields) contain pointers to the next b250 item.
    // every line, in the for loop, MAY progress the pointer by 1, if that b250 was used for that row (all are used for the 
    // fields, but only those info subfields defined in the INFO names of a particular line are used in that line).
            
    // now reconstruct the lines, one line at a time

    buf_alloc (vb, &vb->txt_data, (vb->vb_data_size + 1), 1.1, "txt_data", vb->vblock_i); // +1 bc piz_vcf_reconstruct_samples needs room for a temporary redundant '*' in case ht='*'

    for (unsigned vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        // re-construct variant data (fields CHROM to FORMAT, including INFO subfields) into vb->txt_data
        bool has_13 = false;
        bool is_line_included = piz_vcf_reconstruct_fields (vb, vb_line_i, &has_13);

        // transform sample blocks (each block: n_lines x s_samples) into line components (each line: 1 line x ALL_samples)
        if (!flag_drop_genotypes) {
            // note: we always call piz_vcf_reconstruct_genotype_data_line even if !is_line_included, bc we need to advance the iterators
            if (vb->has_genotype_data && !flag_strip && !flag_gt_only)  
                piz_vcf_reconstruct_genotype_data_line (vb, vb_line_i, is_line_included);

            if (is_line_included)  {
                if (vb->phase_type == PHASE_MIXED_PHASED) 
                    piz_vcf_get_phase_data_line (vb, vb_line_i);

                if (vb->has_haplotype_data) 
                    piz_vcf_get_haplotype_data_line (vb, vb_line_i, ht_columns_data);

                piz_vcf_reconstruct_samples (vb, vb_line_i, has_13);
            }
        }
            
        // reset len for next line - no need to alloc as all the lines are the same size?
        vb->line_ht_data.len = vb->line_gt_data.len = vb->line_phase_data.len = 0;
    }

    COPY_TIMER(vb->profile.piz_reconstruct_vb);
}

static void piz_vcf_uncompress_all_sections (VBlockVCF *vb)
{
    // The VB is read from disk in zfile_vcf_read_one_vb(), in the I/O thread, and is decompressed here in the 
    // Compute thread, with the exception of dictionaries that are processed by the I/O thread
    // Order of sections in a V2 VB:
    // 1. SEC_VB_HEADER - its data is the haplotype index
    // 2. (the dictionaries were here in the file, but they are omitted from vb->z_data)
    // 3. b250 of CHROM to FORMAT fields
    // 4. SEC_VCF_INFO_SF_B250 - All INFO subfield data
    // 5. All sample data: multiple sections per sample block:
    //    5a. SEC_VCF_GT_DATA - genotype data
    //    5b. SEC_VCF_PHASE_DATA - phase data
    //    5c. SEC_HT_DATA (1 section) or SEC_HAPLOTYPE_GTSHARK (5 sections) - haplotype data

    ARRAY (const unsigned, section_index, vb->z_section_headers);

    SectionHeaderVbHeaderVCF *header = (SectionHeaderVbHeaderVCF *)(vb->z_data.data + section_index[0]);
    vb->first_line              = BGEN32 (header->first_line);
    vb->lines.len               = BGEN32 (header->num_lines);
    vb->phase_type              = (PhaseType)header->phase_type;
    vb->has_genotype_data       = header->has_genotype_data;
    vb->num_haplotypes_per_line = BGEN32 (header->num_haplotypes_per_line);
    vb->has_haplotype_data      = vb->num_haplotypes_per_line > 0;
    vb->num_sample_blocks       = BGEN32 (header->num_sample_blocks);
    vb->num_samples_per_block   = BGEN32 (header->num_samples_per_block);
    vb->ploidy                  = BGEN16 (header->ploidy);
    vb->max_gt_line_len         = BGEN32 (header->max_gt_line_len);
    vb->vb_data_size            = BGEN32 (header->vb_data_size);

    // this can if 1. VCF has no samples or 2. num_samples was not re-written to genozip header (for example if we were writing to stdout)
    if (!global_vcf_num_samples) 
        global_vcf_num_samples = BGEN32 (header->num_samples);
    else {
        ASSERT (global_vcf_num_samples == BGEN32 (header->num_samples), "Error: Expecting variant block to have %u samples, but it has %u", global_vcf_num_samples, BGEN32 (header->num_samples));
    }

    // if the user filtered out all samples, we don't need to even uncompress them
    if (flag_samples && !global_vcf_num_displayed_samples) {
        vb->has_genotype_data = false;
        vb->has_haplotype_data = false;
    }

    // in case of --split, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every vcf component
    if (flag_split) 
        vb->vblock_i = BGEN32 (header->h.vblock_i);
    
    // unsqueeze permutation index - if this VCF has samples, AND this vb has any haplotype data
    if (global_vcf_num_samples && vb->num_haplotypes_per_line && !flag_drop_genotypes) {

       zfile_uncompress_section ((VBlockP)vb, &vb->z_data.data[section_index[0]], &vb->haplotype_permutation_index_squeezed, 
                                 "haplotype_permutation_index_squeezed", SEC_VB_HEADER);

        buf_alloc (vb, &vb->haplotype_permutation_index, vb->num_haplotypes_per_line * sizeof(uint32_t), 0, 
                    "haplotype_permutation_index", vb->first_line);

        unsqueeze (vb,
                   (unsigned *)vb->haplotype_permutation_index.data, 
                   (uint8_t *)vb->haplotype_permutation_index_squeezed.data, 
                   BGEN16 (header->haplotype_index_checksum),
                   vb->num_haplotypes_per_line);
    }

    unsigned section_i=1;

    // uncompress the 8 fields (CHROM to FORMAT)    
    piz_uncompress_fields ((VBlockP)vb, section_index, &section_i);

    for (uint8_t sf_i=0; sf_i < vb->num_info_subfields ; sf_i++) {
        
        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[section_i++]);

        if (flag_strip) continue; // we don't need to create INFO subfield ctx if --strip

        MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, vb->dict_id_to_did_i_map, &vb->num_dict_ids, &vb->num_info_subfields, 
                                                  header->dict_id, SEC_VCF_INFO_SF_DICT);

        zfile_uncompress_section ((VBlockP)vb, header, &ctx->b250, "mtf_ctx.b250", SEC_VCF_INFO_SF_B250);    
    }

    if (is_v5_or_above) {
        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[section_i++]);
        
        zfile_uncompress_section ((VBlockP)vb, header, &vb->id_numeric_data, "id_numeric_data", SEC_NUMERIC_ID_DATA);
    }

    if (flag_drop_genotypes) return; // if --drop-genotypes was requested - no need to uncompress the following sections

    // we allocate memory for the Buffer arrays only once the first time this VBlockVCF
    // is used. Subsequent blocks reusing the memory will have the same number of samples (by VCF spec)
    // BUG: this won't work if we're doing mutiple unrelated VCF on the command line
    if (vb->has_genotype_data && !vb->genotype_sections_data) 
        vb->genotype_sections_data  = (Buffer *)calloc (vb->num_sample_blocks, sizeof (Buffer));

    if (vb->phase_type == PHASE_MIXED_PHASED && !vb->phase_sections_data) 
        vb->phase_sections_data     = (Buffer *)calloc (vb->num_sample_blocks, sizeof (Buffer));
    
    if (vb->num_haplotypes_per_line && !vb->haplotype_sections_data) 
        vb->haplotype_sections_data = (Buffer *)calloc (vb->num_sample_blocks, sizeof (Buffer));

    // get data for sample blocks - each block *may* have up to 3 file sections - genotype, phase and haplotype

    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        // skip uncompressing this sample block if it is excluded by --samples
        if (! *ENT(bool, vb->is_sb_included, sb_i)) {
            section_i += (vb->has_genotype_data ? 1 : 0) + 
                         (vb->phase_type == PHASE_MIXED_PHASED) + 
                         (vb->has_haplotype_data ? (header->is_gtshark ? 5 : 1) : 0); // just advance section_i
            continue;
        }

        unsigned num_samples_in_sb = vb_vcf_num_samples_in_sb (vb, sb_i);

        // if genotype data exists, it appears first
        if (vb->has_genotype_data) {
            if (!flag_strip && !flag_gt_only) 
                zfile_uncompress_section ((VBlockP)vb, vb->z_data.data + section_index[section_i], &vb->genotype_sections_data[sb_i], "genotype_sections_data", SEC_VCF_GT_DATA);
            section_i++;
        }                

        // next, comes phase data
        if (vb->phase_type == PHASE_MIXED_PHASED) {
            
            zfile_uncompress_section ((VBlockP)vb, vb->z_data.data + section_index[section_i++], &vb->phase_sections_data[sb_i], "phase_sections_data", SEC_VCF_PHASE_DATA);
            
            unsigned expected_size = vb->lines.len * num_samples_in_sb;
            ASSERT (vb->phase_sections_data[sb_i].len == expected_size, 
                    "Error: unexpected size of phase_sections_data[%u]: expecting %u but got %u", 
                    sb_i, expected_size, (uint32_t)vb->phase_sections_data[sb_i].len)
        }

        // finally, comes haplotype data
        if (vb->has_haplotype_data) {
            
            if (!header->is_gtshark) {
                zfile_uncompress_section ((VBlockP)vb, vb->z_data.data + section_index[section_i++], &vb->haplotype_sections_data[sb_i], 
                                          "haplotype_sections_data", SEC_HT_DATA );
            }
            else { // gtshark

                zfile_uncompress_section ((VBlockP)vb, vb->z_data.data + section_index[section_i++], &vb->gtshark_exceptions_line_i, 
                                          "gtshark_exceptions_line_i", SEC_HT_GTSHARK_X_LINE);
                vb->gtshark_exceptions_line_i.len /= sizeof (uint32_t);

                zfile_uncompress_section ((VBlockP)vb, vb->z_data.data + section_index[section_i++], &vb->gtshark_exceptions_ht_i, 
                                          "gtshark_exceptions_ht_i", SEC_HT_GTSHARK_X_HTI);
                vb->gtshark_exceptions_ht_i.len /= sizeof (uint16_t);

                zfile_uncompress_section ((VBlockP)vb, vb->z_data.data + section_index[section_i++], &vb->gtshark_exceptions_allele, 
                                          "gtshark_exceptions_allele", SEC_HT_GTSHARK_X_ALLELE);

                zfile_uncompress_section ((VBlockP)vb, vb->z_data.data + section_index[section_i++], &vb->gtshark_db_db_data, 
                                          "gtshark_db_db_data", SEC_HT_GTSHARK_DB_DB);

                zfile_uncompress_section ((VBlockP)vb, vb->z_data.data + section_index[section_i++], &vb->gtshark_db_gt_data, 
                                          "gtshark_db_gt_data", SEC_HT_GTSHARK_DB_GT);

                gtshark_uncompress_haplotype_data (vb, sb_i);
            }

            unsigned expected_size = (uint32_t)vb->lines.len * num_samples_in_sb * vb->ploidy;
            ASSERT (vb->haplotype_sections_data[sb_i].len == expected_size, 
                    "Error: unexpected size of haplotype_sections_data[%u]: expecting %u but got %u", 
                    sb_i, expected_size, (uint32_t)vb->haplotype_sections_data[sb_i].len)
        }
    }
}

// this is the compute thread entry point. It receives all data of a variant block and processes it
// in memory to the uncompressed format. This thread then terminates the I/O thread writes the output.
void piz_vcf_uncompress_one_vb (VBlock *vb_)
{
    START_TIMER;

    VBlockVCF *vb = (VBlockVCF *)vb_;

    if (is_v2_or_above) {
        piz_vcf_uncompress_all_sections (vb);

        // combine all the sections of a variant block to regenerate the variant_data, haplotype_data,
        // genotype_data and phase_data for each row of the variant block
        piz_vcf_reconstruct_vb (vb);
    }

    // v1 c
    else {
        void v1_piz_vcf_uncompress_all_sections (VBlockVCFP vb); // forwwrd declaration - these are included at the end of this file
        v1_piz_vcf_uncompress_all_sections (vb);

        void v1_piz_vcf_reconstruct_vb (VBlockVCFP vb);
        v1_piz_vcf_reconstruct_vb (vb);
    }

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway

    COPY_TIMER (vb->profile.compute);
}

#define V1_PIZ // select the piz functions of v1.c
#include "v1_vcf.c"
