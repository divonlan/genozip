// ------------------------------------------------------------------
//   piz.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "zfile.h"
#include "vcffile.h"
#include "vcf_header.h"
#include "segregate.h"
#include "vb.h"
#include "base250.h"
#include "dispatcher.h"
#include "move_to_front.h"
#include "file.h"
#include "endianness.h"
#include "squeeze.h"
#include "piz.h"
#include "sections.h"
#include "random_access.h"

typedef struct {
    unsigned num_subfields;         // number of subfields this FORMAT has
    MtfContext *ctx[MAX_SUBFIELDS]; // pointer to the ctx of each format subfield
} FormatInfo;

// decode the delta-encoded value of the POS field
static inline void piz_decode_pos (VariantBlock *vb, const char *delta_snip, unsigned delta_snip_len,
                                   char *pos_str, unsigned *pos_len) // out
{
    START_TIMER;

    int32_t delta=0;

    // we parse the string ourselves - this is hugely faster than sscanf.
    unsigned negative = *delta_snip == '-'; // 1 or 0

    const char *s; for (s=(delta_snip + negative); *s != '\t' ; s++)
        delta = delta*10 + (*s - '0');

    if (negative) delta = -delta;

    vb->last_pos += delta;
    
    ASSERT (vb->last_pos >= 0, "Error: vb->last_pos=%d is negative", vb->last_pos);

    // create number string without calling slow sprintf

    // create reverse string
    char reverse_pos_str[50];
    uint32_t n = vb->last_pos;

    unsigned len=0; 
    if (n) {
        while (n) {
            reverse_pos_str[len++] = '0' + (n % 10);
            n /= 10;
        }

        // reverse it
        for (unsigned i=0; i < len; i++) pos_str[i] = reverse_pos_str[len-i-1];
        pos_str[len] = '\0';

        *pos_len = len;
    }
    else {  // n=0. note: POS=0 is not compliant with the VCF spcification
        pos_str[0] = 0;
        *pos_len = 1;
    }
    COPY_TIMER(vb->profile.piz_decode_pos);
}


// 1. populates vb->format_info_buf with info about each unique format type in this vb (FormatInfo structure)
// 2. for each line, populates vb->data_lines[line_i].format_mtf_i - an index into vb->format_info_buf
static void piz_get_format_info (VariantBlock *vb)
{    
    START_TIMER;

    MtfContext *format_ctx = &vb->mtf_ctx[FORMAT];

    // get number of subfields for each FORMAT item in dictionary, by traversing the FORMAT dectionary mtf array

    buf_alloc (vb, &vb->format_info_buf, sizeof(FormatInfo) * format_ctx->word_list.len, 1.5, "format_info_buf", 0);
    FormatInfo *format_num_subfields = (FormatInfo *)vb->format_info_buf.data;

    const MtfWord *snip_list = (const MtfWord *)format_ctx->word_list.data;

    for (unsigned format_i=0; format_i < format_ctx->word_list.len; format_i++)
    {
        const char *format_snip = &format_ctx->dict.data[snip_list[format_i].char_index];
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

        format_num_subfields[format_i].num_subfields = num_colons - 1 - format_has_gt_subfield; // if FORMAT has a GT subfield - don't count it
        
        for (unsigned sf_i=0; sf_i < format_num_subfields[format_i].num_subfields; sf_i++) {
            
            // construct dict_id for this format subfield
            
            const char *start = &format_snip[colons[sf_i + format_has_gt_subfield] + 1];
            unsigned len = colons[sf_i + format_has_gt_subfield + 1] - colons[sf_i + format_has_gt_subfield] - 1;

            DictIdType dict_id = dict_id_make (start, len);

            // get the did_i of this subfield. note: did_i can be NIL if the subfield appeared in a FORMAT field
            // in this VB, but never had any value in any sample on any line in this VB
            int did_i = mtf_get_existing_did_i_by_dict_id (vb, dict_id);
            format_num_subfields[format_i].ctx[sf_i] = (did_i != NIL) ? &vb->mtf_ctx[did_i] : NULL;
        }
    }

    // now, get the FORMAT type (format_mtf_i) in each line of the VB, by traversing the FORMAT b250 data
    for (unsigned line_i=0; line_i < vb->num_lines; line_i++) 
        vb->data_lines[line_i].format_mtf_i = mtf_get_next_snip (vb, format_ctx, NULL, NULL, NULL, vb->first_line + line_i);

    // reset format_ctx reader iterator fields, as we are going to traverse FORMAT again when reconstructing the lines
    mtf_init_iterator (format_ctx);

    COPY_TIMER (vb->profile.piz_get_format_info)
}

// for each unique type of INFO fields (each one containing multiple names), create a unique mapping
// info field node index (i.e. in b250) -> list of names, lengths and the context of the subfields
static void piz_map_iname_subfields (VariantBlock *vb)
{
    // terminology: we call a list of INFO subfield names, an "iname". An iname looks something like
    // this: "I1=I2=I3=". Each iname consists of info subfields. These fields are not unique to this
    // iname and can appear in other inames. The INFO field contains the iname, and values of the subfields.
    // iname _mapper maps these subfields. This function creates an iname_mapper for every unique iname.

    const MtfContext *info_ctx = &vb->mtf_ctx[INFO];
    vb->iname_mapper_buf.len = info_ctx->word_list.len;
    buf_alloc (vb, &vb->iname_mapper_buf, sizeof (SubfieldInfoMapperPiz) * vb->iname_mapper_buf.len,
               1, "iname_mapper_buf", 0);
    buf_zero (&vb->iname_mapper_buf);

    SubfieldInfoMapperPiz *all_iname_mappers = (SubfieldInfoMapperPiz*)vb->iname_mapper_buf.data;

    const MtfWord *all_inames = (const MtfWord *)info_ctx->word_list.data;

    for (unsigned iname_i=0; iname_i < vb->iname_mapper_buf.len; iname_i++) {

        const char *iname = (const char *)&info_ctx->dict.data[all_inames[iname_i].char_index]; // e.g. "I1=I2=I3=" - pointer into the INFO dictionary
        unsigned iname_len = all_inames[iname_i].snip_len; 
        SubfieldInfoMapperPiz *iname_mapper = &all_iname_mappers[iname_i]; // iname_mapper of this specific set of names "I1=I2=I3="

        // get INFO subfield snips - which are the values of the INFO subfield, where the names are
        // in the INFO snip in the format "info1=info2=info3="). 
        DictIdType dict_id;
        iname_mapper->num_subfields = 0;

        // traverse the subfields of one iname. E.g. if the iname is "I1=I2=I3=" then we traverse I1, I2, I3
        for (unsigned i=0; i < iname_len; i++) {
            
            iname_mapper->names[iname_mapper->num_subfields] = &iname[i];

            // traverse the iname, and get the dict_id for each subfield name (using only the first 8 characers)
            dict_id.num = 0;
            unsigned j=0; 
            while (iname[i] != '=' && iname[i] != '\t') { // value-less INFO names can be terminated by the end-of-word \t in the dictionary
                if (j < DICT_ID_LEN) dict_id.id[j] = iname[i]; // scan the whole name, but copy only the first 8 bytes to dict_id
                i++, j++;
            }
            dict_id = dict_id_info_subfield (dict_id);

            iname_mapper->name_lens[iname_mapper->num_subfields] = j + (iname[i] == '='); // including the '='

            int did_i = mtf_get_existing_did_i_by_dict_id (vb, dict_id); // it will be NIL if this is an INFO name without values            
            if (did_i != NIL) {
                iname_mapper->ctx[iname_mapper->num_subfields] = &vb->mtf_ctx[did_i];

                ASSERT (iname_mapper->ctx[iname_mapper->num_subfields]->dict_id.num == dict_id.num, "Error: unexpected dict_id. iname_mapper->ctx->dict_id=%.*s dict_id=%.*s", DICT_ID_LEN, 
                        iname_mapper->ctx[iname_mapper->num_subfields]->dict_id.id, DICT_ID_LEN, dict_id.id);
            }

            iname_mapper->num_subfields++;
        }
    }
}

static void piz_get_variant_data_line (VariantBlock *vb, unsigned vb_line_i)
{
    START_TIMER;

    const char *snip[NUM_VCF_B250S]; // snip (pointer into dictionary) and snip_len of each field in this line
    uint32_t snip_len[NUM_VCF_B250S];
    memset (snip_len, 0, sizeof(snip_len));
    memset (snip, 0, sizeof(snip));

    const char *info_sf_value_snip[MAX_SUBFIELDS]; // snip (pointer into dictionary) and snip_len of each field in this line
    uint32_t info_sf_value_snip_len[MAX_SUBFIELDS];

    // get mtf_i and variant data length
    unsigned line_len = 0;
    char pos_str[50];
    SubfieldInfoMapperPiz *iname_mapper = NULL;

    // extract snips and calculate length of variant data
    for (VcfFields f=CHROM; f <= FORMAT; f++) {

        // if the VB doesn't have FORMAT at all - skip it
        if (f==FORMAT && !vb->has_genotype_data && !vb->has_haplotype_data && vb->mtf_ctx[f].dict_section_type != SEC_FORMAT_DICT) continue;

        uint32_t index = mtf_get_next_snip (vb, &vb->mtf_ctx[f], NULL, &snip[f], &snip_len[f], vb->first_line + vb_line_i);

        // reconstruct pos from delta
        if (f == POS) {
            piz_decode_pos (vb, snip[POS], snip_len[POS], pos_str, &snip_len[POS]); 
            snip[POS] = pos_str;
        }

        // add the INFO subfield values
        else if (f == INFO) {
            ASSERT (index >= 0 && index < vb->iname_mapper_buf.len, 
                    "Error: iname_mapper index out of range: index=%d, vb->iname_mapper_buf.len=%u", index, vb->iname_mapper_buf.len);

            iname_mapper = &((SubfieldInfoMapperPiz *)vb->iname_mapper_buf.data)[index];
            for (unsigned sf_i = 0; sf_i < iname_mapper->num_subfields; sf_i++) {
                                
                if (!iname_mapper->ctx[sf_i]) continue; // a name without values

                mtf_get_next_snip (vb, iname_mapper->ctx[sf_i], NULL, &info_sf_value_snip[sf_i], &info_sf_value_snip_len[sf_i], vb->first_line + vb_line_i);

                line_len += info_sf_value_snip_len[sf_i];
            }

            // add the ; between name=value pairs in the INFO data (e.g. "name1=info1;name2=info2")
            line_len += iname_mapper->num_subfields - 1;
        }

        // add the field (for INFO - the names, values are added already ^ )
        line_len += snip_len[f] + 1; // add \t or \n separator
    }
    
    buf_alloc (vb, &vb->line_variant_data, line_len, 1.5, "line_variant_data", 0);

    // reconstrut the line
    for (VcfFields f=CHROM; f <= FORMAT; f++) {

        // info subfield eg "info1=value1;info2=value2" - "info1=", "info2=" are the name snips
        // while "value1" and "value2" are the value snips - we merge them here
        if (f == INFO) {
            for (unsigned sf_i=0; sf_i < iname_mapper->num_subfields ; sf_i++) {

                buf_add (&vb->line_variant_data, iname_mapper->names[sf_i], iname_mapper->name_lens[sf_i]); // name inc. '=' e.g. "Info1="
                
                if (iname_mapper->ctx[sf_i])  // some info names can be without values, in which case there will be no ctx
                    buf_add (&vb->line_variant_data, info_sf_value_snip[sf_i], info_sf_value_snip_len[sf_i]); // value e.g "value1"
    
                if (sf_i != iname_mapper->num_subfields-1)
                    buf_add (&vb->line_variant_data, ";", 1); // seperator between each two name=value pairs e.g "name1=value;name2=value2"
            }
        }

        // other, non-INFO fields
        else if (snip_len[f])  // FORMAT can have snip_len=0, in which case its the end of the line and no \t either
            buf_add (&vb->line_variant_data, snip[f], snip_len[f]);

        if (f != INFO || snip_len[FORMAT]) // add a tab after the field EXCEPT for an INFO before an empty FORMAT
            buf_add (&vb->line_variant_data, (f == FORMAT ? "\n" : "\t"), 1); // \n at end of line, \t between other fields
    }

    COPY_TIMER(vb->profile.piz_get_variant_data_line);
}

// initialize vb->sample_iterator to the first line in the gt data for each sample (column) 
static void piz_initialize_sample_iterators (VariantBlock *vb)
{
    START_TIMER;
    
    buf_alloc (vb, &vb->sample_iterator, sizeof(SnipIterator) * global_num_samples, 1, "sample_iterator", 0);
    SnipIterator *sample_iterator = (SnipIterator *)vb->sample_iterator.data; // an array of SnipIterator
    
    FormatInfo *format_num_subfields = (FormatInfo *)vb->format_info_buf.data;

    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned num_samples_in_sb = vb_num_samples_in_sb (vb, sb_i);
        
        const uint8_t *next = (const uint8_t *)vb->genotype_sections_data[sb_i].data;
        const uint8_t *after = next + vb->genotype_sections_data[sb_i].len;

        unsigned sample_after = sb_i * vb->num_samples_per_block + num_samples_in_sb;
        
        unsigned sample_i = sb_i * vb->num_samples_per_block; 
        for (;sample_i < sample_after && next < after; sample_i++) {
            
            sample_iterator[sample_i].next_b250 = next; // line=0 of each sample_i (column)
            sample_iterator[sample_i].prev_word_index = 1;

            // now skip all remaining genotypes in this column, arriving at the beginning of the next column
            // (gt data is stored transposed - i.e. column by column)
            for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {
                
                FormatInfo *line_format_info = &format_num_subfields[vb->data_lines[line_i].format_mtf_i];
                uint32_t num_subfields = line_format_info->num_subfields;
                
                for (unsigned sf=0; sf < num_subfields; sf++) 
                    next += base250_len (next, line_format_info->ctx[sf] ? line_format_info->ctx[sf]->encoding : B250_ENC_NONE); // if this format has no non-GT subfields, it will not have a ctx 
            }
        }

        // sanity checks to see we read the correct amount of genotypes
        ASSERT (sample_i == sample_after, "Error: expected to find %u genotypes in sb_i=%u of variant_block_i=%u, but found only %u",
                vb->num_lines * num_samples_in_sb, sb_i, vb->variant_block_i, vb->num_lines * (sample_i - sb_i * vb->num_samples_per_block));

        ASSERT (next == after, "Error: unused data remains in buffer after processing genotype data for sb_i=%u of variant_block_i=%u (%u lines x %u samples)",
                sb_i, vb->variant_block_i, vb->num_lines, num_samples_in_sb);
    }

    COPY_TIMER (vb->profile.piz_initialize_sample_iterators)
}

// convert genotype data from sample block format of indices in base-250 to line format
// of tab-separated genotype data string, each string being a colon-seperated list of subfields, 
// the subfields being defined in the FORMAT of this line
static void piz_get_genotype_data_line (VariantBlock *vb, unsigned vb_line_i)
{
    START_TIMER;

    DataLine *dl = &vb->data_lines[vb_line_i];

    SnipIterator *sample_iterator = (SnipIterator *)vb->sample_iterator.data; // for convenience

    const FormatInfo *format_num_subfields = (const FormatInfo *)vb->format_info_buf.data;
    const FormatInfo *line_format_info = &format_num_subfields[dl->format_mtf_i];

    char *next = vb->line_gt_data.data;
    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned first_sample = sb_i * vb->num_samples_per_block;
        unsigned num_samples_in_sb = vb_num_samples_in_sb (vb, sb_i);

        for (unsigned sample_i=first_sample; 
             sample_i < first_sample + num_samples_in_sb; 
             sample_i++) {

            const char *snip = NULL; // will be set to a pointer into a dictionary
            
            for (unsigned sf_i=0; sf_i < line_format_info->num_subfields; sf_i++) {

                MtfContext *sf_ctx = line_format_info->ctx[sf_i];

                ASSERT (sf_ctx || *sample_iterator[sample_i].next_b250 == BASE250_MISSING_SF, 
                        "Error: line_format_info->ctx[sf_i=%u] for line %u sample %u (both counting from 1) is NULL, indicating that this subfield has no value in the vb in any sample or any line. And yet, it does...", 
                        sf_i, vb_line_i + vb->first_line, sample_i+1);

                // add a colon before, if needed
                if (snip) *(next++) = ':'; // this works for empty "" snip too

                unsigned snip_len;
                mtf_get_next_snip (vb, sf_ctx, &sample_iterator[sample_i], &snip, &snip_len, vb->first_line + vb_line_i);

                if (snip && snip_len) { // it can be a valid empty subfield if snip="" and snip_len=0
                    memcpy (next, snip, snip_len); 
                    next += snip_len;
                }
            }

            // if we ended with a : - remove it
            next -= (next[-1] == ':');

            // add sample terminator - \t
            *(next++) = '\t';

            // safety
            ASSERT (next <= vb->line_gt_data.data + vb->line_gt_data.size, 
                    "Error: line_gt_data buffer overflow. variant_block_i=%u line_i=%u sb_i=%u sample_i=%u",
                    vb->variant_block_i, vb_line_i + vb->first_line, sb_i, sample_i);
        } // for sample
    } // for sample block
    
    // change last terminator to a \n
    next[-1] = '\n';

    vb->line_gt_data.len = next - vb->line_gt_data.data;

    dl->has_genotype_data = (vb->line_gt_data.len > global_num_samples); // not all just \t

    COPY_TIMER(vb->profile.piz_get_genotype_data_line);
}

static void piz_get_phase_data_line (VariantBlock *vb, unsigned vb_line_i)
{
    START_TIMER;

    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned num_samples_in_sb = vb_num_samples_in_sb (vb, sb_i);

        memcpy (&vb->line_phase_data.data[sb_i * vb->num_samples_per_block],
                &vb->phase_sections_data[sb_i].data[vb_line_i * num_samples_in_sb], 
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
        unsigned sb_ht_i = row * vb->num_lines;              // index within haplotype block 

        ht_columns_data[ht_i] = &vb->haplotype_sections_data[sb_i].data[sb_ht_i];
    }

    // provide 7 extra zero-columns for the convenience of the permuting loop (supporting 32bit and 64bit assignments)
    if (!buf_is_allocated (&vb->column_of_zeros)) {
        // a constant array of zero, preserved between VBs and concatenated files - as global_max_lines_per_vb is not changed
        buf_alloc (vb, &vb->column_of_zeros, global_max_lines_per_vb, 1, "column_of_zeros", 0);
        buf_zero (&vb->column_of_zeros);
    } 

    for (unsigned ht_i=vb->num_haplotypes_per_line; ht_i < vb->num_haplotypes_per_line + 7; ht_i++)
        ht_columns_data[ht_i] = vb->column_of_zeros.data;

    return ht_columns_data;
}

// build haplotype for a line - reversing the permutation and the transposal.
static void piz_get_haplotype_data_line (VariantBlock *vb, unsigned vb_line_i, const char **ht_columns_data)
{
    START_TIMER;

    // this loop can consume up to 50% of the entire decompress compute time (tested with 1KGP data)
    // note: we do memory assignment 32 bit at time (its about 10% faster than byte-by-byte)
    // TO DO: this is LITTLE ENDIAN, need an #ifdef to support BIG ENDIAN
    // TO DO: we can also #ifdef between 32bit and 64bit compilers and do 8 at a time for the latter
    unsigned *next = (unsigned *)vb->line_ht_data.data;
    for (unsigned ht_i=0; ht_i < vb->num_haplotypes_per_line; ht_i += 4) 
        *(next++) = ((unsigned)(unsigned char)ht_columns_data[ht_i    ][vb_line_i]      ) |  // this is LITTLE ENDIAN order
                    ((unsigned)(unsigned char)ht_columns_data[ht_i + 1][vb_line_i] << 8 ) |
                    ((unsigned)(unsigned char)ht_columns_data[ht_i + 2][vb_line_i] << 16) |
                    ((unsigned)(unsigned char)ht_columns_data[ht_i + 3][vb_line_i] << 24) ;  // no worries if num_haplotypes_per_line is not a multiple of 4 - we have extra columns of zero

    // check if this row has now haplotype data (no GT field) despite some other rows in the VB having data
    DataLine *dl = &vb->data_lines[vb_line_i];
    dl->has_haplotype_data = vb->line_ht_data.data[0] != '-'; // either the entire line is '-' or there is no '-' in the line

    COPY_TIMER(vb->profile.piz_get_haplotype_data_line);
}

// merge line components (variant, haplotype, genotype, phase) back into a line
static void piz_merge_line(VariantBlock *vb, unsigned vb_line_i)
{
    START_TIMER;

    DataLine *dl = &vb->data_lines[vb_line_i]; 

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
    buf_alloc (vb, &dl->line, (dl->line.len + 2), 1.1, "dl->line", vb->first_line + vb_line_i); // +1 for string terminator, +1 for temporary additonal phase in case of '*'

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
            ASSERT (phase=='/' || phase=='|' || phase=='1' || phase=='*', "Error: invalid phase character '%c' line_i=%u sample_i=%u", phase, vb_line_i + vb->first_line, sample_i+1);

            for (unsigned p=0; p < vb->ploidy ; p++) {
                
                unsigned char ht = vb->line_ht_data.data[sample_i * vb->ploidy + p];
                if (ht == '.') // unknown 
                    *(next++) = ht;

                else if (ht == '*')  // missing haplotype - delete previous phase
                    *(--next) = 0;

                else { // allele 0 to 99
                    unsigned allele = ht - '0'; // allele 0->99 represented by ascii 48->147
                    ASSERT (allele <= 99, "Error: allele out of range: %u line_i=%u sample_i=%u", allele, vb->first_line + vb_line_i, sample_i+1);
                    
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

        // add tab separator after sample data
        if (dl->has_haplotype_data || dl->has_genotype_data) *(next++) = '\t';
    }

    // trim trailing tabs due to missing data
    for (; next >= dl->line.data+2 && next[-2] == '\t'; next--); // after this loop, next points to the first tab after the last non-tab character

    next[-1] = '\n'; // replace the last tab with a newline
    next[0]  = '\0'; // end of string
    
    // sanity check (the actual can be smaller in a line with missing samples)
    ASSERT (next - dl->line.data <= dl->line.len, "Error: unexpected line size in line_i=%u: calculated=%u, actual=%u", 
            vb->first_line + vb_line_i, dl->line.len, (unsigned)(next - dl->line.data));

    dl->line.len = next - dl->line.data; // update line len to actual, which will be smaller in case of missing samples

    COPY_TIMER (vb->profile.piz_merge_line);
}

// combine all the sections of a variant block to regenerate the variant_data, haplotype_data,
// genotype_data and phase_data for each row of the variant block
static void piz_reconstruct_line_components (VariantBlock *vb)
{
    START_TIMER;

    if (!vb->data_lines) 
        vb->data_lines = calloc (global_max_lines_per_vb, sizeof (DataLine));

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
    
    // initialize genotype stuff
    if (vb->has_genotype_data) {
        
        // get info about the different types of FORMAT in this vb (vb->format_info_buf)
        // as well as which is used for each line (dl->format_mtf_i)
        piz_get_format_info (vb);

        // initialize vb->sample_iterator to the first line in the gt data for each sample (column) 
        piz_initialize_sample_iterators(vb);

        buf_alloc (vb, &vb->line_gt_data, vb->max_gt_line_len, 1, "line_gt_data", 0);
    }

    // this arrays (for fields) and iname_mapper->next (for info subfields)  contain pointers to the next b250 item.
    // every line, in the for loop, MAY progress the pointer by 1, if that b250 was used for that row (all are used for the 
    // fields, but only those info subfields defined in the INFO names of a particular line are used in that line).
            
    // create mapping for info subfields
    piz_map_iname_subfields (vb);

    // now reconstruct the lines, one line at a time
    for (unsigned vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {

        // re-construct variant data (fields CHROM to FORMAT, including INFO subfields) into vb->line_variant_data
        piz_get_variant_data_line (vb, vb_line_i);

        // transform sample blocks (each block: n_lines x s_samples) into line components (each line: 1 line x ALL_samples)
        if (vb->has_genotype_data)  
            piz_get_genotype_data_line (vb, vb_line_i);

        if (vb->phase_type == PHASE_MIXED_PHASED) 
            piz_get_phase_data_line (vb, vb_line_i);

        if (vb->has_haplotype_data) 
            piz_get_haplotype_data_line (vb, vb_line_i, ht_columns_data);

        piz_merge_line (vb, vb_line_i);
        
        // reset len for next line - no need to alloc as all the lines are the same size?
        vb->line_ht_data.len = vb->line_gt_data.len = vb->line_phase_data.len = 0;
        buf_free (&vb->line_variant_data);
    }

    COPY_TIMER(vb->profile.piz_reconstruct_line_components);
}

static void piz_uncompress_all_sections (VariantBlock *vb)
{
    // The VB is read from disk in zfile_read_one_vb(), in the I/O thread, and is decompressed here in the 
    // Compute thread, with the exception of dictionaries that are processed by the I/O thread
    // Order of sections in a V2 VB:
    // 1. SEC_VB_HEADER - its data is the haplotype index
    // 2. (the dictionaries were here in the file orecn disk, but they are omitted from vb->z_data)
    // 3. SEC_INFO_SUBFIELD_B250 - All INFO subfield data
    // 4. All sample data: up 3 sections per sample block:
    //    4a. SEC_GENOTYPE_DATA - genotype data
    //    4b. SEC_PHASE_DATA - phase data
    //    4c. SEC_HAPLOTYPE_DATA - haplotype data

    unsigned *section_index = (unsigned *)vb->z_section_headers.data;

    SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)(vb->z_data.data + section_index[0]);
    vb->first_line              = BGEN32 (header->first_line);
    vb->num_lines               = BGEN32 (header->num_lines);
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
    if (!global_num_samples) 
        global_num_samples = BGEN32 (header->num_samples);
    else {
        ASSERT (global_num_samples == BGEN32 (header->num_samples), "Error: Expecting variant block to have %u samples, but it has %u", global_num_samples, BGEN32 (header->num_samples));
    }

    // in case of --split, the variant_block_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every vcf component
    if (flag_split) 
        vb->variant_block_i = BGEN32 (header->h.variant_block_i);
    
    // unsqueeze permutation index - if this VCF has samples, AND this vb has any haplotype data
    if (global_num_samples && vb->num_haplotypes_per_line) {

       zfile_uncompress_section (vb, &vb->z_data.data[section_index[0]], &vb->haplotype_permutation_index_squeezed, 
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
    for (VcfFields f=CHROM; f <= FORMAT; f++) {

        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[section_i++]);

        zfile_uncompress_section (vb, header, &vb->mtf_ctx[f].b250, "mtf_ctx.b250", SEC_CHROM_B250 + f*2);
        vb->mtf_ctx[f].encoding = header->encoding;
    }

    for (unsigned sf_i=0; sf_i < vb->num_info_subfields ; sf_i++) {
        
        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[section_i++]);

        MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, &vb->num_dict_ids, &vb->num_info_subfields, header->dict_id, 
                                                  SEC_INFO_SUBFIELD_DICT);

        zfile_uncompress_section (vb, header, &ctx->b250, "mtf_ctx.b250", SEC_INFO_SUBFIELD_B250);    
        ctx->encoding = header->encoding;    
    }

    // we allocate memory for the Buffer arrays only once the first time this VariantBlock
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

        unsigned num_samples_in_sb = (sb_i == vb->num_sample_blocks-1 ? global_num_samples % vb->num_samples_per_block : vb->num_samples_per_block);

        // if genotype data exists, it appears first
        if (vb->has_genotype_data) {
            zfile_uncompress_section (vb, vb->z_data.data + section_index[section_i++], &vb->genotype_sections_data[sb_i], "genotype_sections_data", SEC_GENOTYPE_DATA);

            // all genotype dictionaries are 8bit - for now (set during compression in seg_decide_encodings())
            for (unsigned did_i=0; did_i < MAX_DICTS; did_i++)
                if (vb->mtf_ctx[did_i].b250_section_type == SEC_GENOTYPE_DATA)
                    vb->mtf_ctx[did_i].encoding = B250_ENC_8;
        }
        
        // next, comes phase data
        if (vb->phase_type == PHASE_MIXED_PHASED) {
            
            zfile_uncompress_section (vb, vb->z_data.data + section_index[section_i++], &vb->phase_sections_data[sb_i], "phase_sections_data", SEC_PHASE_DATA);
            
            unsigned expected_size = vb->num_lines * num_samples_in_sb;
            ASSERT (vb->phase_sections_data[sb_i].len==expected_size, 
                    "Error: unexpected size of phase_sections_data[%u]: expecting %u but got %u", sb_i, expected_size, vb->phase_sections_data[sb_i].len)
        }

        // finally, comes haplotype data
        if (vb->has_haplotype_data) {
            
            zfile_uncompress_section (vb, vb->z_data.data + section_index[section_i++], &vb->haplotype_sections_data[sb_i], "haplotype_sections_data", SEC_HAPLOTYPE_DATA);
            
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

    if (vb->z_file->genozip_version > 1) {
        piz_uncompress_all_sections (vb);

        // combine all the sections of a variant block to regenerate the variant_data, haplotype_data,
        // genotype_data and phase_data for each row of the variant block
        piz_reconstruct_line_components (vb);
    }

    // v1 compatability
    else {
        void v1_piz_uncompress_all_sections (VariantBlockP vb); // forwwrd declaration - these are included at the end of this file
        v1_piz_uncompress_all_sections (vb);

        void v1_piz_reconstruct_line_components (VariantBlockP vb);
        v1_piz_reconstruct_line_components (vb);
    }

    buf_test_overflows(vb);

    COPY_TIMER (vb->profile.compute);

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway
}

// Called by PIZ I/O thread: read all the sections at the end of the file, before starting to process VBs
static int16_t piz_read_global_area (VariantBlock *pseudo_vb, bool need_random_access,
                                     Md5Hash *original_file_digest) // out
{
    File *zfile = pseudo_vb->z_file;

    int16_t data_type = zfile_read_genozip_header (pseudo_vb, original_file_digest);
    if (data_type == MAYBE_V1 || data_type == EOF) return data_type;

    // read dictionaries
    zfile_read_all_dictionaries (pseudo_vb, 0);
   
    // read random access, but only if we are going to need it
    if (need_random_access) {
        zfile_read_one_section (pseudo_vb, &pseudo_vb->z_data, "z_data", sizeof (SectionHeader), SEC_RANDOM_ACCESS);

        zfile_uncompress_section (pseudo_vb, pseudo_vb->z_data.data, &zfile->ra_buf, "ra_buf", SEC_RANDOM_ACCESS);

        zfile->ra_buf.len /= random_access_sizeof_entry();
        BGEN_random_access (&zfile->ra_buf);

        buf_free (&pseudo_vb->z_data);
    }

    file_seek (zfile, 0, SEEK_SET, false);

    return DATA_TYPE_VCF;
}

// returns true is successfully outputted a vcf file
bool piz_dispatcher (const char *z_basename, File *z_file, File *vcf_file, bool test_mode, unsigned max_threads, 
                     bool is_first_vcf_component, bool is_last_file)
{
    // static dispatcher - with flag_split, we use the same dispatcher when unzipping components
    static Dispatcher dispatcher = NULL;
    bool piz_successful = false;
    
    if (flag_split && !sections_has_more_vcf_components (z_file)) return false; // no more components

    if (!dispatcher) 
        dispatcher = dispatcher_init (max_threads, POOL_ID_UNZIP, 0, vcf_file, z_file, test_mode, is_last_file,
                                      !test_mode, // in test mode, we leave it to zip_dispatcher to display the progress indicator
                                      z_basename);
    VariantBlock *pseudo_vb = dispatcher_get_pseudo_vb (dispatcher);

    // read genozip header
    Md5Hash original_file_digest;
    
    // read genozip header and set the data type when reading the first vcf component of in case of --split, 
    static int16_t data_type = EOF; 
    if (is_first_vcf_component) {
        data_type = piz_read_global_area (pseudo_vb, true, &original_file_digest);

        if (data_type != MAYBE_V1)  // genozip v2+ - move cursor past first vcf header
            ASSERT (sections_get_next_header_type(z_file) == SEC_VCF_HEADER, "Error: unable to find VCF Header data in %s", file_printname (z_file));
    }

    if (data_type == EOF) goto finish;

    ASSERT (z_file->genozip_version > 1 || !flag_split || is_first_vcf_component,
            "Error: %s was compressed with genozip v1, splitting out components is not supported. Only the first file was extracted.",
            file_printname (z_file));

    // read and write VCF header. in split mode this also opens vcf_file
    piz_successful = (data_type != MAYBE_V1) ? vcf_header_genozip_to_vcf (pseudo_vb, &original_file_digest)
                                             : v1_vcf_header_genozip_to_vcf (pseudo_vb, &original_file_digest);

    if (!piz_successful) goto finish; // empty file - not an error
    
    vcf_file = pseudo_vb->vcf_file; // update local var - in case vcf file was opened by vcf_header_genozip_to_vcf()

    if (flag_split) 
        dispatcher_resume (dispatcher, vcf_file); // accept more input 

    // this is the dispatcher loop. In each iteration, it can do one of 3 things, in this order of priority:
    // 1. In input is not exhausted, and a compute thread is available - read a variant block and compute it
    // 2. Wait for the first thread (by sequential order) to complete and write data

    bool header_only_file = true; // initialize
    do {
        // PRIORITY 1: In input is not exhausted, and a compute thread is available - read a variant block and compute it
        if (!dispatcher_is_input_exhausted (dispatcher) && dispatcher_has_free_thread (dispatcher)) {

            bool compute = false;
            if (z_file->genozip_version > 1) {
                
                switch (sections_get_next_header_type(z_file)) {
                    case SEC_VB_HEADER:  
                        zfile_read_one_vb (dispatcher_generate_next_vb (dispatcher)); 
                        compute = true;
                        break;

                    case SEC_EOF: 
                        break; 

                    case SEC_VCF_HEADER: // 2nd+ vcf header of a concatenated file
                        if (!flag_split) {
                            vcf_header_genozip_to_vcf (pseudo_vb, NULL); // skip 2nd+ vcf header if concatenating
                            continue;
                        }
                        break; // eof if splitting

                    default: ABORT0 ("Error: unexpected section_type");
                }
            }
            else compute = v1_zfile_read_one_vb (dispatcher_generate_next_vb (dispatcher));  // genozip v1

            if (compute) {
                dispatcher_compute (dispatcher, piz_uncompress_variant_block);
                header_only_file = false;                
            }
            else { // eof
                dispatcher_input_exhausted (dispatcher);

                if (header_only_file)
                    dispatcher_finalize_one_vb (dispatcher, z_file, vcf_file->vcf_data_so_far, 0);
            }
        }

        // PRIORITY 2: Wait for the first thread (by sequential order) to complete and write data
        else if (dispatcher_has_processed_vb (dispatcher, NULL)) {
            VariantBlock *processed_vb = dispatcher_get_processed_vb (dispatcher, NULL); 
    
            vcffile_write_one_variant_block (vcf_file, processed_vb);

            z_file->vcf_data_so_far += processed_vb->vb_data_size; 

            dispatcher_finalize_one_vb (dispatcher, z_file, vcf_file->vcf_data_so_far, 0);
        }

    } while (!dispatcher_is_done (dispatcher));

    // verify file integrity, if the genozip compress was run with --md5
    static Buffer md5_verif_str = EMPTY_BUFFER;
    Md5Hash decompressed_file_digest;
    md5_finalize (&vcf_file->md5_ctx_concat, &decompressed_file_digest); // z_file might be a concatenation - this is the MD5 of the entire concatenation

    ASSERT (md5_is_equal (decompressed_file_digest, original_file_digest) || md5_is_zero (original_file_digest), // v1 files might be without md5
            "File integrity error: MD5 of decompressed file %s is %s, but the original VCF file's was %s", 
            vcf_file->name, md5_display (&decompressed_file_digest, false), md5_display (&original_file_digest, false));

    if (flag_split) file_close (&pseudo_vb->vcf_file, pseudo_vb); // close this component file

finish:
    // in split mode - we continue with the same dispatcher in the next component. otherwise, we finish with it here
    if (!flag_split || !piz_successful) {

        dispatcher_finish (&dispatcher, NULL);

        if (buf_is_allocated (&md5_verif_str)) {
            fprintf (stderr, "%.*s", md5_verif_str.len, md5_verif_str.data);
            buf_free (&md5_verif_str);
        }
    }
    else
        dispatcher_pause (dispatcher);

    return piz_successful;
}

#define V1_PIZ // select the piz functions of v1.c
#include "v1.c"
