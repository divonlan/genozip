// ------------------------------------------------------------------
//   vcf_zip.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <math.h>
#include "vcf_private.h"
#include "buffer.h"
#include "file.h"
#include "zfile.h"
#include "seg.h"
#include "vblock.h"
#include "context.h"
#include "zip.h"
#include "base250.h"
#include "endianness.h"
#include "strings.h"

#define DATA_LINE(i) ENT (ZipDataLineVCF, vb->lines, i)

uint32_t global_vcf_samples_per_block = 0; 

static pthread_mutex_t best_gt_data_compressor_mutex;

void vcf_zip_initialize (void)
{
    pthread_mutex_init (&best_gt_data_compressor_mutex, NULL);
}

void vcf_zip_set_global_samples_per_block (const char *num_samples_str)
{
    unsigned len = strlen (num_samples_str);
    for (unsigned i=0; i < len; i++) 
        ASSERT (IS_DIGIT (num_samples_str[i]), "Error: invalid argument of --sblock: %s. Expecting an integer number between 1 and 65535", num_samples_str);

    global_vcf_samples_per_block = atoi (num_samples_str);

    ASSERT (global_vcf_samples_per_block >= 1 && global_vcf_samples_per_block <= 65535, "Error: invalid argument of --sblock: %s. Expecting a number between 1 and 65535", num_samples_str);
}

#define SBL(line_i,sb_i) ((line_i) * vb->num_sample_blocks + (sb_i))

static unsigned vcf_zip_get_genotype_vb_start_len (VBlockVCF *vb)
{
    buf_alloc (vb, &vb->genotype_section_lens_buf, sizeof(unsigned) * vb->num_sample_blocks, 1, "section_lens_buf", 0);
    unsigned section_0_len = 0; // all sections are the same length except the last that might be shorter

    // offsets into genotype data of individual lines
    buf_alloc (vb, &vb->gt_sb_line_starts_buf, 
               vb->lines.len * vb->num_sample_blocks * sizeof(uint32_t*), 
               0, "gt_sb_line_starts_buf", vb->vblock_i);
    ARRAY (uint32_t *, gt_sb_line_starts, vb->gt_sb_line_starts_buf);
    
    // each entry is the length of a single line in a sample block
    buf_alloc (vb, &vb->gt_sb_line_lengths_buf, 
               vb->lines.len * vb->num_sample_blocks * sizeof(unsigned), 
               0, "gt_sb_line_lengths_buf", vb->vblock_i);
    ARRAY (unsigned, gt_sb_line_lengths, vb->gt_sb_line_lengths_buf);
    
    // calculate offsets and lengths of genotype data of each sample block
    for (uint32_t line_i=0; line_i < (uint32_t)vb->lines.len; line_i++) {

        uint32_t *gt_data  = (uint32_t*)GENOTYPE_DATA(vb, DATA_LINE (line_i));        
        unsigned format_mtf_i = DATA_LINE (line_i)->format_mtf_i;
        SubfieldMapper *format_mapper = ENT (SubfieldMapper, vb->format_mapper_buf, format_mtf_i);
        
        for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

            unsigned num_samples_in_sb = vcf_vb_num_samples_in_sb (vb, sb_i);

            gt_sb_line_starts[SBL(line_i, sb_i)] = 
                &gt_data[global_vcf_samples_per_block * sb_i * format_mapper->num_subfields];

            unsigned num_subfields_in_sample_line = format_mapper->num_subfields * num_samples_in_sb; // number of uint32_t
            gt_sb_line_lengths[SBL(line_i, sb_i)] = num_subfields_in_sample_line;

            if (!sb_i) section_0_len += num_subfields_in_sample_line;
        }
    }

    return section_0_len; // in subfields (uint32_t)
}

// split genotype data to sample groups, within a sample group genotypes are separated by a tab
static void vcf_zip_generate_genotype_one_section (VBlockVCF *vb, unsigned sb_i)
{
    START_TIMER;

    // build sample block genetype data
    ARRAY (uint8_t, dst_next, vb->genotype_one_section_data);
    
    // move the GT items from the line data to the permuted data - with each 
    // sample block of gt data containing the data in TRANSPOSED order - i.e. first
    // the gt data for all the variants for sample 1, then all of samples 2 etc.
    uint32_t num_samples_in_sb = vcf_vb_num_samples_in_sb (vb, sb_i);
    for (uint32_t sample_i=0; sample_i < num_samples_in_sb; sample_i++) {

        if (flag_show_gt_nodes) fprintf (stderr, "sample=%u (vb_i=%u sb_i=%u):\n", sb_i * global_vcf_samples_per_block + sample_i + 1, vb->vblock_i, sb_i);

        for (uint32_t line_i=0; line_i < (uint32_t)vb->lines.len; line_i++) {

            if (flag_show_gt_nodes) fprintf (stderr, "  L%u: ", line_i);

            ZipDataLineVCF *dl = DATA_LINE (line_i);

            SubfieldMapper *format_mapper = ENT (SubfieldMapper, vb->format_mapper_buf, dl->format_mtf_i);

            ARRAY (uint32_t *, sb_lines, vb->gt_sb_line_starts_buf);
            uint32_t *this_line = sb_lines[SBL(line_i, sb_i)];
            
            // lookup word indices in the global dictionary for all the subfields
            const uint8_t *dst_start = dst_next;

            int num_subfields = format_mapper->num_subfields;
            ASSERT (num_subfields >= 0, "Error: format_mapper->num_subfields=%d", format_mapper->num_subfields);

            // if this VB has subfields in some line, but not in this line, then we have filled it in vcf_seg_complete_missing_lines(), 
            // therefore we have 1 fake subfield
            if (vb->has_genotype_data > 0 && num_subfields==0) num_subfields = 1;

            for (unsigned sf=0; sf < format_mapper->num_subfields; sf++) { // iterate on the order as in the line
            
                uint32_t node_index = this_line[format_mapper->num_subfields * sample_i + sf];

                if (node_index <= WORD_INDEX_MAX_INDEX) { // normal index

                    Context *ctx = MAPPER_CTX (format_mapper, sf);
                    MtfNode *node = mtf_node_vb (ctx, node_index, NULL, NULL);
                    Base250 index = node->word_index;

                    if (flag_show_gt_nodes) fprintf (stderr, "%.*s:%u ", DICT_ID_LEN, ctx->name, index.n);

                    base250_copy (dst_next, index);
                    dst_next += base250_len (index.encoded.numerals);
                }
                else if (node_index == WORD_INDEX_MISSING_SF) {
                    *(dst_next++) = BASE250_MISSING_SF;
                }
                else {  // node_index == WORD_INDEX_EMPTY_SF
                    *(dst_next++) = BASE250_EMPTY_SF;
                }
            }

            vb->genotype_one_section_data.len += dst_next - dst_start;
            
            if (flag_show_gt_nodes) fprintf (stderr, "\n");
        }
    }

    COPY_TIMER (vb->profile.zip_generate_genotype_sections)
}

// split phase data to sample groups, in each group is a string of / | or -
static void vcf_zip_generate_phase_sections (VBlockVCF *vb)
{   
    START_TIMER;

    // we allocate memory for the Buffer array only once the first time this VBlockVCF
    // is used. Subsequent blocks reusing the memory will have the same number of samples (by VCF spec)
    if (!vb->phase_sections_data) 
        vb->phase_sections_data = (Buffer *)calloc (vb->num_sample_blocks, sizeof(Buffer)); // allocate once, never free

    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned num_samples_in_sb = vcf_vb_num_samples_in_sb (vb, sb_i); 
    
        // allocate memory for phase data for each sample block - one character per sample
        buf_alloc (vb, &vb->phase_sections_data[sb_i], vb->lines.len * num_samples_in_sb, 
                   0, "phase_sections_data", vb->vblock_i);

        // build sample block genetype data
        char *next = vb->phase_sections_data[sb_i].data;
        
        for (uint32_t line_i=0; line_i < (uint32_t)vb->lines.len; line_i++) {
            
            ZipDataLineVCF *dl = DATA_LINE (line_i);
            if (dl->phase_type == PHASE_MIXED_PHASED) 
                memcpy (next, &PHASE_DATA(vb,dl)[sb_i * num_samples_in_sb], num_samples_in_sb);
            else
                memset (next, (char)dl->phase_type, num_samples_in_sb);

            next += num_samples_in_sb;
        }
        vb->phase_sections_data[sb_i].len = num_samples_in_sb * (uint32_t)vb->lines.len;
    }
          
    // add back the phase data bytes that weren't actually "saved"
    COPY_TIMER (vb->profile.vcf_zip_generate_phase_sections)
}

typedef struct {
    int num_alt_alleles;
    unsigned index_in_original_line;
    unsigned index_in_sorted_line;
} HaploTypeSortHelperIndex;

static int sort_by_alt_allele_comparator(const void *p, const void *q)  
{ 
    int l = ((HaploTypeSortHelperIndex *)p)->num_alt_alleles; 
    int r = ((HaploTypeSortHelperIndex *)q)->num_alt_alleles;  
    return (l - r); 
}

static int sort_by_original_index_comparator(const void *p, const void *q)  
{ 
    int l = ((HaploTypeSortHelperIndex *)p)->index_in_original_line; 
    int r = ((HaploTypeSortHelperIndex *)q)->index_in_original_line;  
    return (l - r); 
}

static HaploTypeSortHelperIndex *vcf_zip_construct_ht_permutation_helper_index (VBlockVCF *vb)
{
    START_TIMER; 

    buf_alloc (vb, &vb->helper_index_buf, vb->num_haplotypes_per_line * sizeof(HaploTypeSortHelperIndex), 0,
               "helper_index_buf", vb->vblock_i);
    buf_zero (&vb->helper_index_buf);
    ARRAY (HaploTypeSortHelperIndex, helper_index, vb->helper_index_buf);

    // build index array 
    for (unsigned ht_i=0; ht_i < vb->num_haplotypes_per_line; ht_i++) 
        helper_index[ht_i].index_in_original_line = ht_i;

    for (uint32_t line_i=0; line_i < (uint32_t)vb->lines.len; line_i++) {
        for (unsigned ht_i=0; ht_i < vb->num_haplotypes_per_line; ht_i++) {

            // we count as alt alleles : 1 - 99 (ascii 49 to 147)
            //             ref alleles : 0 . (unknown) - (missing) * (ploidy padding)

            char *haplotype_data = HAPLOTYPE_DATA (vb, DATA_LINE (line_i));
            char one_ht = haplotype_data[ht_i];
            if (one_ht >= '1')
                helper_index[ht_i].num_alt_alleles++;
        }
    }
    COPY_TIMER (vb->profile.count_alt_alleles);

    return helper_index;
}

// sort haplogroups by alt allele count within the variant group, create an index for it, and split
// it to sample groups. for each sample a haplotype is just a string of 1 and 0 etc (could be other alleles too)
static void vcf_zip_generate_haplotype_sections (VBlockVCF *vb)
{
    START_TIMER;

    // we allocate memory for the Buffer array only once the first time this VBlockVCF
    // is used. Subsequent blocks reusing the memory will have the same number of samples (by VCF spec)
    if (!vb->haplotype_sections_data) 
        vb->haplotype_sections_data = (Buffer *)calloc (vb->num_sample_blocks, sizeof(Buffer)); // allocate once, never free
    
    // create a permutation index for the whole variant block, and permuted haplotypes for each sample block        
    buf_alloc (vb, &vb->haplotype_permutation_index, vb->num_haplotypes_per_line * sizeof(uint32_t), 
               0, "haplotype_permutation_index", vb->vblock_i);

    HaploTypeSortHelperIndex *helper_index = vcf_zip_construct_ht_permutation_helper_index (vb);

    // set dl->haplotype_ptr for all lines (for effeciency in the time loop below)
    for (unsigned line_i=0; line_i < vb->lines.len; line_i++) 
        DATA_LINE (line_i)->haplotype_ptr = HAPLOTYPE_DATA (vb, DATA_LINE (line_i));

    // now build per-sample-block haplotype array, picking haplotypes by the order of the helper index array
    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned num_haplotypes_in_sample_block = 
            vb->ploidy * vcf_vb_num_samples_in_sb (vb, sb_i); 

        unsigned helper_index_sb_i = sb_i * vb->num_samples_per_block * vb->ploidy;

        // sort the portion of the helper index related to this sample block. We sort by number of alt alleles.
        if (!flag_gtshark) // gtshark does better without sorting
            qsort (&helper_index[helper_index_sb_i], num_haplotypes_in_sample_block, sizeof (HaploTypeSortHelperIndex), sort_by_alt_allele_comparator);

        // allocate memory for haplotype data for each sample block - one character per haplotype
        buf_alloc (vb, &vb->haplotype_sections_data[sb_i], vb->lines.len * num_haplotypes_in_sample_block, 
                   0, "haplotype_sections_data", vb->vblock_i);

        // build sample block haplptype data - 
        // -- using the helper index to access the haplotypes in sorted order
        // -- transposing the array
        char *next = vb->haplotype_sections_data[sb_i].data;
        
        {   // this loop, tested with 1KGP data, takes up to 1/5 of total compute time, so its highly optimized
            START_TIMER;

            for (unsigned ht_i=0; ht_i < num_haplotypes_in_sample_block; ht_i++) {
                
                unsigned haplotype_data_char_i = helper_index[helper_index_sb_i + ht_i].index_in_original_line;
                const char **ht_data_ptr = &DATA_LINE (0)->haplotype_ptr; // this pointer moves sizeof(ZipDataLineVCF) bytes each iteration - i.e. to the exact same field in the next line

                for (unsigned line_i=0; line_i < vb->lines.len; line_i++, ht_data_ptr += sizeof(ZipDataLineVCF)/sizeof(char*)) 
                    *(next++) = (*ht_data_ptr)[haplotype_data_char_i];
            }
            COPY_TIMER (vb->profile.sample_haplotype_data);
        }
        vb->haplotype_sections_data[sb_i].len = num_haplotypes_in_sample_block * (uint32_t)vb->lines.len;
    }

    // final step - build the reverse index that will allow access by the original index to the sorted array
    // this will be included in the genozip file
    for (unsigned ht_i=0; ht_i < vb->num_haplotypes_per_line ; ht_i++)
        helper_index[ht_i].index_in_sorted_line = ht_i;

    // sort array back to its original order
    qsort (helper_index, vb->num_haplotypes_per_line, sizeof (HaploTypeSortHelperIndex), sort_by_original_index_comparator);

    // construct file index
    ARRAY (unsigned, hp_index, vb->haplotype_permutation_index);
    for (unsigned ht_i=0; ht_i < vb->num_haplotypes_per_line ; ht_i++)
        hp_index[ht_i] = helper_index[ht_i].index_in_sorted_line;

    buf_free (&vb->helper_index_buf);

    COPY_TIMER (vb->profile.vcf_zip_generate_haplotype_sections);
}

static CompressionAlg vcf_zip_get_best_gt_compressor (VBlock *vb, Buffer *test_data)
{
    static CompressionAlg best_gt_data_compressor = COMP_UNKNOWN;
    static Buffer compressed = EMPTY_BUFFER; // no thread issues as protected my mutex

    // get best compression algorithm for gt data - lzma or bzlib - their performance varies considerably with
    // the type of data - with either winning by a big margin
    pthread_mutex_lock (&best_gt_data_compressor_mutex);    

    if (best_gt_data_compressor != COMP_UNKNOWN) goto finish; // answer already known

    #define TEST_BLOCK_SIZE 100000
    buf_alloc (vb, &compressed, TEST_BLOCK_SIZE+1000, 1, "compressed_data_test", 0);

    uint32_t uncompressed_len = MIN (test_data->len, TEST_BLOCK_SIZE);

    uint32_t bzlib_comp_len = compressed.size;
    comp_compress_bzlib (vb, test_data->data, uncompressed_len, NULL, compressed.data, &bzlib_comp_len, false);
    
    uint32_t lzma_comp_len = compressed.size;
    comp_compress_lzma (vb, test_data->data, uncompressed_len, NULL, compressed.data, &lzma_comp_len, false);
    
    if      (bzlib_comp_len < uncompressed_len && bzlib_comp_len < lzma_comp_len) best_gt_data_compressor = COMP_BZ2;
    else if (lzma_comp_len  < uncompressed_len && lzma_comp_len < bzlib_comp_len) best_gt_data_compressor = COMP_LZMA;
    else                                                                          best_gt_data_compressor = COMP_PLN;

    buf_free (&compressed);

finish:
    pthread_mutex_unlock (&best_gt_data_compressor_mutex);
    return best_gt_data_compressor;
}

void vcf_zip_generate_ht_gt_compress_vb_header (VBlockP vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;
    
    if (vb->has_haplotype_data)
        vcf_seg_complete_missing_lines (vb);

    // if block has haplotypes - handle them now
    if (vb->has_haplotype_data)
        vcf_zip_generate_haplotype_sections (vb); 

    // if block has genetype data - calculate starts, lengths and allocate memory
    vb->max_genotype_section_len = vb->has_genotype_data ? vcf_zip_get_genotype_vb_start_len (vb) : 0; // number of b250s in one gt section matrix

    // if block has phase data - handle it
    if (vb->phase_type == PHASE_MIXED_PHASED) 
        vcf_zip_generate_phase_sections (vb);

    //unsigned variant_data_header_pos = vb->z_data.len;
    vcf_zfile_compress_vb_header (vb_); // variant data header + ht index    
}

// this function receives all lines of a vblock and processes them
// in memory to the compressed format. This thread then terminates the I/O thread writes the output.
void vcf_zip_compress_one_vb (VBlockP vb_)
{ 
    VBlockVCF *vb = (VBlockVCF *)vb_;
    
    // compress the sample data - genotype, haplotype and phase sections. genotype data is generated here too.
    CompressionAlg gt_data_alg;
    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {
        
        if (vb->has_genotype_data) {

            // in the worst case scenario, each subfield is represnted by 4 bytes in Base250
            buf_alloc (vb, &vb->genotype_one_section_data, vb->max_genotype_section_len * MAX_BASE250_NUMERALS, 1, "genotype_one_section_data", sb_i);

            // we compress each section at a time to save memory
            vcf_zip_generate_genotype_one_section (vb, sb_i); 

            gt_data_alg = vcf_zip_get_best_gt_compressor (vb_, &vb->genotype_one_section_data);

            COMPRESS_DATA_SECTION (SEC_VCF_GT_DATA, genotype_one_section_data, char, gt_data_alg, false); // gt data

            buf_free (&vb->genotype_one_section_data);
        }

        if (vb->phase_type == PHASE_MIXED_PHASED)
            COMPRESS_DATA_SECTION (SEC_VCF_PHASE_DATA, phase_sections_data[sb_i], char, COMP_BZ2, true); // optional - phase data

        if (vb->has_haplotype_data) {
            if (!flag_gtshark)
                COMPRESS_DATA_SECTION (SEC_VCF_HT_DATA, haplotype_sections_data[sb_i], char, COMP_BZ2, false) // ht data
            else 
                vcf_zfile_compress_haplotype_data_gtshark (vb_, &vb->haplotype_sections_data[sb_i], sb_i);
        }
    }
}

