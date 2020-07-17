// ------------------------------------------------------------------
//   vcf_private.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef VCF_PRIVATE_INCLUDED
#define VCF_PRIVATE_INCLUDED

#include "vblock.h"
#include "vcf.h"

typedef enum { PHASE_UNKNOWN      = '-',
               PHASE_HAPLO        = '1',
               PHASE_PHASED       = '|',
               PHASE_NOT_PHASED   = '/',
               PHASE_MIXED_PHASED = '+'    } PhaseType;

#define GENOTYPE_DATA(vb,dl)  ((dl)->genotype_data_spillover  ? &(vb)->txt_data_spillover.data[(dl)->genotype_data_start] \
                                                              : &(vb)->txt_data.data[(dl)->genotype_data_start])
#define HAPLOTYPE_DATA(vb,dl) ((dl)->haplotype_data_spillover ? &(vb)->txt_data_spillover.data[(dl)->haplotype_data_start] \
                                                              : &(vb)->txt_data.data[(dl)->haplotype_data_start])
#define PHASE_DATA(vb,dl)     ((dl)->phase_data_spillover     ? &(vb)->txt_data_spillover.data[(dl)->phase_data_start] \
                                                              : &(vb)->txt_data.data[(dl)->phase_data_start])
// IMPORTANT: if changing fields in DataLine, also update vb_release_vb
typedef struct {
    // the following 3 are indeces, lens into txt_data or txt_data_spillover. 
    bool genotype_data_spillover, haplotype_data_spillover, phase_data_spillover;
    uint32_t genotype_data_start, haplotype_data_start, phase_data_start;
    uint32_t genotype_data_len, haplotype_data_len, phase_data_len;

    const char *haplotype_ptr; // this is HAPLOTYPE_DATA - set for effeciency AFTER segregation is complete and there will be no more reallocs of txt_data_spillover

    PhaseType phase_type;    // phase type of this line
    
    bool has_haplotype_data; // FORMAT field contains GT
    bool has_genotype_data;  // FORMAT field contains subfields other than GT

    uint32_t format_mtf_i;   // the mtf_i into contexts[VCF_FORMAT].mtf and also format_mapper_buf that applies to this line. Data on the fields is in vb->format_mapper_buf[dl.format_mtf_i]
} ZipDataLineVCF;

// IMPORTANT: if changing fields in DataLine, also update vb_release_vb
typedef struct {
    bool has_haplotype_data; // FORMAT field contains GT
    bool has_genotype_data;  // FORMAT field contains subfields other than GT

    uint32_t format_mtf_i;   // the mtf_i into contexts[VCF_FORMAT].mtf and also format_mapper_buf that applies to this line. Data on the fields is in vb->format_mapper_buf[dl.format_mtf_i]
} PizDataLineVCF;

// IMPORTANT: if changing fields in VBlockVCF, also update vb_release_vb
typedef struct VBlockVCF {

    VBLOCK_COMMON_FIELDS

    // charactaristics of the data
    uint16_t ploidy;
    
    uint32_t num_sample_blocks;
    uint32_t num_samples_per_block; // except last sample block that may have less
    uint32_t num_haplotypes_per_line;
    uint32_t max_genotype_section_len; // number of b250s in one gt section matrix
    bool has_genotype_data;    // if any variant has genotype data, then the block is considered to have it
    bool has_haplotype_data;   // ditto for haplotype data
    PhaseType phase_type;      // phase type of this variant block
    
    // working memory for segregate - we segregate a line components into these buffers, and when done
    // we copy it back to DataLine - the buffers overlaying the line field
    Buffer line_gt_data;       // \t separated genotype data for the line. last one has \t too. no \0. exists if any variant in the variant blck has FORMAT other that "GT"
    Buffer line_ht_data;       // length=ploidy*num_samples. exists if the GT subfield exists in any variant in the variant block
    Buffer line_phase_data;    // used only if phase is mixed. length=num_samples. exists if haplotype data exists and ploidy>=2
    uint32_t max_gt_line_len;  // length of longest gt line in this vb after segregation 
    
    // values stored during segging the INFO field and used when finalizing it - attempting to represent AC
    // as a function of AN and AF if possible
    const char *ac, *an, *af;  
    uint32_t ac_len, af_len, an_len;
    bool is_an_before_ac, is_af_before_ac;

    // section data - ready to compress
    Buffer haplotype_permutation_index;
    Buffer haplotype_permutation_index_squeezed; // used by piz to unsqueeze the index and zfile to compress it
    Buffer optimized_gl_dict;         // GL dictionary data after extra optimization

    // these are Buffer arrays of size vb->num_sample_blocks allocated once when used for the first time and never freed,
    // unless a subsequent VCF file has more sample blocks, in which case they are first destroyed by vb_cleanup_memory().
    // Each buffer is a string of the data as written to the GENOZIP file

    Buffer *haplotype_sections_data;  // ZIP & PIZ: this is the haplotype character for each haplotype in the transposed sample block
    Buffer *phase_sections_data;      // ZIP & PIZ: this is the phase character for each genotype in the sample block
    Buffer *genotype_sections_data;   // PIZ only:  each entry is a sample block, scanned columns first, each cell containing num_subfields indices (in base250 - 1 to 5 bytes each) into the subfield dictionaries
    Buffer is_sb_included;            // PIZ only:  array of bool indicating for each sample block whether it is included, based on --samples 
    Buffer genotype_one_section_data; // ZIP only:  for zip we need only one section data
    
    Buffer gt_sb_line_starts_buf,     // used by vcf_zip_get_genotype_vb_start_len 
           gt_sb_line_lengths_buf,
           genotype_section_lens_buf; 

    Buffer helper_index_buf;          // used by zip_do_haplotypes
 
    Buffer ht_columns_data;           // used by piz_get_ht_permutation_lookups

    Buffer sample_iterator;           // an array of SnipIterator - one for each sample. used for iterate on gt samples to get one snip at a time 
     
    Buffer column_of_zeros;           // used by vcf_piz_get_ht_columns_data

    // dictionaries stuff 
    Buffer format_mapper_buf;         // ZIP only: an array of type SubfieldMapper - one entry per entry in vb->contexts[VCF_FORMAT].mtf   

    // stuff related to compressing haplotype data with gtshark
    Buffer gtshark_db_db_data;        // ZIP & PIZ
    Buffer gtshark_db_gt_data;        // ZIP & PIZ
    Buffer gtshark_exceptions_line_i; // ZIP & PIZ: uint32_t list of vb_line_i that have any allele >= '3'
    Buffer gtshark_exceptions_ht_i;   // ZIP & PIZ: delta-encoded (within the line) list of ht_i. For each exception line, there's the list of its ht_i's followed by a 0.
    Buffer gtshark_exceptions_allele; // ZIP & PIZ: each index (including terminating 0) corresponding to the index in exception_ht_i_offset
    Buffer gtshark_vcf_data;          // PIZ only
} VBlockVCF;

typedef VBlockVCF *VBlockVCFP;

extern unsigned vcf_vb_num_samples_in_sb (const VBlockVCF *vb, unsigned sb_i);
extern uint32_t global_vcf_samples_per_block;
extern void vcf_seg_complete_missing_lines (VBlockVCFP vb);

// Samples stuff
extern void samples_digest_vcf_header (Buffer *vcf_header_buf);
extern char *vcf_samples_is_included;
#define samples_am_i_included(sample_i) (!flag_samples || ((bool)(vcf_samples_is_included[sample_i]))) // macro for speed - this is called in the critical loop of reconstructing samples
extern bool samples_get_is_sb_included (uint32_t num_samples_per_block, uint32_t sb_i);

// Squeeze stuff
extern unsigned squeeze_len (unsigned int len);
extern void squeeze (VBlockVCFP vb, uint8_t *dst, uint16_t *squeezed_checksum, const unsigned *src,unsigned src_len);
extern void unsqueeze (VBlockVCFP vb, unsigned *normal, const uint8_t *squeezed, uint16_t squeezed_checksum,unsigned normal_len);

// ZFILE stuff
extern void vcf_zfile_compress_haplotype_data_gtshark (VBlockP vb, ConstBufferP haplotype_sections_data, unsigned sb_i);

// GTshark stuff
extern void gtshark_compress_haplotype_data (VBlockVCFP vb, ConstBufferP section_data, unsigned sb_i);
extern void gtshark_uncompress_haplotype_data (VBlockVCFP vb, unsigned sb_i);

#endif

