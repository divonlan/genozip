// ------------------------------------------------------------------
//   vcf_private.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef VCF_PRIVATE_INCLUDED
#define VCF_PRIVATE_INCLUDED

#include "vblock.h"
#include "vcf.h"

// IMPORTANT: if changing fields in DataLine, also update vb_release_vb
typedef struct {
    bool has_haplotype_data; // FORMAT field contains GT
    bool has_genotype_data;  // FORMAT field contains subfields other than GT

    WordIndex format_node_i;  // the node_index into contexts[VCF_FORMAT].nodes and also format_mapper_buf that applies to this line. Data on the fields is in vb->format_mapper_buf[dl.format_node_i]
} ZipDataLineVCF;

// IMPORTANT: if changing fields in VBlockVCF, also update vb_release_vb
typedef struct VBlockVCF {

    VBLOCK_COMMON_FIELDS

    // charactaristics of the data
    uint16_t ploidy;           // ZIP only
        
    // values stored during segging the INFO field and used when finalizing it - attempting to represent AC
    // as a function of AN and AF if possible
    const char *ac, *an, *af;  
    uint32_t ac_len, af_len, an_len;
    bool is_an_before_ac, is_af_before_ac, has_basecounts;
    
    Context *sf_ctx;
    enum { USE_SF_UNKNOWN, USE_SF_YES, USE_SF_NO } use_special_sf;
    Buffer last_sf, sf_snip; // INFO/SF data as it appears in the snip being constructed

    // used for segging FORMAT/GT
    uint32_t gt_prev_ploidy;
    char gt_prev_phase;

    // dictionaries stuff 
    Buffer format_mapper_buf;      // ZIP only: an array of type Container - one entry per entry in vb->contexts[VCF_FORMAT].nodes   

    // used by HT matrix codec 
    uint32_t num_haplotypes_per_line; 
    Context *ht_matrix_ctx; 
    
    // used by CODEC_HAPM (for VCF haplotype matrix) 
    Context *hapmat_index_ctx; 
    Buffer hapmat_helper_index_buf; // ZIP: used by codec_hapmat_count_alt_alleles 
    Buffer hapmat_columns_data;     // used by codec_hapmat_piz_get_one_line 
    Buffer hapmat_column_of_zeros;  // used by codec_hapmat_piz_calculate_columns   
    Buffer hapmat_one_array;        // one line or column 
    
    // used by CODEC_GTSHARK 
    Context *gtshark_gt_ctx, *gtshark_db_ctx, *gtshark_ex_ctx;
} VBlockVCF;

typedef VBlockVCF *VBlockVCFP;

// Samples stuff
extern uint32_t vcf_num_samples; // ZIP
extern uint32_t vcf_num_displayed_samples; // genocat with --samples
extern void samples_digest_vcf_header (Buffer *vcf_header_buf);
extern char *vcf_samples_is_included;
#define samples_am_i_included(sample_i) (!flag.samples || ((bool)(vcf_samples_is_included[sample_i]))) // macro for speed - this is called in the critical loop of reconstructing samples

#endif

