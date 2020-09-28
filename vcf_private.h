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

    WordIndex format_mtf_i;  // the mtf_i into contexts[VCF_FORMAT].mtf and also format_mapper_buf that applies to this line. Data on the fields is in vb->format_mapper_buf[dl.format_mtf_i]
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
    bool is_an_before_ac, is_af_before_ac;

    // used for segging FORMAT/GT
    uint32_t gt_prev_ploidy;
    char gt_prev_phase;

    // dictionaries stuff 
    Buffer format_mapper_buf;         // ZIP only: an array of type Structured - one entry per entry in vb->contexts[VCF_FORMAT].mtf   
} VBlockVCF;

typedef VBlockVCF *VBlockVCFP;

// Samples stuff
extern void samples_digest_vcf_header (Buffer *vcf_header_buf);
extern char *vcf_samples_is_included;
#define samples_am_i_included(sample_i) (!flag_samples || ((bool)(vcf_samples_is_included[sample_i]))) // macro for speed - this is called in the critical loop of reconstructing samples

#endif

