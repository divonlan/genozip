// ------------------------------------------------------------------
//   vcf_zfile.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "vcf_private.h"
#include "sections.h"
#include "endianness.h"
#include "crypt.h"
#include "zfile.h"

// called by ZIP compute thread - called from vcf_zip_compress_one_vb()
void vcf_compress_vb_header (VBlock *vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    uint32_t my_squeeze_len = squeeze_len (vb->num_haplotypes_per_line);
    uint32_t sizeof_header = sizeof (SectionHeaderVbHeaderVCF);

    SectionHeaderVbHeaderVCF vb_header;
    memset (&vb_header, 0, sizeof(SectionHeaderVbHeaderVCF)); // safety
    
    vb_header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    vb_header.h.section_type          = SEC_VB_HEADER;
    vb_header.h.data_uncompressed_len = BGEN32 (my_squeeze_len);
    vb_header.h.compressed_offset     = BGEN32 (sizeof_header);
    vb_header.h.vblock_i              = BGEN32 (vb->vblock_i);
    vb_header.h.sec_compression_alg   = COMP_NONE; // in versions 2-4 it was (needlessly) compressed with bzlib
    vb_header.num_lines               = BGEN32 ((uint32_t)vb->lines.len);
    vb_header.longest_line_len        = BGEN32 (vb->longest_line_len);
    vb_header.md5_hash_so_far         = vb->md5_hash_so_far;
    vb_header.phase_type              = (char)vb->phase_type; 
    vb_header.has_genotype_data       = vb->has_genotype_data;
    vb_header.num_samples             = BGEN32 (global_vcf_num_samples);
    vb_header.num_haplotypes_per_line = BGEN32 (vb->num_haplotypes_per_line);
    vb_header.num_samples_per_block   = BGEN32 (vb->num_samples_per_block);
    vb_header.num_sample_blocks       = BGEN32 (vb->num_sample_blocks);
    vb_header.ploidy                  = BGEN16 (vb->ploidy);
    vb_header.vb_data_size            = BGEN32 (vb->vb_data_size);
    vb_header.max_gt_line_len         = BGEN32 (vb->max_gt_line_len);
    vb_header.is_gtshark              = flag_gtshark;

    // create squeezed index - IF we have haplotype data AND more than one haplotype per line (i.e. my_squeeze_len > 0)
    if (my_squeeze_len) {
        buf_alloc (vb, &vb->haplotype_permutation_index_squeezed, my_squeeze_len, 1, "haplotype_permutation_index_squeezed", 0);

        uint16_t haplotype_index_checksum;
        squeeze (vb, (uint8_t *)vb->haplotype_permutation_index_squeezed.data, &haplotype_index_checksum,
                (unsigned *)vb->haplotype_permutation_index.data, 
                vb->num_haplotypes_per_line);

        vb_header.haplotype_index_checksum = BGEN16(haplotype_index_checksum);
    }
    
    // compress section into z_data - to be eventually written to disk by the I/O thread
    vb->z_data.name = "z_data"; // this is the first allocation of z_data - comp_compress requires that it is pre-named
    vb->z_data.param = vb->vblock_i;
    comp_compress (vb_, &vb->z_data, false, (SectionHeader*)&vb_header, 
                   my_squeeze_len ? vb->haplotype_permutation_index_squeezed.data : NULL, NULL);
}

// ZIP only: called by the I/O thread in the sequential order of VBs: updating of the already compressed
// variant data section (compressed by the compute thread in vcf_compress_vb_header) just before writing it to disk
// note: this updates the z_data in memory (not on disk)
void vcf_update_compressed_vb_header (VBlock *vb, uint32_t txt_first_line_i)
{
    SectionHeaderVbHeaderVCF *vb_header = (SectionHeaderVbHeaderVCF *)vb->z_data.data;
    vb_header->z_data_bytes = BGEN32 ((uint32_t)vb->z_data.len);
    vb_header->first_line   = BGEN32 (txt_first_line_i);

    if (flag_show_vblocks) {
        fprintf (stderr, "vb_i=%u first_line=%u num_lines=%u uncomprssed=%u compressed=%u num_sample_blocks=%u ht_per_line=%u max_gt_line_len=%u\n",
                 vb->vblock_i, txt_first_line_i, BGEN32 (vb_header->num_lines), 
                 BGEN32 (vb_header->vb_data_size), BGEN32 (vb_header->z_data_bytes), 
                 BGEN32 (vb_header->num_sample_blocks), BGEN32 (vb_header->num_haplotypes_per_line), 
                 BGEN32 (vb_header->max_gt_line_len));
    }
    // now we can finally encrypt the header - if needed
    if (crypt_have_password())  
        crypt_do (vb, (uint8_t*)vb_header, BGEN32 (vb_header->h.compressed_offset),
                  BGEN32 (vb_header->h.vblock_i), vb_header->h.section_type, true);
}

void vcf_zfile_compress_haplotype_data_gtshark (VBlock *vb_, const Buffer *haplotype_sections_data, unsigned sb_i)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;
    
    gtshark_compress_haplotype_data (vb, haplotype_sections_data, sb_i); // populates vb->gtshark_*
    
    vb->gtshark_exceptions_line_i.len *= sizeof (uint32_t);
    vb->gtshark_exceptions_ht_i.len   *= sizeof (uint16_t);

    zfile_compress_section_data     (vb_, SEC_VCF_HT_GTSHARK, &vb->gtshark_exceptions_line_i);
    zfile_compress_section_data     (vb_, SEC_VCF_HT_GTSHARK, &vb->gtshark_exceptions_ht_i);
    zfile_compress_section_data     (vb_, SEC_VCF_HT_GTSHARK, &vb->gtshark_exceptions_allele);
    zfile_compress_section_data_alg (vb_, SEC_VCF_HT_GTSHARK, &vb->gtshark_db_db_data, 0, 0, COMP_NONE);
    zfile_compress_section_data_alg (vb_, SEC_VCF_HT_GTSHARK, &vb->gtshark_db_gt_data, 0, 0, COMP_NONE);

    // free buffers - they will be needed by the next section
    buf_free (&vb->gtshark_exceptions_line_i);
    buf_free (&vb->gtshark_exceptions_ht_i);
    buf_free (&vb->gtshark_exceptions_allele);
    buf_free (&vb->gtshark_db_db_data);
    buf_free (&vb->gtshark_db_gt_data);
}
