// ------------------------------------------------------------------
//   vb.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef VB_INCLUDED
#define VB_INCLUDED

#include "genozip.h"
#include "buffer.h"
#include "profiler.h"
#include "aes.h"
#include "move_to_front.h"

typedef enum { PHASE_UNKNOWN      = '-',
               PHASE_HAPLO        = '1',
               PHASE_PHASED       = '|',
               PHASE_NOT_PHASED   = '/',
               PHASE_MIXED_PHASED = '+'    } PhaseType;

// IMPORTANT: if changing fields in DataLine, also update vb_release_vb
typedef struct {

    uint32_t line_i;

    // initially, data from vcf file line, later segregated to components "stored" in the overlay buffers below
    Buffer line;             

    // the following 4 buffers are overlay buffers onto line, so that they don't consume memory
    Buffer variant_data;     // string terminated by a newline. len includes the newline.
    Buffer genotype_data;    // \t separated genotype data for the line. last one has \t too. no \0. exists if any variant in the variant blck has FORMAT other that "GT"
    Buffer haplotype_data;   // length=ploidy*num_samples. exists if the GT subfield exists in any variant in the variant block

    Buffer phase_data;       // used only if phase is mixed. length=num_samples. exists if haplotype data exists and ploidy>=2
    PhaseType phase_type;    // phase type of this line
    
    bool has_haplotype_data; // FORMAT field contains GT
    bool has_genotype_data;  // FORMAT field contains subfields other than GT

    unsigned num_subfields;
    unsigned sf_i[MAX_SUBFIELDS]; // array in the order it appears in FORMAT - each an index into vb->mtf_ctx[]
} DataLine;

// IMPORTANT: if changing fields in VariantBlock, also update vb_release_vb
typedef struct variant_block_ {

    unsigned id;               // id of vb within the vb pool
    PoolId pool_id; 

    FileP vcf_file, z_file;  // pointers to objects that span multiple VBs

    // memory management
    Buffer buffer_list;        // a buffer containing an array of pointers to all buffers allocated for this VB (either by the I/O thread or its compute thread)

    bool ready_to_dispatch;    // line data is read, and dispatcher can dispatch this VB to a compute thread
    bool is_processed;         // thread completed processing this VB - it is ready for outputting
    bool in_use;               // this vb is in use
        
    DataLine *data_lines;      // if allocated, this array is of length global_max_lines_per_vb. for ZIP this is determined by the number of samples, and for PIZ by the SectionHeaderVCFHeader.max_lines_per_vb
    uint32_t num_lines;        // number of lines in this variant block
    uint32_t first_line;       // line number in VCF file (counting from 1), of this variant block
    uint32_t variant_block_i;  // number of variant block within VCF file

    // tracking execution
    uint32_t vb_data_size;     // size of variant block as it appears in the source file
    uint32_t max_gt_line_len;  // length of longest gt line in this vb after segregation

    ProfilerRec profile;

    // charactaristics of the data
    uint16_t ploidy;
    uint32_t num_sample_blocks;
    uint32_t num_samples_per_block; // except last sample block that may have less
    uint32_t num_haplotypes_per_line;
    bool has_genotype_data;    // if any variant has genotype data, then the block is considered to have it
    bool has_haplotype_data;   // ditto for haplotype data
    PhaseType phase_type;      // phase type of this variant block

    // chrom and pos
    char chrom[MAX_CHROM_LEN]; // a null-terminated ID of the chromosome
    int64_t min_pos, max_pos;  // minimum and maximum POS values in this VB. -1 if unknown
    uint64_t last_pos;         // value of POS field of the previous line, to do delta encoding
    bool is_sorted_by_pos;     // true if it this variant block is sorted by POS    

    // working memory for segregate - we segregate a line components into these buffers, and when done
    // we copy it back to DataLine - the buffers overlaying the line field
    Buffer line_variant_data;  // string terminated by a newline. len includes the newline.
    Buffer line_gt_data;       // \t separated genotype data for the line. last one has \t too. no \0. exists if any variant in the variant blck has FORMAT other that "GT"
    Buffer line_ht_data;       // length=ploidy*num_samples. exists if the GT subfield exists in any variant in the variant block
    Buffer line_phase_data;    // used only if phase is mixed. length=num_samples. exists if haplotype data exists and ploidy>=2

    // crypto stuff
    Buffer spiced_pw;  // used by crypt_generate_aes_key()
    uint8_t aes_round_key[240];// for 256 bit aes
    uint8_t aes_iv[AES_BLOCKLEN];
    int bi;

    // section data - ready to compress
    Buffer variant_data_section_data;    // all fields until FORMAT, newline-separated, \0-termianted. .len includes the terminating \0
    Buffer haplotype_permutation_index;
    Buffer optimized_gl_dict;  // GL dictionary data after extra optimization

    // these are Buffer arrays of size vb->num_sample_blocks allocated once when used for the first time and never freed. 
    // Subsequent variant blocks that re-use the memory have the same number of samples by VCF spec. Each buffer is a string of the data as written to the GENOZIP file

    // Note: The sample blocks for phase and genotype data are subsequent blocks of num_samples_per_block samples.
    // in contrast the sample blocks of haplotypes are num_samples_per_block*ploidy haplotypes as permuted. 
    // e.g. genotypes and phase in sample_block_i=1 don't necessary correspond to haplotypes in sample_block_i=1.

    Buffer *haplotype_sections_data;  // this is the haplotype character for each haplotype in the transposed sample group
    Buffer *phase_sections_data;      // this is the phase character for each genotype in the sample group
    Buffer *genotype_sections_data;   // this is for piz - each entry is a sample block, scanned columns first, each cell containing num_subfields indices (in base250 - 1 to 5 bytes each) into the subfield dictionaries
    Buffer genotype_one_section_data; // for zip we need only one section data

    // compresssed file data 
    Buffer z_data;                    // all headers and section data as read from disk

    int16_t z_next_header_i;          // next header of this VB to be encrypted or decrypted

    Buffer z_section_headers;         // (used by piz) an array of unsigned offsets of section headers within z_data

    Buffer gt_sb_line_starts_buf,     // used by zip_get_genotype_vb_start_len 
           gt_sb_line_lengths_buf,
           genotype_section_lens_buf; 

    Buffer helper_index_buf;          // used by zip_do_haplotypes

    Buffer vardata_header_buf;        // used by zfile_compress_variant_data

    Buffer compressed;                // used by various zfile functions
 
    Buffer ht_columns_data;           // used by piz_get_ht_permutation_lookups

    Buffer next_gt_in_sample;         // used for reconstructing genotype data by piz

    Buffer subfields_start_buf;       // these 3 are used by piz_reconstruct_line_components
    Buffer subfields_len_buf;
    Buffer num_subfields_buf;

    Buffer column_of_zeros;           // used by piz_get_ht_columns_data

    // subfields stuff 
    unsigned num_subfields;
    MtfContext mtf_ctx[MAX_SUBFIELDS];

    // Information content stats - how many bytes does this section have more than the corresponding part of the vcf file    
    int add_bytes[NUM_SEC_TYPES];                
    uint32_t vcf_section_bytes[NUM_SEC_TYPES];  // how many bytes did each section have in the original vcf file - should add up to the file size
    uint32_t z_section_bytes[NUM_SEC_TYPES];    // how many bytes does each section type have (including headers) in the genozip file - should add up to the file size

#   define NUM_COMPRESS_BUFS 4                  // bzlib2 compress requires 4 and decompress requires 2
    Buffer compress_bufs[NUM_COMPRESS_BUFS];    // memory allocation for compressor so it doesn't do its own malloc/free

} VariantBlock;

extern void vb_cleanup_memory (PoolId pool_id);
extern VariantBlock *vb_get_vb (PoolId pool_id, FileP vcf_file, FileP z_file, unsigned variant_block_i);
extern unsigned vb_num_samples_in_sb (const VariantBlock *vb, unsigned sb_i);
extern unsigned vb_num_sections(VariantBlock *vb);
extern void vb_release_vb (VariantBlock **vb_p);

typedef struct {
    unsigned num_vbs;
    VariantBlock vb[]; // variable length
} VariantBlockPool;
extern void vb_create_pool (PoolId pool_id, unsigned num_vbs);
extern VariantBlockPool *vb_get_pool (PoolId pool_id);

#endif