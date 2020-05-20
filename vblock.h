// ------------------------------------------------------------------
//   vblock.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef VBLOCK_INCLUDED
#define VBLOCK_INCLUDED

#include "genozip.h"
#include "buffer.h"
#include "profiler.h"
#include "aes.h"
#include "move_to_front.h"

extern void vb_vcf_release_vb(VBlockP), vb_sam_release_vb(VBlockP), vb_fast_release_vb(VBlockP), vb_gff3_release_vb(VBlockP);
extern void vb_vcf_destroy_vb(VBlockP), vb_sam_destroy_vb(VBlockP), vb_gff3_destroy_vb(VBlockP);
extern void vb_sam_initialize_vb(VBlockP);
extern void vb_vcf_cleanup_memory(VBlockP);

// PIZ only: can appear in did_i of an INFO subfield mapping, indicating that this INFO has an added
// ":#" indicating that the original VCF line had a Windows-style \r\n ending
#define DID_I_HAS_13 254 
#ifndef DID_I_NONE // also defined in move_to_front.h
#define DID_I_NONE   255
#endif
typedef struct SubfieldMapper {
    uint8_t num_subfields;        // (uint8_t)NIL if this mapper is not defined
    uint8_t did_i[MAX_SUBFIELDS]; // array in the order the subfields appears in FORMAT or INFO - each an index into vb->mtf_ctx[]
} SubfieldMapper;

#define MAPPER_CTX(mapper,sf) (((mapper)->did_i[(sf)] != (uint8_t)NIL) ? &vb->mtf_ctx[(mapper)->did_i[(sf)]] : NULL)

#define NUM_COMPRESS_BUFS 7   // bzlib2 compress requires 4 and decompress requires 2 ; lzma compress requires 7 and decompress 1

// IMPORTANT: if changing fields in VBlockVCF, also update vb_release_vb
#define VBLOCK_COMMON_FIELDS \
    uint32_t vblock_i;         /* number of variant block within VCF file */\
    int id;                    /* id of vb within the vb pool (-1 is the external vb) */\
    DataType data_type;        /* type of this VB */\
    \
    /* memory management  */\
    Buffer buffer_list;        /* a buffer containing an array of pointers to all buffers allocated for this VB (either by the I/O thread or its compute thread) */\
    \
    bool ready_to_dispatch;    /* line data is read, and dispatcher can dispatch this VB to a compute thread */\
    bool is_processed;         /* thread completed processing this VB - it is ready for outputting */\
    bool in_use;               /* this vb is in use */\
    \
    Buffer lines;              /* An array of *DataLine* - the lines in this VB */\
    uint32_t first_line;       /* PIZ only: line number in VCF file (counting from 1), of this variant block */\
    uint32_t num_lines_at_1_3, num_lines_at_2_3; /* ZIP VB=1 the number of lines segmented when 1/3 + 2/3 of estimate was reached  */\
    \
    /* tracking execution */\
    uint64_t vb_position_txt_file; /* position of this VB's data in the plain text file (i.e after decompression if the txt_file is compressed) */\
    int32_t vb_data_size;      /* ZIP: actual size of txt read file file ; PIZ: expected size of decompressed txt. Might be different than original if --optimize is used. */\
    uint32_t vb_data_read_size;/* ZIP only: amount of data read in txtfile_read_block() (either plain VCF or gz or bz2) for this VB */\
    uint32_t longest_line_len; /* length of longest line of text line in this vb */\
    uint32_t line_i;           /* ZIP: current line in VB (0-based) being segmented PIZ: current line in txt file */\
    \
    ProfilerRec profile; \
    \
    /* random access, chrom, pos */ \
    Buffer ra_buf;             /* ZIP only: array of RAEntry - copied to z_file at the end of each vb compression, then written as a SEC_RANDOM_ACCESS section at the end of the genozip file */\
    int32_t chrom_node_index;  /* ZIP: index into ra_buf: used by random_access_update_chrom/random_access_update_pos to sync between them */\
    int32_t last_pos;          /* value of POS field of the previous line, to do delta encoding - we do delta encoding even across chromosome changes */\
    \
    /* regions & filters */ \
    Buffer region_ra_intersection_matrix;  /* PIZ: a byte matrix - each row represents an ra in this vb, and each column is a region specieid in the command. the cell contains 1 if this ra intersects with this region */\
    \
    /* crypto stuff */\
    Buffer spiced_pw;  /* used by crypt_generate_aes_key() */\
    uint8_t aes_round_key[240];/* for 256 bit aes */\
    uint8_t aes_iv[AES_BLOCKLEN]; \
    int bi; \
    \
    /* file data */\
    Buffer z_data;                    /* all headers and section data as read from disk */\
    \
    Buffer txt_data;                  /* ZIP only: txt_data as read from disk - either the VCF header (in evb) or the VB data lines */\
    uint32_t txt_data_next_offset;    /* we re-use txt_data memory to overlay stuff in segregate */\
    Buffer txt_data_spillover;        /* when re-using txt_data, if it is too small, we spill over to this buffer */\
    \
    int16_t z_next_header_i;          /* next header of this VB to be encrypted or decrypted */\
    \
    Buffer z_section_headers;         /* PIZ only: an array of unsigned offsets of section headers within z_data */\
    \
    Buffer compressed;                /* used by various zfile functions */\
    \
    /* dictionaries stuff - we use them for 1. subfields with genotype data, 2. fields 1-9 of the VCF file 3. infos within the info field */\
    uint32_t num_dict_ids;            /* total number of dictionaries of all types */\
    MtfContext mtf_ctx[MAX_DICTS];    \
    uint8_t dict_id_to_did_i_map[65536];       /* map for quick look up of did_i from dict_id */\
    \
    /* Information content stats - how many bytes does this section have more than the corresponding part of the vcf file */\
    Buffer show_headers_buf;                   /* ZIP only: we collect header info, if --show-headers is requested, during compress, but show it only when the vb is written so that it appears in the same order as written to disk */\
    Buffer show_b250_buf;                      /* ZIP only: for collecting b250 during generate - so we can print at onces without threads interspersing */\
    Buffer section_list_buf;                   /* ZIP only: all the sections non-dictionary created in this vb. we collect them as the vb is processed, and add them to the zfile list in correct order of VBs. */\
    \
    Buffer compress_bufs[NUM_COMPRESS_BUFS];   /* memory allocation for compressor so it doesn't do its own malloc/free */

typedef struct VBlock {
    VBLOCK_COMMON_FIELDS
} VBlock;

extern void vb_cleanup_memory(void);
extern VBlock *vb_get_vb (unsigned vblock_i);
extern void vb_external_vb_initialize(void);
extern void vb_release_vb (VBlock *vb);

typedef struct {
    unsigned num_vbs; // length of array of pointers to VBlock
    unsigned num_allocated_vbs; // number of VBlocks allocated ( <= num_vbs )
    VBlock *vb[]; // variable length
} VBlockPool;
extern void vb_create_pool (unsigned num_vbs);
extern VBlockPool *vb_get_pool(void);

//-------------------------------
// VCF stuff
//-------------------------------

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

    uint32_t format_mtf_i;   // the mtf_i into mtf_ctx[VCF_FORMAT].mtf and also format_mapper_buf that applies to this line. Data on the fields is in vb->format_mapper_buf[dl.format_mtf_i]
    uint32_t info_mtf_i;     // the mtf_i into mtx_ctx[VCF_INFO].mtf and also iname_mapper_buf that applies to this line. Data on the infos is in  vb->iname_mapper_buf[dl.info_mtf_i]. either SubfieldInfoMapperPiz or SubfieldInfoZip
} ZipDataLineVCF;

// IMPORTANT: if changing fields in DataLine, also update vb_release_vb
typedef struct {
    bool has_haplotype_data; // FORMAT field contains GT
    bool has_genotype_data;  // FORMAT field contains subfields other than GT

    uint32_t format_mtf_i;   // the mtf_i into mtf_ctx[VCF_FORMAT].mtf and also format_mapper_buf that applies to this line. Data on the fields is in vb->format_mapper_buf[dl.format_mtf_i]
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
    
    Buffer gt_sb_line_starts_buf,     // used by zip_vcf_get_genotype_vb_start_len 
           gt_sb_line_lengths_buf,
           genotype_section_lens_buf; 

    Buffer helper_index_buf;          // used by zip_do_haplotypes
 
    Buffer ht_columns_data;           // used by piz_get_ht_permutation_lookups

    Buffer sample_iterator;           // an array of SnipIterator - one for each sample. used for iterate on gt samples to get one snip at a time 
     
    Buffer column_of_zeros;           // used by piz_vcf_get_ht_columns_data

    // dictionaries stuff 
    uint8_t num_info_subfields;       // e.g. if one inames is I1=I2=I3 and another one is I2=I3=I4= then we have two inames
                                      // entries in the mapper, which have we have num_info_subfields=4 (I1,I2,I3,I4) between them    
    Buffer iname_mapper_buf;          // ZIP only: an array of type SubfieldMapper - one entry per entry in vb->mtf_ctx[VCF_INFO].mtf
    uint8_t num_format_subfields;     // number of format subfields in this VB. num_subfields <= num_dict_ids-9.
    Buffer format_mapper_buf;         // ZIP only: an array of type SubfieldMapper - one entry per entry in vb->mtf_ctx[VCF_FORMAT].mtf   

    // stuff related to compressing haplotype data with gtshark
    Buffer gtshark_db_db_data;        // ZIP & PIZ
    Buffer gtshark_db_gt_data;        // ZIP & PIZ
    Buffer gtshark_exceptions_line_i; // ZIP & PIZ: uint32_t list of vb_line_i that have any allele >= '3'
    Buffer gtshark_exceptions_ht_i;   // ZIP & PIZ: delta-encoded (within the line) list of ht_i. For each exception line, there's the list of its ht_i's followed by a 0.
    Buffer gtshark_exceptions_allele; // ZIP & PIZ: each index (including terminating 0) corresponding to the index in exception_ht_i_offset
    Buffer gtshark_vcf_data;          // PIZ only

    // backward compatibility with genozip v1 
    Buffer v1_variant_data_section_data;  // all fields until FORMAT, newline-separated, \0-termianted. .len includes the terminating \0 (used for decompressed V1 files)
    Buffer v1_subfields_start_buf;        // v1 only: these 3 are used by piz_vcf_reconstruct_vb
    Buffer v1_subfields_len_buf;
    Buffer v1_num_subfields_buf;
} VBlockVCF;

extern unsigned vb_vcf_num_samples_in_sb (const VBlockVCF *vb, unsigned sb_i);

//-------------------------------
// SAM stuff
//-------------------------------

// IMPORTANT: if changing fields in DataLine, also update vb_release_vb
typedef struct {
    // the following 4 are indeces, lens into txt_data 
    uint32_t seq_data_start, qual_data_start, e2_data_start, u2_data_start, bd_data_start, bi_data_start; // start within vb->txt_data
    uint32_t seq_data_len, qual_data_len, e2_data_len, u2_data_len, bd_data_len, bi_data_len;             // length within vb->txt_data
    uint32_t seq_len;        // actual sequence length determined from any or or of: CIGAR, SEQ, QUAL. If more than one contains the length, they must all agree
} ZipDataLineSAM;

// IMPORTANT: if changing fields in VBlockSAM, also update vb_sam_release_vb and vb_sam_destroy_vb
typedef struct VBlockSAM {

    VBLOCK_COMMON_FIELDS

    SubfieldMapper qname_mapper;         // ZIP & PIZ

    uint32_t rname_index_minus_1;        // ZIP & PIZ: RNAME node index of previous line
    uint32_t rname_index_minus_2;        // ZIP & PIZ: RNAME node index of line before the previous line
    uint32_t rname_index_minus_3;        // ZIP & PIZ: RNAME node index of line before that

    const char *last_tlen_abs;           // ZIP & PIZ: last tlen segmented - pointer into vb->txt_data.data (ZIP) or TLEN dictionary (PIZ)        
    uint32_t last_tlen_abs_len;          // absolute value of last tlen, its sign, and the string length of last_tlen_abs
    bool last_tlen_is_positive;

    int32_t last_pnext_delta;            // last delta calculated for PNEXT

    // PIZ-only stuff
    Buffer optional_mapper_buf;          // PIZ: an array of type PizSubfieldMapper - one entry per entry in vb->mtf_ctx[SAM_OPTIONAL].mtf
} VBlockSAM;

//-----------------------------------------
// FASTQ & FASTA STUFF
//-----------------------------------------

// IMPORTANT: if changing fields in DataLine, also update vb_release_vb
typedef struct {
    uint32_t seq_data_start, qual_data_start; // start within vb->txt_data
    uint32_t seq_len;                         // length within vb->txt_data (in case of FASTQ, this length is also applies to quality, per FASTQ spec)
} ZipDataLineFAST;

// IMPORTANT: if changing fields in VBlockFASTQQ, also update vb_fastq_release_vb and vb_fastq_destroy_vb
typedef struct VBlockFAST { // for FASTA and FASTQ

    VBLOCK_COMMON_FIELDS

    // shared fields FASTA and FASTQ
    SubfieldMapper desc_mapper; // ZIP & PIZ
    enum { GS_READ, GS_TEST, GS_UNCOMPRESS } grep_stages;   // PIZ: tell zfile_is_skip_section what to skip

    // FASTA-only fields
    bool fasta_prev_vb_last_line_was_grepped; // PIZ: whether previous VB's last line was grepped successfully. 

    // note: last_line is initialized to FASTA_LINE_SEQ (=0) so that a ; line as the first line of the VB
    // is interpreted as a description, not a comment
    FastaFields last_line; // ZIP

} VBlockFAST;

//-----------------------------------------
// GFF3 STUFF
//-----------------------------------------

// IMPORTANT: if changing fields in DataLine, also update vb_release_vb
typedef struct {
    uint32_t attrs_mtf_i;     // the mtf_i into mtx_ctx[VCF_INFO].mtf and also iname_mapper_buf that applies to this line. Data on the infos is in  vb->iname_mapper_buf[dl.info_mtf_i]. either SubfieldInfoMapperPiz or SubfieldInfoZip
} ZipDataLineGFF3;

// IMPORTANT: if changing fields in VBlockFASTQQ, also update vb_fastq_release_vb and vb_fastq_destroy_vb
typedef struct VBlockGFF3 {

    VBLOCK_COMMON_FIELDS

    uint32_t last_id;             // used for detla'ing the ID subfield

    uint8_t num_info_subfields;   // e.g. if one inames is I1=I2=I3 and another one is I2=I3=I4= then we have two inames
                                  // entries in the mapper, which have we have num_info_subfields=4 (I1,I2,I3,I4) between them    
    Buffer iname_mapper_buf;      // ZIP only: an array of type SubfieldMapper - one entry per entry in vb->mtf_ctx[VCF_INFO].mtf

} VBlockGFF3;

#endif

//-----------------------------------------
// VBlock utilities
//-----------------------------------------

// NOT thread safe, use only in execution-terminating messages
extern const char *err_vb_pos (void *vb);
