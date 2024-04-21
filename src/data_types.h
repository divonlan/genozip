
// ------------------------------------------------------------------
//   data_types.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "dict_id.h"
#include "genozip.h"
#include "digest.h"
#include "sections.h"
#include "flags.h"
#include "reference.h"
#include "txtfile.h"
#include "container.h"
#include "dict_id_gen.h"

// data type files
#include "vcf.h"
#include "sam.h"
#include "fastq.h"
#include "fasta.h"
#include "gff.h"
#include "me23.h"
#include "generic.h"
#include "locs.h"
#include "bed.h"

#ifndef MAX // can't use MAX_ in a global variable
#define MAX(a, b) (((a) > (b)) ? (a) : (b) )
#endif

typedef TRANSLATOR_FUNC ((*PizTranslator));
#define MAX_NUM_TRANSLATORS MAX (NUM_SAM_TRANS,      \
                            MAX (NUM_ME23_TRANS,     \
                                 NUM_LOCS_TRANS))

typedef enum { 
    HDR_NONE,   // no header for this data type
    HDR_OK,     // header is optional
    HDR_OK_0,   // header is optional for comp_i=0, and must not exist for comp_i > 0
    HDR_MUST,   // header must exist for every component
    HDR_MUST_0  // header must exist for comp_i=0 and must not exist for comp_i > 0
} TxtHeaderRequirement;

typedef struct DataTypeProperties {
    
    // ZIP properties and functions
    rom name;
    bool is_binary; 
    bool vb_end_nl;                             // the last character of every VB.txt_data is a newline
    bool uses_reference;                        // if --reference is given, it is needed for reconstruction
    DataType txt_type, bin_type;                // the textual and binary equivalents
    unsigned line_height;                       // how many actual txt file lines are in one seg line (seg lines are counted in lines.len). drop_lines in container_reconstruct needs to match the maximum.
    unsigned (*sizeof_vb)(DataType dt);
    unsigned (*sizeof_zip_dataline)(void);

    // TXT HEADER stuff
    TxtHeaderRequirement txt_header_required;
    rom hdr_contigs;                            // format of contigs in the txtheader
    char txt_header_1st_char;                   // first character in each line in the text file header (-1 if TXT_HEADER_IS_ALLOWED is false)

    #define HEADER_NEED_MORE         -1
    #define HEADER_DATA_TYPE_CHANGED -2         // when reading header, data type of txt/z_file is changed
    int32_t (*is_header_done) (bool is_eof);    // ZIP: header length if header read is complete, HEADER_NEED_MORE if not complete yet + sets lines.len
    int32_t (*unconsumed) (VBlockP, uint32_t first_i, int32_t *i);  // called by main thread called by txtfile_get_unconsumed_to_pass_to_next_vb to get the length of unconsumed txt to pass to next vb. returns -1 if first_i is too high and it needs more data.
    bool (*inspect_txt_header) (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags); // called by main thread to verify the txt header. returns false if this txt file should be skipped
    bool (*is_data_type) (STRp(header), bool *need_more); // ZIP: check for signature of the data type, to rescue generic files. If adding, also add to test.sh:test_redirected

    // ZIP callbacks - per txt file
    void (*zip_initialize)(void);               // called by main thread after the txt header is read
    void (*zip_after_segconf)(void);            // called by main thread after segconf has completed
    void (*zip_after_vbs)(void);                // called by main thread after all VBs of all components are exhausted and before compressing global sections
    void (*zip_finalize)(bool is_last_user_txt_file); // called by main thread after each txt file compressing is done

    // ZIP callbacks - per z_file
    void (*zip_free_end_of_z)(void);            // called by main thread after writing global area

    // ZIP callbacks - per VB
    void (*zip_init_vb)(VBlockP);               // called by main thread after reading txt of one vb into vb->txt_data
    void (*zip_after_compute)(VBlockP);         // called by main thread after completing compute thread of VB
    bool (*zip_dts_flag)(int dts);              // called to set FlagsGenozipHeader.dt_specific
    void (*zip_set_vb_header_specific)(VBlockP, SectionHeaderVbHeaderP); // main thread from zfile_compress_vb_header: set "specific" fields of VB Header
    void (*zip_set_txt_header_flags)(struct FlagsTxtHeader *f); // main thread from txtheader_compress: set data-type-specific fields of TXT Header
    void (*seg_initialize)(VBlockP);            // called by Compute thread at the beginning of Seg
    rom (*seg_txt_line)(VBlockP, rom field_start_line, uint32_t remaining_txt_len, bool *has_13);  // Called by Compute thread to Seg one line
    rom (*assseg_line)(VBlockP);                // called by ASSSEG to show errornous line before terminating . 
    bool (*seg_is_big)(ConstVBlockP, DictId, DictId); // returns true if dict_id is expected to have a large number of values 
    bool (*seg_is_small)(ConstVBlockP, DictId); // returns true if dict_id is known to have less than ~ 30K values (i.e. its context needs only a small hash table)
    void (*seg_finalize)(VBlockP);              // called by Compute thread at the end of Seg
    void (*zip_custom_merge)(VBlockP);          // called by ctx_merge_in_vb_ctx for a vb-specific non-context merge task
    bool seg_modifies;                          // true if Seg reuses txt_data, not including user requested modifications such as --optimize (in which case it should calculate vb->debug_line_hash)
    void (*zip_after_compress)(VBlockP);        // called by Compute thread after compressing contexts

    void (*generate_recon_plan)(void);
    void (*stats_reallocate)(void);             // called by stats_get_compressed_sizes - reallocate compressed sizes between contexts for stats reporting. total must not change.
    void (*zip_genozip_header)(SectionHeaderGenozipHeaderP); // called when writing the Genozip header

    // PIZ callbacks
    void (*piz_genozip_header)(ConstSectionHeaderGenozipHeaderP ); // called when reading the Genozip header
    void (*piz_after_global_area)(void);        // called once per z_file, after reading the global area
    bool (*piz_preprocess) (Dispatcher dispatcher); // called for preprocessing
    void (*piz_header_init)(void);              // called at the beginning of reconstruction of every txt header that is reconstructed
    bool (*piz_initialize)(CompIType comp_i);   // called at the beginning of each output txt_file (after the global sections have been read already)
    void (*piz_finalize)(bool is_last_z_file);  // called by main thread after each txt file reconstruction is done
    void (*piz_before_read)(VBlockP);           // called by main thread after reading VB_HEADER and before reading local/b250 sections from z_file
    bool (*piz_init_vb)(VBlockP, ConstSectionHeaderVbHeaderP);  // called by main thread after all sections of a VB have been read, before dispatching the compute thread
    void (*piz_recon_init)(VBlockP);            // called by the compute thread after uncompressing and before reconstructing
    void (*piz_init_line)(VBlockP);             // called before starting reconstruction of each line
    void (*piz_after_recon)(VBlockP);           // called by the compute thread after reconstruction from piz_reconstruct_one_vb. order of VBs is arbitrary
    void (*piz_process_recon)(VBlockP);         // called by the main thread after VB reconstruction, usually in order of VBs as in recon_plan, but out-of-order if --test with no writer

    void (*piz_after_preproc)(VBlockP);         // called by the main thread after joinig a preprocessing thread
    void (*piz_xtra_line_data)(VBlockP);        // called by container_verify_line_integrity to show extra line data for defective line

    Skip is_skip_section;
    void (*reconstruct_seq)(VBlockP, rom , unsigned, ReconType);
    CONTAINER_FILTER_FUNC ((*container_filter));// called for containers as defined in the container
    CONTAINER_CALLBACK ((*container_cb));       // called after of every repeat in containers with .callback=true
    CONTAINER_ITEM_CALLBACK ((*con_item_cb));   // called after reconstruction of an item, if CI1_ITEM_CB is specified

    // Special reconstruction functions corresponding to SNIP_SPECIAL snips
    unsigned num_special;
    PizSpecialReconstructor special[MAX_NUM_SPECIAL];

    // Container item translators, invoked after reconstruction, for translating from one data type to another
    unsigned num_translators;
    PizTranslator translator[MAX_NUM_TRANSLATORS];
    
    // misc properties and functions
    rom line_name;                              // the header displayed in --stats
    rom dtype_names[4];                         // the dictionary type names
} DataTypeProperties;

#define usz(type) ((unsigned)sizeof(type))
#define DATA_TYPE_PROPERTIES { \
/*    name         is_bin \n-end use_ref txt_type     bin_type    ht sizeof_vb       sizeof_zip_dataline   txt_headr   hdr_contigs     1st  is_header_done          unconsumed          inspect_txt_header      is_data_type  zip_initialize          zip_after_segconf  zip_after_vbs,     zip_finalize         zip_free_end_of_z      zip_init_vb        zip_after_compute          zip_dts_flag          zip_set_vb_header_specific        zip_set_txt_header_flags         seg_initialize          seg_txt_line          assseg_line          seg_is_big            seg_is_small          seg_finalize          zip_custom_merge seg_modifies zip_after_compress        generate_recon_plan          stats_reallocate        zip_genozip_header        piz_genozip_header        piz_after_global_area   piz_preprocess                          piz_header_init        piz_initialize          piz_finalize         piz_before_read        piz_init_vb          piz_recon_init              piz_init_line       piz_after_recon             piz_process_recon         piz_after_preproc      piz_xtra_line_data      is_skip_section             reconstruct_seq            container_filter       container_cb              con_item_cb          num_special          special          num_trans        translators         line_name        dtype_names                             */ \
    { "REFERENCE", false, false, false,  DT_REF,      DT_NONE,    1, fasta_vb_size,  fasta_vb_zip_dl_size, HDR_NONE,   NULL,           -1,  NULL,                   fasta_unconsumed,   NULL,                   is_ref,       ref_make_ref_init,      NULL,              NULL,              ref_make_finalize,   NULL,                  NULL,              ref_make_after_compute,    NULL,                 NULL,                             NULL,                            ref_make_seg_initialize,fasta_seg_txt_line,   NULL,                NULL,                 fasta_seg_is_small,   fasta_seg_finalize,   NULL,            false,       ref_make_create_range,    NULL,                        NULL,                   ref_make_genozip_header,  NULL,                     NULL,                   NULL,                                   NULL,                  NULL,                   NULL,                NULL,                  NULL,                NULL,                       NULL,               NULL,                       NULL,                     NULL,                  NULL,                   NULL,                       NULL,                      NULL,                  0,                        0,                   0,                   {},              0,               {},                 "line",          { "FIELD", "DESC",   "ERROR!"         } }, \
    { "VCF",       false, true,  true,   DT_VCF,      DT_NONE,    1, vcf_vb_size,    vcf_vb_zip_dl_size,   HDR_MUST,   VCF_CONTIG_FMT, '#', NULL,                   NULL,               vcf_inspect_txt_header, is_vcf,       vcf_zip_initialize,     NULL,              vcf_zip_after_vbs, vcf_zip_finalize,    vcf_header_finalize,   vcf_zip_init_vb,   vcf_zip_after_compute,     NULL,                 vcf_zip_set_vb_header_specific,   vcf_zip_set_txt_header_flags,    vcf_seg_initialize,     vcf_seg_txt_line,     NULL,                vcf_seg_is_big,       vcf_seg_is_small,     vcf_seg_finalize,     NULL,            false,       vcf_zip_after_compress,   NULL,                        NULL,                   vcf_zip_genozip_header,   vcf_piz_genozip_header,   NULL,                   NULL,                                   vcf_piz_header_init,   NULL,                   vcf_piz_finalize,    NULL,                  vcf_piz_init_vb,     vcf_piz_recon_init,         vcf_reset_line,     NULL,                       NULL,                     NULL,                  NULL,                   vcf_piz_is_skip_section,    NULL,                      vcf_piz_filter,        vcf_piz_container_cb,     vcf_piz_con_item_cb, NUM_VCF_SPECIAL,     VCF_SPECIAL,     0,               {},                 "variant",       { "FIELD", "INFO",   "FORMAT", "BOTH" } }, \
    { "SAM",       false, true,  true,   DT_SAM,      DT_BAM,     1, sam_vb_size,    sam_vb_zip_dl_size,   HDR_OK_0,   SAM_CONTIG_FMT, '@', NULL,                   NULL,               sam_header_inspect,     is_sam,       sam_zip_initialize,     sam_set_sag_type,  sam_zip_after_vbs, sam_zip_finalize,    sam_zip_free_end_of_z, sam_zip_init_vb,   sam_zip_after_compute,     sam_zip_dts_flag,     sam_zip_set_vb_header_specific,   NULL,                            sam_seg_initialize,     sam_seg_txt_line,     NULL,                sam_seg_is_big,       sam_seg_is_small,     sam_seg_finalize,     sam_deep_merge,  false,       sam_zip_after_compress,   sam_zip_generate_recon_plan, sam_stats_reallocate,   sam_zip_genozip_header,   sam_piz_genozip_header,   sam_piz_load_sags,      sam_piz_dispatch_one_load_sag_vb,       sam_piz_header_init,   sam_piz_initialize,     sam_piz_finalize,    NULL,                  sam_piz_init_vb,     sam_piz_recon_init,         sam_reset_line,     sam_piz_after_recon,        sam_piz_process_recon,    sam_piz_after_preproc, sam_piz_xtra_line_data, sam_piz_is_skip_section,    sam_reconstruct_SEQ_vs_ref,sam_piz_filter,        sam_piz_container_cb,     sam_piz_con_item_cb, NUM_SAM_SPECIAL,     SAM_SPECIAL,     NUM_SAM_TRANS,   SAM_TRANSLATORS,    "alignment",     { "FIELD", "QNAME",  "OPTION"         } }, \
    { "FASTQ",     false, true,  true,   DT_FASTQ,    DT_NONE,    4, fastq_vb_size,  fastq_vb_zip_dl_size, HDR_NONE,   NULL,           -1,  NULL,                   fastq_unconsumed,   NULL,                   is_fastq,     fastq_zip_initialize,   NULL,              NULL,              fastq_zip_finalize,  NULL,                  fastq_zip_init_vb, fastq_zip_after_compute,   NULL,                 NULL,                             fastq_zip_set_txt_header_flags,  fastq_seg_initialize,   fastq_seg_txt_line,   fastq_assseg_line,   NULL,                 fastq_seg_is_small,   fastq_seg_finalize,   NULL,            false,       NULL,                     NULL,                        NULL,                   fastq_zip_genozip_header, fastq_piz_genozip_header, NULL,                   NULL,                                   fastq_piz_header_init, fastq_piz_initialize,   NULL,                fastq_piz_before_read, fastq_piz_init_vb,   NULL,                       fastq_reset_line,   NULL,                       fastq_piz_process_recon,  NULL,                  NULL,                   fastq_piz_is_skip_section,  fastq_recon_aligned_SEQ,   fastq_piz_filter,      fastq_piz_container_cb,   0,                   NUM_FASTQ_SPECIAL,   FASTQ_SPECIAL,   0,               {},                 "read",          { "FIELD", "QNAME",  "AUX"            } }, \
    { "FASTA",     false, false, false,  DT_FASTA,    DT_NONE,    1, fasta_vb_size,  fasta_vb_zip_dl_size, HDR_NONE,   NULL,           -1,  NULL,                   fasta_unconsumed,   NULL,                   is_fasta,     fasta_zip_initialize,   NULL,              NULL,              NULL,                NULL,                  NULL,              fasta_zip_after_compute,   NULL,                 fasta_zip_set_vb_header_specific, NULL,                            fasta_seg_initialize,   fasta_seg_txt_line,   NULL,                fasta_seg_is_big,     fasta_seg_is_small,   fasta_seg_finalize,   NULL,            false,       NULL,                     NULL,                        NULL,                   NULL,                     NULL,                     NULL,                   NULL,                                   NULL,                  fasta_piz_initialize,   NULL,                NULL,                  fasta_piz_init_vb,   NULL,                       NULL,               NULL,                       NULL,                     NULL,                  NULL,                   fasta_piz_is_skip_section,  NULL,                      fasta_piz_filter,      NULL,                     0,                   NUM_FASTA_SPECIAL,   FASTA_SPECIAL,   0,               {},                 "line",          { "FIELD", "DESC",   "ERROR!"         } }, \
    { "GFF",       false, true,  false,  DT_GFF,      DT_NONE,    1, gff_vb_size,    0,                    HDR_OK,     NULL,           '#', NULL,                   gff_unconsumed,     gff_header_inspect,     is_gff,       gff_zip_initialize,     NULL,              NULL,              NULL,                NULL,                  NULL,              NULL,                      NULL,                 NULL,                             NULL,                            gff_seg_initialize,     gff_seg_txt_line,     NULL,                gff_seg_is_big,       gff_seg_is_small,     gff_seg_finalize,     NULL,            false,       NULL,                     NULL,                        NULL,                   NULL,                     NULL,                     NULL,                   NULL,                                   NULL,                  NULL,                   NULL,                NULL,                  gff_piz_init_vb,     NULL,                       gff_reset_line,     NULL,                       NULL,                     NULL,                  NULL,                   NULL,                       NULL,                      gff_piz_filter,        gff_piz_container_cb,     0,                   NUM_GFF_SPECIAL,     GFF_SPECIAL,     0,               {},                 "sequence",      { "FIELD", "ATTRS",  "ENST"           } }, \
    { "23ANDME",   false, true,  false,  DT_ME23,     DT_NONE,    1, 0,              0,                    HDR_MUST,   NULL,           '#', NULL,                   NULL,               me23_header_inspect,    is_me23,      NULL,                   NULL,              NULL,              NULL,                NULL,                  NULL,              NULL,                      NULL,                 NULL,                             NULL,                            me23_seg_initialize,    me23_seg_txt_line,    NULL,                NULL,                 me23_seg_is_small,    me23_seg_finalize,    NULL,            false,       NULL,                     NULL,                        NULL,                   NULL,                     NULL,                     NULL,                   NULL,                                   NULL,                  NULL,                   NULL,                NULL,                  NULL,                NULL,                       NULL,               NULL,                       NULL,                     NULL,                  NULL,                   NULL,                       NULL,                      NULL,                  NULL,                     0,                   0,                   {},              NUM_ME23_TRANS,  ME23_TRANSLATORS,   "SNP",           { "FIELD", "ERROR!", "ERROR!"         } }, \
    { "BAM",       true,  false, true,   DT_SAM,      DT_NONE,    0, sam_vb_size,    sam_vb_zip_dl_size,   HDR_MUST_0, SAM_CONTIG_FMT, -1,  bam_is_header_done,     bam_unconsumed,     sam_header_inspect,     is_bam,       sam_zip_initialize,     sam_set_sag_type,  sam_zip_after_vbs, sam_zip_finalize,    sam_zip_free_end_of_z, sam_zip_init_vb,   sam_zip_after_compute,     sam_zip_dts_flag,     sam_zip_set_vb_header_specific,   NULL,                            bam_seg_initialize,     bam_seg_txt_line,     bam_assseg_line,     sam_seg_is_big,       sam_seg_is_small,     sam_seg_finalize,     sam_deep_merge,  true,        sam_zip_after_compress,   sam_zip_generate_recon_plan, sam_stats_reallocate,   sam_zip_genozip_header,   sam_piz_genozip_header,   sam_piz_load_sags,      sam_piz_dispatch_one_load_sag_vb,       sam_piz_header_init,   sam_piz_initialize,     sam_piz_finalize,    NULL,                  sam_piz_init_vb,     sam_piz_recon_init,         sam_reset_line,     sam_piz_after_recon,        sam_piz_process_recon,    sam_piz_after_preproc, sam_piz_xtra_line_data, NULL,                       NULL,                      sam_piz_filter,        sam_piz_container_cb,     0/*cb only in SAM*/, NUM_SAM_SPECIAL,     SAM_SPECIAL,     NUM_SAM_TRANS,   SAM_TRANSLATORS,    "alignment",     { "FIELD", "DESC",   "OPTION"         } }, \
    { "BCF",       false, true,  true,   DT_VCF,      DT_NONE,    1, vcf_vb_size,    vcf_vb_zip_dl_size,   HDR_MUST,   VCF_CONTIG_FMT, '#', NULL,                   NULL,               vcf_inspect_txt_header, NULL,         vcf_zip_initialize,     NULL,              vcf_zip_after_vbs, vcf_zip_finalize,    vcf_header_finalize,   vcf_zip_init_vb,   vcf_zip_after_compute,     NULL,                 vcf_zip_set_vb_header_specific,   vcf_zip_set_txt_header_flags,    vcf_seg_initialize,     vcf_seg_txt_line,     NULL,                vcf_seg_is_big,       vcf_seg_is_small,     vcf_seg_finalize,     NULL,            false,       vcf_zip_after_compress,   NULL,                        NULL,                   vcf_zip_genozip_header,   vcf_piz_genozip_header,   NULL,                   NULL,                                   vcf_piz_header_init,   NULL,                   vcf_piz_finalize,    NULL,                  vcf_piz_init_vb,     vcf_piz_recon_init,         vcf_reset_line,     NULL,                       NULL,                     NULL,                  NULL,                   vcf_piz_is_skip_section,    NULL,                      vcf_piz_filter,        vcf_piz_container_cb,     vcf_piz_con_item_cb, NUM_VCF_SPECIAL,     VCF_SPECIAL,     0,               {},                 "variant",       { "FIELD", "INFO",   "FORMAT", "BOTH" } }, \
    { "GENERIC",   true,  false, false,  DT_GENERIC,  DT_GENERIC, 0, 0,              0,                    HDR_OK,     NULL,           -1,  generic_is_header_done, generic_unconsumed, NULL,                   NULL,         NULL,                   NULL,              NULL,              NULL,                NULL,                  NULL,              NULL,                      NULL,                 NULL,                             NULL,                            generic_seg_initialize, generic_seg_txt_line, generic_assseg_line, NULL,                 generic_seg_is_small, generic_seg_finalize, NULL,            false,       NULL,                     NULL,                        NULL,                   NULL,                     NULL,                     NULL,                   NULL,                                   NULL,                  NULL,                   NULL,                NULL,                  NULL,                NULL,                       NULL,               NULL,                       NULL,                     NULL,                  NULL,                   NULL,                       NULL,                      NULL,                  NULL,                     0,                   NUM_GENERIC_SPECIAL, GENERIC_SPECIAL, 0,               {},                 "N/A",           { "FIELD", "ERROR!", "ERROR!"         } }, \
    { "PHYLIP" }, /* PHYLIP: 9.0.11 - 15.0.41 */ \
    { "CHAIN"  }, /* CHAIN:  11.0.8 - 15.0.41 */ \
    { "KRAKEN" }, /* KRAKEN: 12.0.2 - 15.0.41 */ \
    { "LOCS",      true,  false, false,  DT_NONE,     DT_LOCS,    0, 0,              0,                    HDR_MUST,   NULL,           -1,  locs_is_header_done,    locs_unconsumed,    NULL,                   NULL,         NULL,                   NULL,              NULL,              NULL,                NULL,                  NULL,              NULL,                      NULL,                 NULL,                             NULL,                            locs_seg_initialize,    locs_seg_txt_line,    NULL,                NULL,                 locs_seg_is_small,    locs_seg_finalize,    NULL,            false,       NULL,                     NULL,                        NULL,                   NULL,                     NULL,                     NULL,                   NULL,                                   NULL,                  NULL,                   NULL,                NULL,                  NULL,                NULL,                       NULL,               NULL,                       NULL,                     NULL,                  NULL,                   NULL,                       NULL,                      NULL,                  NULL,                     0,                   NUM_LOCS_SPECIAL,   LOCS_SPECIAL,   NUM_LOCS_TRANS,  LOCS_TRANSLATORS,   "cluster" ,      { "FIELD", "ERROR!", "ERROR!"         } }, \
    { "BED",       false, true,  false,  DT_BED,      DT_NONE,    1, 0,              0,                    HDR_OK,     NULL,           -1,  bed_is_header_done,     NULL,               NULL,                   is_bed,       bed_zip_initialize,     NULL,              NULL,              NULL,                NULL,                  NULL,              NULL,                      NULL,                 NULL,                             NULL,                            bed_seg_initialize,     bed_seg_txt_line,     NULL,                NULL,                 bed_seg_is_small,     bed_seg_finalize,     NULL,            false,       NULL,                     NULL,                        NULL,                   NULL,                     NULL,                     NULL,                   NULL,                                   NULL,                  NULL,                   NULL,                NULL,                  NULL,                NULL,                       NULL,               NULL,                       NULL,                     NULL,                  NULL,                   NULL,                       NULL,                      NULL,                  NULL,                     0,                   NUM_BED_SPECIAL,    BED_SPECIAL,    0,               {},                 "line",          { "FIELD", "ERROR!", "ERROR!"         } }, \
}  
#define DATA_TYPE_FUNCTIONS_DEFAULT /* only applicable to (some) functions */ \
    { "DEFAULT",   false, false, false,  DT_NONE,     DT_NONE,    0, def_vb_size,    0,                    0,          0,              0,   def_is_header_done,     def_unconsumed,     0,                      NULL,         0,                      NULL,              0,                 0,                   0,                     0,                 0,                         0,                    NULL,                             NULL,                            0,                      0,                    textual_assseg_line, NULL,                 0,                    0,                    NULL,            false,       NULL,                     0,                           NULL,                   NULL,                     NULL,                     NULL,                   NULL,                                   0,                     0,                      NULL,                0,                     0,                   0,                          NULL,               0,                          NULL,                     NULL,                  NULL,                   0,                          0,                         default_piz_filter,    0,                        0,                   0,                  {},             0,               {},                 "",              { "FIELD", "DTYPE1", "DTYPE2"         } }

extern DataTypeProperties dt_props[NUM_DATATYPES], dt_props_def;

#define DT_(src,prop) ((src->data_type == DT_NONE) ? dt_props_def.prop : dt_props[src->data_type].prop)
#define DTP(prop)  DT_((vb),prop)
#define DTPZ(prop) DT_(z_file,prop) 
#define DTPT(prop) DT_(txt_file,prop)

#define ASSERT_DT_FUNC(src,prop) /* use this to assert the a function exists if it is mandatory */ \
    ASSERT (DT_(src,prop) || dt_props_def.prop, "undefined callback function " #src "->" #prop " for %s", dt_name(src->data_type))

// expession evaluates to def_value if neither the src nor the default function exists
#define DT_FUNC_OPTIONAL(src,prop,def_value) (src->data_type == DT_NONE || (!DT_(src,prop) && !dt_props_def.prop)) ? def_value : (DT_(src,prop) ? DT_(src,prop) : dt_props_def.prop) // caller adds function arguments here

// expession evaluates to 0 if neither the src nor the default function exists
// !!WARNING!!: to use DT_FUNC in an expression it MUST be enclosed in parathesis: (DT_FUNC(vb,func)(args))
#define DT_FUNC(src,prop) DT_FUNC_OPTIONAL (src, prop, 0)

typedef struct DataTypeFields {
    const unsigned num_fields;
    const Did pos, chrom; // the fields, or DID_NONE if this data type doesn't have them
    const Did qname; // some ID of the line (in addition to CHROM and POS) to be printed by ASSPIZ in case of an error
    DictId eol, toplevel;
    struct { 
        DictId dict_id; // these names go into the dictionary names on disk. to preserve backward compatibility, they should not be changed (names are not longer than 8=DICT_ID_LEN as the code assumes it)
        STR (tag_name);
    } predefined[MAX_NUM_PREDEFINED]; // consumed by ctx_initialize_predefined_ctxs 
} DataTypeFields;
#define CHROM (Did)0 // chrom is always the first field

#define TOPLEVEL "TOPLEVEL"

#define DATA_TYPE_FIELDS { \
/* num_fields        pos              chrom           qname            eol             toplevel             predefined        */ \
  {NUM_REF_FIELDS,   DID_NONE,        CHROM,          DID_NONE,        {},             {},                  REF_PREDEFINED    }, \
  {NUM_VCF_FIELDS,   VCF_POS,         CHROM,          DID_NONE,        {_VCF_EOL},     {_VCF_TOPLEVEL},     VCF_PREDEFINED    }, \
  {NUM_SAM_FIELDS,   SAM_POS,         CHROM,          SAM_QNAME,       {_SAM_EOL},     {_SAM_TOPLEVEL},     SAM_PREDEFINED    }, \
  {NUM_FASTQ_FIELDS, DID_NONE,        CHROM,          FASTQ_QNAME,     {_FASTQ_E1L},   {_FASTQ_TOPLEVEL},   FASTQ_PREDEFINED  }, \
  {NUM_FASTA_FIELDS, DID_NONE,        CHROM,          DID_NONE,        {_FASTA_EOL},   {_FASTA_TOPLEVEL},   FASTA_PREDEFINED  }, \
  {NUM_GFF_FIELDS,   GFF_START,       CHROM,          DID_NONE,        {_GFF_EOL},     {_GFF_TOPLEVEL},     GFF_PREDEFINED    }, \
  {NUM_ME23_FIELDS,  ME23_POS,        CHROM,          DID_NONE,        {_ME23_EOL},    {_ME23_TOPLEVEL},    ME23_PREDEFINED   }, \
  {NUM_SAM_FIELDS,   SAM_POS,         CHROM,          SAM_QNAME,       {_SAM_EOL},     {_SAM_TOP2BAM},      SAM_PREDEFINED    }, \
  {NUM_VCF_FIELDS,   VCF_POS,         CHROM,          DID_NONE,        {_VCF_EOL},     {_VCF_TOPLEVEL},     VCF_PREDEFINED    }, \
  {NUM_GNRIC_FIELDS, DID_NONE,        DID_NONE,       DID_NONE,        {},             {_GNRIC_TOPLEVEL},   GNRIC_PREDEFINED  }, \
  { /* obsolete PHYLIP*/ },                                                                                                      \
  { /* obsolete CHAIN */ },                                                                                                      \
  { /* obsolete KRKAEN*/ },                                                                                                      \
  {NUM_LOCS_FIELDS,  DID_NONE,        DID_NONE,       DID_NONE,        {},             {_LOCS_TOPLEVEL},    LOCS_PREDEFINED   }, \
  {NUM_BED_FIELDS,   BED_START,       CHROM,          DID_NONE,        {_BED_EOL},     {_BED_TOPLEVEL},     BED_PREDEFINED    }, \
}

extern DataTypeFields dt_fields[NUM_DATATYPES];
#define DTF(prop)  (dt_fields[vb->      data_type].prop)
#define DTFZ(prop) (dt_fields[z_file->  data_type].prop)
#define DTFT(prop) (dt_fields[txt_file->data_type].prop)

// list of ctx who's local data is compressed via a callback function
#define LOCAL_GET_LINE_CALLBACKS {          \
    SAM_LOCAL_GET_LINE_CALLBACKS(DT_SAM)    \
    SAM_LOCAL_GET_LINE_CALLBACKS(DT_BAM)    \
    FASTQ_LOCAL_GET_LINE_CALLBACKS          \
    FASTA_LOCAL_GET_LINE_CALLBACKS          \
}

// aliases - these are used only in PIZ, as a way for multiple dict_id's to get access to the same data storage
// on the ZIP side, the data is just placed directly in the destination ctx
#define DICT_ID_ALIASES {       \
    VCF_DICT_ID_ALIASES         \
    SAM_DICT_ID_ALIASES(DT_SAM) \
    SAM_DICT_ID_ALIASES(DT_BAM) \
    FASTQ_DICT_ID_ALIASES       \
    GFF_DICT_ID_ALIASES         \
}

typedef struct DtTranslation { 
    DataType src_z_non_bin_dt; 
    uint8_t src_z_is_binary;  // same type as z_file.flags
    DataType dst_txt_dt;
    DictId toplevel;
    float factor;             // specifies the size of vb.txt_data allocated in PIZ vs the size recorded in SectionHeaderVbHeader.recon_size.
                              // It is fully allocated in advance and cannot be extended.
    TXTHEADER_TRANSLATOR ((*txtheader_translator));
    bool trans_containers;    // whether to invoke translators in containers
    bool is_src_dt;           // this translation dst is actually the source file dt that was compressed (of which we have the MD5)
    bool (*is_translation)(VBlockP vb); // if non-NULL - called to determine if this VB should be translated
} DtTranslation;

#define TRANSLATIONS { /* note: translation i is not significant. translations can be added or removed. */ \
    /*                           non-bin-dt  binary   dst_txt_dt toplevel             factor  txtheader_transl.   trans_con. src_dt is_trans.             */ \
    /* SAM to BAM           */ { DT_SAM,     false,   DT_BAM,    { _SAM_TOP2BAM },       1,   sam_header_sam2bam,  true,     false, NULL                  }, /* BAM is expected to be smaller, edge cases will be handled by buf_alloc in sam_reset_line */ \
    /* Deep SAM to FASTQ    */ { DT_SAM,     false,   DT_FASTQ,  { _SAM_TOP2NONE },      1,   txtheader_sam2fq,    false,    false, NULL                  }, /* Deep file, and writting FASTQ-only data. SAM/BAM is partially reconstructed but not written */ \
    /* Binary SAM to SAM    */ { DT_SAM,     true,    DT_SAM,    { _SAM_TOPLEVEL },      2,   sam_header_bam2sam,  false,    false, NULL                  }, /* Since 15.0.28, segconf reports actual factor which can vary from 1 to 4. "factor" here is a fallback value. if factor is not sufficient we will buf_alloc in sam_reset_line */ \
    /* Binary SAM to BAM    */ { DT_SAM,     true,    DT_BAM,    { _SAM_TOP2BAM },       1,   NULL,                true,     true,  NULL                  }, /* BAM to BAM */ \
    /* Deep Bin. SAM to FQ  */ { DT_SAM,     true,    DT_FASTQ,  { _SAM_TOP2NONE },      1.5, txtheader_sam2fq,    false,    false, NULL                  }, /* Deep file, and writting FASTQ-only data. SAM/BAM is partially reconstructed but not written */ \
    /* 23andMe to VCF       */ { DT_ME23,    false,   DT_VCF,    { _ME23_TOP2VCF },      4,   txtheader_me232vcf,  true,     false, NULL                  }, \
    /* LOCS reconstruction  */ { DT_LOCS,    true,    DT_LOCS,   { _LOCS_TOPLEVEL },     1,   NULL,                true,     true,  NULL                  }, \
}

extern void dt_initialize (void);
extern const DtTranslation dt_get_translation (VBlockP vb);
extern DataType dt_get_txt_dt (DataType dt);
extern rom dt_name (DataType data_type);
extern rom z_dt_name (void);

