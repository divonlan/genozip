// ------------------------------------------------------------------
//   segconf.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <stdarg.h>
#include "genozip.h"
#include "vblock.h"
#include "file.h"
#include "txtfile.h"
#include "seg.h"
#include "segconf.h"
#include "strings.h"
#include "codec.h"
#include "arch.h"
#include "bgzf.h"
#include "tip.h"

SegConf segconf = {}; // system-wide global
static VBlockP segconf_vb = NULL;

// figure out if file is sorted or not
void segconf_test_sorted (VBlockP vb, WordIndex prev_line_chrom, PosType32 pos, PosType32 prev_line_pos)
{
    // evidence of not being sorted: our CHROM is the same as the previous line, but POS has decreased
    if (segconf.is_sorted && prev_line_chrom == vb->chrom_node_index && prev_line_pos > pos)
        segconf.is_sorted = false;

    // evidence of not being sorted: our CHROM is different than previous line, but we encountered it before
    if (segconf.is_sorted && (prev_line_chrom != NODE_INDEX_NONE) && (prev_line_chrom != vb->chrom_node_index) && 
        *B32 (CTX(CHROM)->counts, vb->chrom_node_index) > 1) // 1 if it has been segged on this line for the first time
        segconf.is_sorted = false;
    
    // evidence of being sorted: same RNAME, increasing POS
    if (prev_line_chrom == vb->chrom_node_index && prev_line_pos <= pos)
        segconf.evidence_of_sorted = true;
}

// mark contexts as used, for calculation of vb_size
void segconf_mark_as_used (VBlockP vb, unsigned num_ctxs, ...)
{
    ASSERTNOTZERO (segconf.running);

    va_list args;
    va_start (args, num_ctxs);

    for (unsigned i=0; i < num_ctxs; i++) {
        Did did_i = (Did)va_arg (args, int);
        CTX(did_i)->local.len32 = 1;
    } 
    
    va_end (args);
}

static void segconf_set_vb_size (VBlockP vb, uint64_t curr_vb_size)
{
    #define VBLOCK_MEMORY_MIN_DYN   (16  MB) // VB memory - min/max when set in segconf_calculate
    #define VBLOCK_MEMORY_MAX_DYN   (512 MB) 
    #define VBLOCK_MEMORY_BEST      (512 MB) // VB memory with --best 
    #define VBLOCK_MEMORY_LOW_MEM   (16  MB) // VB memory - default with --low-memort
    #define VBLOCK_MEMORY_MAKE_REF  (1   MB) // VB memory with --make-reference - reference data 
    #define VBLOCK_MEMORY_GENERIC   (16  MB) // VB memory for the generic data type
    #define VBLOCK_MEMORY_MIN_SMALL (4   MB) // minimum VB memory for small files

    unsigned num_used_contexts=0;

    segconf.vb_size = curr_vb_size;

    if (segconf.vb_size) {
        // already set from previous components of this z_file - do nothing (in particular, FASTQ PAIR_2 must have the same vb_size as PAIR_1)
        // note: for 2nd+ components, we may set other aspects of segconf, but not vb_size
    }

    // if user requested explicit vblock - use it
    else if (flag.vblock) {
        int vblock_len = strlen (flag.vblock);
        
        // case: normal usage - specifying megabytes within the permitted range
        if (vblock_len < 2 || flag.vblock[vblock_len-1] != 'B') { 
            int64_t mem_size_mb;
            ASSINP (str_get_int_range64 (flag.vblock, vblock_len, MIN_VBLOCK_MEMORY, MAX_VBLOCK_MEMORY, &mem_size_mb), 
                    "invalid argument of --vblock: \"%s\". Expecting an integer between 1 and %u. The file will be read and processed in blocks of this number of megabytes.",
                    flag.vblock, MAX_VBLOCK_MEMORY);

            segconf.vb_size = (uint64_t)mem_size_mb MB;
        }

        // case: developer option - a number of bytes eg "100000B"
        else {
            int64_t bytes;
            ASSINP (str_get_int_range64 (flag.vblock, vblock_len-1, ABSOLUTE_MIN_VBLOCK_MEMORY, ABSOLUTE_MAX_VBLOCK_MEMORY, &bytes),  
                    "invalid argument of --vblock: \"%s\". Expecting an integer eg 100000B between %"PRIu64" and %"PRIu64, 
                    flag.vblock, ABSOLUTE_MIN_VBLOCK_MEMORY, ABSOLUTE_MAX_VBLOCK_MEMORY);

            segconf.vb_size = bytes;
        }
    }

    // set memory if --make_reference and user didn't specify --vblock
    else if (flag.make_reference) 
        segconf.vb_size = VBLOCK_MEMORY_MAKE_REF;

    else if (TXT_DT(GENERIC)) 
        segconf.vb_size = VBLOCK_MEMORY_GENERIC;

    // if we failed to calculate an estimated size or file is very small - use default
    else if (txtfile_get_seggable_size() < VBLOCK_MEMORY_GENERIC)
        segconf.vb_size = VBLOCK_MEMORY_GENERIC;

    else { 
        // count number of contexts used
        for_ctx_that (ctx->b250.len || ctx->local.len)
            num_used_contexts++;
            
        uint32_t vcf_samples = (TXT_DT(VCF) || TXT_DT(BCF)) ? vcf_header_get_num_samples() : 0;
        
        // formula - 1MB for each contexts, 128K for each VCF sample
        uint64_t bytes = ((uint64_t)num_used_contexts MB) + (vcf_samples << 17);

        // higher minimum memory for long reads in sorted SAM - enables CPU scaling
        if (segconf.is_long_reads && segconf.is_sorted) {
            uint64_t min_memory = global_max_threads <= 8  ? VBLOCK_MEMORY_MIN_DYN // eg a personal computer
                                : global_max_threads <= 20 ? 64 MB
                                : global_max_threads <= 35 ? 128 MB
                                :                            256 MB;
            // actual memory setting VBLOCK_MEMORY_MIN_DYN to VBLOCK_MEMORY_MAX_DYN
            segconf.vb_size = MIN_(MAX_(bytes, min_memory), VBLOCK_MEMORY_MAX_DYN);
        }

        // larger VBs allow high core scaling due to less contention in merge
        else {
            segconf.vb_size = MIN_(MAX_(bytes, VBLOCK_MEMORY_MIN_DYN), VBLOCK_MEMORY_MAX_DYN);

            if (global_max_threads > 36)
                segconf.vb_size = MIN_(segconf.vb_size * 2, ABSOLUTE_MAX_VBLOCK_MEMORY);
        
            else if (global_max_threads > 20)
                segconf.vb_size = MIN_(segconf.vb_size * 1.5, ABSOLUTE_MAX_VBLOCK_MEMORY);
        }

        int64_t est_seggable_size = txtfile_get_seggable_size();

        // for small files - reduce VB size, to take advantage of all cores (subject to a minimum VB size)
        if (!flag.best && est_seggable_size && global_max_threads > 1) {
            uint64_t new_vb_size = MAX_(VBLOCK_MEMORY_MIN_SMALL, est_seggable_size * 1.2 / global_max_threads);

            // in case of drastic reduction of large-sample VCF file, provide tip
            if (vcf_samples > 10 && new_vb_size < segconf.vb_size / 2) 
                TIP0 ("Using --best can significantly improve compression for this particular file");

            segconf.vb_size = MIN_(segconf.vb_size, new_vb_size);
        } 
        
        // on Windows (inc. WSL) and Mac - which tend to have less memory in typical configurations, warn if we need a lot
        // (note: if user sets --vblock, we won't get here)
        if (flag.is_windows || flag.is_mac || flag.is_wsl) {
            segconf.vb_size = MIN_(segconf.vb_size, 32 MB); // limit to 32MB per VB unless users says otherwise to protect OS UI interactivity 

            int concurrent_vbs = 1 + (est_seggable_size ? MIN_(1+ est_seggable_size / segconf.vb_size, global_max_threads) : global_max_threads);

            if (segconf.vb_size * concurrent_vbs > MEMORY_WARNING_THREASHOLD)
                WARN_ONCE ("\nWARNING: For this file, Genozip selected an optimal setting which consumes a lot of RAM:\n"
                           "%u threads, each processing %u MB of input data at a time (and using working memory too)\n"
                           "To reduce RAM consumption, use --low-memory. To silense this warning use --quiet.\n",
                           global_max_threads, (uint32_t)(segconf.vb_size >> 20));
        }

        segconf.vb_size = ROUNDUP1M (segconf.vb_size);
    }
    
    if (flag.best && !flag.vblock)
        segconf.vb_size = MAX_(segconf.vb_size, VBLOCK_MEMORY_BEST);
    
    if (flag.low_memory && !flag.vblock)
        segconf.vb_size = MIN_(segconf.vb_size, VBLOCK_MEMORY_LOW_MEM);

    if (flag.show_memory && num_used_contexts) 
        iprintf ("\nvblock size set to %u MB (num_used_contexts=%u%s)\n", 
                 (unsigned)(segconf.vb_size >> 20), num_used_contexts, 
                 cond_int (Z_DT(VCF) || Z_DT(BCF), " num_vcf_samples=", vcf_header_get_num_samples()));
}

// this function is called to set is_long_reads, and may be also called while running segconf before is_long_reads is set
bool segconf_is_long_reads(void) 
{ 
    return TECH(PACBIO) || TECH(ONP) || 
           segconf.longest_seq_len > MAX_SHORT_READ_LEN ||
           flag.debug_LONG;
}

static bool segconf_no_calculate (void)
{
    return (Z_DT(FASTQ) && flag.pair == PAIR_R2) || // FASTQ: no recalculating for 2nd pair 
           ((Z_DT(SAM) || Z_DT(BAM)) && flag.deep && flag.zip_comp_i >= SAM_COMP_FQ01); // --deep: no recalculating for second (or more) FASTQ file
}

void segconf_initialize (void)
{
    if (segconf_no_calculate()) return;

    // note: in Deep, we re-segconf the first FASTQ component, but don't remove the SAM segconf data    
    if (!((Z_DT(SAM) || Z_DT(BAM)) && flag.deep && flag.zip_comp_i == SAM_COMP_FQ00))
        segconf = (SegConf){
            .vb_size           = z_file->num_txts_so_far ? segconf.vb_size : 0, // components after first inherit vb_size from first
            .is_sorted         = true,      // initialize optimistically
            
            // SAM stuff
            .is_collated       = true,      // initialize optimistically
            .MAPQ_has_single_value = true,  // initialize optimistically
            .NM_after_MD       = true,      // initialize optimistically
            .nM_after_MD       = true,      // initialize optimistically
            .sam_is_unmapped   = true,      // we will reset this if finding a line with POS>0
            .SA_HtoS           = unknown,
            .sam_XG_inc_S      = unknown,
            .sam_has_BWA_XA_Z  = unknown,
            .CY_con_snip_len   = sizeof (segconf.CY_con_snip),
            .QT_con_snip_len   = sizeof (segconf.QT_con_snip),
            .CB_con_snip_len   = sizeof (segconf.CB_con_snip),

            // FASTQ stuff
            .deep_qtype        = QNONE,

            // FASTA stuff
            .fasta_has_contigs = true, // initialize optimistically
        };

    // Deep - 1st FASTQ file
    else {
        // reset flavor data so we can recalcualte for FASTQ (might be related but different flavor)
        memset (segconf.qname_flavor, 0, sizeof (segconf.qname_flavor));
        memset (segconf.qname_flavor_rediscovered, 0, sizeof (segconf.qname_flavor_rediscovered));
        memset (segconf.qname_line0, 0, sizeof (segconf.qname_line0));
    }
    
    mutex_initialize (segconf.PL_mux_by_DP_mutex);
}

static void segconf_show_has (void)
{
    iprint0 ("Fields recorded in segconf.has:\n");

    bool found = false;
    for_zctx_that (segconf.has[zctx->did_i]) {
        if (TXT_DT(VCF) || TXT_DT(BCF))
            iprintf ("%s/%s ", dict_id_display_type (txt_file->data_type, zctx->dict_id), zctx->tag_name);
        else
            iprintf ("%s ", zctx->tag_name);
        
        found = true;
    }

    iprint0 (found ? "\n" : "None\n");

    exit_ok();
}

// ZIP: Seg a small sample of data of the beginning of the data, to pre-calculate several Seg configuration parameters
void segconf_calculate (void)
{
    // check for components that don't need segconf
    if (segconf_no_calculate()) return; 
    
    if (TXT_DT(GENERIC) ||              // no need for a segconf test VB in generic files
        (flag.biopsy && flag.force)) {  // no need for segconf for a biopsy if --force (allows for biopsy of defective files)
        segconf_set_vb_size (NULL, segconf.vb_size);
        return;
    }

    segconf.running = true;

    uint64_t save_vb_size = segconf.vb_size;
    segconf_vb = vb_initialize_nonpool_vb (VB_ID_SEGCONF, txt_file->data_type, "segconf");
    #define vb segconf_vb

    // note: in case of BZ2, needs to be big enough to overcome the block nature of BZ2 (64K block -> 200-800K text) to get a reasonable size estimate
    uint32_t vb_sizes[] = { 300000, 1500000, 5000000 };
    
    for (int s = (txt_file->codec == CODEC_BZ2); s < ARRAY_LEN(vb_sizes) && !Ltxt; s++) {
        segconf.vb_size = vb_sizes[s];
        txtfile_read_vblock (vb);
        if (txt_file->header_only) break;
    }

    if (!Ltxt) {
        ASSERTW (txt_file->header_only, "FYI: Segconf didn't run because either there is not even one full line, OR first line is larger than %u", vb_sizes[ARRAY_LEN(vb_sizes)-1]);
    
        segconf_set_vb_size (vb, save_vb_size);
        goto done; // cannot find a single line - vb_size set to default and other segconf fields remain default, or previous file's setting
    }
    
    // make a copy of txt_data as seg may modify it
    static Buffer txt_data_copy = {};
    buf_copy (evb, &txt_data_copy, &vb->txt_data, char, 0, 0, "txt_data_copy");

    // segment this VB
    ctx_clone (vb);

    int32_t save_luft_reject_bytes = Z_DT(VCF) ? vcf_vb_get_reject_bytes (vb) : 0;

    SAVE_FLAGS;
    flag.show_alleles = flag.show_digest = flag.show_codec = flag.show_hash = flag.show_reference = false;
    flag.quiet = true;
    flag.show_vblocks = NULL;

    seg_all_data_lines (vb);      
    SAVE_FLAG (aligner_available); // might have been set in sam_seg_finalize_segconf
    RESTORE_FLAGS;
    RESTORE_FLAG (aligner_available);

    if (flag.show_segconf_has) 
        segconf_show_has();

    segconf_set_vb_size (vb, save_vb_size);

    // case: long reads. we can't use the genozip aligner.
    if (segconf.is_long_reads)
        flag.aligner_available = false; // note: original value of flag will be restore when zip of this file is complete (see main_genozip)

    // return the data to txt_file->unconsumed_txt - squeeze it in before the passed-up data
    buf_insert (evb, txt_file->unconsumed_txt, char, 0, txt_data_copy.data, txt_data_copy.len, "txt_file->unconsumed_txt");
    buf_destroy (txt_data_copy);

    if (txt_file->codec == CODEC_BGZF)
        bgzf_return_segconf_blocks (vb); // return BGZF used by the segconf VB to the unconsumed BGZF blocks

    // in case of generated component data - undo
    vb->gencomp_lines.len = 0;

    if (Z_DT(VCF) || Z_DT(BCF))
        txt_file->reject_bytes += save_luft_reject_bytes; // return reject bytes to txt_file, to be reassigned to VB

    // require compressing with a reference when using --best with SAM/BAM/FASTQ, with some exceptions
    bool best_requires_ref = ((VB_DT(SAM) || VB_DT(BAM)) && !(segconf.is_long_reads && segconf.sam_is_unmapped))
                          || (VB_DT(FASTQ) && !segconf.is_long_reads);

    if (segconf.is_long_reads && ((VB_DT(FASTQ) && !flag.deep) || segconf.sam_is_unmapped)) {
        flag.reference = REF_NONE;
        flag.aligner_available = false;
    }

    ASSINP (!flag.best                    // no --best specified
         || flag.force                    // --force overrides
         || !best_requires_ref            // --best doesn't require ref
         || IS_REF_EXTERNAL || IS_REF_EXT_STORE   // reference is provided
         || flag.zip_no_z_file,           // we're not creating a compressed format
            "Using --best on a %s file also requires using --reference. Override with --force.", dt_name(vb->data_type));

done:
    vb_destroy_vb (&vb);

    // restore (used for --optimize-DESC / --add-line-numbers)
    txt_file->num_lines = 0;
    segconf.running = false;
    #undef vb
}

rom segconf_sam_mapper_name (void)
{
    rom mapper_name[] = SAM_MAPPER_NAME;
    ASSERT0 (ARRAY_LEN(mapper_name) == NUM_MAPPERS, "Invalid SAM_MAPPER_NAME array length - perhaps missing commas between strings?");

    return (segconf.sam_mapper >= 0 && segconf.sam_mapper < NUM_MAPPERS) ? mapper_name[segconf.sam_mapper] : "INVALID_MAPPER";
}

rom segconf_tech_name (void)
{
    rom tech_name[] = TECH_NAME;
    ASSERT0 (ARRAY_LEN(tech_name) == NUM_TECHS, "Invalid TECH_NAME array length - perhaps missing commas between strings?");

    return (segconf.tech >= 0 && segconf.tech < NUM_TECHS) ? tech_name[segconf.tech] : "INVALID_MAPPER";
}

