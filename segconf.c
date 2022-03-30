// ------------------------------------------------------------------
//   segconf.c
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "vblock.h"
#include "file.h"
#include "txtfile.h"
#include "seg.h"
#include "segconf.h"
#include "strings.h"
#include "codec.h"
#include "arch.h"

SegConf segconf = {}; // system-wide global

static void segconf_set_vb_size (ConstVBlockP vb, uint64_t curr_vb_size)
{
    #define VBLOCK_MEMORY_MIN_DYN  (16   << 20) // VB memory - min/max when set in segconf_calculate
    #define VBLOCK_MEMORY_MAX_DYN  (512  << 20) 
    #define VBLOCK_MEMORY_BEST     (512  << 20) // VB memory with --best 
    #define VBLOCK_MEMORY_MAKE_REF (1    << 20) // VB memory with --make-reference - reference data 
    #define VBLOCK_MEMORY_GENERIC  (16   << 20) // VB memory for the generic data type
    
    segconf.vb_size = curr_vb_size;

    if (segconf.vb_size) {
        // already set from previous components of this z_file - do nothing (in particular, FASTQ PAIR_2 must have the same vb_size as PAIR_1)
        // note: for 2nd+ components, we may set other aspects of segconf, but not vb_size
    }

    // if user requested explicit vblock - use it
    else if (flag.vblock) {
        int vblock_len = strlen(flag.vblock);
        
        // case: normal usage - specifying megabytes within the permitted range
        if (vblock_len < 2 || flag.vblock[vblock_len-1] != 'B') { 
            int64_t mem_size_mb;
            ASSINP (str_get_int_range64 (flag.vblock, vblock_len, MIN_VBLOCK_MEMORY, MAX_VBLOCK_MEMORY, &mem_size_mb), 
                    "invalid argument of --vblock: \"%s\". Expecting an integer between 1 and %u. The file will be read and processed in blocks of this number of megabytes.",
                    flag.vblock, MAX_VBLOCK_MEMORY);

            segconf.vb_size = (uint64_t)mem_size_mb << 20;
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

    else if (TXT_DT(DT_GENERIC)) 
        segconf.vb_size = VBLOCK_MEMORY_GENERIC;

    // if we failed to calculate an estimated size or file is very small - use default
    else if (txtfile_get_seggable_size() < VBLOCK_MEMORY_GENERIC)
        segconf.vb_size = VBLOCK_MEMORY_GENERIC;

    else { 
        // count number of contexts used
        unsigned num_used_contexts=0;
        for (DidIType did_i=0; did_i < vb->num_contexts ; did_i++)
            if (CTX(did_i)->b250.len || CTX(did_i)->local.len)
                num_used_contexts++;
            
        // formula - 1MB for each contexts, 128K for each VCF sample
        uint64_t bytes = ((uint64_t)num_used_contexts << 20) + 
                            (vcf_header_get_num_samples() << 17 /* 0 if not vcf */);

        uint64_t min_memory = !segconf.sam_is_sorted     ? VBLOCK_MEMORY_MIN_DYN
                            : !segconf.is_long_reads   ? VBLOCK_MEMORY_MIN_DYN
                            : arch_get_num_cores() <= 8  ? VBLOCK_MEMORY_MIN_DYN // eg a personal computer
                            : arch_get_num_cores() <= 20 ? (128 << 20)           // higher minimum memory for long reads in sorted SAM - enables CPU scaling
                            :                              (256 << 20);

        // actual memory setting VBLOCK_MEMORY_MIN_DYN to VBLOCK_MEMORY_MAX_DYN
        segconf.vb_size = MIN_(MAX_(bytes, min_memory), VBLOCK_MEMORY_MAX_DYN);
        
        int64_t est_seggable_size = txtfile_get_seggable_size();
        if (est_seggable_size) segconf.vb_size = MIN_(segconf.vb_size, est_seggable_size * 1.5);

        if (flag.show_memory)
            iprintf ("\nDyamically set vblock_memory to %u MB (num_contexts=%u num_vcf_samples=%u)\n", 
                        (unsigned)(segconf.vb_size >> 20), num_used_contexts, vcf_header_get_num_samples());

        // on Windows (inc. WSL2) and Mac - which tend to have less memory in typical configurations, warn if we need a lot
        // (note: if user sets --vblock, we won't get here)
        if (flag.is_windows || flag.is_mac || strstr (arch_get_os(), "microsoft-standard") /* WSL2 */) {
            segconf.vb_size = MIN_(segconf.vb_size, 32 << 20); // limit to 32MB per VB unless users says otherwise to protect OS UI interactivity 

            uint64_t concurrent_vbs = 1 + (txt_file->disk_size ? MIN_(1+ txt_file->disk_size / segconf.vb_size, global_max_threads)
                                                               : global_max_threads);

            ASSERTW (segconf.vb_size * concurrent_vbs < MEMORY_WARNING_THREASHOLD,
                    "\nWARNING: For this file, Genozip selected an optimal setting which consumes a lot of RAM:\n"
                    "%u threads, each processing %u MB of input data at a time (and using working memory too)\n"
                    "To reduce RAM consumption, you may use:\n"
                    "   --threads to set the number of threads (affects speed)\n"
                    "   --vblock to set the amount of input data (in MB) a thread processes (affects compression ratio)\n"
                    "   --quiet to silence this warning",
                    global_max_threads, (uint32_t)(segconf.vb_size >> 20));
        }
    }
    
    if (flag.best && !flag.vblock)
        segconf.vb_size = MAX_(segconf.vb_size, VBLOCK_MEMORY_BEST);
}

// this function is called to set is_long_reads, and may be also called while running segconf before is_long_reads is set
bool segconf_is_long_reads(void) 
{ 
    return segconf.tech == TECH_PACBIO    || 
           segconf.tech == TECH_ONP       || 
           segconf.longest_seq_len > 2000 ||
           flag.debug_LONG;
}


static bool segconf_no_calculate (void)
{
    return (Z_DT(DT_FASTQ) && flag.pair == PAIR_READ_2); // FASTQ: no recalculating for 2nd pair 
}

void segconf_initialize (void)
{
    if (segconf_no_calculate()) return;

    uint64_t save_vb_size = segconf.vb_size;
    segconf = (SegConf){}; // reset for new component

    if (z_file->num_txts_so_far > 0)
        segconf.vb_size = save_vb_size; // components after first inherit vb_size from first

    mutex_initialize (segconf.PL_mux_by_DP_mutex);
}

// ZIP: Seg a small sample of data of the beginning of the data, to pre-calculate several Seg configuration parameters
void segconf_calculate (void)
{
    if (segconf_no_calculate()) return;

    if (TXT_DT(DT_GENERIC)) {                                     // nothing to calculate in generic files    
        segconf.vb_size = VBLOCK_MEMORY_GENERIC;
        return;
    }

    segconf.running = true;

    uint64_t save_vb_size = segconf.vb_size;
    VBlockP vb = vb_initialize_nonpool_vb (VB_ID_SEGCONF, z_file->data_type, "segconf");

    // note: in case of BZ2, needs to be big enough to overcome the block nature of BZ2 (64K block -> 200-800K text) to get a reasonable size estimate
    uint32_t vb_sizes[] = { 300000, 1500000, 5000000 };
    
    for (int s = (txt_file->codec == CODEC_BZ2); s < ARRAY_LEN(vb_sizes) && !vb->txt_data.len; s++) {
        segconf.vb_size = vb_sizes[s];
        txtfile_read_vblock (vb);
        if (txt_file->header_only) break;
    }

    if (!vb->txt_data.len) {
        ASSERTW (txt_file->header_only, "Segconf didn't run because first line is larger than %u", vb_sizes[ARRAY_LEN(vb_sizes)-1]);
    
        segconf_set_vb_size (vb, save_vb_size);
        goto done; // cannot find a single line - vb_size set to default and other segconf fields remain default, or previous file's setting
    }
    
    // make a copy of txt_data as seg may modify it
    static Buffer txt_data_copy = {};
    buf_copy (evb, &txt_data_copy, &vb->txt_data, char, 0, 0, "txt_data_copy");

    // segment this VB
    ctx_clone (vb);

    int32_t save_luft_reject_bytes = Z_DT(DT_VCF) ? vcf_vb_get_reject_bytes (vb) : 0;

    SAVE_FLAGS;
    flag.show_alleles = flag.show_digest = flag.show_codec = flag.show_hash =
    flag.show_reference = flag.show_vblocks = false;
    flag.quiet = true;

    seg_all_data_lines (vb);        
    RESTORE_FLAGS;

    segconf_set_vb_size (vb, save_vb_size);
    segconf.is_long_reads = segconf_is_long_reads();

    // return the data to txt_file->unconsumed_txt - squeeze it in before the passed-up data
    buf_alloc (evb, &txt_file->unconsumed_txt, txt_data_copy.len, 0, char, 0, "txt_file->unconsumed_txt");
    memmove (&txt_file->unconsumed_txt.data[txt_data_copy.len], txt_file->unconsumed_txt.data, txt_file->unconsumed_txt.len);
    memcpy (txt_file->unconsumed_txt.data, txt_data_copy.data, txt_data_copy.len);
    txt_file->unconsumed_txt.len += txt_data_copy.len;
    buf_destroy (txt_data_copy);

    // in case of generated component data - undo
    vb->gencomp_lines.len = 0;

    if (Z_DT(DT_VCF))
        txt_file->reject_bytes += save_luft_reject_bytes; // return reject bytes to txt_file, to be reassigned to VB

    // require compressing with a reference when using --best with SAM/BAM/FASTQ, with some exceptions
    bool best_requires_ref = ((VB_DT(DT_SAM) || VB_DT(DT_BAM)) && !(segconf.is_long_reads && segconf.sam_is_unmapped))
                          || (VB_DT(DT_FASTQ) && !segconf.is_long_reads);

    if (segconf.is_long_reads && (VB_DT(DT_FASTQ) || segconf.sam_is_unmapped)) {
        flag.reference = REF_NONE;
        flag.aligner_available = false;
    }

    ASSINP (flag.best != 1                // no --best specified, or specified --best=NO_REF
         || !best_requires_ref            // --best doesn't require ref
         || flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE // reference is provided
         || flag.seg_only || flag.biopsy, // we're not creating a compressed format
            "Using --best on a %s file also requires using --reference. You can override this requirement with --best=NO_REF", dt_name(vb->data_type));

done:
    vb_destroy_vb (&vb);

    // restore (used for --optimize-DESC / --add-line-numbers)
    txt_file->num_lines = 0;
    segconf.running = false;
}

void segconf_update_qual (STRp (qual))
{
    if (segconf.nontrivial_qual) return; // already established

    for (uint32_t i=1; i < qual_len; i++) 
        if (qual[i] != qual[0]) {
            segconf.nontrivial_qual = true;
            break;
        }
}
