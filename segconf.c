// ------------------------------------------------------------------
//   segconf.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "vblock.h"
#include "file.h"
#include "txtfile.h"
#include "seg.h"
#include "segconf.h"
#include "strings.h"
#include "codec.h"

SegConf segconf = {}; // system-wide global

static void segconf_set_vb_size (const VBlock *vb)
{
    #define MAX_VBLOCK_MEMORY      2048         // in MB
    #define VBLOCK_MEMORY_MIN_DYN  (16   << 20) // VB memory - min/max when set in segconf_calculate
    #define VBLOCK_MEMORY_MAX_DYN  (512  << 20) 
    #define VBLOCK_MEMORY_FAST     (16   << 20) // VB memory with --fast
    #define VBLOCK_MEMORY_MAKE_REF (1    << 20) // VB memory with --make-reference - reference data 
    #define VBLOCK_MEMORY_GENERIC  (16   << 20) // VB memory for the generic data type

    if (segconf.vb_size) {
        // already set from previous components of this z_file - do nothing (in particular, FASTQ PAIR_2 must have the same vb_size as PAIR_1)
    }

    // if user requested explicit vblock - use it
    else if (flag.vblock) {
        int64_t mem_size_mb;
        ASSINP (str_get_int_range64 (flag.vblock, 0, 1, MAX_VBLOCK_MEMORY, &mem_size_mb), 
                "invalid argument of --vblock: \"%s\". Expecting an integer between 1 and %u. The file will be read and processed in blocks of this number of megabytes.",
                flag.vblock, MAX_VBLOCK_MEMORY);

        segconf.vb_size = (uint64_t)mem_size_mb << 20;
    }

    // set memory if --fast and user didn't specify --vblock
    else if (flag.fast) 
        segconf.vb_size = VBLOCK_MEMORY_FAST;

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

        // actual memory setting VBLOCK_MEMORY_MIN_DYN to VBLOCK_MEMORY_MAX_DYN
        segconf.vb_size = MIN_(MAX_(bytes, VBLOCK_MEMORY_MIN_DYN), VBLOCK_MEMORY_MAX_DYN);
        
        int64_t est_seggable_size = txtfile_get_seggable_size();
        if (est_seggable_size) segconf.vb_size = MIN_(segconf.vb_size, est_seggable_size * 1.5);

        if (flag.show_memory)
            iprintf ("\nDyamically set vblock_memory to %u MB (num_contexts=%u num_vcf_samples=%u)\n", 
                        (unsigned)(segconf.vb_size >> 20), num_used_contexts, vcf_header_get_num_samples());

        // on Windows and Mac - which tend to have less memory in typical configurations, warn if we need a lot
        #if defined _WIN32 || defined APPLE
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
        #endif
    }
}

// ZIP: Seg a small sample of data of the beginning of the data, to pre-calculate several Seg configuration parameters
void segconf_calculate (void)
{
    if (flag.pair && z_file->num_txt_components_so_far % 2 == 0)
        return; // in FASTQ pairs, we only calculate for the first one

    if (z_file->num_txt_components_so_far == 1)
        segconf = (SegConf){}; // reset for new z_file

    segconf.running = true;

    uint64_t save_vb_size = segconf.vb_size;
    segconf.vb_size = txt_file->codec == CODEC_BZ2 ? 1500000 // needs to be big enough to overcome the block nature of BZ2 (64K block -> 200-800K text) to get a reasonable size estimate
                    :                                300000; // big enough of a long-read SAM or FASTQ line
    
    VBlockP vb = vb_initialize_nonpool_vb (VB_ID_SEGCONF, z_file->data_type, "segconf");

    txtfile_read_vblock (vb);
    if (!vb->txt_data.len) {
        segconf_set_vb_size (vb);
        goto done; // cannot find a single line - vb_size set to default and other segconf fields remain default, or previous file's setting
    }
    
    // make a copy of txt_data as seg may modify it
    static Buffer txt_data_copy = {};
    buf_copy (evb, &txt_data_copy, &vb->txt_data, char, 0, 0, "txt_data_copy");

    // segment this VB
    ctx_clone (vb);

    int32_t save_luft_reject_bytes = vb->reject_bytes;

    SAVE_FLAGS;
    flag.show_alleles = flag.show_digest = flag.show_codec = flag.show_hash =
    flag.show_reference = flag.show_vblocks = false;
    flag.quiet = true;

    seg_all_data_lines (vb);        

    RESTORE_FLAGS;

    segconf.vb_size = save_vb_size;
    segconf_set_vb_size (vb);
    
    // return the data to txt_file->unconsumed_txt - squeeze it in before the passed-up data
    buf_alloc (evb, &txt_file->unconsumed_txt, txt_data_copy.len, 0, char, 0, "txt_file->unconsumed_txt");
    memmove (&txt_file->unconsumed_txt.data[txt_data_copy.len], txt_file->unconsumed_txt.data, txt_file->unconsumed_txt.len);
    memcpy (txt_file->unconsumed_txt.data, txt_data_copy.data, txt_data_copy.len);
    txt_file->unconsumed_txt.len += txt_data_copy.len;
    buf_destroy (&txt_data_copy);

    // in case of Luft file - undo
    vb->lo_rejects[0].len = vb->lo_rejects[1].len = 0;
    txt_file->reject_bytes += save_luft_reject_bytes; // return reject bytes to txt_file, to be reassigned to VB

done:
    vb_destroy_vb (&vb);

    // restore (used for --optimize-DESC / --add-line-numbers)
    txt_file->num_lines = 0;
    segconf.running = false;
}
