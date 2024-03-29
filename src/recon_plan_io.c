// ------------------------------------------------------------------
//   recon_plan_io.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <errno.h>
#include "genozip.h"
#include "sections.h"
#include "vblock.h"
#include "context.h"
#include "zfile.h"
#include "endianness.h"
#include "file.h"
#include "strings.h"
#include "flags.h"
#include "buffer.h"
#include "dispatcher.h"
#include "piz.h"
#include "compressor.h"
#include "tip.h"

#define MIN_FRAG_LEN_BITS 17 // part of the file format and cannot be changed
#define FRAG_SIZE_BITS 20    // ZIP. fragment size in log2(bytes). valid values: 17-32. transfered to piz via FlagsReconPlan.frag_len_bits. Not part of the file format - can be changed.
#define frag_len(bits) ((1<<(bits)) / sizeof(ReconPlanItem)) // fragment length in units of ReconPlanItem

rom recon_plan_flavors[8] = PLAN_FLAVOR_NAMES;

void recon_plan_show (FileP file, uint32_t conc_writing_vbs, uint32_t vblock_mb)
{
    iprintf ("\nReconstruction plan: entries=%"PRIu64" conc_writing_vbs=%u x %u MB\n", 
             file->recon_plan.len, conc_writing_vbs, vblock_mb);

    for (uint32_t i=0; i < file->recon_plan.len32; i++) {
        ReconPlanItemP p = B(ReconPlanItem, file->recon_plan, i);
        rom comp = p->vb_i ? comp_name(sections_vb_header(p->vb_i)->comp_i) : NULL;
        rom flav = recon_plan_flavors[p->flavor];

        rom no_recon = IS_ZIP                                                                         ? ""
                     : (p->vb_i && !writer_does_vb_need_recon (p->vb_i))                              ? " no_recon"
                     : (p->flavor == PLAN_TXTHEADER && !writer_does_txtheader_need_recon (p->comp_i)) ? " no_recon" 
                     :                                                                                  "";

        rom no_write = IS_ZIP                                                                         ? ""
                     : (p->vb_i && !writer_does_vb_need_write (p->vb_i))                              ? " no_write"
                     : (p->flavor == PLAN_TXTHEADER && !writer_does_txtheader_need_write (p->comp_i)) ? " no_write" 
                     :                                                                                  "";

        switch (p->flavor) {
            case PLAN_RANGE      : iprintf ("%-10s vb=%s/%u\tstart_line=%u\tnum_lines=%u%s%s\n", flav, comp, p->vb_i, p->start_line, p->num_lines, no_recon, no_write); break;
            case PLAN_FULL_VB    : iprintf ("%-10s vb=%s/%u%s%s%s\n", flav, comp, p->vb_i, cond_int (IS_PIZ, "\t\t\tnum_lines=", p->num_lines), no_recon, no_write); break;
            case PLAN_TXTHEADER  : iprintf ("%-10s %s%s%s\n", flav, comp_name(p->comp_i), no_recon, no_write); break;
            case PLAN_REMOVE_ME  : iprintf ("%-10s\n", flav); break;
            case PLAN_DOWNSAMPLE : iprintf ("%-10s vb=%s/%u\tnum_lines=%u%s%s\n", flav, comp, p->vb_i, p->num_lines, no_recon, no_write); break;
            case PLAN_END_OF_VB  : iprintf ("%-10s vb=%s/%u%s%s\n", flav, comp, p->vb_i, no_recon, no_write); break;
            case PLAN_INTERLEAVE : iprintf ("%-10s vb1=%s/%u\tvb2=%s/%u", flav, comp, p->vb_i,
                                            comp_name(sections_vb_header(p->vb2_i)->comp_i), p->vb2_i); 
                                   iprintf (command==PIZ ? "\tnum_lines=%u%s%s\n" : "\n", p->num_lines, no_recon, no_write); // num_lines populated by writer in PIZ
                                   break;
            default              : ABORT ("Unknown flavor %u", p->flavor);
        }
    }               

    if (IS_PIZ && flag.to_stdout)
        TIP0 ("This is the reconstruction plan with flag.no_writer=true because decompressing to stdout is suppressed\n"
              "by --show-plan. To see the recon plan actually used for reconstruction, use --output to send file elsewhere.\n");
}

// -------------------------------------------------------------------------------
// convert ReconPlanItem.start_line between absolute line numbers and deltas
// -------------------------------------------------------------------------------

// ZIP main thread
static void recon_plan_deltify (void)
{
    START_TIMER;

    ARRAY (ReconPlanItem, plan, txt_file->recon_plan);

    VBIType max_vb_i = 0; // Note: this might be more than z_file->num_vbs, as in SAM recon_plan is compressed before PRIM/DEPN components are segged
    for (uint32_t i=0; i < plan_len; i++)
        if (plan[i].vb_i > max_vb_i) max_vb_i = plan[i].vb_i;

    ASSERTNOTINUSE (evb->codec_bufs[0]);
    ARRAY_alloc (uint32_t, next_line, max_vb_i+1, true, evb->codec_bufs[0], evb, "codec_bufs[0]");

    for (uint32_t i=0; i < plan_len; i++)
        if (plan[i].flavor == PLAN_RANGE) {
            uint32_t save_start_line = plan[i].start_line;
            plan[i].start_line -= next_line[plan[i].vb_i];
            next_line[plan[i].vb_i] = save_start_line + plan[i].num_lines;
        }

    buf_free (evb->codec_bufs[0]);

    COPY_TIMER_EVB (recon_plan_deltify);
}

// PIZ main thread
static void recon_plan_dedeltify (void)
{
    ARRAY (ReconPlanItem, plan, evb->scratch);

    // count vbs
    uint32_t max_vb = 0;
    for (uint32_t i=0; i < plan_len; i++)
        max_vb = MAX_(plan[i].vb_i, max_vb);

    ASSERTNOTINUSE (evb->codec_bufs[0]);
    ARRAY_alloc (uint32_t, next_line, max_vb+1, true, evb->codec_bufs[0], evb, "codec_bufs[0]");

    for (uint32_t i=0; i < plan_len; i++)
        if (plan[i].flavor == PLAN_RANGE) {
            plan[i].start_line += next_line[plan[i].vb_i];
            next_line[plan[i].vb_i] = plan[i].start_line + plan[i].num_lines;
        }

    buf_free (evb->codec_bufs[0]);
}

// -------------------------------------
// ZIP: Compress and output recon_plan
// -------------------------------------

// constants used for all fragments
static Codec frag_codec = CODEC_UNKNOWN;
static uint32_t conc_writing_vbs = 0, vblock_mb = 0;

// compress the fragment - either an entire recon_plan, or divide it to fragments if large to allow multi-threaded
// compression and decompression
static void recon_plan_prepare_for_compress (VBlockP vb)
{
    #define FRAG_LEN_ZIP frag_len(FRAG_SIZE_BITS)

    uint32_t frag_i = vb->vblock_i - 1;
    if (frag_i * FRAG_LEN_ZIP >= txt_file->recon_plan.len && vb->vblock_i > 1) return; // don't dispatch (but we do, if its the first VB of an empty recon_plan)

    vb->fragment_start = (char *)B(ReconPlanItem, txt_file->recon_plan, frag_i * FRAG_LEN_ZIP);
    vb->fragment_len   = MIN_(FRAG_LEN_ZIP, txt_file->recon_plan.len - frag_i * FRAG_LEN_ZIP) * sizeof (ReconPlanItem);
    vb->dispatch       = READY_TO_COMPUTE;
}

static void recon_plan_compress_one_fragment (VBlockP vb)
{
    START_TIMER;

    // BGEN
    uint32_t *data32 = (uint32_t *)vb->fragment_start;
    for (uint32_t i=0; i < vb->fragment_len / sizeof (uint32_t); i++, data32++)
        *data32 = BGEN32 (*data32);

    SectionHeaderReconPlan header = (SectionHeaderReconPlan){
        .magic                 = BGEN32 (GENOZIP_MAGIC),
        .section_type          = SEC_RECON_PLAN,
        .data_uncompressed_len = BGEN32 (vb->fragment_len),
        .codec                 = frag_codec,
        .flags.recon_plan.frag_len_bits = FRAG_SIZE_BITS - MIN_FRAG_LEN_BITS, 
        .conc_writing_vbs      = BGEN32 (conc_writing_vbs),
        .vblock_mb             = BGEN32 ((uint32_t)(segconf.vb_size >> 20))
    };

    if (flag.show_time) codec_show_time (vb, st_name(SEC_RECON_PLAN), NULL, frag_codec);

    vb->comp_i = 0; // goes into SectionEnt.comp_i (this is correct for SAM)
    comp_compress (vb, NULL, &vb->z_data, &header, vb->fragment_start, NO_CALLBACK, "SEC_RECON_PLAN");

    COPY_TIMER (recon_plan_compress_one_fragment)    

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}

// called by main thread in zip_write_global_area
void recon_plan_compress (uint32_t my_conc_writing_vbs)
{
    START_TIMER;

    if (flag.show_recon_plan && txt_file->recon_plan.len)
        recon_plan_show (txt_file, my_conc_writing_vbs, (uint32_t)(segconf.vb_size >> 20));

    // replace ReconPlanItem.start_line with delta vs last line of same ReconPlanItem.vb_i
    recon_plan_deltify();

    // assign codec
    // txt_file->recon_plan.len *= sizeof (ReconPlanItem);
    // frag_codec = codec_assign_best_codec (evb, NULL, &txt_file->recon_plan, SEC_RECON_PLAN);
    // txt_file->recon_plan.len /= sizeof (ReconPlanItem);
    
    // better not to do expensive codec assignment in the main thread. Also, CODEC_BSC is always selected 
    // by codec_assign_best_codec, in fast/normal/best modes, no need to test
    frag_codec = CODEC_BSC; 

    conc_writing_vbs = my_conc_writing_vbs;

    // divvy up recon_plan to fragments of about ~1MB and compress in parallel
    dispatcher_fan_out_task ("compress_recon_plan", NULL, 0, "Writing reconstruction plan...", false, false, false, 0, 20000, true,
                             recon_plan_prepare_for_compress, 
                             recon_plan_compress_one_fragment, 
                             zfile_output_processed_vb);

    buf_free (txt_file->recon_plan);

    COPY_TIMER_EVB (recon_plan_compress);
}

// -------------------------------------
// PIZ: Read and decompress recon_plan
// -------------------------------------
static Section next_sec = NULL; 

// PIZ main thread
static void recon_plan_read_one_vb (VBlockP vb)
{
    if (next_sec->st != SEC_RECON_PLAN)
        return; // we're done - no more SEC_RECON_PLAN sections 

    zfile_read_section (z_file, vb, next_sec->vblock_i, &vb->z_data, "z_data", SEC_RECON_PLAN, next_sec);    
    SectionHeaderReconPlanP header =  B1ST(SectionHeaderReconPlan, vb->z_data);

    uint64_t max_frag_size = frag_len(header->flags.recon_plan.frag_len_bits + MIN_FRAG_LEN_BITS) * sizeof (ReconPlanItem); // in bytes

    // first vb - set constants that are the same for all VBs, and allocate memory
    if (vb->vblock_i == 1) {
        conc_writing_vbs = BGEN32 (header->conc_writing_vbs);
        vblock_mb        = BGEN32 (header->vblock_mb);

        ASSERTNOTINUSE (evb->scratch);

        // count fragments
        evb->scratch.count = 0; // number of fragments so far
        for (Section sec=next_sec; sec->st == SEC_RECON_PLAN; sec++) evb->scratch.count++;

        // allocate maximal memory 
        buf_alloc (evb, &evb->scratch, 0, max_frag_size * evb->scratch.count, char, 0, "scratch");
        buf_set_shared (&evb->scratch);
    }

    vb->fragment_len   = BGEN32 (header->data_uncompressed_len);
    vb->fragment_start = Bc (evb->scratch, (vb->vblock_i-1) * max_frag_size);

    // overlay vb->scratch on the part of evb->scratch belonging to this thread for storing UNcompressed data
    buf_overlay_partial (vb, &vb->scratch, &evb->scratch, BNUM(evb->scratch, vb->fragment_start), "scratch");

    evb->scratch.len += vb->fragment_len / sizeof (ReconPlanItem);

    ASSERT (evb->scratch.len * sizeof (ReconPlanItem) <= evb->scratch.size, "%u X recon_plan.len=%"PRIu64" exceeds allocated size=%"PRIu64" (num_fragments=%u max_frag_size=%"PRIu64")", 
            (int)sizeof (ReconPlanItem), evb->scratch.len, (uint64_t)evb->scratch.size, (unsigned)evb->scratch.count, max_frag_size);

    next_sec++;
    vb->dispatch = READY_TO_COMPUTE;
}

// entry point of compute thread of recon_plan decompression
static void recon_plan_uncompress_one_vb (VBlockP vb)
{
    SectionHeaderDictionaryP header = (SectionHeaderDictionaryP)vb->z_data.data;
    uint32_t uncomp_len = BGEN32 (header->data_uncompressed_len); // in bytes

    ASSERT (vb->fragment_start + uncomp_len <= evb->scratch.data + evb->scratch.size, 
            "Buffer overflow when uncompressing recon_plan vb_i=%u", vb->vblock_i);

    zfile_uncompress_section (vb, header, &vb->scratch, NULL, 0, SEC_RECON_PLAN); // NULL name prevents buf_alloc

    vb->scratch.len = uncomp_len / sizeof (uint32_t);
    BGEN_u32_buf (&vb->scratch, 0);

    buf_free (vb->scratch); // un-overlay

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}

// PIZ main thread. Reads, uncompresses and reassembles recon_plan into evb->scratch
void recon_plan_uncompress (Section sec, uint32_t *out_conc_writing_vbs, uint32_t *out_vblock_mb)
{
    ASSERT (sec && sec->st == SEC_RECON_PLAN, "sec=%p is not SEC_RECON_PLAN", sec);

    next_sec = sec;

    dispatcher_fan_out_task ("uncompress_recon_plan", NULL, 0, 0, true, flag.test, false, 0, 1000, true,
                             recon_plan_read_one_vb, recon_plan_uncompress_one_vb, NO_CALLBACK);
 
    // assign outs
    if (out_conc_writing_vbs) *out_conc_writing_vbs = conc_writing_vbs;
    if (out_vblock_mb) *out_vblock_mb = vblock_mb;

    if (VER(14)) // we started to deltify in v14
        recon_plan_dedeltify();
}
