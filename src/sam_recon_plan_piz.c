// ------------------------------------------------------------------
//   sam_recon_plan_piz.c
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <errno.h>
#include "zfile.h"
#include "buffer.h"
#include "dispatcher.h"
#include "compressor.h"
#include "writer_private.h"
#include "file.h"
#include "tip.h"
#include "strings.h"

#define MIN_FRAG_LEN_BITS 17 // part of the file format and cannot be changed
#define FRAG_SIZE_BITS 20    // ZIP. fragment size in log2(bytes). valid values: 17-32. transfered to piz via FlagsReconPlan.frag_len_bits. Not part of the file format - can be changed.
#define frag_len(bits) ((1<<(bits)) / sizeof(ReconPlanItem)) // fragment length in units of ReconPlanItem

StrText display_plan_item (ReconPlanItemP p)
{
    StrText s;
    
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

    VbInfo *v = (IS_PIZ && p->vb_i) ? VBINFO(p->vb_i) : NULL;
    uint32_t num_dropped_lines = (v && v->is_dropped) ? bits_num_set_bits (v->is_dropped) : 0;
    rom dropped_lines = cond_int (num_dropped_lines, " dropped_lines=", num_dropped_lines);

    switch (p->flavor) {
        case PLAN_VB_PLAN    : snprintf (s.s, sizeof(s), "%-10s vb=%.12s/%u%.10s", flav, comp, p->vb_i, cond_int (IS_PIZ, "\tnum_lines=", p->num_lines)); break; // note: num_lines including integrated PRIM and DEPN lines
        case PLAN_RANGE      : snprintf (s.s, sizeof(s), "%-10s vb=%.12s/%u\tstart_line=%u\tnum_lines=%u%s%s", flav, comp, p->vb_i, p->start_line, p->num_lines, no_recon, no_write); break;
        case PLAN_FULL_VB    : snprintf (s.s, sizeof(s), "%-10s vb=%.12s/%u%.10s%.10s%s%s", flav, comp, p->vb_i, cond_int (IS_PIZ, "\tnum_lines=", p->num_lines), dropped_lines, no_recon, no_write); break;
        case PLAN_TXTHEADER  : snprintf (s.s, sizeof(s), "%-10s %.12s%s%s", flav, comp_name(p->comp_i), no_recon, no_write); break;
        case PLAN_REMOVE_ME  : snprintf (s.s, sizeof(s), "%-10s vb=%.12s/%u\tnum_lines=%u", flav, comp, p->vb_i, p->num_lines); break;
        case PLAN_DOWNSAMPLE : snprintf (s.s, sizeof(s), "%-10s vb=%.12s/%u\tnum_lines=%u%s%s", flav, comp, p->vb_i, p->num_lines, no_recon, no_write); break;
        case PLAN_END_OF_VB  : snprintf (s.s, sizeof(s), "%-10s vb=%.12s/%u%.10s%s%s", flav, comp, p->vb_i, dropped_lines, no_recon, no_write); break;
        case PLAN_INTERLEAVE : snprintf (s.s, sizeof(s), "%-10s vb1=%.8s/%u\tvb2=%.8s/%u%.10s%.10s%.9s%.9s", flav, comp, p->vb_i,
                                         comp_name(sections_vb_header(p->vb2_i)->comp_i), p->vb2_i,
                                         cond_int (IS_PIZ, "\tnum_lines=", p->num_lines), dropped_lines, no_recon, no_write); break;
        default              : ABORT ("Unknown flavor %u", p->flavor);
    }  

    return s;
}

void recon_plan_show (int vb_i, // -1 means for entire file, 0 means for txt header
                      int64_t start, 
                      int64_t len) // -1 = until plan's end
{
    iprintf ("\nReconstruction plan: %s%s entries=%"PRIu64" conc_writing_vbs=%u vb_1_size=%s\n", 
             cond_int (vb_i >= 1, "for MAIN vblock_i=", vb_i),
             vb_i==0 ? "for TXT_HEADER" : "",
             z_file->recon_plan.len, z_file->max_conc_writing_vbs, str_size (segconf.vb_size).s);

    if (len == -1) 
        len = z_file->recon_plan.len - start;
    
    for (int64_t i=start; i < start + len; i++) 
        iprintf ("%s\n", display_plan_item (B(ReconPlanItem, z_file->recon_plan, i)).s);

    if (IS_PIZ && flag.to_stdout)
        TIP0 ("This is the reconstruction plan with flag.no_writer=true because decompressing to stdout is suppressed\n"
              "by --show-plan. To see the recon plan actually used for reconstruction, use --output to send file elsewhere.\n");
}

// ---------------------------------------------------------------------------------------------------------------
// PIZ: Read and decompress recon_plan (used for DVCF v12-15.0.42, and for SAM gencomp v14-15.0.63)
// ---------------------------------------------------------------------------------------------------------------

static Section next_sec = NULL; 

// PIZ main thread: convert ReconPlanItem.start_line: absolute-line-numbers ⇔ deltas
static void recon_plan_dedeltify (void)
{
    ARRAY (ReconPlanItem, plan, z_file->recon_plan);

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
        z_file->max_conc_writing_vbs = BGEN32 (header->conc_writing_vbs);
        segconf.vb_size = BGEN32 (header->vblock_mb) MB; // VB size of vb=1 as proxy for a typical VB size. Note that other VBs might be larger or smaller.

        // add memory to recon_plan
        buf_alloc (wvb, &z_file->recon_plan, max_frag_size * sections_get_recon_plan(NULL), 0, char, 0, "z_file->recon_plan");
        buf_set_shared (&z_file->recon_plan);
    }

    vb->fragment_len   = BGEN32 (header->data_uncompressed_len);
    vb->fragment_start = Bc (z_file->recon_plan, (vb->vblock_i-1) * max_frag_size + sizeof (ReconPlanItem)/*PLAN_TXTHEADER*/);

    // overlay vb->scratch on the part of z_file->recon_plan belonging to this thread for storing UNcompressed data
    buf_overlay_partial (vb, &vb->scratch, &z_file->recon_plan, BNUM(z_file->recon_plan, vb->fragment_start), "scratch");

    z_file->recon_plan.len += vb->fragment_len / sizeof (ReconPlanItem);

    ASSERT (z_file->recon_plan.len * sizeof (ReconPlanItem) <= z_file->recon_plan.size, "%u X recon_plan.len=%"PRIu64" exceeds allocated size=%"PRIu64" (num_fragments=%u max_frag_size=%"PRIu64")", 
            (int)sizeof (ReconPlanItem), z_file->recon_plan.len, (uint64_t)z_file->recon_plan.size, sections_get_recon_plan (NULL), max_frag_size);

    next_sec++;
    vb->dispatch = READY_TO_COMPUTE;
}

// entry point of compute thread of recon_plan decompression
static void recon_plan_uncompress_one_vb (VBlockP vb)
{
    SectionHeaderDictionaryP header = (SectionHeaderDictionaryP)vb->z_data.data;
    uint32_t uncomp_len = BGEN32 (header->data_uncompressed_len); // in bytes

    ASSERT (vb->fragment_start + uncomp_len <= z_file->recon_plan.data + z_file->recon_plan.size, 
            "Buffer overflow when uncompressing recon_plan vb_i=%u", vb->vblock_i);

    zfile_uncompress_section (vb, header, &vb->scratch, NULL, 0, SEC_RECON_PLAN); // NULL name prevents buf_alloc

    vb->scratch.len = uncomp_len / sizeof (uint32_t);
    BGEN_u32_buf (&vb->scratch, 0);

    buf_destroy (vb->scratch); // un-overlay

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}

// PIZ main thread. 
// case: v14-15.0.63 (and since v12 for DVCF): a SAM/BAM with gencomp had a SEC_RECON_PLAN section:
// Reads, uncompresses and reassembles recon_plan appending z_file->recon_plan
void recon_plan_add_prescribed_by_recon_plan_section (void)
{
    Section recon_plan_sec;
    if (!sections_get_recon_plan (&recon_plan_sec)) {
        writer_add_trivial_plan (SAM_COMP_MAIN, PLAN_FULL_VB); // SEC_RECON_PLAN was not written if file is only full MAIN VBs
        return;
    }

    next_sec = recon_plan_sec;

    dispatcher_fan_out_task ("uncompress_recon_plan", NULL, 0, 0, true, flag.test, false, 0, 1000, true,
                            recon_plan_read_one_vb, recon_plan_uncompress_one_vb, NO_CALLBACK);

    recon_plan_dedeltify();

    // add num_lines that is left out in the file to improve compression of the SEC_RECON_PLAN section
    // and keep only VBs that need reconstruction
    ARRAY (ReconPlanItem, plan, z_file->recon_plan);
    ReconPlanItem *dst = plan + 1; // skip PLAN_TXTHEADER

    // keep only VBs that needs_recon
    for (uint32_t i=1; i < plan_len; i++) { // 1 to skip PLAN_TXTHEADER
        VbInfo *v = VBINFO(plan[i].vb_i);

        if (v->needs_recon) {
            if (plan[i].flavor == PLAN_FULL_VB) 
                plan[i].num_lines = v->num_lines; // note: in the file format num_lines=0 - improves compression of RECON_PLAN section

            *dst++ = plan[i];
        }
    }

    z_file->recon_plan.len = dst - plan;
}

