// ------------------------------------------------------------------
//   sam_gc.c
//   Copyright (C) 2022-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "sam_private.h"
#include "sections.h"
#include "codec.h"
#include "compressor.h"

typedef struct __attribute__ ((__packed__)) {
    uint32_t line_i;
    uint32_t line_len : 31;
    uint32_t is_prim  : 1;
} GencompLineIEntry;

typedef struct {
    uint32_t vb_i;
    uint32_t num_lines;
    uint64_t txt_len;
} SamGcVbInfo;

#define SAM_GC_UPDATE_PRIM 0xfffffffe
#define SAM_GC_UPDATE_DEPN 0xffffffff

// returns true if this line is a Primary or Dependent line of a supplementary/secondary group - 
// and should be moved to a generated component. Called by compute thread in seg of Normal VB.
bool sam_seg_is_gc_line (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(alignment), STRps(aux), bool is_bam)
{
    STR(NH); STR(HI); 
    SamComponentType gc_type = 0;;

    // dependent - always has secondary or supplementary flag
    if (dl->FLAG.bits.secondary || dl->FLAG.bits.supplementary) 
        gc_type = CT_SA_DEPN;

    // primary - no special flag, identify by SA or (NH>1 && HI==1)
    else if (sam_seg_get_aux ("SA:Z", STRas(aux), NULL, is_bam) ||
               (((NH = sam_seg_get_aux ("NH:i", STRas(aux), &NH_len, is_bam)) && (NH_len > 1 || *NH != '1')) &&
               (((HI = sam_seg_get_aux ("HI:i", STRas(aux), &HI_len, is_bam)) && (HI_len ==1 || *HI == '1')))))
        gc_type = CT_SA_PRIM;

    else
        return false;

    buf_add_more (VB, &vb->gencomp[gc_type-1], STRa(alignment), "gencomp");

    // store location where this gc line should be inserted    
    buf_alloc (VB, &vb->gc_line_info, 1, 100, GencompLineIEntry, 2, "gc_line_info");
    NEXTENT (GencompLineIEntry, vb->gc_line_info) = (GencompLineIEntry){ .line_i   = vb->line_i, 
                                                                         .line_len = alignment_len,
                                                                         .is_prim  = (gc_type == CT_SA_PRIM) };
    vb->line_i--;
    vb->recon_size -= alignment_len;

    return true;
}

// called main thread after_compute, in order of VBs
void sam_zip_gc_after_compute (VBlockSAMP vb)
{
    if (flag.gencomp_num) return;

    BufferP recon_plan = &txt_file->recon_plan;
    ARRAY (GencompLineIEntry, gc_lines, vb->gc_line_info);

    // we create a "Full VB" recon plan if there are no lines in a Supplementary/Secondary group in this
    // VB. Other VBs might have them. If no VB has them, we will get rid of the recon_plan in sam_zip_compress_recon_plan
    if (!gc_lines_len) {
        buf_alloc (evb, recon_plan, 1, 100000, ReconPlanItem, 2, "txt_file->recon_plan");

        NEXTENT (ReconPlanItem, *recon_plan) = (ReconPlanItem){
            .full_vb.plan_type = PLAN_FULL_VB,
            .full_vb.vb_i      = vb->vblock_i
        };
    }

    // case: we removed some gencomp lines from this VB, we now create a recon plan to insert them back.
    // the actual line location in the gencomp vb will be updated later
    else {
        // write flush gencomp lines to the gencomp files
        if (vb->gencomp[CT_SA_PRIM-1].len) gencomp_append_file (VB, CT_SA_PRIM, "PRIM");
        if (vb->gencomp[CT_SA_DEPN-1].len) gencomp_append_file (VB, CT_SA_DEPN, "DEPN");

        buf_alloc (evb, recon_plan, gc_lines_len * 2 + 2, 100000, ReconPlanItem, 2, "txt_file->recon_plan");
        buf_alloc (evb, &txt_file->line_info[0], gc_lines_len, 100000, uint32_t, 2, "txt_file->line_info");
        buf_alloc (evb, &txt_file->line_info[1], gc_lines_len, 100000, uint32_t, 2, "txt_file->line_info");

        uint32_t normal_line_i=0;
        for (uint32_t gc_i=0 ; gc_i < gc_lines_len; gc_i++) {
            GencompLineIEntry gc_line_i = gc_lines[gc_i];
            
            // insert normal lines before the next gc line
            if (gc_line_i.line_i > normal_line_i) {

                NEXTENT (ReconPlanItem, *recon_plan) = (ReconPlanItem){
                    .range.vb_i       = vb->vblock_i,
                    .range.start_line = normal_line_i,
                    .range.num_lines  = gc_line_i.line_i - normal_line_i
                };

                normal_line_i = gc_line_i.line_i; // next normal line will be after this gc line and possibly subsequent ones that are marked as coming before it
            }

            // insert gc lines - a bunch of them that have the same component and are consecutive
            // vb_i within the gencomp components and start_line will be updated later as we don't know them yet
            NEXTENT (ReconPlanItem, *recon_plan) = 
                (ReconPlanItem){ .range.vb_i = gc_line_i.is_prim ? SAM_GC_UPDATE_PRIM : SAM_GC_UPDATE_DEPN };
        
            // store line lengths, to be used later to calculate vb_info
            NEXTENT (uint32_t, txt_file->line_info[!gc_line_i.is_prim]) = gc_line_i.line_len;
        }
        
        // insert final normal lines
        if (normal_line_i < vb->lines.len)
            NEXTENT (ReconPlanItem, *recon_plan) = (ReconPlanItem){
                .range.vb_i       = vb->vblock_i,
                .range.start_line = normal_line_i,
                .range.num_lines  = vb->lines.len - normal_line_i
            };

        // insert end-of-VB
        NEXTENT (ReconPlanItem, txt_file->recon_plan) = (ReconPlanItem){ 
            .end_of_vb.vb_i      = vb->vblock_i,
            .end_of_vb.plan_type = PLAN_END_OF_VB,
        }; 
    }
}

// callback function of compress_recon_plan, called from zip_one_file
static bool sam_zip_recon_plan_full_vb_only (void)
{
    ARRAY (ReconPlanItem, recon_plan, txt_file->recon_plan);

    for (uint32_t i=0; i < recon_plan_len; i++)
        if (recon_plan[i].x.plan_type != PLAN_FULL_VB) return false;

    return true;
}

// predict (precisely!) how txtfile will split the vb->gencomp data into VBs
static void sam_zip_gc_calc_vb_info (void)
{
    uint64_t max_memory_per_vb = txtfile_max_memory_per_vb();
    Buffer *line_info = txt_file->line_info;
    uint32_t next_vb_i = txt_file->num_vbs + 1;

    for (int gc_i=0; gc_i < 2; gc_i++) {

        uint64_t txt_len = 0;
        uint32_t num_lines=0;
        
        for (uint32_t i=0; i < line_info[gc_i].len; i++) {
            uint32_t line_len = *ENT (uint32_t, line_info[gc_i], i);
            
            if (txt_len + line_len > max_memory_per_vb) {
                buf_alloc (evb, &z_file->vb_info[gc_i], 1, 20, SamGcVbInfo, 2, "z_file->vb_info");
                NEXTENT (SamGcVbInfo, z_file->vb_info[gc_i]) = 
                    (SamGcVbInfo) {.vb_i = next_vb_i, .num_lines = num_lines, .txt_len = txt_len };

                num_lines = 0;
                next_vb_i++;
                txt_len = 0;
            }

            txt_len += line_len;
            num_lines++;
        }

        if (txt_len) {
            buf_alloc (evb, &z_file->vb_info[gc_i], 1, 20, SamGcVbInfo, 2, "z_file->vb_info");
            NEXTENT (SamGcVbInfo, z_file->vb_info[gc_i]) = 
                (SamGcVbInfo) {.vb_i = next_vb_i, .num_lines = num_lines, .txt_len = txt_len };
            next_vb_i++;
        }
    }
} 

// called by main thread after reading a VB - callback of zip_read_one_vb
void sam_zip_verify_gc_vb (VBlockP vb)
{
    if (!flag.gencomp_num || !vb->txt_data.len) return;

    ARRAY (SamGcVbInfo, vb_info, z_file->vb_info[flag.gencomp_num-1]);

    ASSERT (vb_info_len, "No vb_info for gencomp_num=%u, while verifying vb=%u", flag.gencomp_num, vb->vblock_i);

    int i = vb->vblock_i - vb_info[0].vb_i; 

    ASSERT (i < vb_info_len, "Cannot find vb_i=%u of gencomp_num=%u in vb_info", vb->vblock_i, flag.gencomp_num);

    ASSERT (vb_info[i].txt_len == vb->txt_data.len, "Expected vb_i=%u of gencomp_num=%u to have txt_len=%"PRIu64" but it has txt_data.len=%"PRIu64,
            vb->vblock_i, flag.gencomp_num, vb_info[i].txt_len, vb->txt_data.len);
}

static void sam_zip_recon_plan_add_gc_lines (void)
{
    sam_zip_gc_calc_vb_info();

    SamGcVbInfo *vb_info[2] = { FIRSTENT (SamGcVbInfo, z_file->vb_info[CT_SA_PRIM-1]), 
                                FIRSTENT (SamGcVbInfo, z_file->vb_info[CT_SA_DEPN-1]) };
    uint32_t gc_vb_line_i[2] = {};
    
    for (uint32_t i=0; i < txt_file->recon_plan.len; i++) {
        ReconPlanItem *pi = ENT (ReconPlanItem, txt_file->recon_plan, i);

        if (pi->x.plan_type >= MIN_PLAN_TYPE || pi->range.vb_i < SAM_GC_UPDATE_PRIM) continue; // not a gc item

        bool is_depn = (pi->range.vb_i == SAM_GC_UPDATE_DEPN);

        ASSERTNOTNULL (vb_info[is_depn]);
        ASSERT (vb_info[is_depn] < AFTERENT (SamGcVbInfo, z_file->vb_info[is_depn]),
                "vb_info[%u] out if bounds", is_depn);
                
        *pi = (ReconPlanItem) {
            .range.vb_i       = vb_info[is_depn]->vb_i,
            .range.start_line = gc_vb_line_i[is_depn]++,
            .range.num_lines  = 1
        };

        vb_info[is_depn]->num_lines--;
        if (!vb_info[is_depn]->num_lines) {
            
            // add PLAN_END_OF_VB after last line of the VB
            *INSERTENTAFTER (ReconPlanItem, txt_file->recon_plan, i) = (ReconPlanItem){
                .end_of_vb.plan_type = PLAN_END_OF_VB,
                .end_of_vb.vb_i = vb_info[is_depn]->vb_i,
            };
            
            i++; // skip inserted PLAN_END_OF_VB
            vb_info[is_depn]++;
            gc_vb_line_i[is_depn] = 0;
        }
    }
}

// merge consecutive recon plan entries of the same VB to reduce the recon plan size
static void sam_zip_recon_plan_optimize_gc_lines (void)
{
    ARRAY (ReconPlanItem, recon_plan, txt_file->recon_plan);

    uint32_t last_item = 0;
    for (uint32_t i=1; i < recon_plan_len; i++) 
        if (recon_plan[i]        .range.num_lines < MIN_PLAN_TYPE
         && recon_plan[last_item].range.num_lines < MIN_PLAN_TYPE
         && recon_plan[i].range.vb_i == recon_plan[last_item].range.vb_i) 
            recon_plan[last_item].range.num_lines++;
        
        else 
            recon_plan[++last_item] = recon_plan[i];

    txt_file->recon_plan.len = last_item + 1;
}

void sam_zip_compress_recon_plan (void)
{
    // note: the normal component has a recon plan that pulls for generated components too, 
    // genenerated componets have an empty plan - causing writer to not write them despite being reconstructed
    if (!flag.gencomp_num) {

        // if we have no gencomp components, we don't need a recon plan
        if (sam_zip_recon_plan_full_vb_only()) return;

        // update recon plan with gencomp lines
        sam_zip_recon_plan_add_gc_lines();

        // merge consecutive GC lines from the same GC & VB
        sam_zip_recon_plan_optimize_gc_lines();

        if (flag.show_recon_plan)
            linesorter_show_recon_plan (txt_file, false, 3, (uint32_t)(segconf.vb_size >> 20));
    }

    // get best codec for the reconstruction plan data
    Codec codec = codec_assign_best_codec (evb, 0, &txt_file->recon_plan, SEC_RECON_PLAN);
    if (codec == CODEC_UNKNOWN) codec = CODEC_NONE; // too small for compression

    // prepare section header and compress
    SectionHeaderReconPlan header = (SectionHeaderReconPlan){
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_RECON_PLAN,
        .h.compressed_offset     = BGEN32 (sizeof(SectionHeaderReconPlan)),
        .h.data_uncompressed_len = BGEN32 (txt_file->recon_plan.len * sizeof (ReconPlanItem)),
        .h.codec                 = codec,
        .h.flags.recon_plan.luft = false,
        .conc_writing_vbs        = 3, // VB being reconstructed, one of each PRIM and DEPN generated components at a time
        .vblock_mb               = BGEN32 ((uint32_t)(segconf.vb_size >> 20))
    };

    txt_file->recon_plan.len *= 3; // each ReconPlanItem is 3xuint32_t
    BGEN_u32_buf (&txt_file->recon_plan, 0);

    comp_compress (evb, &evb->z_data, (SectionHeader*)&header, txt_file->recon_plan.data, NULL);
}