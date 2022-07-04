// ------------------------------------------------------------------
//   vcf_ps_pid.c
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "vcf_private.h"
#include "seg.h"
#include "context.h"
#include "file.h"
#include "strings.h"
#include "dict_id.h"
#include "reconstruct.h"
#include "piz.h"
#include "lookback.h"
#include "stats.h"

// ------------------------
// FORMAT/PS and FORMAT/PID
// ------------------------

// PS/PID for a particular sample often skips variants, so we need to lookback more than a single line
#define MAX_PS_PID_LOOKBACK_LINES 16
#define PS_PID_LOOKBACK_LINES (flag.fast ? 4 : flag.best ? 16 : 10) // may be modified without affecting backward compatability

static char ps_lookback_snips[MAX_PS_PID_LOOKBACK_LINES][32], ps_pra_snip[200];
static unsigned ps_lookback_snip_lens[MAX_PS_PID_LOOKBACK_LINES], ps_pra_snip_len;

void vcf_samples_zip_initialize_PS_PID (void)
{
    // FORMAT/PS related stuff (0=lookback 1 line, 1=lookback 2 lines etc)
    for (int i=0; i < PS_PID_LOOKBACK_LINES; i++)
        seg_prepare_snip_other_chari (SNIP_LOOKBACK, _VCF_LOOKBACK, 'T'+i, ps_lookback_snip, i);

    SmallContainer con_PS_pos_ref_alt = {
        .repeats   = 1,
        .nitems_lo = 3,
        .items     = { { .dict_id = (DictId)_FORMAT_PSpos, .separator = "_"},
                       { .dict_id = (DictId)_FORMAT_PSref, .separator = "_"},
                       { .dict_id = (DictId)_FORMAT_PSalt                  } } };                       

    ps_pra_snip_len = sizeof (ps_pra_snip);
    container_prepare_snip ((ContainerP)&con_PS_pos_ref_alt, 0, 0, ps_pra_snip, &ps_pra_snip_len); 
}

// initialize VCF_LOOKBACK used by PS and PID
void vcf_samples_seg_initialize_LOOKBACK (VBlockVCFP vb)
{
    ContextP lookback_ctx = CTX(VCF_LOOKBACK);

    // note: we don't actually Seg anything into VCF_LOOKBACK as we get the lookback from the number of samples
    // this just carries the lb_size as an empty local section
    lookback_ctx->flags.store    = STORE_INT;
    lookback_ctx->ltype          = LT_DYN_INT;
    lookback_ctx->local_param    = true;
    lookback_ctx->local.prm8[0]  = lookback_size_to_local_param (PS_PID_LOOKBACK_LINES * vcf_num_samples + 1); // 1+ number of lookback values
    lookback_ctx->local_always   = (lookback_ctx->local.prm8[0] != 0); // no need for a SEC_LOCAL section if the parameter is 0 (which is the default anyway)
    lookback_ctx->is_initialized = true;

    CTX(VCF_SAMPLES)->flags.store = STORE_INDEX; // last_value is number of samples (=con.repeats)

    lookback_init (VB, lookback_ctx, CTX(FORMAT_PS),  STORE_INT); // lookback_ctx->local.param must be set before
    lookback_init (VB, lookback_ctx, CTX(FORMAT_PID), STORE_INT);
}

void vcf_samples_seg_initialize_PS_PID (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    // initialize all-the-same contexts for the REF and ALT container items of PS_POS_REF_ALT
    ctx_create_node (VB, FORMAT_PSref, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_COPY_REForALT, '0' }), 3);
    ctx_create_node (VB, FORMAT_PSalt, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_COPY_REForALT, '1' }), 3);
    ctx_create_node (VB, FORMAT_PSpos, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_COPYPOS,       '0' }), 3);

    // analyze PS or PID and determine its type
    if (ctx->did_i == FORMAT_PS && str_is_int (STRa(value)))   // this format appears only for PS, not PID
        ctx->ps_type = PS_POS;  // eg "73218731"
    
    else {
        str_split (value, value_len, 3, '_', item, true);
        if (n_items && str_is_int (STRi(item, 0))) 
            ctx->ps_type = PS_POS_REF_ALT; // eg "18182014_G_A"
        else
            ctx->ps_type = PS_UNKNOWN;
    }
}

void vcf_samples_seg_finalize_PS_PID (VBlockVCFP vb)
{
    // remove PSpos, PSref, PSalt if not needed
    if (CTX(FORMAT_PID)->ps_type || CTX(FORMAT_PS)->ps_type) {
        if (!ctx_get_count (VB, CTX(FORMAT_PSpos), 0)) ctx_free_context(CTX(FORMAT_PSpos), FORMAT_PSpos);
        if (!ctx_get_count (VB, CTX(FORMAT_PSref), 0)) ctx_free_context(CTX(FORMAT_PSref), FORMAT_PSref);
        if (!ctx_get_count (VB, CTX(FORMAT_PSalt), 0)) ctx_free_context(CTX(FORMAT_PSalt), FORMAT_PSalt);
    }

    // consolidate to the context actually used 
    if (CTX(FORMAT_PID)->nodes.len || CTX(FORMAT_PID)->ol_nodes.len)
        ctx_consolidate_stats (VB, FORMAT_PID, 5, VCF_LOOKBACK, FORMAT_PSref, FORMAT_PSalt, FORMAT_PSpos, DID_EOL);
    else
        ctx_consolidate_stats (VB, FORMAT_PS, 5, VCF_LOOKBACK, FORMAT_PSref, FORMAT_PSalt, FORMAT_PSpos, DID_EOL);
}

static inline bool vcf_seg_FORMAT_PS_PID_is_same_alt1 (VBlockVCFP vb, STRp(alt))
{
    if (alt_len == vb->main_alt_len)
        return !memcmp (alt, vb->main_refalt + vb->main_ref_len + 1, alt_len);

    // case: possibly multi-alt - compare to first ALT (ALT1)
    else {
        str_split (vb->main_refalt + vb->main_ref_len + 1, vb->main_alt_len, MAX_ALLELES, ',', alt, false);
        return str_issame_(STRa(alt), STRi(alt,0));
    }
}

// returns number ([1,PS_PID_LOOKBACK_LINES]) of lines back, or 0 if none
static inline unsigned vcf_seg_FORMAT_PS_PID_test_lookback (VBlockVCFP vb, ContextP ctx, STRp(value), uint32_t lookback)
{
    for (int lb_lines=1; lb_lines <= PS_PID_LOOKBACK_LINES; lb_lines++)
        if (lookback_is_same_txt (VB, VCF_LOOKBACK, ctx, lb_lines * lookback, STRa(value)))
            return lb_lines;

    return 0;
}

static inline bool vcf_seg_FORMAT_PS_PID_ps_matches_pid (VBlockVCFP vb, STRp(ps))
{
    rom pid = last_txt (VB, FORMAT_PID);
    unsigned pid_len = vb->last_txt_len (FORMAT_PID);

    return (pid_len > ps_len) && (pid[ps_len] == '_') && str_is_numeric (STRa(ps)) && !memcmp (ps, pid, ps_len);
}

// <ID=PS,Number=1,Type=Integer,Description="Phasing set (typically the position of the first variant in the set)">
// <ID=PS,Number=1,Type=Integer,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
// <ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
// Encountered formats: 1. PID="18182014_G_A", PS="18182014" ; 2. PS="18182014_G_A", no PID
// Value is the same as POS,REF,ALT of this line, or the same as PS/PID of this sample in one of the few previous lines
void vcf_seg_FORMAT_PS_PID (VBlockVCFP vb, ZipDataLineVCF *dl, ContextP ctx, STRp(value))
{
    int64_t ps_value;
    uint32_t lookback, lb_lines;    
    bool is_missing = (value_len==1 && *value=='.');

    // set ps_type if not already set
    if (!ctx->ps_type && !is_missing) 
        vcf_samples_seg_initialize_PS_PID (vb, ctx, STRa(value)); // set global segconf parameter - no harm even if multiple threads will set concurrently as regardless of the winner, the value is legitimate

    // case: '.' value - seg normally
    if (is_missing) goto fallback;

    // case: this is PS and we also have PID on this line - they are usually the same POS (so no need to lookback)
    else if (ctx->did_i == FORMAT_PS && CTX(FORMAT_PID)->ps_type == PS_POS_REF_ALT && 
             ctx_encountered (VB, FORMAT_PID) && vcf_seg_FORMAT_PS_PID_ps_matches_pid (vb, STRa(value)))
        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_PS_BY_PID }), 2, ctx, value_len);

    // case: this line is in the same Phase Set as the previous line
    else if ((lookback = (uint32_t)CTX(VCF_SAMPLES)->last_value.i/*last_line_num_samples*/) 
          && (lb_lines = vcf_seg_FORMAT_PS_PID_test_lookback (vb, ctx, STRa(value), lookback))) 
        seg_by_ctx (VB, STRi(ps_lookback_snip, lb_lines-1), ctx, value_len); 

    // case: not the same as previous line - seg according to ps_type
    else {
        // PS_POS (only applicable to PS, not PID): delta vs POS
        if (ctx->ps_type == PS_POS && vb->line_coords == DC_PRIMARY &&
            str_get_int (STRa(value), &ps_value)) {

            SNIPi2 (SNIP_SPECIAL, VCF_SPECIAL_COPYPOS, ps_value - dl->pos[0]);
            seg_by_ctx (VB, STRa(snip), ctx, value_len);
        }

        // PS_POS_REF_ALT: copy POS, REF, ALT1 in a container
        else if (ctx->ps_type == PS_POS_REF_ALT && vb->line_coords == DC_PRIMARY) {
            str_split (value, value_len, 3, '_', item, true);

            if (n_items && 
                str_issame_(STRi(item,0), last_txt(VB, VCF_POS), vb->last_txt_len(VCF_POS)) &&
                str_issame_(STRi(item,1), vb->main_refalt, vb->main_ref_len) &&
                vcf_seg_FORMAT_PS_PID_is_same_alt1(vb, STRi(item,2))) {

                // replace PS with COPY_POS, COPY_REF, COPY_ALT container 
                seg_by_ctx (VB, STRa(ps_pra_snip), ctx, value_len); 
                
                // items PSpos, PSref and PSalt of the container are all_the_same - initialized in vcf_samples_seg_initialize
                ctx_increment_count (VB, CTX(FORMAT_PSpos), 0);
                ctx_increment_count (VB, CTX(FORMAT_PSref), 0);
                ctx_increment_count (VB, CTX(FORMAT_PSalt), 0);
            }
            else   
                goto fallback;
        }
        else 
            fallback:            
            seg_by_ctx (VB, STRa(value), ctx, value_len); // segging once during segconf is enough - to create a context        
    }

    lookback_insert_txt (VB, VCF_LOOKBACK, ctx->did_i, STRa(value));

    seg_set_last_txt (VB, ctx, STRa(value));
}

// called if value is defined in FORMAT but missing at the end of the sample string
void vcf_seg_FORMAT_PS_PID_missing_value (VBlockVCFP vb, ContextP ctx, rom end_of_sample)
{
    lookback_insert_txt (VB, VCF_LOOKBACK, ctx->did_i, end_of_sample, 0); 

    // special case: a missing PS which follows a PID='.' - we generate a SPECIAL which then uses 
    // ctx->value_is_missing to achieve the same effect as WORD_INDEX_MISSING, so that b250 is have near-all SNIP_SPECIAL.
    // note: we DONT generate a SPECIAL if PS is '.'
    if (ctx->did_i == FORMAT_PS && ctx_encountered (VB, FORMAT_PID) && 
        (CTX(FORMAT_PID)->ps_type == PS_POS_REF_ALT || CTX(FORMAT_PID)->ps_type == PS_UNKNOWN) &&    
        vb->last_txt_len(FORMAT_PID)==1 && *last_txt(VB, FORMAT_PID)=='.') {

        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_PS_BY_PID }), 2, ctx, 0); 
    }
    else
        seg_by_ctx (VB, NULL, 0, ctx, 0); // generates WORD_INDEX_MISSING
}

//---------
// PIZ side
//---------

// called by compute thread, after uncompress, before reconstruct
void vcf_piz_initialize_ps_pid (VBlockP vb)
{
    if (CTX(FORMAT_PS)->dict.len)  // this file has PS
        lookback_init (vb, CTX(VCF_LOOKBACK), CTX (FORMAT_PS), STORE_INT);

    if (CTX(FORMAT_PID)->dict.len) // this file has PID
        lookback_init (vb, CTX(VCF_LOOKBACK), CTX (FORMAT_PID), STORE_INT);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_PS_by_PID)
{
    if (reconstruct) {
        STR(pid);
        reconstruct_peek (vb, CTX(FORMAT_PID), pSTRa(pid)); // note: we can't use last_txt, because PS might be reconstructed before PID, as its peeked by GT
        
        // if this SPECIAL was used with PID='.'
        if (pid_len == 1 && *pid == '.') 
            ctx->value_is_missing = true;

        else {
            unsigned ps_len = (char*)memchr (pid, '_', pid_len) - pid;
            ASSPIZ (ps_len < 15, "Failed to reconstructed PS from invalid PID: \"%.*s\"", pid_len, pid);

            RECONSTRUCT (pid, ps_len);
        }
    }

    return NO_NEW_VALUE;
}
