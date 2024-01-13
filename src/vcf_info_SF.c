// ------------------------------------------------------------------
//   vcf_info_SF.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "piz.h"
#include "context.h"
#include "strings.h"
#include "dict_id.h"
#include "reconstruct.h"

#define adjustment ctx->last_delta

static void vcf_seg_INFO_SF_seg (VBlockP vb); // forward

// INFO/SF contains a comma-seperated list of the 0-based index of the samples that are NOT '.'
// example: "SF=0,1,15,22,40,51,59,78,88,89,90,112,124,140,147,155,156,164,168,183,189,197,211,215,216,217,222,239,244,256,269,270,277,281,290,291,299,323,338,340,348" 
// Algorithm: SF is segged either as an as-is string, or as a SPECIAL that includes the index of all the non-'.' samples. 
// if sf.SF_by_GT=YES, we use SNIP_SPECIAL and we validate the correctness during vcf_seg_FORMAT_GT -
// if it is wrong we set sf.SF_by_GT=NO. The assumption is that normally, it is either true for all lines or false.
bool vcf_seg_INFO_SF_init (VBlockVCFP vb, ContextP ctx, STRp(sf))
{
    SEGCONF_RECORD_WIDTH (SF, sf_len);

    switch (ctx->sf.SF_by_GT) {
        case no: 
            return true; // "special" is suppressed - caller should go ahead and seg normally

        case unknown: 
            ctx->sf.SF_by_GT = yes; // first call to this function, after finding that we have an INFO/SF field - set and fall through

        case yes: 
            seg_set_last_txt (VB, ctx, STRa(sf));

            vb_add_to_deferred_q (VB, ctx, vcf_seg_INFO_SF_seg, vb->idx_SF);

            ctx->sf.next = 0; // relative to last_txt
            adjustment = 0;      
            
            // snip being contructed 
            buf_alloc (vb, &ctx->deferred_snip, 0, sf_len + 20, char, 2, "contexts->deferred_snip"); // initial value - we will increase if needed
            BNXTc (ctx->deferred_snip) = SNIP_SPECIAL;
            BNXTc (ctx->deferred_snip) = VCF_SPECIAL_SF;

            return false; // caller should not seg as we already did

        default:
            ABOSEG ("invalid sf.SF_by_GT=%d", ctx->sf.SF_by_GT);
    }
}

// verify next number on the list of the SF field is sample_i (called from vcf_seg_FORMAT_GT)
void vcf_seg_INFO_SF_one_sample (VBlockVCFP vb)
{
    decl_ctx(INFO_SF);
    STRlast (sf, INFO_SF);

    // case: no more SF values left to compare - we ignore this sample
    while (ctx->sf.next < sf_len) {

        buf_alloc (vb, &ctx->deferred_snip, 10, 0, char, 2, "contexts->deferred_snip");

        rom sf_one_value = &sf[ctx->sf.next]; 
        char *after;
        int32_t value = strtol (sf_one_value, &after, 10); // terminated with , ; \t (not \n bc we know there are samples)

        int32_t adjusted_sample_i = (int32_t)(vb->sample_i + adjustment); // adjustment is the number of values in SF that are not in samples

        // case: badly formatted SF field
        if (*after != ',' && *after != ';' && *after != '\t') {
            CTX(INFO_SF)->sf.SF_by_GT = no; // failed - turn off for the remainder of this vb
            break;
        }
        
        // case: value exists in SF and samples
        else if (value == adjusted_sample_i) {
            BNXTc (ctx->deferred_snip) = ',';
            ctx->sf.next = after - sf + 1; // +1 to skip comma
            break;
        }

        // case: value in SF file doesn't appear in samples - keep the value in the snip
        else if (value < adjusted_sample_i) {
            ctx->deferred_snip.len += str_int (value, BAFTc (ctx->deferred_snip));
            BNXTc (ctx->deferred_snip) = ',';
            adjustment++;
            ctx->sf.next = after - sf + 1; // +1 to skip comma
            // continue and read the next value
        }

        // case: value in SF is larger than current sample - don't advance iterator - perhaps future sample will cover it
        else { // value > adjusted_sample_i
            BNXTc (ctx->deferred_snip) = '~'; // skipped sample
            break; 
        }
    }
}

static void vcf_seg_INFO_SF_seg (VBlockP vb_)
{   
    VBlockVCFP vb = (VBlockVCFP)vb_;
    decl_ctx(INFO_SF);
    STRlast (sf, INFO_SF);

    // case: SF data remains after all samples - copy it
    int32_t remaining_len = (int32_t)ctx->last_txt.len - (int32_t)ctx->sf.next; // -1 if all done, because we skipped a non-existing comma
    if (remaining_len > 0) {
        buf_add_more (VB, &ctx->deferred_snip, &sf[ctx->sf.next], remaining_len, "contexts->deferred_snip");
        BNXTc (ctx->deferred_snip) = ','; // buf_add_more allocates one character extra
    }

    if (CTX(INFO_SF)->sf.SF_by_GT == yes) 
        seg_by_ctx (VB, STRb(ctx->deferred_snip), ctx, sf_len);
    
    else if (CTX(INFO_SF)->sf.SF_by_GT == no)
        seg_by_ctx (VB, STRa(sf), ctx, sf_len);

    buf_free (ctx->deferred_snip);
}

#define sample_i ctx->last_value.i
#define snip_i   ctx->deferred_snip.param

// leave space for reconstructing SF - actual reconstruction will be in vcf_piz_container_cb
SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_INFO_SF)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    if (reconstruct) {
        if (!flag.drop_genotypes && !flag.samples) {
            adjustment = sample_i = snip_i = 0;

            // temporary place for SF
            buf_alloc (vb, &ctx->insertion, 0, 5 * vcf_header_get_num_samples(), char, 1, "contexts->insertion"); // initial estimate, we may further grow it later
            ctx->insertion.len = 0;

            // save snip for later (note: the SNIP_SPECIAL+code are already removed)
            buf_add_moreS (vb, &ctx->deferred_snip, snip, "contexts->deferred_snip");

            vcf_piz_defer_to_after_samples (SF);
        }
        else
            RECONSTRUCT ("N/A", 3);
    }

    return NO_NEW_VALUE;
}

// While reconstructing the GT fields of the samples - calculate the INFO/SF field
void vcf_piz_GT_cb_calc_INFO_SF (VBlockVCFP vb, unsigned rep, char *recon, int32_t recon_len)
{
    decl_ctx (INFO_SF);

    if (rep != 0 || !CTX(INFO_SF)->recon_insertion) return; // we only look at the first ht in a sample, and only if its not '.'/'%'

    if (*recon == '.' || *recon == '%') { // . can be written as % in vcf_seg_FORMAT_GT
        sample_i++;
        return;
    }

    ARRAY (const char, snip, ctx->deferred_snip);

    while (snip_i < snip_len) {

        buf_alloc (vb, &ctx->insertion, 12, 0, char, 2, "insertion"); // sufficient for int32 + ','

        int32_t adjusted_sample_i = (int32_t)(sample_i + adjustment); // last_delta is the number of values in SF that are not in samples

        // case: we derive the value from sample_i
        if (snip[snip_i] == ',') {
            snip_i++;

            buf_add_int_as_text (&ctx->insertion, adjusted_sample_i);

            if (snip_i < snip_len) // add comma if not done yet
                BNXTc (ctx->insertion) = ',';

            break;
        }

        // case: sample was not used to derive any value
        if (snip[snip_i] == '~') {
            snip_i++;
            break;
        }

        // case: value quoted verbatim - output and continue searching for , or ~ for this sample
        else {
            char *comma = strchr (&snip[snip_i], ',');
            unsigned len = comma - &snip[snip_i];
            bool add_comma = (snip_i + len + 1 < snip_len); // add comma if not last

            buf_add (&ctx->insertion, &snip[snip_i], len + add_comma);
            snip_i += len + 1;

            adjustment++;
        }
    }

    sample_i++;
}

// Upon completing the line - insert the calculated INFO/SF field to its place
void vcf_piz_insert_INFO_SF (VBlockVCFP vb)
{
    decl_ctx(INFO_SF);

    if (!ctx->recon_insertion) return;

    ARRAY (const char, snip, ctx->deferred_snip);

    // if there are some items remaining in the snip (values that don't appear in samples) - copy them
    if (snip_i < snip_len) {
        unsigned remaining_len = snip_len - snip_i - 1; // all except the final comma
        buf_add_more ((VBlockP)vb, &ctx->insertion, &snip[snip_i], remaining_len, "contexts->insertion"); 
    }

    // make room for the SF txt and copy it to its final location
    vcf_piz_insert_field (vb, ctx, STRb(ctx->insertion), segconf.wid_SF.width);

    buf_free (ctx->deferred_snip);
    buf_free (ctx->insertion);
}
