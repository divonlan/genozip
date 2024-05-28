// ------------------------------------------------------------------
//   vcf_format_GT.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"

//------------------
// ZIP
//------------------

// increase ploidy of the previous lines, if higher ploidy was encountered
static void vcf_seg_FORMAT_GT_increase_ploidy (VBlockVCFP vb, ContextP ht_ctx, unsigned new_ploidy)
{
    buf_alloc (vb, &ht_ctx->local, 0, new_ploidy * vcf_num_samples * vb->lines.len, 
               char, CTX_GROWTH, CTX_TAG_LOCAL);

    uint32_t num_samples = ht_ctx->HT_n_lines * vcf_num_samples + vb->sample_i; // all samples in previous lines + previous samples in current line
    char *ht_data = B1STc (ht_ctx->local);

    // copy the haplotypes backwards (to avoid overlap), padding with '*' (which are NOT counted in .repeats of the GT container)
    for (int sam_i = num_samples-1; sam_i >= 0; sam_i--) {

        int ht_i=new_ploidy-1 ; for (; ht_i >= vb->ploidy; ht_i--) 
            ht_data[sam_i * new_ploidy + ht_i] = '*'; 

        for (; ht_i >= 0; ht_i--)
            ht_data[sam_i * new_ploidy + ht_i] = ht_data[sam_i * vb->ploidy + ht_i];
    }

    vb->ploidy = new_ploidy;
    ht_ctx->ht_per_line = vb->ploidy * vcf_num_samples;
}

WordIndex vcf_seg_FORMAT_GT (VBlockVCFP vb, ContextP ctx, ZipDataLineVCF *dl, STRp(gt), bool has_ps, bool has_null_dp)
{
    ContextP ht_ctx = CTX(FORMAT_GT_HT);

    // the GT field is represented as a Container, with a single item repeating as required by poidy, and the separator 
    // determined by the phase
    MiniContainer con = { .repeats   = 1, 
                          .nitems_lo = 1, 
                          .callback  = true, // see vcf_piz_container_cb
                          .drop_final_repsep     = true, 
                          .items[0].dict_id.num  = _FORMAT_GT_HT };

    unsigned save_gt_len = gt_len;

    // update repeats according to ploidy, and separator according to phase
    for (unsigned i=1; i < gt_len-1; i++)
        if (gt[i] == '|' || gt[i] == '/') {
            con.repeats++;
            con.repsep[0] = gt[i];
        }

    ASSVCF (con.repeats <= VCF_MAX_PLOIDY, "ploidy=%u exceeds the maximum of %u", con.repeats, VCF_MAX_PLOIDY);
    
    // if the ploidy of this line is bigger than the ploidy of the data in this VB so far, then
    // we have to increase ploidy of all the haplotypes read in in this VB so far. This can happen for example in 
    // the X chromosome if initial samples are male with ploidy=1 and then a female sample with ploidy=2
    if (vb->ploidy && con.repeats > vb->ploidy) 
        vcf_seg_FORMAT_GT_increase_ploidy (vb, ht_ctx, con.repeats);

    if (!vb->ploidy) {
        vb->ploidy = con.repeats; // very first sample in the vb
        ht_ctx->ht_per_line = vb->ploidy * vcf_num_samples;
    }

    buf_alloc (vb, &ht_ctx->local, vb->ploidy, ht_ctx->ht_per_line * vb->lines.len, char, CTX_GROWTH, CTX_TAG_LOCAL);

    // note - ploidy of this sample might be smaller than vb->ploidy (eg a male sample in an X chromosesome that was preceded by a female sample, or "." sample)
    Allele *ht_data = B(Allele, ht_ctx->local, ht_ctx->HT_n_lines * ht_ctx->ht_per_line + vb->ploidy * vb->sample_i);

    // number of allele other than REF (i.e. ht >= 1)
    // note: up to 15.0.35, dosage was -1 if any allele was not '0' or '1', or if ploidy exceeded 2
    int64_t dosage=0; 
    
    for (unsigned ht_i=0; ht_i < con.repeats; ht_i++) {

        Allele ht = *(gt++); 
        gt_len--;

        ASSVCF (IS_DIGIT(ht) || ht == '.', 
                "invalid VCF file - expecting an allele in a sample to be a number 0-99 or . , but seeing %c (ht_i=%u)", ht, ht_i);

        // single-digit allele numbers
        ht_data[ht_i] = ht;

        // calculate dosage contribution of this ht (to be used in vcf_seg_get_mux_channel_i)
        if (ht != '0' && ht != '.') dosage++; 
    
        if (ht != '.') CTX(INFO_AN)->an.count_ht++;

        if (!gt_len) break;

        // handle 2-digit allele numbers
        if (ht != '.' && IS_DIGIT (*gt)) {
            unsigned allele = 10 * (ht-'0') + (*(gt++) - '0');
            gt_len--;

            // make sure there isn't a 3rd digit
            ASSVCF (!gt_len || !IS_DIGIT (*gt), "VCF file sample %u - genozip currently supports only alleles up to 99", vb->sample_i+1);

            ht_data[ht_i] = '0' + allele; // use ascii 48->147
        }
        
        // read and verify phase
        if (con.repeats > 1 && ht_i < con.repeats-1) {
            
            char phase = *(gt++);
            gt_len--;

            ASSVCF (phase != ' ', "invalid VCF file - expecting a tab or newline after sample %u but seeing a space", vb->sample_i+1);
            ASSVCF (phase == con.repsep[0], "invalid VCF file -  unable to parse sample %u: expecting a %c but seeing %c", vb->sample_i+1, con.repsep[0], phase);
        }
    } // for characters in a sample

    ctx->gt.actual_last_ploidy = con.repeats; // set before increasing con.repeats to vb->ploidy

    // if the ploidy of the sample is lower than vb->ploidy, set missing ht as '-' (which will cause deletion of themselves and their separator)
    // and set the ploidy to vb->ploidy - to avoid increase in entroy of GT.b250
    if (con.repeats != vb->ploidy) {
        
        for (unsigned ht_i=con.repeats; ht_i < vb->ploidy; ht_i++) 
            ht_data[ht_i] = '-'; // unlike '*', we DO count '-' in .repeats (so that we can have the same number of repeats = lower entroy in GT.b250)

        con.repeats = vb->ploidy;
        if (!con.repsep[0]) con.repsep[0] = ctx->gt.prev_phase; // this happens in case if a 1-ploid sample
    }

    // case DP='.' - we predict that GT=./.
    bool no_duplicate=false;
    if (has_null_dp) { // segconf.use_null_DP_method=true && line has DP && it is '.'
        // case: prediction is correct - re-write GT as "0/0" or "0|0" (use previous phase)
        if (con.repeats==2 && ht_data[0]=='.' && con.repsep[0]=='/' && ht_data[1]=='.') {
            ht_data[0] = ht_data[1] = '0';
            if (ctx->gt.prev_phase) con.repsep[0] = ctx->gt.prev_phase;
        }

        // case: prediction is incorrect tell piz to NOT change the HTs and phase according to the prediction.
        // Modifting con here adds entropy to the GT context, hopefully the prediction is usually correct and this rarely happens 
        else {
            con.items[0].separator[1] = CI1_ITEM_PRIVATE; // override (hopefully almost never used)
            no_duplicate = true; // prevent seg_duplicate_last - we shouldn't copy the previous sample and the next sample shouldn't copy us
        }
    }

    // case: phase is predictable from has_ps and ht_data[0]
    // note: a case where this prediction fails is with alleles > 1 eg "1|2". In GIAB, this never have a PS, but can be | or / 
    else if (CTX(FORMAT_PS)->ps_type && con.repeats==2 &&
            ((ht_data[0] != '.' && (has_ps == (con.repsep[0]=='|'))) || // ht!=. --> predicted to be as has_ps says
             (ht_data[0] == '.' && con.repsep[0]=='/'))) {              // ht==. --> predicted to be /

        // generate the phase from 
        con.repsep[0] = '&'; // re-write from prediction
        con.callback = true; // vcf_piz_filter will re-write repsep based on prediction
    }

    // if this sample is a "./." - replace it with "%|%" or "%/%" according to the previous sample's phase -  
    // so that the container is likely identical and we reduce GT.b250 entropy. Reason: many tools
    // (including bcftools merge) produce "./." for missing samples even if all other samples are phased
    else if (ht_data[0]=='.' && con.repeats==2 && ht_data[1]=='.' && con.repsep[0]=='/') {
        con.repsep[0] = ctx->gt.prev_phase ? ctx->gt.prev_phase : '|'; // '|' is arbitrary
        ht_data[0] = ht_data[1] = '%';
    }

    // in case we have INFO/SF, we verify that it is indeed the list of samples for which the first ht is not '.'
    if (CTX(INFO_SF)->sf.SF_by_GT == yes && ht_data[0] != '.') 
        vcf_seg_INFO_SF_one_sample (vb);

    ctx_set_last_value (VB, ctx, dosage); // to be used in vcf_seg_get_mux_channel_i

    ASSVCF (!gt_len, "Invalid GT data in sample_i=%u", vb->sample_i+1);
    
    // shortcut if we have the same ploidy and phase as previous GT (saves re-genetrating base64 in container_seg)
    if (con.repeats == ctx->gt.prev_ploidy && con.repsep[0] == ctx->gt.prev_phase && !no_duplicate) 
        return seg_duplicate_last (VB, ctx, save_gt_len);

    else {
        ctx->gt.prev_ploidy = no_duplicate ? 0 : con.repeats; // if no_duplicate - 0 to prevent next sample from duplicating this one
        ctx->gt.prev_phase  = con.repsep[0];
        return container_seg (vb, ctx, (ContainerP)&con, 0, 0, save_gt_len); 
    }
}

void vcf_seg_FORMAT_GT_finalize_line (VBlockVCFP vb, uint32_t line_n_samples)
{
    ContextP ht_ctx = CTX(FORMAT_GT_HT);
    if (!ht_ctx->use_HT_matrix) return;

    // some real-world files encountered have too-short lines due to human errors. we pad them
    if (line_n_samples < vcf_num_samples) {
        char *ht_data = Bc(ht_ctx->local, ht_ctx->HT_n_lines * vb->ploidy * vcf_num_samples + vb->ploidy * line_n_samples);
        uint32_t num_missing = vb->ploidy * (vcf_num_samples - line_n_samples); 
        memset (ht_data, '*', num_missing);
    }

    ht_ctx->HT_n_lines++;
    ht_ctx->local.len32 = ht_ctx->HT_n_lines * ht_ctx->ht_per_line; // not just "+=ht_per_line", so this also works if ploidy has increased
}

//------------------
// PIZ
//------------------

// return dosage derived from GT - equivalent to dosage calculation in vcf_seg_FORMAT_GT
int vcf_piz_GT_get_last_dosage (VBlockP vb)
{
    STRlast(gt, FORMAT_GT);
    int dosage = 0;

    if (z_file->max_ploidy_for_mux)  // since 15.0.36
        for (unsigned i=0; i < gt_len; i++) {
            if (gt[i] == '/' || gt[i] == '|') continue;
            if (gt[i] >= '1' && gt[i] <= '9') dosage++;
            if (i < gt_len-1 && IS_DIGIT(gt[i+1])) i++;
        }
    
    else { // up to 15.0.35
        for (unsigned i=1; i < gt_len; i += 2)
            if (gt[i] != '/' && gt[i] != '|') return -1; // we have an allele >= 10 - dosage is -1

        for (unsigned i=0; i < gt_len; i += 2) {
            if (gt[i] != '0' && gt[i] != '1') return -1; // dosage only defined if all allele are 0 or 1
            dosage += (gt[i] - '0');
        }
    }

    return dosage;
}

static inline bool vcf_piz_is_in_FORMAT (VBlockP vb, rom tag/*2-chars*/)
{
    STRlast (format, VCF_FORMAT);

    // check if a :PS: or :PS\t exists (using PS as an example)
    for (int i=3; i < format_len-1; i++) // start from 3 - skip "GT:"
        if (format[i]==tag[0] && format[i+1]==tag[1] && format[i-1]==':' && (format[i+2]==':' || format[i+2]=='\t'))
            return true;

    return false; // not found
}

// rewrite '&' to the predicted phase 
void vcf_piz_FORMAT_GT_rewrite_predicted_phase (VBlockP vb, char *recon, uint32_t recon_len)
{
    // prediction is '/' by default. it is '|' only if: 
    // 1. ht[0] is not '.' AND 2. we have PS on this line AND 3. PS is not (empty or '.')
    char prediction = '/';

    if (*recon != '.' && vcf_piz_is_in_FORMAT (vb, "PS")) {

        // case: if its a non-trivial PS - not empty or '.' - the phase is |, otherwise it's /
        STR(ps);
        reconstruct_peek (vb, CTX(FORMAT_PS), pSTRa(ps));
        if (ps_len && !(ps_len==1 && *ps=='.')) prediction = '|';
    }
        
    // we only use prediction in in diploid GTs, and genozip allows only single or double-digit alleles
    if      (recon[1] == '&') recon[1] = prediction;
    else if (recon[2] == '&') recon[2] = prediction;
    else ABORT_PIZ ("Cannot find '&' predictable phase in recon: \"%.*s\"", recon_len, recon);
}

// called after 2nd HT in a diploid GT, if has_null_DP alg is activated 
void vcf_piz_GT_cb_null_GT_if_null_DP (VBlockP vb , char *recon)
{
    if (recon[0]=='0' && recon[-2]=='0' && recon[-3]=='\t' && vcf_piz_is_in_FORMAT (vb, "DP")) { 

        STR(DP);
        reconstruct_peek (VB, CTX(FORMAT_DP), pSTRa(DP));

        if (IS_PERIOD (DP)) {
            recon[-2] = recon[0] = '.';
            recon[-1] = '/';
        }
    }
}

