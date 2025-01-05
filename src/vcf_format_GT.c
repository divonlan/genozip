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

// splits a string with up to (max_items-1) separators (doesn't need to be nul-terminated) to up to or exactly max_items integers
// returns the actual number of items, or 0 is unsuccessful
static uint32_t str_split_gt_do (VBlockVCFP vb, STRp(gt), 
                                 int32_t *items/*-1 means '.'*/, char *phase) // out
{
    rom c = gt;
    rom after = &gt[gt_len];
    SAFE_NUL (after);

    uint32_t item_i;
    for (item_i=0; item_i < VCF_MAX_PLOIDY && c < after; item_i++) {
        if (*c == '.') {
            items[item_i] = -1;
            c++;
        }
        
        else {
            items[item_i] = 0;
            while (IS_DIGIT(*c))
                items[item_i] = items[item_i] * 10 + (*c++ - '0');
        }

        if (item_i < VCF_MAX_PLOIDY-1 && (*c == '/' || *c == '|')) {
            if (!item_i) 
                *phase = *c;
            else 
                ASSVCF (*c == *phase, "Genozip limitation: cannot handle multiploid GT=\"%.*s\" field with both / and  |", STRf(gt));
            c++;
        }

        else if (item_i < VCF_MAX_PLOIDY-1 && *c != '/' && *c != '|' && *c != 0) 
            break; // fail - not integer
    }

    ASSSEG (c == after, "Invalid GT value \"%*.s\"", STRf(gt));

    SAFE_RESTORE;
    return item_i;
}
#define str_split_gt(gt,gt_len,name,phase)          \
    int32_t name[VCF_MAX_PLOIDY];                   \
    char phase=0;                                   \
    int n_##name##s;                                \
    if ((gt_len) == 3 && ((gt)[1] == '/' || (gt)[1] == '|')) { /* shortcut for most common cases */\
        n_##name##s = 2;                            \
        name[0] = (gt)[0]=='.' ? -1 : ((gt)[0]-'0');\
        name[1] = (gt)[2]=='.' ? -1 : ((gt)[2]-'0');\
        phase   = (gt)[1];                          \
    }                                               \
    else if (gt_len == 1) {                         \
        n_##name##s = 1;                            \
        name[0] = (gt)[0]=='.' ? -1 : ((gt)[0]-'0');\
    }                                               \
    else                                            \
        n_##name##s = str_split_gt_do (vb, (gt), (gt_len), name, &phase); 

// increase ploidy of the previous lines, if higher ploidy was encountered
static void vcf_seg_FORMAT_GT_increase_ploidy (VBlockVCFP vb, ContextP ht_ctx, unsigned new_ploidy)
{
    if (GT_USES_PBWT) {
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
    }

    vb->ploidy = new_ploidy;
    ht_ctx->ht_per_line = vb->ploidy * vcf_num_samples;
}

static void vcf_seg_analyze_GT (VBlockVCFP vb, int n_hts, int32_t *ht)
{
    // number of allele other than REF (i.e. ht >= 1)
    // note: up to 15.0.35, dosage was -1 if any allele was not '0' or '1', or if ploidy exceeded 2
    int64_t dosage=0; 
    for (int ht_i=0; ht_i < n_hts; ht_i++)
        if (ht[ht_i] >= 0) {
            CTX(INFO_AN)->an.count_ht++;
            if (ht[ht_i] >= 1) dosage++;
        }
    
    ctx_set_last_value (VB, CTX(FORMAT_GT), dosage); // to be used in vcf_seg_get_mux_channel_i

    // in case we have INFO/SF, we verify that it is indeed the list of samples for which the first ht is not '.'
    if (CTX(INFO_SF)->sf.SF_by_GT == yes && ht[0] >= 0) 
        vcf_seg_INFO_SF_one_sample (vb);
}

void vcf_seg_analyze_copied_GT (VBlockVCFP vb, STRp(gt))
{
    START_TIMER;

    str_split_gt (gt, gt_len, ht, phase);
    vcf_seg_analyze_GT (vb, n_hts, ht);

    COPY_TIMER (vcf_seg_analyze_copied_GT);
}

void vcf_seg_FORMAT_GT (VBlockVCFP vb, ContextP ctx, ZipDataLineVCF *dl, STRp(gt), ContextP *ctxs, STRps(sf))
{
    set_last_txt (FORMAT_GT, gt);

    ContextP ht_ctx = CTX(FORMAT_GT_HT);

    str_split_gt (gt, gt_len, h, phase);

    // the GT field is represented as a Container, with a single item repeating as required by poidy, and the separator 
    // determined by the phase
    MiniContainer con = { .repeats              = n_hs, 
                          .nitems_lo            = 1, 
                          .callback             = true, // see vcf_piz_container_cb
                          .drop_final_repsep    = true, 
                          .items[0].dict_id.num = _FORMAT_GT_HT,
                          .repsep[0]            = phase };

    unsigned save_gt_len = gt_len;
    
    // if the ploidy of this line is bigger than the ploidy of the data in this VB so far, then
    // we have to increase ploidy of all the haplotypes read in in this VB so far. This can happen for example in 
    // the X chromosome if initial samples are male with ploidy=1 and then a female sample with ploidy=2
    if (vb->ploidy && n_hs > vb->ploidy) 
        vcf_seg_FORMAT_GT_increase_ploidy (vb, ht_ctx, n_hs);

    if (!vb->ploidy) {
        vb->ploidy = n_hs; // very first sample in the vb
        ht_ctx->ht_per_line = n_hs * vcf_num_samples;
    }

    vcf_seg_analyze_GT (vb, n_hs, h);

    // store first 2 hts
    ctx->gt.ht[0] =                h[0]==-1?'.' : ('0' + (uint8_t)(MIN_(h[0], NUM_SMALL_ALLELES-1)));
    ctx->gt.ht[1] = (n_hs < 2)?0 : h[1]==-1?'.' : ('0' + (uint8_t)(MIN_(h[1], NUM_SMALL_ALLELES-1)));
    
    // case: copy_sample method: we seg GT as a snip in samples that are not copied ()
    if (segconf.vcf_sample_copy) 
        seg_mux_by_is_prev_sample_copied (vb, dl, ctx, &vb->mux_GT, STRa(gt)); 

    else if (!GT_USES_PBWT)
        seg_by_ctx (VB, STRa(gt), ctx, gt_len);

    else {
        buf_alloc (vb, &ht_ctx->local, vb->ploidy, ht_ctx->ht_per_line * vb->lines.len, char, CTX_GROWTH, CTX_TAG_LOCAL);

        // note - ploidy of this sample might be smaller than vb->ploidy (eg a male sample in an X chromosesome that was preceded by a female sample, or "." sample)
        Allele *ht_data = B(Allele, ht_ctx->local, ht_ctx->HT_n_lines * ht_ctx->ht_per_line + vb->ploidy * vb->sample_i);

        for (int ht_i=0; ht_i < n_hs; ht_i++) 
            if (h[ht_i] == -1)
                ht_data[ht_i] = '.';

            else if (h[ht_i] < NUM_SMALL_ALLELES)            
                ht_data[ht_i] = '0' + (uint8_t)h[ht_i]; // use ascii 48->255,0->36 (rely on uint8_t arithmetic to round robin)
            
            else { // big allele: '&' in ht_data and actual allele value in FORMAT_GT_HT_BIG
                ht_data[ht_i] = '&';
                seg_integer (VB, CTX(FORMAT_GT_HT_BIG), h[ht_i] - NUM_SMALL_ALLELES, false, 0);
            }

        ctx->gt.actual_last_ploidy = n_hs; // set before increasing con.repeats to vb->ploidy

        // if the ploidy of the sample is lower than vb->ploidy, set missing ht as '-' (which will cause deletion of themselves and their separator)
        // and set the ploidy to vb->ploidy - to avoid increase in entroy of GT.b250
        if (n_hs != vb->ploidy) {
            for (int ht_i=n_hs; ht_i < vb->ploidy; ht_i++) 
                ht_data[ht_i] = '-'; // unlike '*', we DO count '-' in .repeats (so that we can have the same number of repeats = lower entroy in GT.b250)

            con.repeats = vb->ploidy;
            if (!con.repsep[0]) con.repsep[0] = ctx->gt.prev_phase; // this happens in case if a 1-ploid sample
        }

        // case DP='.' - we predict that GT=./.
        bool no_duplicate=false;

        if (segconf.use_null_DP_method && vcf_seg_sample_has_null_value (FORMAT_DP, ctxs, STRas(sf))) { // segconf.use_null_DP_method=true && line has DP && it is '.'
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
                ((ht_data[0] != '.' && (vcf_seg_sample_has_PS (vb, ctxs, STRas(sf)) == (con.repsep[0]=='|'))) || // ht!=. --> predicted to be as has_ps says
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
        
        // shortcut if we have the same ploidy and phase as previous GT (saves re-genetrating base64 in container_seg)
        if (con.repeats == ctx->gt.prev_ploidy && con.repsep[0] == ctx->gt.prev_phase && !no_duplicate) 
            seg_duplicate_last (VB, ctx, save_gt_len);

        else {
            ctx->gt.prev_ploidy = no_duplicate ? 0 : con.repeats; // if no_duplicate - 0 to prevent next sample from duplicating this one
            ctx->gt.prev_phase  = con.repsep[0];
            container_seg (vb, ctx, (ContainerP)&con, 0, 0, save_gt_len); 
        }
    }
}

void vcf_seg_FORMAT_GT_finalize_line (VBlockVCFP vb, uint32_t line_n_samples)
{
    if (!GT_USES_PBWT) return;

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
int vcf_piz_GT_get_last_dosage (VBlockVCFP vb)
{
    STRlast(gt, FORMAT_GT);
    int dosage = 0;

    if (z_file->max_ploidy_for_mux)  // since 15.0.36
        for (unsigned i=0; i < gt_len; i++) {
            if (gt[i] == '/' || gt[i] == '|') continue;
            if (gt[i] >= '1' && gt[i] <= '9') dosage++; // inspect first digit of allele
            while (i < gt_len-1 && IS_DIGIT(gt[i+1])) i++; // skip remaining digits
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

static inline bool vcf_piz_is_in_FORMAT (VBlockVCFP vb, rom tag/*2-chars*/)
{
    STRlast (format, VCF_FORMAT);

    // check if a :PS: or :PS\t exists (using PS as an example)
    for (int i=3; i < format_len-1; i++) // start from 3 - skip "GT:"
        if (format[i]==tag[0] && format[i+1]==tag[1] && format[i-1]==':' && (format[i+2]==':' || format[i+2]=='\t'))
            return true;

    return false; // not found
}

// rewrite '&' to the predicted phase 
void vcf_piz_FORMAT_GT_rewrite_predicted_phase (VBlockVCFP vb, char *recon, uint32_t recon_len)
{
    // prediction is '/' by default. it is '|' only if: 
    // 1. ht[0] is not '.' AND 2. we have PS on this line AND 3. PS is not (empty or '.')
    char prediction = '/';

    if (*recon != '.' && vcf_piz_is_in_FORMAT (vb, "PS")) {

        // case: if its a non-trivial PS - not empty or '.' - the phase is |, otherwise it's /
        STR(ps);
        reconstruct_peek (VB, CTX(FORMAT_PS), pSTRa(ps));
        if (ps_len && !(ps_len==1 && *ps=='.')) prediction = '|';
    }
        
    // we only use prediction in in diploid GTs, and genozip allows only single or double-digit alleles
    if      (recon[1] == '&') recon[1] = prediction;
    else if (recon[2] == '&') recon[2] = prediction;
    else ABORT_PIZ ("Cannot find '&' predictable phase in recon: \"%.*s\"", recon_len, recon);
}

// called after 2nd HT in a diploid GT, if has_null_DP alg is activated 
void vcf_piz_GT_cb_null_GT_if_null_DP (VBlockVCFP vb , char *recon)
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

void vcf_piz_GT_update_other_fields (VBlockVCFP vb, rom recon) // note: only first char of recon is used (i.e. ht=0)
{
    CTX(INFO_AN)->an.count_ht += (*recon != '.'); // used for reconstructing INFO/AN

    if (CTX(INFO_SF)->deferred_snip.len > 0/*have INFO_SF*/) 
        vcf_piz_GT_cb_calc_INFO_SF (vb, recon);
}
