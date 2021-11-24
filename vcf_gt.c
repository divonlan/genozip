// ------------------------------------------------------------------
//   vcf_gt.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "vcf_private.h"
#include "seg.h"
#include "context.h"
#include "file.h"
#include "strings.h"
#include "dict_id.h"
#include "reconstruct.h"
#include "piz.h"

//------------------
// ZIP
//------------------

// complete haplotypes of lines that don't have GT, if any line in the vblock does have GT.
// In this case, the haplotype matrix must include the lines without GT too
void vcf_seg_FORMAT_GT_complete_missing_lines (VBlockVCF *vb)
{
    buf_alloc (vb, &CTX(FORMAT_GT_HT)->local, 0, vb->lines.len * vb->ht_per_line, char, CTX_GROWTH, "contexts->local");

    for (vb->line_i=0; vb->line_i < (uint32_t)vb->lines.len; vb->line_i++) {

        if (CTX(FORMAT_GT_HT) && !DATA_LINE (vb->line_i)->has_haplotype_data) {
            char *ht_data = ENT (char, CTX(FORMAT_GT_HT)->local, vb->line_i * vb->ht_per_line);
            memset (ht_data, '*', vb->ht_per_line);

            // NOTE: we DONT set dl->has_haplotype_data to true bc downstream we still
            // count this row as having no GT field when analyzing gt data
        }
    }

    CTX(FORMAT_GT_HT)->local.len = vb->lines.len * vb->ht_per_line;
}

// increase ploidy of the previous lines, if higher ploidy was encountered
static void vcf_seg_FORMAT_GT_increase_ploidy (VBlockVCF *vb, unsigned new_ploidy)
{
    buf_alloc (vb, &CTX(FORMAT_GT_HT)->local, 0, new_ploidy * vcf_num_samples * vb->line_i, 
               char, CTX_GROWTH, "contexts->local");

    uint32_t num_samples = vb->line_i * vcf_num_samples + vb->sample_i; // all samples in previous lines + previous samples in current line
    char *ht_data = FIRSTENT (char, CTX(FORMAT_GT_HT)->local);

    // copy the haplotypes backwards (to avoid overlap), padding with '*' (which are NOT counted in .repeats of the GT container)
    for (int sam_i = num_samples-1; sam_i >= 0; sam_i--) {

        int ht_i=new_ploidy-1 ; for (; ht_i >= vb->ploidy; ht_i--) 
            ht_data[sam_i * new_ploidy + ht_i] = '*'; 

        for (; ht_i >= 0; ht_i--)
            ht_data[sam_i * new_ploidy + ht_i] = ht_data[sam_i * vb->ploidy + ht_i];
    }

    vb->ploidy = new_ploidy;
    vb->ht_per_line = vb->ploidy * vcf_num_samples;
}

WordIndex vcf_seg_FORMAT_GT (VBlockVCFP vb, ContextP ctx, ZipDataLineVCF *dl, STRp(cell), bool has_ps, bool has_null_dp)
{
    // the GT field is represented as a Container, with a single item repeating as required by poidy, and the separator 
    // determined by the phase
    MiniContainer gt = { .repeats   = 1, 
                         .nitems_lo = 1, 
                         .drop_final_repeat_sep = true, 
                         .callback = (vb->use_special_sf == USE_SF_YES || segconf.use_null_DP_method),
                         .items = { { .dict_id = (DictId)_FORMAT_GT_HT } },
                       };

    unsigned save_cell_len = cell_len;

    // update repeats according to ploidy, and separator according to phase
    for (unsigned i=1; i < cell_len-1; i++)
        if (cell[i] == '|' || cell[i] == '/') {
            gt.repeats++;
            gt.repsep[0] = cell[i];
        }

    ASSVCF (gt.repeats <= VCF_MAX_PLOIDY, "ploidy=%u exceeds the maximum of %u", gt.repeats, VCF_MAX_PLOIDY);
    
    // if the ploidy of this line is bigger than the ploidy of the data in this VB so far, then
    // we have to increase ploidy of all the haplotypes read in in this VB so far. This can happen for example in 
    // the X chromosome if initial samples are male with ploidy=1 and then a female sample with ploidy=2
    if (vb->ploidy && gt.repeats > vb->ploidy) 
        vcf_seg_FORMAT_GT_increase_ploidy (vb, gt.repeats);

    if (!vb->ploidy) {
        vb->ploidy = gt.repeats; // very first sample in the vb
        vb->ht_per_line = vb->ploidy * vcf_num_samples;
    }

    buf_alloc (vb, &CTX(FORMAT_GT_HT)->local, vb->ploidy, vb->ht_per_line * vb->lines.len, char, CTX_GROWTH, "contexts->local");

    // note - ploidy of this sample might be smaller than vb->ploidy (eg a male sample in an X chromosesome that was preceded by a female sample, or "." sample)
    Allele *ht_data = ENT (Allele, CTX(FORMAT_GT_HT)->local, vb->line_i * vb->ht_per_line + vb->ploidy * vb->sample_i);

    int64_t dosage=0; // sum of allele values
    for (unsigned ht_i=0; ht_i < gt.repeats; ht_i++) {

        Allele ht = *(cell++); 
        cell_len--;

        ASSVCF (IS_DIGIT(ht) || ht == '.', 
                "invalid VCF file - expecting an allele in a sample to be a number 0-99 or . , but seeing %c (ht_i=%u)", ht, ht_i);

        // single-digit allele numbers
        ht_data[ht_i] = ht;

        // calculate dosage contribution of this ht (to be used in vcf_seg_FORMAT_mux_by_dosage)
        // note: only set for ploidy=1 or 2, and only if GT has 0 or 1, so values are: 0/0->0 ; 0/1->1 ; 1/0->1 ; 1/1->2 ; other->-1
        if (dosage >= 0 && gt.repeats <= 2 && (ht == '0' || ht == '1'))
            dosage += ht - '0'; // dosage only works if alleles are 0 or 1
        else
            dosage = -1; // no dosage
    
        if (!cell_len) break;

        // handle 2-digit allele numbers
        if (ht != '.' && IS_DIGIT (*cell)) {
            unsigned allele = 10 * (ht-'0') + (*(cell++) - '0');
            cell_len--;

            // make sure there isn't a 3rd digit
            ASSVCF (!cell_len || !IS_DIGIT (*cell), "VCF file sample %u - genozip currently supports only alleles up to 99", vb->sample_i+1);

            ht_data[ht_i] = '0' + allele; // use ascii 48->147

            dosage = -1; // no dosage (since allele is not 0 or 1)
        }

        // read and verify phase
        if (gt.repeats > 1 && ht_i < gt.repeats-1) {
            
            char phase = *(cell++);
            cell_len--;

            ASSVCF (phase != ' ', "invalid VCF file - expecting a tab or newline after sample %u but seeing a space", vb->sample_i+1);
            ASSVCF (phase == gt.repsep[0], "invalid VCF file -  unable to parse sample %u: expecting a %c but seeing %c", vb->sample_i+1, gt.repsep[0], phase);
        }
    } // for characters in a sample

    // if the ploidy of the sample is lower than vb->ploidy, set missing ht as '-' (which will cause deletion of themselves and their separator)
    // and set the ploidy to vb->ploidy - to avoid increase in entroy of GT.b250
    if (gt.repeats != vb->ploidy) {
        
        for (unsigned ht_i=gt.repeats; ht_i < vb->ploidy; ht_i++) 
            ht_data[ht_i] = '-'; // unlike '*', we DO count '-' in .repeats (so that we can have the same number of repeats = lower entroy in GT.b250)

        gt.repeats = vb->ploidy;
        if (!gt.repsep[0]) gt.repsep[0] = vb->gt_prev_phase; // this happens in case if a 1-ploid sample
    }

    // case DP='.' - we predict that GT=./.
    bool no_duplicate=false;
    if (has_null_dp) { // segconf.use_null_DP_method=true && line has DP && it is '.'
        // case: prediction is correct - re-write GT as "0/0" or "0|0" (use previous phase)
        if (gt.repeats==2 && ht_data[0]=='.' && gt.repsep[0]=='/' && ht_data[1]=='.') {
            ht_data[0] = ht_data[1] = '0';
            if (vb->gt_prev_phase) gt.repsep[0] = vb->gt_prev_phase;
        }

        // case: prediction is incorrect tell piz to NOT change the HTs and phase according to the prediction.
        // Modifting gt here adds entropy to the GT context, hopefully the prediction is usually correct and this rarely happens 
        else {
            gt.items[0].separator[1] = CI1_ITEM_PRIVATE; // override (hopefully almost never used)
            no_duplicate = true; // prevent seg_duplicate_last - we shouldn't copy the previous sample and the next sample shouldn't copy us
        }
    }

    // case: phase is predictable from has_ps and ht_data[0]
    // note: a case where this prediction fails is with alleles > 1 eg "1|2". In GIAB, this never have a PS, but can be | or / 
    else if (segconf.ps_type && gt.repeats==2 &&
            ((ht_data[0] != '.' && (has_ps == (gt.repsep[0]=='|'))) || // ht!=. --> predicted to be as has_ps says
             (ht_data[0] == '.' && gt.repsep[0]=='/'))) {              // ht==. --> predicted to be /

        // generate the phase from 
        gt.repsep[0] = '&'; // re-write from prediction
        gt.callback = true; // vcf_piz_filter will re-write repsep based on prediction
    }

    // if this sample is a "./." - replace it with "%|%" or "%/%" according to the previous sample's phase -  
    // so that the gt container is likely identical and we reduce GT.b250 entropy. Reason: many tools
    // (including bcftools merge) produce "./." for missing samples even if all other samples are phased
    else if (ht_data[0]=='.' && gt.repeats==2 && ht_data[1]=='.' && gt.repsep[0]=='/') {
        gt.repsep[0] = vb->gt_prev_phase ? vb->gt_prev_phase : '|'; // '|' is arbitrary
        ht_data[0] = ht_data[1] = '%';
    }

    // in case we have INFO/SF, we verify that it is indeed the list of samples for which the first ht is not '.'
    if (vb->use_special_sf == USE_SF_YES && ht_data[0] != '.') 
        vcf_seg_INFO_SF_one_sample (vb);

    ctx_set_last_value (VB, ctx, dosage); // to be used in vcf_seg_FORMAT_mux_by_dosage
    
    if (segconf.running) 
        segconf.count_dosage[dosage >= 0 && dosage <= 2]++;
    
    ASSVCF (!cell_len, "Invalid GT data in sample_i=%u", vb->sample_i+1);
    
    // shortcut if we have the same ploidy and phase as previous GT (saves re-genetrating base64 in container_seg)
    if (gt.repeats == vb->gt_prev_ploidy && gt.repsep[0] == vb->gt_prev_phase && !no_duplicate) 
        return seg_duplicate_last (VB, ctx, save_cell_len);

    else {
        vb->gt_prev_ploidy = no_duplicate ? 0 : gt.repeats; // if no_duplicate - 0 to prevent next sample from duplicating this one
        vb->gt_prev_phase  = gt.repsep[0];
        return container_seg (vb, ctx, (ContainerP)&gt, 0, 0, save_cell_len); 
    }
}

//------------------
// PIZ
//------------------

static inline bool vcf_piz_is_in_FORMAT (VBlockP vb, const char *tag/*2-chars*/)
{
    // get FORMAT field
    const char *format = last_txt (vb, VCF_FORMAT);
    uint32_t format_len = vb->last_txt_len (VCF_FORMAT);

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
    else ASSPIZ (false, "Cannot find '&' predictable phase in recon: \"%.*s\"", recon_len, recon);
}

// called after 2nd HT in a diploid GT, if has_null_DP alg is activated 
void vcf_piz_GT_cb_null_GT_if_null_DP (VBlockP vb , char *recon)
{
    if (recon[0]=='0' && recon[-2]=='0' && recon[-3]=='\t' && vcf_piz_is_in_FORMAT (vb, "DP")) { 

        STR(DP);
        reconstruct_peek (VB, CTX(FORMAT_DP), pSTRa(DP));

        if (DP_len == 1 && *DP == '.') {
            recon[-2] = recon[0] = '.';
            recon[-1] = '/';
        }
    }
}

// Lift-over translator assigned to a FORMAT/GT item, IF it is bi-allelic and we have a ALT<>REF switch. No limitations on ploidy.
// We switch 0<>1. If its unphased (only /) - we list the 0s first, then the 1s
TRANSLATOR_FUNC (vcf_piz_luft_GT)
{
    // validate. make sure this is a bi-allelic genotype (no ploidy limitation)
    for (uint32_t i=0; i < recon_len; i += 2)  
        if (recon[i] != '0' && recon[i] != '1' && recon[i] != '.') return false;

    for (uint32_t i=1; i < recon_len; i += 2)  
        if (recon[i] != '/' && recon[i] != '|') return false;

    VB_VCF->gt_prev_ploidy = (recon_len+1) / 2; // consumed by vcf_piz_luft_PLOIDY

    if (validate_only) return true;

    // exchange 0 <> 1
    for (uint32_t i=0; i < recon_len; i += 2)
        if      (recon[i] == '0') recon[i] = '1';
        else if (recon[i] == '1') recon[i] = '0';

    return true;    
}
