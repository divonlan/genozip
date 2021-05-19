// ------------------------------------------------------------------
//   vcf_samples.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "vcf_private.h"
#include "seg.h"
#include "context.h"
#include "optimize.h"
#include "file.h"
#include "strings.h"
#include "dict_id.h"
#include "codec.h"
#include "reconstruct.h"

static Container con_FORMAT_AD={}, con_FORMAT_ADALL={}, con_FORMAT_ADF={}, con_FORMAT_ADR={}, 
                 con_FORMAT_F1R2={}, con_FORMAT_F2R1={};

void vcf_format_zip_initialize (void) 
{
    if (!con_FORMAT_AD.repeats) { // not yet initialized (in previous files)
        con_FORMAT_AD    = seg_initialize_container_array (dict_id_FORMAT_AD,    false);    
        con_FORMAT_ADALL = seg_initialize_container_array (dict_id_FORMAT_ADALL, false);    
        con_FORMAT_ADF   = seg_initialize_container_array (dict_id_FORMAT_ADF,   false);    
        con_FORMAT_ADR   = seg_initialize_container_array (dict_id_FORMAT_ADR,   false);    
        con_FORMAT_F1R2  = seg_initialize_container_array (dict_id_FORMAT_F1R2,  false);    
        con_FORMAT_F2R1  = seg_initialize_container_array (dict_id_FORMAT_F2R1,  false);    
    }
}

// used for DP, GQ, A0D and otheres - store in transposed matrix in local 
static inline WordIndex vcf_seg_FORMAT_transposed (VBlockVCF *vb, Context *ctx, const char *cell, unsigned cell_len, unsigned add_bytes)
{
    ctx->ltype = LT_UINT32_TR;
    ctx->flags.store = STORE_INT;

    buf_alloc (vb, &ctx->local, 1, vb->lines.len * vcf_num_samples, uint32_t, 1, "contexts->local");

    if (cell_len == 1 && cell[0] == '.') {
        NEXTENT (uint32_t, ctx->local) = 0xffffffff;
    }
    else {
        ASSSEG (str_get_int (cell, cell_len, &ctx->last_value.i) && ctx->last_value.i >= 0 && ctx->last_value.i <= 0xfffffffe, 
                cell, "While compressing %s expecting an integer in the range [0, 0xfffffffe] or a '.', but found: %.*s", 
                ctx->name, cell_len, cell);

        NEXTENT (uint32_t, ctx->local) = (uint32_t)ctx->last_value.i;
    }

    // add a LOOKUP to b250
    seg_by_ctx (vb, (char []){ SNIP_LOOKUP }, 1, ctx, add_bytes);

    return 0;
}

// a comma-separated array - each element goes into its own item context, single repeat (somewhat similar to compound, but 
// intended for simple arrays - just comma separators, no delta between lines or optimizations)
#define MAX_AD_ARRAY_ITEMS 36 // same as MAX_COMPOUND_COMPONENTS
static WordIndex vcf_seg_FORMAT_A_R_G (VBlockVCF *vb, Context *ctx, Container con /* by value */, const char *value, int value_len)
{   
    int64_t sum = 0; // initialize
    ctx->flags.store  = STORE_INT; // tell container_reconstruct_do to set last_value of container to the sum of its items
    
    bool all_ints = true; 
    for (con.nitems_lo=0; con.nitems_lo < MAX_AD_ARRAY_ITEMS && value_len > 0; con.nitems_lo++) { // value_len will be -1 after last number

        Context *ctx_item = ctx_get_ctx (vb, con.items[con.nitems_lo].dict_id);
        ctx_item->flags.store = STORE_INT;
        ctx_item->st_did_i = ctx->did_i;

        const char *snip = value;
        for (; value_len && *value != ','; value++, value_len--) {};

        unsigned number_len = (unsigned)(value - snip);

        all_ints &= str_get_int_range64 (snip, number_len, 0, 1000000000LL, &ctx_item->last_value.i); // must be non-negative ints
        if (all_ints) 
            sum += ctx_item->last_value.i;

        if (con.nitems_lo == MAX_AD_ARRAY_ITEMS-1) // final permitted repeat - take entire remaining string
            number_len += value_len;

        if (value_len > 0) con.items[con.nitems_lo].seperator[0] = ','; 
        
        if (con.nitems_lo == 0)
            // the first value is usually somewhat related to the overall sample depth, therefore values within 
            // a sample are expected to be correlated - so we store it transposed
            vcf_seg_FORMAT_transposed (vb, ctx_item, snip, number_len, number_len + (value_len>0 /* has comma */));
        else
            seg_by_ctx (vb, snip, number_len, ctx_item, number_len + (value_len>0 /* has comma */));
        
        value_len--; // skip comma
        value++;
    }

    if (all_ints) // a valid value for downstream fields to use
        ctx_set_last_value (vb, ctx, sum);

    ctx->last_txt_len = con.nitems_lo; // Seg only: for use by other fields, eg vcf_seg_FORMAT_F2R1

    return container_seg_by_ctx (vb, ctx, (ContainerP)&con, 0, 0, 0);
}

//----------
// FORMAT/PS
// ---------

static inline WordIndex vcf_seg_FORMAT_PS (VBlockVCF *vb, Context *ctx, const char *cell, unsigned cell_len)
{
    ctx->flags.store = STORE_INT;
    ctx->no_stons = true;

    int64_t ps_value=0;
    if (str_get_int (cell, cell_len, &ps_value) && ps_value == ctx->last_value.i) // same as previous line
        return seg_by_ctx (vb, ((char []){ SNIP_SELF_DELTA, '0' }), 2, ctx, cell_len);

    return seg_delta_vs_other (vb, ctx, &vb->contexts[VCF_POS], cell, cell_len, 1000);
}

//----------
// FORMAT/DP
// ---------

static inline WordIndex vcf_seg_FORMAT_DP (VBlockVCF *vb, Context *ctx, const char *cell, unsigned cell_len)
{
    Context *other_ctx;

    // case - we have FORMAT/AD - calculate delta vs the sum of AD components
    if (ctx_has_value_in_line (vb, dict_id_FORMAT_AD, &other_ctx))
        return seg_delta_vs_other (vb, ctx, other_ctx, cell, cell_len, -1);

    // case: there is only one sample there is an INFO/DP too, we store a delta 
    else if (vcf_num_samples == 1 && ctx_has_value_in_line (vb, dict_id_INFO_DP, &other_ctx)) 
        return seg_delta_vs_other (vb, ctx, other_ctx, cell, cell_len, -1);

    // case: no FORMAT/AD and no INFO/DP - store in transposed matrix
    else 
        return vcf_seg_FORMAT_transposed (vb, ctx, cell, cell_len, cell_len); // this handles DP that is an integer or '.'
}

//----------
// FORMAT/DS
// ---------

// the DS (allele DoSage) value is usually close to or exactly the sum of '1' alleles in GT. we store it as a delta from that,
// along with the floating point format to allow exact reconstruction
static inline WordIndex vcf_seg_FORMAT_DS (VBlockVCF *vb, Context *ctx, const char *cell, unsigned cell_len)
{
    int64_t dosage = vb->gt_ctx->last_value.i; // dosage store here by vcf_seg_FORMAT_GT
    double ds_val;

    if (dosage < 0 || (ds_val = str_get_positive_float (cell, cell_len)) < 0) 
        return seg_by_ctx (vb, cell, cell_len, ctx, cell_len);

    char snip[30] = { SNIP_SPECIAL, VCF_SPECIAL_DS }; 
    unsigned snip_len = 2 + str_get_float_format (cell, cell_len, &snip[2]);
    snip[snip_len++] = ' ';
    snip_len += str_int ((int64_t)((ds_val - dosage) * 1000000), &snip[snip_len]);

    return seg_by_ctx (vb, snip, snip_len, ctx, cell_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_FORMAT_DS)
{
    if (!reconstruct) goto done;

    VBlockVCFP vcf_vb = (VBlockVCFP)vb;

    if (!vcf_vb->gt_ctx) vcf_vb->gt_ctx = ctx_get_existing_ctx (vb, dict_id_FORMAT_GT);

    // we are guaranteed that if we have a special snip, then all values are either '0' or '1';
    char *gt = last_txt (vb, vcf_vb->gt_ctx->did_i);
    unsigned dosage=0;
    for (unsigned i=0; i < vcf_vb->gt_ctx->last_txt_len; i+=2) 
        dosage += gt[i]-'0';

    char float_format[10];
    int32_t val;
    sscanf (snip, "%s %d", float_format, &val); // snip looks like eg: "%5.3f 50000"

    bufprintf (vb, &vb->txt_data, float_format, (double)val / 1000000 + dosage);

done:
    return false; // no new value
}

//----------
// FORMAT/GT
// ---------

// complete haplotypes of lines that don't have GT, if any line in the vblock does have GT.
// In this case, the haplotype matrix must include the lines without GT too
void vcf_seg_FORMAT_GT_complete_missing_lines (VBlockVCF *vb)
{
    for (vb->line_i=0; vb->line_i < (uint32_t)vb->lines.len; vb->line_i++) {

        if (vb->ht_matrix_ctx && !DATA_LINE (vb->line_i)->has_haplotype_data) {
            char *ht_data = ENT (char, vb->ht_matrix_ctx->local, vb->line_i * vb->ht_per_line);
            memset (ht_data, '*', vb->ht_per_line);

            // NOTE: we DONT set dl->has_haplotype_data to true bc downstream we still
            // count this row as having no GT field when analyzing gt data
        }
    }

    vb->ht_matrix_ctx->local.len = vb->lines.len * vb->ht_per_line;
}

// increase ploidy of the previous lines, if higher ploidy was encountered
static void vcf_seg_FORMAT_GT_increase_ploidy (VBlockVCF *vb, unsigned new_ploidy, unsigned sample_i, uint32_t max_new_size)
{
    // protect against highly unlikely case that we don't have enough consumed txt data to store increased-ploidy ht data 
    ASSVCF (new_ploidy * vb->line_i * vcf_num_samples <= max_new_size, 
            "haplotype data overflow due to increased ploidy on line %u", vb->line_i);

    uint32_t num_samples = vb->line_i * vcf_num_samples + sample_i; // all samples in previous lines + previous samples in current line
    char *ht_data = FIRSTENT (char, vb->ht_matrix_ctx->local);

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

static inline WordIndex vcf_seg_FORMAT_GT (VBlockVCF *vb, Context *ctx, ZipDataLineVCF *dl, const char *cell, unsigned cell_len, unsigned sample_i)
{
    vb->gt_ctx = ctx;

    // the GT field is represented as a Container, with a single item repeating as required by poidy, and the seperator 
    // determined by the phase
    MiniContainer gt = { .repeats = 1, 
                         .nitems_lo = 1, 
                         .drop_final_repeat_sep = true, 
                         .callback = (vb->use_special_sf == USE_SF_YES),
                         .items = { { .dict_id = (DictId)dict_id_FORMAT_GT_HT } },
                       };

    unsigned save_cell_len = cell_len;

    // update repeats according to ploidy, and separator according to phase
    for (unsigned i=1; i<cell_len-1; i++)
        if (cell[i] == '|' || cell[i] == '/') {
            gt.repeats++;
            gt.repsep[0] = cell[i];
        }

    ASSVCF (gt.repeats <= VCF_MAX_PLOIDY, "ploidy=%u exceeds the maximum of %u", gt.repeats, VCF_MAX_PLOIDY);
    
    // if the ploidy of this line is bigger than the ploidy of the data in this VB so far, then
    // we have to increase ploidy of all the haplotypes read in in this VB so far. This can happen for example in 
    // the X chromosome if initial samples are male with ploidy=1 and then a female sample with ploidy=2
    if (vb->ploidy && gt.repeats > vb->ploidy) 
        vcf_seg_FORMAT_GT_increase_ploidy (vb, gt.repeats, sample_i, ENTNUM (vb->txt_data, cell));

    if (!vb->ploidy) {
        vb->ploidy = gt.repeats; // very first sample in the vb
        vb->ht_per_line = vb->ploidy * vcf_num_samples;
    }

    if (vb->ht_matrix_ctx->local.type != BUF_OVERLAY) { // first time
        // we overlay on the txt to save memory. since the HT data is by definition a subset of txt, we only overwrite txt
        // areas after we have already consumed them
        buf_set_overlayable (&vb->txt_data);
        buf_overlay ((VBlockP)vb, &vb->ht_matrix_ctx->local, &vb->txt_data, "contexts->local");
    }

    // note - ploidy of this sample might be smaller than vb->ploidy (eg a male sample in an X chromosesome that was preceded by a female sample, or "." sample)
    Allele *ht_data = ENT (Allele, vb->ht_matrix_ctx->local, vb->line_i * vb->ht_per_line + vb->ploidy * sample_i);

    int64_t dosage=0; // sum of allele values
    for (unsigned ht_i=0; ht_i < gt.repeats; ht_i++) {

        Allele ht = *(cell++); 
        cell_len--;

        ASSVCF (IS_DIGIT(ht) || ht == '.', 
                "invalid VCF file - expecting an allele in a sample to be a number 0-99 or . , but seeing %c (ht_i=%u)", ht, ht_i);

        // single-digit allele numbers
        ht_data[ht_i] = ht;

        // calculate dosage contribution of this ht (to be used in vcf_seg_FORMAT_DS)
        if (dosage >= 0 && (ht == '0' || ht == '1'))
            dosage += ht - '0'; // dosage only works if alleles are 0 or 1
        else
            dosage = -1; // no dosage

        if (!cell_len) break;

        // handle 2-digit allele numbers
        if (ht != '.' && IS_DIGIT (*cell)) {
            unsigned allele = 10 * (ht-'0') + (*(cell++) - '0');
            cell_len--;

            // make sure there isn't a 3rd digit
            ASSVCF (!cell_len || !IS_DIGIT (*cell), "VCF file sample %u - genozip currently supports only alleles up to 99", sample_i+1);

            ht_data[ht_i] = '0' + allele; // use ascii 48->147

            dosage = -1; // no dosage (since allele is not 0 or 1)
        }

        // read and verify phase
        if (gt.repeats > 1 && ht_i < gt.repeats-1) {
            
            char phase = *(cell++);
            cell_len--;

            ASSVCF (phase != ' ', "invalid VCF file - expecting a tab or newline after sample %u but seeing a space", sample_i+1);
            ASSVCF (phase == gt.repsep[0], "invalid VCF file -  unable to parse sample %u: expecting a %c but seeing %c", sample_i+1, gt.repsep[0], phase);
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

    // if this sample is a "./." - replace it with "%|%" or "%/%" according to the previous sample's phase -  
    // so that the gt container is likely identical and we reduce GT.b250 entropy. Reason: many tools
    // (including bcftools merge) produce "./." for missing samples even if all other samples are phased
    if (ht_data[0]=='.' && gt.repeats==2 && ht_data[1]=='.' && gt.repsep[0]=='/') {
        gt.repsep[0] = vb->gt_prev_phase ? vb->gt_prev_phase : '|'; // '|' is arbitrary
        ht_data[0] = ht_data[1] = '%';
    }

    // in case we have INFO/SF, we verify that it is indeed the list of samples for which the first ht is not '.'
    if (vb->use_special_sf == USE_SF_YES && ht_data[0] != '.') 
        vcf_seg_INFO_SF_one_sample (vb, sample_i);

    ctx->last_value.i = dosage; // to be used in vcf_seg_FORMAT_DS

    ASSVCF (!cell_len, "Invalid GT data in sample_i=%u", sample_i);

    // shortcut if we have the same ploidy and phase as previous GT (saves re-genetrating base64 in container_seg_by_ctx)
    if (gt.repeats == vb->gt_prev_ploidy && gt.repsep[0] == vb->gt_prev_phase) {

        buf_alloc (vb, &ctx->b250, 1, vb->lines.len, uint32_t, CTX_GROWTH, "contexts->b250");

        WordIndex node_index = *LASTENT (uint32_t, ctx->b250);
        NEXTENT (uint32_t, ctx->b250) = node_index;
        ctx->txt_len += save_cell_len;

        return node_index;
    }
    else {
        vb->gt_prev_ploidy = gt.repeats;
        vb->gt_prev_phase  = gt.repsep[0];
        return container_seg_by_ctx (vb, ctx, (ContainerP)&gt, 0, 0, save_cell_len); 
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

    if (validate_only) return true;

    // exchange 0 <> 1
    for (uint32_t i=0; i < recon_len; i += 2)
        if      (recon[i] == '0') recon[i] = '1';
        else if (recon[i] == '1') recon[i] = '0';

    return true;    
}

//------------
// FORMAT/F2R1
//------------

// F2R1 is expected to be (AD-F1R2) - if it indeed is, we use a special snip
static WordIndex vcf_seg_FORMAT_F2R1 (VBlockVCF *vb, Context *f2r1_ctx, const char *value, int value_len)
{
    Context *ad_ctx, *f1r2_ctx;

    if (!ctx_has_value_in_line (vb, dict_id_FORMAT_AD, &ad_ctx) ||
        !ctx_has_value_in_line (vb, dict_id_FORMAT_F1R2, &f1r2_ctx))
        return vcf_seg_FORMAT_A_R_G (vb, f2r1_ctx, con_FORMAT_F2R1, value, value_len); //can't "goto" into scope of variable length arrays

    unsigned num_items = ad_ctx->last_txt_len;
    const char *f2r1_items[num_items];
    unsigned f2r1_item_lens[num_items];

    if (!num_items || num_items != f1r2_ctx->last_txt_len)  // not the same positive number of items
        goto fallback; 

    if (!str_split (value, value_len, num_items, ',', f2r1_items, f2r1_item_lens, true, NULL))
        goto fallback; // F2R1 doesn't have the correct number of items;

    int64_t f2r1_item_value;
    for (unsigned i=0; i < num_items; i++) {
        if (!str_get_int (f2r1_items[i], f2r1_item_lens[i], &f2r1_item_value)) 
            goto fallback; // not an integer
            
        Context *ad_item_ctx   = ctx_get_ctx (vb, con_FORMAT_AD.items[i].dict_id);
        Context *f1r2_item_ctx = ctx_get_ctx (vb, con_FORMAT_F1R2.items[i].dict_id);
        
        if (f2r1_item_value != ad_item_ctx->last_value.i - f1r2_item_ctx->last_value.i)
            goto fallback; // not expected value
    }

    char special[2] = { SNIP_SPECIAL, VCF_SPECIAL_F2R1 };
    return seg_by_ctx (vb, special, sizeof(special), f2r1_ctx, value_len);

fallback:
    return vcf_seg_FORMAT_A_R_G (vb, f2r1_ctx, con_FORMAT_F2R1, value, value_len);
}

// F2R1 = AD - F1R2
SPECIAL_RECONSTRUCTOR (vcf_piz_special_FORMAT_F2R1)
{
    if (!reconstruct) goto done;

    // number of items in F2R1 is the same as AD - count AD's commas 
    Context *ad_ctx = ctx_get_existing_ctx (vb, dict_id_FORMAT_AD);
    unsigned num_items = str_count_char (last_txt (vb, ad_ctx->did_i), ad_ctx->last_txt_len, ',') + 1;

    for (unsigned i=0; i < num_items; i++) {
        DictId ad_item_dict_id   = { .id = {'A', base36(i), 'D'} }; // see seg_initialize_container_array
        DictId f1r2_item_dict_id = { .id = {'F', base36(i), '1', 'R', '2'} };
        
        Context *ad_item_ctx   = ctx_get_existing_ctx (vb, ad_item_dict_id); 
        Context *f1r2_item_ctx = ctx_get_existing_ctx (vb, f1r2_item_dict_id); 

        char f2r1_item[30];
        unsigned f2r1_item_len = str_int (ad_item_ctx->last_value.i - f1r2_item_ctx->last_value.i, f2r1_item);
        
        RECONSTRUCT (f2r1_item, f2r1_item_len);
        if (i < num_items-1) RECONSTRUCT1 (',');
    }

done:
    return false; // no new value
}

//------------------------------------------------
// FORMAT and INFO - subfields with G, R, R2 types
//------------------------------------------------

// Lift-over ALT<>REF switch translator for bi-allelic multi-value fields: 
// three cases: (1) R1,R2->R2,R1 (2) Ra1,Rb1,Ra2,Rb2->Ra2,Rb2,Ra1,Rb1 (3) G11,G12,G22->G22,G12,G11
// returns true if successful (return value used only if validate_only)
static bool vcf_piz_luft_switch_first_last (VBlockP vb, ContextP ctx, char *recon, int32_t recon_len, 
                                           unsigned num_items, char field_type, bool validate_only)
{
    char copy[recon_len];
    memcpy (copy, recon, recon_len);
    
    const char *items[num_items];
    unsigned item_lens[num_items];
    
    // if item format is inconsistent with VCF header - we won't translate it
    if (!str_split (copy, recon_len, num_items, ',', items, item_lens, true, NULL))
        return false;

    if (validate_only) return true;
    
    vb->txt_data.len -= recon_len;
    
    if (num_items==2 || num_items == 3) {
        RECONSTRUCT_SEP (items[num_items-1], item_lens[num_items-1], ',');

        if (num_items==3)
            RECONSTRUCT_SEP (items[1], item_lens[1], ',');

        RECONSTRUCT (items[0], item_lens[0]);
    }
    else if (num_items == 4) { // Ra1,Rb1,Ra2,Rb2 -> Ra2,Rb2,Ra1,Rb1
        RECONSTRUCT_SEP (items[2], item_lens[2], ',');
        RECONSTRUCT_SEP (items[3], item_lens[3], ',');
        RECONSTRUCT_SEP (items[0], item_lens[0], ',');
        RECONSTRUCT     (items[1], item_lens[1]);
    }
    
    return true;
}

// Lift-over translator assigned to a Number=R item, IF it is bi-allelic and we have a ALT<>REF switch.
// 'R' : We switch between the two comma-separated values.
// 'R2': We switch between the two PAIRS of comma-separated values.
// 'G' : We have 3 values which represent the genotypes REF/REF,REF/ALT,ALT/ALT We switch between the 1st and 3rd value.
// returns true if successful 
TRANSLATOR_FUNC (vcf_piz_luft_R)  { return vcf_piz_luft_switch_first_last (vb, ctx, recon, recon_len, 2, 'R', validate_only); }
TRANSLATOR_FUNC (vcf_piz_luft_R2) { return vcf_piz_luft_switch_first_last (vb, ctx, recon, recon_len, 4, '.', validate_only); }
TRANSLATOR_FUNC (vcf_piz_luft_G)  { return vcf_piz_luft_switch_first_last (vb, ctx, recon, recon_len, 3, 'G', validate_only); }

//------------------------------------------------------------------------
// Validate that ALL subfields in ALL samples can luft-translate as needed
//------------------------------------------------------------------------

// if any context fails luft-translation, returns that context, or if all is good, returns NULL
static inline Context *vcf_seg_validate_luft_trans_one_sample (VBlockVCF *vb, ContextP *ctxs, uint32_t num_items, uint32_t sample_i,
                                                               char *sample, unsigned sample_len)
{
    char *items[num_items]; unsigned item_lens[num_items];
    ASSVCF ((num_items = str_split (sample, sample_len, num_items, ':', (const char **)items, item_lens, false, NULL)), // possibly less than declared
            "Invalid sample_i=%u: \"%.*s\"", sample_i, sample_len, sample);

    for (unsigned i=0; i < num_items; i++) 
        if (needs_translation (ctxs[i]) && item_lens[i]) {
            if ((vb->line_coords == DC_LUFT && !vcf_lo_seg_lift_back_to_primary (vb, ctxs[i], items[i], item_lens[i], NULL, NULL)) ||
                (vb->line_coords == DC_PRIMARY && !(DT_FUNC(vb, translator)[ctxs[i]->luft_trans]((VBlockP)vb, ctxs[i], items[i], item_lens[i], true))))
                return ctxs[i];  // failed translation
        }

    return NULL; // all good
}

// If ALL subfields in ALL samples can luft-translate as required: 1.sets ctx->line_is_luft_trans for all contexts 2.lifted-back if this is a LUFT lne
// if NOT: ctx->line_is_luft_trans=false for all contexts, line is rejects (LO_FORMAT), and keeps samples in their original LUFT or PRIMARY coordinates.
static inline void vcf_seg_validate_luft_trans_all_samples (VBlockVCF *vb, uint32_t num_items, ContextP *ctxs, 
                                                            int32_t len, char *samples_start,
                                                            const char *backup_luft_samples, uint32_t backup_luft_samples_len)
{
    const char *field_start, *next_field = samples_start;
    unsigned field_len=0;

    // initialize optimistically. we will roll back and set to false if ANY subfield in ANY sample fails to translate, and re-seg all samples
    for (unsigned sf_i=0; sf_i < num_items; sf_i++)
        ctxs[sf_i]->line_is_luft_trans = needs_translation (ctxs[sf_i]); 

    // 0 or more samples
    unsigned sample_i=1;
    for (char separator=0 ; separator != '\n'; sample_i++) {

        field_start = next_field;
        next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_SEP, GN_IGNORE, &field_len, &separator, NULL, "sample-subfield");
        ASSVCF (field_len, "Error: invalid VCF file - expecting sample data for sample # %u, but found a tab character", sample_i);

        Context *failed_ctx = vcf_seg_validate_luft_trans_one_sample (vb, ctxs, num_items, sample_i, (char *)field_start, field_len);
        if (failed_ctx) { // some context doesn't luft-translate as required
            vcf_lo_seg_rollback_and_reject (vb, LO_FORMAT, failed_ctx); // rollback liftover

            // make all contexts untranslateable in this line
            for (unsigned i=0; i < num_items; i++)  // iterate on the order as in the line
                ctxs[i]->line_is_luft_trans = false;

            // if this is an untranslatable LUFT-only line, recover the original LUFT-coordinates samples
            if (vb->line_coords == DC_LUFT) 
                memcpy (samples_start, backup_luft_samples, backup_luft_samples_len);
        }
    }
}

// ----------
// One sample
// ----------

// returns length of the snip, not including the separator
static inline unsigned seg_snip_len_tnc (const char *snip, bool *has_13)
{
    const char *s; for (s=snip ; *s != '\t' && *s != ':' && *s != '\n'; s++);
    
    *has_13 = (*s == '\n' && s > snip && s[-1] == '\r'); // check if we have a Windows-style \r\n line ending

    return s - snip - *has_13;
}

// if any context fails luft-translation, returns that context, or if all is good, returns NULL
static inline void vcf_seg_one_sample (VBlockVCF *vb, ZipDataLineVCF *dl, ContextP *ctxs, ContainerP samples, uint32_t sample_i,
                                       char *cell,           // beginning of sample, also beginning of first "cell" (subfield)
                                       unsigned sample_len,  // not including the \t or \n 
                                       unsigned *num_colons,
                                       bool *has_13)
{
    Context *other_ctx;
    bool end_of_sample = !sample_len;
    uint32_t num_items = con_nitems (*samples);

    for (unsigned sf=0; sf < num_items; sf++) { // iterate on the order as in the line

        // move next to the beginning of the subfield data, if there is any
        unsigned cell_len = end_of_sample ? 0 : seg_snip_len_tnc (cell, has_13);
        DictId dict_id = samples->items[sf].dict_id;
        Context *ctx = ctxs[sf];

        WordIndex node_index;
        unsigned opt_snip_len = 0;
        char opt_snip[OPTIMIZE_MAX_SNIP_LEN];

#       define SEG_OPTIMIZED do { \
            node_index = seg_by_ctx (vb, opt_snip, opt_snip_len, ctx, opt_snip_len); \
            vb->vb_data_size -= (int)cell_len - (int)opt_snip_len; \
        } while (0)
        
        // two cases of cell_len=2:  
        // node_index is WORD_INDEX_MISSING_SF if its a premature end of sample ; WORD_INDEX_EMPTY_SF if an empty subfield
        if (!cell || !cell_len)
            node_index = seg_by_ctx (vb, end_of_sample ? NULL : cell, cell_len, ctx, cell_len); 

        // note: cannot use switch bc dict_id_* are variables, not constants

        // ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        else if (dict_id.num == dict_id_FORMAT_GT)
            node_index = vcf_seg_FORMAT_GT (vb, ctx, dl, cell, cell_len, sample_i);

        // ## Allele DoSage
        // Also: ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        else if (dict_id.num == dict_id_FORMAT_DS && // DS special only works if we also have GT
                 samples->items[0].dict_id.num == dict_id_FORMAT_GT)
            node_index = vcf_seg_FORMAT_DS (vb, ctx, cell, cell_len);

        // ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">       
        else if (flag.optimize_PL && dict_id.num == dict_id_FORMAT_PL && 
            optimize_vcf_pl (cell, cell_len, opt_snip, &opt_snip_len)) 
            SEG_OPTIMIZED;

        else if (flag.optimize_GL && dict_id.num == dict_id_FORMAT_GL &&
            optimize_vector_2_sig_dig (cell, cell_len, opt_snip, &opt_snip_len))
            SEG_OPTIMIZED;
    
        else if (flag.optimize_GP && dict_id.num == dict_id_FORMAT_GP && 
            optimize_vector_2_sig_dig (cell, cell_len, opt_snip, &opt_snip_len))
            SEG_OPTIMIZED;

        // ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">       
        else if (dict_id.num == dict_id_FORMAT_PL) // not optimized
            node_index = seg_array ((VBlockP)vb, ctx, ctx->did_i, cell, cell_len, ',', 0, false);

        // case: PS ("Phase Set") - might be the same as POS (for example, if set by Whatshap: https://whatshap.readthedocs.io/en/latest/guide.html#features-and-limitations)
        // or might be the same as the previous line
        else if (dict_id.num == dict_id_FORMAT_PS) 
            node_index = vcf_seg_FORMAT_PS (vb, ctx, cell, cell_len);

        // ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        else if (dict_id.num == dict_id_FORMAT_GQ) 
            node_index = vcf_seg_FORMAT_transposed (vb, ctx, cell, cell_len, cell_len);
            
        // ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
        else if (dict_id.num == dict_id_FORMAT_DP) 
            node_index = vcf_seg_FORMAT_DP (vb, ctx, cell, cell_len);
            
        // ##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
        // case: MIN_DP - it is slightly smaller and usually equal to DP - we store MIN_DP as the delta DP-MIN_DP
        // note: the delta is vs. the DP field that preceeds MIN_DP - we take the DP as 0 there is no DP that preceeds
        else if (dict_id.num == dict_id_FORMAT_MIN_DP && ctx_has_value_in_line (vb, dict_id_FORMAT_DP, &other_ctx)) 
            node_index = seg_delta_vs_other (vb, ctx, other_ctx, cell, cell_len, -1);

        // ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">  
        else if (dict_id.num == dict_id_FORMAT_AD) 
            node_index = vcf_seg_FORMAT_A_R_G (vb, ctx, con_FORMAT_AD, cell, cell_len);

        else if (dict_id.num == dict_id_FORMAT_ADALL) 
            node_index = vcf_seg_FORMAT_A_R_G (vb, ctx, con_FORMAT_ADALL, cell, cell_len);

        else if (dict_id.num == dict_id_FORMAT_ADF) 
            node_index = vcf_seg_FORMAT_A_R_G (vb, ctx, con_FORMAT_ADF, cell, cell_len);

        else if (dict_id.num == dict_id_FORMAT_ADR) 
            node_index = vcf_seg_FORMAT_A_R_G (vb, ctx, con_FORMAT_ADR, cell, cell_len);

        // ##FORMAT=<ID=F1R2,Number=R,Type=Integer,Description="Count of reads in F1R2 pair orientation supporting each allele">
        else if (dict_id.num == dict_id_FORMAT_F1R2) 
            node_index = vcf_seg_FORMAT_A_R_G (vb, ctx, con_FORMAT_F1R2, cell, cell_len);

        // ##FORMAT=<ID=F2R1,Number=R,Type=Integer,Description="Count of reads in F2R1 pair orientation supporting each allele">
        else if (dict_id.num == dict_id_FORMAT_F2R1) 
            node_index = vcf_seg_FORMAT_F2R1 (vb, ctx, cell, cell_len);

        else // default
            node_index = seg_by_ctx (vb, cell, cell_len, ctx, cell_len);

        if (node_index != WORD_INDEX_MISSING_SF) 
            cell += cell_len + 1 + *has_13; // skip separator too

        end_of_sample = end_of_sample || cell[-1] != ':'; // a \t or \n encountered

        if (num_colons) *num_colons += !end_of_sample;

        ASSVCF (!end_of_sample || !cell || cell[-1] == '\t' || cell[-1] == '\n', 
                "Error in vcf_seg_one_sample (sample_i=%u) - end_of_sample and yet separator is '%c' (ASCII %u) is not \\t or \\n (cell[-3]='%c' (%u) cell[-2]='%c' (%u) *has_13=%u cell_len=%u node_index=%d)",
                sample_i, cell[-1], (uint8_t)cell[-1], cell[-3], (uint8_t)cell[-3], cell[-2], (uint8_t)cell[-2], *has_13, cell_len, node_index);
    }

    ASSVCF (end_of_sample, "Sample %u has more subfields than specified in the FORMAT field (%u)", sample_i+1, num_items);
}

//------------
// All samples
//------------

const char *vcf_seg_samples (VBlockVCF *vb, ZipDataLineVCF *dl, int32_t *len, char *next_field, bool *has_13,
                             const char *backup_luft_samples, uint32_t backup_luft_samples_len)
{
    // Container for samples - we have:
    // - repeats as the number of samples in the line (<= vcf_num_samples)
    // - num_items as the number of FORMAT subfields (inc. GT)

    Container samples = *ENT (Container, vb->format_mapper_buf, dl->format_node_i); // make a copy of the template
    ContextP *ctxs = ENT (ContextP, vb->format_contexts, dl->format_node_i * MAX_SUBFIELDS);
    uint32_t num_items = con_nitems (samples);

    // check that all subfields in all samples can be luft-translated as required, or make this a LUFT-only / PRIMARY-only line.
    // Also, if the data is in LUFT coordinates and is indeed translatable, then this lifts-back the samples to PRIMARY coordinates
    if (z_dual_coords && LO_IS_OK (last_ostatus))
        vcf_seg_validate_luft_trans_all_samples (vb, num_items, ctxs, *len, next_field, backup_luft_samples, backup_luft_samples_len);

    const char *field_start;
    unsigned field_len=0, num_colons=0;

    // 0 or more samples
    for (char separator=0 ; separator != '\n'; samples.repeats++) {

        field_start = next_field;
        next_field = (char *)seg_get_next_item (vb, field_start, len, GN_SEP, GN_SEP, GN_IGNORE, &field_len, &separator, has_13, "sample-subfield");

        ASSVCF (field_len, "Error: invalid VCF file - expecting sample data for sample # %u, but found a tab character", 
                samples.repeats+1);

        vcf_seg_one_sample (vb, dl, ctxs, &samples, samples.repeats, (char *)field_start, field_len, &num_colons, has_13);

        ASSVCF (samples.repeats < vcf_num_samples || separator == '\n',
                "invalid VCF file - expecting a newline after the last sample (sample #%u)", vcf_num_samples);
    }

    ASSVCF (samples.repeats <= vcf_num_samples, "according the VCF header, there should be %u sample%s per line, but this line has %u samples - that's too many",
            vcf_num_samples, vcf_num_samples==1 ? "" : "s", samples.repeats);

    // in some real-world files I encountered have too-short lines due to human errors. we pad them
    if (samples.repeats < vcf_num_samples) {
        WARN_ONCE ("FYI: the number of samples in variant CHROM=%.*s POS=%"PRId64" is %u, different than the VCF column header line which has %u samples",
                   vb->chrom_name_len, vb->chrom_name, vb->last_int (VCF_POS), samples.repeats, vcf_num_samples);

        if (dl->has_haplotype_data) {
            char *ht_data = ENT (char, vb->ht_matrix_ctx->local, vb->line_i * vb->ploidy * vcf_num_samples + vb->ploidy * samples.repeats);
            unsigned num_missing = vb->ploidy * (vcf_num_samples - samples.repeats); 
            memset (ht_data, '*', num_missing);
        }
    }
    
    // assign all translators. note: we either have translators for all translatable items, or none at all.
    for (uint32_t i=0; i < num_items; i++)
        if (ctxs[i]->line_is_luft_trans)
            samples.items[i].translator = ctxs[i]->luft_trans;

    container_seg_by_ctx (vb, &vb->contexts[VCF_SAMPLES], &samples, 0, 0, samples.repeats + num_colons); // account for : and \t \r \n separators

    if (vb->ht_matrix_ctx)
        vb->ht_matrix_ctx->local.len = (vb->line_i+1) * vb->ht_per_line;
 
    return next_field;
}

