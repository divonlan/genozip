// ------------------------------------------------------------------
//   vcf_seg.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <math.h>
#include "vcf_private.h"
#include "seg.h"
#include "context.h"
#include "random_access.h"
#include "optimize.h"
#include "file.h"
#include "strings.h"
#include "zip.h"
#include "dict_id.h"
#include "codec.h"
#include "reference.h"

#define DATA_LINE(i) ENT (ZipDataLineVCF, vb->lines, i)

static void vcf_seg_complete_missing_lines (VBlockVCF *vb);

// called from seg_all_data_lines
void vcf_seg_initialize (VBlock *vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    vb->contexts[VCF_CHROM] .no_stons   = true; // needs b250 node_index for random access
    vb->contexts[VCF_FORMAT].no_stons   = true;
    vb->contexts[VCF_INFO]  .no_stons   = true;
    vb->contexts[VCF_TOPLEVEL].no_stons = true; // keep in b250 so it can be eliminated as all_the_same

    ctx_get_ctx (vb, dict_id_FORMAT_GT)->no_stons = true; // we store the GT matrix in local, so cannot accomodate singletons

    // room for already existing FORMATs from previous VBs
    vb->format_mapper_buf.len = vb->contexts[VCF_FORMAT].ol_nodes.len;
    buf_alloc (vb, &vb->format_mapper_buf, vb->format_mapper_buf.len * sizeof (Container), 1.2, "format_mapper_buf");
    buf_zero (&vb->format_mapper_buf);

    if (flag.gtshark) codec_gtshark_comp_init (vb_);
    else              codec_hapmat_comp_init  (vb_);
}             

void vcf_seg_finalize (VBlockP vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;
    
    if (vb->ht_matrix_ctx) 
        vcf_seg_complete_missing_lines (vb);

    // top level snip
    SmallContainer top_level = { 
        .repeats   = vb->lines.len,
        .is_toplevel = true,
        .nitems_lo = 10,
        .items     = { { .dict_id = (DictId)dict_id_fields[VCF_CHROM],   .seperator = "\t" },
                       { .dict_id = (DictId)dict_id_fields[VCF_POS],     .seperator = "\t" },
                       { .dict_id = (DictId)dict_id_fields[VCF_ID],      .seperator = "\t" },
                       { .dict_id = (DictId)dict_id_fields[VCF_REFALT],  .seperator = "\t" },
                       { .dict_id = (DictId)dict_id_fields[VCF_QUAL],    .seperator = "\t" },
                       { .dict_id = (DictId)dict_id_fields[VCF_FILTER],  .seperator = "\t" },
                       { .dict_id = (DictId)dict_id_fields[VCF_INFO],    .seperator = "\t" },
                       { .dict_id = (DictId)dict_id_fields[VCF_FORMAT],  .seperator = "\t" },
                       { .dict_id = (DictId)dict_id_fields[VCF_SAMPLES], .seperator = ""   },
                       { .dict_id = (DictId)dict_id_fields[VCF_EOL],     .seperator = ""   } }
    };
    
    container_seg_by_ctx (vb_, &vb->contexts[VCF_TOPLEVEL], (ContainerP)&top_level, 0, 0, 0);

    if (flag.show_alleles && vb->ht_matrix_ctx) {
        printf ("After segmenting (lines=%u samples=%u ploidy=%u len=%u):\n", (uint32_t)vb->lines.len, vcf_num_samples, vb->ploidy, (unsigned)vb->ht_matrix_ctx->local.len);

        for (uint32_t line_i=0; line_i < vb->lines.len; line_i++)
            printf ("Line %-2u: %.*s\n", line_i, vb->num_haplotypes_per_line, ENT (char, vb->ht_matrix_ctx->local, line_i * vb->num_haplotypes_per_line));
    }
}

// optimize REF and ALT, for simple one-character REF/ALT (i.e. mostly a SNP or no-variant)
static void vcf_seg_optimize_ref_alt (VBlockP vb, const char *start_line, char vcf_ref, char vcf_alt)
{
    char new_ref=0, new_alt=0;

    // if we have a reference, we use it
    if (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE) {
        PosType pos = vb->contexts[VCF_POS].last_value.i;

        RefLock lock;
        Range *range = ref_seg_get_locked_range (vb, pos, 1, start_line, &lock);
        uint32_t index_within_range = pos - range->first_pos;

        ref_assert_nucleotide_available (range, pos);
        char ref = ref_get_nucleotide (range, index_within_range);

        if (vcf_ref == ref) new_ref = '-'; // this should always be the case...
        if (vcf_alt == ref) new_alt = '-'; 

        if (flag.reference == REF_EXT_STORE)
            bit_array_set (&range->is_set, index_within_range);

        ref_unlock (lock);
    }

    // replace the most common SNP with +
    // based on counting simple SNPs from from chr22 of 1000 genome project phase 3:
    // G to: A=239681 C=46244 T=44084
    // C to: T=238728 G=46508 A=43685
    // A to: G=111967 C=30006 T=26335
    // T to: C=111539 G=29504 A=25599

    if      (vcf_alt == 'A' && vcf_ref == 'G') new_alt = '+';
    else if (vcf_alt == 'C' && vcf_ref == 'T') new_alt = '+';
    else if (vcf_alt == 'G' && vcf_ref == 'A') new_alt = '+';
    else if (vcf_alt == 'T' && vcf_ref == 'C') new_alt = '+';

    // if anything was done, we create a "special" snip
    if (new_ref || new_alt) {
        char refalt_special[4] = { SNIP_SPECIAL, VCF_SPECIAL_REFALT };
        refalt_special[2] = new_ref ? new_ref : vcf_ref;
        refalt_special[3] = new_alt ? new_alt : vcf_alt;

        seg_by_did_i (vb, refalt_special, sizeof(refalt_special), VCF_REFALT, 4 /* 2 characters and 2 tabs */);
    }
    // if not - just the normal snip
    else {
        char refalt_normal[3] = { 0, '\t', 0 };
        refalt_normal[0] = vcf_ref;
        refalt_normal[2] = vcf_alt;

        seg_by_did_i (vb, refalt_normal, sizeof(refalt_normal), VCF_REFALT, 4 /* 2 characters and 2 tabs */);
    }
}

// traverses the FORMAT field, gets ID of subfield, and moves to the next subfield
static DictId vcf_seg_get_format_subfield (const char **str, uint32_t *len) // remaining length of line 
{
    unsigned i=0; for (; i < *len && (*str)[i] != ':' && (*str)[i] != '\t' && (*str)[i] != '\n'; i++);
    
    DictId dict_id;
    // case: normal field - starts with a letter or another character in the range
    if ((*str)[0] >= 64 && (*str)[0] <= 127) 
        dict_id = dict_id_make (*str, i, DTYPE_VCF_FORMAT);
    
    // case: unusual field - starts with an out-range character, eg a digit - prefix with @ so its a legal FORMAT dict_id
    else {
        SAFE_ASSIGN (*str - 1, '@');
        dict_id = dict_id_make (*str-1, i+1, DTYPE_VCF_FORMAT);
        SAFE_RESTORE;
    }

    *str += i+1;
    *len -= i+1;
    return dict_id; 
}

static void vcf_seg_format_field (VBlockVCF *vb, ZipDataLineVCF *dl, const char *field_start, int field_len)
{
    const char *str = field_start;
    int len = field_len;

    if (!vcf_num_samples) {
        seg_by_did_i_ex (vb, field_start, field_len, VCF_FORMAT, field_len + 1 /* \n */, NULL);
        return; // if we're not expecting any samples, no need to analyze the FORMAT field
    }

    ASSSEG0 (field_len >= 2, field_start, "missing or invalid FORMAT field");

    Container format_mapper = (Container){ 
        .drop_final_repeat_sep = true,
        .filter_items          = true,
        .filter_repeats        = true,
        .repsep                = "\t"
    };

    dl->has_haplotype_data = (str[0] == 'G' && str[1] == 'T' && (str[2] == ':' || field_len==2)); // GT field in FORMAT columns - must always appear first per VCF spec (if at appears)
    dl->has_genotype_data  = (field_len > 2 || (!dl->has_haplotype_data && field_len > 0));

    bool last_item = false;
    do {
        ASSSEG (con_nitems (format_mapper) < MAX_SUBFIELDS, field_start,
                "FORMAT field has too many subfields, the maximum allowed is %u: \"%.*s\"",  
                MAX_SUBFIELDS, field_len, field_start);

        DictId dict_id = vcf_seg_get_format_subfield (&str, (unsigned *)&len);
        last_item = (str[-1] == '\t' || str[-1] == '\n');

        format_mapper.items[con_nitems (format_mapper)] = (ContainerItem) {
            .dict_id   = dict_id,
            .seperator = { last_item ? 0 : ':' },
        };
        con_inc_nitems (format_mapper);

        ASSSEG (dict_id_is_vcf_format_sf (dict_id), field_start,
                "string %.*s in the FORMAT field \"%.*s\" is not a legal subfield", 
                DICT_ID_LEN, dict_id.id, field_len, field_start);

    } 
    while (!last_item && len > 0);
    
    bool is_new;
    char snip[field_len+2];
    snip[0] = SNIP_SPECIAL;
    snip[1] = VCF_SPECIAL_FORMAT;
    memcpy (&snip[2], field_start, field_len);

    uint32_t node_index = seg_by_did_i_ex (vb, snip, field_len+2, VCF_FORMAT, field_len + 1 /* \t or \n */, &is_new);

    dl->format_node_i = node_index;

    if (is_new) {
        ASSERTE (node_index == vb->format_mapper_buf.len, 
                 "node_index=%u different than vb->format_mapper_buf.len=%u", node_index, (uint32_t)vb->format_mapper_buf.len);

        buf_alloc (vb, &vb->format_mapper_buf, (++vb->format_mapper_buf.len) * sizeof (Container), 2, "format_mapper_buf");
    }    

    ContainerP con = ENT (Container, vb->format_mapper_buf, node_index);
    if (is_new || !con_nitems (*con)) // assign if not already assigned. 
        *con = format_mapper; 
}

// return true if caller still needs to seg 
static bool vcf_seg_INFO_DP (VBlockVCF *vb, const char *value, int value_len)
{
    ContextP ctx_dp = ctx_get_ctx (vb, dict_id_INFO_DP);

    // also tried delta vs DP4, but it made it worse
    if (vb->has_basecounts) {
        seg_delta_vs_other ((VBlockP)vb, ctx_dp, ctx_get_existing_ctx (vb, dict_id_INFO_BaseCounts), 
                            value, value_len, -1);
        return false; // caller needn't seg
    }
    else {
        // store last_value of INFO/DP field in case we have FORMAT/DP as well (used in vcf_seg_one_sample)
        ctx_dp->last_value.i = atoi (value);
        return true; // caller should seg
    }
}

// return true if caller still needs to seg 
static bool vcf_seg_INFO_BaseCounts (VBlockVCF *vb, const char *value, int value_len)
{
    Context *ctx_basecounts= ctx_get_ctx (vb, dict_id_INFO_BaseCounts);
    Context *ctx_refalt = &vb->contexts[VCF_REFALT];

    if (ctx_refalt->last_txt_len != 3) return true; // not simple two bases - caller should seg

    char *str = (char *)value;
    ctx_basecounts->last_value.i = 0;

    uint32_t counts[4], sorted_counts[4] = {}; // corresponds to A, C, G, T

    SAFE_ASSIGN (&value[value_len], 0);
    for (unsigned i=0; i < 4; i++) {
        counts[i] = strtoul (str, &str, 10);
        str++; // skip comma seperator
        ctx_basecounts->last_value.i += counts[i];
    }
    SAFE_RESTORE;

    if (str - value != value_len + 1 /* +1 due to final str++ */) return true; // invalid BaseCounts data - caller should seg

    const char *refalt = ENT (const char, vb->txt_data, ctx_refalt->last_txt);

    unsigned ref_i = acgt_encode[(int)refalt[0]];
    unsigned alt_i = acgt_encode[(int)refalt[2]];

    bool used[4] = {};
    sorted_counts[0] = counts[ref_i]; // first - the count of the REF base
    sorted_counts[1] = counts[alt_i]; // second - the count of the ALT base
    used[ref_i] = used[alt_i] = true;

    // finally - the other two cases in their original order (usually these are 0)
    for (unsigned sc_i=2; sc_i <= 3; sc_i++)
        for (unsigned c_i=0; c_i <= 3; c_i++)
            if (!used[c_i]) { // found a non-zero count
                sorted_counts[sc_i] = counts[c_i];
                used[c_i] = true;
                break;
            }

    char snip[2 + value_len + 1]; // +1 for \0
    sprintf (snip, "%c%c%u,%u,%u,%u", SNIP_SPECIAL, VCF_SPECIAL_BaseCounts, 
             sorted_counts[0], sorted_counts[1], sorted_counts[2], sorted_counts[3]);

    seg_by_ctx ((VBlockP)vb, snip, value_len+2, ctx_basecounts, value_len, NULL); 
    
    ctx_basecounts->flags.store = STORE_INT;
    vb->has_basecounts = true;

    return false; // we already segged - caller needn't seg
}

// INFO fields with a format originating from the VEP software, eg
// vep=T|intergenic_variant|MODIFIER|||Intergenic||||||||||||1|||SNV||||||||||||||||||||||||
static inline void vcf_seg_INFO_CSQ (VBlock *vb, DictId dict_id, const char *field, unsigned field_len)
{
    Container con = { .repsep = { ',' }, 
                      .drop_final_repeat_sep = true,
                      .keep_empty_item_sep   = true }; // don't delete the | before an empty item

    Context *vep_ctx = ctx_get_ctx (vb, dict_id);
    Context *sf_ctxs[MAX_SUBFIELDS] = {};

    uint32_t item_i=0;
    const char *item_start = field;
    for (uint32_t i=0; i < field_len+1; i++) {
        if (i == field_len || field[i] == ',' || field[i] == '|') { // end of item
            if (item_i == con_nitems(con)) {
                ASSSEG (!con.repeats, field, 
                        "expecting all repeats of %s to have the same number of items, %u, as the first repeat, but repeat %u (0-based) has more: %.*s", 
                        vep_ctx->name, con_nitems(con), con.repeats, field_len, field);

                ASSSEG (item_i < MIN (126, MAX_SUBFIELDS), field, "exceeded the max number of %s items=%u", 
                        vep_ctx->name, MIN (126, MAX_SUBFIELDS)); // the 126 constraint is just the context naming scheme
                
                char name[8];
                sprintf (name, "%c%c_%.3s", item_i < 63 ? '_' : '`', '@' + (item_i % 63), vep_ctx->name);
                DictId dict_id = dict_id_make (name, 6, DTYPE_VCF_INFO);

                sf_ctxs[item_i] = ctx_get_ctx (vb, dict_id);
                sf_ctxs[item_i]->st_did_i = vep_ctx->did_i;

                con.items[item_i] = (ContainerItem){ .dict_id   = dict_id, 
                                                     .seperator = { field[i]=='|' ? '|' : 0 } };
                con_inc_nitems (con);                                     
            }

            unsigned item_len = &field[i] - item_start;
            seg_by_ctx (vb, item_start, item_len, sf_ctxs[item_i], item_len + (i != field_len), NULL);

            item_i++;
            item_start = &field[i+1];

            if (field[i] != '|') { // end of repeat
                ASSSEG (!con.repeats || item_i == con_nitems(con), field, 
                        "expecting all repeats of %s to have the same number of items, %u, as the first repeat, but repeat %u (0-based) has only %u items: %.*s", 
                        vep_ctx->name, con_nitems(con), con.repeats, item_i, field_len, field);
            
                con.repeats++;
                item_i=0;
            }
        }
    }

    container_seg_by_ctx (vb, vep_ctx, &con, 0, 0, 0);

    // TODO: recover and seg just the field as-is, rather than throw an error, if its not the CSQ format we expect
    // (need to roll back all changes to subfields)
}

static inline bool vcf_seg_test_svlen (VBlockVCF *vb, const char *svlen_str, unsigned svlen_str_len)
{
    int64_t svlen;
    if (!str_get_int (svlen_str, svlen_str_len, &svlen)) return false;

    int64_t last_delta = vb->contexts[VCF_POS].last_delta; // INFO_END is an alias of POS - so the last delta would be between END and POS
    return last_delta == -svlen;
}

static bool vcf_seg_special_info_subfields (VBlockP vb_, DictId dict_id, 
                                            const char **this_value, unsigned *this_value_len, char *optimized_snip)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;
    unsigned optimized_snip_len;

    // Optimize VQSLOD
    if (flag.optimize_VQSLOD && (dict_id.num == dict_id_INFO_VQSLOD) &&
        optimize_float_2_sig_dig (*this_value, *this_value_len, 0, optimized_snip, &optimized_snip_len)) {
        
        vb->vb_data_size -= (int)(*this_value_len) - (int)optimized_snip_len;
        *this_value = optimized_snip;
        *this_value_len = optimized_snip_len;
    }

    // POS and END share the same delta stream - the next POS will be a delta vs this END)
    else if (dict_id.num == dict_id_INFO_END) {
        seg_pos_field ((VBlockP)vb, VCF_POS, VCF_POS, true, *this_value, *this_value_len, 0, *this_value_len); // END is an alias of POS
        return false; // do not add to dictionary/b250 - we already did it
    }

    // if SVLEN is negative, it is expected to be minus the delta between END and POS
    else if (dict_id.num == dict_id_INFO_SVLEN && vcf_seg_test_svlen (vb, *this_value, *this_value_len)) {
        Context *ctx = ctx_get_ctx (vb, dict_id_INFO_SVLEN);
        seg_by_ctx ((VBlockP)vb, (char [2]){ SNIP_SPECIAL, VCF_SPECIAL_SVLEN }, 2, ctx, *this_value_len, NULL);
        return false; // do not add to dictionary/b250 - we already did it
    }

    else if (dict_id.num == dict_id_INFO_BaseCounts)
        return vcf_seg_INFO_BaseCounts (vb, *this_value, *this_value_len);
    
    else if (dict_id.num == dict_id_INFO_DP) 
        return vcf_seg_INFO_DP (vb, *this_value, *this_value_len);

    // in case of AC, AN and AF - we store the values, and we postpone handling AC to the finalization
    else if (dict_id.num == dict_id_INFO_AC) { 
        vb->ac = *this_value; 
        vb->ac_len = *this_value_len; 
        return false; 
    } 
    
    else if (dict_id.num == dict_id_INFO_AN) { 
        vb->an = *this_value; 
        vb->an_len = *this_value_len;
        vb->is_an_before_ac = (vb->ac == NULL);
        Context *ctx = ctx_get_ctx(vb, dict_id);
        ctx->flags.store = STORE_INT;
        ctx->no_stons = true;
    }

    else if (dict_id.num == dict_id_INFO_AF) { 
        vb->af = *this_value; 
        vb->af_len = *this_value_len;
        vb->is_af_before_ac = (vb->ac == NULL);     
        Context *ctx = ctx_get_ctx(vb, dict_id);
        ctx->flags.store = STORE_FLOAT;
        ctx->no_stons = true;
    }

    else if (dict_id.num == dict_id_INFO_CSQ) {
        vcf_seg_INFO_CSQ (vb_, dict_id, *this_value, *this_value_len);
        return false; // caller shouldn't seg because we already did
    }

    //else if (dict_id.num == dict_id_INFO_SF) {
    //   tried: seg_array - with or without delta
    //          storing as bitmap (aligned to 64bits boundary per line) - set bits for SFs used
    //          all resulted in inferior compression to the default
    //}

    else if (dict_id.num == dict_id_INFO_vep ||
             dict_id.num == dict_id_INFO_DP_HIST ||
             dict_id.num == dict_id_INFO_GQ_HIST ||
             dict_id.num == dict_id_INFO_AGE_HISTOGRAM_HET ||
             dict_id.num == dict_id_INFO_AGE_HISTOGRAM_HOM) {

        Context *ctx = ctx_get_ctx (vb, dict_id);
        seg_array (vb_, ctx, ctx->did_i, *this_value, *this_value_len, ',', '|', false);
        return false;
    }

    else if (dict_id.num == dict_id_INFO_DP4) {
        // note: I tried this for dict_id_INFO_SF too, but it made it worse (both with/without delta)
        Context *ctx = ctx_get_ctx (vb, dict_id);
        seg_array (vb_, ctx, ctx->did_i, *this_value, *this_value_len, ',', 0, false);
        return false;
    }

    // finalization of this INFO field
    else if (!dict_id.num) { 

        // if we have AC - we try to represent it as derived from AN and AF, and if not possible, we just
        // store it normally. Note: we can normally do this accurately because AC is an integer. If we tried
        // to instead represent AF as a function of AC and AN we would run into floating point issues
        if (vb->ac) {
            bool special = false;
            if (vb->an && vb->af) {

                char *after_af; double af = strtod (vb->af, &after_af); 
                char *after_an; int    an = strtol (vb->an, &after_an, 10); 
                char *after_ac; int    ac = strtol (vb->ac, &after_ac, 10); 

                if (vb->an + vb->an_len == after_an && // conversion to a number consumed the entire snip
                    vb->af + vb->af_len == after_af && 
                    vb->ac + vb->ac_len == after_ac && 
                    (int)round (af * an) == ac) { 

                    special = true;

                    char special_snip[4] = { SNIP_SPECIAL,
                                             VCF_SPECIAL_AC,
                                             '0' + (char)vb->is_an_before_ac,
                                             '0' + (char)vb->is_af_before_ac };

                    Context *ctx = ctx_get_ctx (vb, dict_id_INFO_AC);
                    ctx->no_stons = true;
                    seg_by_ctx ((VBlockP)vb, special_snip, sizeof(special_snip), ctx, vb->ac_len, NULL);
                }
            }

            if (!special) {
                seg_by_dict_id (vb, vb->ac, vb->ac_len, dict_id_INFO_AC, vb->ac_len);
            }
        }
    }

    return true; // procedue with adding to dictionary/b250
}

static inline WordIndex vcf_seg_FORMAT_PS (VBlockVCF *vb, Context *ctx, const char *cell, unsigned cell_len)
{
    ctx->flags.store = STORE_INT;
    ctx->no_stons = true;

    int64_t ps_value=0;
    if (str_get_int (cell, cell_len, &ps_value) && ps_value == ctx->last_value.i) // same as previous line
        return seg_by_ctx ((VBlockP)vb, (char []){ SNIP_SELF_DELTA, '0' }, 2, ctx, cell_len, NULL);

    return seg_delta_vs_other ((VBlockP)vb, ctx, &vb->contexts[VCF_POS], cell, cell_len, 1000);
}

// used for DP and GQ - store in transposed matrix in local 
static inline WordIndex vcf_seg_FORMAT_transposed (VBlockVCF *vb, Context *ctx, const char *cell, unsigned cell_len, unsigned add_bytes)
{
    ctx->ltype = LT_UINT32_TR;
    ctx->flags.store = STORE_INT;

    buf_alloc_more (vb, &ctx->local, 1, vb->lines.len * vcf_num_samples, uint32_t, 1, "contexts->local");

    if (cell_len == 1 && cell[0] == '.') {
        NEXTENT (uint32_t, ctx->local) = 0xffffffff;
    }
    else {
        ASSERTE (str_get_int (cell, cell_len, &ctx->last_value.i) && ctx->last_value.i >= 0 && ctx->last_value.i <= 0xfffffffe, 
                 "expecting an integer in the range [0, 0xfffffffe] or a '.', but found: %.*s", cell_len, cell);

        NEXTENT (uint32_t, ctx->local) = (uint32_t)ctx->last_value.i;
    }

    // add a LOOKUP to b250
    static const char lookup[1] = { SNIP_LOOKUP };
    seg_by_ctx ((VBlockP)vb, lookup, 1, ctx, add_bytes, NULL);

    return 0;
}

static inline WordIndex vcf_seg_FORMAT_DP (VBlockVCF *vb, Context *ctx, const char *cell, unsigned cell_len,
                                           bool ad_has_sum)
{
    // case: if there is only one sample there is an INFO/DP too, and it is the same - we store a delta of 0 
    Context *info_dp_ctx;
    if (vcf_num_samples == 1 && (info_dp_ctx = ctx_get_existing_ctx (vb, dict_id_INFO_DP))) 
        return seg_delta_vs_other ((VBlockP)vb, ctx, info_dp_ctx, cell, cell_len, -1);

    // case: no valid AD in this sample - store in transposed matrix
    if (!ad_has_sum)
        return vcf_seg_FORMAT_transposed (vb, ctx, cell, cell_len, cell_len); // this handles DP that is an integer or '.'

    return seg_delta_vs_other ((VBlockP)vb, ctx, ctx_get_existing_ctx (vb, dict_id_FORMAT_AD), cell, cell_len, -1);
}

// the DS (allele DoSage) value is usually close to or exactly the sum of '1' alleles in GT. we store it as a delta from that,
// along with the floating point format to allow exact reconstruction
static inline WordIndex vcf_seg_FORMAT_DS (VBlockVCF *vb, Context *ctx, const char *cell, unsigned cell_len)
{
    int64_t dosage = ctx_get_ctx (vb, (DictId)dict_id_FORMAT_GT)->last_value.i; // dosage store here by vcf_seg_FORMAT_GT
    double ds_val;

    if (dosage < 0 || (ds_val = str_get_positive_float (cell, cell_len)) < 0) 
        return seg_by_ctx ((VBlockP)vb, cell, cell_len, ctx, cell_len, NULL);

    char snip[30] = { SNIP_SPECIAL, VCF_SPECIAL_DS }; 
    unsigned snip_len = 2 + str_get_float_format (cell, cell_len, &snip[2]);
    snip[snip_len++] = ' ';
    snip_len += str_int ((int64_t)((ds_val - dosage) * 1000000), &snip[snip_len]);

    return seg_by_ctx ((VBlockP)vb, snip, snip_len, ctx, cell_len, NULL);
}

// increase ploidy of the previous lines, if higher ploidy was encountered
static void vcf_seg_increase_ploidy (VBlockVCF *vb, unsigned new_ploidy, unsigned sample_i, uint32_t max_new_size)
{
    // protect against highly unlikely case that we don't have enough consumed txt data to store increased-ploidy ht data 
    ASSERTE (new_ploidy * vb->line_i * vcf_num_samples <= max_new_size, 
             "haplotype data overflow due to increased ploidy on line %u", vb->line_i);

    uint32_t num_samples = vb->line_i * vcf_num_samples + sample_i; // all samples in previous lines + previous samples in current line
    char *ht_data = vb->ht_matrix_ctx->local.data;

    // copy the haplotypes backwards (to avoid overlap), padding with '*' (which are NOT counted in .repeats of the GT container)
    for (int sam_i = num_samples-1; sam_i >= 0; sam_i--) {

        int ht_i=new_ploidy-1 ; for (; ht_i >= vb->ploidy; ht_i--) 
            ht_data[sam_i * new_ploidy + ht_i] = '*'; 

        for (; ht_i >= 0; ht_i--)
            ht_data[sam_i * new_ploidy + ht_i] = ht_data[sam_i * vb->ploidy + ht_i];
    }

    vb->ploidy = new_ploidy;
    vb->num_haplotypes_per_line = vb->ploidy * vcf_num_samples;
}

static inline WordIndex vcf_seg_FORMAT_GT (VBlockVCF *vb, Context *ctx, ZipDataLineVCF *dl, const char *cell, unsigned cell_len, unsigned sample_i)
{
    // the GT field is represented as a Container, with a single item repeating as required by poidy, and the seperator 
    // determined by the phase
    MiniContainer gt = { .repeats = 1, 
                         .nitems_lo = 1, 
                         .drop_final_repeat_sep = true, 
                         .items = { { .dict_id = (DictId)dict_id_FORMAT_GT_HT } } };

    unsigned save_cell_len = cell_len;

    // update repeats according to ploidy, and separator according to phase
    for (unsigned i=1; i<cell_len-1; i++)
        if (cell[i] == '|' || cell[i] == '/') {
            gt.repeats++;
            gt.repsep[0] = cell[i];
        }

    ASSSEG (gt.repeats <= VCF_MAX_PLOIDY, cell, "ploidy=%u exceeds the maximum of %u", gt.repeats, VCF_MAX_PLOIDY);
    
    // if the ploidy of this line is bigger than the ploidy of the data in this VB so far, then
    // we have to increase ploidy of all the haplotypes read in in this VB so far. This can happen for example in 
    // the X chromosome if initial samples are male with ploidy=1 and then a female sample with ploidy=2
    if (vb->ploidy && gt.repeats > vb->ploidy) 
        vcf_seg_increase_ploidy (vb, gt.repeats, sample_i, (uint32_t)(cell - vb->txt_data.data));

    if (!vb->ploidy) {
        vb->ploidy = gt.repeats; // very first sample in the vb
        vb->num_haplotypes_per_line = vb->ploidy * vcf_num_samples;
    }

    if (sample_i == 0 && vb->line_i == 0) {
        // we overlay on the txt to save memory. since the HT data is by definition a subset of txt, we only overwrite txt
        // areas after we have already consumed them
        buf_set_overlayable (&vb->txt_data);
        buf_overlay ((VBlockP)vb, &vb->ht_matrix_ctx->local, &vb->txt_data, "contexts->local");
    }

    // note - ploidy of this sample might be smaller than vb->ploidy (eg a male sample in an X chromosesome that was preceded by a female sample, or "." sample)
    uint8_t *ht_data = ENT (uint8_t, vb->ht_matrix_ctx->local, vb->line_i * vb->num_haplotypes_per_line + vb->ploidy * sample_i);

    int64_t dosage=0; // sum of allele values
    for (unsigned ht_i=0; ht_i < gt.repeats; ht_i++) {

        uint8_t ht = *(cell++); 
        cell_len--;

        ASSSEG (IS_DIGIT(ht) || ht == '.', cell,
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
            ASSSEG (!cell_len || !IS_DIGIT (*cell), cell, "VCF file sample %u - genozip currently supports only alleles up to 99", sample_i+1);

            ht_data[ht_i] = '0' + allele; // use ascii 48->147

            dosage = -1; // no dosage (since allele is not 0 or 1)
        }

        // read and verify phase
        if (gt.repeats > 1 && ht_i < gt.repeats-1) {
            
            char phase = *(cell++);
            cell_len--;

            ASSSEG (phase != ' ', cell, "invalid VCF file - expecting a tab or newline after sample %u but seeing a space", sample_i+1);
            ASSSEG (phase == gt.repsep[0], cell, "invalid VCF file -  unable to parse sample %u: expecting a %c but seeing %c", sample_i+1, gt.repsep[0], phase);
        }
    } // for characters in a sample

    // if the ploidy of the sample is lower than vb->ploidy, set missing ht as '\b' (backspace) (which will cause deletion of themselves and their separator)
    // and set the ploidy to vb->ploidy - to avoid increase in entroy of GT
    if (gt.repeats != vb->ploidy) {
        
        for (unsigned ht_i=gt.repeats; ht_i < vb->ploidy; ht_i++) 
            ht_data[ht_i] = '\b'; // unlike '*', we DO count '\b' in .repeats

        gt.repeats = vb->ploidy;
        gt.repsep[0] = vb->gt_prev_phase;
    }

    ctx->last_value.i = dosage; // to be used in vcf_seg_FORMAT_DS

    ASSSEG (!cell_len, cell, "Invalid GT data in sample_i=%u", sample_i);

    // shortcut if we have the same ploidy and phase as previous GT (saves re-genetrating base64 in container_seg_by_ctx)
    if (gt.repeats == vb->gt_prev_ploidy && gt.repsep[0] == vb->gt_prev_phase) {

        buf_alloc_more (vb, &ctx->b250, 1, vb->lines.len, uint32_t, CTX_GROWTH, "contexts->b250");

        WordIndex node_index = *LASTENT (uint32_t, ctx->b250);
        NEXTENT (uint32_t, ctx->b250) = node_index;
        ctx->txt_len += save_cell_len;

        return node_index;
    }
    else {
        vb->gt_prev_ploidy = gt.repeats;
        vb->gt_prev_phase  = gt.repsep[0];
        return container_seg_by_ctx ((VBlockP)vb, ctx, (ContainerP)&gt, 0, 0, save_cell_len); 
    }
}

// returns length of the snip, not including the separator
static inline unsigned seg_snip_len_tnc (const char *snip, bool *has_13)
{
    const char *s; for (s=snip ; *s != '\t' && *s != ':' && *s != '\n'; s++);
    
    *has_13 = (*s == '\n' && s > snip && s[-1] == '\r'); // check if we have a Windows-style \r\n line ending

    return s - snip - *has_13;
}

// a comma-separated array - each element goes into its own item context, single repeat (somewhat similar to compound, but 
// intended for simple arrays - just comma separators, no delta between lines or optimizations)
#define MAX_AD_ARRAY_ITEMS 36
static WordIndex vcf_seg_FORMAT_AD (VBlockVCF *vb, Context *ctx, const char *value, int value_len,
                                    bool *ad_has_sum) // optional out
{   
    Container con = seg_initialize_container_array ((VBlockP)vb, ctx->dict_id, false);    

    if (ad_has_sum) {
        ctx->last_value.i = 0; // initialize
        ctx->flags.store = STORE_INT;
        *ad_has_sum = true;
    }
    
    for (con.nitems_lo=0; con.nitems_lo < MAX_AD_ARRAY_ITEMS && value_len > 0; con.nitems_lo++) { // value_len will be -1 after last number

        const char *snip = value;
        for (; value_len && *value != ','; value++, value_len--) {};

        unsigned number_len = (unsigned)(value - snip);

        if (ad_has_sum && *ad_has_sum) {
            int64_t number;
            if (str_get_int (snip, number_len, &number))
                ctx->last_value.i += number;
            else
                *ad_has_sum = false; // not a valid number - we cannot return a sum
        }

        if (con.nitems_lo == MAX_AD_ARRAY_ITEMS-1) // final permitted repeat - take entire remaining string
            number_len += value_len;

        if (value_len > 0) con.items[con.nitems_lo].seperator[0] = ','; 
        
        Context *ctx_item = ctx_get_ctx (vb, con.items[con.nitems_lo].dict_id);
        ctx_item->flags.store = STORE_INT;
        ctx_item->st_did_i = ctx->did_i;

        if (con.nitems_lo == 0)
            // the first value is usually somewhat related to the overall sample depth, therefore values within 
            // a sample are expected to be correlated - so we store it transposed
            vcf_seg_FORMAT_transposed (vb, ctx_item, snip, number_len, number_len + (value_len>0 /* has comma */));
        else
            seg_by_ctx ((VBlockP)vb, snip, number_len, ctx_item, number_len + (value_len>0 /* has comma */), NULL);
        
        value_len--; // skip comma
        value++;
    }

    return container_seg_by_ctx ((VBlockP)vb, ctx, (ContainerP)&con, 0, 0, 0);
}

static void vcf_seg_one_sample (VBlockVCF *vb, ZipDataLineVCF *dl, ContainerP samples, uint32_t sample_i,
                                const char *cell, // beginning of sample, also beginning of first "cell" (subfield)
                                unsigned sample_len,  // not including the \t or \n 
                                bool is_vcf_string,
                                unsigned *num_colons,
                                bool *has_13)
{
    Context *dp_ctx = NULL;
    bool end_of_sample = !sample_len;
    bool ad_has_sum = false; // field AD has a valid last_value

    uint32_t num_items = con_nitems (*samples);
    for (unsigned sf=0; sf < num_items; sf++) { // iterate on the order as in the line

        // move next to the beginning of the subfield data, if there is any
        unsigned cell_len = end_of_sample ? 0 : seg_snip_len_tnc (cell, has_13);
        DictId dict_id = samples->items[sf].dict_id;
        Context *ctx = ctx_get_ctx (vb, dict_id);

        // DP here means we can use it for delta of MIN_DP
        if (cell && ctx->dict_id.num == dict_id_FORMAT_DP) dp_ctx = ctx;

        WordIndex node_index;
        unsigned optimized_snip_len;
        char optimized_snip[OPTIMIZE_MAX_SNIP_LEN];

#       define EVAL_OPTIMIZED { \
            node_index = seg_by_ctx ((VBlockP)vb, optimized_snip, optimized_snip_len, ctx, cell_len, NULL); \
            vb->vb_data_size -= (int)cell_len - (int)optimized_snip_len; \
        }

        if (!cell || !cell_len)
            node_index = seg_by_ctx ((VBlockP)vb, cell, cell_len, ctx, cell_len, NULL);

        // note: cannot use switch bc dict_id_* are variables, not constants

        else if (dict_id.num == dict_id_FORMAT_GT)
            node_index = vcf_seg_FORMAT_GT (vb, ctx, dl, cell, cell_len, sample_i);

        else if (dict_id.num == dict_id_FORMAT_DS && // DS special only works if we also have GT
                 samples->items[0].dict_id.num == dict_id_FORMAT_GT)
            node_index = vcf_seg_FORMAT_DS (vb, ctx, cell, cell_len);

        else if (flag.optimize_PL && dict_id.num == dict_id_FORMAT_PL && 
            optimize_vcf_pl (cell, cell_len, optimized_snip, &optimized_snip_len)) 
            EVAL_OPTIMIZED

        else if (flag.optimize_GL && dict_id.num == dict_id_FORMAT_GL &&
            optimize_vector_2_sig_dig (cell, cell_len, optimized_snip, &optimized_snip_len))
            EVAL_OPTIMIZED
    
        else if (flag.optimize_GP && dict_id.num == dict_id_FORMAT_GP && 
            optimize_vector_2_sig_dig (cell, cell_len, optimized_snip, &optimized_snip_len))
            EVAL_OPTIMIZED

        else if (dict_id.num == dict_id_FORMAT_PL) // not optimized
            node_index = seg_array ((VBlockP)vb, ctx, ctx->did_i, cell, cell_len, ',', 0, false);

        // case: PS ("Phase Set") - might be the same as POS (for example, if set by Whatshap: https://whatshap.readthedocs.io/en/latest/guide.html#features-and-limitations)
        // or might be the same as the previous line
        else if (dict_id.num == dict_id_FORMAT_PS) 
            node_index = vcf_seg_FORMAT_PS (vb, ctx, cell, cell_len);

        else if (dict_id.num == dict_id_FORMAT_GQ) 
            node_index = vcf_seg_FORMAT_transposed (vb, ctx, cell, cell_len, cell_len);
            
        else if (dict_id.num == dict_id_FORMAT_DP) 
            node_index = vcf_seg_FORMAT_DP (vb, ctx, cell, cell_len, ad_has_sum);
            
        // case: MIN_DP - it is slightly smaller and usually equal to DP - we store MIN_DP as the delta DP-MIN_DP
        // note: the delta is vs. the DP field that preceeds MIN_DP - we take the DP as 0 there is no DP that preceeds
        else if (dict_id.num == dict_id_FORMAT_MIN_DP) 
            node_index = seg_delta_vs_other ((VBlockP)vb, ctx, dp_ctx, cell, cell_len, -1);

        else if (dict_id.num == dict_id_FORMAT_AD  || dict_id.num == dict_id_FORMAT_ADALL || 
                 dict_id.num == dict_id_FORMAT_ADF || dict_id.num == dict_id_FORMAT_ADR) 
            node_index = vcf_seg_FORMAT_AD (vb, ctx, cell, cell_len, 
                                            dict_id.num == dict_id_FORMAT_AD ? &ad_has_sum : NULL);

        else // default
            node_index = seg_by_ctx ((VBlockP)vb, cell, cell_len, ctx, cell_len, NULL);
        
        if (node_index != WORD_INDEX_MISSING_SF) 
            cell += cell_len + 1 + *has_13; // skip separator too

        end_of_sample = end_of_sample || cell[-1] != ':'; // a \t or \n encountered

        if (num_colons) *num_colons += !end_of_sample;

        ASSSEG (!end_of_sample || !cell || cell[-1] == '\t' || cell[-1] == '\n', &cell[-1], 
                "Error in vcf_seg_one_sample - end_of_sample and yet separator is %c (ASCII %u) is not \\t or \\n",
                cell[-1], cell[-1]);
    }
    ASSSEG0 (end_of_sample, cell, "More FORMAT subfields data than expected by the specification in the FORMAT field");
}

static const char *vcf_seg_samples (VBlockVCF *vb, ZipDataLineVCF *dl, int32_t *len, const char *next_field, 
                                    bool *has_13)
{
    // Container for samples - we have:
    // - repeats as the number of samples in the line (<= vcf_num_samples)
    // - num_items as the number of FORMAT subfields (inc. GT)

    Container samples = *ENT (Container, vb->format_mapper_buf, dl->format_node_i); // make a copy of the template

    const char *field_start;
    unsigned field_len=0, num_colons=0;

    // 0 or more samples
    for (char separator=0 ; separator != '\n'; samples.repeats++) {

        field_start = next_field;
        next_field = seg_get_next_item (vb, field_start, len, true, true, false, &field_len, &separator, has_13, "sample-subfield");

        ASSSEG (field_len, field_start, "Error: invalid VCF file - expecting sample data for sample # %u, but found a tab character", 
                samples.repeats+1);

        vcf_seg_one_sample (vb, dl, &samples, samples.repeats, field_start, field_len, true, &num_colons, has_13);

        ASSSEG (samples.repeats < vcf_num_samples || separator == '\n', next_field,
                "invalid VCF file - expecting a newline after the last sample (sample #%u)", vcf_num_samples);
    }

    // in some real-world files I encountered have too-short lines due to human errors. we pad them
    if (samples.repeats < vcf_num_samples) {
        ASSERTW (false, "Warning: the number of samples in vb->line_i=%u is %u, different than the VCF column header line which has %u samples",
                 vb->line_i, samples.repeats, vcf_num_samples);

        if (dl->has_haplotype_data) {
            char *ht_data = ENT (char, vb->ht_matrix_ctx->local, vb->line_i * vb->ploidy * vcf_num_samples + vb->ploidy * samples.repeats);
            memset (ht_data, '*', vb->ploidy * (vcf_num_samples - samples.repeats));
        }
    }
    
    container_seg_by_ctx ((VBlockP)vb, &vb->contexts[VCF_SAMPLES], &samples, 0, 0, samples.repeats + num_colons); // account for : and \t \r \n separators

    if (vb->ht_matrix_ctx)
        vb->ht_matrix_ctx->local.len = (vb->line_i+1) * vb->num_haplotypes_per_line;
 
    return next_field;
}

// complete haplotypes of lines that don't have GT, if any line in the vblock does have GT.
// In this case, the haplotype matrix must include the lines without GT too
static void vcf_seg_complete_missing_lines (VBlockVCF *vb)
{
    for (vb->line_i=0; vb->line_i < (uint32_t)vb->lines.len; vb->line_i++) {

        if (vb->ht_matrix_ctx && !DATA_LINE (vb->line_i)->has_haplotype_data) {
            char *ht_data = ENT (char, vb->ht_matrix_ctx->local, vb->line_i * vb->num_haplotypes_per_line);
            memset (ht_data, '*', vb->num_haplotypes_per_line);

            // NOTE: we DONT set dl->has_haplotype_data to true bc downstream we still
            // count this row as having no GT field when analyzing gt data
        }
    }

    vb->ht_matrix_ctx->local.len = vb->lines.len * vb->num_haplotypes_per_line;
}

/* segment a VCF line into its fields:
   fields CHROM->FORMAT are normal contexts
   all samples go into the SAMPLES context, which is a Container
   Haplotype and phase data are stored in a separate buffers + a SNIP_SPECIAL in the GT context 
*/
const char *vcf_seg_txt_line (VBlock *vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    VBlockVCF *vb = (VBlockVCF *)vb_;
    ZipDataLineVCF *dl = DATA_LINE (vb->line_i);
    vb->ac = vb->an = vb->af = NULL;
    vb->has_basecounts = false;

    const char *next_field=field_start_line, *field_start;
    unsigned field_len=0;
    char separator;

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    GET_NEXT_ITEM ("CHROM");
    seg_chrom_field (vb_, field_start, field_len);

    GET_NEXT_ITEM ("POS");
    seg_pos_field (vb_, VCF_POS, VCF_POS, false, field_start, field_len, 0, field_len+1);
    
    // POS <= 0 not expected in a VCF file
    ASSERTW (vb->contexts[VCF_POS].last_value.i > 0, "Warning: invalid POS=%"PRId64" value in vb_i=%u vb_line_i=%u: line will be compressed, but not indexed", 
             vb->contexts[VCF_POS].last_value.i, vb->vblock_i, vb->line_i);
             
    random_access_update_pos (vb_, VCF_POS);

    GET_NEXT_ITEM ("ID");
    seg_id_field (vb_, (DictId)dict_id_fields[VCF_ID], field_start, field_len, true);

    // REF + ALT
    // note: we treat REF+\t+ALT as a single field because REF and ALT are highly corrected, in the case of SNPs:
    // e.g. GG has a probability of 0 and GC has a higher probability than GA.
    GET_NEXT_ITEM ("REF");
    
    unsigned alt_len=0;
    const char *alt_start = next_field;
    next_field = seg_get_next_item (vb, alt_start, &len, false, true, false, &alt_len, &separator, NULL, "ALT");

    // optimize ref/alt in the common case of single-character
    if (field_len == 1 && alt_len == 1) 
        vcf_seg_optimize_ref_alt (vb_, field_start_line, *field_start, *alt_start);
    else
        seg_by_did_i (vb, field_start, field_len + alt_len + 1, VCF_REFALT, field_len + alt_len + 2);

    vb->contexts[VCF_REFALT].last_txt = field_start - vb->txt_data.data; // used by vcf_seg_INFO_BaseCounts
    vb->contexts[VCF_REFALT].last_txt_len = field_len + alt_len + 1;

    SEG_NEXT_ITEM (VCF_QUAL);
    SEG_NEXT_ITEM (VCF_FILTER);

    // INFO
    if (vcf_num_samples)
        GET_NEXT_ITEM (DTF(names)[VCF_INFO]) // pointer to string to allow pointer comparison 
    else
        GET_MAYBE_LAST_ITEM (DTF(names)[VCF_INFO]); // may or may not have a FORMAT field

    seg_info_field (vb_, vcf_seg_special_info_subfields, field_start, field_len);

    bool has_samples = true;
    if (separator != '\n') { // has a FORMAT field

        // FORMAT
        GET_MAYBE_LAST_ITEM ("FORMAT");
        vcf_seg_format_field (vb, dl, field_start, field_len);

        ASSSEG0 (separator == '\n' || dl->has_genotype_data || dl->has_haplotype_data, field_start,
                 "expecting line to end as it has no genotype or haplotype data, but it is not");

        if (separator != '\n') // has samples
            next_field = vcf_seg_samples (vb, dl, &len, next_field, has_13); // All sample columns
        else {
            has_samples = false;
            seg_by_did_i (vb, NULL, 0, VCF_SAMPLES, 0); // case no samples: WORD_INDEX_MISSING_SF
        }
    }

    // case no format or samples
    else {
        has_samples = false;
        seg_by_did_i (vb, NULL, 0, VCF_FORMAT, 0); 
        seg_by_did_i (vb, NULL, 0, VCF_SAMPLES, 0);
    }

    ASSERTW (has_samples || !vcf_num_samples, "Warning: vb->line_i=%u has no samples", vb->line_i);

    SEG_EOL (VCF_EOL, false);

    return next_field;
}
