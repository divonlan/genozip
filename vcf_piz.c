// ------------------------------------------------------------------
//   vcf_piz.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <math.h>
#include "vcf_private.h"
#include "zfile.h"
#include "txtfile.h"
#include "seg.h"
#include "context.h"
#include "file.h"
#include "endianness.h"
#include "sections.h"
#include "random_access.h"
#include "dict_id.h"
#include "strings.h"
#include "codec.h"
#include "reference.h"
#include "reconstruct.h"

// returns true if section is to be skipped reading / uncompressing
bool vcf_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    if (flag.drop_genotypes && // note: if all samples are filtered out with --samples then flag.drop_genotypes=true (set in vcf_samples_analyze_field_name_line)
        (dict_id.num == dict_id_fields[VCF_FORMAT] || dict_id.num == dict_id_fields[VCF_SAMPLES] || dict_id_is_vcf_format_sf (dict_id)))
        return true;

    if (flag.gt_only && dict_id_is_vcf_format_sf (dict_id) 
        && dict_id.num != dict_id_FORMAT_GT
        && dict_id.num != dict_id_FORMAT_GT_HT
        && dict_id.num != dict_id_FORMAT_GT_HT_INDEX
        && dict_id.num != dict_id_FORMAT_GT_SHARK_DB
        && dict_id.num != dict_id_FORMAT_GT_SHARK_GT
        && dict_id.num != dict_id_FORMAT_GT_SHARK_EX
        && dict_id.num != dict_id_PBWT_RUNS
        && dict_id.num != dict_id_PBWT_FGRC)
        return true;

    return false;
}

// filter is called before reconstruction of a repeat or an item, and returns false if item should 
// not be reconstructed. contexts are not consumed.
CONTAINER_FILTER_FUNC (vcf_piz_filter)
{
    if (flag.gt_only && dict_id.num == dict_id_fields[VCF_SAMPLES] && item >= 0
        && con->items[item].dict_id.num != dict_id_FORMAT_GT) 
        return false; 

    // in dual-coordinates files - get the oSTATUS at the beginning of each line
    else if (item >= 0 && con->items[item].dict_id.num == dict_id_fields[VCF_oSTATUS]) {
        if (z_file->z_flags.dual_coords)
            *reconstruct = false; // set last_index, without reconstructing
        else
            return false; // filter out entirely without consuming (non-DC files have no oStatus data)
    }

    // for dual-coordinates genozip files - in which rows capable of liftover have both 
    // INFO/LIFTOVER and INFO/LIFTBACK subfields - we select which one to filter out based on flag.luft
    else if (item >= 0 && con->items[item].dict_id.num == dict_id_INFO_LIFTOVER && flag.luft)
        return false;
        
    else if (item >= 0 && con->items[item].dict_id.num == dict_id_INFO_LIFTBACK && !flag.luft)
        return false;
        
    return true;    
}

// FORMAT - obey GT-only and drop-genotypes ; load haplotype line
SPECIAL_RECONSTRUCTOR (vcf_piz_special_FORMAT)
{
    VBlockVCF *vb_vcf = (VBlockVCF *)vb;

    if (flag.drop_genotypes) goto done;

    bool has_GT = (snip_len>=2 && snip[0]=='G' && snip[1] == 'T' && (snip_len==2 || snip[2] == ':'));
    
    if (reconstruct) {
        if (flag.gt_only) {
            if (has_GT)
                RECONSTRUCT ("GT\t", 3);
        }
        else 
            RECONSTRUCT (snip, snip_len);
    }

    // initialize haplotype stuff
    if (has_GT && !vb_vcf->ht_matrix_ctx) {

        ASSERT ((vb_vcf->ht_matrix_ctx = ctx_get_existing_ctx (vb, dict_id_FORMAT_GT_HT)), 
                "vb_i=%u: cannot find GT_HT data", vb->vblock_i);

        // will exist in case of use of CODEC_HAPMAT but not CODEC_GTSHARK        
        vb_vcf->hapmat_index_ctx = ctx_get_existing_ctx (vb, dict_id_FORMAT_GT_HT_INDEX);
        
        if (vb_vcf->hapmat_index_ctx)
            codec_hapmat_piz_calculate_columns (vb);
    }
done:
    return false; // no new value
}

// REFALT - reconstruct from reference and/or common SNPs
SPECIAL_RECONSTRUCTOR (vcf_piz_special_REFALT)
{
    if (!reconstruct) goto done;

    ASSERT (snip_len==2, "expecting snip_len=2 but seeing %u", snip_len);

    // snip is 3 characters - REF, \t, ALT
    char ref_alt[3] = { 0, '\t', 0 };
    char ref_value = 0;
    
    if (snip[0] == '-' || snip[1] == '-') { 
        PosType pos = vb->last_int (VCF_POS);

        const Range *range = ref_piz_get_range (vb, pos, 1);
        ASSERT (range, "failed to find range for chrom='%s' pos=%"PRId64, vb->chrom_name, pos);
        
        uint32_t idx = pos - range->first_pos;
        ASSERT (ref_is_nucleotide_set (range, idx), "reference is not set: chrom=%.*s pos=%"PRId64, range->chrom_name_len, range->chrom_name, pos);
        ref_value = ref_get_nucleotide (range, idx);
    }

    // recover ref
    if (snip[0] == '-') 
        ref_alt[0] = ref_value;
    else 
        ref_alt[0] = snip[0];

    // recover alt
    if (snip[1] == '+') { // the alt has the most common value for a SNP
        if      (ref_alt[0] == 'A') ref_alt[2] = 'G';
        else if (ref_alt[0] == 'C') ref_alt[2] = 'T';
        else if (ref_alt[0] == 'G') ref_alt[2] = 'A';
        else if (ref_alt[0] == 'T') ref_alt[2] = 'C';
    }
    else if (snip[1] == '-')  // the alt has the reference value
        ref_alt[2] = ref_value;

    else // the alt is specified verbatim
        ref_alt[2] = snip[1];

    RECONSTRUCT (ref_alt, sizeof (ref_alt));

done:
    return false; // no new value
}   

SPECIAL_RECONSTRUCTOR (vcf_piz_special_AC)
{
    if (!reconstruct) goto done;
    
    bool is_an_before_ac = (bool)(snip[0] - '0');
    bool is_af_before_ac = (bool)(snip[1] - '0');

    Context *ctx_an = ctx_get_existing_ctx (vb, dict_id_INFO_AN);
    Context *ctx_af = ctx_get_existing_ctx (vb, dict_id_INFO_AF);

    uint32_t an = is_an_before_ac ? ctx_an->last_value.i : ctx_peek_next_int (vb, ctx_an);
    double   af = is_af_before_ac ? ctx_af->last_value.f : ctx_peek_next_float (vb, ctx_af);

    char ac_str[30];
    unsigned ac_str_len = str_int ((int64_t)round(an * af), ac_str);    

    RECONSTRUCT (ac_str, ac_str_len);

done:
    return false; // no new value
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_DS)
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

SPECIAL_RECONSTRUCTOR (vcf_piz_special_BaseCounts)
{
    Context *ctx_refalt = &vb->contexts[VCF_REFALT];
    ASSERT (ctx_refalt->last_txt_len == 3, "Expecting ctx_refalt->last_txt_len=%u to be 3", ctx_refalt->last_txt_len);
    const char *refalt = last_txt (vb, VCF_REFALT);

    uint32_t counts[4], sorted_counts[4] = {}; // counts of A, C, G, T

    new_value->i = 0;
    char *str = (char *)snip;

    for (unsigned i=0; i < 4; i++) {
        sorted_counts[i] = strtoul (str, &str, 10);
        str++; // skip comma seperator
        new_value->i += sorted_counts[i];
    }

    if (!reconstruct) goto done; // just return the new value

    ASSERT (str - snip == snip_len + 1, "expecting (str-snip)=%d == (snip_len+1)=%u", (int)(str - snip), snip_len+1);

    unsigned ref_i = acgt_encode[(int)refalt[0]];
    unsigned alt_i = acgt_encode[(int)refalt[2]];
    
    counts[ref_i] = sorted_counts[0];
    counts[alt_i] = sorted_counts[1];
    
    unsigned sc_i=2;
    for (unsigned i=0; i <= 3; i++)
        if (ref_i != i && alt_i != i) counts[i] = sorted_counts[sc_i++];

    bufprintf (vb, &vb->txt_data, "%u,%u,%u,%u", counts[0], counts[1], counts[2], counts[3]);

done:
    return true; // has new value
}

#define adjustment vcf_vb->sf_ctx->last_delta
#define sample_i   vcf_vb->sf_ctx->last_value.i
#define snip_i     vcf_vb->sf_snip.param

// leave space for reconstructing SF - actual reconstruction will be in vcf_piz_container_cb
SPECIAL_RECONSTRUCTOR (vcf_piz_special_SF)
{
    VBlockVCFP vcf_vb = (VBlockVCFP)vb;

    if (reconstruct) {
        vcf_vb->sf_ctx = ctx;
        adjustment    = 0;
        sample_i      = 0; 
        snip_i        = 0;

        // temporary place for SF
        buf_alloc_old (vb, &vcf_vb->sf_txt, 5 * vcf_header_get_num_samples(), 1, "sf_txt" ); // initial estimate, we may further grow it later
        vcf_vb->sf_txt.len = 0;

        // copy snip to sf_snip (note: the SNIP_SPECIAL+code are already removed)
        buf_alloc_old (vb, &vcf_vb->sf_snip, snip_len, 2, "sf_snip");
        vcf_vb->sf_snip.len = snip_len; 
        memcpy (vcf_vb->sf_snip.data, snip, snip_len); 
    }

    return false; // no new value
}

// the case where SVLEN is minus the delta between END and POS
SPECIAL_RECONSTRUCTOR (vcf_piz_special_SVLEN)
{
    if (!reconstruct) goto done;

    int64_t value = -vb->contexts[VCF_POS].last_delta; // END is a alias of POS - they share the same data stream - so last_delta would be the delta between END and POS
    char str[30];
    unsigned str_len = str_int (value, str);
    RECONSTRUCT (str, str_len);

done:
    return false; // no new value
}

// While reconstructing the GT fields of the samples - calculate the INFO/SF field
static void inline vcf_piz_GT_cb_calc_INFO_SF (VBlockVCFP vcf_vb, unsigned rep, char *reconstructed, int32_t reconstructed_len)
{
    if (rep != 0) return; // we only look at the first ht in a sample, and only if its not '.'/'%'

    if (*reconstructed == '.' || *reconstructed == '%') { // . can be written as % in vcf_seg_FORMAT_GT
        sample_i++;
        return;
    }

    ARRAY (const char, sf_snip, vcf_vb->sf_snip);

    while (snip_i < sf_snip_len) {

        buf_alloc (vcf_vb, &vcf_vb->sf_txt, 12, 0, char, 2, "sf_txt"); // sufficient for int32 + ','

        int32_t adjusted_sample_i = (int32_t)(sample_i + adjustment); // last_delta is the number of values in SF that are not in samples

        // case: we derive the value from sample_i
        if (sf_snip[snip_i] == ',') {
            snip_i++;

            buf_add_int ((VBlockP)vcf_vb, &vcf_vb->sf_txt, adjusted_sample_i);

            if (snip_i < sf_snip_len) // add comma if not done yet
                NEXTENT (char, vcf_vb->sf_txt) = ',';

            break;
        }

        // case: sample was not used to derive any value
        if (sf_snip[snip_i] == '~') {
            snip_i++;
            break;
        }

        // case: value quoted verbatim - output and continue searching for , or ~ for this sample
        else {
            char *comma = strchr (&sf_snip[snip_i], ',');
            unsigned len = comma - &sf_snip[snip_i];
            bool add_comma = (snip_i + len + 1 < sf_snip_len); // add comma if not last

            buf_add (&vcf_vb->sf_txt, &sf_snip[snip_i], len + add_comma);
            snip_i += len + 1;

            adjustment++;
        }
    }

    sample_i++;
}

// Upon completing the line - insert the calculated INFO/SF field to its place
static void inline vcf_piz_TOPLEVEL_cb_insert_INFO_SF (VBlockVCFP vcf_vb)
{
    ARRAY (const char, sf_snip, vcf_vb->sf_snip);

    // if there are some items remaining in the snip (values that don't appear in samples) - copy them
    if (snip_i < sf_snip_len) {
        unsigned remaining_len = sf_snip_len - snip_i - 1; // all except the final comma
        buf_add_more (vcf_vb, &vcf_vb->sf_txt, &sf_snip[snip_i], remaining_len, "sf_txt"); 
    }

    // make room for the SF txt and copy it to its final location
    char *sf_txt = last_txt (vcf_vb, vcf_vb->sf_ctx->did_i);
    memmove (sf_txt + vcf_vb->sf_txt.len, sf_txt, AFTERENT (char, vcf_vb->txt_data) - sf_txt); // make room
    memcpy (sf_txt, vcf_vb->sf_txt.data, vcf_vb->sf_txt.len); // copy

    vcf_vb->txt_data.len += vcf_vb->sf_txt.len;

    buf_free (&vcf_vb->sf_snip);
    buf_free (&vcf_vb->sf_txt);
}

static void inline vcf_piz_SAMPLES_subset_samples (VBlockVCFP vb, unsigned rep, int32_t reconstructed_len)
{
    if (!samples_am_i_included (rep))
        vb->txt_data.len -= reconstructed_len;
}

// ---------------------------------------------------------------------------------
// Translators and Special functions for Luft (secondary coordinates) reconstruction
// Invoked by genocat --luft
// ---------------------------------------------------------------------------------

// Translator called for CHROM (main field) if genocat --luft. if LO_OK, replaces CHROM with oCHROM, and saves CHROM in vb->liftover.
TRANSLATOR_FUNC (vcf_piz_luft_CHROM)
{
    if (vb->is_rejects_vb) return 0;

    // note: we have oSTATUS snips for all lines, but oCHROM, oPOS, oREF snips only for the lines with LO_OK
    if (vb->last_index (VCF_oSTATUS) != LO_OK) return 0;

    // save the primary-coord chrom in liftback (this is part of INFO/LIFTBACK used in vcf_piz_special_LIFTBACK)
    Buffer *liftback = &((VBlockVCF *)vb)->liftover;
    liftback->len = 0;
    buf_alloc_old (vb, liftback, reconstructed_len + 100, 0, "liftover"); // enough for POS (9 chars), and usually also REFALT
    bufprintf (vb, liftback, "%.*s,", reconstructed_len, reconstructed);

    // delete the primary-coord chrom
    vb->txt_data.len -= reconstructed_len;

    reconstruct_from_ctx (vb, VCF_oCHROM, 0, true);
    return 0;    
}

// Translator called for POS (main field) if genocat --luft. if LO_OK, replaces POS with oPOS, and saves POS in vb->liftover.
// (if not LO_OK, line will be dropped by vcf_piz_container_cb)
TRANSLATOR_FUNC (vcf_piz_luft_POS)
{
    if (vb->is_rejects_vb || vb->last_index (VCF_oSTATUS) != LO_OK) return 0;

    // save the primary-coord chrom in liftback (this is part of INFO/LIFTBACK used in vcf_piz_special_LIFTBACK)
    bufprintf (vb, &((VBlockVCF *)vb)->liftover, "%.*s,", reconstructed_len, reconstructed);

    // delete the primary-coord pos
    vb->txt_data.len -= reconstructed_len;

    reconstruct_from_ctx (vb, VCF_oPOS, 0, true);
    
    return 0;    
}

// Translator called for REFALT (main field) if genocat --luft. if LO_OK, replaces REFALT with oREF + calculated oALT
// and saves REF in vb->liftover (actually liftback).
// -f REF != oREF and oREF=ALT, swaps REF and ALT, inc. in nagative strand. vb->liftover updated.
// This is the only form of REF change supported by genozip currently, other changes will 
// be set by Seg as LO_UNSUPPORTED. 
TRANSLATOR_FUNC (vcf_piz_luft_REFALT)
{
    if (vb->is_rejects_vb || vb->last_index (VCF_oSTATUS) != LO_OK) return 0;

    // save the primary-coord chrom in liftback (this is part of INFO/LIFTBACK used in vcf_piz_special_LIFTBACK)
    Buffer *liftback = &((VBlockVCF *)vb)->liftover;
    buf_add_more (vb, liftback, reconstructed, reconstructed_len, "liftover");

    // split REFALT
    const char *ref_alt[2]; unsigned ref_alt_lens[2];
    ASSERT (str_split (reconstructed, reconstructed_len, 2, '\t', ref_alt, ref_alt_lens), 
            "expecting one tab in the REFALT snip: \"%.*s\"", reconstructed_len, reconstructed);

    liftback->len -= ref_alt_lens[1] + 1; // remove \tALT from INFO/LIFTBACK data (but data is still in the buffer...)

    // reconstruct oREF, overwriting REFALT
    uint64_t txt_data_len_with_refalt = vb->txt_data.len;    
    vb->txt_data.len -= reconstructed_len;

    reconstruct_from_ctx (vb, VCF_oREF, '\t', true); // reconstructs special snip using vcf_piz_special_oREF

    // case: REF unchanged - oALT is ALT - just restore the length
    if (vb->last_int (VCF_oREF) == 0) // oREF value stored in vcf_piz_special_oREF()
        vb->txt_data.len = txt_data_len_with_refalt;

// TODO - the old ALT to be rev complemented in strand=-

    // case: REF and ALT were switched - oALT is REF
    else 
// TODO - the new ALT to be rev complemented in strand=-
        RECONSTRUCT (ref_alt[0], ref_alt_lens[0]);

    return 0;    
}

// used for:
// 1. Primary VCF - reconstruct oREF in INFO/LIFTOVER (reconstruction invoked from LIFTOVER container)
// 2. Luft VCF    - reconstruct oREF in primary REF field (reconstruction invoked by vcf_piz_luft_REFALT translator)
SPECIAL_RECONSTRUCTOR (vcf_piz_special_oREF)
{
    // make a copy of refalt, as we will be overwriting it (since call from from vcf_piz_luft_REFALT)
    unsigned ref_alt_str_len = vb->last_txt_len (VCF_REFALT);

    char ref_alt_str[ref_alt_str_len];
    memcpy (ref_alt_str, last_txt (vb, VCF_REFALT), ref_alt_str_len);
    const char *ref_alt[2]; unsigned ref_alt_len[2];
    ASSERT (str_split (ref_alt_str, ref_alt_str_len, 2, '\t', ref_alt, ref_alt_len),
            "expecting one tab in the REFALT snip: \"%.*s\"", ref_alt_str_len, ref_alt_str);

    if (snip_len==1 && *snip=='0')      // ALT
        RECONSTRUCT (ref_alt[1], ref_alt_len[1]);

    else if (snip_len==1 && *snip=='1') // REF
        RECONSTRUCT (ref_alt[0], ref_alt_len[0]);

    else
        ABORT ("Invalid OREF special snip: %.*s", snip_len, snip);

    new_value->i = *snip - '0';
    return true; // new_value was set
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_LIFTREJT)
{
    ctx_get_snip_by_word_index (&vb->contexts[VCF_oSTATUS], vb->last_index(VCF_oSTATUS), &snip, &snip_len);
    RECONSTRUCT (snip, snip_len);
    return false; // new_value was not set
}

// reconstruct LIFTBACK - from CHROM,POS,REF previously reconstructed and stored in vb->liftover, and oSTRAND, oALTRULE, oSTATUS
SPECIAL_RECONSTRUCTOR (vcf_piz_special_LIFTBACK)
{
    ARRAY (char, liftover, ((VBlockVCF *)vb)->liftover);
    
    RECONSTRUCT_SEP (liftover, liftover_len, ',');            
    reconstruct_from_ctx (vb, VCF_oSTRAND, ',', true);
    reconstruct_from_ctx (vb, VCF_oALTRULE,  0, true);

    return false; // new_value was not set
}

// Callback at the end of reconstructing a VCF line with --luft, i.e. with luft coordinates:
// we drop lines that failed liftovers. Seperately, these lines will be preserved in the VCF header.
void vcf_piz_TOPLEVEL_cb_drop_line_if_bad_oSTATUS_or_no_header (VBlock *vb)
{
    // conditions for dropping a line in --luft
    if ((!vb->is_rejects_vb && vb->last_index (VCF_oSTATUS) != LO_OK) || // drop primary lines that are rejected (note: for sorted files, rejects are already excluded from the reconstruction plan in sorter_zip_merge_vb_do)
        ( vb->is_rejects_vb && (flag.no_header || flag.header_one)))     // drop rejects in --no-header or --header-one 
        vb->dont_show_curr_line = true;
}

TRANSLATOR_FUNC (vcf_piz_luft_AC)
{
    
    return 0;    
}

TRANSLATOR_FUNC (vcf_piz_luft_AF)
{
    
    return 0;    
}

TRANSLATOR_FUNC (vcf_piz_luft_AD)
{
    
    return 0;    
}

TRANSLATOR_FUNC (vcf_piz_luft_END)
{
    
    return 0;    
}

TRANSLATOR_FUNC (vcf_piz_luft_GT)
{
    
    return 0;    
}

TRANSLATOR_FUNC (vcf_piz_luft_GL)
{
    return 0;    
}

// ------------------------------------
// callback called after reconstruction
// ------------------------------------

CONTAINER_CALLBACK (vcf_piz_container_cb)
{
    VBlockVCFP vcf_vb = (VBlockVCFP)vb;
    bool have_INFO_SF = vcf_vb->sf_snip.len > 0;

    // case: we have an INFO/SF field and we reconstructed and we reconstructed a repeat (i.e. one ht) GT field of a sample 
    if (dict_id.num == dict_id_FORMAT_GT && have_INFO_SF) 
        vcf_piz_GT_cb_calc_INFO_SF (vcf_vb, rep, reconstructed, reconstructed_len);

    else if (dict_id.num == dict_id_fields[VCF_TOPLEVEL]) {

        // case: we have an INFO/SF field and we reconstructed one VCF line
        if (have_INFO_SF) 
            vcf_piz_TOPLEVEL_cb_insert_INFO_SF (vcf_vb); // cleans up allocations - call even if line will be dropped due oSTATUS

        // case: we are reconstructing with --luft and we reconstructed one VCF line
        if (flag.luft) 
            vcf_piz_TOPLEVEL_cb_drop_line_if_bad_oSTATUS_or_no_header (vb);
    }

    else if (dict_id.num == dict_id_fields[VCF_SAMPLES] && flag.samples) 
        vcf_piz_SAMPLES_subset_samples (vcf_vb, rep, reconstructed_len);
}


#undef sample_i
#undef adjustment
#undef snip_i

