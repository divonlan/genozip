// ------------------------------------------------------------------
//   vcf_modify.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <math.h>
#include "vcf_private.h"
#include "zip.h"
#include "chrom.h"
#include "container.h"
#include "context.h"

static inline uint32_t vcf_get_n_repeats_genotypes (VBlockVCFP vb)
{
    uint32_t n_alts = N_ALTS;

    // common cases of G. note: we can extend to more cases, but never reduce, to not break back comp
    return vb->ploidy==2 && n_alts==1 ? 3 
         : vb->ploidy==2 && n_alts==2 ? 6 
         : vb->ploidy==2 && n_alts==3 ? 10 
         : vb->ploidy==2 && n_alts==4 ? 15 
         : vb->ploidy==2 && n_alts==5 ? 21 
         : vb->ploidy==2 && n_alts==6 ? 28 // note on ploidy=2 values: 28-21=7, which is one more than 21-15=6 etc. 
         : vb->ploidy==1              ? (n_alts + 1) 
         :                              0;
}

// get number of comma-separated repeats expected in an INFO/FORMAT field
static inline uint32_t vcf_get_n_repeats (VBlockVCFP vb, ContextP ctx)
{
    uint32_t n_alts = N_ALTS;

    return (ctx->header_info.vcf.Number >= 1)                  ? ctx->header_info.vcf.Number // note: VCF_QUAL.header_info is set vcf_zip_initialize
         : (ctx->header_info.vcf.Number == NUMBER_A)           ? n_alts // 0 if n_alts is not set yet
         : (ctx->header_info.vcf.Number == NUMBER_R && n_alts) ? (n_alts + 1)
         : (ctx->header_info.vcf.Number == NUMBER_G && 
            dict_id_is_vcf_format_sf (ctx->dict_id))           ? vcf_get_n_repeats_genotypes (vb) // note: for INFO fields, we don't know the ploidy yet
         : /* fallback */                                        0; 
}

void vcf_segconf_finalize_optimizations (VBlockVCFP vb)
{
    segconf.optimize[VCF_INFO] = !CTX(VCF_INFO)->flags.all_the_same; 
    segconf.optimize[VCF_QUAL] = !CTX(VCF_QUAL)->flags.all_the_same; 
    
    // conditional optimizations    
    segconf.optimize[FORMAT_GP]  = segconf.has[FORMAT_GP] && (segconf.FMT_GP_content == GP_phred || segconf.FMT_GP_content == GP_probabilities);
    segconf.optimize[FORMAT_PP]  = segconf.has[FORMAT_PP] && !segconf.vcf_is_pindel; // note: in Pindel, PP is used for something else
    segconf.optimize[FORMAT_QR]  = segconf.has[FORMAT_QR] && segconf.vcf_is_freebayes; 

    // unconditional optimizations
    segconf.optimize[FORMAT_GL]  = segconf.has[FORMAT_GL];
    segconf.optimize[FORMAT_PL]  = segconf.has[FORMAT_PL];
    segconf.optimize[FORMAT_SPL] = segconf.has[FORMAT_SPL];
    segconf.optimize[FORMAT_PRI] = segconf.has[FORMAT_PRI];
    segconf.optimize[FORMAT_GQ]  = segconf.has[FORMAT_GQ];

    // optimize all floats. note: segconf_calculate guarantees us that all segconf vctx and zctx are the same at this point
    for (Did did_i=VCF_FIRST_OPTIONAL_DID; did_i < vb->num_contexts; did_i++) {
        decl_zctx (did_i);
        if (segconf.has[did_i] && zctx->header_info.vcf.Type == VCF_Float)
            segconf.optimize[did_i] = true;
    }
}

// sorts by alphabetical order, except END that is always first.
static SORTER (INFO_sorter)
{ 
    InfoItemP ina = (InfoItemP)a;
    InfoItemP inb = (InfoItemP)b;
    
    // END comes first (as eg vcf_INFO_SVLEN_prediction depends on POS.last_delta and vcf_seg_melt_ADJLEFT relies on it too)
    if (str_issame_(STRa(ina->name), "END=", 4)) return -1;
    if (str_issame_(STRa(inb->name), "END=", 4)) return 1;
    
    return strncmp (ina->name, inb->name, MIN_(ina->name_len, inb->name_len));
}

// ZIP main thread
void vcf_zip_add_line_numbers_init_vb (VBlockVCFP vb)
{
    // note: normally we update these in zip_update_txt_counters to avoid counting \n of the entire VB upfront, but we need to count in this case
    vb->first_line = txt_file->num_lines + 1; 
    vb->lines.len  = str_count_char (STRb(vb->txt_data), '\n');

    txt_file->num_lines += vb->lines.len;
}

void vcf_add_line_numbers_seg_initialize (VBlockVCFP vb)
{
    // container just for adding a prefix to the delta-encoded line number (the container is all_the_same)
    MiniContainer line_number_container = {
        .repeats   = 1,
        .nitems_lo = 1,
        .items     = { { .dict_id = { _VCF_LINE_NUM } } }
    };

    // create a b250 and dict entry for VCF_LINE_NUM, VCF_ID - these become "all_the_same" so no need to seg them explicitly hereinafter        
    const char prefix[] = { CON_PX_SEP, // has prefix 
                            CON_PX_SEP, // end of (empty) container-wide prefix
                            'L', 'N', '=', CON_PX_SEP };  // NOTE: if changing prefix, update LN_PREFIX_LEN

    container_seg (vb, CTX(VCF_ID), (Container *)&line_number_container, prefix, sizeof prefix, 0); 
    ctx_decrement_count (VB, CTX(VCF_ID), 0);
    
    ctx_set_no_stons (VB, VCF_ID, VCF_LINE_NUM, DID_EOL);
    ctx_consolidate_stats (VB, VCF_ID, VCF_LINE_NUM, DID_EOL);
}

void vcf_seg_line_number_ID (VBlockVCFP vb, STRp(id))
{
    uint32_t snip_len = 3 + str_get_uint_textual_len(vb->first_line + vb->line_i);

    if (vb->line_i == 0)
        seg_integer_as_snip_do (VB, CTX(VCF_LINE_NUM), vb->first_line, snip_len + 1); // +1 for \t
    else
        seg_by_did (VB, (char[]){ SNIP_SELF_DELTA, '1' }, 2, VCF_LINE_NUM, snip_len + 1);
}

static char *vcf_add_line_numbers (VBlockVCFP vb, STRp(id), char *next)
{
    if (flag.add_line_numbers) {
        char *save_next = next;

        next = mempcpy (next, "LN=", 3); 
        next += str_int (vb->first_line + vb->line_i, next);

        CTX(VCF_ID)->txt_shrinkage += (int32_t)id_len - (int32_t)(next - save_next);
    }
    else 
        next = mempcpy (next, id, id_len);

    *next++ = '\t';
    return next;
}

// rounds floating point number to 3 significant digits. Used for SAM too.
// Note that remainder of fraction's digits are truncated rather than rounded.
char *optimize_float_3_sig_dig (VBlockP vb, ContextP ctx, STRp(snip), char *out)
{
    START_TIMER;

    rom save_snip = snip, after = snip + snip_len;
    char *save_out = out;

    if (*snip == '-') 
        *out++ = *snip++; // copy minus sign intact

    // copy int and count int digits
    int sig_digits=0; // significant digits
    if (*snip == '0') // 0.ddd - the 0 is not counted as a signifcant digit 
        *out++ = *snip++;
    else
        while (IS_DIGIT(*snip)) {
            *out++ = *snip++;
            sig_digits++;
        }

    if (snip == after) goto done; // number is an integer

    // expecting a . after the integer part
    if (*snip++ != '.' || snip == after) // no '.' or only '.'
        goto fallback; // bad float format

    // case: we're done with only the integer part - no need for fractions
    if (sig_digits >= 3) {
        // expecting fractional part to be only digits (esp. no e if integer part has more than 1 digit)
        if (!str_is_numeric (snip, after - snip))   
            goto fallback; 
        
        goto done; // outputted integer only 
    }

    *out++ = '.';

    // if no significant integer, all leading 0s in fraction are not significant either - just copy them
    if (!sig_digits)
        while (snip != after && *snip == '0')
            *out++ = *snip++;

    // complete the required 3 significant digits from the fraction
    for (; sig_digits < 3 && snip != after && IS_DIGIT(*snip); sig_digits++)
        *out++ = *snip++;

    // skip the rest of the fraction's digits
    while (snip != after && IS_DIGIT(*snip)) snip++;

    // case: number was a normal fraction
    if (snip == after) {
        // remove trailing 0s from the fraction, and possible the '.' too
        while (out[-1] == '0') out--;
        if    (out[-1] == '.') out--; 

        goto done;
    }

    // verify that the next character (if any) is 'e' or 'E'
    if (*snip != 'e' && *snip != 'E') 
        goto fallback; // invalid format

    // special case: drop e+00 or E+00
    if (after-snip == 4 && snip[1]=='+' && snip[2]=='0' && snip[3]=='0')
        goto done;

    // case: number is a floating point in mantissa/exponent format
    out = mempcpy (out, snip, after - snip); // copy the e/E and all the follows

done:
    ctx->txt_shrinkage += (int32_t)snip_len - (int32_t)(out - save_out);
    COPY_TIMER (optimize_float_3_sig_dig);
    return out;

fallback:
    COPY_TIMER (optimize_float_3_sig_dig);
    return mempcpy (save_out, save_snip, snip_len); 
}

static char *vcf_optimize_float_or_not (VBlockVCFP vb, ContextP ctx, STRp(value), char *next)
{
    if (str_is_1char (value, '.')) goto fallback;
    
    if (value && ctx && ctx->header_info.vcf.Type == VCF_Float && ctx->dict_id.num && segconf.optimize[ctx->did_i]) {
        if (ctx->header_info.vcf.Number == 1)
            next = optimize_float_3_sig_dig (VB, ctx, STRa(value), next);

        else {
            uint32_t max_items = vcf_get_n_repeats (vb, ctx);
            str_split (value, value_len, max_items, ',', item, (max_items>0)); 
            if (!n_items) goto fallback;

            for (uint32_t i=0; i < n_items; i++) {
                next = optimize_float_3_sig_dig (VB, ctx, STRi(item,i), next);
                *next++ = ',';
            }
            next--; // remove final comma
        }

        return next;
    }

fallback:
    return mempcpy (next, value, value_len);
}

static char *vcf_optimize_INFO (VBlockVCFP vb, char *next, STRp(INFO))
{
    START_TIMER;

    vcf_parse_info_subfields (vb, STRa(INFO));

    // case: --optimize-sort: sort the info fields in a consistent orderto reduce word count in the INFO dict. 
    if (segconf.optimize[VCF_INFO] && ii_buf.len32 > 1) 
        qsort (B1ST(InfoItem, ii_buf), ii_buf.len32, sizeof(InfoItem), INFO_sorter);

    // INFO subfield optimizations
    for_buf (InfoItem, ii, ii_buf) {
        next = mempcpy (next, ii->name, ii->name_len); // including '=' if there is one
        next = vcf_optimize_float_or_not (vb, ii->ctx, STRa(ii->value), next);
        *next++ = ';';
    }
    next--; // remove final ';'

    ii_buf.len32 = 0; // reset as seg will re-parse

    COPY_TIMER (vcf_optimize_INFO);
    return next;
} 

// convert an array of probabilities (values∈[0,1]) to an array of integer phred scores capped at 60
static char *vcf_convert_probabilites_to_phred (VBlockVCFP vb, ContextP ctx, STRp(snip), char *next)
{
    START_TIMER;

    char *save_next = next;

    str_split_floats (snip, snip_len, 0, ',', prob, false, 0);
    if (!n_probs) goto error;

    // verify that the field consists of probabilities adding up approximately to 1
    double sum=0;
    for (unsigned i=0; i < n_probs; i++) {
        if (probs[i] < 0 || probs[i] > 1) goto error;
        sum += probs[i];
    }

    if (sum < 0.98 || sum > 1.02) goto error; // allow an epsilon due to rounding errors

    for (unsigned i=0; i < n_probs; i++) {
        int64_t phred = (probs[i] > 1e-60) ? (int64_t)(((-log10 (probs[i])) * 10)+0.5) : 60; // round to the nearest int, capped at 60

        next += str_int (phred, next);
        *next++ = ',';
    }
    next--; // remove final ','

    ctx->txt_shrinkage += (int32_t)snip_len - (int32_t)(next - save_next);
    
    COPY_TIMER (vcf_convert_probabilites_to_phred);
    return next;

error:
    COPY_TIMER (vcf_convert_probabilites_to_phred);
    return mempcpy (save_next, snip, snip_len); // format error - don't optimize
}

// convert an array of likelihoods (non-positive values that are log10(probability)) to an array of integer phred scores capped at 60
static char *vcf_convert_likelihoods_to_phred (VBlockVCFP vb, ContextP ctx, STRp(snip), char *next)
{
    START_TIMER;

    char *save_next = next;

    str_split_floats (snip, snip_len, 0, ',', lhood, false, 0);
    if (!n_lhoods) goto error;

    // verify that the field consists of likelihoods 
    for (unsigned i=0; i < n_lhoods; i++) 
        if (lhoods[i] > 0) goto error;

    for (unsigned i=0; i < n_lhoods; i++) {
        int64_t phred = MIN_(60, (int64_t)(((-lhoods[i]) * 10)+0.5)); // round to the nearest int, capped at 60

        next += str_int (phred, next);
        *next++ = ',';
    }
    next--; // remove final ','

    ctx->txt_shrinkage += (int32_t)snip_len - (int32_t)(next - save_next);
    COPY_TIMER (vcf_convert_likelihoods_to_phred);
    return next;

error:
    COPY_TIMER (vcf_convert_likelihoods_to_phred);
    return mempcpy (save_next, snip, snip_len); // format error - don't optimize
}

// converts an array of phred scores (possibly floats) to integers capped at 60
static char *vcf_phred_optimize (VBlockVCFP vb, ContextP ctx, STRp(snip), char *next)
{
    START_TIMER;

    char *save_next = next;
    if (str_is_1char (snip, '.')) fallback: {
        COPY_TIMER (vcf_phred_optimize);
        return mempcpy (save_next, snip, snip_len); // format error - don't optimize
    }

    uint32_t max_items = vcf_get_n_repeats (vb, ctx);

    if (ctx->header_info.vcf.Type == VCF_Integer) {
        str_split (snip, snip_len, max_items, ',', item, (max_items>0)); 
        if (!n_items) goto fallback; // not an array of ints, or too long of an array

        for (unsigned i=0; i < n_items; i++) {
            int64_t val;
            if (!str_get_int (STRi(item,i), &val)) goto fallback;
            
            next = (val < 60) ? mempcpy (next, items[i], item_lens[i])
                              : mempcpy (next, "60", 2);
            *next++ = ',';
        }
    }

    else { // Float or Unknown
        str_split_floats (snip, snip_len, max_items, ',', item, (max_items>0), '.'); 
        if (!n_items) goto fallback; // not an array of floats, or too long of an array

        for (unsigned i=0; i < n_items; i++) {
            if (isnan (items[i]))  // a '.'
                *next++ = '.';
                
            else {
                int64_t new_phred = MIN_(60, (int64_t)(items[i] + 0.5));
                next += str_int (new_phred, next);
            }
            *next++ = ',';
        }
    }

    next--; // remove final ','

    ctx->txt_shrinkage += (int32_t)snip_len - (int32_t)(next - save_next);
    COPY_TIMER (vcf_phred_optimize);
    return next;
}

static char *vcf_optimize_QUAL (VBlockVCFP vb, char *next, STRp(qual))
{
    START_TIMER;

    decl_ctx(VCF_QUAL);
    
    if (segconf.optimize[VCF_QUAL]) {
        // case: QUAL is expected to be approximately the same as GP[0] - optimize like we optimize GP, so it continues to compress well 
        if (segconf.FMT_GP_content == GP_phred && segconf.vcf_QUAL_method == VCF_QUAL_by_GP)
            next = vcf_phred_optimize (vb, ctx, STRa(qual), next);
        else
            next = optimize_float_3_sig_dig (VB, ctx, STRa(qual), next);
    }
    else
        next = mempcpy (next, qual, qual_len);

    *next++ = '\t';

    COPY_TIMER (vcf_optimize_QUAL);
    return next;
} 

static void vcf_modify_format (STRc (format), bool GL_to_PL, bool GP_to_PP)
{
    if (!GL_to_PL && !GP_to_PP) return;

    str_split (format, format_len, 0, ':', sf, false);

    for (int i=0; i < n_sfs; i++)
        if ((GL_to_PL && str_issame_(STRi(sf,i), "GL", 2)) ||
            (GP_to_PP && str_issame_(STRi(sf,i), "GP", 2)))
            *(char *)&sfs[i][0] = 'P';
}

static char *vcf_optimize_samples (VBlockVCFP vb, STRp(format), rom samples, rom after, char *next, bool *has_13)
{
    START_TIMER;

    str_split (format, format_len, 0, ':', fmt,  false);

    ContextP ctxs[n_fmts];
    for (int i=0; i < n_fmts; i++)
        ctxs[i] = ctx_get_ctx_tag (vb, dict_id_make (STRi(fmt,i), DTYPE_2), fmts[i], fmt_lens[i]); 

    // 0 or more samples (note: we don't use str_split, because samples could be very numerous)
    while (samples[-1] != '\n' && samples[-1] != '\r') {
        rom sample = samples;
        unsigned sample_len;
        int remaining = after - samples;
        char sep;
        samples = (char *)seg_get_next_item (VB, samples, &remaining, GN_SEP, GN_SEP, GN_IGNORE, &sample_len, &sep, has_13, "sample-subfield");

        str_split (sample, sample_len, n_fmts, ':', sf, false);
        if (!n_sfs) { // possibly a split error, don't analyze here, seg will, just keep as is
            next = mempcpy (next, sample, sample_len);
            *next++ = '\t';
        }
        
        else {
            for (uint32_t i=0; i < n_sfs; i++) {
                #define H2(c1,c2)    (c1 | (c2<<8))
                #define H3(c1,c2,c3) (c1 | (c2<<8) | (c3<<16))

                #define CASE2(f,c1,c2)    case H2(c1,c2):    if (!segconf.optimize[FORMAT_##f]) break; 
                #define CASE3(f,c1,c2,c3) case H3(c1,c2,c3): if (!segconf.optimize[FORMAT_##f]) break; 

                if (fmt_lens[i] == 2)
                    switch (H2(fmts[i][0], fmts[i][1])) {
                        case H2('G','T'): vb->ploidy = (sf_lens[i]==1)?1 : (sf_lens[i]==3 && (sfs[i][1] == '/' || sfs[i][1] == '|'))?2 : 0; goto fallback; // set ploidy to 1, 2 or leave 0 if otherwise

                        CASE2(GL,'G','L') // GL=-7.61618,-0.447624,-0.193264 ➔ PL=60,4,2
                            next = vcf_convert_likelihoods_to_phred (vb, ctxs[i], STRi(sf, i), next); break; 

                        CASE2(GP,'G','P') 
                            if (segconf.FMT_GP_content == GP_phred)
                                next = vcf_phred_optimize (vb, ctxs[i], STRi(sf, i), next);  // GP is already Phred, just optimize it
                            else if (segconf.FMT_GP_content == GP_probabilities)
                                next = vcf_convert_probabilites_to_phred (vb, ctxs[i], STRi(sf, i), next); // convert GP (probabilities) to PP (phred values). PP was introduced in VCF v4.3. GP=0.996,0.004,0.000 ➔ PP=0,24,60
                            else
                                goto fallback;
                            break; 

                        CASE2(PL,'P','L')  next = vcf_phred_optimize (vb, ctxs[i], STRi(sf, i), next); break; 
                        CASE2(PP,'P','P')  next = vcf_phred_optimize (vb, ctxs[i], STRi(sf, i), next); break; 
                        CASE2(GQ,'G','Q')  next = vcf_phred_optimize (vb, ctxs[i], STRi(sf, i), next); break; 
                                                
                        default: goto fallback;
                    }

                else if (fmt_lens[i] == 3)
                    switch (H3(fmts[i][0], fmts[i][1], fmts[i][2])) {
                        CASE3(SPL,'S','P','L') next = vcf_phred_optimize (vb, ctxs[i], STRi(sf, i), next); break; 
                        CASE3(PRI,'P','R','I') next = vcf_phred_optimize (vb, ctxs[i], STRi(sf, i), next); break; 
                        default: goto fallback;
                    }

                else fallback: // length other than 2,3
                    next = vcf_optimize_float_or_not (vb, ctxs[i], STRi(sf, i), next);

                *next++ = ':';
            }
            next[-1] = '\t'; // replace final colon with a tab
        }
    }

    next--; // remove final tab

    COPY_TIMER (vcf_optimize_samples);
    return next;
}

rom vcf_zip_modify (VBlockP vb_, rom line_start, uint32_t remaining)
{
    START_TIMER;

    VBlockVCFP vb = (VBlockVCFP)vb_;
    rom next_field=line_start, field_start;
    int32_t len = remaining;
    unsigned field_len=0;
    char separator;
    bool has_13B=false, *has_13=&has_13B;
    rom after = memchr (line_start, '\n', remaining) + 1;

    buf_alloc (vb, &vb->optimized_line, 0, (after - line_start) * 2 + 10 * vcf_num_samples + 1000, char, 0, "optimized_line"); // should be plenty for the types of modifications we have so far.        

    GET_NEXT_ITEM (VCF_CHROM);
    GET_NEXT_ITEM (VCF_POS);
    char *next = mempcpy (B1STc(vb->optimized_line), line_start, next_field - line_start); // initialize to exact copy of fields 1-10

    GET_NEXT_ITEM (VCF_ID);
    
    // add line numbers instead of ID, if needed
    next = vcf_add_line_numbers (vb, field_start, field_len, next);
    
    rom ref = next_field;
    GET_NEXT_ITEM (VCF_REF);
    GET_NEXT_ITEM (VCF_ALT);
    next = mempcpy (next, ref, next_field - ref); // including \t
    N_ALTS = str_count_char (field_start, field_len, ',') + 1;

    GET_NEXT_ITEM (VCF_QUAL);

    // optimize QUAL if needed
    next = vcf_optimize_QUAL (vb, next, field_start, field_len);

    GET_NEXT_ITEM (VCF_FILTER);
    next = mempcpy (next, field_start, next_field - field_start); // including \t

    if (vcf_num_samples) { GET_NEXT_ITEM       (VCF_INFO); }
    else                 { GET_MAYBE_LAST_ITEM (VCF_INFO); } // may or may not have a FORMAT field

    // optimize INFO if needed
    next = vcf_optimize_INFO (vb, next, field_start, field_len); 

    // handle FORMAT and samples
    if (separator != '\n') { // has a FORMAT field
        *next++ = '\t';

        if (vcf_num_samples) { GET_MAYBE_LAST_ITEM (VCF_FORMAT); } // possibly no samples this line
        else                 { GET_LAST_ITEM       (VCF_FORMAT); }

        char *format = next;
        next = mempcpy (next, field_start, field_len);
        uint32_t format_len = next - format;
        
        if (separator != '\n') { // has samples
            *next++ = '\t';

            // case: samples modifications
            if (segconf.optimize[FORMAT_GL] || segconf.optimize[FORMAT_GP] || 
                segconf.optimize[FORMAT_PL] || segconf.optimize[FORMAT_PP] || segconf.optimize[FORMAT_PRI]) 
                next = vcf_optimize_samples (vb, field_start, field_len, next_field, after, next, has_13);
            
            // case: no mod - just copy all samples (excluding \r and \n)
            else {
                *has_13 = (after[-2] == '\r');
                next = mempcpy (next, next_field, after - next_field - 1 - *has_13);
            }
        }

        // now that we have completed optimizing the subfields (which required the original FORMAT)
        // we can chage GP->PP and/or GL->PL in format string
        vcf_modify_format (STRa(format), segconf.optimize[FORMAT_GL], segconf.optimize[FORMAT_GP] && segconf.FMT_GP_content == GP_probabilities);
    }

    *next++ = '\n'; // note: no \r even if there was one

    CTX(VCF_EOL)->txt_shrinkage += *has_13;

    vb->optimized_line.len32 = BNUM(vb->optimized_line, next);
    
    vb->ploidy = N_ALTS = 0; // reset

    COPY_TIMER (vcf_zip_modify);
    return after;
}
