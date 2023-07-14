// ------------------------------------------------------------------
//   vcf_mastermind.c
//   Copyright (C) 2022-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "piz.h"
#include "context.h"
#include "reconstruct.h"
#include "url.h"

sSTRl (hgvsg_con_snip, 64);
sSTRl (copy_chrom_snip, 32);
static Did hgvsg_did_accession, hgvsg_did_hgvs;

sSTRl (copy_gene_snip, 32);

void vcf_mastermind_zip_initialize (void)
{
    // init HGVSG
    DictId hgvsg_dict_id_accession = (DictId)DICT_ID_MAKE1_6("H0GVSG");
    DictId hgvsg_dict_id_hgvs = (DictId)DICT_ID_MAKE1_6("H1GVSG");
    
    hgvsg_did_accession = ctx_add_new_zf_ctx_from_txtheader ("HGVSG_ACCESSION", 15, hgvsg_dict_id_accession, 0)->did_i;
    hgvsg_did_hgvs = ctx_add_new_zf_ctx_from_txtheader ("HGVSG_HGVS", 10, hgvsg_dict_id_hgvs,  0)->did_i;
    
    SmallContainer con = { .repeats   = 1,
                           .nitems_lo = 2,
                           .items     = { { .dict_id = hgvsg_dict_id_accession, .separator[0] = ':' },
                                          { .dict_id = hgvsg_dict_id_hgvs                           } } };

    container_prepare_snip ((ContainerP)&con, 0, 0, qSTRa(hgvsg_con_snip));

    seg_prepare_snip_other (SNIP_COPY, _VCF_CHROM, false, 0, copy_chrom_snip);
    seg_prepare_snip_other (SNIP_COPY, _INFO_GENE, false, 0, copy_gene_snip);
}

void vcf_mastermind_seg_initialize (VBlockVCFP vb)
{
    ctx_consolidate_stats (VB, INFO_HGVSG, hgvsg_did_accession, hgvsg_did_hgvs, DID_EOL);
    vcf_seg_hgvs_consolidate_stats (vb, INFO_HGVSG);
}

// <ID=HGVSG,Number=1,Type=String,Description="HGVS genomic notation for this variant">
// Example: HGVSG=NC_000001.10:g.69511A>C
void vcf_seg_mastermind_HGVSG (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    str_split (value, value_len, 2, ':', item, true);
    
    if (n_items == 2) {
        // CHROM
        if (str_issame_(STRi(item,0), STRa(vb->chrom_name)))
            seg_by_did (VB, STRa(copy_chrom_snip), hgvsg_did_accession, item_lens[0]);
        else
            seg_by_did (VB, STRi(item,0),          hgvsg_did_accession, item_lens[0]);

        // HGVS
        vcf_seg_INFO_HGVS (VB, CTX(hgvsg_did_hgvs), STRi(item,1), 0);
    
        // container
        seg_by_ctx (VB, STRa (hgvsg_con_snip), ctx, 1); // 1 for ':' separator
    }

    else
        seg_by_ctx (VB, STRa (value), ctx, value_len); 

    set_last_txt (ctx->did_i, value);
}

static bool vcf_seg_INFO_MMID3_gene (VBlockP vb, ContextP ctx, STRp(mmid3_gene), uint32_t rep)
{
    STRlast (gene, INFO_GENE);

    if (ctx_encountered (VB, INFO_GENE) && str_issame (gene, mmid3_gene)) 
        seg_by_ctx (VB, STRa(copy_gene_snip), ctx, mmid3_gene_len);
    else
        seg_by_ctx (VB, STRa(mmid3_gene), ctx, mmid3_gene_len);

    return true;
}

// <ID=MMID3,Number=.,Type=String,Description="Mastermind variant identifiers, as gene:key, for MMCNT3">
// Example: MMID3=OR4F5:G136S
void vcf_seg_INFO_MMID3 (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    static MediumContainer con = { 
        .repsep[0] = ',',
        .drop_final_repsep = true,
        .nitems_lo = 2,
        .items     = { { .dict_id.num = DICT_ID_MAKE1_6("M0MID3"), .separator[0] = ':' },
                       { .dict_id.num = DICT_ID_MAKE1_6("M1MID3")                      } } 
    };

    SegCallback callbacks[2] = { vcf_seg_INFO_MMID3_gene, 0 };

    seg_array_of_struct (VB, ctx, con, STRa(value), callbacks, value_len);
}

// <ID=MMURI3,Number=1,Type=String,Description="Mastermind search URI for articles including other DNA-level variants resulting in the same amino acid change">
// Example: MMURI3=https://mastermind.genomenon.com/detail?mutation=NC_000001.10%3Ag.69511A%3EG&ref=cvr
void vcf_seg_INFO_MMURI3 (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    if (ctx_encountered (VB, INFO_HGVSG)) {
        STRlast (hgvsg, INFO_HGVSG);
        char escaped_hgvsg[hgvsg_len * 3 + 1];

        rom mut = memchr (value, '=', value_len);
        if (!mut) goto fallback;
        mut++;

        rom after = memchr (mut, '&', value_len - (mut-value));
        if (!after) goto fallback;

        SAFE_NULT(hgvsg);
        url_esc_non_valid_chars_(hgvsg, escaped_hgvsg, false);
        SAFE_RESTORE;

        if (str_issame_(escaped_hgvsg, strlen(escaped_hgvsg), mut, after - mut)) {
            char snip[2 + mut-value + 1 + value+value_len-after];
            char *next = snip;

            *next++ = SNIP_SPECIAL;
            *next++ = VCF_SPECIAL_MMURI;
            next = mempcpy (next, value, mut-value); // substring before '='
            *next++ = 255; // means: reconstruct escaped HGVSG
            next = mempcpy (next, after, value+value_len - after); // substring following '&'

            seg_by_ctx (VB, snip, next - snip, ctx, value_len); 
        }
        else 
            goto fallback;
    }    
    else fallback: 
        seg_by_ctx (VB, STRa (value), ctx, value_len); 
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MMURI)
{
    STRlast (hgvsg, INFO_HGVSG);
    char escaped_hgvsg[hgvsg_len * 3 + 1];

    SAFE_NULT(hgvsg);
    url_esc_non_valid_chars_(hgvsg, escaped_hgvsg, false);
    SAFE_RESTORE;

    rom placeholder = memchr (snip, 255, snip_len);
    ASSPIZ0 (placeholder, "placeholder not found");

    RECONSTRUCT (snip, placeholder - snip);
    RECONSTRUCT (escaped_hgvsg, strlen (escaped_hgvsg));
    RECONSTRUCT (placeholder + 1, snip + snip_len - (placeholder + 1));

    return NO_NEW_VALUE;
}
