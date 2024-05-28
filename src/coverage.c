// ------------------------------------------------------------------
//   coverage.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "coverage.h"
#include "txtfile.h"
#include "contigs.h"
#include "sam_private.h"

void coverage_initialize (VBlockP vb)
{
    vb->coverage.len = CTX(CHROM)->word_list.len + NUM_COVER_TYPES;
    buf_alloc (vb, &vb->coverage, 0, vb->coverage.len, uint64_t, 1, "coverage");
    buf_zero (&vb->coverage); // zero every VB even if already allocated from prev VB

    vb->read_count.len = CTX(CHROM)->word_list.len + NUM_COVER_TYPES;
    buf_alloc (vb, &vb->read_count, 0, vb->read_count.len, uint64_t, 1, "read_count");
    buf_zero (&vb->read_count); // zero every VB even if already allocated from prev VB
    
    vb->unmapped_read_count.len = CTX(CHROM)->word_list.len;
    buf_alloc (vb, &vb->unmapped_read_count, 0, vb->unmapped_read_count.len, uint64_t, 1, "unmapped_read_count");
    buf_zero (&vb->unmapped_read_count); // zero every VB even if already allocated from prev VB
}

// PIZ: main thread: piz_process_recon callback: called in order of VBs (bc can't be called with --test)
void coverage_add_one_vb (VBlockP vb)
{
    buf_add_vector (evb, uint64_t, txt_file->coverage, vb->coverage, "txt_file->coverage");
    buf_add_vector (evb, uint64_t, txt_file->read_count, vb->read_count, "txt_file->read_count");
    buf_add_vector (evb, uint64_t, txt_file->unmapped_read_count, vb->unmapped_read_count, "txt_file->unmapped_read_count");
}

static ConstContigPkgP coverage_get_contigs (rom option_name)
{
    ConstContigPkgP contigs = ((Z_DT(SAM) || Z_DT(BAM)) && sam_hdr_contigs) ? sam_hdr_contigs : ref_get_ctgs (gref);

    ASSINP (!Z_DT(FASTQ) || contigs->contigs.len, "%s: %s for FASTQ only works on files compressed with a reference", z_name, option_name);

    return contigs;
}

typedef struct { WordIndex chrom; STR (chrom_name); uint64_t coverage, read_count; PosType64 LN; } CovDis;

static __attribute__((unused)) DESCENDING_SORTER (sort_by_LN, CovDis, LN)

// output of genocat --coverage, called from piz_one_txt_file
void coverage_show_coverage (void)
{
    ASSERTNOTEMPTY (txt_file->coverage);

    ConstContigPkgP contigs = coverage_get_contigs ("--coverage");

    unsigned chr_width = (flag.show_coverage == COV_ALL ? 26 : 13);

    if (flag.show_coverage != COV_ONE)
        iprintf (is_info_stream_terminal ? "\n--coverage for: %s\n%-*s  %-8s  %-11s  %-10s  %-5s   %s\n" 
                                         : "\n--coverage for: %s\n%*s\t%s\t%s\t%s\t%s\t%s\n", 
                 z_name, chr_width, "Contig", "LN", Z_DT(FASTQ) ? "Reads" : "Alignments", "Bases", "Bases", "Depth");

    txt_file->coverage.len -= NUM_COVER_TYPES; // real contigs only
    ARRAY (uint64_t, coverage, txt_file->coverage);
    ARRAY (uint64_t, read_count, txt_file->read_count);

    uint64_t *coverage_special   = &coverage[coverage_len];
    uint64_t *read_count_special = &read_count[coverage_len];

    // calculate CVR_TOTAL - all bases in the file
    for (uint64_t i=0; i < coverage_len + NUM_COVER_TYPES - 1; i++) {
        coverage_special  [CVR_TOTAL] += coverage[i];
        read_count_special[CVR_TOTAL] += read_count[i];
    }

    // calculate CVR_PRIMARY - all bases in the contig (i.e. excluding unmapped, secondary, supplementary, duplicate, soft-clipped)
    bool has_short_name_contigs = false;
    for (uint64_t i=0; i < coverage_len; i++) {
        coverage_special  [CVR_ALL_CONTIGS] += coverage[i];
        read_count_special[CVR_ALL_CONTIGS] += read_count[i];

        has_short_name_contigs |= ({ STR(chrom_name) ; ctx_get_snip_by_word_index (ZCTX(CHROM), i, chrom_name) ; chrom_name_len <= 5; });
    }

    // case: if we have no short-name contigs, show all contigs 
    if (!has_short_name_contigs) flag.show_coverage = COV_ALL; 

    // sort by contig LN
    buf_alloc (evb, &evb->scratch, 0, coverage_len, CovDis, 0, "scratch");
    for (uint64_t i=0; i < coverage_len; i++) 
        if (coverage[i]) {
            STR (chrom_name);
            ctx_get_snip_by_word_index (ZCTX(CHROM), i, chrom_name);

            BNXT (CovDis, evb->scratch) = (CovDis){
                .chrom          = i,
                .chrom_name     = chrom_name,
                .chrom_name_len = chrom_name_len,
                .coverage       = coverage[i],
                .read_count     = read_count[i],
                .LN             = IS_ASTERISK (chrom_name) ? 0 : contigs_get_LN (contigs, i),
            };
        }

    // no sorting for now - keep in order of header or reference (to do: sort by numeric component chrom number)
    // qsort (STRb(evb->scratch)), sizeof (CovDis), sort_by_LN);

    PosType64 genome_nbases = 0;

    for_buf (CovDis, ctg, evb->scratch) {
        if (flag.show_coverage == COV_ALL || (flag.show_coverage == COV_CHROM && ctg->chrom_name_len <= 5))
            iprintf (is_info_stream_terminal ? "%-*s  %-8s  %-11s  %-10s  %-4.1f%%  %6.2f\n" : "%*s\t%s\t%s\t%s\t%4.1f\t%6.2f\n", 
                     chr_width, ctg->chrom_name, str_bases(ctg->LN).s, str_int_commas (ctg->read_count).s, str_bases(ctg->coverage).s, 
                     100.0 * (float)ctg->coverage / (float)coverage_special[CVR_TOTAL], 
                     ctg->LN ? (float)ctg->coverage / (float)ctg->LN : 0);

        else {
            coverage_special[CVR_OTHER_CONTIGS]   += ctg->coverage; // other non-chromosome contigs
            read_count_special[CVR_OTHER_CONTIGS] += ctg->read_count;
        }

        genome_nbases += ctg->LN;
    }

    char all_coverage[7] = "0";
    if (genome_nbases) // avoid division by 0
        snprintf (all_coverage, sizeof (all_coverage), "%*.2f", (int)sizeof(all_coverage)-1, (float)coverage_special[CVR_ALL_CONTIGS] / (float)genome_nbases);

    if (flag.show_coverage == COV_ONE) 
        iprintf ("%s:\t%s\n", z_name, all_coverage);

    else {
        char *cvr_names[NUM_COVER_TYPES] = { "Other contigs", "All contigs", "Soft clip", "Unmapped", "Secondary", "Supplementary", "Failed filters", "Duplicate", "TOTAL"};

        for (uint64_t i=0; i < NUM_COVER_TYPES; i++) {

            if (coverage_special[i])
                iprintf (is_info_stream_terminal ? "%-*s  %-8s  %-11s  %-10s  %-4.1f%%  %s\n" : "%*s\t%s\t%s\t%s\t%4.1f%%\t%s\n", 
                        chr_width, cvr_names[i], "",
                        (i == CVR_SOFT_CLIP ? "" : str_int_commas (read_count_special[i]).s),
                        str_bases(coverage_special[i]).s,
                        100.0 * (float)coverage_special[i] / (float)coverage_special[CVR_TOTAL],
                        i == CVR_ALL_CONTIGS ? all_coverage : "");
            
            if (i == CVR_OTHER_CONTIGS && is_info_stream_terminal)
                iprint0 ("-----\n");
        }
    }
}

// output of genocat --idxstats - designed to identical to samtools idxstats
void coverage_show_idxstats (void)
{
    ASSERTNOTEMPTY (txt_file->coverage);

    ConstContigPkgP contigs = coverage_get_contigs ("--idxstats");

    txt_file->read_count.len -= NUM_COVER_TYPES; // real contigs only
    ARRAY (uint64_t, read_count, txt_file->read_count);
    ARRAY (uint64_t, unmapped_read_count, txt_file->unmapped_read_count);

    for (uint64_t i=0; i < read_count_len; i++) {
        STR(chrom_name);
        ctx_get_snip_by_word_index (ZCTX(CHROM), i, chrom_name);

        PosType64 len = IS_ASTERISK (chrom_name) ? 0 : contigs_get_LN (contigs, i);

        iprintf ("%.*s\t%"PRIu64"\t%"PRId64"\t%"PRIu64"\n", chrom_name_len, chrom_name, len, read_count[i], unmapped_read_count[i]);
    }

    // FASTQ (but not SAM) unmapped reads
    uint64_t unmapped_unknown = BAFT(uint64_t, txt_file->read_count)[CVR_UNMAPPED]; 
    if (unmapped_unknown)
        iprintf ("*\t0\t0\t%"PRIu64"\n", unmapped_unknown);

    fflush (info_stream); // in case output is redirected
}
