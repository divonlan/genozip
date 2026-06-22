// ------------------------------------------------------------------
//   fastq_parse.c
//   Copyright (C) 2026-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "fastq_private.h"

#define NUM_SPLIT_CTXS 6

// global, used to report exceptions
#define MAX_NUM_LINKERS   2 // if increasing, also define more NONBIO_LINK dids in fastq.h and sam.h
#define MIN_LINKER_LEN    8
#define MAX_LINKER_LEN   33 // longest known linker is Parse v.1/2 - 30 bases

static int n_linkers, umi_len, next_in_master;
bool error_free; // sequencing errors have been corrected  
static struct Linker { // 36 bytes
    uint16_t index_in_read_seq;
    uint8_t seq_len;
    char seq[MAX_LINKER_LEN];
} linker[MAX_NUM_LINKERS];

static Did bc_did[2];

rom fastq_nonbio_type_name (void)
{
    switch (segconf.nonbio_type) {
        case NONBIO_NONE   : return "None";
        case NONBIO_Parse  : return "Parse";
        case NONBIO_10xGen : return "10xGenomics";
        default            : return "Invalid_Nonbio";
    }
}

int fastq_nonbio_get_n_linkers (void)
{
    return n_linkers;
}

// report unrecognized nonbio
StrText fastq_nonbio_get_linker_for_stats (unsigned linker_i)
{
    // note: we use non-ascii commas ﹐ to avoid splitting by comma
    StrText s = {};
    if (!n_linkers || segconf.nonbio_type) return s;

    if (linker_i < n_linkers)
        snprintf (s.s, sizeof(s)-1, "linker=(index:%u len:%u),", 
                  linker[linker_i].index_in_read_seq, linker[linker_i].seq_len);

    else // UMI appears is one of QNAME's barcodes
        for (int bc_i=0; bc_i <= 1; bc_i++) // at most one of them will be set
            if (bc_did[bc_i] != DID_NONE)
                snprintf (s.s, sizeof(s)-1, "umi=(did:%.16s len:%u),", ZCTX(bc_did[bc_i])->tag_name, umi_len);

    return s;
}

static char linker_histo_to_base (const int histo[5], int max_histo)
{
    if      (max_histo == histo[0]) return 'A';
    else if (max_histo == histo[1]) return 'C';
    else if (max_histo == histo[2]) return 'G';
    else                            return 'T';    
}

static void fastq_segconf_find_linkers (VBlockP vb, int seq_len)
{
    int linker_histo[seq_len][4];
    memset (linker_histo, 0, seq_len * 4 * sizeof (int));

    for_line {
        if (DATA_LINE(line_i)->seq.len < seq_len) continue; // ignore reads that are too short

        rom seq  = Btxt(DATA_LINE(line_i)->seq.index);
        rom qual = Btxt(DATA_LINE(line_i)->qual.index);

        for (int i=0; i < seq_len; i++)
            if (qual[i] >= 'A')  // only consider high quality bases
                linker_histo[i][acgt_encode(seq[i])]++;
    }
    
    LineWord l={};
    n_linkers = next_in_master = 0;
    memset (linker, 0, sizeof (linker));

    // find linkers: at least MIN_LINKER_LEN consecutive bases identical in 99% of bases considered
    int max_base_pc[seq_len], max_base_count[seq_len];

    for (int i=0; i < seq_len; i++) {
        max_base_count[i] = MAX_(MAX_(MAX_(linker_histo[i][0], linker_histo[i][1]), linker_histo[i][2]), linker_histo[i][3]);
        int total_base_count = linker_histo[i][0] + linker_histo[i][1] + linker_histo[i][2] + linker_histo[i][3]; // ignoring [4] which is 'N' etc
        max_base_pc[i] = percent (max_base_count[i], total_base_count);
    }

    for (int i=0; i < seq_len; i++) {
        bool is_linker_base = max_base_pc[i] >= 80 || // 80% is a linker
                             (max_base_pc[i] >= 70 && ((i && max_base_pc[i-1] >= 80) || (i < seq_len-1 && max_base_pc[i+1] >= 80))); // we allow 70% on condition that one of the two neighors is 80%

        // case: start of new linker
        if (!l.len && is_linker_base) 
            l = (LineWord){ .index = i, .len = 1 };

        // case: continuation of linker
        else if (l.len && is_linker_base)
            l.len++;

        // case: end of linker - record it if we have room
        else if (l.len && !is_linker_base) {
            if (l.len >= MIN_LINKER_LEN && l.len <= MAX_LINKER_LEN) {
                    
                linker[n_linkers] = (struct Linker){ .index_in_read_seq    = l.index, 
                                                     .seq_len         = l.len  };

                for (int i=0; i < l.len; i++) {
                    linker[n_linkers].seq[i] = linker_histo_to_base (linker_histo[l.index + i], max_base_count[l.index + i]);

                    if (max_base_pc[l.index + i] < 100) 
                        error_free = false; // not all reads agree on the base
                }
                
                if (flag.show_segconf_has)
                    iprintf ("linker[%u]: index_in_read_seq=%u len=%u linker_seq=%.*s\n", n_linkers, 
                             linker[n_linkers].index_in_read_seq, linker[n_linkers].seq_len, STRf(linker[n_linkers].seq));
                
                if (++n_linkers == MAX_NUM_LINKERS) break;
            }
            
            l = (LineWord){}; // reset and continue looking for more linkers
        }
    }    
}

static void fastq_segconf_find_umi_by_qname (VBlockP vb, int seq_len)
{
    bc_did[0] = qf_get_barcode_did(1);
    bc_did[1] = qf_get_barcode_did(2);
    umi_len = 0;

    // test up to 5 lines - that's enough
    for (LineIType line_i=0; (bc_did[0] != DID_NONE || bc_did[1] != DID_NONE) && line_i < MIN_(5, vb->lines.len32); line_i++) {
        rom seq = Btxt(DATA_LINE(line_i)->seq.index);

        for (int bc_i=0; bc_i < 2; bc_i++) 
            if (bc_did[bc_i] != DID_NONE) {
                STR(bc_seq);
                qname_get_barcode (STRtxt (DATA_LINE(line_i)->qname), bc_i + 1, pSTRa(bc_seq));

                if (bc_seq_len > seq_len || str_count_mismatches (bc_seq, seq, bc_seq_len) > 1) // allow one mismatch due to a sequencing error
                    bc_did[bc_i] = DID_NONE;

                else {
                    if (!umi_len) 
                        umi_len = bc_seq_len;
                    else if (umi_len != bc_seq_len)
                        bc_did[bc_i] = DID_NONE; // failed condition that all UMIs must be the same length
                }
            }
    }
}

static SORTER (sort_seq_len)
{
    return ASCENDING_RAW(*(uint32_t *)a, *(uint32_t *)b);
}

static uint32_t get_5_percentile_seq_len (VBlockP vb)
{
    // get array of seq_len of all lines sorted 
    ARRAY_alloc (uint32_t, lens, vb->lines.len, false, vb->scratch, vb, "scratch");
    for (LineIType i=0; i < lens_len; i++)
        lens[i] = DATA_LINE(i)->seq.len;

    qsort (lens, lens_len, sizeof (uint32_t), sort_seq_len);

    // return the %5 shortest (i.e. near the shortest, but not quite)
    uint32_t lens_5_pc = lens[lens_len / 20];

    buf_free (vb->scratch);
    return lens_5_pc;
}

void fastq_segconf_check_if_Parse (VBlockP vb)
{
    // seq_len is the sequence length for >95% of the reads are at least this long
    int seq_len = get_5_percentile_seq_len (vb);
    error_free = true; // initialize optimistically - linkers (and therefore, we assume, barcodes) are free of sequencing errors

    // find linkers (each base is constant across large majority of reads - allow sequencing errors)
    fastq_segconf_find_linkers (vb, seq_len);

    if (!n_linkers) return; // not nonbio

    // find UMI
    fastq_segconf_find_umi_by_qname (vb, seq_len);

    // check for SPLiT-seq / Parse Biosciences: 
    // Parse v3/v4 (58 bases): UMI(10) + BC3(8) + LINKER(12, ATGAGGGGTCAG) + BC2(8) + LINKER(12, TCCAACCACCTC) + BC1(8)
    // SLPiT-seq (86 bases):   UMI(10) + BC3(8) + LINKER(30, GTGGCCGATGTTTCGCATCGGCGTACGACT) + BC2(8) + LINKER(22, ATCCACGTGCTTGAGACTGTGG) + BC1(8). See: https://splitcode.readthedocs.io/en/latest/tutorials_splitseq.html

    if ((!umi_len || umi_len == 10) 
     && n_linkers == 2 
     && linker[0].index_in_read_seq == 18 // 10 bp UMI + 8 bp first barcode 
     && linker[1].index_in_read_seq - (linker[0].index_in_read_seq + linker[0].seq_len) == 8 // 2nd barcode is 8 bp
     && seq_len - (linker[1].index_in_read_seq + linker[1].seq_len) >= 8) { // 3rd barcode is 8 bp - we might have biological data after it 

        segconf.nonbio_type = NONBIO_Parse;

        #define X segconf.split_seq
        X.barcode_index[0] = 10; // fixed UMI length
        X.linker_index[0]  = X.barcode_index[0] + 8; // barcode length is always 8
        X.linker_len[0]    = linker[0].seq_len;
        X.barcode_index[1] = linker[0].index_in_read_seq + linker[0].seq_len;
        X.linker_index[1]  = X.barcode_index[1] + 8;
        X.linker_len[1]    = linker[1].seq_len;
        X.barcode_index[2] = linker[1].index_in_read_seq + linker[1].seq_len;
        X.non_bio_len      = X.barcode_index[2] + 8;

        if      (bc_did[0] != DID_NONE) X.umi_did_i = bc_did[0];
        else if (bc_did[1] != DID_NONE) X.umi_did_i = bc_did[1];

        // prepare container
        const Container_7 con = {
            .nitems_lo = 7,
            .repeats   = 1,
            .items = {
                { .dict_id.num = _FASTQ_NONBIO_UMI,   .separator = { CI0_FIXED_0_PAD, 10              } }, // UMI
                { .dict_id.num = _FASTQ_NONBIO_BC0,   .separator = { CI0_FIXED_0_PAD, 8               } }, // barcode
                { .dict_id.num = _FASTQ_NONBIO_LINK0, .separator = { CI0_FIXED_0_PAD, X.linker_len[0] } }, // linker
                { .dict_id.num = _FASTQ_NONBIO_BC1,   .separator = { CI0_FIXED_0_PAD, 8               } }, // barcode
                { .dict_id.num = _FASTQ_NONBIO_LINK1, .separator = { CI0_FIXED_0_PAD, X.linker_len[1] } }, // linker
                { .dict_id.num = _FASTQ_NONBIO_BC2,   .separator = { CI0_FIXED_0_PAD, 8               } }, // barcode
                { .dict_id.num = _FASTQ_NONBIO_EXCESS                                                   }, // excess (possibly none)
            }
        };

        segconf.nonbio_con_snip_len = sizeof (segconf.nonbio_con_snip);
        container_prepare_snip (&con, 0, 0, qSTRa(segconf.nonbio_con_snip));

        if (X.umi_did_i != DID_NONE) {
            segconf.copy_qname_umi_snip_len = sizeof (segconf.copy_qname_umi_snip);
            seg_prepare_snip_other (SNIP_COPY, ZCTX(X.umi_did_i)->dict_id, 0, 0, segconf.copy_qname_umi_snip);
        }
        #undef X
    }

    // case: unrecognized non-biological: report first 6 QNAMEs (if not already reported)
    else if (n_linkers && !segconf.nonbio_type && z_file) {
        // report first QNAMEs, unless already reported, linkers will be reported in stats_calc_hash_occ
        if (!z_file->n_1st_flav_qnames[QNAME1]) {
            z_file->n_1st_flav_qnames[QNAME1] = MIN_(vb->lines.len, NUM_COLLECTED_WORDS);
            
            for (int line_i=0; line_i < z_file->n_1st_flav_qnames[QNAME1]; line_i++) {
                rom qname = Btxt (DATA_LINE(line_i)->qname.index);
                uint32_t qname_len = DATA_LINE(line_i)->qname.len;
                memcpy (z_file->unk_flav_qnames[QNAME][line_i], qname, MIN_(qname_len, UNK_QNANE_LEN));
            }
        }
    }

    // add linkers as word_index=0 in their respective dictionaries
    if (segconf.nonbio_type) {
        for (int i=0; i < n_linkers; i++)
            ctx_populate_zf_ctx (FASTQ_NONBIO_LINK0 + i, STRa(linker[i].seq), 0);
    }
}

void fastq_parse_seg_initialize (VBlockFASTQP vb)
{
    ctx_consolidate_stats (VB, FASTQ_SQBITMAP, FASTQ_NONBIO, FASTQ_NONBIO_BC0, FASTQ_NONBIO_BC1, FASTQ_NONBIO_BC2,
                           FASTQ_NONBIO_UMI, FASTQ_NONBIO_LINK0, FASTQ_NONBIO_LINK1, FASTQ_NONBIO_EXCESS, DID_EOL); 
    
    // barcodes and the umi are expected to repeat many times in the file
    if (error_free)
        ctx_set_no_stons (VB, FASTQ_NONBIO_BC0, FASTQ_NONBIO_BC1, FASTQ_NONBIO_BC2, DID_EOL);

    // UINT8 required by seg_diff
    ctx_set_ltype (VB, LT_UINT8, FASTQ_NONBIO_LINK0, FASTQ_NONBIO_LINK1, DID_EOL);

    // UMIs are a completely random sequence
    ctx_set_ltype (VB, LT_BLOB, FASTQ_NONBIO_UMI, DID_EOL);
    CTX(FASTQ_NONBIO_UMI)->lcodec = CODEC_RANB;
    CTX(FASTQ_NONBIO_UMI)->lcodec_hard_coded = true;

} 

bool fastq_parse_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ𐤐  dl, STRp(seq), PosType64 gpos_R1, bool is_forward_R1,
                          bool *is_excess_aligned)
{
    #define X segconf.split_seq
    if (seq_len < X.non_bio_len) return false;

    // UMI - copy from QNAME if possible
    if (X.umi_did_i != DID_NONE && ctx_encountered_in_line (VB, X.umi_did_i) &&
        str_issame_(STRlst(X.umi_did_i), seq, 10))

        seg_by_did (VB, STRa(segconf.copy_qname_umi_snip), FASTQ_NONBIO_UMI, 10);
    else 
        seg_add_to_local_fixed (VB, CTX(FASTQ_NONBIO_UMI), seq, 10, LOOKUP_WITH_LENGTH, 10);

    // seg barcodes
    for (int i=0; i < 3; i++)
        seg_by_did (VB, &seq[X.barcode_index[i]], 8, FASTQ_NONBIO_BC0 + i, 8);

    // seg linkers
    for (int i=0; i < 2; i++)
        if (error_free)
            seg_by_ctx (VB, &seq[X.linker_index[i]], X.linker_len[i], CTX(FASTQ_NONBIO_LINK0 + i), X.linker_len[i]);
        else
            seg_diff (VB, CTX(FASTQ_NONBIO_LINK0 + i), vs_ACGTN_DICT, NULL, 
                      &seq[X.linker_index[i]], X.linker_len[i], true, X.linker_len[i]);

    // handle the biological portion - the "excess" based beyond the non-biological
    STRinc (seq, X.non_bio_len);
    MappingType mapping_type;
    *is_excess_aligned = false; // initialize

    // add excess: aligner if possible
    if (seq_len && flag.aligner_available &&
        (mapping_type = aligner_seg_seq (VB, STRa(seq), (vb->R1_vb_i > 0), gpos_R1, is_forward_R1))) {
        
        // lookup from SQBITMAP - reconstructed by fastq_recon_aligned_SEQ
        STRl(snip,48)=48;    
        seg_prepare_snip_other (SNIP_OTHER_LOOKUP, _FASTQ_SQBITMAP, false, 0, snip);

        if (MT(PERFECT) || MT(PERFECT_SPLICED)) snip[snip_len++] = ALIGNED_PERFECT;
        if (MT(SPLICED) || MT(PERFECT_SPLICED)) snip[snip_len++] = ALIGNED_SPLICED; 
        snip_len += str_int (seq_len, &snip[snip_len]);
        
        seg_by_did (VB, STRa(snip), FASTQ_NONBIO_EXCESS, 0);
        CTX(FASTQ_SQBITMAP)->txt_len += seq_len; // account in SQBITMAP - like other aligned reads
        
        *is_excess_aligned = true;
    }

    // add excess: verbatim if not aligned
    else if (seq_len) {
        seg_lookup_with_length_other (VB, CTX(FASTQ_NONBIO_EXCESS), _FASTQ_NONREF, seq_len, 0);
        seg_add_to_local_fixed (VB, CTX(FASTQ_NONREF), STRa(seq), LOOKUP_NONE, seq_len);

        vb->num_verbatim++;
    }

    // add excess: no excess
    else
        seg_by_did (VB, "", 0, FASTQ_NONBIO_EXCESS, 0); // all-the-same if no excess in the file.

    seg_by_did (VB, STRa(segconf.nonbio_con_snip), FASTQ_NONBIO, 0);

    return true;
    #undef X
}

