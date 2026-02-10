// ------------------------------------------------------------------
//   sam_bam_assist.c
//   Copyright (C) 2024-2026 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <math.h>
#include "fastq_private.h"
#include "dispatcher.h"
#include "txtfile.h"
#include "qname.h"
#include "progress.h"
#include "mgzip.h"
#include "sections.h"
#include "aligner.h"
#include "refhash.h"
#include "codec.h"
#include "huffman.h"
#include "contigs.h"
#include "libdeflate_1.19/libdeflate.h"

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// ZIP data formats
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

typedef struct {
    uint64_t first_ent, first_aln;
} BamAssIndexer;

// DATA FORMAT IN bamass_alns:
// 4 or 5 bytes: gpos - depending on BamAssEnt.is_long_gpos). 
// (var length)  nico-compressed textual cigar (H and P omitted, X and = collapsed to M, reversed if SAM_FLAG.rev_comp). 

Buffer bamass_ents             = {}; // ZIP: of type BamAssEnt (entry per usable alignment: persist between fastq z_files)
Buffer bamass_heads            = {}; // ZIP: hash table (entry per hash value that is a subset of the qname_hash bits) - head of linked list for each hash value
static Buffer bamass_alns      = {}; // ZIP: variable-length gpos and cigar data described above
static Buffer bamass_ents_bufs = {}; // ZIP: array of BufferP - one per VB: grabbed vb_bamass_ents
static Buffer bamass_alns_bufs = {}; // ZIP: array of BufferP - one per VB: grabbed vb_bamass_alns
static Buffer save_cigar_huff  = {}; // ZIP: save CIGAR huffman between z_files

static double bamass_trims_ratio[3]; // ZIP: bamass compression ratio of SEQ of vb=1,2,3: in vb=1,2,3 we use the 3 bamass_trims options and subsequent VBs use the results to decide which to use
static int generate_vbs_per_linker_vb;

//********************************************************
// PART 1: read SAM/BAM file and populate z_file->bamass_*
//********************************************************

#define BAMASS_INC_STATS(reason) ({ if (flag_show_deep && !segconf_running) vb->deep_stats[reason]++; })

#define UNUSABLE(reason) ({ BAMASS_INC_STATS(reason); *usable = false; goto done; })

#define vb_bamass_ents  vb->codec_bufs[0] // ZIP: BAM-Assist entries
#define vb_bamass_alns  vb->codec_bufs[1] // ZIP: BAM-Assist alignments
#define vb_bamass_cigar CTX(FASTQ_CIGAR)->bamass_cigar // ZIP: bamass cigar of this line (modified binary cigar) 

#define cigar_L_flank B1ST(BamCigarOp, vb_bamass_cigar)
#define cigar_R_flank BLST(BamCigarOp, vb_bamass_cigar)

static void bamass_zip_display_reasons (void)
{
    iprintf ("\n\n%s Alignments breakdown by Bam-Assist usability:\n", dt_name (txt_file->data_type));

    for (int i=0; i < BA_AFTER; i++) 
        if (z_file->deep_stats[i])
            iprintf ("%-14.14s: %"PRIu64" (%.1f%%)\n", (rom[])BAMASS_STATS_NAMES[i], z_file->deep_stats[i], 100.0 * (double)z_file->deep_stats[i] / (double)z_file->deep_stats[0]);

    iprintf ("\nRAM consumption:\nbamass_index: %s\nbamass_ents: %s\nbamass_alns: %s\n",
             str_size (bamass_heads.len * sizeof (uint32_t)).s, 
             str_size (bamass_ents.len * sizeof (BamAssEnt)).s, 
             str_size (bamass_alns.len).s);
}

static void inline bamass_alloc_vb_ents (VBlockFASTQP vb, uint32_t n_cigar_op)
{
    uint32_t est_n_lines  = segconf.line_len ? (vb->txt_data.len / segconf.line_len) : 10000;

    // case: initial allocation
    if (!vb_bamass_ents.len) {
        buf_alloc (vb, &vb_bamass_ents,  0, est_n_lines, BamAssEnt, 1.1, "bamass_ents"); 
        buf_alloc (vb, &vb_bamass_alns,  0, 6 * est_n_lines, uint8_t, 1.1, "bamass_alns"); 
    }

    // add space for one more (usually not needed)
    buf_alloc (vb, &vb_bamass_ents, 1, 0, BamAssEnt, 1.15, NULL);
    buf_alloc (vb, &vb_bamass_alns, nico_max_comp_len (FASTQ_CIGAR, n_cigar_op) + 8, 0, uint8_t, 1.15, NULL); 
}

static void inline bamass_alloc_z_ents (void)
{
    uint64_t est_num_alns = sam_deep_calc_hash_bits(); // note: est_num_alns is a bit of an over-estimate as not all lines are usable

    // pointers into bamass_ents - initialize to NO_NEXT
    bamass_heads.can_be_big = true;
    buf_alloc_exact_255 (evb, bamass_heads, ((uint64_t)1 << num_hash_bits), uint32_t, "global.bamass_heads"); 

    double est_num_vbs = MAX_(1, ceil((double)txt_file->est_seggable_size / (double)segconf.vb_size));
    buf_alloc_zero (evb, &bamass_ents_bufs, 0, est_num_vbs+1, BufferP, 0, "global.bamass_ents_bufs");
    buf_alloc_zero (evb, &bamass_alns_bufs, 0, est_num_vbs+1, BufferP, 0, "global.bamass_alns_bufs");

    if (flag_show_deep) 
        printf ("num_hash_bits=%u est_num_lines*1.15=%"PRIu64" bamass_heads.len=%"PRIu64" bamass_ents.len=%"PRIu64"\n", 
                num_hash_bits, est_num_alns, bamass_heads.len, bamass_ents.len);
}

static rom bamass_get_one_bam_aln (VBlockFASTQP vb, BAMAlignmentFixedP aln, uint32_t remaining_txt_len, 
                                   bool *usable, pSTRp(qname), SamFlags *sam_flags, ConstContigP *hdr_contig,
                                   PosType32 *pos, uint32_t *n_cigar_op, BamCigarOpP *cigar, pSTRp(seq)) // out   
{
    START_TIMER;

    uint32_t aln_len = LTEN32(aln->block_size) + 4;

    // a non-sensical block_size might indicate an false-positive identification of a BAM alignment in bam_unconsumed
    ASSERT (IN_RANGX (aln_len, sizeof (BAMAlignmentFixed), remaining_txt_len), 
            "%s: aln_len=%u is out of range - too small, or goes beyond end of txt data: remaining_txt_len=%u",
            LN_NAME, aln_len, remaining_txt_len);

    int32_t hdr_contig_wi = LTEN32 (aln->ref_id);
    if (!IN_RANGE(hdr_contig_wi, 0, sam_hdr_contigs->contigs.len32) || aln->pos == -1 || aln->n_cigar_op == 0) 
        UNUSABLE(BA_UNMAPPED);

    // consensus sequences don't match any FASTQ sequence
    if (aln->l_read_name > 4 && !memcmp (aln->read_name, "cons", 4))
        UNUSABLE(BA_CONSENSUS);
        
    if (aln->l_seq == 0) UNUSABLE(BA_NO_SEQ);

    *hdr_contig = B(Contig, sam_hdr_contigs->contigs, hdr_contig_wi);

    if ((*hdr_contig)->ref_index == WORD_INDEX_NONE) 
        UNUSABLE(BA_NOT_IN_REF);

    *sam_flags = (SamFlags){ .value = LTEN16(aln->flag) };
    if (sam_flags->secondary || sam_flags->supplementary) UNUSABLE(BA_NOT_PRIMARY);
    if (sam_flags->unmapped) UNUSABLE(BA_UNMAPPED);
    
    // handle segconf
    if (segconf_running) {
        if (aln->next_pos >= 0 || (sam_flags->is_last && !sam_flags->is_first)) 
            segconf.is_paired = true;
    }

    *qname_len = aln->l_read_name - 1;  // -1 bc without \0
    *qname     = aln->read_name;
    *pos       = LTEN32 (aln->pos) + 1; // POS in BAM is 0-based

    *cigar = (BamCigarOpP)((bytes)(aln+1) + aln->l_read_name);
    *n_cigar_op = LTEN16(aln->n_cigar_op);

    uint32_t l_seq = LTEN32(aln->l_seq);
     
    if ((*cigar)[0].n == l_seq && (*cigar)[0].op == BC_S) {
        // to do: look for CG:Z, and if present, replace *cigar, *n_cigar_op with it (also enforce that n_cigar_op<1MB as required nico_compress_one_cigar)
    }

    bytes bam_seq = (bytes)(*cigar + *n_cigar_op);

    ASSERTNOTINUSE (vb->scratch);
    buf_alloc (vb, &vb->scratch, 0, l_seq + 2, char, 0, "scratch");
    bam_seq_to_sam (VB, bam_seq, l_seq, false, false, &vb->scratch, false);

    *seq       = B1STc(vb->scratch);
    *seq_len   = l_seq;
    buf_free (vb->scratch);

done:
    COPY_TIMER (bamass_get_one_bam_aln);
    return (rom)aln + aln_len;
}

static rom bamass_get_one_sam_aln (VBlockFASTQP vb, rom alignment, uint32_t remaining_txt_len, 
                                   bool *usable, pSTRp(qname), SamFlags *sam_flags, ConstContigP *hdr_contig,
                                   PosType32 *pos, uint32_t *n_cigar_op, BamCigarOpP *cigar, pSTRp(seq)) // out   
{
    START_TIMER;

    str_split_by_tab (alignment, remaining_txt_len, 11, NULL, true, true, true); // also advances alignment to next line

    if (IS_ASTERISKi(fld,RNAME) || str_is_1chari(fld,POS,'0') || IS_ASTERISKi(fld,CIGAR))
        UNUSABLE(BA_UNMAPPED);
    
    if (IS_ASTERISKi(fld,SEQ)) UNUSABLE(BA_NO_SEQ);

    // consensus sequences don't match any FASTQ sequence
    if (fld_lens[QNAME] > 4 && !memcmp (flds[QNAME], "cons", 4))
        UNUSABLE(BA_CONSENSUS);

    ASSERT (str_get_uint16 (STRi(fld,FLAG), &sam_flags->value), "%s: invalid FLAG=\"%.*s\"", LN_NAME, STRfi(fld,FLAG));

    if (sam_flags->secondary || sam_flags->supplementary) UNUSABLE(BA_NOT_PRIMARY);
    if (sam_flags->unmapped) UNUSABLE(BA_UNMAPPED);

    // handle segconf
    if (segconf_running) {
        if (!str_is_1chari(fld,PNEXT,'0') || (sam_flags->is_last && !sam_flags->is_first)) 
            segconf.is_paired = true;
    }

    sam_cigar_textual_to_binary (VB, STRi(fld,CIGAR), &vb_bamass_cigar, "bamass_cigar");
    *cigar = cigar_L_flank;
    *n_cigar_op = vb_bamass_cigar.len32;

    if ((*cigar)[0].n == fld_lens[SEQ] && (*cigar)[0].op == BC_S) {
        // to do: look for CG:Z, and if present, replace *cigar, *n_cigar_op with it
    }

    *qname_len = fld_lens[QNAME];
    *qname     = flds[QNAME];

    WordIndex header_ctg_wi = contigs_get_by_name (sam_hdr_contigs, STRi(fld,RNAME));
    if (header_ctg_wi == WORD_INDEX_NONE) UNUSABLE (BA_NOT_IN_REF);

    *hdr_contig = B(Contig, sam_hdr_contigs->contigs, header_ctg_wi);
    
    if ((*hdr_contig)->ref_index == WORD_INDEX_NONE) 
        UNUSABLE(BA_NOT_IN_REF);

    ASSERT (str_get_uint32 (STRi(fld,POS), (uint32_t *)pos), "%s: invalid POS=\"%.*s\"", LN_NAME, STRfi(fld,POS));

    *seq       = flds[SEQ];
    *seq_len   = fld_lens[SEQ];

done:
    COPY_TIMER (bamass_get_one_sam_aln);
    return alignment; // next alignment
}

// compute thread entry point
static void bamass_generate_bamass_ents (VBlockP vb_)
{
    VBlockFASTQP vb = (VBlockFASTQP)vb_;

    START_TIMER; 

    // if the txt file is compressed with BGZF, we uncompress now, in the compute thread
    if (TXT_IS_BGZF) 
        mgzip_uncompress_vb (VB, CODEC_BGZF);    // some of the blocks might already have been decompressed while reading - we decompress the remaining

    rom next  = B1STtxt;
    rom after = BAFTtxt;

    bool attempted_rediscover = false;

    while (next < after) {
        STR0(qname); STR0(seq); // all textual
        bool usable=true/*optimistic*/;
        SamFlags sam_flags;
        PosType32 pos;
        BamCigarOpP cigar = NULL;
        uint32_t n_cigar_op = 0;
        rom this = next;

        vb->scratch.len = 0;
        vb_bamass_cigar.len32 = 0;
        vb_bamass_cigar.name = "bamass_cigar";
        
        BAMASS_INC_STATS(BA_TOTAL);

        ConstContigP hdr_contig = NULL;
        next = IS_BAM_ZIP ? bamass_get_one_bam_aln (vb, (BAMAlignmentFixedP)this, after - this, &usable, pSTRa(qname), &sam_flags, &hdr_contig, &pos, &n_cigar_op, &cigar, pSTRa(seq))
                          : bamass_get_one_sam_aln (vb, this,                     after - this, &usable, pSTRa(qname), &sam_flags, &hdr_contig, &pos, &n_cigar_op, &cigar, pSTRa(seq));

        // case: alignment can participate in bamass
        if (usable) {
            RangeP r = ref_get_range_by_ref_index (VB, hdr_contig->ref_index); // note: bamass_get_one_*_aln verify that ref_index is valid 
            ASSERT (r, "Failed get to get range for ref_index=%d", hdr_contig->ref_index); // not expected, as we already verified ref_index
            
            // replace X, = with M ; N with D ; remove H, P ; reverse if rev_comp ; remove n=0 ops 
            sam_prepare_deep_cigar (VB, cigar, n_cigar_op, sam_flags.rev_comp);

            // gpos of lowest-gpos base consumed in the reference. note: if rev_comp, this correspond to the *last* ref-consumed base of the FASTQ SEQ 
            PosType64 gpos = r->gpos + (pos-1); 
            
            BAMASS_INC_STATS(BA_USABLE);

            if (segconf_running) {
                nico_chew_one_cigar (FASTQ_CIGAR, CIG(CTX(SAM_CIGAR)->deep_cigar));
                
                segconf.total_usable_len += (next - this); 
                vb->lines.count++; 

                DO_ONCE 
                    qname_segconf_discover_flavor (VB, QNAME1, STRa(qname));

                else if (!attempted_rediscover && qname_test_flavor (STRa(qname), QNAME1, segconf.qname_flavor[QNAME1], true) != QTR_SUCCESS) {
                    qname_segconf_rediscover_flavor (VB, QNAME1, STRa(qname));
                    attempted_rediscover = true;
                }

                goto next_iter;
            }

            bamass_alloc_vb_ents (vb, n_cigar_op);

            thool is_last_for_hash = segconf.deep_paired_qname ? sam_flags.is_last : unknown;

            START_TIMER;
            // set BamAssEnt in pre-linking VB format
            BamAssEnt e = { 
                .vb_qname_hash   = deep_qname_hash (VB, QNAME1, STRa(qname), is_last_for_hash, NULL),
                .seq_hash        = deep_seq_hash (VB, STRa(seq), sam_flags.rev_comp),
                .aln_lo          = vb_bamass_alns.len32, // note: within a VB, definitely aln_lens < 4GB
                .vb_is_forward   = !sam_flags.rev_comp,
                .vb_is_long_gpos = (gpos > 0xffffffffULL),
            };
            COPY_TIMER (bamass_generate_bamass_ents_hash);

            BNXT (BamAssEnt, vb_bamass_ents) = e;

            uint8_t *next = BAFT8 (vb_bamass_alns);

            // gpos - 4 or 5 bytes 
            if (e.vb_is_long_gpos) *(uint40_t *)next = U64to40(gpos);
            else                   *(unaligned_uint32_t *)next = gpos;  
            next += e.vb_is_long_gpos ? sizeof (uint40_t) : sizeof (uint32_t);

            next += nico_compress_cigar (VB, FASTQ_CIGAR, CIG(CTX(SAM_CIGAR)->deep_cigar), next, vb_bamass_alns.size - BNUM(vb_bamass_alns, next));

            vb_bamass_alns.len32 = BNUM (vb_bamass_alns, next); 
            
            vb->lines.count++; // count only usable lines

            if (flag.show_deep == SHOW_DEEP_ALL || 
                (flag.show_deep == SHOW_DEEP_ONE_HASH && flag.debug_deep_hash.qname == e.vb_qname_hash && flag.debug_deep_hash.seq == e.seq_hash)) {
                
                rom rname = Bc(ZCTX(CHROM)->dict, hdr_contig->char_index); // note: might be different that reference contig name
                uint32_t rname_len = hdr_contig->snip_len;

                iprintf ("\n%s_inspect/%u/%u Recording SAM alignment hash={%016"PRIx64",%08x}\nQNAME=\"%.*s\" (is_last=%s)\nSEQ=\"%.*s\"\nRNAME=\"%.*s\" POS=%u BAMASS_CIGAR=\"%s\"\n",
                         dt_name(txt_file->data_type), vb->vblock_i, vb->line_i, e.vb_qname_hash, e.seq_hash, STRf(qname), YN(is_last_for_hash), STRf(seq),
                         STRf(rname), pos, dis_binary_cigar (VB, cigar, n_cigar_op, &vb->codec_bufs[4]).s);
            }
        }

    next_iter:
        vb->line_i++;
        buf_free (vb_bamass_cigar);
    }

    buf_free (vb->scratch); // stores textual_seq in bam

    // tell dispatcher this thread is done and can be joined.
    vb_set_is_processed (VB); 
    COPY_TIMER (bamass_generate_bamass_ents);
}

static void bamass_segconf (void)
{
    START_TIMER;

    SAVE_FLAGS;
    flag.show_reference = false;
    flag.show_vblocks = NULL;

    segconf.vb_size = 4 MB; // this is quite large, because we are generating cigar stats

    VBlockP vb = vb_initialize_nonpool_vb (VB_ID_SEGCONF, txt_file->data_type, "segconf");

    txtfile_read_vblock (vb);

    nico_start_chewing (FASTQ_CIGAR);

    segconf.running = true;
    bamass_generate_bamass_ents (vb);

    segconf.line_len = vb->lines.count ? (segconf.total_usable_len / vb->lines.count) : 0; 

    segconf.deep_paired_qname = segconf.is_paired; // note: true means that the BAM may have 2 primary alignments with the same canonical read name. we need to set regardless of the number of FASTQs we are compressing, so that each primary alignment has a good chance of a unique deep_hash(qname) 

    segconf.sam_qname_line0 = segconf.qname_line0[QNAME1]; // survives to fastq in deep/bamass

    segconf.sam_is_unmapped = false; // always false, as we discovered usable (i.e. mapped) aligments

    segconf_set_use_insertion_ctxs();

    // note: we (sometimes) replace the I with an S and not the other way around, bc in fastq_bamass_retrieve_ent we might extend the S
    nico_produce_compressor (FASTQ_CIGAR, true);
    
    // return segconf data 
    buf_insert (evb, txt_file->unconsumed_txt, char, 0, vb->txt_data.data, vb->txt_data.len32, "txt_file->unconsumed_txt");

    if (TXT_IS_BGZF) 
        mgzip_return_segconf_blocks (vb); // return BGZF used by the segconf VB to the unconsumed BGZF blocks

    vb_destroy_vb (&vb);
    txt_file->num_lines = 0;
    segconf.running = false;

    RESTORE_FLAGS;

    COPY_TIMER_EVB (bamass_segconf);
}

static void bamass_read_one_vb (VBlockP vb)
{
    START_TIMER;

    txtfile_read_vblock (vb);

    if (Ltxt) {
        vb->preprocessing = true;
        vb->dispatch = READY_TO_COMPUTE;
    }

    else 
        vb->dispatch = DATA_EXHAUSTED;

    COPY_TIMER_EVB (bamass_read_one_vb);
}

static void bamass_append_z_ents (VBlockP vb_)
{
    VBlockFASTQP vb = (VBlockFASTQP)vb_; 
    START_TIMER;

    // we can only handle 4G entries, because hash entry and "next" fields are 32 bit.
    int64_t excess = bamass_ents.len + vb_bamass_ents.len - MAX_BAMASS_ENTS;  
    if (excess > 0) {
        WARN_ONCE ("The number of usable alignments in %s exceeds %u - the maximum supported for Bam-Assist compression. This can affect the compression ratio, but is otherwise unharmful. %s", 
                    txt_name, MAX_BAMASS_ENTS, report_support());

        vb_bamass_ents.len -= excess;
        vb->deep_stats[BA_USABLE] -= excess;
        vb->deep_stats[BA_OVERFLOW] += excess;

        if (!vb_bamass_ents.len) goto done;
    }

    // if this VB alns go beyond maximum, discard entire VB (we don't know how to take partial alns)
    if (bamass_alns.len + vb_bamass_alns.len > MAX_BAMASS_ALNS_LEN) {
        vb->deep_stats[BA_OVERFLOW] += vb->deep_stats[BA_USABLE];
        vb->deep_stats[BA_USABLE]   = 0;
        goto done;
    }
    
    // grab the vb ent/alns buffers and store them for now. bamass_link_entries will copy them to the global buffer.
    // note: bamass_ents_bufs is an array of BufferP, not Buffer! Reason is that we don't yet know
    // how many VBs we will have, and hence might realloc the bamass_ents_bufs - but it is not possible to 
    // realloc an array of Buffers as it messes up the buf_list. This way each vb ents/alns Buffer is 
    // fixed in memory throughout its life.
    vb_bamass_ents.param = bamass_ents.len; // start entry in global array
    bamass_ents.len += vb_bamass_ents.len;  // just count length for now, we will allocate later
    buf_alloc_zero (NULL, &bamass_ents_bufs, 0, vb->vblock_i+1, BufferP, 0, NULL);
    MAXIMIZE (bamass_ents_bufs.len32, vb->vblock_i+1);
    *B(BufferP, bamass_ents_bufs, vb->vblock_i) = CALLOC (sizeof (Buffer)); 
    buf_grab (evb, **B(BufferP, bamass_ents_bufs, vb->vblock_i), NULL, vb_bamass_ents);

    vb_bamass_alns.param = bamass_alns.len; 
    bamass_alns.len += vb_bamass_alns.len;
    buf_alloc_zero (NULL, &bamass_alns_bufs, 0, vb->vblock_i+1, BufferP, 0, NULL);
    MAXIMIZE (bamass_alns_bufs.len32, vb->vblock_i+1);
    *B(BufferP, bamass_alns_bufs, vb->vblock_i) = CALLOC (sizeof (Buffer)); 
    buf_grab (evb, **B(BufferP, bamass_alns_bufs, vb->vblock_i), NULL, vb_bamass_alns);

    if (flag.show_deep) {
        #define ONCE_EVERY 10000000 // print once every 10M ents or so
        static int64_t prev_print_ents_len = -ONCE_EVERY;
        
        if (prev_print_ents_len / ONCE_EVERY != bamass_ents.len / ONCE_EVERY) { 
            if (prev_print_ents_len < 0) iprint_newline();

            iprintf ("bamass entries: %s\tbamass_heads: %s\tbamass_ents: %s\tbamass_alns: %s\ttotal: %s\n", 
                    str_int_commas (bamass_ents.len).s, 
                    str_size (bamass_heads.len * sizeof (uint32_t)).s, // doesn't change after initialization
                    str_size (bamass_ents.len * sizeof (BamAssEnt)).s, 
                    str_size (bamass_alns.len).s,
                    str_size (bamass_heads.len * sizeof (uint32_t) + bamass_ents.len * sizeof (BamAssEnt) + bamass_alns.len).s);

            prev_print_ents_len = bamass_ents.len;
        }
    }

done:
    for (BamassStats i=0; i < BA_AFTER; i++)
        z_file->deep_stats[i] += vb->deep_stats[i];

    COPY_TIMER_EVB (bamass_append_z_ents);
}

static void bamass_prepare_to_link_entries (VBlockP vb)
{
    if (vb->vblock_i < bamass_ents_bufs.len32) {
        vb->preprocessing = true;
        vb->dispatch = READY_TO_COMPUTE;
    }
    
    else
        vb->dispatch = DATA_EXHAUSTED;
}

// compute thread entry: index and update ents generated by one vb 
static void bamass_link_entries (VBIType vb_i)
{
    // initialize ents, alns to vb arrays
    ConstBufferP vb_ents_buf = *B(BufferP, bamass_ents_bufs, vb_i);
    const BamAssEnt *vb_start_ents = B1ST(BamAssEnt, *vb_ents_buf);
    BamAssEntP start_ents = B(BamAssEnt, bamass_ents, vb_ents_buf->param);
    memcpy (start_ents, vb_start_ents, vb_ents_buf->len * sizeof (BamAssEnt));

    ConstBufferP vb_alns_buf = *B(BufferP, bamass_alns_bufs, vb_i);
    rom vb_start_alns = B1STc(*vb_alns_buf);
    char *start_alns = Bc(bamass_alns, vb_alns_buf->param);
    memcpy (start_alns, vb_start_alns, vb_alns_buf->len);

    uint64_t mask = bitmask64 (num_hash_bits); 
    uint32_t *heads = B1ST32 (bamass_heads);
    BamAssEntP first_ent = B1ST(BamAssEnt, bamass_ents);

    // add ents to the linked lists: each thread has exclusive access to the ents in its fragment, 
    // and access to bamass_heads which contains the heads of each linked list is managed by atomic ops 
    for (BamAssEntP ent=start_ents, after_ent = start_ents + vb_ents_buf->len; ent < after_ent; ent++) {
        
        // make our ent the new head of the linked list and set old_head to the previous head
        uint32_t hash     = ent->vb_qname_hash & mask; 
        uint32_t new_head = ent - first_ent;
        uint32_t old_head = heads[hash]; // modified by each failed called to __atomic_compare_exchange_n until we successfully assign the new head
        while (!__atomic_compare_exchange_n (&heads[hash], &old_head, new_head, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {};

        // we overwrite the 10 LSb of qname_hash: 2 bits for flags, and 8 for z_aln_hi ; the 54 MSb of qname_hash remain intact 
        ent->z_is_forward   = ent->vb_is_forward;
        ent->z_is_long_gpos = ent->vb_is_long_gpos;

        // update our ent (our thread has exclusive access to it)
        ent->next = old_head; // previous head (or NO_NEXT) is now the 2nd element in the linked list
        uint64_t z_aln = vb_alns_buf->param + ent->aln_lo;
        ASSERTINRANGX (z_aln, 0, MAX_BAMASS_ALNS_LEN);

        ent->aln_lo   = (uint32_t)z_aln;
        ent->z_aln_hi = z_aln >> 32; // 8bits that also overwrite LSbs of qname_hash
    }
}

static void bamass_link_entries_several_vbs (VBlockP vb)
{
    START_TIMER;

    VBIType first_vb_i = generate_vbs_per_linker_vb * (vb->vblock_i-1) + 1;
    VBIType after_vb_i = MIN_(first_vb_i + generate_vbs_per_linker_vb, bamass_ents_bufs.len32);

    for (VBIType vb_i=first_vb_i; vb_i < after_vb_i; vb_i++)
        bamass_link_entries (vb_i);

    vb_set_is_processed (vb); 
    COPY_TIMER (bamass_link_entries); 
}

static void bamass_after_link_entries (VBlockP vb)
{
    START_TIMER;

    VBIType first_vb_i = generate_vbs_per_linker_vb * (vb->vblock_i-1) + 1;
    VBIType after_vb_i = MIN_(first_vb_i + generate_vbs_per_linker_vb, bamass_ents_bufs.len32);

    // free grabbed buffers
    for (VBIType vb_i=first_vb_i; vb_i < after_vb_i; vb_i++) {
        buf_destroy (**B(BufferP, bamass_ents_bufs, vb_i));
        FREE (*B(BufferP, bamass_ents_bufs, vb_i));

        buf_destroy (**B(BufferP, bamass_alns_bufs, vb_i));
        FREE (*B(BufferP, bamass_alns_bufs, vb_i));
    }

    COPY_TIMER (bamass_after_link_entries);
}

// ZIP: runs during FASTQ / FASTA zip_initialize if --bamass specified
void fastq_bamass_populate (void)
{
    START_TIMER;

    // 2nd+ file or pair - bamass is already populated, just restore the cigar huffman
    if (file_i) {
        buf_move (evb, ZCTX(FASTQ_CIGAR)->huffman, NULL, save_cigar_huff);
        return;
    }

    ASSINP0 (*p_genome_nbases <= MAX_POS40, "Reference file too big: --bamass can handle reference files up to 1 trillion basepairs");

    const rom task_name = "bam_assist";

    FileP save_txt_file = txt_file;
    txt_file = file_open_txt_read (flag.bam_assist);
    ASSERTNOTNULL (txt_file);

    ASSERT (TXT_DT(BAM) || TXT_DT(SAM), "--bamass requires an argument that is a SAM, BAM or CRAM file, but \"%s\" isn't one", flag.bam_assist);

    // verify reference matches SAM header and populates contigs
    txtfile_read_header (true); 
    sam_header_inspect (evb, &evb->txt_data, (struct FlagsTxtHeader){}); 

    ASSINP (sam_hdr_contigs, "%s is missing a SAM header, or SAM header is lacking @SQ records", txt_name);

    buf_free (evb->txt_data);   // discard the header
    buf_free (evb->gz_blocks);

    flag.preprocessing = PREPROC_RUNNING;

    bamass_segconf();

    bamass_alloc_z_ents();

    segconf_set_vb_size (NULL, 0); // fastq will reset this to its own vb_size later

    // step one: create ents
    { START_TIMER;
    dispatcher_fan_out_task (task_name, 
                             save_txt_file->basename, // not txt_file-> bc it will be closed in a sec, while the progress component will continue to the main zip fan_out 
                             0, "Inspecting BAM...",   
                             false, false, flag.xthreads, 0, 5000, true,
                             bamass_read_one_vb, 
                             bamass_generate_bamass_ents, 
                             bamass_append_z_ents);
    COPY_TIMER_EVB (bamass_generate); }

    if (flag_show_deep)
        bamass_zip_display_reasons();

    // note: we now know the exact lengths, so we can allocate
    bamass_ents.can_be_big = true;
    buf_alloc (evb, &bamass_ents, 0, bamass_ents.len, BamAssEnt, 0, "global.bamass_ents");

    bamass_alns.can_be_big = true;
    buf_alloc (evb, &bamass_alns, 0, bamass_alns.len, uint8_t, 0, "global.bamass_alns"); 

    generate_vbs_per_linker_vb = ceil ((double)bamass_ents_bufs.len32 / (double)global_max_threads);

    // step two: create linked lists
    { START_TIMER;
    dispatcher_fan_out_task (task_name, 
                             save_txt_file->basename, // not txt_file-> bc it will be closed in a sec, while the progress component will continue to the main zip fan_out 
                             0, "Inspecting BAM...",   
                             false, false, flag.xthreads, 0, 5000, true,
                             bamass_prepare_to_link_entries, 
                             bamass_link_entries_several_vbs, 
                             bamass_after_link_entries);
    COPY_TIMER_EVB (bamass_link); }

    progress_erase(); // erase "Inspecting BAM..."
    
    file_close (&txt_file);
    segconf.vb_size = 0; // so fastq can recalculate to its liking
    flag.preprocessing = NO_PREPROC;
    
    txt_file = save_txt_file;

    // return unused memory to libc
    buf_destroy (bamass_ents_bufs);
    buf_destroy (bamass_alns_bufs);
    buf_low_level_release_memory_back_to_kernel();
    
    COPY_TIMER_EVB (fastq_bamass_populate);
}

//*******************************************************
// PART 2: seg FASTQ: find alignment in bamass_*
//*******************************************************

void fastq_bamass_seg_initialize (VBlockFASTQP vb)
{
    sam_seg_SEQ_initialize (VB);

    // to decide how to handle trims and S in CIGAR (3 options: make them all M, all S, or M unless flank is S)    
    // we test in vb=1,2,3 and subsequent VBs use the best of the 3 if the results are already known
    // or TRIM_IS_M if not known yet.
    __atomic_thread_fence (__ATOMIC_ACQUIRE); // make sure all bamass_trims_ratio entries are visible in this thread if already set in another thread

    // take a snapshop in time as these values can be changed by other threads
    double r[3];
    memcpy (r, bamass_trims_ratio, sizeof(r)); 

    CTX(FASTQ_CIGAR)->bamass_trims = 
        IN_RANGX(vb->vblock_i, 1, 3)                                    ? vb->vblock_i-1  // vb=1 tests M, vb=2 tests S, vb=3 tests MS
      : (!r[TRIM_IS_MS] && r[TRIM_IS_M] >= r[TRIM_IS_S])                ? TRIM_IS_M       // MS is not known : judge by M and S (note that if both are 0, default is M) 
      : (!r[TRIM_IS_MS] && r[TRIM_IS_M] < r[TRIM_IS_S])                 ? TRIM_IS_S
      : (r[TRIM_IS_MS] && (!r[TRIM_IS_M] || !r[TRIM_IS_S]))             ? TRIM_IS_M       // MS is known but either M or S are not: default to M
      : (r[TRIM_IS_M] >= r[TRIM_IS_S] && r[TRIM_IS_M] >= r[TRIM_IS_MS]) ? TRIM_IS_M       // M is known to be best
      : (r[TRIM_IS_S] >= r[TRIM_IS_M] && r[TRIM_IS_S] >= r[TRIM_IS_MS]) ? TRIM_IS_S       // S is known to be best
      :                                                                   TRIM_IS_MS;     // MS is known to be best

    segconf.bamass_trims = CTX(FASTQ_CIGAR)->bamass_trims; // note: unstable until all three are tested
}

void fastq_bamass_seg_finalize (VBlockFASTQP vb)
{
    bits_truncate ((BitsP)&CTX(SAM_SQBITMAP)->local, CTX(SAM_SQBITMAP)->next_local); // remove unused bits due to perfect mathcing of reference
}

void fastq_bamass_zip_comp_cb (VBlockFASTQP vb, ContextP ctx, SectionType st, uint32_t comp_len)
{
    if (vb->vblock_i <= 3 && (st == SEC_B250 || st == SEC_LOCAL) && ctx->st_did_i == FASTQ_SQBITMAP) 
        vb->seq_comp_len += comp_len; 
}

void fastq_bamass_zip_after_compress (VBlockFASTQP vb)
{
    if (vb->vblock_i <= 3) {
        int tested = vb->vblock_i-1;
        bamass_trims_ratio[tested] = (double)CTX(FASTQ_SQBITMAP)->txt_len / (double)vb->seq_comp_len;
    
        if (flag.show_deep)
            iprintf ("TRIM_IS_%-2s SEQ compression ratio: %2.2f\n", 
                     (rom[])TRIM_IS_NAMES[tested], bamass_trims_ratio[tested]);
    }
}

void fastq_bamass_consider_stopping_aligner (VBlockFASTQP vb)
{
    // in case the vast majority of the non-bamass reads are contamination (cannot be aligned to the reference)
    // stop using the aligner from now on. note: aligner is very slow on reads that fail alignment. 
    double threashold = flag.fast?0.80 : flag.best?0.99 : 0.95;
    double fail_rate;
    uint32_t num_attempted_alignment = vb->num_verbatim + vb->num_aligned; 
    if (num_attempted_alignment > 100 && // at least 100 reads in this VB had alignment attempted
        (fail_rate = (double)vb->num_verbatim / (double)num_attempted_alignment) > threashold) { // too many reads failed alignment
        
        store_release (flag.aligner_available, false); 
        if (flag.show_deep)
            iprintf ("\n%s: FYI: Stopped using aligner from now on, because its fail-rate for non-bamass reads is %.2f%%\n", VB_NAME, fail_rate);
    }
}

static void cigar_remove_flanking_0_ops (BufferP cigar)
{
    // remove flanking n=0 ops (possibly added M or S)
    if (B1ST(BamCigarOp, *cigar)->n == 0)
        buf_remove (*cigar, BamCigarOp, 0, 1);

    if (cigar->len32 && BLST(BamCigarOp, *cigar)->n == 0)
        cigar->len32--;
}

static inline StrTextMegaLong display_bamass_cigar (VBlockP vb, thool is_fwd)
{
    ARRAY (BamCigarOp, cigar, vb_bamass_cigar);
    
    // NOTE: destructive is is_fwd != unknown !!! caller should take care
    if (is_fwd != unknown) { // recover SAM cigar
        cigar_remove_flanking_0_ops (&vb_bamass_cigar);

        if (!is_fwd) {
            for (uint32_t op_i=0; op_i < vb_bamass_cigar.len32 / 2; op_i++) 
                SWAP(cigar[op_i], cigar[vb_bamass_cigar.len32 - op_i - 1]);
        }
    }
    
    StrTextMegaLong s = dis_binary_cigar (vb, cigar_L_flank, vb_bamass_cigar.len32, &vb->codec_bufs[6]);

    return s;
}

// sets dl->sam_seq_len, vb->ref_consumed, vb->gpos and bamass_cigar
void fastq_bamass_retrieve_ent (VBlockP vb, // note: doesn't require FASTQ VB, can ran on evb
                                const BamAssEnt *e, 
                                bool get_cigar, // if true, cigar is written to vb_bamass_cigar
                                bool *is_fwd, PosType64 *gpos, uint32_t *seq_len, uint32_t *ref_consumed, uint32_t *ref_and_seq_consumed, uint32_t *insertions) // optional out 
{
    START_TIMER;
    
    uint8_t *next = B8(bamass_alns, (uint64_t)e->aln_lo | ((uint64_t)e->z_aln_hi << 32));
    
    PosType64 my_gpos = e->z_is_long_gpos ? U40to64(*(uint40_t *)next) : *(unaligned_uint32_t *)next;
    next += e->z_is_long_gpos ? sizeof (uint40_t) : sizeof (uint32_t);

    next += nico_uncompress_cigar (VB, FASTQ_CIGAR, next, &vb_bamass_cigar, "bamass_cigar");

    uint32_t sam_seq_len=0, sam_ref_consumed=0, sam_ref_and_seq_consumed=0, sam_insertions=0; 

    if (seq_len || ref_consumed) 
        for_cigar (vb_bamass_cigar) {
            case BC_M : sam_seq_len += op->n; sam_ref_consumed += op->n; sam_ref_and_seq_consumed += op->n; break;
            case BC_I : sam_insertions   += op->n; // fallthrough
            case BC_S : sam_seq_len      += op->n; break;
            case BC_D : sam_ref_consumed += op->n; break;
            // note: bamass cigar does not have N,H,X,=,P
        }

    uint32_t my_seq_len              = sam_seq_len;
    uint32_t my_ref_consumed         = sam_ref_consumed;
    uint32_t my_ref_and_seq_consumed = sam_ref_and_seq_consumed;

    // case: add 0S or 0M to each side that doesn't already have an S/M - to be used later to account for trimming
    if (get_cigar) {
        // case TRIM_IS_S: add a flanking S if one doesn't already exist
        switch (CTX(FASTQ_CIGAR)->bamass_trims) {
            case TRIM_IS_S: // trims are treated as S
                if (cigar_L_flank->op != BC_S)
                    buf_insert (vb, vb_bamass_cigar, BamCigarOp, 0, &(BamCigarOp){ .op = BC_S }, 1, "bamass_cigar");

                if (cigar_R_flank->op != BC_S)
                    buf_append_one (vb_bamass_cigar, (BamCigarOp){ .op = BC_S });
                break;

            case TRIM_IS_MS: // trims extend existing flanking S or M, and if flank is not S or M, it is M 
                if (cigar_L_flank->op != BC_M && cigar_L_flank->op != BC_S)
                    buf_insert (vb, vb_bamass_cigar, BamCigarOp, 0, &(BamCigarOp){ .op = BC_M }, 1, "bamass_cigar");

                if (cigar_R_flank->op != BC_M && cigar_R_flank->op != BC_S)
                    buf_append_one (vb_bamass_cigar, (BamCigarOp){ .op = BC_M });
                break;

            case TRIM_IS_M: // S flanks are turned into M
                if (cigar_L_flank->op == BC_S) {
                    cigar_L_flank->op = BC_M;
                    if (e->z_is_forward) my_gpos -= cigar_L_flank->n; 
                    my_ref_consumed += cigar_L_flank->n;
                    my_ref_and_seq_consumed += cigar_L_flank->n;

                    // case: two first ops are M - merge them
                    if (vb_bamass_cigar.len32 > 1 && (cigar_L_flank+1)->op == BC_M) { 
                        (cigar_L_flank+1)->n += cigar_L_flank->n;
                        buf_remove (vb_bamass_cigar, BamCigarOp, 0, 1);
                    }        
                }

                if (cigar_L_flank->op != BC_M) // test if originally M or perhaps just converted S->M
                    buf_insert (vb, vb_bamass_cigar, BamCigarOp, 0, &(BamCigarOp){ .op = BC_M }, 1, "bamass_cigar");

                if (cigar_R_flank->op == BC_S) {
                    cigar_R_flank->op = BC_M;
                    if (!e->z_is_forward) my_gpos -= cigar_R_flank->n; 
                    my_ref_consumed += cigar_R_flank->n;
                    my_ref_and_seq_consumed += cigar_R_flank->n;

                    // case: two last ops are M - merge them
                    if (vb_bamass_cigar.len32 > 1 && (cigar_R_flank-1)->op == BC_M) {
                        (cigar_R_flank-1)->n += cigar_R_flank->n;
                        vb_bamass_cigar.len32--;
                    }
                }

                if (cigar_R_flank->op != BC_M)
                    buf_append_one (vb_bamass_cigar, (BamCigarOp){ .op = BC_M });
                
                break;
        }
    }

    else // not get_cigar
        buf_free (vb_bamass_cigar);

    if (is_fwd)               *is_fwd               = e->z_is_forward;
    if (gpos)                 *gpos                 = my_gpos;
    if (seq_len)              *seq_len              = my_seq_len;
    if (ref_consumed)         *ref_consumed         = my_ref_consumed;
    if (ref_and_seq_consumed) *ref_and_seq_consumed = my_ref_and_seq_consumed;
    if (insertions)           *insertions           = sam_insertions;

    COPY_TIMER (fastq_bamass_retrieve_ent);
}

// displays bamass_ents AFTER linking
StrTextLong bamass_dis_ent (VBlockP vb, const BamAssEnt *e, uint64_t qname_hash)
{
    StrTextLong s;
    uint32_t s_len = 0;

    // save bamass_cigar
    ASSERTNOTINUSE (vb->codec_bufs[5]);
    if (vb_bamass_cigar.len) 
        buf_move (vb, vb->codec_bufs[5], "codec_bufs[5]", vb_bamass_cigar);

    PosType64 gpos;
    fastq_bamass_retrieve_ent (vb, e, true, NULL, &gpos, NULL, NULL, NULL, NULL);

    PosType32 pos;
    WordIndex ref_index = ref_contig_get_by_gpos (gpos, 0, &pos, false);

    STR(rname);
    rname = ref_contigs_get_name (ref_index, &rname_len);

    // sanity
    ASSERT (qname_hash >> 10 == e->z_qname_hash_hi, "qname_hash>>10=%016"PRIx64" ≠ %016"PRIx64, qname_hash >> 10, (uint64_t)e->z_qname_hash_hi);

    SNPRINTF (s, "hash={%016"PRIx64",%08x} fwd=%u gpos=%-10"PRIu64" sam={ %-5.*s %-10u %-20.256s }", 
              qname_hash, e->seq_hash, e->z_is_forward, gpos, MIN_(rname_len,32), rname, pos,
              display_bamass_cigar(vb, e->z_is_forward).s); // destroys vb_bamass_cigar (so we saved it)

    buf_free (vb_bamass_cigar);

    // restore
    if (vb->codec_bufs[5].len) 
        buf_move (vb, vb_bamass_cigar, "codec_bufs[5]", vb->codec_bufs[5]);

    return s;
}

// called by main thread after each FASTQ file (or pair of files) is done compressing with bamass
void fastq_bamass_zip_finalize (bool is_last_fastq) 
{
    if (flag_show_deep) {
        iprintf ("\nZIP: %s%s reads breakdown by bammassability:\n", z_dt_name_faf(), IS_R2 ? " (R1+R2)" : "");
        
        for (int i=BA_AFTER+1; i < NUM_DEEP_STATS_ZIP; i++) 
            if (z_file->deep_stats[i])
                iprintf ("%-11.11s: %"PRIu64" (%.1f%%)\n", (rom[])DEEP_STATS_NAMES_ZIP[i], z_file->deep_stats[i], 100.0 * (double)z_file->deep_stats[i] / (double)z_file->deep_stats[NDP_FQ_READS]);
    }

    // if we have more fastq files with the same bamass - save cigar huffman as it is a different z_file
    if (!is_last_fastq) {
        buf_move (evb, save_cigar_huff, NULL, ZCTX(FASTQ_CIGAR)->huffman);
        return;
    }

    if (flag_show_deep) 
        iprintf ("\nZIP: BAM bamassable alignments=%"PRIu64" consumed=%"PRIu64" not_consumed=%"PRIu64"\n", 
                 bamass_ents.len, global_num_consumed, bamass_ents.len);

    if (!flag.let_OS_cleanup_on_exit) {
        buf_destroy (bamass_ents);
        buf_destroy (bamass_alns);
        buf_destroy (bamass_heads);
    }
}

DeepStatsZip fastq_seg_find_bamass (VBlockFASTQP vb, ZipDataLineFASTQ *dl, DeepHash *deep_hash, STRp(seq),  
                                    BamAssEnt **matching_ent) // out
{
    START_TIMER;
    #define RETURN(x) ({ COPY_TIMER (fastq_seg_find_bamass); return (x); })
    
    ARRAY (BamAssEnt, ents, bamass_ents);

    for (uint32_t ent_i = *B32(bamass_heads, deep_hash->qname & bitmask64 (num_hash_bits)); 
         ent_i != NO_NEXT; 
         ent_i = ents[ent_i].next) {

        BamAssEnt *e = &ents[ent_i];

        if (e->next != NO_NEXT) __builtin_prefetch (&ents[e->next]); // prefetch next entry

        if (e->z_qname_hash_hi != deep_hash->qname >> 10) // test if high 54 bits are the same ; low 10 bits already known to be the same as on the same linked list (num_hash_bits >= 10)
            continue; // possible, because hash uses less bits that e->hash.qname
            
        fastq_bamass_retrieve_ent (VB, e, true, &vb->is_forward, &vb->gpos, &dl->sam_seq_len, &vb->ref_consumed, 
                                   &vb->ref_and_seq_consumed, &vb->insertions);

        // change: fastq_bamass_retrieve_ent changed S->M to account for trimming, and now the alignment falls of the genome
        if (vb->gpos < 0 || vb->gpos + vb->ref_consumed > *p_genome_nbases)
            continue;

        if (dl->sam_seq_len > seq_len || // definitely not a match if FASTQ seq_len is shorter than SAM's   
            (dl->sam_seq_len < seq_len && !segconf.deep_has_trimmed)) // not a match if FASTQ seq_len is too long, and we don't allow trimming
             continue;     

        // case: QNAME matches: see if we can find a subsequence of SEQ that matches
        if (!fastq_deep_seg_find_subseq (vb, STRa(seq), dl->sam_seq_len, e->seq_hash, segconf.deep_has_trimmed_left, &dl->sam_seq_offset))
            continue; 

        // update to subsequence found (for messages)
        deep_hash->seq = e->seq_hash;
        
        *matching_ent = e;

        // represent left and right trims as enlarging the S or M ops in the cigar (flanking S/M op added in fastq_bamass_retrieve_ent)
        uint32_t left_trim  = dl->sam_seq_offset;
        uint32_t right_trim = seq_len - dl->sam_seq_len - dl->sam_seq_offset;
        
        cigar_L_flank->n += left_trim;
        cigar_R_flank->n += right_trim;

        if (cigar_L_flank->op == BC_M) {
            if (vb->is_forward) vb->gpos -= left_trim; 
            vb->ref_consumed += left_trim;
            vb->ref_and_seq_consumed += left_trim;
        }

        if (cigar_R_flank->op == BC_M) {
            if (!vb->is_forward) vb->gpos -= right_trim; 
            vb->ref_consumed += right_trim;
            vb->ref_and_seq_consumed += right_trim;
        }

        RETURN ((dl->sam_seq_len == seq_len) ? NDP_DEEPABLE : NDP_DEEPABLE_TRIM);
    }

    RETURN (NDP_NO_MATCH);
    
    #undef RETURN
}

// true if next seq-consuming op is I
static inline bool next_op_is_I (BamCigarOpP cigar, BamCigarOpP after_cigar)
{
    for (cigar++; cigar < after_cigar; cigar++) 
        switch (cigar->op) {
            case BC_I:            return true;
            case BC_M: case BC_S: return false; // note: bamass doesn't have X,=,H,P,N
            case BC_D:            continue;     // skip D
        }
    
    return false; // there are no further seq-consuming ops
}

MappingType fastq_bamass_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(seq), bool is_pair_2, PosType64 pair_gpos, bool pair_is_forward)
{
    declare_seq_contexts;   
    START_TIMER;

    BitsP bitmap = (BitsP)&bitmap_ctx->local;

    ASSERTW (seq_len < 100000 || segconf_running || segconf.is_long_reads, 
             "%s: Warning: fastq_bamass_seg_SEQ: seq_len=%u is suspiciously high and might indicate a bug", LN_NAME, seq_len);

    RefLock lock = IS_REF_EXTERNAL ? REFLOCK_NONE : ref_lock (vb->gpos, vb->ref_consumed);

    buf_alloc_bits (vb, &bitmap_ctx->local, vb->ref_and_seq_consumed, vb->lines.len32 / 16 * segconf.std_seq_len,
                    SET/*initialize to "no mismatches"*/, CTX_GROWTH, CTX_TAG_LOCAL); 

    for (int i=0; i < 4; i++) 
        buf_alloc (vb, &seqmis_ctx[i].local, vb->ref_and_seq_consumed, 0, char, CTX_GROWTH, CTX_TAG_LOCAL); 

    if (segconf.use_insertion_ctxs)    
        for (int i=0; i < 4; i++) 
            buf_alloc (vb, &seqins_ctx[i].local, vb->insertions, 0, char, CTX_GROWTH, CTX_TAG_LOCAL); 

    buf_alloc (vb, &nonref_ctx->local, seq_len + 3, 0, uint8_t, CTX_GROWTH, CTX_TAG_LOCAL); 

    bitmap_ctx->local_num_words++;

    // get reference sequence
    ASSERTNOTINUSE (vb->scratch);
    buf_alloc_exact (vb, vb->scratch, vb->ref_consumed, char, "scratch"); 
    rom ref = ref_get_textual_seq (vb->gpos, B1STc(vb->scratch), vb->ref_consumed, !vb->is_forward);
    
    rom after_ref = ref + vb->ref_consumed;
    rom after_seq = seq + seq_len;

    if (IS_REF_EXT_STORE) 
        bits_set_region (ref_get_genome_is_set(), vb->gpos, vb->ref_consumed); // we will need this ref to reconstruct

    aligner_seg_gpos_and_fwd (VB, vb->gpos, vb->is_forward, is_pair_2, pair_gpos, pair_is_forward);

    bool perfect=true;
    
    for_cigar (vb_bamass_cigar) {
        case BC_M: {
            uint32_t n = op->n;

            ASSSEG (seq + n <= after_seq, "seq out of range: n=%u seq_len=%u cigar=%s", n, seq_len, display_bamass_cigar(VB, unknown).s);
            ASSSEG (ref + n <= after_ref, "ref out of range: n=%u ref_consumed=%u cigar=%s", n, vb->ref_consumed, display_bamass_cigar(VB, unknown).s);

            uint32_t bit_i = bitmap_ctx->next_local; // copy to automatic variable for performance
            
            while (n) {
                // case: our seq is identical to the reference at this site
                if (*seq == *ref) 
                    bit_i++; // bit remains set. 

                // case: ref is set to a different value - we store our value in seqmis_ctx
                else {
                    uint8_t ref_base_2bit = nuke_encode(*ref);
                            
                    BNXTc (seqmis_ctx[ref_base_2bit].local) = *seq;
                    bits_clear (bitmap, bit_i++);
                    perfect = false;
                } 

                n--; seq++; ref++;
            }
            
            bitmap_ctx->next_local = bit_i;
            break; 
        }
        case BC_I: case BC_S:
            ASSSEG (seq + op->n <= after_seq, "seq out of range: n=%u seq_len=%u cigar=%s", op->n, seq_len, display_bamass_cigar(VB, unknown).s);

            if (op->op == BC_I && segconf.use_insertion_ctxs) {
                // mux by base after insertion in SEQ. 
                // note: the byte after SEQ, \n or \r, mapped to 0, so no need to test explicitly
                // note: if the next seq-consuming op is also I (e.g. in RNA: 621I2611N310I), map to 0, as we won't yet know seq[n] yet during reconstruction 
                int ins_ctx_i = (next_op_is_I (op, after_op) ? 0 : acgt_encode[(int8_t)seq[op->n]]);

                buf_add (&seqins_ctx[ins_ctx_i].local, seq, op->n);
            }
            else 
                buf_add (&nonref_ctx->local, seq, op->n);

            seq += op->n;
            break;

        case BC_D: 
            ASSSEG (ref + op->n <= after_ref, "ref out of range: n=%u ref_consumed=%u cigar=%s", 
                    op->n, vb->ref_consumed, display_bamass_cigar(VB, unknown).s);

            ref += op->n;
            break;

        default: // bamass cigar doesn't have H,P,X,=,N
            ABOSEG ("Invalid CIGAR op=%u", op->op);        
    }

    ref_unlock (&lock); // does nothing if REFLOCK_NONE
    buf_free (vb->scratch);

    // an error in the first can indicate that CIGAR is inconsistent with sequence, and in the 2nd it is a bug
    ASSSEG (seq == after_seq, "seq_len=%u but consumed only %d: cigar=%s", seq_len, (int)seq_len - (int)(after_seq - seq), display_bamass_cigar(VB, unknown).s);
    ASSSEG (ref == after_ref, "ref_consumed=%u but consumed only %d: cigar=%s", vb->ref_consumed, (int)vb->ref_consumed - (int)(after_ref - ref), display_bamass_cigar(VB, unknown).s);

    COPY_TIMER (fastq_bamass_seg_SEQ);

    if (perfect) {
        bitmap_ctx->next_local -= vb->ref_and_seq_consumed; // note: we truncate bitmap, if needed, in fastq_seg_finalize
        return MAPPING_PERFECT; // no mismatches
    }

    else
        return MAPPING_ALIGNED;
}

void fastq_bamass_seg_CIGAR (VBlockFASTQP vb)
{
    START_TIMER
    decl_ctx (FASTQ_CIGAR);

    cigar_remove_flanking_0_ops (&vb_bamass_cigar);

    BamCigarOpP largest_op = NULL;
    uint32_t save_largest_n = 0;

    // seq_len is calculable from QNAME or std_seq_len - improve cigar compression by dropping largest N   
    // (note: cigar is guaranteed by sam_prepare_deep_cigar to not contain any n=0 op) 
    if ((segconf.seq_len_dict_id.num  && vb->seq_len == ECTX(segconf.seq_len_dict_id)->last_value.i)
     || (!segconf.seq_len_dict_id.num && vb->seq_len == segconf.std_seq_len)) {

        // find largest seq-consuming op
        for_buf (BamCigarOp, op, vb_bamass_cigar) 
            if ((!largest_op || op->n > largest_op->n) && op->op != BC_D) // op is M, I or S
                largest_op = op;

        save_largest_n = largest_op->n;
        largest_op->n = 0;
    }

    ASSERTNOTINUSE (vb->scratch);
    sam_cigar_binary_to_textual (VB, CIG(vb_bamass_cigar), false, &vb->scratch);

    if (largest_op)
        largest_op->n = save_largest_n; // restore

    // case: long CIGAR - to local
    if (vb_bamass_cigar.len32 >= 3) 
        seg_add_to_local_string (VB, ctx, STRb(vb->scratch), LOOKUP_SIMPLE, 0);

    // case: short cigar
    else 
        seg_by_ctx (VB, STRb(vb->scratch), ctx, 0); 

    buf_free (vb->scratch);

    COPY_TIMER(fastq_bamass_seg_CIGAR);
}

//*******************************************************
// PART 3: piz FASTQ segged with bamass
//*******************************************************

static inline void update_by_op (VBlockFASTQP vb, ConstBamCigarOpP op)
{
    switch (op->op) {
        case BC_M : vb->ref_consumed += op->n; vb->seq_len += op->n; vb->ref_and_seq_consumed += op->n; break;
        case BC_D : vb->ref_consumed += op->n; break;
        case BC_I : case BC_S : vb->seq_len += op->n; break;
        default   : ABORT_PIZ ("invalid cigar_op=%u", op->op); // note: H,P,=,X,N are never included in bamass cigars
    }   
}
// sets vb_bamass_cigar, vb->seq_len, vb->ref_consumed, vb->ref_and_seq_consumed
static void fastq_bamass_recon_cigar (VBlockFASTQP vb)
{
    char *recon = BAFTtxt;
    reconstruct_from_ctx (vb, FASTQ_CIGAR, 0, true); // reconstruct textual cigar

    sam_cigar_textual_to_binary (VB, recon, BAFTtxt - recon, &vb_bamass_cigar, "bamass_cigar");
    Ltxt = BNUMtxt (recon); // remove reconstructed textual cigar

    vb->seq_len = 0; // possibly set earlier by QNAME length

    BamCigarOpP zero_n_op = NULL;

    for_buf (BamCigarOp, op, vb_bamass_cigar) {
        if (!op->n) zero_n_op = op;
        update_by_op (vb, op);
    }

    // recover n=0 to actual value
    if (zero_n_op) {
        uint32_t real_seq_len = segconf.seq_len_dict_id.num ? reconstruct_peek_by_dict_id (VB, segconf.seq_len_dict_id, 0, 0).i // peek, since length can come from either line1 or line3
                                                            : segconf.std_seq_len;
        zero_n_op->n = real_seq_len - vb->seq_len;
        update_by_op (vb, zero_n_op);
    }
}

SPECIAL_RECONSTRUCTOR_DT (fastq_special_SEQ_by_bamass)
{
    START_TIMER;

    VBlockFASTQP vb = (VBlockFASTQP)vb_;
    declare_seq_contexts;

#ifdef DEBUG // this function is a performance hotspot, therefore this code is only a compilation, not runtime, option
    #define verify_is_set(ref) ASSPIZ (!flag.debug || IS_REF_EXTERNAL || ((ref)-start_ref < bitmap_ctx->piz_is_set.len32 && *Bc(bitmap_ctx->piz_is_set, (ref)-start_ref) == 1), \
                                       "Expecting POS=%d + base_i=%u to have is_set=1", (PosType32)vb->contexts[SAM_POS].last_value.i, (int)((ref)-start_ref))
    #define verify_is_set_n(ref, n) ({ for (uint32_t i=0; i < (n); i++) verify_is_set ((ref) + i); })
#else
    #define verify_is_set(ref)
    #define verify_is_set_n(ref, n)
#endif

    if (!bitmap_ctx->is_loaded) return NO_NEW_VALUE; // if case we need to skip the SEQ field (for the entire VB)

    bool is_perfect = snip[0] == '1';

    ASSERT0 (is_perfect || bitmap_ctx->local.len32 || nonref_ctx->local.len32, "No SEQ data, perhaps sections were skipped?");

    bitmap_ctx->last_value.i = bitmap_ctx->next_local; // for SEQ, we use last_value for storing the beginning of the sequence

    // set vb_bamass_cigar, vb->seq_len, vb->ref_consumed
    fastq_bamass_recon_cigar (vb);

    if (flag.collect_coverage) {
        fastq_update_coverage_aligned (vb);
        goto done;
    }

    if (vb->R1_vb_i) // R2 
        fastq_piz_R1_test_aligned (vb); // set r1_is_aligned

    // get gpos and is_forward
    aligner_recon_get_gpos_and_fwd (VB, vb->R1_vb_i > 0, &vb->gpos, &vb->is_forward);

    // get reference - the needed ref_consumed bases - revcomped if needed
    ASSERTNOTINUSE (vb->scratch);
    buf_alloc_exact (vb, vb->scratch, vb->ref_consumed, char, "scratch"); 
    rom start_ref = ref_get_textual_seq (vb->gpos, B1STc(vb->scratch), vb->ref_consumed, !vb->is_forward); 
    rom ref = start_ref;

    char *nonref      = Bc(nonref_ctx->local, nonref_ctx->next_local); // possibly, this VB has no nonref (i.e. everything is ref), in which case nonref would be an invalid pointer. That's ok, as it will not be accessed.
    char *recon       = BAFTtxt;
    rom start_recon   = recon;
    bool deferred_ins = false;

    if (!IS_REF_EXTERNAL && flag.debug) {
        ref_get_is_set_bytemap (VB, vb->gpos, vb->ref_consumed, !vb->is_forward, &bitmap_ctx->piz_is_set, "piz_is_set"); // we test is_set only in debug
        verify_is_set_n (ref, vb->ref_consumed);
    }

    for_cigar (vb_bamass_cigar) {
        case BC_M : // this might have been one or more X/= in the original SAM cigar
            memcpy (recon, ref, op->n);

            if (is_perfect) {
                recon += op->n;
                ref   += op->n;
            }

            else 
                for (uint32_t i=0; i < op->n; i++) {
                    uint32_t matches = bits_get_run ((BitsP)&bitmap_ctx->local, bitmap_ctx->next_local, op->n - i);

                    bool has_mismatch = (i + matches < op->n);
                    recon += matches;
                    ref   += matches;
                    bitmap_ctx->next_local += matches + has_mismatch;
                    i += matches - 1 + has_mismatch;

                    // case: we have a mismatch: overwrite base copied from reference with mismatch
                    if (has_mismatch) {
                        uint8_t ref_base_2bit = nuke_encode(*ref++);
                        ContextP mis_ctx = &seqmis_ctx[ref_base_2bit];
                        *recon++ = *Bc(mis_ctx->local, mis_ctx->next_local++);
                    }
                }
            break;

        case BC_I : 
            if (segconf.use_insertion_ctxs) {
                recon += op->n; // skip for now, we will fill in the gaps later
                deferred_ins = true;
                break;
            }
            // fallthrough

        case BC_S : // this might have been an I in the original SAM cigar
            if (op->n == 1)  // shortcut for most common case
                *recon++ = *nonref++;

            else { // larger op.n - we are better off with memcpy
                recon = mempcpy (recon, nonref, op->n);
                nonref += op->n;
            }
            break;

        case BC_D : // this might have been an N in the original SAM cigar
            ref += op->n;
            break;
    }

    ASSPIZ (recon - start_recon == vb->seq_len,  "expecting seq_consumed(%u) == vb->seq_len(%u)", (int)(recon - start_recon), vb->seq_len);
    ASSPIZ (ref - start_ref == vb->ref_consumed, "expecting ref_consumed(%u) == vb->ref_consumed(%u)", (int)(ref - start_ref), vb->ref_consumed);

    // case: insertions were muxed by the base after - we can fill them in now
    if (deferred_ins) {
        SAFE_ASSIGN (recon, 0); // in case last insertion goes to the end of SEQ - it will use this "base" for muxing
        recon -= vb->seq_len;   // go back to the start of the SEQ reconstruction

        for_cigar (vb_bamass_cigar) {        
            case BC_I: {
                int ins_ctx_i = next_op_is_I (op, after_op) ? 0 : acgt_encode[(int8_t)recon[op->n]];
                ContextP ins_ctx = &seqins_ctx[ins_ctx_i]; // mux by base after insertion in SEQ. 
                
                ASSPIZ (ins_ctx->next_local + op->n <= ins_ctx->local.len32, NEXT_ERRFMT " op_i=%u cigar=\"%s\"", 
                        __FUNCTION__, ins_ctx->tag_name, ins_ctx->next_local, op->n, ins_ctx->local.len32, BNUM(vb_bamass_cigar,op), display_bamass_cigar(VB, unknown).s); 
       
                recon = mempcpy (recon, Bc(ins_ctx->local, ins_ctx->next_local), op->n);
                
                ins_ctx->next_local += op->n;
                break;
            }

            case BC_M: case BC_S: 
                recon += op->n;
                break;

            case BC_D: break; 
        }

        SAFE_RESTORE;
    }
    
    Ltxt = BNUMtxt (recon);
    nonref_ctx->next_local = BNUM (nonref_ctx->local, nonref);

done:
    buf_free (vb->scratch); 
    buf_free (vb_bamass_cigar);
    if (flag.debug) buf_free (bitmap_ctx->piz_is_set); 

    COPY_TIMER (fastq_special_SEQ_by_bamass);
    return NO_NEW_VALUE;
}
