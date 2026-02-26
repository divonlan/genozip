// ------------------------------------------------------------------
//   bai.c
//   Copyright (C) 2026-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <errno.h>
#include "bai.h"
#include "file.h"
#include "sam_friend.h"
#include "seg.h"
#include "mgzip.h"
#include "endianness.h"
#include "codec.h"
#include "profiler.h"
#include "reconstruct.h"

// note: to see the BINs in a BAM file: show_BIN.sh file.bam

#define BAI_MAGIC "BAI\1"
#define TBI_MAGIC "TBI\1"

#define MAX_BAI_CONTIG_LEN (512 MB - 1) // sam spec §5.1.1

// Bins: 
// Level Range   Bin_ID     Num_Bins
// 0     512 MB  0	        1
// 1     64 MB   1–8        8
// 2     8 MB    9–72       64
// 3     1 MB    73–584	    512
// 4     128 KB  585–4680   4096
// 5     16 KB   4681–37448 32768
#define NUM_BINS  37451 // per spec 
#define BIN_STATS 37450 // pseudo-bin for holding contig stats 

static bool bai_is_aborted;

FileType index_type = UNKNOWN_FILE_TYPE; // BAI or TBI
#define IS_BAI (index_type == BAI)
#define IS_TBI (index_type == TBI)

typedef union {
    struct { uint64_t off_in_bb : 16; uint64_t bb_i_in_vb    : 48; }; // phase 1: index by bb_i in VB
    struct { uint64_t xx1       : 16; uint64_t bb_off_in_vb  : 48; }; // phase 2: index by bb offset in VB
    struct { uint64_t xx2       : 16; uint64_t bb_off_in_txt : 48; }; // phase 3: index by bb offset in txt file (.bam, .sam.gz, .vcf.gz)
    uint64_t value;
} VirtualFileOffset;

typedef struct {
    VirtualFileOffset start, after; // a [VFO₁ — VFO₂) range
} VirtualFileRange;

#define VOFF_FMT "%"PRIu64"⁀%u"
#define VOFF_RNG_FMT "↕=["VOFF_FMT" — "VOFF_FMT")"
#define VOFF_VAL(x) (uint64_t)x.bb_i_in_vb, (uint32_t)x.off_in_bb
#define VOFF_RNG_VAL(x) (uint64_t)(x).start.bb_i_in_vb, (uint32_t)(x).start.off_in_bb, (uint64_t)(x).after.bb_i_in_vb, (uint32_t)(x).after.off_in_bb

// a chunk is a consecutive group of alignments, which all have BIN=bin, or BIN=(a bin fully contained in bin)
typedef struct { // 40 bytes
    BaiBinType bin;
    uint16_t rname_lo; 
    uint32_t num_lines;     // number of alignment in this chunk
    union { 
    VirtualFileOffset nxtbb;// the next bb after "after", except if after.off_in_bb==0, in which case its the same as after (note: field off_in_bb is not used)
    struct { uint64_t rname_hi : 16; uint64_t xx3 : 48; }; // utilize the unused off_in_bb bits
    };
    VirtualFileRange;       // VFO of first alignment in chunk and alignment after chunk (even if beyond end of data)
} BaiChunk, *BaiChunkP;
#define CHUNK_RNAME(c) ((((WordIndex)(c)->rname_hi) << 16) | (WordIndex)(c)->rname_lo)

// VFO of alignment which is the first alignment to overlap a certain 16 KBp window
typedef struct { // 16 bytes
    WordIndex rname;
    uint32_t window;        // 0-based pos in the rname contig, divisible by 16384
    VirtualFileOffset start;
} BaiLinear;

typedef struct { // stats for one rnames
    WordIndex rname;
    uint64_t n_mapped;      // number of alignments with UNMAPPED flag not set
    uint64_t n_unmapped;    // number of alignments with UNMAPPED set (placed if RNAME, and unplaced if no RNAME)
    VirtualFileOffset ref_beg, ref_end; // [VFO₁ — VFO₂) range of reads of this contig (whether mapped or not)
} BaiStats;

typedef struct __attribute__((packed, aligned(1))) {
    char magic[4];      // "TBI\1"
    uint32_t n_ref;     // number of contig names in names
    uint32_t format;    // 0: generic; 1: SAM; 2: VCF           (VCF: 2)
    uint32_t col_rname; // Column for the rname                 (VCF: 1)  
    uint32_t col_pos;   // Column for the start of a region     (VCF: 2)
    uint32_t col_end;   // Column for the end of a region       (VCF: 0)
    uint32_t meta;      // Leading character for comment lines  (VCF: '#')
    uint32_t skip;      // # lines to skip at the beginning     (VCF: 0)
    uint32_t names_len; // length of names array (inc. \0s)
    char names[];       // concatenated contig names, each terminated with \0
} TbiHeader;

static int32_t num_bai_contigs; // note: signed to avoid unsigned comparisons with rname=-1

// calculate the shortest bin (i.e. highest bin value) which includes the entire range [first_pos,last_pos].
// code adapted from SAM specification §5.3 and Tabix specification
// NOTE: sam_cigar_special_CIGAR relies on this calculation in PIZ, so if ever changing, need to take care of backcomp.
BaiBinType bai_reg2bin (PosType32 first_pos, PosType32 last_pos) // 1-based closed interval [first_pos,last_pos]
{
    first_pos--; last_pos--; // make it 0-based
      
    if (first_pos>>14 == last_pos>>14) return ((1<<15)-1)/7 + (first_pos>>14);
    if (first_pos>>17 == last_pos>>17) return ((1<<12)-1)/7 + (first_pos>>17);
    if (first_pos>>20 == last_pos>>20) return ((1<<9 )-1)/7 + (first_pos>>20);
    if (first_pos>>23 == last_pos>>23) return ((1<<6 )-1)/7 + (first_pos>>23);
    if (first_pos>>26 == last_pos>>26) return ((1<<3 )-1)/7 + (first_pos>>26);
    return 0;
}

bool bai_is_native_indexing_supported (DataType dt)
{
    return dt == DT_BAM || dt == DT_SAM || dt == DT_VCF;
}

void bai_set_show_bai (rom optarg)
{
    if (!optarg)
        flag.show_bai = is_genozip ? SHOW_BAI_UNSORTED : SHOW_BAI_CHUNKS; // default mode

    else {
        flag.show_bai = (optarg && !strcmp (optarg,"sort"))     ? SHOW_BAI_SORT 
                      : (optarg && !strcmp (optarg,"chunks"))   ? SHOW_BAI_CHUNKS 
                      : (optarg && !strcmp (optarg,"raw"))      ? SHOW_BAI_RAW 
                      : (optarg && !strcmp (optarg,"unsorted")) ? SHOW_BAI_UNSORTED
                      : (optarg && !strcmp (optarg,"linear"))   ? SHOW_BAI_LINEAR
                      :                                           SHOW_BAI_NONE;

        if ((is_genounzip && !(flag.show_bai == SHOW_BAI_CHUNKS || flag.show_bai == SHOW_BAI_RAW || flag.show_bai == SHOW_BAI_LINEAR)) // genounzip
        ||  (is_genozip   && !(flag.show_bai == SHOW_BAI_SORT   || flag.show_bai == SHOW_BAI_UNSORTED || flag.show_bai == SHOW_BAI_CHUNKS))) // genozip --show-bai
            ABORTINP ("unrecognized --show-bai mode: \"%s\"", optarg);
    }
}

// PIZ: can be called at any time by any thread to cancel bai creation
static void bai_abort (rom reason, STRp(qname), bool debug_bai_only)
{
    // only one thread gets to display abortion reason
    if (!__atomic_test_and_set (&bai_is_aborted, __ATOMIC_RELAXED)) {
        if (!debug_bai_only || flag.debug_bai) {
            if (qname) 
                WARN ("FYI: Not creating a %s file. Reason: %s. Offending %s: %.*s ", 
                      IS_BAI?".bai":".tbi", reason, DTPZ (line_name), STRf(qname));
            else
                WARN ("FYI: Not creating a %s file. Reason: %s.", IS_BAI?".bai":".tbi", reason);
        }
    }

    flag.make_bai = false; // can propagate to other threads at its leisure
} 
#define BAI_ABORT(reason) { bai_abort (reason, STRa(qname), false); return; }

// PIZ main thread: called after reconstructing SAM header
void bai_initialize (rom txt_header/*NUL-terminated*/, uint32_t num_header_contigs, PosType64 len_of_longest_contig)
{
    ASSERTNOTNULL (txt_file); // must be called after txt_file initialized, eg from inspect_txt_header

    if (Z_DT(SAM)) { // SAM or BAM
        num_bai_contigs = num_header_contigs; // note: in VCF, TBI contigs come from the data, not VCF header
        index_type = BAI;
    }
    else if (Z_DT(VCF)) {
        num_bai_contigs = ZCTX(VCF_CHROM)->word_list.len32; // initial, to be updated later based on the actual number of CHROMs in the file
        index_type = TBI;
    }
    else
        ABORT ("file type %s not supported for native indexing yet", z_dt_name());

    // give a reason if BAI is not created due to the details of the BAM file, but not if user is not likely to expect it anyway
    if (txt_file->effective_codec == CODEC_BGZF && // .sam.gz or .bam
        !txt_file->redirected && txt_file->name && !flag.no_writer && !flag.no_index) {
        
        if (flag.bgzf == BGZF_BY_ZFILE)  // too complicated to calculate BAI if bgzf=exact, because a BGZF block can span two VBs, and also original .bai is still good
            bai_abort ("using --bgzf=exact", 0, 0, false);

        else if (Z_DT(SAM) && !num_bai_contigs)
            bai_abort ("BAM header has no SQ records", 0, 0, false);

        else if (Z_DT(SAM) && len_of_longest_contig > MAX_BAI_CONTIG_LEN)
            bai_abort ("BAM contains a contig longer than 512 Mbps", 0, 0, false);

        else if (Z_DT(SAM) && !segconf.is_sorted) 
            bai_abort ("BAM is not sorted (1)", 0, 0, false); 
        
        else if (Z_DT(VCF) && strstr (txt_header, "##FORMAT=<ID=LEN,"))
            // TO DO: files with FORMAT/LEN - LEN needs to be considered (bug 1229). For now we abort.
            bai_abort ("Genozip limitation: cannot index a VCF with FORMAT/LEN. Please contact "EMAIL_SUPPORT, 0, 0, false);

        else {
            flag.make_bai = true;
            bai_is_aborted = false;

            buf_alloc_exact_zero (evb, txt_file->bai_stats, num_bai_contigs, BaiStats, "txt_file->bai_stats");
            
            buf_set_promiscuous (&txt_file->bai_chunks, "txt_file->bai_chunks");   // allow other threads to enlarge
            buf_alloc (evb, &txt_file->bai_chunks, 10000, 0, BaiChunk,  0, NULL);  // initial allocation

            buf_set_promiscuous (&txt_file->bai_linear, "txt_file->bai_linear");   
            buf_alloc (evb, &txt_file->bai_linear, 10000, 0, BaiLinear, 0, NULL); 
        }
    }
}

static inline rom bai_get_line_vcf (rom next, rom after_vb,
                                    WordIndex *rname, PosType32 *pos, 
                                    BaiBinType *bin, uint32_t *ref_consumed, pSTRp(id))
{
    typedef enum { RNAME, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT } VcfFields __attribute__((unused)); // quick way to define constants
    extern uint32_t vcf_num_samples; 

    // get the fields efficiently (10x faster than str_split_by_tab) taking advantage of the data already tested during seg and verified in recon
    mSTR (fld, 8);
    for (int i=0; i < 7 + (vcf_num_samples > 0); i++) { // first 7 columns always end with \t, INFO ends with \t iff file has samples
        flds[i] = next;
        while (*next != '\t') next++;
        fld_lens[i] = next - flds[i];
        next++; // past \t
    }

    // INFO is the last field if no samples
    if (!vcf_num_samples) {
        flds[INFO] = next;
        while (*next != '\n' && *next != '\r') next++;
        fld_lens[INFO] = next - flds[INFO];

        next--; // move back so next while loop will find \n regardless on whether we are on a \r or \n
    }

    while (*next != '\n') next++; // move past samples

    *id     = flds[ID];
    *id_len = fld_lens[ID];
    *rname  = ctx_get_word_index_by_snip (VCF_CHROM, STRi(fld, RNAME));

    ASSERT (*rname != WORD_INDEX_NONE, "CHROM=%.*s not found in CHROM dict", STRfi(fld,RNAME));

    ASSERT (str_get_int_range32 (STRfld (POS), 1, MAX_POS32, pos),
            "Invalid POS \"%.*s\": expecting an integer [1,%d]", STRfi (fld, POS), (int)MAX_POS32);
    (*pos)--; // convert to 0-based.

    // ref_consumed is as defined in the rlen field in the VCF 4.5 spec §6.3.1
    // if END exists, it is the same as rlen. starting VCF 4.5, per §1.6.1 rlen/END are calculatable as:
    // "maximum end reference position (1-based) of: VCF spec §6.3.1
    // - the position of the final base of the REF allele
    // - the end position corresponding to the SVLEN of a symbolic SV allele,
    // - end positions calculated from FORMAT LEN for the <*> symbolic allele."
    
    *ref_consumed = fld_lens[REF]; // at least this much

    SAFE_NUL (&flds[INFO][fld_lens[INFO]]);
    rom info;

    // If there is an END (=rlen), use it
    PosType32 end, svlen;
    if (     (info = strstr (flds[INFO], "END="))   && (info[-1]=='\t' || info[-1]==';') && (end   = atoi (info + STRLEN("END="))) && end-1 >= *pos)
        MAXIMIZE (*ref_consumed, 1+ ((end-1/*0-based*/) - *pos)); // both start and end position are counted

    // if no END, check for SVLEN
    else if ((info = strstr (flds[INFO], "SVLEN=")) && (info[-1]=='\t' || info[-1]==';') && (svlen = atoi (info + STRLEN("SVLEN="))))
        MAXIMIZE (*ref_consumed, ABS(svlen)); // absolute value per VCF spec 4.5 §3

    // TO DO maximize also with FORMAT/LEN (bug 1229)

    SAFE_RESTORE;

    PosType32 last_pos = (*pos + *ref_consumed - 1); 
    *bin = bai_reg2bin (*pos + 1, last_pos + 1); 

    return next + 1; 
}

static inline rom bai_get_line_sam (rom next, rom after_vb,
                                    SamFlags *sam_flags, WordIndex *rname, PosType32 *pos, 
                                    BaiBinType *bin, uint32_t *ref_consumed, pSTRp(qname))
{
    // get the fields efficiently (10x faster than str_split_by_tab) taking advantage of the data already tested during seg and verified in recon
    mSTR (fld, 6);
    for (int i=0; i < 6; i++) {
        flds[i] = next;
        while (*next != '\t') next++;
        fld_lens[i] = next - flds[i];
        next++; // past \t
    }

    *qname     = flds[QNAME];
    *qname_len = fld_lens[QNAME];
    
    if (IS_ASTERISKi(fld,RNAME))
        *rname = -1; // "*" in BAM terms
    else {
        *rname = sam_get_contig_by_name (STRi(fld,RNAME)); 
        ASSERT (*rname != WORD_INDEX_NONE, "RNAME=%.*s not found in SAM header. Use --no-bai", STRfi(fld,RNAME)); // should never happen - tested in seg
    }

    ASSERT (str_get_int_range32 (STRfld (POS), 0, MAX_POS_SAM, pos),
            "Invalid POS \"%.*s\": expecting an integer [0,%d]", STRfi (fld, POS), (int)MAX_POS_SAM);
    (*pos)--; // convert to 0-based. note: "0" (unmapped) becomes -1

    ASSERT (str_get_int_range16 (STRfld (FLAG), 0, SAM_MAX_FLAG, &sam_flags->value), 
            "Invalid FLAG field: \"%.*s\"", STRfi (fld, FLAG));

    *ref_consumed = sam_cigar_get_ref_consumed (STRi(fld,CIGAR), false, true);

    PosType32 last_pos = sam_flags->unmapped ? *pos : (*pos + *ref_consumed - 1); // note: it is possible to have RNAME/POS set and FLAG.unmapped. The SAM spec calls this "placed unmapped"; 
    *bin = bai_reg2bin (*pos + 1, last_pos + 1); 

    while (*next != '\n') next++;
    return next + 1; 
}

static inline rom bai_get_line_bam (rom next_field, 
                                    SamFlags *sam_flags, WordIndex *rname, PosType32 *pos, 
                                    BaiBinType *bin, uint32_t *ref_consumed, pSTRp(qname))
{
    // get the data we need from the alignment: rname, pos and bin
    uint32_t block_size = NEXT_UINT32;
    rom after_aln = next_field + block_size;

    *rname              = NEXT_UINT32;
    *pos                = NEXT_UINT32;  // 0-based
    *qname_len          = NEXT_UINT8 - 1; // -1 to exclude the \0
    next_field++;                       // skip MAPQ
    *bin                = NEXT_UINT16;
    uint16_t n_cigar_op = NEXT_UINT16;
    sam_flags->value    = NEXT_UINT16; 
    
    next_field += 4 * sizeof(uint32_t); // skip l_seq, RNEXT, PNEXT, TLEN
    
    *qname = next_field;
    next_field += *qname_len + 1;       // skip read_name and \0

    ASSERT0 (next_field + sizeof (BamCigarOp) * n_cigar_op < after_aln, "CIGAR: txt_data overflow");

    *ref_consumed = sam_cigar_get_ref_consumed (next_field, n_cigar_op, true, true);

    return after_aln; // skip SEQ, QUAL and all AUX fields
}

static VirtualFileOffset bai_increment_vfo (VBlockP vb, VirtualFileOffset vfo, uint64_t inc)
{
    const BgzfBlockPiz *curr_blk = B(BgzfBlockPiz, vb->gz_blocks, vfo.bb_i_in_vb);

    while (inc) {
        ASSERT (vfo.bb_i_in_vb < vb->gz_blocks.len, "%s: went beyond end gz_blocks.len=%u", VB_NAME, vb->gz_blocks.len32);

        uint64_t advance = MIN_(inc, curr_blk->txt_size - vfo.off_in_bb);
        
        vfo.off_in_bb += advance;
        inc -= advance;

        if (vfo.off_in_bb == curr_blk->txt_size) {
            curr_blk++;

            vfo.bb_i_in_vb++;
            vfo.off_in_bb = 0;
        }
    }

    return vfo;
}

static inline BaiBinType get_parent (BaiBinType bin)
{
    return (bin - 1) / 8;
}

// BGZF compute thread: called before BGZF-compressing a BGZF VB 
void bai_calculate_one_vb (VBlockP vb)
{
    START_TIMER;

    // initial allocations
    buf_alloc (vb, &vb->bai_stats,  0, 10,   BaiStats,  0, "bai_stats"); 
    buf_alloc (vb, &vb->bai_linear, 0, 100,  BaiLinear, 0, "bai_linear"); 
    buf_alloc (vb, &vb->bai_chunks, 0, 1000, BaiChunk,  0, "bai_chunks"); 

    // note: since we are not in --bgzf=exact, the BGZF vb is guaranteed to contain only full alignments
    
    rom line = B1STtxt;
    rom after_vb = BAFTtxt;
    
    WordIndex prev_rname=-2, rname=0; // invalid prev_rname, to force new contig in linear index 
    PosType32 prev_pos=0, pos=0;      // zero-based pos
    VirtualFileOffset line_vfo={};
    PosType32 next_win=0; // BAI linear index: A 0-based pos value divisible by 16384

    BaiChunk *bin_tracker[NUM_BINS];   // array of pointers to first chunk of each BIN in current bb_i
    uint64_t bb_i_in_vb = (uint64_t)-1; // current bb_i
    
    vb->line_i=0;
    while (line < after_vb) {
        uint32_t ref_consumed;
        SamFlags sam_flags = {};
        BaiBinType bin; // bin ∈ [0,37449] 
        STR(qname);     // note: in VCF, we used the ID field as a pseudo-qname

        START_TIMER;
        rom next_line = OUT_DT(BAM) ? bai_get_line_bam (line,           &sam_flags, &rname, &pos, &bin, &ref_consumed, pSTRa(qname))
                      : OUT_DT(SAM) ? bai_get_line_sam (line, after_vb, &sam_flags, &rname, &pos, &bin, &ref_consumed, pSTRa(qname))
                      : /* VCF */     bai_get_line_vcf (line, after_vb,             &rname, &pos, &bin, &ref_consumed, pSTRa(qname));
        COPY_TIMER (bai_get_line);

        VirtualFileOffset next_line_vfo = bai_increment_vfo (vb, line_vfo, next_line - line);

        if (vb->line_i == 0) {
            vb->first_rname = rname;
            vb->first_pos   = pos;
        }

        // initialize bb bin tracker in this is a new bgzf block or new rname
        if (line_vfo.bb_i_in_vb != bb_i_in_vb || rname != prev_rname) {
            memset (bin_tracker, 0, sizeof (bin_tracker));
            bb_i_in_vb = line_vfo.bb_i_in_vb;
        }

        // errors only relevant to SAM/BAMs
        if (IS_BAI) { 
            // signs of a corrupt file - qname is also like not valid - don't attempt to print it
            if (pos < -1)                                 bai_abort ("POS out of range", 0, 0, false);
            if (ref_consumed == (uint32_t)-1)/*bad cigar*/bai_abort ("Invalid CIGAR", 0, 0, false);
            if (!IN_RANGE(rname, -1, num_bai_contigs))    bai_abort ("RNAME out of range", 0, 0, false); // comparison as signed captures missing RNAME==-1

            if (rname == -1 && !sam_flags.unmapped)       BAI_ABORT ("An alignment has UNMAPPED flag off, but RNAME is missing");
            if (pos == -1 && !sam_flags.unmapped)         BAI_ABORT ("An alignment has UNMAPPED flag off, but POS is missing");
            if (rname == -1 && pos != -1)                 BAI_ABORT ("An alignment has RNAME, but POS is missing");
            if (pos == -1 && rname != -1)                 BAI_ABORT ("An alignment has POS, but RNAME is missing");
            if (vb->line_i && ( ((uint32_t)rname < (uint32_t)prev_rname) || (pos < prev_pos && rname == prev_rname) ) ) // unsigned comparison bc missing rname (0xffffffff) expected at end of file
                                                          BAI_ABORT ("BAM is not sorted (2)");
        }

        // errors only relevant to VCF
        else {
            // in VCF, rnames are created in CHROM dict in the order in which VBs are merged into zctx in Seg, 
            // which is not necessary the order they appear in the file - so we can test sortness only related to POS, not RNAME
            if (vb->line_i && pos < prev_pos && rname == prev_rname)
                BAI_ABORT ("VCF is not sorted (2)");
        }
        
        if (!IN_RANGE(bin, 0, NUM_BINS-1)) bai_abort ("BIN out of range", 0, 0, false);
        
        // the bb following the last bb overlapping an alignment in this chunk (used for splicing)
        #define next_bb (VirtualFileOffset){ .bb_i_in_vb = next_line_vfo.bb_i_in_vb + (next_line_vfo.off_in_bb > 0) } 
        
        // case: add alignement to existing chunk 
        if (bin_tracker[bin]) {
            bin_tracker[bin]->num_lines++;
            bin_tracker[bin]->after = next_line_vfo;
            bin_tracker[bin]->nxtbb = next_bb; 
        }
        
        // case: bin changed, start new chunk
        else {
            BaiChunk new_chunk = { .rname_lo  = rname & 0xffff, 
                                   .bin       = bin, 
                                   .num_lines = 1, 
                                   .start     = line_vfo, 
                                   .after     = next_line_vfo, 
                                   .nxtbb     = next_bb };
            new_chunk.rname_hi = rname >> 16; // (first 16 bites of .nxtbb) any way to put it in this ^ initialization causes clang warning (gcc is fine)

            buf_append_one (vb->bai_chunks, new_chunk); // careful, rname_lo overlaps
            bin_tracker[bin] = BLST(BaiChunk, vb->bai_chunks);
        }
        
        // if this is a new contig (inc. first line in VB)
        if (rname != prev_rname) {
            next_win = -1; // reset linear index

            // stats for this rname (inc. if rname=-1)            
            buf_append_one (vb->bai_stats, ((BaiStats){ .rname=rname, .ref_beg=line_vfo, .ref_end=next_line_vfo })); 
        }

        // another alignment on same contig - expand range
        else
            BLST(BaiStats, vb->bai_stats)->ref_end = next_line_vfo;

        if (!sam_flags.unmapped) { // always true for VCF
            BLST(BaiStats, vb->bai_stats)->n_mapped++; // stats for this rname
    
            #define POS2WINDOW(pos) ((pos) & ~0x3FFFUL)

            // traverse all 16 Kpbs windows overlapping this alignment and not yet in the linear array
            for (PosType32 win=POS2WINDOW(pos); win <= POS2WINDOW(pos + ref_consumed - 1); win += 16 KB) 
                // case: this alignment overlaps the next 16 Kb window or one even further down
                if (win >= next_win) { // note: >= and not == bc possibly some windows skipped if no alignment overlaps them
                    buf_append_one (vb->bai_linear, ((BaiLinear){ .rname=rname, .window=win, .start=line_vfo }));

                    next_win = win + 16 KB;
                }
        }

        else  // unmapped (never for VCF)
            BLST(BaiStats, vb->bai_stats)->n_unmapped++; // note: this is "unmapped placed" or "unmapped unplaced" depending on whether rname==-1

        if (flag.debug_bai) { // debug goes to info channel, main output to stdout
            iprintf ("%s: rname=%d bin=%u parent_bin=%u pos=%d ref_consumed=%u "VOFF_RNG_FMT,
                     VB_NAME, rname, bin, get_parent(bin), pos, ref_consumed, VOFF_VAL(line_vfo), VOFF_VAL(next_line_vfo));
    
            iprint_newline();
        }

        line        = next_line;
        line_vfo    = next_line_vfo;
        prev_rname = rname;
        prev_pos   = pos;
        vb->line_i++;
    }

    vb->last_rname = rname;
    vb->last_pos   = pos;

    // sanity
    ASSERT (line == after_vb, "%s: expecting next_field==after_vb but (next_field-after_vb)=%d", 
            VB_NAME, (int)(line - after_vb));

    ASSERT (line_vfo.bb_i_in_vb == vb->gz_blocks.len && line_vfo.off_in_bb == 0, 
            "%s: expecting line_vfo="VOFF_FMT" but it is "VOFF_FMT, 
            VB_NAME, vb->gz_blocks.len, 0, VOFF_VAL(line_vfo));

    COPY_TIMER (bai_calculate_one_vb);
}

// Writer thread: after BGZF-compressing one BGZF VB: For all VFOs, change from bb_i_in_vb to bb_off_in_vb
static void bai_update_vfo_to_bgzf_offsets (VBlockP vb)
{
    ARRAY (BgzfBlockPiz, bgzf_blocks, vb->gz_blocks);

    // note: we have one extra entry at the end of bgzf_blocks with the compressed_index of the start of the next VB (see bgzf_compress_vb)
    for_buf (BaiChunk, ent, vb->bai_chunks) {
        ent->start.bb_off_in_vb = bgzf_blocks[ent->start.bb_i_in_vb].compressed_index;
        ent->after.bb_off_in_vb = bgzf_blocks[ent->after.bb_i_in_vb].compressed_index; 
        ent->nxtbb.bb_off_in_vb = bgzf_blocks[ent->nxtbb.bb_i_in_vb].compressed_index;
    }

    for_buf (BaiLinear, ent, vb->bai_linear) 
        ent->start.bb_off_in_vb = bgzf_blocks[ent->start.bb_i_in_vb].compressed_index;

    for_buf (BaiStats, ent, vb->bai_stats) {
        ent->ref_beg.bb_off_in_vb = bgzf_blocks[ent->ref_beg.bb_i_in_vb].compressed_index;
        ent->ref_end.bb_off_in_vb = bgzf_blocks[ent->ref_end.bb_i_in_vb].compressed_index;
    }
}

// Writer thread: after BGZF-compressing one BGZF VB
static void bai_finalize_one_vb_update_stats (VBlockP vb)
{
    // add this VBs alignment stats to txt_file buffer
    ARRAY (BaiStats, stats, txt_file->bai_stats);

    for_buf (BaiStats, vb_stats, vb->bai_stats) 
        if (vb_stats->rname == -1)
            txt_file->bai_stats.count += vb_stats->n_unmapped; // unmapped-unplaced (global)
        
        else {            
            vb_stats->ref_beg.bb_off_in_txt += txt_file->disk_so_far;
            vb_stats->ref_end.bb_off_in_txt += txt_file->disk_so_far;

            // first VB covering this rname
            if (!stats[vb_stats->rname].n_mapped) {
                stats[vb_stats->rname].ref_beg = vb_stats->ref_beg;
                stats[vb_stats->rname].ref_end = vb_stats->ref_end;
            }

            // non-first VB expanding the range of this rname
            else
                stats[vb_stats->rname].ref_end = vb_stats->ref_end;

            stats[vb_stats->rname].n_mapped   += vb_stats->n_mapped;   // mapped (this rname)
            stats[vb_stats->rname].n_unmapped += vb_stats->n_unmapped; // unmapped-placed (this rname)
        }
}

// Writer thread: after BGZF-compressing one BGZF VB
static void bai_finalize_one_vb_merge_chunks (VBlockP vb)
{
    // VFOs: update bb_off_in_vb to bb_off_in_txt 
    for_buf (BaiChunk, chunk, vb->bai_chunks) {
        chunk->start.bb_off_in_txt += txt_file->disk_so_far;
        chunk->after.bb_off_in_txt += txt_file->disk_so_far;
        chunk->nxtbb.bb_off_in_txt += txt_file->disk_so_far;
    }

    // merge first vb chunk into last txt_file->bai_chunks chunk if same bin
    uint64_t start_chunk=0;
    if (txt_file->bai_chunks.len && vb->bai_chunks.len) {
        BaiChunk *prev = BLST(BaiChunk, txt_file->bai_chunks);
        BaiChunk *new  = B1ST(BaiChunk, vb->bai_chunks);

        if (CHUNK_RNAME(prev) == CHUNK_RNAME(new) && prev->bin == new->bin) {
            prev->after = new->after;
            prev->num_lines += new->num_lines;
            start_chunk = 1;
        }
    }

    // append the new chunks (but not the merged one)    
    buf_append (NULL, txt_file->bai_chunks, BaiChunk, B(BaiChunk, vb->bai_chunks, start_chunk), vb->bai_chunks.len - start_chunk, NULL);
}

// Writer thread: after BGZF-compressing one BGZF VB
static void bai_finalize_one_vb_merge_linear (VBlockP vb)
{
    // VFOs: update bb_off_in_vb to bb_off_in_txt 
    for_buf (BaiLinear, ent, vb->bai_linear)
        ent->start.bb_off_in_txt += txt_file->disk_so_far;
    
    BaiLinear *prev = txt_file->bai_linear.len ? BLST(BaiLinear, txt_file->bai_linear) : NULL;
    const BaiLinear *new = vb->bai_linear.len  ? B1ST(BaiLinear, vb->bai_linear) : NULL,
                    *first_new = new, *after_new = BAFT(BaiLinear, vb->bai_linear);

    // skip new windows that are previously recorded (i.e. same or smaller than prev window) 
    uint64_t start_linear=0;
    if (new && txt_file->bai_linear.len) {        
        while (new < after_new && new->rname == prev->rname && new->window <= prev->window) 
            new++;

        start_linear = (new - first_new);
    }

    // append the new linear entries (but not the skipped ones)    
    buf_append (NULL, txt_file->bai_linear, BaiLinear, B(BaiLinear, vb->bai_linear, start_linear), vb->bai_linear.len - start_linear, NULL);
}

static bool bai_verify_sortedness_between_vbs (VBlockP vb)
{
    // verify file is sorted between VBs (bai_calculate_one_vb verified this with each VB)
    static uint32_t last_rname_last_vb=0; // unsigned: unmapped reads with rname=-1 are expected to be last
    static PosType32 last_pos_last_vb=0;

    if (txt_file->bai_chunks.len && // not first VB
        (  (last_rname_last_vb >  vb->first_rname)
        || (last_rname_last_vb != vb->first_rname && IS_TBI) // in VCF, CHROM WordIndex is not necessarily order of appearance in file (rather in order of merger of VBs, unless in VCF header)
        || (last_pos_last_vb > vb->first_pos && last_rname_last_vb == vb->first_rname))) {

        bai_abort ("file is not sorted (3)", 0, 0, false);
        return false;
    }

    last_rname_last_vb = vb->last_rname;
    last_pos_last_vb   = vb->last_pos;
    return true;
} 

// Writer thread: after BGZF-compressing one BGZF VB, but before writing it to disk. Called in order of writing.
// Add VB BAI data to txt_file->bai_*.
void bai_finalize_one_vb (VBlockP vb)
{
    if (!bai_verify_sortedness_between_vbs (vb)) return;

    bai_update_vfo_to_bgzf_offsets (vb);

    bai_finalize_one_vb_update_stats (vb);

    bai_finalize_one_vb_merge_chunks (vb);

    bai_finalize_one_vb_merge_linear (vb);
}

// PIZ main thread: merge two chunks of same (rname,bin) - if the first.after == second.start (e.g. if coming from consecutive VBs)
static void bai_splice_touching_chunks (void)
{
    BaiChunk *next_out = B1ST(BaiChunk, txt_file->bai_chunks), *first = next_out; // pointer aliasing - update in-place
    bool modified = false;

    for_buf (BaiChunk, chunk, txt_file->bai_chunks) {
        if (next_out != first                                                    
        &&  chunk->start.value == (next_out-1)->after.value // touching
        &&  chunk->  bin       == (next_out-1)->bin
        &&  CHUNK_RNAME(chunk) == CHUNK_RNAME(next_out-1)) {

            // extend previous chunk until end of this chunk
            (next_out-1)->after      = chunk->after;
            (next_out-1)->nxtbb      = chunk->nxtbb;
            (next_out-1)->num_lines += chunk->num_lines;
            
            modified = true;
        }

        else {
            if (modified) *next_out = *chunk;
            next_out++;
        }
    }

    txt_file->bai_chunks.len = next_out - first;
}

// PIZ main thread: merge two chunks of same (rname,bin) - in the same or in the consecutive bb
static void bai_splice_chunks (BaiChunkP rname_chunks,        // in/out updates in-place
                               uint64_t *n_chunks_in_rname,   // in/out
                               BaiChunkP *first_chunk_in_bin, // out array
                               uint32_t *n_bins_in_rname,     // out
                               uint32_t *n_chunks_in_bin)     // out array
{
    BaiChunk *first = rname_chunks, *after = rname_chunks + *n_chunks_in_rname, 
             *out = B(BaiChunk, txt_file->bai_chunks, txt_file->bai_chunks.next);
    
    *n_bins_in_rname = 0;

    for (BaiChunkP chunk=first; chunk < after; chunk++) {

        // case: two consecutive chunks with the end or next bb of the first = the start bb of the second
        if (chunk != first 
        &&  (  chunk->start.bb_off_in_txt == (out-1)->after.bb_off_in_txt 
            || chunk->start.bb_off_in_txt == (out-1)->nxtbb.bb_off_in_txt)
        &&  chunk->bin == (out-1)->bin) {

            (out-1)->after      = chunk->after;
            (out-1)->nxtbb      = chunk->nxtbb;
            (out-1)->num_lines += chunk->num_lines;

            (*n_chunks_in_rname)--;
        }

        else {
            n_chunks_in_bin[chunk->bin]++;

            if (n_chunks_in_bin[chunk->bin] == 1) {
                first_chunk_in_bin[chunk->bin] = out;
                (*n_bins_in_rname)++;
            }

            *out++ = *chunk;
        }        
    }

    txt_file->bai_chunks.next = BNUM64 (txt_file->bai_chunks, out);
}

// PIZ main thread
static ASCENDING_ASCENDING_SORTER (sort_chunks_by_bin_and_start, BaiChunk, bin, start.value); // short by bin, then by start.value

// returns total number of bins across all rnames
static uint64_t bai_sort_chunks_by_bin_start (uint64_t *n_chunk_in_rname, BaiChunkP *first_chunk_in_rname)
{
    BaiChunkP chunk = B1ST(BaiChunk, txt_file->bai_chunks), after=BAFT(BaiChunk, txt_file->bai_chunks);
    uint64_t total_bins = 0;

    for (WordIndex rname=0; rname < num_bai_contigs; rname++) {

        *first_chunk_in_rname = chunk;

        // count how many chunks in this rname
        while (chunk < after && CHUNK_RNAME(chunk) == rname) chunk++;

        *n_chunk_in_rname = chunk - *first_chunk_in_rname;

        if (*n_chunk_in_rname) { 
            // sort chunks of this rname by bin value and the start
            qsort (*first_chunk_in_rname, *n_chunk_in_rname, sizeof (BaiChunk), sort_chunks_by_bin_and_start);

            // count bins now while data is warm in cache
            BaiBinType last_bin = BIND_NONE;
            for (chunk=*first_chunk_in_rname; chunk < after && CHUNK_RNAME(chunk) == rname; chunk++)
                if (chunk->bin != last_bin) {
                    total_bins++;
                    last_bin = chunk->bin;
                }
        }

        n_chunk_in_rname++;
        first_chunk_in_rname++;
    }

    return total_bins;
}

// PIZ main thread
static char *bai_write_next_bin (WordIndex rname, uint32_t bin, char *out, const BaiChunk *chunk, uint32_t n_chunks, char *after_txt)
{
    const BaiChunk *first=chunk, *after=chunk + n_chunks;

    ASSERT0 (out + 2 * sizeof(uint32_t) <= after_txt, "txt_data overflow");

    PUT_UINT32 (out, chunk->bin); 
    char *n_chunks_p = out + sizeof(uint32_t); 
    out += 2 * sizeof(uint32_t);

    while (chunk < after) {
        ASSERT (chunk->bin == bin && CHUNK_RNAME(chunk) == rname, "chunk_i=%u: expecting (rname,bin)=(%u,%u) but seeing (%u,%u)", 
                (int)(chunk - first), rname, bin, CHUNK_RNAME(chunk), chunk->bin);
                
        if (flag.debug_bai) 
            iprintf ("writing chunk: rname=%u bin=%u chunk="VOFF_RNG_FMT"\n", rname, bin, VOFF_RNG_VAL(*chunk));
        
        ASSERT0 (out + sizeof (VirtualFileRange) <= after_txt, "txt_data overflow");
        
        out = mempcpy (out, &chunk->start, sizeof (VirtualFileRange)); // TO DO: Endianity
        chunk++;
    }
    
    PUT_UINT32 (n_chunks_p, n_chunks);
    
    return out;
}

// PIZ main thread
static char *bai_write_pseudo_bin (WordIndex rname, char *out, char *after)
{
    ASSERT0 (out + 40 <= after, "txt_data overflow");

    const BaiStats *stats = B(BaiStats, txt_file->bai_stats, rname);

    PUT_UINT32 (out, BIN_STATS);            out += sizeof (uint32_t);
    PUT_UINT32 (out, 2);                    out += sizeof (uint32_t); // number of "chunks" in this pseudo bin
    PUT_UINT64 (out, stats->ref_beg.value); out += sizeof (uint64_t);
    PUT_UINT64 (out, stats->ref_end.value); out += sizeof (uint64_t);
    PUT_UINT64 (out, stats->n_mapped);      out += sizeof (uint64_t);
    PUT_UINT64 (out, stats->n_unmapped);    out += sizeof (uint64_t);
        
    return out;
}

static uint64_t bai_count_intervals (void)
{
    uint64_t count=0;
    const BaiLinear *last = BLST(BaiLinear, txt_file->bai_linear);

    for_buf (const BaiLinear, lin, txt_file->bai_linear) {
        ASSERT (lin->window % (16 KB) == 0 && lin->window < 2 GB, "invalid window=%u (rname=%d)", lin->window, lin->rname); // sanity check while we're at it

        // for each rname, at intervals up to an including last window of this rname
        if (lin == last || lin->rname != (lin+1)->rname)
            count += 1 + (lin->window / (16 KB));
    }
    return count;
}

// PIZ main thread: Generate intervals from BaiLinear data for one rname
static uint32_t bai_get_intervals_one_rname (WordIndex rname, VirtualFileOffset *intervals)
{
    const BaiLinear *first = B(BaiLinear, txt_file->bai_linear, txt_file->bai_linear.next), 
                    *next, *after = BAFT(BaiLinear, txt_file->bai_linear); // 16-byte units

    uint32_t n_intervals = 0; // counting intervals of this rname

    if (flag.show_bai == SHOW_BAI_LINEAR)
        printf ("Linear 🧬=%u\n", rname); // note: main output goes to stdout, while debug, show-time etc go to info_channel 

    for (next=first; next < after && next->rname == rname; next++) 

        // make linear entry, and the same entry for all previous windows that had no alignments
        while (n_intervals <= next->window / (16 KB)) {
            intervals[n_intervals++] = next->start;
            
            if (flag.show_bai == SHOW_BAI_LINEAR)
                printf ("🧬=%u.%u ➰="VOFF_FMT"\n", rname, n_intervals, VOFF_VAL(next->start));
        }            

    txt_file->bai_linear.next += next - first;

    return n_intervals;
}

static ASCENDING_SORTER (chunk_sort_by_start, BaiChunk, start.value)
static void bai_show_chunks (void)
{
    FILE *stream = IS_PIZ ? info_stream // genounzip / genocat --show-bai
                 :          stdout;     // genozip --show-bai

    if (is_genounzip) { // in genounzip, make a copy to not ruin the original
        ASSERTNOTINUSE (evb->scratch);
        buf_copy (evb, &evb->scratch, &txt_file->bai_chunks, BaiChunk, 0, 0, "scratch");
    }

    BufferP chunks_buf = is_genozip ? &evb->bai_chunks : &evb->scratch;

    // sort chunks by start of range (assuming: also stays sorted by rname, given the file is sorted)
    qsort (B1ST(BaiChunk, *chunks_buf), chunks_buf->len, sizeof (BaiChunk), chunk_sort_by_start);

    uint64_t highest_after = B1ST(BaiChunk, *chunks_buf)->start.value;

    fprintf (stream, "🧬=contig 🗑 =bin\tⓅ=parent 📑=n_lines   ↕=[VFO₁ — VFO₂)\n");

    for_buf (BaiChunk, chunk, *chunks_buf) {
        bool overlap = chunk->start.value != highest_after;

        fprintf (stream, "🧬=%u\t🗑 =%-5u\tⓅ=%-4u%s\t" VOFF_RNG_FMT "%s", 
                 CHUNK_RNAME(chunk), chunk->bin, get_parent (chunk->bin), 
                 cond_int (is_genounzip, "\t📑=", chunk->num_lines),
                 VOFF_RNG_VAL(*chunk), (overlap ? "    \tOVERLAPPING\n" : "\n"));

        MAXIMIZE (highest_after, chunk->after.value);
    }

    buf_free (evb->scratch);
}

// in VCF, the CHROM word_indices are according to context merge order (unless in the VCF header).
// now we find a mapping between CHROM word_index and RNAME (the order they appear in the file = order in reference file)
static void bai_map_rname_wi (void)
{
    #define wi_to_rname_buf evb->codec_bufs[5]
    #define rname_to_wi_buf evb->codec_bufs[6]

    ARRAY_alloc (WordIndex, wi_to_rname, ZCTX(VCF_CHROM)->word_list.len, false, wi_to_rname_buf, evb, "evb->codec_bufs[5]");
    ARRAY_alloc (WordIndex, rname_to_wi, ZCTX(VCF_CHROM)->word_list.len, false, rname_to_wi_buf, evb, "evb->codec_bufs[6]");
    memset (wi_to_rname, 0xff, wi_to_rname_len * sizeof (WordIndex));
    memset (rname_to_wi, 0xff, rname_to_wi_len * sizeof (WordIndex));

    // bai_linear is in the order of the txt file. use it to assign rname (in order of appearance in file) to word_index
    WordIndex next_rname = 0, last_wi = WORD_INDEX_NONE;
    for_buf (BaiLinear, ent, txt_file->bai_linear) {
        // first encounter with a new rname
        if (ent->rname != last_wi) { // rname in bai_linear are still actually word_index
            wi_to_rname[ent->rname] = next_rname;
            rname_to_wi[next_rname] = ent->rname;
            last_wi = ent->rname;

            next_rname++;
        }

        ent->rname = wi_to_rname[ent->rname]; // update - now it is truely the rname (not CHROM word_index)
    }
    rname_to_wi_buf.len = num_bai_contigs = next_rname; // update to actual number of CHROMs in this VCF file

    // update rname in chunks
    for_buf (BaiChunk, ent, txt_file->bai_chunks) {
        WordIndex wi = CHUNK_RNAME(ent);
        ent->rname_lo = wi_to_rname[wi] & 0xffff;
        ent->rname_hi = wi_to_rname[wi] >> 16;
    }

    // re-sort stats by rname instead of wi
    ASSERTNOTINUSE (evb->scratch);
    buf_copy (evb, &evb->scratch, &txt_file->bai_stats, BaiStats, 0, 0, "scratch");

    txt_file->bai_stats.len = num_bai_contigs;
    for_buf2 (BaiStats, dst, rname, txt_file->bai_stats)
        *dst = *B(BaiStats, evb->scratch, rname_to_wi[rname]);

    buf_free (evb->scratch);
}

static void bai_prepare_tbi_header (void)
{
    ARRAY (const CtxWord, words, ZCTX(VCF_CHROM)->word_list);
    ARRAY (const char,    dict,  ZCTX(VCF_CHROM)->dict);
    ARRAY (const WordIndex, rname_to_wi, evb->codec_bufs[6]);

    // find num_bai_contigs and names_len: all CHROM snips, until the last wi which appears in the file
    uint32_t names_len=0;
    for (WordIndex rname=0; rname < num_bai_contigs; rname++) 
        names_len += words[rname_to_wi[rname]].snip_len + 1; // +1 for \0s
  
    ASSERTNOTINUSE (evb->txt_data);
    ARRAY_alloc (char, tbi_header, sizeof (TbiHeader) + names_len, false, evb->txt_data, evb, "txt_data");
    *(TbiHeader *)tbi_header = (TbiHeader){
        .magic     = TBI_MAGIC,
        .n_ref     = LTEN32 (num_bai_contigs),
        .format    = LTEN32 (2),
        .col_rname = LTEN32 (1),
        .col_pos   = LTEN32 (2),
        .meta      = LTEN32 ('#'),
        .names_len = LTEN32 (names_len)
    };
    
    char *next = ((TbiHeader *)tbi_header)->names, *first=next;
    for (WordIndex rname=0; rname < num_bai_contigs; rname++) {
        CtxWord word = words[rname_to_wi[rname]];
        next = mempcpy (next, &dict[word.char_index], word.snip_len + 1); // including \0
    }

    // sanity
    ASSERT (next - first == names_len, "names miscalculation: names.len=%u ≠ names_len=%u", (int)(next-first), names_len);
    ASSERT (BNUM(evb->txt_data, next) == evb->txt_data.len, "wrong tbi header len: actual=%u ≠ expected=%u", BNUM(evb->txt_data, next), evb->txt_data.len32);
}

// PIZ main thread
void bai_write (void)
{
    START_TIMER;
    VBlockP vb = evb; // so we can use RECON macros
    vb->data_type = DT_GNRIC; // to avoid error in piz_dis_qname that can be called by RECONSTRUCT->ASSPIZ

    char filename[strlen(txt_file->name) + 5];
    sprintf (filename, "%s.%s", txt_file->name, IS_BAI ? "bai" : "tbi");

    file_remove (filename, true);

    FILE *fp = fopen (filename, "wb");
    ASSRET (fp,, "WARNING: Failed to write %s. fopen: %s", filename, strerror (errno));

    ARRAY_alloc (uint64_t, n_chunks_in_rname,      num_bai_contigs,  true,  evb->codec_bufs[0], evb, "codec_bufs[0]");
    ARRAY_alloc (BaiChunk *, first_chunk_in_rname, num_bai_contigs,  true,  evb->codec_bufs[1], evb, "codec_bufs[1]");
    ARRAY_alloc (VirtualFileOffset, intervals,     (2 GB) / (16 KB), true,  evb->codec_bufs[2], evb, "codec_bufs[2]"); // max POS is 2 GB per SAM spec, so this is the maximum theoritcal intervals per rname
    ARRAY_alloc (BaiChunkP, first_chunk_in_bin,    NUM_BINS,         false, evb->codec_bufs[3], evb, "codec_bufs[3]");
    ARRAY_alloc (uint32_t, n_chunks_in_bin,        NUM_BINS,         false, evb->codec_bufs[4], evb, "codec_bufs[4]");

    // Show raw: chunks by order of VFO, before splicing
    if (flag.show_bai == SHOW_BAI_RAW) {
        bai_splice_touching_chunks(); // minimal splicing - just remove artifacts of VB boundaries
        bai_show_chunks();
    }

    if (IS_TBI) { // write TBI header into evb->txt_data
        bai_map_rname_wi();
        bai_prepare_tbi_header(); 
    }

    // after this, chunks are sorted by rname, then bin, then start. 
    uint64_t total_bins = bai_sort_chunks_by_bin_start (n_chunks_in_rname, first_chunk_in_rname);

    uint64_t size_upperbound = (IS_BAI ? 8/* magic, n_ref */ : evb->txt_data.len)    // header
                             + num_bai_contigs * (sizeof (uint32_t) * 2 + 40)        // n_bins, n_intervals and 40-byte stats pseudo-bin for each contig 
                             + total_bins * sizeof (uint32_t) * 2                    // bin, n_chunks
                             + txt_file->bai_chunks.len * sizeof (VirtualFileRange)  // chunks (possibly less, after splicing)
                             + bai_count_intervals() * sizeof (VirtualFileOffset)    // intervals                             
                             + sizeof (uint64_t);                                    // n_no_coord

    buf_alloc (evb, &evb->txt_data, 0, size_upperbound, char, 0, "txt_data");
 
    if (IS_BAI) { // write BAI header
        RECONSTRUCT(BAI_MAGIC, 4);
        RECONSTRUCT_BIN32 (num_bai_contigs);
    }

    for (WordIndex rname=0; rname < num_bai_contigs; rname++) {

        if (n_chunks_in_rname[rname]) {
            memset (first_chunk_in_bin, 0, NUM_BINS * sizeof (BaiChunkP)); // initialize start & n_chunks of bins of this rname
            memset (n_chunks_in_bin,    0, NUM_BINS * sizeof (uint32_t));

            // splice same-bin chunks in same- or adjacent bb's (even if the bins are not adjacent) 
            uint32_t n_bins_in_rname;
            bai_splice_chunks (first_chunk_in_rname[rname], &n_chunks_in_rname[rname], first_chunk_in_bin, &n_bins_in_rname, n_chunks_in_bin); 

            // note: unlike samtools, we don't merge non-same-bin chunks: samtools merges some small chunks into  adjacent chunks of their parent bin. 
            // samtools's algorithm for deciding which to merge is neither documented nor obvious from studying the BAI files, and has neglibile benefit in terms of BAI file size.

            char *after_txt = BAFTtxt + sizeof (uint32_t)/*n_bins*/ + n_bins_in_rname * 8/*bin headers*/ + n_chunks_in_rname[rname] * sizeof (VirtualFileRange) + 40/*stats "bin"*/;

            char *next = BAFTtxt + sizeof(uint32_t); // leave room for n_bins
            
            for (BaiBinType bin=0; bin < NUM_BINS; bin++) 
                if (n_chunks_in_bin[bin])
                    next = bai_write_next_bin (rname, bin, next, first_chunk_in_bin[bin], n_chunks_in_bin[bin], after_txt);

            next = bai_write_pseudo_bin (rname, next, after_txt);

            ASSERT (next == after_txt, "unexpected amount of chunk data: remaining=%d\n", (int)(BAFTc(evb->txt_data) - next));

            RECONSTRUCT_BIN32 (n_bins_in_rname ? (1 + n_bins_in_rname) : 0); // 1 extra bin for stats

            Ltxt = BNUMtxt (after_txt);
        }

        else  // empty rname 
            RECONSTRUCT_BIN32 (0);
        
        uint32_t n_intervals = bai_get_intervals_one_rname (rname, intervals); // modifies txt_file->bai_linear in-place
        RECONSTRUCT_BIN32 (n_intervals);

        if (n_intervals) RECONSTRUCT (intervals, n_intervals * sizeof (VirtualFileOffset));
    }

    RECONSTRUCT_BIN64 (IS_BAI ? txt_file->bai_stats.count : 0); // n_no_coord: unmapped unplaced (global) (always 0 for VCF)

    txt_file->bai_chunks.len = txt_file->bai_chunks.next; // update after all the splicing

    // Show chunks: final BAI, after splicing
    if (flag.show_bai == SHOW_BAI_CHUNKS) 
        bai_show_chunks();
    
    codec_free_all (evb); // free evb->codec_bufs

    ASSERT (txt_file->bai_linear.next == txt_file->bai_linear.len, 
            "txt_file->bai_linear not perfectly consumed: len=%u next=%u", txt_file->bai_linear.len32, (int)txt_file->bai_linear.next);

    if (IS_BAI)
        ASSGOTO (fwrite (vb->txt_data.data, vb->txt_data.len, 1, fp) == 1, 
                 "WARNING: Failed to write %s. fwrite: %s", filename, strerror (errno));

    // TBI needs to be BGZF-compressed
    else {
        bgzf_compress_tbi();

        ASSGOTO (fwrite (vb->comp_txt_data.data, vb->comp_txt_data.len, 1, fp) == 1, 
                 "WARNING: Failed to write %s. fwrite: %s", filename, strerror (errno));
        
        buf_destroy (vb->comp_txt_data);
    }
    
    ASSGOTO (!fclose (fp), "WARNING: Failed to write %s. fclose: %s", filename, strerror (errno));
    
    buf_free (vb->txt_data);
    flag.make_bai = false; 

    COPY_TIMER_EVB (bai_write);
    return;

error:
    file_remove (filename, true);
}

//-----------------------------------------------------
// code for genozip --show-bai/tbi file.bai or file.tbi
//-----------------------------------------------------

// structure of one bin, per SAM spec §5.2
typedef struct { 
    uint32_t bin, n_chunks; 
    union {
        // normal bin
        VirtualFileRange chunks[0];

        struct { // stats pseudo-bin
            VirtualFileRange contig;
            uint64_t n_mapped, n_unmapped;
        };
    };
} BinStruct, *BinStructP, __attribute__((aligned(1))) *UnalignedBinStructP;

static void bai_show_one_bin (BinStructP bin, rom rname_s)
{
    // optional pseudo-bin: documented for BAI, undocumented for TBI, yet appears in tabix-produced tbi files
    if (bin->bin == BIN_STATS) {
        ASSINP0 (bin->n_chunks == 2, "expecting pseudo-bin n_chunks==2");

        printf ("🧬=%s\t🗑 =%u x 1\t%s=%"PRIu64"%s  "VOFF_RNG_FMT,
                rname_s, bin->bin, 
                IS_BAI ? "✅✅✅" : "⌬⌬⌬", bin->n_mapped, 
                cond_int (IS_BAI, " ❌❌❌=", bin->n_unmapped), 
                VOFF_RNG_VAL(bin->contig));
        DO_ONCE if (!flag.quiet) 
            printf ("   ⇐ %s", 
                    IS_BAI ? "reads placed on this contig: mapped, unmapped, range" : "variants in this contig: number and range");

        printf ("\n");
    }

    else for (uint32_t chunk_i=0; chunk_i < bin->n_chunks; chunk_i++) {
        printf ("🧬=%s\t🗑 =%u x %u\tⓅ=%-4u\t" VOFF_RNG_FMT,  
                rname_s, bin->bin, bin->n_chunks, get_parent (bin->bin), VOFF_RNG_VAL(bin->chunks[chunk_i]));

        DO_ONCE {      if (!flag.quiet) printf ("   ⇐ a chunk: a [VFO₁ — VFO₂) range between two V̲irtual F̲ile O̲ffsets"); }
        else { DO_ONCE if (!flag.quiet) printf ("   ⇐ VFO ≐ BGZF_blk-offset-in-file⁀offset-in-uncomp-BGZF_blk"); }
        printf ("\n");
    }
}

// read BAI or TBI file data into evb->txt_data
static void bai_show_read_file (rom filename)
{
    ASSINP (file_exists (filename), "Error: file not found: %s", filename);
    
    uint64_t file_size = file_get_size (filename);
    ASSINP (file_size, "Error: file is empty: %s", filename);
        
    txt_file = file_open_txt_read (filename);
    ASSERTNOTNULL (txt_file);

    ASSINP (txt_file->type == BAI || txt_file->type == TBI, "%s: Expecting a file with a .bai or .tbi extension", filename);
    index_type = txt_file->type;

    // read the entire file (and uncompress it if needed) into evb->txt_data. note: TBI is BGZF-compressed and BAI is not.
    segconf.vb_size = MAX_(MIN_VBLOCK_MEMORY MB, MIN_(MAX_VBLOCK_MEMORY MB, (IS_TBI ? 8 : 1.1) * file_size)); // should be large enough for reasonable BGZF-compressed TBI files
    txtfile_read_vblock (evb); 

    ASSINP (txt_file->no_more_blocks, "%s file %s size=%"PRIu64" is beyond genozip's limits", 
            IS_BAI?"BAI":"TBI", filename, file_size);

    file_close (&txt_file);
}

// execute command --show-bai
static ASCENDING_SORTER_INDIRECT (bin_sorter, BinStruct, bin)
void bai_show (rom filename)
{
    START_TIMER;

    #define GET_BAI32 ({ bai += sizeof (uint32_t); GET_UINT32(bai-4); })
    #define GET_BAI64 ({ bai += sizeof (uint64_t); GET_UINT64(bai-8); })
    
    profiler_new_z_file(); // initialize profiler

    bai_show_read_file (filename);

    rom bai = B1STc(evb->txt_data);

    if ((IS_BAI && memcmp (bai, BAI_MAGIC, STRLEN(BAI_MAGIC)))
    ||  (IS_TBI && memcmp (bai, TBI_MAGIC, STRLEN(TBI_MAGIC))))
        ABORTINP ("%s: bad magic: not a %s file", filename, IS_BAI?"BAI":"TBI");

    // handle TBI header
    if (IS_TBI) {
        TbiHeader *h = (TbiHeader *)bai;
        
        ASSINP (h->format==2 && h->col_rname==1 && h->col_pos==2 && h->col_end==0 && h->meta=='#' && h->skip==0,
                "Not VCF TBI file: format=%u col_rname=%u col_pos=%u col_end=%u meta=%c skip=%u", 
                h->format, h->col_rname, h->col_pos, h->col_end, h->meta, h->skip);

        num_bai_contigs = h->n_ref;
        buf_alloc_exact_zero (evb, evb->tbi_contigs, num_bai_contigs, rom, "tbi_contigs"); // not on stack, because there could be lots on contigs

        bai = h->names;
        for_buf (rom, contig_name, evb->tbi_contigs) {
            *contig_name = bai;
            bai += strlen (bai) + 1; // past name and \0
        }   
    } 

    // handle BAI header
    else { 
        bai += 4; // past magic
        num_bai_contigs = GET_BAI32;
    }

    bool show_chunks = (flag.show_bai == SHOW_BAI_CHUNKS);
    if (!show_chunks) {
        printf ("🧬🧬🧬=%u", num_bai_contigs);

        DO_ONCE if (!flag.quiet) 
            printf ("   ⇐ number of contigs in %s (--quiet to remove comments)", 
                    IS_BAI ? "BAM header" : "VCF file");
        printf ("\n");

    }
    else 
        buf_alloc (evb, &evb->bai_chunks, 0, evb->txt_data.len / sizeof (VirtualFileRange), BaiChunk, 0, "bai_chunks"); // theoretical max size

    #define rname_str (IS_BAI ? str_int_s(rname).s : *B(rom, evb->tbi_contigs, rname))

    for (uint32_t rname=0; rname < num_bai_contigs; rname++) {
        uint32_t n_bins = GET_BAI32;

        if (!show_chunks) {
            printf ("R-Tree 🧬=%s\t🗑🗑🗑 =%u", rname_str, n_bins);
            { DO_ONCE if (!flag.quiet) printf ("   ⇐ number of bins in this contig"); }
            printf ("\n");
        }

        // get all all bins of this rname, moving bai to after this rname data
        UnalignedBinStructP bins[n_bins];
        for (int bin_i=0; bin_i < n_bins; bin_i++) {
            bins[bin_i] = (BinStructP)bai; 
            bai += 2 * sizeof(uint32_t) + 2 * sizeof(uint64_t) * bins[bin_i]->n_chunks; // skip to after bin
        }

        // sort bins of this rname by bin value (technically: sorting array of pointers to the start of each bin's data)
        if (flag.show_bai == SHOW_BAI_SORT)
            qsort (bins, n_bins, sizeof (BinStructP), bin_sorter);

        if (!show_chunks) 
            for (int bin_i=0; bin_i < n_bins; bin_i++) 
                bai_show_one_bin (bins[bin_i], rname_str);

        else for (int bin_i=0; bin_i < n_bins; bin_i++) { // show_chunks
            UnalignedBinStructP bin = bins[bin_i];

            if (bin->bin == BIN_STATS) continue;

            for (uint32_t chunk_i=0; chunk_i < bin->n_chunks; chunk_i++) 
                BNXT(BaiChunk, evb->bai_chunks) = (BaiChunk){
                    .rname_lo = rname & 0xffff,
                    .rname_hi = rname >> 16,
                    .bin      = bin->bin,
                    .start    = bin->chunks[chunk_i].start,
                    .after    = bin->chunks[chunk_i].after     
                };
        }

        uint32_t n_intervals = GET_BAI32; // number of linear intervals
        
        if (show_chunks)
            bai += n_intervals * sizeof(uint64_t); // move to after intervals

        else {
            printf ("Linear 🧬=%s ➰➰➰=%u", rname_str, n_intervals);
            DO_ONCE if (!flag.quiet) printf ("   ⇐ number of linear index intervals in this contig");
            printf ("\n");

            for (uint32_t interval_i=0; interval_i < n_intervals; interval_i++) {
                VirtualFileOffset ioffset = { .value = GET_BAI64 };  
                
                printf ("🧬=%s.%u ➰="VOFF_FMT, rname_str, interval_i, VOFF_VAL(ioffset));
                DO_ONCE if (!flag.quiet) printf ("   ⇐ VFO of first alignment in this interval");
                printf ("\n");
            }
        }
    } // for rname
    
    if (show_chunks)
        bai_show_chunks();
    
    if (bai == BAFTc(evb->txt_data) - sizeof(uint64_t)) {
        uint64_t n_unplaced_unmapped = GET_BAI64; // optional
    
        if (!show_chunks && IS_BAI) { // n_unplaced_unmapped might exist in TBI, but is not used for VCF
            printf ("❌❌❌=%"PRIu64, n_unplaced_unmapped);
            { DO_ONCE if (!flag.quiet) printf ("   ⇐ number of unplaced-unmapped reads in the BAM file"); }
            printf ("\n");
        }
    }

    ASSINP (bai == BAFTc(evb->txt_data), "Unexpected end of bai: after-bai=%d words", (int)(BAFTc(evb->txt_data) - bai));

    COPY_TIMER_EVB(show_bai);

    if (flag.show_time)
        profiler_add_evb_and_print_report();
}
