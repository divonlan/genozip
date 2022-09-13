// ------------------------------------------------------------------
//   codec_pbwt.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

// The codec implements a modified version of the PBWT algorithm. The permutation logic is loosely based on the logic in
// Durbin R. Efficient haplotype matching and storage using the positional Burrows-Wheeler transform (PBWT). Bioinformatics. 2014 

#include "genozip.h"
#include "codec.h"
#include "buffer.h"
#include "file.h"
#include "vblock.h"
#include "dict_id.h"
#include "reconstruct.h"
#include "strings.h"
#include "compressor.h"
#include "profiler.h"
#include "context.h"
#include "endianness.h"
#include "piz.h"

typedef struct {
    Allele allele;      // copied from line on which permutation is applied
    uint32_t index;     // index into ht_one_line from which this allele was copied
} PermEnt;

typedef struct {
    BufferP runs;       // array of uint32_t - alternating bg ('0') / fg runs
    BufferP fgrc;       // array PbwtFgRunCount describing the fg runs
    Allele run_allele;  // current allele for which we're constructing a run 
    PermEnt *perm;      // a permutation - expressed as indices into the haplotype line 
    PermEnt *temp;      // working memory
} PbwtState;

// this struct is part of the file format
typedef struct __attribute__ ((__packed__)) {
    uint32_t fg_allele : 8;
    uint32_t count     : 24;   // number of consecutive runs in PbwtState->runs with fg_allele 
} PbwtFgRunCount;

#define htputc(c) fprintf (info_stream, IS_NON_WS_PRINTABLE(c) ? "%c" : "\\x%02x", (c))

void codec_pbwt_display_ht_matrix (VBlockP vb, uint32_t max_rows)
{
    ARRAY (Allele, ht, vb->ht_matrix_ctx->local);
    uint32_t cols = vb->ht_per_line;
    uint32_t rows = ht_len / vb->ht_per_line;
    max_rows = max_rows ? MIN_(rows, max_rows) : rows;

    iprint0 ("\nht_matrix:\n");
    for (uint32_t r=0; r < max_rows; r++) {
        for (uint32_t c=0; c < cols; c++) 
            iputc (ht[r * cols + c]);
        iprint0 ("\n");
    }
} 

static void show_runs (const PbwtState *state)
{
    fprintf (info_stream,  "RUNS: ");
    
    PbwtFgRunCount rc = {};
    PbwtFgRunCount *next_allele = B1ST (PbwtFgRunCount, *state->fgrc);

    for (unsigned i=0, fg=0; i < state->runs->len; i++, fg=!fg) {

        if (fg && !rc.count) {
            ASSERT0 (next_allele < BAFT (PbwtFgRunCount, *state->fgrc), "premature end of alleles buffer");
            rc = *next_allele++;
        }
    
        iputc ('\''); 
        htputc (fg ? rc.fg_allele : '0');
        fprintf (info_stream, "':%u ", *B32 (*state->runs, i));

        if (fg) rc.count--;
    }
    iprint0 ("\n");
}

#define SHOW(msg, cmd) do { \
    fprintf (info_stream, msg " %-2u: ", line_i); \
    for (uint32_t i=0; i < vb->ht_per_line; i++) cmd; \
    iprint0 ("\n"); } while (0) // flush

#define show_line    if (flag.show_alleles) SHOW ("LINE", (htputc (*B(Allele, vb->ht_matrix_ctx->local, line_i * vb->ht_per_line + i))))
#define show_perm(s) if (flag.show_alleles) SHOW ("PERM", (fprintf (info_stream, "%d ", (s)->perm[i].index)));       

static PbwtState codec_pbwt_initialize_state (VBlockP vb, BufferP runs, BufferP fgrc)
{
    buf_alloc (vb, &vb->codec_bufs[0], 0, vb->ht_per_line * 2, PermEnt, 1, "codec_bufs");
    buf_zero (&vb->codec_bufs[0]); // re-zero every time
    ARRAY (PermEnt, state_data, vb->codec_bufs[0]);

    PbwtState state = {
        .runs = runs,
        .fgrc = fgrc,
        .perm = &state_data[0],               // size: ht_per_line X PermEnt
        .temp = &state_data[vb->ht_per_line], // size: ht_per_line X PermEnt
    };
        
    return state;
} 

// update the permutation for the next row: we re-sort it to make indices containing the same allele grouped
// first '0', then '1' etc - but the order within each of these allele groups remains as in the current row's premutation 
// (this is why we traverse the permuted line rather than the ht_matrix line)
static void inline codec_pbwt_calculate_permutation (PbwtState *state, const Allele *line, uint32_t line_len, bool is_first_line)
{
    // populate permutation index - by re-ordering according to previous line's alleles
    if (!is_first_line) {

        // all valid allele values - in the order we want them in the permutation
        #define NUM_VALIDS 104
        static const Allele valids[NUM_VALIDS] = { 
            48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,
            75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,
            106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,
            129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,
            '.', '*', '%', '-'
        };
        
        // mark which allele values exist in prev line
        bool has_allele[256] = {};
        for (uint32_t ht_i=0; ht_i < line_len; ht_i++) 
            has_allele[state->perm[ht_i].allele] = true;

        // re-order permutation - first taking the indices of the '0' alleles, then '1', '2' etc - but keeping
        // the order within each allele as it was 
        uint32_t temp_i=0;        
        for (Allele i=0; i < NUM_VALIDS; i++) {
            Allele bg = valids[i];
        
            if (has_allele[bg]) 
                for (uint32_t ht_i=0; ht_i < line_len; ht_i++) 
                    if (state->perm[ht_i].allele == bg) 
                        state->temp[temp_i++].index = state->perm[ht_i].index;
        }

        SWAP (state->perm, state->temp);
    }

    else // first line - initialize permutation to identity
        for (uint32_t ht_i=0; ht_i < line_len; ht_i++) 
            state->perm[ht_i].index = ht_i;

    // ZIP: populate alleles
    if (IS_ZIP) 
        for (uint32_t ht_i=0; ht_i < line_len; ht_i++)
            state->perm[ht_i].allele = line[state->perm[ht_i].index];
}

// -------------
// ZIP side
// -------------

// called by vcf_seg_finalize to create all contexts - must be done before merge
void codec_pbwt_seg_init (VBlockP vb, Context *runs_ctx, Context *fgrc_ctx)
{
    vb->runs_ctx = runs_ctx;
    vb->fgrc_ctx = fgrc_ctx;
    
    vb->ht_matrix_ctx->lcodec = CODEC_PBWT; // this will trigger codec_pbwt_compress even though the section is not written to the file

    // create a contexts (if one doesn't exist already) for all alleles observed in this VB. 
    // Note: we can't do this in codec_pbwt_compress because it must be done before merge, do be copied to z_file->contexts

    // this context will contain alternate run lengths of the background allele and any forward allele
    runs_ctx->ltype         = LT_UINT32;
    runs_ctx->local_dep     = DEP_L1; // RUNS.local is generated when PBWT_HT_MATRIX.local is compressed
    
    // this context is used to determine which forward allele is in each forward run in RUNS: it contains
    // an array of PbwtFgRunCount
    fgrc_ctx->ltype         = LT_UINT32;
    fgrc_ctx->lsubcodec_piz = CODEC_PBWT;
    fgrc_ctx->local_dep     = DEP_L2; // FGRC.local is generated when PBWT_HT_MATRIX.local is compressed, and *must* be after RUNS in z_file
}

// updates FGRC - our list of foreground run alleles. Each entry represents "count" consecutive runs of this same "fg_allele"
// returns true if a background run needs to be added to account for a missing background run between two foreground runs of different alleles
static bool inline codec_pbwt_udpate_fgrc (PbwtState *state, Allele allele_of_completed_run)
{
    // case: new run a foreground allele as previous run was a background - update foreground run count
    if (allele_of_completed_run == '0') {
        PbwtFgRunCount *rc;

        // case: previous foreground run was same allele as this - increment the count
        if (  state->fgrc->len && // we already have foreground runs
                state->run_allele == (rc = BLST (PbwtFgRunCount, *state->fgrc))->fg_allele) { 
            rc->count++; // one more consecutive foreground run of this allele
            ASSERT (rc->count, "Number of consecutive forward runs with allele=%c exceeded the maximum of %u",
                    state->run_allele, 0xffffff); // error is count round-robined back to 0
        }
        // case: the previous fg run (if any) was of a different allele - start RunCount entry
        else
            BNXT (PbwtFgRunCount, *state->fgrc) = (PbwtFgRunCount){ .fg_allele = state->run_allele, .count = 1 };

        return false; // no need to add a background run
    }

    // case: this run is fg, and prev run was fg too (or this is the first run: run_allele==0) and too - skip the background run
    else if (state->run_allele != '0') {
        BNXT (PbwtFgRunCount, *state->fgrc) = (PbwtFgRunCount){ .fg_allele = state->run_allele, .count = 1 }; // new foreground
        return true; // add a 0 background run to account for the missing background run between two different foreground alleles
    }

    // case: run completed is a foreground and new run is a background - no need for a FGRC entry for background runs
    else 
        return false; // no need for an extra empty background run
}

// add data to the run-length encoding (may be called with the full ht matrix or just pieces of it)
static void inline codec_pbwt_run_len_encode (PbwtState *state, uint32_t line_len, bool backwards) 
{ 
    for (uint32_t ht_i=0; ht_i < line_len; ) {		

        #define ALLELE(ht_i) state->perm[backwards ? line_len - (ht_i) -1 : (ht_i)].allele

        uint32_t run_len=0;
        for (; ht_i < line_len && ALLELE(ht_i) == state->run_allele; ht_i++) 
            run_len++;

        // case: we have more of the existing run - extend the current run length
        if (run_len) 
            *BLST32 (*state->runs) += run_len;
        
        // case: we need to start new run 
        else {
            Allele allele_of_completed_run = state->run_allele;
            state->run_allele = ALLELE(ht_i); // new run allele            

            if (codec_pbwt_udpate_fgrc (state, allele_of_completed_run)) 
                BNXT32 (*state->runs) = 0; // entry for missing background run (length 0) between runs of 2 different foreground alleles

            BNXT32 (*state->runs) = 0;     // initial next run - its length is 0 so far
        }
    }
}

// this function is first called to compress the ht_matrix_ctx.
// but it creates no compressed data for ht_matrix_ctx - instead it generates one or more PBWT sections (one for each allele
// in the VB) with CODEC_PBWT. Since these sections have a higher did_i, this function will call again compressing these sections,
// and this time it will simply copy the data to z_data.
COMPRESS (codec_pbwt_compress)
{
    START_TIMER;

    PbwtState state = codec_pbwt_initialize_state (vb, &vb->runs_ctx->local, &vb->fgrc_ctx->local); 
 
    buf_alloc (vb, &vb->runs_ctx->local, 0, MAX_(vb->ht_per_line, vb->ht_matrix_ctx->local.len / 5 ), char, CTX_GROWTH, "contexts->local"); // initial allocation
    buf_alloc (vb, &vb->fgrc_ctx->local, 0, MAX_(vb->ht_per_line, vb->ht_matrix_ctx->local.len / 30), char, CTX_GROWTH, "contexts->local");
        
    ARRAY (Allele, ht_data, vb->ht_matrix_ctx->local);
    
    uint32_t num_lines = ht_data_len / vb->ht_per_line;

    for (uint32_t line_i=0; line_i < num_lines; line_i++) {

        codec_pbwt_calculate_permutation (&state, &ht_data[line_i * vb->ht_per_line], vb->ht_per_line, line_i==0);

        // grow local if needed (unlikely) to the worst case scenario - all ht foreground, no two consecutive are similar -> 2xlen runs, half of them fg runs
        buf_alloc (vb, &vb->runs_ctx->local, 2 * vb->ht_per_line, 0, uint32_t, CTX_GROWTH, "contexts->local"); 
        buf_alloc (vb, &vb->fgrc_ctx->local,     vb->ht_per_line, 0, uint32_t, CTX_GROWTH, "contexts->local"); 

        show_line; show_perm(&state); 
        bool backward_permuted_ht_line = line_i % 2; // even rows are forward, odd are backward - better run length encoding 

        codec_pbwt_run_len_encode (&state, vb->ht_per_line, backward_permuted_ht_line);
    }
  
    if (flag.show_alleles) show_runs (&state);   
    
    // add ht_matrix_ctx.len to the end of fgrc_ctx.local (this should really be in the section header, but we don't
    // want to change SectionHeaderCtx (now in genozip v11)
    buf_alloc (vb, &vb->fgrc_ctx->local, 2, 0, uint32_t, 1, "contexts->local");
    BNXT32 (vb->fgrc_ctx->local) = vb->ht_matrix_ctx->local.len & 0xffffffffULL; // 32 LSb
    BNXT32 (vb->fgrc_ctx->local) = vb->ht_matrix_ctx->local.len >> 32;           // 32 MSb

    buf_free (vb->codec_bufs[0]); // allocated in codec_pbwt_initialize_state

    // note: we created the data in PBWT contexts - no section should be created for the ht_matrix context it in the file
    buf_free (vb->ht_matrix_ctx->local); 
    *compressed_len = 0;

    COPY_TIMER_COMPRESS (compressor_pbwt);

    return true;
}

// ----------
// PIZ side
// ----------

static void codec_pbwt_decode_init_ht_matrix (VBlockP vb, const uint32_t *rc_data, uint32_t *rc_data_len)
{
    vb->ht_matrix_ctx = ECTX (_PBWT_HT_MATRIX);

    // extract the HT matrix size from the last 2 words (64 bits) of rc_data
    *rc_data_len -= 2;

    uint64_t uncompressed_len = (uint64_t)rc_data[*rc_data_len] | ((uint64_t)rc_data[*rc_data_len + 1] << 32);
    
    ASSERT (vb->lines.len && uncompressed_len, 
            "%s: Expecting num_lines=%u and uncompressed_len=%"PRIu64" to be >0", VB_NAME, vb->lines.len32, uncompressed_len);

    buf_alloc (vb, &vb->ht_matrix_ctx->local, 0, uncompressed_len, char, 1, "contexts->local");

    vb->ht_matrix_ctx->local.len = uncompressed_len;
    vb->ht_matrix_ctx->lcodec    = CODEC_PBWT;
    vb->ht_matrix_ctx->ltype     = LT_CODEC; // reconstruction will go to codec_pbwt_reconstruct as defined in codec_args for CODEC_PBWT

    vb->ht_per_line = uncompressed_len / vb->lines.len; // note: the matrix always includes data for every line, even for lines without a GT field (it would be '*')
}

#define next param // decode use param as "next";

static inline void pbwt_decode_one_line (VBlockP vb, PbwtState *state, uint32_t line_i)
{
    ASSPIZ0 (state->runs->len, "No runs found");

    uint32_t *runs = B32 (*state->runs, state->runs->next);
    uint32_t *start_runs = runs;

    if (!state->run_allele) state->run_allele = '0'; // first run in VB - always start with background

    // de-permute one line onto the ht_matrix
    Allele *ht_one_line = B(Allele, vb->ht_matrix_ctx->local, line_i * vb->ht_per_line);

    codec_pbwt_calculate_permutation (state, ht_one_line, vb->ht_per_line, line_i==0);
    
    bool backwards = line_i % 2;

    for (uint32_t ht_i=0; ht_i < vb->ht_per_line; ) { 
        uint32_t run_len = *runs;

        uint32_t line_part_of_run = MIN_(run_len, vb->ht_per_line - ht_i); // the run could be shared with the next line
        
        for (uint32_t i=0; i < line_part_of_run; i++, ht_i++) {
            uint32_t oriented_ht_i = backwards ? vb->ht_per_line - ht_i - 1 : ht_i;
            ht_one_line[state->perm[oriented_ht_i].index] = state->perm[oriented_ht_i].allele = state->run_allele;
        }

        // case: we consumed only part of the run - leave the remaining part for the next line
        if (line_part_of_run < run_len) *runs -= line_part_of_run; 

        // case: run full consumed - move to the next run 
        if (line_part_of_run == run_len) {
            runs++; 
            
            // case: this run was background - next is foreground - pick the right foreground
            if (state->run_allele == '0') {
                PbwtFgRunCount *rc = B(PbwtFgRunCount, *state->fgrc, state->fgrc->next);
                state->run_allele = rc->fg_allele;
                rc->count--; // consume one
                if (!rc->count) state->fgrc->next++; // next run is the last run with this fg, afterwards, its a new fg
            }

            // case: this run was foreground - next it background
            else
                state->run_allele = '0';
        }
    }

    show_perm(state); show_line; 

    state->runs->next += runs - start_runs;
    ASSERT (state->runs->next <= state->runs->len, "%s: state.runs->next=%u is out of range", 
            VB_NAME, (uint32_t)state->runs->next);
}

// this function is called for the FORMAT_PBWT_FGRC section - after the FORMAT_PBWT_RUNS was already decompressed
UNCOMPRESS (codec_pbwt_uncompress)
{
    START_TIMER;

    uint32_t *rc_data = (uint32_t *)compressed;
    uint32_t rc_data_len = compressed_len / sizeof (PbwtFgRunCount);

    for (uint32_t i=0; i < rc_data_len; i++) 
        rc_data[i] = BGEN32 (rc_data[i]);

    // retrieve uncompressed_len stored at the end of the compressed data, initiatlize ht_matrix_ctx
    codec_pbwt_decode_init_ht_matrix (vb, rc_data, &rc_data_len);

    Context *runs_ctx = ECTX (_PBWT_RUNS); // might be different did_i for different data types
    ASSERT (runs_ctx, "%s: \"%s\": Cannot find context for PBWT_RUNS", VB_NAME, name);

    PbwtState state = codec_pbwt_initialize_state (vb, &runs_ctx->local, &vb->scratch); // background allele is provided to us in param ; this is a subcodec, so vb->scratch.data == compressed

    if (flag.show_alleles) show_runs (&state);

    state.runs->next = state.fgrc->next = 0; 

    // generate ht_matrix
    for (uint32_t line_i=0; line_i < vb->lines.len32; line_i++) 
        pbwt_decode_one_line (vb, &state, line_i); 	

    buf_free (vb->codec_bufs[0]); // free state data  
    buf_free (vb->scratch);    // free PBWT_ALLELES data

    COPY_TIMER (compressor_pbwt);
}

#undef next

CODEC_RECONSTRUCT (codec_pbwt_reconstruct)
{
    // find next allele - skipping unused spots ('*')
    Allele ht = '*';
    do { 
        ht = NEXTLOCAL (Allele, vb->ht_matrix_ctx);
    } while (ht == '*' && vb->ht_matrix_ctx->next_local < vb->ht_matrix_ctx->local.len);

    if (vb->drop_curr_line) return;

    if (IS_DIGIT(ht) || ht == '.') 
        RECONSTRUCT1 (ht);
    
    else if (ht == '-') 
        vb->txt_data.len--; // ploidy padding (starting 10.0.2) - appears as the 2nd+ HT - counted in GT.repeats: don't reconstruct anything, just remove previous phase character

    // % means "replace the previous phase and insert .". this is meant to reduce entroy for the GT.b250
    // in the case we have phased data, but missing samples appear as "./.". With this, we seg ".|%" instead of "./."
    // and hence the repsep and hence the entire GT container are the same beteen "./." and eg "1|0" samples.
    else if (ht == '%') {
        if (*BLSTtxt == '|' || *BLSTtxt == '/') { // second % in "%|%" 
            vb->txt_data.len -= 1;  // remove | or /
            RECONSTRUCT ("/.", 2);
        }
        else
            RECONSTRUCT1 ('.'); // first % in "%|%" 
    }

    else  // allele 10 to 99 (ascii 58 to 147)
        RECONSTRUCT_INT (ht - '0');
}
