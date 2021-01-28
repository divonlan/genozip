// ------------------------------------------------------------------
//   codec_pbwt.c
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// The codec implements a modified version of the PBWT algorithm. The permutation logic is loosely based on the logic in
// Durbin R. Efficient haplotype matching and storage using the positional Burrows-Wheeler transform (PBWT). Bioinformatics. 2014 

#include "genozip.h"
#include "codec.h"
#include "buffer.h"
#include "file.h"
#include "vblock.h"
#include "dict_id.h"
#include "reconstruct.h"
#include "vcf_private.h"
#include "strings.h"
#include "compressor.h"
#include "profiler.h"
#include "context.h"

typedef struct {
    uint8_t allele;     // copied from line on which permutation is applied
    uint32_t index;     // index into ht_one_line from which this allele was copied
} PermEnt;

typedef struct {
    Buffer *runs;       // array of uint32_t - alternating bg ('0') / fg runs
    Buffer *fgrc;       // array PbwtFgRunCount describing the fg runs
    uint8_t run_allele; // current allele for which we're constructing a run 
    PermEnt *perm;      // a permutation - expressed as indices into the haplotype line 
    PermEnt *temp;      // working memory
} PbwtState;

// this struct is part of the file format
typedef struct __attribute__ ((__packed__)) {
    uint32_t fg_allele : 8;
    uint32_t count     : 24;   // number of consecutive runs in PbwtState->runs with fg_allele 
} PbwtFgRunCount;

#define htputc(c) fprintf (info_stream, IS_NON_WS_PRINTABLE(c) ? "%c" : "\\x%02x", (c))

static void show_runs (const PbwtState *state)
{
    fprintf (info_stream,  "RUNS: ");
    
    PbwtFgRunCount rc = {};
    PbwtFgRunCount *next_allele = FIRSTENT (PbwtFgRunCount, *state->fgrc);

    for (unsigned i=0, fg=0; i < state->runs->len; i++, fg=!fg) {

        if (fg && !rc.count) {
            ASSERTE0 (next_allele < AFTERENT (PbwtFgRunCount, *state->fgrc), "premature end of alleles buffer");
            rc = *next_allele++;
        }
    
        iputc ('\''); 
        htputc (fg ? rc.fg_allele : '0');
        fprintf (info_stream, "':%u ", *ENT (uint32_t, *state->runs, i));

        if (fg) rc.count--;
    }
    iprint0 ("\n");
}

#define SHOW(msg, cmd) do { \
    fprintf (info_stream, msg " %-2u: ", line_i); \
    for (uint32_t i=0; i < vb->ht_per_line; i++) cmd; \
    iprint0 ("\n"); } while (0) // flush

#define show_line    if (flag.show_alleles) SHOW ("LINE", (htputc (*ENT (uint8_t, vb->ht_matrix_ctx->local, line_i * vb->ht_per_line + i))))
#define show_perm(s) if (flag.show_alleles) SHOW ("PERM", (fprintf (info_stream, "%d ", (s)->perm[i].index)));       

static PbwtState codec_pbwt_initialize_state (VBlockVCFP vb, Buffer *runs, Buffer *fgrc)
{
    buf_alloc (vb, &vb->codec_bufs[0], vb->ht_per_line * 2 * sizeof (uint32_t), 1, "codec_bufs");
    buf_zero (&vb->codec_bufs[0]); // re-zero every time
    ARRAY (PermEnt, state_data, vb->codec_bufs[0]);

    PbwtState state = {
        .runs = runs,
        .fgrc = fgrc,
        .perm = &state_data[0],               // size: ht_per_line X uint32_t
        .temp = &state_data[vb->ht_per_line], // size: ht_per_line X uint32_t
    };
        
    return state;
} 

// update the permutation for the next row: we re-sort it to make indices containing the same allele grouped
// first '0', then '1' etc - but the order within each of these allele groups remains as in the current row's premutation 
// (this is why we traverse the permuted line rather than the ht_matrix line)
static void inline codec_pbwt_calculate_permutation (PbwtState *state, const uint8_t *line, uint32_t line_len, bool is_first_line)
{
    uint32_t temp_i=0;        

    // populate permutation index - by re-ordering according to previous line's alleles
    if (!is_first_line) {

        // all valid allele values - in the order we want them in the permutation
        #define NUM_VALIDS 104
        static const uint8_t valids[NUM_VALIDS] = { 
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
        for (uint8_t i=0; i < NUM_VALIDS; i++) {
            uint8_t bg = valids[i];
        
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
    if (command == ZIP) 
        for (uint32_t ht_i=0; ht_i < line_len; ht_i++)
            state->perm[ht_i].allele = line[state->perm[ht_i].index];
}

// -------------
// ZIP side
// -------------

static inline DictId codec_pbwt_allele_dict_id (uint8_t allele)
{
    return dict_id_make ((char[]){ '@', allele, '-','P','B','W','T' }, 7, DTYPE_VCF_FORMAT);
}

// called by vcf_seg_finalize to create all contexts - must be done before merge
void codec_pbwt_comp_init (VBlock *vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    vb->ht_matrix_ctx->lcodec = CODEC_PBWT; // this will trigger codec_pbwt_compress even though the section is not written to the file

    DidIType st_did_i = ctx_get_existing_ctx (vb, dict_id_FORMAT_GT)->did_i;

    // create a contexts (if one doesn't exist already) for all alleles observed in this VB. 
    // Note: we can't do this in codec_pbwt_compress because it must be done before merge, do be copied to z_file->contexts

    // this context will contain alternate run lengths of the background allele and any forward allele
    Context *runs_ctx = ctx_get_ctx (vb, dict_id_PBWT_RUNS);
    runs_ctx->st_did_i      = st_did_i;   // in --stats, consoliate stats into GT
    runs_ctx->ltype         = LT_UINT32;
    
    // this context is used to determine which forward allele is in each forward run in RUNS: it contains
    // an array of PbwtFgRunCount
    Context *allele_ctx = ctx_get_ctx (vb, dict_id_PBWT_FGRC);
    allele_ctx->st_did_i      = st_did_i;   // in --stats, consoliate stats into GT
    allele_ctx->ltype         = LT_UINT32;
    allele_ctx->lsubcodec_piz = CODEC_PBWT;
}

// updates FGRC - our list of foreground run alleles. Each entry represents "count" consecutive runs of this same "fg_allele"
// returns true if a background run needs to be added to account for a missing background run between two foreground runs of different alleles
static bool inline codec_pbwt_udpate_fgrc (PbwtState *state, uint8_t allele_of_completed_run)
{
    // case: new run a foreground allele as previous run was a background - update foreground run count
    if (allele_of_completed_run == '0') {
        PbwtFgRunCount *rc;

        // case: previous foreground run was same allele as this - increment the count
        if (  state->fgrc->len && // we already have foreground runs
                state->run_allele == (rc = LASTENT (PbwtFgRunCount, *state->fgrc))->fg_allele) { 
            rc->count++; // one more consecutive foreground run of this allele
            ASSERTE (rc->count, "Number of consecutive forward runs with allele=%c exceeded the maximum of %u",
                    state->run_allele, 0xffffff); // error is count round-robined back to 0
        }
        // case: the previous fg run (if any) was of a different allele - start RunCount entry
        else
            NEXTENT (PbwtFgRunCount, *state->fgrc) = (PbwtFgRunCount){ .fg_allele = state->run_allele, .count = 1 };

        return false; // no need to add a background run
    }

    // case: this run is fg, and prev run was fg too (or this is the first run: run_allele==0) and too - skip the background run
    else if (state->run_allele != '0') {
        NEXTENT (PbwtFgRunCount, *state->fgrc) = (PbwtFgRunCount){ .fg_allele = state->run_allele, .count = 1 }; // new foreground
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
            *LASTENT (uint32_t, *state->runs) += run_len;
        
        // case: we need to start new run 
        else {
            uint8_t allele_of_completed_run = state->run_allele;
            state->run_allele = ALLELE(ht_i); // new run allele            

            if (codec_pbwt_udpate_fgrc (state, allele_of_completed_run)) 
                NEXTENT (uint32_t, *state->runs) = 0; // entry for missing background run (length 0) between runs of 2 different foreground alleles

            NEXTENT (uint32_t, *state->runs) = 0;     // initial next run - its length is 0 so far
        }
    }
}

// this function is first called to compress the ht_matrix_ctx, as we set its codec to CODEC_PBWT in codec_gtshark_comp_init.
// but it creates no compressed data for ht_matrix_ctx - instead it generates one or more PBWT sections (one for each allele
// in the VB) with CODEC_PBWT. Since these sections have a higher did_i, this function will call again compressing these sections,
// and this time it will simply copy the data to z_data.
bool codec_pbwt_compress (VBlock *vb_, 
                          SectionHeader *header,    
                          const char *uncompressed, // option 1 - compress contiguous data
                          uint32_t *uncompressed_len, 
                          LocalGetLineCB callback,  // option 2 - not supported
                          char *compressed, uint32_t *compressed_len /* in/out */, 
                          bool soft_fail)           // soft fail not supported
{
    START_TIMER;

    VBlockVCF *vb = (VBlockVCF *)vb_;

    Context *runs_ctx = ctx_get_existing_ctx (vb, dict_id_PBWT_RUNS);
    Context *fgrc_ctx = ctx_get_existing_ctx (vb, dict_id_PBWT_FGRC);

    PbwtState state = codec_pbwt_initialize_state (vb, &runs_ctx->local, &fgrc_ctx->local); 
 
    buf_alloc (vb, &runs_ctx->local, MAX (vb->ht_per_line, vb->ht_matrix_ctx->local.len / 5 ), CTX_GROWTH, "contexts->local"); // initial allocation
    buf_alloc (vb, &fgrc_ctx->local, MAX (vb->ht_per_line, vb->ht_matrix_ctx->local.len / 30), CTX_GROWTH, "contexts->local");
        
    uint8_t *ht_data = FIRSTENT (uint8_t, vb->ht_matrix_ctx->local);
    
    for (uint32_t line_i=0; line_i < vb->lines.len; line_i++) {

        codec_pbwt_calculate_permutation (&state, &ht_data[line_i * vb->ht_per_line], vb->ht_per_line, line_i==0);

        // grow local if needed (unlikely) to the worst case scenario - all ht foreground, no two consecutive are similar -> 2xlen runs, half of them fg runs
        buf_alloc_more (vb, &runs_ctx->local, 2 * vb->ht_per_line, 0, uint32_t, CTX_GROWTH, "contexts->local"); 
        buf_alloc_more (vb, &fgrc_ctx->local,     vb->ht_per_line, 0, uint32_t, CTX_GROWTH, "contexts->local"); 

        show_line; show_perm(&state); 
        bool backward_permuted_ht_line = line_i % 2; // even rows are forward, odd are backward - better run length encoding 

        codec_pbwt_run_len_encode (&state, vb->ht_per_line, backward_permuted_ht_line);
    }
  
    if (flag.show_alleles) show_runs (&state);   
    
    // add ht_matrix_ctx.len to the end of allele_ctx.local (this should really be in the section header, but we don't
    // want to change SectionHeaderCtx (now in genozip v11)
    buf_alloc_more (vb, &fgrc_ctx->local, 2, 0, uint32_t, 1, "contexts->local");
    NEXTENT (uint32_t, fgrc_ctx->local) = vb->ht_matrix_ctx->local.len & 0xffffffffULL; // 32 LSb
    NEXTENT (uint32_t, fgrc_ctx->local) = vb->ht_matrix_ctx->local.len >> 32;           // 32 MSb

    BGEN_u32_buf (&runs_ctx->local, NULL);
    BGEN_u32_buf (&fgrc_ctx->local, NULL);

    buf_free (&vb->codec_bufs[0]); // allocated in codec_pbwt_initialize_state

    // the allele sections are further compressed with the best simple codec 
    // note: this simple codec (not CODEC_PBWT) will be the codec stored in zf_ctx->lcodec
    PAUSE_TIMER; //  don't include sub-codec compressor - it accounts for itself
    codec_assign_best_codec (vb_, runs_ctx, NULL, SEC_LOCAL);
    codec_assign_best_codec (vb_, fgrc_ctx, NULL, SEC_LOCAL);
    RESUME_TIMER (compressor_pbwt);

    // note: we created the data in PBWT contexts - no section should be created for the ht_matrix context it in the file
    buf_free (&vb->ht_matrix_ctx->local); 
    *compressed_len = 0;

    COPY_TIMER (compressor_pbwt);

    return true;
}

// ----------
// PIZ side
// ----------

static void codec_pbwt_decode_init_ht_matrix (VBlockVCF *vb, const uint32_t *rc_data, uint32_t *rc_data_len)
{
    vb->ht_matrix_ctx = ctx_get_ctx (vb, dict_id_FORMAT_GT_HT); // create new context - it doesn't exist in the genozip file

    // extract the HT matrix size from the last 2 words (64 bits) of rc_data
    *rc_data_len -= 2;

    uint64_t uncompressed_len = (uint64_t)rc_data[*rc_data_len] | ((uint64_t)rc_data[*rc_data_len + 1] << 32);
    
    ASSERTE (vb->lines.len && uncompressed_len && vcf_header_get_num_samples(), 
                "Expecting num_lines=%u, num_samples=%u and uncompressed_len=%"PRIu64" to be >0", (uint32_t)vb->lines.len, vcf_header_get_num_samples(), uncompressed_len);

    buf_alloc (vb, &vb->ht_matrix_ctx->local, uncompressed_len, 1, "contexts->local");

    vb->ht_matrix_ctx->local.len = uncompressed_len;
    vb->ht_matrix_ctx->lcodec    = CODEC_PBWT;
    vb->ht_matrix_ctx->ltype     = LT_CODEC; // reconstruction will go to codec_pbwt_reconstruct as defined in codec_args for CODEC_PBWT

    vb->ht_per_line = uncompressed_len / vb->lines.len; // note: the matrix always includes data for every line, even for lines without a GT field (it would be '*')
    vb->ploidy      = vb->ht_per_line / vcf_header_get_num_samples();
}

#define next param // decode use param as "next";

static inline void pbwt_decode_one_line (VBlockVCFP vb, PbwtState *state, uint32_t line_i)
{
    uint32_t *runs = ENT (uint32_t, *state->runs, state->runs->next);
    uint32_t *start_runs = runs;

    if (!state->run_allele) state->run_allele = '0'; // first run in VB - always start with background

    // de-permute one line onto the ht_matrix
    uint8_t *ht_one_line = ENT (uint8_t, vb->ht_matrix_ctx->local, line_i * vb->ht_per_line);

    codec_pbwt_calculate_permutation (state, ht_one_line, vb->ht_per_line, line_i==0);
    
    bool backwards = line_i % 2;

    for (uint32_t ht_i=0; ht_i < vb->ht_per_line; ) { 
        uint32_t run_len = *runs;

        uint32_t line_part_of_run = MIN (run_len, vb->ht_per_line - ht_i); // the run could be shared with the next line
        
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
                PbwtFgRunCount *rc = ENT (PbwtFgRunCount, *state->fgrc, state->fgrc->next);
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
    ASSERTE (state->runs->next <= state->runs->len, "state.runs->next=%u is out of range", (uint32_t)state->runs->next);
}


// this function is called for the PBWT_ALLELES section - after the PBWT_RUNS was already decompressed
void codec_pbwt_uncompress (VBlock *vb_, Codec codec, uint8_t unused /* param */,
                            const char *compressed, uint32_t compressed_len,
                            Buffer *uncompressed_buf, uint64_t uncompressed_len,
                            Codec sub_codec)
{
    START_TIMER;

    VBlockVCF *vb = (VBlockVCF *)vb_;
    uint32_t *rc_data = (uint32_t *)compressed;
    uint32_t rc_data_len = compressed_len / sizeof (PbwtFgRunCount);

    for (uint32_t i=0; i < rc_data_len; i++) 
        rc_data[i] = BGEN32 (rc_data[i]);

    // retrieve uncompressed_len stored at the end of the compressed data, initiatlize ht_matrix_ctx
    codec_pbwt_decode_init_ht_matrix (vb, rc_data, &rc_data_len);

    Context *runs_ctx = ctx_get_existing_ctx (vb, dict_id_PBWT_RUNS);
    
    PbwtState state = codec_pbwt_initialize_state (vb, &runs_ctx->local, &vb->compressed); // background allele is provided to us in param ; this is a subcodec, so vb->compressed.data == compressed

    if (flag.show_alleles) show_runs (&state);

    state.runs->next = state.fgrc->next = 0; 

    // generate ht_matrix
    for (uint32_t line_i=0; line_i < vb->lines.len; line_i++) 
        pbwt_decode_one_line (vb, &state, line_i); 	

    buf_free (&vb->codec_bufs[0]); // free state data  
    buf_free (&vb->compressed);    // free PBWT_ALLELES data

    COPY_TIMER (compressor_pbwt);
}

#undef next

void codec_pbwt_reconstruct (VBlock *vb_, Codec codec, Context *ctx)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    // find next allele - skipping unused spots ('*')
    uint8_t ht = '*';
    do { 
        ht = NEXTLOCAL (uint8_t, vb->ht_matrix_ctx);
    } while (ht == '*' && vb->ht_matrix_ctx->next_local < vb->ht_matrix_ctx->local.len);

    if (vb->dont_show_curr_line) return;

    if (IS_DIGIT(ht) || ht == '.') 
        RECONSTRUCT1 (ht);
    
    else if (ht == '-') 
        vb->txt_data.len--; // ploidy padding (starting 10.0.2) - appears as the 2nd+ HT - counted in GT.repeats: don't reconstruct anything, just remove previous phase character

    // % means "replace the previous phase and insert .". this is meant to reduce entroy for the GT.b250
    // in the case we have phased data, but missing samples appear as "./.". With this, we seg ".|%" instead of "./."
    // and hence the repsep and hence the entire GT container are the same beteen "./." and eg "1|0" samples.
    else if (ht == '%') {
        if (*LASTENT (char, vb->txt_data) == '|' || *LASTENT (char, vb->txt_data) == '/') { // second % in "%|%" 
            vb->txt_data.len -= 1;  // remove | or /
            RECONSTRUCT ("/.", 2);
        }
        else
            RECONSTRUCT1 ('.'); // first % in "%|%" 
    }

    else { // allele 10 to 99 (ascii 58 to 147)
        RECONSTRUCT_INT (ht - '0');
    }
}
