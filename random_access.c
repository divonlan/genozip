// ------------------------------------------------------------------
//   random_access.h
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#include "buffer.h"
#include "random_access.h"
#include "vblock.h"
#include "file.h"
#include "endianness.h"
#include "regions.h"
#include "mutex.h"
#include "zfile.h"
#include "strings.h"
#include "zfile.h"
#include "sections.h"
#include "codec.h"

static Mutex ra_mutex[2] = {};

#define RA_UNKNOWN_CHROM_SKIP_POS 1

// --------------------
// misc stuff
// --------------------

void random_access_initialize(void)
{    
    mutex_initialize (ra_mutex[0]);
    mutex_initialize (ra_mutex[1]);
}

// Called by PIZ main thread (piz_read_global_area) and ZIP main thread (zip_write_global_area)
static void BGEN_random_access (BufferP ra_buf)
{
    for_buf (RAEntry, ra, *ra_buf) {
        ra->vblock_i    = BGEN32 (ra->vblock_i);
        ra->chrom_index = BGEN32 (ra->chrom_index);
        ra->min_pos     = BGEN64 (ra->min_pos);
        ra->max_pos     = BGEN64 (ra->max_pos);
    }
}

static void random_access_show_index (ConstBufferP ra_buf, bool from_zip, Did chrom_did_i, rom msg)
{
    iprintf ("\n%s:\n", msg);
    
    Context *ctx = ZCTX (chrom_did_i);

    for_buf (const RAEntry, ra, *ra_buf) {
        
        STR(chrom_snip);
        if (from_zip) {
            if (ra->chrom_index != WORD_INDEX_NONE) {
                CtxNode *chrom_node = ctx_get_node_by_word_index (ctx, ra->chrom_index);
                chrom_snip     = Bc (ctx->dict, chrom_node->char_index);
                chrom_snip_len = chrom_node->snip_len;
            }
            else {
                chrom_snip     = "CONT_FROM_PREV_VB";
                chrom_snip_len = strlen (chrom_snip);
            }
        }
        else {
            CtxWord *chrom_word = B(CtxWord, ctx->word_list, ra->chrom_index);
            chrom_snip = Bc (ctx->dict, chrom_word->index);
            chrom_snip_len = chrom_word->len;
        }
        iprintf ("vb_i=%u chrom='%.*s' (chrom_word_index=%d) min_pos=%"PRId64" max_pos=%"PRId64"\n",
                 ra->vblock_i, chrom_snip_len, chrom_snip, ra->chrom_index, ra->min_pos, ra->max_pos);
    }
}

// --------------------
// ZIP stuff
// --------------------

// ZIP only
void random_access_alloc_ra_buf (VBlockP vb, int ra_i, int32_t chrom_node_index)
{
    vb->ra_buf[ra_i].len = MAX_(vb->ra_buf[ra_i].len, chrom_node_index + 2);
    buf_alloc_zero (vb, &vb->ra_buf[ra_i], 0, MAX_(vb->ra_buf[ra_i].len, 500), RAEntry, 2, "ra_buf");
}

// ZIP only: called from Seg when the CHROM changed - this might be a new chrom, or
// might be an exiting chrom (for example, in an unsorted SAM or VCF). we maitain one ra field per chrom per vb
// chrom_node_index==WORD_INDEX_NONE if we don't yet know what chrom this is and we will update it later (see fasta_seg_txt_line)
void random_access_update_chrom (VBlockP vb, int ra_i, WordIndex chrom_node_index, STRp (chrom_name))
{
    // note: when FASTA calls this for a sequence that started in the previous vb, and hence chrom is unknown, chrom_node_index==-1.
    ASSERT (chrom_node_index >= -1, "chrom_node_index=%d in vb_i=%u", chrom_node_index, vb->vblock_i);

    // if this is an "unavailable" chrom ("*" in SAM, "." in VCF) we don't store it and signal not to store POS either
    if (chrom_name_len == 1 && (*chrom_name == '*' || *chrom_name == '.')) {
        vb->ra_buf[ra_i].param = RA_UNKNOWN_CHROM_SKIP_POS;
        vb->chrom_name       = chrom_name;
        vb->chrom_name_len   = chrom_name_len;
        return;
    }

    // allocate
    random_access_alloc_ra_buf (vb, ra_i, chrom_node_index);

    RAEntry *ra_ent = B(RAEntry, vb->ra_buf[ra_i], chrom_node_index + 1); // chrom_node_index=-1 goes into entry 0 etc
    if (!ra_ent->vblock_i) { // first occurance of this chrom
        ra_ent->chrom_index = chrom_node_index;
        ra_ent->vblock_i    = vb->vblock_i;
    }
}

// ZIP only: update the pos in the existing chrom entry
void random_access_update_pos (VBlockP vb, int ra_i, Did did_i_pos)
{
    // if the last chrom was unknown ("*"), we do nothing with pos either
    if (vb->ra_buf[ra_i].param == RA_UNKNOWN_CHROM_SKIP_POS) {
        vb->ra_buf[ra_i].param = 0; // reset
        return;
    }
    
    PosType this_pos = vb->last_int (did_i_pos);

    if (this_pos <= 0) return; // ignore pos<=0 (in SAM, 0 means unmapped POS)

    RAEntry *ra_ent = B(RAEntry, vb->ra_buf[ra_i], vb->chrom_node_index + 1); // chrom_node_index=-1 goes into entry 0 etc

    if (!ra_ent->min_pos) // first line this chrom is encountered in this vb - initialize pos
        ra_ent->min_pos = ra_ent->max_pos = this_pos;

    else if (this_pos < ra_ent->min_pos) ra_ent->min_pos = this_pos; 
    
    else if (this_pos > ra_ent->max_pos) ra_ent->max_pos = this_pos;
}

// ZIP: update increment reference pos
void random_access_increment_last_pos (VBlockP vb, int ra_i, PosType increment)
{
    ASSERTISALLOCED (vb->ra_buf[ra_i]);

    RAEntry *ra_ent = B(RAEntry, vb->ra_buf[ra_i], vb->chrom_node_index + 1); // chrom_node_index=-1 goes into entry 0 etc

    if (!ra_ent->min_pos) ra_ent->min_pos = 1; 
    ra_ent->max_pos += increment;
}

// ZIP: update last reference pos
void random_access_update_last_pos (VBlockP vb, int ra_i, PosType last_pos)
{
    ASSERTISALLOCED (vb->ra_buf[ra_i]);
    
    RAEntry *ra_ent = B(RAEntry, vb->ra_buf[ra_i], vb->chrom_node_index + 1); // chrom_node_index=-1 goes into entry 0 etc
    if (last_pos > ra_ent->max_pos) ra_ent->max_pos = last_pos;
}

void random_access_update_first_last_pos (VBlockP vb, int ra_i, WordIndex chrom_node_index, STRp (first_pos), STRp (last_pos))
{
    ASSERTISALLOCED (vb->ra_buf[ra_i]);

    PosType first_pos_value, last_pos_value;    
    if (!str_get_int (STRa(first_pos), &first_pos_value)) return; // fail silently
    if (!str_get_int (STRa(last_pos),  &last_pos_value )) return; 
    
    RAEntry *ra_ent = B(RAEntry, vb->ra_buf[ra_i], chrom_node_index + 1); // chrom_node_index=-1 goes into entry 0 etc
    if (first_pos_value < ra_ent->min_pos) ra_ent->min_pos = first_pos_value;
    if (last_pos_value  > ra_ent->max_pos) ra_ent->max_pos = last_pos_value;
}

void random_access_update_to_entire_chrom (VBlockP vb, int ra_i, PosType first_pos_of_chrom, PosType last_pos_of_chrom)
{
    ASSERTISALLOCED (vb->ra_buf[ra_i]);
    
    RAEntry *ra_ent = B(RAEntry, vb->ra_buf[ra_i], vb->chrom_node_index + 1); // chrom_node_index=-1 goes into entry 0 etc
    ra_ent->min_pos = first_pos_of_chrom;
    ra_ent->max_pos = last_pos_of_chrom;
}

// called by ZIP compute thread, while holding the z_file mutex: merge in the VB's ra_buf in the global z_file one
// note: the order of the merge is not necessarily the sequential order of VBs
void random_access_merge_in_vb (VBlockP vb, int ra_i)
{
    if (!vb->ra_buf[ra_i].len) return; // nothing to merge

    BufferP z_buf  = ra_i==0 ? &z_file->ra_buf : &z_file->ra_buf_luft;
    ARRAY (const RAEntry, src_ra, vb->ra_buf[ra_i]);
     
    mutex_lock (ra_mutex[ra_i]);

    START_TIMER; // not including mutex wait time

    buf_alloc (evb, z_buf, 0, z_buf->len + src_ra_len, RAEntry, 2, "z_file->ra_buf"); 

    Context *chrom_ctx = CTX(ra_i==0 ? DTF(prim_chrom) : DTF(luft_chrom));
    ASSERT0 (chrom_ctx, "cannot find chrom_ctx");

    for (unsigned i=0; i < src_ra_len; i++) {
        
        if (!src_ra[i].vblock_i) continue; // not used in this VB

        RAEntry *dst_ra = &BNXT (RAEntry, *z_buf);

        dst_ra->vblock_i = vb->vblock_i;
        dst_ra->min_pos  = src_ra[i].min_pos;
        dst_ra->max_pos  = src_ra[i].max_pos;

        if (src_ra[i].chrom_index != WORD_INDEX_NONE) {
            CtxNode *chrom_node = ctx_node_vb (chrom_ctx, (WordIndex)src_ra[i].chrom_index, NULL, NULL);

            //if (chrom_node->word_index.n == WORD_INDEX_NONE) { // this contig was canceled by seg_rollback
            if (chrom_node->word_index == WORD_INDEX_NONE) { // this contig was canceled by seg_rollback
                z_buf->len--;
                continue;
            }
            
            // dst_ra->chrom_index = chrom_node->word_index.n; // note: in the VB we store the node index, while in zfile we store tha word index
            dst_ra->chrom_index = chrom_node->word_index; // note: in the VB we store the node index, while in zfile we store tha word index
        }
        else 
            dst_ra->chrom_index = WORD_INDEX_NONE; // to be updated in random_access_finalize_entries()
    }

    COPY_TIMER(random_access_merge_in_vb);

    mutex_unlock (ra_mutex[ra_i]);
}

// ZIP
static BufferP ra_buf_being_sorted;
int random_access_sort_by_vb_i (const void *a_, const void *b_)
{
    ARRAY (RAEntry, ra, *ra_buf_being_sorted);

    int32_t a = *(int32_t *)a_;
    int32_t b = *(int32_t *)b_;
    
    int vb_i_diff = (int)ra[a].vblock_i - (int)ra[b].vblock_i;

    if (vb_i_diff)
        return vb_i_diff;   // from lower to higher vb_i
    else
        return a - b; // for same vb_i, keep order as it currently is in array (order will be from lower to higher address)
}

// ZIP (main thread) sort RA, update overflowing chroms, create and merge in evb ra
void random_access_finalize_entries (BufferP ra_buf)
{
    START_TIMER;

    if (!ra_buf->len) return; // no random access

    // build an index into ra_buf that we will sort. we need that, because for same-vb entries we need to 
    // maintain their current order - sorted by the index
    int32_t *sorter = MALLOC (ra_buf->len * sizeof (int32_t));

    for (int32_t i=0; i < ra_buf->len; i++) sorter[i] = i;

    // at this point, VBs in RA might be out of order, but all entries of a specific VB are grouped together in order
    // we sort so that the VBs are in order, and the order within a VB is unchanged
    ra_buf_being_sorted = ra_buf;
    qsort (sorter, ra_buf->len, sizeof (uint32_t), random_access_sort_by_vb_i);

    // use sorter to consturct a sorted RA
    static Buffer sorted_ra_buf = EMPTY_BUFFER; // must be static because its added to buf_list
    buf_alloc (evb, &sorted_ra_buf, 0, ra_buf->len, RAEntry, 1, ra_buf->name);
    sorted_ra_buf.len = ra_buf->len;

    for (uint32_t i=0; i < ra_buf->len32; i++) 
        *B(RAEntry, sorted_ra_buf, i) = *B(RAEntry, *ra_buf, sorter[i]);

    // replace ra_buf with sorted one
    buf_destroy (*ra_buf);
    buf_move (evb, ra_buf, evb, &sorted_ra_buf);

    FREE (sorter);

    // now that the VBs are in order, we can updated the "CONTINUED FROM PREVIOUS VB"(=WORD_INDEX_NONE) ra's to their final values
    if (Z_DT(FASTA) || Z_DT(REF)) 
        for_buf2 (RAEntry, ra, i, *ra_buf) 
            if (ra->chrom_index == WORD_INDEX_NONE) {
                // we expect this ra to be the first in its VB, and the previous ra to be of the previous VB
                ASSERT (i && (ra->vblock_i == (ra-1)->vblock_i+1), "corrupt ra[%u]: chrom_index=WORD_INDEX_NONE but vb_i=%u and (ra-1)->vb_i=%u (expecting it to be %u)",
                        i, ra->vblock_i, (ra-1)->vblock_i, ra->vblock_i-1);

                ra->chrom_index = (ra-1)->chrom_index;
                ra->min_pos    += (ra-1)->max_pos;
                ra->max_pos    += (ra-1)->max_pos;
            } 

    COPY_TIMER_VB (evb, random_access_finalize_entries);
}

Codec random_access_compress (ConstBufferP ra_buf_, SectionType sec_type, Codec codec, int ra_i, rom msg)
{
    START_TIMER;

    if (!ra_buf_->len) goto done; // no random access

    BufferP ra_buf = (BufferP)ra_buf_; // we will restore everything so that Const is honoured

    if (msg) random_access_show_index (ra_buf, true, ra_i ? DTFZ(luft_chrom) : DTFZ(prim_chrom), msg);
    
    BGEN_random_access (ra_buf); // make ra_buf into big endian

    SectionFlags sec_flags = {};
    if (ra_i) sec_flags.random_access.luft = true;

    // assigned codec to lens buffer (lens_ctx.local)
    ra_buf->len *= sizeof (RAEntry); // change len to count bytes
    
    if (codec == CODEC_UNKNOWN) codec = codec_assign_best_codec (evb, NULL, ra_buf, sec_type);
    if (codec == CODEC_UNKNOWN) codec = CODEC_NONE; // really small

    zfile_compress_section_data_ex (evb, NULL, sec_type, ra_buf, 0,0, codec, sec_flags, NULL); // ra data compresses better with LZMA than BZLIB

    // restore so we leave it intact as promised by Const
    ra_buf->len /= sizeof (RAEntry); // restore
    BGEN_random_access (ra_buf); // make ra_buf into big endian

done:
    COPY_TIMER_VB (evb, random_access_compress);
    return codec;
}

// --------------------
// PIZ stuff
// --------------------

// PIZ: binary search for first entry of vb_i, using the fact that in random_access_finalize_entries (in ZIP) we sorted the entries by vb_i
static const RAEntry *random_access_get_first_ra_of_vb_do (VBIType vb_i, const RAEntry *first_ra, const RAEntry *last_ra)
{
    if (first_ra > last_ra) return NULL; // vb_i not found

    const RAEntry *mid_ra = first_ra + ((last_ra - first_ra) / 2);
    bool first_ra_in_vb = (mid_ra == B1ST (RAEntry, z_file->ra_buf)) || (mid_ra->vblock_i != (mid_ra-1)->vblock_i);

    if (mid_ra->vblock_i < vb_i) 
        return random_access_get_first_ra_of_vb_do (vb_i, mid_ra+1, last_ra);

    if (mid_ra->vblock_i > vb_i || !first_ra_in_vb)
        return random_access_get_first_ra_of_vb_do (vb_i, first_ra, mid_ra-1);

    return mid_ra;
}
#define random_access_get_first_ra_of_vb(vb_i) random_access_get_first_ra_of_vb_do (vb_i, B1ST (RAEntry, z_file->ra_buf), BLST (RAEntry, z_file->ra_buf))

bool random_access_has_filter (void)
{
    return flag.regions &&      // --regions specified
           z_file->ra_buf.len;  // this file has RA data (eg not unaligned SAM/BAM)
}

// PIZ main thread (called from writer_init_vb_info): check if for the given VB,
// the ranges in random access (from the file) overlap with the ranges in regions (from the command line --regions)
bool random_access_is_vb_included (VBIType vb_i)
{
    if (!random_access_has_filter())  
        return true; // all VBs are included

    const RAEntry *ra = random_access_get_first_ra_of_vb (vb_i);
    
    // case: an entire VB without RA data while some other VBs do have. For example - a sorted SAM where unaligned reads are pushed to the end of the file
    if (!ra) return false; // don't include this VB

    for (uint32_t ra_i=0; ra->vblock_i == vb_i; ra_i++, ra++) 
        if (regions_get_ra_intersection (ra->chrom_index, ra->min_pos, ra->max_pos))
            return true; // vb is included

    return false; 
}

// FASTA PIZ
bool random_access_does_last_chrom_continue_in_next_vb (VBIType vb_i)
{
    const RAEntry *ra = random_access_get_first_ra_of_vb (vb_i+1);
    if (!ra) return false; // vb_i is the last vb (there is no vb_i+1) - so it doesn't continue in the next vb...

    return ra->chrom_index == (ra-1)->chrom_index;
}

// FASTA+FASTQ PIZ - number of chroms in this VB, excluding one that started before
uint32_t random_access_num_chroms_start_in_this_vb (VBIType vb_i)
{
    const RAEntry *ra = random_access_get_first_ra_of_vb (vb_i);
    ASSERT (ra, "no ra for vb_i%u", vb_i);

    // count the first RA if its its the first RA in the file OR it is the first RA with this chrom 
    // (NOT a continuation of the chrom of the last RA of the previous VB)
    int32_t count = (ra==B1ST (RAEntry, z_file->ra_buf)) ? 1 : (ra-1)->chrom_index != ra->chrom_index;

    for (ra=ra+1; ra < BAFT (RAEntry, z_file->ra_buf) && ra->vblock_i == vb_i; ra++) 
        count++; // count all additional chroms starting in this VB

    return count;
}


// FASTA PIZ main thread: check if all contigs have the same max pos, and return it
uint32_t random_access_verify_all_contigs_same_length (void)
{
    static Buffer max_lens_buf = EMPTY_BUFFER;
    const Context *ctx = ZCTX(CHROM);
    buf_alloc (evb, &max_lens_buf, 0, ctx->word_list.len, PosType, 1, "max_lens");
    buf_zero (&max_lens_buf);

    ARRAY (PosType, max_lens, max_lens_buf);

    PosType max_of_maxes=0;
    for_buf (const RAEntry, ra, z_file->ra_buf) {
        max_lens[ra->chrom_index] = MAX_(ra->max_pos, max_lens[ra->chrom_index]);
        max_of_maxes = MAX_(max_of_maxes, max_lens[ra->chrom_index]);
    }

    for (uint32_t contig_i=0; contig_i < ctx->word_list.len32; contig_i++) 
        ASSINP (max_lens[contig_i] == max_of_maxes, "file %s cannot be displayed in Phylip format because not all contigs are the same length: contig %s has length=%"PRIu64" which is shorter than other contigs of length=%"PRIu64,
                z_name, ctx_get_words_snip (ctx, contig_i), max_lens[contig_i], max_of_maxes);

    buf_free (max_lens_buf);

    return max_of_maxes;
}

// PIZ: read SEC_RANDOM_ACCESS. If --luft - read the luft section instead.
void random_access_load_ra_section (SectionType sec_type, Did chrom_did_i, BufferP ra_buf, rom buf_name, rom show_index_msg)
{
    Section ra_sl = sections_first_sec (sec_type, true);
    if (!ra_sl) return; // section doesn't exist

    // if we're pizzing as lift-over, get the next section which has the Luft random access info (zfile_get_global_section will verify)
    if (flag.luft) {
        ra_sl++; 
        if (ra_sl->st != sec_type) return; // no luft over section
    }
    
    zfile_get_global_section (SectionHeader, ra_sl, ra_buf, buf_name);
    
    if (ra_buf->len) {
        ra_buf->len /= sizeof (RAEntry);
        BGEN_random_access (ra_buf);
    }
    
    if (show_index_msg) {
        random_access_show_index (ra_buf, false, chrom_did_i, show_index_msg);
        if (is_genocat) exit_ok(); // in genocat --show-index, we only show the index, not the data
    }
}

// returns true if successful
void random_access_get_ra_info (VBIType vblock_i, WordIndex *chrom_index, PosType *min_pos, PosType *max_pos)
{
    const RAEntry *ra = random_access_get_first_ra_of_vb (vblock_i);
    ASSERT (ra, "vblock_i=%u has no RAEntry", vblock_i);

    *chrom_index = ra->chrom_index;
    *min_pos     = ra->min_pos;
    *max_pos     = ra->max_pos;
}