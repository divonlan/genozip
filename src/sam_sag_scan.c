// ------------------------------------------------------------------
//   sam_sag_scan.c
//   Copyright (C) 2022-2022 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

// In BAM files with SAG_BY_FLAG, absent an SA, NH or CC tag, we are unable to distiguish one of the alignments of
// being primary during, and unable to tell whether a non-secondary, non-supplementary alignment has dependent alignments or not.
// 
// This module scans the txt_file for depn lines, creating a sort-uniqed array of hash(QNAME) in z_file->sag_depn_index
// for each dependent line - during segconf, and a function for determining if a depn exists by looking up in this array during seg

#include "genozip.h"
#include "sam_private.h"
#include "sections.h"
#include "dispatcher.h"
#include "profiler.h"
#include "txtfile.h"
#include "bgzf.h"
#include "libdeflate/libdeflate.h"

static VBlockP real_evb;

#define vb_depn_index  vb->scratch
#define z_qname_index  evb->z_data

typedef struct { uint32_t hash; VBIType vb_i; } QnameIndexEnt;

static rom scan_index_one_line (VBlockSAMP vb, rom alignment, uint32_t remaining_txt_len)   
{
    STR(qname);
    SamFlags flag;
    uint32_t alignment_len;

    if (IS_BAM_ZIP) {
        alignment_len = GET_UINT32 (alignment) + 4;

        // a non-sensical block_size might indicate an false-positive identification of a BAM alignment in bam_unconsumed
        ASSERT (alignment_len >= sizeof (BAMAlignmentFixed) && alignment_len <= remaining_txt_len, 
                "%s: alignment_len=%u is out of range - too small, or goes beyond end of txt data: remaining_txt_len=%u",
                LN_NAME, alignment_len, remaining_txt_len);

        flag      = (SamFlags){ .value = GET_UINT16 (&alignment[18]) };
        qname_len = (uint8_t)alignment[12] - 1; // -1 bc without \0
        qname     = &alignment[36];
    }
    
    else { // SAM
        rom nl = memchr (alignment, '\n', remaining_txt_len);
        ASSERT (nl, "%s: last line in this VB has no newline", VB_NAME);

        alignment_len = nl - alignment + 1;

        rom tab = memchr (alignment, '\t', alignment_len);
        ASSERT (tab, "%s: line has no \\t after QNAME", LN_NAME);

        qname = alignment;
        qname_len = tab - qname;
        
        rom flag_str = tab+1;
        tab = memchr (flag_str, '\t', alignment_len - (qname_len+1));
        ASSERT (tab, "%s: line has no \\t after FLAG", VB_NAME);

        ASSERT (str_get_int_range16 (flag_str, tab - flag_str, 0, SAM_MAX_FLAG, &flag.value), "%s: invalid FLAG=%.*s", VB_NAME, (int)(tab - flag_str), flag_str);
    }

    uint32_t hash = QNAME_HASH (qname, qname_len, flag.is_last);

    if (vb->preprocessing) {
        // depn index
        if (!flag.unmapped && (flag.supplementary || flag.secondary)) {
            buf_alloc (vb, &vb_depn_index, 1, 50000, uint32_t, 2, "scratch");
            BNXT32 (vb_depn_index) = hash;
        }

        // qname index
        if (!flag.unmapped) {
            buf_alloc (vb, &vb->qname_count, 1, 5000000, QnameIndexEnt, 2, "qname_count");
            BNXT(QnameIndexEnt, vb->qname_count) = (QnameIndexEnt){ .hash=hash, .vb_i=vb->vblock_i };
        }
    }
    else {
        if (!flag.unmapped) {
            buf_alloc (vb, &vb->qname_count, 1, 0, QnameCount, 2, "qname_count");
            BNXT(QnameCount, vb->qname_count) = (QnameCount){ .hash=hash, .count=1 };
        }
    }
    return alignment + alignment_len;
}

//--------------------------------------------------------------------------
// Scanning during pre-processing - for SAM_BY_FLAG - so we can know whether
// a non-depn line is a prim that has depn's in other VBs
//--------------------------------------------------------------------------

static void scan_read_one_vb (VBlockP vb)
{
    txtfile_read_vblock (vb);

    if (vb->txt_data.len) {
        vb->preprocessing = true;
        vb->dispatch = READY_TO_COMPUTE;
        txt_file->num_vbs_dispatched++;
    }

    else 
        vb->dispatch = DATA_EXHAUSTED;
}

static SORTER (depn_index_sorter)
{
    return ASCENDING_RAW (*(uint32_t *)a, *(uint32_t *)b);
}

static void scan_sort_unique_depn_index (Buffer *buf)
{
    if (!buf->len) return;

    ARRAY (uint32_t, arr, *buf);

    // sort-uniq depn index    
    qsort (arr, arr_len, sizeof(uint32_t), depn_index_sorter);

    // unique
    uint32_t new_len=1;
    
    for (uint32_t i=1; i < arr_len; i++)
        if (arr[i] != arr[i-1]) 
            arr[new_len++] = arr[i];
    
    buf->len32 = new_len;
}

static ASCENDING_SORTER (qname_index_sorter, QnameIndexEnt, hash)

static void scan_sort_unique_qname_index (Buffer *buf)
{
    if (!buf->len) return;

    ARRAY (QnameIndexEnt, arr, *buf);

    // sort-uniq depn index    
    qsort (arr, arr_len, sizeof(QnameIndexEnt), qname_index_sorter);

    // unique
    uint32_t new_len=1;
    
    for (uint32_t i=1; i < arr_len; i++)
        if (arr[i].hash != arr[i-1].hash) 
            arr[new_len++] = arr[i];

        else if (arr[i].vb_i != arr[new_len-1].vb_i)
            arr[new_len-1].vb_i = 0; // qname hash exists in more than 1 vb
    
    buf->len32 = new_len;
}

static void scan_index_qnames_preprocessing (VBlockP vb)
{
    START_TIMER; 

    // if the txt file is compressed with BGZF, we uncompress now, in the compute thread
    if (txt_file->codec == CODEC_BGZF) 
        bgzf_uncompress_vb (vb);    // some of the blocks might already have been decompressed while reading - we decompress the remaining

    rom next  = B1STtxt;
    rom after = BAFTtxt;

    while (next < after)
        next = scan_index_one_line (VB_SAM, next, after - next);

    scan_sort_unique_depn_index (&vb_depn_index);
    scan_sort_unique_qname_index (&VB_SAM->qname_count);

    // tell dispatcher this thread is done and can be joined.
    vb_set_is_processed (vb); 
    
    COPY_TIMER (scan_index_qnames_preprocessing);
}

static void scan_append_index (VBlockP vb)
{
    buf_add_buf (real_evb, &z_file->sag_depn_index, &vb_depn_index, uint32_t, "z_file->sag_depn_index");
    buf_add_buf (evb, &z_qname_index, &VB_SAM->qname_count, QnameIndexEnt, "z_qname_index");
}

// remove depn for which qnames_hash instances (prim and depn) are contained in a single VB
static void scan_remove_single_vb_depns (void)
{
    if (!z_file->sag_depn_index.len) return;

    ARRAY (uint32_t, depn, z_file->sag_depn_index);
    ARRAY (QnameIndexEnt, qname, z_qname_index);

    uint64_t d=0, new_depn_len=0;
    for (uint64_t q=0; q < qname_len; q++) {
        if (qname[q].hash == depn[d]) {
            if (qname[q].vb_i == 0) // appears in more than 1 VB - we'll keep it
                depn[new_depn_len++] = depn[d];
            d++;
        } 
    }

    z_file->sag_depn_index.len = new_depn_len;
}

// ZIP: runs during segconf: scan the entire file and build z_file->sag_depn_index - a uniq-sorted list of hash(QNAME) of all depn
// alignments in the file
void sam_sag_by_flag_scan_for_depn (void)
{
    const rom task_name = "scan_for_depn";

    VBlockP scan_vb = vb_initialize_nonpool_vb (VB_ID_DEPN_SCAN, DT_BAM, task_name);
    VBlockP real_evb = evb;
    evb = scan_vb;

    File *save_txt_file = txt_file;
    txt_file = file_open (save_txt_file->name, READ, TXT_FILE, DT_BAM);
    
    txtfile_read_header (true); // reads into evb->txt_data and evb->lines.len
    buf_free (evb->txt_data);   // discard the header
    
    dispatcher_fan_out_task (task_name, txt_file->basename, 
                             0, "Preprocessing...", // this is not the same as the preprocessing that happens in PIZ - loading sags'
                             false, false, flag.xthreads, 0, 5000,
                             scan_read_one_vb, 
                             scan_index_qnames_preprocessing, 
                             scan_append_index);

    // most of the sort/uniq work was already done by the compute threads, so this final pass should be fast
    scan_sort_unique_depn_index (&z_file->sag_depn_index);
    scan_sort_unique_qname_index (&z_qname_index);

    // keep only depns that have a qname that appears in 2 or VBs. If all qnames are in a single VB, we will use saggy.
    scan_remove_single_vb_depns ();

    file_close (&txt_file, false, true);

    evb = real_evb;
    vb_destroy_vb (&scan_vb);
    
    txt_file = save_txt_file;
}

// for SAG_BY_FLAG: returns true if qname_hash represets a sag that has alignments in more than one VB
static bool sam_sag_is_multi_vb_SAG_BY_FLAG (int64_t first, int64_t last, uint32_t qname_hash)
{
    if (last < first) return false; // qname_hash doesn't exist in the array

    int64_t mid = (first + last) / 2;
    uint32_t mid_hash = *B32 (z_file->sag_depn_index, mid);

    if (mid_hash == qname_hash) return true; // found qname_hash in the array

    if (mid_hash < qname_hash) return sam_sag_is_multi_vb_SAG_BY_FLAG (mid+1, last, qname_hash);
    else                       return sam_sag_is_multi_vb_SAG_BY_FLAG (first, mid-1, qname_hash);
}

//---------------------------------------------------------------------------
// Scanning during sam_seg_initialize - for non-SAG_by_FLAG files - so we can 
// know whether all alignments of a sag are contained in a VB (in which case
// no need to move to the prim/depn components) - BAM only
//---------------------------------------------------------------------------

static ASCENDING_SORTER (qname_count_sorter, QnameCount, hash)

void scan_index_qnames_seg (VBlockSAMP vb)
{
    rom next  = B1STtxt;
    rom after = BAFTtxt;

    while (next < after)
        next = scan_index_one_line (vb, next, after - next);

    ARRAY (QnameCount, counts, vb->qname_count);

    if (!counts_len) return; // all lines are unmapped
    
    // sort-uniq counts    
    qsort (counts, counts_len, sizeof(QnameCount), qname_count_sorter);

    // unique
    uint32_t new_len=1;
    for (uint32_t i=1; i < counts_len; i++)
        if (counts[i].hash != counts[i-1].hash) 
            counts[new_len++] = counts[i];

        else 
            counts[new_len-1].count++; 

    // keep only those with 2 or more counts, as NH=1 is not a candidate for gencomp anyway
    uint32_t new_len2=0;
    
    for (uint32_t i=0; i < new_len; i++)
        if (counts[i].count >= 2) 
            counts[new_len2++] = counts[i];
    
    vb->qname_count.len32 = new_len2;
}

// returns the number of alignments in this VB with qname_hash
static uint32_t sam_sag_get_qname_counts_this_vb (QnameCount *counts, uint32_t counts_len, int32_t first, int32_t last, uint32_t qname_hash)
{
    if (last < first) return 0; // qname_hash doesn't exist in the array

    int32_t mid = (first + last) / 2;
    uint32_t mid_hash = counts[mid].hash;

    if (mid_hash == qname_hash) return counts[mid].count; // found qname_hash in the array

    if (mid_hash < qname_hash) return sam_sag_get_qname_counts_this_vb (counts, counts_len, mid+1, last, qname_hash);
    else                       return sam_sag_get_qname_counts_this_vb (counts, counts_len, first, mid-1, qname_hash);
}

//---------------------------------------------------
// Shared - called from sam_seg_is_gc_line during seg
//---------------------------------------------------

// ZIP: true if this prim or depn sag line *might* have saggies in other VBs
bool sam_might_have_saggies_in_other_VBs (VBlockSAMP vb, ZipDataLineSAM *dl, int32_t n_alns)
{
    if (!vb->qname_count.len32) return true; // we didn't count qnames, so we don't have proof that there aren't any saggies in other VBs

    uint32_t qname_hash = QNAME_HASH (Btxt(dl->QNAME.index), dl->QNAME.len, dl->FLAG.is_last);

    if (IS_SAG_CC) 
        return true;   // we don't do qname accounting for SAG_BY_CC

    else if (IS_SAG_FLAG)
        return sam_sag_is_multi_vb_SAG_BY_FLAG (0, z_file->sag_depn_index.len - 1, qname_hash);

    else {
        if (!n_alns) n_alns = str_count_char (STRauxZ (SA_Z, true), ';'); // happens iff depn of SAG_BY_SA

        return !vb->qname_count.len32 &&
               n_alns != sam_sag_get_qname_counts_this_vb (B1ST (QnameCount, vb->qname_count), vb->qname_count.len32, 0, vb->qname_count.len - 1, qname_hash);
    }
}
