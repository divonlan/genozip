// ------------------------------------------------------------------
//   sam_sag_scan.c
//   Copyright (C) 2022-2026 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

// In BAM files with SAG_BY_FLAG, absent an SA, NH or CC tag, we are unable to distiguish one of the alignments of
// being primary during, and unable to tell whether a non-secondary, non-supplementary alignment has dependent alignments or not.
// 
// This module scans the txt_file for depn lines, creating a sort-uniqed array of hash(QNAME) in z_file->sag_depn_index
// for each dependent line - during segconf, and a function for determining if a depn exists by looking up in this array during seg

#include "sam_private.h"
#include "dispatcher.h"
#include "txtfile.h"
#include "qname.h"
#include "progress.h"
#include "mgzip.h"
#include "sorter.h"

static VBlockP real_evb;
static VBlockP scan_vb;

#define vb_depn_index  vb->scratch // QNAME hash of all sec and supp alignments in this VB
#define z_qname_index  evb->z_data // QNAME hash of all mapped alignments in this file

typedef struct { uint32_t hash; VBIType vb_i; } QnameIndexEnt;
static rom scan_index_one_line (VBlockSAMP vb, rom alignment, uint32_t remaining_txt_len)   
{
    STR(qname);
    SamFlags flag;
    uint32_t alignment_len;

    if (IS_BAM_ZIP) {
        alignment_len = GET_UINT32 (alignment) + 4;

        // a non-sensical block_size might indicate an false-positive identification of a BAM alignment in bam_unconsumed
        ASSERT (IN_RANGX (alignment_len, sizeof (BAMAlignmentFixed), remaining_txt_len), 
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

    uint32_t hash = qname_calc_hash (QNAME1, COMP_NONE, qname, qname_len, flag.is_last, false, CRC32, NULL);

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

    if (Ltxt) {
        vb->preprocessing = true;
        vb->dispatch = READY_TO_COMPUTE;
    }

    else 
        vb->dispatch = DATA_EXHAUSTED;
}

static SORTER (depn_index_sorter)
{
    return ASCENDING_RAW (*(uint32_t *)a, *(uint32_t *)b);
}

// sort and uniq list of qname hashes that appear in this VB
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

// sort and uniq QNAME hashes, and also mark if they appear in more than one VB
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
    if (TXT_IS_BGZF) 
        mgzip_uncompress_vb (vb, CODEC_BGZF);    // some of the blocks might already have been decompressed while reading - we decompress the remaining

    rom next  = B1STtxt;
    rom after = BAFTtxt;

    #define SET_HAS_SA \
        ({ segconf.has[OPTION_SA_Z] = true; /* likely, bot not 100% sure, that an SAZ string indicates an SA:Z optional field */ \
           DO_ONCE if (flag.show_scan) \
               iprintf ("Scan: SA:Z found in vb=%u", vb->vblock_i); })

    // in case of force-gencomp, we also try to detect SA_Z in the first 16 VBs
    if (flag.force_gencomp && !segconf.has[OPTION_SA_Z]) {
        if (IS_BAM_ZIP) {
            for (rom c=next; c < after-4 ; c++) 
                if (c[0] == 'S' && c[1]== 'A' && c[2] == 'Z' && IS_ALPHANUMERIC(c[3])
                    && sam_zip_is_valid_SA (&c[3], after - &c[3], true)) {
                    SET_HAS_SA;
                    break;
                }
        }
        else {
            SAFE_NUL (after);
            rom sa = strstr (next, "SA:Z:");
            SAFE_RESTORE;

            if (sa && sam_zip_is_valid_SA (sa + 5, after - (sa+5), false)) 
                SET_HAS_SA;
        }
    }

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
    buf_append_buf (real_evb, &z_file->sag_depn_index, &vb_depn_index, uint32_t, NULL);
    buf_append_buf (evb, &z_qname_index, &VB_SAM->qname_count, QnameIndexEnt, "z_qname_index");
}

// remove depn for which qnames_hash instances (prim and depn) are contained in a single VB
static void scan_remove_single_vb_depns (void)
{
    START_TIMER;

    if (!z_file->sag_depn_index.len) return;

    ARRAY (uint32_t, depn, z_file->sag_depn_index);

    uint64_t d=0, new_depn_len=0;
    for_buf (QnameIndexEnt, q, z_qname_index)
        if (q->hash == depn[d]) {
            if (q->vb_i == 0) // appears in more than 1 VB - we'll keep it
                depn[new_depn_len++] = depn[d];
            d++;
        } 

    z_file->sag_depn_index.len = new_depn_len;

    if (flag.show_scan) {
        iprintf ("Scan: Unique mapped QNAMExfirst in file: %"PRIu64"\n", z_qname_index.len);
        iprintf ("Scan: Unique mapped QNAMExfirst that appear in multiple VBs: %"PRIu64"\n", new_depn_len);
    }

    COPY_TIMER_EVB (scan_remove_single_vb_depns);
}

// ZIP: runs during segconf if SAG_BY_FLAG is suspected: scan the entire file and build 
// z_file->sag_depn_index - a uniq-sorted list of hash(QNAME) of all depn alignments in the file
void sam_sag_by_flag_scan_for_depn (void)
{
    START_TIMER;

    const rom task_name = "scan_for_depn";

    buf_alloc (evb, &z_file->sag_depn_index, 0, 100, uint32_t, 0, "z_file->sag_depn_index");

    scan_vb = vb_initialize_nonpool_vb (VB_ID_SCAN_VB, DT_BAM, task_name);
    VBlockP real_evb = evb;
    evb = scan_vb;

    FileP save_txt_file = txt_file;
    txt_file = file_open_txt_read (save_txt_file->name);
    ASSERTNOTNULL (txt_file);
    
    TEMP_FLAG(biopsy, false);
    SAVE_FLAG(biopsy_line); flag.biopsy_line.line_i = NO_LINE;
    txtfile_read_header (true); // reads into evb->txt_data and evb->lines.len
    buf_free (evb->txt_data);   // discard the header
    
    dispatcher_fan_out_task (task_name, 
                             save_txt_file->basename, // not txt_file-> bc it will be closed in a sec, while the progress component will continue to the main zip fan_out 
                             0, "Preprocessing...",   // this is not the same as the preprocessing that happens in PIZ - loading sags'
                             false, false, flag.xthreads, 0, 5000, true,
                             scan_read_one_vb, 
                             scan_index_qnames_preprocessing, 
                             scan_append_index);

    RESTORE_FLAG (biopsy);
    RESTORE_FLAG (biopsy_line);

    progress_erase(); // erase "Preprocessing..."
    
    // most of the sort/uniq work was already done by the compute threads, so this final pass should be fast
    // TO DO: replace sort with combine sorted arrays 
    { START_TIMER;  
    scan_sort_unique_depn_index (&z_file->sag_depn_index); 
    COPY_TIMER_EVB (sam_sag_by_flag_scan_sag_depn_index); }
    
    { START_TIMER;  
    scan_sort_unique_qname_index (&z_qname_index);   
    COPY_TIMER_EVB (sam_sag_by_flag_scan_sort_qname_index); }

    // keep only depns that have a qname that appears in 2 or VBs. If all qnames are in a single VB, we will use saggy.
    scan_remove_single_vb_depns();

    file_close (&txt_file);

    evb = real_evb;
    vb_destroy_vb (&scan_vb);
    
    txt_file = save_txt_file;

    COPY_TIMER_EVB (sam_sag_by_flag_scan_for_depn);
}

//---------------------------------------------------------------------------
// Scanning during sam_seg_initialize - for non-SAG_by_FLAG files - so we can 
// know whether all alignments of a sag are contained in a VB (in which case
// no need to move to the prim/depn components) 
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

//---------------------------------------------------
// Shared - called from sam_seg_is_gc_line during seg
//---------------------------------------------------

// for SAG_BY_FLAG: returns non-NULL if qname_hash represets a sag that has alignments in more than one VB
static BINARY_SEARCHER_INTEGRAL (sam_sag_is_multi_vb_SAG_BY_FLAG, uint32_t, true, IfNotExact_ReturnNULL)

static BINARY_SEARCHER (sam_sag_get_qname_counts_this_vb, QnameCount, uint32_t, hash, true, IfNotExact_ReturnNULL)

// ZIP: true if this prim or depn sag line *might* have saggies in other VBs
bool sam_might_have_saggies_in_other_VBs (VBlockSAMP vb, ZipDataLineSAMP dl, int32_t n_alns)
{
    if (!vb->qname_count.len32) return true; // we didn't count qnames, so we don't have proof that there aren't any saggies in other VBs

    uint32_t qname_hash = qname_calc_hash (QNAME1, COMP_NONE, STRqname(dl), dl->FLAG.is_last, false, CRC32, NULL);

    if (IS_SAG_CC) 
        return true;   // we don't do qname accounting for SAG_BY_CC

    else if (IS_SAG_FLAG) 
        return binary_search (sam_sag_is_multi_vb_SAG_BY_FLAG, uint32_t, z_file->sag_depn_index, qname_hash) != NULL;

    else {
        if (!n_alns) n_alns = str_count_char (STRauxZ (SA_Z, true), ';'); // n_alns=0 only possible in SAG_BY_SA
        QnameCount *count = binary_search (sam_sag_get_qname_counts_this_vb, QnameCount, vb->qname_count, qname_hash);

        return !count || n_alns != count->count;
    }
}
