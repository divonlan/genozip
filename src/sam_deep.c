// ------------------------------------------------------------------
//   deep_sam.c
//   Copyright (C) 2023-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "deep.h"
#include "bits.h"
#include "piz.h"
#include "htscodecs/rANS_static4x16.h"

// --------
// ZIP side
// --------

// Called during zip_finalize of the SAM component for a Deep compression
void sam_deep_zip_finalize (void)
{
    // Build vb_start_deep_line: first deep_line (0-based SAM-wide line_i, but not counting SUPP/SEC lines, monochar reads, SEQ.len=0) of each SAM VB
    buf_alloc (evb, &z_file->vb_start_deep_line, 0, z_file->num_vbs, uint64_t, 0, "z_file->vb_start_deep_line");

    uint64_t count=0;
    for (VBIType vb_i=1; vb_i <= z_file->num_vbs; vb_i++) {
        BNXT64(z_file->vb_start_deep_line) = count;
        count += sections_vb_header (vb_i, HARD_FAIL)->num_deep_lines;
    }
}

// SAM seg: calculate hash values into dl->deep_*

#define MAX_ENTRIES 0xfffffffe // our current implementation is limited to 4G reads

// ZIP compute thread while segging SEQ
void sam_deep_set_QNAME_hash (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qname))
{
    if (IS_DEEPABLE(dl->FLAG)) {
        dl->deep_hash.qname = deep_qname_hash (STRa(qname));
        vb->lines.count++; // counts primary lines
    }
}

// ZIP compute thread while segging SEQ: set hash of forward SEQ (i.e. as it would appear in the FASTQ file) 
void sam_deep_set_SEQ_hash (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(textual_seq))
{
    if (!IS_DEEPABLE(dl->FLAG) || (!dl->SEQ.len && textual_seq_len <= 1)) return; // note: SEQ.len=0 also for unmapped reads, but we do want to include them

    // we don't Deep monochar sequences as there is too much contention
    if (str_is_monochar (STRa(textual_seq))) 
        dl->monochar_seq = true;

    else
        dl->deep_hash.seq = deep_seq_hash (VB, STRa(textual_seq), dl->FLAG.rev_comp);
}

// ZIP compute thread while segging QUAL: hash of forward QUAL (i.e. as it would appear in the FASTQ file) for Deep 
void sam_deep_set_QUAL_hash (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qual))
{
    if (!IS_DEEPABLE(dl->FLAG) || !vb->has_qual || segconf.has_bqsr || dl->monochar_seq) return; // dl->deep_qual_hash remains 0

    dl->deep_hash.qual = deep_qual_hash (VB, STRa(qual), dl->FLAG.rev_comp); 

    if (flag.debug_deep == 2 && deephash_issame (dl->deep_hash, flag.debug_deep_hash)) 
        printf ("%s Found deep_hash=%u,%u,%u\nQNAME=\"%.*s\"\nSEQ=\"%.*s\"\nQUAL=\"%.*s\"\n",
                LN_NAME, DEEPHASHf(dl->deep_hash), STRfw(dl->QNAME), STRfb(vb->textual_seq), STRfw(dl->QUAL));
}

// ZIP compute thread: callback from ctx_merge_in_vb_ctx during merge: add VB's deep_hash to z_file->deep_hash/deep_ents
void sam_deep_merge (VBlockP vb_)
{
    if (!flag.deep) return;

    VBlockSAMP vb = (VBlockSAMP)vb_;

    START_TIMER;

    // initialize
    if (!z_file->deep_hash.len) {
        double est_num_vbs = MAX_(1, (double)txtfile_get_seggable_size() / (double)vb->txt_data.len * 1.05/* additonal PRIM VBs */);
        uint64_t est_num_lines = MIN_(MAX_ENTRIES, (uint64_t)(est_num_vbs * (double)vb->lines.len * 1.15)); // inflate the estimate by 15% - to reduce hash contention, and to reduce realloc chances for z_file->deep_ents

        // number of bits of hash table size (maximum 31)
        int clzll = (int)__builtin_clzll (est_num_lines * 2); 
        num_hash_bits = MIN_(31, 64 - clzll); // num hash bits

        z_file->deep_hash.can_be_big = true;
        buf_alloc_exact_255 (evb, z_file->deep_hash, ((uint64_t)1 << num_hash_bits), uint32_t, "z_file->deep_hash"); // pointer into deep_ents
        
        z_file->deep_ents.can_be_big = true;
        buf_alloc_exact_255 (evb, z_file->deep_ents, MAX_(est_num_lines, vb->lines.count), ZipZDeep, "z_file->deep_ents");

        if (flag.debug_deep) 
            printf ("num_hash_bits=%u est_num_lines*1.15=%"PRIu64" clzll=%d deep_hash.len=%"PRIu64" deep_ents.len=%"PRIu64"\n", 
                    num_hash_bits, est_num_lines, clzll, z_file->deep_hash.len, z_file->deep_ents.len);

        z_file->deep_ents.len = 0;
    }

    buf_alloc_255 (evb, &z_file->deep_ents, vb->lines.count/*num primary lines in VB*/, 0, ZipZDeep, CTX_GROWTH, "z_file->deep_ents"); 

    uint32_t ent_i = z_file->deep_ents.len32;

    ARRAY (uint32_t, hash_table, z_file->deep_hash);
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    uint32_t mask = bitmask32 (num_hash_bits);
printf ("xxx mask=%u\n", mask);
    uint32_t deep_line_i = 0; // line_i within the VB, but counting only primary lines that have SEQ.len > 0 and are not monochar
    for_buf2 (ZipDataLineSAM, dl, line_i, vb->lines) {
        if (!IS_DEEPABLE (dl->FLAG) || !dl->SEQ.len || dl->monochar_seq) continue; // case: secondary or supplementary line, or no sequence

        // we can only handle 4G entries, because the hash table as well as num_hash_bits is limited 32 bit.
        if (ent_i == 0xffffffff) {
            WARN_ONCE ("The number of alignments in %s exceeds the maximum supported for optimal Deep compression. The compression will not be as good as it could have. Please report this to " EMAIL_SUPPORT, txt_name);
            break;
        }

        uint32_t hash = dl->deep_hash.seq & mask; 
        
        // case: first entry for this hash value
        if (hash_table[hash] == NO_NEXT) 
            hash_table[hash] = ent_i; // first entry on linked list

        // case: add entry add end of linked list
        else {
            uint32_t linked_ent_i = hash_table[hash];
            while (deep_ents[linked_ent_i].next != NO_NEXT) {
if (deephash_issame(deep_ents[linked_ent_i].hash, dl->deep_hash)) printf ("xxx11 dup ent: old_hash=%u,%u,%u old vb,line=%u,%u new_hash=%u,%u,%u new vb,line=%u,%u\n", 
DEEPHASHf(deep_ents[linked_ent_i].hash), deep_ents[linked_ent_i].vb_i, deep_ents[linked_ent_i].line_i , DEEPHASHf(dl->deep_hash), vb->vblock_i, line_i);                
                linked_ent_i = deep_ents[linked_ent_i].next;
            }

if (deephash_issame(deep_ents[linked_ent_i].hash, dl->deep_hash)) printf ("xxx22 dup ent: old_hash=%u,%u,%u old vb,line=%u,%u new_hash=%u,%u,%u new vb,line=%u,%u\n", 
DEEPHASHf(deep_ents[linked_ent_i].hash), deep_ents[linked_ent_i].vb_i, deep_ents[linked_ent_i].line_i , DEEPHASHf(dl->deep_hash), vb->vblock_i, line_i);                

            deep_ents[linked_ent_i].next = ent_i; // extend linked list
        }

        // create new entry
        deep_ents[ent_i++] = (ZipZDeep){ .next   = NO_NEXT, 
                                         .hash   = dl->deep_hash,
                                         .vb_i   = vb->vblock_i,
                                         .line_i = deep_line_i++ };
    }

    z_file->deep_ents.len32 = ent_i;

    COPY_TIMER(sam_deep_merge);
} 

//-----------------------------------------------------------------------------
// SAM PIZ: step 1: output QNAME,SEQ,QUAL into vb->deep_ents / vb->deep_index
//-----------------------------------------------------------------------------

void sam_piz_deep_init_vb (VBlockSAMP vb, const SectionHeaderVbHeader *header)
{
    vb->rans_compress_bound_longest_seq_len = rans_compress_bound_4x16 (vb->longest_seq_len, X_NOSZ);
}

// strive to alloc exactly twice: 128K initialy, and then a good guess of the total
static void inline sam_piz_alloc_deep_ents (VBlockSAMP vb, uint32_t size)
{
    // case: initial allocation
    if (!vb->deep_ents.len) {
        buf_alloc (vb, &vb->deep_ents, size, 128 KB, uint8_t, 0, "deep_ents"); 
        buf_alloc (vb, &vb->deep_index, 0, vb->lines.len, uint32_t, 0, "deep_index"); 
    }

    // case: have enough memory, no need to allocate
    else if (vb->deep_ents.len + size <= vb->deep_ents.size)
        return;

    // case: we already have stats from a few lines - allocate now for the entire VB
    else {
        uint32_t avg_per_line_so_far = vb->deep_index.len ? (vb->deep_ents.len32 / vb->deep_index.len32) : 0;
        buf_alloc (vb, &vb->deep_ents, size, avg_per_line_so_far * vb->lines.len, uint8_t, 1.1, "deep_ents"); 
    }
}

// copy QNAME as is to deep_ents
static void sam_piz_deep_add_qname (VBlockSAMP vb, STRp(qname))
{
    sam_piz_alloc_deep_ents (vb, 2 + qname_len);

    BNXT32 (vb->deep_index) = vb->deep_ents.len32; // index in deep_ents of data of current deepable line
    vb->deep_ents.len32++; // // skip the PizZDeepFlags byte

    if (!segconf.deep_no_qname) {
        uint8_t *next = BAFT8 (vb->deep_ents);

        ASSPIZ (qname_len <= 254, "QNAME.len=%u, but per BAM specfication it is limited to 254 characters. QNAME=\"%.*s\"", qname_len, STRf(qname));
        *next = qname_len;

        memcpy (next + 1, qname, qname_len);

        vb->deep_ents.len32 += qname_len + 1;
    }
}

// compress QUAL or SEQ - adding a 1B or 4B length
static bool sam_piz_deep_compress (VBlockSAMP vb, STRp(data), bool is_seq)
{
    // reverse if needed
    if (last_flags.rev_comp) {
        ASSERTNOTINUSE (vb->scratch);
        buf_alloc_exact (vb, vb->scratch, data_len, char, "scratch");

        if (is_seq) str_revcomp (B1STc(vb->scratch), STRa(data));
        else        str_reverse (B1STc(vb->scratch), STRa(data));
    
        data = B1STc(vb->scratch);
    }

    sam_piz_alloc_deep_ents (vb, 4 + vb->rans_compress_bound_longest_seq_len);

    uint8_t *next = BAFT8 (vb->deep_ents);

    uint32_t comp_len = vb->rans_compress_bound_longest_seq_len;
    
    bool len_is_1_byte = (data_len < (is_seq ? 1024 : 512)); // just a guess for now
    uint8_t *guess_addr = next + (len_is_1_byte ?  1 : 4);   // we don't know comp_len yet, so just a guess
    ASSPIZ (rans_compress_to_4x16 (evb, (uint8_t *)STRa(data), guess_addr, &comp_len, X_NOSZ) && comp_len,
            "%s: Failed to compress data: data_len=%u", LN_NAME, data_len);

    // we guessed wrong that we need 1 byte for the length - actually need 4 
    if (len_is_1_byte && comp_len > 255) {
        memmove (guess_addr+3, guess_addr, comp_len); // make room for another 3 bytes
        len_is_1_byte = false;
    }

    if (len_is_1_byte)  // 1B length
        *next = comp_len;

    else  // 4B length
        memcpy (next, &comp_len, 4); // memcpy bc non-aligned

    vb->deep_ents.len32 += comp_len + (len_is_1_byte ? 1 : 4);

    if (last_flags.rev_comp)   
        buf_free (vb->scratch);

    return !len_is_1_byte; // true long length
}

// pack SEQ into 2-bit and into deep_ents or compress if it has non-ACGT (eg N)
void sam_piz_deep_add_seq (VBlockSAMP vb, STRp(seq))
{
    vb->seq_is_monochar = str_is_monochar (STRa(seq));

    // case: not deepable bc no SEQ or monochar - rollback
    if (vb->seq_missing || vb->seq_is_monochar) {
        vb->deep_ents.len32 = *BLST32 (vb->deep_index);
        vb->deep_index.len32--;
        return;
    }

    PizZDeepFlags *deep_flags = B(PizZDeepFlags, vb->deep_ents, *BLST32(vb->deep_index));
    *deep_flags = (PizZDeepFlags){}; // initialize

    if (seq_len <= 255) 
        BNXT8 (vb->deep_ents) = seq_len;
    
    else {
        deep_flags->is_long_seq = true;
        memcpy (BAFT8 (vb->deep_ents), &seq_len, sizeof (uint32_t)); // unaligned
        vb->deep_ents.len32 += sizeof (uint32_t);
    }

    // if only A,C,G,T - pack in 2-bits
    if (str_is_only_ACGT (STRa(seq), NULL)) {
        sam_piz_alloc_deep_ents (vb, 4 + seq_len + 2*sizeof(uint64_t)/*padding for Bits*/); 

        uint32_t pack_index = ROUNDUP8 (vb->deep_ents.len32); // 64b word-align as required by Bits
        Bits pack = bits_init (seq_len * 2, B8(vb->deep_ents, pack_index), vb->deep_ents.size - pack_index, false);
        sam_seq_pack (vb, &pack, 0, STRa(seq), false, last_flags.rev_comp, HARD_FAIL);

        uint32_t bytes = roundup_bits2bytes (pack.nbits);

        if (BAFT8 (vb->deep_ents) != (uint8_t*)pack.words)
            memmove (BAFT8 (vb->deep_ents), pack.words, bytes); // move back to (unaligned) place

        vb->deep_ents.len32 += bytes;
    }

    // if has an N or other IUPACs - compress
    else {
        deep_flags->is_seq_compressed = true;
     
        deep_flags->is_long_seq_comp = sam_piz_deep_compress (vb, STRa(seq), true);
    }
}

// compress QUAL into deep_ents
void sam_piz_deep_add_qual (VBlockSAMP vb, STRp(qual))
{
    if (vb->seq_missing || vb->seq_is_monochar || VB_SAM->qual_missing) return; // not a deepable alignment or no QUAL

    PizZDeepFlags *deep_flags = B(PizZDeepFlags, vb->deep_ents, *BLST32(vb->deep_index));

    deep_flags->is_long_qual_comp = sam_piz_deep_compress (vb, STRa(qual), false);
}


// called for top-level SAM (but not BAM) items in case file was compressed with --deep
CONTAINER_ITEM_CALLBACK (sam_piz_con_item_cb)
{
    if (!IS_DEEPABLE (last_flags)) return;

    switch (con_item->did_i_small) {
        case SAM_QNAME    : 
        case SAM_QNAMESA  : sam_piz_deep_add_qname (VB_SAM, STRa(recon)); break;
        case SAM_SQBITMAP : sam_piz_deep_add_seq   (VB_SAM, STRa(recon)); break;
        case SAM_QUAL     : // note: QUAL is in toplevel in MAIN QUALSA is toplevel in PRIM. DEPN VBs don't have this callback.
        case SAM_QUALSA   : sam_piz_deep_add_qual  (VB_SAM, STRa(recon)); break; 
        default           : ASSPIZ (false, "Invalid con_item=%s", dis_dict_id (con_item->dict_id).s);
    }
}

// called for QNAME in BAM (not SAM)
TRANSLATOR_FUNC (sam_piz_sam2bam_QNAME)
{
    if (flag.deep && IS_DEEPABLE (last_flags))
        sam_piz_deep_add_qname (VB_SAM, STRa(recon));

    return 0;
}

//-----------------------------------------------------------------------------
// SAM PIZ: step 2: move VB deep* buffers to z_file
//-----------------------------------------------------------------------------

// PIZ: called by the main thread for each reconstructed VB, in the order of VBs
void sam_piz_deep_grab_deep_ents (VBlockSAMP vb)
{
    static VBIType num_sam_vbs=0; // updated when pizzing vb_i==1 for each z_file

    if (vb->vblock_i == 1) {
        num_sam_vbs = sections_get_num_vbs (SAM_COMP_MAIN) + sections_get_num_vbs (SAM_COMP_PRIM) + sections_get_num_vbs (SAM_COMP_DEPN);

        buf_alloc_exact_zero (evb, z_file->vb_start_deep_line, num_sam_vbs, uint64_t, "z_file->vb_start_deep_line");
        buf_alloc_exact_zero (evb, z_file->deep_ents,          num_sam_vbs, Buffer,   "z_file->deep_ents"); 
        buf_alloc_exact_zero (evb, z_file->deep_index,         num_sam_vbs, Buffer,   "z_file->deep_index");        
    }

    // set start for next VB - after lines of current VB (start for vb_i=1 remains 0)
    if (vb->vblock_i < num_sam_vbs) 
        VB_start_deep_line (vb->vblock_i + 1) = VB_start_deep_line(vb->vblock_i) + vb->deep_index.len; 

    buf_grab (evb, *B(Buffer, z_file->deep_ents,  vb->vblock_i-1), "z_file->deep_ents",  vb->deep_ents);
    buf_grab (evb, *B(Buffer, z_file->deep_index, vb->vblock_i-1), "z_file->deep_index", vb->deep_index);    
}
