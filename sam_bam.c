// ------------------------------------------------------------------
//   sam_bam.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "md5.h"
#include "buffer.h"
#include "vblock.h"
#include "txtfile.h"
#include "file.h"
#include "endianness.h"
#include "sam_private.h"
#include "seg.h"
#include "strings.h"
#include "random_access.h"
#include "dict_id.h"
#include "codec.h"

// loading a Little Endian uint32_t from an unaligned buffer
#define GET_UINT16(p) (((uint8_t*)p)[0] | (((uint8_t*)p)[1] << 8))
#define GET_UINT32(p) (((uint8_t*)p)[0] | (((uint8_t*)p)[1] << 8) | (((uint8_t*)p)[2] << 16) | (((uint8_t*)p)[3] << 24))

// getting integers from the BAM data
#define NEXT_UINT8  *((const uint8_t *)next_field++)
#define NEXT_UINT16 GET_UINT16 (next_field); next_field += sizeof (uint16_t);
#define NEXT_UINT32 GET_UINT32 (next_field); next_field += sizeof (uint32_t);

#define READ_UINT32(vb,u) \
    const char *u##_ = bam_read_exact(vb, 4, false); \
    uint32_t u = GET_UINT32 (u##_);

#define BAM_MAGIC "BAM\1" // first 4 characters of a BAM file

// copy header ref data during reading header
static uint32_t header_n_ref=0;
static char **header_ref_contigs = NULL;
static uint32_t *header_ref_contigs_len = NULL;

// read exactly num_bytes from the txt file, append txt_data
static char *bam_read_exact (VBlockP vb, uint32_t num_bytes, bool allow_eof)
{
    buf_alloc (vb, &vb->txt_data, vb->txt_data.len + num_bytes, 2, "txt_data", 0); // increase size if needed
    char *start = AFTERENT (char, vb->txt_data);

    int32_t bytes_read = txtfile_read_block (start, num_bytes);
    
    if (!bytes_read && allow_eof) return NULL;

    ASSERT (bytes_read == (int32_t)num_bytes, "Error in bam_read_exact: needed to read %u bytes, but read only %d",
            num_bytes, bytes_read);

    vb->txt_data.len += bytes_read;
    return start;
}

// ZIP I/O thread: returns the hash of the header
Md5Hash bam_read_txt_header (bool is_first_txt, bool header_required_unused, char first_char_unused)
{
    START_TIMER;

    const bool has_ref = (flag_reference == REF_EXTERNAL || flag_reference == REF_EXT_STORE);

    buf_alloc (evb, &evb->txt_data, 100000, 1, "txt_data", 0); // usally, this is more than enough, but bam_read_exact will enlarge if needed

    // interpret the BAM header by the spec on page 15 here: https://samtools.github.io/hts-specs/SAMv1.pdf  
    ASSERT (!memcmp (bam_read_exact (evb, 4, false), BAM_MAGIC, 4), // magic
            "Error in bam_read_txt_header: %s doesn't have a BAM magic - it doesn't seem to be a BAM file", txt_name);

    // read sam header (text) and count the lines
    READ_UINT32 (evb, l_text); 
    const char *text = bam_read_exact (evb, l_text, false); 

    for (int i=0; i < l_text; i++) 
        if (text[i] == '\n') evb->lines.len++;   

    READ_UINT32 (evb, n_ref);

    ASSERT (!has_ref || n_ref == ref_num_loaded_contigs(), "Error: contig mismatch: file %s has %u contigs in its header, but reference file %s has %u contigs",
            txt_name, n_ref, ref_filename, ref_num_loaded_contigs());
    
    // free allocation from previous file
    if (header_ref_contigs) {
        for (unsigned i=0; i < header_n_ref; i++) 
            FREE (header_ref_contigs[i]);    
        FREE (header_ref_contigs);
        FREE (header_ref_contigs_len);
    }
    // allocate storage for contigs for use in seg
    header_n_ref = n_ref;
    header_ref_contigs = MALLOC (n_ref * sizeof (char *));
    header_ref_contigs_len = MALLOC (n_ref * sizeof (uint32_t));

    for (uint32_t i=0; i < n_ref; i++) {
        READ_UINT32 (evb, l_name);
        const char *chrom_name = bam_read_exact (evb, l_name, false);
        READ_UINT32 (evb, l_ref);

        // check that the contigs specified in the header are consistent with the reference given in --reference/--REFERENCE
        if (has_ref) 
            ref_contigs_verify_identical_chrom (chrom_name, l_name-1, l_ref, i);
         
        // case: store contigs for use in seg. for vb_1 to enter into context
        header_ref_contigs[i] = MALLOC (l_name);
        header_ref_contigs_len[i] = l_name-1;
        memcpy (header_ref_contigs[i], chrom_name, l_name);
    }
    
    txt_file->txt_data_so_far_single = evb->txt_data.len;

    // md5 header - with logic related to is_first
    txtfile_update_md5 (evb->txt_data.data, evb->txt_data.len, !is_first_txt);
    Md5Hash header_md5 = md5_snapshot (&z_file->md5_ctx_single);

    COPY_TIMER_VB (evb, txtfile_read_header); // use same profiler id as standard

    return header_md5;
}

// PIZ I/O thread: make the txt header either SAM or BAM according to flag_reconstruct_binary, and regardless of the source file
void bam_prepare_txt_header (Buffer *txt)
{
    bool is_bam_header = (txt->len >= 4) && !memcmp (txt->data, BAM_MAGIC, 4);

    // case: we need to convert a BAM header to a SAM header
    if (is_bam_header && !flag_reconstruct_binary) {
        uint32_t l_text = GET_UINT32 (ENT (char, *txt, 4));
        memcpy (txt->data, ENT (char, *txt, 8), l_text);
        txt->len = l_text;
    }

    // case: we need to convert a SAM header to a BAM header
    else if (!is_bam_header && flag_reconstruct_binary) {

    }
}

// read vblock (instead of txtfile_read_vblock, since BAM is a binary file)
void bam_read_vblock (VBlock *vb)
{
    START_TIMER;

    int64_t pos_before = file_tell (txt_file);

    buf_alloc (vb, &vb->txt_data, global_max_memory_per_vb, 1, "txt_data", vb->vblock_i);    
    
    uint32_t block_size = txt_file->unconsumed_txt.param; // maybe passed from previous VB
    if (block_size) {
        *FIRSTENT (uint32_t, vb->txt_data) = LTEN32 (block_size); // encode back in its original Little Endian
        vb->txt_data.len += sizeof (uint32_t);
    }

    txt_file->unconsumed_txt.param = 0; // nothing yet to pass to the next_field VB
    
    while (vb->txt_data.len < global_max_memory_per_vb - sizeof (uint32_t) /* for block_size */) {  

        char *start = AFTERENT (char, vb->txt_data);

        bool block_size_is_read = false;
        if (!block_size) { // block size not already read in previous VB - read it now
            const uint8_t *p = (const uint8_t *)bam_read_exact(vb, 4, true); 
            if (!p) break; // EOF

            block_size = GET_UINT32 (p); // little endian, unaligned word
            block_size_is_read = true;

            ASSERT0 (block_size, "Error in bam_read_vblock: found block_size=0");
        }

        // case: we don't have enough space for this block - pass block_size to next_field VB
        if (global_max_memory_per_vb - vb->txt_data.len < block_size) {
            txt_file->unconsumed_txt.param = block_size; // pass on to next_field VB
            vb->txt_data.len -= sizeof (uint32_t);
            
            if (flag_md5) // snapshot at end of this VB's data - before adding the block_size
                vb->md5_hash_so_far = md5_snapshot (flag_bind ? &z_file->md5_ctx_bound : &z_file->md5_ctx_single);

            if (block_size_is_read) // add MD5 of little endian (i.e. as in the file) block_size
                txtfile_update_md5 (start, sizeof (uint32_t), false);

            break; 
        }

        (void)!bam_read_exact (vb, block_size, false);
            
        // note: we md_udpate after every block, rather on the complete data (vb or txt header) when its done
        // because this way the OS read buffers / disk cache get pre-filled in parallel to our md5
        // include the block_size itself if we read it in this VB, but not if it was passed from previous VB
        txtfile_update_md5 (start, AFTERENT (char, vb->txt_data) - start, false);

        block_size=0;
        vb->lines.len++;
    }

    // case: we completed full blocks are are not passing on a block_size to the next_field VB - we take the snapshot now
    if (!txt_file->unconsumed_txt.param && flag_md5) 
        vb->md5_hash_so_far = md5_snapshot (flag_bind ? &z_file->md5_ctx_bound : &z_file->md5_ctx_single);

    vb->vb_position_txt_file = txt_file->txt_data_so_far_single;

    txt_file->txt_data_so_far_single += vb->txt_data.len;
    vb->vb_data_size = vb->txt_data.len; // initial value. it may change if --optimize is used.
    
    vb->vb_data_read_size = file_tell (txt_file) - pos_before; // apporx. bgz compressed bytes read (measuring bytes uncompressed by zlib, but might still reside in zlib's output buffer)

    COPY_TIMER (txtfile_read_vblock); // use same profiler id as standard
}

void bam_seg_initialize (VBlock *vb)
{
    sam_seg_initialize (vb);

    // copy contigs from header in vb=1 for RNAME dictionary
    if (vb->vblock_i == 1) {
        for (unsigned i=0; i < header_n_ref; i++) 
            mtf_evaluate_snip_seg (vb, &vb->contexts[SAM_RNAME], header_ref_contigs[i], header_ref_contigs_len[i], NULL);

        vb->contexts[SAM_RNAME].inst |= CTX_INST_NO_VB1_SORT; // keep the contigs in the order as in BAM header
    }

}

static inline void bam_seg_bin (VBlock *vb, uint16_t bin)
{
    Context *ctx = &vb->contexts[BAM_BIN];
    char snip[20] = { SNIP_SELF_DELTA };

    PosType delta = (int64_t)bin - ctx->last_value.i;
    unsigned snip_len = 1 + str_int (delta, &snip[1]);
    ctx->flags |= CTX_FL_STORE_INT;
    ctx->last_value.i = (int64_t)bin;

    seg_by_ctx (vb, snip, snip_len, ctx, sizeof (uint16_t), NULL);
}

static inline void bam_seg_ref_id (VBlockP vb, DidIType did_i, int32_t ref_id, int32_t compare_to_ref_i)
{
    ASSERT (ref_id >= -1 && ref_id < (int32_t)header_n_ref, "Error in bam_seg_ref_id: vb=%u: encountered ref_id=%d but header has only %u contigs",
            vb->vblock_i, ref_id, header_n_ref);

    // get snip and snip_len
    DECLARE_SNIP;
    if (ref_id >= 0) {
        if (ref_id == compare_to_ref_i) {
            snip = "=";
            snip_len = 1;
        }
        else {
            snip = header_ref_contigs[ref_id];
            snip_len = header_ref_contigs_len[ref_id];
        }
    }
    else { 
        snip = "*";
        snip_len = 1;
    }

    WordIndex chrom_index = seg_by_did_i (vb, snip, snip_len, did_i, sizeof (int32_t));
        
    if (did_i==CHROM) 
        random_access_update_chrom (vb, chrom_index, snip, snip_len);
}

static inline void bam_rewrite_cigar (VBlockSAM *vb, uint16_t n_cigar_op, const uint32_t *cigar)
{
    if (!n_cigar_op) {
        buf_alloc (vb, &vb->textual_cigar, 2, 2, "textual_cigar", 0);
        NEXTENT (char, vb->textual_cigar) = '*';
        goto finish;
    }

    // calculate length
    unsigned len=0;
    for (uint16_t i=0; i < n_cigar_op; i++) {
        uint32_t op_len = LTEN32 (cigar[i]) >> 4;
        if      (op_len < 10)      len += 2; // 1 for the op, 1 for the number
        else if (op_len < 100)     len += 3; 
        else if (op_len < 1000)    len += 4; 
        else if (op_len < 10000)   len += 5; 
        else if (op_len < 100000)  len += 6; 
        else if (op_len < 1000000) len += 7;
        else ABORT ("op_len=%u too long in vb=%u: ", vb->vblock_i, op_len); 
    }

    buf_alloc (vb, &vb->textual_cigar, len + 1 /* for \0 */, 2, "textual_cigar", 0);

    for (uint16_t i=0; i < n_cigar_op; i++) {
        uint32_t subcigar = LTEN32 (cigar[i]);
        uint32_t op_len = subcigar >> 4;

        vb->textual_cigar.len += str_int (op_len, AFTERENT (char, vb->textual_cigar));

        static const char op_table[16] = "MIDNSHP=X";
        NEXTENT (char, vb->textual_cigar) = op_table[subcigar & 0xf];
    }

finish:
    *AFTERENT (char, vb->textual_cigar) = 0; // null terminate
    vb->last_cigar = FIRSTENT (char, vb->textual_cigar);
}

static void bam_seg_cigar_field (VBlockSAM *vb, ZipDataLineSAM *dl, uint32_t n_cigar_op)
{
    char cigar_snip[vb->textual_cigar.len + 20];
    cigar_snip[0] = SNIP_SPECIAL;
    cigar_snip[1] = SAM_SPECIAL_CIGAR;
    unsigned cigar_snip_len=2;

    // case: SEQ is "*" - we add a '-' to the CIGAR
    if (!dl->seq_len) cigar_snip[cigar_snip_len++] = '-';

    // case: CIGAR is "*" - we add the length to CIGAR eg "151*"
    if (*vb->last_cigar == '*')  // CIGAR is not available
        cigar_snip_len += str_int (MAX (dl->seq_len, 1), &cigar_snip[cigar_snip_len]); // if seq_len=0, then we add "1" because we have the * in seq and qual (and consistency with sam_seg_cigar_field)

    memcpy (&cigar_snip[cigar_snip_len], vb->textual_cigar.data, vb->textual_cigar.len);
    cigar_snip_len += vb->textual_cigar.len;

    seg_by_did_i (vb, cigar_snip, cigar_snip_len, SAM_CIGAR, 
                  n_cigar_op * sizeof (uint32_t) /* cigar */ + sizeof (uint16_t) /* n_cigar_op */ );
}

static inline void bam_rewrite_seq (VBlockSAM *vb, ZipDataLineSAM *dl, const char *next_field)
{
    buf_alloc (vb, &vb->textual_seq, dl->seq_len+1 /* for last half-byte */, 1.5, "textual_seq", 0);

    if (!dl->seq_len) {
        NEXTENT (char, vb->textual_seq) = '*';
        return;        
    }

    char *next = FIRSTENT (char, vb->textual_seq);

    for (uint32_t i=0; i < (dl->seq_len+1) / 2; i++) {
        static const char base_codes[16] = "=ACMGRSVTWYHKDBN";

        *next++ = base_codes[*(uint8_t*)next_field >> 4];
        *next++ = base_codes[*(uint8_t*)next_field & 0xf];
        next_field++;
    }

    vb->textual_seq.len = dl->seq_len;

    ASSERTW (!(dl->seq_len % 2) || (*AFTERENT(char, vb->textual_seq)=='='), 
            "Warning in bam_rewrite_seq vb=%u: expecting the unused lower 4 bits of last seq byte in an odd-length seq_len=%u to be 0, but its not. This will cause an incorrect MD5",
             vb->vblock_i, dl->seq_len);
}

// Rewrite the QUAL field - add +33 to Phred scores to make them ASCII
static inline void bam_rewrite_qual (char *qual, uint32_t qual_len)
{
    for (uint32_t i=0; i < qual_len; i++)
        qual[i] += 33;
}

static inline const char *bam_rewrite_one_optional_number (VBlockSAM *vb, const char *next_field, uint8_t type)
{
    static const char special_float[2] = { SNIP_SPECIAL, SAM_SPECIAL_FLOAT };

    switch (type) {
        case 'c': { uint8_t  n = NEXT_UINT8 ; vb->textual_opt.len += str_int ((int8_t )n, AFTERENT (char, vb->textual_opt)); break; }
        case 'C': { uint8_t  n = NEXT_UINT8 ; vb->textual_opt.len += str_int (         n, AFTERENT (char, vb->textual_opt)); break; }
        case 's': { uint16_t n = NEXT_UINT16; vb->textual_opt.len += str_int ((int16_t)n, AFTERENT (char, vb->textual_opt)); break; }
        case 'S': { uint16_t n = NEXT_UINT16; vb->textual_opt.len += str_int (         n, AFTERENT (char, vb->textual_opt)); break; }
        case 'i': { uint32_t n = NEXT_UINT32; vb->textual_opt.len += str_int ((int32_t)n, AFTERENT (char, vb->textual_opt)); break; }
        case 'F': buf_add (&vb->textual_opt, special_float, 2); // add this prefix and then fall through - the float casted to uint32
        case 'I': { uint32_t n = NEXT_UINT32; vb->textual_opt.len += str_int (         n, AFTERENT (char, vb->textual_opt)); break; }
        default: ABORT ("Error in bam_rewrite_one_optional_number: enrecognized Optional field type '%c' (ASCII %u) in vb=%u", 
                        type, type, vb->vblock_i);
    }    

    return next_field;
} 

const char *bam_get_one_optional (VBlockSAM *vb, const char *next_field,
                                  const char **tag, char *type, const char **value, unsigned *value_len) // out
{
    *tag  = next_field;
    *type = next_field[2];
    next_field += 3;

    if (*type == 'Z' || *type == 'H') {
        *value = next_field;
        *value_len = strlen (*value);
        return next_field + *value_len + 1; // +1 for \0
    }
    
    else if (*type == 'A') {
        *value = next_field;
        *value_len = 1;
        return next_field + 1;
    } 

    uint32_t max_len = (*type == 'B') ? (12 * GET_UINT32 (next_field) + 10) : // worst case scenario for item: "-1000000000,"
                       30;
    buf_alloc (vb, &vb->textual_opt, max_len, 2, "textual_opt", 0); // a rather inefficient max in case of arrays, to do: tighten the calculation

    if (*type != 'B')
        next_field = bam_rewrite_one_optional_number (vb, next_field, *type);

    else { // 'B'
        char subtype = *next_field++; // type of elements of array
        uint32_t num_elements = NEXT_UINT32;

        NEXTENT (uint8_t, vb->textual_opt) = subtype;
        
        for (uint32_t i=0; i < num_elements; i++) {
            NEXTENT (char, vb->textual_opt) = ',';
            next_field = bam_rewrite_one_optional_number (vb, next_field, subtype);
        }
    }    

    *value     = vb->textual_opt.data;
    *value_len = vb->textual_opt.len;

    return next_field;
}

const char *bam_seg_txt_line (VBlock *vb_, const char *alignment /* BAM terminology for one line */, bool *has_13_unused)   
{
    VBlockSAM *vb = (VBlockSAM *)vb_;
    ZipDataLineSAM *dl = DATA_LINE (vb->line_i);
    const char *next_field = alignment;
    char snip[50];
    unsigned snip_len;

    uint32_t alignment_size = NEXT_UINT32;

    // ref_id (RNAME)
    int32_t ref_id = (int32_t)NEXT_UINT32; // corresponding to CHROMs in the BAM header
    bam_seg_ref_id (vb_, SAM_RNAME, ref_id, -1);

    // POS
    PosType this_pos = 1 + NEXT_UINT32; // pos in BAM is 0 based, -1 for unknown? 
    seg_pos_field (vb_, SAM_POS, SAM_POS, false, 0, 0, this_pos, sizeof (uint32_t));
    ASSERT (!(ref_id==-1 && this_pos), "Error in %s vb=%u: RNAME=\"*\" - expecting POS to be 0 but it is %"PRId64, txt_name, vb->vblock_i, this_pos);

    uint8_t l_read_name = NEXT_UINT8; // QNAME length

    // MAPQ
    snip_len = str_int (NEXT_UINT8, snip);
    seg_by_did_i (vb, snip, snip_len, SAM_MAPQ, 1);

    // BIN
    uint16_t bin = NEXT_UINT16;
    bam_seg_bin (vb_, bin);

    uint16_t n_cigar_op = NEXT_UINT16;

    // FLAG
    uint16_t flag = NEXT_UINT16;
    snip_len = str_int (flag, snip);
    seg_by_did_i (vb, snip, snip_len, SAM_FLAG, 2);
    
    // l_seq (=seq_len)
    dl->seq_len = NEXT_UINT32;

    // next_ref_id (RNEXT)
    int32_t next_ref_id = (int32_t)NEXT_UINT32; // corresponding to CHROMs in the BAM header
    bam_seg_ref_id (vb_, SAM_RNEXT, next_ref_id, ref_id);

    // PNEXT
    PosType pnext = 1 + NEXT_UINT32; // pos in BAM is 0 based, -1 for unknown? 
    seg_pos_field (vb_, SAM_PNEXT, SAM_POS, false, 0, 0, pnext, sizeof (uint32_t));

    // TLEN
    int32_t tlen_value = (int32_t)NEXT_UINT32;
    sam_seg_tlen_field (vb, 0, 0, (int64_t)tlen_value, vb->contexts[SAM_PNEXT].last_delta, dl->seq_len);

    // QNAME
    seg_compound_field ((VBlockP)vb, &vb->contexts[SAM_QNAME], next_field, l_read_name-1, false, 0, 2 /* account for \0 and l_read_name */);
    next_field += l_read_name; // inc. \0

    // CIGAR
    bam_rewrite_cigar (vb, n_cigar_op, (uint32_t*)next_field); // re-write BAM format CIGAR as SAM textual format in vb->textual_cigar
    sam_analyze_cigar (vb->textual_cigar.data, vb->textual_cigar.len, NULL, &vb->ref_consumed, &vb->ref_and_seq_consumed);
    next_field += n_cigar_op * sizeof (uint32_t);

    // SEQ - calculate diff vs. reference (self or imported)
    bam_rewrite_seq (vb, dl, next_field);
    sam_seg_seq_field (vb, vb->textual_seq.data, vb->textual_seq.len, this_pos, vb->last_cigar, 0, vb->textual_seq.len, 
                       vb->last_cigar, (dl->seq_len+1)/2 + sizeof (uint32_t) /* account for l_seq and seq fields */);
    next_field += (dl->seq_len+1)/2; 

    // QUAL
    if (dl->seq_len) {
        bam_rewrite_qual ((char*)next_field, dl->seq_len); // add 33 to Phred scores to make them ASCII
        sam_seg_qual_field (vb, dl, next_field, dl->seq_len, dl->seq_len /* account for qual field */ );
    }
    else {
        *(char *)alignment = '*'; // overwrite as we need it somewhere in txt_data
        sam_seg_qual_field (vb, dl, alignment, 1, 0);
    }
    next_field += dl->seq_len; 

    // finally we can seg the textual CIGAR now (including if n_cigar_op=0)
    bam_seg_cigar_field (vb, dl, n_cigar_op);

    // OPTIONAL fields - up to MAX_SUBFIELDS of them
    next_field = sam_seg_optional_all (vb, dl, next_field, 0,0,0, alignment + alignment_size + sizeof (uint32_t));
    
    buf_free (&vb->textual_cigar);
    buf_free (&vb->textual_seq);

    return next_field;
}
