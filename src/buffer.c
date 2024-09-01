// ------------------------------------------------------------------
//   buffer.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "buffer.h"
#include "file.h"

// writes a buffer to a file, return true if successful
// note: this is designed to run in any there, so it cannot create any buffers in evb
bool buf_dump_to_file (rom filename, ConstBufferP buf, unsigned buf_word_width, bool including_control_region, 
                       bool no_dirs, bool verbose, bool do_gzip)
{
    RETURNW (buf->type == BUF_REGULAR, false, 
             "FYI: failed to dump buffer.type=%s name=%s while putting %s", 
             buf_type_name (buf), buf->name ? buf->name : "(null)", filename);

    int fn_len = strlen(filename);
    char update_filename[fn_len + 10];
    strcpy (update_filename, filename);

    if (no_dirs) {
        for (unsigned i=0; i < fn_len; i++)
            if (filename[i] == '/' || (flag.is_windows && (filename[i] == '\\' || (filename[i] == ':' && i!=1))))
                update_filename[i] = '-';
        filename = update_filename;
    }

    bool success;
    if (including_control_region) {
        ASSERT (*(uint64_t *)(buf->memory)           == UNDERFLOW_TRAP, "dumping to %s: buffer has underflowed", filename);
        ASSERT (*(uint64_t *)(buf->data + buf->size) == OVERFLOW_TRAP,  "dumping to %s: buffer has overflowed",  filename);

        success = file_put_data (update_filename, buf->memory, buf_mem_size (buf), 0);
    }
    else
        success = file_put_data (update_filename, buf->data, buf->len * buf_word_width, 0);

    if (success && do_gzip) 
        file_gzip (update_filename); // updates filename if successful

    if (success && verbose) iprintf ("\nDumped file %s\n", update_filename);

    return success;
}

// copy data - possibly within the same buffer
void buf_copy_do (VBlockP dst_vb, BufferP dst, ConstBufferP src, 
                  uint64_t bytes_per_entry, // how many bytes are counted by a unit of .len
                  uint64_t src_start_entry, uint64_t max_entries,  // if 0 copies the entire buffer 
                  FUNCLINE,
                  rom dst_name) // dst buffer settings, or take from src if 0
{
    ASSERTNOTNULL (src);
    ASSERTNOTNULL (dst);

    ASSERT (src->data, "called from %s:%u: src->data is NULL", func, code_line);
    
    ASSERT (!max_entries || src_start_entry < src->len, 
            "buf_copy of %s called from %s:%u: src_start_entry=%"PRIu64" is larger than src->len=%"PRIu64, buf_desc(src).s, func, code_line, src_start_entry, src->len);

    uint64_t num_entries = src->len - MIN_(src_start_entry, src->len);
    if (max_entries && max_entries < num_entries) num_entries = max_entries;
    // up to 15.0.63 this was buggy: uint64_t num_entries = max_entries ? MIN_(max_entries, src->len - src_start_entry) : src->len - src_start_entry;

    if (!bytes_per_entry) bytes_per_entry=1;
    
    if (num_entries) {
        buf_alloc_(dst_vb, dst, 0, num_entries * bytes_per_entry, 1, 1, dst_name ? dst_name : src->name, func, code_line); 

        if (dst != src || src_start_entry >= num_entries)
            memcpy (dst->data, &src->data[src_start_entry * bytes_per_entry], num_entries * bytes_per_entry);
        else 
            memmove (dst->data, &src->data[src_start_entry * bytes_per_entry], num_entries * bytes_per_entry); // need memmove for overlapping areas
    }

    dst->len = num_entries;  
}   

// removes a section from the buffer
void buf_remove_do (BufferP buf, unsigned sizeof_item, uint64_t remove_start, uint64_t remove_len)
{
    if (!remove_len) return;

    ASSERT (remove_start + remove_len <= buf->len, "Out of range: remove_start=%"PRIu64" + remove_len=%"PRIu64" > buf->len=%"PRIu64,
            remove_start, remove_len, buf->len);

    if (remove_len != buf->len) { // skip in common case of deleting entire buffer 
        uint64_t remove_start_byte = remove_start * sizeof_item;
        uint64_t remove_bytes      = remove_len   * sizeof_item;
        memmove (buf->data + remove_start_byte, 
                 buf->data + remove_start_byte + remove_bytes, 
                 (buf->len * sizeof_item) - (remove_start_byte + remove_bytes));
    }

    buf->len -= remove_len;
}

void buf_add (BufferP buf, STRp(data))
{
    if (!data_len) return; // don't test for space (and "data" pointer) if length is 0

    ASSERT (buf_has_space (buf, data_len),
            "buf_add: buffer %s is out of space: len=%u size=%u data_len=%u",
            buf_desc (buf).s, (uint32_t)buf->len, (uint32_t)buf->size, data_len);
    buf_add_do (buf, data, data_len);
}

void buf_insert_do (VBlockP vb, BufferP buf, unsigned width, uint64_t insert_at, const void *new_data, uint64_t new_data_len, rom name, FUNCLINE) 
{ 
    if (!new_data_len) return;

    buf_alloc_(vb ? vb : buf->vb, buf, new_data_len + (width==1)/*room for \0 or separator if char*/, 0, width, CTX_GROWTH, name, func, code_line); 

    if (insert_at != buf->len) {
        ASSERT (insert_at < buf->len, "called from %s:%u: expecting insert_at=%"PRIu64" <= buf->len=%"PRIu64" in buf=%s", 
                func, code_line, insert_at, buf->len, buf_desc(buf).s);

        memmove (&buf->data[(insert_at + new_data_len) * width], &buf->data[insert_at * width], (buf->len - insert_at) * width);
    }

    memcpy (&buf->data[insert_at * width], new_data, new_data_len * width);   
    buf->len += new_data_len; 

    ASSERT (BOVERFLOW(buf) == OVERFLOW_TRAP, "buffer overflow: %s", buf_desc(buf).s);
} 

void buf_append_string (VBlockP vb, BufferP buf, rom str) 
{ 
    uint64_t len = strlen (str); 
    ASSERT (len < 10000000, "len=%"PRIu64" too long, looks like a bug", len);

    buf_add_more (vb, buf, str, len, buf->name ? buf->name : "string_buf"); // allocates one char extra
    *BAFTc (*buf) = '\0'; // string terminator without increasing buf->len
}

// swaps buffers' content without affecting buffer list
void buf_swap (BufferP buf1, BufferP buf2)
{
    ASSERT (buf1->vb == buf2->vb && 
            buf1->type == BUF_REGULAR && buf2->type == BUF_REGULAR &&
            !buf1->shared && !buf2->shared,
            "buf1->vb != buf2->vb or not REGULAR or shared. buf1=%s buf2=%s", buf_desc (buf1).s, buf_desc (buf2).s);

    SWAP (buf1->memory,   buf2->memory);
    SWAP (buf1->data,     buf2->data);
    SWAP (buf1->len,      buf2->len);
    SWAP (buf1->param,    buf2->param);
    SWAPbits (buf1->size, buf2->size);
}

void buf_print (BufferP buf, bool add_newline)
{
    for (uint64_t i=0; i < buf->len; i++) 
        fputc (buf->data[i], info_stream);  // safer than printf %.*s ?

    iprint0 (add_newline ? "\n" : "");
}

// iterator on a buffer containing newline-terminated lines
// false means continue iterating, true means stop
char *buf_foreach_line (BufferP buf,
                        bool reverse, // iterate backwards
                        TxtIteratorCallback callback, 
                        void *cb_param1, void *cb_param2, unsigned cb_param3, // passed as-is to callback
                        int64_t *line_len) // out
{
    if (line_len) *line_len = 0;

    if (!buf->len) return NULL;

    char *firstbuf = buf->data;
    char *afterbuf = BAFTc (*buf);

    char *first = !reverse ? firstbuf : 0;
    char *after = !reverse ? 0 : afterbuf;

    while (1) {
            
        // get one line - searching forward or backwards
        if (!reverse) {
            for (after=first ; after < afterbuf && *after != '\n' ; after++);
            after++; // skip newline
        }
        else {
            for (first=after-2 /* skip final \n */; first >= firstbuf && *first != '\n'; first--);
            first++; // after detected \n or at start of line
        }

        if (!reverse && after > afterbuf) return NULL; // we don't call callback if after>afterbuf - beyond end of line
            
        if (callback (first, after - first, cb_param1, cb_param2, cb_param3)) {
            if (line_len) *line_len = after - first;
            return first;
        }

        if (reverse && first == firstbuf) return NULL; // beginning of line - we called the cb

        if (!reverse) first=after;
        else          after=first;
    }

    return 0; // never reaches here
}   

//---------------------
// Bits stuff
//---------------------

BitsP buf_alloc_bits_do (VBlockP vb, BufferP buf, uint64_t nbits, BitsInitType init_to, float grow_at_least_factor, rom name, FUNCLINE)
{
    ASSERT0 (buf->type == BUF_UNALLOCATED || buf->type == BUF_REGULAR, "buf needs to be BUF_UNALLOCATED or BUF_REGULAR");

    uint64_t nwords = roundup_bits2words64 (nbits);
    uint64_t old_nbits = buf->nbits;
    
    buf_alloc_(vb, buf, 0, nwords, sizeof(uint64_t), grow_at_least_factor, name, func, code_line);

    buf->nbits  = nbits;   
    buf->nwords = nwords; 

    bits_clear_excess_bits_in_top_word ((BitsP)buf, false);

    if (init_to == CLEAR) {
        if (!old_nbits) bits_clear_all ((BitsP)buf);
        else if (nbits > old_nbits) bits_clear_region ((BitsP)buf, old_nbits, nbits - old_nbits);           
    } 
    if (init_to == SET) {
        if (!old_nbits) bits_set_all ((BitsP)buf);
        else if (nbits > old_nbits) bits_set_region ((BitsP)buf, old_nbits, nbits - old_nbits);           
    } 

    return (BitsP)buf;
}

BitsP buf_overlay_bits_do (VBlockP vb,
                           BufferP top_buf, BufferP bottom_buf,  
                           uint64_t start_byte_in_bottom_buf,
                           uint64_t nbits,
                           FUNCLINE, rom name)
{
    uint64_t nwords = roundup_bits2words64 (nbits);

    buf_overlay_do (evb, top_buf, bottom_buf, start_byte_in_bottom_buf, false, func, code_line, name);

    top_buf->nbits  = nbits;
    top_buf->nwords = nwords;
    return (BitsP)top_buf;
}

// convert a Buffer from a z_file section whose len is in char to a bits
Bits *buf_zfile_buf_to_bits (BufferP buf, uint64_t nbits)
{
    ASSERT (roundup_bits2bytes (nbits) <= buf->len, "nbits=%"PRId64" indicating a length of at least %"PRId64", but buf->len=%"PRId64,
            nbits, roundup_bits2bytes (nbits), buf->len);

    Bits *bits = (BitsP)buf;
    bits->nbits  = nbits;
    bits->nwords = roundup_bits2words64 (bits->nbits);

    ASSERT (roundup_bits2bytes64 (nbits) <= buf->size, "buffer to small: buf->size=%"PRId64" but bits has %"PRId64" words and hence requires %"PRId64" bytes",
            (uint64_t)buf->size, bits->nwords, bits->nwords * sizeof(uint64_t));

    LTEN_bits (bits);

    bits_clear_excess_bits_in_top_word (bits, false);

    return bits;
}

void buf_add_bit (BufferP buf, int64_t new_bit) 
{
    Bits *bar = (BitsP)buf;

    ASSERT (bar->nbits < buf->size * 8, "no room in Buffer %s to extend the bitmap", buf->name);
    bar->nbits++;     
    if (bar->nbits % 64 == 1) { // starting a new word                
        bar->nwords++;
        bar->words[bar->nwords-1] = new_bit; // LSb is as requested, other 63 bits are 0
    } 
    else
        bits_assign (bar, bar->nbits-1, new_bit);  
}
 
uint64_t buf_extend_bits (BufferP buf, int64_t num_new_bits) 
{
    Bits *bar = (BitsP)buf;

    ASSERT (bar->nbits + num_new_bits <= buf->size * 8, "Error in %s:%u: no room in Buffer %s to extend the bitmap: nbits=%"PRIu64", num_new_bits=%"PRId64", buf->size=%"PRIu64, 
            __FUNCLINE, buf->name, bar->nbits, num_new_bits, (uint64_t)buf->size);
    
    uint64_t next_bit = bar->nbits;

    bar->nbits += num_new_bits;     
    bar->nwords = roundup_bits2words64 (bar->nbits);
    bits_clear_excess_bits_in_top_word (bar, true);

    return next_bit;
}
 
//---------------------
// Endianity stuff
//---------------------

void interlace_d8_buf       (BufferP buf, LocalType *lt) { for_buf (int8_t,  num, *buf) *num =        (INTERLACE(int8_t,  *num)); }
void BGEN_interlace_d16_buf (BufferP buf, LocalType *lt) { for_buf (int16_t, num, *buf) *num = BGEN16 (INTERLACE(int16_t, *num)); }
void BGEN_interlace_d32_buf (BufferP buf, LocalType *lt) { for_buf (int32_t, num, *buf) *num = BGEN32 (INTERLACE(int32_t, *num)); }
void BGEN_interlace_d64_buf (BufferP buf, LocalType *lt) { for_buf (int64_t, num, *buf) *num = BGEN64 (INTERLACE(int64_t, *num)); }
void LTEN_interlace_d16_buf (BufferP buf, LocalType *lt) { for_buf (int16_t, num, *buf) *num = LTEN16 (INTERLACE(int16_t, *num)); }
void LTEN_interlace_d32_buf (BufferP buf, LocalType *lt) { for_buf (int32_t, num, *buf) *num = LTEN32 (INTERLACE(int32_t, *num)); }
void LTEN_interlace_d64_buf (BufferP buf, LocalType *lt) { for_buf (int64_t, num, *buf) *num = LTEN64 (INTERLACE(int64_t, *num)); }

void BGEN_u8_buf  (BufferP buf, LocalType *lt) {}
void BGEN_u16_buf (BufferP buf, LocalType *lt) { if ( flag.is_lten) for_buf (uint16_t, num, *buf) *num = BGEN16 (*num); }
void BGEN_u32_buf (BufferP buf, LocalType *lt) { if ( flag.is_lten) for_buf (uint32_t, num, *buf) *num = BGEN32 (*num); }
void BGEN_u64_buf (BufferP buf, LocalType *lt) { if ( flag.is_lten) for_buf (uint64_t, num, *buf) *num = BGEN64 (*num); }
void LTEN_u16_buf (BufferP buf, LocalType *lt) { if (!flag.is_lten) for_buf (uint16_t, num, *buf) *num = LTEN16 (*num); }
void LTEN_u32_buf (BufferP buf, LocalType *lt) { if (!flag.is_lten) for_buf (uint32_t, num, *buf) *num = LTEN32 (*num); }
void LTEN_u64_buf (BufferP buf, LocalType *lt) { if (!flag.is_lten) for_buf (uint64_t, num, *buf) *num = LTEN64 (*num); }

// number of columns is trasmitted in the count, except if this is a matrix of VCF samples, in which case param=0 and we take 
// the number of columns to be the number of samples in the VCF header
static inline uint32_t BGEN_transpose_num_cols (ConstBufferP buf)
{
    uint32_t cols = buf->n_cols; // cols and rows in terms of the target non-transposed matrix (0 if vcf_num_samples)

    if (!cols) cols = vcf_header_get_num_samples(); 
    ASSERT0 (cols, "vcf_header_get_num_samples=0");
    
    return cols;
}

void BGEN_transpose_u8_buf (BufferP buf, LocalType *lt)
{
    if (!buf->len) return;

    uint32_t cols = BGEN_transpose_num_cols (buf);
    uint32_t rows = buf->len / cols;

    buf_alloc (buf->vb, &buf->vb->scratch, 0, buf->len, uint8_t, 1, "scratch");
    ARRAY (uint8_t, target, buf->vb->scratch);
    ARRAY (uint8_t, transposed, *buf);

    for (uint32_t c=0; c < cols; c++) 
        for (uint32_t r=0; r < rows; r++) 
            target[r * cols + c] = transposed[c * rows + r];

    buf->vb->scratch.len = buf->len;
    buf_copy (buf->vb, buf, &buf->vb->scratch, uint8_t, 0, 0, CTX_TAG_LOCAL); // copy and not move, so we can keep local's memory for next vb

    buf_free (buf->vb->scratch);

    if (lt) *lt = LT_UINT8; // no longer transposed
}

void BGEN_transpose_u16_buf (BufferP buf, LocalType *lt)
{
    if (!buf->len) return;

    uint32_t cols = BGEN_transpose_num_cols (buf);
    uint32_t rows = buf->len / cols;

    buf_alloc (buf->vb, &buf->vb->scratch, 0, buf->len, uint16_t, 1, "scratch");
    ARRAY (uint16_t, target, buf->vb->scratch);
    ARRAY (uint16_t, transposed, *buf);

    for (uint32_t c=0; c < cols; c++) 
        for (uint32_t r=0; r < rows; r++) 
            target[r * cols + c] = BGEN16 (transposed[c * rows + r]);

    buf->vb->scratch.len = buf->len;
    buf_copy (buf->vb, buf, &buf->vb->scratch, uint16_t, 0, 0, CTX_TAG_LOCAL); // copy and not move, so we can keep local's memory for next vb

    buf_free (buf->vb->scratch);

    *lt = LT_UINT16; // no longer transposed
}

void BGEN_transpose_u32_buf (BufferP buf, LocalType *lt)
{
    if (!buf->len) return;

    uint32_t cols = BGEN_transpose_num_cols (buf);
    uint32_t rows = buf->len / cols;

    buf_alloc (buf->vb, &buf->vb->scratch, 0, buf->len, uint32_t, 1, "scratch");
    ARRAY (uint32_t, target, buf->vb->scratch);
    ARRAY (uint32_t, transposed, *buf);

    for (uint32_t c=0; c < cols; c++) 
        for (uint32_t r=0; r < rows; r++) 
            target[r * cols + c] = BGEN32 (transposed[c * rows + r]);

    buf->vb->scratch.len = buf->len;
    buf_copy (buf->vb, buf, &buf->vb->scratch, uint32_t, 0, 0, CTX_TAG_LOCAL); // copy and not move, so we can keep local's memory for next vb

    buf_free (buf->vb->scratch);

    *lt = LT_UINT32; // no longer transposed
}

void BGEN_deinterlace_d8_buf (BufferP buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint8_t unum = *B8 (*buf, i);
        *B(int8_t, *buf, i) = DEINTERLACE(int8_t,unum); 
    }
}

void BGEN_deinterlace_d16_buf (BufferP buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint16_t num_big_en = *B16 (*buf, i);
        uint16_t unum = BGEN16 (num_big_en);
        *B(int16_t, *buf, i) = DEINTERLACE(int16_t,unum); 
    }
}

void BGEN_deinterlace_d32_buf (BufferP buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint32_t num_big_en = *B32 ( *buf, i);
        uint32_t unum = BGEN32 (num_big_en);
        *B(int32_t, *buf, i) = DEINTERLACE(int32_t,unum); 
    }
}

void BGEN_deinterlace_d64_buf (BufferP buf, LocalType *lt)
{
    for (uint64_t i=0; i < buf->len; i++) {
        uint64_t num_big_en = *B64 (*buf, i);
        uint64_t unum = BGEN64 (num_big_en);
        *B(int64_t, *buf, i) = DEINTERLACE(int64_t,unum); 
    }
}
