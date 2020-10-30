// ------------------------------------------------------------------
//   sam_zip.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "sam_private.h"
#include "reference.h"
#include "seg.h"
#include "context.h"
#include "file.h"
#include "random_access.h"
#include "endianness.h"
#include "strings.h"
#include "zip.h"
#include "optimize.h"
#include "dict_id.h"
#include "codec.h"
#include "aligner.h"
#include "container.h"

#define IS_BAM (txt_file->data_type==DT_BAM)

// called by I/O thread (zip_initialize callback)
void sam_zip_initialize (void)
{
    // if there is no external reference provided, then we create our internal one, and store it
    // (if an external reference IS provided, the user can decide whether to store it or not, with --store-reference)
    if (flag_reference == REF_NONE) flag_reference = REF_INTERNAL;

    // in case of internal reference, we need to initialize. in case of --reference, it was initialized by ref_load_external_reference()
    if (!flag_reference || flag_reference == REF_INTERNAL) ref_initialize_ranges (RT_DENOVO); // it will be REF_INTERNAL if this is the 2nd+ non-conatenated file

    // evb buffers must be alloced by I/O threads, since other threads cannot modify evb's buf_list
    random_access_alloc_ra_buf (evb, 0);
}

bool sam_inspect_txt_header (BufferP txt_header)
{
    if (!(flag_reference == REF_EXTERNAL || flag_reference == REF_EXT_STORE)) return true; // we're not using a reference - all is good

    char *line = txt_header->data;

    while (1) {
        line = strchr (line, '@');
        if (!line) break;

        if (line[1] == 'S' && line[2] == 'Q') { // this test will always be in the buffer - possible in the overflow area
            char *newline = strchr (line, '\n');
            *newline = 0;

            // line looks something like: @SQ\tSN:chr10\tLN:135534747 (no strict order between SN and LN and there could be more parameters)
            char *chrom_name   = strstr (line, "SN:"); 
            char *last_pos_str = strstr (line, "LN:"); 

            if (chrom_name && last_pos_str) {
                unsigned chrom_name_len = strcspn (&chrom_name[3], "\t\n\r");
                PosType last_pos = (PosType)strtoull (&last_pos_str[3], NULL, 10);
                ref_contigs_verify_identical_chrom (&chrom_name[3], chrom_name_len, last_pos, WORD_INDEX_NONE);
            }

            line = newline+1;

            *newline = '\n'; // restore
        }
        else 
            line++;
    }

    return true;
}

// callback function for compress to get data of one line (called by codec_bz2_compress)
void sam_zip_qual (VBlock *vb, uint32_t vb_line_i, 
                                        char **line_qual_data, uint32_t *line_qual_len, // out
                                        char **line_u2_data,   uint32_t *line_u2_len,
                                        uint32_t maximum_len) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_qual_len  = MIN (maximum_len, dl->qual_data_len);
    
    maximum_len -= *line_qual_len; 
    *line_u2_len    = MIN (maximum_len, dl->u2_data_len);

    if (!line_qual_data) return; // only lengths were requested

    *line_qual_data = ENT (char, vb->txt_data, dl->qual_data_start);
    *line_u2_data   = dl->u2_data_start ? ENT (char, vb->txt_data, dl->u2_data_start) : NULL;

    // if QUAL is just "*" (i.e. unavailable) replace it by " " because '*' is a legal PHRED quality value that will confuse PIZ
    if (dl->qual_data_len == 1 && (*line_qual_data)[0] == '*') 
        *line_qual_data = " "; // pointer to static string

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    else if (flag_optimize_QUAL) {
        optimize_phred_quality_string (*line_qual_data, *line_qual_len);
        if (*line_u2_data) optimize_phred_quality_string (*line_u2_data, *line_u2_len);
    }
}

// callback function for compress to get BD_BI data of one line: this is an
// interlaced line containing a character from BD followed by a character from BI - since these two fields are correlated
// note: if only one of BD or BI exists, the missing data in the interlaced string will be 0 (this should is not expected to ever happen)
void sam_zip_bd_bi (VBlock *vb_, uint32_t vb_line_i, 
                    char **line_data, uint32_t *line_len,  // out 
                    char **unused1,  uint32_t *unused2,
                    uint32_t maximum_len)
{
    VBlockSAM *vb = (VBlockSAM *)vb_;
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
    
    const char *bd = dl->bdbi_data_start[0] ? ENT (char, vb->txt_data, dl->bdbi_data_start[0]) : NULL;
    const char *bi = dl->bdbi_data_start[1] ? ENT (char, vb->txt_data, dl->bdbi_data_start[1]) : NULL;
    
    if (!bd && !bi) return; // no BD or BI on this line

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_len  = MIN (maximum_len, dl->seq_len * 2);

    if (!line_data) return; // only length was requested

    buf_alloc (vb, &vb->bd_bi_line, dl->seq_len * 2, 2, "bd_bi_line", 0);

    // calculate character-wise delta
    for (unsigned i=0; i < dl->seq_len; i++) {
        *ENT (uint8_t, vb->bd_bi_line, i*2    ) = bd ? bd[i] : 0;
        *ENT (uint8_t, vb->bd_bi_line, i*2 + 1) = bi ? bi[i] - (bd ? bd[i] : 0) : 0;
    }

    *line_data = FIRSTENT (char, vb->bd_bi_line);
}   

void sam_seg_initialize (VBlock *vb)
{
    vb->contexts[SAM_RNAME].inst     = CTX_INST_NO_STONS; // needs b250 node_index for random access
    vb->contexts[SAM_SQBITMAP].ltype = LT_BITMAP;
    vb->contexts[SAM_TLEN].flags     = CTX_FL_STORE_INT;
    vb->contexts[SAM_OPTIONAL].flags = CTX_FL_CONTAINER;
    vb->contexts[SAM_STRAND].ltype   = LT_BITMAP;
    vb->contexts[SAM_GPOS].ltype     = LT_UINT32;
    vb->contexts[SAM_GPOS].flags     = CTX_FL_STORE_INT;

    // when reconstructing BAM, we output the word_index instead of the string
    vb->contexts[SAM_RNAME].flags    = CTX_FL_STORE_INDEX;
    vb->contexts[SAM_RNEXT].flags    = CTX_FL_STORE_INDEX;

    // in --stats, consolidate stats into SQBITMAP
    vb->contexts[SAM_GPOS].st_did_i = vb->contexts[SAM_STRAND].st_did_i =
    vb->contexts[SAM_NONREF].st_did_i = vb->contexts[SAM_NONREF_X].st_did_i = SAM_SQBITMAP;

    codec_acgt_comp_init (vb);
}

void sam_seg_finalize (VBlockP vb)
{
    // for qual data - select domqual compression if possible, or fallback 
    if (!codec_domq_comp_init (vb, sam_zip_qual)) {
        vb->contexts[SAM_QUAL].ltype  = LT_SEQUENCE; 
        vb->contexts[SAM_QUAL].inst   = 0; // don't inherit from previous file 
    }

    // top level snip - reconstruction as SAM
    Container top_level_sam = { 
        .repeats   = vb->lines.len,
        .flags     = CONTAINER_TOPLEVEL,
        .num_items = 13,
        .items     = { { (DictId)dict_id_fields[SAM_QNAME],    DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[SAM_FLAG],     DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[SAM_RNAME],    DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[SAM_POS],      DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[SAM_MAPQ],     DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[SAM_CIGAR],    DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[SAM_RNEXT],    DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[SAM_PNEXT],    DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[SAM_TLEN],     DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[SAM_SQBITMAP], DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[SAM_QUAL],     DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[SAM_OPTIONAL], DID_I_NONE       },
                       { (DictId)dict_id_fields[SAM_EOL],      DID_I_NONE       } }
    };
    seg_container_by_ctx (vb, &vb->contexts[SAM_TOPLEVEL], &top_level_sam, 0, 0, 0);

    // top level snip - reconstruction as BAM
    // strategy: we start by reconstructing the variable-length fields first (after a prefix that sets them in place) 
    // - read_name, cigar, seq and qual - and then go back and fill in the fixed-location fields
    // Translation (a feature of Container): items reconstruct their data and then call a translation function to translate it to the desired format
    Container top_level_bam = { 
        .repeats   = vb->lines.len,
        .flags     = CONTAINER_TOPLEVEL,
        .num_items = 13,
        .items     = { { (DictId)dict_id_fields[SAM_QNAME],    DID_I_NONE, { CI_NUL_TERMINATE } }, // Translate - move QNAME forward in txt_data to its correct place + nul-terminate
                       { (DictId)dict_id_fields[SAM_CIGAR],    DID_I_NONE, "", TRS_SAM2BAM_CIGAR    }, // Translate - translate textual to BAM CIGAR format + reconstruct l_read_name, n_cigar_op, l_seq, block_size
                       { (DictId)dict_id_fields[SAM_SQBITMAP], DID_I_NONE, "", TRS_SAM2BAM_SEQ      }, // Translate - textual format to BAM format
                       { (DictId)dict_id_fields[SAM_QUAL],     DID_I_NONE, "", TRS_SAM2BAM_QUAL     }, // Translate - textual format to BAM format, save txt_data.len, and move it back to the position of ref_id
                       { (DictId)dict_id_fields[SAM_RNAME],    DID_I_NONE, { CI_DONT_RECONSTRUCT }, TRS_SAM2BAM_RNAME    }, // Translate - output word_index instead of string
                       { (DictId)dict_id_fields[SAM_POS],      DID_I_NONE, { CI_DONT_RECONSTRUCT }, TRS_SAM2BAM_POS      }, // Translate - output little endian POS-1
                       { (DictId)dict_id_fields[SAM_MAPQ],     DID_I_NONE, { CI_DONT_RECONSTRUCT }, TRS_U8               }, // Translate - textual to binary number
                       { (DictId)dict_id_fields[SAM_BAM_BIN],  DID_I_NONE, { CI_DONT_RECONSTRUCT }, TRS_LTEN_U16         }, // Translate - textual to binary number
                       { (DictId)dict_id_fields[SAM_FLAG],     DID_I_NONE, { CI_DONT_RECONSTRUCT }, TRS_LTEN_U16         }, // Translate - textual to binary number
                       { (DictId)dict_id_fields[SAM_RNEXT],    DID_I_NONE, { CI_DONT_RECONSTRUCT }, TRS_SAM2BAM_RNAME    }, // Translate - output word_index instead of string
                       { (DictId)dict_id_fields[SAM_PNEXT],    DID_I_NONE, { CI_DONT_RECONSTRUCT }, TRS_SAM2BAM_POS      }, // Translate - output little endian POS-1
                       { (DictId)dict_id_fields[SAM_TLEN],     DID_I_NONE, { CI_DONT_RECONSTRUCT }, TRS_LTEN_I32         }, // Translate - textual to binary number
                       { (DictId)dict_id_fields[SAM_OPTIONAL], DID_I_NONE, { CI_DONT_RECONSTRUCT }, TRS_SAM2BAM_OPTIONAL } } // transforms prefixes MX:i->MXC and moves txt_data.len forward to after qual
    };

    // 36 characters (of 0) will be written first, before the QNAME. We will override them after.
    static const char bam_line_prefix[37] = { [36] = SNIP_CONTAINER }; 

    seg_container_by_ctx (vb, &vb->contexts[SAM_TOP2BAM], &top_level_bam, bam_line_prefix, sizeof(bam_line_prefix), 
                          IS_BAM ? sizeof (uint32_t) * vb->lines.len : 0); // if BAM, account for block_size

    // top level snip - reconstruction as FASTQ
    Container top_level_fastq = { 
        .repeats   = vb->lines.len,
        .flags     = CONTAINER_TOPLEVEL | CONTAINER_FILTER_REPEATS, // filter - drop non-primary chimeric reads and reads without QUAL data
        .num_items = 3,
        .items     = { { (DictId)dict_id_fields[SAM_QNAME],    DID_I_NONE, "\n" },
                       { (DictId)dict_id_fields[SAM_SQBITMAP], DID_I_NONE, "\n", TRS_SAM2FASTQ_SEQ  }, // Translate - reverse complement if FLAGS & 0x10
                       { (DictId)dict_id_fields[SAM_QUAL],     DID_I_NONE, "\n" } // Translate - textual format to BAM format, save txt_data.len, and move it back to the position of ref_id
                     }
    };

    // use a prefix to at the + line    
    static const char fastq_line_prefix[6] = { SNIP_CONTAINER, SNIP_CONTAINER, SNIP_CONTAINER, '+', '\n', SNIP_CONTAINER };

    seg_container_by_ctx (vb, &vb->contexts[SAM_TOP2FQ], &top_level_fastq, fastq_line_prefix, sizeof(fastq_line_prefix), 0);
}

// TLEN - 3 cases: 
// 1. if a non-zero value that is the negative of the previous line - a SNIP_DELTA & "-" (= value negation)
// 2. else, tlen>0 and pnext_pos_delta>0 and seq_len>0 tlen is stored as SNIP_SPECIAL & tlen-pnext_pos_delta-seq_len
// 3. else, stored as is
void sam_seg_tlen_field (VBlockSAM *vb, 
                         const char *tlen, unsigned tlen_len, // option 1
                         int64_t tlen_value, // option 2
                         PosType pnext_pos_delta, int32_t cigar_seq_len)
{
    Context *ctx = &vb->contexts[SAM_TLEN];

    if (tlen) { // option 1
        ASSSEG (tlen_len, tlen, "%s: empty TLEN", global_cmd);

        bool is_int = str_get_int (tlen, tlen_len, &tlen_value); // note: tlen_value remains 0 if not a valid integer
        ASSSEG (is_int, tlen, "%s: expecting TLEN to be an integer, but found \"%.*s\"", global_cmd, tlen_len, tlen);
    }

    unsigned add_bytes = IS_BAM ? sizeof (uint32_t) : tlen_len + 1;

    // case 1
    if (tlen_value && tlen_value == -ctx->last_value.i) {
        char snip_delta[2] = { SNIP_SELF_DELTA, '-' };
        seg_by_ctx ((VBlockP)vb, snip_delta, 2, ctx, add_bytes, NULL);
    }
    // case 2:
    else if (tlen_value > 0 && pnext_pos_delta > 0 && cigar_seq_len > 0) {
        char tlen_by_calc[30] = { SNIP_SPECIAL, SAM_SPECIAL_TLEN };
        unsigned tlen_by_calc_len = str_int (tlen_value - pnext_pos_delta - (int64_t)cigar_seq_len, &tlen_by_calc[2]);
        seg_by_ctx ((VBlockP)vb, tlen_by_calc, tlen_by_calc_len + 2, ctx, add_bytes, NULL);
    }
    // case default: add as is (option 1)
    else if (tlen)
        seg_by_ctx ((VBlockP)vb, tlen, tlen_len, ctx, add_bytes, NULL);

    // case default: add as is (option 2)
    else {
        char snip[20];
        unsigned snip_len = str_int (tlen_value, snip);
        seg_by_ctx ((VBlockP)vb, snip, snip_len, ctx, add_bytes, NULL);
    }

    ctx->last_value.i = tlen_value;
}

// returns length of string ending with separator, or -1 if separator was not found
static inline int sam_seg_get_next_subitem (const char *str, int str_len, char separator)
{
    for (int i=0; i < str_len; i++) {
        if (str[i] == separator) return i;
        if (str[i] == ',' || str[i] == ';') return -1; // wrong separator encountered
    }
    return -1;
}

#define DO_SSF(ssf,sep) \
        ssf = &field[i]; \
        ssf##_len = sam_seg_get_next_subitem (&field[i], field_len-i, sep); \
        if (ssf##_len == -1) goto error; /* bad format */ \
        i += ssf##_len + 1; /* skip snip and separator */        

#define DEC_SSF(ssf) const char *ssf; \
                     int ssf##_len;

// Creates a bitmap from seq data - exactly one bit per base that is mapped to the reference (e.g. not for INSERT bases)
// - Normal SEQ: tracking CIGAR, we compare the sequence to the reference, indicating in a SAM_SQBITMAP whether this
//   base in the same as the reference or not. In case of REF_INTERNAL, if the base is not already in the reference, we add it.
//   bases that differ from the reference are stored in SAM_NONREF
// - Edge case: no POS (i.e. unaligned read) - we just store the sequence in SAM_NONREF
// - Edge case: no CIGAR (it is "*") - we just treat it as an M and compare to the reference
// - Edge case: no SEQ (it is "*") - we '*' in SAM_NONREF and indicate "different from reference" in the bitmap. We store a
//   single entry, regardless of the number of entries indicated by CIGAR
//
// Best explanation of CIGAR operations: https://davetang.org/wiki/tiki-index.php?page=SAM
void sam_seg_seq_field (VBlockSAM *vb, const char *seq, uint32_t seq_len, PosType pos, const char *cigar, 
                        unsigned recursion_level, uint32_t level_0_seq_len, const char *level_0_cigar, unsigned add_bytes)
{
    START_TIMER;

    Context *bitmap_ctx = &vb->contexts[SAM_SQBITMAP];
    Context *nonref_ctx = &vb->contexts[SAM_NONREF];

    ASSERT0 (recursion_level < 4, "Error in sam_seg_seq_field: excess recursion"); // this would mean a read of about 4M bases... in 2020, this looks unlikely

    bitmap_ctx->txt_len += add_bytes; // byte counts for --show-sections

    // for unaligned lines, if we have refhash loaded, use the aligner instead of CIGAR-based segmenting
    if (!pos && flag_ref_use_aligner) {
        aligner_seg_seq ((VBlockP)vb, bitmap_ctx, seq, seq_len);
        goto align_nonref_local;
    }

    BitArray *bitmap = buf_get_bitarray (&bitmap_ctx->local);

    if (!recursion_level) {
        // allocate bitmap - provide name only if buffer is not allocated, to avoid re-writing param which would overwrite num_of_bits that overlays it
        buf_alloc (vb, &bitmap_ctx->local, MAX (bitmap_ctx->local.len + roundup_bits2bytes64 (seq_len), vb->lines.len * (seq_len+5) / 8), CTX_GROWTH, 
                buf_is_allocated (&bitmap_ctx->local) ? NULL : "context->local", 0); 
        
        buf_alloc (vb, &nonref_ctx->local, MAX (nonref_ctx->local.len + seq_len + 3, vb->lines.len * seq_len / 4), CTX_GROWTH, "context->local", nonref_ctx->did_i); 

        buf_extend_bits (&bitmap_ctx->local, vb->ref_and_seq_consumed);
    }

    // we can't compare to the reference if POS is 0: we store the seqeuence in SEQ_NOREF without an indication in the bitmap
    if (!pos) {
        buf_add (&nonref_ctx->local, seq, seq_len);
        goto align_nonref_local; 
    }

    if (seq[0] == '*') goto done; // we already handled a missing seq (SEQ="*") by adding a '-' to CIGAR - no data added here

    RefLock lock;
    Range *range = ref_seg_get_locked_range ((VBlockP)vb, pos, vb->ref_consumed, seq, &lock);

    // Cases where we don't consider the refernce and just copy the seq as-is
    if (!range || // 1. (denovo:) case reference range is NULL as the hash entry for this range is unfortunately already occupied by another range
                  // 2. (loaded:) case contig doesn't exist in the reference
        (cigar[0] == '*' && cigar[1] == 0)) { // 3. in case there's no CIGAR. The sequence is not aligned to the reference even if we have RNAME and POS (and its length can exceed the reference contig)

        buf_add (&nonref_ctx->local, seq, seq_len);
        
        bit_array_clear_region (bitmap, bitmap_ctx->next_local, vb->ref_and_seq_consumed); // note: vb->ref_and_seq_consumed==0 if cigar="*"
        bitmap_ctx->next_local += vb->ref_and_seq_consumed;

        random_access_update_last_pos ((VBlockP)vb, pos + vb->ref_consumed - 1);

        if (range) ref_unlock (lock);
        
        // note: in case of a missing range (which can be the first range in this seq, or a subsequent range), we zero the entire remaining bitmap.
        // this is because, absent a reference, we don't know how much ref is consumed by this missing range.
        goto align_nonref_local; 
    }    

    uint32_t pos_index     = pos - range->first_pos;
    uint32_t next_ref      = pos_index;
    const char *next_cigar = cigar;

    uint32_t i=0;
    int subcigar_len=0;
    char cigar_op;

    uint32_t ref_len_this_level = (flag_reference == REF_INTERNAL ? MIN (vb->ref_consumed, range->last_pos - pos + 1)
                                                                  : vb->ref_consumed); // possibly going around the end of the chromosome in case of a circular chromosome                                   

    uint32_t range_len = (range->last_pos - range->first_pos + 1);
    
    while (i < seq_len || next_ref < pos_index + ref_len_this_level) {

        ASSERT0 (i <= seq_len && next_ref <= pos_index + ref_len_this_level, "Error in sam_seg_seq_field: i or next_ref are out of range");

        subcigar_len = strtod (next_cigar, (char **)&next_cigar); // get number and advance next_cigar
        
        cigar_op = *(next_cigar++);

        if (cigar_op == 'M' || cigar_op == '=' || cigar_op == 'X') { // alignment match or sequence match or mismatch

            ASSERT (subcigar_len > 0 && subcigar_len <= (seq_len - i), 
                    "Error in sam_seg_seq_field: CIGAR %s implies seq_len longer than actual seq_len=%u (recursion_level=%u level0: cigar=%s seq_len=%u)", 
                    cigar, seq_len, recursion_level, level_0_cigar, level_0_seq_len);

            uint32_t bit_i = bitmap_ctx->next_local; // copy to automatic variable for performance
            uint32_t start_i = i;
            while (subcigar_len && next_ref < pos_index + ref_len_this_level) {

                // when we have an X we don't enter it into our internal ref, and we wait for a read with a = or M for that site,
                // as we assume that on average, more reads will have the reference base, leading to better compression
            
                bool normal_base = IS_NUCLEOTIDE (seq[i]);

                // circle around to beginning of chrom if out of range (can only happen with external reference, expected only with circular chromosomes) 
                uint32_t actual_next_ref = next_ref % range_len; 

                // case: we have not yet set a value for this site - we set it now. note: in ZIP, is_set means that the site
                // will be needed for pizzing. With REF_INTERNAL, this is equivalent to saying we have set the ref value for the site
                if (flag_reference == REF_INTERNAL && range && normal_base 
                    && !ref_is_nucleotide_set (range, actual_next_ref)) { 
                    
                    ref_set_nucleotide (range, actual_next_ref, seq[i]);
                    bit_array_set (&range->is_set, actual_next_ref); // we will need this ref to reconstruct
                    bit_array_set (bitmap, bit_i); bit_i++; // cannot increment inside the macro
                }

                // case our seq is identical to the reference at this site
                else if (range && normal_base && seq[i] == ref_get_nucleotide (range, actual_next_ref)) {
                    bit_array_set (bitmap, bit_i); bit_i++;

                    if (flag_reference == REF_EXT_STORE)
                        bit_array_set (&range->is_set, actual_next_ref); // we will need this ref to reconstruct
                }
                
                // case: ref is set to a different value - we store our value in nonref_ctx
                else {
                    NEXTENT (char, nonref_ctx->local) = seq[i];
                    bit_array_clear (bitmap, bit_i); bit_i++;
                } 

                subcigar_len--;
                next_ref++;
                i++;
            }
            vb->ref_and_seq_consumed -= (i - start_i); // update in case a range in a subsequent recursion level is missing and we need to clear the bitmap
            bitmap_ctx->next_local = bit_i;
        } // end if 'M', '=', 'X'

        // for Insertion or Soft clipping - this SEQ segment doesn't align with the reference - we leave it as is 
        else if (cigar_op == 'I' || cigar_op == 'S') {

            ASSSEG (subcigar_len > 0 && subcigar_len <= (seq_len - i), seq,
                    "Error in sam_seg_seq_field: CIGAR %s implies seq_len longer than actual seq_len=%u", cigar, seq_len);

            buf_add (&nonref_ctx->local, &seq[i], subcigar_len);
            i += subcigar_len;
            subcigar_len = 0;
        }

        // for Deletion or Skipping - we move the next_ref ahead
        else if (cigar_op == 'D' || cigar_op == 'N') {
            unsigned ref_consumed = (flag_reference == REF_INTERNAL ? MIN (subcigar_len, range_len - next_ref)
                                                                    : subcigar_len);
            next_ref     += ref_consumed;
            subcigar_len -= ref_consumed;
        }

        // Hard clippping (H) or padding (P) we do nothing
        else if (cigar_op == 'H' || cigar_op == 'P') {}

        else {
            ASSSEG (cigar_op, vb->last_cigar, "Error in sam_seg_seq_field: End of CIGAR reached but we still have %u reference and %u sequence bases to consume"
                    "(cigar=%s pos=%"PRId64" recursion_level=%u level_0_cigar=%s level_0_seq_len=%u) (vb->ref_consumed=%d next_ref=%u pos_index=%u ref_len_this_level=%u subcigar_len=%u range=[%.*s %"PRId64"-%"PRId64"])",
                    pos_index + ref_len_this_level - next_ref, seq_len-i,   cigar, pos, recursion_level, level_0_cigar, level_0_seq_len,
                    vb->ref_consumed, next_ref, pos_index, ref_len_this_level, subcigar_len, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos);        

            ASSSEG (false, vb->last_cigar, "Error in sam_seg_seq_field: Invalid CIGAR op: '%c' (ASCII %u)", cigar_op, cigar_op);        
        }

        // case: we're at the end of the reference AND we want more of it
        if (next_ref == pos_index + ref_len_this_level && subcigar_len) break;
    }

    if (range) ref_unlock (lock);       

    uint32_t this_seq_last_pos = pos + (next_ref - pos_index) - 1;

    // in REF_INTERNAL, the sequence can flow over to the next range as each range is 1M bases. this cannot happen
    // in REF_EXTERNAL as each range is the entire contig
    ASSERT (flag_reference == REF_INTERNAL || i == seq_len, "Error in sam_seg_seq_field: expecting i(%u) == seq_len(%u) pos=%"PRId64" range=[%.*s %"PRId64"-%"PRId64"] (cigar=%s recursion_level=%u level0: cigar=%s seq_len=%u)", 
            i, seq_len, pos, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos, cigar, recursion_level, level_0_cigar, level_0_seq_len);

    // case: we have reached the end of the current reference range, but we still have sequence left - 
    // call recursively with remaining sequence and next reference range 
    if (i < seq_len) {

        ASSSEG (this_seq_last_pos <= MAX_POS, cigar, "%s: Error: POS=%"PRId64" and the consumed reference implied by CIGAR=\"%s\", exceeding MAX_POS=%"PRId64
                " (next_ref=%u pos_index=%u ref_len_this_level=%u subcigar_len=%u range=[%.*s %"PRId64"-%"PRId64"])",
                global_cmd, pos, cigar, MAX_POS, next_ref, pos_index, ref_len_this_level, subcigar_len, 
                range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos);

        vb->ref_consumed -= ref_len_this_level;

        char updated_cigar[100];
        if (subcigar_len) sprintf (updated_cigar, "%u%c%s", subcigar_len, cigar_op, next_cigar);

        sam_seg_seq_field (vb, seq + i, seq_len - i, range->last_pos + 1, subcigar_len ? updated_cigar : next_cigar, recursion_level + 1, level_0_seq_len, level_0_cigar, 0);
    }
    else { // update RA of the VB with last pos of this line as implied by the CIGAR string
        if (this_seq_last_pos <= range->last_pos) // always the case in INTERNAL and non-circular EXTERNAL 
            random_access_update_last_pos ((VBlockP)vb, this_seq_last_pos);

        else  // we circled back to the beginning for the chromosome - i.e. this VB RA is the entire chromosome
            random_access_update_to_entire_chrom ((VBlockP)vb, range->first_pos, range->last_pos);
    }
align_nonref_local: {
    // we align nonref_ctx->local to a 4-character boundary. this is because CODEC_ACGT squeezes every 4 characters into a byte,
    // before compressing it with LZMA. In sorted SAM, we want subsequent identical sequences to have the same byte alignment
    // so that LZMA can catch their identicality.
    uint64_t add_chars = (4 - (nonref_ctx->local.len & 3)) & 3;
    if (add_chars) buf_add (&nonref_ctx->local, "AAA", add_chars); // add 1 to 3 As
}
done:
    COPY_TIMER (sam_seg_seq_field);
}

static void sam_seg_SA_or_OA_field (VBlockSAM *vb, DictId subfield_dict_id, 
                                    const char *field, unsigned field_len, const char *field_name)
{
    // OA and SA format is: (rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+ . in OA - NM is optional (but its , is not)
    // Example SA:Z:chr13,52863337,-,56S25M70S,0,0;chr6,145915118,+,97S24M30S,0,0;chr18,64524943,-,13S22M116S,0,0;chr7,56198174,-,20M131S,0,0;chr7,87594501,+,34S20M97S,0,0;chr4,12193416,+,58S19M74S,0,0;
    // See: https://samtools.github.io/hts-specs/SAMtags.pdf
    static const Container container_SA_OA = {
        .repeats     = 0, 
        .num_items   = 6, 
        .flags       = 0,
        .repsep      = {0,0},
        .items       = { { .dict_id = {.id="@RNAME" }, .seperator = {','}, .did_i = DID_I_NONE},  // we don't mix with primary as primary is often sorted, and mixing will ruin its b250 compression
                         { .dict_id = {.id="@POS"   }, .seperator = {','}, .did_i = DID_I_NONE},  // we don't mix with primary as these are local-stored random numbers anyway - no advantage for mixing, and it would obscure the stats
                         { .dict_id = {.id="@STRAND"}, .seperator = {','}, .did_i = DID_I_NONE},
                         { .dict_id = {.id={'C'&0x3f,'I','G','A','R'}}, .seperator = {','}, .did_i = DID_I_NONE}, // we mix with primary - CIGAR tends to be a rather large dictionary, so better not have two copies of it
                         { .dict_id = {.id="@MAPQ"  }, .seperator = {','}, .did_i = DID_I_NONE},  // we don't mix with primary as primary often has a small number of values, and mixing will ruin its b250 compression
                         { .dict_id = {.id="NM:i"   }, .seperator = {';'}, .did_i = DID_I_NONE} } // we mix together with the NM option field
    };

    DEC_SSF(rname); DEC_SSF(pos); DEC_SSF(strand); DEC_SSF(cigar); DEC_SSF(mapq); DEC_SSF(nm); 

    Container sa_oa = container_SA_OA;

    for (uint32_t i=0; i < field_len; sa_oa.repeats++) {

        ASSSEG (sa_oa.repeats <= CONTAINER_MAX_REPEATS, field, "Error in sam_seg_SA_or_OA_field - exceeded maximum repeats allowed (%lu) while parsing %s",
                CONTAINER_MAX_REPEATS, err_dict_id (subfield_dict_id));

        DO_SSF (rname,  ','); // these also do sanity checks
        DO_SSF (pos,    ','); 
        DO_SSF (strand, ','); 
        DO_SSF (cigar,  ','); 
        DO_SSF (mapq,   ','); 
        DO_SSF (nm,     ';'); 

        // sanity checks before adding to any dictionary
        if (strand_len != 1 || (strand[0] != '+' && strand[0] != '-')) goto error; // invalid format
        
        PosType pos_value = seg_scan_pos_snip ((VBlockP)vb, pos, pos_len, true);
        if (pos_value < 0) goto error;

        seg_by_dict_id (vb, rname,  rname_len,  container_SA_OA.items[0].dict_id, 1 + rname_len);
        seg_by_dict_id (vb, strand, strand_len, container_SA_OA.items[2].dict_id, 1 + strand_len);
        seg_by_dict_id (vb, cigar,  cigar_len,  container_SA_OA.items[3].dict_id, 1 + cigar_len);
        seg_by_dict_id (vb, mapq,   mapq_len,   container_SA_OA.items[4].dict_id, 1 + mapq_len);
        seg_by_dict_id (vb, nm,     nm_len,     container_SA_OA.items[5].dict_id, 1 + nm_len);
        
        Context *pos_ctx = mtf_get_ctx (vb, container_SA_OA.items[1].dict_id);
        pos_ctx->ltype   = LT_UINT32;
        seg_add_to_local_uint32 ((VBlockP)vb, pos_ctx, pos_value, 1 + pos_len);
    }

    seg_container_by_dict_id (vb, subfield_dict_id, &sa_oa, 1 /* 1 for \t in SAM and \0 in BAM */);
    
    return;

error:
    // if the error occurred on on the first repeat - this file probably has a different
    // format - we just store as a normal subfield
    // if it occurred on the 2nd+ subfield, after the 1st one was fine - we reject the file
    ASSSEG (!sa_oa.repeats, field, "Invalid format in repeat #%u of field %s. snip: %.*s",
            sa_oa.repeats+1, err_dict_id (subfield_dict_id), field_len, field);

    seg_by_dict_id (vb, field, sizeof(field_len), subfield_dict_id, field_len + 1 /* 1 for \t in SAM and \0 in BAM */); 
}

static void sam_seg_XA_field (VBlockSAM *vb, const char *field, unsigned field_len)
{
    // XA format is: (chr,pos,CIGAR,NM;)*  pos starts with +- which is strand
    // Example XA:Z:chr9,-60942781,150M,0;chr9,-42212061,150M,0;chr9,-61218415,150M,0;chr9,+66963977,150M,1;
    // See: http://bio-bwa.sourceforge.net/bwa.shtml
    static const Container container_XA = {
        .repeats     = 0, 
        .num_items   = 5, 
        .flags       = 0,
        .repsep      = {0,0},
        .items       = { { .dict_id = {.id="@RNAME"  }, .seperator = {','}, .did_i = DID_I_NONE },
                         { .dict_id = {.id="@STRAND" }, .seperator = { 0 }, .did_i = DID_I_NONE },
                         { .dict_id = {.id="@POS"    }, .seperator = {','}, .did_i = DID_I_NONE },
                         { .dict_id = {.id={'C'&0x3f,'I','G','A','R'}}, .seperator = {','}, .did_i = DID_I_NONE},
                         { .dict_id = {.id="NM:i"    }, .seperator = {';'}, .did_i = DID_I_NONE } }     
    };

    Container xa = container_XA;

    DEC_SSF(rname); DEC_SSF(pos); DEC_SSF(cigar); DEC_SSF(nm); 

    for (uint32_t i=0; i < field_len; xa.repeats++) {

        ASSSEG (xa.repeats <= CONTAINER_MAX_REPEATS, field, "Error in sam_seg_XA_field - exceeded maximum repeats allowed (%lu) while parsing XA",
                CONTAINER_MAX_REPEATS);

        DO_SSF (rname,  ','); 
        DO_SSF (pos,    ','); 
        DO_SSF (cigar,  ','); 
        DO_SSF (nm,     ';'); 

        // sanity checks before adding to any dictionary
        if (pos_len < 2 || (pos[0] != '+' && pos[0] != '-')) goto error; // invalid format - expecting pos to begin with the strand

        PosType pos_value = seg_scan_pos_snip ((VBlockP)vb, &pos[1], pos_len-1, true);
        if (pos_value < 0) goto error;

        seg_by_dict_id (vb, rname,  rname_len, container_XA.items[0].dict_id, 1 + rname_len);
        seg_by_dict_id (vb, pos,    1,         container_XA.items[1].dict_id, 1); // strand is first character of pos
        seg_by_dict_id (vb, cigar,  cigar_len, container_XA.items[3].dict_id, 1 + cigar_len);
        seg_by_dict_id (vb, nm,     nm_len,    container_XA.items[4].dict_id, 1 + nm_len);
        
        Context *pos_ctx = mtf_get_ctx (vb, container_XA.items[2].dict_id);
        pos_ctx->ltype  = LT_UINT32;

        seg_add_to_local_uint32 ((VBlockP)vb, pos_ctx, pos_value, pos_len); // +1 for seperator, -1 for strand
    }

    seg_container_by_dict_id (vb, dict_id_OPTION_XA, &xa, 1 /* 1 for \t in SAM and \0 in BAM */);
    return;

error:
    // if the error occurred on on the first repeat - this file probably has a different
    // format - we just store as a normal subfield
    // if it occurred on the 2nd+ subfield, after the 1st one was fine - we reject the file
    ASSSEG (!xa.repeats, field, "Invalid format in repeat #%u of field XA. snip: %.*s", xa.repeats+1, field_len, field);

    seg_by_dict_id (vb, field, field_len, dict_id_OPTION_XA, field_len + 1 /* 1 for \t in SAM and \0 in BAM */); 
}

uint32_t sam_seg_get_seq_len_by_MD_field (const char *md_str, unsigned md_str_len)
{
    uint32_t result=0, curr_num=0;

    for (unsigned i=0; i < md_str_len; i++) {   
        if (IS_DIGIT (md_str[i])) 
            curr_num = curr_num * 10 + (md_str[i] - '0');

        else {
            result += curr_num + 1; // number terminates here + one character
            curr_num = 0;
        }
    }

    result += curr_num; // in case the string ends with a number

    return result;
}

// in the case where sequence length as calculated from the MD is the same as that calculated
// from the CIGAR/SEQ/QUAL (note: this is required by the SAM spec but nevertheless genozip doesn't require it):
// MD is shortened to replace the last number with a *, since it can be calculated knowing the length. The result is that
// multiple MD values collapse to one, e.g. "MD:Z:119C30" and "MD:Z:119C31" both become "MD:Z:119C*" hence improving compression.
// In the case where the MD is simply a number "151" and drop it altogether and keep just an empty string.
static inline bool sam_seg_get_shortened_MD (const char *md_str, unsigned md_str_len, uint32_t seq_len,
                                             char *new_md_str, unsigned *new_md_str_len)
{
    uint32_t seq_len_by_md = sam_seg_get_seq_len_by_MD_field (md_str, md_str_len);

    if (seq_len_by_md != seq_len) return false;  // MD string doesn't comply with SAM spec and is therefore not changed
    
    // case - MD ends with a number eg "119C31" - we replace it with prefix+"119C". if its all digits then just prefix
    if (IS_DIGIT (md_str[md_str_len-1])) {

        int i=md_str_len-1; for (; i>=0; i--)
            if (!IS_DIGIT (md_str[i])) break;

        new_md_str[0] = SNIP_SPECIAL;
        new_md_str[1] = SAM_SPECIAL_MD;
        if (i >= 0) memcpy (&new_md_str[2], md_str, i+1);
        
        *new_md_str_len = i+3;
        return true;
    }

    return false; // MD doesn't end with a number and is hence unchanged (this normally doesn't occur as the MD would finish with 0)
}

// AS and XS are values (at least as set by BWA) at most the seq_len, and AS is often equal to it. we modify
// it to be new_value=(value-seq_len) 
static inline void sam_seg_AS_field (VBlockSAM *vb, ZipDataLineSAM *dl, DictId dict_id, 
                                     const char *snip, unsigned snip_len, unsigned add_bytes)
{
    bool positive_delta = true;

    // verify that its a unsigned number
    for (unsigned i=0; i < snip_len; i++)
        if (!IS_DIGIT (snip[i])) positive_delta = false;

    int32_t as;
    if (positive_delta) {
        as = atoi (snip); // type i is signed 32 bit by SAM specification
        if (dl->seq_len < as) positive_delta=false;
    }

    // if possible, store a special snip with the positive delta
    if (positive_delta) {
        char new_snip[20] = { SNIP_SPECIAL, SAM_SPECIAL_AS };
        unsigned delta_len = str_int (dl->seq_len-as, &new_snip[2]);

        seg_by_dict_id (vb, new_snip, delta_len+2, dict_id, add_bytes); 
    }

    // not possible - just store unmodified
    else
        seg_by_dict_id (vb, snip, snip_len, dict_id, add_bytes); 
}

// optimization for Ion Torrent flow signal (ZM) - negative values become zero, positives are rounded to the nearest 10
static void sam_optimize_ZM (const char **snip, unsigned *snip_len, char *new_str)
{
    char *after;
    int number = strtoul (*snip, &after, 10);

    if ((unsigned)(after - *snip) > 0) {
        if (number >= 0) number = ((number + 5) / 10) * 10;
        else             number = 0;

        *snip_len = str_int (number, new_str);
        *snip = new_str;
    }    
}

static inline ContainerItemTransform sam_seg_optional_transform (char type)
{
    switch (type) {
        case 'c': return TRS_I8;
        case 'C': return TRS_U8;
        case 's': return TRS_LTEN_I16;
        case 'S': return TRS_LTEN_U16;
        case 'i': return TRS_LTEN_I32;
        case 'I': return TRS_LTEN_U32;
        case 'f': return TRS_SAM2BAM_FLOAT;
        default : return TRS_NONE;
    }
}

static inline unsigned sam_seg_optional_add_bytes (char type, unsigned value_len, bool is_bam)
{
    if (is_bam)
        switch (type) {
            case 'c': case 'C': case 'A': return 1;
            case 's': case 'S':           return 2;
            case 'i': case 'I': case 'f': return 4;
            case 'Z': case 'H':           return value_len + 1; // +1 for \0
            default : return 0;
        }
    else // SAM
        return value_len + 1; // +1 for \t
}

static inline char sam_seg_bam_type_to_sam_type (char type)
{
    return (type=='c' || type=='C' || type=='s' || type=='S' || type=='I') ? 'i' : type;
}

// process an optional subfield, that looks something like MX:Z:abcdefg. We use "MX" for the field name, and
// the data is abcdefg. The full name "MX:Z:" is stored as part of the OPTIONAL dictionary entry
static DictId sam_seg_optional_field (VBlockSAM *vb, ZipDataLineSAM *dl, bool is_bam, 
                                      const char *tag, char type, const char *value, unsigned value_len)
{
    char dict_name[4] = { tag[0], tag[1], ':', sam_seg_bam_type_to_sam_type (type) };
    DictId dict_id = sam_dict_id_optnl_sf (dict_id_make (dict_name, 4));

    unsigned add_bytes = sam_seg_optional_add_bytes (type, value_len, is_bam);

    if (dict_id.num == dict_id_OPTION_SA || dict_id.num == dict_id_OPTION_OA)
        sam_seg_SA_or_OA_field (vb, dict_id, value, value_len, dict_id.num == dict_id_OPTION_SA ? "SA" : "OA");

    else if (dict_id.num == dict_id_OPTION_XA) 
        sam_seg_XA_field (vb, value, value_len);

    // fields containing CIGAR format data
    else if (dict_id.num == dict_id_OPTION_MC || dict_id.num == dict_id_OPTION_OC) 
        seg_by_did_i (vb, value, value_len, SAM_CIGAR, add_bytes)

    // MD's logical length is normally the same as seq_len, we use this to optimize it.
    // In the common case that it is just a number equal the seq_len, we replace it with an empty string.
    else if (dict_id.num == dict_id_OPTION_MD) {
        // if MD value can be derived from the seq_len, we don't need to store - store just an empty string

#define MAX_SAM_MD_LEN 1000 // maximum length of MD that is shortened.
        char new_md[MAX_SAM_MD_LEN];
        unsigned new_md_len = 0;
        bool md_is_special  = (value_len-2 <= MAX_SAM_MD_LEN);

        if (md_is_special) 
            md_is_special = sam_seg_get_shortened_MD (value, value_len, dl->seq_len, new_md, &new_md_len);

        // not sure which of these two is better....
        seg_by_dict_id (vb,                                 
                        md_is_special ? new_md : value, 
                        md_is_special ? new_md_len : value_len,
                        dict_id, add_bytes);
    }

    // BD and BI set by older versions of GATK's BQSR is expected to be seq_len (seen empircally, documentation is lacking)
    else if ((dict_id.num == dict_id_OPTION_BD || dict_id.num == dict_id_OPTION_BI) && value_len == dl->seq_len) {
        
        bool is_bi = (dict_id.num == dict_id_OPTION_BI);
        dl->bdbi_data_start[is_bi] = value - vb->txt_data.data;

        Context *ctx = mtf_get_ctx (vb, dict_id_OPTION_BD_BI);
        ctx->txt_len += add_bytes; 
        ctx->ltype   = LT_SEQUENCE;

        if (!dl->bdbi_data_start[!is_bi]) // the first of BD and BI increments local.len, so it is incremented even if just one of BD/BI appears
            ctx->local.len += value_len * 2;

        // we can't use local for singletons in BD or BI as next_local is used by sam_piz_special_BD_BI to point into BD_BI
        Context *this_ctx  = mtf_get_ctx (vb, dict_id);
        this_ctx->inst     = CTX_INST_NO_STONS; 
        this_ctx->st_did_i = ctx->did_i; 

        const char special_snip[2] = { SNIP_SPECIAL, SAM_SPECIAL_BDBI };
        seg_by_dict_id (vb, special_snip, 2, dict_id, 0);
    }

    // AS is a value (at least as set by BWA) at most the seq_len, and often equal to it. we modify
    // it to be new_AS=(AS-seq_len) 
    else if (dict_id.num == dict_id_OPTION_AS) 
        sam_seg_AS_field (vb, dl, dict_id, value, value_len, add_bytes);
    
    // mc:i: (output of bamsormadup? - mc in small letters) appears to a pos value usually close to POS.
    // we encode as a delta.
    else if (dict_id.num == dict_id_OPTION_mc) {
        uint8_t mc_did_i = mtf_get_ctx (vb, dict_id)->did_i;

        seg_pos_field ((VBlockP)vb, mc_did_i, SAM_POS, true, value, value_len, 0, add_bytes);
    }

    // E2 - SEQ data (note: E2 doesn't have a context - it shares with SEQ)
    else if (dict_id.num == dict_id_OPTION_E2) { 
        ASSERT (value_len == dl->seq_len, 
                "Error in %s: Expecting E2 data to be of length %u as indicated by CIGAR, but it is %u. E2=%.*s",
                txt_name, dl->seq_len, value_len, value_len, value);

        PosType this_pos = vb->contexts[SAM_POS].last_value.i;
        sam_seg_seq_field (vb, (char *)value, value_len, this_pos, vb->last_cigar, 0, value_len, // remove const bc SEQ data is actually going to be modified
                           vb->last_cigar, add_bytes); 
    }

    // U2 - QUAL data (note: U2 doesn't have a context - it shares with QUAL)
    else if (dict_id.num == dict_id_OPTION_U2) {
        ASSERT (value_len == dl->seq_len, 
                "Error in %s: Expecting U2 data to be of length %u as indicated by CIGAR, but it is %u. E2=%.*s",
                txt_name, dl->seq_len, value_len, value_len, value);

        dl->u2_data_start = value - vb->txt_data.data;
        dl->u2_data_len   = value_len;
        vb->contexts[SAM_QUAL].txt_len   += add_bytes;
        vb->contexts[SAM_QUAL].local.len += value_len;
    }

    // Numeric array array
    else if (type == 'B') {
        SegOptimize optimize = NULL;

        if (flag_optimize_ZM && dict_id.num == dict_id_OPTION_ZM && value_len > 3 && value[0] == 's')  // XM:B:s,
            optimize = sam_optimize_ZM;

        ContainerItemTransform transform = sam_seg_optional_transform (value[0]); // instructions on how to transform array items if reconstructing as BAM (value[0] is the subtype of the array)

        uint32_t repeats = seg_array_field ((VBlockP)vb, dict_id, value, value_len, !IS_BAM, transform, optimize);

        // add bytes here in case of BAM - all to main field
        if (IS_BAM) {
            unsigned add_bytes_per_repeat = sam_seg_optional_add_bytes (value[0], 0, true);
            mtf_get_ctx (vb, dict_id)->txt_len += add_bytes_per_repeat * (repeats-1) + 4 // don't include 1st repeat which is the type, +4 for count
                                                  + 1; // type (but not tag and 'B' which are part of OPTIONAL)
        }
    }

    // All other subfields - have their own dictionary
    else        
        seg_by_dict_id (vb, value, value_len, dict_id, add_bytes); 

    return dict_id;
}

const char *sam_get_one_optional (VBlockSAM *vb, const char *next_field, int32_t len, char *separator_p, bool *has_13, 
                                  const char **tag, char *type, const char **value, unsigned *value_len) // out
{
    unsigned field_len;
    const char *field_start;

    char separator;
    GET_MAYBE_LAST_ITEM ("OPTIONAL-subfield"); 

    ASSSEG0 (field_len, field_start, "Error: line invalidly ends with a tab");

    ASSSEG (field_len >= 6 && field_start[2] == ':' && field_start[4] == ':', field_start, "Error: invalid optional field format: %.*s",
            field_len, field_start);

    *tag         = field_start;
    *type        = field_start[3];
    *value       = field_start + 5;
    *value_len   = field_len - 5;
    *separator_p = separator;

    return next_field;
}

const char *sam_seg_optional_all (VBlockSAM *vb, ZipDataLineSAM *dl, const char *next_field,
                                  int32_t len, bool *has_13, char separator, // sam only
                                  const char *after_field) // bam only 
{
    const bool is_bam = IS_BAM;
    Container con = { .repeats=1 };
    char prefixes[MAX_SUBFIELDS * 6 + 2]; // each name is 5 characters per SAM specification, eg "MC:Z:" followed by SNIP_CONTAINER ; +2 for the initial SNIP_CONTAINER
    prefixes[0] = prefixes[1] = SNIP_CONTAINER; // initial SNIP_CONTAINER follow by seperator of empty Container-wide prefix
    unsigned prefixes_len=2;
    const char *value, *tag;
    char type;
    unsigned value_len;

    while (is_bam ? (next_field < after_field) : (separator != '\n')) {
    
        next_field = is_bam ? bam_get_one_optional (vb, next_field,                          &tag, &type, &value, &value_len)
                            : sam_get_one_optional (vb, next_field, len, &separator, has_13, &tag, &type, &value, &value_len);

        con.items[con.num_items++] = (ContainerItem) {
            .dict_id      = sam_seg_optional_field (vb, dl, is_bam, tag, type, value, value_len),
            .transform    = sam_seg_optional_transform (type), // how to transform the field if reconstructing to BAM
            .seperator[0] = '\t',
            .seperator[1] = 0,
            .did_i        = DID_I_NONE, // seg always puts NONE, PIZ changes it
        };

        ASSSEG (con.num_items <= MAX_SUBFIELDS, value, "Error: too many optional fields, limit is %u", MAX_SUBFIELDS);

        // in the optional field prefix (unlike array type), all integer types become 'i'.
        char prefix_type = sam_seg_bam_type_to_sam_type (type);

        char prefix[6] = { tag[0], tag[1], ':', prefix_type, ':', SNIP_CONTAINER}; 
        memcpy (&prefixes[prefixes_len], prefix, 6);
        prefixes_len += 6;

        if (is_bam) buf_free (&vb->textual_opt);
    }

    if (con.num_items) {
        con.items[con.num_items-1].seperator[0] = 0; // last Optional field has no tab
        seg_container_by_ctx ((VBlockP)vb, &vb->contexts[SAM_OPTIONAL], &con, prefixes, prefixes_len, (is_bam ? 3 : 5) * con.num_items); // account for : SAM: "MX:i:" BAM: "MXi"
    }
    else
        seg_by_did_i (vb, NULL, 0, SAM_OPTIONAL, 0); // NULL means MISSING Container item - will cause deletion of previous separator (\t)

    return next_field;        
}

static void sam_seg_cigar_field (VBlockSAM *vb, ZipDataLineSAM *dl, unsigned last_cigar_len,
                                 const char *seq,  uint32_t seq_data_len, 
                                 const char *qual, uint32_t qual_data_len)
{
    bool qual_is_available = qual_data_len != 1 || *qual != '*';
    bool seq_is_available  = seq_data_len != 1  || *seq  != '*';

    ASSSEG (!(seq_is_available && *seq=='*'), seq, "seq=%.*s (seq_len=%u), but expecting a missing seq to be \"*\" only (1 character)", seq_data_len, seq, seq_data_len);

    char cigar_snip[last_cigar_len + 50];
    cigar_snip[0] = SNIP_SPECIAL;
    cigar_snip[1] = SAM_SPECIAL_CIGAR;
    unsigned cigar_snip_len=2;

    // case: SEQ is "*" - we add a '-' to the CIGAR
    if (!seq_is_available) cigar_snip[cigar_snip_len++] = '-';

    // case: CIGAR is "*" - we get the dl->seq_len directly from SEQ or QUAL, and add the length to CIGAR eg "151*"
    if (!dl->seq_len) { // CIGAR is not available
        ASSSEG (!seq_data_len || !qual_is_available || seq_data_len==dl->qual_data_len, seq,
                "Bad line: SEQ length is %u, QUAL length is %u, unexpectedly differ. SEQ=%.*s QUAL=%.*s", 
                seq_data_len, dl->qual_data_len, seq_data_len, seq, dl->qual_data_len, qual);    

        dl->seq_len = MAX (seq_data_len, dl->qual_data_len); // one or both might be not available and hence =1

        cigar_snip_len += str_int (dl->seq_len, &cigar_snip[cigar_snip_len]);
    } 
    else { // CIGAR is available - just check the seq and qual lengths
        ASSSEG (!seq_is_available || seq_data_len == dl->seq_len, seq,
                "Bad line: according to CIGAR, expecting SEQ length to be %u but it is %u. SEQ=%.*s", 
                dl->seq_len, seq_data_len, seq_data_len, seq);

        ASSSEG (!qual_is_available || qual_data_len == dl->seq_len, qual,
                "Bad line: according to CIGAR, expecting QUAL length to be %u but it is %u. QUAL=%.*s", 
                dl->seq_len, dl->qual_data_len, dl->qual_data_len, qual);    
    }

    memcpy (&cigar_snip[cigar_snip_len], vb->last_cigar, last_cigar_len);
    cigar_snip_len += last_cigar_len;

    seg_by_did_i (vb, cigar_snip, cigar_snip_len, SAM_CIGAR, last_cigar_len+1); // +1 for \t
}

void sam_seg_qual_field (VBlockSAM *vb, ZipDataLineSAM *dl, const char *qual, uint32_t qual_data_len, unsigned add_bytes)
{
    dl->qual_data_start = qual - vb->txt_data.data;
    dl->qual_data_len   = qual_data_len;

    Context *qual_ctx = &vb->contexts[SAM_QUAL];
    qual_ctx->local.len += dl->qual_data_len;
    qual_ctx->txt_len   += add_bytes;
}

// test function called from main_load_reference -> txtfile_test_data: returns true if this line as pos=0 (i.e. unaligned)
bool sam_zip_is_unaligned_line (const char *line, int len)
{
    VBlock *vb = evb;
    const char *field_start, *next_field=line;
    char separator;
    unsigned field_len;

    GET_NEXT_ITEM ("QNAME");
    GET_NEXT_ITEM ("FLAG");
    GET_NEXT_ITEM ("RNAME");
    GET_NEXT_ITEM ("POS");

    return (field_len == 1 && *field_start == '0');
}

const char *sam_seg_txt_line (VBlock *vb_, const char *field_start_line, bool *has_13)     // index in vb->txt_data where this line starts
{
    VBlockSAM *vb = (VBlockSAM *)vb_;
    ZipDataLineSAM *dl = DATA_LINE (vb->line_i);

    const char *next_field=field_start_line, *field_start;
    unsigned field_len=0;
    char separator;

    int32_t len = AFTERENT (char, vb->txt_data) - field_start_line;

    // QNAME - We break down the QNAME into subfields separated by / and/or : - these are vendor-defined strings. Examples:
    // Illumina: <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> for example "A00488:61:HMLGNDSXX:4:1101:15374:1031" see here: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
    // PacBio BAM: {movieName}/{holeNumber}/{qStart}_{qEnd} see here: https://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
    GET_NEXT_ITEM ("QNAME");
    seg_compound_field ((VBlockP)vb, &vb->contexts[SAM_QNAME], field_start, field_len, false, 0, 1 /* \n */);

    SEG_NEXT_ITEM (SAM_FLAG);

    GET_NEXT_ITEM ("RNAME");
    seg_chrom_field (vb_, field_start, field_len);
    bool rname_is_missing = (*field_start == '*' && field_len == 1);

    GET_NEXT_ITEM ("POS");
    PosType this_pos = seg_pos_field (vb_, SAM_POS, SAM_POS, false, field_start, field_len, 0, field_len+1);
    ASSSEG (!(rname_is_missing && this_pos), field_start, "Error: RNAME=\"*\" - expecting POS to be 0 but it is %"PRId64, this_pos);

    random_access_update_pos (vb_, SAM_POS);

    SEG_NEXT_ITEM (SAM_MAPQ);

    // CIGAR - we wait to get more info from SEQ and QUAL
    GET_NEXT_ITEM ("CIGAR");
    sam_analyze_cigar (field_start, field_len, &dl->seq_len, &vb->ref_consumed, &vb->ref_and_seq_consumed);
    vb->last_cigar = field_start;
    unsigned last_cigar_len = field_len;
    ((char *)vb->last_cigar)[field_len] = 0; // null-terminate CIGAR string

    SEG_NEXT_ITEM (SAM_RNEXT);
    
    GET_NEXT_ITEM ("PNEXT");
    seg_pos_field (vb_, SAM_PNEXT, SAM_POS, false, field_start, field_len, 0, field_len+1);

    GET_NEXT_ITEM ("TLEN");
    sam_seg_tlen_field (vb, field_start, field_len, 0, vb->contexts[SAM_PNEXT].last_delta, dl->seq_len);

    GET_NEXT_ITEM ("SEQ");

    // calculate diff vs. reference (self or imported)
    sam_seg_seq_field (vb, field_start, field_len, this_pos, vb->last_cigar, 0, field_len, vb->last_cigar, field_len+1);
    const char *seq = field_start;
    uint32_t seq_data_len = field_len;

    GET_MAYBE_LAST_ITEM ("QUAL");
    sam_seg_qual_field (vb, dl, field_start, field_len, field_len + 1); 

    // finally we can seg CIGAR now
    sam_seg_cigar_field (vb, dl, last_cigar_len, seq, seq_data_len, field_start, field_len);
    
    // OPTIONAL fields - up to MAX_SUBFIELDS of them
    next_field = sam_seg_optional_all (vb, dl, next_field, len, has_13, separator, 0);

    SEG_EOL (SAM_EOL, false); /* last field accounted for \n */

    return next_field;
}
