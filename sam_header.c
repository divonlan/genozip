// ------------------------------------------------------------------
//   sam_header.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// ----------------------
// SAM / BAM header stuff
// ----------------------

#include "sam_private.h"
#include "file.h"
#include "random_access.h"
#include "reference.h"
#include "zfile.h"

#define HDRLEN evb->txt_data.len
#define HDRSKIP(n) if (HDRLEN < next + n) goto incomplete_header; next += n
#define HDR32 (next + 4 <= HDRLEN ? GET_UINT32 (&evb->txt_data.data[next]) : 0) ; if (HDRLEN < next + 4) goto incomplete_header; next += 4;

// call a callback for each SQ line (contig). Note: callback function is the same as ref_contigs_iterate
void sam_iterate_SQ_lines (const char *txt_header, // nul-terminated string
                           RefContigsIteratorCallback callback, void *callback_param)
{
    const char *line = txt_header;

    while (1) {
        line = strchr (line, '@');
        if (!line) break;

        if (line[1] == 'S' && line[2] == 'Q') { // this test will always be in the buffer - possible in the overflow area
            char *newline = strchr (line, '\n');
            *newline = 0;

            // line looks something like: @SQ\tSN:chr10\tLN:135534747 (no strict order between SN and LN and there could be more parameters)
            const char *chrom_name   = strstr (line, "SN:"); 
            const char *last_pos_str = strstr (line, "LN:"); 

            if (chrom_name && last_pos_str) {
                unsigned chrom_name_len = strcspn (&chrom_name[3], "\t\n\r");
                PosType last_pos = (PosType)strtoull (&last_pos_str[3], NULL, 10);

                ASSERT (last_pos <= MAX_POS_SAM, "Error: @SQ record in header contains LN:%"PRId64" which is beyond the maximum permitted by the SAM spec of %"PRId64,
                        last_pos, MAX_POS_SAM);

                callback (&chrom_name[3], chrom_name_len, last_pos, callback_param);
            }

            line = newline+1;

            *newline = '\n'; // restore
        }
        else 
            line++;
    }
}

// call a callback for each SQ line (contig). Note: callback function is the same as ref_contigs_iterate
static void bam_iterate_SQ_lines (const char *txt_header, // binary BAM header
                                  RefContigsIteratorCallback callback, void *callback_param)
{
    uint32_t next=0;

    // skip magic, l_text, text
    HDRSKIP(4);
    uint32_t l_text = HDR32;
    HDRSKIP(l_text);

    uint32_t n_ref = HDR32;

    for (uint32_t ref_i=0; ref_i < n_ref; ref_i++) {
        uint32_t l_name = HDR32;
        const char *name = &txt_header[next];
        HDRSKIP (l_name);
        uint32_t l_ref = HDR32;

        callback (name, l_name-1, l_ref, callback_param);
    }
    return;

incomplete_header:
    ABORT0 ("Error in bam_iterate_SQ_lines: incomplete BAM header");
}   

void sam_header_get_contigs (ConstBufferP *contigs_dict, ConstBufferP *contigs)
{
    if (!header_contigs.len) return; // SQ-less SAM

    *contigs      = &header_contigs;
    *contigs_dict = &header_contigs_dict;
}

void sam_header_finalize (void)
{
    buf_free (&header_contigs);
    buf_free (&header_contigs_dict);
}

static void sam_header_add_contig (const char *chrom_name, unsigned chrom_name_len, PosType last_pos, void *callback_param)
{
    WordIndex ref_chrom=WORD_INDEX_NONE;

    if (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE) 
        ref_chrom = ref_contigs_ref_chrom_from_header_chrom (chrom_name, chrom_name_len, last_pos, header_contigs.len);

    // add to header_contigs
    NEXTENT (RefContig, header_contigs) = (RefContig){ 
        .max_pos     = last_pos, 
        .char_index  = header_contigs_dict.len, 
        .snip_len    = chrom_name_len,
        .chrom_index = ref_chrom
    };

    // add to header_contigs_dict
    buf_add (&header_contigs_dict, chrom_name, chrom_name_len);
    NEXTENT (char, header_contigs_dict) = 0; // nul-termiante
}

typedef struct { uint32_t num_contigs, dict_len; } NumRangesCbParam;

static void sam_header_get_num_ranges_cb (const char *chrom_name, unsigned chrom_name_len, PosType last_pos, void *cb_param)
{
    // note: for simplicity, we consider contigs to be in the range [0, last_pos] even though POS=0 doesn't exist - last_pos+1 loci
    ((NumRangesCbParam*)cb_param)->num_contigs++;
    ((NumRangesCbParam*)cb_param)->dict_len += chrom_name_len + 1;
}

// constructs header_contigs, and in ZIP, also initialzes refererence ranges and random_access
bool sam_header_inspect (BufferP txt_header)
{    
    if (command == ZIP) {
        // if there is no external reference provided, then we create our internal one, and store it
        // (if an external reference IS provided, the user can decide whether to store it or not, with --reference vs --REFERENCE)
        if (flag.reference == REF_NONE) flag.reference = REF_INTERNAL;

        // in case of internal reference, we need to initialize. in case of --reference, it was initialized by ref_load_external_reference()
        if (!flag.reference || flag.reference == REF_INTERNAL) ref_initialize_ranges (RT_DENOVO); // it will be REF_INTERNAL if this is the 2nd+ non-conatenated file

        // evb buffers must be alloced by I/O threads, since other threads cannot modify evb's buf_list
        random_access_alloc_ra_buf (evb, 0);
    }

    if (!IS_BAM) *AFTERENT (char, *txt_header) = 0; // nul-terminate as required by sam_iterate_SQ_lines
    
    // initialize ranges array - one range per contig
    NumRangesCbParam ranges_data = {}; 
    (IS_BAM ? bam_iterate_SQ_lines : sam_iterate_SQ_lines)(txt_header->data, sam_header_get_num_ranges_cb, &ranges_data);

    // note: In BAM, a contig-less header is considered an unaligned BAM in which all RNAME=*
    // In SAM, a contig-less header just means we don't know the contigs and they are defined in RNAME 
    has_header_contigs = IS_BAM || (ranges_data.num_contigs > 0);
    
    buf_alloc (evb, &header_contigs, ranges_data.num_contigs * sizeof (RefContig), 1, "header_contigs"); 
    buf_alloc (evb, &header_contigs_dict, ranges_data.dict_len, 1, "header_contigs_dict"); 
    
    (IS_BAM ? bam_iterate_SQ_lines : sam_iterate_SQ_lines) (txt_header->data, sam_header_add_contig, NULL);

    // If we have a header with SQ lines in SAM, and always in BAM, all RNAME values must be defined in the header
    flag.const_chroms = (command == ZIP) && (IS_BAM || header_contigs.len); 

    return true;
}

// returns header length if header read is complete + sets lines.len, -1 not complete yet 
// note: usually a BAM header fits into a single 512KB READ BUFFER, so this function is called only twice (without and then with data).
// callback from DataTypeProperties.is_header_done
int32_t bam_is_header_done (void)
{
    uint32_t next=0;

    HDRSKIP(4); // magic
    ASSERT (!memcmp (evb->txt_data.data, BAM_MAGIC, 4), // magic
            "Error in bam_read_txt_header: %s doesn't have a BAM magic - it doesn't seem to be a BAM file", txt_name);

    // sam header text
    uint32_t l_text = HDR32;
    const char *text = ENT (const char, evb->txt_data, next);
    HDRSKIP(l_text);

    uint32_t header_n_ref = HDR32;

    for (uint32_t ref_i=0; ref_i < header_n_ref; ref_i++) {
        uint32_t l_name = HDR32;
        HDRSKIP (l_name + 4); // skip name and l_ref
    }

    // we have the entire header - count the text lines in the SAM header
    evb->lines.len = 0;
    for (unsigned i=0; i < l_text; i++)
        if (text[i] == '\n') evb->lines.len++;

    return next; // return BAM header length

incomplete_header:
    return -1;
}   

