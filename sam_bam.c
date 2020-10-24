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

static char *bam_read_exact (uint32_t num_bytes)
{
    char *start = AFTERENT (char, evb->txt_data);
    
    buf_alloc (evb, &evb->txt_data, evb->txt_data.len + num_bytes, 2, "txt_data", 0); // increase size if needed

    int32_t bytes_read = txtfile_read_block (start, num_bytes);
    
    ASSERT (bytes_read == (int32_t)num_bytes, "Error in bam_read_exact: needed to read %u bytes, but read only %d",
            num_bytes, bytes_read);

    evb->txt_data.len += bytes_read;
    return start;
}

// ZIP: returns the hash of the header
Md5Hash bam_read_txt_header (bool is_first_txt, bool header_required_unused, char first_char_unused)
{
    START_TIMER;

    // loading an integer from an unaligned buffer
    #define READ_UINT32(u) char *u##_ = bam_read_exact(4); \
                           uint32_t u = (((uint8_t*)u##_)[0] | (((uint8_t*)u##_)[1] << 8) | (((uint8_t*)u##_)[2] << 16) | (((uint8_t*)u##_)[2] << 24));

    buf_alloc (evb, &evb->txt_data, 20000, 1, "txt_data", 0); // this should be more than enough

    // interpret the BAM header by the spec on page 15 here: https://samtools.github.io/hts-specs/SAMv1.pdf  
    ASSERT (!memcmp (bam_read_exact (4), "BAM\1", 4), // magic
            "Error in bam_read_txt_header: %s doesn't have a BAM magic - it doesn't seem to be a BAM file", txt_name);

    // read sam header (text) and count the lines
    READ_UINT32 (l_text); 
    const char *text = bam_read_exact (l_text); 

    for (int i=0; i < l_text; i++) 
        if (text[i] == '\n') evb->lines.len++;   

    READ_UINT32 (n_ref);

    for (uint32_t i=0; i < n_ref; i++) {
        READ_UINT32 (l_name);
        const char *name = bam_read_exact (l_name);
        READ_UINT32 (l_ref);

        // check that the contigs specified in the header are consistent with the reference given in --reference/--REFERENCE
        ref_contigs_verify_identical_chrom (name, l_name-1, l_ref);

        // TODO: we need to build a translation table between ref contig index and bam contig index or something like that
    }
    
    txt_file->txt_data_so_far_single = evb->txt_data.len;

    // md5 header - with logic related to is_first
    txtfile_update_md5 (evb->txt_data.data, evb->txt_data.len, !is_first_txt);
    Md5Hash header_md5 = md5_snapshot (&z_file->md5_ctx_single);

    COPY_TIMER_VB (evb, txtfile_read_header);

    return header_md5;
}
