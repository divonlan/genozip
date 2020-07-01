// ------------------------------------------------------------------
//   sam_shared.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "sam_private.h"
#include "strings.h"
#include "file.h"

// table to determine whether a character is a valid cigar_op as defined in https://samtools.github.io/hts-specs/SAMv1.pdf
// bits - 1=digit 2=consumes_query 4=consumes_ref 8=error
const uint8_t cigar_lookup[256] = {
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, // ASCII 0-15
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, // ASCII 16-31
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, // ASCII 32-47  '*' consumes both
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 6, 8, 8, // ASCII 48-63  '0'-'9' are digits, '=' consumes both
    8, 8, 8, 8, 4, 8, 8, 8, 0, 2, 8, 8, 8, 6, 4, 8, // ASCII 64-79  'D', 'N' - referece only, 'I' - query only, 'H' - none, 'M' - both
    0, 8, 8, 2, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, // ASCII 80-95  'P' - none, 'S' - query only, 'X' - both
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, // ASCII 96-111  
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, // ASCII 112-127 
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, // ASCII 128-255
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 };

unsigned sam_vb_size (void) { return sizeof (VBlockSAM); }
unsigned sam_vb_zip_dl_size (void) { return sizeof (ZipDataLineSAM); }

void sam_vb_release_vb (VBlockSAM *vb)
{
    memset (&vb->qname_mapper, 0, sizeof (vb->qname_mapper));
    buf_free (&vb->optional_mapper_buf);
    vb->last_cigar = NULL;
    vb->ref_consumed = 0;
}

void sam_vb_destroy_vb (VBlockSAM *vb)
{
    buf_destroy (&vb->optional_mapper_buf);
}

// calculate the expected length of SEQ and QUAL from the CIGAR string
// A CIGAR looks something like: "109S19M23S", See: https://samtools.github.io/hts-specs/SAMv1.pdf 
void sam_analyze_cigar (const char *cigar, unsigned cigar_len, 
                        unsigned *seq_consumed, unsigned *ref_consumed, unsigned *seq_and_ref) // optional outs
{
    if (seq_consumed) *seq_consumed = 0;
    if (ref_consumed) *ref_consumed = 0;
    if (seq_and_ref)  *seq_and_ref  = 0;

    // ZIP case: if the CIGAR is "*", we later get the length from SEQ and store it as eg "151*". 
    // In PIZ it will be eg "151*" or "1*" if both SEQ and QUAL are "*"
    if (cigar_len == 1 && cigar[0] == '*') return;

    // PIZ case: CIGAR string starts with '-' (indicating missing SEQ) - just skip the '-' for now
    if (*cigar == '-') {
        cigar++;
        cigar_len--;
    }
    
    unsigned n=0;
    for (unsigned i=0; i < cigar_len; i++) {

        char c = cigar[i];
        char lookup = cigar_lookup[(uint8_t)c];

        ASSERT (lookup != CIGAR_INVALID, "Error: Invalid CIGAR in %s: invalid operation %c. CIGAR=%.*s", 
                txt_name, cigar[i], cigar_len, cigar);

        if (lookup == CIGAR_DIGIT) 
            n = n*10 + (c - '0');
        
        else {
            ASSERT (n, "Error: Invalid CIGAR in %s: operation %c not preceded by a number. CIGAR=%.*s", 
                    txt_name, c, cigar_len, cigar);
            
            if ((lookup & CIGAR_CONSUMES_QUERY)     && seq_consumed) *seq_consumed += n;
            if ((lookup & CIGAR_CONSUMES_REFERENCE) && ref_consumed) *ref_consumed += n;
            if ((lookup & CIGAR_CONSUMES_QUERY) && (lookup & CIGAR_CONSUMES_REFERENCE) && seq_and_ref) *seq_and_ref += n;

            n = 0;
        }
    }                          

    ASSERT (!n, "Error: Invalid CIGAR in %s: expecting it to end with an operation character. CIGAR=%.*s", 
            txt_name, cigar_len, cigar);

    ASSERT (!seq_consumed || *seq_consumed, "Error: Invalid CIGAR in %s: CIGAR implies 0-length SEQ. CIGAR=%.*s", txt_name, cigar_len, cigar);
}
