// ------------------------------------------------------------------
//   sam_shared.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "sam_private.h"
#include "strings.h"
#include "file.h"

unsigned sam_vb_size (void) { return sizeof (VBlockSAM); }
unsigned sam_vb_zip_dl_size (void) { return sizeof (ZipDataLineSAM); }

void sam_vb_release_vb (VBlockSAM *vb)
{
    memset (&vb->qname_mapper, 0, sizeof (vb->qname_mapper));
    buf_free (&vb->optional_mapper_buf);
}

void sam_vb_destroy_vb (VBlockSAM *vb)
{
    buf_destroy (&vb->optional_mapper_buf);
}

// calculate the expected length of SEQ and QUAL from the CIGAR string
// A CIGAR looks something like: "109S19M23S", See: https://samtools.github.io/hts-specs/SAMv1.pdf 
uint32_t sam_seq_len_from_cigar (const char *cigar, unsigned cigar_len)
{
    // ZIP case: if the CIGAR is "*", we later get the length from SEQ and store it as eg "151*". 
    // In PIZ it will be eg "151*" or "1*" if both SEQ and QUAL are "*"
    if (cigar_len == 1 && cigar[0] == '*') return 0;

    unsigned seq_len=0, n=0;

    for (unsigned i=0; i < cigar_len; i++) {
        char c = cigar[i];
        if (IS_DIGIT (c)) 
            n = n*10 + (c - '0');
        
        else if (c=='M' || c=='I' || c=='S' || c=='=' || c=='X' || c=='*') { // these "consume" sequences
            ASSERT (n, "Error: Invalid CIGAR in %s: operation %c not preceded by a number. CIGAR=%.*s", 
                    txt_name, c, cigar_len, cigar);
            seq_len += n;
            n = 0;
        }

        else if (c=='D' || c=='N' || c=='H' || c=='P') { // these don't consume sequence
            ASSERT (n, "Error: Invalid CIGAR in %s: operation %c not preceded by a number. CIGAR=%.*s", 
                    txt_name, c, cigar_len, cigar);
            n = 0;
        }

        else ABORT ("Error: Invalid CIGAR in %s: invalid operation %c. CIGAR=%.*s", 
                    txt_name, c, cigar_len, cigar);
    }                          

    ASSERT (!n, "Error: Invalid CIGAR in %s: expecting it to end with an operation character. CIGAR=%.*s", 
            txt_name, cigar_len, cigar);

    ASSERT (seq_len, "Error: Invalid CIGAR in %s: CIGAR implies 0-length SEQ. CIGAR=%.*s", 
            txt_name, cigar_len, cigar);

    return seq_len;
}
