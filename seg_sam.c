// ------------------------------------------------------------------
//   seg_sam.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "move_to_front.h"
#include "header.h"
#include "file.h"

#define DATA_LINE(vb,i) (&((ZipDataLineSAM *)((vb)->data_lines))[(i)])

static inline uint32_t seg_sam_one_field (VBlockSAM *vb, const char *str, unsigned len, unsigned vb_line_i, SamFields f, 
                                          bool *is_new)  // optional out
{
    return seg_one_field ((VBlockP)vb, str, len, vb_line_i, f, SEC_SAM_QNAME_B250 + f*2, is_new);
}                                          
  

// We break down the QNAME into subfields separated by / and/or : - these are vendor-defined strings.
// Up to MAX_SUBFIELDS subfields are permitted.
// each subfield is stored in its own directory called QNAME_n where n is the subfield number starting from 0. 
// The separators are made into a string we call "template" that is stored in the QNAME directory - we anticipate that 
// usually all lines have the same format, but we allow lines to have different formats.
// If a subfield is an integer, we delta-encode it.
// QNAME formats:
// Illumina: <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> for example "A00488:61:HMLGNDSXX:4:1101:15374:1031" see here: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
// PacBio BAM: {movieName}/{holeNumber}/{qStart}_{qEnd} see here: https://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
static void seg_sam_qname_field (VBlockSAM *vb, const char *qname, unsigned qname_len, unsigned vb_line_i)
{
    static DictIdType sf_dict_id = { .id = { 'Q' | 0xc0 ,'N','A','M','E','0','0', 0} }; // 0xc0 because it is a QNAME subfield (see dict_id.h)

    const char *snip = qname;
    unsigned snip_len = 0;
    unsigned sf_i = 0;
    char template[MAX_SUBFIELDS];
    MtfNode *node;

    // add each subfield to its dictionary QNAME00 through QNAME62
    for (unsigned i=0; i <= qname_len; i++) { // one more than qname_len - to finalize the last subfield

        if (i==qname_len || qname[i]==':' || qname[i]=='/') { // a subfield ended - separator between subfields
            
            // process the subfield that just ended
            MtfContext *sf_ctx;

            if (vb->qname_mapper.num_subfields == sf_i) { // new subfield in this VB (sf_ctx might exist from previous VBs)
                sf_dict_id.id[5] = sf_i / 10 + '0';
                sf_dict_id.id[6] = sf_i % 10 + '0';

                sf_ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, &vb->num_dict_ids, NULL, sf_dict_id, SEC_SAM_QNAME_SF_DICT);
                vb->qname_mapper.did_i[sf_i] = sf_ctx->did_i;
                vb->qname_mapper.num_subfields++;
            }
            else 
                sf_ctx = MAPPER_CTX (&vb->qname_mapper, sf_i);

            ASSERT0 (sf_ctx, "Error in seg_sam_qname_field: sf_ctx is NULL");

            // allocate memory if needed
            buf_alloc (vb, &sf_ctx->mtf_i, MIN (vb->num_lines, sf_ctx->mtf_i.len + 1) * sizeof (uint32_t),
                       1.5, "mtf_ctx->mtf_i", SEC_SAM_QNAME_SF_DICT);

            NEXTENT (uint32_t, sf_ctx->mtf_i) = mtf_evaluate_snip_seg ((VBlockP)vb, sf_ctx, snip, snip_len, &node, NULL);

            // finalize this subfield and get ready for reading the next one
            if (i < qname_len) {
                template[sf_i] = qname[i];
                snip = &qname[i+1];
                snip_len = 0;
            }
            sf_i++;
        }
        else snip_len++;
    }

    // if template is empty, make it "*"
    if (sf_i==1) template[0] = '*';

    // add template to the QNAME dictionary (note: template may be of length zero if qname has no / or :)
    MtfContext *qname_ctx = &vb->mtf_ctx[SAM_QNAME];
    NEXTENT (uint32_t, qname_ctx->mtf_i) = mtf_evaluate_snip_seg ((VBlockP)vb, qname_ctx, template, MAX (1, sf_i-1), &node, NULL);

    // byte counts for --show-sections 
    vb->txt_section_bytes[SEC_SAM_QNAME_B250]    += sf_i; // sf_i has 1 for each separator including the terminating \t
    vb->txt_section_bytes[SEC_SAM_QNAME_SF_B250] += qname_len - (sf_i-1); // the entire field except for the / and : separators
}

// process an optional subfield, that looks something like MX:Z:abcdefg. We use "MX" for the field name, and
// the data is abcdefg. The full name "MX:Z:" is stored as part of the OPTIONAL dictionary entry
static void seg_sam_optional_field (VBlockSAM *vb, const char *field, unsigned field_len, char separator, unsigned vb_line_i)
{
    ASSERT (field_len >= 6 && field[2] == ':' && field[4] == ':', "Error in %s: invalid optional field format: %.*s",
            txt_name, field_len, field);

    DictIdType dict_id = dict_id_sam_optnl_sf (dict_id_make (field, 2));

    MtfNode *node;
    MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, &vb->num_dict_ids, NULL, dict_id, SEC_SAM_OPTNL_SF_DICT);

    // allocate memory if needed
    buf_alloc (vb, &ctx->mtf_i, MIN (vb->num_lines, ctx->mtf_i.len + 1) * sizeof (uint32_t),
                1.5, "mtf_ctx->mtf_i", SEC_SAM_OPTNL_SF_DICT);

    NEXTENT (uint32_t, ctx->mtf_i) = mtf_evaluate_snip_seg ((VBlockP)vb, ctx, &field[5], field_len-5, &node, NULL);

    // account for the subfield value and \t, but not \n which is already accounted for in SEC_SAM_OPTIONAL_B250
    vb->txt_section_bytes[SEC_SAM_OPTNL_SF_B250] += (field_len-5) + (separator == '\t'); 
}

// calculate the expected length of SEQ and QUAL from the CIGAR string
// A CIGAR looks something like: "109S19M23S", See: https://samtools.github.io/hts-specs/SAMv1.pdf 
uint32_t seg_sam_seq_len_from_cigar (const char *cigar, unsigned cigar_len)
{
    // ZIP case: if the CIGAR is "*", we later get the length from SEQ and store it as eg "151*". 
    // In PIZ it will be eg "151*" or "1*" if both SEQ and QUAL are "*"
    if (cigar_len == 1 && cigar[0] == '*') return 0;

    unsigned seq_len=0, n=0;

    for (unsigned i=0; i < cigar_len; i++) {
        char c = cigar[i];
        if (c >= '0' && c <= '9') 
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

static void seg_sam_seq_qual_fields (VBlockSAM *vb, ZipDataLineSAM *dl, unsigned seq_len, unsigned qual_len, unsigned vb_line_i)
{
    const char *seq  = &vb->txt_data.data[dl->seq_data_start];
    const char *qual = &vb->txt_data.data[dl->qual_data_start];

    bool seq_is_available  = seq_len  != 1 || *seq  != '*';
    bool qual_is_available = qual_len != 1 || *qual != '*';

    // handle the case where CIGAR is "*" and we get the dl->seq_len directly from SEQ or QUAL
    if (!dl->seq_len) { // CIGAR is not available
        ASSERT (!seq_is_available || !qual_is_available || seq_len==qual_len, 
                "Bad line in %s: SEQ length is %u, QUAL length is %u, unexpectedly differ. SEQ=%.*s QUAL=%.*s", 
                txt_name, seq_len, qual_len, seq_len, seq, qual_len, qual);    

        dl->seq_len = MAX (seq_len, qual_len); // one or both might be not available and hence =1
    
        char new_cigar[20];
        sprintf (new_cigar, "%u*", dl->seq_len); // eg "151*". if both SEQ and QUAL are unavailable it will be "1*"
        unsigned new_cigar_len = strlen (new_cigar);

        seg_sam_one_field (vb, new_cigar, new_cigar_len, vb_line_i, SAM_CIGAR, NULL); 
        vb->txt_section_bytes[SEC_SAM_CIGAR_B250] -= (new_cigar_len-1); // adjust - actual SAM length was only one (seg_sam_one_field added new_cigar_len)
    } 
    else { // CIGAR is available - just check the seq and qual lengths
        ASSERT (!seq_is_available || seq_len == dl->seq_len, 
                "Bad line in %s: according to CIGAR, expecting SEQ length to be %u but it is %u. SEQ=%.*s", 
                txt_name, dl->seq_len, seq_len, seq_len, seq);

        ASSERT (!qual_is_available || qual_len == dl->seq_len, 
                "Bad line in %s: according to CIGAR, expecting QUAL length to be %u but it is %u. QUAL=%.*s", 
                txt_name, dl->seq_len, qual_len, qual_len, qual);    
    }

    // byte counts for --show-sections 
    vb->txt_section_bytes[SEC_SAM_SEQ_DATA]  += seq_len  + 1; // +1 for terminating \t
    vb->txt_section_bytes[SEC_SAM_QUAL_DATA] += qual_len + 1; 
}

/* split the data line into sections - 
   1. variant data - each of 9 fields (CHROM to FORMAT) is a section
   2. genotype data (except the GT subfield) - is a section
   3. haplotype data (the GT subfield) - a string of eg 1 or 0 (no whitespace) - ordered by the permutation order
   4. phase data - only if MIXED phase - a string of | and / - one per sample
*/
const char *seg_sam_data_line (VBlock *vb_,   
                               const char *field_start_line,     // index in vb->txt_data where this line starts
                               uint32_t vb_line_i) // line within this vb (starting from 0)
{
    VBlockSAM *vb = (VBlockSAM *)vb_;
    ZipDataLineSAM *dl = DATA_LINE (vb, vb_line_i);

    const char *next_field, *field_start;
    unsigned field_len=0, seq_len=0, qual_len=0;
    char separator;
    bool has_13 = false; // does this line end in Windows-style \r\n rather than Unix-style \n

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    // QNAME
    field_start = field_start_line;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "QNAME");
    seg_sam_qname_field (vb, field_start, field_len, vb_line_i);

    // FLAG
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "FLAG");
    seg_sam_one_field (vb, field_start, field_len, vb_line_i, SAM_FLAG, NULL);

    // RNAME
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "RNAME");
    seg_sam_one_field (vb, field_start, field_len, vb_line_i, SAM_RNAME, NULL);

    // POS - delta encode
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "POS");
    vb->last_pos = seg_pos_field (vb_, vb->last_pos, SAM_POS, SEC_SAM_POS_B250, field_start, field_len, vb_line_i);

    // MAPQ
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "MAPQ");
    seg_sam_one_field (vb, field_start, field_len, vb_line_i, SAM_MAPQ, NULL);

    // CIGAR - if CIGAR is "*" and we wait to get the length from SEQ or QUAL
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "CIGAR");
    dl->seq_len = seg_sam_seq_len_from_cigar (field_start, field_len);
    if (dl->seq_len) seg_sam_one_field (vb, field_start, field_len, vb_line_i, SAM_CIGAR, NULL); // not "*" - all good!

    // RNEXT
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "RNEXT");
    seg_sam_one_field (vb, field_start, field_len, vb_line_i, SAM_RNEXT, NULL);

    // PNEXT - encode as delta from POS (except if its unavailable - encode as "*")
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "PNEXT");
    if (field_len == 1 && *field_start == '0') // PNEXT is "unavailable" per SAM specification
        seg_sam_one_field (vb, "*", 1, vb_line_i, SAM_PNEXT, NULL); 
    else
        seg_pos_field (vb_, vb->last_pos, SAM_PNEXT, SEC_SAM_PNEXT_B250, field_start, field_len, vb_line_i);

    // TLEN
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "TLEN");
    seg_sam_one_field (vb, field_start, field_len, vb_line_i, SAM_TLEN, NULL);

    // SEQ & QUAL
    dl->seq_data_start = next_field - vb->txt_data.data;
    next_field = seg_get_next_item (next_field, &len, false, true, false, vb_line_i, &seq_len, &separator, &has_13, "SEQ");

    dl->qual_data_start = next_field - vb->txt_data.data;
    next_field = seg_get_next_item (next_field, &len, true, true, false, vb_line_i, &qual_len, &separator, &has_13, "QUAL");

    seg_sam_seq_qual_fields (vb, dl, seq_len, qual_len, vb_line_i); // also updates cigar if it is "*"

    // OPTIONAL fields - up to MAX_SUBFIELDS of them
    char oname[MAX_SUBFIELDS * 5 + 1]; // each name is 5 characters per SAM specification, eg "MC:Z:" ; +1 for # in case of \r
    unsigned oname_len;

    for (oname_len=0; oname_len < MAX_SUBFIELDS*5 && separator != '\n'; oname_len += 5) {
        field_start = next_field;
        next_field = seg_get_next_item (field_start, &len, true, true, false, vb_line_i, &field_len, &separator, &has_13, "OPTIONAL-subfield");
        seg_sam_optional_field (vb, field_start, field_len, separator, vb_line_i);
        memcpy (&oname[oname_len], field_start, 5);
    }
    ASSERT (separator=='\n', "Error in %s: too many optional fields, limit is %u", txt_name, MAX_SUBFIELDS);

    // treat a Windows style \r (as part of a \r\n line termination) as part of the oname - add a # to represent it
    if (has_13) oname[oname_len++] = '#';

    // if oname is empty, we make it "*" (and adjust size accounting)
    if (!oname_len) oname[oname_len++] = '*';

    seg_sam_one_field (vb, oname, oname_len, vb_line_i, SAM_OPTIONAL, NULL);

    // if we have no real OPTIONAL, seg_sam_one_field adds +1 for the character (correct for '#' as it stands for \r, but incorrect for '*') and +1 for the wrongly presumsed \t
    if (oname_len == 1) 
        vb->txt_section_bytes[SEC_SAM_OPTIONAL_B250] -= (1 + (oname[0]=='*'));  

    return next_field;
}
