// ------------------------------------------------------------------
//   seg_sam.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
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

                sf_ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, &vb->num_dict_ids, &vb->qname_mapper.num_subfields,
                                              sf_dict_id, SEC_SAM_QNAME_SF_DICT);
                vb->qname_mapper.did_i[sf_i] = sf_ctx->did_i;
            }
            else sf_ctx = MAPPER_CTX (&vb->qname_mapper, sf_i);

            ASSERT0 (sf_ctx, "Error in seg_sam_qname_field: sf_ctx is NULL");

            // allocate memory if needed
            buf_alloc (vb, &sf_ctx->mtf_i, MIN (vb->num_lines, sf_ctx->mtf_i.len + 1) * sizeof (uint32_t),
                       1.5, "mtf_ctx->mtf_i", SEC_SAM_QNAME_SF_DICT);

            *NEXTENT (uint32_t, &sf_ctx->mtf_i) = mtf_evaluate_snip_seg ((VBlockP)vb, sf_ctx, snip, snip_len, &node, NULL);

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

    // add template to the QNAME dictionary (note: template may be of length zero if qname has no / or :)
    MtfContext *qname_ctx = &vb->mtf_ctx[SAM_QNAME];
    *NEXTENT (uint32_t, &qname_ctx->mtf_i) = mtf_evaluate_snip_seg ((VBlockP)vb, qname_ctx, template, sf_i-1, &node, NULL);
}

// process an optional subfield, that looks something like MX:Z:abcdefg. We use "MX" for the field name, and
// the data is abcdefg. The full name "MX:Z:" is stored as part of the OPTIONAL dictionary entry
static void seg_sam_optional_field (VBlockSAM *vb, const char *field, unsigned field_len, unsigned vb_line_i)
{
    ASSERT (field_len >= 6 && field[2] == ':' && field[4] == ':', "Error in %s: invalid optional field format: %.*s",
            txt_name, field_len, field);

    DictIdType dict_id = dict_id_sam_optnl_sf (dict_id_make (field, 2));

    MtfNode *node;
    MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, &vb->num_dict_ids, NULL, dict_id, SEC_SAM_OPTNL_SF_DICT);

    // allocate memory if needed
    buf_alloc (vb, &ctx->mtf_i, MIN (vb->num_lines, ctx->mtf_i.len + 1) * sizeof (uint32_t),
                1.5, "mtf_ctx->mtf_i", SEC_SAM_OPTNL_SF_DICT);

    *NEXTENT (uint32_t, &ctx->mtf_i) = mtf_evaluate_snip_seg ((VBlockP)vb, ctx, &field[5], field_len-5, &node, NULL);
}

static void seg_sam_rname_rnext_fields (VBlockSAM *vb, 
                                        const char *rname, unsigned rname_len, 
                                        const char *rnext, unsigned rnext_len,
                                        unsigned vb_line_i)
{
#   define MAX_RNAME_NX 1000
    char rname_nx[MAX_RNAME_NX];
    unsigned rname_nx_len = rname_len + rnext_len + 1; // +1 for the \t
    
    ASSERT (rname_nx_len <= MAX_RNAME_NX, 
            "Error in %s: the length of RNAME and RNEXT fields, combined, cannot exceed %u characters. RNAME=%.*s RNEXT=%*.s",
            txt_name, MAX_RNAME_NX-1, rname_len, rname, rnext_len, rnext);

    memcpy (rname_nx, rname, rname_len);
    rname_nx[rname_len] = '\t';
    memcpy (&rname_nx[rname_len+1], rnext, rnext_len);

    seg_sam_one_field (vb, rname_nx, rname_nx_len, vb_line_i, SAM_RNAMENX, NULL);
}   

// calculate the expected length of SEQ and QUAL from the CIGAR string
uint32_t seg_sam_seq_len_from_cigar (const char *cigar, unsigned cigar_len)
{
    // if the CIGAR is '*', we get the length from SEQ and store it as eg "151*"
    if (cigar_len == 1 && cigar[0] == '*') return 0;

    unsigned seq_len=0, n=0;

    for (unsigned i=0; i < cigar_len; i++) {
        char c = cigar[i];
        if (c >= '0' && c <= '9') 
            n = n*10 + (c - '0');
        
        else if (c=='M' || c=='I' || c=='S' || c=='=' || c=='X') {
            ASSERT (n, "Error: Invalid CIGAR in %s: operation %c not preceded by a number. CIGAR=%.*s", 
                    txt_name, c, cigar_len, cigar);
            seq_len += n;
            n = 0;
        }

        else if (c=='D' || c=='N' || c=='H' || c=='P') {
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
    unsigned field_len=0;
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

    // RNAME - wait to store together with RNEXT, because they are highly correlated
    const char *rname_field_start = next_field;
    unsigned rname_field_len;
    next_field = seg_get_next_item (rname_field_start, &len, false, true, false, vb_line_i, &rname_field_len, &separator, &has_13, "RNAME");

    // POS - delta encode
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "POS");
    vb->last_pos = seg_pos_field (vb_, vb->last_pos, SAM_POS, SEC_SAM_POS_B250, field_start, field_len, vb_line_i);

    // MAPQ
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "MAPQ");
    seg_sam_one_field (vb, field_start, field_len, vb_line_i, SAM_MAPQ, NULL);

    // CIGAR - if CIGAR is "*" and we wait to get the length from SEQ
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "CIGAR");
    dl->seq_len = seg_sam_seq_len_from_cigar (field_start, field_len);
    if (dl->seq_len) seg_sam_one_field (vb, field_start, field_len, vb_line_i, SAM_CIGAR, NULL); // not "*" - all good!

    // RNEXT - store together with RNAME as they are highly correlated
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "RNEXT");
    seg_sam_rname_rnext_fields (vb, rname_field_start, rname_field_len, field_start, field_len, vb_line_i);

    // PNEXT - encode as delta from POS
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "PNEXT");
    seg_pos_field (vb_, vb->last_pos, SAM_PNEXT, SEC_SAM_PNEXT_B250, field_start, field_len, vb_line_i);

    // TLEN
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "TLEN");
    seg_sam_one_field (vb, field_start, field_len, vb_line_i, SAM_TLEN, NULL);

    // SEQ
    dl->seq_data_start = next_field - vb->txt_data.data;
    next_field = seg_get_next_item (next_field, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "SEQ");
    
    // handle the case where CIGAR is "*" and we get the seq_len directly from SEQ
    char new_cigar[20];
    if (!dl->seq_len) { // CIGAR is "*"
        sprintf (new_cigar, "%u*", field_len); // eg "151*"
        seg_sam_one_field (vb, new_cigar, strlen (new_cigar), vb_line_i, SAM_CIGAR, NULL); 
        dl->seq_len = field_len;
    }

    ASSERT (field_len == dl->seq_len, "Bad line in %s: according to CIGAR, expecting SEQ length to be %u but it is %u",
            txt_name, dl->seq_len, field_len);

    // QUAL
    dl->qual_data_start = next_field - vb->txt_data.data;
    next_field = seg_get_next_item (next_field, &len, true, true, false, vb_line_i, &field_len, &separator, &has_13, "QUAL");
    ASSERT (field_len == dl->seq_len, "Bad line in %s: according to CIGAR, expecting QUAL length to be %u but it is %u",
            txt_name, dl->seq_len, field_len);

    // OPTIONAL fields - up to MAX_SUBFIELDS of them
    char oname[MAX_SUBFIELDS * 5 + 1]; // each name is 5 characters per SAM specification, eg "MC:Z:" ; +1 for # in case of \r
    unsigned oname_len;

    for (oname_len=0; oname_len < MAX_SUBFIELDS*5 && separator != '\n'; oname_len += 5) {
        field_start = next_field;
        next_field = seg_get_next_item (field_start, &len, true, true, false, vb_line_i, &field_len, &separator, &has_13, "SEQ");
        seg_sam_optional_field (vb, field_start, field_len, vb_line_i);
        memcpy (&oname[oname_len], field_start, 5);
    }
    ASSERT (separator=='\n', "Error in %s: too many optional fields, limit is %u", txt_name, MAX_SUBFIELDS);

    // treat a Windows style \r (as part of a \r\n line termination) as part of the oname - add a # to represent it
    if (has_13) {
        oname[oname_len++] = '#';
        vb->txt_section_bytes[SEC_STATS_HT_SEPERATOR]--; // the \r in case of Windows \r\n line ending (WHY IS THIS?)
    }

    // we record the oname, even if its an empty string 
    seg_sam_one_field (vb, oname, oname_len, vb_line_i, SAM_OPTIONAL, NULL);

    return next_field;
}
