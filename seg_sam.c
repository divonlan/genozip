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
#include "random_access.h"

#define DATA_LINE(vb,i) (&((ZipDataLineSAM *)((vb)->data_lines))[(i)])

void seg_sam_initialize (VBlockSAM *vb)
{
    buf_alloc (vb, &vb->random_pos_data, 1000000, 1, "random_pos_data", vb->vblock_i); // initial ~1MB allocation

    vb->last_pos = 0;
    vb->last_rname_node_index = (uint32_t)-1;

    vb->mc_did_i = DID_I_NONE;
}

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

                sf_ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, vb->dict_id_to_did_i_map, &vb->num_dict_ids, NULL, sf_dict_id, SEC_SAM_QNAME_SF_DICT);
                vb->qname_mapper.did_i[sf_i] = sf_ctx->did_i;
                vb->qname_mapper.num_subfields++;
            }
            else 
                sf_ctx = MAPPER_CTX (&vb->qname_mapper, sf_i);

            ASSERT0 (sf_ctx, "Error in seg_sam_qname_field: sf_ctx is NULL");

            // allocate memory if needed
            buf_alloc (vb, &sf_ctx->mtf_i, MAX (vb->num_lines, sf_ctx->mtf_i.len + 1) * sizeof (uint32_t),
                       CTX_GROWTH, "mtf_ctx->mtf_i", sf_ctx->did_i);

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

// returns length of string ending with separator, or -1 if separator was not found
static inline int seg_sam_get_next_subitem (const char *str, int str_len, char separator)
{
    for (int i=0; i < str_len; i++) {
        if (str[i] == separator) return i;
        if (str[i] == ',' || str[i] == ';') return -1; // wrong separator encountered
    }
    return -1;
}

#define DO_SSF(ssf,sep) \
        ssf = &field[i]; \
        ssf##_len = seg_sam_get_next_subitem (&field[i], field_len-i, sep); \
        if (ssf##_len == -1) goto error; /* bad format */ \
        i += ssf##_len + 1; /* skip snip and separator */        

#define DEC_SSF(ssf) const char *ssf; \
                     int ssf##_len;

static unsigned seg_sam_sa_or_oa_field (VBlockSAM *vb, const char *field, unsigned field_len, unsigned vb_line_i)
{
    // OA and SA format is: (rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+ . in OA - NM is optional (but its , is not)
    // Example SA:Z:chr13,52863337,-,56S25M70S,0,0;chr6,145915118,+,97S24M30S,0,0;chr18,64524943,-,13S22M116S,0,0;chr7,56198174,-,20M131S,0,0;chr7,87594501,+,34S20M97S,0,0;chr4,12193416,+,58S19M74S,0,0;
    // See: https://samtools.github.io/hts-specs/SAMtags.pdf
    DEC_SSF(rname); DEC_SSF(pos); DEC_SSF(strand); DEC_SSF(cigar); DEC_SSF(mapq); DEC_SSF(nm); 
    int repeats = 0;

    for (unsigned i=0; i < field_len; repeats++) {
        DO_SSF (rname,  ','); 
        DO_SSF (pos,    ','); 
        DO_SSF (strand, ','); 
        DO_SSF (cigar,  ','); 
        DO_SSF (mapq,   ','); 
        DO_SSF (nm,     ';'); 

        // sanity checks before adding to any dictionary
        if (strand_len != 1 || (strand[0] != '+' && strand[0] != '-')) goto error; // invalid format

        // pos - add to the random pos data together with all other random pos data (originating from POS, PNEXT etc).
        buf_alloc_more (vb, &vb->random_pos_data, pos_len+1, vb->num_lines, char, 2);
        buf_add (&vb->random_pos_data, pos, pos_len); 
        buf_add (&vb->random_pos_data, "\t", 1); 
        vb->txt_section_bytes[SEC_SAM_RAND_POS_DATA] += pos_len + 1; 

        // nm : we store in the same dictionary as the Optional subfield NM
        seg_one_subfield ((VBlockP)vb, nm, nm_len, vb_line_i, (DictIdType)dict_id_OPTION_NM, 
                          SEC_SAM_OPTNL_SF_B250, nm_len+1);

        // strand : we store in in the STRAND dictionary (used by SA, OA, XA)
        seg_one_subfield ((VBlockP)vb, strand, 1, vb_line_i, (DictIdType)dict_id_OPTION_STRAND, 
                          SEC_SAM_OPTNL_SF_B250, 2);

        // rname, cigar and mapq: we store in the same dictionary as the primary fields
        seg_sam_one_field (vb, rname, rname_len, vb_line_i, SAM_RNAME, NULL);
        seg_sam_one_field (vb, cigar, cigar_len, vb_line_i, SAM_CIGAR, NULL);
        seg_sam_one_field (vb, mapq,  mapq_len, vb_line_i, SAM_MAPQ, NULL);
    }

    return repeats;

error:
    // if the error occurred on on the first repeat - this file probably has a different
    // format - we return 0 and the subfield will be stored as a normal subfield
    // if it occurred on the 2nd+ subfield, after the 1st one was fine - we reject the file
    ASSERT (!repeats, "Invalid format in repeat #%u of field SA or OA, vb_i=%u vb_line_i=%u. data: %.*s",
            repeats+1, vb->vblock_i, vb_line_i, field_len, field);
    return 0;
}

static unsigned seg_sam_xa_field (VBlockSAM *vb, const char *field, unsigned field_len, unsigned vb_line_i)
{
    // XA format is: (chr,pos,CIGAR,NM;)*  pos starts with +- which is strand
    // Example XA:Z:chr9,-60942781,150M,0;chr9,-42212061,150M,0;chr9,-61218415,150M,0;chr9,+66963977,150M,1;
    // See: http://bio-bwa.sourceforge.net/bwa.shtml
    DEC_SSF(rname); DEC_SSF(pos); DEC_SSF(cigar); DEC_SSF(nm); 
    int repeats = 0;

    for (unsigned i=0; i < field_len; repeats++) {
        DO_SSF (rname,  ','); 
        DO_SSF (pos,    ','); 
        DO_SSF (cigar,  ','); 
        DO_SSF (nm,     ';'); 

        // sanity checks before adding to any dictionary
        if (pos_len < 2 || (pos[0] != '+' && pos[0] != '-')) goto error; // invalid format - expecting pos to begin with the strand

        // pos - add to the pos data together with all other pos data (POS, PNEXT etc),
        // strand - add to separate STRAND dictionary, to not adversely affect compression of POS.
        // there is no advantage to storing the strand together with pos as they are not correlated.
        buf_alloc_more (vb, &vb->random_pos_data, pos_len, 0, char, 2);
        buf_add (&vb->random_pos_data, pos+1, pos_len-1); 
        buf_add (&vb->random_pos_data, "\t", 1); 
        vb->txt_section_bytes[SEC_SAM_RAND_POS_DATA] += pos_len; // one less for the strand, one more for the , 

        // strand
        seg_one_subfield ((VBlockP)vb, pos, 1, vb_line_i, (DictIdType)dict_id_OPTION_STRAND, 
                          SEC_SAM_OPTNL_SF_B250, 1);

        // nm : we store in the same dictionary as the Optional subfield NM
        seg_one_subfield ((VBlockP)vb, nm, nm_len, vb_line_i, (DictIdType)dict_id_OPTION_NM, 
                          SEC_SAM_OPTNL_SF_B250, nm_len+1);

        // rname and cigar: we store in the same dictionary as the primary fields
        seg_sam_one_field (vb, rname, rname_len, vb_line_i, SAM_RNAME, NULL);
        seg_sam_one_field (vb, cigar, cigar_len, vb_line_i, SAM_CIGAR, NULL);
    }

    return repeats;

error:
    // if the error occurred on on the first repeat - this file probably has a different
    // format - we return 0 and the subfield will be stored as a normal subfield
    // if it occurred on the 2nd+ subfield, after the 1st one was fine - we reject the file
    ASSERT (!repeats, "Invalid format in repeat #%u of field XA (0-based), vb_i=%u vb_line_i=%u. data: %.*s",
            repeats, vb->vblock_i, vb_line_i, field_len, field);
    return 0;
}

// process an optional subfield, that looks something like MX:Z:abcdefg. We use "MX" for the field name, and
// the data is abcdefg. The full name "MX:Z:" is stored as part of the OPTIONAL dictionary entry
static void seg_sam_optional_field (VBlockSAM *vb, const char *field, unsigned field_len, 
                                    char separator, unsigned vb_line_i)
{
    ASSERT (field_len, "Error in %s: line invalidly ends with a tab. vb_i=%u vb_line_i=%u", 
            txt_name, vb->vblock_i, vb_line_i);

    ASSERT (field_len >= 6 && field[2] == ':' && field[4] == ':', "Error in %s: invalid optional field format: %.*s",
            txt_name, field_len, field);

    DictIdType dict_id = dict_id_sam_optnl_sf (dict_id_make (field, 2));

    // some fields have "repeats" - multiple instances of the same format of data separated by a ;
    unsigned repeats = 0; 
    if (dict_id.num == dict_id_OPTION_SA || dict_id.num == dict_id_OPTION_OA)
        repeats = seg_sam_sa_or_oa_field (vb, &field[5], field_len-5, vb_line_i);

    else if (dict_id.num == dict_id_OPTION_XA && field[3] == 'Z') 
        repeats = seg_sam_xa_field (vb, &field[5], field_len-5, vb_line_i);

    // special subfield - in the subfield, store just the (char)-1 + the number of repeats.
    // the data itself was already stored in the subsubfields in seg_sam_*_field
    if (repeats) {
        char repeats_str[20];
        sprintf (repeats_str, "%c%u", -1, repeats); // (char)-1 to indicate repeats (invalid char per SAM specification, so it cannot appear organically)
        seg_one_subfield ((VBlock *)vb, repeats_str, strlen (repeats_str), vb_line_i, dict_id, SEC_SAM_OPTNL_SF_B250, 1 /* \t */); 
    }

    else { // non-repeating fields

        // fields containing CIGAR format data
        if (dict_id.num == dict_id_OPTION_MC || dict_id.num == dict_id_OPTION_OC) 
            seg_sam_one_field (vb, &field[5], field_len-5, vb_line_i, SAM_CIGAR, NULL);

        // mc:i: (output of bamsormadup? - mc in small letters) appears to a pos value usually close to POS.
        // we encode as a delta.
        else if (dict_id.num == dict_id_OPTION_mc && field[3] == 'i') {
            if (vb->mc_did_i == DID_I_NONE)
                vb->mc_did_i = mtf_get_ctx_by_dict_id (vb->mtf_ctx, vb->dict_id_to_did_i_map, &vb->num_dict_ids, NULL, dict_id, SEC_SAM_OPTNL_SF_DICT)->did_i;

            seg_pos_field ((VBlockP)vb, vb->last_pos, vb->mc_did_i, SEC_SAM_OPTNL_SF_B250, &field[5], field_len-5, vb_line_i, "mc");
        }

        // E2 - SEQ data (note: E2 doesn't have a dictionary)
        else if (dict_id.num == dict_id_OPTION_E2) { 
            ZipDataLineSAM *dl = DATA_LINE (vb, vb_line_i);
            ASSERT (field_len-5 == dl->seq_len, 
                    "Error in %s: Expecting E2 data to be of length %u as indicated by CIGAR, but it is %u. E2=%.*s",
                    txt_name, dl->seq_len, field_len-5, field_len-5, &field[5]);

            dl->e2_data_start = &field[5] - vb->txt_data.data;
            vb->txt_section_bytes[SEC_SAM_SEQ_DATA] += dl->seq_len + 1; // +1 for \t
        }
        // U2 - QUAL data (note: U2 doesn't have a dictionary)
        else if (dict_id.num == dict_id_OPTION_U2) {
            ZipDataLineSAM *dl = DATA_LINE (vb, vb_line_i);
            ASSERT (field_len-5 == dl->seq_len, 
                    "Error in %s: Expecting U2 data to be of length %u as indicated by CIGAR, but it is %u. E2=%.*s",
                    txt_name, dl->seq_len, field_len-5, field_len-5, &field[5]);

            dl->u2_data_start = &field[5] - vb->txt_data.data;
            vb->txt_section_bytes[SEC_SAM_QUAL_DATA] += dl->seq_len + 1; // +1 for \t
        }
        // All other subfields - have their own dictionary
        else        
            seg_one_subfield ((VBlock *)vb, &field[5], field_len-5, vb_line_i, dict_id, SEC_SAM_OPTNL_SF_B250,
                             (field_len-5) + 1); // +1 for \t
    }

    if (separator == '\n') 
        vb->txt_section_bytes[SEC_SAM_OPTNL_SF_B250]--; // \n is already accounted for in SEC_SAM_OPTIONAL_B250
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
    uint32_t rname_node_index = seg_sam_one_field (vb, field_start, field_len, vb_line_i, SAM_RNAME, NULL);

    random_access_update_chrom (vb_, vb_line_i, rname_node_index);

    // POS - two options:
    // 1. if RNAME is the same as the previous line - store the delta in SAM_POS dictionary
    // 2. If its a different RNAME - add to vb->random_pos_data (not a dictionary)
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "POS");
    
    if (rname_node_index != vb->last_rname_node_index) { // different rname than prev line - store full pos in random
        buf_alloc_more (vb, &vb->random_pos_data, field_len+1, 0, char, 2);
        buf_add (&vb->random_pos_data, field_start, field_len+1); // including the \t
        vb->txt_section_bytes[SEC_SAM_RAND_POS_DATA] += field_len+1;
        vb->last_pos = seg_pos_snip_to_int (field_start, vb_line_i, "POS");
    }
    else // same rname - do a delta
        vb->last_pos = seg_pos_field (vb_, vb->last_pos, SAM_POS, SEC_SAM_POS_B250, field_start, field_len, vb_line_i, "POS");

    random_access_update_pos (vb_, vb->last_pos);

    vb->last_rname_node_index = rname_node_index;

    // MAPQ
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "MAPQ");
    seg_sam_one_field (vb, field_start, field_len, vb_line_i, SAM_MAPQ, NULL);

    // CIGAR - if CIGAR is "*" and we wait to get the length from SEQ or QUAL
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "CIGAR");
    dl->seq_len = seg_sam_seq_len_from_cigar (field_start, field_len);
    if (dl->seq_len) seg_sam_one_field (vb, field_start, field_len, vb_line_i, SAM_CIGAR, NULL); // not "*" - all good!

    // RNEXT - add to RNAME dictionary
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "RNEXT");
    uint32_t rnext_node_index = seg_sam_one_field (vb, field_start, field_len, vb_line_i, SAM_RNAME, NULL);
    
    bool rnext_same_as_rname = (field_len == 1 && (field_start[0] == '=' || field_start[0] == '*')) || 
                                (rnext_node_index == rname_node_index);
    // PNEXT - 2.5 options:
    // 1. if RNEXT is the same as RNAME or unknown (i.e. either it is "=" or "*" or the string is the same) - store the delta
    //    vs. POS in SAM_PNEXT dictionary
    //    1A. If 1, but PNEXT is "0" and hence "unavailable" per SAM specification - store "*" in the SAM_PNEXT dictionary
    // 2. If RNEXT and RNAME are different - add to vb->random_pos_data (not a dictionary)
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "PNEXT");
        
    if (rnext_same_as_rname) { // RNAME and RNEXT are the same - store delta in dictionary
        if (field_len == 1 && *field_start == '0') // PNEXT is "unavailable"
            seg_sam_one_field (vb, "*", 1, vb_line_i, SAM_PNEXT, NULL); 
        else
            seg_pos_field (vb_, vb->last_pos, SAM_PNEXT, SEC_SAM_PNEXT_B250, field_start, field_len, vb_line_i, "PNEXT");
    }

    else { // RNAME and RNEXT differ - store in random_pos
        buf_alloc_more (vb, &vb->random_pos_data, field_len+1, 0, char, 2);
        buf_add (&vb->random_pos_data, field_start, field_len+1); // including the \t
        vb->txt_section_bytes[SEC_SAM_RAND_POS_DATA] += field_len+1;
    }

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

    vb->longest_line_len = MAX (vb->longest_line_len, (next_field - field_start_line));

    return next_field;
}
