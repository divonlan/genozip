// ------------------------------------------------------------------
//   segregate.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "segregate.h"
#include "vb.h"
#include "move_to_front.h"
#include "vcf_header.h"
#include "endianness.h"
#include "random_access.h"

// returns the node index
static int32_t seg_one_field (VariantBlock *vb, const char *str, unsigned len, unsigned vcf_line_i, VcfFields f, 
                              bool *is_new) // optional out
{
    uint32_t *this_field_section = (uint32_t *)vb->mtf_ctx[f].mtf_i.data;

    MtfContext *ctx = &vb->mtf_ctx[f];
    MtfNode *node;
    int32_t node_index = mtf_evaluate_snip (vb, ctx, str, len, &node, is_new);
    this_field_section[vb->mtf_ctx[f].mtf_i.len++] = node_index;

    vb->vcf_section_bytes[SEC_CHROM_B250 + f*2] += len + 1;

    return node_index;
}

// returns true if this line has the same chrom as this VB, or if it is the first line
static void seg_chrom_field (VariantBlock *vb, const char *chrom_str, unsigned chrom_str_len, unsigned vcf_line_i)
{
    int32_t chrom_node_index = seg_one_field (vb, chrom_str, chrom_str_len, vcf_line_i, CHROM, NULL);
    uint32_t vb_line_i = vcf_line_i - vb->first_line;
    
    // case: first vb line or change in chrom - start new entry
    if (!vb_line_i || chrom_node_index != random_access_get_last_chrom_node_index(vb)) 
        random_access_new_entry (vb, vb_line_i, chrom_node_index);
}

static void seg_pos_field (VariantBlock *vb, const char *pos_str, unsigned pos_len, unsigned vcf_line_i)
{
    // scan by ourselves - hugely faster the sscanf
    int64_t this_pos_64=0; // long long so we can test for overflow
    const char *s; for (s=pos_str; *s != '\t'; s++)
        this_pos_64 = this_pos_64 * 10 + (*s - '0');

    int32_t this_pos = (int32_t)this_pos_64;

    ASSERT (this_pos_64 >= 1 && this_pos_64 <= 0x7fffffff, 
            "Error: Invalid POS in line %u - value should be between 1 and %u, but found %u", vcf_line_i, 0x7fffffff, this_pos);
    
    int32_t pos_delta = this_pos - vb->last_pos;
    
    // print our string without expensive sprintf
    char pos_delta_str[50], reverse_pos_delta_str[50];
    
    bool negative = (pos_delta < 0);
    if (negative) pos_delta = -pos_delta;

    // create reverse string
    unsigned len=0; 
    if (pos_delta) { 
        while (pos_delta) {
            reverse_pos_delta_str[len++] = '0' + (pos_delta % 10);
            pos_delta /= 10;
        }
        if (negative) reverse_pos_delta_str[len++] = '-';

        // reverse it
        for (unsigned i=0; i < len; i++) pos_delta_str[i] = reverse_pos_delta_str[len-i-1];
    }
    else { //pos_delta==0
        pos_delta_str[0] = '0';
        len = 1;
    }

    seg_one_field (vb, pos_delta_str, len, vcf_line_i, POS, NULL);

    vb->vcf_section_bytes[SEC_POS_B250] += pos_len - len; // re-do the calculation - seg_one_field doesn't do it good in our case

    vb->last_pos = this_pos;

    random_access_update_last_entry (vb, this_pos);
}

// traverses the FORMAT field, gets ID of subfield, and moves to the next subfield
DictIdType seg_get_format_subfield (const char **str, uint32_t *len, // remaining length of line
                                    unsigned vcf_line_i) 
{
    DictIdType subfield = { 0 };

    unsigned i=0; for (; i < *len; i++) {
        if ((*str)[i] != ':' && (*str)[i] != '\t' && (*str)[i] != '\n') { // note: in vcf files, we will see \t at end of FORMAT, but in variant data sections of piz files, we will see a \n
            
            if (i < DICT_ID_LEN) // note, a subfield might be longer that DICT_ID_LEN, but we copy only the first DICT_ID_LEN characters
                subfield.id[i] = (*str)[i];
        } 
        else break;
    }

    *str += i+1;
    *len -= i+1;
    return subfield; // never reach here, just to avoid a compiler warning
}

static void seg_format_field (VariantBlock *vb, DataLine *dl, 
                              const char *field_start, int field_len, 
                              unsigned vcf_line_i)
{
    const char *str = field_start;
    int len = field_len;
    SubfieldMapperZip format_mapper;
    memset (&format_mapper, 0, sizeof (format_mapper));

    if (!global_num_samples) return; // if we're not expecting any samples, we ignore the FORMAT field

    ASSERT (field_len, "Error: missing FORMAT field in line %u", vcf_line_i);

    if (field_len >= 2 && str[0] == 'G' && str[1] == 'T' && (field_len == 2 || str[2] == ':')) // GT field in FORMAT columns - must always appear first per VCF spec (if at appears)
        dl->has_haplotype_data = true; 

    // we have genotype data, iff FORMAT is not "GT" or ""
    if (field_len > 2 || (!dl->has_haplotype_data && field_len > 0)) {
        dl->has_genotype_data = true; // no GT subfield, or subfields in addition to GT

        if (dl->has_haplotype_data) {
            str +=3; // skip over GT
            len -= 3;
        }

        do {
            ASSERT (format_mapper.num_subfields < MAX_SUBFIELDS, 
                    "Error: FORMAT field in line %u has too many subfields, the maximum allowed is %u (excluding GT)", vcf_line_i, MAX_SUBFIELDS);

            DictIdType subfield = seg_get_format_subfield (&str, (unsigned *)&len, vcf_line_i);

            ASSERT (dict_id_is_gtdata_subfield (subfield), 
                    "Error: string %.*s in the FORMAT field of line=%u is not a legal subfield", DICT_ID_LEN, subfield.id, vcf_line_i);

            MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, &vb->num_dict_ids, &vb->num_format_subfields, subfield, SEC_FRMT_SUBFIELD_DICT);
            
            format_mapper.did_i[format_mapper.num_subfields++] = ctx ? ctx->did_i : (uint8_t)NIL;
        } 
        while (str[-1] != '\t' && str[-1] != '\n' && len > 0);
    }

    if (dl->has_genotype_data)
        buf_alloc (vb, &vb->line_gt_data, format_mapper.num_subfields * global_num_samples * sizeof(uint32_t), 1, "line_gt_data", vcf_line_i); 

    bool is_new;
    int32_t node_index = seg_one_field (vb, field_start, field_len, vcf_line_i, FORMAT, &is_new);

    dl->format_mtf_i = node_index;

    // if this is a new format - add mapper
    if (is_new) {
        ASSERT (node_index == vb->format_mapper_buf.len, 
                "Error: node_index=%u different than vb->format_mapper_buf.len=%u", node_index, vb->format_mapper_buf.len);

        vb->format_mapper_buf.len++;
        buf_alloc (vb, &vb->format_mapper_buf, vb->format_mapper_buf.len * sizeof (SubfieldMapperZip), 2, "format_mapper_buf", 0);
    }

    // it is possible that the mapper is not set yet even though not new - if the node is from a previous VB and
    // we have not yet encountered in node in this VB
    ((SubfieldMapperZip *)vb->format_mapper_buf.data)[node_index] = format_mapper;
}

static void seg_info_field (VariantBlock *vb, DataLine *dl, const char *info_str, unsigned info_len, unsigned vcf_line_i)
{
    #define MAX_INFO_NAMES_LEN 1000 // max len of just the names string, without the data eg "INFO1=INFO2=INFO3="
    char iname[MAX_INFO_NAMES_LEN];
    unsigned iname_len = 0;
    const char *this_name = info_str;
    unsigned this_name_len = 0;
    const char *this_value;
    unsigned this_value_len=0;
    unsigned sf_i=0;

    // count infos
    SubfieldMapperZip iname_mapper;
    memset (&iname_mapper, 0, sizeof (iname_mapper));

    for (unsigned i=0; i < info_len; i++) 
        if (info_str[i] == '=') iname_mapper.num_subfields++;
    
    // get name / value pairs - and insert values to the "name" dictionary
    bool reading_name = true;
    for (unsigned i=0; i < info_len + 1; i++) {

        char c = (i==info_len) ? ';' : info_str[i]; // add an artificial ; at the end of the INFO data

        if (reading_name) {
            iname[iname_len++] = c; // info names inc. the =. the = terminats each name
            ASSERT (iname_len <= MAX_INFO_NAMES_LEN, "Error: INFO field too long in line %u", vcf_line_i);

            if (c == '=') {  // end of name

                ASSERT (this_name_len > 0, "Error: INFO field in line %u, contains a = without a preceding subfield name", vcf_line_i);

                ASSERT (this_name[0] >= 64 && this_name[0] <= 127, 
                        "Error: INFO field in line %u, contains a name %.*s starting with an illegal character", vcf_line_i, this_name_len, this_name);

                reading_name = false; 
                this_value = &info_str[i+1]; 
                this_value_len = 0;
            }
            else if (c == ';' && i==info_len) { // name without value, can happen for example if the field is "." only
                iname_len--; // remove ;
                continue;
            }
            
            else
                this_name_len++; // don't count the = or ; in the len
        }
        else {
            if (c == ';') { // end of value

                ASSERT (this_value_len > 0, 
                        "Error: INFO field in line %u, subfield %.*s, does not contain a value", vcf_line_i, this_name_len, this_name);

                // find (or create) an MTF context (= a dictionary) for this name
                DictIdType dict_id = dict_id_info_subfield (dict_id_make (this_name, this_name_len));

                // find which DictId (did_i) this subfield belongs to
                MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, &vb->num_dict_ids, &vb->num_info_subfields, dict_id, SEC_INFO_SUBFIELD_DICT);
                iname_mapper.did_i[sf_i] = ctx ? ctx->did_i : (uint8_t)NIL;

                // allocate memory if needed (check before calling buf_alloc - we're in a tight loop)
                Buffer *mtf_i_buf = &ctx->mtf_i;
                if (!buf_is_allocated(mtf_i_buf) || (mtf_i_buf->len + 1) * sizeof(uint32_t) > mtf_i_buf->size)
                    buf_alloc (vb, mtf_i_buf, MIN (global_max_lines_per_vb, mtf_i_buf->len + 1) * sizeof (uint32_t),
                               1.15, "mtf_ctx->mtf_i", ctx->dict_section_type);

                MtfNode *sf_node;

                ((uint32_t *)mtf_i_buf->data)[mtf_i_buf->len++] = 
                    mtf_evaluate_snip (vb, ctx, this_value, this_value_len, &sf_node, NULL);

                vb->vcf_section_bytes[SEC_INFO_SUBFIELD_B250] += (this_value_len+1); // including the separator (; or \n)

                reading_name = true;  // end of value - move to the next time
                this_name = &info_str[i+1]; // move to next field in info string
                this_name_len = 0;
                sf_i++;
            }
            else  
                this_value_len++;
        }
    }

    // now insert the info names - a snip is a string that looks like: "INFO1=INFO2=INFO3="
    // 1. find it's mtf_i (and add to dictionary if a new name)
    // 2. place mtf_i in INFO section of this VB
    uint32_t *info_field_mtf_i = (uint32_t *)vb->mtf_ctx[INFO].mtf_i.data;
    MtfContext *info_ctx = &vb->mtf_ctx[INFO];
    MtfNode *node;
    bool is_new;
    uint32_t node_index = mtf_evaluate_snip (vb, info_ctx, iname, iname_len, &node, &is_new);
    info_field_mtf_i[vb->mtf_ctx[INFO].mtf_i.len++] = node_index;

    // if this is a totally new subfield (first time in this file) - make a new SubfieldMapperZip for it.
    if (is_new) {   
        ASSERT (node_index == vb->iname_mapper_buf.len, "Error: node_index=%u different than vb->iname_mapper_buf.len=%u", node_index, vb->iname_mapper_buf.len);
    
        vb->iname_mapper_buf.len++;
        buf_alloc (vb, &vb->iname_mapper_buf, MAX (100, vb->iname_mapper_buf.len) * sizeof (SubfieldMapperZip), 1.5, "iname_mapper_buf", 0);
    }

    // it is possible that the iname_mapper is not set yet even though not new - if the node is from a previous VB and
    // we have not yet encountered in node in this VB
    ((SubfieldMapperZip *)vb->iname_mapper_buf.data)[node_index] = iname_mapper;

    dl->info_mtf_i = node_index;

    vb->vcf_section_bytes[SEC_INFO_B250] += iname_len + !sf_i; // this includes all the = (the \t is included with the values, except if there are none)
}

static void seg_increase_ploidy_one_line (VariantBlock *vb, Buffer *line, unsigned new_ploidy, unsigned num_samples)
{
    // copy the haplotypes backwards (to avoid overlap), padding with '*'
    for (int sam_i = num_samples-1; sam_i >= 0; sam_i--) {

        int ht_i=new_ploidy-1 ; for (; ht_i >= vb->ploidy; ht_i--) 
            line->data[sam_i * new_ploidy + ht_i] = '*';

        for (; ht_i >= 0; ht_i--)
            line->data[sam_i * new_ploidy + ht_i] = line->data[sam_i * vb->ploidy + ht_i];
    
        // note: no change in phase type of previous rows, even if they were PHASE_TYPE_HAPLO
    }
}

// increase ploidy of the previous lines, if higher ploidy was encountered
static void seg_increase_ploidy (VariantBlock *vb, unsigned new_ploidy, unsigned vcf_line_i, unsigned sample_i)
{
    unsigned vb_line_i = vcf_line_i - vb->first_line;

    // increase ploidy of all previous lines. at this point, dl->haplotype_data is overlaid dl->line_data
    for (unsigned i=0; i < vb_line_i; i++) {
        DataLine *dl = &vb->data_lines[i];

        dl->haplotype_data.len = global_num_samples * new_ploidy; // increase len

        // note: line_ht_data is overlaid at the end of dl->line
        unsigned new_total_len = dl->haplotype_data.len + dl->genotype_data.len + dl->phase_data.len;            
        buf_alloc (vb, &dl->line, new_total_len, 1, "dl->line", vcf_line_i);

        // incease ploidy if this row has haplotypes already. if not, will add it at the end.
        if (dl->has_haplotype_data)
            seg_increase_ploidy_one_line (vb, &dl->haplotype_data, new_ploidy, global_num_samples);
    }

    // increase ploidy in all previous samples of this line
    if (sample_i) seg_increase_ploidy_one_line (vb, &vb->line_ht_data, new_ploidy, sample_i);

    vb->ploidy = new_ploidy;
}

static void seg_haplotype_area (VariantBlock *vb, DataLine *dl, const char *str, unsigned len, unsigned vcf_line_i, unsigned sample_i,
                                bool is_vcf_string)
{
    // check ploidy
    unsigned ploidy=1;
    for (unsigned i=1; i<len-1; i++)
        if (str[i] == '|' || str[i] == '/') ploidy++;

    if (is_vcf_string) {
        vb->vcf_section_bytes[SEC_PHASE_DATA]     += ploidy-1;
        vb->vcf_section_bytes[SEC_HAPLOTYPE_DATA] += ploidy;
    }

    ASSERT (ploidy <= MAX_PLOIDY, "Error: ploidy=%u exceeds the maximum of %u in line %u", ploidy, MAX_PLOIDY, vcf_line_i);
    
    // if the ploidy of this line is bigger than the ploidy of the data in this VB so far, then
    // we have to increase ploidy of all the haplotypes read in in this VB so far. This can happen for example in 
    // the X chromosome if initial samples are male with ploidy=1 and then a female sample with ploidy=2
    if (vb->ploidy && ploidy > vb->ploidy) 
        seg_increase_ploidy (vb, ploidy, vcf_line_i, sample_i);

    if (!vb->ploidy) vb->ploidy = ploidy; // very first sample in the vb

    // note - ploidy of this sample might be smaller than vb->ploidy (eg a male sample in an X chromosesome that was preceded by a female sample)
    if (sample_i == 0) {
        buf_alloc (vb, &vb->line_ht_data, vb->ploidy * global_num_samples, 1, "line_ht_data", vcf_line_i);
        dl->phase_type = (vb->ploidy==1 ? PHASE_HAPLO : PHASE_UNKNOWN);
    }

    char *ht_data = &vb->line_ht_data.data[vb->ploidy * sample_i];

    PhaseType ht0_phase_type = PHASE_UNKNOWN;
    for (unsigned ht_i=0; ht_i < ploidy; ht_i++) {

        char ht = *(str++); 
        len--;

        ASSERT ((ht >= '0' && ht <= '9') || ht == '.' || ht == '*',
                "Error: invalid VCF file - line %u - expecting an allele in a sample to be a number 0-9 or . , but seeing %c", vcf_line_i, *str);

        // single-digit allele numbers
        ht_data[ht_i] = ht;

        if (!len) break;

        // handle 2-digit allele numbers
        if (ht != '.' && *str >= '0' && *str <= '9') {
            unsigned allele = 10 * (ht-'0') + (*(str++) - '0');
            len--;

            // make sure there isn't a 3rd digit
            ASSERT (!len || *str < '0' || * str > '9', "Error: VCF file - line %u sample %u - genozip currently supports only alleles up to 99", vcf_line_i, sample_i+1);

            ht_data[ht_i] = '0' + allele; // use ascii 48->147

            vb->vcf_section_bytes[SEC_HAPLOTYPE_DATA]++;
        }

        // get phase
        if (ploidy > 1 && ht_i < ploidy-1) {
            
            PhaseType cell_phase_type = (PhaseType)*(str++);
            len--;

            ASSERT (cell_phase_type != ' ', "Error: invalid VCF file - line %u - expecting a tab or newline after sample %u but seeing a space", vcf_line_i, sample_i+1);
            ASSERT (cell_phase_type == '|' || cell_phase_type == '/', "Error: invalid VCF file - line %u - unable to parse sample %u: expecting a | or / but seeing %c", vcf_line_i, sample_i+1, cell_phase_type);

            // deal with phase - only at the first separator eg 1|1|0|1
            if (ht_i==0) { 
                if (cell_phase_type == dl->phase_type) {} // do nothing

                else if ((cell_phase_type == '|' || cell_phase_type == '/') && 
                         (dl->phase_type == PHASE_UNKNOWN || dl->phase_type == PHASE_HAPLO))
                    dl->phase_type = cell_phase_type;

                else if ((dl->phase_type == '|' && cell_phase_type == '/') || 
                         (dl->phase_type == '/' && cell_phase_type == '|')) {
                    dl->phase_type = PHASE_MIXED_PHASED;
                
                    buf_alloc (vb, &vb->line_phase_data, global_num_samples, 1, "line_phase_data", vcf_line_i);

                    // fill in the so-far uniform phase (which is different than the current one)
                    memset (vb->line_phase_data.data, (cell_phase_type == '|' ? '/' : '|'), sample_i);
                    vb->line_phase_data.data[sample_i] = (char)cell_phase_type;
                }
                else if (dl->phase_type == PHASE_MIXED_PHASED)
                    vb->line_phase_data.data[sample_i] = (char)cell_phase_type;

                ht0_phase_type = cell_phase_type;
            }
            // subsequent seperators - only make sure they're consistent
            else
                ASSERT (cell_phase_type==ht0_phase_type, "Error: invalid VCF file - line %u - unable to parse sample %u: inconsistent phasing symbol '|' '/'", vcf_line_i, sample_i+1);
        }
    } // for characters in a sample

    // if the ploidy of the sample is lower than vb->ploidy, set missing ht as '*' ("ploidy padding")
    if (ploidy != vb->ploidy) {
        
        for (unsigned ht_i=ploidy; ht_i < vb->ploidy; ht_i++) 
            ht_data[ht_i] = '*';
    }

    if (ploidy==1 && vb->ploidy > 1 && dl->phase_type == PHASE_MIXED_PHASED)
        vb->line_phase_data.data[sample_i] = (char)PHASE_HAPLO;
}

// returns length of the snip, not including the separator
static inline unsigned seg_snip_len_tnc (const char *snip)
{
    const char *s; for (s=snip ; *s != '\t' && *s != ':' && *s != '\n'; s++);
    return s - snip;
}

static void seg_genotype_area (VariantBlock *vb, DataLine *dl, 
                               const char *cell_gt_data, 
                               unsigned cell_gt_data_len,  // not including the \t or \n 
                               unsigned vcf_line_i,
                               bool is_vcf_string)
{
    uint32_t *next = (uint32_t *)(&vb->line_gt_data.data[vb->line_gt_data.len]);
    SubfieldMapperZip *format_mapper = &((SubfieldMapperZip *)vb->format_mapper_buf.data)[dl->format_mtf_i];

    bool end_of_cell = !cell_gt_data_len;
    for (unsigned sf=0; sf < format_mapper->num_subfields; sf++) { // iterate on the order as in the line

        // move next to the beginning of the subfield data, if there is any
        unsigned len = end_of_cell ? 0 : seg_snip_len_tnc (cell_gt_data);

        MtfNode *node;
        int32_t node_index = mtf_evaluate_snip (vb, MAPPER_CTX (format_mapper, sf), cell_gt_data, len, &node, NULL);
        *(next++) = node_index;

        if (node_index != WORD_INDEX_MISSING_SF) // don't skip the \t if we might have more missing subfields
            cell_gt_data += len + 1; // skip separator too

        end_of_cell = end_of_cell || cell_gt_data[-1] != ':'; // a \t or \n encountered
    }
    ASSERT0 (end_of_cell, "Error: invalid reading of genotype data");

    vb->line_gt_data.len += format_mapper->num_subfields * sizeof(uint32_t);  // len is number of bytes

    if (is_vcf_string)
        // size including : (if we have both ht and gt), but not including \t which goes into SEC_STATS_HT_SEPERATOR
        vb->vcf_section_bytes[SEC_GENOTYPE_DATA] += cell_gt_data_len + (dl->has_haplotype_data && dl->has_genotype_data);
}

static inline const char *seg_get_next_item (const char *str, unsigned *str_len, bool allow_newline, bool allow_tab, bool allow_colon, unsigned vcf_line_i, // line in vcf file,
                                             unsigned *len, char *separator, const char *item_name) // out
{
    unsigned i=0; for (; i < *str_len; i++)
        if ((allow_tab     && str[i] == '\t') ||
            (allow_colon   && str[i] == ':')  ||
            (allow_newline && str[i] == '\n')) {
                *len = i;
                *separator = str[i];
                *str_len -= i+1;
                return str + i+1; // beyond the separator
        }
        else if ((!allow_tab     && str[i] == '\t') ||  // note: a colon with allow_colon=false is not an error, its just part of the string rather than being a separator
                 (!allow_newline && str[i] == '\n')) break;
            
    ASSERT (*str_len, "Error: missing %s field in line %u", item_name, vcf_line_i);

    ASSERT (str[i] != '\t' || strcmp (item_name, vcf_field_names[INFO]), 
           "Error: while segmenting %s in line %u: expecting a NEWLINE after the INFO field, because this VCF file has no samples (individuals) declared in the header line",
            item_name, vcf_line_i);

    ABORT ("Error: while segmenting %s in line %u: expecting a %s %s %s after \"%.*s\"", 
            item_name, vcf_line_i, 
            allow_newline ? "NEWLINE" : "", allow_tab ? "TAB" : "", allow_colon ? "\":\"" : "", 
            MIN (i, 1000), str);

    return 0; // avoid compiler warning - never reaches here
}

// in some real-world files I encountered have too-short lines due to human errors. we pad them
static void seg_add_samples_missing_in_line (VariantBlock *vb, DataLine *dl, unsigned *gt_line_len, 
                                             unsigned num_samples, unsigned vcf_line_i)
{
    ASSERTW (false, "Warning: the number of samples in line %u is %u, different than the VCF column header line which has %u samples",
             vcf_line_i, num_samples, global_num_samples);

    for (; num_samples < global_num_samples; num_samples++) {

        if (dl->has_haplotype_data) {
            // '*' (haplotype padding) with ploidy 1
            seg_haplotype_area (vb, dl, "*", 1, vcf_line_i, num_samples, false);
        }

        if (dl->has_genotype_data) {
            seg_genotype_area (vb, dl, NULL, 0, vcf_line_i, false);
            (*gt_line_len)++; // adding the WORD_INDEX_MISSING_SF
        }
    }
}

/* split the data line into sections - 
   1. variant data - each of 9 fields (CHROM to FORMAT) is a section
   2. genotype data (except the GT subfield) - is a section
   3. haplotype data (the GT subfield) - a string of eg 1 or 0 (no whitespace) - ordered by the permutation order
   4. phase data - only if MIXED phase - a string of | and / - one per sample
*/
static void seg_data_line (VariantBlock *vb, /* may be NULL if testing */
                           DataLine *dl, 
                           unsigned vcf_line_i) // line in original VCF file
{
    dl->phase_type = PHASE_UNKNOWN;

    unsigned sample_i = 0;
    unsigned gt_line_len=0;

    const char *field_start, *next_field;
    unsigned field_len;
    char separator;

    unsigned len = dl->line.len;

    // CHROM
    field_start = dl->line.data;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vcf_line_i, &field_len, &separator, "CHROM");
    seg_chrom_field (vb, field_start, field_len, vcf_line_i);

    // POS
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vcf_line_i, &field_len, &separator, "POS");
    seg_pos_field (vb, field_start, field_len, vcf_line_i);;

    // ID
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vcf_line_i, &field_len, &separator, "ID");
    seg_one_field (vb, field_start, field_len, vcf_line_i, ID, NULL);

    // REF + ALT
    // note: we treat REF+\t+ALT as a single field because REF and ALT are highly corrected, in the case of SNPs:
    // e.g. GG has a probability of 0 and GC has a higher probability than GA.
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vcf_line_i, &field_len, &separator, "REF");

    unsigned alt_len;
    next_field = seg_get_next_item (next_field, &len, false, true, false, vcf_line_i, &alt_len, &separator, "ALT");
    seg_one_field (vb, field_start, field_len+alt_len+1, vcf_line_i, REFALT, NULL);

    // QUAL
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vcf_line_i, &field_len, &separator, "QUAL");
    seg_one_field (vb, field_start, field_len, vcf_line_i, QUAL, NULL);

    // FILTER
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vcf_line_i, &field_len, &separator, "FILTER");
    seg_one_field (vb, field_start, field_len, vcf_line_i, FILTER, NULL);

    // INFO
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, global_num_samples==0, global_num_samples>0, false, vcf_line_i, &field_len, &separator, "INFO");
    seg_info_field (vb, dl, field_start, field_len, vcf_line_i);

    if (separator != '\n') {

        // FORMAT
        field_start = next_field;
        next_field = seg_get_next_item (field_start, &len, true, true, false, vcf_line_i, &field_len, &separator, "FORMAT");
        seg_format_field (vb, dl, field_start, field_len, vcf_line_i);

        ASSERT (separator == '\n' || dl->has_genotype_data || dl->has_haplotype_data, 
                "Error: expecting line %u to end as it has no genotype or haplotype data, but it is not", vcf_line_i);

        // 0 or more samples
        while (separator != '\n') {
            
            // get haplotype data
            bool has_genotype_data = dl->has_genotype_data;
            if (dl->has_haplotype_data) { // FORMAT declares GT, we may have it or not

                field_start = next_field;
                next_field = seg_get_next_item (field_start, &len, true, true, dl->has_genotype_data, vcf_line_i, &field_len, &separator, "GT");
                seg_haplotype_area (vb, dl, field_start, field_len, vcf_line_i, sample_i, true);

                if (separator != ':' && has_genotype_data) {
                    // missing genotype data despite being declared in FORMAT - this is permitted by VCF spec
                    has_genotype_data = false;
                    seg_genotype_area (vb, dl, NULL, 0, vcf_line_i, false);
                    gt_line_len++; // adding the WORD_INDEX_MISSING_SF
                }
            }

            if (has_genotype_data) { // FORMAT declares other subfields, we may have them or not
                field_start = next_field;
                next_field = seg_get_next_item (field_start, &len, true, true, false, vcf_line_i, &field_len, &separator, "Non-GT");
                ASSERT (field_len, "Error: invalid VCF file - expecting sample data for sample # %u on line %u, but found a tab character", 
                        sample_i+1, vcf_line_i);

                seg_genotype_area (vb, dl, field_start, field_len, vcf_line_i, true);
                gt_line_len += field_len + 1; // including the \t or \n
            }

            sample_i++;

            vb->vcf_section_bytes[SEC_STATS_HT_SEPERATOR]++; // the \t or \n following a sample

            ASSERT (sample_i < global_num_samples || separator == '\n',
                    "Error: invalid VCF file - expecting a newline after the last sample (sample #%u) on line %u", global_num_samples, vcf_line_i);
        }
    }

    // in some real-world files I encountered have too-short lines due to human errors. we pad them
    if (sample_i < global_num_samples) 
        seg_add_samples_missing_in_line (vb, dl, &gt_line_len, sample_i, vcf_line_i);

    // update lengths
    if (dl->has_haplotype_data) {
        vb->line_ht_data.len = global_num_samples * vb->ploidy;

        if (dl->phase_type == PHASE_MIXED_PHASED) 
            vb->line_phase_data.len = global_num_samples;
    } else 
        vb->line_ht_data.len = 0;

    vb->max_gt_line_len = MAX (vb->max_gt_line_len, gt_line_len);

    // now, overlay the data over the line memory so we can re-use our working buffers vb->line_* for the next line
    unsigned total_len = vb->line_gt_data.len + vb->line_phase_data.len + vb->line_ht_data.len;
    buf_alloc (vb, &dl->line, total_len, 1, "dl->line", vcf_line_i);

    unsigned offset_in_line = 0;

    // note we put haplotype_data at the end to ease realloc in case of ploidy increase
    if (dl->has_genotype_data) 
        buf_overlay (&dl->genotype_data, &dl->line, &vb->line_gt_data, &offset_in_line, "dl->genotype_data", vcf_line_i);

    if (dl->has_haplotype_data && dl->phase_type == PHASE_MIXED_PHASED)
        buf_overlay (&dl->phase_data, &dl->line, &vb->line_phase_data, &offset_in_line, "dl->phase_data", vcf_line_i);    

    if (dl->has_haplotype_data) {
        buf_overlay (&dl->haplotype_data, &dl->line, &vb->line_ht_data, &offset_in_line, "dl->haplotype_data", vcf_line_i);
        
        if (flag_show_alleles)
            printf ("%.*s\n", dl->haplotype_data.len, dl->haplotype_data.data);
    }

    buf_free (&vb->line_ht_data);
    buf_free (&vb->line_phase_data);
    buf_free (&vb->line_gt_data);
}

static void seg_update_vb_from_dl (VariantBlock *vb, DataLine *dl)
{
    // update block data
    vb->has_genotype_data = vb->has_genotype_data || dl->has_genotype_data;

    vb->has_haplotype_data = vb->has_haplotype_data || dl->has_haplotype_data;
    
    if (vb->phase_type == PHASE_UNKNOWN) 
        vb->phase_type = dl->phase_type;
    
    else if ((vb->phase_type == PHASE_PHASED     && dl->phase_type == PHASE_NOT_PHASED) ||
             (vb->phase_type == PHASE_NOT_PHASED && dl->phase_type == PHASE_PHASED) ||
              dl->phase_type  == PHASE_MIXED_PHASED)
        vb->phase_type = PHASE_MIXED_PHASED;    
}

// complete haplotype and genotypes of lines that don't have them, if we found out that some
// other line has them, and hence all lines in the vb must have them
void seg_complete_missing_lines (VariantBlock *vb)
{
    vb->num_haplotypes_per_line = vb->ploidy * global_num_samples;
    
    for (unsigned vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {

        DataLine *dl = &vb->data_lines[vb_line_i];                                

        if (vb->has_haplotype_data && !dl->has_haplotype_data) {
            // realloc line which is the buffer on which haplotype data is overlaid
            buf_alloc (vb, &dl->line, dl->line.len + vb->num_haplotypes_per_line, 1, dl->line.name, dl->line.param);

            // overlay the haplotype buffer overlaid at the end of the line buffer
            buf_overlay (&dl->haplotype_data, &dl->line, NULL, &dl->line.len, "dl->haplotype_data", vb->first_line + vb_line_i);

            memset (dl->haplotype_data.data, '-', vb->num_haplotypes_per_line); // '-' means missing haplotype - note: ascii 45 (haplotype values start at ascii 48)
            dl->haplotype_data.len = vb->num_haplotypes_per_line;

            // NOTE: we DONT set dl->has_haplotype_data to true bc downstream we still
            // count this row as having no GT field when analyzing gt data
        }

        if (vb->has_genotype_data && !dl->has_genotype_data) {
            // realloc line which is the buffer on which haplotype data is overlaid
            buf_alloc (vb, &dl->line, dl->line.len + global_num_samples * sizeof(uint32_t), 1, dl->line.name, dl->line.param);

            // overlay the haplotype buffer overlaid at the end of the line buffer
            buf_overlay (&dl->genotype_data, &dl->line, NULL, &dl->line.len, "dl->genotype_data", vb->first_line + vb_line_i);
            for (unsigned i=0; i < global_num_samples; i++) 
                ((uint32_t*)dl->genotype_data.data)[i] = WORD_INDEX_MISSING_SF;

            dl->genotype_data.len += global_num_samples * sizeof(uint32_t);
        }
    }
}


// split each lines in this variant block to its components
void seg_all_data_lines (VariantBlock *vb, Buffer *lines_orig /* for testing */)
{
    START_TIMER;
    
    vb->num_dict_ids = MAX (FORMAT+1, vb->num_dict_ids); // first 8 mtf_ctx are reserved for the VCF fields (up to FORMAT) (vb->num_dict_ids might be already higher due to previous VBs)

    for (VcfFields f=CHROM; f <= FORMAT; f++) {
        buf_alloc (vb, &vb->mtf_ctx[f].mtf_i, vb->num_lines * sizeof (uint32_t), 1, "mtf_ctx.mtf_i", f);
     
        strcpy ((char *)vb->mtf_ctx[f].dict_id.id, vcf_field_names[f]); // length of all is <= 7 so fits in id
        vb->mtf_ctx[f].dict_id = dict_id_vardata_field (vb->mtf_ctx[f].dict_id);

        vb->mtf_ctx[f].b250_section_type = SEC_CHROM_B250 + f*2;
        vb->mtf_ctx[f].dict_section_type = SEC_CHROM_DICT + f*2;
    }

    for (unsigned vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {

        //printf ("vb_line_i=%u\n", vb_line_i);
        DataLine *dl = &vb->data_lines[vb_line_i];

        if (lines_orig) buf_copy (vb, &lines_orig[vb_line_i], &dl->line, 1, 0, dl->line.len+1, "lines_orig", vb->variant_block_i); // if testing
        seg_data_line (vb, dl, vb->first_line + vb_line_i);
        seg_update_vb_from_dl (vb, dl);
    }

    if (/*vb->has_genotype_data || */vb->has_haplotype_data)
        seg_complete_missing_lines(vb);

    COPY_TIMER(vb->profile.seg_all_data_lines);
}