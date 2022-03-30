// ------------------------------------------------------------------
//   vcf_tags.c
//   Copyright (C) 2021-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt
//
// Handles renaming INFO and FORMAT tags in DVCF. Tags can be defined explicitly in the VCF header or implicitly in the data lines

#include "genozip.h"
#include "buffer.h"
#include "vcf_private.h"
#include "dict_id.h"
#include "context.h"
#include "file.h"
#include "strings.h"

#define DROP_PREFIX "DROP_"
#define DROP_PREFIX_LEN (sizeof DROP_PREFIX-1)

#define foreach_dtype(d) for (DictIdType one_dtype=((d&DTYPE_VCF_INFO) ? DTYPE_VCF_INFO : DTYPE_VCF_FORMAT) ; one_dtype <= ((d&DTYPE_VCF_FORMAT) ? DTYPE_VCF_FORMAT : DTYPE_VCF_INFO); one_dtype++)
 
static Buffer command_line_tags = EMPTY_BUFFER;

// compares between tags - returns 0 if they are the same 
static inline int tags_cmp (STRp (tag_name_a), DictIdType dtype_a, STRp (tag_name_b), DictIdType dtype_b)
{
    return (tag_name_a_len != tag_name_b_len) ? ((int)tag_name_a_len - (int)tag_name_b_len) // sort by length first, as it is a good differentiator and fast to test
         : (dtype_a != dtype_b)               ? (dtype_a - dtype_b) 
         :                                      memcmp (tag_name_a, STRa (tag_name_b));                                   
}

// ------------------------------------------------------------
// ingesting command line options --dvcf-rename and --dvcf-drop
// ------------------------------------------------------------

// get tag_name and dtype from "FORMAT/ADF"
static void vcf_tags_cmdline_type_and_name_do (STRp(s), rom option,
                                               DictIdType *dtype, STRp (*tag_name)) // out
{
        str_split (s, s_len, 2, '/', left, false);
        ASSINP (n_lefts > 0, "%s: invalid format of tag name: \"%.*s\". Expecting something like \"INFO/AC\" or just \"AC\". See " WEBSITE_DVCF,
                option, s_len, s);

        *tag_name = lefts[n_lefts-1]; 
        *tag_name_len = left_lens[n_lefts-1];

        ASSINP (*tag_name_len > 0, "%s: empty tag name. See: " WEBSITE_DVCF, option);    
        
        // get tag type
        if (n_lefts == 1)                                *dtype = DTYPE_VCF_INFO | DTYPE_VCF_FORMAT;
        else if (str_issame_(STRi(left,0), "INFO", 4))   *dtype = DTYPE_VCF_INFO;
        else if (str_issame_(STRi(left,0), "FORMAT", 6)) *dtype = DTYPE_VCF_FORMAT;
        else ABORTINP ("%s: Invalid tag type \"%.*s\", expecting INFO or FORMAT. See: " WEBSITE_DVCF, option, left_lens[0], lefts[0]);
}
#define vcf_tags_cmdline_type_and_name(name,i,option) \
    STR (tag_name); DictIdType dtype;\
    vcf_tags_cmdline_type_and_name_do (name##s[i], name##_lens[i], (option), &dtype, pSTRa(tag_name))


// add tag from command line (--chain)  
// adds to last tag in buffer, or creates a new tag if different than last one
static void vcf_tags_cmdline_add_attr (rom option, DictIdType dtype, STRp(attr_str), STRp(tag_name), STRp(dest), VcfTagSource source)
{
    // get attr type
    RenameAttr attr=0;
    if      (str_issame_(STRa(attr_str), cSTR("STRAND"))) attr = RA_STRAND;
    else if (str_issame_(STRa(attr_str), cSTR("REFALT"))) attr = RA_REFALT;
    else if (str_issame_(STRa(attr_str), cSTR("TLAFER"))) attr = RA_TLAFER;
    else if (str_issame_(STRa(attr_str), cSTR("ALWAYS"))) attr = RA_ALWAYS;
    else ABORTINP ("%s: Unrecognized attr: \"%.*s\". See: " WEBSITE_DVCF, option, STRf(attr_str));

    ASSINP (dest_len < MAX_TAG_LEN, "tag \"%.*s\" is beyond Genozip's length limit of %u", STRf(dest), MAX_TAG_LEN); // leave room for \0

    Buffer *tags = z_file ? &z_file->apriori_tags : &command_line_tags; // command line tags are added before z_file is created

    // check if tag already exists - searching backwards
    Tag *tag = NULL;
    for (int i=tags->len-1; i >= 0; i--) {
        Tag *candidate = B(Tag, *tags, i);
        if (!tags_cmp (STRa(tag_name), dtype, STRa(candidate->tag_name), candidate->dtype)) {
            tag = candidate;
            break;
        }
    }

    // case: tag not found - create a new one
    if (!tag) { 
        buf_alloc_zero (evb, tags, 1, 100, Tag, 2, z_file ? "z_file->apriori_tags" : "command_line_tags"); 
        
        tag = &BNXT (Tag, *tags);
        STRcpy (tag->tag_name, tag_name);
        tag->dtype  = dtype;
        tag->source = source; 
        STRcpyi(tag->dest, attr, dest); 
    } 
    
    // case: tag exists, but attribute is not set - set it now (note: if attribute is already set, we don't overwrite it)
    else {
        if (!tag->dest_lens[attr]) 
            STRcpyi(tag->dest, attr, dest);

        tag->source = MAX_(tag->source, source); // the more authortative prevails 
    }
}

// add symmetrical attrs to the destination tags
void vcf_tags_cmdline_add_symmetrical (rom option, DictIdType dtype, STRp(attr_str), STRp(tag_name), STRp(dest_tag_name))
{
    // use the existing tag, if one exists
    uint64_t save_len = command_line_tags.len ;

    for (int tag_i=command_line_tags.len-1; tag_i >= 0; tag_i--) { // search backwards as more likely towards the end
        Tag *tag = B(Tag, command_line_tags, tag_i);
        if (!tags_cmp (STRa(tag_name), dtype, STRa (tag->tag_name), tag->dtype)) {
            command_line_tags.len = tag_i+1; // set len, so vcf_tags_cmdline_add_attr will use existing tag
            break;
        }
    }

    // add to new or existing tag
    vcf_tags_cmdline_add_attr (option, dtype, STRa(attr_str), STRa(tag_name), STRa(dest_tag_name), TAG_CMDLINE_DST);
    command_line_tags.len = MAX_(command_line_tags.len, save_len); // restore
}

// eg "FORMAT/ADF:STRAND>ADR|REFALT>ADF2,..."
void vcf_tags_cmdline_rename_option(void)
{
    if (!flag.dvcf_rename) return; // nothing to do

    str_split (flag.dvcf_rename, strlen (flag.dvcf_rename), 0, ',', unit, false);

    for (unsigned unit_i=0; unit_i < n_units; unit_i++) {

        str_split (units[unit_i], unit_lens[unit_i], 2, ':', top, true);
        ASSINP0 (n_tops == 2, "--dvcf-rename: invalid format. See: " WEBSITE_DVCF);

        // get tag_name and dtype from "FORMAT/ADF"
        vcf_tags_cmdline_type_and_name (top, 0, "--dvcf-rename");

        // split to intstructions from "STRAND>ADR|REFALT>DROP_ADF"
        str_split (tops[1], top_lens[1], 0, '|', attr, false);
        ASSINP0 (attr_lens[n_attrs-1] > 0, "--dvcf-rename: missing attrs. See: " WEBSITE_DVCF);    

        // add to tag
        foreach_dtype(dtype) {

            for (unsigned i=0; i < n_attrs; i++) {
                // parse one attr: "REFALT>DROP_ADF"
                str_split (attrs[i], attr_lens[i], 2, '>', item, true);
                ASSINP (n_items == 2 && item_lens[1] > 0, "--dvcf-rename: invalid attr \"%.*s\". See: " WEBSITE_DVCF, attr_lens[i], attrs[i]);

                if (items[0][item_lens[0]-1] == '\\') item_lens[0]--; // remove escape char for escaping >
                    
                vcf_tags_cmdline_add_attr ("--dvcf-rename", one_dtype, STRi(item,0), STRa(tag_name), STRi(item,1), TAG_CMDLINE);
            }

            // add symmetrical attrs to the destination tags
            for (unsigned i=0; i < n_attrs; i++) {
                str_split (attrs[i], attr_lens[i], 2, '>', item, true);
                if (items[0][item_lens[0]-1] == '\\') item_lens[0]--; // remove escape char for escaping >

                vcf_tags_cmdline_add_symmetrical ("--dvcf-rename", one_dtype, STRi(item,0), STRi(item,1), STRa(tag_name));
            }
        }
    }
}

// eg "FORMAT/ADF:STRAND|REFALT,..."
void vcf_tags_cmdline_drop_option (void)
{
    if (!flag.dvcf_drop) return; // nothing to do

    str_split (flag.dvcf_drop, strlen (flag.dvcf_drop), 0, ',', unit, false);

    for (unsigned unit_i=0; unit_i < n_units; unit_i++) {

        str_split (units[unit_i], unit_lens[unit_i], 2, ':', top, true);
        ASSINP0 (n_tops == 2, "--dvcf-drop: invalid format. See: " WEBSITE_DVCF);

        // get tag_name and dtype from "FORMAT/ADF"
        vcf_tags_cmdline_type_and_name (top, 0, "--dvcf-drop");

        // split to intstructions from "STRAND|REFALT"
        str_split (tops[1], top_lens[1], 0, '|', attr, false);
        ASSINP0 (attr_lens[n_attrs-1] > 0, "--dvcf-drop: missing attrs. See: " WEBSITE_DVCF);    

        // destination is DROP_*
        unsigned drop_tag_len = DROP_PREFIX_LEN + tag_name_len;
        char drop_tag[drop_tag_len];
        memcpy (drop_tag, DROP_PREFIX, DROP_PREFIX_LEN);
        memcpy (drop_tag + DROP_PREFIX_LEN, tag_name, tag_name_len);

        foreach_dtype(dtype) {
            // add to tag
            for (unsigned attr_i=0; attr_i < n_attrs; attr_i++) 
                vcf_tags_cmdline_add_attr ("--dvcf-drop", one_dtype, STRi(attr,attr_i), STRa(tag_name), STRa(drop_tag), TAG_CMDLINE);

            // add symmetrical attrs to DROP_tag
            for (unsigned attr_i=0; attr_i < n_attrs; attr_i++) 
                vcf_tags_cmdline_add_attr ("--dvcf-drop", one_dtype, STRi(attr, attr_i), STRa(drop_tag), STRa(tag_name), TAG_CMDLINE_DST);
        }
    }
}

// -------------------------------------------
// prcessing tags while reading the VCF header
// -------------------------------------------

static void vcf_tags_show_rename_tags (void)
{
    if (!z_file || !z_file->apriori_tags.len) {
        iprint0 ("\nThere are no renamed tags defined in the VCF meta-information lines or --dvcf-rename or --dvcf-drop options\n");
        return;
    }

    iprint0 ("\nRenamed tags defined in the VCF meta-information lines or --dvcf-rename or --dvcf-drop options:\n");
    for (unsigned i=0; i < z_file->apriori_tags.len; i++) {
        Tag *tag = B(Tag, z_file->apriori_tags, i);
        rom src_names[] = TAG_SOURCE_NAMES;

        iprintf ("%s/%.*s:\tSrc=%s\tNumber=%.*s\tType=%.*s\t", 
                 DTPT(dtype_names)[tag->dtype], STRf(tag->tag_name), src_names[tag->source], STRf(tag->number), STRf(tag->type));

        if (tag->rendalg_len)
            iprintf (HK_RENDALG_ATTR"%.*s\t", STRf(tag->rendalg));

        for (unsigned attr_i=0; attr_i < NUM_RENAME_ATTRS; attr_i++) 
            if (tag->dest_lens[attr_i]) 
                iprintf ("%s%.*s\t", vcf_header_rename_attrs[attr_i], tag->dest_lens[attr_i], tag->dests[attr_i]);
        
        iprint0 ("\n");
    }

    iprint0 ("\n");
}

// sort 
static SORTER (tags_sorter)
{ 
    Tag *tag_a = (Tag *)a, *tag_b = (Tag *)b;

    return tags_cmp (STRa(tag_a->tag_name), tag_a->dtype, STRa(tag_b->tag_name), tag_b->dtype);
}

// sort, but place TAG_NO_SRC and TAG_GENOZIP at the end
static SORTER (tags_sorter_demote_defaults) 
{ 
    Tag *tag_a = (Tag *)a, *tag_b = (Tag *)b;

    return (tag_a->source != tag_b->source) ? DESCENDING (Tag, source) // A
                                            : tags_cmp (STRa(tag_a->tag_name), tag_a->dtype, STRa(tag_b->tag_name), tag_b->dtype);
}

static Tag *vcf_tags_get_apriori_tag_do (DictIdType dtype, STRp(tag_name), int first, int last)
{
    if (first > last) return NULL;
    
    int mid = (first + last) / 2;
    Tag *ent = B(Tag, z_file->apriori_tags, mid);

    // first compare by dtype, and then by tag_name
    int cmp = tags_cmp (STRa(tag_name), dtype, STRa(ent->tag_name), ent->dtype);
    if (!cmp) return ent;

    if (cmp > 0) return vcf_tags_get_apriori_tag_do (dtype, STRa(tag_name), mid+1, last);
    else         return vcf_tags_get_apriori_tag_do (dtype, STRa(tag_name), first, mid-1);
}

static Tag *vcf_tags_get_apriori_tag (DictIdType dtype, STRp(tag_name))
{
    // binary search the sorted tags (param is the length of the sorted portion)
    SAFE_NUL (&tag_name[tag_name_len]);    
    Tag *tag = vcf_tags_get_apriori_tag_do (dtype, STRa(tag_name), 0, (int)z_file->apriori_tags.param-1); 
    if (tag) goto done;

    // linear search the rest
    for (uint64_t i=z_file->apriori_tags.param; i < z_file->apriori_tags.len; i++) {
        tag = B(Tag, z_file->apriori_tags, i);
        if (!tags_cmp(STRa(tag_name), dtype, STRa (tag->tag_name), tag->dtype)) goto done;
    }

    tag = NULL; // not found

done:
    SAFE_RESTORE;
    return tag;
}

// add system defaults tag renaming, unless tag specified in the command line (with any attribute)
typedef struct { STR(tag_name); DictIdType dtype; STR(refalt); STR(strand); STR(tlafer); STR(always); } DefaultTag;

static void vcf_tags_add_system_default (DefaultTag new_tag)
{
    Tag *tag = vcf_tags_get_apriori_tag_do (new_tag.dtype, STRa(new_tag.tag_name), 0, (int)z_file->apriori_tags.param-1); 
    
    if (tag) return; // command line overrides system defaults

    buf_alloc_zero (evb, &z_file->apriori_tags, 1, 100, Tag, 2, "z_file->apriori_tags"); 
    tag = &BNXT (Tag, z_file->apriori_tags);
    
    tag->dtype  = new_tag.dtype;
    tag->source = TAG_GENOZIP;
    STRcpy (tag->tag_name, new_tag.tag_name);
    STRcpyi(tag->dest, RA_REFALT, new_tag.refalt);
    STRcpyi(tag->dest, RA_STRAND, new_tag.strand);
    STRcpyi(tag->dest, RA_TLAFER, new_tag.tlafer);
    STRcpyi(tag->dest, RA_ALWAYS, new_tag.always);
}

// called from vcf_inspect_txt_header_zip: we populate its apriori tags from the command line
void vcf_tags_populate_tags_from_command_line (void)
{
    // sort command line tags, as these will be searched when adding the header tags (command line tags never change after sorted here, so we only sort in the first call to this function)
    DO_ONCE qsort (STRb(command_line_tags), sizeof (Tag), tags_sorter); 

    if (command_line_tags.len && !z_file->apriori_tags.len) {
        buf_copy (evb, &z_file->apriori_tags, &command_line_tags, Tag, 0, 0, "z_file->apriori_tags");
        z_file->apriori_tags.count = command_line_tags.len; // sorted length
    }

    // add some system-default values (the user can override these with --dvcf-rename and --dvcf-drop)
    vcf_tags_add_system_default ((DefaultTag){ cSTR("MAX_AF"),       DTYPE_VCF_INFO,   .refalt=cSTR("DROP_MAX_AF") });
    vcf_tags_add_system_default ((DefaultTag){ cSTR("DROP_MAX_AF"),  DTYPE_VCF_INFO,   .refalt=cSTR("MAX_AF") });
    vcf_tags_add_system_default ((DefaultTag){ cSTR("CLNHGVS"),      DTYPE_VCF_INFO,   .always=cSTR("DROP_CLNHGVS") });
    vcf_tags_add_system_default ((DefaultTag){ cSTR("DROP_CLNHGVS"), DTYPE_VCF_INFO,   .always=cSTR("CLNHGVS") });
    vcf_tags_add_system_default ((DefaultTag){ cSTR("ADF"),          DTYPE_VCF_FORMAT, .strand=cSTR("ADR") });
    vcf_tags_add_system_default ((DefaultTag){ cSTR("ADR"),          DTYPE_VCF_FORMAT, .strand=cSTR("ADF") });
    vcf_tags_add_system_default ((DefaultTag){ cSTR("F1R2"),         DTYPE_VCF_FORMAT, .strand=cSTR("F2R1") });
    vcf_tags_add_system_default ((DefaultTag){ cSTR("F2R1"),         DTYPE_VCF_FORMAT, .strand=cSTR("F1R2") });
    vcf_tags_add_system_default ((DefaultTag){ cSTR("REF_F1R2"),     DTYPE_VCF_FORMAT, .strand=cSTR("REF_F2R1"), .refalt = cSTR("ALT_F1R2"), .tlafer = cSTR("ALT_F2R1") });
    vcf_tags_add_system_default ((DefaultTag){ cSTR("REF_F2R1"),     DTYPE_VCF_FORMAT, .strand=cSTR("REF_F1R2"), .refalt = cSTR("ALT_F2R1"), .tlafer = cSTR("ALT_F1R2") });
    vcf_tags_add_system_default ((DefaultTag){ cSTR("ALT_F1R2"),     DTYPE_VCF_FORMAT, .strand=cSTR("ALT_F2R1"), .refalt = cSTR("REF_F1R2"), .tlafer = cSTR("REF_F2R1") });
    vcf_tags_add_system_default ((DefaultTag){ cSTR("ALT_F2R1"),     DTYPE_VCF_FORMAT, .strand=cSTR("ALT_F1R2"), .refalt = cSTR("REF_F2R1"), .tlafer = cSTR("REF_F1R2") });

    // VarScan tags
    vcf_tags_add_system_default ((DefaultTag){ cSTR("RDF"),          DTYPE_VCF_FORMAT, .strand=cSTR("RDR") }); 
    vcf_tags_add_system_default ((DefaultTag){ cSTR("RDR"),          DTYPE_VCF_FORMAT, .strand=cSTR("RDF") }); 
    vcf_tags_add_system_default ((DefaultTag){ cSTR("ABQ"),          DTYPE_VCF_FORMAT, .refalt=cSTR("RBQ") }); 
    vcf_tags_add_system_default ((DefaultTag){ cSTR("RBQ"),          DTYPE_VCF_FORMAT, .refalt=cSTR("ABQ") }); 
    // note: in VarScan, we should also have AD<>RD ADF<>RDF ADR<>RDR in case of Ref/Alt switch but
    // AD, ADF and ADR conflict the standard 'R' tags with the same need (Bug 453)

    // sort again
    qsort (STRb(z_file->apriori_tags), sizeof (Tag), tags_sorter); 
    z_file->apriori_tags.count = z_file->apriori_tags.len;
}

// iterator for getting tags that are in apriori_tags, but not in VCF header. 
Tag *vcf_tags_get_next_missing_tag (Tag *tag)
{
    for (tag = tag ? tag+1 : B1ST (Tag, z_file->apriori_tags);
         tag < BAFT (Tag, z_file->apriori_tags);
         tag++) 
        // cases we add the tag to the header: it is defined in the command line, or the symmetrical tag is in the header
        if (tag->source == TAG_CMDLINE || tag->source == TAG_CMDLINE_DST || tag->source == TAG_HEADER_DST) return tag;

    return NULL;
}

// called from vcf_inspect_txt_header_zip: sort apriori tags (command line tags + header tags)
void vcf_tags_finalize_tags_from_vcf_header (void)
{
    if (!z_file->apriori_tags.len) return; // nothing to do
 
    // sort, placing the Genozip defaults that are not in header or command line, last
    qsort (STRb(z_file->apriori_tags), sizeof (Tag), tags_sorter_demote_defaults); 

    ARRAY (Tag, tags, z_file->apriori_tags);

    // verify attribute consistency
    for (int i=tags_len-1; i >= 0; i--) {

        // case: tag is a Genozip-default tag, but doesn't appear in command line or header - we remove it
        // as we can't use it even if tag appears in data - because to cross-render, it needs the Rename attributes in the header
        if (tags[i].source == TAG_GENOZIP)
            z_file->apriori_tags.len--;

        #define HAS(x) (tags[i].dest_lens[RA_##x] > 0)

        // verify that source is assigned
        ASSINP (tags[i].source, "Tag %s/%.*s has no source", DTPT(dtype_names)[tags[i].dtype], STRf(tags[i].tag_name));

        // verify that if there is both STRAND and REFALT, there is also TLAFER
        ASSINP (!(HAS(REFALT) && HAS(STRAND)) || HAS(TLAFER), 
                "Tag %s/%.*s has both REFALT and STRAND renaming - it is required to have TLAFER renaming as well", 
                DTPT(dtype_names)[tags[i].dtype], STRf(tags[i].tag_name)) ;

        // verify that if there is TLAFER both there are also both STRAND and REFALT
        ASSINP (!HAS(TLAFER) || (HAS(REFALT) && HAS(STRAND)), 
                "Tag %s/%.*s has TLAFER - it is required to both REFALT and STRAND renaming as well", 
                DTPT(dtype_names)[tags[i].dtype], STRf(tags[i].tag_name));

        ASSINP (!HAS(ALWAYS) || (!HAS(REFALT) && !HAS(STRAND) && !HAS(TLAFER)),
                "Tag %s/%.*s has ALWAYS renaming, and therefore it cannot also have REFALT, STRAND or TLAFER renaming", 
                DTPT(dtype_names)[tags[i].dtype], STRf(tags[i].tag_name));
        #undef HAS
    }

    // sort again, without taking source into account
    qsort (STRb(z_file->apriori_tags), sizeof (Tag), tags_sorter); 
    z_file->apriori_tags.count = z_file->apriori_tags.len; // sorted length

    if (flag.show_rename_tags) vcf_tags_show_rename_tags ();
}

// called when reading VCF header: sets the attr to dest only if it doesn't yet exists, and returns the final value
// (the command line value prevails in case of conflicting values)
bool vcf_tags_add_attr_from_header (DictIdType dtype, STRp(tag_name), RenameAttr attr, 
                                    STRp (number), STRp (type), STRp (rendalg),  // Number, Types, RendAlg attributes of this tag
                                    pSTRp(dest), // in/out
                                    bool recursive)
{
    // allocate 2, one for symmetrical, so that we don't need to allocate in the recursive call (which would invalidate tag and dest)
    if (!recursive)
        buf_alloc_zero (evb, &z_file->apriori_tags, 2, 100, Tag, 2, "z_file->apriori_tags"); 

    Tag *tag = vcf_tags_get_apriori_tag (dtype, STRa(tag_name));

    // if tag doesn't exist - add one if we have a value to assign    
    if (!tag && *dest) {
        ASSINP (*dest_len < MAX_TAG_LEN, "attribute %s=%.*s appe    aring in a ##%s line with ID=%.*s is beyond Genozip's length limit of %u", 
                vcf_header_rename_attrs[attr], STRf(*dest), DTPT(dtype_names)[dtype], STRf (tag_name), MAX_TAG_LEN); // leave room for \0

        tag = &BNXT (Tag, z_file->apriori_tags);
    }

    // case: attribute not specified in command line - use attribute from VCF header
    if (*dest_len && !tag->dest_lens[attr]) {
        STRcpyi(tag->dest, attr, *dest);

        if (!tag->tag_name_len) 
            STRcpy (tag->tag_name, tag_name);
    }

    // return pre-existing attribute in tag, or newly assigned one
    bool different_than_header = false;
    if (tag) {
        different_than_header = !str_issame_(STRa(*dest), STRi(tag->dest, attr));
        
        if (!recursive) {
            *dest     = tag->dests[attr];
            *dest_len = tag->dest_lens[attr];
        }

        // initialize Number, Type and RendAlg, source upon first encounter in header
        if (tag->source != TAG_HEADER) {
            
            #define CHECK_CONSISTENCY(x,MAX_X,NAME_X) \
                ASSINP (x##_len <= MAX_X, NAME_X"=%.*s appearing in a ##%s line with ID=%.*s is beyond Genozip's length limit of %u",  \
                        STRf(x), DTPT(dtype_names)[dtype], STRf(tag_name), MAX_X); \
                \
                ASSINP (!x##_len || !tag->x##_len || str_issame (x, tag->x), \
                        "Tag %s/%.*s has "NAME_X"=%.*s, inconsistent with its Rename target which has "NAME_X"=%.*s", \
                        DTPT(dtype_names)[dtype], STRf(tag_name), STRf(x), STRf(tag->x)); \
                \
                if (!tag->x##_len && x##_len) { \
                    memcpy (tag->x, x, x##_len);\
                    tag->x##_len = x##_len; \
                }
            CHECK_CONSISTENCY (number, MAX_NUMBER_LEN, "Number");
            CHECK_CONSISTENCY (type, MAX_NUMBER_LEN, "Type");
            CHECK_CONSISTENCY (rendalg, MAX_RENDALG_LEN, HK_RENDALG);

            tag->dtype = dtype; // if it was BOTH, it is now set to INFO or FORMAT according to the header
            tag->source = !recursive ? TAG_HEADER
                        : (tag->source == TAG_NO_SRC || tag->source == TAG_GENOZIP) ? TAG_HEADER_DST
                        : tag->source;
        }
    }

    if (!recursive && *dest_len && tag_name_len) 
        vcf_tags_add_attr_from_header (dtype, STRa(*dest), attr, STRa(number), STRa (type), STRa (rendalg), pSTRa(tag_name), true);

    return different_than_header;
}

// -----------------------------------------------------------------------
// Seg stuff - operates on actual tags observed in data stored in vb->tags
// -----------------------------------------------------------------------

// called by Seg to add a tag encountered in the data
void vcf_tags_add_tag (VBlockVCFP vb, Context *ctx, DictIdType dtype, STRp(tag_name))
{
    // check if tag was already added, considering a multiple tags can be mapped to the same ctx (bc they have the same dict_id)
    if (ctx->tag_i != -1) {
        Tag *old_tag = B(Tag, vb->tags, ctx->tag_i);
        
        // case that tag pointed to from ctx, is indeed the tag we want to add - done!
        if (str_issame (old_tag->tag_name, tag_name)) return; // already added

        // since ctx->tag_i has a value, and it is not our tag - check if any of the other tags are our tag
        for (unsigned i=0; i < vb->tags.len; i++) {
            old_tag = B(Tag, vb->tags, i);
            if (str_issame (old_tag->tag_name, tag_name)) return; // already added
        }
    }

    buf_alloc_zero (vb, &vb->tags, 1, 100, Tag, 2, "tags");

    ctx->tag_i = vb->tags.len;

    Tag *tag = &BNXT (Tag, vb->tags);

    // case: we address this tag in the command line (--chain) or VCF header attributes (DVCF files)
    Tag *apriori_tag = vcf_tags_get_apriori_tag (dtype, STRa(tag_name));
    
    if (apriori_tag)
        *tag = *apriori_tag; // an apriori tag (from --dvcf-rename / --dvcf-drop / VCF header) takes precedence over the Genozip defaults

    else {
        tag->dtype = dtype; // note: it is possible that ctx->dtype than the type implied by ctx->dict_id bc the context might cover several tags with the same dict_id
        STRcpy (tag->tag_name, tag_name); 
    }
}

static const Tag *vcf_tags_rename_get_tag (VBlockVCFP vb, ConstContextP ctx, STRp(tag_name))
{
    const Tag *tag = B(const Tag, vb->tags, ctx->tag_i);

    // tag[i], pointed to from the context is usually the right tag, but sometimes it might be the wrong tag
    // if two or more tags have the same dict_id and hence the same context. in this case, we search for the tag
    if (!str_issame (tag_name, tag->tag_name)) {
        tag = NULL;
        for (unsigned t=0; t < vb->tags.len; t++) {
            const Tag *candidate = B(Tag, vb->tags, t);
            if (str_issame (candidate->tag_name, tag_name)) {
                tag = candidate;
                break;
            }
        }
        ASSVCF (tag, "cannot find tag=\"%.*s\" in vb->tags(len=%"PRIu64")", STRf(tag_name), vb->tags.len);
    }

    return tag;
}

// checks if any of the tags in the container have renaming, given the ostatus and xstrand of the current line,
// if yes - renamed string is in "renamed" and length is returned; if no - returns 0
unsigned vcf_tags_rename (VBlockVCFP vb, 
                          unsigned num_tags,
                          const ContextPBlock ctxs, rom sf_names[], const unsigned sf_name_lens[], // option 1 - FORMAT 
                          const InfoItem *ii,   // option 2 - INFO - an array of InfoItem
                          char *renamed)  // out - allocated by caller
{
    ASSERT (vb->tags.len >= num_tags, "Expecting vb->tags.len=%"PRIu64" >= num_tags=%u", vb->tags.len, num_tags);

    #define HAS(x) (tags[i]->dest_lens[RA_##x] > 0)

    bool is_xstrand = vb->last_index (VCF_oXSTRAND);
    bool is_switch = LO_IS_OK_SWITCH (last_ostatus);
    bool is_format = !ii;

    const Tag *tags[num_tags];

    // check first whether renaming is needed at all for this variant before expensive string-copying
    bool needs_renaming=false;
    for (unsigned i=0; i < num_tags; i++) {
        tags[i] = vcf_tags_rename_get_tag (vb, is_format ? ctxs[i]         : ii[i].ctx,
                                               is_format ? sf_names[i]     : ii[i].name,
                                               is_format ? sf_name_lens[i] : ii[i].name_len - !!ii[i].value); // without '=' for valueful items
        
        needs_renaming |=                             HAS (ALWAYS)
                       || (is_switch               && HAS (REFALT))
                       || (is_xstrand              && HAS (STRAND))
                       || (is_switch && is_xstrand && HAS (TLAFER));
    }

    if (!needs_renaming) return 0;

    // created renamed string
    char *p = renamed;

    if (ii) {
        p[0] = p[1] = CON_PX_SEP;
        p += 2;
    }

    for (unsigned i=0; i < num_tags; i++) {

        #define COPY_TAG(x) do { \
            memcpy (p, tags[i]->dests[RA_##x], tags[i]->dest_lens[RA_##x]); \
            p += tags[i]->dest_lens[RA_##x]; \
        } while (0)

        // concatenated (possibly renamed) tag
        if      (                           HAS (ALWAYS)) COPY_TAG (ALWAYS);
        else if (is_switch && is_xstrand && HAS (TLAFER)) COPY_TAG (TLAFER);
        else if (is_switch               && HAS (REFALT)) COPY_TAG (REFALT);
        else if (is_xstrand              && HAS (STRAND)) COPY_TAG (STRAND);
        else { memcpy (p, tags[i]->tag_name, tags[i]->tag_name_len); p += tags[i]->tag_name_len; }                                             
            
        // concatenate separator
        if (is_format) {
            if (i < num_tags-1) *(p++) = ':';
        }
        
        else {
            if (ii[i].value) *(p++) = '=';
            *(p++) = CON_PX_SEP;
        }
    }

    return p - renamed; // renamed length
    #undef HAS
}

