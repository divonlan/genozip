// ------------------------------------------------------------------
//   header_vcf.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "zfile.h"
#include "txtfile.h"
#include "crypt.h"
#include "dispatcher.h"
#include "txtfile.h"
#include "txtheader.h"
#include "contigs.h"
#include "stats.h"
#include "tip.h"
#include "arch.h"

// VCF standard keys
#define HK_GENOZIP_CMD   "##genozip_command="
#define TAG_SOURCE       Source=\""GENOZIP_URL"\""
#define HK_INFO          "##INFO="
#define HK_FORMAT        "##FORMAT="
#define HK_CONTIG        "##contig="

// Globals
static VcfVersion vcf_version;
uint32_t vcf_num_samples = 0; // number of samples in the file
static uint32_t vcf_num_hdr_contigs = 0;
static uint64_t vcf_num_hdr_nbases = 0;
static uint32_t vcf_num_displayed_samples = 0; // PIZ only: number of samples to be displayed - might be less that vcf_num_samples if --samples is used
static Buffer vcf_field_name_line = {};  // header line of first VCF file read - use to compare to subsequent files to make sure they have the same header during bound

static bool vcf_has_file_format_vcf = false; 
bool vcf_header_get_has_fileformat (void) { return vcf_has_file_format_vcf; }

#define FI_LEN (is_info ? 7 : 9) // length of ##FORMAT= or ##INFO=

// referring to sample strings from the --samples command line option
static Buffer cmd_samples_buf = {}; // an array of (char *)
static bool cmd_is_negative_samples = false;

// referring to samples in the vcf file
char *vcf_samples_is_included;         // a bytemap indicating for each sample if it is included
static char **vcf_sample_names;        // an array of char * to nul-terminated names of samples 
static char *vcf_sample_names_data;    // vcf_sample_names point into here

#define LINEIS(s) (line_len > (sizeof s - 1) && !memcmp (line, (s), sizeof s - 1)) 

#define SUBST_LABEL(old,new) ({ buf_add_more (NULL, new_txt_header, (new), (sizeof new - 1), NULL); \
                                buf_add_more (NULL, new_txt_header, line + (sizeof old - 1), line_len - (sizeof old - 1), NULL); })

typedef struct { STR (key); STR (value); } Attr;

int32_t vcf_is_header_done (bool is_eof)
{
    // verify that this isn't by any chance a BCF...
    ASSINP0 (!is_bcf (STRb(evb->txt_data), NULL), "This is BCF file, but Genozip got confused because the file's name implies it is a VCF. Solution: re-run and add the command line option \"--input bcf\"");

    return def_is_header_done (is_eof);
}

void vcf_piz_header_init (void)
{
    vcf_num_samples = vcf_num_hdr_contigs = 0;
    vcf_num_hdr_nbases = 0;
    buf_free (vcf_field_name_line);

    // note: we don't re-initialize vcf_num_displayed_samples - this is calculated only once
}

void vcf_header_finalize (void)
{
    vcf_num_samples = vcf_num_hdr_contigs = 0;
    vcf_num_hdr_nbases = 0;

    buf_destroy (vcf_field_name_line);
}

static void vcf_header_subset_samples (BufferP vcf_header);
static bool vcf_header_set_globals (rom filename, BufferP vcf_header, FailType soft_fail);
static void vcf_header_trim_field_name_line (BufferP vcf_header);

// returns the length of the first line, if it starts with ##fileformat, or 0
// note: per VCF spec, the ##fileformat line must be the first in the file
static inline void vcf_header_set_vcf_version (ConstBufferP txt_header)
{
    ARRAY (char, line, *txt_header);
    if (!LINEIS ("##fileformat")) return;

    rom newline = memchr (line, '\n', line_len);
    unsigned len = newline ? newline - line + 1 : 0;

    if (len >= 20 && str_isprefix_(line, len, _S("##fileformat=VCFv4."))) 
        switch (line[19]) {
            case '1' : vcf_version = VCF_v4_1; break;
            case '2' : vcf_version = VCF_v4_2; break;
            case '3' : vcf_version = VCF_v4_3; break;
            case '4' : vcf_version = VCF_v4_4; break;
            case '5' : vcf_version = VCF_v4_5; break;
            default: {}
        }             
}

static inline bool is_field_name_line (STRp(line))
{
    // compare to "#CHROM" - not to entire VCF_FIELD_NAMES as it may be a line maybe separated by spaces instead of tabs
    return LINEIS ("#CHROM");
}

// returns the length of the last line, if it starts with #CHROM, or 0
static bool vcf_header_get_last_line_cb (STRp(line), void *unused1, void *unused2, unsigned unused3)
{ return true; }
static unsigned vcf_header_get_last_line (BufferP txt_header, char **line_p)
{
    int64_t line_len;
    char *line = txtfile_foreach_line (txt_header, true, vcf_header_get_last_line_cb, 0, 0, 0, &line_len);
    if (line_p) *line_p = line;

    return is_field_name_line (STRa(line)) ? (unsigned)line_len : 0;
}

static void vcf_header_add_genozip_command (VBlockP txt_header_vb, BufferP txt_header)
{
    // the command line length is unbound, careful not to put it in a bufprintf
    buf_append_string (txt_header_vb, txt_header, HK_GENOZIP_CMD"\"");
    buf_append_string (txt_header_vb, txt_header, flags_command_line());
    bufprintf (txt_header_vb, txt_header, "\" %s\n", str_time().s);
}

static void vcf_header_get_attribute (STRp(line), unsigned key_len, STRp(attr), bool remove_quotes, bool enforce,// in
                                      pSTRp(snip)) // out
{
    SAFE_NUL (&line[line_len]);
    rom start = strstr (line + key_len, attr);
    SAFE_RESTORE;

    if (!start) goto missing; // attr doesn't exist or malformed line

    start += attr_len;

    rom comma   = memchr (start, ',', &line[line_len] - start);
    rom bracket = memchr (start, '>', &line[line_len] - start);

    if (!comma && !bracket) goto missing; // attr isn't terminated with a comma or >
    rom after = !comma ? bracket : !bracket ? comma : MIN_(comma,bracket);
    
    *snip_len = after - start ;
    *snip = *snip_len ? start : NULL;

    if (remove_quotes && *snip) {
        ASSINP ((*snip)[0]=='"' && (*snip)[*snip_len-1]=='"', "Expecting value an attribute %.*s to be enclosed in quotes. Line=\"%.*s\"", STRf(attr), STRf(line));
        (*snip)++;
        (*snip_len) -= 2;
    }
    return;

missing:
    *snip = NULL;
    *snip_len = 0;
    ASSINP (!enforce, "missing attr %s in header line %.*s", attr, line_len, line);
}

// ZIP: VCF main component (rejects components txt headers don't have contigs)
static void vcf_header_consume_contig (STRp (contig_name), PosType64 *LN)
{
    WordIndex ref_index = WORD_INDEX_NONE;

    // case: we have a reference - verify length
    if (flag.reference & REF_ZIP_LOADED)
        ref_index = ref_contigs_ref_chrom_from_header_chrom (gref, STRa(contig_name), LN);

    // also populates chrom2ref_map if external reference is loaded
    ctx_populate_zf_ctx (VCF_CHROM, STRa(contig_name), ref_index); 

    if (flag.show_txt_contigs) 
        iprintf ("\"%.*s\" LN=%"PRId64" ref_index=%d\n", STRf(contig_name), *LN, ref_index);
}

// in case of --optimize: add PL and PP (if converting from GL or GP) - these might be dup lines, that's ok
// bug 1079
// static void vcf_header_add_FORMAT_lines (BufferP txt_header)
// {
//     txt_header->len -= vcf_field_name_line.len; // remove field name line

//     #define FORMAT_PP "##FORMAT=<ID=PP,Number=G,Type=Integer,Description=\"Phred-scaled genotype posterior probabilities rounded to the closest integer\">\n"
//     buf_add_moreC (evb, txt_header, FORMAT_PP, "txt_data");
//     z_file->header_size += STRLEN(FORMAT_PP); // header_size has grown

//     #define FORMAT_PL "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods rounded to the closest integer\">\n"
//     buf_add_moreC (evb, txt_header, FORMAT_PL, "txt_data");
//     z_file->header_size += STRLEN(FORMAT_PL); // header_size has grown

//     // add back the field name (#CHROM) line
//     buf_add_buf (evb, txt_header, &vcf_field_name_line, char, "txt_data");
// }

static void vcf_header_zip_create_contig (STRp(line))
{
    // parse eg "##contig=<ID=NKLS02001838.1,length=29167>" 
    SAFE_NUL (&line[line_len]);

    STR(contig_name);
    vcf_header_get_attribute (STRa(line), STRLEN (HK_CONTIG), cSTR("ID="),     false, true,  pSTRa(contig_name));

    STR(length_str);
    vcf_header_get_attribute (STRa(line), STRLEN (HK_CONTIG), cSTR("length="), false, false, pSTRa(length_str));
    PosType64 LN = 0;
    str_get_int_range64 (STRa(length_str), 1, 1000000000000ULL, &LN); // length stays 0 if length_str_len=0

    vcf_header_consume_contig (STRa (contig_name), &LN); // also verifies / updates length
    vcf_num_hdr_contigs++;
    vcf_num_hdr_nbases += LN;

    SAFE_RESTORE;
}

static void vcf_header_zip_create_ctx (STRp(line), bool is_info)
{
    SAFE_NUL (&line[line_len]);
    int key_len = is_info ? STRLEN (HK_INFO) : STRLEN (HK_FORMAT);

    STR(tag_name);
    vcf_header_get_attribute (STRa(line), key_len, cSTR("ID="), false, false, pSTRa(tag_name));
    
    STR(number); // expecting: an integer string, "A", "R", "G" or "." 
    vcf_header_get_attribute (STRa(line), key_len, cSTR("Number="), false, false, pSTRa(number));
    
    STR(type); // expecting:  "Integer", "Float", "Flag", "Character", or "String"
    vcf_header_get_attribute (STRa(line), key_len, cSTR("Type="), false, false, pSTRa(type));

    if (!tag_name || !number || !type)
        goto done; // ignore this header line if one of the needed attributes is missing;

    DictId dict_id = dict_id_make (STRa(tag_name), is_info ? DTYPE_1 : DTYPE_2);

    ContextP zctx = ctx_add_new_zf_ctx_at_init (STRa(tag_name), dict_id);
    if (!zctx) zctx = ctx_get_existing_zctx (dict_id); // pre-defined context

    int64_t number_int;    
    if      (str_issame_(STRa(number), "A", 1)) zctx->header_info.vcf.Number = NUMBER_A;
    else if (str_issame_(STRa(number), "G", 1)) zctx->header_info.vcf.Number = NUMBER_G;
    else if (str_issame_(STRa(number), "R", 1)) zctx->header_info.vcf.Number = NUMBER_R;
    else if (str_issame_(STRa(number), ".", 1)) zctx->header_info.vcf.Number = NUMBER_VAR;
    else if (str_get_int (STRa(number), &number_int)) 
        zctx->header_info.vcf.Number = (number_int >= 1 && number_int <= 127) ? number_int : NUMBER_LARGE; 
    else goto done; // unrecgonized Number - remains NOT_IN_HEADER

    if      (str_issame_(STRa(type), _S("Integer")))   zctx->header_info.vcf.Type = VCF_Integer;
    else if (str_issame_(STRa(type), _S("Float")))     zctx->header_info.vcf.Type = VCF_Float;
    else if (str_issame_(STRa(type), _S("Flag")))      zctx->header_info.vcf.Type = VCF_Flag;
    else if (str_issame_(STRa(type), _S("Character"))) zctx->header_info.vcf.Type = VCF_Character;
    else if (str_issame_(STRa(type), _S("String")))    zctx->header_info.vcf.Type = VCF_String;
    else zctx->header_info.vcf.Number = NOT_IN_HEADER; // unrecognized Type

    if (zctx->did_i == FORMAT_GP) {
        if (strstr (line, "Phred") || strstr (line, "phred")) 
            segconf.FMT_GP_content = GP_phred;
        else if (strstr (line, "Estimated Posterior Probabilities"))
            segconf.FMT_GP_content = GP_probabilities;
    }

done:
    SAFE_RESTORE;
}

static bool vcf_header_zip_parse_one_line (STRp(line), void *unused1, void *unused2, unsigned unused3)
{
    bool is_info=false;

    // populate contigs
    if (LINEIS (HK_CONTIG)) 
        vcf_header_zip_create_contig (STRa(line));

    // create contexts that appear in header 
    else if ((is_info = LINEIS(HK_INFO)) || LINEIS(HK_FORMAT)) 
        vcf_header_zip_create_ctx (STRa(line), is_info);

    // keywords that designate programs
    #define ISPROG(s) \
    else if (LINEIS (s)) stats_add_one_program (&line[STRLEN(s)], line_len-STRLEN(s)-1) // without the \n

    // identify programs
    ISPROG("##source=");
    ISPROG("##imputation=");
    ISPROG("##phasing=");
    ISPROG("##annotator=");

    return false; // continue iterating
}

static bool vcf_inspect_txt_header_zip (BufferP txt_header)
{
    if (!vcf_header_set_globals (txt_file->name, txt_header, true)) return false; // samples are different than a previous concatented file

    // populate contigs, identify programs, identify float fields for optimization    
    txtfile_foreach_line (txt_header, false, vcf_header_zip_parse_one_line, 0, 0, 0, 0);
    
    SAFE_NULB (*txt_header);
    #define IF_IN_SOURCE(signature, segcf) if (stats_is_in_programs (signature)) segconf.segcf = true
    #define IF_IN_HEADER(signature, segcf, program) ({ if ((sig = strstr (txt_header->data, (signature)))) {                    \
                                                           segconf.segcf = true;                                                \
                                                           if (program[0]) stats_add_one_program (program, strlen (program)); } })
    
    // when adding here, also add to stats_output_file_metadata()
    rom sig;
    IF_IN_SOURCE ("infiniumFinalReportConverter", vcf_is_infinium);
    IF_IN_SOURCE ("VarScan", vcf_is_varscan);
    IF_IN_SOURCE ("dbSNP", vcf_is_dbSNP);
    IF_IN_SOURCE ("COSMIC", vcf_is_cosmic);
    IF_IN_SOURCE ("ClinVar", vcf_is_clinvar);
    IF_IN_SOURCE ("IsaacVariantCaller", vcf_is_isaac);
    IF_IN_SOURCE ("starling", vcf_is_isaac);
    IF_IN_SOURCE ("Platypus", vcf_is_platypus); // https://github.com/andyrimmer/Platypus
    IF_IN_SOURCE ("GenerateSVCandidates", vcf_is_manta); // https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md
    IF_IN_SOURCE ("svaba", vcf_is_svaba);
    IF_IN_SOURCE ("pbsv", vcf_is_pbsv);
    IF_IN_SOURCE ("callMom", vcf_is_callmom);
    IF_IN_SOURCE ("MELT", vcf_is_melt);
    IF_IN_SOURCE ("GLIMPSE_phase", vcf_is_GLIMPSE_phase);
    IF_IN_SOURCE ("Gencove", vcf_is_gencove);
    IF_IN_SOURCE ("freeBayes", vcf_is_freebayes);
    IF_IN_HEADER ("GenotypeGVCFs", vcf_is_gatk_gvcf, "GenotypeGVCFs");
    IF_IN_HEADER ("CombineGVCFs", vcf_is_gatk_gvcf, "CombineGVCFs");
    if (segconf.vcf_is_gatk_gvcf) segconf.vcf_is_gvcf = true;
    if (segconf.vcf_is_isaac) IF_IN_HEADER ("gvcf", vcf_is_gvcf, "");
    IF_IN_HEADER ("beagle", vcf_is_beagle, "beagle");
    IF_IN_HEADER ("Pindel", vcf_is_pindel, "Pindel");
    IF_IN_HEADER ("caveman", vcf_is_caveman, "CaVEMan");
    IF_IN_HEADER ("gnomAD", vcf_is_gnomad, "gnomAD");
    IF_IN_HEADER ("ExAC", vcf_is_exac, ""); // "" = not a program
    IF_IN_HEADER ("Mastermind", vcf_is_mastermind, "Mastermind");
    IF_IN_HEADER ("VcfAnnotate.pl", vcf_is_vagrent, "VAGrENT");
    IF_IN_HEADER ("ICGC", vcf_is_icgc, "ICGC");
    IF_IN_HEADER ("Illumina GenCall", vcf_illum_gtyping, "Illumina GenCall");
    IF_IN_HEADER ("Log R Ratio", vcf_illum_gtyping, "Illumina Genotyping");
    IF_IN_HEADER ("Number of cases used to estimate genetic effect", vcf_is_gwas, "GWAS_1.0"); // v1.0
    IF_IN_HEADER ("##trait", vcf_is_gwas, "GWAS_1.2"); /*v1.2*/
    IF_IN_HEADER ("DRAGEN", vcf_is_dragen, "DRAGEN");
    IF_IN_HEADER ("DeepVariant", vcf_is_deep_variant, "DeepVariant");
    IF_IN_HEADER ("##hailversion", vcf_is_hail, "Hail");
    IF_IN_HEADER ("Giggle", vcf_is_giggle, "Giggle"); // appears in description of fields such as GT and GQ
    
    if (segconf.vcf_is_svaba) 
        stats_remove_data_after_program_name ("svaba"); // remove the command line that appears after "svaba"

    segconf.vcf_is_sv = segconf.vcf_is_svaba || segconf.vcf_is_manta || segconf.vcf_is_pbsv;

    #define VEP_SIGNATURE "VEP. Format: "
    IF_IN_HEADER (VEP_SIGNATURE, vcf_is_vep, "VEP");
    if (segconf.vcf_is_vep) {
        rom vep = sig + STRLEN(VEP_SIGNATURE);
        rom quote = strpbrk (vep, "\"\n");
        if (quote && *quote == '"') { // cursory verification of field format
            unsigned vep_len = quote - vep;
            segconf.vcf_vep_spec = CALLOC (vep_len + 1); 
            memcpy (segconf.vcf_vep_spec, vep, vep_len); // consumed and freed by by vcf_vep_zip_initialize
        }
    }
    
    bool has_PROBE = !!strstr (txt_header->data, "##INFO=<ID=PROBE_A");
    SAFE_RESTORE;

    if (!flag.reference && segconf.vcf_illum_gtyping && !flag.seg_only && has_PROBE)
        TIP ("Compressing an Illumina Genotyping VCF file using a reference file can reduce the compressed file's size by 20%%.\n"
             "Use: \"%s --reference <ref-file> %s\". ref-file may be a FASTA file or a .ref.genozip file.\n",
             arch_get_argv0(), txt_file->name);
    
    // (bug 1079) add PP and/or PL INFO lines if needed
    // if (flag.optimize)
    //     vcf_header_add_FORMAT_lines (txt_header);

    // set vcf_version
    vcf_header_set_vcf_version (txt_header);

    return true; // all good
}

static bool vcf_inspect_txt_header_piz (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags)
{
    vcf_header_set_globals (z_file->name, txt_header, false);

    if (flag.genocat_no_reconstruct) return true;

    // remove #CHROM line (it is saved in vcf_field_name_line by vcf_header_set_globals()) - or everything
    // if we ultimately want the #CHROM line. We will add back the #CHROM line later.
    if (flag.header_one) 
        txt_header->len = 0; 
    else
        txt_header->len -= vcf_field_name_line.len; 

    // if genocat --samples is used, update vcf_header and vcf_num_displayed_samples
    // note: we subset reject samples (in the header) as well - we analyze only in the rejects section
    if (flag.samples) vcf_header_subset_samples (&vcf_field_name_line);
    else              vcf_num_displayed_samples = vcf_num_samples;

    // add genozip command line
    if (!flag.header_one && is_genocat && !flag.genocat_no_reconstruct && !flag.no_pg && flag.piz_txt_modified) 
        vcf_header_add_genozip_command (txt_header_vb, txt_header);

    if (flag.drop_genotypes) 
        vcf_header_trim_field_name_line (&vcf_field_name_line); // drop FORMAT and sample names

    // add the (perhaps modified) field name (#CHROM) line
    buf_add_more (txt_header_vb, txt_header, vcf_field_name_line.data, vcf_field_name_line.len, "txt_data");

    return true; // all good
}

bool vcf_inspect_txt_header (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags)
{
    return (IS_ZIP) ? vcf_inspect_txt_header_zip (txt_header)
                            : vcf_inspect_txt_header_piz (txt_header_vb, txt_header, txt_header_flags);
}

static bool vcf_header_set_globals (rom filename, BufferP vcf_header, FailType soft_fail)
{
    // check for ##fileformat=VCF
    vcf_has_file_format_vcf = str_isprefix_(STRb(*vcf_header), _S("##fileformat=VCF"));

    // count tabs in last line which should be the field header line
    unsigned tab_count = 0;
    for (int i=vcf_header->len-1; i >= 0; i--) {
        
        #define THIS (*Bc(*vcf_header, i))
        #define PREV ((i>=1) ? *Bc(*vcf_header, i-1) : 0)

        if (THIS == '\t')
            tab_count++;
        
        // some times files have spaces instead of \t 
        else if (THIS == ' ') {
            tab_count++;
            while (PREV==' ') i--; // skip run of spaces
        }

        // if this is the beginning of field header line 
        else if (THIS == '#' && (i==0 || PREV == '\n' || PREV == '\r')) {

            ASSINP0 (*Bc(*vcf_header, i+1) != '#', "Error: Missing VCF field/samples header line");  

            // note: we don't memcmp to the entire VCF_FIELD_NAMES, bc some files have spaces instead of tabs
            ASSINP (vcf_header->len - i > STRLEN(VCF_FIELD_NAMES) && !memcmp ("#CHROM", Bc(*vcf_header, i), 6), 
                    "Error: Invalid VCF field/samples header line, found (partial): \"%.*s\"", (int)MIN_(vcf_header->len - i, STRLEN(VCF_FIELD_NAMES)), Bc(*vcf_header, i));

            // ZIP: if first vcf file ; PIZ: everytime - copy the header to the global
            if (!buf_is_alloc (&vcf_field_name_line) || IS_PIZ) 
                buf_copy (evb, &vcf_field_name_line, vcf_header, char, i, vcf_header->len - i, "vcf_field_name_line");

            // count samples
            vcf_num_samples = (tab_count >= 9) ? tab_count-8 : 0; 
            
            // note: a VCF file without samples may or may not have a "FORMAT" in the header, i.e. tab_count==7 or 8 (8 or 9 fields).
            // however, even if it has a FORMAT in the header, it won't have a FORMAT column in the data

            ASSINP (tab_count >= 7, "Error: Invalid VCF field/samples header line - it contains only %d fields, expecting at least 8", tab_count+1);

            return true; 
        }
    }

    ABORT ("Error: invalid VCF file - it does not contain a field header line; tab_count=%u", tab_count+1);
}

// genocat: remove FORMAT and sample names from the vcf header line, in case of --drop-genotypes
// TODO - remove FORMAT lines from header itself
static void vcf_header_trim_field_name_line (BufferP vcf_header)
{
    char *line;
    unsigned line_len = vcf_header_get_last_line (vcf_header, &line);

    if (is_field_name_line (STRa(line))) {
        vcf_header->len = strstr (line, "INFO") + 5 - B1STc (*vcf_header);
        *BLSTc (*vcf_header) = '\n';
    }                  
}

uint32_t vcf_header_get_num_samples (void)
{
    if Z_DT(VCF)
        return vcf_num_samples;
    else
        return 0;
}

uint32_t vcf_header_get_num_contigs (void)
{
    if Z_DT(VCF)
        return vcf_num_hdr_contigs;
    else
        return 0;
}

uint64_t vcf_header_get_nbases (void)
{
    if Z_DT(VCF)
        return vcf_num_hdr_nbases;
    else
        return 0;
}

// -------------
// samples stuff
// -------------

// processes the vcf header sample line according to the --samples option, removing samples that are not required
// and building a bytemap PIZ filter the samples
static void vcf_header_subset_samples (BufferP vcf_field_name_line)
{
    // accept a sample from the vcf file's samples as consistent with the --samples requested
    #define samples_accept(sample_str) { \
        buf_append_string (evb, vcf_field_name_line, sample_str); \
        buf_append_string (evb, vcf_field_name_line, "\t"); \
        vcf_num_displayed_samples++; \
    }

    ARRAY (char, line, *vcf_field_name_line);

    flag.samples = is_field_name_line (STRa(line));
    RETURNW (flag.samples,, "Warning: found non-standard VCF sample header line. Ingoring --samples : \n%.*s", (int)line_len, line);

    int32_t num_samples=-8;
    for (unsigned i=0; i < line_len; i++)
        if (line[i] == '\t') num_samples++;
        
    ASSINP (num_samples >= 1, "Cannot use --samples with %s, as this file has no samples", z_name);
    
    vcf_samples_is_included = MALLOC (num_samples);
    memset (vcf_samples_is_included, cmd_is_negative_samples, num_samples); // 0 if not included unless list says so (positive) and vice versa

    unsigned vcf_names_start_index = (line + sizeof VCF_FIELD_NAMES_LONG-1 + 1) - B1STc (*vcf_field_name_line);
    unsigned vcf_names_data_len = vcf_field_name_line->len - vcf_names_start_index;
    vcf_sample_names_data = MALLOC (vcf_names_data_len);
    memcpy (vcf_sample_names_data, &vcf_field_name_line->data[vcf_names_start_index], vcf_names_data_len);
    vcf_sample_names_data[vcf_names_data_len-1] = '\t'; // change last separator from \n to \t

    vcf_sample_names = MALLOC (num_samples * sizeof (char *));

    vcf_field_name_line->len = vcf_names_start_index;
    char *next_token = vcf_sample_names_data;
    
    // go through the vcf file's samples and add those that are consistent with the --samples requested
    vcf_num_displayed_samples = 0;

    buf_copy (evb, &evb->codec_bufs[0], &cmd_samples_buf, char*, 0, 0, 0);
    ARRAY (rom , snames, evb->codec_bufs[0]);

    int fixed_len = (snames_len == 1 && atoi (snames[0]) >= 1) ? atoi (snames[0]) : 0;

    for (unsigned i=0; i < num_samples; i++) {
        vcf_sample_names[i] = strtok_r (next_token, "\t", &next_token);
        bool handled = false;

        // case: --samples <number>
        if (fixed_len && i < fixed_len) {
            vcf_samples_is_included[i] = true;
            samples_accept (vcf_sample_names[i]);
            (*(uint64_t *)&snames_len) = 0;
            continue;
        }
        
        for (unsigned s=0; s < snames_len; s++) 
            if (!strcmp (vcf_sample_names[i], snames[s])) { // found
                vcf_samples_is_included[i] = !cmd_is_negative_samples;
                if (!cmd_is_negative_samples) 
                    samples_accept (vcf_sample_names[i]);

                // remove this sample from the --samples list as we've found it already
                memcpy (&snames[s], &snames[s+1], (snames_len-s-1) * sizeof (char *));
                (*(uint64_t *)&snames_len)--; // override const
                handled = true;
                break;
            }
            
        if (!handled && cmd_is_negative_samples) 
            samples_accept (vcf_sample_names[i]);
    }

    *BLSTc (*vcf_field_name_line) = '\n'; // change last separator from \t

    // warn about any --samples items that were not found in the vcf file (all the ones that still remain in the buffer)
    for (unsigned s=0; s < snames_len; s++) 
        ASSERTW (false, "Warning: requested sample '%s' is not found in the VCF file, ignoring it", snames[s]);

    // if the user filtered out all samples, its equivalent of drop_genotypes
    if (!vcf_num_displayed_samples) flag.drop_genotypes = true;

    buf_free (evb->codec_bufs[0]);
}

// called from genozip.c for processing the --samples flag
void vcf_samples_add  (rom samples_str)
{
    ASSERTNOTNULL (samples_str);

    // --samples 0 is interpreted as --drop-genotypes
    if (samples_str[0]=='0' && !samples_str[1] && !cmd_samples_buf.len) {
        flag.drop_genotypes = true;
        flag.samples = false;
    }

    bool is_negated = samples_str[0] == '^';

    bool is_conflicting_negation = (cmd_samples_buf.len && (cmd_is_negative_samples != is_negated));
    ASSINP0 (!is_conflicting_negation, "Error: inconsistent negation - all samples listed must either be negated or not");

    cmd_is_negative_samples = is_negated;

    // make a copy of the string and leave the original one for error message. 
    // we don't free this memory as chrom fields in regions will be pointing to it
    char *next_region_token = MALLOC (strlen (samples_str)+1); // heap memory, as cmd_samples_buf elements point into this
    strcpy (next_region_token, samples_str + is_negated); // drop the ^ if there is one

    while (1) {
        char *one_sample = strtok_r (next_region_token, ",", &next_region_token);
        if (!one_sample) break;

        bool is_duplicate = false;
        for (unsigned s=0; s < cmd_samples_buf.len32; s++)
            if (!strcmp (one_sample, *B(char *, cmd_samples_buf, s))) {
                is_duplicate = true;
                break;
            }
        if (is_duplicate) continue; // skip duplicates "genocat -s sample1,sample2,sample1"

        buf_alloc (evb, &cmd_samples_buf, 1, 100, char*, 2, "cmd_samples_buf");

        BNXT (char *, cmd_samples_buf) = one_sample;
    }
}

VcfVersion vcf_header_get_version (void) { return vcf_version; }