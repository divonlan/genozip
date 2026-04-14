// ------------------------------------------------------------------
//   sections_show.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "sections.h"
#include "vblock.h"
#include "file.h"
#include "mgzip.h"
#include "codec.h"
#include "license.h"
#include "crypt.h"
#include "zfile.h"

extern const struct {rom name; uint32_t header_size; } abouts[NUM_SEC_TYPES];

typedef struct SectionEnt SectionEntModifiable;

rom st_name (SectionType sec_type)
{
    static char invalid[24]; // not thread safe, but not expected except in error situations
    if (sec_type < SEC_NONE || sec_type >= NUM_SEC_TYPES) {
        snprintf (invalid, sizeof (invalid), "INVALID(%d)", sec_type);
        return invalid;
    }

    return (sec_type == SEC_NONE) ? "SEC_NONE" : type_name (sec_type, &abouts[sec_type].name , ARRAY_LEN(abouts));
}

StrText comp_name_(CompIType comp_i)
{
    static int max_comps_by_dt[NUM_DATATYPES] = { [DT_SAM]=3/*except if --deep*/, [DT_BAM]=3, [DT_FASTQ]=2 };
    
    static rom comp_names[NUM_DATATYPES][5]   = { [DT_SAM] = SAM_COMP_NAMES, [DT_BAM] = SAM_COMP_NAMES,
                                                  [DT_FASTQ] = FASTQ_COMP_NAMES };
    
    DataType dt = z_file ? z_file->data_type : DT_NONE;
    StrText s;

    if (!z_file) // can happen if another thread is busy exiting the process
        strcpy (s.s, "EMEM");
    
    else if (!max_comps_by_dt[dt] && comp_i==0)
        strcpy (s.s, "MAIN");

    else if (comp_i==COMP_NONE)
        strcpy (s.s, "NONE");

    else if ((dt==DT_SAM || dt==DT_BAM) && comp_i >= SAM_COMP_FQ00) 
        snprintf (s.s, sizeof (s.s), "FQ%02u", comp_i - SAM_COMP_FQ00); 

    else if (comp_i >= max_comps_by_dt[dt]) 
        snprintf (s.s, sizeof (s.s), "inv-comp-%d", comp_i);

    else 
        strcpy (s.s, comp_names[dt][comp_i]);

    return s;
}

static StrText comp_name_ex (CompIType comp_i, SectionType st)
{
    if (ST(RECON_PLAN) || ST(TXT_HEADER)) {
        StrText s;
        snprintf (s.s, sizeof (s.s), "%s%02u",ST(RECON_PLAN) ? "RP" : "TX", comp_i);
        return s;
    }
        
    else {
        if (IS_COMP_SEC(st)) return comp_name_(comp_i);
        if (ST(REFERENCE))   return (StrText){ "REFR" };
        if (ST(REF_IS_SET))  return (StrText){ "REFR" };
        if (ST(REF_HASH))    return (StrText){ "REFH" };
        if (ST(DICT))        return (StrText){ "DICT" };
        else                 return (StrText){ "----" };
    }
}

StrText vb_name (VBlockP vb)
{
    StrText s;
    if (vb && vb->pool == POOL_BGZF)
        snprintf (s.s, sizeof (s.s), "BGZF/%u", vb->vblock_i);
    else if (vb && vb->vblock_i)
        snprintf (s.s, sizeof (s.s), "%.10s/%u", comp_name (vb->comp_i), vb->vblock_i);
    else if (vb && segconf_running)
        strcpy (s.s, "SEGCONF");
    else
        strcpy (s.s, "MAIN");

    return s;
}

StrText line_name (VBlockP vb)
{
    StrText s;
    snprintf (s.s, sizeof (s.s), "%.10s/%u%s", vb_name(vb).s, vb->line_i, vb->preprocessing ? "(preproc)" : "");
    return s;
}

rom lt_name (LocalType lt)
{
    if (lt >= 0 && lt < NUM_LOCAL_TYPES) 
        return lt_desc[lt].name;
    else
        return "INVALID_LT";
}

rom store_type_name (StoreType store)
{
    switch (store) {
        case STORE_NONE  : return "NONE";
        case STORE_INT   : return "INT";
        case STORE_FLOAT : return "FLOAT";
        case STORE_INDEX : return "INDEX";
        default          : return "InvalidStoreType";
    }
}

rom make_ref_size_name (MakeRefSize size)
{
    static rom names[NUM_MAKE_REF] = MAKE_REF_SIZE_NAMES;

    return IN_RANGE(size, 0, NUM_MAKE_REF) ? names[size] : "InvalidMakeRefSize";
}

// called to parse the optional argument to --show-headers. we accept eg "REFERENCE" or "SEC_REFERENCE" or even "ref"
void sections_set_show_headers (rom arg)
{
    flag.show_headers = true;

    if (!arg) {
        memset (flag.show_sec_headers, true, NUM_SEC_TYPES); // all section types
        return; 
    }
    
    uint32_t arg_len = strlen(arg); 

    if (str_get_uint32 (STRa(arg), (uint32_t*)&flag.show_header_section_i)) {
        memset (flag.show_sec_headers, true, NUM_SEC_TYPES); // all section types
        return;
    }

    char upper_arg[arg_len+1];
    str_toupper_(arg, upper_arg, arg_len+1); // arg is case insesitive (compare uppercase)

    // consider multiple comma-separated section names and dicts
    str_split (upper_arg, arg_len, NUM_SEC_TYPES, ',', sec_name, false);

    for (int i=0; i < n_sec_names; i++) {
        SectionType candidate = SEC_NONE; 
        unsigned candidate_len = 100000; // arbitrary large number
        *(char *)&sec_names[i][sec_name_lens[i]] = 0; // nul-terminate

        for (SectionType st=0; st < NUM_SEC_TYPES; st++) 
            // choose the shortest section for which arg is a case-insensitive substring (eg "dict" will result in SEC_DICT not SEC_DICT_ID_ALIASES)
            if (strstr (abouts[st].name, sec_names[i])) { // found if arg is a substring of the section name
                unsigned len = strlen (abouts[st].name);
                if (len < candidate_len) { // best candidate so far
                    candidate     = st;
                    candidate_len = len;
                }
            }
        
        // case: this is a section name
        if (candidate != SEC_NONE) 
            flag.show_sec_headers[candidate] = true; 

        // case: if not a section name, we assume it is a dict name (only one is possible)
        else 
            flag.show_header_dict_name = arg;
    }
}


FlagStr sections_dis_flags (SectionFlags f, SectionType st, DataType dt/*data type of component*/,
                            bool is_r2) // relevant only if dt=FASTQ
{
    rom dts[NUM_DATATYPES]  = { [DT_FASTQ]=(!VER(15) ? " v14_is_paired=" : 0), [DT_SAM]=" is_ref_internal=" };
    rom dts2[NUM_DATATYPES] = { [DT_FASTQ]=" is_bamass=", [DT_SAM]=" is_deep=" };

    FlagStr str = {};
    rom paired_label = is_r2 ? " pair_assisted=" : " pair_identical=";

    switch (st) {
        case SEC_GENOZIP_HEADER:
            snprintf (str.s, sizeof (str.s), "aligner=%u%.24s%.24s txt_is_bin=%u %s=%u not_md5=%u has_gencomp=%u has_taxid=%u",
                     f.genozip_header.aligner, 
                     cond_int (dt != DT_NONE && dts[dt],  dts[dt],  f.genozip_header.dt_specific), 
                     cond_int (dt != DT_NONE && dts2[dt], dts2[dt], f.genozip_header.dt_specific2), 
                     f.genozip_header.txt_is_bin, 
                     VER(15) ? "has_digest" : "v14_bgzf", f.genozip_header.has_digest, 
                     f.genozip_header.not_md5, f.genozip_header.has_gencomp,f.genozip_header.has_taxid);
            break;

        case SEC_TXT_HEADER: {
            char extra[64] = {};
            if ((dt==DT_SAM || dt==DT_BAM || dt==DT_FASTQ) && VER(15)) snprintf (extra, sizeof (extra), " pair=%s", pair_type_name (f.txt_header.pair));
            snprintf (str.s, sizeof (str.s), "no_gz_ext=%u %s", f.txt_header.no_gz_ext, extra);
            break;
        }

        case SEC_VB_HEADER:
            switch (dt) {
                case DT_VCF: 
                    snprintf (str.s, sizeof (str.s), "null_DP=%u", f.vb_header.vcf.use_null_DP_method);
                    break;
                
                case DT_SAM: case DT_BAM:
                    if (VER(13) && !VER(14)) snprintf (str.s, sizeof (str.s), "sorted=%u collated=%u", f.vb_header.sam.v13_is_sorted, f.vb_header.sam.v13_is_collated);
                    break;
                
                case DT_GFF:
                    if (VER(15)) snprintf (str.s, sizeof (str.s), "embedded_fasta=%u", f.vb_header.gff.embedded_fasta);
                    break;
                
                default:
                    str.s[0] = 0;
            }
            break;

        case SEC_GZ_ISIZES:
            snprintf (str.s, sizeof (str.s), "library=%s level=%u OLD_has_eof=%u", bgzf_library_name(f.mgzip.library, false), f.mgzip.level, f.mgzip.OLD_has_eof_block); 
            break;

        case SEC_LOCAL:
            snprintf (str.s, sizeof (str.s), "store=%-5s per_ln=%u delta=%u%.22s spl_custom=%u specific=%u",
                     store_type_name(f.ctx.store), f.ctx.store_per_line, f.ctx.store_delta, 
                     cond_int (dt == DT_FASTQ, paired_label, f.ctx.paired),
                     f.ctx.spl_custom, f.ctx.ctx_specific_flag); // note: we don't print ctx_specific as its not currently used
            break;

        case SEC_B250:
            snprintf (str.s, sizeof (str.s), "store=%-5s per_ln=%u delta=%u%.22s spl_custom=%u specific=%u same=%u",
                     store_type_name(f.ctx.store), f.ctx.store_per_line, f.ctx.store_delta,
                     cond_int (dt == DT_FASTQ, paired_label, f.ctx.paired),
                     f.ctx.spl_custom, f.ctx.ctx_specific_flag, f.ctx.all_the_same); 
            break;
 
        case SEC_RANDOM_ACCESS:
        case SEC_RECON_PLAN:
            snprintf (str.s, sizeof (str.s), "frag_len_bits=%u", f.recon_plan.frag_len_bits);
            break;

        case SEC_DICT:
            snprintf (str.s, sizeof (str.s), "deep_sam=%u deep_fq=%u all_the_same_wi=%u", 
                     f.dictionary.deep_sam, f.dictionary.deep_fastq, f.dictionary.all_the_same_wi);
            break;

        default: 
            str.s[0] = 0;
    }

    return str;
}

static bool is_show_header (SectionType st, uint32_t sec_i, DictId sec_dict_id/*DICT_ID_NONE if not known*/)
{
    // case: specific section number requested
    if (flag.show_header_section_i != -1)
        return flag.show_header_section_i == sec_i;

    // case: specific dict requested
    if (flag.show_header_dict_name) {
        if (!IS_DICTED_SEC(st)) return false;
        
        return !sec_dict_id.num || // not known yet 
               !strcmp (flag.show_header_dict_name, dis_dict_id (sec_dict_id).s); // the dict user requested
    }

    // case: --debug-read-ctxs - show all dicted sections
    if (flag.debug_read_ctxs && IS_DICTED_SEC(st)) 
        return true;

    // loading reference: show reference sections only if requested explicitly
    if (flag_loading_auxiliary)
        return flag.show_sec_headers[st] && // requested this section
               memchr (flag.show_sec_headers, 0, NUM_SEC_TYPES); // but not requested all

    else 
        return flag.show_sec_headers[st];
}

void sections_show_header (ConstSectionHeaderP header, 
                           VBlockP vb,       // optional if output to buffer (allowed only in ZIP)
                           CompIType comp_i, 
                           uint64_t offset, char rw)
{
    #define DT(x) ((dt) == DT_##x)

    ASSERTW0 (!vb || IS_ZIP, "sections_show_header: expecting vb=NULL in PIZ"); // because we dump the show_headers_buf to the terminal only in ZIP

    SectionType st = header->section_type;
    
    if (!is_show_header (st, header->section_i, sections_get_dict_id (header)))
        return; // we don't need to show this section
              
    bool is_dict_offset = (header->section_type == SEC_DICT && rw == 'W'); // at the point calling this function in zip, SEC_DICT offsets are not finalized yet and are relative to the beginning of the dictionary area in the genozip file
    bool v12 = (IS_ZIP || VER(12));
    bool v14 = (IS_ZIP || VER(14));
    bool v15 = (IS_ZIP || VER(15));

    char str[4096];
    #define PRINT { if (vb) buf_append_string (vb, &vb->show_headers_buf, str); else iprintf ("%s", str); } 
    
    rom SEC_TAB = isatty(1) ? "            ++  " : " "; // single line if not terminal - eg for grepping

    snprintf (str, sizeof (str), "%c %s%-*"PRIu64" %-19s %-4.4s %.8s%.15s vb=%-3u %s=%*u txt_len=%-7u z_len=%-7u enc_len=%-7u ",
              rw, 
              is_dict_offset ? "~" : "", 9-is_dict_offset, offset, 
              st_name(header->section_type), 
              codec_name (header->codec), 
              cond_str (header->section_type == SEC_LOCAL && header->sub_codec, "sub=", codec_name (header->sub_codec)),
              cond_int (header->section_type == SEC_DICT, "dict_helper=", header->dict_helper),
              BGEN32 (header->vblock_i), 
              v15 ? "z_digest" : "comp_off",
              v15 ? -10 : -4,
              v15 ? BGEN32 (header->z_digest) : BGEN32 (header->v14_compressed_offset), 
              BGEN32 (header->data_uncompressed_len), 
              BGEN32 (header->data_compressed_len), 
              BGEN32 (header->data_encrypted_len));
    PRINT;

    SectionFlags f = header->flags;
    DataType dt    = (flag.deep && comp_i >= SAM_COMP_FQ00) ? DT_FASTQ : z_file->data_type;
    bool is_r2     = (dt == DT_FASTQ && comp_i == FQ_COMP_R2) ||
                     (flag.deep && flag.pair && comp_i % 2 == 0); // SAM_COMP_FQ01=4, SAM_COMO_FQ03=6...

    switch (st) {

    case SEC_GENOZIP_HEADER: {
        SectionHeaderGenozipHeaderP h = (SectionHeaderGenozipHeaderP)header;
        z_file->z_flags.not_md5 = h->flags.genozip_header.not_md5; // needed by digest_display_ex
        char dt_specific[REF_FILENAME_LEN + 2 KB] = "";
        dt = BGEN16 (h->data_type); // for GENOZIP_HEADER, go by what header declares
        if (dt >= NUM_DATATYPES) dt = DT_NONE;

        if ((DT(VCF) || DT(BCF)) && v14)
            snprintf (dt_specific, sizeof (dt_specific), "%ssegconf=(has_RGQ=%s,del_svlen_is_neg=%s,FMT_GQ_method=%s,FMT_DP_method=%s,INFO_DP_method=%s,mate_id_chars=%s,allow_cp_smp=%s) width=(AC=%u,AN=%u,MLEAC=%u,DP=%u,QD=%u,SF=%u,AS_SB_TABLE=%u,QUAL=%u,BaseCounts=%u,DPB=%u) max_ploidy_for_mux=%u\n", 
                      SEC_TAB, TF(h->vcf.segconf_has_RGQ), TF(h->vcf.segconf_del_svlen_is_neg), FMT_GQ_method_name (h->vcf.segconf_GQ_method), FMT_DP_method_name (h->vcf.segconf_FMT_DP_method), INFO_DP_method_name (h->vcf.segconf_INF_DP_method), 
                      (rom[])MATEID_METHOD_NAMES[h->vcf.segconf_MATEID_method], TF(h->vcf.segconf_sample_copy),
                      h->vcf.width.AC, h->vcf.width.AN, h->vcf.width.MLEAC, h->vcf.width.DP, h->vcf.width.QD, h->vcf.width.SF, h->vcf.width.AS_SB_TABLE, h->vcf.width.QUAL, h->vcf.width.BaseCounts, h->vcf.width.DPB,
                      h->vcf.max_ploidy_for_mux);

        else if ((DT(SAM) || DT(BAM)) && v14)
            snprintf (dt_specific, sizeof (dt_specific), "%ssegconf=(sorted=%u,collated=%u,std_seq_len=%u,seq_len_to_cm=%u,s1_to_cm_32=%u,ms_type=%u,has_MD_or_NM=%u,bisulfite=%u,MD_NM_by_unconverted=%u,predict_meth=%u,is_paired=%u,sag_type=%s,sag_has_AS=%u,pysam_qual=%u,cellranger=%u,SA_HtoS=%u,seq_len_dict_id=%s,deep_qname1=%u,deep_qname2=%u,deep_no_qual=%u,use_ins_ctx=%u,MAPQ_use_xq=%u,has_MQ=%u,SA_CIGAR_abb=%u,SA_NM_by_X=%u,CIGAR_has_eqx=%u,is_ileaved=%u,xcons_std_seq_len_M100=%u,%uXsam_factor=%u)\n", 
                      SEC_TAB, h->sam.segconf_is_sorted, h->sam.segconf_is_collated, BGEN32(h->sam.segconf_std_seq_len), 
                      h->sam.segconf_seq_len_cm, h->sam.segconf_s1_to_cm_32, 
                      h->sam.segconf_ms_type, h->sam.segconf_has_MD_or_NM, 
                      h->sam.segconf_bisulfite, h->sam.segconf_MD_NM_by_un, h->sam.segconf_predict_meth, 
                      h->sam.segconf_is_paired, sag_type_name(h->sam.segconf_sag_type), h->sam.segconf_sag_has_AS, 
                      h->sam.segconf_pysam_qual, h->sam.segconf_10xGen, h->sam.segconf_SA_HtoS, dis_dict_id(h->sam.segconf_seq_len_dict_id).s,
                      h->sam.segconf_deep_qname1, h->sam.segconf_deep_qname2, h->sam.segconf_deep_no_qual, 
                      h->sam.segconf_use_ins_ctxs, h->sam.segconf_MAPQ_use_xq, h->sam.segconf_has_MQ, h->sam.segconf_SA_CIGAR_abb, h->sam.segconf_SA_NM_by_X, 
                      h->sam.segconf_CIGAR_has_eqx, h->sam.segconf_is_ileaved, h->sam.xcons_std_seq_len_M100,
                      SAM_FACTOR_MULT, h->sam.segconf_sam_factor);

        else if (DT(REF)) {
            if (v15) snprintf (dt_specific, sizeof (dt_specific), "%sgenome_digest=%s%.30s%.30s%.30s%.30s\n", 
                               SEC_TAB, digest_display (h->genome_digest).s,
                               cond_int (VER2(15,81), " bases_per_hash=", h->ref.bases_per_hash),
                               cond_int (VER2(15,81), " bits_per_hash_out=", h->ref.bits_per_hash_out),
                               cond_int (VER2(15,81), " gpos_bytes=", h->ref.gpos_bytes),
                               cond_str (VER2(15,81), " make_ref_size=", make_ref_size_name(h->ref.make_ref_size)));
            else     snprintf (dt_specific, sizeof (dt_specific), "%sfasta_md5=%s\n", SEC_TAB, digest_display (h->v14_REF_fasta_md5).s);
        }

        else if (DT(FASTQ) && v14) 
            snprintf (dt_specific, sizeof (dt_specific), "%sFASTQ_v13_digest_bound=%s segconf=(seq_len_dict_id=%s,fa_as_fq=%s,is_ileaved=%s,std_seq_len=%u,use_ins_ctxs=%u)\n", 
                      SEC_TAB, digest_is_zero(h->FASTQ_v13_digest_bound) ? "N/A" : digest_display (h->FASTQ_v13_digest_bound).s, 
                      dis_dict_id(h->fastq.segconf_seq_len_dict_id).s, TF(h->fastq.segconf_fa_as_fq), TF(h->fastq.segconf_is_ileaved), BGEN32(h->fastq.segconf_std_seq_len),
                      h->fastq.segconf_use_ins_ctxs);

        snprintf (str, sizeof (str), "\n%sver=%.10s%.10s modified=%u lic=%.20s private=%u enc=%.10s dt=%.10s usize=%"PRIu64" lines=%"PRIu64" secs=%u txts=%u %.10s%.20s%.6s\n" 
                                     "%s%s %s=\"%.*s\" %s:%s%s\n"
                                     "%s" // dt_specific, if there is any
                                     "%screated=\"%.*s\"\n",
                  SEC_TAB, STRver_(h->genozip_version, h->genozip_minor_ver).s, 
                  cond_str(VER2(15,81), " ref_ver=", STRver_(h->ref_genozip_ver, BGEN16(h->ref_genozip_minor_ver)).s),
                  h->is_modified/*15.0.60*/, lic_type_name (h->lic_type)/*15.0.59*/, 
                  h->private_file, encryption_name (h->encryption_type), dt_name (dt), 
                  BGEN64 (h->recon_size), BGEN64 (h->num_lines_bound), BGEN32 (h->num_sections), h->num_txt_files,
                  cond_int(!VER2(15,65), "vb_size=", BGEN16(h->old_vb_size)),
                  cond_int(VER2(15,65), "segconf_vb_size=", BGEN32(h->segconf_vb_size)),
                  cond_int ((DT(SAM) || DT(BAM)) && v15, " conc_writing_vbs=", BGEN16(h->sam.conc_writing_vbs)), 
                  SEC_TAB, sections_dis_flags (f, st, dt, is_r2).s,
                  DT(REF) ? "fasta" : "ref", REF_FILENAME_LEN, h->ref_filename, 
                  DT(REF) ? "refhash_digest" : "ref_genome_digest", DT(REF) ? digest_display(h->refhash_digest).s : digest_display_ex_ (h->ref_genome_digest, DD_NORMAL, VER2(15,81) ? h->ref_genome_digest_alg : DIGEST_UNKNOWN, !VER2(15,81)/*guess alg for older files*/).s,
                  cond_str (!DT(REF) && VER2(15,81), " ref_genome_digest_alg=", digest_alg_name (h->ref_genome_digest_alg)),
                  dt_specific, 
                  SEC_TAB, FILE_METADATA_LEN, h->created);
        break;
    }

    case SEC_TXT_HEADER: {
        SectionHeaderTxtHeaderP h = (SectionHeaderTxtHeaderP)header;
        if (!v15)
            snprintf (str, sizeof (str), "\n%stxt_data_size=%"PRIu64" txt_header_size=%"PRIu64" lines=%"PRIu64" max_lines_per_vb=%u digest=%s digest_header=%s\n" 
                      "%ssrc_codec=%s OLD_gz_size_3LSB=%u %s txt_filename=\"%.*s\"\n",
                      SEC_TAB, BGEN64 (h->txt_data_size), v12 ? BGEN64 (h->txt_header_size) : 0, BGEN64 (h->txt_num_lines), BGEN32 (h->max_lines_per_vb), 
                      digest_display (h->digest).s, digest_display (h->digest_header).s, 
                      SEC_TAB, codec_name (h->src_codec), GET_UINT24(h->OLD_gz_size_3LSB), 
                      sections_dis_flags (f, st, dt, is_r2).s, TXT_FILENAME_LEN, h->txt_filename);
        else
            snprintf (str, sizeof (str), "\n%stxt_data_size=%"PRIu64" txt_header_size=%"PRIu64" lines=%"PRIu64" max_lines_per_vb=%u digest=%s digest_header=%s\n" 
                      "%ssrc_codec=%s OLD_gz_size_3LSB=%u %s txt_filename=\"%.*s\" flav_prop=(has_seq_len,consensus,mated,cnn,tokenized)=[[%u,%u,%u,'%s',%u],[%u,%u,%u,'%s',%u],[%u,%u,%u,'%s',%u]]\n",
                      SEC_TAB, BGEN64 (h->txt_data_size), v12 ? BGEN64 (h->txt_header_size) : 0, BGEN64 (h->txt_num_lines), BGEN32 (h->max_lines_per_vb), 
                      digest_display (h->digest).s, digest_display (h->digest_header).s, 
                      SEC_TAB, codec_name (h->src_codec), GET_UINT24(h->OLD_gz_size_3LSB),
                      sections_dis_flags (f, st, dt, is_r2).s, TXT_FILENAME_LEN, h->txt_filename,
                      h->flav_prop[0].has_seq_len, h->flav_prop[0].is_consensus, h->flav_prop[0].is_mated, char_to_printable((char[])CNN_TO_CHAR[h->flav_prop[0].cnn]).s, h->flav_prop[0].is_tokenized,
                      h->flav_prop[1].has_seq_len, h->flav_prop[1].is_consensus, h->flav_prop[1].is_mated, char_to_printable((char[])CNN_TO_CHAR[h->flav_prop[1].cnn]).s, h->flav_prop[1].is_tokenized,
                      h->flav_prop[2].has_seq_len, h->flav_prop[2].is_consensus, h->flav_prop[2].is_mated, char_to_printable((char[])CNN_TO_CHAR[h->flav_prop[2].cnn]).s, h->flav_prop[2].is_tokenized);

        break;
    }

    case SEC_VB_HEADER: {
        SectionHeaderVbHeaderP h = (SectionHeaderVbHeaderP)header;
        if Z_DT(VCF) 
            snprintf (str, sizeof (str), 
                      "\n%srecon_size=%u longest_line=%u HT_n_lines=%u z_data_bytes=%u digest=%s %s\n",
                      SEC_TAB, BGEN32 (h->recon_size), BGEN32 (h->longest_line_len), BGEN32 (h->vcf_HT_n_lines),
                      BGEN32 (h->z_data_bytes), digest_display (h->digest).s, sections_dis_flags (f, st, dt, is_r2).s);
        else if (Z_DT(SAM) && comp_i == SAM_COMP_MAIN)
            snprintf (str, sizeof (str), 
                      "\n%srecon_size=%-8u longest_line=%u longest_seq=%u z_data_bytes=%u digest=%s %s\n", 
                      SEC_TAB, BGEN32 (h->recon_size),  
                      BGEN32 (h->longest_line_len), BGEN32(h->longest_seq_len),
                      BGEN32 (h->z_data_bytes), digest_display (h->digest).s, 
                      sections_dis_flags (f, st, dt, 0).s);
        else if (Z_DT(SAM) && comp_i == SAM_COMP_PRIM)
            snprintf (str, sizeof (str), 
                      "\n%srecon_size=%-8u longest_line=%u longest_seq=%u z_data_bytes=%u digest=%s PRIM=(seq=%u comp_qual=%u qname=%u num_alns=%u first_grp_i=%u %s=%u) %s\n", 
                      SEC_TAB, BGEN32 (h->recon_size),  
                      BGEN32 (h->longest_line_len), BGEN32(h->longest_seq_len),
                      BGEN32 (h->z_data_bytes), digest_display (h->digest).s, 
                      v14 ? BGEN32 (h->sam_prim_seq_len)          : 0,
                      v14 ? BGEN32 (h->sam_prim_comp_qual_len)    : 0,
                      v14 ? BGEN32 (h->sam_prim_comp_qname_len)   : 0,
                      v14 ? BGEN32 (h->sam_prim_num_sag_alns)     : 0,
                      v14 ? BGEN32 (h->sam_prim_first_grp_i)      : 0,
                      IS_SAG_SA?"comp_cigars" : IS_SAG_SOLO?"solo_data" : "unused",
                      v14 ? BGEN32 (h->sam_prim_comp_cigars_len)  : 0,
                      sections_dis_flags (f, st, dt, 0).s);
        else if (Z_DT(SAM) && comp_i == SAM_COMP_DEPN)
            snprintf (str, sizeof (str), 
                      "\n%srecon_size=%-8u longest_line=%u longest_seq=%u z_data_bytes=%u digest=%s %sDEPN\n", 
                      SEC_TAB, BGEN32 (h->recon_size),  
                      BGEN32 (h->longest_line_len), BGEN32(h->longest_seq_len),
                      BGEN32 (h->z_data_bytes), digest_display (h->digest).s, 
                      sections_dis_flags (f, st, dt, 0).s);
        else if (Z_DT(FASTQ))
            snprintf (str, sizeof (str), 
                      "\n%srecon_size=%u longest_line=%u longest_seq=%u z_data_bytes=%u digest=%s %s\n",
                      SEC_TAB, BGEN32 (h->recon_size),BGEN32 (h->longest_line_len), BGEN32(h->longest_seq_len),
                      BGEN32 (h->z_data_bytes), digest_display (h->digest).s, sections_dis_flags (f, st, dt, is_r2).s);
        else
            snprintf (str, sizeof (str), 
                      "\n%srecon_size=%u longest_line=%u z_data_bytes=%u digest=%s %s\n",
                      SEC_TAB, BGEN32 (h->recon_size),BGEN32 (h->longest_line_len), 
                      BGEN32 (h->z_data_bytes), digest_display (h->digest).s, sections_dis_flags (f, st, dt, is_r2).s);

        break;
    }

    case SEC_REFERENCE:
    case SEC_REF_IS_SET: {
        SectionHeaderReferenceP h = (SectionHeaderReferenceP)header;
        snprintf (str, sizeof (str), "pos=%-9"PRIu64" gpos=%-9"PRIu64" num_bases=%-6u chrom_word_index=%-4d\n",
                  BGEN64 (h->pos), BGEN64 (h->gpos), BGEN32 (h->num_bases), BGEN32 (h->chrom_word_index)); 
        break;
    }
    
    case SEC_REF_HASH: {
        SectionHeaderRefHashP h = (SectionHeaderRefHashP)header;
        if (VER2(15,81))
            snprintf (str, sizeof (str), "first_ent=%"PRIu64"\n", BGEN64 (h->first_ent)); 
        else
            snprintf (str, sizeof (str), "num_layers=%u layer_i=%u layer_bits=%u start_in_layer=%u\n",
                      h->OLD.num_layers, h->OLD.layer_i, h->OLD.layer_bits, BGEN32 (h->OLD.start_in_layer)); 
        break;
    }
    
    case SEC_REF_CONTIGS: {
        snprintf (str, sizeof (str), "sequential_ref_index=%u\n", header->flags.ref_contigs.sequential_ref_index);
        break;
    }
    
    case SEC_RECON_PLAN: {
        SectionHeaderReconPlanP h = (SectionHeaderReconPlanP)header;
        snprintf (str, sizeof (str), "conc_writing_vbs=%u %s\n", BGEN32 (h->conc_writing_vbs), sections_dis_flags (f, st, dt, 0).s); 
        break;
    }
        
    case SEC_GZ_DIGESTS: {
        SectionHeaderGzDigestsP h = (SectionHeaderGzDigestsP)header;
        snprintf (str, sizeof (str), "gz_file_size=%"PRIu64" effective_codec=%s %s\n", 
                  BGEN64(h->gz_file_size), codec_name (h->effective_codec), sections_dis_flags (f, st, dt, 0).s); 
        break;
    }


    case SEC_GZ_ISIZES:
    case SEC_RANDOM_ACCESS: {
        snprintf (str, sizeof (str), "%s%s\n", SEC_TAB, sections_dis_flags (f, st, dt, 0).s); 
        break;
    }
    
    case SEC_B250: {
        SectionHeaderCtxP h = (SectionHeaderCtxP)header;
        snprintf (str, sizeof (str), "%s/%-8s\tb250_size=%s param=%u %s\n",
                  dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s,  
                  h->b250_size==B250_BYTES_1?"1" : h->b250_size==B250_BYTES_2?"2" : h->b250_size==B250_BYTES_3?"3" : h->b250_size==B250_BYTES_4?"4" : h->b250_size==B250_VARL?"VARL" : "INVALID",
                  h->param, sections_dis_flags (f, st, dt, is_r2).s);
        break;
    }

    case SEC_LOCAL: {
        SectionHeaderCtxP h = (SectionHeaderCtxP)header;
        snprintf (str, sizeof (str), "%s/%-8s\tltype=%s param=%u%s %s\n",
                  dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s, lt_name (h->ltype), h->param, 
                  cond_str (lt_max(h->ltype), " nothing_char=", char_to_printable (h->nothing_char).s), // nothing_char only defined for integer ltypes
                  sections_dis_flags (f, st, dt, is_r2).s);
        break;
    }

    case SEC_DICT: {
        SectionHeaderDictionaryP h = (SectionHeaderDictionaryP)header;
        snprintf (str, sizeof (str), "%s/%-8s\tnum_snips=%u %s\n", dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s, BGEN32 (h->num_snips), 
                  sections_dis_flags (f, st, dt, 0).s); 
        break;
    }

    case SEC_COUNTS: {
        SectionHeaderCountsP h = (SectionHeaderCountsP)header;
        snprintf (str, sizeof (str), "  %s/%-8s param=%"PRId64"\t\n", dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s, h->nodes_param); 
        break;
    }

    case SEC_SUBDICTS: {
        SectionHeaderSubDictsP h = (SectionHeaderSubDictsP)header;
        snprintf (str, sizeof (str), "  %s/%-8s param=%"PRId64"\t\n", dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s, h->param); 
        break;
    }

    case SEC_HUFFMAN: {
        SectionHeaderHuffmanP h = (SectionHeaderHuffmanP)header;
        snprintf (str, sizeof (str), "  %s/%-8s \t\n", dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s); 
        break;
    }


    default: 
        str[0] = '\n'; str[1] = 0; 
    }

    // if not going directly to the terminal, replace non-final newlines with \t, to allow grep etc
    if (!isatty(1)) str_replace_letter (str, strlen(str)-1, '\n', '\t');

    PRINT;
}

// called from main()
void noreturn genocat_show_headers (rom z_filename)
{
    TEMP_FLAG(quiet, false);
    z_file = file_open_z_read (z_filename);    
    RESTORE_FLAG (quiet);
    ASSERTNOTNULL (z_file);
    
    SectionHeaderGenozipHeader header;

    TEMP_FLAG (show_headers, false); // avoid showing the genozip header twice 
    TEMP_FLAG (genocat_no_reconstruct, 1);
    zfile_read_genozip_header (&header, HARD_FAIL); // also sets z_file->genozip_ver
    RESTORE_FLAG (show_headers);
    RESTORE_FLAG (genocat_no_reconstruct);

    // normal --show-headers - go by the section list
    if (!flag.force) {
        for_buf2 (SectionEnt, sec, sec_i, z_file->section_list)
            if (flag.show_headers && is_show_header (sec->st, sec_i, DICT_ID_NONE)) {                
                header = zfile_read_section_header (evb, sec, SEC_NONE).genozip_header; // we assign the largest of the SectionHeader* types
                header.section_i = sec_i; // note: replaces magic, 32 bit only. nonsense if sec is not in z_file->section_list.

                sections_show_header ((SectionHeaderP)&header, NULL, sec->comp_i, sec->offset, 'R');
            }        
    }

    // --show-headers --force - search for actual headers in case of file corruption
    else {
        file_seek (z_file, 0, SET, READ, HARD_FAIL);
        uint64_t gap; // gap before section
        uint64_t accumulated_gap = 0;
        SectionEntModifiable sec = { .st = SEC_GENOZIP_HEADER/*largest header*/ };

        for (int sec_i=0; zfile_advance_to_next_header (&sec.offset, &gap); sec_i++) {
            if (gap || accumulated_gap) 
                iprintf ("ERROR: unexpected of %"PRIu64" bytes before next section\n", gap + accumulated_gap);
            
            header = zfile_read_section_header (evb, &sec, SEC_NONE).genozip_header; // we assign the largest of the SectionHeader* types
            if (header.section_type < 0 || header.section_type >= NUM_SEC_TYPES) { // not true section - magic matches by chance
                sec.offset += 4;
                accumulated_gap += 4;
                continue;
            }
            
            int header_size = st_header_size (header.section_type);
            if (header.data_encrypted_len) header_size = ROUNDUP16(header_size);

            // up to v14 we verify v14_compressed_offset
            if (!VER(15) && header_size != BGEN32 (header.v14_compressed_offset)) {
                sec.offset += 4;
                accumulated_gap += 4;
                continue;
            }

            if (flag.show_sec_headers[header.section_type])
                iprintf ("%5u ", sec_i);
            
            header.section_i = sec_i;
            sections_show_header ((SectionHeaderP)&header, NULL, sec.comp_i, sec.offset, 'R');

            sec.offset += header_size + BGEN32 (header.data_compressed_len);
            accumulated_gap = 0;
        }

        if (flag.show_headers) {
            if (gap) iprintf ("ERROR: unexpected gap of %"PRIu64" bytes before Footer\n", gap);

            if ((sec.offset = zfile_read_genozip_header_get_offset (true)))
                iprintf ("R %9"PRIu64" FOOTER              genozip_header_offset=%"PRIu64"\n", 
                        z_file->disk_size - sizeof (SectionFooterGenozipHeader), sec.offset);
            else
                iprint0 ("ERROR: no valid Footer\n");
        }
    }

    exit_ok;
}

void sections_show_section_list (DataType dt, BufferP section_list, SectionType only_this_st/*SEC_NONE if all*/)
{    
    for_buf (SectionEnt, s, *section_list)
        if (only_this_st != SEC_NONE && s->st != only_this_st)
            continue;

        else if (IS_B250(s) || IS_LOCAL(s) || s->st == SEC_DICT) {
            DataType my_dt = (s->st != SEC_DICT && flag.deep && s->comp_i >= SAM_COMP_FQ00) ? DT_FASTQ : dt;

            bool is_r2 = (dt == DT_FASTQ && s->comp_i == FQ_COMP_R2) ||
                         (flag.deep && flag.pair && s->comp_i % 2 == 0); // SAM_COMP_FQ01=4, SAM_COMO_FQ03=6...

            iprintf ("%5u %-20.20s %s%s%-8.8s\tvb=%s/%-4u offset=%-10"PRIu64"  size=%-6u  %s\n", 
                     BNUM(*section_list, s), st_name(s->st), 
                     s->dict_id.num ? dtype_name_z(s->dict_id) :"     ", 
                     s->dict_id.num ? "/" : "", 
                     s->dict_id.num ? dis_dict_id (s->dict_id).s : "", 
                     comp_name_ex (s->comp_i, s->st).s, s->vblock_i, s->offset, s->size, 
                     sections_dis_flags (s->flags, s->st, my_dt, is_r2).s);
        }
        
        else if (IS_VB_HEADER(s))
            iprintf ("%5u %-20.20s\t\t\tvb=%s/%-4u offset=%-10"PRIu64"  size=%-6u  num_lines=%-8u%s\n", 
                     BNUM(*section_list, s), st_name(s->st), comp_name (s->comp_i), 
                     s->vblock_i, s->offset, s->size, s->num_lines, 
                     sections_dis_flags (s->flags, s->st, dt, 0).s);

        else if (IS_FRAG_SEC(s->st) || s->st == SEC_GZ_ISIZES || s->st == SEC_GZ_DIGESTS)
            iprintf ("%5u %-20.20s\t\t\tvb=%s/%-4u offset=%-10"PRIu64"  size=%-6u  %s\n", 
                     BNUM(*section_list, s), st_name(s->st), 
                     comp_name_ex (s->comp_i, s->st).s, s->vblock_i, s->offset, s->size, sections_dis_flags (s->flags, s->st, dt, 0).s);

        else if (IS_DICTED_SEC(s->st))
            iprintf ("%5u %-20.20s %s%s%-8.8s\t             offset=%-10"PRIu64"  size=%-6u  %s\n", 
                     BNUM(*section_list, s), st_name(s->st),
                     s->dict_id.num ? dtype_name_z(s->dict_id) :"     ", 
                     s->dict_id.num ? "/" : "", 
                     s->dict_id.num ? dis_dict_id (s->dict_id).s : "", 
                     s->offset, s->size, sections_dis_flags (s->flags, s->st, dt, 0).s);
        else
            iprintf ("%5u %-20.20s\t\t\t             offset=%-10"PRIu64"  size=%-6u  %s\n", 
                     BNUM(*section_list, s), st_name(s->st),
                     s->offset, s->size, sections_dis_flags (s->flags, s->st, dt, 0).s);
}

void sections_show_gheader (ConstSectionHeaderGenozipHeaderP header)
{
    bool v15 = (IS_ZIP || VER(15));

    if (flag_loading_auxiliary) return; // don't show gheaders of an auxiliary file
    
    DataType dt = BGEN16 (header->data_type);
    ASSERT (dt < NUM_DATATYPES, "Invalid data_type=%u", dt);
    
    if (header) {
        iprintf ("Contents of the SEC_GENOZIP_HEADER section (output of --show-gheader) of %s:\n", z_name);
        iprintf ("  genozip_version: %s\n",         STRver_(header->genozip_version, header->genozip_minor_ver).s); // note: minor version always 0 before 15.0.28
        iprintf ("  data_type: %s\n",               dt_name (dt));
        iprintf ("  encryption_type: %s\n",         encryption_name (header->encryption_type)); 
        iprintf ("  recon_size: %s\n",              str_int_commas (BGEN64 (header->recon_size)).s);
        iprintf ("  num_lines_bound: %"PRIu64"\n",  BGEN64 (header->num_lines_bound));
        iprintf ("  num_sections: %u\n",            z_file->section_list.len32);
        iprintf ("  num_txt_files: %u\n",           header->num_txt_files);
        iprintf ("  modified_by_zip: %s\n", TF(header->is_modified));
        if (dt == DT_REF)
            iprintf ("  %s: %s\n", (v15 ? "genome_digest" : "REF_fasta_md5"), digest_display (header->genome_digest).s);
        iprintf ("  created: %*s\n",                -FILE_METADATA_LEN, header->created);
        iprintf ("  license_hash: %s\n",            digest_display (header->license_hash).s);
        iprintf ("  private_file: %s\n",            TF(header->private_file));
        if (header->ref_filename[0]) {
            iprintf ("  reference filename: %s\n",  header->ref_filename);
            iprintf ("  reference file hash: %s\n", digest_display (header->ref_genome_digest).s);
        }
        iprintf ("  flags: %s\n",                   sections_dis_flags (header->flags, SEC_GENOZIP_HEADER, dt, 0).s);
    }

    iprint0 ("  sections:\n");

    sections_show_section_list (dt, &z_file->section_list, SEC_NONE);
}
