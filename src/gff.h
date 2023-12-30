// ------------------------------------------------------------------
//   gff.h
//   Copyright (C) 2020-2024 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

#define DTYPE_GFF_ATTR DTYPE_1
#define DTYPE_GFF_ENST DTYPE_2

#pragma GENDICT_PREFIX GFF

// -----------------------------------------------------------------------------------------------------------
// Common contexts of FASTA and GFF - these MUST be first; in same order; same dict_id
// -----------------------------------------------------------------------------------------------------------

#pragma GENDICT GFF_SEQID=DTYPE_FIELD=SEQID             // must be 1st as this is GFF's CHROM. We don't use FASTA_CONTIG for gff-embedded FASTAs
#pragma GENDICT GFF_FASTA_LINEMETA=DTYPE_FIELD=LINEMETA
#pragma GENDICT GFF_EOL=DTYPE_FIELD=EOL
#pragma GENDICT GFF_FASTA_DESC=DTYPE_FIELD=DESC
#pragma GENDICT GFF_COMMENT=DTYPE_FIELD=COMMENT
#pragma GENDICT GFF_FASTA_NONREF=DTYPE_FIELD=NONREF 
#pragma GENDICT GFF_FASTA_NONREF_X=DTYPE_FIELD=NONREF_X
#pragma GENDICT GFF_TOPLEVEL=DTYPE_FIELD=TOPLEVEL
#pragma GENDICT GFF_FASTA_TAXID=DTYPE_FIELD=TAXID
#pragma GENDICT GFF_DEBUG_LINES=DTYPE_FIELD=DBGLINES    // used by --debug-lines

// -----------------------------------------------------------------------------------------------------------
// End of common contexts of FASTA and GFF
// -----------------------------------------------------------------------------------------------------------

#pragma GENDICT GFF_SOURCE=DTYPE_FIELD=SOURCE
#pragma GENDICT GFF_TYPE=DTYPE_FIELD=TYPE
#pragma GENDICT GFF_START=DTYPE_FIELD=START
#pragma GENDICT GFF_END=DTYPE_FIELD=END
#pragma GENDICT GFF_SCORE=DTYPE_FIELD=SCORE
#pragma GENDICT GFF_STRAND=DTYPE_FIELD=STRAND
#pragma GENDICT GFF_PHASE=DTYPE_FIELD=PHASE
#pragma GENDICT GFF_ATTRS=DTYPE_FIELD=ATTRS

// standard GFF3 attributes defined in https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
#pragma GENDICT ATTR_ID=DTYPE_1=ID
#pragma GENDICT ATTR_Name=DTYPE_1=Name
#pragma GENDICT ATTR_Alias=DTYPE_1=Alias
#pragma GENDICT ATTR_Parent=DTYPE_1=Parent
#pragma GENDICT ATTR_ParentItem=DTYPE_1=parentItm           // an item of the Parent array (defined explicity and cannot be changed, as defined in alias list, including of old files)
#pragma GENDICT ATTR_Target=DTYPE_1=Target
#pragma GENDICT ATTR_Target_ID=DTYPE_1=T1gtID
#pragma GENDICT ATTR_Target_POS=DTYPE_1=T2gtPOS
#pragma GENDICT ATTR_Target_STRAND=DTYPE_1=T3gtSTRAND
#pragma GENDICT ATTR_Gap=DTYPE_1=Gap
#pragma GENDICT ATTR_Derives_from=DTYPE_1=Derives_from
#pragma GENDICT ATTR_Note=DTYPE_1=Note
#pragma GENDICT ATTR_Dbxref=DTYPE_1=Dbxref
#pragma GENDICT ATTR_db_xref=DTYPE_1=db_xref
#pragma GENDICT ATTR_DBXdb=DTYPE_1=DBXdb
#pragma GENDICT ATTR_DBXid=DTYPE_1=DBxid
#pragma GENDICT ATTR_Ontology_term=DTYPE_1=Ontology_term
#pragma GENDICT ATTR_Is_circular=DTYPE_1=Is_circular

// standard GVF fields 
#pragma GENDICT ATTR_Variant_seq=DTYPE_1=Variant_seq
#pragma GENDICT ATTR_Reference_seq=DTYPE_1=Reference_seq
#pragma GENDICT ATTR_Variant_freq=DTYPE_1=Variant_freq

// fields added in the GVFs of GRCh37/38
#pragma GENDICT ATTR_ancestral_allele=DTYPE_1=ancestral_allele
#pragma GENDICT ATTR_Variant_effect=DTYPE_1=Variant_effect
#pragma GENDICT ATTR_sift_prediction=DTYPE_1=sift_prediction
#pragma GENDICT ATTR_polyphen_prediction=DTYPE_1=polyphen_prediction
#pragma GENDICT ATTR_variant_peptide=DTYPE_1=variant_peptide

// GTF fields
#pragma GENDICT ATTR_gene_id=DTYPE_1=gene_id                // eg "ENSMUSG00000051951.5"
#pragma GENDICT ATTR_gene_type=DTYPE_1=gene_type            // eg "snRNA"
#pragma GENDICT ATTR_gene_status=DTYPE_1=gene_status        // eg "KNOWN"
#pragma GENDICT ATTR_gene_name=DTYPE_1=gene_name            // eg "Xkr4"
#pragma GENDICT ATTR_transcript_id=DTYPE_1=transcript_id    // eg "ENSMUST00000159265.1"
#pragma GENDICT ATTR_transcript_name=DTYPE_1=transcript_name// eg "Xkr4-202"
#pragma GENDICT ATTR_transcript_name_gene=DTYPE_1=t0ranscript_name
#pragma GENDICT ATTR_transcript_name_num=DTYPE_1=t1ranscript_name

#pragma GENDICT ATTR_protein_id=DTYPE_1=protein_id          // eg "ENSMUSP00000070648.4"
#pragma GENDICT ATTR_ccds_id=DTYPE_1=ccds_id                // eg "CCDS47836"

#pragma GENDICT ATTR_exon_id=DTYPE_1=exon_id                // eg "ENSMUSE00000485541.3"
#pragma GENDICT ATTR_exon_number=DTYPE_1=exon_number        // eg "1"

// Gencode GTF: https://www.gencodegenes.org/pages/data_format.html
#pragma GENDICT ATTR_havana_gene=DTYPE_1=havana_gene        // eg "OTTMUSG00000026353.2" 
#pragma GENDICT ATTR_havana_transcript=DTYPE_1=havana_transcript    // eg "OTTMUST00000065166.1"
#pragma GENDICT ATTR_ccdsid=DTYPE_1=ccdsid                  // eg "CCDS14803.1"

// MGKit (Metagenomics): https://mgkit.readthedocs.io/en/latest/gff.html#gff-specs
#pragma GENDICT ATTR_db=DTYPE_1=db                          // identifies the database used to make the gene_id prediction: any string, like UNIPROT-SP, UNIPROT-TR, NCBI-NT
#pragma GENDICT ATTR_taxon_db=DTYPE_1=taxon_db              // identifies the database used to make the taxon_id prediction: any string, like UNIPROT-SP, UNIPROT-TR, NCBI-NT
#pragma GENDICT ATTR_dbq=DTYPE_1=dbq                        // identifies the quality of the database, used when filtering annotations: integer
#pragma GENDICT ATTR_taxon_id=DTYPE_1=taxon_id              // identifies the annotation taxon, NCBI taxonomy is used: integer
#pragma GENDICT ATTR_uid=DTYPE_1=uid                        // unique identifier for the annotation, any string is accepted but a value is assigned by using uuid.uuid4(): string
#pragma GENDICT ATTR_cov=DTYPE_1=cov                        // coverage for the annotation over all samples, keys ending with _cov indicates coverage for each sample : integer
#pragma GENDICT ATTR_exp_syn=DTYPE_1=exp_syn                // expected number of synonymous changes for the annotation
#pragma GENDICT ATTR_exp_nonsyn=DTYPE_1=exp_nonsyn          // expected number of non-synonymous changes for the annotation
#pragma GENDICT ATTR_taxon_name=DTYPE_1=taxon_name          // name of the taxon : string
#pragma GENDICT ATTR_lineage=DTYPE_1=lineage                // taxon lineage : string
#pragma GENDICT ATTR_EC=DTYPE_1=EC                          // list of EC numbers associated to the annotation : comma separated values

// prodigal (prokaryote gene prediction): https://github.com/hyattpd/prodigal/wiki
#pragma GENDICT ATTR_partial=DTYPE_1=partial
#pragma GENDICT ATTR_start_type=DTYPE_1=start_type
#pragma GENDICT ATTR_rbs_motif=DTYPE_1=rbs_motif
#pragma GENDICT ATTR_rbs_spacer=DTYPE_1=rbs_spacer
#pragma GENDICT ATTR_gc_cont=DTYPE_1=gc_cont
#pragma GENDICT ATTR_conf=DTYPE_1=conf
#pragma GENDICT ATTR_score=DTYPE_1=score
#pragma GENDICT ATTR_cscore=DTYPE_1=cscore
#pragma GENDICT ATTR_sscore=DTYPE_1=sscore
#pragma GENDICT ATTR_rscore=DTYPE_1=rscore
#pragma GENDICT ATTR_uscore=DTYPE_1=uscore
#pragma GENDICT ATTR_tscore=DTYPE_1=tscore

// other fields
#pragma GENDICT ATTR_chr=DTYPE_1=chr

#pragma GENDICT ENSTid=DTYPE_2=ENSTid 

// ZIP Stuff
extern void gff_zip_initialize (void);
extern bool is_gff (STRp(header), bool *need_more);
extern bool gff_header_inspect (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags);
extern int32_t gff_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern rom gff_seg_txt_line (VBlockP vb_, rom field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern void gff_seg_initialize (VBlockP vb_);
extern void gff_seg_finalize (VBlockP vb);
extern bool gff_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern bool gff_seg_is_big (ConstVBlockP vb, DictId dict_id, DictId st_dict_id);
extern void gff_reset_line (VBlockP vb);

// VBlock stuff
extern unsigned gff_vb_size (DataType dt);

// PIZ stuff
extern bool gff_piz_init_vb (VBlockP vb, ConstSectionHeaderVbHeaderP header, uint32_t *txt_data_so_far_single_0_increment);
extern CONTAINER_FILTER_FUNC (gff_piz_filter);
extern CONTAINER_CALLBACK (gff_piz_container_cb);

#define GFF_SPECIAL {  gff_piz_special_exon_number }
SPECIAL (GFF, 0,  exon_number, gff_piz_special_exon_number);
#define NUM_GFF_SPECIAL 1

#define GFF_DICT_ID_ALIASES                                              \
    /*        type        alias                   maps to this ctx   */  \
    { DT_GFF, ALIAS_CTX,  _ATTR_Variant_seq,      _ATTR_Reference_seq }, \
    { DT_GFF, ALIAS_CTX,  _ATTR_ancestral_allele, _ATTR_Reference_seq }, 
