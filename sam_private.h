// ------------------------------------------------------------------
//   sam_private.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef SAM_PRIVATE_INCLUDED
#define SAM_PRIVATE_INCLUDED

#include "sam.h"
#include "vblock.h"
#include "reference.h"
#include "contigs.h"

#define DTYPE_QNAME        DTYPE_1
#define DTYPE_SAM_OPTIONAL DTYPE_2

#pragma GENDICT_PREFIX SAM

// Fields
#pragma GENDICT SAM_RNAME=DTYPE_FIELD=RNAME // RNAME must be first
#pragma GENDICT SAM_QNAME=DTYPE_FIELD=QNAME
#pragma GENDICT SAM_FLAG=DTYPE_FIELD=FLAG
#pragma GENDICT SAM_POS=DTYPE_FIELD=POS
#pragma GENDICT SAM_MAPQ=DTYPE_FIELD=MAPQ
#pragma GENDICT SAM_CIGAR=DTYPE_FIELD=CIGAR
#pragma GENDICT SAM_RNEXT=DTYPE_FIELD=RNEXT
#pragma GENDICT SAM_PNEXT=DTYPE_FIELD=PNEXT
#pragma GENDICT SAM_TLEN=DTYPE_FIELD=TLEN
#pragma GENDICT SAM_OPTIONAL=DTYPE_FIELD=OPTIONAL
#pragma GENDICT SAM_SQBITMAP=DTYPE_FIELD=SQBITMAP
#pragma GENDICT SAM_NONREF=DTYPE_FIELD=NONREF    // these 4 fields must be in this order, right after SAM_SQBITMAP
#pragma GENDICT SAM_NONREF_X=DTYPE_FIELD=NONREF_X
#pragma GENDICT SAM_GPOS=DTYPE_FIELD=GPOS
#pragma GENDICT SAM_STRAND=DTYPE_FIELD=STRAND
#pragma GENDICT SAM_QUAL=DTYPE_FIELD=QUAL 
#pragma GENDICT SAM_DOMQRUNS=DTYPE_FIELD=DOMQRUNS // must be right after SAM_QUAL
#pragma GENDICT SAM_EOL=DTYPE_FIELD=EOL
#pragma GENDICT SAM_BAM_BIN=DTYPE_FIELD=BAM_BIN
#pragma GENDICT SAM_TOPLEVEL=DTYPE_FIELD=TOPLEVEL // must be called TOPLEVEL
#pragma GENDICT SAM_TOP2BAM=DTYPE_FIELD=TOP2BAM
#pragma GENDICT SAM_TOP2FQ=DTYPE_FIELD=TOP2FQ
#pragma GENDICT SAM_TOP2FQEX=DTYPE_FIELD=TOP2FQEX
#pragma GENDICT SAM_TAXID=DTYPE_FIELD=TAXID

#pragma GENDICT OPTION_AM=DTYPE_2=AM:i
#pragma GENDICT OPTION_AS=DTYPE_2=AS:i
#pragma GENDICT OPTION_CC=DTYPE_2=CC:Z
#pragma GENDICT OPTION_CM=DTYPE_2=CM:i
#pragma GENDICT OPTION_E2=DTYPE_2=E2:Z
#pragma GENDICT OPTION_2NONREF=DTYPE_2=N2ONREF // these 4 fields must be in this order, right after OPTION_E2
#pragma GENDICT OPTION_N2ONREFX=DTYPE_2=n2ONREFX
#pragma GENDICT OPTION_2GPOS=DTYPE_FIELD=G2POS
#pragma GENDICT OPTION_S2TRAND=DTYPE_2=S2TRAND
#pragma GENDICT OPTION_FI=DTYPE_2=FI:i
#pragma GENDICT OPTION_H0=DTYPE_2=H0:i
#pragma GENDICT OPTION_H1=DTYPE_2=H1:i
#pragma GENDICT OPTION_H2=DTYPE_2=H2:i
#pragma GENDICT OPTION_LB=DTYPE_2=LB:Z
#pragma GENDICT OPTION_MC=DTYPE_2=MC:Z
#pragma GENDICT OPTION_MD=DTYPE_2=MD:Z
#pragma GENDICT OPTION_MQ=DTYPE_2=MQ:i
#pragma GENDICT OPTION_NH=DTYPE_2=NH:i
#pragma GENDICT OPTION_NM=DTYPE_2=NM:i
#pragma GENDICT OPTION_OA=DTYPE_2=OA:Z
#pragma GENDICT OPTION_OC=DTYPE_2=OC:Z
#pragma GENDICT OPTION_PG=DTYPE_2=PG:Z
#pragma GENDICT OPTION_PQ=DTYPE_2=PQ:i
#pragma GENDICT OPTION_PU=DTYPE_2=PU:Z
#pragma GENDICT OPTION_RG=DTYPE_2=RG:Z
#pragma GENDICT OPTION_SA=DTYPE_2=SA:Z
#pragma GENDICT OPTION_SM=DTYPE_2=SM:i
#pragma GENDICT OPTION_TC=DTYPE_2=TC:i
#pragma GENDICT OPTION_UQ=DTYPE_2=UQ:i
#pragma GENDICT OPTION_U2=DTYPE_2=U2:Z
#pragma GENDICT OPTION_D2OMQRUN=DTYPE_2=D2OMQRUN // must be right after OPTION_U2
                
// Ion Torrent flow signal array
#pragma GENDICT OPTION_ZM=DTYPE_2=ZM:B

// bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
#pragma GENDICT OPTION_X0=DTYPE_2=X0:i 
#pragma GENDICT OPTION_X1=DTYPE_2=X1:i 
#pragma GENDICT OPTION_XA=DTYPE_2=XA:Z 
#pragma GENDICT OPTION_XN=DTYPE_2=XN:i 
#pragma GENDICT OPTION_XM=DTYPE_2=XM:i 
#pragma GENDICT OPTION_XO=DTYPE_2=XO:i
#pragma GENDICT OPTION_XG=DTYPE_2=XG:i 
#pragma GENDICT OPTION_XS=DTYPE_2=XS:i 
#pragma GENDICT OPTION_XE=DTYPE_2=XE:i

// subfields used for XA,SA,OA
#pragma GENDICT OPTION_XA_RNAME=DTYPE_2=X0A_RNAME 
#pragma GENDICT OPTION_XA_STRAND=DTYPE_2=X1A_STRAND 
#pragma GENDICT OPTION_XA_POS=DTYPE_2=X2A_POS 
#pragma GENDICT OPTION_XA_CIGAR=DTYPE_2=X3A_CIGAR 
#pragma GENDICT OPTION_XA_NM=DTYPE_2=X4A_NM 
#pragma GENDICT OPTION_XA_STRAND_POS=DTYPE_2=X5A_STRAND_POS 

#pragma GENDICT OPTION_OA_RNAME=DTYPE_2=O0A_RNAME 
#pragma GENDICT OPTION_OA_STRAND=DTYPE_2=O1A_STRAND 
#pragma GENDICT OPTION_OA_POS=DTYPE_2=O2A_POS 
#pragma GENDICT OPTION_OA_CIGAR=DTYPE_2=O3A_CIGAR 
#pragma GENDICT OPTION_OA_NM=DTYPE_2=O4A_NM 
#pragma GENDICT OPTION_OA_MAPQ=DTYPE_2=O5A_MAPQ 

#pragma GENDICT OPTION_SA_RNAME=DTYPE_2=S0A_RNAME 
#pragma GENDICT OPTION_SA_STRAND=DTYPE_2=S1A_STRAND 
#pragma GENDICT OPTION_SA_POS=DTYPE_2=S2A_POS 
#pragma GENDICT OPTION_SA_CIGAR=DTYPE_2=S3A_CIGAR 
#pragma GENDICT OPTION_SA_NM=DTYPE_2=S4A_NM 
#pragma GENDICT OPTION_SA_MAPQ=DTYPE_2=S5A_MAPQ 

// biobambam tags
#pragma GENDICT OPTION_mc=DTYPE_2=mc:i
#pragma GENDICT OPTION_ms=DTYPE_2=ms:i

// added by GATK's BQSR (Base Quality Score Recalibration)
#pragma GENDICT OPTION_BD=DTYPE_2=BD:Z // not used in newer versions of GATK
#pragma GENDICT OPTION_BI=DTYPE_2=BI:Z // not used in newer versions of GATK
#pragma GENDICT OPTION_BD_BI=DTYPE_2=BD_BI

// our private dictionary for + or 0 strands
#pragma GENDICT OPTION_STRAND=DTYPE_2=@STRAND
#pragma GENDICT OPTION_RNAME=DTYPE_2=@RNAME
#pragma GENDICT OPTION_POS=DTYPE_2=@POS
#pragma GENDICT OPTION_CIGAR=DTYPE_2=@CIGAR
#pragma GENDICT OPTION_MAPQ=DTYPE_2=@MAPQ
#pragma GENDICT OPTION_TX=DTYPE_2=TX:i

#define BAM_MAGIC "BAM\1" // first 4 characters of a BAM file

#define IS_BAM (command==ZIP ? txt_file->data_type==DT_BAM \
                             : (z_file->data_type==DT_SAM && z_file->z_flags.txt_is_bin))

// as defined in https://samtools.github.io/hts-specs/SAMv1.pdf 1.4.2
#define SAM_FLAG_MULTI_SEGMENTS 0x0001
#define SAM_FLAG_IS_ALIGNED     0x0002
#define SAM_FLAG_UNMAPPED       0x0004
#define SAM_FLAG_NEXT_UNMAPPED  0x0008
#define SAM_FLAG_REV_COMP       0x0010
#define SAM_FLAG_NEXT_REV_COMP  0x0020
#define SAM_FLAG_IS_FIRST       0x0040
#define SAM_FLAG_IS_LAST        0x0080
#define SAM_FLAG_SECONDARY      0x0100
#define SAM_FLAG_FAILED_FILTERS 0x0200
#define SAM_FLAG_DUPLICATE      0x0400
#define SAM_FLAG_SUPPLEMENTARY  0x0800

typedef struct {
    uint32_t qual_data_start, u2_data_start, bdbi_data_start[2]; // start within vb->txt_data
    uint32_t qual_data_len, u2_data_len; // length within vb->txt_data
    uint32_t seq_len;        // actual sequence length determined from any or or of: CIGAR, SEQ, QUAL. If more than one contains the length, they must all agree
} ZipDataLineSAM;

typedef struct VBlockSAM {
    VBLOCK_COMMON_FIELDS
    const char *last_cigar;        // ZIP/PIZ: last CIGAR
    Buffer textual_cigar;          // ZIP: Seg of BAM
    Buffer textual_seq;            // ZIP: Seg of BAM
    Buffer textual_opt;            // ZIP: Seg of BAM
    uint32_t ref_consumed;         // ZIP/PIZ: how many bp of reference are consumed according to the last_cigar
    uint32_t ref_and_seq_consumed; // ZIP: how many bp in the last seq consumes both ref and seq, according to CIGAR
    Buffer bd_bi_line;             // ZIP: interlaced BD and BI data for one line

    uint32_t soft_clip;            // PIZ: number of bases that were soft-clipped in this line
    
    // data used in genocat --show-sex
    WordIndex x_index, y_index, a_index;    // word index of the X, Y and chr1 chromosomes
    uint64_t x_bases, y_bases, a_bases;     // counters of number chromosome X, Y and chr1 bases
} VBlockSAM;

// fixed-field part of a BAM alignment, see https://samtools.github.io/hts-specs/SAMv1.pdf
typedef struct __attribute__ ((__packed__)) {
    uint32_t block_size;
    int32_t ref_id;
    int32_t pos;
    uint8_t l_read_name;
    uint8_t mapq;
    uint16_t bin;
    uint16_t n_cigar_op;
    uint16_t flag;
    uint32_t l_seq;
    int32_t next_ref_id;
    int32_t next_pos;  
    int32_t tlen;
    char read_name[]; // char[l_read_name]
} BAMAlignmentFixed;

typedef VBlockSAM *VBlockSAMP;

#define DATA_LINE(i) ENT (ZipDataLineSAM, vb->lines, i)

#define MAX_POS_SAM ((PosType)0x7fffffff)

#define CIGAR_DIGIT              1
#define CIGAR_CONSUMES_QUERY     2
#define CIGAR_CONSUMES_REFERENCE 4
extern const uint8_t cigar_lookup_sam[256];
extern const uint8_t cigar_lookup_bam[16];

// loading a Little Endian uint32_t from an unaligned buffer
#define GET_UINT16(p) (((uint8_t*)p)[0] | (((uint8_t*)p)[1] << 8))
#define GET_UINT32(p) (((uint8_t*)p)[0] | (((uint8_t*)p)[1] << 8) | (((uint8_t*)p)[2] << 16) | (((uint8_t*)p)[3] << 24))

// getting integers from the BAM data
#define NEXT_UINT8  *((const uint8_t *)next_field++)
#define NEXT_UINT16 GET_UINT16 (next_field); next_field += sizeof (uint16_t);
#define NEXT_UINT32 GET_UINT32 (next_field); next_field += sizeof (uint32_t);

#define dict_id_is_sam_qname_sf dict_id_is_type_1
#define dict_id_sam_qname_sf dict_id_type_1

extern void sam_seg_qname_field (VBlockSAM *vb, STRp(qname), unsigned add_additional_bytes);
extern void sam_seg_rname_rnext (VBlockP vb, DidIType did_i, STRp (chrom), unsigned add_bytes);
extern void sam_analyze_cigar (VBlockSAMP vb, STRp(cigar), unsigned *seq_consumed, unsigned *ref_consumed, unsigned *seq_and_ref, unsigned *coverage);
extern void sam_seg_tlen_field (VBlockSAM *vb, STRp(tlen), int64_t tlen_value, PosType pnext_pos_delta, int32_t cigar_seq_len);
extern void sam_seg_qual_field (VBlockSAM *vb, ZipDataLineSAM *dl, const char *qual, uint32_t qual_data_len, unsigned add_bytes);
extern void sam_seg_seq_field (VBlockSAM *vb, DidIType bitmap_did, STRp(seq), PosType pos, const char *cigar, unsigned recursion_level, uint32_t level_0_seq_len, const char *level_0_cigar, unsigned add_bytes);
extern const char *sam_seg_optional_all (VBlockSAM *vb, ZipDataLineSAM *dl, const char *next_field, int32_t len, bool *has_13, char separator, const char *after_field);
extern const char *bam_get_one_optional (VBlockSAM *vb, const char *next_field, const char **tag, char *type, const char **value, unsigned *value_len);
extern uint16_t bam_reg2bin (int32_t first_pos, int32_t last_pos);
extern void bam_seg_bin (VBlockSAM *vb, uint16_t bin, uint16_t flag, PosType this_pos);
extern void sam_seg_verify_rname_pos (VBlockP vb, const char *p_into_txt, PosType this_pos);

#endif
