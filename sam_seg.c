// ------------------------------------------------------------------
//   sam_seg.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "sam_private.h"
#include "reference.h"
#include "seg.h"
#include "context.h"
#include "file.h"
#include "random_access.h"
#include "endianness.h"
#include "strings.h"
#include "zip.h"
#include "optimize.h"
#include "dict_id.h"
#include "codec.h"
#include "aligner.h"
#include "container.h"
#include "stats.h"
#include "txtheader.h"
#include "kraken.h"

static const StoreType optional_field_store_flag[256] = {
    ['c']=STORE_INT, ['C']=STORE_INT, 
    ['s']=STORE_INT, ['S']=STORE_INT,
    ['i']=STORE_INT, ['I']=STORE_INT,
    ['f']=STORE_FLOAT
};

static const char optional_sep_by_type[2][256] = { { // compressing from SAM
        ['c']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['C']=CI_NATIVE_NEXT | CI_TRANS_NOR, // reconstruct number and \t separator is SAM, and don't reconstruct anything if BAM (reconstruction will be done by translator)
        ['s']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['S']=CI_NATIVE_NEXT | CI_TRANS_NOR, // -"-
        ['i']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['I']=CI_NATIVE_NEXT | CI_TRANS_NOR, // -"-
        ['f']=CI_NATIVE_NEXT | CI_TRANS_NOR,                                      // compressing SAM - a float is stored as text, and when piz with translate to BAM - is not reconstructed, instead - translated
        ['Z']=CI_NATIVE_NEXT | CI_TRANS_NUL, ['H']=CI_NATIVE_NEXT | CI_TRANS_NUL, // reconstruct text and then \t seperator if SAM and \0 if BAM 
        ['A']=CI_NATIVE_NEXT,                                                     // reconstruct character and then \t seperator if SAM and no seperator for BAM
        ['B']=CI_NATIVE_NEXT                                                      // reconstruct array and then \t seperator if SAM and no seperator for BAM
}, 
{ // compressing from BAM
        ['c']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['C']=CI_NATIVE_NEXT | CI_TRANS_NOR, // reconstruct number and \t separator is SAM, and don't reconstruct anything if BAM (reconstruction will be done by translator)
        ['s']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['S']=CI_NATIVE_NEXT | CI_TRANS_NOR, // -"-
        ['i']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['I']=CI_NATIVE_NEXT | CI_TRANS_NOR, // -"-
        ['f']=CI_NATIVE_NEXT,                                                     // compressing SAM - a float is stored as a SPECIAL, and the special reconstructor handles the SAM and BAM reconstructing
        ['Z']=CI_NATIVE_NEXT | CI_TRANS_NUL, ['H']=CI_NATIVE_NEXT | CI_TRANS_NUL, // reconstruct text and then \t seperator if SAM and \0 if BAM 
        ['A']=CI_NATIVE_NEXT,                                                     // reconstruct character and then \t seperator if SAM and no seperator for BAM
        ['B']=CI_NATIVE_NEXT                                                      // reconstruct array and then \t seperator if SAM and no seperator for BAM
} };

static char taxid_redirection_snip[100];
static unsigned taxid_redirection_snip_len;

// ----------------------
// Compressor callbacks
// ----------------------

// callback function for compress to get data of one line (called by codec_bz2_compress)
void sam_zip_qual (VBlock *vb, uint64_t vb_line_i, char **line_qual_data, uint32_t *line_qual_len, uint32_t maximum_len) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_qual_len  = MIN (maximum_len, dl->qual_data_len);
    
    if (!line_qual_data) return; // only lengths were requested

    *line_qual_data = ENT (char, vb->txt_data, dl->qual_data_start);

    // if QUAL is just "*" (i.e. unavailable) replace it by " " because '*' is a legal PHRED quality value that will confuse PIZ
    if (dl->qual_data_len == 1 && (*line_qual_data)[0] == '*') 
        *line_qual_data = " "; // pointer to static string

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    else if (flag.optimize_QUAL) 
        optimize_phred_quality_string (*line_qual_data, *line_qual_len);
}

// callback function for compress to get data of one line
void sam_zip_u2 (VBlock *vb, uint64_t vb_line_i, char **line_u2_data,  uint32_t *line_u2_len, uint32_t maximum_len) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    *line_u2_len = MIN (maximum_len, dl->u2_data_len);

    if (!line_u2_data) return; // only lengths were requested

    *line_u2_data = ENT (char, vb->txt_data, dl->u2_data_start);

    if (flag.optimize_QUAL)
        optimize_phred_quality_string (*line_u2_data, *line_u2_len);
}

// callback function for compress to get BD_BI data of one line: this is an
// interlaced line containing a character from BD followed by a character from BI - since these two fields are correlated
// note: if only one of BD or BI exists, the missing data in the interlaced string will be 0 (this should is not expected to ever happen)
void sam_zip_bd_bi (VBlock *vb_, uint64_t vb_line_i, 
                    char **line_data, uint32_t *line_len,  // out 
                    uint32_t maximum_len)
{
    VBlockSAM *vb = (VBlockSAM *)vb_;
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
    
    const char *bd = dl->bdbi_data_start[0] ? ENT (char, vb->txt_data, dl->bdbi_data_start[0]) : NULL;
    const char *bi = dl->bdbi_data_start[1] ? ENT (char, vb->txt_data, dl->bdbi_data_start[1]) : NULL;
    
    if (!bd && !bi) return; // no BD or BI on this line

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_len  = MIN (maximum_len, dl->seq_len * 2);

    if (!line_data) return; // only length was requested

    buf_alloc (vb, &vb->bd_bi_line, 0, dl->seq_len * 2, uint8_t, 2, "bd_bi_line");

    // calculate character-wise delta
    for (unsigned i=0; i < dl->seq_len; i++) {
        *ENT (uint8_t, vb->bd_bi_line, i*2    ) = bd ? bd[i] : 0;
        *ENT (uint8_t, vb->bd_bi_line, i*2 + 1) = bi ? bi[i] - (bd ? bd[i] : 0) : 0;
    }

    *line_data = FIRSTENT (char, vb->bd_bi_line);
}   

// called by zfile_compress_genozip_header to set FlagsGenozipHeader.dt_specific
bool sam_zip_dts_flag (void)
{
    return flag.reference == REF_INTERNAL;
}

// ----------------------
// Seg stuff
// ----------------------

void sam_zip_initialize (void)
{
    taxid_redirection_snip_len = sizeof (taxid_redirection_snip);
    seg_prepare_snip_other (SNIP_REDIRECTION, _SAM_TAXID, false, 0, 
                            taxid_redirection_snip, &taxid_redirection_snip_len);
}

void sam_seg_initialize (VBlock *vb)
{
    START_TIMER;

    // note: all numeric fields needs STORE_INT to be reconstructable to BAM (possibly already set)
    // via the translators set in the SAM_TOP2BAM Container
    CTX(SAM_SQBITMAP)->ltype    = LT_BITMAP;
    CTX(SAM_SQBITMAP)->local_always = true;
    CTX(SAM_TLEN)->flags.store  = STORE_INT;
    CTX(SAM_MAPQ)->flags.store  = STORE_INT;
    CTX(SAM_FLAG)->flags.store  = STORE_INT;
    CTX(SAM_RNAME)->flags.store = STORE_INDEX; // since v12
    CTX(SAM_POS)->flags.store   = STORE_INT;   // since v12
    CTX(SAM_PNEXT)->flags.store = STORE_INT;
    CTX(SAM_STRAND)->ltype      = LT_BITMAP;
    CTX(SAM_GPOS)->ltype        = LT_UINT32;
    CTX(SAM_GPOS)->flags.store  = STORE_INT;

    CTX(SAM_TOPLEVEL)->no_stons = true; // keep in b250 so it can be eliminated as all_the_same
    CTX(SAM_TOP2BAM)->no_stons  = true;
    CTX(SAM_TOP2FQ)->no_stons   = true;
    CTX(SAM_TOP2FQEX)->no_stons = true;

    Context *rname_ctx = CTX(SAM_RNAME);
    Context *rnext_ctx = CTX(SAM_RNEXT);

    rname_ctx->flags.store = rnext_ctx->flags.store = STORE_INDEX; // when reconstructing BAM, we output the word_index instead of the string
    rname_ctx->no_stons = rnext_ctx->no_stons = true;  // BAM reconstruction needs RNAME, RNEXT word indices. also needed for random access.

    // in --stats, consolidate stats 
    stats_set_consolidation (vb, SAM_SQBITMAP, 4, SAM_NONREF, SAM_NONREF_X, SAM_GPOS, SAM_STRAND);
    stats_set_consolidation (vb, OPTION_E2, 4, OPTION_2NONREF, OPTION_N2ONREFX, OPTION_2GPOS, OPTION_S2TRAND);

    codec_acgt_comp_init (vb);

    if (kraken_is_loaded) {
        CTX(SAM_TAXID)->flags.store    = STORE_INT;
        CTX(SAM_TAXID)->no_stons       = true; // must be no_stons the SEC_COUNTS data needs to mirror the dictionary words
        CTX(SAM_TAXID)->counts_section = true; 
    }

    COPY_TIMER (seg_initialize);
}

void sam_seg_finalize (VBlockP vb)
{
    // for qual data - select domqual compression if possible, or fallback 
    if (!codec_domq_comp_init (vb, SAM_QUAL, sam_zip_qual)) 
        CTX(SAM_QUAL)->ltype  = LT_SEQUENCE; 

    if (!codec_domq_comp_init (vb, OPTION_U2, sam_zip_u2)) 
        CTX(OPTION_U2)->ltype  = LT_SEQUENCE; 

    // top level snip - reconstruction as SAM
    SmallContainer top_level_sam = { 
        .repeats   = vb->lines.len,
        .is_toplevel = true,
        .callback  = true,
        .nitems_lo = 13,
        .items     = { { .dict_id = { _SAM_QNAME },    .seperator = "\t" },
                       { .dict_id = { _SAM_FLAG },     .seperator = "\t" },
                       { .dict_id = { _SAM_RNAME },    .seperator = "\t" },
                       { .dict_id = { _SAM_POS },      .seperator = "\t" },
                       { .dict_id = { _SAM_MAPQ },     .seperator = "\t" },
                       { .dict_id = { _SAM_CIGAR },    .seperator = "\t" },
                       { .dict_id = { _SAM_RNEXT },    .seperator = "\t" },
                       { .dict_id = { _SAM_PNEXT },    .seperator = "\t" },
                       { .dict_id = { _SAM_TLEN },     .seperator = "\t" },
                       { .dict_id = { _SAM_SQBITMAP }, .seperator = "\t" },
                       { .dict_id = { _SAM_QUAL },     .seperator = "\t" },
                       { .dict_id = { _SAM_OPTIONAL },                   },
                       { .dict_id = { _SAM_EOL },                        } }
    };
    container_seg (vb, CTX(SAM_TOPLEVEL), (ContainerP)&top_level_sam, 0, 0, 0);

    // top level snip - reconstruction as BAM
    // strategy: we start by reconstructing the variable-length fields first (after a prefix that sets them in place) 
    // - read_name, cigar, seq and qual - and then go back and fill in the fixed-location fields
    // Translation (a feature of Container): items reconstruct their data and then call a translation function to translate it to the desired format
    SmallContainer top_level_bam = { 
        .repeats   = vb->lines.len,
        .is_toplevel = true,
        .callback  = true,
        .nitems_lo = 13,
        .items     = { { .dict_id = { _SAM_RNAME },    .seperator = { CI_TRANS_NOR                    }, SAM2BAM_RNAME    }, // Translate - output word_index instead of string
                       { .dict_id = { _SAM_POS },      .seperator = { CI_TRANS_NOR | CI_TRANS_MOVE, 1 }, SAM2BAM_POS      }, // Translate - output little endian POS-1
                       { .dict_id = { _SAM_MAPQ },     .seperator = { CI_TRANS_NOR                    }, SAM2BAM_U8       }, // Translate - textual to binary number
                       { .dict_id = { _SAM_BAM_BIN },  .seperator = { CI_TRANS_NOR | CI_TRANS_MOVE, 2 }, SAM2BAM_LTEN_U16 }, // Translate - textual to binary number
                       { .dict_id = { _SAM_FLAG },     .seperator = { CI_TRANS_NOR | CI_TRANS_MOVE, 4 }, SAM2BAM_LTEN_U16 }, // Translate - textual to binary number
                       { .dict_id = { _SAM_RNEXT },    .seperator = { CI_TRANS_NOR                    }, SAM2BAM_RNAME    }, // Translate - output word_index instead of string
                       { .dict_id = { _SAM_PNEXT },    .seperator = { CI_TRANS_NOR | CI_TRANS_MOVE, 4 }, SAM2BAM_POS      }, // Translate - output little endian POS-1
                       { .dict_id = { _SAM_QNAME },    .seperator = { CI_TRANS_NUL                    }                   }, // normal 
                       { .dict_id = { _SAM_CIGAR },    .seperator = ""                                                    }, // handle in special reconstructor - translate textual to BAM CIGAR format + reconstruct l_read_name, n_cigar_op, l_seq
                       { .dict_id = { _SAM_TLEN },     .seperator = { CI_TRANS_NOR                    }, SAM2BAM_TLEN     }, // must be after CIGAR bc sam_piz_special_TLEN needs vb->seq_num
                       { .dict_id = { _SAM_SQBITMAP }, .seperator = "",                                  SAM2BAM_SEQ      }, // Translate - textual format to BAM format
                       { .dict_id = { _SAM_QUAL },     .seperator = "",                                  SAM2BAM_QUAL     }, // Translate - textual format to BAM format, set block_size
                       { .dict_id = { _SAM_OPTIONAL }, .seperator = { CI_TRANS_NOR                    }                   }, // up to v11, this had the SAM2BAM_OPTIONAL translator
                     }
    };

    // no container wide-prefix, skip l_name with a 4-character prefix
    static const char bam_line_prefix[] = { CON_PREFIX_SEP, // has prefix 
                                            CON_PREFIX_SEP, // end of (empty) container-wide prefix
                                            ' ',' ',' ',' ', CON_PREFIX_SEP }; // first item prefix - 4 spaces (place holder for block_size)

    container_seg (vb, CTX(SAM_TOP2BAM), (ContainerP)&top_level_bam, bam_line_prefix, sizeof(bam_line_prefix), 
                          IS_BAM ? sizeof (uint32_t) * vb->lines.len : 0); // if BAM, account for block_size

    // top level snip - reconstruction as FASTQ
    SmallContainer top_level_fastq = { 
        .repeats     = vb->lines.len,
        .is_toplevel = true,
        .callback    = true,  // drop non-primary chimeric reads and reads without QUAL data
        .nitems_lo   = 7,
        .items       = { { .dict_id = { _SAM_QNAME },    .seperator = "\n"                 }, 
                         { .dict_id = { _SAM_RNAME },    .seperator = { CI_TRANS_NOR }     }, // needed for reconstructing seq 
                         { .dict_id = { _SAM_POS },      .seperator = { CI_TRANS_NOR }     }, // needed for reconstructing seq
                         { .dict_id = { _SAM_FLAG },     .seperator = { CI_TRANS_NOR }, .translator = SAM2FASTQ_FLAG }, // need to know if seq is reverse complemented & if it is R2 ; reconstructs "1" for R1 and "2" for R2
                         { .dict_id = { _SAM_CIGAR },    .seperator = { CI_TRANS_NOR }     }, // needed for reconstructing seq
                         { .dict_id = { _SAM_SQBITMAP }, .seperator = "\n",             .translator =SAM2FASTQ_SEQ  }, 
                         { .dict_id = { _SAM_QUAL },     .seperator = "\n",             .translator =SAM2FASTQ_QUAL }, // also moves fastq "line" to R2 (paired file) if needed
                       }
    };

    // add a '@' to the description line, use a prefix to add the + line    
    static const char fastq_line_prefix[] = { CON_PREFIX_SEP, CON_PREFIX_SEP, '@', CON_PREFIX_SEP, CON_PREFIX_SEP, CON_PREFIX_SEP, 
                                              CON_PREFIX_SEP, CON_PREFIX_SEP, CON_PREFIX_SEP, '+', '\n', CON_PREFIX_SEP };

    container_seg (vb, CTX(SAM_TOP2FQ), (ContainerP)&top_level_fastq, fastq_line_prefix, sizeof(fastq_line_prefix), 0);

    // top level snip - reconstruction as FASTQ "extended" - with all the SAM fields in the description line
    SmallContainer top_level_fastq_ext = { 
        .repeats     = vb->lines.len,
        .is_toplevel = true,
        .callback    = true,  // drop non-primary chimeric reads and reads without QUAL data
        .nitems_lo   = 12,
        .items       = { { .dict_id = { _SAM_QNAME },    .seperator = "\t" }, 
                         { .dict_id = { _SAM_FLAG },     .seperator = "\t", .translator = SAM2FASTQ_FLAG }, // need to know if seq is reverse complemented & if it is R2 ; reconstructs "1" for R1 and "2" for R2
                         { .dict_id = { _SAM_RNAME },    .seperator = "\t" },
                         { .dict_id = { _SAM_POS },      .seperator = "\t" },
                         { .dict_id = { _SAM_MAPQ },     .seperator = "\t" },
                         { .dict_id = { _SAM_CIGAR },    .seperator = "\t" },
                         { .dict_id = { _SAM_RNEXT },    .seperator = "\t" },
                         { .dict_id = { _SAM_PNEXT },    .seperator = "\t" },
                         { .dict_id = { _SAM_TLEN },     .seperator = "\t" },
                         { .dict_id = { _SAM_OPTIONAL }, .seperator = "\n" },
                         { .dict_id = { _SAM_SQBITMAP }, .seperator = "\n", .translator =SAM2FASTQ_SEQ  }, 
                         { .dict_id = { _SAM_QUAL },     .seperator = "\n", .translator =SAM2FASTQ_QUAL }, // also moves fastq "line" to R2 (paired file) if needed
                       }
    };

    // add a '@' to the description line, use a prefix to add the + line    
    static const char fastq_ext_line_prefix[] = { 
        CON_PREFIX_SEP, CON_PREFIX_SEP, 
        '@', CON_PREFIX_SEP, 
        'F','L','A','G',':', CON_PREFIX_SEP, 
        'R','N','A','M','E',':', CON_PREFIX_SEP, 
        'P','O','S',':', CON_PREFIX_SEP, 
        'M','A','P','Q',':', CON_PREFIX_SEP, 
        'C','I','G','A','R',':', CON_PREFIX_SEP, 
        'R','N','E','X','T',':', CON_PREFIX_SEP, 
        'P','N','E','X','T',':', CON_PREFIX_SEP, 
        'T','L','E','N',':', CON_PREFIX_SEP, 
        CON_PREFIX_SEP,              // optional
        CON_PREFIX_SEP,              // SEQ
        '+', '\n', CON_PREFIX_SEP }; // QUAL

    container_seg (vb, CTX(SAM_TOP2FQEX), (ContainerP)&top_level_fastq_ext, 
                   fastq_ext_line_prefix, sizeof(fastq_ext_line_prefix), 0);
}

bool sam_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return 
        // typically small 
        dict_id.num == _SAM_TOPLEVEL  ||
        dict_id.num == _SAM_TOP2BAM   ||
        dict_id.num == _SAM_TOP2FQ    ||
        dict_id.num == _SAM_TOP2FQEX  ||
        dict_id.num == _SAM_FLAG      ||
        dict_id.num == _SAM_MAPQ      ||
        dict_id.num == _OPTION_MAPQ   ||
        dict_id.num == _SAM_QNAME     ||
        dict_id.num == _SAM_OPTIONAL  ||
        dict_id.num == _SAM_EOL       ||
        dict_id.num == _SAM_TAXID     ||

        // standard tags, see here: https://samtools.github.io/hts-specs/SAMtags.pdf
        dict_id.num == _OPTION_AM ||
        dict_id.num == _OPTION_AS ||
        dict_id.num == _OPTION_CM ||
        dict_id.num == _OPTION_LB ||
        dict_id.num == _OPTION_FI ||
        dict_id.num == _OPTION_H0 ||
        dict_id.num == _OPTION_H1 ||
        dict_id.num == _OPTION_H2 ||
        dict_id.num == _OPTION_MQ ||
        dict_id.num == _OPTION_NH ||
        dict_id.num == _OPTION_NM ||
        dict_id.num == _OPTION_OC ||
        dict_id.num == _OPTION_PG ||
        dict_id.num == _OPTION_PQ ||
        dict_id.num == _OPTION_PU ||
        dict_id.num == _OPTION_RG ||
        dict_id.num == _OPTION_SA ||
        dict_id.num == _OPTION_SM ||
        dict_id.num == _OPTION_TC ||
        dict_id.num == _OPTION_UQ ||
        
        // bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
        dict_id.num == _OPTION_X0 ||
        dict_id.num == _OPTION_X1 ||
        dict_id.num == _OPTION_XA ||
        dict_id.num == _OPTION_XN ||
        dict_id.num == _OPTION_XM ||
        dict_id.num == _OPTION_XO ||
        dict_id.num == _OPTION_XG ||
        dict_id.num == _OPTION_XS ||
        dict_id.num == _OPTION_XE ||
        
        // biobambam tags        
        dict_id.num == _OPTION_STRAND     ||            

        // typically smallish - a few thousands
        dict_id.num == _SAM_RNAME ||
        dict_id.num == _SAM_RNEXT ||
        dict_id.num == _OPTION_RNAME      ||
        dict_id.num == _OPTION_CC;
}

void sam_seg_verify_rname_pos (VBlock *vb, const char *p_into_txt, PosType this_pos)
{
    const Buffer *header_contigs = txtheader_get_contigs();

    if (flag.reference == REF_INTERNAL && (!header_contigs /* SQ-less SAM */ || !header_contigs->len /* SQ-less BAM */)) return;
    if (!this_pos) return; // unaligned
    
    PosType max_pos;
    if (header_contigs) {
        // since this SAM file has a header, all RNAMEs must be listed in it
        ASSSEG (vb->chrom_node_index < header_contigs->len, p_into_txt, "RNAME \"%.*s\" does not have an SQ record in the header", vb->chrom_name_len, vb->chrom_name);
        max_pos = ENT (RefContig, *header_contigs, vb->chrom_node_index)->max_pos;
    }
    else {
        WordIndex ref_contig = ref_contig_get_by_chrom (vb, gref, vb->chrom_node_index, vb->chrom_name, vb->chrom_name_len, &max_pos); // possibly an alt contig
        if (ref_contig == WORD_INDEX_NONE) {
            WARN_ONCE ("FYI: RNAME \"%.*s\" (and possibly others) is missing in the reference file. No harm.", 
                       vb->chrom_name_len, vb->chrom_name);
            return; // the sequence will be segged as unaligned
        }
    }

    ASSINP (this_pos <= max_pos, "%s: Error POS=%"PRId64" is beyond the size of \"%.*s\" which is %"PRId64". In vb=%u line_i=%"PRIu64" chrom_node_index=%d", 
            txt_name, this_pos, vb->chrom_name_len, vb->chrom_name, max_pos, vb->vblock_i, vb->line_i, vb->chrom_node_index);
}

// TLEN - 3 cases: 
// 1. if a non-zero value that is the negative of the previous line - a SNIP_DELTA & "-" (= value negation)
// 2. else, tlen>0 and pnext_pos_delta>0 and seq_len>0 tlen is stored as SNIP_SPECIAL & tlen-pnext_pos_delta-seq_len
// 3. else, stored as is
void sam_seg_tlen_field (VBlockSAM *vb, 
                         const char *tlen, unsigned tlen_len, // option 1
                         int64_t tlen_value, // option 2
                         PosType pnext_pos_delta, int32_t cigar_seq_len)
{
    Context *ctx = CTX(SAM_TLEN);

    if (tlen) { // option 1
        ASSSEG0 (tlen_len, tlen, "empty TLEN");

        bool is_int = str_get_int (tlen, tlen_len, &tlen_value); // note: tlen_value remains 0 if not a valid integer
        ASSSEG (is_int, tlen, "expecting TLEN to be an integer, but found \"%.*s\"", tlen_len, tlen);
    }

    unsigned add_bytes = IS_BAM ? sizeof (uint32_t) : tlen_len + 1;

    // case 1
    if (tlen_value && tlen_value == -ctx->last_value.i) {
        char snip_delta[2] = { SNIP_SELF_DELTA, '-' };
        seg_by_ctx (vb, snip_delta, 2, ctx, add_bytes);
    }
    // case 2:
    else if (tlen_value > 0 && pnext_pos_delta > 0 && cigar_seq_len > 0) {
        char tlen_by_calc[30] = { SNIP_SPECIAL, SAM_SPECIAL_TLEN };
        unsigned tlen_by_calc_len = str_int (tlen_value - pnext_pos_delta - (int64_t)cigar_seq_len, &tlen_by_calc[2]);
        seg_by_ctx (vb, tlen_by_calc, tlen_by_calc_len + 2, ctx, add_bytes);
    }
    // case default: add as is (option 1)
    else if (tlen)
        seg_by_ctx (vb, tlen, tlen_len, ctx, add_bytes);

    // case default: add as is (option 2)
    else {
        char snip[20];
        unsigned snip_len = str_int (tlen_value, snip);
        seg_by_ctx (vb, snip, snip_len, ctx, add_bytes);
    }

    ctx->last_value.i = tlen_value;
}

// Creates a bitmap from seq data - exactly one bit per base that is mapped to the reference (e.g. not for INSERT bases)
// - Normal SEQ: tracking CIGAR, we compare the sequence to the reference, indicating in a SAM_SQBITMAP whether this
//   base in the same as the reference or not. In case of REF_INTERNAL, if the base is not already in the reference, we add it.
//   bases that differ from the reference are stored in SAM_NONREF
// - Edge case: no POS (i.e. unaligned read) - we just store the sequence in SAM_NONREF
// - Edge case: no CIGAR (it is "*") - we just treat it as an M and compare to the reference
// - Edge case: no SEQ (it is "*") - we '*' in SAM_NONREF and indicate "different from reference" in the bitmap. We store a
//   single entry, regardless of the number of entries indicated by CIGAR
//
// Best explanation of CIGAR operations: https://davetang.org/wiki/tiki-index.php?page=SAM
void sam_seg_seq_field (VBlockSAM *vb, DidIType bitmap_did, const char *seq, uint32_t seq_len, PosType pos, const char *cigar, 
                        unsigned recursion_level, uint32_t level_0_seq_len, const char *level_0_cigar, unsigned add_bytes)
{
    START_TIMER;

    Context *bitmap_ctx = CTX(bitmap_did);
    Context *nonref_ctx = bitmap_ctx + 1;

    ASSERT (recursion_level < 500, "excess recursion recursion_level=%u seq_len=%u level_0_cigar=%s", // a large number of recursion calls can happen if CIGAR=9M10910N86M3274690N30M1S as observed with the STAR aligner https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
            recursion_level, seq_len, level_0_cigar);

    bitmap_ctx->txt_len += add_bytes; // byte counts for --show-sections

    // for unaligned lines, if we have refhash loaded, use the aligner instead of CIGAR-based segmenting
    if (!pos && flag.ref_use_aligner) {
        aligner_seg_seq ((VBlockP)vb, bitmap_ctx, seq, seq_len);
        goto align_nonref_local;
    }

    BitArray *bitmap = buf_get_bitarray (&bitmap_ctx->local);

    if (!recursion_level) {

        ASSERTW (seq_len < 1000000, "Warning: sam_seg_seq_field: seq_len=%u is suspeciously high and might indicate a bug", seq_len);

        buf_alloc (vb, &bitmap_ctx->local, roundup_bits2bytes64 (vb->ref_and_seq_consumed), vb->lines.len * (vb->ref_and_seq_consumed+5) / 8, uint8_t, CTX_GROWTH, "contexts->local"); 
        buf_extend_bits (&bitmap_ctx->local, vb->ref_and_seq_consumed);

        buf_alloc (vb, &nonref_ctx->local, seq_len + 3, vb->lines.len * seq_len / 4, uint8_t, CTX_GROWTH, "contexts->local"); 
    }

    // we can't compare to the reference if it is unaligned: we store the seqeuence in nonref without an indication in the bitmap
    if (!pos || (vb->chrom_name_len==1 && vb->chrom_name[0]=='*')) {
        buf_add (&nonref_ctx->local, seq, seq_len);
        goto align_nonref_local; 
    }

    if (seq[0] == '*') goto done; // we already handled a missing seq (SEQ="*") by adding a '-' to CIGAR - no data added here

    RefLock lock;
    Range *range = ref_seg_get_locked_range ((VBlockP)vb, gref, vb->chrom_node_index, vb->chrom_name, vb->chrom_name_len, pos, vb->ref_consumed, seq, &lock);

    // Cases where we don't consider the refernce and just copy the seq as-is
    if (!range || // 1. (denovo:) this contig defined in @SQ went beyond the maximum genome size of 4B and is thus ignored
                  // 2. (loaded:) case contig doesn't exist in the reference
        (cigar[0] == '*' && cigar[1] == 0)) { // case: there's no CIGAR. The sequence is not aligned to the reference even if we have RNAME and POS (and its length can exceed the reference contig)

        buf_add (&nonref_ctx->local, seq, seq_len);
        
        bit_array_clear_region (bitmap, bitmap_ctx->next_local, vb->ref_and_seq_consumed); // note: vb->ref_and_seq_consumed==0 if cigar="*"
        bitmap_ctx->next_local += vb->ref_and_seq_consumed;

        random_access_update_last_pos ((VBlockP)vb, DC_PRIMARY, pos + vb->ref_consumed - 1);

        if (range) ref_unlock (gref, lock);
        
        // note: in case of a missing range (which can be the first range in this seq, or a subsequent range), we zero the entire remaining bitmap.
        // this is because, absent a reference, we don't know how much ref is consumed by this missing range.
        goto align_nonref_local; 
    }    

    uint32_t pos_index     = pos - range->first_pos;
    uint32_t next_ref      = pos_index;
    const char *next_cigar = cigar;

    uint32_t i=0;
    int subcigar_len=0;
    char cigar_op;

    uint32_t ref_len_this_level = (flag.reference == REF_INTERNAL ? MIN (vb->ref_consumed, range->last_pos - pos + 1)
                                                                  : vb->ref_consumed); // possibly going around the end of the chromosome in case of a circular chromosome                                   

    uint32_t range_len = (range->last_pos - range->first_pos + 1);
    
    while (i < seq_len || next_ref < pos_index + ref_len_this_level) {

        ASSERT0 (i <= seq_len && next_ref <= pos_index + ref_len_this_level, "i or next_ref are out of range");

        subcigar_len = strtod (next_cigar, (char **)&next_cigar); // get number and advance next_cigar
        
        cigar_op = *(next_cigar++);

        if (cigar_op == 'M' || cigar_op == '=' || cigar_op == 'X') { // alignment match or sequence match or mismatch

            ASSERT (subcigar_len > 0 && subcigar_len <= (seq_len - i), 
                    "CIGAR %s implies seq_len longer than actual seq_len=%u (recursion_level=%u level0: cigar=%s seq_len=%u)", 
                    cigar, seq_len, recursion_level, level_0_cigar, level_0_seq_len);

            uint32_t bit_i = bitmap_ctx->next_local; // copy to automatic variable for performance
            uint32_t start_i = i;
            while (subcigar_len && next_ref < pos_index + ref_len_this_level) {

                // when we have an X we don't enter it into our internal ref, and we wait for a read with a = or M for that site,
                // as we assume that on average, more reads will have the reference base, leading to better compression
            
                bool normal_base = IS_NUCLEOTIDE (seq[i]);

                // circle around to beginning of chrom if out of range (can only happen with external reference, expected only with circular chromosomes) 
                uint32_t actual_next_ref = next_ref % range_len; 

                // case: we have not yet set a value for this site - we set it now. note: in ZIP, is_set means that the site
                // will be needed for pizzing. With REF_INTERNAL, this is equivalent to saying we have set the ref value for the site
                if (flag.reference == REF_INTERNAL && range && normal_base 
                    && !ref_is_nucleotide_set (range, actual_next_ref)) { 
                    
                    ref_set_nucleotide (range, actual_next_ref, seq[i]);
                    bit_array_set (&range->is_set, actual_next_ref); // we will need this ref to reconstruct
                    bit_array_set (bitmap, bit_i); bit_i++; // cannot increment inside the macro
                }

                // case our seq is identical to the reference at this site
                else if (range && normal_base && seq[i] == ref_base_by_idx (range, actual_next_ref)) {
                    bit_array_set (bitmap, bit_i); bit_i++;

                    if (flag.reference == REF_EXT_STORE)
                        bit_array_set (&range->is_set, actual_next_ref); // we will need this ref to reconstruct
                }
                
                // case: ref is set to a different value - we store our value in nonref_ctx
                else {
                    NEXTENT (char, nonref_ctx->local) = seq[i];
                    bit_array_clear (bitmap, bit_i); bit_i++;
                } 

                subcigar_len--;
                next_ref++;
                i++;
            }
            vb->ref_and_seq_consumed -= (i - start_i); // update in case a range in a subsequent recursion level is missing and we need to clear the bitmap
            bitmap_ctx->next_local = bit_i;
        } // end if 'M', '=', 'X'

        // for Insertion or Soft clipping - this SEQ segment doesn't align with the reference - we leave it as is 
        else if (cigar_op == 'I' || cigar_op == 'S') {

            ASSSEG (subcigar_len > 0 && subcigar_len <= (seq_len - i), seq,
                    "CIGAR %s implies seq_len longer than actual seq_len=%u", cigar, seq_len);

            buf_add (&nonref_ctx->local, &seq[i], subcigar_len);
            i += subcigar_len;
            subcigar_len = 0;
        }

        // for Deletion or Skipping - we move the next_ref ahead
        else if (cigar_op == 'D' || cigar_op == 'N') {
            unsigned ref_consumed = (flag.reference == REF_INTERNAL ? MIN (subcigar_len, range_len - next_ref)
                                                                    : subcigar_len);
            next_ref     += ref_consumed;
            subcigar_len -= ref_consumed;
        }

        // Hard clippping (H) or padding (P) - nothing much to do
        else if (cigar_op == 'H' || cigar_op == 'P') {
            subcigar_len = 0;
        }

        else {
            ASSSEG (cigar_op, vb->last_cigar, "Error in sam_seg_seq_field: End of CIGAR reached but we still have %u reference and %u sequence bases to consume"
                    "(cigar=%s pos=%"PRId64" recursion_level=%u level_0_cigar=%s level_0_seq_len=%u) (vb->ref_consumed=%d next_ref=%u pos_index=%u ref_len_this_level=%u subcigar_len=%u range=[%.*s %"PRId64"-%"PRId64"])",
                    pos_index + ref_len_this_level - next_ref, seq_len-i,   cigar, pos, recursion_level, level_0_cigar, level_0_seq_len,
                    vb->ref_consumed, next_ref, pos_index, ref_len_this_level, subcigar_len, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos);        

            ASSSEG (false, vb->last_cigar, "Invalid CIGAR op: '%c' (ASCII %u)", cigar_op, cigar_op);        
        }

        // case: we're at the end of the reference AND we want more of it
        if (next_ref == pos_index + ref_len_this_level && subcigar_len) break;
    }

    if (range) ref_unlock (gref, lock);       

    uint32_t this_seq_last_pos = pos + (next_ref - pos_index) - 1;

    // in REF_INTERNAL, the sequence can flow over to the next range as each range is 1M bases. this cannot happen
    // in REF_EXTERNAL as each range is the entire contig
    ASSERT (flag.reference == REF_INTERNAL || i == seq_len, "expecting i(%u) == seq_len(%u) pos=%"PRId64" range=[%.*s %"PRId64"-%"PRId64"] (cigar=%s recursion_level=%u level0: cigar=%s seq_len=%u)", 
            i, seq_len, pos, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos, cigar, recursion_level, level_0_cigar, level_0_seq_len);

    // case: we have reached the end of the current reference range, but we still have sequence left - 
    // call recursively with remaining sequence and next reference range 
    if (i < seq_len) {

        ASSSEG (this_seq_last_pos <= MAX_POS, cigar, "POS=%"PRId64" and the consumed reference implied by CIGAR=\"%s\", exceeding MAX_POS=%"PRId64
                " (next_ref=%u pos_index=%u ref_len_this_level=%u subcigar_len=%u range=[%.*s %"PRId64"-%"PRId64"])",
                pos, cigar, MAX_POS, next_ref, pos_index, ref_len_this_level, subcigar_len, 
                range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos);

        vb->ref_consumed -= ref_len_this_level;

        char updated_cigar[100];
        if (subcigar_len) sprintf (updated_cigar, "%u%c%s", subcigar_len, cigar_op, next_cigar);

        sam_seg_seq_field (vb, bitmap_did, seq + i, seq_len - i, range->last_pos + 1, subcigar_len ? updated_cigar : next_cigar, recursion_level + 1, level_0_seq_len, level_0_cigar, 0);
    }
    else { // update RA of the VB with last pos of this line as implied by the CIGAR string
        if (this_seq_last_pos <= range->last_pos) // always the case in INTERNAL and non-circular EXTERNAL 
            random_access_update_last_pos ((VBlockP)vb, DC_PRIMARY, this_seq_last_pos);

        else  // we circled back to the beginning for the chromosome - i.e. this VB RA is the entire chromosome
            random_access_update_to_entire_chrom ((VBlockP)vb, DC_PRIMARY, range->first_pos, range->last_pos);
    }
align_nonref_local: {
    // we align nonref_ctx->local to a 4-character boundary. this is because CODEC_ACGT squeezes every 4 characters into a byte,
    // before compressing it with LZMA. In sorted SAM, we want subsequent identical sequences to have the same byte alignment
    // so that LZMA can catch their identicality.
    uint64_t add_chars = (4 - (nonref_ctx->local.len & 3)) & 3;
    if (add_chars) buf_add (&nonref_ctx->local, "AAA", add_chars); // add 1 to 3 As
}
done:
    COPY_TIMER (sam_seg_seq_field);
}

// returns length of string ending with separator, or -1 if separator was not found
static inline int sam_seg_get_next_subitem (const char *str, int str_len, char separator)
{
    for (int i=0; i < str_len; i++) {
        if (str[i] == separator) return i;
        if (str[i] == ',' || str[i] == ';') return -1; // wrong separator encountered
    }
    return -1;
}

#define DO_SSF(ssf,sep) \
        ssf = &field[i]; \
        ssf##_len = sam_seg_get_next_subitem (&field[i], field_len-i, sep); \
        if (ssf##_len == -1) goto error; /* bad format */ \
        i += ssf##_len + 1; /* skip snip and separator */        

#define DEC_SSF(ssf) const char *ssf; \
                     int ssf##_len;   \
                     Context *ssf##_ctx = ctx_get_ctx (vb, con.items[item_i++].dict_id); \
                     ssf##_ctx->st_did_i = ctx->did_i; 

static void sam_seg_SA_or_OA_field (VBlockSAM *vb, DictId subfield_dict_id, 
                                    const char *field, unsigned field_len, const char *field_name)
{
    // OA and SA format is: (rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+ . in OA - NM is optional (but its , is not)
    // Example SA:Z:chr13,52863337,-,56S25M70S,0,0;chr6,145915118,+,97S24M30S,0,0;chr18,64524943,-,13S22M116S,0,0;chr7,56198174,-,20M131S,0,0;chr7,87594501,+,34S20M97S,0,0;chr4,12193416,+,58S19M74S,0,0;
    // See: https://samtools.github.io/hts-specs/SAMtags.pdf
    // note: even though SA, OA, XA contain similar fields amongst each other and similar to the primary fields,
    // the values of subsequent lines tend to be similar for each one of them seperately, so we maintain separate contexts
    #define CONTAINER_SA_OA(s) {   \
        .repeats     = 0,          \
        .nitems_lo   = 6,          \
        .repsep      = {0,0},      \
        .items       = { { .dict_id = {.id=s "0ARNAME" }, .seperator = {','} },  \
                         { .dict_id = {.id=s "1APOS"   }, .seperator = {','} },  \
                         { .dict_id = {.id=s "2ASTRAN" }, .seperator = {','} },  \
                         { .dict_id = {.id=s "3ACIGAR" }, .seperator = {','} },  \
                         { .dict_id = {.id=s "4AMAPQ"  }, .seperator = {','} },  \
                         { .dict_id = {.id=s "5ANM"    }, .seperator = {';'} } } \
    }
    static const SmallContainer container_SA = CONTAINER_SA_OA("S"), container_OA = CONTAINER_SA_OA("O");

    SmallContainer con = (subfield_dict_id.num == _OPTION_SA) ? container_SA : container_OA; // make a copy
    Context *ctx = ctx_get_ctx (vb, subfield_dict_id);

    unsigned item_i=0;
    DEC_SSF(rname); DEC_SSF(pos); DEC_SSF(strand); DEC_SSF(cigar); DEC_SSF(mapq); DEC_SSF(nm); 

    for (uint32_t i=0; i < field_len; con.repeats++) {

        ASSSEG (con.repeats <= CONTAINER_MAX_REPEATS, field, "exceeded maximum repeats allowed (%u) while parsing %s",
                CONTAINER_MAX_REPEATS, dis_dict_id (subfield_dict_id).s);

        DO_SSF (rname,  ','); // these also do sanity checks
        DO_SSF (pos,    ','); 
        DO_SSF (strand, ','); 
        DO_SSF (cigar,  ','); 
        DO_SSF (mapq,   ','); 
        DO_SSF (nm,     ';'); 

        // sanity checks before adding to any dictionary
        if (strand_len != 1 || (strand[0] != '+' && strand[0] != '-')) goto error; // invalid format
        
        seg_by_ctx (vb, rname,  rname_len,  rname_ctx,  1 + rname_len);
        seg_by_ctx (vb, strand, strand_len, strand_ctx, 1 + strand_len);
        seg_by_ctx (vb, cigar,  cigar_len,  cigar_ctx,  1 + cigar_len);
        seg_by_ctx (vb, mapq,   mapq_len,   mapq_ctx,   1 + mapq_len);
        seg_by_ctx (vb, nm,     nm_len,     nm_ctx,     1 + nm_len);
        
        Context *pos_ctx = ctx_get_ctx (vb, con.items[1].dict_id);
        seg_pos_field ((VBlockP)vb, pos_ctx->did_i, pos_ctx->did_i, 0, 0, pos, pos_len, 0, 1 + pos_len);
    }

    container_seg_by_dict_id (vb, subfield_dict_id, (ContainerP)&con, 1 /* 1 for \t in SAM and \0 in BAM */);
    
    return;

error:
    // if the error occurred on on the first repeat - this file probably has a different
    // format - we just store as a normal subfield
    // if it occurred on the 2nd+ subfield, after the 1st one was fine - we reject the file
    ASSSEG (!con.repeats, field, "Invalid format in repeat #%u of field %s. snip: %.*s",
            con.repeats+1, dis_dict_id (subfield_dict_id).s, field_len, field);

    seg_by_dict_id (vb, field, field_len, subfield_dict_id, field_len + 1 /* 1 for \t in SAM and \0 in BAM */); 
}

static void sam_seg_XA_field (VBlockSAM *vb, const char *field, unsigned field_len)
{
    // XA format is: (chr,pos,CIGAR,NM;)*  pos starts with +- which is strand
    // Example XA:Z:chr9,-60942781,150M,0;chr9,-42212061,150M,0;chr9,-61218415,150M,0;chr9,+66963977,150M,1;
    // See: http://bio-bwa.sourceforge.net/bwa.shtml
    static const SmallContainer container_XA = {
        .repeats     = 0, 
        .nitems_lo   = 5, 
        .repsep      = {0,0},
        .items       = { { .dict_id = {.id="X0ARNAME" }, .seperator = {','} }, // note: optional fields are DTYPE_2, in which short ids are left as-is, so we can skip dict_id_make
                         { .dict_id = {.id="X1ASTRAN" }, .seperator = { 0 } },
                         { .dict_id = {.id="X2APOS"   }, .seperator = {','} },
                         { .dict_id = {.id="X3ACIGAR" }, .seperator = {','} }, // we don't mix the primary as the primary has a SNIP_SPECIAL
                         { .dict_id = {.id="X4ANM"    }, .seperator = {';'} } }     
    };

    Context *ctx = CTX(OPTION_XA);

    SmallContainer con = container_XA;

    unsigned item_i=0;
    DEC_SSF(rname); DEC_SSF(strand); DEC_SSF(pos); DEC_SSF(cigar); DEC_SSF(nm); 

    for (uint32_t i=0; i < field_len; con.repeats++) {

        ASSSEG (con.repeats <= CONTAINER_MAX_REPEATS, field, "exceeded maximum repeats allowed (%u) while parsing XA",
                CONTAINER_MAX_REPEATS);

        DO_SSF (rname,  ','); 
        DO_SSF (pos,    ','); // includes strand
        DO_SSF (cigar,  ','); 
        DO_SSF (nm,     ';'); 

        // split the pos string, eg "-10000" to strand "-" and pos "10000"
        if (pos_len < 2 || (pos[0] != '+' && pos[0] != '-')) goto error; // invalid format - expecting pos to begin with the strand
        strand = pos++;
        pos_len--;
        strand_len = 1;

        seg_by_ctx (vb, rname,  rname_len,  rname_ctx,  1 + rname_len);
        seg_by_ctx (vb, strand, strand_len, strand_ctx, strand_len); // strand is first character of pos - no separator
        seg_by_ctx (vb, cigar,  cigar_len,  cigar_ctx,  1 + cigar_len);
        seg_by_ctx (vb, nm,     nm_len,     nm_ctx,     1 + nm_len);
        
        seg_integer_or_not ((VBlockP)vb, pos_ctx, pos, pos_len, 1+pos_len);
    }

    container_seg_by_dict_id (vb, _OPTION_XA, (ContainerP)&con, 1 /* 1 for \t in SAM and \0 in BAM */);
    return;

error:
    // if the error occurred on on the first repeat - this file probably has a different
    // format - we just store as a normal subfield
    // if it occurred on the 2nd+ subfield, after the 1st one was fine - we reject the file
    ASSSEG (!con.repeats, field, "Invalid format in repeat #%u of field XA. snip: %.*s", con.repeats+1, field_len, field);

    seg_by_dict_id (vb, field, field_len, _OPTION_XA, field_len + 1 /* 1 for \t in SAM and \0 in BAM */); 
}

uint32_t sam_seg_get_seq_len_by_MD_field (const char *md_str, unsigned md_str_len)
{
    uint32_t result=0, curr_num=0;

    for (unsigned i=0; i < md_str_len; i++) {   
        if (IS_DIGIT (md_str[i])) 
            curr_num = curr_num * 10 + (md_str[i] - '0');

        else {
            result += curr_num + 1; // number terminates here + one character
            curr_num = 0;
        }
    }

    result += curr_num; // in case the string ends with a number

    return result;
}

// in the case where sequence length as calculated from the MD is the same as that calculated
// from the CIGAR/SEQ/QUAL (note: this is required by the SAM spec but nevertheless genozip doesn't require it):
// MD is shortened to replace the last number with a *, since it can be calculated knowing the length. The result is that
// multiple MD values collapse to one, e.g. "MD:Z:119C30" and "MD:Z:119C31" both become "MD:Z:119C*" hence improving compression.
// In the case where the MD is simply a number "151" and drop it altogether and keep just an empty string.
static inline bool sam_seg_get_shortened_MD (const char *md_str, unsigned md_str_len, uint32_t seq_len,
                                             char *new_md_str, unsigned *new_md_str_len)
{
    uint32_t seq_len_by_md = sam_seg_get_seq_len_by_MD_field (md_str, md_str_len);

    if (seq_len_by_md != seq_len) return false;  // MD string doesn't comply with SAM spec and is therefore not changed
    
    // case - MD ends with a number eg "119C31" - we replace it with prefix+"119C". if its all digits then just prefix
    if (IS_DIGIT (md_str[md_str_len-1])) {

        int i=md_str_len-1; for (; i>=0; i--)
            if (!IS_DIGIT (md_str[i])) break;

        new_md_str[0] = SNIP_SPECIAL;
        new_md_str[1] = SAM_SPECIAL_MD;
        if (i >= 0) memcpy (&new_md_str[2], md_str, i+1);
        
        *new_md_str_len = i+3;
        return true;
    }

    return false; // MD doesn't end with a number and is hence unchanged (this normally doesn't occur as the MD would finish with 0)
}

// AS and XS are values (at least as set by BWA) at most the seq_len, and AS is often equal to it. we modify
// it to be new_value=(value-seq_len) 
static inline void sam_seg_AS_field (VBlockSAM *vb, ZipDataLineSAM *dl, DictId dict_id, 
                                     const char *snip, unsigned snip_len, unsigned add_bytes)
{
    bool positive_delta = true;

    // verify that its a unsigned number
    for (unsigned i=0; i < snip_len; i++)
        if (!IS_DIGIT (snip[i])) positive_delta = false;

    int32_t as;
    if (positive_delta) {
        as = atoi (snip); // type i is signed 32 bit by SAM specification
        if (dl->seq_len < as) positive_delta=false;
    }

    // if possible, store a special snip with the positive delta
    if (positive_delta) {
        char new_snip[20] = { SNIP_SPECIAL, SAM_SPECIAL_AS };
        unsigned delta_len = str_int (dl->seq_len-as, &new_snip[2]);

        seg_by_dict_id (vb, new_snip, delta_len+2, dict_id, add_bytes); 
    }

    // not possible - just store unmodified
    else
        seg_by_dict_id (vb, snip, snip_len, dict_id, add_bytes); 
}

// mc:i: (output of bamsormadup and other biobambam tools - mc in small letters) 
// appears to be a pos value usually close to PNEXT, but it is -1 is POS=PNEXT.
static inline void sam_seg_mc_field (VBlockSAM *vb, DictId dict_id, 
                                     const char *snip, unsigned snip_len, unsigned add_bytes)
{
    uint8_t mc_did_i = ctx_get_ctx (vb, dict_id)->did_i;
    
    // if snip is "-1", store as simple snip
    if (snip_len == 2 && snip[0] == '-' && snip[1] == '1')
        seg_by_did_i (vb, snip, snip_len, mc_did_i, add_bytes);
    
    // delta vs PNEXT
    else
        seg_pos_field ((VBlockP)vb, mc_did_i, SAM_PNEXT, SPF_BAD_SNIPS_TOO, 0, snip, snip_len, 0, add_bytes);
}

// optimization for Ion Torrent flow signal (ZM) - negative values become zero, positives are rounded to the nearest 10
static void sam_optimize_ZM (const char **snip, unsigned *snip_len, char *new_str)
{
    char *after;
    int number = strtoul (*snip, &after, 10);

    if ((unsigned)(after - *snip) > 0) {
        if (number >= 0) number = ((number + 5) / 10) * 10;
        else             number = 0;

        *snip_len = str_int (number, new_str);
        *snip = new_str;
    }    
}

static inline TranslatorId optional_field_translator (char type)
{
    switch (type) {
        case 'c' : return SAM2BAM_I8;
        case 'C' : return SAM2BAM_U8;
        case 's' : return SAM2BAM_LTEN_I16;
        case 'S' : return SAM2BAM_LTEN_U16;
        case 'i' : return SAM2BAM_LTEN_I32;
        case 'I' : return SAM2BAM_LTEN_U32;
        case 'f' : return IS_BAM ? 0 : SAM2BAM_FLOAT; // when reconstucting BAM->BAM we use SPECIAL_FLOAT rather than a translator
        default  : return 0;
    }
}

static inline unsigned sam_seg_optional_add_bytes (char type, unsigned value_len, bool is_bam)
{
    if (is_bam)
        switch (type) {
            case 'c': case 'C': case 'A': return 1;
            case 's': case 'S':           return 2;
            case 'i': case 'I': case 'f': return 4;
            case 'Z': case 'H':           return value_len + 1; // +1 for \0
            default : return 0;
        }
    else // SAM
        return value_len + 1; // +1 for \t
}

static inline char sam_seg_bam_type_to_sam_type (char type)
{
    return (type=='c' || type=='C' || type=='s' || type=='S' || type=='I') ? 'i' : type;
}

// an array - all elements go into a single item context, multiple repeats
static void sam_seg_array_field (VBlock *vb, DictId dict_id, const char *value, unsigned value_len)
{   
    // get optimization function, if there is one
    SegOptimize optimize = NULL;
    if (flag.optimize_ZM && dict_id.num == _OPTION_ZM && value_len > 3 && value[0] == 's')  // XM:B:s,
        optimize = sam_optimize_ZM;

    // prepare array container - a single item, with number of repeats of array element. array type is stored as a prefix
    Context *container_ctx     = ctx_get_ctx (vb, dict_id);

    SmallContainer con = { .nitems_lo = 2, 
                           .drop_final_item_sep_of_final_repeat = true, // TODO - get rid of this flag and move to making the seperators to be repeat seperators as they should have been, using drop_final_repeat_sep and obsoleting this flag 
                           .repsep    = {0,0}, 
                           .items     = { { .translator = SAM2BAM_ARRAY_SELF  },  // item[0] is translator-only item - to translate the Container itself in case of reconstructing BAM 
                                          { .seperator  = {0, ','}            } } // item[1] is actual array item
                         };
    
    char prefixes[] = { CON_PREFIX_SEP, value[0], ',', CON_PREFIX_SEP }; // prefix contains type eg "i,"
    
    const char *str = value + 2;      // remove type and comma
    int str_len = (int)value_len - 2; // must be int, not unsigned, for the for loop

    // prepare context where array elements will go in
    char arr_dict_id_str[8]   = "XX_ARRAY";
    arr_dict_id_str[0]        = FLIP_CASE (dict_id.id[0]);
    arr_dict_id_str[1]        = FLIP_CASE (dict_id.id[1]);
    con.items[1].dict_id      = dict_id_make (arr_dict_id_str, 8, DTYPE_SAM_OPTIONAL);
    con.items[1].translator   = optional_field_translator ((uint8_t)value[0]); // instructions on how to transform array items if reconstructing as BAM (value[0] is the subtype of the array)
    con.items[1].seperator[0] = optional_sep_by_type[IS_BAM][(uint8_t)value[0]];
    
    Context *element_ctx      = ctx_get_ctx (vb, con.items[1].dict_id);
    element_ctx->st_did_i     = container_ctx->did_i;
    element_ctx->flags.store  = optional_field_store_flag[(uint8_t)value[0]];

    for (con.repeats=0; con.repeats < CONTAINER_MAX_REPEATS && str_len > 0; con.repeats++) { // str_len will be -1 after last number

        const char *snip = str;
        for (; str_len && *str != ','; str++, str_len--) {};

        unsigned number_len = (unsigned)(str - snip);
        unsigned snip_len   = number_len; // might be changed by optimize
             
        char new_number_str[30];
        if (optimize && snip_len < 25)
            optimize (&snip, &snip_len, new_number_str);

        seg_by_ctx (vb, snip, snip_len, element_ctx, IS_BAM ? 0 : number_len+1);
        
        str_len--; // skip comma
        str++;
    }

    ASSSEG (con.repeats < CONTAINER_MAX_REPEATS, value, "array has too many elements, more than %u", CONTAINER_MAX_REPEATS);

    // add bytes here in case of BAM - all to main field
    unsigned container_add_bytes=0;
    if (IS_BAM) {
        unsigned add_bytes_per_repeat = sam_seg_optional_add_bytes (value[0], 0, true);
        container_add_bytes = add_bytes_per_repeat * con.repeats + 4 /* count */ + 1 /* type */ ;
    }
    else 
        container_add_bytes = 2; // type - eg "i,"

    container_seg (vb, container_ctx, (ContainerP)&con, prefixes, sizeof(prefixes), container_add_bytes);
}

// process an optional subfield, that looks something like MX:Z:abcdefg. We use "MX" for the field name, and
// the data is abcdefg. The full name "MX:Z:" is stored as part of the OPTIONAL dictionary entry
static DictId sam_seg_optional_field (VBlockSAM *vb, ZipDataLineSAM *dl, bool is_bam, 
                                      const char *tag, char bam_type, const char *value, unsigned value_len)
{
    char sam_type = sam_seg_bam_type_to_sam_type (bam_type);
    char dict_name[4] = { tag[0], tag[1], ':', sam_type };
    DictId dict_id = dict_id_make (dict_name, 4, DTYPE_SAM_OPTIONAL);

    unsigned add_bytes = sam_seg_optional_add_bytes (bam_type, value_len, is_bam);

    if (dict_id.num == _OPTION_SA || dict_id.num == _OPTION_OA)
        sam_seg_SA_or_OA_field (vb, dict_id, value, value_len, dict_id.num == _OPTION_SA ? "SA" : "OA");

    else if (dict_id.num == _OPTION_XA) 
        sam_seg_XA_field (vb, value, value_len);

    // fields containing CIGAR format data - aliases of _OPTION_CIGAR (not the main CIGAR field that all snips have SNIP_SPECIAL)
    // MC: "Mate Cigar", added by eg https://manpages.debian.org/unstable/biobambam2/bamsort.1.en.html  
    else if (dict_id.num == _OPTION_MC || dict_id.num == _OPTION_OC) 
        seg_by_dict_id (vb, value, value_len, _OPTION_CIGAR, add_bytes); 

    // MD's logical length is normally the same as seq_len, we use this to optimize it.
    // In the common case that it is just a number equal the seq_len, we replace it with an empty string.
    else if (dict_id.num == _OPTION_MD) {
        // if MD value can be derived from the seq_len, we don't need to store - store just an empty string

#define MAX_SAM_MD_LEN 1000 // maximum length of MD that is shortened.
        char new_md[MAX_SAM_MD_LEN];
        unsigned new_md_len = 0;
        bool md_is_special  = (value_len-2 <= MAX_SAM_MD_LEN);

        if (md_is_special) 
            md_is_special = sam_seg_get_shortened_MD (value, value_len, dl->seq_len, new_md, &new_md_len);

        // not sure which of these two is better....
        seg_by_dict_id (vb,                                 
                        md_is_special ? new_md : value, 
                        md_is_special ? new_md_len : value_len,
                        dict_id, add_bytes);
    }

    // BD and BI set by older versions of GATK's BQSR is expected to be seq_len (seen empircally, documentation is lacking)
    else if ((dict_id.num == _OPTION_BD || dict_id.num == _OPTION_BI) && value_len == dl->seq_len) {
        
        bool is_bi = (dict_id.num == _OPTION_BI);
        dl->bdbi_data_start[is_bi] = value - vb->txt_data.data;

        Context *ctx = CTX(OPTION_BD_BI);
        ctx->txt_len += add_bytes; 
        ctx->ltype   = LT_SEQUENCE;

        if (!dl->bdbi_data_start[!is_bi]) // the first of BD and BI increments local.len, so it is incremented even if just one of BD/BI appears
            ctx->local.len += value_len * 2;

        // we can't use local for singletons in BD or BI as next_local is used by sam_piz_special_BD_BI to point into BD_BI
        Context *this_ctx  = ctx_get_ctx (vb, dict_id);
        this_ctx->no_stons = true; 
        this_ctx->st_did_i = ctx->did_i; 

        const char special_snip[2] = { SNIP_SPECIAL, SAM_SPECIAL_BDBI };
        seg_by_dict_id (vb, special_snip, 2, dict_id, 0);
    }

    // AS is a value (at least as set by BWA) at most the seq_len, and often equal to it. we modify
    // it to be new_AS=(AS-seq_len) 
    else if (dict_id.num == _OPTION_AS) 
        sam_seg_AS_field (vb, dl, dict_id, value, value_len, add_bytes);
    
    // mc:i: (output of bamsormadup and other biobambam tools - mc in small letters) 
    // appears to be a pos value usually close to PNEXT, but it is -1 is POS=PNEXT.
    else if (dict_id.num == _OPTION_mc) 
        sam_seg_mc_field (vb, dict_id, value, value_len, add_bytes);

    // TX:i: - we seg this as a primary field SAM_TAX_ID
    else if (dict_id.num == _OPTION_TX) 
        seg_by_dict_id (vb, taxid_redirection_snip, taxid_redirection_snip_len, dict_id, add_bytes); 

    // E2 - SEQ data. Currently broken. To do: fix.
/*    else if (dict_id.num == _OPTION_E2) {
        ASSSEG0 (dl->seq_len, value, "E2 tag without a SEQ"); 
        ASSINP (value_len == dl->seq_len, 
                "Error in %s: Expecting E2 data to be of length %u as indicated by CIGAR, but it is %u. E2=%.*s",
                txt_name, dl->seq_len, value_len, value_len, value);

        PosType this_pos = vb->last_int(SAM_POS);

        sam_seg_seq_field (vb, OPTION_E2, (char *)value, value_len, this_pos, vb->last_cigar, 0, value_len, // remove const bc SEQ data is actually going to be modified
                           vb->last_cigar, add_bytes); 
    }
*/
    // U2 - QUAL data (note: U2 doesn't have a context - it shares with QUAL)
    else if (dict_id.num == _OPTION_U2) {
        ASSSEG0 (dl->seq_len, value, "U2 tag without a SEQ"); 
        ASSINP (value_len == dl->seq_len, 
                "Error in %s: Expecting U2 data to be of length %u as indicated by CIGAR, but it is %u. E2=%.*s",
                txt_name, dl->seq_len, value_len, value_len, value);

        dl->u2_data_start = value - vb->txt_data.data;
        dl->u2_data_len   = value_len;
        CTX(OPTION_U2)->txt_len   += add_bytes;
        CTX(OPTION_U2)->local.len += value_len;
    }

    // Numeric array array
    else if (bam_type == 'B') 
        sam_seg_array_field ((VBlockP)vb, dict_id, value, value_len);

    // All other subfields - normal snips in their own dictionary
    else        
        seg_by_dict_id (vb, value, value_len, dict_id, add_bytes); 

    // integer and float fields need to be STORE_INT/FLOAT to be reconstructable as BAM
    if (optional_field_store_flag[(uint8_t)sam_type]) {
        Context *ctx;
        if ((ctx = ECTX (dict_id)))
            ctx->flags.store = optional_field_store_flag[(uint8_t)sam_type];
    }
 
    return dict_id;
}

const char *sam_get_one_optional (VBlockSAM *vb, const char *next_field, int32_t len, char *separator_p, bool *has_13, 
                                  const char **tag, char *type, const char **value, unsigned *value_len) // out
{
    unsigned field_len;
    const char *field_start;

    char separator;
    GET_MAYBE_LAST_ITEM (subfield); 

    ASSSEG0 (field_len, field_start, "line invalidly ends with a tab");

    ASSSEG (field_len >= 6 && field_start[2] == ':' && field_start[4] == ':', field_start, "invalid optional field format: %.*s",
            field_len, field_start);

    *tag         = field_start;
    *type        = field_start[3];
    *value       = field_start + 5;
    *value_len   = field_len - 5;
    *separator_p = separator;

    return next_field;
}

static const char *sam_seg_get_kraken (VBlockSAM *vb, const char *next_field, bool *has_kraken, 
                                       char *taxid_str, const char **tag, char *type,  // out
                                       const char **value, unsigned *value_len, // out
                                       bool is_bam)
{
    *tag        = "TX"; // genozip introduced tag (=taxid)
    *type       = 'i';
    *value      = taxid_str;
    *has_kraken = false;
    *value_len  = kraken_seg_taxid_do ((VBlockP)vb, SAM_TAXID, last_txt (vb, SAM_QNAME), vb->last_txt_len (SAM_QNAME),
                                       taxid_str, true);

    vb->recon_size += is_bam ? 7 : (*value_len + 6); // txt modified

    return next_field; // unmodified
}

const char *sam_seg_optional_all (VBlockSAM *vb, ZipDataLineSAM *dl, const char *next_field,
                                  int32_t len, bool *has_13, char separator, // sam only
                                  const char *after_field) // bam only 
{
    const bool is_bam = IS_BAM;
    Container con = { .repeats=1 };
    char prefixes[MAX_FIELDS * 6 + 3]; // each name is 5 characters per SAM specification, eg "MC:Z:" followed by CON_PREFIX_SEP ; +3 for the initial CON_PREFIX_SEP
    prefixes[0] = prefixes[1] = prefixes[2] = CON_PREFIX_SEP; // initial CON_PREFIX_SEP follow by seperator of empty Container-wide prefix followed by seperator for empty prefix for translator-only item[0]
    unsigned prefixes_len=3;
    const char *value, *tag;
    char type;
    unsigned value_len;
    char added_value[20];

    // item[0] is translator-only item - to translate the Container itself in case of reconstructing BAM 
    con.items[con_nitems(con)] = (ContainerItem){ .translator = SAM2BAM_OPTIONAL_SELF }; 
    con_inc_nitems (con);
    
    bool has_kraken = kraken_is_loaded;

    while ((is_bam ? (next_field < after_field) : (separator != '\n'))
         || has_kraken) {

        next_field = has_kraken ? sam_seg_get_kraken   (vb, next_field, &has_kraken, added_value, &tag, &type, &value, &value_len, is_bam)  
                   : is_bam     ? bam_get_one_optional (vb, next_field,                           &tag, &type, &value, &value_len)
                   :              sam_get_one_optional (vb, next_field, len, &separator, has_13,  &tag, &type, &value, &value_len);

        con.items[con_nitems(con)] = (ContainerItem) {
            .dict_id    = sam_seg_optional_field (vb, dl, is_bam, tag, type, value, value_len),
            .translator = optional_field_translator ((uint8_t)type), // how to transform the field if reconstructing to BAM
            .seperator  = { optional_sep_by_type[is_bam][(uint8_t)type], '\t' },
        };
        con_inc_nitems (con);

        ASSSEG (con_nitems(con) <= MAX_FIELDS, value, "too many optional fields, limit is %u", MAX_FIELDS);

        // in the optional field prefix (unlike array type), all integer types become 'i'.
        char prefix_type = sam_seg_bam_type_to_sam_type (type);

        char prefix[6] = { tag[0], tag[1], ':', prefix_type, ':', CON_PREFIX_SEP}; 
        memcpy (&prefixes[prefixes_len], prefix, 6);
        prefixes_len += 6;

        if (is_bam) buf_free (&vb->textual_opt);
    }

    uint32_t num_items = con_nitems(con);
    if (num_items > 1) { // we have Optional fields, not just the translator item
        if (con.items[num_items-1].seperator[0] & 0x80) // is a flag
            con.items[num_items-1].seperator[0] &= ~(CI_NATIVE_NEXT & ~(uint8_t)0x80); // last Optional field has no tab
        con.items[num_items-1].seperator[1] = 0;
        container_seg (vb, CTX(SAM_OPTIONAL), &con, prefixes, prefixes_len, (is_bam ? 3 : 5) * (num_items-1)); // account for : SAM: "MX:i:" BAM: "MXi"
    }
    else
        // NULL means MISSING Container item (of the toplevel container) - will cause container_reconstruct_do of 
        // the toplevel container to delete of previous separator (\t)
        container_seg (vb, CTX(SAM_OPTIONAL), 0, 0, 0, 0); 

    return next_field;        
}

static void sam_seg_cigar_field (VBlockSAM *vb, ZipDataLineSAM *dl, unsigned last_cigar_len,
                                 const char *seq,  uint32_t seq_data_len, 
                                 const char *qual, uint32_t qual_data_len)
{
    bool qual_is_available = (qual_data_len != 1 || *qual != '*');
    bool seq_is_available  = (seq_data_len  != 1 || *seq  != '*');

    ASSSEG (!(seq_is_available && *seq=='*'), seq, "seq=%.*s (seq_len=%u), but expecting a missing seq to be \"*\" only (1 character)", seq_data_len, seq, seq_data_len);

    char cigar_snip[last_cigar_len + 50];
    cigar_snip[0] = SNIP_SPECIAL;
    cigar_snip[1] = SAM_SPECIAL_CIGAR;
    unsigned cigar_snip_len=2;

    // case: SEQ is "*" - we add a '-' to the CIGAR
    if (!seq_is_available) cigar_snip[cigar_snip_len++] = '-';

    // case: CIGAR is "*" - we get the dl->seq_len directly from SEQ or QUAL, and add the length to CIGAR eg "151*"
    if (!dl->seq_len) { // CIGAR is not available
        ASSSEG (!seq_data_len || !qual_is_available || seq_data_len==dl->qual_data_len, seq,
                "Bad line: SEQ length is %u, QUAL length is %u, unexpectedly differ. SEQ=%.*s QUAL=%.*s", 
                seq_data_len, dl->qual_data_len, seq_data_len, seq, dl->qual_data_len, qual);    

        dl->seq_len = MAX (seq_data_len, dl->qual_data_len); // one or both might be not available and hence =1

        cigar_snip_len += str_int (dl->seq_len, &cigar_snip[cigar_snip_len]);
    } 
    else { // CIGAR is available - just check the seq and qual lengths
        ASSSEG (!seq_is_available || seq_data_len == dl->seq_len, seq,
                "Bad line: according to CIGAR, expecting SEQ length to be %u but it is %u. SEQ=%.*s", 
                dl->seq_len, seq_data_len, seq_data_len, seq);

        ASSSEG (!qual_is_available || qual_data_len == dl->seq_len, qual,
                "Bad line: according to CIGAR, expecting QUAL length to be %u but it is %u. QUAL=%.*s", 
                dl->seq_len, dl->qual_data_len, dl->qual_data_len, qual);    
    }

    memcpy (&cigar_snip[cigar_snip_len], vb->last_cigar, last_cigar_len);
    cigar_snip_len += last_cigar_len;

    seg_by_did_i (vb, cigar_snip, cigar_snip_len, SAM_CIGAR, last_cigar_len+1); // +1 for \t
}

void sam_seg_qual_field (VBlockSAM *vb, ZipDataLineSAM *dl, const char *qual, uint32_t qual_data_len, unsigned add_bytes)
{
    dl->qual_data_start = qual - vb->txt_data.data;
    dl->qual_data_len   = qual_data_len;

    Context *qual_ctx = CTX(SAM_QUAL);
    qual_ctx->local.len += dl->qual_data_len;
    qual_ctx->txt_len   += add_bytes;
}

// test function called from main_load_reference -> txtfile_test_data: returns true if this line as pos=0 (i.e. unaligned)
bool sam_zip_is_unaligned_line (const char *line, int len)
{
    VBlock *vb = evb;
    const char *field_start, *next_field=line;
    char separator;
    unsigned field_len;

    GET_NEXT_ITEM (SAM_QNAME);
    GET_NEXT_ITEM (SAM_FLAG);
    GET_NEXT_ITEM (SAM_RNAME);
    GET_NEXT_ITEM (SAM_POS);

    return (field_len == 1 && *field_start == '0');
}

const char *sam_seg_txt_line (VBlock *vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    VBlockSAM *vb = (VBlockSAM *)vb_;
    ZipDataLineSAM *dl = DATA_LINE (vb->line_i);

    const char *next_field=field_start_line, *field_start;
    unsigned field_len=0;
    char separator;

    int32_t len = AFTERENT (char, vb->txt_data) - field_start_line;

    // QNAME - We break down the QNAME into subfields separated by / and/or : - these are vendor-defined strings. Examples:
    // Illumina: <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> for example "A00488:61:HMLGNDSXX:4:1101:15374:1031" see here: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
    // PacBio BAM: {movieName}/{holeNumber}/{qStart}_{qEnd} see here: https://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
    // BGI: E100020409L1C001R0030000234 (E100020409=Flow cell serial number, L1=Lane 1, C001R003=column 1 row 3, 0000234=Tile) Also see: https://github.com/IMB-Computational-Genomics-Lab/BGIvsIllumina_scRNASeq
    GET_NEXT_ITEM (SAM_QNAME);
    seg_compound_field (vb_, CTX(SAM_QNAME), field_start, field_len, sep_without_space, 0, 1 /* \n */);
    CTX(SAM_QNAME)->last_txt_index = ENTNUM (vb->txt_data, field_start); // store for kraken
    vb->last_txt_len (SAM_QNAME) = field_len;

    SEG_NEXT_ITEM (SAM_FLAG);
    int64_t flag;
    ASSSEG (str_get_int (field_start, field_len, &flag), field_start, "invalid FLAG field: %.*s", field_len, field_start);

    GET_NEXT_ITEM (SAM_RNAME);
    seg_chrom_field (vb_, field_start, field_len);

    // note: pos can have a value even if RNAME="*" - this happens if a SAM with a RNAME that is not in the header is converted to BAM with samtools
    GET_NEXT_ITEM (SAM_POS);
    PosType this_pos = seg_pos_field (vb_, SAM_POS, SAM_POS, 0, 0, field_start, field_len, 0, field_len+1);
    sam_seg_verify_rname_pos (vb_, SAM_RNAME_str, this_pos);
    
    random_access_update_pos (vb_, DC_PRIMARY, SAM_POS);

    SEG_NEXT_ITEM (SAM_MAPQ);

    // CIGAR - we wait to get more info from SEQ and QUAL
    GET_NEXT_ITEM (SAM_CIGAR);
    sam_analyze_cigar (vb, field_start, field_len, &dl->seq_len, &vb->ref_consumed, &vb->ref_and_seq_consumed, NULL);
    vb->last_cigar = field_start;
    unsigned last_cigar_len = field_len;
    ((char *)vb->last_cigar)[field_len] = 0; // nul-terminate CIGAR string

    SEG_NEXT_ITEM (SAM_RNEXT);
    
    GET_NEXT_ITEM (SAM_PNEXT);
    seg_pos_field (vb_, SAM_PNEXT, SAM_POS, 0, 0, field_start, field_len, 0, field_len+1);

    GET_NEXT_ITEM (SAM_TLEN);
    sam_seg_tlen_field (vb, field_start, field_len, 0, CTX(SAM_PNEXT)->last_delta, dl->seq_len);

    GET_NEXT_ITEM (SAM_SEQ);

    ASSSEG (dl->seq_len == field_len || vb->last_cigar[0] == '*' || field_start[0] == '*', field_start, 
            "seq_len implied by CIGAR=%s is %u, but actual SEQ length is %u, SEQ=%.*s", 
            vb->last_cigar, dl->seq_len, field_len, field_len, field_start);

    // calculate diff vs. reference (denovo or loaded)
    uint32_t save_ref_and_seq_consumed = vb->ref_and_seq_consumed; // save in case we have E2
    sam_seg_seq_field (vb, SAM_SQBITMAP, field_start, field_len, this_pos, vb->last_cigar, 0, field_len, vb->last_cigar, field_len+1);
    vb->ref_and_seq_consumed = save_ref_and_seq_consumed; // restore

    GET_MAYBE_LAST_ITEM (SAM_QUAL);
    sam_seg_qual_field (vb, dl, field_start, field_len, field_len + 1); 

    // finally we can seg CIGAR now
    sam_seg_cigar_field (vb, dl, last_cigar_len, SAM_SEQ_str, SAM_SEQ_len, field_start, field_len);
    
    // add BIN so this file can be reconstructed as BAM
    bam_seg_bin (vb, 0, (uint16_t)flag, this_pos);

    // OPTIONAL fields - up to MAX_FIELDS of them
    next_field = sam_seg_optional_all (vb, dl, next_field, len, has_13, separator, 0);

    SEG_EOL (SAM_EOL, false); /* last field accounted for \n */

    return next_field;
}
