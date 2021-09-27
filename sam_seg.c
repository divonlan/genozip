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
#include "kraken.h"
#include "segconf.h"
#include "contigs.h"
#include "chrom.h"


// callback function for compress to get data of one line (called by codec_bz2_compress)
void sam_zip_qual (VBlock *vb, uint64_t vb_line_i, char **line_qual_data, uint32_t *line_qual_len, uint32_t maximum_len) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_qual_len  = MIN_(maximum_len, dl->qual_data_len);
    
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
    bool has_hdr_contigs = sam_hdr_contigs && sam_hdr_contigs->contigs.len;
    
    // we can use the aligner for unaligned reads IF it is loaded and we have no header contigs 
    segconf.sam_use_aligner = flag.aligner_available && !has_hdr_contigs; 

    // Copy header contigs to RNAME and RNEXT upon first component. This is in the order of the
    // header, as required by BAM (it encodes ref_id based on header order). Note, subsequent
    // bound files are required to have the same contigs or at least a contiguous subset starting a contig 0.
    // BAM always has header contigs (might be 0 of them, for an unaligned file), while SAM is allowed to be header-less
    if (z_file->num_txt_components_so_far == 1 && has_hdr_contigs) { 
        ctx_populate_zf_ctx_from_contigs (gref, SAM_RNAME, sam_hdr_contigs); 
        ctx_populate_zf_ctx_from_contigs (gref, SAM_RNEXT, sam_hdr_contigs);
    }

    // with REF_EXTERNAL and unaligned data, we don't know which chroms are seen (bc unlike REF_EXT_STORE, we don't use is_set), so
    // we just copy all reference contigs. this are not need for decompression, just for --coverage/--sex/--idxstats
    if (z_file->num_txt_components_so_far == 1 && segconf.sam_use_aligner && flag.reference == REF_EXTERNAL) 
        ctx_populate_zf_ctx_from_contigs (gref, SAM_RNAME, ref_get_ctgs (gref)); 

    sam_aux_zip_initialize();
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

    CTX(OPTION_BI)->no_stons    = CTX(OPTION_BD)->no_stons = true; // we can't use local for singletons in BD or BI as next_local is used by sam_piz_special_BD_BI to point into BD_BI
    CTX(OPTION_BI)->st_did_i    = CTX(OPTION_BD)->st_did_i = OPTION_BD_BI; 
    CTX(OPTION_BD_BI)->ltype    = LT_SEQUENCE;
    CTX(OPTION_NM)->flags.store = STORE_INT;

    Context *rname_ctx = CTX(SAM_RNAME);
    Context *rnext_ctx = CTX(SAM_RNEXT);

    rname_ctx->flags.store = rnext_ctx->flags.store = STORE_INDEX; // when reconstructing BAM, we output the word_index instead of the string
    rname_ctx->no_stons = rnext_ctx->no_stons = true;  // BAM reconstruction needs RNAME, RNEXT word indices. also needed for random access.

    // in --stats, consolidate stats 
    stats_set_consolidation (vb, SAM_SQBITMAP, 4, SAM_NONREF, SAM_NONREF_X, SAM_GPOS, SAM_STRAND);
    stats_set_consolidation (vb, OPTION_E2, 4, OPTION_2NONREF, OPTION_N2ONREFX, OPTION_2GPOS, OPTION_S2TRAND);
    stats_set_consolidation (vb, OPTION_SA, 6, OPTION_SA_RNAME, OPTION_SA_POS, OPTION_SA_STRAND, OPTION_SA_CIGAR, OPTION_SA_MAPQ, OPTION_SA_NM);
    stats_set_consolidation (vb, OPTION_OA, 6, OPTION_OA_RNAME, OPTION_OA_POS, OPTION_OA_STRAND, OPTION_OA_CIGAR, OPTION_OA_MAPQ, OPTION_OA_NM);
    stats_set_consolidation (vb, OPTION_XA, 6, OPTION_XA_RNAME, OPTION_XA_POS, OPTION_XA_STRAND, OPTION_XA_CIGAR, OPTION_XA_NM, OPTION_XA_STRAND_POS);

    CTX(OPTION_SA_CIGAR)->no_stons = CTX(OPTION_OA_CIGAR)->no_stons = CTX(OPTION_XA_CIGAR)->no_stons = true;
    CTX(OPTION_SA_CIGAR)->ltype    = CTX(OPTION_OA_CIGAR)->ltype    = CTX(OPTION_XA_CIGAR)->ltype    = LT_TEXT;

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

    if (!codec_domq_comp_init (vb, OPTION_U2, sam_zip_U2)) 
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

static bool sam_seg_get_MD (const char *txt /* points to SEQ field */, uint32_t txt_len, pSTRp(md))
{
    const char *after = &txt[txt_len];
    SAFE_ASSIGN(after, '\n'); // for safety

    unsigned column = 9; // SEQ
    *md = NULL;

    for (;*txt != '\n'; txt++) 
        if (*txt == '\t') {
            if (! *md) { 
                column++;
                // case: start of MD
                if (column >= 11 && after - txt > 6 && !memcmp (txt+1, "MD:Z:", 5)) 
                    *md = txt+6; // skip \tMD:Z:
            }
            // case: \t at end of MD
            else break;
        }

    // next currently points either to the \n or \t at after of the MD field
    if (*md) *md_len = txt - *md;

    SAFE_RESTORE;
    
    return *md != NULL;
}

void sam_seg_verify_rname_pos (VBlock *vb, const char *p_into_txt, PosType this_pos)
{
    if (segconf.running) return;

    if (flag.reference == REF_INTERNAL && (!sam_hdr_contigs /* SQ-less SAM */ || !sam_hdr_contigs->contigs.len /* SQ-less BAM */)) return;
    if (!this_pos) return; // unaligned
    
    PosType LN;
    if (sam_hdr_contigs) {
        // since this SAM file has a header, all RNAMEs must be listed in it (the header contigs appear first in CTX(RNAME), see sam_zip_initialize
        ASSSEG (vb->chrom_node_index < sam_hdr_contigs->contigs.len, p_into_txt, "RNAME \"%.*s\" does not have an SQ record in the header", vb->chrom_name_len, vb->chrom_name);
        LN = contigs_get_LN (sam_hdr_contigs, vb->chrom_node_index);
    }
    else { // headerless SAM
        WordIndex ref_index = chrom_2ref_seg_get (gref, vb, vb->chrom_node_index); // possibly an alt contig
        if (ref_index == WORD_INDEX_NONE) {
            WARN_ONCE ("FYI: RNAME \"%.*s\" (and possibly others) is missing in the reference file. No harm.", 
                       vb->chrom_name_len, vb->chrom_name);
            return; // the sequence will be segged as unaligned
        }
        LN = contigs_get_LN (ref_get_ctgs (gref), ref_index);
    }

    ASSINP (this_pos <= LN, "%s: Error POS=%"PRId64" is beyond the size of \"%.*s\" which is %"PRId64". In vb=%u line_i=%"PRIu64" chrom_node_index=%d", 
            txt_name, this_pos, vb->chrom_name_len, vb->chrom_name, LN, vb->vblock_i, vb->line_i, vb->chrom_node_index);
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
        seg_by_ctx (VB, snip_delta, 2, ctx, add_bytes);
    }
    // case 2:
    else if (tlen_value > 0 && pnext_pos_delta > 0 && cigar_seq_len > 0) {
        char tlen_by_calc[30] = { SNIP_SPECIAL, SAM_SPECIAL_TLEN };
        unsigned tlen_by_calc_len = str_int (tlen_value - pnext_pos_delta - (int64_t)cigar_seq_len, &tlen_by_calc[2]);
        seg_by_ctx (VB, tlen_by_calc, tlen_by_calc_len + 2, ctx, add_bytes);
    }
    // case default: add as is (option 1)
    else if (tlen)
        seg_by_ctx (VB, tlen, tlen_len, ctx, add_bytes);

    // case default: add as is (option 2)
    else {
        char snip[20];
        unsigned snip_len = str_int (tlen_value, snip);
        seg_by_ctx (VB, snip, snip_len, ctx, add_bytes);
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
void sam_seg_seq_field (VBlockSAM *vb, DidIType bitmap_did, STRp(seq), const PosType pos, const char *cigar, 
                        unsigned recursion_level, uint32_t level_0_seq_len, const char *level_0_cigar, unsigned add_bytes)
{
    START_TIMER;

    Context *bitmap_ctx = CTX(bitmap_did);
    Context *nonref_ctx = bitmap_ctx + 1;

    ASSERT (recursion_level < 500, "excess recursion recursion_level=%u seq_len=%u level_0_cigar=%s", // a large number of recursion calls can happen if CIGAR=9M10910N86M3274690N30M1S as observed with the STAR aligner https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
            recursion_level, seq_len, level_0_cigar);

    bitmap_ctx->txt_len += add_bytes; // byte counts for --show-sections

    // for unaligned lines, if we have refhash loaded, use the aligner instead of CIGAR-based segmenting
    if (!pos && segconf.sam_use_aligner) {
        aligner_seg_seq (VB, bitmap_ctx, seq, seq_len);
        goto align_nonref_local;
    }

    BitArray *bitmap = buf_get_bitarray (&bitmap_ctx->local);
    uint32_t bitmap_start = bitmap_ctx->next_local;        

    if (!recursion_level) {

        ASSERTW (seq_len < 1000000, "Warning: sam_seg_seq_field: seq_len=%u is suspiciously high and might indicate a bug", seq_len);

        buf_alloc (vb, &bitmap_ctx->local, roundup_bits2bytes64 (vb->ref_and_seq_consumed), vb->lines.len * (vb->ref_and_seq_consumed+5) / 8, uint8_t, CTX_GROWTH, "contexts->local"); 
        buf_extend_bits (&bitmap_ctx->local, vb->ref_and_seq_consumed); 
        bit_array_set_region (bitmap, bitmap_ctx->next_local, vb->ref_and_seq_consumed); // we initiaze all the bits to "set", and clear as needed.

        buf_alloc (vb, &nonref_ctx->local, seq_len + 3, vb->lines.len * seq_len / 4, uint8_t, CTX_GROWTH, "contexts->local"); 
    
        if (segconf.running) return; // case segconf: we created the contexts for segconf_set_vb_size accounting. that's enough - actually segging will will mark is_set and break sam_md_analyze.
    }

    // we can't compare to the reference if it is unaligned: we store the seqeuence in nonref without an indication in the bitmap
    if (!pos || (vb->chrom_name_len==1 && vb->chrom_name[0]=='*')) {
        buf_add (&nonref_ctx->local, seq, seq_len);
        goto align_nonref_local; 
    }

    if (seq[0] == '*') {
        vb->md_verified = false;    
        goto done; // we already handled a missing seq (SEQ="*") by adding a '-' to CIGAR - no data added here
    }

    bool no_cigar = cigar[0] == '*' && cigar[1] == 0; // there's no CIGAR. 

    RefLock lock;
    Range *range = no_cigar ? NULL : ref_seg_get_locked_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, vb->ref_consumed, WORD_INDEX_NONE, seq, &lock);

    // Cases where we don't consider the refernce and just copy the seq as-is
    if (!range) { // 1. (denovo:) this contig defined in @SQ went beyond the maximum genome size of 4B and is thus ignored
                  // 2. (loaded:) case contig doesn't exist in the reference
                  // 3. no cigar - the sequence is not aligned to the reference even if we have RNAME and POS (and its length can exceed the reference contig)
        buf_add (&nonref_ctx->local, seq, seq_len);
        
        bit_array_clear_region (bitmap, bitmap_ctx->next_local, vb->ref_and_seq_consumed); // note: vb->ref_and_seq_consumed==0 if cigar="*"
        bitmap_ctx->next_local += vb->ref_and_seq_consumed;

        random_access_update_last_pos (VB, DC_PRIMARY, pos + vb->ref_consumed - 1);

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

    uint32_t ref_len_this_level = (flag.reference == REF_INTERNAL ? MIN_(vb->ref_consumed, range->last_pos - pos + 1)
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
                    bit_i++; // bit remains set. cannot increment inside the macro
                }

                // case our seq is identical to the reference at this site
                else if (range && normal_base && seq[i] == ref_base_by_idx (range, actual_next_ref)) {
                    bit_i++; // bit remains set. 

                    if (flag.reference == REF_EXT_STORE)
                        bit_array_set (&range->is_set, actual_next_ref); // we will need this ref to reconstruct
                }
                
                // case: ref is set to a different value - we store our value in nonref_ctx
                else {
                    NEXTENT (char, nonref_ctx->local) = seq[i];
                    bit_array_clear (bitmap, bit_i); bit_i++;
                    vb->mismatch_bases++;
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
            unsigned ref_consumed = (flag.reference == REF_INTERNAL ? MIN_(subcigar_len, range_len - next_ref)
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

        char updated_cigar[strlen (next_cigar) + 20];
        if (subcigar_len) sprintf (updated_cigar, "%u%c%s", subcigar_len, cigar_op, next_cigar);

        sam_seg_seq_field (vb, bitmap_did, seq + i, seq_len - i, range->last_pos + 1, subcigar_len ? updated_cigar : next_cigar, recursion_level + 1, level_0_seq_len, level_0_cigar, 0);
    }
    else { // update RA of the VB with last pos of this line as implied by the CIGAR string
        if (this_seq_last_pos <= range->last_pos) // always the case in INTERNAL and non-circular EXTERNAL 
            random_access_update_last_pos (VB, DC_PRIMARY, this_seq_last_pos);

        else  // we circled back to the beginning for the chromosome - i.e. this VB RA is the entire chromosome
            random_access_update_to_entire_chrom (VB, DC_PRIMARY, range->first_pos, range->last_pos);
    }

    // final verification step - does MD:Z correctly reflect matches and mismatches of M/X/=
    if (!recursion_level) {
        BitArray *M_is_ref = buf_get_bitarray (&vb->md_M_is_ref);

        bool bitmap_matches_MD = vb->md_verified && 
                                !bit_array_manhattan_distance (M_is_ref, 0, bitmap, bitmap_start, M_is_ref->nbits);

        if (flag.show_wrong_md && vb->md_verified && !bitmap_matches_MD) {
            iprintf ("vb=%u line=%"PRIu64" RNAME=%.*s POS=%"PRId64" CIGAR=%s MD=%.*s SEQ=%.*s\n", 
                    vb->vblock_i, vb->line_i, STRf(vb->chrom_name), pos, vb->last_cigar, vb->last_txt_len(OPTION_MD), last_txt(vb, OPTION_MD), STRf(seq));
            bit_array_print_substr ("SEQ match to ref", bitmap, bitmap_start, M_is_ref->nbits, info_stream);
            bit_array_print_substr ("MD implied match", M_is_ref, 0, M_is_ref->nbits, info_stream); 
        }

        vb->md_verified = bitmap_matches_MD;
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

static const char *sam_seg_get_kraken (VBlockSAM *vb, const char *next_field, bool *has_kraken, 
                                       char *taxid_str, const char **tag, char *type,  // out
                                       const char **value, unsigned *value_len, // out
                                       bool is_bam)
{
    *tag        = "TX"; // genozip introduced tag (=taxid)
    *type       = 'i';
    *value      = taxid_str;
    *has_kraken = false;
    *value_len  = kraken_seg_taxid_do (VB, SAM_TAXID, last_txt (vb, SAM_QNAME), vb->last_txt_len (SAM_QNAME),
                                       taxid_str, true);

    vb->recon_size += is_bam ? 7 : (*value_len + 6); // txt modified

    return next_field; // unmodified
}


static const char *sam_get_one_optional (VBlockSAM *vb, const char *next_field, int32_t len, char *separator_p, bool *has_13, 
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

        dl->seq_len = MAX_(seq_data_len, dl->qual_data_len); // one or both might be not available and hence =1

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

    seg_by_did_i (VB, cigar_snip, cigar_snip_len, SAM_CIGAR, last_cigar_len+1); // +1 for \t
}

void sam_seg_qual_field (VBlockSAM *vb, ZipDataLineSAM *dl, const char *qual, uint32_t qual_data_len, unsigned add_bytes)
{
    dl->qual_data_start = qual - vb->txt_data.data;
    dl->qual_data_len   = qual_data_len;

    Context *qual_ctx = CTX(SAM_QUAL);
    qual_ctx->local.len += dl->qual_data_len;
    qual_ctx->txt_len   += add_bytes;
}

void sam_seg_qname_field (VBlockSAM *vb, const char *qname, uint32_t qname_len, unsigned add_additional_bytes)
{
    // I tried: SNIP_COPY on collation (with segconf.is_collated) but QNAME.b250 worsening (from all_the_same) outweighs the b250 compression 
    // improvement of the container items
    ContextP ctx = CTX(SAM_QNAME);
    seg_compound_field (VB, ctx, qname, qname_len, sep_without_space, 0, add_additional_bytes);

    seg_set_last_txt (VB, ctx, STRa(qname), STORE_NONE); // store for kraken
}

void sam_seg_rname_rnext (VBlockP vb, DidIType did_i, STRp (chrom), unsigned add_bytes)
{
    bool is_new;
    chrom_seg_ex (VB, did_i, STRa(chrom), 0, NULL, add_bytes, !IS_BAM, &is_new);

    // don't allow adding chroms to a BAM file or a SAM that has SQ lines in the header, but we do allow to add to a headerless SAM.
    ASSSEG (!is_new || !sam_hdr_contigs || segconf.running || (chrom_len==1 && (*chrom=='*' || *chrom=='=')),
            chrom, "contig '%.*s' appears in file, but is missing in the %s header", STRf(chrom), dt_name (vb->data_type));
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

    int32_t len = REMAINING (vb->txt_data, field_start_line);

    // QNAME - We break down the QNAME into subfields separated by / and/or : - these are vendor-defined strings. Examples:
    // Illumina: <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> for example "A00488:61:HMLGNDSXX:4:1101:15374:1031" see here: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
    // PacBio BAM: {movieName}/{holeNumber}/{qStart}_{qEnd} see here: https://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
    // BGI: E100020409L1C001R0030000234 (E100020409=Flow cell serial number, L1=Lane 1, C001R003=column 1 row 3, 0000234=Tile) Also see: https://github.com/IMB-Computational-Genomics-Lab/BGIvsIllumina_scRNASeq
    GET_NEXT_ITEM (SAM_QNAME);
    sam_seg_qname_field (vb, STRdid(SAM_QNAME), 1);

    SEG_NEXT_ITEM (SAM_FLAG);
    int64_t flag;
    ASSSEG (str_get_int (STRdid(SAM_FLAG), &flag), field_start, "invalid FLAG field: %.*s", field_len, field_start);

    GET_NEXT_ITEM (SAM_RNAME);
    sam_seg_rname_rnext (VB, SAM_RNAME, STRdid(SAM_RNAME), SAM_RNAME_len+1);

    // note: pos can have a value even if RNAME="*" - this happens if a SAM with a RNAME that is not in the header is converted to BAM with samtools
    GET_NEXT_ITEM (SAM_POS);
    PosType this_pos = seg_pos_field (VB, SAM_POS, SAM_POS, 0, 0, STRdid(SAM_POS), 0, SAM_POS_len+1);
    
    if (SAM_RNAME_len != 1 || *SAM_RNAME_str != '*')
        sam_seg_verify_rname_pos (VB, SAM_RNAME_str, this_pos);
    
    random_access_update_pos (VB, DC_PRIMARY, SAM_POS);

    SEG_NEXT_ITEM (SAM_MAPQ);

    // CIGAR - we wait to get more info from SEQ and QUAL
    GET_NEXT_ITEM (SAM_CIGAR);
    sam_analyze_cigar (vb, STRdid(SAM_CIGAR), &dl->seq_len);
    vb->last_cigar = SAM_CIGAR_str;
    unsigned last_cigar_len = SAM_CIGAR_len;
    ((char *)vb->last_cigar)[SAM_CIGAR_len] = 0; // nul-terminate CIGAR string

    GET_NEXT_ITEM (SAM_RNEXT);
    sam_seg_rname_rnext (VB, SAM_RNEXT, STRdid(SAM_RNEXT), SAM_RNEXT_len+1);
    
    GET_NEXT_ITEM (SAM_PNEXT);
    seg_pos_field (VB, SAM_PNEXT, SAM_POS, 0, 0, STRdid(SAM_PNEXT), 0, SAM_PNEXT_len+1);

    GET_NEXT_ITEM (SAM_TLEN);
    sam_seg_tlen_field (vb, STRdid(SAM_TLEN), 0, CTX(SAM_PNEXT)->last_delta, dl->seq_len);

    // we search forward for MD:Z now, as we will need it for SEQ if it exists
    if (segconf.has_MD && !segconf.running) {
        STR(md); 
        if (sam_seg_get_MD (next_field, remaining_txt_len, pSTRa(md)))
            sam_md_analyze (vb, STRa(md), this_pos, vb->last_cigar);
    }

    GET_NEXT_ITEM (SAM_SEQ);

    ASSSEG (dl->seq_len == field_len || vb->last_cigar[0] == '*' || SAM_SEQ_str[0] == '*', SAM_SEQ_str, 
            "seq_len implied by CIGAR=%s is %u, but actual SEQ length is %u, SEQ=%.*s", 
            vb->last_cigar, dl->seq_len, SAM_SEQ_len, SAM_SEQ_len, SAM_SEQ_str);

    // calculate diff vs. reference (denovo or loaded)
    uint32_t save_ref_and_seq_consumed = vb->ref_and_seq_consumed; // save in case we have E2
    sam_seg_seq_field (vb, SAM_SQBITMAP, STRdid(SAM_SEQ), this_pos, vb->last_cigar, 0, field_len, vb->last_cigar, SAM_SEQ_len+1);
    vb->ref_and_seq_consumed = save_ref_and_seq_consumed; // restore

    GET_MAYBE_LAST_ITEM (SAM_QUAL);
    sam_seg_qual_field (vb, dl, STRdid(SAM_QUAL), SAM_QUAL_len + 1); 

    // finally we can seg CIGAR now
    sam_seg_cigar_field (vb, dl, last_cigar_len, STRdid(SAM_SEQ), STRdid(SAM_QUAL));
    
    // add BIN so this file can be reconstructed as BAM
    bam_seg_bin (vb, 0, (uint16_t)flag, this_pos);

    // OPTIONAL fields - up to MAX_FIELDS of them
    next_field = sam_seg_optional_all (vb, dl, next_field, len, has_13, separator, 0);

    SEG_EOL (SAM_EOL, false); /* last field accounted for \n */

    return next_field;
}
