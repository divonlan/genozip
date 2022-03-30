// ------------------------------------------------------------------
//   sam_seg.c
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
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
#include "dict_id.h"
#include "codec.h"
#include "container.h"
#include "stats.h"
#include "kraken.h"
#include "segconf.h"
#include "contigs.h"
#include "chrom.h"
#include "qname.h"
#include "lookback.h"
#include "gencomp.h"
#include "libdeflate/libdeflate.h"

typedef enum { QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, AUX } SamFields __attribute__((unused)); // quick way to define constants

char taxid_redirection_snip[100], xa_strand_pos_snip[100], XS_snip[30], XM_snip[30], MC_buddy_snip[30], 
     MQ_buddy_snip[30], MAPQ_buddy_snip[30], QUAL_buddy_snip[30], XA_lookback_snip[30],
     AS_buddy_snip[30], YS_buddy_snip[30], POS_buddy_snip[100], PNEXT_buddy_snip[100];
unsigned taxid_redirection_snip_len, xa_strand_pos_snip_len, XS_snip_len, XM_snip_len, MC_buddy_snip_len,
     MQ_buddy_snip_len, MAPQ_buddy_snip_len, QUAL_buddy_snip_len, XA_lookback_snip_len,
     AS_buddy_snip_len, YS_buddy_snip_len, POS_buddy_snip_len, PNEXT_buddy_snip_len;
WordIndex xa_lookback_strand_word_index = WORD_INDEX_NONE, xa_lookback_rname_word_index = WORD_INDEX_NONE;

// called by zfile_compress_genozip_header to set FlagsGenozipHeader.dt_specific
bool sam_zip_dts_flag (void)
{
    return flag.reference == REF_INTERNAL;
}

// ----------------------
// Seg stuff
// ----------------------

void sam_zip_free_end_of_z (void)
{
    sam_header_finalize(); 
}

// main thread, called for each component
void sam_zip_initialize (void)
{
    bool has_hdr_contigs = sam_hdr_contigs && sam_hdr_contigs->contigs.len;
    
    // we can use the aligner for unaligned reads IF it is loaded and we have no header contigs 
    segconf.sam_use_aligner = flag.aligner_available && !has_hdr_contigs; 

    // Copy header contigs to RNAME and RNEXT upon first component. This is in the order of the
    // header, as required by BAM (it encodes ref_id based on header order). Note, subsequent
    // bound files are required to have the same contigs or at least a contiguous subset starting a contig 0.
    // BAM always has header contigs (might be 0 of them, for an unaligned file), while SAM is allowed to be header-less
    if (z_file->num_txts_so_far == 1 && has_hdr_contigs) { 
        ctx_populate_zf_ctx_from_contigs (gref, SAM_RNAME, sam_hdr_contigs); 
        ctx_populate_zf_ctx_from_contigs (gref, SAM_RNEXT, sam_hdr_contigs);
    }

    // with REF_EXTERNAL and unaligned data, we don't know which chroms are seen (bc unlike REF_EXT_STORE, we don't use is_set), so
    // we just copy all reference contigs. this are not needed for decompression, just for --coverage/--sex/--idxstats
    if (z_file->num_txts_so_far == 1 && segconf.sam_use_aligner && flag.reference == REF_EXTERNAL) 
        ctx_populate_zf_ctx_from_contigs (gref, SAM_RNAME, ref_get_ctgs (gref)); 

    seg_prepare_snip_other (SNIP_REDIRECTION, _SAM_TAXID, false, 0, taxid_redirection_snip);

    static SmallContainer xa_strand_pos_con = { .repeats=1, .nitems_lo=2, .items = { { { _OPTION_XA_STRAND } }, { { _OPTION_XA_POS } } } };
    xa_strand_pos_snip_len = sizeof (xa_strand_pos_snip);
    container_prepare_snip ((ConstContainerP)&xa_strand_pos_con, 0, 0, xa_strand_pos_snip, &xa_strand_pos_snip_len); 

    seg_prepare_snip_other (SNIP_COPY, _OPTION_AS_i, false, 0, XS_snip);
    seg_prepare_snip_other (SNIP_COPY, _OPTION_NM_i, false, 0, XM_snip);

    qname_zip_initialize (SAM_QNAME);
    
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _SAM_PNEXT,   false, 0, POS_buddy_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _SAM_POS,     false, 0, PNEXT_buddy_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _SAM_MAPQ,    false, 0, MQ_buddy_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _SAM_QUAL,    false, 0, QUAL_buddy_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _OPTION_MQ_i, false, 0, MAPQ_buddy_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _OPTION_YS_i, false, 0, AS_buddy_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _OPTION_AS_i, false, 0, YS_buddy_snip);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY_MC, MC_buddy_snip, _SAM_CIGAR);

    seg_prepare_snip_other (SNIP_LOOKBACK, (DictId)_OPTION_XA_LOOKBACK, false, 0, XA_lookback_snip);

    sam_sa_prim_initialize_ingest(); // the PRIM component is compressed (out-of-band) at the same time as MAIN
    gencomp_initialize (SAM_COMP_PRIM, GCT_OOB); 
    gencomp_initialize (SAM_COMP_DEPN, GCT_DEPN); 
}

// called after each file
void sam_zip_finalize (void)
{
    gencomp_destroy();
}

void sam_zip_init_vb (VBlockP vb)
{
    vb->chrom_node_index = NODE_INDEX_NONE;    
    VB_SAM->check_for_gc = (!segconf.running && segconf.sam_is_sorted && vb->comp_i == SAM_COMP_MAIN);

    sam_gc_zip_init_vb (vb);
}

// called compute thread after compress, order of VBs is arbitrary
void sam_zip_after_compress (VBlockP vb)
{
    // Only the MAIN component produces gencomp lines, however we are processing VBs in order, so out-of-band VBs
    // need to be sent too, just to advance the serializing mutex
    if (vb->comp_i == SAM_COMP_MAIN || vb->comp_i == SAM_COMP_PRIM) 
        gencomp_absorb_vb_gencomp_lines (vb);
}

// called main thread, in order of VBs
void sam_zip_after_compute (VBlockP vb)
{
    if (vb->comp_i == SAM_COMP_MAIN) 
        sam_zip_gc_after_compute_main (VB_SAM);

    else if (vb->comp_i == SAM_COMP_PRIM) 
        gencomp_sam_prim_vb_has_been_ingested (vb);
}

// main thread: writing data-type specific fields to genozip header
void sam_zip_genozip_header (SectionHeaderGenozipHeader *header)
{
    header->sam.segconf_seq_len = BGEN32 (segconf.sam_seq_len); // introduced in v14
}

static void sam_seg_initialize_0X (VBlockP vb, DidIType lookback_did_i, DidIType rname_did_i, DidIType strand_did_i, DidIType pos_did_i, DidIType cigar_did_i)
{
    ContextP rname_ctx    = CTX(rname_did_i);
    ContextP strand_ctx   = CTX(strand_did_i);
    ContextP pos_ctx      = CTX(pos_did_i);
    ContextP cigar_ctx    = CTX(cigar_did_i);

    // note: we need to allocate lookback even if reps_per_line=0, lest an XA shows up despite not being in segconf
    if (lookback_did_i != DID_I_NONE) {
        ContextP lookback_ctx = CTX(lookback_did_i); // invalid if lookback_did_i=DID_I_NONE, that's ok
    
        rname_ctx->no_stons  = true;  // as we store by index
        strand_ctx->no_stons = true;
        lookback_ctx->flags.store        = STORE_INT;
        lookback_ctx->dynamic_size_local = true;
        lookback_ctx->local_param        = true;
        lookback_ctx->local.param        = lookback_size_to_local_param (1024);
        lookback_ctx->local_always       = (lookback_ctx->local.param != 0); // no need for a SEC_LOCAL section if the parameter is 0 (which is the default anyway)

        lookback_init (vb, lookback_ctx, rname_ctx,  STORE_INDEX); // lookback_ctx->local.param must be set before
        lookback_init (vb, lookback_ctx, strand_ctx, STORE_INDEX);
        lookback_init (vb, lookback_ctx, pos_ctx,    STORE_INT);
    }

    // create strand nodes (nodes will be deleted in sam_seg_finalize if not used)
    ctx_create_node (vb, strand_did_i, cSTR("-"));
    ctx_create_node (vb, strand_did_i, cSTR("+"));
    CTX(strand_did_i)->no_vb1_sort = true; // keep them in this ^ order

    cigar_ctx->no_stons = true; // as we use local to store long CIGARs in sam_seg_0A_cigar_cb
}

void sam_seg_initialize (VBlockP vb)
{
    START_TIMER;

    // all numeric fields need STORE_INT / STORE_FLOAT to be reconstructable to BAM (possibly already set)
    // via the translators set in the SAM_TOP2BAM Container
    DidIType numeric_fields[] = { SAM_TLEN, SAM_MAPQ, SAM_FLAG, SAM_POS, SAM_PNEXT, SAM_GPOS,
    /* non-fallback fields */     OPTION_NM_i, OPTION_AS_i, OPTION_MQ_i, OPTION_XS_i, OPTION_XM_i, OPTION_mc_i, OPTION_ms_i, OPTION_Z5_i, OPTION_tx_i,
                                  OPTION_YS_i };
    for (int i=0; i < ARRAY_LEN(numeric_fields); i++)
        CTX(numeric_fields[i])->flags.store = STORE_INT;

    CTX(SAM_SQBITMAP)->ltype    = LT_BITMAP;
    CTX(SAM_RNAME)->flags.store = STORE_INDEX; // since v12
    CTX(SAM_STRAND)->ltype      = LT_BITMAP;
    CTX(SAM_GPOS)->ltype        = LT_UINT32;
    CTX(SAM_BUDDY)->dynamic_size_local = true;
    CTX(SAM_QNAME)->no_stons    = true;        // no singletons, bc sam_piz_filter uses PEEK_SNIP
    CTX(SAM_QUAL) ->flags.store = STORE_INT;   // since v13 - store QUAL_score for buddy ms:i

    // MAPQ is uint8_t by BAM specification
    CTX(SAM_MAPQ)->ltype = CTX(OPTION_SA_MAPQ)->ltype = CTX(OPTION_OA_MAPQ)->ltype = LT_UINT8;
    
    if (segconf.sam_is_collated) 
        CTX(SAM_POS)->flags.store_delta = true; // since v12.0.41

    // a bug that existed 12.0.41-13.0.1 (bug 367): we stored buddy in machine endianty instead of BGEN32.
    // we use local.param=1 to indicate to reconstruct_set_buddy that this bug is now fixed.
    CTX(SAM_BUDDY)->local.param = 1;
    CTX(SAM_BUDDY)->local_param = true;

    // DEPN stuff
    if (sam_is_depn_vb) {
        CTX(SAM_SAGROUP)->ltype = sizeof(SAGroup)==4 ? LT_UINT32 : LT_UINT64;
        CTX(SAM_SAALN)->ltype   = LT_UINT16; // index of alignment with SA Group
        CTX(SAM_SEQSA)->ltype   = LT_BITMAP; // bitwise-xor of prim vs depn sequence
        CTX(SAM_QUALSA)->ltype  = LT_INT8;   // diff of prim vs depn qual

        buf_add_to_buffer_list_(vb, &CTX(SAM_SEQSA)->local, "local"); // because first alloc might be an anonymous realloc in bit_array_concat 
    }

    // PRIM stuff
    else if (sam_is_prim_vb) { 
        CTX(SAM_SAALN)->ltype              = LT_UINT16; // index of alignment with SA Group
        CTX(OPTION_SA_Z)->ltype            = LT_UINT8;    // we store num_alns in local
        CTX(SAM_QNAMESA)->st_did_i         = SAM_QNAME;   // consolidate stats
        
        // set store, consumed by sam_piz_prim_add_Alns
        CTX(OPTION_SA_RNAME)->flags.store  = STORE_INDEX;
        CTX(OPTION_SA_STRAND)->flags.store = STORE_INDEX; // 1 for - and 0 for +
        CTX(OPTION_SA_POS)->flags.store    = STORE_INT;
        CTX(OPTION_SA_MAPQ)->flags.store   = STORE_INT;
        CTX(OPTION_SA_NM)->flags.store     = STORE_INT;
    }

    // we sometimes copy from buddy alignments (mates etc), but only in sorted files
    if (segconf.sam_is_sorted) { 
        CTX(SAM_QNAME  )->flags.store_per_line = true; // 12.0.41
        CTX(SAM_FLAG   )->flags.store_per_line = true;
        CTX(SAM_POS    )->flags.store_per_line = true;
        CTX(SAM_PNEXT  )->flags.store_per_line = true;
        CTX(SAM_CIGAR  )->flags.store_per_line = true;
        CTX(OPTION_MC_Z)->flags.store_per_line = true;
        CTX(SAM_MAPQ   )->flags.store_per_line = true;
        CTX(OPTION_MQ_i)->flags.store_per_line = true; // 13.0.0
    }

    if (segconf.sam_bowtie2) {
        CTX(OPTION_AS_i)->flags.store_per_line = true; // 13.0.7
        CTX(OPTION_YS_i)->flags.store_per_line = true; // 13.0.7
    }

    if (segconf.sam_buddy_RG)
        CTX(OPTION_RG_Z)->flags.store_per_line = true;

    CTX(OPTION_MC_Z)->no_stons = true; // we're offloading to local ourselves
    CTX(SAM_CIGAR)->no_stons   = true;
    CTX(OPTION_MD_Z)->no_stons = true;
    if (segconf.has_bsseeker2) 
        CTX(OPTION_XM_Z)->no_stons = CTX(OPTION_XG_Z)->no_stons = true;

    CTX(OPTION_BI_Z)->no_stons = CTX(OPTION_BD_Z)->no_stons = true; // we can't use local for singletons in BD or BI as next_local is used by sam_piz_special_BD_BI to point into BD_BI
    CTX(OPTION_BD_BI)->ltype   = LT_SEQUENCE;
    
    // set dynamic_size_local to allow use of seg_integer
    CTX(OPTION_NM_i)->dynamic_size_local = true;
    CTX(OPTION_XM_i)->dynamic_size_local = true;
    CTX(OPTION_ms_i)->dynamic_size_local = true;
    
    Context *rname_ctx = CTX(SAM_RNAME);
    Context *rnext_ctx = CTX(SAM_RNEXT);

    rname_ctx->flags.store = rnext_ctx->flags.store = STORE_INDEX; // when reconstructing BAM, we output the word_index instead of the string
    rname_ctx->no_stons = rnext_ctx->no_stons = true;  // BAM reconstruction needs RNAME, RNEXT word indices. also needed for random access.

    // in --stats, consolidate stats 
    stats_set_consolidation (vb, SAM_SQBITMAP, 5, SAM_NONREF, SAM_NONREF_X, SAM_GPOS, SAM_STRAND, SAM_SEQSA);
    stats_set_consolidation (vb, SAM_QUAL,     2, SAM_DOMQRUNS, SAM_QUALSA);
    stats_set_consolidation (vb, OPTION_E2_Z,  4, OPTION_2NONREF, OPTION_N2ONREFX, OPTION_2GPOS, OPTION_S2TRAND);
    stats_set_consolidation (vb, OPTION_SA_Z,  6, OPTION_SA_RNAME, OPTION_SA_POS, OPTION_SA_STRAND, OPTION_SA_CIGAR, OPTION_SA_MAPQ, OPTION_SA_NM);
    stats_set_consolidation (vb, OPTION_OA_Z,  6, OPTION_OA_RNAME, OPTION_OA_POS, OPTION_OA_STRAND, OPTION_OA_CIGAR, OPTION_OA_MAPQ, OPTION_OA_NM);
    stats_set_consolidation (vb, OPTION_XA_Z,  7, OPTION_XA_RNAME, OPTION_XA_POS, OPTION_XA_STRAND, OPTION_XA_CIGAR, OPTION_XA_NM, OPTION_XA_STRAND_POS, OPTION_XA_LOOKBACK);
    stats_set_consolidation (vb, SAM_AUX,      1, SAM_MC_Z);
    stats_set_consolidation (vb, SAM_QNAME,    1, SAM_BUDDY);
    stats_set_consolidation (vb, OPTION_BD_BI, 2, OPTION_BI_Z, OPTION_BD_Z);

    codec_acgt_comp_init (vb);

    if (kraken_is_loaded) {
        CTX(SAM_TAXID)->flags.store    = STORE_INT;
        CTX(SAM_TAXID)->no_stons       = true; // must be no_stons the SEC_COUNTS data needs to mirror the dictionary words
        CTX(SAM_TAXID)->counts_section = true; 
    }

    qname_seg_initialize (vb, SAM_QNAME);
    sam_seg_QUAL_initialize (VB_SAM);
    sam_seg_SEQ_initialize (VB_SAM);

    if (segconf.running) {
        segconf.sam_is_sorted = segconf.sam_is_collated = segconf.MAPQ_has_single_value = true; // initialize optimistically
        segconf.sam_is_unmapped = true;  // we will reset this if finding a line with POS>0
        segconf.qname_flavor  = 0; // unknown
    }

    // initial allocations based on segconf data
    if (!segconf.running && segconf.sam_is_sorted) 
        buf_alloc_255 (vb, &VB_SAM->qname_hash, 0, (1 << BUDDY_HASH_BITS), int32_t, 1, "qname_hash");    

    sam_seg_initialize_0X (vb, (segconf.sam_is_sorted ? OPTION_XA_LOOKBACK : DID_I_NONE), OPTION_XA_RNAME, OPTION_XA_STRAND, OPTION_XA_POS, OPTION_XA_CIGAR);
    sam_seg_initialize_0X (vb, DID_I_NONE, OPTION_SA_RNAME, OPTION_SA_STRAND, OPTION_SA_POS, OPTION_SA_CIGAR);
    sam_seg_initialize_0X (vb, DID_I_NONE, OPTION_OA_RNAME, OPTION_OA_STRAND, OPTION_OA_POS, OPTION_OA_CIGAR);

    // create an "all the same" node for SAM_MC_Z
    ctx_create_node (vb, SAM_MC_Z, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_CONSUME_MC_Z }, 2);
    
    COPY_TIMER (seg_initialize);
}

void sam_seg_finalize (VBlockP vb)
{
    vb->flags.sam.is_collated = segconf.sam_is_collated;
    vb->flags.sam.is_sorted   = segconf.sam_is_sorted;

    // We always include the SQBITMAP local section, except if no lines
    if (vb->lines.len)
        CTX(SAM_SQBITMAP)->local_always = true;

    // assign the QUAL codec
    codec_assign_best_qual_codec (vb, SAM_QUAL, sam_zip_qual, VB_SAM->qual_codec_no_longr);

    if (!codec_domq_comp_init (vb, OPTION_U2_Z, sam_zip_U2)) 
        CTX(OPTION_U2_Z)->ltype  = LT_SEQUENCE; 

    // determine if sam_piz_sam2bam_SEQ ought to store vb->textual_seq
    CTX(SAM_SQBITMAP)->flags.no_textual_seq = CTX(SAM_QUAL)->lcodec != CODEC_LONGR
                                           && !segconf.has_bsseeker2;

    uint64_t qname_dict_id = (sam_is_prim_vb ? _SAM_QNAMESA : _SAM_QNAME);
    uint64_t seq_dict_id   = (sam_is_prim_vb ? _SAM_SEQSA   : _SAM_SQBITMAP);
    uint64_t qual_dict_id  = (sam_is_prim_vb ? _SAM_QUALSA  : _SAM_QUAL);
    
    // top level snip - reconstruction as SAM
    SmallContainer top_level_sam = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .callback     = true,
        .filter_items = true,
        .nitems_lo    = 13,
        .items        = { { .dict_id = { qname_dict_id }, .separator = "\t" },
                          { .dict_id = { _SAM_FLAG     }, .separator = "\t" },
                          { .dict_id = { _SAM_RNAME    }, .separator = "\t" },
                          { .dict_id = { _SAM_POS      }, .separator = "\t" },
                          { .dict_id = { _SAM_MAPQ     }, .separator = "\t" },
                          { .dict_id = { _SAM_CIGAR    }, .separator = "\t" },
                          { .dict_id = { _SAM_RNEXT    }, .separator = "\t" },
                          { .dict_id = { _SAM_PNEXT    }, .separator = "\t" },
                          { .dict_id = { _SAM_TLEN     }, .separator = "\t" },
                          { .dict_id = { seq_dict_id   }, .separator = "\t" },
                          { .dict_id = { qual_dict_id  }, .separator = "\t" },
                          { .dict_id = { _SAM_AUX      },                   },
                          { .dict_id = { _SAM_EOL      },                   } 
                        }
    };
    container_seg (vb, CTX(SAM_TOPLEVEL), (ContainerP)&top_level_sam, 0, 0, 0);

    // top level snip - reconstruction as BAM
    // strategy: we start by reconstructing the variable-length fields first (after a prefix that sets them in place) 
    // - read_name, cigar, seq and qual - and then go back and fill in the fixed-location fields
    // Translation (a feature of Container): items reconstruct their data and then call a translation function to translate it to the desired format
    SmallContainer top_level_bam = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .callback     = true,
        .filter_items = true,
        .nitems_lo    = 14,
        .items        = { { .dict_id = { _SAM_RNAME    }, .separator = { CI0_TRANS_NOR                     }, SAM2BAM_RNAME    }, // Translate - output word_index instead of string
                          { .dict_id = { _SAM_POS      }, .separator = { CI0_TRANS_NOR | CI0_TRANS_MOVE, 1 }, SAM2BAM_POS      }, // Translate - output little endian POS-1
                          { .dict_id = { _SAM_MAPQ     }, .separator = { CI0_TRANS_NOR                     }, SAM2BAM_U8       }, // Translate - textual to binary number
                          { .dict_id = { _SAM_BAM_BIN  }, .separator = { CI0_TRANS_NOR | CI0_TRANS_MOVE, 2 }, SAM2BAM_LTEN_U16 }, // Translate - textual to binary number
                          { .dict_id = { _SAM_FLAG     }, .separator = { CI0_TRANS_NOR | CI0_TRANS_MOVE, 4 }, SAM2BAM_LTEN_U16 }, // Translate - textual to binary number
                          { .dict_id = { _SAM_RNEXT    }, .separator = { CI0_TRANS_NOR                     }, SAM2BAM_RNAME    }, // Translate - output word_index instead of string
                          { .dict_id = { _SAM_PNEXT    }, .separator = { CI0_TRANS_NOR | CI0_TRANS_MOVE, 4 }, SAM2BAM_POS      }, // Translate - output little endian POS-1
                          { .dict_id = { qname_dict_id }, .separator = { CI0_TRANS_NUL                     }                   }, // normal 
                          { .dict_id = { _SAM_CIGAR    }, .separator = ""                                                      }, // handle in special reconstructor - translate textual to BAM CIGAR format + reconstruct l_read_name, n_cigar_op, l_seq
                          { .dict_id = { _SAM_TLEN     }, .separator = { CI0_TRANS_NOR                     }, SAM2BAM_TLEN     }, // must be after CIGAR bc sam_piz_special_TLEN_old needs vb->seq_num
                          { .dict_id = { seq_dict_id   }, .separator = "",                                    SAM2BAM_SEQ      }, // Translate - textual format to BAM format
                          { .dict_id = { qual_dict_id  }, .separator = "",                                    SAM2BAM_QUAL     }, // Translate - textual format to BAM format, set block_size
                          { .dict_id = { _SAM_AUX      }, .separator = { CI0_TRANS_NOR                     }                   }, // up to v11, this had the SAM2BAM_AUX translator
                        }
    };

    // no container wide-prefix, skip l_name with a 4-character prefix
    static const char bam_line_prefix[] = { CON_PX_SEP, // has prefix 
                                            CON_PX_SEP, // end of (empty) container-wide prefix
                                            ' ',' ',' ',' ', CON_PX_SEP }; // first item prefix - 4 spaces (place holder for block_size)

    container_seg (vb, CTX(SAM_TOP2BAM), (ContainerP)&top_level_bam, bam_line_prefix, sizeof(bam_line_prefix), 
                   IS_BAM_ZIP ? sizeof (uint32_t) * vb->lines.len : 0); // if BAM, account for block_size

    // top level snip - reconstruction as FASTQ
    SmallContainer top_level_fastq = { 
        .repeats     = vb->lines.len,
        .is_toplevel = true,
        .callback    = true,  // drop supplementary alignments and alignments without QUAL data
        .nitems_lo   = 9,
        .items       = { { .dict_id = { qname_dict_id }, .separator = "\n"                  }, 
                         { .dict_id = { _SAM_RNAME    }, .separator = { CI0_TRANS_NOR }     }, // needed for reconstructing seq 
                         { .dict_id = { _SAM_POS      }, .separator = { CI0_TRANS_NOR }     }, // needed for reconstructing seq
                         { .dict_id = { _SAM_PNEXT    }, .separator = { CI0_TRANS_NOR }     }, // needed for reconstructing POS (in case of BUDDY)
                         { .dict_id = { _SAM_FLAG     }, .separator = { CI0_TRANS_NOR }, .translator = SAM2FASTQ_FLAG }, // need to know if seq is reverse complemented & if it is R2 ; reconstructs "1" for R1 and "2" for R2
                         { .dict_id = { _SAM_CIGAR    }, .separator = { CI0_TRANS_NOR }     }, // needed for reconstructing seq
                         { .dict_id = { _SAM_MC_Z     }, .separator = { CI0_TRANS_NOR }     }, // consumes OPTION_MC_Z if its on this line, might be needed for buddy CIGAR
                         { .dict_id = { seq_dict_id   }, .separator = "\n",              .translator = SAM2FASTQ_SEQ  }, 
                         { .dict_id = { qual_dict_id  }, .separator = "\n",              .translator = SAM2FASTQ_QUAL }, // also moves fastq "line" to R2 (paired file) if needed
                       }
    };

    // add a '@' to the description line, use a prefix to add the + line    
    static const char fastq_line_prefix[] = { CON_PX_SEP, CON_PX_SEP, '@', CON_PX_SEP, CON_PX_SEP, CON_PX_SEP, 
                                              CON_PX_SEP, CON_PX_SEP, CON_PX_SEP, CON_PX_SEP, CON_PX_SEP, '+', '\n', CON_PX_SEP };

    container_seg (vb, CTX(SAM_TOP2FQ), (ContainerP)&top_level_fastq, fastq_line_prefix, sizeof(fastq_line_prefix), 0);

    // top level snip - reconstruction as FASTQ "extended" - with all the SAM fields in the description line
    SmallContainer top_level_fastq_ext = { 
        .repeats     = vb->lines.len,
        .is_toplevel = true,
        .callback    = true,  // drop non-primary chimeric reads and reads without QUAL data
        .nitems_lo   = 12,
        .items       = { { .dict_id = { qname_dict_id }, .separator = "\t"                               }, 
                         { .dict_id = { _SAM_FLAG     }, .separator = "\t", .translator = SAM2FASTQ_FLAG }, // need to know if seq is reverse complemented & if it is R2 ; reconstructs "1" for R1 and "2" for R2
                         { .dict_id = { _SAM_RNAME    }, .separator = "\t"                               },
                         { .dict_id = { _SAM_POS      }, .separator = "\t"                               },
                         { .dict_id = { _SAM_MAPQ     }, .separator = "\t"                               },
                         { .dict_id = { _SAM_CIGAR    }, .separator = "\t"                               },
                         { .dict_id = { _SAM_RNEXT    }, .separator = "\t"                               },
                         { .dict_id = { _SAM_PNEXT    }, .separator = "\t"                               },
                         { .dict_id = { _SAM_TLEN     }, .separator = "\t"                               },
                         { .dict_id = { _SAM_AUX      }, .separator = "\n"                               },
                         { .dict_id = { seq_dict_id   }, .separator = "\n", .translator =SAM2FASTQ_SEQ   }, 
                         { .dict_id = { qual_dict_id  }, .separator = "\n", .translator =SAM2FASTQ_QUAL  }, // also moves fastq "line" to R2 (paired file) if needed
                       }
    };

    // add a '@' to the description line, use a prefix to add the + line    
    static const char fastq_ext_line_prefix[] = { 
        CON_PX_SEP, CON_PX_SEP, 
        '@', CON_PX_SEP, 
        'F','L','A','G',':', CON_PX_SEP, 
        'R','N','A','M','E',':', CON_PX_SEP, 
        'P','O','S',':', CON_PX_SEP, 
        'M','A','P','Q',':', CON_PX_SEP, 
        'C','I','G','A','R',':', CON_PX_SEP,
        'R','N','E','X','T',':', CON_PX_SEP, 
        'P','N','E','X','T',':', CON_PX_SEP, 
        'T','L','E','N',':', CON_PX_SEP, 
        CON_PX_SEP,              // optional
        CON_PX_SEP,              // SEQ
        '+', '\n', CON_PX_SEP }; // QUAL

    container_seg (vb, CTX(SAM_TOP2FQEX), (ContainerP)&top_level_fastq_ext, 
                   fastq_ext_line_prefix, sizeof(fastq_ext_line_prefix), 0);

    // finalize Seg configuration parameters    
    if (segconf.running) {
        segconf.sam_buddy_RG  = segconf.sam_is_sorted && (CTX(OPTION_RG_Z)->nodes.len >= BUDDY_MIN_RG_COUNT); // it only makes sense to buddy RGs in a sorted (not collated) file, and where there are "a lot" of RGs or else we would be increasing rather than decreasing entropy
        segconf.sam_cigar_len = 1 + ((segconf.sam_cigar_len-1) / vb->lines.len); // set to the average CIGAR len (rounded up)
        segconf.sam_seq_len   = (uint32_t)(0.5 + (double)segconf.sam_seq_len / (double)vb->lines.len); // set average seq_len - rounded to the nearest 

        if (segconf.is_long_reads) {

            if (segconf.has[OPTION_ms_i])
                segconf.sam_ms_type = ms_MINIMAP2; // definitely not biobambam's MateBaseScore as long reads don't have mates
        
            segconf.sam_is_collated = false; // long reads are never paired-end
        }
        else {
            if (segconf.has[OPTION_ms_i] && segconf.sam_is_paired)
                // TODO: test not strong enough - minimap2 may be used for short reads too
                segconf.sam_ms_type = ms_BIOBAMBAM;
        }
    }

    // get rid of the 0A strand contexts (nodes added in sam_seg_initialize_0X) if we ended up not have using them
    // (note: we use SA_STRAND in PRIM lines even if there is no SA:Z field)
    if (!CTX(OPTION_XA_STRAND)->b250.len) ctx_free_context (CTX(OPTION_XA_STRAND), OPTION_XA_STRAND);
    if (!CTX(OPTION_OA_STRAND)->b250.len) ctx_free_context (CTX(OPTION_OA_STRAND), OPTION_OA_STRAND);
    if (!CTX(OPTION_SA_STRAND)->b250.len) ctx_free_context (CTX(OPTION_SA_STRAND), OPTION_SA_STRAND);

    // case: we have XA but didn't use lookback (because its a non-sorted file) - create an all-the-same XA_LOOKBACK dict
    if (CTX(OPTION_XA_Z)->b250.len && !CTX(OPTION_XA_Z)->local.len)
        ctx_create_node (vb, OPTION_XA_LOOKBACK, "0", 1);

    // Primary VB - ingest VB data into z_file->sa_*
    if (sam_is_prim_vb)
        sam_zip_prim_ingest_vb (VB_SAM);
}

// main thread: called after all VBs, before compressing global sections
void sam_zip_after_vbs (void)
{
    // shorten unused RNAME / RNEXT dictionary strings to "" (dict pre-populated in sam_zip_initialize)
    // note: we prevent shortening of RNAMEs that appear in SA:Z fields of primary alignments, but not
    // in the main RNAME field, by incrementing their count in sam_seg_prim_add_sa_group
    ctx_shorten_unused_dict_words (SAM_RNAME);
    ctx_shorten_unused_dict_words (SAM_RNEXT);
}

bool sam_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return 
        // typically small 
        dict_id.num == _SAM_TOPLEVEL     ||
        dict_id.num == _SAM_TOP2BAM      ||
        dict_id.num == _SAM_TOP2FQ       ||
        dict_id.num == _SAM_TOP2FQEX     ||
        dict_id.num == _SAM_FLAG         ||
        dict_id.num == _SAM_MAPQ         ||
        dict_id.num == _SAM_QNAME        ||
        dict_id.num == _SAM_RNAME        ||
        dict_id.num == _SAM_RNEXT        ||
        dict_id.num == _SAM_AUX          ||
        dict_id.num == _SAM_EOL          ||
        dict_id.num == _SAM_TAXID        ||
        dict_id.num == _SAM_BUDDY        ||

        // standard tags, see here: https://samtools.github.io/hts-specs/SAMtags.pdf
        dict_id.num == _OPTION_AM_i      ||
        dict_id.num == _OPTION_AS_i      ||
        dict_id.num == _OPTION_CC_Z      ||
        dict_id.num == _OPTION_CM_i      ||
        dict_id.num == _OPTION_LB_Z      ||
        dict_id.num == _OPTION_FI_i      ||
        dict_id.num == _OPTION_H0_i      ||
        dict_id.num == _OPTION_H1_i      ||
        dict_id.num == _OPTION_H2_i      ||
        dict_id.num == _OPTION_HI_i      ||
        dict_id.num == _OPTION_MQ_i      ||
        dict_id.num == _OPTION_NH_i      ||
        dict_id.num == _OPTION_NM_i      ||
        dict_id.num == _OPTION_OA_Z      ||
        dict_id.num == _OPTION_OA_RNAME  ||
        dict_id.num == _OPTION_OA_NM     ||
        dict_id.num == _OPTION_OA_STRAND ||
        dict_id.num == _OPTION_OA_MAPQ   || 
        dict_id.num == _OPTION_OA_CIGAR  || 
        dict_id.num == _OPTION_OC_Z      ||
        dict_id.num == _OPTION_PG_Z      ||
        dict_id.num == _OPTION_PQ_i      ||
        dict_id.num == _OPTION_PU_Z      ||
        dict_id.num == _OPTION_RG_Z      ||
        dict_id.num == _OPTION_SA_Z      ||
        dict_id.num == _OPTION_SA_RNAME  ||
        dict_id.num == _OPTION_SA_NM     ||
        dict_id.num == _OPTION_SA_STRAND ||
        dict_id.num == _OPTION_SA_MAPQ   || 
        dict_id.num == _OPTION_SA_CIGAR  || 
        dict_id.num == _OPTION_SM_i      ||
        dict_id.num == _OPTION_TC_i      ||
        dict_id.num == _OPTION_UQ_i      ||
        
        // bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
        dict_id.num == _OPTION_X0_i      ||
        dict_id.num == _OPTION_X1_i      ||
        dict_id.num == _OPTION_XA_Z      ||
        dict_id.num == _OPTION_XA_RNAME  ||
        dict_id.num == _OPTION_XA_NM     ||
        dict_id.num == _OPTION_XA_STRAND ||
        dict_id.num == _OPTION_XA_CIGAR  || 
        dict_id.num == _OPTION_XN_i      ||
        dict_id.num == _OPTION_XM_i      ||
        dict_id.num == _OPTION_XO_i      ||
        dict_id.num == _OPTION_XG_i      ||
        dict_id.num == _OPTION_XS_i      ||
        dict_id.num == _OPTION_XE_i;
}

static void sam_get_one_optional (VBlockSAMP vb, STRp(aux),
                                  rom *tag, char *type, char *array_subtype, pSTRp(value)) // out
{
    ASSSEG (aux_len >= 6 && aux[2] == ':' && aux[4] == ':', aux, "invalid optional field format: %.*s", aux_len, aux);

    *tag         = aux;
    *type        = aux[3];
    
    if (*type == 'B') {
        *array_subtype = aux[5];
        *value         = aux + 7;
        *value_len     = aux_len - 7;
    }
    else {
        *array_subtype = 0;
        *value         = aux + 5;
        *value_len     = aux_len - 5;
    }

    ASSSEG0 (*value_len < 10000000, aux, "Invalid array field format"); // check that *value_len didn't underflow beneath 0
}

// returns aux field if it exists or NULL if it doesn't
// In BAM, with a numeric field, result is returned in numeric, otherwise in the return value + value_len
static rom sam_seg_get_aux (VBlockSAMP vb, rom name, STRps (aux), uint32_t *value_len, ValueType *numeric, bool is_bam)
{
    if ((int32_t)n_auxs < 1) return NULL; // this line has no AUX fields (possibly negatuve)

    for (uint32_t f=0; f < n_auxs; f++) {
        rom aux = auxs[f];
        if (aux_lens[f] > (is_bam ? 3 : 5) && 
            name[0] == aux[0] && 
            name[1] == aux[1] && 
            name[3] == (is_bam ? sam_seg_bam_type_to_sam_type(aux[2]) : aux[3])) {

            STR0 (my_value);
            rom tag;
            char sam_type=0, bam_type=0, array_subtype=0;

            if (is_bam) bam_get_one_optional (vb, STRi(aux,f), &tag, &bam_type, &array_subtype, pSTRa(my_value), numeric);
            else        sam_get_one_optional (vb, STRi(aux,f), &tag, &sam_type, &array_subtype, pSTRa(my_value));

            if (value_len) *value_len = my_value_len;
            return (is_bam && !my_value) ? &aux[3] : my_value; // in BAM, if NULL, value is in numeric
        }
    }

    if (value_len) *value_len = 0;
    return NULL;
}

rom sam_seg_get_aux_str (VBlockSAMP vb, rom name, STRps (aux), uint32_t *value_len, bool is_bam)
{
    ValueType numeric __attribute__((unused)); 
    return sam_seg_get_aux (vb, name, STRas(aux), value_len, &numeric, is_bam);
}

// returns the length of the field, or 0 if it is not found
uint32_t sam_seg_get_aux_int (VBlockSAMP vb, rom name, STRps (aux), 
                              int64_t *number, // modified only if integer is parsed
                              bool is_bam)
{
    ValueType numeric; 
    uint32_t value_len;
    rom value = sam_seg_get_aux (vb, name, STRas(aux), &value_len, &numeric, is_bam);
    if (value && is_bam) {
        *number = numeric.i;
        return value_len;
    }
    else if (value && str_get_int (STRa(value), number))
        return value_len;
    else 
        return 0; // this line doesn't have this field, or (SAM only) the field is not a valid integer
}

void sam_seg_verify_RNAME_POS (VBlockSAMP vb, rom p_into_txt, PosType this_pos)
{
    if (segconf.running) return;

    if (flag.reference == REF_INTERNAL && (!sam_hdr_contigs /* SQ-less SAM */ || !sam_hdr_contigs->contigs.len /* SQ-less BAM */)) return;
    if (!this_pos || (vb->chrom_name_len==1 && *vb->chrom_name=='*')) return; // unaligned
    
    PosType LN;
    if (sam_hdr_contigs) {
        // since this SAM file has a header, all RNAMEs must be listed in it (the header contigs appear first in CTX(RNAME), see sam_zip_initialize
        ASSSEG (vb->chrom_node_index < sam_hdr_contigs->contigs.len, p_into_txt, "RNAME \"%.*s\" does not have an SQ record in the header", vb->chrom_name_len, vb->chrom_name);
        LN = contigs_get_LN (sam_hdr_contigs, vb->chrom_node_index);
    }
    else { // headerless SAM
        WordIndex ref_index = chrom_2ref_seg_get (gref, VB, vb->chrom_node_index); // possibly an alt contig
        if (ref_index == WORD_INDEX_NONE) {
            WARN_ONCE ("FYI: RNAME \"%.*s\" (and possibly others) is missing in the reference file. This might impact the compression ratio.", 
                       vb->chrom_name_len, vb->chrom_name);
            return; // the sequence will be segged as unaligned
        }
        LN = contigs_get_LN (ref_get_ctgs (gref), ref_index);
    }

    ASSINP (this_pos <= LN, "%s: Error POS=%"PRId64" is beyond the size of \"%.*s\" which is %"PRId64". In vb=%u line_i=%"PRIu64" chrom_node_index=%d", 
            txt_name, this_pos, vb->chrom_name_len, vb->chrom_name, LN, vb->vblock_i, vb->line_i, vb->chrom_node_index);
}

static void sam_seg_get_kraken (VBlockSAMP vb, bool *has_kraken, 
                                char *taxid_str, rom *tag, char *sam_type, 
                                pSTRp (value), ValueType *numeric, // out - value for SAM, numeric for BAM
                                bool is_bam)
{
    *tag        = "tx"; // genozip introduced tag (=taxid)
    *sam_type   = 'i';
    *value      = taxid_str;
    *has_kraken = false;
    *value_len  = kraken_seg_taxid_do (VB, SAM_TAXID, last_txt (vb, SAM_QNAME), vb->last_txt_len (SAM_QNAME),
                                       taxid_str, true);
    *numeric = CTX(SAM_TAXID)->last_value;

    vb->recon_size += is_bam ? 7 : (*value_len + 6); // txt modified
}

void sam_seg_aux_all (VBlockSAMP vb, ZipDataLineSAM *dl, STRps(aux))
{
    const bool is_bam = IS_BAM_ZIP;
    Container con = { .repeats=1 };
    char prefixes[MAX_FIELDS * 6 + 3]; // each name is 5 characters per SAM specification, eg "MC:Z:" followed by CON_PX_SEP ; +3 for the initial CON_PX_SEP
    prefixes[0] = prefixes[1] = prefixes[2] = CON_PX_SEP; // initial CON_PX_SEP follow by separator of empty Container-wide prefix followed by separator for empty prefix for translator-only item[0]
    unsigned prefixes_len=3;

    // item[0] is translator-only item - to translate the Container itself in case of reconstructing BAM 
    con.items[con_nitems(con)] = (ContainerItem){ .translator = SAM2BAM_AUX_SELF }; 
    con_inc_nitems (con);
    
    bool has_kraken = kraken_is_loaded;

    for (uint32_t f=0 ; f < n_auxs + has_kraken; f++) {

        STR0(value);
        ValueType numeric = {};
        rom tag;
        char sam_type=0, bam_type=0, array_subtype=0;
        char taxid_str[20];

        if (f == n_auxs) sam_seg_get_kraken   (vb, &has_kraken, taxid_str,  &tag, &bam_type, pSTRa(value), &numeric, is_bam);
        else if (is_bam) bam_get_one_optional (vb, STRi(aux,f), &tag, &bam_type, &array_subtype, pSTRa(value), &numeric);
        else             sam_get_one_optional (vb, STRi(aux,f), &tag, &sam_type, &array_subtype, pSTRa(value));

        if (sam_type == 'i') {
            ASSERT (str_get_int (STRa(value), &numeric.i), "Expecting integer value for auxiliary field %c%c but found \"%.*s\". vb=%u line=%"PRIu64,
                    tag[0], tag[1], value_len, value, vb->vblock_i, vb->line_i);
            value=0;
        }
            
        if (!bam_type) bam_type = sam_seg_sam_type_to_bam_type (sam_type, numeric.i);
        if (!sam_type) sam_type = sam_seg_bam_type_to_sam_type (bam_type);
        
        ASSERT (bam_type, "value %"PRId64" of field %c%c is out of range of the BAM specification: [%d-%u]. vb=%u line=%"PRIu64, 
                numeric.i, tag[0], tag[1], -0x80000000, 0x7fffffff, vb->vblock_i, vb->line_i);

        con.items[con_nitems(con)] = (ContainerItem) {
            .dict_id    = sam_seg_aux_field (vb, dl, is_bam, tag, bam_type, array_subtype, STRa(value), numeric),
            .translator = aux_field_translator ((uint8_t)bam_type), // how to transform the field if reconstructing to BAM
            .separator  = { aux_sep_by_type[is_bam][(uint8_t)bam_type], '\t' },
        };
        con_inc_nitems (con);

        ASSSEG (con_nitems(con) <= MAX_FIELDS, value, "too many optional fields, limit is %u", MAX_FIELDS);

        // note: in the optional field prefix (unlike array type), all integer types become 'i'.
        char prefix[6] = { tag[0], tag[1], ':', sam_type, ':', CON_PX_SEP}; 
        memcpy (&prefixes[prefixes_len], prefix, 6);
        prefixes_len += 6;
    }

    uint32_t num_items = con_nitems(con);
    if (num_items > 1) { // we have Aux fields, not just the translator item
        if (con.items[num_items-1].separator[0] & 0x80) // is a flag
            con.items[num_items-1].separator[0] &= ~(CI0_NATIVE_NEXT & ~(uint8_t)0x80); // last Aux field has no tab
        con.items[num_items-1].separator[1] = 0;
        container_seg (vb, CTX(SAM_AUX), &con, prefixes, prefixes_len, (is_bam ? 3 : 5) * (num_items-1)); // account for : SAM: "MX:i:" BAM: "MXi"
    }
    else
        // NULL means MISSING Container item (of the toplevel container) - will cause container_reconstruct_do of 
        // the toplevel container to delete of previous separator (\t)
        container_seg (vb, CTX(SAM_AUX), 0, 0, 0, 0); 
}

void sam_seg_QNAME (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qname), unsigned add_additional_bytes)
{
    // in segconf, identify if this file is collated (each QNAME appears in two or more consecutive lines)
    if (segconf.running) {
        if (segconf.sam_is_collated) { // still collated, try to find evidence to the contrary
            #define last_is_new ctx_specific // QNAME is different than previous line's QNAME  
            bool is_new = (qname_len != vb->last_txt_len (SAM_QNAME) || memcmp (qname, last_txt(vb, SAM_QNAME), qname_len));
            if (is_new && CTX(SAM_QNAME)->last_is_new)
                segconf.sam_is_collated = false; // two new QNAMEs in a row = not collated
            CTX(SAM_QNAME)->last_is_new = is_new;
        }
        
        if (vb->line_i==0)
            qname_segconf_discover_flavor (VB, SAM_QNAME, STRa(qname));
    }

    ContextP ctx = CTX(SAM_QNAME);
    bool segged_against_buddy = false;

    // I tried: SNIP_COPY on collation (with segconf.sam_is_collated) but QNAME.b250 worsening (from all_the_same) outweighs the b250 compression 
    // improvement of the container items
    if (segconf.sam_is_sorted && !segconf.running) {
        uint32_t hash = adler32 (1971, qname, qname_len) & ((1 << BUDDY_HASH_BITS)-1);
        int32_t *hash_ent = B(int32_t, vb->qname_hash, hash);
        int32_t buddy_line_i = *hash_ent;
        
        ZipDataLineSAM *buddy_dl = DATA_LINE (buddy_line_i); // possibly an invalid address if prev_line_i=-1, that's ok

        // case: we found a previous line with an identical QNAME - that line is now our buddy
        if (buddy_line_i != -1 && buddy_dl->QNAME.len == qname_len && 
            !memcmp (qname, Bc (vb->txt_data, buddy_dl->QNAME.index), qname_len)) {

            seg_add_to_local_resizable (VB, CTX(SAM_BUDDY), vb->line_i - buddy_line_i, 0); // add buddy (delta) >= 1 . 

            // case: Seg against buddy - the second SNIP_COPY_BUDDY means that reconstruction of this snip will also set buddy
            // note: if DEPN with SA Group we could Seg against buddy or against SA Group. Segging against buddy
            // note: in PRIM segging against buddy here is used for loading SA Groups, for reconstruction SAM_QNAMESA is used
            // has the advantage that all the other buddied fields (RNEXT, TLEN...) also benefit from the buddy.
            seg_by_ctx (VB, (char[]){ SNIP_COPY_BUDDY, SNIP_COPY_BUDDY }, 2, ctx, qname_len + add_additional_bytes); // seg QNAME as copy-from-buddy (an extra SNIP_COPY_BUDDY indicates that reconstruct_from_buddy should set buddy_line_i here)
            
            vb->buddy_line_i = buddy_line_i; 
            segged_against_buddy = true;
        }

        // update hash, possibly override existing value - that's ok, we prefer a line higher up (smaller delta)
        // exception: if current read is Supplamentary and buddy is not, we keep it as is so it buddies with its mate and its additional supps too
        if (buddy_line_i == -1 || buddy_dl->FLAG.bits.supplementary || buddy_dl->FLAG.bits.secondary || 
            !(dl->FLAG.bits.supplementary || dl->FLAG.bits.secondary))
            *hash_ent = vb->line_i; 
    }

    // case: DEPN with SA Group: seg against SA Group (unless already segged against buddy)
    if (sam_is_depn_vb && vb->sa_grp && !segged_against_buddy)
        sam_seg_against_sa_group (vb, ctx, qname_len + add_additional_bytes);

    // case: PRIM: seg against SA Group - store in SAM_QNAMESA - Reconstruct will take from here in PRIM per Toplevel container
    else if (sam_is_prim_vb) {
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_PRIM_QNAME }, 2, SAM_QNAMESA, 0); // consumed when reconstructing PRIM vb
        if (!segged_against_buddy) qname_seg (VB, ctx, STRa(qname), add_additional_bytes); // consumed when loading SA Groups
    }

    // MAIN, DEPN without SA Group - seg normally if not already segged against buddy
    // Note: In PRIM, this is consumed by sam_piz_load_SA_Groups, not reconstruct
    else if ((sam_is_main_vb || (sam_is_depn_vb && !vb->sa_grp)) && !segged_against_buddy)
        qname_seg (VB, ctx, STRa(qname), add_additional_bytes);
}

// We seg against a previous buddy line's MQ if one exists, but not if this is a single-MAPQ-value file
void sam_seg_MAPQ (VBlockSAMP vb, ZipDataLineSAM *dl, unsigned add_bytes)
{
    ContextP ctx = CTX(SAM_MAPQ);

    if (segconf.running && dl->MAPQ) {
        if (!segconf.MAPQ_value) 
            segconf.MAPQ_value = dl->MAPQ;
        else if (segconf.MAPQ_value != dl->MAPQ) 
            segconf.MAPQ_has_single_value = false;
    }

    ctx_set_last_value (VB, ctx, (int64_t)dl->MAPQ);

    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // an invalid pointer if buddy_line_i is -1

    if (sam_seg_has_SA_Group(vb)) {
        sam_seg_against_sa_group (vb, ctx, add_bytes);

        // in PRIM, we also seg it as the first SA alignment (used for PIZ to load alignments to memory, not used for reconstructing SA)
        if (sam_is_prim_vb) {
            ContextP sa_mapq_ctx = CTX(OPTION_SA_MAPQ);
            seg_add_to_local_nonresizeable (vb, sa_mapq_ctx, dl->MAPQ, false, 0);

            // count MAPQ field contribution to OPTION_SA_MAPQ, so sam_stats_reallocate can allocate the z_data between MAPQ and SA:Z
            sa_mapq_ctx->counts.count += add_bytes; 
        }
    }

    else if (!segconf.running && vb->buddy_line_i != -1 && dl->MAPQ && !segconf.MAPQ_has_single_value && dl->MAPQ == buddy_dl->MQ)
        seg_by_ctx (VB, STRa(MAPQ_buddy_snip), ctx, add_bytes); // copy MQ from earlier-line buddy 

    else 
        seg_add_to_local_nonresizeable (vb, CTX(SAM_MAPQ), dl->MAPQ, true, add_bytes);
}

void sam_seg_RNAME_RNEXT (VBlockSAMP vb, DidIType did_i, STRp (chrom), unsigned add_bytes)
{
    if (did_i == SAM_RNEXT || !sam_seg_has_SA_Group(vb)) {
        bool is_new;
        chrom_seg_ex (VB, did_i, STRa(chrom), 0, NULL, add_bytes, !IS_BAM_ZIP, &is_new);

        // don't allow adding chroms to a BAM file or a SAM that has SQ lines in the header, but we do allow to add to a headerless SAM.
        ASSSEG (!is_new || !sam_hdr_contigs || segconf.running || (chrom_len==1 && (*chrom=='*' || *chrom=='=')),
                chrom, "contig '%.*s' appears in file, but is missing in the %s header", STRf(chrom), dt_name (vb->data_type));
    }

    // case: RNAME - PRIM or DEPN vb - seg against SA group with alignments
    // Note: in DEPN, rname already verified in sam_sa_seg_depn_find_sagroup to be as in SA alignment
    else {
        sam_seg_against_sa_group (vb, CTX(SAM_RNAME), add_bytes);

        if (sam_is_prim_vb) {
            // in PRIM, we also seg it as the first SA alignment (used for PIZ to load alignments to memory, not used for reconstructing SA)
            seg_by_did_i (VB, STRa(chrom), OPTION_SA_RNAME, 0);

            // protect rname from removal by ctx_shorten_unused_dict_words if it is unused in the main RNAME field. 
            ctx_protect_from_removal (VB, CTX(SAM_RNAME), vb->chrom_node_index); 

            // count RNAME field contribution to OPTION_SA_RNAME, so sam_stats_reallocate can allocate the z_data between RNAME and SA:Z
            CTX(OPTION_SA_RNAME)->counts.count += add_bytes; 
        }

        // note: vb->chrom_node_index already set in bam/sam_seg_txt_line
        STRset (vb->chrom_name, chrom);

        random_access_update_chrom (VB, 0, vb->chrom_node_index, chrom, chrom_len); 
    }
}

// test function called from main_load_reference -> txtfile_test_data: returns true if this line as pos=0 (i.e. unaligned)
bool sam_zip_is_unaligned_line (rom line, int len)
{
    VBlockP vb = evb;
    rom field_start, next_field=line;
    char separator;
    unsigned field_len;

    GET_NEXT_ITEM (SAM_QNAME);
    GET_NEXT_ITEM (SAM_FLAG);
    GET_NEXT_ITEM (SAM_RNAME);
    GET_NEXT_ITEM (SAM_POS);

    return (field_len == 1 && *field_start == '0');
}

rom sam_seg_txt_line (VBlockP vb_, rom next_line, uint32_t remaining_txt_len, bool *has_13)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    sam_reset_line (VB);

    WordIndex prev_line_chrom = vb->chrom_node_index;
    PosType prev_line_pos = vb->last_int (SAM_POS);

    ZipDataLineSAM *dl = DATA_LINE (vb->line_i);
        
    str_split_by_tab (next_line, remaining_txt_len, MAX_FIELDS+AUX, has_13); // also advances next_line to next line

    dl->SEQ.index = BNUMtxt (flds[SEQ]); // SEQ.len to be determined by sam_cigar_analyze
    dl->QNAME = TXTWORDi(fld, QNAME);
    dl->QUAL  = TXTWORDi(fld, QUAL);
    
    if (fld_lens[QUAL]==1 && flds[QUAL][0]=='*') 
        vb->qual_missing = QUAL_MISSING_STANDARD;

    // lazy way to get vb->chrom* (inc. if --match-chrom), rollback later if seg RNAME is not needed
    if (vb->check_for_gc || !sam_is_main_vb)
        seg_create_rollback_point (VB, NULL, 1, SAM_RNAME);

    chrom_seg_ex (VB, SAM_RNAME, STRfld(RNAME), 0, NULL, fld_lens[RNAME]+1, true, NULL);

    // xxx // for a header-full SAM, we expected to be in the SAM header, needed by sam_seg_is_gc_line and 
    // // sam_sa_seg_depn_find_sagroup (that only run in headerful files). Otherwise it is set in chrom_seg_ex.
    // if (!segconf.running && sam_hdr_contigs && !(fld_lens[RNAME]==1 && *flds[RNAME]=='*')) {
    //     vb->chrom_node_index = ctx_get_ol_node_index_by_snip (VB, CTX(SAM_RNAME), STRfld(RNAME));
    //     vb->chrom_name       = flds[RNAME];
    //     vb->chrom_name_len   = fld_lens[RNAME];
    //     ASSSEG (vb->chrom_node_index != WORD_INDEX_NONE, flds[RNAME], "RNAME=\"%.*s\", not found in SAM header", vb->chrom_name_len, vb->chrom_name);
    // }

    ASSSEG (str_get_int (STRfld(POS), &dl->POS), 
            flds[POS], "Invalid POS \"%.*s\": expecting an integer", fld_lens[POS], flds[POS]);

    ASSSEG (str_get_int_range8 (STRfld(MAPQ), 0, 255, &dl->MAPQ), 
            flds[MAPQ], "Invalid MAPQ \"%.*s\": expecting an integer [0,255]", fld_lens[MAPQ], flds[MAPQ]);

    ASSSEG (str_get_int_range16 (STRfld(FLAG), 0, SAM_MAX_FLAG, &dl->FLAG.value), flds[FLAG], "invalid FLAG field: \"%.*s\"", fld_lens[FLAG], flds[FLAG]);

    // analyze (but don't seg yet) CIGAR
    uint32_t seq_len;
    sam_cigar_analyze (vb, STRfld(CIGAR), false, &seq_len);
    dl->SEQ.len = seq_len; // do it this way to avoid compiler warning

    ASSSEG (dl->SEQ.len == fld_lens[SEQ] || (fld_lens[CIGAR] == 1 && *flds[CIGAR] == '*') || flds[SEQ][0] == '*', flds[SEQ], 
            "seq_len implied by CIGAR=%.*s is %u, but actual SEQ length is %u, SEQ=%.*s", 
            fld_lens[CIGAR], flds[CIGAR], dl->SEQ.len, fld_lens[SEQ], fld_lens[SEQ], flds[SEQ]);

    // if this is a secondary / supplamentary read (aka Dependent) or a read that has an associated sec/sup read (aka Primary) - move
    // the line to the appropriate component and skip it here (no segging done yet)
    if (vb->check_for_gc && sam_seg_is_gc_line (vb, dl, flds[0], next_line - flds[0], STRasi(fld,AUX), false)) {
        seg_rollback (VB); // cancelling segging of RNAME
        goto done;
    }

    vb->last_cigar = flds[CIGAR];
    SAFE_NUL (&vb->last_cigar[fld_lens[CIGAR]]); // nul-terminate CIGAR string

    if (!sam_is_main_vb) {
        sam_seg_sa_group_stuff (vb, dl , STRasi(fld,AUX), STRfld(CIGAR), flds[SEQ], false);

        // re-seg rname, against SA group
        seg_rollback (VB);
        sam_seg_RNAME_RNEXT (vb, SAM_RNAME, STRfld(RNAME), fld_lens[RNAME]+1);
    }

    sam_seg_QNAME (vb, dl, STRfld(QNAME), 1);

    sam_seg_FLAG (vb, dl, fld_lens[FLAG]+1);
    

    // note: pos can have a value even if RNAME="*" - this happens if a SAM with a RNAME that is not in the header is converted to BAM with samtools
    sam_seg_POS (vb, dl, prev_line_chrom, fld_lens[POS]+1);

    if (fld_lens[RNAME] != 1 || *flds[RNAME] != '*')
        sam_seg_verify_RNAME_POS (vb, flds[RNAME], dl->POS);

    sam_seg_MAPQ (vb, dl, fld_lens[MAPQ]+1);

    sam_seg_RNAME_RNEXT (vb, SAM_RNEXT, STRfld(RNEXT), fld_lens[RNEXT]+1);
    
    sam_seg_PNEXT (vb, dl, STRfld(PNEXT), 0, prev_line_pos, fld_lens[PNEXT]+1);

    // we search forward for MD:Z now, XG:Z as we will need it for SEQ if it exists
    STR(MD);
    if (segconf.has[OPTION_MD_Z] && !segconf.running &&
        (MD = sam_seg_get_aux_str (vb, "MD:Z", STRasi(fld,AUX), &MD_len, false))) 
        sam_seg_MD_Z_analyze (vb, STRa(MD), dl->POS, vb->last_cigar);

    STR(XG);
    if (segconf.has_bsseeker2 &&
        (XG = sam_seg_get_aux_str (vb, "XG:Z", STRasi(fld,AUX), &XG_len, false))) 
        sam_seg_XG_Z_analyze (vb, dl, STRa(XG), dl->POS);

    seg_set_last_txt (VB, CTX(SAM_SQBITMAP), STRfld(SEQ));

    // calculate diff vs. reference (denovo or loaded)
    sam_seg_SEQ (vb, SAM_SQBITMAP, STRfld(SEQ), dl->POS, vb->last_cigar, dl->FLAG.bits.rev_comp, vb->ref_consumed, vb->ref_and_seq_consumed, 
                 0, fld_lens[SEQ], vb->last_cigar, fld_lens[SEQ]+1);
    
    sam_seg_QUAL (vb, dl, STRfld(QUAL), fld_lens[QUAL] + 1); 
    
    ASSSEG (str_is_in_range (flds[QUAL], fld_lens[QUAL], 33, 126), flds[QUAL], "Invalid QUAL - it contains non-Phred characters: \"%.*s\"", 
            fld_lens[QUAL], flds[QUAL]);

    if (fld_lens[SEQ] != fld_lens[QUAL])
        vb->qual_codec_no_longr = true; // we cannot compress QUAL with CODEC_LONGR in this case

    // finally we can seg CIGAR now
    sam_cigar_seg_textual (vb, dl, fld_lens[CIGAR], STRfld(SEQ), STRfld(QUAL));
    
    // add BIN so this file can be reconstructed as BAM
    bam_seg_BIN (vb, dl, 0, false);

    // AUX fields - up to MAX_FIELDS of them
    sam_seg_aux_all (vb, dl, n_flds-AUX, &flds[AUX], &fld_lens[AUX]);

    // finally, we can seg TLEN now, after MC:Z, if it exists
    bool is_rname_rnext_same = (fld_lens[RNEXT]==1 && *flds[RNEXT]=='=') || 
                               (fld_lens[RNEXT]==fld_lens[RNAME] && !memcmp (flds[RNEXT], flds[RNAME], fld_lens[RNAME]));

    sam_seg_TLEN (vb, dl, STRfld(TLEN), 0, is_rname_rnext_same);

    SEG_EOL (SAM_EOL, false); /* last field accounted for \n */
    SAFE_RESTORE; // restore \t after CIGAR

done:
    return next_line;
}
