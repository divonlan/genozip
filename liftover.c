// ------------------------------------------------------------------
//   liftover.c
//   Copyright (C) 2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// Dual-coordinates genozip/genocat file flow:
//
// genozip -C source.txt       --> dual-coord.genozip - component_i=0 contains all lines, each with LIFTOVER or LIFTREJT 
//                                                                    LIFTREFD are also written to the rejects file
//                                                      component_i=1 is the rejects file containing LIFTREFD lines (again)
//
// genocat dual-crd.genozip    --> primary.txt -        component_i=0 is reconstructed - lines have LIFTOVER or LIFTREJT 
//                                                                    lines are IN ORDER
//                                                      component_i=1 is skipped.
//
// genozip primary.txt         --> dual-coord.genozip - component_i=0 contains all lines, each with LIFTOVER or LIFTREJT 
//                                                                    LIFTREFD are also written to the rejects file
//                                                      component_i=1 is the rejects file containing LIFTREFD lines (again)
//
// genocat -v dual-crd.genozip --> luft.txt -           component_i=1 with LIFTREJT is reconstructed first and becomes part of the header
//                                                      component_i=0 is reconstructed - dropping LIFTREJT lines and lifting over LIFTOVER->LIFTBACK 
//                                                                    lines are OUT OF ORDER due to change in coordinates 
//
// genozip luft.txt            --> dual-coord.genozip - LIFTREJT header lines are sent to vblock_i=1 via unconsumed_txt
//                                                      component_i=0 contains all lines first all LIFTREJT lines followed by all LIFTOVER lines  
//                                                                    LIFTREFD lines are also written to the rejects file
//                                                      component_i=1 is the rejects file containing LIFTREFD lines (again)
//
// genocat dual-crd.genozip    --> primary.txt          component_i=0 is reconstructed - first LIFTREJT then LIFTOVER lines
//                                                                    lines are OUT OF ORDER (LIFTREJT first)
//                                                      component_i=1 is skipped.
//

#include <errno.h>
#include <math.h>
#include "liftover.h"
#include "base250.h"
#include "context.h"
#include "flags.h"
#include "chain.h"
#include "data_types.h"
#include "file.h"
#include "buffer.h"
#include "vblock.h"
#include "hash.h"
#include "container.h"
#include "seg.h"
#include "dict_id.h"
#include "codec.h"
#include "strings.h"
#include "random_access.h"

enum { OCHROM_OFFSET, OPOS_OFFSET, OREF_OFFSET, OSTRAND_OFFSET, OALTRULE_OFFSET, OSTATUS_OFFSET };

const char *liftover_status_names[] = LIFTOVER_STATUS_NAMES;

#define WORD_INDEX_MISSING -2 // a value in liftover_chrom2chainsrc

// constant snips initialized at the beginning of the first file zip
static char liftover_snip[200];
static unsigned liftover_snip_len;
static char liftback_snip_id;

// ---------------
// ZIP & SEG stuff
// ---------------

// ZIP: called by the main thread from *_zip_initialize 
// 
void liftover_zip_initialize (DidIType dst_contig_did_i, char special_liftback_snip_id)
{
    // case: --chain : create SRCNAME context from liftover contigs and dict
    if (chain_is_loaded)
        chain_copy_dst_contigs_to_z_file (dst_contig_did_i);

    // prepare (constant) snips. note: we need these snips both for --chain and when zipping dual-coord files

    SmallContainer liftover_con = {
        .repeats     = 1,
        .nitems_lo   = 5,
        .items       = { { .dict_id = (DictId)dict_id_fields[VCF_oCHROM],  .seperator = ","  },
                         { .dict_id = (DictId)dict_id_fields[VCF_oPOS],    .seperator = ","  },
                         { .dict_id = (DictId)dict_id_fields[VCF_oREF],    .seperator = ","  },
                         { .dict_id = (DictId)dict_id_fields[VCF_oSTRAND], .seperator = ","  },
                         { .dict_id = (DictId)dict_id_fields[VCF_oALTRULE]                   } }
    };

    liftover_snip_len = sizeof (liftover_snip);
    container_prepare_snip ((Container*)&liftover_con, 0, 0, liftover_snip, &liftover_snip_len);

    liftback_snip_id = special_liftback_snip_id;
}

// SEG: map coordinates primary->luft
static LiftOverStatus liftover_get_liftover_coords (VBlockP vb, Buffer *liftover_chrom2chainsrc, 
                                                    WordIndex *dst_contig_index, PosType *dst_1pos) // out
{
    // extend liftover_chrom2chainsrc if needed
    if (vb->chrom_node_index >= liftover_chrom2chainsrc->len) { 
        buf_alloc (vb, liftover_chrom2chainsrc, 0, MAX (vb->chrom_node_index+1, 100), WordIndex, 0, "liftover_chrom2chainsrc");

        // initialize new entries allocated to WORD_INDEX_MISSING
        uint32_t size = liftover_chrom2chainsrc->size / sizeof (WordIndex);
        for (uint32_t i=liftover_chrom2chainsrc->param; i < size; i++)
            *ENT (WordIndex, *liftover_chrom2chainsrc, i) = WORD_INDEX_MISSING;

        liftover_chrom2chainsrc->param = size; // param holds the number of entries initialized >= len
        liftover_chrom2chainsrc->len   = vb->chrom_node_index + 1;
    }

    ARRAY (WordIndex, map, *liftover_chrom2chainsrc);

    // if chrom is not yet mapped to src_contig, map it now
    if (map[vb->chrom_node_index] == WORD_INDEX_MISSING) 
        map[vb->chrom_node_index] = chain_get_src_contig_index (vb->chrom_name, vb->chrom_name_len);

    *dst_1pos = 0; // initialize
    if (map[vb->chrom_node_index] == WORD_INDEX_NONE) {
        *dst_contig_index = WORD_INDEX_NONE;
        *dst_1pos = 0;
        return LO_NO_CHROM;
    }
    else 
        return chain_get_liftover_coords (map[vb->chrom_node_index], vb->last_int(DTF(pos)), dst_contig_index, dst_1pos);
}                                

// --------------------------------------------------------------
// Segging a NON-dual-coordinates file called when genozip --chain
// --------------------------------------------------------------

static char liftover_seg_get_oREF (VBlockP vb, WordIndex ochrom, PosType opos)
{
    Range *range = ref_seg_get_locked_range (vb, ochrom, opos, 1, ENT (char, vb->txt_data, vb->line_start), NULL);
    uint32_t index_within_range = opos - range->first_pos;

    ref_assert_nucleotide_available (range, opos);
    char ref = ref_get_nucleotide (range, index_within_range);

    return ref;
}

// Create either (LIFTOVER and LIFTBACK) or LIFTREJT records when aligning a NON-dual-coordinates file (using --chain)
LiftOverStatus liftover_seg_add_INFO_LIFT_fields (VBlockP vb, DidIType ochrom_did_i, char orefalt_special_snip_id,
                                                  DictId liftover_dict_id, DictId liftback_dict_id, DictId liftrejt_dict_id,
                                                  ZipDataLine *dl)
{
    LiftOverStatus ostatus = liftover_get_liftover_coords (vb, &vb->liftover, &dl->chrom_index[1], &dl->pos[1]);
    
    seg_by_did_i (vb, liftover_status_names[ostatus], strlen (liftover_status_names[ostatus]), ochrom_did_i + OSTATUS_OFFSET, 0);
    vb->last_index (ochrom_did_i + OSTATUS_OFFSET) = ostatus;

    // we add oPOS, oCHROM, oREF and oSTRAND only in case status is OK, they will be consumed by the liftover_snip container
    if (ostatus == LO_OK) {
        unsigned opos_len = str_int_len (dl->pos[1]);
        seg_pos_field (vb, VCF_oPOS, VCF_oPOS, false, true, 0, 0, dl->pos[1], opos_len);

        buf_alloc (vb, &vb->contexts[ochrom_did_i].b250, 1, vb->lines.len, WordIndex, CTX_GROWTH, "contexts->b250");
        NEXTENT (WordIndex, vb->contexts[ochrom_did_i].b250) = dl->chrom_index[1];
        
        unsigned ochrom_len = ENT (CtxNode, vb->contexts[ochrom_did_i].ol_nodes, dl->chrom_index[1])->snip_len;
        vb->contexts[ochrom_did_i].txt_len += ochrom_len;

// TODO - check with reference if REF changed
        char oref = liftover_seg_get_oREF (vb, dl->chrom_index[1], dl->pos[1]);

// TO DO handle orientation change
        unsigned ref_len = vb->last_txt_len(VCF_oREF); // vcf_seg_txt_line uses this location to store ref_len (not oREF!)
        unsigned alt_len = vb->last_txt_len(VCF_REFALT) - ref_len - 1;

        bool oref_is_ref = (ref_len == 1 && *last_txt (vb, VCF_REFALT) == oref);
        unsigned oref_len = oref_is_ref ? ref_len : alt_len; 
if (!oref_is_ref) fprintf (stderr, "xxx changed ref: before: %.*s after: %c\n", ref_len, last_txt (vb, VCF_REFALT), oref);   

// TODO - get strand - change of strand withour change in REF is no change

        char oref_special[3] = { SNIP_SPECIAL, orefalt_special_snip_id, '0' + oref_is_ref };
        seg_by_did_i (vb, oref_special, 3, ochrom_did_i + OREF_OFFSET, oref_len); 

        seg_by_did_i (vb, "+", 1, ochrom_did_i + OSTRAND_OFFSET, 1);
        #define ostrand_len 1

        seg_by_did_i (vb, "N", 1, ochrom_did_i + OALTRULE_OFFSET, 1);
        #define oaltrule_len 1

        seg_by_dict_id (vb, liftover_snip, liftover_snip_len, liftover_dict_id, 4); // account for 4 commas
        seg_by_dict_id (vb, ((char[]){ SNIP_SPECIAL, liftback_snip_id }), 2, liftback_dict_id, 0); // don't account as this container is not reconstructed by default

        // we modified the txt by adding these 5 fields to INFO/LIFTOVER. We account for them now, and we will account for the INFO name etc in seg_info_field
        vb->vb_data_size += opos_len + ochrom_len + oref_len + ostrand_len + oaltrule_len + 4 /* commas */;
    }
    else {
        dl->chrom_index[1] = WORD_INDEX_NONE;
        dl->pos[1] = 0;

        unsigned ostatus_len = strlen (liftover_status_names[ostatus]);
        seg_by_dict_id (vb, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_LIFTREJT }), 2, liftrejt_dict_id, ostatus_len);

        vb->vb_data_size += ostatus_len; // we modified the txt by adding this string to INFO/LIFTREJT - update
    }

    return ostatus;
}

// -----------------------------------------------------------------------
// Segging LIFTOVER / LIFTBACK / LIFTREJT records of dual-coordinates files
// Liftover record: CHROM,POS,REF,STRAND,ALTRULE
// Liftback record: CHROM,POS,REF,STRAND,ALTRULE
// Liftrejt snip:   REJECTION_REASON
// -----------------------------------------------------------------------

// parse INFO_LIFTOVER/BACK to its componets. note that it is destructive - it replaces the ','s with 0
static void liftover_seg_parse_record (VBlockP vb, char *value, int value_len, const char *field_name, 
                                       const char **chrom, unsigned *chrom_len,
                                       PosType *pos,       unsigned *pos_len,
                                       const char **ref,   unsigned *ref_len,
                                       char *strand,
                                       char *altrule)
{
    const char *strs[5];
    unsigned str_lens[5];
    ASSSEG (str_split (value, value_len, 5, ',', strs, str_lens, 0), value, "Invalid %s field: \"%.*s\"", field_name, value_len, value);

    *chrom     = strs[0]; 
    *chrom_len = str_lens[0];
    *pos_len   = str_lens[1];
    *ref       = strs[2];
    *ref_len   = str_lens[2];
    *strand    = strs[3][0];
    *altrule   = strs[4][0];

    ASSSEG (str_get_int_range64 (strs[1], *pos_len, 0, MAX_POS, pos), value, "Invalid POS value in liftover record: \"%s\"", strs[1]);
}

// parse the of LIFTOVER record, seg all the o* fields, and add the liftover and liftback snips
void liftover_seg_LIFTOVER (VBlockP vb, DictId liftover_dict_id, DictId liftback_dict_id,
                            DidIType ochrom_did_i, char orefalt_special_snip_id,
                            const char *ref, unsigned ref_len, // primary REF
                            const char *alt, unsigned alt_len, // optional - primary ALT
                            char *value, int value_len,
                            ZipDataLine *dl)
{
    ASSINP (!chain_is_loaded, "%s: --chain cannot be used with this file because it is already a dual-coordinates VCF file - it contains variants with the INFO/"INFO_LIFTOVER" subfield", txt_name);

    // parse record
    const char *ochrom, *oref; char ostrand, oaltrule; PosType opos; unsigned ochrom_len, opos_len, oref_len;
    liftover_seg_parse_record (vb, value, value_len, "INFO/LIFTOVER",
                               &ochrom, &ochrom_len, &opos, &opos_len, &oref, &oref_len, &ostrand, &oaltrule); 

    // oCHROM
    dl->chrom_index[1] = seg_by_did_i (vb, ochrom, ochrom_len, ochrom_did_i + OCHROM_OFFSET, ochrom_len);

    // oPOS
    dl->pos[1] = seg_pos_field (vb, ochrom_did_i+1, ochrom_did_i + OPOS_OFFSET, false, true, 0, 0, opos, opos_len);

    // oREF
    bool oref_is_ref =        (ref_len==oref_len && str_case_compare (ref, oref, oref_len, NULL)); // alleles are case insensitive per VCF spec
    bool oref_is_alt = alt && (alt_len==oref_len && str_case_compare (alt, oref, oref_len, NULL));
    
    // currently, genozip only supports unchanged REF or switched REF/ALT, variants with other type REF changes should have been rejected
    ASSSEG (oref_is_ref || oref_is_alt, value, "expecting OREF=\"%.*s\" to be the same as either REF=\"%.*s\" or ALT=\"%.*s\"",
            oref_len, oref, ref_len, ref, alt_len, alt ? alt : "");

    char oref_special[3] = { SNIP_SPECIAL, orefalt_special_snip_id, '0' + oref_is_ref };
    seg_by_did_i (vb, oref_special, 3, ochrom_did_i + OREF_OFFSET, oref_len); 

    // oSTRAND, oALTRULE and oSTATUS
    seg_by_did_i (vb, &ostrand,  1, ochrom_did_i + OSTRAND_OFFSET,  1);
    seg_by_did_i (vb, &oaltrule, 1, ochrom_did_i + OALTRULE_OFFSET, 1);
    seg_by_did_i (vb, liftover_status_names[LO_OK], strlen (liftover_status_names[LO_OK]), ochrom_did_i + OSTATUS_OFFSET, 0); // 0 bc doesn't reconstruct

    // LIFTOVER and LIFTBACK container snips (INFO container filter will determine which is reconstructed)
    seg_by_dict_id (vb, liftover_snip, liftover_snip_len, liftover_dict_id, 4); // account for 4 commas
    seg_by_dict_id (vb, ((char[]){ SNIP_SPECIAL, liftback_snip_id }), 2, liftback_dict_id, 0); // don't account as this container is not reconstructed by default

    vb->last_index (ochrom_did_i + OSTATUS_OFFSET) = LO_OK;
}

// parse the of LIFTBACK record, seg all the primary coordinate fields, generate oREF, and add the liftover and liftback snips
void liftover_seg_LIFTBACK (VBlockP vb, DictId liftover_dict_id, DictId liftback_dict_id,
                            DidIType ochrom_did_i, DidIType pos_did_i, char orefalt_special_snip_id,
                            const char *oref, unsigned oref_len, 
                            const char *oalt, unsigned oalt_len, // optional
                            void (*seg_ref_alt_cb)(VBlockP vb, const char *ref, unsigned ref_len, const char *alt, unsigned alt_len), // call back to seg primary REF and ALT
                            char *value, int value_len,
                            ZipDataLine *dl)
{
    ASSINP (!chain_is_loaded, "%s: --chain cannot be used with this file because it is already a dual-coordinates VCF file - it contains variants with the INFO/"INFO_LIFTBACK" subfield", txt_name);

    // parse record
    const char *chrom, *ref; char strand, altrule; PosType pos; unsigned pos_len, ref_len, chrom_len;
    liftover_seg_parse_record (vb, value, value_len, "INFO/LIFTBACK", 
                               &chrom, &chrom_len, &pos, &pos_len, &ref, &ref_len, &strand, &altrule); 

    // CHROM
    dl->chrom_index[0] = seg_chrom_field (vb, chrom, chrom_len);

    // POS
    dl->pos[0] = seg_pos_field (vb, pos_did_i, pos_did_i, false, true, 0, 0, pos, pos_len+1);
            
    random_access_update_pos (vb, pos_did_i);

    // REFALT & OREF
    bool ref_is_oref =         (ref_len==oref_len && str_case_compare (ref, oref, ref_len, NULL)); // alleles are case insensitive per VCF spec
    bool ref_is_oalt = oalt && (ref_len==oalt_len && str_case_compare (ref, oalt, ref_len, NULL));
    
    // currently, genozip only supports unchanged REF or switched REF/ALT, variants with other type REF changes should have been rejected
    ASSSEG (ref_is_oref || ref_is_oalt, value, "expecting REF=\"%.*s\" to be the same as either OREF=\"%.*s\" or OALT=\"%.*s\"",
            ref_len, ref, oref_len, oref, oalt_len, oalt ? oalt : "");

    // REFALT (primary coords) - to be interpreted the REFALT translator that will translate it to oREF and oALT using oref_special function
    if (ref_is_oref) 
        seg_ref_alt_cb (vb, oref, oref_len, oalt, oalt_len); // also accounts for 2 tabs
    else  // switch ref and alt
        seg_ref_alt_cb (vb, oalt, oalt_len, oref, oref_len); // also accounts for 2 tabs

    // oREF (we couldn't seg it in the main field, as it refers to REF)
    char oref_special[3] = { SNIP_SPECIAL, orefalt_special_snip_id, '0' + ref_is_oref }; // should be the same as in vcf_seg_initialize()
    seg_by_did_i (vb, oref_special, 3, ochrom_did_i + OREF_OFFSET, ref_is_oref ? oref_len : oalt_len); 

    // oSTRAND, oALTRULE and oSTATUS 
    seg_by_did_i (vb, &strand,  1, ochrom_did_i + OSTRAND_OFFSET,  1); // +1 because displayed by default (in LIFTOVER record)
    seg_by_did_i (vb, &altrule, 1, ochrom_did_i + OALTRULE_OFFSET, 1);
    seg_by_did_i (vb, liftover_status_names[LO_OK], strlen (liftover_status_names[LO_OK]), ochrom_did_i + OSTATUS_OFFSET, 0); // not displayed by default

    // LIFTOVER and LIFTBACK container snips (INFO container filter will determine which is reconstructed)
    seg_by_dict_id (vb, liftover_snip, liftover_snip_len, liftover_dict_id, 4); // account for 4 commas
    seg_by_dict_id (vb, ((char[]){ SNIP_SPECIAL, liftback_snip_id }), 2, liftback_dict_id, 0); // don't account as this container is not reconstructed by default

    vb->last_index (ochrom_did_i + OSTATUS_OFFSET) = LO_OK;
}

void liftover_seg_LIFTREJT (VBlockP vb, DictId dict_id, DidIType ochrom_did_i, const char *value, int value_len,
                            ZipDataLine *dl)
{
    ASSINP (!chain_is_loaded, "%s: --chain cannot be used with this file because it is already a dual-coordinates VCF file - it contains variants with the INFO/"INFO_LIFTREJT" subfield", txt_name);

    dl->chrom_index[1] = WORD_INDEX_NONE;
    dl->pos[1] = 0;

    // get status from string
    SAFE_NUL (&value[value_len]);
    
    for (unsigned i=0; i < NUM_LO_STATUSES; i++)
        if (!strcmp (value, liftover_status_names[i])) { // found

            unsigned status_len = strlen (value);
            seg_by_did_i (vb, value, status_len, ochrom_did_i + OSTATUS_OFFSET, status_len);
            seg_by_dict_id (vb, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_LIFTREJT }), 2, dict_id_INFO_LIFTREJT, 0); // 0 bc no comma and name is accounted for by INFO

            SAFE_RESTORE;

            vb->last_index (ochrom_did_i + OSTATUS_OFFSET) = i;

            return;
        }

    ASSSEG (false, value, "Invalid oSTATUS name: %.*s", value_len, value);
}


// ----------------------------------------------------------------
// Rejects file stuff - file contains data in the native txt format
// ----------------------------------------------------------------

// ZIP: called when inspecting the txtheader to add header data, and after each VB to add rejected line
void liftover_append_rejects_file (VBlockP vb)
{
    ASSERTNOTNULL (z_file);

    // create rejects file if not already open
    if (!z_file->rejects_file) {
        z_file->rejects_file_name = malloc (strlen (z_file->name) + 20);
        sprintf (z_file->rejects_file_name, "%s.rejects%s", z_file->name, file_plain_ext_by_dt (z_file->data_type));
        
        z_file->rejects_file = fopen (z_file->rejects_file_name, "wb+");
        ASSERT (z_file->rejects_file, "fopen() failed to open %s: %s", z_file->rejects_file_name, strerror (errno));
    }

    ASSERT0 (z_file->rejects_file, "liftover rejects file is not open");

    fwrite (vb->liftover_rejects.data, 1, vb->liftover_rejects.len, z_file->rejects_file);
    z_file->rejects_disk_size += vb->liftover_rejects.len;

    buf_free (&vb->liftover_rejects);
}

