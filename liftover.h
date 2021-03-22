// ------------------------------------------------------------------
//   liftover.h
//   Copyright (C) 2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef LIFTOVER_INCLUDED
#define LIFTOVER_INCLUDED

#include "genozip.h"

// this is part of the file format (in the oSTATUS context) and CANNOT BE CHANGED (but can be extended)
typedef enum {
    // Algorithm interim statuses not outputed to the genozip file
    LO_UNKNOWN      = -2, // Status of this Liftover line not known yet
    LO_NONE         = -1, // not a Liftover item

    // Rejection reasons defined in the "Dual coordinate VCF" spec 
    LO_OK           = 0, 
    LO_NO_CHROM     = 1,  // chain file doesn't contain a qName which is the primary chrom as 
    LO_NO_MAPPING   = 2,  // chain file doesn't contain a mapping for the primary POS
    LO_REF_SPLIT    = 3,  // A multi-base REF does not map to one contiguous sequence in the laft reference
    LO_END_MAPPING  = 4,  // INFO/END is on a different alignment than POS
    LO_UNSUPPORTED  = 5,  // REF changed in a way that genozip can't handle
    NUM_LO_STATUSES
} LiftOverStatus;
extern const char *liftover_status_names[];

#define LIFTOVER_STATUS_NAMES { "OK", "NO_CHROM", "NO_MAPPING", "REF_SPLIT", "END_MAPPING", "UNSUPPORTED", "LINE_TOO_LONG" }
 
extern void liftover_copy_data_from_chain_file (void);
extern void liftover_zip_initialize (DidIType tname_did_i, char special_liftback_snip_id);
extern const char *liftover_get_laft_contig (uint32_t contig_i);

// Seg
extern void liftover_append_rejects_file (VBlockP vb);
extern LiftOverStatus liftover_seg_add_INFO_LIFT_fields (VBlockP vb, DidIType ochrom_did_i, char orefalt_special_snip_id, DictId liftover_dict_id, DictId liftback_dict_id, DictId liftrejd_dict_id, ZipDataLineP dl);
extern void liftover_seg_LIFTOVER (VBlockP vb, DictId liftover_dict_id, DictId liftback_dict_id, DidIType ochrom_did_i, char orefalt_special_snip_id, const char *ref, unsigned ref_len, const char *alt, unsigned alt_len, char *value, int value_len, ZipDataLineP dl);
extern void liftover_seg_LIFTBACK (VBlockP vb, DictId liftover_dict_id, DictId liftback_dict_id, DidIType ochrom_did_i, DidIType pos_did_i, char orefalt_special_snip_id, const char *oref, unsigned oref_len, const char *oalt, unsigned oalt_len, void (*seg_ref_alt_cb)(VBlockP, const char *, unsigned, const char *, unsigned), char *value, int value_len, ZipDataLineP dl);
extern void liftover_seg_LIFTREJD (VBlockP vb, DictId dict_id, DidIType ochrom_did_i, const char *value, int value_len, ZipDataLineP dl);

// PIZ: section list stuff
extern void liftover_section_list_remove_rejects (void);        // called if no --laft
extern void liftover_section_list_move_rejects_to_front (void); // called if --laft


#endif
