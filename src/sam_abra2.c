// ------------------------------------------------------------------
//   sam_abra2.c
//   Copyright (C) 2024-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#include "container.h"


#define _YA_RNAME  DICT_ID_MAKE2_8("Y0A_RNAM")  
#define _YA_POS    DICT_ID_MAKE2_7("Y1A_POS")    
#define _YA_CIGAR  DICT_ID_MAKE2_8("Y2A_CIGR")  

#define _YO_RNAME  DICT_ID_MAKE2_8("Y0O_RNAM")  
#define _YO_POS    DICT_ID_MAKE2_7("Y1O_POS")    
#define _YO_STRAND DICT_ID_MAKE2_8("Y2O_SRND") 
#define _YO_CIGAR  DICT_ID_MAKE2_8("Y3O_CIGR")  

sSTRl (YO_NA_con_snip, 100);

void sam_abra2_zip_initialize (void)
{
    MiniContainer container_YO_NA = { .nitems_lo = 1, .repeats = 1, .items[0].dict_id = { _YO_POS } };

    container_prepare_snip ((ContainerP)&container_YO_NA, 
                            (char[]){ CON_PX_SEP, 'N', '/', 'A', ':', CON_PX_SEP }, 6, 
                            qSTRa (YO_NA_con_snip));
}

void sam_abra2_seg_initialize (VBlockSAMP vb)
{
    ctx_get_ctx (vb, _YA_CIGAR)->ltype = LT_STRING;
    ctx_get_ctx (vb, _YO_CIGAR)->ltype = LT_STRING;
}

// Contig alignment info: "rname:pos:cigar"
void sam_seg_ABRA2_YA_Z (VBlockSAMP vb, STRp(field), unsigned add_bytes)
{
    static const MediumContainer container_YA = { .nitems_lo = 3,   
                                                  .repeats   = 1,       
                                                  .items     = { { .dict_id = { _YA_RNAME }, .separator = {':'} }, // note: no need to alias dict to RNAME, bc its expected to mostly be "copy from RNAME", not contigs 
                                                                 { .dict_id = { _YA_POS   },   .separator = {':'} },  
                                                                 { .dict_id = { _YA_CIGAR },                    } } };

    SegCallback callbacks[6] = { sam_seg_0A_rname_cb, sam_seg_0A_pos_cb, sam_seg_0A_cigar_cb };
     
    seg_struct (VB, CTX(OPTION_YA_Z), container_YA, STRa(field), callbacks, add_bytes, true);
}

// "rname:pos:orientation:cigar" OR "N/A:pos" (TO DO: find a way to mux between them)
void sam_seg_ABRA2_YO_Z (VBlockSAMP vb, STRp(field), unsigned add_bytes)
{
    decl_ctx (OPTION_YO_Z);

    if (field_len > 4 && field[0]=='N' && field[1]=='/' && field[2]=='A' && field[3]==':') {
        sam_seg_0A_pos_cb (VB, ctx_get_ctx (vb, _YO_POS), field+4, field_len-4, 0);
        seg_by_ctx (VB, STRa(YO_NA_con_snip), ctx, add_bytes - (field_len-4));
    }

    else {
        static const MediumContainer container_YO = { .nitems_lo = 4,   
                                                      .repeats   = 1,       
                                                      .items     = { { .dict_id = { _YO_RNAME  },  .separator = {':'} }, // note: no need to alias dict to RNAME, bc its expected to mostly be "copy from RNAME", not contigs 
                                                                     { .dict_id = { _YO_POS    },    .separator = {':'} },  
                                                                     { .dict_id = { _YO_STRAND }, .separator = {':'} },  
                                                                     { .dict_id = { _YO_CIGAR  },                     } } };

        SegCallback callbacks[6] = { sam_seg_0A_rname_cb, sam_seg_0A_pos_cb, 0, sam_seg_0A_cigar_cb };
        
        seg_struct (VB, ctx, container_YO, STRa(field), callbacks, add_bytes, true);
    }
}
