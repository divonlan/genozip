// ------------------------------------------------------------------
//   txtfile.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef TXTFILE_INCLUDED
#define TXTFILE_INCLUDED

#include "genozip.h"

//extern bool txtfile_get_line (VBlockP vb, unsigned line_i_in_file, bool skip_md5_vcf_header, BufferP line, const char *buf_name);

extern void txtfile_read_vcf_header (bool is_first_vcf);
extern void txtfile_read_variant_block (VBlockP vb);
extern unsigned txtfile_write_to_disk (ConstBufferP buf);
extern void txtfile_estimate_txt_data_size (VBlockP vb);

extern void txtfile_write_one_vblock_vcf (VBlockVCFP vb);
extern void txtfile_write_one_vblock_sam (VBlockSAMP vb);

#endif