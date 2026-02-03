// ------------------------------------------------------------------
//   dict_io.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

extern void dict_io_read_all_dictionaries (void);
extern void dict_io_compress_dictionaries (void);

extern StrTextMegaLong str_snip_ex (DataType dt, STRp(snip), bool add_quote);
#define str_snip str_snip_ex (DT_NONE, snip, snip_##len, true).s
extern void dict_io_print (FILE *fp, DictId dict_id, STRp(data), bool with_word_index, bool add_quotation_marks, bool add_newline, bool remove_non_contigs);
extern void dict_io_show_singletons (VBlockP vb, ContextP ctx);
