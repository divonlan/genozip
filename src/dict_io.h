// ------------------------------------------------------------------
//   dict_io.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

extern void dict_io_read_all_dictionaries (void);
extern void dict_io_compress_dictionaries (void);

extern StrTextMegaLong str_snip_ex (STRp(snip), bool add_quote);
#define str_snip str_snip_ex (snip, snip_##len, true).s
extern void dict_io_print (FILE *fp, STRp(data), bool with_word_index, bool add_quotation_marks, bool add_newline, bool remove_equal_asterisk);
extern void dict_io_show_singletons (VBlockP vb, ContextP ctx);
