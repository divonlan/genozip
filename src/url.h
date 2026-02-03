// ------------------------------------------------------------------
//   url.h
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include <stdio.h>
#include "genozip.h"
#include "stream.h"
#include "buffer.h"

// remote file stream
extern FILE *url_open_remote_file  (StreamP parent_stream, rom url);
extern void url_reset_if_remote_file_stream (StreamP maybe_remote_file_stream);
extern void url_close_remote_file_stream (FILE **copy_of_input_pipe);
extern void url_disconnect_from_remote_file_stream (FILE **copy_of_input_pipe);

// access URLs
extern rom url_get_status (rom url, thool *is_file_exists, int64_t *file_size);
extern int32_t url_read_string (rom url, STRc(STRc), bool blocking, bool follow_redirects, rom show_errors);
extern bool url_get_redirect (rom url, STRc(redirect_url), StreamP *redirect_stream);

// string operations
extern bool url_is_url (rom filename);

extern char *url_esc_non_valid_chars_(rom in, char *out, bool esc_all_or_none);
static inline char *url_esc_non_valid_chars (rom in) { return url_esc_non_valid_chars_ (in, NULL, false); } // on heap
static inline char *url_esc_all_or_none (rom in)     { return url_esc_non_valid_chars_ (in, NULL, true ); } // on heap

extern StrTextLong url_esc_non_valid_charsS (rom in); // for short strings - on stack
