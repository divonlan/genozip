// ------------------------------------------------------------------
//   url.h
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef URL_INCLUDED
#define URL_INCLUDED

#include <stdio.h>
#include "genozip.h"
#include "stream.h"

extern bool url_is_url (const char *filename);

extern const char *url_get_status (const char *url, bool *is_file_exists, int64_t *file_size);

extern FILE *url_open (StreamP parent_stream, const char *url);
extern void url_reset_if_curl (StreamP maybe_curl_stream);

extern int32_t url_read_string (const char *url, char *data, uint32_t data_size);

extern void url_kill_curl (void);

extern char *url_esc_non_valid_chars (const char *in);

#endif
