// ------------------------------------------------------------------
//   url.h
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include <stdio.h>
#include "genozip.h"
#include "stream.h"

extern bool url_is_url (rom filename);

extern rom url_get_status (rom url, bool *is_file_exists, int64_t *file_size);

extern FILE *url_open (StreamP parent_stream, rom url);
extern void url_reset_if_curl (StreamP maybe_curl_stream);

extern int32_t url_read_string (rom url, char *data, uint32_t data_size);

extern void url_kill_curl (void);

extern char *url_esc_non_valid_chars (rom in);
