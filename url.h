// ------------------------------------------------------------------
//   url.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef URL_INCLUDED
#define URL_INCLUDED

#include <stdio.h>
#include "genozip.h"
#include "stream.h"

extern bool url_is_url (const char *filename);

extern const char *url_get_status (const char *url, bool *file_exists, int64_t *file_size);

extern FILE *url_open (StreamP parent_stream, const char *url);

extern void url_kill_curl (void);

#endif
