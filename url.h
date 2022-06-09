// ------------------------------------------------------------------
//   url.h
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include <stdio.h>
#include "genozip.h"
#include "stream.h"
#include "buffer.h"

extern bool url_is_url (rom filename);

extern rom url_get_status (rom url, bool *is_file_exists, int64_t *file_size);

extern FILE *url_open (StreamP parent_stream, rom url);
extern void url_reset_if_curl (StreamP maybe_curl_stream);

extern int32_t url_read_string (rom url, char *data, uint32_t data_size);

extern void url_get_redirect (rom url, STRc(redirect_url));

extern void url_kill_curl (void);

extern char *url_esc_non_valid_chars_(rom in, char *out);
static inline char *url_esc_non_valid_chars (rom in) { return url_esc_non_valid_chars_ (in, NULL); } // on heap

typedef struct { char s[256]; } UrlStr;
static inline UrlStr url_esc_non_valid_charsS (rom in) // for short strings - on stack
{
    rom esc = url_esc_non_valid_chars_(in, NULL); // note: might be longer than UrlStr
    
    UrlStr out;
    int out_len = MIN_(sizeof (out.s)-1, strlen(esc)); // trim if needed
    memcpy (out.s, esc, out_len);
    out.s[out_len] = 0;
    
    FREE (esc);
    return out;
}
