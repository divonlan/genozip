// ------------------------------------------------------------------
//   user_message.c
//   Copyright (C) 2023-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "file.h"
#include "user_message.h"
#include "zfile.h"
#include "license.h"

static Buffer user_message = {};

void user_message_init (rom filename)
{
    file_get_file (evb, filename, &user_message, "user_message", 0, VERIFY_UTF8, false);
}

bool has_user_msg (void)
{
    return user_message.len > 0;
}

// ZIP main thread: write to global section
void user_message_compress (void)
{
    if (user_message.len) 
        zfile_compress_section_data (evb, SEC_USER_MESSAGE, &user_message);
}

// PIZ main thread
void user_message_display (void)
{
    if (flag.quiet || flag.test) return;

    Section sec = sections_first_sec (SEC_USER_MESSAGE, SOFT_FAIL);
    if (!sec) return; // no user message in this file

    zfile_get_global_section (SectionHeader, sec, &user_message, "user_message");

    iprintf ("%.*s", STRfb(user_message));
}
