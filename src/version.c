// ------------------------------------------------------------------
//   version.c
//   Copyright (C) 2022-2023 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <pthread.h>

#include "genozip.h"
#include "url.h"
#include "version.h"
#include "flags.h"
#include "strings.h"
#include "arch.h"

static rom latest_version = NULL;
static bool thread_running = false;
static pthread_t thread_id;
static char redirect_url[128];

// thread entry
static void *version_background_test_for_newer_do (void *unused)
{
    thread_running = true;
    pthread_setcanceltype (PTHREAD_CANCEL_ASYNCHRONOUS, NULL); // thread can be canceled at any time

    if (url_get_redirect (GITHUB_LATEST_RELEASE, redirect_url, sizeof redirect_url)) {
        latest_version = strrchr (redirect_url, '-'); // url looks something like: "https://github.com/divonlan/genozip/releases/tag/genozip-13.0.18"
        if (latest_version) {
            latest_version++;
            if (!flag.debug_latest && !strcmp (latest_version, GENOZIP_CODE_VERSION))
                latest_version = NULL; // same as current version
        } 
    }

    thread_running = false;
    return NULL;
}

// tests for a newer version on github in a background threads, and sets the url_has_newer_version global variable
void version_background_test_for_newer (void)
{
    // case: current version ends with "-1\0" - its the next major release dev version - already the newest
    if (!flag.debug_latest && GENOZIP_CODE_VERSION[sizeof GENOZIP_CODE_VERSION-3] == '-') return; 

    pthread_create (&thread_id, NULL, version_background_test_for_newer_do, NULL); // ignore errors
}

void version_print_notice_if_has_newer (void)
{
    // case: Genozip finished its work while thread is still running - kill it
    if (thread_running) 
        pthread_cancel (thread_id);
    
    // case: thread completed, and there is a new version
    else if (latest_version) {
        iprintf ("\nA newer & better version of Genozip is available - version %s. You are currently running version %s\n", 
                 latest_version, GENOZIP_CODE_VERSION);

        if (is_info_stream_terminal) {
            // note: distribution names are defined in the Makefile
            if (!strcmp (arch_get_distribution(), "InstallForge")) 
                iprintf ("You can install the latest version from here: %s\n", GITHUB_WINDOWS_INSTALLER);
            
            else if (!strcmp (arch_get_distribution(), "linux-x86_64"))
                iprintf ("You can download the latest version from here: %s\n", GITHUB_LINUX_TARBALL);

#ifndef _WIN32
            else if (!strcmp (arch_get_distribution(), "conda") &&
                     str_query_user_yn ("Do you want to upgrade now?", QDEF_YES))
                system ("conda update genozip");
#endif
            else 
                iprintf ("Installation instructions: %s\n", WEBSITE_INSTALLING);
        }
        flag.no_tip = true; // printed - no more tips
    }
}

