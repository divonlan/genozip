// ------------------------------------------------------------------
//   version.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <pthread.h>
#include <errno.h>
#include <libgen.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include "genozip.h"
#include "url.h"
#include "version.h"
#include "flags.h"
#include "strings.h"
#include "arch.h"
#include "website.h"
#include "file.h"

static rom latest_version = NULL;
static bool thread_running = false;
static pthread_t thread_id;
static char redirect_url[128];

// version of the genozip executable running
int code_version_major (void)
{
    return (GENOZIP_CODE_VERSION[0]-'0') * 10 + (GENOZIP_CODE_VERSION[1]-'0');
}

int code_version_minor (void)
{
    rom ver = GENOZIP_CODE_VERSION;
    return atoi (strrchr (ver, '.') + 1);
}

StrText version_str (void)
{
    StrText s={};

    if (IS_ZIP || !z_file)  
        sprintf (s.s, "version=%s", GENOZIP_CODE_VERSION);
    else if (!VER2(15,28))
        sprintf (s.s, "code_version=%s file_version=%u", GENOZIP_CODE_VERSION, z_file->genozip_version);
    else        
        sprintf (s.s, "code_version=%s file_version=%u.0.%u", GENOZIP_CODE_VERSION, z_file->genozip_version, z_file->genozip_minor_ver);

    return s;
}

rom genozip_update_msg (void)
{
    return "Please update Genozip to the latest version. See: " WEBSITE_INSTALLING;
}

// thread entry
static void *version_background_test_for_newer_do (void *unused)
{
    thread_running = true;
    pthread_setcanceltype (PTHREAD_CANCEL_ASYNCHRONOUS, NULL); // thread can be canceled at any time

    if (url_get_redirect (GITHUB_LATEST_RELEASE, redirect_url, sizeof redirect_url)) {
        latest_version = strrchr (redirect_url, '-'); // url looks something like: "https://github.com/divonlan/genozip/releases/tag/genozip-13.0.18"
        if (latest_version) {
            latest_version++;

            char copy[100];
            strcpy (copy, GENOZIP_CODE_VERSION); // bc str_split_ints doesn't work on string literals

            str_split_ints (latest_version, strlen (latest_version), 3, '.', latest, true);
            str_split_ints (copy, strlen (copy), 3, '.', this, true);

            if (!flag.debug_latest && n_latests == 3 && n_thiss == 3 && latests[0] == thiss[0] && latests[2] <= thiss[2])
                latest_version = NULL; // same or newer than current version
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

#ifndef _WIN32
static void my_system (rom command)
{
    #pragma GCC diagnostic push 
    #pragma GCC diagnostic ignored "-Wpragmas"         // avoid warning if "-Wuse-after-free" is not defined in this version of gcc
    #pragma GCC diagnostic ignored "-Wunused-result"   // avoid compiler warning due to not checking system() return value

    int ret = system (command);
    ASSERT (ret != -1, "Failed to create child process for: \"%s\": %s", command, strerror (ret));
    ASSERT (!ret, "command \"%s\" failed: status=%d", command, ret);

    #pragma GCC diagnostic pop
}
#endif

#ifdef __linux__
static bool is_updatable (void)
{
    // check that executable is "genozip" (not genozip-debug or renamed to anything else). note: even in piz with --test, executable is "genozip"
    rom exe_base = basename (arch_get_executable().s); // note: basename might modify string
    if (strcmp (exe_base, "genozip"))
        return false; 

    // check that executable directory is writeable
    char *dir_name = dirname (arch_get_executable().s);
    if (access (dir_name, W_OK)) return false;

    return true;
}

static void udpate_do (void)
{
    // cd to /tmp
    rom save_cwd = getcwd (0, 0);
    ASSERT (save_cwd, "getcwd failed: %s", strerror (errno));

    ASSERT (!chdir ("/tmp"), "chdir(\"/tmp\") failed: %s", strerror (errno));

    // download
    my_system ("rm -f " TARBALL_NAME);
    my_system ("wget --quiet " GITHUB_LINUX_TARBALL);

    // untar
    char cmd1[100 + STRLEN(TARBALL_NAME)];
    sprintf (cmd1, "tar xf %s", TARBALL_NAME);
    my_system (cmd1);

    // mv to destination 
    rom dir_name = dirname (arch_get_executable().s);
    char cmd2[strlen(dir_name) + 100];

    rom exe_names[] = { "genozip", "genounzip", "genocat", "genols" };
    for (int exe_i=0; exe_i < ARRAY_LEN(exe_names); exe_i++) {    
        // note: mv better than cp because it preserves the hard links ; mv better than rename because it works across filesystems
        sprintf (cmd2, "mv genozip-linux-x86_64/%s %s", exe_names[exe_i], dir_name);
        my_system (cmd2);
    }

    // cleanup
    sprintf (cmd1, "rm -Rf genozip-linux-x86_64 %s", TARBALL_NAME); // note: "genozip-linux-x86_64" is defined in the Makefile
    my_system (cmd1);

    // revert current directory
    ASSERT (!chdir (save_cwd), "chdir(\"%s\") failed: %s", save_cwd, strerror (errno));
    FREE (save_cwd);

    iprint0 ("Genozip has been updated.\n");
}

#elif defined _WIN32
// static bool is_updatable (void)
// {
//     // check that executable is "genozip.exe" (not genozip-debug.exe or renamed to anything else). note: even in piz with --test, executable is "genozip.exe"
//     rom exe_base = basename (arch_get_executable().s); // note: basename might modify string
//     if (strcmp (exe_base, "genozip.exe"))
//         return false; 

//     return true;
// }

// static void udpate_do (void)
// {
//     rom dir_name = dirname (arch_get_executable().s);
//     int dir_name_len = strlen (dir_name);

//     char src[dir_name_len + 100], dst[dir_name_len + 100], old[dir_name_len + 100];

//     sprintf (src, "%s/genozip.exe", dir_name);
//     sprintf (old, "%s/genozip-old.exe", dir_name); // cannot delete running executable
//     DeleteFile (old); // perhaps from previous upgrade - ignore errors
    
//     ASSERT (MoveFile (src, old), "Failed to MoveFile %s to %s: %s", src, old, str_win_error());

//     // DOESNT WORK: any attempt to create an exe file - either by downloading it, or by downloading
//     // to another extension and then CopyFile or MoveFile to an .exe - no error is returned, but exe is not created.
//     sprintf (dst, "%s/%s", dir_name, WINDOWS_UPDATE_NAME);
//     DeleteFile (dst); 
//     ASSERT (URLDownloadToFile (NULL, GITHUB_WINDOWS_UPDATE, dst, 0, NULL) == S_OK,
//             "Failed to download %s", GITHUB_WINDOWS_UPDATE);

//     ASSERT (MoveFile (dst, src), "Failed to MoveFile %s to %s: %s", dst, src, str_win_error());

//     sprintf (dst, "%s/LICENSE.txt", dir_name);
//     ASSERT (URLDownloadToFile (NULL, GITHUB_LICENSE_TXT, dst, 0, NULL) == S_OK,
//             "Failed to download %s", GITHUB_LICENSE_TXT);

//     rom exe_names[] = { "genounzip", "genocat", "genols" };
//     for (int exe_i=0; exe_i < ARRAY_LEN (exe_names); exe_i++) {
//         // rename current file to -old
//         sprintf (old, "%s/%s-old.exe", dir_name, exe_names[exe_i]); // cannot delete running executable
//         DeleteFile (old); // perhaps from previous upgrade - ignore errors
//         ASSERT (MoveFile (src, old), "Failed to MoveFile %s to %s: %s", src, old, str_win_error());
        
//         sprintf (dst, "%s/genounzip.exe", dir_name);
//         ASSERT (CopyFile (src, dst, false), "Failed to CopyFile %s to %s: %s", src, dst, str_win_error());
//     }    

//     iprint0 ("Genozip has been updated.\n");
// }
#endif

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
            if (0) {}
#ifdef _WIN32
            // else if (is_updatable() &&                      
            //     str_query_user_yn ("Do you want to update Genozip now?", QDEF_YES)) 
            //     udpate_do();
            
#elif defined __linux__            
            else if (strcmp (arch_get_distribution(), "conda") && // anything but conda. note: distribution names are defined in the Makefile
                is_updatable() &&                      
                str_query_user_yn ("Do you want to update Genozip now?", QDEF_YES)) 
                udpate_do();
#endif            

#ifndef _WIN32
            else if (!strcmp (arch_get_distribution(), "conda") &&
                     str_query_user_yn ("Do you want to update Genozip now?", QDEF_YES)) {
                my_system ("conda update genozip");
            }
#endif
            else if (!strcmp (arch_get_distribution(), "InstallForge")) 
                iprintf ("You can install the latest version from here: %s\n", GITHUB_WINDOWS_INSTALLER);

            else 
                iprintf ("Installation instructions: %s\n", WEBSITE_INSTALLING);
        }
        flag.no_tip = true; // printed - no more tips
    }
}

