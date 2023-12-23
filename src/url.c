// ------------------------------------------------------------------
//   url.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifndef _WIN32
#include <sys/socket.h>
#include <fcntl.h>
#endif

#define Z_LARGE64
#ifdef __APPLE__
    #define off64_t __int64_t
#endif
#include "genozip.h"
#include "url.h"
#include "stream.h"
#include "file.h"
#include "strings.h"

// our instance of curl for url_open - only one at a time is permitted
static StreamP curl = NULL;

bool url_is_url (rom filename)
{
    // in a URL, immediately following a sequence of letters, we're expecting "://". check this.
    for (rom c = filename; *c && c < (filename + 8); c++) // relevant URI schemes are http, https, ftp, file... 8 charactersnto be on the safe side (see: https://en.wikipedia.org/wiki/List_of_URI_schemes)
        if (!IS_LETTER (*c)) 
            return c > filename && c[0]==':' && c[1]=='/' && c[2]=='/' && c[3];

    return false; // all letters, no ://
}

#define CURL_RESPONSE_LEN 4096

// returns error string if curl itself (not server) failed, or NULL if successful
static void url_do_curl (rom url, char *stdout_data, unsigned *stdout_len/*in/out*/,
                         char *error, unsigned *error_len,
                         bool blocking, bool follow_redirects) 
{
    // our own instance of curl - to not conflict with url_open
    StreamP curl = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 0, 0,
                                  "To read from a URL", // reason in case of failure to execute curl
                                  "curl", 
                                  flag.is_windows  ? "--ssl-no-revoke" : SKIP_ARG, 
                                  follow_redirects ? "--location"      : SKIP_ARG,
                                  url, NULL); 

    FILE *data_stream  = stream_from_stream_stdout (curl);
    FILE *error_stream = stream_from_stream_stdout (curl);

#ifndef _WIN32
    // set non-blocking (unfortunately _pipes in Windows are always blocking)
    if (!blocking)
        ASSERT (!fcntl(fileno(data_stream), F_SETFL, fcntl (fileno(data_stream), F_GETFL) | O_NONBLOCK), "fcntl failed: %s", strerror (errno));
#endif

    #define RETRIES 70
    for (int i=0; i < RETRIES; i++) {
        int ret = fread (stdout_data, 1, *stdout_len - 1, data_stream); // -1 to leave room for \0
        if (ret == 0) {
            if (errno == EAGAIN && i < RETRIES-1) { // doesn't happen on Windows
                usleep (100000); // 100ms
                continue;
            }

            strcpy (error, "Failed to read (possibly no connection)");
            *error_len = strlen (error);
            *stdout_len = 0;

            stream_close (&curl, STREAM_DONT_WAIT_FOR_PROCESS);
            return;
        }
    
        *stdout_len = (unsigned)ret;
        break;
    }

    stdout_data[*stdout_len] = '\0'; // terminate string

    *error_len = fread (error, 1, CURL_RESPONSE_LEN-1, error_stream); 
    error[*error_len] = '\0'; // terminate string

    int curl_exit_code = stream_close (&curl, STREAM_DONT_WAIT_FOR_PROCESS); // Don't wait for process - Google Forms call hangs if waiting
    if (!curl_exit_code) 
        return; // curl itself is good - we may have or not an error in "error" from the server or in case of no connection

#ifdef _WIN32
    if (curl_exit_code == ENFILE) return; // we didn't read all the data on the pipe - that's ok
#endif

    // case: for non-HTTP urls (eg ftp:// file://) or for HTTP urls where the error occurred before connecting
    // to the webserver (eg bad url) the error comes in stderr, and curl exit code is non-0.
    // curl message format looks like this: "curl: (3) URL using bad/illegal format or missing URL"

    // get message itself
    char *msg_start = strstr (error, ") ");
    if (msg_start) { 
    
        msg_start += 2;

        char *msg_after = strchr (msg_start, '\r'); // in case of Windows-style end-of-line
        if (!msg_after) msg_after = strchr (msg_start, '\n'); // try again, with \n

        if (msg_start && msg_after) {
            memcpy (error, msg_start, msg_after - msg_start);
            error[msg_after - msg_start] = '\0';
        }
    }
    else 
        strcpy (error, strerror (curl_exit_code));

    *error_len = strlen (error);
}

static int url_do_curl_head (rom url,
                             char *stdout_data, unsigned *stdout_len,
                             char *stderr_data, unsigned *stderr_len)
{
    ASSERT (stream_is_exec_in_path ("curl"),
            "Failed to open URL %s because curl was not found in the execution path", url);

    // our own instance of curl - to not conflict with url_open
    StreamP curl = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 0, 0,
                                  "To compress files from a URL",
                                  "curl", 
                                  flag.is_windows ? "--ssl-no-revoke" : SKIP_ARG,                          
                                  "--head", url, NULL); // not silent - we want to collect errors

    *stdout_len = fread (stdout_data, 1, CURL_RESPONSE_LEN-1, stream_from_stream_stdout (curl));
    stdout_data[*stdout_len] = '\0'; // terminate string

    *stderr_len = fread (stderr_data, 1, CURL_RESPONSE_LEN-1, stream_from_stream_stderr (curl));
    stderr_data[*stderr_len] = '\0'; // terminate string

    int curl_exit_code = stream_close (&curl, STREAM_WAIT_FOR_PROCESS);

    return curl_exit_code;
}

// for a url, returns whether that URL exists, and if possible, its file_size, or -1 if its not available
// note that the file_size availability is at the discretion of the web/ftp site. 
// in case of an error, returns the error string : caller should FREE() the error string
rom url_get_status (rom url, thool *is_file_exists, int64_t *file_size)
{
    *is_file_exists = unknown;
    *file_size = 0;
    char response[CURL_RESPONSE_LEN];
    char *error = MALLOC (CURL_RESPONSE_LEN);
    unsigned response_len, error_len;

    // run 'curl --head url'
    int curl_exit_code = url_do_curl_head (url, response, &response_len, error, &error_len);
    
    // case: for non-HTTP urls (eg ftp:// file://) or for HTTP urls where the error occurred before connecting
    // to the webserver (eg bad url) the error comes in stderr, and curl exit code is non-0.
    // curl message format looks like this: "curl: (3) URL using bad/illegal format or missing URL"
    if (curl_exit_code) {
        // get message itself
        char *msg_start = strstr (error, ") ");
        if (msg_start) msg_start += 2;

        char *msg_after = strchr (msg_start, '\r'); // in case of Windows-style end-of-line
        if (!msg_after) msg_after = strchr (msg_start, '\n'); // try again, with \n

        if (msg_start && msg_after) {
            memcpy (error, msg_start, msg_after - msg_start);
            error[msg_after - msg_start] = '\0';
        }
        return error;
    }

    // Get the first line of text - since curl completed successfully, we expect this to always exist
    char *first_eol = strchr (response, '\r'); // in case of Windows-style end-of-line
    if (!first_eol) first_eol = strchr (response, '\n'); // try again, with \n

    // case: we got exit_code=0 (OK) but no response. this happens, for example, in some FTP servers when the URL is a directory
    if (!first_eol) {
        strcpy (error, "Server did not respond to curl -I. Please check the URL");
        return error;
    }

    // for HTTP urls, the status (success or error) is the first line. For now, we will copy the first line regardless of url type
    strncpy (error, response, first_eol - response);
    error[first_eol - response] = '\0';

    // case HTTP: check status
    if (strstr (error, "http") || strstr (error, "HTTP")) {
        if (strstr (error, "200")) 
            *is_file_exists = true;
        else if (strstr (error, "404")) // note: some servers respond with 403 ("Forbidden") to a --head request, but allow downloading the file
            *is_file_exists = false;
        else
            return error;
    } 
        
    rom len_start = NULL;
    if      ((len_start = strstr (response, "content-length:"))) len_start += STRLEN("content-length:");
    else if ((len_start = strstr (response, "Content-Length:"))) len_start += STRLEN("Content-Length:");

    // Case: we got the file length - file exists even if we didn't get an HTTP status (eg because URL is not http)
    if (len_start) {
        *file_size = strtoull (len_start, NULL, 10);
        *is_file_exists = true; // if we got its length, we know it exists (for non-HTTP, this is the way we check existence)
    }
    else *file_size = -1; // file possibly exists (if this is HTTP and response was 200), but length is unknown

    FREE (error);
    return NULL; // no error
}

// reads a string response from a URL, returns a nul-terminated string and the number of characters (excluding \0), or -1 if failed
int32_t url_read_string (rom url, STRc(data), bool blocking, bool follow_redirects)
{
    char *response, local_response[CURL_RESPONSE_LEN], error[CURL_RESPONSE_LEN];
    unsigned response_len=0, error_len=0;

    int url_len = strlen (url);
    for (int i=0, in_arg=false; i < url_len ; i++) {
        ASSERT (!in_arg || IS_VALID_URL_CHAR(url[i]) || url[i]=='&' || url[i]=='=' || url[i]=='%', 
                "Invalid url character [%u]=%c(%u). url=\"%s\"", i, url[i], (unsigned char)url[i], url); 
        if (url[i] == '?') in_arg = true;
    }

    response     = data ? data     : local_response;
    response_len = data ? data_len : CURL_RESPONSE_LEN;
    url_do_curl (url, response, &response_len, error, &error_len, blocking, follow_redirects);

    if (error_len && !response_len) return -1; // failure

    if (response_len && strstr (response, "Bad Request"))
        return -2;

    return data ? response_len : 0;
}

bool url_get_redirect (rom url, STRc(redirect_url))
{
    if (!stream_is_exec_in_path ("curl")) return false;

    // our own instance of curl - to not conflict with url_open
    StreamP curl = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 0, 0,
                                  "To get a URL's redirect", // reason in case of failure to execute curl
                                  "curl", url, "--location", "--silent", "--output", "/dev/null",
                                  flag.is_windows ? "--ssl-no-revoke" : SKIP_ARG,                          
                                  "--write-out", "%{url_effective}", NULL); 

    redirect_url_len = fread (redirect_url, 1, redirect_url_len - 1, stream_from_stream_stdout (curl));
    redirect_url[redirect_url_len] = '\0'; // terminate string

    stream_close (&curl, STREAM_WAIT_FOR_PROCESS);

    return true;
}

// returns a FILE* which streams the content of a URL 
// Note: FILE* returned is a *copy* of the FILE* in the curl stream - it should not be FCLOSEd
// directly, rather call url_kill_curl or url_disconnect_from_curl
FILE *url_open (StreamP parent_stream, rom url)
{
    ASSERT0 (!curl, "Error url_open failed because curl is already running");

    char str5[6] = "";
    strncpy (str5, url, 5);
    bool is_file = str_case_compare (str5, "file:", NULL); // wget doesn't support file://

    // check if wget exists, it is better than curl in flakey connections
    bool has_wget = !flag.is_windows && stream_is_exec_in_path ("wget");

    ASSERT (has_wget || stream_is_exec_in_path ("curl"),
            "Failed to open URL %s because %s not found in the execution path", 
            url, (flag.is_windows || is_file) ? "curl was" : "wget or curl were");

    if (has_wget && !is_file) // note: stream_create doesn't support soft_fail (yet)
        curl = stream_create (parent_stream, DEFAULT_PIPE_SIZE, 0, 0, 0, 0, 0,
                              "To compress files from a URL", "wget", "--tries=16", "--quiet", "--waitretry=3", "--output-document=/dev/stdout", 
                              url, NULL);
    else // curl
        curl = stream_create (parent_stream, DEFAULT_PIPE_SIZE, 0, 0, 0, 0, 0,
                              "To compress files from a URL", "curl", "--silent", 
                              flag.is_windows ? "--ssl-no-revoke" : SKIP_ARG,                          
                              url, NULL);

    return stream_from_stream_stdout (curl);
}

void url_reset_if_curl (StreamP maybe_curl_stream)
{
    if (curl == maybe_curl_stream) curl = NULL;
}

// kill curl - used in case of an
void url_kill_curl (FILE **copy_of_input_pipe)
{
    if (!curl) return; // nothing to do
    
    if (copy_of_input_pipe) *copy_of_input_pipe = NULL;
    stream_close (&curl, STREAM_KILL_PROCESS);
}

// close stream without killing curl process - used if stream is used by forked sub-process
void url_disconnect_from_curl (FILE **copy_of_input_pipe)
{
    if (!curl) return; // nothing to do
    
    if (copy_of_input_pipe) *copy_of_input_pipe = NULL;
    stream_close (&curl, STREAM_DONT_WAIT_FOR_PROCESS);

}

// make a string into a a string containing only valid url characters, eg "first last" -> "first%20last"
char *url_esc_non_valid_chars_(rom in, char *out/*malloced if NULL*/, bool esc_all_or_none)
{
    if (!out) out = MALLOC (strlen(in) * 3 + 1);
    char *next = out;
    rom save_in = in;
    bool any_invalid=false;

    for (; *in; in++) {
        if (IS_VALID_URL_CHAR(*in)) {
            if (esc_all_or_none) {
                sprintf (next, "%%%02X", (unsigned char)*in);
                next +=3;
            }
            else
                *(next++) = *in;
        }

        else {
            sprintf (next, "%%%02X", (unsigned char)*in);
            next +=3;
            any_invalid = true;
        }
    }
    *next = '\0';

    if (esc_all_or_none && !any_invalid)
        strcpy (out, save_in);

    return out;
}

UrlStr url_esc_non_valid_charsS (rom in) // for short strings - on stack
{
    rom esc = url_esc_non_valid_chars_(in, NULL, false); // note: might be longer than UrlStr
    
    UrlStr out;
    int out_len = MIN_(sizeof (out.s)-1, strlen(esc)); // trim if needed - possibly resulting in an invalid URL!
    memcpy (out.s, esc, out_len);
    out.s[out_len] = 0;
    
    FREE (esc);
    return out;
}
