// ------------------------------------------------------------------
//   url.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
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
#include "arch.h"

#define CURL_RESPONSE_LEN 4096

static StreamP remote_file_stream = NULL;

bool url_is_url (rom filename)
{
    // in a URL, immediately following a sequence of letters, we're expecting "://". check this.
    for (rom c = filename; *c && c < (filename + 8); c++) // relevant URI schemes are http, https, ftp, file... 8 charactersnto be on the safe side (see: https://en.wikipedia.org/wiki/List_of_URI_schemes)
        if (!IS_LETTER (*c)) 
            return c > filename && c[0]==':' && c[1]=='/' && c[2]=='/' && c[3];

    return false; // all letters, no ://
}

static inline bool is_wget (StreamP stream) { return stream_exec_is (stream, "wget"); }

static StreamP url_open (StreamP parent_stream, rom url, bool head_only)
{
    char str5[6] = "";
    strncpy (str5, url, 5);
    bool is_file = str_case_compare (str5, "file:", NULL); 

    // wget is better than curl in flakey connections
    if (!is_file &&                       // wget doesn't support file://
        !(head_only && curl_available())  // wget --spider doesn't follow redirects, so for header_only, we prefer curl if its available
        && wget_available())              // note: wget is not supported for Windows (see wget_available)
        
        return stream_create (parent_stream, 
                              head_only ? 0 : DEFAULT_PIPE_SIZE, // in wget, header arrives in error channel and data channel is empty
                              DEFAULT_PIPE_SIZE, 0, 0, 0, 0,
                              "To compress files from a URL", "wget", "--tries=16", 
                              "--quiet", 
                              "--waitretry=3", 
                              head_only ? "--server-response"           : SKIP_ARG,                         
                              head_only ? "--spider"                    : SKIP_ARG,                         
                              head_only ? "--output-document=/dev/null" : "--output-document=/dev/stdout", 
                              url, NULL);
    
    else if (curl_available()) 
        return stream_create (parent_stream, DEFAULT_PIPE_SIZE, 0, 0, 0, 0, 0,
                              "To compress files from a URL", "curl", "--silent", "--location",
                              flag.is_windows ? "--ssl-no-revoke" : SKIP_ARG,
                              head_only       ? "--head"          : SKIP_ARG,                         
                              url, NULL);

    else 
        ABORT ("Failed to open URL %s because %s not found in the execution path", 
                url, (flag.is_windows || is_file) ? "curl was" : "wget or curl were");
}

static bool url_get_head (rom url, qSTRp(data), qSTRp(error),
                          StreamP *url_stream) // set while stream is open, so other threads can kill it if needed
{
    store_release (*url_stream, (StreamP)NULL);

    store_release (*url_stream, url_open (NULL, url, true));
    
    bool wget = is_wget (*url_stream);

    FILE *f;
    if (wget) {
        *data_len  = 0;
        *error_len = (f=stream_from_stream_stderr (*url_stream)) ? fread (error, 1, *error_len - 1, f)   : 0;
    }
    else {
        *data_len  = (f=stream_from_stream_stdout (*url_stream)) ? fread (data,  1, *data_len  - 1, f) : 0;
        *error_len = (f=stream_from_stream_stderr (*url_stream)) ? fread (error, 1, *error_len - 1, f) : 0;
    }

    int exit_code = stream_close (url_stream, STREAM_WAIT_FOR_PROCESS);

    if ((flag.is_windows && exit_code == ENFILE) || // we didn't read all the data on the pipe - that's ok
        (wget && exit_code == 1)) // not really an error
        exit_code = 0; 

    if (wget && exit_code && ! *error_len) {
        // see: https://www.gnu.org/software/wget/manual/html_node/Exit-Status.html
        static rom wget_exit_errors[] = { 
            [1]="Generic error code",
            [2]="Parse error",
            [3]="File I/O error",
            [4]="Network failure",
            [5]="SSL verification failure",
            [6]="Authentication failure",
            [7]="Protocol error",
            [8]="Server error" 
        };

        // mimic curls' format: "curl: (3) URL using bad/illegal format or missing URL\n"
        if (exit_code < ARRAY_LEN(wget_exit_errors))
            *error_len = snprintf (error, CURL_RESPONSE_LEN, "wget: (%u) %s\n", exit_code, wget_exit_errors[exit_code]);
        else
            *error_len = snprintf (error, CURL_RESPONSE_LEN, "wget: (%u) Unknown error\n", exit_code);
    }

    if (!wget && exit_code && ! *error_len) { 
        // see: https://curl.se/libcurl/c/libcurl-errors.html
        static rom some_curl_exit_errors[] = { // curl has 101 error codes... we display common ones, and refer to the manual for the rest
            [1]="Unsupported protocol",
            [3]="URL malformat",
            [5]="Couldn't resolve proxy",
            [6]="Couldn't resolve host",
            [7]="Couldn't connect",
            [9]="Remote access denied",
            [10]="FTP access denied",
            [12]="FTP ACCEPT timeout", 
            [22]="HTTP error",
        };

        if (exit_code < ARRAY_LEN(some_curl_exit_errors) && some_curl_exit_errors[exit_code])
            *error_len = snprintf (error, CURL_RESPONSE_LEN, "curl: (%u) %s\n", exit_code, some_curl_exit_errors[exit_code]);
        else
            *error_len = snprintf (error, CURL_RESPONSE_LEN, "curl: (%u) see: https://curl.se/libcurl/c/libcurl-errors.html\n", exit_code);
    }

    // in wget, header arrives in the error channel - move it to data if successful
    if (wget && !exit_code) {        
        memcpy (data, error, *error_len); 
        *data_len = *error_len;
        *error_len = 0;
    }

    data[*data_len] = error[*error_len] = '\0'; // terminate strings

    return !exit_code; // true if successful 
}

bool url_get_redirect (rom url, STRc(redirect_url),
                       StreamP *redirect_stream) // optional out: used to free the stream if the running thread is canceled. 
{
    if (!curl_available() && !wget_available()) return false;

    STRlic (response, CURL_RESPONSE_LEN);
    STRlic (error, CURL_RESPONSE_LEN);

    url_get_head (url, qSTRa(response), qSTRa(error), redirect_stream);

    // we're looking for a line that looks like this (URL is terminated by some white space)
    // Location: https://github.com/divonlan/genozip/releases/tag/genozip-15.0.57
    #define LOCATION "Location: "
    char *location = strstr (response, LOCATION);
    if (!location) return false;
    location += STRLEN (LOCATION);

    uint32_t location_len = strcspn (location, " \n\r\t");
    if (location_len >= redirect_url_len) return false;

    memcpy (redirect_url, location, location_len); 
    redirect_url[location_len] = 0;

    return true;
}

// for a url, returns whether that URL exists, and if possible, its file_size, or -1 if its not available
// note that the file_size availability is at the discretion of the web/ftp site. 
// in case of an error, returns the error string : caller should FREE() the error string
rom url_get_status (rom url, thool *is_file_exists, int64_t *file_size)
{
    *is_file_exists = unknown;
    *file_size = 0;
    STRlic (response, CURL_RESPONSE_LEN);
    char *error = MALLOC (CURL_RESPONSE_LEN);
    uint32_t error_len = CURL_RESPONSE_LEN;
    
    // case: for non-HTTP urls (eg ftp:// file://) or for HTTP urls where the error occurred before connecting
    // to the webserver (eg bad url) the error comes in stderr, and curl exit code is non-0.
    // curl message format looks like this: "curl: (3) URL using bad/illegal format or missing URL"
    StreamP url_stream = 0;
    if (!url_get_head (url, qSTRa(response), qSTRa(error), &url_stream)) {
        // get message itself
        char *msg_after=0, *msg_start = strstr (error, ") ");
        if (msg_start) {
            msg_start += 2;
            
            msg_after = strchr (msg_start, '\r'); // in case of Windows-style end-of-line
            if (!msg_after) msg_after = strchr (msg_start, '\n'); // try again, with \n
        
            if (msg_after) {
                memcpy (error, msg_start, msg_after - msg_start);
                error[msg_after - msg_start] = '\0';
            }
        }
        
        *is_file_exists = false;
        return error;
    }

    // Get the first line of text - since curl/wget completed successfully, we expect this to always exist
    char *first_eol = strchr (response, '\r'); // in case of Windows-style end-of-line
    if (!first_eol) first_eol = strchr (response, '\n'); // try again, with \n

    // case: we got exit_code=0 (OK) but no response. this happens, for example, in some FTP servers when the URL is a directory
    if (!first_eol) {
        strcpy (error, "Server did not respond to request for the URL head. Please check the URL");
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

// returns error string if curl/wget itself (not server) failed, or NULL if successful
static void url_read_string_do (rom url, qSTRp(data), qSTRp(error), bool blocking, bool follow_redirects) 
{
    StreamP url_stream = url_open (NULL, url, false);
    FILE *data_stream  = stream_from_stream_stdout (url_stream);
    FILE *error_stream = stream_from_stream_stdout (url_stream);

#ifndef _WIN32
    // set non-blocking (unfortunately _pipes in Windows are always blocking)
    if (!blocking)
        ASSERT (!fcntl(fileno(data_stream), F_SETFL, fcntl (fileno(data_stream), F_GETFL) | O_NONBLOCK), "fcntl failed: %s", strerror (errno));
#endif

    #define RETRIES 70
    for (int i=0; i < RETRIES; i++) {
        int ret = fread (data, 1, *data_len - 1, data_stream); // -1 to leave room for \0
        if (ret == 0) {
            if (errno == EAGAIN && i < RETRIES-1) { // doesn't happen on Windows
                usleep (300000); // 300ms
                continue;
            }

            strcpy (error, "Failed to read (possibly no connection)");
            *error_len = strlen (error);
            *data_len = 0;

            stream_close (&url_stream, STREAM_DONT_WAIT_FOR_PROCESS);
            return;
        }
    
        *data_len = (unsigned)ret;
        break;
    }

    data[*data_len] = '\0'; // terminate string

    *error_len = fread (error, 1, CURL_RESPONSE_LEN-1, error_stream); 
    error[*error_len] = '\0'; // terminate string

    int exit_code = stream_close (&url_stream, STREAM_DONT_WAIT_FOR_PROCESS); // Don't wait for process - Google Forms call hangs if waiting
    if (!exit_code) 
        return; // curl/wget itself is good - we may have or not an error in "error" from the server or in case of no connection

#ifdef _WIN32
    if (exit_code == ENFILE) return; // we didn't read all the data on the pipe - that's ok
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
        strcpy (error, strerror (exit_code));

    *error_len = strlen (error);
}

// reads a string response from a URL, returns a nul-terminated string and the number of characters (excluding \0), or -1 if failed
int32_t url_read_string (rom url, STRc(data), bool blocking, bool follow_redirects, rom show_errors)
{
    rom action = show_errors ? show_errors : "url_read_string";

    if (!wget_available() && !curl_available()) {
        if (flag.debug_submit || show_errors) 
            fprintf (stderr, "\nError in %s: neither curl nor wget are available\n", action);
        
        return -3; // failure
    }

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
    url_read_string_do (url, qSTRa(response), qSTRa(error), blocking, follow_redirects);

    if (error_len && !response_len) {
        if (flag.debug_submit || show_errors) 
            fprintf (stderr, "\nError in %s: %.*s\n", action, STRf(error));

        return -1; // failure
    }

    if (response_len && strstr (response, "Bad Request")) {
        if (flag.debug_submit || show_errors) 
            fprintf (stderr, "\nError in %s: Bad Request\n", action);

        return -2;
    }

    return data ? response_len : 0;
}

// returns a FILE* which streams the content of a URL 
// Note: FILE* returned is a *copy* of the FILE* in the url_stream - it should not be FCLOSEd
// directly, rather call url_close_remote_file_stream or url_disconnect_from_remote_file_stream
// Note: only one open URL file is possible at any time, as it goes into the global "curl".
//       This is reserved for reading a remote file and should not be used for other purposes.
FILE *url_open_remote_file (StreamP parent_stream, rom url)
{
    ASSERTISNULL (remote_file_stream);
    remote_file_stream = url_open (parent_stream, url, false);

    return stream_from_stream_stdout (remote_file_stream);
}

void url_reset_if_remote_file_stream (StreamP maybe_remote_file_stream)
{
    if (remote_file_stream == maybe_remote_file_stream) remote_file_stream = NULL;
}

// kill curl/wget 
void url_close_remote_file_stream (FILE **copy_of_input_pipe)
{
    if (!remote_file_stream) return; // nothing to do
    
    if (copy_of_input_pipe) *copy_of_input_pipe = NULL;

    stream_close (&remote_file_stream, STREAM_KILL_PROCESS);
}

// close stream without killing curl/wget process - used if stream is used by forked sub-process
void url_disconnect_from_remote_file_stream (FILE **copy_of_input_pipe)
{
    if (!remote_file_stream) return; // nothing to do
    
    if (copy_of_input_pipe) *copy_of_input_pipe = NULL;
    stream_close (&remote_file_stream, STREAM_DONT_WAIT_FOR_PROCESS);
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
                snprintf (next, 4, "%%%02X", (unsigned char)*in); // 4 bc inc. room for \0 required by snprintf
                next += 3;
            }
            else
                *(next++) = *in;
        }

        else {
            snprintf (next, 4, "%%%02X", (unsigned char)*in);
            next += 3;
            any_invalid = true;
        }
    }
    *next = '\0';

    if (esc_all_or_none && !any_invalid)
        strcpy (out, save_in);

    return out;
}

StrTextLong url_esc_non_valid_charsS (rom in) // for short strings - on stack
{
    rom esc = url_esc_non_valid_chars_(in, NULL, false); // note: might be longer than StrTextLong
    
    StrTextLong out;
    int out_len = MIN_(sizeof (out.s)-1, strlen(esc)); // trim if needed - possibly resulting in an invalid URL!
    memcpy (out.s, esc, out_len);
    out.s[out_len] = 0;
    
    FREE (esc);
    return out;
}
