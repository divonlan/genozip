// ------------------------------------------------------------------
//   file.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#define Z_LARGE64
#ifdef __APPLE__
    #define off64_t __int64_t
#endif
#include "genozip.h"
#include "url.h"
#include "stream.h"
#include "file.h"
#include "strings.h"

// our instance of curl - only one at a time is permitted
static StreamP curl = NULL;

bool url_is_url (const char *filename)
{
    return !!strstr (filename, "://");
}

#define CURL_RESPONSE_LEN 4096

// returns error string if curl itself (not server) failed, or NULL if successful
static void url_do_curl (const char *url, bool head,
                         char *stdout_data, unsigned *stdout_len,
                         char *error, unsigned *error_len,
                         const char *reason) // optional text in case curl execution fails
{
    curl = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 0,
                          reason ? reason : "To read from a URL",
                          "curl", url, 
                          head ? "--head" : NULL,
                          NULL); // not silent - we want to collect errors

    *stdout_len = read (fileno (stream_from_stream_stdout (curl)), stdout_data, CURL_RESPONSE_LEN-1);
    stdout_data[*stdout_len] = '\0'; // terminate string

    *error_len = read (fileno (stream_from_stream_stderr (curl)), error, CURL_RESPONSE_LEN-1);
    error[*error_len] = '\0'; // terminate string

    int curl_exit_code = stream_close (&curl, STREAM_WAIT_FOR_PROCESS);
    if (!curl_exit_code) 
        return; // curl itself is good - we may have or not an error in "error" from the server

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

static int url_do_curl_head (const char *url,
                             char *stdout_data, unsigned *stdout_len,
                             char *stderr_data, unsigned *stderr_len)
{
    curl = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 0,
                          "To compress files from a URL",
                          "curl", "--head", url, NULL); // not silent - we want to collect errors

    *stdout_len = fread (stdout_data, 1, CURL_RESPONSE_LEN-1, stream_from_stream_stdout (curl));
    stdout_data[*stdout_len] = '\0'; // terminate string

    *stderr_len = fread (stderr_data, 1, CURL_RESPONSE_LEN-1, stream_from_stream_stderr (curl));
    stderr_data[*stderr_len] = '\0'; // terminate string

    int curl_exit_code = stream_close (&curl, STREAM_WAIT_FOR_PROCESS);

    return curl_exit_code;
}

// for a url, returns whether that URL exists, and if possible, its file_size, or -1 if its not available
// note that the file_size availability is at the discretion of the web/ftp site. 
// in case of an error, returns the error string
const char *url_get_status (const char *url, bool *file_exists, int64_t *file_size)
{
    *file_exists = false;
    *file_size = 0;
    char response[CURL_RESPONSE_LEN];
    char *error = malloc (CURL_RESPONSE_LEN);
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
            *file_exists = true;
        else
            return error;
    } 
        
    const char *len_start = NULL;
    if      ((len_start = strstr (response, "content-length:"))) len_start += strlen ("content-length:");
    else if ((len_start = strstr (response, "Content-Length:"))) len_start += strlen ("Content-Length:");

    // Case: we got the file length - file exists even if we didn't get an HTTP status (eg because URL is not http)
    if (len_start) {
        *file_size = strtoull (len_start, NULL, 10);
        *file_exists = true; // if we got its length, we know it exists (for non-HTTP, this is the way we check existence)
    }
    else *file_size = -1; // file possibly exists (if this is HTTP and response was 200), but length is unknown

    return NULL; // no error
}


// reads a string response from a URL, returns a null-terminated string and the number of characters (excluding \0)
uint32_t url_read_string (const char *url, char *data, uint32_t data_size,
                          const char *reason) // optional text in case curl execution fails
{
    char response[CURL_RESPONSE_LEN], error[CURL_RESPONSE_LEN];
    unsigned response_len, error_len;

    url_do_curl (url, false, response, &response_len, error, &error_len, reason);

    ASSERT (!error_len || response_len, "Error accessing the Internet: %s", error);

    if (data) {
        unsigned len = MIN (data_size-1, response_len);
        memcpy (data, response, len);
        data[len] = '\0';
        return len;
    }
    else return 0;
}

// returns a FILE* which streams the content of a URL
FILE *url_open (StreamP parent_stream, const char *url)
{
    ASSERT0 (!curl, "Error url_open failed because curl is already running");

    curl = stream_create (parent_stream, global_max_memory_per_vb, 0, 0, 0, 0, 
                          "To compress files from a URL", "curl", "--silent", url, NULL);
    return stream_from_stream_stdout (curl);
}


// kill curl - used in case of an
void url_kill_curl (void)
{
    if (!curl) return; // nothing to do
    
    stream_close (&curl, STREAM_KILL_PROCESS);
}

// make a string into a a string containing only valid url characters, eg "first last" -> "first%20last"
char *url_esc_non_valid_chars (const char *in)
{
    char *out = malloc (strlen(in) * 3 + 1);
    char *next = out;

    for (; *in; in++) {
        if (IS_VALID_URL_CHAR(*in)) 
            *(next++) = *in;

        else {
            sprintf (next, "%%%02X", *in);
            next +=3;
        }
    }
    *next = '\0';

    return out;
}