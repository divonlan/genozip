// ------------------------------------------------------------------
//   strings.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "strings.h"
#include "dict_id.h"
#ifndef WIN32
#include <sys/ioctl.h>
#endif

void str_to_lowercase (char *s)
{
    // lowercase argv[0] to allow case-insensitive comparison in Windows
    for (; *s; s++) 
        if (IS_CLETTER (*s))
            *s += 'a' - 'A';
}

char *str_size (int64_t size, char *str /* out */)
{
    if      (size > (1LL << 40)) sprintf (str, "%3.1lf TB", ((double)size) / (double)(1LL << 40));
    else if (size > (1LL << 30)) sprintf (str, "%3.1lf GB", ((double)size) / (double)(1LL << 30));
    else if (size > (1LL << 20)) sprintf (str, "%3.1lf MB", ((double)size) / (double)(1LL << 20));
    else if (size > (1LL << 10)) sprintf (str, "%3.1lf KB", ((double)size) / (double)(1LL << 10));
    else                         sprintf (str, "%3d B"    ,     (int)size)                       ;

    return str; // for convenience so caller can use in printf directly
}

char *str_int (int64_t n, char *str /* out */, unsigned *len)
{
    *len=0;

    if (n==0) {
        str[0] = '0';
        *len=1;
    }

    else {
        bool is_negative = (n<0);
        if (is_negative) n = -n;

        char rev[50] = {}; // "initialize" to avoid compiler warning
        while (n) {
            rev[(*len)++] = '0' + n % 10;
            n /= 10;
        }
        // now reverse it
        for (int i=0; i < (*len); i++) str[i + is_negative] = rev[(*len)-i-1];

        if (is_negative) {
            str[0] = '-';
            (*len)++;
        }
    }

    str[(*len)] = '\0'; // string terminator
    return str;
}

char *str_uint_commas (int64_t n, char *str /* out */)
{
    unsigned len = 0, orig_len=0;

    if (n==0) {
        str[0] = '0';
        len = 1;
    }
    else {
        char rev[50] = {}; // "initialize" to avoid compiler warning
        while (n) {
            if (orig_len && orig_len % 3 == 0) rev[len++] = ',';    
            rev[len++] = '0' + n % 10;
            orig_len++;
            n /= 10;
        }
        // now reverse it
        for (int i=0; i < len; i++) str[i] = rev[len-i-1];
    }

    str[len] = '\0'; // string terminator
    return str;
}

#define POINTER_STR_LEN 19
char *str_pointer (const void *p, char *str /* POINTER_STR_LEN bytes allocated by caller*/)
{
#ifdef _MSC_VER
    sprintf (str, "0x%I64x", (uint64_t)p);
#else
    sprintf (str, "0x%"PRIx64, (uint64_t)p);
#endif
    return str;
}

const char *type_name (unsigned item, 
                       const char * const *name, // the address in which a pointer to name is found, if item is in range
                       unsigned num_names)
{
    if (item > num_names) {
        static char str[50];
        sprintf (str, "%d (out of range)", item);
        return str;
    }
    
    return *name;    
}

int str_print_text (const char **text, unsigned num_lines,
                    const char *wrapped_line_prefix, 
                    const char *newline_separator, 
                    unsigned line_width /* 0=calcuate optimal */)
{                       
    if (!line_width) {
#ifdef _WIN32
        line_width = 120; // default width of cmd window in Windows 10
#else
        // in Linux and Mac, we can get the actual terminal width
        struct winsize w;
        ioctl(0, TIOCGWINSZ, &w);
        line_width = MAX (40, w.ws_col); // our wrapper cannot work with to small line widths 
#endif
    }

    for (unsigned i=0; i < num_lines; i++)  {
        const char *line = text[i];
        unsigned line_len = strlen (line);
        
        // print line with line wraps
        bool wrapped = false;
        while (line_len + (wrapped ? strlen (wrapped_line_prefix) : 0) > line_width) {
            int c; for (c=line_width-1 - (wrapped ? strlen (wrapped_line_prefix) : 0); 
                        c>=0 && (IS_LETTER(line[c]) || IS_DIGIT (line[c])); // wrap lines at - and | too, so we can break very long regex strings like in genocat
                        c--); // find 
            printf ("%s%.*s\n", wrapped ? wrapped_line_prefix : "", c, line);
            line += c + (line[c]==' '); // skip space too
            line_len -= c + (line[c]==' ');
            wrapped = true;
        }
        printf ("%s%s%s", wrapped ? wrapped_line_prefix : "", line, newline_separator);
    }
    return 0;
}

// receives a user response, a default "Y" or "N" (or NULL) and modifies the response to be "Y" or "N"
bool str_verify_y_n (char *response, unsigned response_size, const char *y_or_n)
{
    ASSERT0 (!y_or_n || (strlen (y_or_n)==1 && (y_or_n[0]=='Y' || y_or_n[0]=='N')), 
             "Error: y_or_n needs to be NULL, \"Y\" or \"N\"");

    // default is N (or no default) and first character of the user's response is y or Y
    if ((!y_or_n || y_or_n[0] == 'N') && (response[0] == 'y' || response[0] == 'Y')) response[0] = 'Y'; 
    
    // default is Y (or no default) and first character of the user's response is n or N
    else if ((!y_or_n || y_or_n[0] == 'Y') && (response[0] == 'n' || response[0] == 'N')) response[0] = 'N';

    // we have a default, and the response of user is not opposite of the default - return default
    else if (y_or_n) response[0] = y_or_n[0]; 

    // we don't have a default - we request the user to respond again
    else return false;

    response[1] = 0;

    return true; // always accept the response
}

bool str_verify_not_empty (char *response, unsigned response_size, const char *unused)
{ 
    unsigned len = strlen (response);

    return !(len==0 || (len==1 && response[0]=='\r')); // not \n or \r\n only
}

void str_query_user (const char *query, char *response, unsigned response_size, 
                     ResponseVerifier verifier, const char *verifier_param)
{
    do {
        fprintf (stderr, "%s", query);

        unsigned bytes = read (STDIN_FILENO, response, response_size); 
        response[bytes-1] = '\0'; // string terminator instead of the newline

    } while (verifier && !verifier (response, response_size, verifier_param));
}
