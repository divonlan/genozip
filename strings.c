// ------------------------------------------------------------------
//   strings.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <time.h>
#include "genozip.h"
#include "strings.h"
#include "flags.h"
#ifndef WIN32
#include <sys/ioctl.h>
#else
#include <windows.h>
#endif

char *str_tolower (const char *in, char *out /* out allocated by caller - can be the same as in */)
{
    char *startout = out;

    for (; *in; in++, out++) 
        *out = (*in >= 'A' && *in <= 'Z') ? *in + 32 : *in;
    
    return startout;
}

char *str_toupper (const char *in, char *out /* out allocated by caller - can be the same as in */)
{
    char *startout = out;

    for (; *in; in++, out++) 
        *out = (*in >= 'a' && *in <= 'z') ? *in - 32 : *in;
    
    return startout;
}

StrText char_to_printable (char c) 
{
    switch ((uint8_t)c) {
        case 32 ... 126 : return (StrText) { .s = {c, 0} };   // printable ASCII
        case '\t'       : 
        case '\n'       : 
        case '\r'       : return (StrText) { .s = {' ', 0} }; // whitespace
        default         : { // unprintable - output eg \xf 
            StrText p;
            sprintf (p.s, "\\x%x", (uint8_t)c);
            return p;
        }
    }
}

// replaces \t, \n, \r, \b with "\t" etc, replaces unprintables with '?'. caller should allocate out. returns out.
// out should be allocated by caller to (in_len*2 + 1), out is null-terminated
char *str_to_single_line_printable (const char *in, unsigned in_len, char *out)
{
    char *start = out;

    for (unsigned i=0; i < in_len; i++)
        switch (in[i]) {
            case '\t' : *out++ = '\\'; *out++ = 't'; break;
            case '\n' : *out++ = '\\'; *out++ = 'n'; break;
            case '\r' : *out++ = '\\'; *out++ = 'r'; break;
            case '\b' : *out++ = '\\'; *out++ = 'b'; break;
            case -128 ... 7: case 11 ... 12: case 14 ... 31: 
                        *out++ = '?';                break;
            default:    *out++ = in[i];
        }
    
    *out = 0;
    return start;
}

StrText str_size (uint64_t size)
{
    StrText s;

    if      (size >= (1LL << 50)) sprintf (s.s, "%3.1lf PB", ((double)size) / (double)(1LL << 50));
    else if (size >= (1LL << 40)) sprintf (s.s, "%3.1lf TB", ((double)size) / (double)(1LL << 40));
    else if (size >= (1LL << 30)) sprintf (s.s, "%3.1lf GB", ((double)size) / (double)(1LL << 30));
    else if (size >= (1LL << 20)) sprintf (s.s, "%3.1lf MB", ((double)size) / (double)(1LL << 20));
    else if (size >= (1LL << 10)) sprintf (s.s, "%3.1lf KB", ((double)size) / (double)(1LL << 10));
    else if (size >  0          ) sprintf (s.s, "%3d B"    ,     (int)size)                       ;
    else                          sprintf (s.s, "-"                       )                       ;

    return s;
}

StrText str_bases (uint64_t num_bases)
{
    StrText s;

    if      (num_bases >= 1000000000000000) sprintf (s.s, "%.1lf Pb", ((double)num_bases) / 1000000000000000.0);
    else if (num_bases >= 1000000000000)    sprintf (s.s, "%.1lf Tb", ((double)num_bases) / 1000000000000.0);
    else if (num_bases >= 1000000000)       sprintf (s.s, "%.1lf Gb", ((double)num_bases) / 1000000000.0);
    else if (num_bases >= 1000000)          sprintf (s.s, "%.1lf Mb", ((double)num_bases) / 1000000.0);
    else if (num_bases >= 1000)             sprintf (s.s, "%.1lf Kb", ((double)num_bases) / 1000.0);
    else if (num_bases >  0   )             sprintf (s.s, "%u b",    (unsigned)num_bases);
    else                                    sprintf (s.s, "-");

    return s;
}

// returns length
unsigned str_int (int64_t n, char *str /* out */)
{
    unsigned len=0;

    if (n==0) {
        str[0] = '0';
        len=1;
    }

    else {
        bool is_negative = (n<0);
        if (is_negative) n = -n;

        char rev[50] = {}; // "initialize" to avoid compiler warning
        while (n) {
            rev[len++] = '0' + n % 10;
            n /= 10;
        }
        // now reverse it
        for (int i=0; i < len; i++) str[i + is_negative] = rev[len-i-1];

        if (is_negative) {
            str[0] = '-';
            len++;
        }
    }

    str[len] = '\0'; // string terminator
    return len;
}

StrText str_int_s (int64_t n)
{
    StrText s;
    str_int (n, s.s);
    return s;
}

// similar to strtoull, except it rejects numbers that are shorter than str_len, or that their reconstruction would be different
// the the original string.
// returns true if successfully parsed an integer of the full length of the string
bool str_get_int (const char *str, unsigned str_len, 
                  int64_t *value) // out - modified only if str is an integer
{
    int64_t out = 0;

     // edge cases - "" ; "-" ; "030" ; "-0" ; "-030" - false, because if we store these as an integer, we can't reconstruct them
    if (!str_len || 
        (str_len == 1 && str[0] == '-') || 
        (str_len >= 2 && str[0] == '0') || 
        (str_len >= 2 && str[0] == '-' && str[1] == '0')) return false;

    unsigned negative = (str[0] == '-');

    for (unsigned i=negative; i < str_len; i++) {
        if (!IS_DIGIT(str[i])) return false;

        int64_t prev_out = out;
        out = (out * 10) + (str[i] - '0');

        if (out < prev_out) return false; // number overflowed beyond maximum int64_t
    }

    if (negative) out *= -1;

    if (value) *value = out; // update only if successful
    return true;
}

#define str_get_int_range_type(func_num,type) \
bool str_get_int_range##func_num (const char *str, unsigned str_len, type min_val, type max_val, type *value /* optional */) \
{                                                                                     \
    if (!str) return false;                                                           \
                                                                                      \
    int64_t value64;                                                                  \
    if (!str_get_int (str, str_len ? str_len : strlen (str), &value64)) return false; \
    if (value) *value = (type)value64;                                                \
                                                                                      \
    return (type)value64 >= min_val && (type)value64 <= max_val;                      \
}
str_get_int_range_type(8,uint8_t)   // unsigned
str_get_int_range_type(16,uint16_t) // unsigned
str_get_int_range_type(32,int32_t)  // signed
str_get_int_range_type(64,int64_t)  // signed

// get a positive hexadecimal integer, may have leading zeros eg 00FFF
bool str_get_int_hex (const char *str, unsigned str_len, 
                      uint64_t *value) // out - modified only if str is an integer
{
    uint64_t out = 0;

    for (unsigned i=0; i < str_len; i++) {
        uint64_t prev_out = out;
        char c = str[i];
        
        if      (c >= '0' && c <= '9') out = (out * 16) + (c - '0');
        else if (c >= 'A' && c <= 'F') out = (out * 16) + (c - 'A' + 10);
        else if (c >= 'a' && c <= 'f') out = (out * 16) + (c - 'a' + 10);
        else return false;

        if (out < prev_out) return false; // number overflowed beyond maximum uint64_t
    }

    if (value) *value = out; // update only if successful
    return true;
}

// positive integer, may be hex prefixed with 0x
#define str_get_int_range_allow_hex_bits(bits) \
bool str_get_int_range_allow_hex##bits (const char *str, unsigned str_len, uint##bits##_t min_val, uint##bits##_t max_val, uint##bits##_t *value) \
{                                                                                             \
    if (!str) return false;                                                                   \
                                                                                              \
    if (str_len >= 3 && str[0]=='0' && str[1]=='x') {                                         \
        uint64_t value64; /* unsigned */                                                      \
        if (!str_get_int_hex (&str[2], str_len ? str_len-2 : strlen (str)-2, &value64)) return false; \
        *value = (uint##bits##_t)value64;                                                     \
    }                                                                                         \
    else {                                                                                    \
        int64_t value64; /* signed bc str_get_int 0 - max is 07fff... */                      \
        if (!str_get_int (str, str_len ? str_len : strlen (str), &value64)) return false;     \
        *value = (uint##bits##_t)value64;                                                     \
    }                                                                                         \
                                                                                              \
    return *value >= min_val && *value <= max_val;                                            \
}
str_get_int_range_allow_hex_bits(8)  // unsigned
str_get_int_range_allow_hex_bits(16) // unsigned
str_get_int_range_allow_hex_bits(32) // unsigned
str_get_int_range_allow_hex_bits(64) // unsigned



StrText str_uint_commas (int64_t n)
{
    StrText s;

    unsigned len = 0, orig_len=0;

    if (n==0) {
        s.s[0] = '0';
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
        for (int i=0; i < len; i++) s.s[i] = rev[len-i-1];
    }

    s.s[len] = '\0'; // string terminator
    return s;
}

// "3.123" -> "%5.3f" ;  returns length of format string
unsigned str_get_float_format (const char *float_str, unsigned float_str_len, char *str /* out */)
{
    unsigned decimal_digits=0;
    for (int i=(int)float_str_len-1; i >= 0; i--)
        if (float_str[i] == '.') {
            decimal_digits = (float_str_len-1) - i;
            break;
        }

    unsigned next=0;
    str[next++] = '%';
    next += str_int (float_str_len, &str[next]);    
    str[next++] = '.';
    next += str_int (decimal_digits, &str[next]);
    str[next++] = 'f';
    str[next] = 0;

    return next;
}

// returns float value, or -1 if not a positive float of at most 9 digits after the decimal point
double str_get_positive_float (const char *float_str, unsigned float_str_len)
{
    bool in_decimals=false;
    unsigned num_decimals=0;
    double val = 0;

    for (unsigned i=0; i < float_str_len; i++) {
        if (float_str[i] == '.' && !in_decimals)
            in_decimals = true;
    
        else if (IS_DIGIT (float_str[i])) {
            val = (val * 10) + (float_str[i] - '0');
            if (in_decimals) num_decimals++;
        }
    
        else return -1; // not a positive float
    }

    static const double pow10[10] = { 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000 };
    if (num_decimals >= sizeof (pow10) / sizeof (pow10[0])) return -1; // too many decimals

    return val / pow10[num_decimals];
}

StrText str_pointer (const void *p)
{
    StrText s;
    sprintf (s.s, "0x%"PRIx64, (uint64_t)p);
    return s;
}

bool str_is_in_range (const char *str, uint32_t str_len, char first_c, char last_c)
{
    for (; str_len ; str++, str_len--)
        if (*str < first_c || *str > last_c) return false;
    return true;
}

// returns true if the strings are case-insensitive similar, and *identical if they are case-sensitive identical
bool str_case_compare (const char *str1, const char *str2, unsigned len, 
                       bool *identical) // optional out
{
    bool my_identical = true; // optimistic (automatic var)

    for (unsigned i=0; i < len; i++)
        if (str1[i] != str2[i]) {
            my_identical = false; // NOT case-senstive identical
            if (UPPER_CASE(str1[i]) != UPPER_CASE(str2[i])) {
                if (identical) *identical = false;
                return false;
            }
        }

    if (identical) *identical = my_identical;
    return true; // case-insenstive identical
}

// splits a string with (num_items-1) separators (doesn't need to be nul-terminated) to up to or exactly num_items
// returns the actual number of items, or 0 is unsuccessful
unsigned str_split (const char *str, unsigned str_len, uint32_t num_items, char sep,
                    const char **items,  // out - array of char* of length num_items - one more than the number of separators
                    unsigned *item_lens, // optional out - corresponding lengths
                    bool exactly,
                    const char *enforce_msg)   // non-NULL if enforcement of length is requested

{
    items[0] = str;
    uint32_t item_i=1;

    for (uint32_t i=0; i < str_len ; i++) 
        if (str[i] == sep) {
            if (item_i == num_items) {
                ASSERT (!enforce_msg, "expecting up to %u %s separators but found more: (100 first) %.*s", 
                        num_items-1, enforce_msg, MIN (str_len, 100), str);
                return 0; // too many separators
            }
            items[item_i++] = &str[i+1];
        }

    if (item_lens) {
        for (uint32_t i=0; i < item_i-1; i++)    
            item_lens[i] = items[i+1] - items[i] - 1;
            
        item_lens[item_i-1] = &str[str_len] - items[item_i-1];
    }

    ASSERT (!exactly || !enforce_msg || item_i == num_items, "Expecting the number of %s to be %u, but it is %u: (100 first) \"%.*s\"", 
            enforce_msg, num_items, str_len ? item_i : 0, MIN (100, str_len), str);
    
    return (!exactly || item_i == num_items) ? item_i : 0; // 0 if requested exactly, but too few separators 
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

void str_print_null_seperated_data (const char *data, unsigned len, bool add_newline, bool remove_equal_asterisk)
{
    for (unsigned i=0; i < len; i++) {
        // in case we are showing chrom data in --list-chroms in SAM - don't show * and =
        if (remove_equal_asterisk && (data[i]=='*' || data[i]=='=') && !data[i+1]) {
            i++;
            continue; // skip character and following separator
        }

        switch (data[i]) {
            case 32 ... 127 : fputc (data[i], info_stream);      break;
            case 0          : fputc (add_newline ? '\n' : ' ', info_stream); break; // snip separator
            case '\t'       : fwrite ("\\t", 1, 2, info_stream); break;
            case '\n'       : fwrite ("\\n", 1, 2, info_stream); break;
            case '\r'       : fwrite ("\\r", 1, 2, info_stream); break;
            default         : iprintf ("\\x%x", data[i]);
        }
        }
}

int str_print_text (const char **text, unsigned num_lines,
                    const char *wrapped_line_prefix, 
                    const char *newline_separator, 
                    unsigned line_width /* 0=calcuate optimal */)
{                       
    ASSERTNOTNULL (text);
    
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
            iprintf ("%s%.*s\n", wrapped ? wrapped_line_prefix : "", c, line);
            line += c + (line[c]==' '); // skip space too
            line_len -= c + (line[c]==' ');
            wrapped = true;
        }
        iprintf ("%s%s%s", wrapped ? wrapped_line_prefix : "", line, newline_separator);
    }
    return 0;
}

// receives a user response, a default "Y" or "N" (or NULL) and modifies the response to be "Y" or "N"
bool str_verify_y_n (char *response, unsigned response_size, const char *y_or_n)
{
    ASSERT0 (!y_or_n || (strlen (y_or_n)==1 && (y_or_n[0]=='Y' || y_or_n[0]=='N')), 
              "y_or_n needs to be NULL, \"Y\" or \"N\"");

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

const char *str_win_error (void)
{
    static char msg[100];
#ifdef _WIN32
    FormatMessageA (FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,                   
                    NULL, GetLastError(), MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), msg, sizeof (msg), NULL);
#endif
    return msg;
}

// current date and time
StrText str_time (void)
{
    StrText s;
    time_t now = time (NULL);
    strftime (s.s, 100, "%Y-%m-%d %H:%M:%S ", localtime (&now));
    strcpy (&s.s[strlen(s.s)], tzname[daylight]);
    return s;
}

