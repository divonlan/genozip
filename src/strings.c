// ------------------------------------------------------------------
//   strings.c
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <time.h>
#include "genozip.h"
#include "strings.h"
#include "flags.h"
#include "context.h"
#ifndef WIN32
#include <sys/ioctl.h>
#else
#include <windows.h>
#endif

bool is_printable[256] = { [9]=1, [10]=1, [13]=1, [32 ... 126]=1 };

uint64_t p10[] = { 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000ULL, 100000000000ULL, 1000000000000ULL, 10000000000000ULL, 100000000000000ULL, 1000000000000000ULL, 100000000000000000ULL, 100000000000000000ULL };

char *str_tolower (rom in, char *out /* out allocated by caller - can be the same as in */)
{
    char *startout = out;

    for (; *in; in++, out++) 
        *out = LOWER_CASE (*in);
    
    *out = 0;

    return startout;
}

char *str_toupper (rom in, char *out /* out allocated by caller - can be the same as in */)
{
    char *startout = out;

    for (; *in; in++, out++) 
        *out = UPPER_CASE (*in);

    *out = 0;

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

char *str_print_snip (STRp(in), char *out) // caller allocated - in_len+20
{
    char *save_out = out;

    if (!in) strcpy (out, "NULL");
    
    else if (!in_len) strcpy (out, "\"\"");
    
    else {
        if (in[0] < NUM_SNIP_CODES) {
            static rom snip_codes[NUM_SNIP_CODES] = SNIP_CODES;
            uint32_t len = strlen (snip_codes[(int)in[0]]);
            strcpy (out, snip_codes[(int)in[0]]);
            out[len] = ' ';
            out += len+1;
            in++;
            in_len--;
        }

        *(out++) = '\"';

        for (uint32_t i=0; i < in_len; i++) 
            switch ((uint8_t)in[i]) {
                case 32 ... 126 :  *(out++) = in[i]; break;   // printable ASCII
                case '\t'       : 
                case '\n'       : 
                case '\r'       :  *(out++) = ' '; break; // whitespace
                default         :  *(out++) = '?'; 
            }

        *(out++) = '\"';
        *(out++) = 0; // nul-terminate
    }

    return save_out;
}

// replaces \t, \n, \r, \b with "\t" etc, replaces unprintables with '?'. caller should allocate out. returns out.
// out should be allocated by caller to (in_len*2 + 1), out is null-terminated
char *str_to_printable (STRp(in), char *out)
{
    char *start = out;

    for (uint32_t i=0; i < in_len; i++)
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
    else if (num_bases >  0   )             sprintf (s.s, "%u b",    (uint32_t)num_bases);
    else                                    sprintf (s.s, "-");

    return s;
}

// returns length
uint32_t str_int_ex (int64_t n, char *str /* out */, bool add_nul_terminator)
{
    uint32_t len=0;

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

    if (add_nul_terminator) str[len] = '\0'; // string terminator
    return len;
}

StrText str_int_s (int64_t n)
{
    StrText s;
    str_int (n, s.s);
    return s;
}

// returns length
uint32_t str_hex_ex (int64_t n, char *str /* out */, bool uppercase, bool add_nul_terminator)
{
    uint32_t len=0;

    if (n==0) {
        str[0] = '0';
        len=1;
    }

    else {
        bool is_negative = (n<0);
        if (is_negative) n = -n;

        char rev[50] = {}; // "initialize" to avoid compiler warning
        while (n) {
            int hexit = n % 16;
            rev[len++] = (hexit < 10) ? ('0' + hexit)
                       : uppercase    ? ('A' + (hexit-10))
                       :                ('a' + (hexit-10));
            n /= 16;
        }
        // now reverse it
        for (int i=0; i < len; i++) str[i + is_negative] = rev[len-i-1];

        if (is_negative) {
            str[0] = '-';
            len++;
        }
    }

    if (add_nul_terminator) str[len] = '\0'; // string terminator
    return len;
}
// similar to strtoull, except it rejects numbers that are shorter than str_len, or that their reconstruction would be different
// the the original string.
// returns true if successfully parsed an integer of the full length of the string
bool str_get_int (STRp(str), 
                  int64_t *value) // out - modified only if str is an integer
{
    int64_t out = 0;

     // edge cases - "" ; "-" ; "030" ; "-0" ; "-030" - false, because if we store these as an integer, we can't reconstruct them
    if (!str_len || 
        (str_len == 1 && str[0] == '-') || 
        (str_len >= 2 && str[0] == '0') || 
        (str_len >= 2 && str[0] == '-' && str[1] == '0')) return false;

    uint32_t negative = (str[0] == '-');

    for (uint32_t i=negative; i < str_len; i++) {
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
bool str_get_int_range##func_num (rom str, uint32_t str_len, int64_t min_val, int64_t max_val, type *value /* optional */) \
{                                                                                     \
    if (!str) return false;                                                           \
                                                                                      \
    int64_t value64;                                                                  \
    if (!str_get_int (str, str_len ? str_len : strlen (str), &value64)) return false; \
    if (value) *value = (type)value64;                                                \
                                                                                      \
    return value64 >= min_val && value64 <= max_val;                                  \
}
str_get_int_range_type(8,uint8_t)   // unsigned
str_get_int_range_type(16,uint16_t) // unsigned
str_get_int_range_type(32,int32_t)  // signed
str_get_int_range_type(64,int64_t)  // signed

bool str_get_uint32 (STRp(str), uint32_t *value)
{
    int64_t v64;
    if (!str_get_int_range64 (STRa(str), 0, 0xffffffff, &v64)) return false;

    *value = v64;
    return true;
}

// get a positive decimal integer, may have leading zeros eg 005
bool str_get_int_dec (STRp(str), 
                      uint64_t *value) // out - modified only if str is an integer
{
    int64_t out = 0;

    for (uint32_t i=0; i < str_len; i++) {
        if (!IS_DIGIT(str[i])) return false;

        uint64_t prev_out = out;
        out = (out * 10) + (str[i] - '0');

        if (out < prev_out) return false; // number overflowed beyond maximum uint64_t
    }

    if (value) *value = out; // update only if successful
    return true;
}

// get a positive hexadecimal integer, may have leading zeros eg 00FFF
bool str_get_int_hex (STRp(str), bool allow_hex, bool allow_HEX,
                      uint64_t *value) // out - modified only if str is an integer
{
    uint64_t out = 0;

    for (uint32_t i=0; i < str_len; i++) {
        uint64_t prev_out = out;
        char c = str[i];
        
        if      (c >= '0' && c <= '9') out = (out << 4) | (c - '0');
        else if (allow_HEX && c >= 'A' && c <= 'F') out = (out << 4) | (c - 'A' + 10);
        else if (allow_hex && c >= 'a' && c <= 'f') out = (out << 4) | (c - 'a' + 10);
        else return false;

        if (out < prev_out) return false; // number overflowed beyond maximum uint64_t
    }

    if (value) *value = out; // update only if successful
    return true;
}

// positive integer, may be hex prefixed with 0x
#define str_get_int_range_allow_hex_bits(bits) \
bool str_get_int_range_allow_hex##bits (rom str, uint32_t str_len, uint##bits##_t min_val, uint##bits##_t max_val, uint##bits##_t *value) \
{                                                                                             \
    if (!str) return false;                                                                   \
                                                                                              \
    if (str_len >= 3 && str[0]=='0' && str[1]=='x') {                                         \
        uint64_t value64; /* unsigned */                                                      \
        if (!str_get_int_hex (&str[2], str_len ? str_len-2 : strlen (str)-2, true, true, &value64)) return false; \
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

// caller should allocate hex_str[data_len*2+1] (or *3 if with_dot). returns nul-terminated string.
rom str_to_hex (bytes data, uint32_t data_len, char *hex_str, bool with_dot)
{
    char *s = hex_str;

    for (int i=0; i < data_len; i++) {
        *s++ = NUM2HEXDIGIT(data[i] >> 4);
        *s++ = NUM2HEXDIGIT(data[i] & 0xf);
        if (with_dot) *s++ = '.';
    }

    if (with_dot) s--; // remove terminating dot

    *s = 0;
    return hex_str;
}

StrText str_hex10 (bytes data, uint32_t data_len)
{
    StrText s;
    str_to_hex (data, MIN_(10, data_len), s.s, false);
    return s;
}


StrText str_int_commas (int64_t n)
{
    StrText s;
    char *c = s.s; 

    uint32_t len = 0, orig_len=0;

    if (n==0) {
        *c++ = '0';
        len = 1;
    }
    else {
        if (n < 0) {
            *c++ = '-';
            n = -n;
        }
    
        char rev[50] = {}; 
        while (n) {
            if (orig_len && orig_len % 3 == 0) rev[len++] = ',';    
            rev[len++] = '0' + n % 10;
            orig_len++;
            n /= 10;
        }
        // now reverse it
        for (int i=0; i < len; i++) *c++ = rev[len-i-1];
    }

    *c = '\0'; // string terminator
    return s;
}

// display comma'd number up and including the limit, and thereafter with K, M etc.
StrText str_uint_commas_limit (uint64_t n, uint64_t limit)
{
    if (n <= limit) return str_int_commas (n);

    StrText s;

    if      (n >= (1LL << 50)) sprintf (s.s, "%3.1lfP", ((double)n) / 1000000000000000.0);
    else if (n >= (1LL << 40)) sprintf (s.s, "%3.1lfT", ((double)n) / 1000000000000.0);
    else if (n >= (1LL << 30)) sprintf (s.s, "%3.1lfG", ((double)n) / 1000000000.0);
    else if (n >= (1LL << 20)) sprintf (s.s, "%3.1lfM", ((double)n) / 1000000.0);
    else if (n >= (1LL << 10)) sprintf (s.s, "%3.1lfK", ((double)n) / 1000.0);
    else if (n >  0          ) sprintf (s.s, "%3d",     (int)n);
    else                       sprintf (s.s, "-");

    return s;
}

uint32_t str_get_uint_textual_len (uint64_t n) 
{ 
    for (int i=0; i < ARRAY_LEN(p10); i++)
        if (n < p10[i]) return i+1;

    ABORT_R ("n=%"PRIu64" too big", n);
}

// returns 32 bit float value and/or format: "3.123" -> "%5.3f" ; false if not a simple float
bool str_get_float (STRp(float_str), 
                    double *value, char format[FLOAT_FORMAT_LEN], uint32_t *format_len) // optional outs (format allocated by caller)
{
    // TODO: add support for %e and %E (matinsa/exponent) formats

    bool in_decimals=false;
    uint32_t num_decimals=0;
    double val = 0;
    bool is_negative = (float_str[0] == '-');

    for (uint32_t i=is_negative; i < float_str_len; i++) {
        if (float_str[i] == '.' && !in_decimals)
            in_decimals = true;
    
        else if (IS_DIGIT (float_str[i])) {
            val = (val * 10) + (float_str[i] - '0');
            if (in_decimals) num_decimals++;
        }
    
        else return false; // can't interpret this string as float
    }

    // calculate format if requested - the format string is in a format expected by reconstruct_from_local_float
    if (format) {

        if (float_str_len > 99) return false; // we support format of float strings up to 99 characters... more than enough
        uint32_t next=0;
        format[next++] = '%';
        if (float_str_len >= 10) format[next++] = '0' + (float_str_len / 10);
        format[next++] = '0' + (float_str_len % 10);
        format[next++] = '.';
        if (num_decimals >= 10) format[next++] = '0' + (num_decimals / 10);
        format[next++] = '0' + (num_decimals % 10);
        format[next++] = 'f';
        format[next] = 0;

        if (format_len) *format_len = next;
    }

    if (value) {
        static const double pow10[16] = { 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0, 10000000.0, 100000000.0, 1000000000.0, 
                                          10000000000.0, 100000000000.0, 1000000000000.0, 10000000000000.0, 100000000000000.0, 1000000000000000.0 };
                                        
        if (num_decimals >= ARRAY_LEN (pow10)) return false; // too many decimals

        *value = (is_negative ? -1 : 1) * (val / pow10[num_decimals]);
    }

    return true;
}

// if value is a float in scientific notation eg 4.31e-03, it is converted to eg 0.00431 and true is returned. otherwise false is returned.
bool str_scientific_to_decimal (STRp(float_str), char *modified, uint32_t *modified_len /* in / out */, double *value)
{
    // short circuit normal floats eg 0.941
    if (float_str_len < 5) return false; // scientific notation has a minimum of 5 characters eg 3e-05 
    
    bool negative = (float_str[0] == '-');
    bool has_decimal = float_str[negative+1] == '.';
    if (has_decimal && float_str_len < 7 + negative) return false; // short circuit common normal numbers eg 0.885
    if (!has_decimal && float_str[negative+1] != 'e' && float_str[negative+1] != 'E') return false; // expecting [-]3. or [-]3e or [-]3E (single digit mantissa)

    SAFE_NUL (&float_str[float_str_len]);  // no "return" until SAFE_RESTORE

        int mantissa_len = strcspn (float_str, "eE");
        if (mantissa_len == float_str_len) goto not_scientific_float; // no e or E
        
        int exp_len = float_str_len - mantissa_len - 1;
        if (exp_len < 3) goto not_scientific_float; // not standard notiation, expecting eg e-02

        char *after;
        int exp = strtod (&float_str[mantissa_len+1], &after); 
        if (after != float_str + float_str_len) goto not_scientific_float;

        double f = atof (float_str);

    SAFE_RESTORE;

    if (exp > -1) return false; // this function currently works only for 0.*** numbers (all digits after the decimal point). TO DO: remove limitation

    int decimal_digits = mantissa_len - negative - has_decimal + (-exp) - 1; // eg. -2.30e-02 --> -0.0230  mantissa_len=5 exp=-2 --> width=7
    
    if (modified) {
        if (*modified_len < decimal_digits + 2 + negative) return false; // not enough room (+1 for \0)
        
        sprintf (modified, "%.*f", decimal_digits, f);
        *modified_len = decimal_digits + 2 + negative;
    }

    if (value) *value = f;

    return true;

not_scientific_float:
    SAFE_RESTORE;
    return false;
}

bool str_is_in_range (rom str, uint32_t str_len, char first_c, char last_c)
{
    for (; str_len ; str++, str_len--)
        if (*str < first_c || *str > last_c) return false;
    return true;
}

// returns true if the strings are case-insensitive similar, and *identical if they are case-sensitive identical
bool str_case_compare (rom str1, rom str2,
                       bool *identical) // optional out
{
    bool my_identical = true; // optimistic (automatic var)
    uint32_t len =strlen (str1);

    if (len != strlen (str2)) goto differ;

    for (uint32_t i=0; i < len; i++)
        if (str1[i] != str2[i]) {
            my_identical = false; // NOT case-senstive identical
            if (UPPER_CASE(str1[i]) != UPPER_CASE(str2[i])) 
                goto differ;
        }

    if (identical) *identical = my_identical;
    return true; // case-insenstive identical

differ:
    if (identical) *identical = false;
    return false;
}

#define ASSSPLIT(condition, format, ...) ({\
    if (!(condition)) { \
        if (enforce_msg) {  \
            progress_newline(); fprintf (stderr, "Error in %s:%u: ", __FUNCLINE); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, SUPPORT); exit_on_error(true); /* same as ASSERT */ \
        } \
        return 0; /* just return 0 if we're not asked to enforce */ \
    } \
})

// splits a string with up to (max_items-1) separators (doesn't need to be nul-terminated) to up to or exactly max_items
// returns the actual number of items, or 0 is unsuccessful
uint32_t str_split_do (STRp(str), 
                       uint32_t max_items,    // optional - if not given, a count of sep is done first
                       char sep,              // separator
                       rom *items,            // out - array of char* of length max_items - one more than the number of separators
                       uint32_t *item_lens,   // optional out - corresponding lengths
                       bool exactly,
                       rom enforce_msg)       // non-NULL if enforcement of length is requested
{
    if (!str) return 0; // note: str!=NULL + str_len==0 results in n_items=1 - one empty item
    
    items[0] = str;
    uint32_t item_i = 1;
    for (uint32_t str_i=0 ; str_i < str_len ; str_i++) 
        if (str[str_i] == sep) {
            ASSSPLIT (item_i < max_items, "expecting up to %u %s separators but found more: (100 first) %.*s", 
                      max_items-1, enforce_msg, MIN_(str_len, 100), str);

            items[item_i++] = &str[str_i+1];
        }

    if (item_lens) {
        for (uint32_t i=0; i < item_i-1; i++)    
            item_lens[i] = items[i+1] - items[i] - 1; 
            
        item_lens[item_i-1] = &str[str_len] - items[item_i-1];
    }

    ASSERT (!exactly || !enforce_msg || item_i == max_items, "Expecting the number of %s to be %u, but it is %u: (100 first) \"%.*s\"", 
            enforce_msg, max_items, str_len ? item_i : 0, MIN_(100, str_len), str);
    
    return (!exactly || item_i == max_items) ? item_i : 0; // 0 if requested exactly, but too few separators 
}

// splits a string by tab, ending with the first \n or \r\n. 
// Returns the address of the byte after the \n if successful, or NULL if not.
rom str_split_by_tab_do (STRp(str), 
                         uint32_t *n_items, // in / out
                         rom *items, uint32_t *item_lens, // out - array of char* of length max_items - one more than the number of separators
                         bool *has_13,      // optional out - true if line is terminated by \r\n instead of \n
                         bool enforce_msg)  // failure handling: if false return NULL ; if true ABORT

{
    ASSERTNOTNULL (str);
    ASSERTNOTZERO (*n_items, "");

    items[0] = str;
    uint32_t item_i = 1;
    uint32_t str_i;
    bool my_has_13;

    for (str_i=0 ; str_i < str_len ; str_i++) {
        char c = str[str_i]; 
        if (c == '\t') {
            ASSSPLIT (item_i < *n_items, "expecting up to %u items but found more", *n_items);

            items[item_i++] = &str[str_i+1];
        }

        else if (c == '\r') {
            ASSSPLIT (str_i+1 < str_len && str[str_i+1] == '\n', "encountered a \\r without a following \\n (str_i=%u)", str_i);
            my_has_13 = true;
            break;
        }

        else if (c == '\n') {
            my_has_13 = false;
            break;
        }
    }

    ASSSPLIT (str_i < str_len, "Line not terminated by newline (str_len=%u)", str_len);

    for (uint32_t i=0; i < item_i-1; i++)    
        item_lens[i] = items[i+1] - items[i] - 1; 
        
    item_lens[item_i-1] = &str[str_i] - items[item_i-1];
    *n_items = item_i;

    if (has_13) *has_13 = my_has_13;

    return &str[str_i + 1 + my_has_13]; // byte after \n
}

// splits a string based on container items (doesn't need to be nul-terminated). 
// returns the number on unskipped items if successful
uint32_t str_split_by_container_do (STRp(str), ConstContainerP con, STRp(con_prefixes),
                                    rom *items,                // out - array of char* of length max_items - one more than the number of separators
                                    uint32_t *item_lens,       // optional out - corresponding lengths
                                    rom enforce_msg)           // non-NULL if enforcement of length is requested
{
    if (!str) return 0; // note: str!=NULL + str_len==0 results in n_items=1 - one empty item

    uint32_t num_items = con_nitems (*con);
    ASSERT0 (num_items, "Container has no items");

    rom *save_items = items;
    rom after_str = &str[str_len], save_str = str;
    rom px = 0, after_px=0;

    // if we have prefixes, set px_i to the first item's prefix, after the container-wide prefix
    if (con_prefixes_len) {
        px = &con_prefixes[1];
        after_px = &con_prefixes[con_prefixes_len];
        while (*px != CON_PX_SEP) px++; // skip over container-wide prefix
        px++; // skip CON_PX_SEP to first item prefix
    }

    for (uint32_t item_i=0; item_i < num_items; item_i++) { 

        char sep = con->items[item_i].separator[0];
        
        // verify item prefix
        if (px < after_px) {
            rom item_px = px;
            while (*px != CON_PX_SEP) px++;

            uint32_t item_px_len = px - item_px;                   
            px++; // skip CON_PX_SEP

            if (item_px_len) {
                ASSERT (sep != CI0_INVISIBLE, "item_i=%u is CI0_INVISIBLE, expecting item_px_len=%u to be 0", item_i, item_px_len);

                if (sep != CI0_SKIP) { // if CI0_SKIP, we ignore the prefix
                    bool prefix_matches = (str + item_px_len <= after_str) && !memcmp (str, item_px, item_px_len);
                    ASSSPLIT (prefix_matches, "prefix mismatch for item_i=%u: expecting \"%.*s\" but seeing \"%.*s\"",
                              item_i-1, item_px_len, item_px, MIN_(item_px_len, (int)(after_str-str)), str);
    
                    str += item_px_len; // advance past prefix
                }
            }
        }

        switch (sep) {

            case CI0_SKIP: 
            case CI0_INVISIBLE:
                *item_lens = 0;
                *items = str; // zero-length item - next item will start from the same str_i
                break;

            case CI0_FIXED_0_PAD: 
                *item_lens = con->items[item_i].separator[1];
                *items = str; // zero-length item - next item will start from the same str_i
                str += *item_lens;
                ASSSPLIT (str <= after_str, "item_i=%u fixed_len=%u goes beyond end of string \"%.*s\"", item_i, *item_lens, str_len, save_str);
                break;
            
            case 0: // no separator - goes to end of string
                *item_lens = after_str - str;
                *items = str; // zero-length item - next item will start from the same str_i
                str = after_str;
                break;

            case CI0_DIGIT: // items goes until a digit or end-of-string is encountered
                *items = str;

                while (str < after_str && !IS_DIGIT(*str)) str++;  

                *item_lens = str - *items;
                break;

            default:
                ASSERT (IS_PRINTABLE(sep), "item_i=%u sep=%u is not a printable character in string \"%.*s\"", item_i, sep, str_len, save_str);

                char sep1 = con->items[item_i].separator[1];
                if (!IS_PRINTABLE (sep1)) sep1 = 0;
                
                *items = str;

                if (!sep1) 
                    while (str < after_str && *str != sep) str++; 
                else 
                    while (str < (after_str-1) && (str[0] != sep || str[1] != sep1)) str++; 

                ASSSPLIT (str < after_str, "item_i=%u reached end of string without finding separator '%c' in string \"%.*s\"", 
                          item_i, sep, str_len, save_str);

                *item_lens = str - *items;
                str += 1 + (sep1 != 0); // skip seperator
        }

        items++; item_lens++; // increment pointers        
    }

    ASSSPLIT (str == after_str, "container consumed only %u of %u characters of string \"%.*s\"", 
              (int)(str-save_str), str_len, str_len, save_str);

    return items - save_items;
}

// remove \r (ASCII Carriage Return) from each lines[] that has it as its final character, but decrementing the matching line_lens
void str_remove_CR_do (uint32_t n_lines, rom *lines, uint32_t *line_lens)
{
    for (uint32_t i=0; i < n_lines; i++)
        if (line_lens[i] >= 1 && lines[i][line_lens[i]-1] == '\r') 
            line_lens[i]--;
}

// replace the last character in each item with \0. these are generated by str_split, so expecting a separator after each string
void str_nul_separate_do (STRps(item))
{
    for (uint32_t i=0; i < n_items; i++)
        ((char**)items)[i][item_lens[i]] = '\0';
}

// generate new string (mem allocated by caller) - copy of in, but removing all ASCII < 33 or > 126 
uint32_t str_remove_whitespace (STRp(in), char *out)
{
    uint32_t out_len = 0;
    for (uint32_t i=0; i < in_len; i++)
        if (in[i] >= 33 && in[i] <= 126) 
            out[out_len++] = in[i];

    return out_len;
}

// in-place removal of flanking whitespace from a null-terminated string
void str_trim (STRe(str))
{
    // remove leading whitespave
    int i=0; for (; i < *str_len; i++)
        if (str[i] != ' ' && str[i] != '\t' && str[i] != '\n' && str[i] != '\r')
            break;

    if (i) {
        *str_len -= i;
        memmove (str, str+i, *str_len + 1); // +1 to move \0 as well
    }

    // remove trailing whitespace
    for (i = *str_len - 1; i >= 0; i--)
        if (str[i] != ' ' && str[i] != '\t' && str[i] != '\n' && str[i] != '\r') 
            break;

    if (i < *str_len - 1) {
        *str_len = i+1;
        str[*str_len] = '\0';
    }
}

// splits a string with up to (max_items-1) separators (doesn't need to be nul-terminated) to up to or exactly max_items integers
// returns the actual number of items, or 0 is unsuccessful
uint32_t str_split_ints_do (STRp(str), uint32_t max_items, char sep, bool exactly,
                            int64_t *items)  // out - array of integers
                       
{
    rom after = &str[str_len];
    SAFE_NUL (after);

    uint32_t item_i;
    for (item_i=0; item_i < max_items && str < after; item_i++, str++) {
        items[item_i] = strtoll (str, (char **)&str, 10);
        if (item_i < max_items-1 && *str != sep) {
            item_i=0;
            break; // fail
        }
    }

    if (str < after) item_i = 0;

    SAFE_RESTORE;

    return (!exactly || item_i == max_items) ? item_i : 0; // 0 if requested exactly, but too few separators 
}

// splits a string with up to (max_items-1) separators (doesn't need to be nul-terminated) to up to or exactly max_items floats
// returns the actual number of items, or 0 is unsuccessful
uint32_t str_split_floats_do (STRp(str), uint32_t max_items, char sep, bool exactly,
                              double *items)  // out - array of floats
                       
{
    rom after = &str[str_len];
    SAFE_NUL (after);

    uint32_t item_i;
    for (item_i=0; item_i < max_items && str < after; item_i++, str++) {
        items[item_i] = strtod (str, (char**)&str);
        if (item_i < max_items-1 && *str != sep) {
            item_i=0;
            break; // fail
        }
    }

    if (str < after) item_i = 0;

    SAFE_RESTORE;

    return (!exactly || item_i == max_items) ? item_i : 0; // 0 if requested exactly, but too few separators 
}

rom type_name (uint32_t item, 
               rom const *name, // the address in which a pointer to name is found, if item is in range
               uint32_t num_names)
{
    if (item > num_names) {
        static char str[50];
        sprintf (str, "%d (out of range)", item);
        return str;
    }
    
    return *name;    
}

void str_print_dict (FILE *fp, STRp(data), bool add_newline, bool remove_equal_asterisk)
{
    WordIndex word_index=0;

    for (uint32_t i=0; i < data_len; i++) {
        // in case we are showing chrom data in --list-chroms in SAM - don't show * and =
        if (remove_equal_asterisk && (data[i]=='*' || data[i]=='=') && !data[i+1]) {
            i++;
            continue; // skip character and following separator
        }

        if (!i || data[i] || data[i-1]) { // skip empty words
            
            if (!i || !data[i-1]) fprintf (fp, "%u%c", word_index, add_newline ? '\t' : '='); 

            switch (data[i]) {
                case 32 ... 127 : fputc (data[i], fp);      break;
                case 0          : fputc (add_newline ? '\n' : ' ', fp); break; // snip separator
                case '\t'       : fwrite ("\\t", 1, 2, fp); break;
                case '\n'       : fwrite ("\\n", 1, 2, fp); break;
                case '\r'       : fwrite ("\\r", 1, 2, fp); break;
                default         : fprintf (fp, "\\x%x", (uint8_t)data[i]);
            }
        }

        if (!data[i]) word_index++;
    }
    fflush (fp);
}

int str_print_text (rom *text, uint32_t num_lines,
                    rom wrapped_line_prefix, rom newline_separator, 
                    rom added_header, // optional
                    uint32_t line_width /* 0=calcuate optimal */)
{                       
    ASSERTNOTNULL (text);
    
    if (!line_width) {
#ifdef _WIN32
        line_width = 120; // default width of cmd window in Windows 10
#else
        // in Linux and Mac, we can get the actual terminal width
        struct winsize w;
        ioctl(0, TIOCGWINSZ, &w);
        line_width = MAX_(40, w.ws_col); // our wrapper cannot work with to small line widths 
#endif
    }

    if (added_header) iprintf ("%s\n", added_header);
    
    for (uint32_t i=0; i < num_lines; i++)  {
        rom line = text[i];
        uint32_t line_len = strlen (line);
        
        // print line with line wraps
        bool wrapped = false;
        while (line_len + (wrapped ? strlen (wrapped_line_prefix) : 0) > line_width) {
            int c; for (c=line_width-1 - (wrapped ? strlen (wrapped_line_prefix) : 0); 
                        c>=0 && (IS_LETTER(line[c]) || IS_DIGIT (line[c]) || line[c]==','); // wrap lines at - and | too, so we can break very long regex strings like in genocat
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
static bool str_verify_y_n (char *response, uint32_t len, rom def_res)
{
    ASSERT0 (!def_res || (strlen(def_res)==1 && (def_res[0]=='Y' || def_res[0]=='N')), 
              "def_res needs to be NULL, \"Y\" or \"N\"");

    if (len > 1 || (!len && !def_res)) return false;
    
    response[0] = !len             ? def_res[0] 
                : response[0]=='n' ? 'N' 
                : response[0]=='y' ? 'Y' 
                :                    response[0];
    
    response[1] = 0;

    return (response[0] == 'N' || response[0] == 'Y'); // return false if invalid response - we request the user to respond again
}

bool str_query_user_yn (rom query, DefAnswerType def_answer)
{
    char query_str[strlen(query)+32];
    sprintf (query_str, "%s (%sy%s or %sn%s) ", query, 
             def_answer==QDEF_YES?"[":"",  def_answer==QDEF_YES?"]":"",  
             def_answer==QDEF_NO ?"[":"",  def_answer==QDEF_NO ?"]":"");

    char y_n[32];
    str_query_user (query_str, y_n, sizeof(y_n), def_answer != QDEF_NONE,
                    str_verify_y_n, def_answer==QDEF_YES ? "Y" : def_answer==QDEF_NO ? "N" : 0);

    return y_n[0] == 'Y';
}

void str_query_user (rom query, char *response, uint32_t response_size, bool allow_empty, 
                     ResponseVerifier verifier, rom verifier_param)
{
    uint32_t len;
    do {
        fprintf (stderr, "%s", query);

        // Linux: in case of non-Latin, fgets doesn't handle backspace well - not completely removing the multi-byte character from the string (bug 593)
        // Windows (MingW): the string is terminated before the first non-Latin character on bash terminal and only ???? on Powershell and Command (bug 594)
        do {
            fgets (response, response_size, stdin); 
            len = strlen (response);
            str_trim (response, &len);
        } while (!len && !allow_empty);

    } while (verifier && !verifier (response, len, verifier_param));
}

rom str_win_error (void)
{
    static char msg[100];
#ifdef _WIN32
    FormatMessageA (FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,                   
                    NULL, GetLastError(), MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), msg, sizeof (msg), NULL);
#endif
    return msg;
}

// C<>G A<>T c<>g a<>t ; IUPACs: R<>Y K<>M B<>V D<>H W<>W S<>S N<>N (+ lowercase); other ASCII 32->126 preserved ; other = 0
const char COMPLEM[256] = "-------------------------------- !\"#$\%&'()*+,-./0123456789:;<=>?@TVGHEFCDIJMLKNOPQYSAUBWXRZ[\\]^_`tvghefcdijmlknopqysaubwxrz{|}~";

// same as COMPLEM[UPPER_CASE(c)]
const char UPPER_COMPLEM[256] = "-------------------------------- !\"#$\%&'()*+,-./0123456789:;<=>?@TVGHEFCDIJMLKNOPQYSAUBWXRZ[\\]^_`TVGHEFCDIJMLKNOPQYSAUBWXRZ{|}~";

// reverse-complements a string in-place
char *str_revcomp_in_out (char *dst_seq, rom src_seq, uint32_t seq_len)
{
    for (uint32_t i=0; i < seq_len / 2; i++) {
        char l_base = src_seq[i];
        char r_base = src_seq[seq_len-1-i];

        dst_seq[i]           = COMPLEM[(uint8_t)r_base];
        dst_seq[seq_len-1-i] = COMPLEM[(uint8_t)l_base];
    }

    if (seq_len % 2) // we have an odd number of bases - now complement the middle one
        dst_seq[seq_len/2] = COMPLEM[(uint8_t)src_seq[seq_len/2]];

    return dst_seq;
}

// reverse-complements a string in-place - complements A,C,G,T characters and leaves others intact. set dst_seq=src_seq to output in-place.
char *str_revcomp_actg (char *dst_seq, rom src_seq, uint32_t seq_len)
{
    for (uint32_t i=0; i < seq_len / 2; i++) {
        char l_base = src_seq[i];
        char r_base = src_seq[seq_len-1-i];

        dst_seq[i]           = complem(r_base);
        dst_seq[seq_len-1-i] = complem(l_base);
    }

    if (seq_len % 2) // we have an odd number of bases - now complement the middle one
        dst_seq[seq_len/2] = complem(src_seq[seq_len/2]);

    return dst_seq;
}

// implementing (slightly modified) memrchr, as it doesn't exist on Windows
rom my_memrchr (rom str, char c, uint32_t str_len)
{
    for (int32_t i=(int32_t)str_len-1; i >= 0; i++)
        if (str[i] == c) return &str[i];

    return NULL;
}

// print duration in human-readable form eg 1h2' or "1 hour 2 minutes"
void str_human_time (unsigned secs, bool compact, char *str /* out */)
{
    unsigned hours = secs / 3600;
    unsigned mins  = (secs % 3600) / 60;
             secs  = secs % 60;

    if (compact)
        sprintf (str, "%uh%u'%u\"", hours, mins, secs);
    else if (hours) 
        sprintf (str, "%u %s %u %s", hours, hours==1 ? "hour" : "hours", mins, mins==1 ? "minute" : "minutes");
    else if (mins)
        sprintf (str, "%u %s %u %s", mins, mins==1 ? "minute" : "minutes", secs, secs==1 ? "second" : "seconds");
    else 
        sprintf (str, "%u %s", secs, secs==1 ? "second" : "seconds");
}

// current date and time
StrTime str_time (void)
{
    StrTime s;
    time_t now = time (NULL);
    strftime (s.s, sizeof (s.s), "%Y-%m-%d %H:%M:%S ", localtime (&now));
    int len = strlen(s.s);

    strncpy (&s.s[len], tzname[daylight], sizeof (s.s) - len);
    return s;
}
