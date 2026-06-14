// ------------------------------------------------------------------
//   strings.c
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <time.h>
#include <math.h>
#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif
#include "genozip.h"
#include "strings.h"
#include "context.h"
#include "file.h"
#ifndef _WIN32
#include <sys/ioctl.h>
#else
#include <windows.h>
#endif

alignas(64) const bool is_printable[256] = { ['\t']=1, ['\n']=1, ['\r']=1, [32 ... 126]=1 };

alignas(64) const bool is_ACGT[256] = { ['A']=true, ['C']=true, ['G']=true, ['T']=true };

alignas(64) const bool is_ACGTN[256] = { ['A']=true, ['C']=true, ['G']=true, ['T']=true, ['N']=true };

// valid characters in a FASTQ sequence
alignas(64) const bool is_fastq_seq[256] = { 
    ['A']=true, ['C']=true, ['D']=true, ['G']=true, ['H']=true, ['K']=true, ['M']=true, ['N']=true, 
    ['R']=true, ['S']=true, ['T']=true, ['V']=true, ['W']=true, ['Y']=true, ['U']=true, ['B']=true 
};

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
        case 32 ... '\\'-1  :
        case '\\'+1 ... 126 : return (StrText) { .s[0] = c   };
        case '\\'           : return (StrText) { .s = "\\\\" };
        case '\t'           : return (StrText) { .s = "\\t"  };
        case '\n'           : return (StrText) { .s = "\\n"  };
        case '\r'           : return (StrText) { .s = "\\r"  };
        case '\b'           : return (StrText) { .s = "\\b"  };
        case 0 ... 7        : return (StrText) { .s[0] = '\\', .s[1] = '0'+c };
        
        default             : { // unprintable - output eg \xf 
            StrText p = {};
            snprintf (p.s, sizeof(p.s), "\\x%x", (uint8_t)c);
            return p;
        }
    }
}

// same as char_to_printable, but escaped for json
StrText char_to_printable_json (char c)  
{
    switch ((uint8_t)c) {
        case 32 ... '"' - 1 :
        case '"'+1 ... '\\'-1:
        case '\\'+1 ... 126 : return (StrText) { .s[0] = c     };
        case '"'            : return (StrText) { .s = "\\\""   };
        case '\\'           : return (StrText) { .s = "\\\\"   };
        case '\t'           : return (StrText) { .s = "\\\\t"  };
        case '\n'           : return (StrText) { .s = "\\\\n"  };
        case '\r'           : return (StrText) { .s = "\\\\r"  };
        case '\b'           : return (StrText) { .s = "\\\\b"  };
        case 0 ... 7        : return (StrText) { .s = { '\\', '\\', '0'+c } };
        
        default             : { // unprintable - output eg \xf 
            StrText p = {};
            snprintf (p.s, sizeof(p.s), "\\x%x", (uint8_t)c);
            return p;
        }
    }
}

// returns length (excluding \0). out should be allocated by caller to (in_len*4 + 1), out is null-terminated
uint32_t str_to_printable (STR𐤐(in), char *restrict s, int s_len)
{
    char *start = s;

    for (uint32_t i=0; i < in_len && s_len > 5; i++) // 5 = 4 characters + nul
        switch (in[i]) {
            case 32 ... '\\'-1: case '\\'+1 ... 126:
                           *s++ = in[i];                   ; s_len -= 1; break;
            case '\t'    : *s++ = '\\'; *s++ = 't'         ; s_len -= 2; break;
            case '\n'    : *s++ = '\\'; *s++ = 'n'         ; s_len -= 2; break;
            case '\r'    : *s++ = '\\'; *s++ = 'r'         ; s_len -= 2; break;
            case '\b'    : *s++ = '\\'; *s++ = 'b'         ; s_len -= 2; break;
            case '\\'    : *s++ = '\\'; *s++ = '\\'        ; s_len -= 2; break;
            case 0 ... 7 : *s++ = '\\'; *s++ = '0' + in[i] ; s_len -= 2; break;
            default      : *s++ = '\\'; *s++ = 'x' ; 
                           *s++ = NUM2HEXDIGIT((uint8_t)in[i] >> 4); 
                           *s++ = NUM2HEXDIGIT((uint8_t)in[i] & 0xf); 
                           s_len -= 4; 
        }
    
    *s = 0;
    return s - start;
}

uint32_t str_to_printable_json (STR𐤐(in), char *restrict s, int s_len)
{
    char *start = s;

    for (uint32_t i=0; i < in_len && s_len > 6; i++) // 6 = 5 characters + nul
        switch (in[i]) {
            case 32 ... '"' - 1 : case '"'+1 ... '\\'-1: case '\\'+1 ... 126:
                           *s++ = in[i];           ; s_len -= 1; break;
            case '\t'    : memcpy (s, "\\\\t", 3)  ; s_len -= 3; break;
            case '\n'    : memcpy (s, "\\\\n", 3)  ; s_len -= 3; break;
            case '\r'    : memcpy (s, "\\\\r", 3)  ; s_len -= 3; break;
            case '\b'    : memcpy (s, "\\\\b", 3)  ; s_len -= 3; break;
            case '\\'    : *s++ = '\\'; *s++ = '\\'; s_len -= 2; break;
            case '"'     : *s++ = '\\'; *s++ = '"' ; s_len -= 2; break;
            case 0 ... 7 : *s++ = '\\'; *s++ = '\\'; *s++ = '0' + in[i] ; s_len -= 3; break;
            default      : *s++ = '\\'; *s++ = '\\'; *s++ = 'x' ; 
                           *s++ = NUM2HEXDIGIT((uint8_t)in[i] >> 4); 
                           *s++ = NUM2HEXDIGIT((uint8_t)in[i] & 0xf); 
                           s_len -= 5; 
        }
    
    *s = 0;
    return s - start;
}

StrText str_size (uint64_t size)
{
    StrText s = {};

    if      (size >= 1 PB) snprintf (s.s, sizeof (s.s), "%3.1lf PB", ((double)size) / (double)(1 PB));
    else if (size >= 1 TB) snprintf (s.s, sizeof (s.s), "%3.1lf TB", ((double)size) / (double)(1 TB));
    else if (size >= 1 GB) snprintf (s.s, sizeof (s.s), "%3.1lf GB", ((double)size) / (double)(1 GB));
    else if (size >= 1 MB) snprintf (s.s, sizeof (s.s), "%3.1lf MB", ((double)size) / (double)(1 MB));
    else if (size >= 1 KB) snprintf (s.s, sizeof (s.s), "%3.1lf KB", ((double)size) / (double)(1 KB));
    else if (size >  0          ) snprintf (s.s, sizeof (s.s), "%3d B"    ,     (int)size)                       ;
    else                          snprintf (s.s, sizeof (s.s), "-"                       )                       ;

    return s;
}

StrText str_bases (uint64_t num_bases)
{
    StrText s = {};

    if      (num_bases >= 1000000000000000) snprintf (s.s, sizeof (s.s), "%.1lf Pb", ((double)num_bases) / 1000000000000000.0);
    else if (num_bases >= 1000000000000)    snprintf (s.s, sizeof (s.s), "%.1lf Tb", ((double)num_bases) / 1000000000000.0);
    else if (num_bases >= 1000000000)       snprintf (s.s, sizeof (s.s), "%.1lf Gb", ((double)num_bases) / 1000000000.0);
    else if (num_bases >= 1000000)          snprintf (s.s, sizeof (s.s), "%.1lf Mb", ((double)num_bases) / 1000000.0);
    else if (num_bases >= 1000)             snprintf (s.s, sizeof (s.s), "%.1lf Kb", ((double)num_bases) / 1000.0);
    else if (num_bases >  0   )             snprintf (s.s, sizeof (s.s), "%u b",    (uint32_t)num_bases);
    else                                    snprintf (s.s, sizeof (s.s), "-");

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
    return len; // excluding the \0
}

StrText str_int_s (int64_t n)
{
    StrText s = {};
    str_int (n, s.s);
    return s;
}

StrText1K str_int_s_(rom label/*max 59 chars*/, int64_t n)
{
    StrText1K s = {};
    int label_len = strlen (label);
    ASSERT (label_len < sizeof(StrText1K)-20/*max len of int64*/ -1/*\0*/, "label \"%s\" it too long (len=%u)", STRa(label)); 

    memcpy (s.s, label, label_len);
    str_int (n, &s.s[label_len]);
    
    return s;
}

StrText1K str_str_s_(rom label, STRp(str))
{
    StrText1K s = {};
    int label_len = strlen (label);

    bool truncated = (label_len + str_len > sizeof(s)-1);

    if (truncated)
        snprintf (s.s, sizeof (s.s), "%.17s%.*s...", label, (int)(sizeof(s)-4-label_len), str);
    else {
        memcpy (s.s, label, label_len);
        memcpy (s.s + label_len, str, str_len);
        *(s.s + label_len + str_len) = 0;
    }

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

// returns true if string is a valid "simple" float (i.e. not exponential notation, no leading zeros, must begin and end with a digit, may be an integer)
bool str_is_simple_float (STRp(str), 
                          uint32_t *decimals) // optional out, set only if true is returned   
{ 
    // enforce: leading zeros not allowed; cannot start or end with a period
    if (!str_len || str[0] == '0' || str[0] == '.' || str[str_len-1] == '.') 
        return false; 
    
    bool period_encountered = false;
    for (uint32_t i=0; i < str_len; i++) {
        if (!IS_DIGIT(str[i]) && str[i] != '.') return false; 
        
        if (str[i] == '.') {
            if (period_encountered) return false; // enforce: at most one period allowed
            if (decimals) *decimals = str_len - i - 1;
            period_encountered = true;
        }
    } 

    if (!period_encountered && decimals)
        *decimals = 0;
        
    return true;
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

    if (negative) out = -out;

    if (value) *value = out; // update only if successful
    return true;
}

#define str_get_int_range_type(func_num,type) \
bool str_get_int_range##func_num (rom𐤐 str, uint32_t str_len, int64_t min_val, int64_t max_val, type *restrict value /* optional */) \
{                                                                                     \
    if (!str) return false;                                                           \
                                                                                      \
    int64_t value64;                                                                  \
    if (!str_get_int (str, str_len ? str_len : strlen (str), &value64)) return false; \
    if (value) *value = (type)value64;                                                \
                                                                                      \
    return IN_RANGX (value64, min_val, max_val);                                      \
}
str_get_int_range_type(8,uint8_t)   // unsigned
str_get_int_range_type(16,uint16_t) // unsigned
str_get_int_range_type(32,int32_t)  // signed
str_get_int_range_type(64,int64_t)  // signed

bool str_get_uint16 (STRp(str), uint16_t *restrict value)
{
    int64_t v64;
    if (!str_get_int_range64 (STRa(str), 0, 0xffffffff, &v64)) return false;

    *value = v64;
    return true;
}

bool str_get_uint32 (STRp(str), uint32_t *restrict value)
{
    int64_t v64;
    if (!str_get_int_range64 (STRa(str), 0, 0xffffffff, &v64)) return false;

    *value = v64;
    return true;
}

// get a positive decimal integer, may have leading zeros eg 005
bool str_get_int_dec (STR𐤐(str), 
                      uint64_t *restrict value) // out - modified only if str is an integer
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
bool str_get_int_hex (STR𐤐(str), bool allow_hex, bool allow_HEX,
                      uint64_t *restrict value) // out - modified only if str is an integer
{
    uint64_t out = 0;

    for (uint32_t i=0; i < str_len; i++) {
        uint64_t prev_out = out;
        char c = str[i];
        
        if      (c >= '0' && c <= '9') out = (out << 4) | (c - '0');
        else if (allow_HEX && c >= 'A' && c <= 'F') out = (out << 4) | (c - 'A' + 10);
        else if (allow_hex && c >= 'a' && c <= 'f') out = (out << 4) | (c - 'a' + 10);
        else    
            return false;

        if (out < prev_out) 
            return false; // number overflowed beyond maximum uint64_t
    }

    if (value) *value = out; // update only if successful
    return true;
}

// positive integer, may be hex prefixed with 0x
#define str_get_int_range_allow_hex_bits(bits) \
bool str_get_int_range_allow_hex##bits (rom𐤐 str, uint32_t str_len, uint##bits##_t min_val, uint##bits##_t max_val, uint##bits##_t *restrict value) \
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
rom str_to_hex_(bytes𐤐 data, uint32_t data_len, char *restrict hex_str, bool with_dot)
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

StrText str_int_commas (int64_t n)
{
    StrText s = {};
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

    StrText s = {};

    if      (n >= 1000000000000000ULL) snprintf (s.s, sizeof (s.s), "%3.1lfP", ((double)n) / 1000000000000000.0);
    else if (n >=    1000000000000ULL) snprintf (s.s, sizeof (s.s), "%3.1lfT", ((double)n) / 1000000000000.0);
    else if (n >=       1000000000ULL) snprintf (s.s, sizeof (s.s), "%3.1lfG", ((double)n) / 1000000000.0);
    else if (n >=          1000000ULL) snprintf (s.s, sizeof (s.s), "%3.1lfM", ((double)n) / 1000000.0);
    else if (n >=                1000) snprintf (s.s, sizeof (s.s), "%3.1lfK", ((double)n) / 1000.0);
    else if (n >                    0) snprintf (s.s, sizeof (s.s), "%3d",     (int)n);
    else                               snprintf (s.s, sizeof (s.s), "-");

    return s;
}

// returns double value and/or format: "3.123" -> "%5.3f" 
bool str_get_float (STR𐤐(float_str), 
                    double *restrict value, char format[FLOAT_FORMAT_LEN], uint32_t *restrict format_len) // optional outs (format allocated by caller)
{
    bool in_decimals=false;
    bool is_negative = (float_str[0] == '-');
    int num_decimals=0, exponent=0;
    double val = 0;

    for (uint32_t i=is_negative; i < float_str_len; i++) {
        if (IS_DIGIT (float_str[i])) {
            val = (val * 10) + (float_str[i] - '0');
            if (in_decimals) num_decimals++;
        }
        
        else if (!in_decimals && float_str[i] == '.')
            in_decimals = true;
    
        // format [eE][+-][0-9][0-9]
        else if ((float_str[i] == 'e' || float_str[i] == 'E') &&
                  i+4 == float_str_len && 
                  (float_str[i+1] == '-' || float_str[i+1] == '+') &&
                  IS_DIGIT (float_str[i+2]) && IS_DIGIT (float_str[i+3])) {
            
            exponent = (float_str[i+2]-'0') * 10 + (float_str[i+3]-'0');
            
            if (float_str[i+1] == '-') exponent = -exponent;
            break;
        }
    
        else return false; // can't interpret this string as float
    }

    // calculate format if requested - the format string is in a format expected by reconstruct_from_local_float
    if (format) {
        if (float_str_len > 99) return false; // we support format of float strings up to 99 characters... more than enough
        uint32_t next=0;
        format[next++] = '%';
        
        if (float_str_len >= 10) {
            format[next++] = '0' + (float_str_len / 10);
            format[next++] = '0' + (float_str_len % 10);
        }
        else
            format[next++] = '0' + float_str_len;

        format[next++] = '.';

        if (num_decimals >= 10) {
            format[next++] = '0' + (num_decimals / 10);
            format[next++] = '0' + (num_decimals % 10);
        }
        else
            format[next++] = '0' + num_decimals;

        format[next++] = exponent ? float_str[float_str_len-4]/*e or E*/ : 'f';
        format[next] = 0;

        if (format_len) *format_len = next;
    }

    if (value) {
        static const double pow10[100] = { 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22, 1e23, 1e24, 1e25, 1e26, 1e27, 1e28, 1e29, 1e30, 1e31, 1e32, 1e33, 1e34, 1e35, 1e36, 1e37, 1e38, 1e39, 1e40, 1e41, 1e42, 1e43, 1e44, 1e45, 1e46, 1e47, 1e48, 1e49, 1e50, 1e51, 1e52, 1e53, 1e54, 1e55, 1e56, 1e57, 1e58, 1e59, 1e60, 1e61, 1e62, 1e63, 1e64, 1e65, 1e66, 1e67, 1e68, 1e69, 1e70, 1e71, 1e72, 1e73, 1e74, 1e75, 1e76, 1e77, 1e78, 1e79, 1e80, 1e81, 1e82, 1e83, 1e84, 1e85, 1e86, 1e87, 1e88, 1e89, 1e90, 1e91, 1e92, 1e93, 1e94, 1e95, 1e96, 1e97, 1e98, 1e99 };
                                        
        if (num_decimals >= ARRAY_LEN (pow10)) return false; // too many decimals

        *value = (is_negative ? -1 : 1) * (val / pow10[num_decimals]);

        if      (exponent > 0) *value *= pow10[exponent];
        else if (exponent < 0) *value /= pow10[-exponent];
    }

    return true;
}

bool str_is_in_range (rom str, uint32_t str_len, char first_c, char last_c)
{
    for (rom after = str + str_len; str < after ; str++)
        if (!IN_RANGX(*str, first_c, last_c)) return false;
    return true;
}

// returns true if the strings are case-insensitive similar, and *identical if they are case-sensitive identical
bool str_case_compare (rom str1, rom str2,
                       bool *identical) // optional out
{
    bool my_identical = true; // optimistic (automatic var)
    uint32_t len = strlen (str1);

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
        if (enforce_msg || flag.debug_split) {   \
            progress_newline(); fprintf (stderr, "Error in %s:%u: ", __FUNCLINE); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "%s", report_support_if_unexpected()); \
            if (enforce_msg) exit_on_error(true); /* same as ASSERT */ \
        } \
        return 0; /* just return 0 if we're not asked to enforce */ \
    } \
})

// splits a string with up to (max_items-1) separators (doesn't need to be nul-terminated) to up to or exactly max_items
// returns the actual number of items, or 0 is unsuccessful
uint32_t str_split_do (STR𐤐(str), 
                       uint32_t max_items,    // optional - if not given, a count of sep is done first
                       char sep,              // separator
                       rom𐤐 *restrict items,  // out - array of char* of length max_items - one more than the number of separators
                       uint32_t *restrict item_lens, // optional out - corresponding lengths
                       bool exactly,
                       rom𐤐 enforce_msg)      // non-NULL if enforcement of length is requested
{
    // IMPORTANT: restrict: since items[] elements do actually alias str, they should never be dereferenced!

    if (!str) return 0; // note: str!=NULL + str_len==0 results in n_items=1 - one empty item
    
    items[0] = str;
    uint32_t item_i = 1;
    for (uint32_t str_i=0 ; str_i < str_len ; str_i++) 
        if (str[str_i] == sep) {
            ASSSPLIT (item_i < max_items, "expecting up to %u %s separators but found more: (100 first) %.*s", 
                      max_items-1, enforce_msg ? enforce_msg : "", MIN_(str_len, 100), str);

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

// get item_i in an array. true if successful, false if array is too short
bool str_item_i (STR𐤐(str), char sep, uint32_t requested_item_i, 𐤐STR𐤐(item))
{
    // IMPORTANT: restrict: since item does actually alias str, it should never be dereferenced!

    if (!str) return false; 

    *item = str; // initialize, also to avoid compile warning
    
    uint32_t item_i = 1;
    for (uint32_t str_i=0 ; str_i < str_len ; str_i++) 
        if (str[str_i] == sep) {
            if (item_i == requested_item_i)
                *item = &str[str_i+1];
            
            else if (item_i == requested_item_i + 1) {
                *item_len = &str[str_i] - *item;
                return true;
            }

            item_i++;
        }

    if (item_i == requested_item_i + 1) {
        *item_len = &str[str_len] - *item;
        return true;
    }

    return false; 
}

// get item i from an array - expected to be a float. returns false if array is too short, or item is not a float
bool str_item_i_int (STRp(str), char sep, uint32_t requested_item_i, int64_t *restrict item)
{
    STR(item_str); // note: aliases str

    if (!str_item_i (STRa(str), sep, requested_item_i, pSTRa(item_str)))
        return false; // array too short

    return str_get_int (STRa(item_str), item);
}

// get item i from an array - expected to be a float. returns false if array is too short, or item is not a float
bool str_item_i_float (STRp(str), char sep, uint32_t requested_item_i, double *restrict item)
{
    STR(item_str); // note: aliases str

    if (!str_item_i (STRa(str), sep, requested_item_i, pSTRa(item_str)))
        return false; // array too short

    return str_get_float (STRa(item_str), item, NULL, NULL); // strtod until 15.0.83
}

// splits a string by tab, ending with the first \n or \r\n. 
// Returns the address of the byte after the \n if successful, or NULL if not.
rom str_split_by_tab_do (STR𐤐(str), 
                         uint32_t *restrict n_flds,  // in / out
                         rom𐤐 *restrict flds, uint32_t *restrict fld_lens, // out - array of char* of length max_flds - one more than the number of separators
                         bool *restrict has_13,      // optional out - true if line is terminated by \r\n instead of \n
                         bool exactly,      // line should have at_least this number of flds (or exactly if allow_more=false)
                         bool ignore_excess,// if true, ignore fields beyond n_flds. if false, its an error if there are more the n_flds
                         bool enforce_msg)
{
    // IMPORTANT: restrict: since flds[] elements do actually alias str, they should never be dereferenced!

    ASSERTNOTNULL (str);
    ASSERTNOTZERO (*n_flds);

    flds[0] = str;
    uint32_t fld_i = 1; 
    uint32_t str_i;
    bool my_has_13;

    for (str_i=0 ; str_i < str_len ; str_i++) {
        char c = str[str_i]; 
        if (c == '\t') {
            if (fld_i >= *n_flds)  // excess
                ASSSPLIT (ignore_excess, "expecting up to %u fields but found more. str_len=%u str(first 1000)=\"%.*s\"", *n_flds, str_len, MIN_(1000, str_len), str);
            else
                flds[fld_i++] = &str[str_i+1];
        }

        else if (c == '\r') {
            ASSSPLIT (str_i+1 < str_len && str[str_i+1] == '\n', "encountered a \\r without a following \\n. str_i=%u str_len=%u str(first 1000)=\"%.*s\" last_20_until_str_i=%s", 
                      str_i, str_len, MIN_(1000, str_len), str, str_to_printable_(&str[str_i-MIN_(str_i,20)+1], MIN_(str_i,20)).s);
            my_has_13 = true;
            break;
        }

        else if (c == '\n') {
            my_has_13 = false;
            break;
        }
    }

    ASSSPLIT (str_i < str_len, "Line not terminated by newline. str_len=%u str(first 1000)=\"%.*s\"", str_len, MIN_(1000, str_len), str);

    for (uint32_t i=0; i < fld_i-1; i++)    
        fld_lens[i] = flds[i+1] - flds[i] - 1; 
        
    fld_lens[fld_i-1] = &str[str_i] - flds[fld_i-1];

    ASSSPLIT (!exactly || fld_i == *n_flds, "expecting %u fields but found more. str=\"%.*s\"", *n_flds, str_len, str);

    *n_flds = fld_i;

    if (has_13) *has_13 = my_has_13;

    return &str[str_i + 1 + my_has_13]; // byte after \n
}

// get up to n_lines from str. ignores subsequent lines.
// returns lines, with their lengths excluding \n and \r
uint32_t str_split_by_lines_do (STR𐤐(str), uint32_t max_lines, rom𐤐 *restrict lines, uint32_t *restrict line_lens)
{
    // IMPORTANT: restrict: since lines[] elements do actually alias str, they should never be dereferenced!

    rom after = str + str_len;
    uint32_t i; for (i=0; i < max_lines && str < after; i++) {
        rom nl = memchr (str, '\n', after - str);
        if (nl) {
            lines[i] = str;
            line_lens[i] = nl - str;
            if (nl != str && nl[-1] == '\r') line_lens[i]--;

            str = nl + 1;
        }
        else 
            break; // we're done
    }

    return i;
}

// splits a string based on container items (doesn't need to be nul-terminated). 
// returns the number on unskipped items if successful
uint32_t str_split_by_container_do (STR𐤐(str), ConstContainer𐤐 con, STR𐤐(con_prefixes),
                                    rom𐤐 *restrict items, // out - array of char* of length max_items - one more than the number of separators
                                    uint32_t *restrict item_lens, // optional out - corresponding lengths
                                    rom𐤐 enforce_msg)     // non-NULL if enforcement of length is requested
{
    // IMPORTANT: restrict: since items[] elements do actually alias str, they should never be dereferenced!

    if (!str) return 0; // note: str!=NULL + str_len==0 results in n_items=1 - one empty item

    uint32_t num_items = con_nitems (*con);
    ASSERT0 (num_items, "Container has no items");

    const rom𐤐 *restrict save_items = items;
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
                    ASSSPLIT (prefix_matches, "prefix mismatch for item_i=%d: expecting \"%.*s\" but seeing \"%.*s\"",
                              item_i-1, item_px_len, item_px, MIN_(item_px_len, (int)(after_str-str)), str);
    
                    str += item_px_len; // advance past prefix
                }
            }
        }

        switch (sep) {
            case 0: long_var_0_pad: 
                // case: no separator and next item has no prefix - goes to end of string
                if (!(px < after_px && IS_PRINTABLE(*px))) { // make sure we don't have an implicit separator - the prefix of the next item
                    *item_lens = after_str - str;
                    *items = str; 
                    str = after_str;
                }

                // case: terminated by prefix of next item
                else {
                    rom c; for (c = px; *c != CON_PX_SEP; c++);
                    uint32_t next_px_len = c - px;
                    
                    *items = str;

                    // increment str to the first character of the next item's prefix 
                    while (str < after_str-next_px_len && memcmp (str, px, next_px_len)) str++; 

                    ASSSPLIT (str < after_str-next_px_len || !memcmp (str, px, next_px_len), 
                            "item_i=%u reached end of string without finding next item prefix \"%.*s\" in string \"%.*s\"", 
                            item_i, next_px_len, px, str_len, save_str);

                    *item_lens = str - *items;
                }
                break;

            case 32 ... 126: case '\t': case '\n': case '\r': { // printable
                char sep1 = con->items[item_i].separator[1];
                if (!IS_PRINTABLE (sep1)) sep1 = 0;
                
                *items = str;

                if (!sep1) 
                    while (str < after_str && *str != sep) str++; 
                else 
                    while (str < (after_str-1) && (str[0] != sep || str[1] != sep1)) str++; 

                ASSSPLIT (str < after_str, "item_i=%u searched remainder of string: \"%.*s\" but could not find separator '%c'. Entire string is \"%.*s\"", 
                        item_i, (int)(after_str - *items), *items, sep, str_len, save_str);

                *item_lens = str - *items;
                str += 1 + (sep1 != 0); // skip seperator

                break;
            }

            case CI0_SKIP: 
            case CI0_INVISIBLE:
                *item_lens = 0;
                *items = str; // zero-length item - next item will start from the same str_i
                break;

            case CI0_VAR_0_PAD:
                // case: e.g. 4-pad but item is 10342 (5 digits)
                if (after_str - str > con->items[item_i].separator[1]) goto long_var_0_pad;
                // fallthrough

            case CI0_FIXED_0_PAD: 
                *item_lens = con->items[item_i].separator[1];
                *items = str; // zero-length item - next item will start from the same str_i
                str += *item_lens;
                ASSSPLIT (str <= after_str, "item_i=%u fixed_len=%u goes beyond end of string \"%.*s\"", item_i, *item_lens, str_len, save_str);
                break;
            
            case CI0_DIGIT: // items goes until a digit or end-of-string is encountered
                *items = str;

                while (str < after_str && !IS_DIGIT(*str)) str++;  

                ASSSPLIT (str < after_str, "item_i=%u reached end of string without finding a DIGIT separator \"%.*s\"", 
                          item_i, str_len, save_str);

                *item_lens = str - *items;
                break;

            case CI0_ACGTN: // items goes until a A,C,G,T or N, but not less than sep[1] characters
                *items = str;

                str += (uint8_t)con->items[item_i].separator[1]; // minimum number of character
                while (str < after_str && !IS_ACGTN(*str)) str++;  

                ASSSPLIT (str < after_str, "item_i=%u reached end of string without finding a ACGTN separator \"%.*s\"", 
                          item_i, str_len, save_str);

                *item_lens = str - *items;
                break;
            
            case CI0_COLONn: 
                sep = ':'; goto multisep;
            
            multisep: {
                int n = (uint8_t)con->items[item_i].separator[1];
                
                *items = str;

                for (int i=0; i < n; i++, str++/*skip sep*/) {
                    while (str < (after_str-1) && (str[0] != sep)) str++; 

                    ASSSPLIT (str < after_str, "item_i=%u reached end of string without finding %u separators '%c' in string \"%.*s\"", 
                              item_i, n, sep, str_len, save_str);
                }

                *item_lens = str - *items -1/*excluding final separator*/;
                break;
            }

            case CI0_LAST_MATCH: { // items goes until LAST match in string of sep1
                char sep1 = con->items[item_i].separator[1];
                *items = str;

                for (str=after_str-1; str >= *items && *str != sep1; str--);
                ASSSPLIT (str >= *items, "item_i=%u reached end of string (LAST_MATCH) without finding separator '%c' in string \"%.*s\"", 
                          item_i, sep, str_len, save_str);

                *item_lens = str - *items;
                str++; // skip separator
                break;
            }

            default: 
                ABORT ("item_i=%u sep=%u is not a printable character in string \"%.*s\"", item_i, sep, str_len, save_str);
        }

        // printf ("item_i=%u item=%.*s\n", item_i, *item_lens, *items);
        
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
void str_nul_separate_do (STR𐤐s(item))
{
    for (uint32_t i=0; i < n_items; i++)
        ((char**)items)[i][item_lens[i]] = '\0';
}

// generate new string (mem allocated by caller) - copy of in, but removing all ASCII < 33 or > 126 
uint32_t str_remove_whitespace (STRp(in), bool also_uppercase, char *out)
{
    uint32_t out_len = 0;
    for (uint32_t i=0; i < in_len; i++)
        if (IS_NON_WS(in[i])) 
            out[out_len++] = also_uppercase ? UPPER_CASE(in[i]) : in[i];

    return out_len;
}

// in-place removal of flanking whitespace from a null-terminated string
void str_trim (qSTRp(str))
{
    // remove leading whitespace
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
uint32_t str_split_ints_do (STR𐤐(str), uint32_t max_items, char sep, bool exactly, int base,
                            int64_t *restrict items)  // out - array of integers
                       
{
    rom after = &str[str_len];
    SAFE_NUL (after); // this doesn't work on string literals

    uint32_t item_i;
    for (item_i=0; item_i < max_items && str < after; item_i++, str++) {
        items[item_i] = (base==16) ? strtoull (str, (char **)&str, base) // allows parsing 64bit (unsigned) hexes
                                   : strtoll (str, (char **)&str, base);

        if (item_i < max_items-1 && *str != sep && *str != 0) {
            item_i=0;
            break; // fail - not integer
        }
    }

    if (str < after) item_i = 0;

    SAFE_RESTORE;

    return (!exactly || item_i == max_items) ? item_i : 0; // 0 if requested exactly, but too few separators 
}

// splits a string with up to (max_items-1) separators (doesn't need to be nul-terminated) to up to or exactly max_items floats
// returns the actual number of items, or 0 is unsuccessful
// WARNING: this function is slow! because strtod is very slow (atof is just as slow).
uint32_t str_split_floats_do (STR𐤐(str), uint32_t max_items, char sep, bool exactly, 
                              char this_char_is_NAN,   // optional (if not 0): if this character is encountered, then the floating point value is NAN
                              double *restrict items)  // out - array of floats                 
{
    rom after = &str[str_len];
    SAFE_NUL (after); // this doesn't work on string literals

    uint32_t item_i;
    for (item_i=0; item_i < max_items && str < after; item_i++, str++) {
        if (this_char_is_NAN && *str == this_char_is_NAN && (str[1] == sep || str[1] == 0)) {
            items[item_i] = NAN;
            str++;    
        }
        
        else {
            rom after_item = memchr (str, sep, after - str);
            if (!after_item) after_item = after;

            if (!str_get_float (str, after_item - str, &items[item_i], 0, 0)) goto error;
            str = after_item;
        }
    }

    if (str < after) item_i = 0;

    SAFE_RESTORE;
    return (!exactly || item_i == max_items) ? item_i : 0; // 0 if requested exactly, but too few separators 

error:
    SAFE_RESTORE;
    return 0;
}

rom type_name (uint32_t item, 
               rom const *name, // the address in which a pointer to name is found, if item is in range
               uint32_t num_names)
{
    if (item > num_names) {
        static char str[50];
        snprintf (str, sizeof(str), "%d (out of range)", item);
        return str;
    }
    
    return *name;    
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
                        c>=0 && (IS_ALPHANUMERIC(line[c]) || line[c]==','); // wrap lines at - and | too, so we can break very long regex strings like in genocat
                        c--); // find 
            iprintf ("%s%.*s\n", wrapped ? wrapped_line_prefix : "", c, line);
            STRinc (line, c + (line[c]==' ')/*skip space too*/);     
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
    int query_str_len = strlen(query)+32; 
    char query_str[query_str_len];
    snprintf (query_str, query_str_len, "%s (%sy%s or %sn%s) ", query, 
             def_answer==QDEF_YES?"[":"",  def_answer==QDEF_YES?"]":"",  
             def_answer==QDEF_NO ?"[":"",  def_answer==QDEF_NO ?"]":"");

    char y_n[32] = {};
    str_query_user (query_str, y_n, sizeof(y_n), def_answer != QDEF_NONE,
                    str_verify_y_n, def_answer==QDEF_YES ? "Y" : def_answer==QDEF_NO ? "N" : 0);

    return y_n[0] == 'Y';
}

void str_query_user (rom query, 
                     char *response, uint32_t response_size, // if response is not empty, it is used as the default
                     bool allow_empty, 
                     ResponseVerifier verifier, rom verifier_param)
{
    uint32_t len = 0;

    char save_response[response_size];
    memcpy (save_response, response, response_size);

    do {
        fprintf (stderr, "%s", query);
        if (save_response[0]) fprintf (stderr, "[%s] ", save_response);

        // Linux: in case of non-Latin, fgets doesn't handle backspace well - not completely removing the multi-byte character from the string (bug 593)
        // Windows (MingW): the string is terminated before the first non-Latin character on bash terminal and only ???? on Powershell and Command (bug 594)
        do {
            if (fgets (response, response_size, stdin)) { // success
                len = strlen (response);
                str_trim (response, &len);
            }
        } while (!len && !allow_empty && !save_response[0]);

        if (!len && save_response[0]) {
            memcpy (response, save_response, response_size);
            len = strlen (response);
        }

    } while (verifier && !verifier (response, len, verifier_param));
}

rom str_win_error_(uint32_t error)
{
    static char msg[100];
#ifdef _WIN32
    FormatMessageA (FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,                   
                    NULL, error, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), msg, sizeof (msg), NULL);
#endif
    return msg;
}

rom str_win_error (void)
{
#ifdef _WIN32
    uint32_t error = GetLastError();
    return str_win_error_(error);
#else
    return "";
#endif
}

// C<>G A<>T c<>g a<>t ; IUPACs: R<>Y K<>M B<>V D<>H W<>W S<>S N<>N (+ lowercase); other ASCII 32->126 preserved ; other = 0
alignas(64) const char COMPLEM[256] = "-------------------------------- !\"#$%&'()*+,-./0123456789:;<=>?@TVGHEFCDIJMLKNOPQYSAUBWXRZ[\\]^_`tvghefcdijmlknopqysaubwxrz{|}~";

// reverse-complements a string (possibly in-place)
char *str_revcomp (char *dst_seq, rom src_seq, uint32_t seq_len)
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

// complements uppercase A,C,G,T only - everything else is left intact
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winitializer-overrides"
alignas(64) const char COMPLEM_ACGT[256] = { 
    0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255, // seq 0 256 | tr "\n" ,
    ['A']='T', ['C']='G', ['G']='C', ['T']='A'
};
#pragma clang diagnostic pop

// reverse-complements a string (possibly in-place) - complements A,C,G,T characters and leaves others intact. set dst_seq=src_seq to output in-place.
char *str_revcomp_ACGT (char *dst_seq, rom src_seq, uint32_t seq_len)
{
    for (uint32_t i=0; i < seq_len / 2; i++) {
        uint8_t l_base = src_seq[i];
        uint8_t r_base = src_seq[seq_len-1-i];

        dst_seq[i]           = COMPLEM_ACGT[r_base];
        dst_seq[seq_len-1-i] = COMPLEM_ACGT[l_base];
    }

    if (seq_len % 2) // we have an odd number of bases - now complement the middle one
        dst_seq[seq_len/2] = COMPLEM_ACGT[(uint8_t)src_seq[seq_len/2]];

    return dst_seq;
}

// implementing memrchr, as it doesn't exist in Windows libc (msvcrt.dll) or Darwin (at least clang)
#if defined _WIN32 || defined __APPLE__
void *memrchr (const void *s, int c/*interpreted as unsigned char*/, size_t n)
{
    for (uint8_t c8=c, *p=((uint8_t *)s)+n-1; p >= (uint8_t *)s; p--)
        if (*p == c8) return p;

    return NULL;
}
#endif

char *memchr2 (const void *p_, char ch1, char ch2, uint32_t count)
{
    rom p = (rom)p_;

    for (rom after = p + count; p < after; p++)
        if (*p == ch1 || *p == ch2) return (char *)p;

    return NULL;
}

// print duration in human-readable form eg 1h2' or "1 hour 2 minutes"
StrText str_human_time (unsigned secs, bool compact)
{
    StrText s;
    int s_len = 0;

    unsigned hours = secs / 3600;
    unsigned mins  = (secs % 3600) / 60;
             secs  = secs % 60;

    if (compact)
        SNPRINTF (s, "%uh%u'%u\\\"", hours, mins, secs); // escape the double-quotes
    else if (hours) 
        SNPRINTF (s, "%u %s %u %s", hours, hours==1 ? "hour" : "hours", mins, mins==1 ? "minute" : "minutes");
    else if (mins)
        SNPRINTF (s, "%u %s %u %s", mins, mins==1 ? "minute" : "minutes", secs, secs==1 ? "second" : "seconds");
    else 
        SNPRINTF (s, "%u %s", secs, secs==1 ? "second" : "seconds");

    return s;
}

// current date and time
StrText1K str_time (void)
{
    StrText1K s = {}; // long, in case of eg Chinese language time zone strings

#ifdef _WIN32
    int s_len = GetDateFormatA (LOCALE_USER_DEFAULT, DATE_SHORTDATE, NULL, NULL, s.s, sizeof(s.s)-1) - 1;
    s.s[s_len++] = ' ';
    s_len += GetTimeFormatA (LOCALE_USER_DEFAULT, TIME_FORCE24HOURFORMAT, NULL, NULL, &s.s[s_len], sizeof(s.s)-s_len-1) - 1;

    s.s[s_len++] = ' ';

    TIME_ZONE_INFORMATION tz_info = {};
    bool is_daylight = (GetTimeZoneInformation (&tz_info) == TIME_ZONE_ID_DAYLIGHT);

    SNPRINTF (s, "%ls", is_daylight ? tz_info.DaylightName : tz_info.StandardName);

#else
    time_t now = time (NULL);
    struct tm lnow;
    localtime_r (&now, &lnow);

    static rom month[12] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };
    snprintf (s.s, sizeof (s.s), "%02u-%s-%04u %02u:%02u:%02u %s", 
             lnow.tm_mday, month[lnow.tm_mon], lnow.tm_year+1900, lnow.tm_hour, lnow.tm_min, lnow.tm_sec, tzname[daylight != 0]);
#endif

    return s;
}

// check if string is UTF-8 based on Unicode version 15 (2022): Table 3-7 in https://www.unicode.org/versions/Unicode15.0.0/UnicodeStandard-15.0.pdf
bool str_is_utf8 (STRp(s))
{
    rom after = s + s_len;

    while (s < after) {
        #define r(i,first,last) ((uint8_t)s[i] >= (first) && (uint8_t)s[i] <= (last))
        if      (               r(0, 0x00, 0x7f)) s += 1;
        else if (s+1 < after && r(0, 0xC2, 0xDF) && r(1, 0x80, 0xBF)) s += 2;
        else if (s+2 < after && r(0, 0xE0, 0xE0) && r(1, 0xA0, 0xBF) && r(2, 0x80, 0xBF)) s += 3;
        else if (s+2 < after && r(0, 0xE1, 0xEC) && r(1, 0x80, 0xBF) && r(2, 0x80, 0xBF)) s += 3;
        else if (s+2 < after && r(0, 0xED, 0xED) && r(1, 0x80, 0x9F) && r(2, 0x80, 0xBF)) s += 3;
        else if (s+2 < after && r(0, 0xEE, 0xEF) && r(1, 0x80, 0xBF) && r(2, 0x80, 0xBF)) s += 3;
        else if (s+3 < after && r(0, 0xF0, 0xF0) && r(1, 0x90, 0xBF) && r(2, 0x80, 0xBF) && r(3, 0x80, 0xBF)) s += 4;
        else if (s+3 < after && r(0, 0xF1, 0xF3) && r(1, 0x80, 0xBF) && r(2, 0x80, 0xBF) && r(3, 0x80, 0xBF)) s += 4;
        else if (s+3 < after && r(0, 0xF4, 0xF4) && r(1, 0x80, 0x8F) && r(2, 0x80, 0xBF) && r(3, 0x80, 0xBF)) s += 4;
        else return false; // not UTF-8
        #undef r
    }

    return true;
}

#ifdef __MINGW64__
void *memmem (const void *haystack, size_t haystack_len, // same as in gcc
              const void *needle,   size_t needle_len)
{
    const void *after_haystack = (rom)haystack + haystack_len - needle_len + 1;
    for (; haystack < after_haystack; haystack=(rom)haystack + 1)
        if (!memcmp (haystack, needle, needle_len)) 
            return (void *)haystack;

    return NULL;
}
#endif

// add packed A,C,G,T data to a byte array. 
// NOTE: bases should have accessible memory for reading 3 bytes before and 3 bytes after (i.e. always ok if bases on in a Buffer)
uint32_t str_pack_bases (uint8_t *restrict packed, STR𐤐(bases), bool revcomp) // no aliasing!
{
    uint32_t num_bytes = roundup_bits2bytes (bases_len * 2);
    uint32_t bases_in_partial_final_byte = bases_len % 4; // 0 if final byte is not partial

    if (!revcomp)
        for (uint32_t byte_i=0; byte_i < num_bytes; byte_i++, bases += 4) {
            uint32_t w = *(unaligned_uint32_t *)bases; // single memory read - possibly reading beyond the end of bases. the garbage bases at the end become garbage bits in the final byte of the packed array

            // bits 1 and 2 of ASCII 'A','C','G','T' are 00, 01, 11, 10 respectively
            uint8_t b0 = (w >> 1 ) & 3;
            uint8_t b1 = (w >> 9 ) & 3;
            uint8_t b2 = (w >> 17) & 3;
            uint8_t b3 = (w >> 25) & 3;

            packed[byte_i] = b0 | (b1 << 2) | (b2 << 4) | (b3 << 6);            
        }

    else {
        bases -= bases_in_partial_final_byte ? (4 - bases_in_partial_final_byte) : 0; // possibly reading before the beginning of bases - the gargabe before becomes garbage bits in the final packed byte - cleared in the end
         
        for (int32_t byte_i=num_bytes - 1; byte_i >= 0; byte_i--, bases += 4) {
            uint32_t w = *(const unaligned_uint32_t *)bases;

            uint8_t b0 = ((w >> 1)  & 3) ^ 2; // ^2 complements: A(00)⇔T(10) ; C(01)⇔G(11)
            uint8_t b1 = ((w >> 9)  & 3) ^ 2; 
            uint8_t b2 = ((w >> 17) & 3) ^ 2; 
            uint8_t b3 = ((w >> 25) & 3) ^ 2; 

            packed[byte_i] = b3 | (b2 << 2) | (b1 << 4) | (b0 << 6);
        }
    }
    
    // clear garbage bits - due to processing of up to 3 characters after the end (if forward) or before the beginning (if reverse) of bases.
    if (bases_in_partial_final_byte)    
        packed[num_bytes-1] &= bitmask8 (2 * bases_in_partial_final_byte);
     
    return num_bytes;
}

// convert a bytes containing a series of 2-bits, a string of A,C,G,T
uint32_t str_unpack_bases (char *restrict dst, bytes𐤐 packed, uint32_t num_bases)
{
    char actg[4] = "ACTG"; // not ACGT!

    uint32_t full_bytes = num_bases / 4;

    for (uint32_t i=0; i < full_bytes; i++, dst += 4) { // 4 bases at a time
        uint8_t b = packed[i];

        uint32_t base0 = actg[(b >> 0) & 3];
        uint32_t base1 = actg[(b >> 2) & 3];
        uint32_t base2 = actg[(b >> 4) & 3];
        uint32_t base3 = actg[(b >> 6) & 3];

#if defined __LITTLE_ENDIAN__
        uint32_t word = base0 | (base1 << 8) | (base2 << 16) | (base3 << 24);
#elif defined __BIG_ENDIAN__
        uint32_t word = base3 | (base2 << 8) | (base1 << 16) | (base0 << 24);
#else
#error  "Neither __BIG_ENDIAN__ nor __LITTLE_ENDIAN__ is defined - is endianness.h included?"
#endif    
        *(unaligned_uint32_t *)dst = word; // single write for 4 bases
    }

    // final partial byte
    uint32_t bytes_consumed = (num_bases + 3) / 4;

    uint8_t b = packed[bytes_consumed-1]; // this way, so read the last bytes again if no remainder, instead of the next byte which might be beyond the end of the array
    
    switch (num_bases % 4) {
        case 3: dst[2] = actg[(b >> 4) & 3]; // fall through
        case 2: dst[1] = actg[(b >> 2) & 3]; // fall through
        case 1: dst[0] = actg[(b >> 0) & 3];        
        default: break;
    }

    return bytes_consumed; // data consumed in bytes
}

// true if all bytes of str are 0
// Note: only worth calling if str_len is usually >32, and there is a reasonable chance of true. Otherwise use str_is_monochar.
bool str_is_zero (STR𐤐(str))
{
#if defined(__AVX2__) // 32 bytes at a time
    uint32_t num_words = str_len / sizeof(__m256i); // 32 byte "words"

    if (num_words) {
        __m256i *data = (__m256i *)str;
        for (uint32_t i=0; i < num_words; i++) {
            __m256i word = _mm256_loadu_si256 (&data[i]); // Load 32 bytes from memory (unaligned)
            if (!_mm256_testz_si256(word, word)) // _mm256_testz_si256 is true if (word & word)==0, i.e. if word==0
                return false;
        }

        // final word (note: there is a 1/32 chance this test is redundant, but this is still less overhead than adding an if statement to avoid this redundancy)
        __m256i word = _mm256_loadu_si256 ((__m256i *)(str + str_len) - 1); // 32B which partially (or rarely, fully) overlap already tested last word
        return _mm256_testz_si256(word, word); 
    }

#else // 8 bytes a time
    uint32_t num_words = str_len / sizeof(uint64_t);

    if (num_words) {
        unaligned_uint64_t *data = (unaligned_uint64_t *)str; // note: buf->data is always word aligned
        for (uint32_t i=0; i < num_words; i++) 
            if (data[i] != 0) 
                return false;

        // final word
        return *(unaligned_uint64_t *)(str + str_len - sizeof(uint64_t)) == 0; // 8B which partially (or rarely, fully) overlap already tested last word
    }
#endif

    // num_words==0 - fallback to byte loop
    else
        return str_is_monochar_(STRa(str), 0);
}