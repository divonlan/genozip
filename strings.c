// ------------------------------------------------------------------
//   strings.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "strings.h"

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

char *str_uint (int64_t n, char *str /* out */, unsigned *len)
{
    *len=0;

    if (n==0) {
        str[0] = '0';
        *len=1;
    }

    else {
        char rev[50] = {}; // "initialize" to avoid compiler warning
        while (n) {
            rev[(*len)++] = '0' + n % 10;
            n /= 10;
        }
        // now reverse it
        for (int i=0; i < (*len); i++) str[i] = rev[(*len)-i-1];
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
