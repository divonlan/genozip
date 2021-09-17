// ------------------------------------------------------------------
//   strings.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef STRINGS_INCLUDED
#define STRINGS_INCLUDED

#include "genozip.h"

#define IS_NUCLEOTIDE(c) ((c) == 'A' || (c) == 'T' || (c) == 'C' || (c) == 'G')
#define IS_DIGIT(c)    ((c)>='0' && (c)<='9')
#define IS_HEXDIGIT(c) (((c)>='0' && (c)<='9') || ((c)>='A' && (c)<='F') || ((c)>='a' && (c)<='f'))
#define IS_CLETTER(c)  ((c)>='A' && (c)<='Z')
#define IS_SLETTER(c)  ((c)>='a' && (c)<='z')
#define IS_LETTER(c) (IS_CLETTER(c) || IS_SLETTER(c))
#define IS_NON_WS_PRINTABLE(c) (((c)>=33) && ((c)<=126))
#define IS_VALID_URL_CHAR(c) (IS_LETTER(c) || IS_DIGIT(c) || c=='-' || c=='_' || c=='.' || c=='~') // characters valid in a URL
#define FLIP_CASE(c) (IS_CLETTER(c) ? ((c)+32) : (IS_SLETTER(c) ? ((c)-32) : (c))) // flips lower <--> upper case
#define UPPER_CASE(c) (IS_SLETTER(c) ? ((c)-32) : (c))
#define LOWER_CASE(c) (IS_CLETTER(c) ? ((c)+32) : (c))

typedef struct { char s[80]; } StrText;

extern StrText char_to_printable (char c);

extern char *str_print_snip_do (const char *in, unsigned in_len, char *out);
#define str_print_snip(in,in_len) ({ struct { char s[in_len+20]; } str; str_print_snip_do (in, in_len, str.s); str; })

extern char *str_to_single_line_printable (const char *in, unsigned in_len, char *out);
extern char *str_tolower (const char *in, char *out /* out allocated by caller - can be the same as in */);
extern char *str_toupper (const char *in, char *out);

static _Bool inline str_issame_ (const char *str1, unsigned str1_len, const char *str2, unsigned str2_len) // true if the same
{
    return (str1_len == str2_len) && !memcmp (str1, str2, str1_len);
}
#define str_issame(str1,str2) str_issame_ (str1, str1##_len, str2, str2##_len)

extern _Bool str_case_compare (const char *str1, const char *str2, _Bool *identical); // similar to stricmp that doesn't exist on all platforms

static inline char *str_tolower_(const char *in, char *out, unsigned len)
{
    for (unsigned i=0; i < len; i++) out[i] = LOWER_CASE (in[i]); 
    return out;
}

static inline char *str_toupper_(const char *in, char *out, unsigned len)
{
    for (unsigned i=0; i < len; i++) out[i] = UPPER_CASE (in[i]);
    return out;
}

static inline char *str_reverse (const char *in, char *out, unsigned len)
{
    for (unsigned i=0; i < len; i++) out[len-1-i] = in[i]; 
    return out;
}

extern const char COMPLEM[], UPPER_COMPLEM[];
extern char *str_revcomp (char *seq, unsigned seq_len);

static inline uint64_t str_count_char (const char *str, uint64_t len, char c)
{
    if (!str) return 0;
    
    uint64_t count=0;
    for (uint64_t i=0; i < len; i++)
        if (str[i] == c) count++;

    return count;
}

extern StrText str_size (uint64_t size);
extern StrText str_bases (uint64_t num_bases);
extern StrText str_uint_commas (int64_t n);
extern StrText str_uint_commas_limit (uint64_t n, uint64_t limit);
extern StrText str_int_s (int64_t n);

// string length of an integer. #include <math.h> if using this.
#define str_int_len(n) ( (n)>=1          ? ((int)log10(n)) + 1     \
                       : (n)>-1 && (n)<1 ? 1                       \
                       :                   ((int)log10(-(n))) + 2)

extern unsigned str_int (int64_t n, char *str /* out */);
extern _Bool str_get_int (const char *str, unsigned str_len, int64_t *value); // note: for a reason beyond me, Docker hub won't compile if its "bool" and not "_Bool"
extern _Bool str_get_int_range8  (const char *str, unsigned str_len, uint8_t  min_val, uint8_t  max_val, uint8_t  *value);
extern _Bool str_get_int_range16 (const char *str, unsigned str_len, uint16_t min_val, uint16_t max_val, uint16_t *value);
extern _Bool str_get_int_range64 (const char *str, unsigned str_len, int64_t  min_val, int64_t  max_val, int64_t  *value);
extern _Bool str_get_int_range32 (const char *str, unsigned str_len, int32_t  min_val, int32_t  max_val, int32_t  *value);

extern _Bool str_get_int_hex (const char *str, unsigned str_len, uint64_t *value); 
extern _Bool str_get_int_range_allow_hex8  (const char *str, unsigned str_len, uint8_t  min_val, uint8_t  max_val, uint8_t  *value);
extern _Bool str_get_int_range_allow_hex16 (const char *str, unsigned str_len, uint16_t min_val, uint16_t max_val, uint16_t *value);
extern _Bool str_get_int_range_allow_hex32 (const char *str, unsigned str_len, uint32_t min_val, uint32_t max_val, uint32_t *value);
extern _Bool str_get_int_range_allow_hex64 (const char *str, unsigned str_len, uint64_t min_val, uint64_t max_val, uint64_t *value);

#define str_is_int(str,str_len) str_get_int ((str), (str_len), NULL)

extern _Bool str_is_in_range (const char *str, uint32_t str_len, char first_c, char last_c);

extern StrText str_pointer (const void *p);
extern StrText str_time (void);

#define FLOAT_FORMAT_LEN 12
extern _Bool str_get_float (const char *float_str, unsigned float_str_len, double *value, char format[FLOAT_FORMAT_LEN], unsigned *format_len);

extern _Bool str_scientific_to_decimal (const char *float_str, unsigned float_str_len, char *modified, unsigned *modified_len, double *value);

extern unsigned str_split_do (const char *str, unsigned str_len, uint32_t max_items, char sep, const char **items, unsigned *item_lens, _Bool exactly, const char *enforce_msg);

// name      : eg "item", macro defines variables "items" (array of pointers), item_lens (array or unsigned), n_items (actual number of items)
// max_items : maximum allowed items, or 0 if not known
#define str_split(str,str_len,max_items,sep,name,exactly) str_split_enforce((str),(str_len),(max_items),(sep),name,(exactly),NULL)

#define str_split_enforce(str,str_len,max_items,sep,name,exactly,enforce) \
    unsigned n_##name##s = (max_items) ? (max_items) : str_count_char ((str), (str_len), (sep)) + 1; /* 0 if str is NULL */ \
    const char *name##s[MAX_(n_##name##s, 1)]; \
    unsigned name##_lens[MAX_(n_##name##s, 1)]; \
    n_##name##s = str_split_do ((str), (str_len), n_##name##s, (sep), name##s, name##_lens, (exactly), (enforce)); 

extern void str_remove_CR_do (unsigned n_lines, const char **lines, unsigned *line_lens);
#define str_remove_CR(name) str_remove_CR_do (n_##name##s, name##s, name##_lens)

extern void str_nul_separate_do (unsigned n_items, const char **items, unsigned *item_lens);
#define str_nul_separate(name) str_nul_separate_do (n_##name##s, name##s, name##_lens)

extern unsigned str_remove_whitespace (const char *in, unsigned in_len, char *out);

extern unsigned str_split_ints_do (const char *str, unsigned str_len, uint32_t max_items, char sep, _Bool exactly, int64_t *items);
#define str_split_ints(str,str_len,max_items,sep,name,exactly) \
    unsigned n_##name##s = (max_items) ? (max_items) : str_count_char ((str), (str_len), (sep)) + 1; \
    int64_t name##s[n_##name##s]; \
    n_##name##s = str_split_ints_do ((str), (str_len), n_##name##s, (sep), (exactly), name##s); 

extern unsigned str_split_floats_do (const char *str, unsigned str_len, uint32_t max_items, char sep, _Bool exactly, double *items);
#define str_split_floats(str,str_len,max_items,sep,name,exactly) \
    unsigned n_##name##s = (max_items) ? (max_items) : str_count_char ((str), (str_len), (sep)) + 1; \
    double name##s[n_##name##s]; \
    n_##name##s = str_split_floats_do ((str), (str_len), n_##name##s, (sep), (exactly), name##s); 

extern const char *type_name (unsigned item, 
                              const char * const *name, // the address in which a pointer to name is found, if item is in range
                              unsigned num_names);

extern void str_print_dict (const char *data, unsigned len, _Bool add_newline, _Bool remove_equal_asterisk);

extern int str_print_text (const char **text, unsigned num_lines,
                           const char *wrapped_line_prefix, 
                           const char *newline_separator, 
                           unsigned line_width /* 0=calcuate optimal */);

typedef _Bool (*ResponseVerifier) (char *response, unsigned response_size, const char *verifier_param);
extern void str_query_user (const char *query, char *response, unsigned response_size, ResponseVerifier verifier, const char *verifier_param);

// ResponseVerifier functions
extern _Bool str_verify_y_n (char *response, unsigned response_size, const char *y_or_n);
extern _Bool str_verify_not_empty (char *response, unsigned response_size, const char *unused);

extern const char *str_win_error (void);

static inline char base36(unsigned n) { return (n < 10) ? ('0' + n) : ('a' + (n-10)); };

#endif
