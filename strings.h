// ------------------------------------------------------------------
//   strings.h
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

#define IS_NUCLEOTIDE(c) ((c) == 'A' || (c) == 'T' || (c) == 'C' || (c) == 'G')
#define IS_DIGIT(c)      ((c)>='0' && (c)<='9')
#define IS_HEXDIGIT(c)   (IS_DIGIT(c) || ((c)>='A' && (c)<='F') || ((c)>='a' && (c)<='f'))
#define IS_HEXDIGITlo(c) (IS_DIGIT(c) || ((c)>='a' && (c)<='f'))
#define IS_HEXDIGITUP(c) (IS_DIGIT(c) || ((c)>='A' && (c)<='F'))
#define IS_CLETTER(c)    ((c)>='A' && (c)<='Z')
#define IS_SLETTER(c)    ((c)>='a' && (c)<='z')
#define IS_WS(c) ((c)=='\t' || (c)=='\r' || (c)=='\n')
#define IS_LETTER(c) (IS_CLETTER(c) || IS_SLETTER(c))
#define IS_NON_WS_PRINTABLE(c) (((c)>=33) && ((c)<=126))
extern bool is_printable[256];
#define IS_PRINTABLE(c) is_printable[(uint8_t)(c)]
#define IS_VALID_URL_CHAR(c) (IS_LETTER(c) || IS_DIGIT(c) || c=='-' || c=='_' || c=='.' || c=='~') // characters valid in a URL
#define FLIP_CASE(c)  (IS_CLETTER(c) ? ((c)+32) : (IS_SLETTER(c) ? ((c)-32) : (c))) // flips lower <--> upper case
#define UPPER_CASE(c) (IS_SLETTER(c) ? ((c)-32) : (c))
#define LOWER_CASE(c) (IS_CLETTER(c) ? ((c)+32) : (c))

#define TF(s) ((s) ? "true" : "false")
#define S(s)  ((s) ? (s) : "(none)")

typedef struct { char s[80]; } StrText;

extern StrText char_to_printable (char c);

extern char *str_print_snip (STRp(in), char *out);

extern char *str_to_printable (STRp(in), char *out);
extern char *str_tolower (const char *in, char *out /* out allocated by caller - can be the same as in */);
extern char *str_toupper (const char *in, char *out);

static bool inline str_issame_ (STRp(str1), STRp(str2)) // true if the same
{
    return (str1_len == str2_len) && !memcmp (str1, str2, str1_len);
}
#define str_issame(str1,str2) str_issame_ (str1, str1##_len, str2, str2##_len)

extern bool str_case_compare (const char *str1, const char *str2, bool *identical); // similar to stricmp that doesn't exist on all platforms

static inline char *str_tolower_(const char *in, char *out, uint32_t len)
{
    for (uint32_t i=0; i < len; i++) out[i] = LOWER_CASE (in[i]); 
    return out;
}

static inline char *str_toupper_(const char *in, char *out, uint32_t len)
{
    for (uint32_t i=0; i < len; i++) out[i] = UPPER_CASE (in[i]);
    return out;
}

static inline char *str_reverse (const char *in, char *out, uint32_t len)
{
    for (uint32_t i=0; i < len; i++) out[len-1-i] = in[i]; 
    return out;
}

extern const char COMPLEM[], UPPER_COMPLEM[];
extern char *str_revcomp_in_out (char *dst_seq, const char *src_seq, uint32_t seq_len);
static inline char *str_revcomp (char *seq, uint32_t seq_len) { return str_revcomp_in_out (seq, seq, seq_len); }

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

extern uint32_t str_int (int64_t n, char *str /* out */);
extern bool str_get_int (STRp(str), int64_t *value); 
extern bool str_get_int_range8  (STRp(str), int64_t min_val, int64_t max_val, uint8_t  *value); // unsigned
extern bool str_get_int_range16 (STRp(str), int64_t min_val, int64_t max_val, uint16_t *value); // unsigned
extern bool str_get_int_range64 (STRp(str), int64_t min_val, int64_t max_val, int64_t  *value); // signed
extern bool str_get_int_range32 (STRp(str), int64_t min_val, int64_t max_val, int32_t  *value); // signed

extern bool str_get_int_dec (STRp(str), uint64_t *value); 
extern bool str_get_int_hex (STRp(str), uint64_t *value); 
extern bool str_get_int_range_allow_hex8  (STRp(str), uint8_t  min_val, uint8_t  max_val, uint8_t  *value);
extern bool str_get_int_range_allow_hex16 (STRp(str), uint16_t min_val, uint16_t max_val, uint16_t *value);
extern bool str_get_int_range_allow_hex32 (STRp(str), uint32_t min_val, uint32_t max_val, uint32_t *value);
extern bool str_get_int_range_allow_hex64 (STRp(str), uint64_t min_val, uint64_t max_val, uint64_t *value);

static bool inline str_is_int(STRp(str))     { return str_get_int (STRa(str), NULL); } // integer - leading zeros not allowed
static bool inline str_is_numeric(STRp(str)) { for (int i=0; i<str_len; i++) if (!IS_DIGIT(str[i])) return false; return true; } // numeric - leading zeros ok
static bool inline str_is_hex(STRp(str))     { for (int i=0; i<str_len; i++) if (!IS_HEXDIGIT(str[i])) return false; return true; } 
static bool inline str_is_hexlo(STRp(str))   { for (int i=0; i<str_len; i++) if (!IS_HEXDIGITlo(str[i])) return false; return true; } 
static bool inline str_is_hexup(STRp(str))   { for (int i=0; i<str_len; i++) if (!IS_HEXDIGITUP(str[i])) return false; return true; } 
static bool inline str_is_no_ws(STRp(str))   { for (int i=0; i<str_len; i++) if (!IS_NON_WS_PRINTABLE(str[i])) return false; return true; } 

extern bool str_is_in_range (STRp(str), char first_c, char last_c);

extern StrText str_pointer (const void *p);
extern StrText str_time (void);

#define FLOAT_FORMAT_LEN 12
extern bool str_get_float (STRp(float_str), double *value, char format[FLOAT_FORMAT_LEN], uint32_t *format_len);

extern bool str_scientific_to_decimal (STRp(float_str), char *modified, uint32_t *modified_len, double *value);

extern uint32_t str_split_do (STRp(str), uint32_t max_items, char sep, const char **items, uint32_t *item_lens, bool exactly, const char *enforce_msg);

#define str_split_enforce(str,str_len,max_items,sep,name,exactly,enforce) \
    uint32_t n_##name##s = (max_items) ? (max_items) : str_count_char ((str), (str_len), (sep)) + 1; /* 0 if str is NULL */ \
    const char *name##s[MAX_(n_##name##s, 1)]; \
    uint32_t name##_lens[MAX_(n_##name##s, 1)]; \
    n_##name##s = str_split_do ((str), (str_len), n_##name##s, (sep), name##s, name##_lens, (exactly), (enforce))

// name      : eg "item", macro defines variables "items" (array of pointers), item_lens (array or uint32_t), n_items (actual number of items)
// max_items : maximum allowed items, or 0 if not known
#define str_split(str,str_len,max_items,sep,name,exactly) str_split_enforce((str),(str_len),(max_items),(sep),name,(exactly),NULL)

extern uint32_t str_split_by_container_do (STRp(str), ConstContainerP con, STRp(con_prefixes), const char **items, uint32_t *item_lens, const char *enforce_msg);

#define str_split_by_container(str,str_len,container,prefix,prefix_len,name) \
    const char *name##s[MAX_(con_nitems(*container), 1)];  \
    uint32_t name##_lens[MAX_(con_nitems(*container), 1)]; \
    uint32_t n_##name##s = str_split_by_container_do ((str), (str_len), (ConstContainerP)(container), (prefix), (prefix_len), name##s, name##_lens, NULL)

extern void str_remove_CR_do (uint32_t n_lines, pSTRp(line));
#define str_remove_CR(name) str_remove_CR_do (n_##name##s, name##s, name##_lens)

extern void str_nul_separate_do (uint32_t n_items, STRps(item));
#define str_nul_separate(name) str_nul_separate_do (n_##name##s, name##s, name##_lens)

extern uint32_t str_remove_whitespace (const char *in, uint32_t in_len, char *out);

extern uint32_t str_split_ints_do (STRp(str), uint32_t max_items, char sep, bool exactly, int64_t *items);
#define str_split_ints(str,str_len,max_items,sep,name,exactly) \
    uint32_t n_##name##s = (max_items) ? (max_items) : str_count_char ((str), (str_len), (sep)) + 1; \
    int64_t name##s[n_##name##s]; \
    n_##name##s = str_split_ints_do ((str), (str_len), n_##name##s, (sep), (exactly), name##s); 

extern uint32_t str_split_floats_do (STRp(str), uint32_t max_items, char sep, bool exactly, double *items);
#define str_split_floats(str,str_len,max_items,sep,name,exactly) \
    uint32_t n_##name##s = (max_items) ? (max_items) : str_count_char ((str), (str_len), (sep)) + 1; \
    double name##s[n_##name##s]; \
    n_##name##s = str_split_floats_do ((str), (str_len), n_##name##s, (sep), (exactly), name##s); 

extern const char *type_name (uint32_t item, 
                              const char * const *name, // the address in which a pointer to name is found, if item is in range
                              uint32_t num_names);

extern void str_print_dict (FILE *fp, STRp(data), bool add_newline, bool remove_equal_asterisk);

extern int str_print_text (const char **text, uint32_t num_lines,
                           const char *wrapped_line_prefix, 
                           const char *newline_separator, 
                           uint32_t line_width /* 0=calcuate optimal */);

typedef bool (*ResponseVerifier) (char *response, uint32_t response_size, const char *verifier_param);
extern void str_query_user (const char *query, char *response, uint32_t response_size, ResponseVerifier verifier, const char *verifier_param);

// ResponseVerifier functions
extern bool str_verify_y_n (char *response, uint32_t response_size, const char *y_or_n);
extern bool str_verify_not_empty (char *response, uint32_t response_size, const char *unused);

extern const char *str_win_error (void);

static inline char base32(uint32_t n) { return (n) < 26 ? ('a' + (n))     // 97->122. 5bits: 1->26    <-- differentiated 5bits, so to work well with ALT_KEY
                                             : (n) < 31 ? ((n)-26 + '[')  // 91->95.  5bits: 27->31
                                             :            '@';          } // 64.      5bits: 0
