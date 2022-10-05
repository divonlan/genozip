// ------------------------------------------------------------------
//   strings.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

#define IS_NUCLEOTIDE(c) ((c) == 'A' || (c) == 'T' || (c) == 'C' || (c) == 'G')
#define IS_DIGIT(c)      ((c)>='0' && (c)<='9')
#define NUM2HEXDIGIT(n)  ((n)<=9 ? '0' + (n) : 'A'+((n)-10))  // converts a number [0,15] to hex digit character
#define HEXDIGIT2NUM(c)  (IS_DIGIT(c) ? ((c)-'0') : ((c)-'A') // converts an uppercase hex digit to a number [0,15]
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

#define IS_ASTERISK(str) (str##_len==1 && *str=='*')
#define IS_EQUAL_SIGN(str) (str##_len==1 && *str=='=')

#define TF(s) ((s) ? "true" : "false")
#define S(s)  ((s) ? (s) : "(none)")

extern StrText char_to_printable (char c);

extern char *str_print_snip (STRp(in), char *out);

extern char *str_to_printable (STRp(in), char *out);
extern char *str_tolower (rom in, char *out /* out allocated by caller - can be the same as in */);
extern char *str_toupper (rom in, char *out);

static bool inline str_issame_ (STRp(str1), STRp(str2)) // true if the same
{
    return (str1_len == str2_len) && !memcmp (str1, str2, str1_len);
}
#define str_issame(str1,str2) str_issame_ (str1, str1##_len, str2, str2##_len)

// same as str_issame, but internally it searches in reverse - good for strings that are known to differ at their end
static bool inline str_issameR_ (STRp(str1), STRp(str2)) // true if the same
{
    if (str1_len != str2_len) return false;

    for (int32_t i=str1_len-1; i >= 0; i--)
        if (str1[i] != str2[i]) return false;

    return true;
}
#define str_issameR(str1,str2) str_issameR_ (str1, str1##_len, str2, str2##_len)

static bool inline str_issame_rev_(STRp(str1), STRp(str2)) // true if the same
{
    if (str1_len != str2_len) return false;
    for (uint32_t i=0; i < str1_len; i++)
        if (str1[i] != str2[str1_len-i-1]) return false;
    return true;
}
#define str_issame_rev(str1,str2) str_issame_rev_ (str1, str1##_len, str2, str2##_len)

// including lowercase, IUPACs
extern const char COMPLEM[], UPPER_COMPLEM[]; 
// only A,C,G,T
static inline char complem(char b) { switch(b) { case'A':return'T' ; case'C':return'G' ; case'G':return'C' ; case'T':return'A' ; default:return b; }} 

static bool inline str_issame_revcomp_(STRp(str1), STRp(str2)) // true if the same
{
    if (str1_len != str2_len) return false;
    for (uint32_t i=0; i < str1_len; i++)
        if (str1[i] != COMPLEM[(uint8_t)str2[str1_len-i-1]]) return false;
    return true;
}
#define str_issame_revcomp(str1,str2) str_issame_revcomp_ (str1, str1##_len, str2, str2##_len)


extern bool str_case_compare (rom str1, rom str2, bool *identical); // similar to stricmp that doesn't exist on all platforms

static inline char *str_tolower_(rom in, char *out, uint32_t len)
{
    for (uint32_t i=0; i < len; i++) out[i] = LOWER_CASE (in[i]); 
    return out;
}

static inline char *str_toupper_(rom in, char *out, uint32_t len)
{
    for (uint32_t i=0; i < len; i++) out[i] = UPPER_CASE (in[i]);
    return out;
}

static inline char *str_reverse (char *dst, rom src, uint32_t len)
{
    for (uint32_t i=0; i < len; i++) dst[len-1-i] = src[i]; 
    return dst;
}

static inline char *str_reverse_in_place (char *str, uint32_t len)
{
    for (uint32_t i=0; i < len/2; i++) SWAP(str[i], str[len-i-1]);
    return str;
}

static inline char *str_replace_letter (STRc(str), char before, char after)
{
    for (uint32_t i=0; i < str_len; i++)
        if (str[i] == before) str[i] = after;

    return str;
}

extern char *str_revcomp_in_out (char *dst_seq, rom src_seq, uint32_t seq_len);
static inline char *str_revcomp (STRc(seq)) { return str_revcomp_in_out (seq, seq, seq_len); }
extern char *str_revcomp_actg (char *dst_seq, rom src_seq, uint32_t seq_len);

// count the number of occurances of a character in a string
static inline uint64_t str_count_char (rom str, uint64_t len, char c)
{
    if (!str) return 0;
    
    uint64_t count=0;
    for (uint64_t i=0; i < len; i++)
        if (str[i] == c) count++;

    return count;
}

// count the number of consecutive occurances of a character
static inline uint64_t str_count_consecutive_char (rom str, uint64_t len, char c)
{
    if (!str) return 0;
    
    uint64_t i=0;
    for (i=0; i < len; i++)
        if (str[i] != c) break;

    return i;
}

extern StrText str_size (uint64_t size);
extern StrText str_bases (uint64_t num_bases);
extern StrText str_int_commas (int64_t n);
extern StrText str_uint_commas_limit (uint64_t n, uint64_t limit);
extern StrText str_int_s (int64_t n);

extern rom str_to_hex (bytes data, uint32_t data_len, char *hex_str, bool with_dot);
extern StrText str_hex10 (bytes data, uint32_t data_len); // up to 10 bytes in hex (21 chars inc. \0)

// string length of an integer. #include <math.h> if using this.
static inline unsigned str_int_len (uint32_t n)
    { return n<10?1 : n<100?2 : n<1000?3 : n<10000?4 : n<100000?5 : n<1000000?6 : n<10000000?7 : n<100000000?8 : n<1000000000?9 : 10; }

extern uint32_t str_int_ex (int64_t n, char *str /* out */, bool add_nul_terminator);
static inline uint32_t str_int (int64_t n, char *str /* out */) { return str_int_ex (n, str, true); }

static inline uint32_t str_int_fast (int64_t n, char *str /* out */) // faster if many of the numbers are expected to be single-digit
{ 
    if (n <= 9) { *str = '0' + n; return 1; }
    else return str_int_ex (n, str, true); 
}

extern uint32_t str_hex_ex (int64_t n, char *str /* out */, bool uppercase, bool add_nul_terminator);
static inline uint32_t str_hex (int64_t n, char *str /* out */, bool uppercase) { return str_hex_ex (n, str, uppercase, true); }

extern bool str_get_int (STRp(str), int64_t *value); 
extern bool str_get_uint32 (STRp(str), uint32_t *value); 
extern bool str_get_int_range8  (STRp(str), int64_t min_val, int64_t max_val, uint8_t  *value); // unsigned
extern bool str_get_int_range16 (STRp(str), int64_t min_val, int64_t max_val, uint16_t *value); // unsigned
extern bool str_get_int_range64 (STRp(str), int64_t min_val, int64_t max_val, int64_t  *value); // signed
extern bool str_get_int_range32 (STRp(str), int64_t min_val, int64_t max_val, int32_t  *value); // signed

extern bool str_get_int_dec (STRp(str), uint64_t *value); 
extern bool str_get_int_hex (STRp(str), bool allow_hex, bool allow_HEX, uint64_t *value); 
extern bool str_get_int_range_allow_hex8  (STRp(str), uint8_t  min_val, uint8_t  max_val, uint8_t  *value);
extern bool str_get_int_range_allow_hex16 (STRp(str), uint16_t min_val, uint16_t max_val, uint16_t *value);
extern bool str_get_int_range_allow_hex32 (STRp(str), uint32_t min_val, uint32_t max_val, uint32_t *value);
extern bool str_get_int_range_allow_hex64 (STRp(str), uint64_t min_val, uint64_t max_val, uint64_t *value);

static bool inline str_is_int(STRp(str))       { return str_get_int (STRa(str), NULL); } // integer - leading zeros not allowed
static bool inline str_is_numeric(STRp(str))   { for (int i=0; i<str_len; i++) if (!IS_DIGIT(str[i])) return false; return true; } // numeric - leading zeros ok
static bool inline str_is_hex(STRp(str))       { for (int i=0; i<str_len; i++) if (!IS_HEXDIGIT(str[i])) return false; return true; } 
static bool inline str_is_hexlo(STRp(str))     { for (int i=0; i<str_len; i++) if (!IS_HEXDIGITlo(str[i])) return false; return true; } 
static bool inline str_is_hexup(STRp(str))     { for (int i=0; i<str_len; i++) if (!IS_HEXDIGITUP(str[i])) return false; return true; } 
static bool inline str_is_no_ws(STRp(str))     { for (int i=0; i<str_len; i++) if (!IS_NON_WS_PRINTABLE(str[i])) return false; return true; } 
static bool inline str_is_only_ACGT(STRp(str), uint32_t *bad_i) { for (int i=0; i<str_len; i++) if (str[i]!='A' && str[i]!='C' && str[i]!='G' && str[i]!='T') { if(bad_i) *bad_i = i; return false; } return true; } 

extern bool str_is_in_range (STRp(str), char first_c, char last_c);

// textual length of a non-negative integer
extern uint32_t str_get_uint_textual_len (uint64_t n);

extern StrText str_time (void);
extern void str_human_time (unsigned secs, bool compact, char *str /* out */);

#define FLOAT_FORMAT_LEN 12
extern bool str_get_float (STRp(float_str), double *value, char format[FLOAT_FORMAT_LEN], uint32_t *format_len);

extern bool str_scientific_to_decimal (STRp(float_str), STRe(modified), double *value);

extern uint32_t str_split_do (STRp(str), uint32_t max_items, char sep, rom *items, uint32_t *item_lens, bool exactly, rom enforce_msg);

#define str_split_enforce(str,str_len,max_items,sep,name,exactly,enforce) \
    uint32_t n_##name##s = (max_items) ? (max_items) : (sep)==0 ? 1 : str_count_char ((str), (str_len), (sep)) + 1; /* 0 if str is NULL */ \
    rom name##s[MAX_(n_##name##s, 1)]; \
    uint32_t name##_lens[MAX_(n_##name##s, 1)]; \
    n_##name##s = str_split_do ((str), (str_len), n_##name##s, (sep), name##s, name##_lens, (exactly), (enforce))

// name      : eg "item", macro defines variables "items" (array of pointers), item_lens (array or uint32_t), n_items (actual number of items)
// max_items : maximum allowed items, or 0 if not known
#define str_split(str,str_len,max_items,sep,name,exactly) str_split_enforce((str),(str_len),(max_items),(sep),name,(exactly),NULL)

extern uint32_t str_split_by_container_do (STRp(str), ConstContainerP con, STRp(con_prefixes), rom *items, uint32_t *item_lens, rom enforce_msg);

#define str_split_by_container(str,str_len,container,prefix,prefix_len,name) \
    rom name##s[MAX_(con_nitems(*container), 1)];  \
    uint32_t name##_lens[MAX_(con_nitems(*container), 1)]; \
    uint32_t n_##name##s = str_split_by_container_do ((str), (str_len), (ConstContainerP)(container), (prefix), (prefix_len), name##s, name##_lens, NULL)

extern rom str_split_by_tab_do (STRp(str), uint32_t *n_items, rom *items, uint32_t *item_lens, bool *has_13, bool enforce);

#define str_split_by_tab(str,max_len,max_items,has_13,enforce) \
    rom flds[(max_items)];   \
    uint32_t fld_lens[(max_items)];  \
    uint32_t n_flds = (max_items);   \
    str = str_split_by_tab_do ((str), (max_len), &n_flds, flds, fld_lens, (has_13), (enforce))
#define STRfld(i) STRi(fld,(i))

extern void str_remove_CR_do (uint32_t n_lines, pSTRp(line));
#define str_remove_CR(name) str_remove_CR_do (n_##name##s, name##s, name##_lens)

extern void str_nul_separate_do (STRps(item));
#define str_nul_separate(name) str_nul_separate_do (n_##name##s, name##s, name##_lens)

extern uint32_t str_remove_whitespace (rom in, uint32_t in_len, char *out);

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

extern rom type_name (uint32_t item, 
                      rom  const *name, // the address in which a pointer to name is found, if item is in range
                      uint32_t num_names);

extern void str_print_dict (FILE *fp, STRp(data), bool add_newline, bool remove_equal_asterisk);

extern int str_print_text (rom *text, uint32_t num_lines,
                           rom wrapped_line_prefix, 
                           rom newline_separator, 
                           uint32_t line_width /* 0=calcuate optimal */);

typedef bool (*ResponseVerifier) (STRc(response), rom verifier_param);
extern void str_query_user (rom query, STRc(response), ResponseVerifier verifier, rom verifier_param);

typedef enum { QDEF_NONE, QDEF_NO, QDEF_YES } DefAnswerType;
extern bool str_query_user_yn (rom query, DefAnswerType def_answer);

// ResponseVerifier functions
extern bool str_verify_not_empty (STRc(response), rom unused);

extern rom my_memrchr (rom str, char c, uint32_t str_len);

extern rom str_win_error (void);

static inline char base32(uint32_t n) { return (n) < 26 ? ('a' + (n))     // 97->122. 5bits: 1->26    <-- differentiated 5bits, so to work well with ALT_KEY
                                             : (n) < 31 ? ((n)-26 + '[')  // 91->95.  5bits: 27->31
                                             :            '@';          } // 64.      5bits: 0

extern uint64_t p10[];