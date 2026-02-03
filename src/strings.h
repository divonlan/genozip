// ------------------------------------------------------------------
//   strings.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

#define IS_DIGIT(c)        IN_RANGX((c), '0', '9')
#define NUM2HEXDIGIT(n)    ((n)<=9 ? '0' + (n) : 'a'+((n)-10))      
#define HEXDIGIT2NUM(c)    (IS_DIGIT(c) ? ((c)-'0') : ((c)-'a'+10)) // converts an lowercase hex digit to a number [0,15]
#define IS_HEXDIGIT(c)     (IS_DIGIT(c) || IN_RANGX((c),'A','F') || IN_RANGX((c),'a','f'))
#define IS_HEXDIGITlo(c)   (IS_DIGIT(c) || IN_RANGX((c),'a','f'))
#define IS_HEXDIGITUP(c)   (IS_DIGIT(c) || IN_RANGX((c),'A','F'))
#define IS_CLETTER(c)      IN_RANGX((c),'A','Z')
#define IS_SLETTER(c)      IN_RANGX((c),'a','z')
#define IS_LETTER(c)       (IS_CLETTER(c) || IS_SLETTER(c))
#define IS_ALPHANUMERIC(c) (IS_LETTER(c)  || IS_DIGIT(c))
#define IS_ACGT(c)         ((c)=='A' || (c)=='T' || (c)=='C' || (c)=='G') // AT often have slightly higher frequency than CG, so testing first
#define IS_ACGTN(c)        (IS_ACGT(c) || (c)=='N')
#define IS_NON_WS(c)       IN_RANGX((c), 33, 126)
#define IS_QUAL_SCORE(c)   IN_RANGX((c), 33, 126)
extern const bool is_printable[256];
#define IS_PRINTABLE(c)    is_printable[(uint8_t)(c)]
extern const bool is_fastq_seq[256];
#define IS_FASTQ_SEQ(c)    is_fastq_seq[(uint8_t)(c)]

#define IS_VALID_URL_CHAR(c) (IS_ALPHANUMERIC(c) || c=='-' || c=='_' || c=='.' || c=='~') // characters valid in a URL
#define FLIP_CASE(c)  (IS_CLETTER(c) ? ((c)+32) : (IS_SLETTER(c) ? ((c)-32) : (c))) // flips lower <--> upper case
#define UPPER_CASE(c) (IS_SLETTER(c) ? ((c)-32) : (c))
#define LOWER_CASE(c) (IS_CLETTER(c) ? ((c)+32) : (c))

#define IS_ASTERISK(str)    str_is_1char(str, '*')
#define IS_ASTERISKi(arr,i) str_is_1chari(arr, (i), '*')
#define IS_PERIOD(str)      str_is_1char(str, '.')
#define IS_SPACE(str)       str_is_1char(str, ' ')
#define IS_EQUAL_SIGN(str)  str_is_1char(str, '=')

#define TF(x) ((x) ? "true" : "false")
#define VX(x) ((x) ? "✓" : "✗")
#define YN(x) ((x)==yes?"Yes" : (x)==no?"No" : "Unknown")
#define S(s)  ((s) ? (s) : "(none)")

extern StrText char_to_printable (char c);

extern uint32_t str_to_printable (STRp(in), char *out, int out_len);
static inline StrTextSuperLong str_to_printable_(STRp(in)) { // for bound-length short texts
    StrTextSuperLong s;
    str_to_printable (STRa(in), s.s, sizeof (s.s));
    return s;
}

extern char *str_tolower (rom in, char *out /* out allocated by caller - can be the same as in */);
extern char *str_toupper (rom in, char *out);

static inline bool str_islower (STRp (str))
{
    for (int i=0; i < str_len; i++)
        if (!IS_SLETTER(str[i])) return false;
    return true;
}

static bool inline str_issame_(STRp(str1), STRp(str2)) // true if the same
{
    return (str1_len == str2_len) && !memcmp (str1, str2, str1_len);
}
#define str_issame(str1,str2) str_issame_ (str1, str1##_len, str2, str2##_len)

static bool inline str_isprefix_(STRp(long_str), STRp(short_str)) // true short_str is a prefix of long_str
{
    return (long_str_len >= short_str_len) && !memcmp (long_str, short_str, short_str_len);
}
#define str_isprefix(long_str,short_str) str_isprefix_ (long_str, long_str##_len, short_str, short_str##_len)

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

#define str_is_1char(str,char) (str##_len==1 && *str==(char))
#define str_is_1chari(arr,i,char) (arr##_lens[i]==1 && *arr##s[i]==(char))

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

extern char *str_revcomp (char *dst_seq, rom src_seq, uint32_t seq_len);
static inline char *str_revcomp_in_place (STRc(seq)) { return str_revcomp (seq, seq, seq_len); }
extern char *str_revcomp_actg (char *dst_seq, rom src_seq, uint32_t seq_len);

// count the number of occurances of a character in a string
static inline uint32_t str_count_char (rom str, uint32_t len, char c)
{
    if (!str) return 0;
    
    uint32_t count=0;
    rom after = str + len;
    while (str < after)
        if (*str++ == c) count++;

    return count;
}

// true if entire string is a single character
static inline bool str_is_monochar (STRp(str))
{
    char mono = str_len ? str[0] : 0;
 
    for (int i=1; i < str_len; i++)
        if (str[i] != mono) return false;
    
    return true;
}

static inline bool str_is_monochar_(STRp(str), char mono)
{
    for (int i=0; i < str_len; i++)
        if (str[i] != mono) return false;
    
    return true;
}

#define is_zero_struct(_struct) str_is_monochar_((rom)&(_struct), sizeof (_struct), 0)

static inline unsigned homopolymer_len (STRp(seq), unsigned start)
{
    char base = seq[start];
    for (unsigned i=start+1; i < seq_len; i++)
        if (seq[i] != base) return i - start;

    return seq_len - start;
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

static inline rom str_a_an (rom s) {
    return (s[0]=='a' || s[0]=='e' || s[0]=='i' || s[0]=='o' || s[0]=='u' ||  
            s[0]=='A' || s[0]=='E' || s[0]=='I' || s[0]=='O' || s[0]=='U') ? "an" : "a";  
}

extern StrText str_size (uint64_t size);
extern StrText str_bases (uint64_t num_bases);
extern StrText str_int_commas (int64_t n);
extern StrText str_uint_commas_limit (uint64_t n, uint64_t limit);
extern StrText str_int_s (int64_t n);
extern StrTextLong str_int_s_(rom label, int64_t n);
#define cond_int(cond, label, n) ((cond) ? str_int_s_((label), (n)).s : "") /* note: n does not evaluate if cond is false! */\

extern StrTextLong str_str_s_(rom label, STRp(str));
#define cond_str(cond, label, str)  ((cond) ? ({ rom str_=(str); str_str_s_((label), str_, strlen (str_)).s; }) : "") /* note: str evaluates once, but only if cond is true */
#define cond_stra(cond, label, str) ((cond) ? str_str_s_((label), str, str_len).s : "") 

extern rom str_to_hex_(bytes data, uint32_t data_lepn, char *hex_str, bool with_dot);
static inline StrText str_to_hex (bytes data, uint32_t data_len) // note: for data_len up to 39 (truncated if longer)
{
    StrText s;
    str_to_hex_(data, MIN_(data_len, sizeof(s)/2 - 1), s.s, false);
    return s;
}

// string length of an integer. #include <math.h> if using this.
static inline unsigned str_int_len (int64_t n_)
{ 
    uint64_t n = ABS(n_); // longest possible number: -9,223,372,036,854,775,808 - 20 characters (without commas)
    return (n_ < 0) +
           (n<10ULL?1 : n<100ULL?2 : n<1000ULL?3 : n<10000ULL?4 : n<100000ULL?5 : n<1000000ULL?6 : n<10000000ULL?7 
          : n<100000000ULL?8 : n<1000000000ULL?9 : n<10000000000ULL?10 : n<100000000000ULL?11 : n<1000000000000ULL?12 
          : n<10000000000000ULL?13 : n<100000000000000ULL?14 : n<1000000000000000ULL?15 : n<10000000000000000ULL?16 
          : n<100000000000000000ULL?17: n<1000000000000000000ULL?18 : 19); 
}

extern uint32_t str_int_ex (int64_t n, char *str /* out */, bool add_nul_terminator);
static inline uint32_t str_int (int64_t n, char *str /* out */) { return str_int_ex (n, str, true); }

static inline uint32_t str_int_fast (int64_t n, char *str /* out */) // faster if many of the numbers are expected to be single-digit
{ 
    if (n <= 9) { *str = '0' + n; return 1; }
    else return str_int_ex (n, str, false); 
}

extern uint32_t str_hex_ex (int64_t n, char *str /* out */, bool uppercase, bool add_nul_terminator);
static inline uint32_t str_hex (int64_t n, char *str /* out */, bool uppercase) { return str_hex_ex (n, str, uppercase, true); }

extern bool str_get_int (STRp(str), int64_t *value); 
extern bool str_get_uint16 (STRp(str), uint16_t *value);
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
extern bool str_is_simple_float (STRp(str), uint32_t *decimals);
static bool inline str_is_hexlo(STRp(str))     { for (int i=0; i<str_len; i++) if (!IS_HEXDIGITlo(str[i]))       return false; return true; } 
static bool inline str_is_hexup(STRp(str))     { for (int i=0; i<str_len; i++) if (!IS_HEXDIGITUP(str[i]))       return false; return true; } 
static bool inline str_is_printable(STRp(str)) { for (int i=0; i<str_len; i++) if (!IS_PRINTABLE(str[i]))        return false; return true; } 
static bool inline str_is_fastq_seq(STRp(str)) { for (int i=0; i<str_len; i++) if (!IS_FASTQ_SEQ(str[i]))        return false; return true; } 
extern bool str_is_utf8 (STRp(str));
static bool inline str_is_no_ws(STRp(str))     { for (int i=0; i<str_len; i++) if (!IS_NON_WS(str[i])) return false; return true; } 
static bool inline str_is_ACGT(STRp(str), uint32_t *bad_i) { for (int i=0; i<str_len; i++) if (!IS_ACGT(str[i])) { if(bad_i) *bad_i = i; return false; } return true; } 
static bool inline str_is_ACGTN(STRp(str))     { for (int i=0; i<str_len; i++) if (!IS_ACGTN(str[i]))            return false; return true; } 

extern bool str_is_in_range (STRp(str), char first_c, char last_c);
static bool inline str_is_upper (STRp(str))    { return str_is_in_range (STRa(str), 'A', 'Z'); }
static bool inline str_is_lower (STRp(str))    { return str_is_in_range (STRa(str), 'a', 'z'); }
static bool inline str_is_numeric(STRp(str))   { return str_is_in_range (STRa(str), '0', '9'); } // numeric - leading zeros ok

extern uint32_t str_pack_bases (uint8_t *packed, STRp(bases), bool revcomp);
extern uint32_t str_unpack_bases (char *bases, bytes packed, uint32_t num_bases);

// textual length of a non-negative integer
extern uint32_t str_get_uint_textual_len (uint64_t n);

extern StrTextLong str_time (void);
extern StrText str_human_time (unsigned secs, bool compact);

#define FLOAT_FORMAT_LEN 12
extern bool str_get_float (STRp(float_str), double *value, char format[FLOAT_FORMAT_LEN], uint32_t *format_len);

extern bool str_scientific_to_decimal (STRp(float_str), qSTRp(modified), double *value);

extern uint32_t str_split_do (STRp(str), uint32_t max_items, char sep, rom *items, uint32_t *item_lens, bool exactly, rom enforce_msg);

#define str_split_enforce(str,str_len,max_items,sep,name,exactly,enforce) \
    uint32_t n_##name##s = (max_items) ? (max_items) : (sep)==0 ? 1 : str_count_char ((str), (str_len), (sep)) + 1; /* 0 if str is NULL */ \
    STR_ARRAY_(name, MAX_(n_##name##s, 1)); \
    n_##name##s = str_split_do ((str), (str_len), n_##name##s, (sep), name##s, name##_lens, (exactly), (enforce))

// name      : eg "item", macro defines variables "items" (array of pointers), item_lens (array or uint32_t), n_items (actual number of items)
// max_items : maximum allowed items, or 0 if not known
#define str_split(str,str_len,max_items,sep,name,exactly) str_split_enforce((str),(str_len),(max_items),(sep),name,(exactly),NULL)

extern uint32_t str_split_by_container_do (STRp(str), ConstContainerP con, STRp(con_prefixes), rom *items, uint32_t *item_lens, rom enforce_msg);

#define str_split_by_container(str,str_len,container,prefix,prefix_len,name,enforce_msg) \
    STR_ARRAY (name, MAX_(con_nitems(*container), 1)) = str_split_by_container_do ((str), (str_len), (ConstContainerP)(container), (prefix), (prefix_len), name##s, name##_lens, (enforce_msg))

extern rom str_split_by_tab_do (STRp(str), uint32_t *n_flds, rom *flds, uint32_t *fld_lens, bool *has_13, bool exactly, bool ignore_excess, bool enforce);
#define str_split_by_tab(str,max_len,max_flds,has_13,exactly,ignore_excess,enforce) \
    STR_ARRAY (fld, (max_flds)) = (max_flds);   \
    str = str_split_by_tab_do ((str), (max_len), &n_flds, flds, fld_lens, (has_13), (exactly), (ignore_excess), (enforce))
#define STRfld(i) STRi(fld,(i))

extern uint32_t str_split_by_lines_do (STRp(str), uint32_t max_lines, rom *lines, uint32_t *line_lens);
#define str_split_by_lines(str,str_len,max_lines) \
    STR_ARRAY (line, (max_lines)) = str_split_by_lines_do ((str), (str_len), max_lines, lines, line_lens)

extern uint32_t str_split_ints_do (STRp(str), uint32_t max_items, char sep, bool exactly, int base, int64_t *items);

#define str_split_ints(str,str_len,max_items,sep,name,exactly) \
    uint32_t n_##name##s = (max_items) ? (max_items) : str_count_char ((str), (str_len), (sep)) + 1; \
    int64_t name##s[n_##name##s]; \
    n_##name##s = str_split_ints_do ((str), (str_len), n_##name##s, (sep), (exactly), 10, name##s); 

#define str_split_hexs(str,str_len,max_items,sep,name,exactly) \
    uint32_t n_##name##s = (max_items) ? (max_items) : str_count_char ((str), (str_len), (sep)) + 1; \
    int64_t name##s[n_##name##s]; \
    n_##name##s = str_split_ints_do ((str), (str_len), n_##name##s, (sep), (exactly), 16, name##s); 

extern uint32_t str_split_floats_do (STRp(str), uint32_t max_items, char sep, bool exactly, char this_char_is_NAN, double *items);
#define str_split_floats(str,str_len,max_items,sep,name,exactly,this_char_is_NAN) \
    uint32_t n_##name##s = (max_items) ? (max_items) : str_count_char ((str), (str_len), (sep)) + 1; \
    double name##s[n_##name##s]; \
    n_##name##s = str_split_floats_do ((str), (str_len), n_##name##s, (sep), (exactly), (this_char_is_NAN), name##s); 

extern bool str_item_i (STRp(str), char sep, uint32_t requested_item_i, pSTRp(item));
extern bool str_item_i_int (STRp(str), char sep, uint32_t requested_item_i, int64_t *item);
extern bool str_item_i_float (STRp(str), char sep, uint32_t requested_item_i, double *item);

extern void str_remove_CR_do (uint32_t n_lines, pSTRp(line));
#define str_remove_CR(name) str_remove_CR_do (n_##name##s, name##s, name##_lens)

extern void str_nul_separate_do (STRps(item));
#define str_nul_separate(name) str_nul_separate_do (n_##name##s, name##s, name##_lens)

extern uint32_t str_remove_whitespace (STRp(in), bool also_uppercase, char *out);
extern void str_trim (qSTRp(str));

extern rom type_name (uint32_t item, 
                      rom  const *name, // the address in which a pointer to name is found, if item is in range
                      uint32_t num_names);

extern int str_print_text (rom *text, uint32_t num_lines, rom wrapped_line_prefix, rom newline_separator,
                           rom added_header, uint32_t line_width /* 0=calcuate optimal */);

typedef bool (*ResponseVerifier) (STRc(response), rom verifier_param);
extern void str_query_user (rom query, STRc(response), bool allow_empty, ResponseVerifier verifier, rom verifier_param);

typedef enum { QDEF_NONE, QDEF_NO, QDEF_YES } DefAnswerType;
extern bool str_query_user_yn (rom query, DefAnswerType def_answer);

extern char *memchr2 (const void *p, char ch1, char ch2, uint32_t count);

// implementing memrchr, as it doesn't exist in Windows libc (msvcrt.dll) or Darwin (at least clang)
#if defined _WIN32 || defined __APPLE__
extern void *memrchr (const void *s, int c, size_t n);
#endif

#ifdef __APPLE__ // is mempcpy available on gcc of darwin? if so, this is ifdef should be just for clang
static inline void *mempcpy(void *restrict dst, const void *restrict src, size_t n)
{
    memcpy (dst, src, n);
    return dst + n;
}
#endif

static inline char *strpcpy(char *restrict dst, const void *restrict src)
{
    int len = strlen (src);
    return mempcpy (dst, src, len);
}

#ifdef __MINGW64__
extern void *memmem (const void *haystack, size_t haystack_len, const void *needle, size_t needle_len);
#endif

extern rom str_win_error_(uint32_t error);
extern rom str_win_error (void);

static inline char base32(uint32_t n) { return (n) < 26 ? ('a' + (n))     // 97->122. 5bits: 1->26    <-- differentiated 5bits, so to work well with ALT_KEY
                                             : (n) < 31 ? ((n)-26 + '[')  // 91->95.  5bits: 27->31
                                             :            '@';          } // 64.      5bits: 0

extern const uint64_t p10[]; // powers of 10

extern uint64_t crc64 (uint64_t crc, bytes data, uint64_t data_len); // implementation is in crc64.c
