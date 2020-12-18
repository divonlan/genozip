// ------------------------------------------------------------------
//   strings.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef STRINGS_INCLUDED
#define STRINGS_INCLUDED

#include "genozip.h"

#define IS_NUCLEOTIDE(c) ((c) == 'A' || (c) == 'T' || (c) == 'C' || (c) == 'G')
#define IS_DIGIT(c)   ((c)>='0' && (c)<='9')
#define IS_CLETTER(c) ((c)>='A' && (c)<='Z')
#define IS_SLETTER(c) ((c)>='a' && (c)<='z')
#define IS_LETTER(c) (IS_CLETTER(c) || IS_SLETTER(c))
#define IS_VALID_URL_CHAR(c) (IS_LETTER(c) || IS_DIGIT(c) || c=='-' || c=='_' || c=='.' || c=='~') // characters valid in a URL
#define FLIP_CASE(c) (IS_CLETTER(c) ? ((c)+32) : (IS_SLETTER(c) ? ((c)-32) : (c))) // flips lower <--> upper case

typedef struct { char s[30]; } StrText;

extern StrText char_to_printable (char c);
extern char *str_tolower (const char *in, char *out /* out allocated by caller - can be the same as in */);

extern StrText str_size (uint64_t size);
extern StrText str_uint_commas (int64_t n);
extern StrText str_int_s (int64_t n);
extern unsigned str_int (int64_t n, char *str /* out */);
extern _Bool str_get_int (const char *str, unsigned str_len, int64_t *value); // note: for a reason beyond me, Docker hub won't compile if its "bool" and not "_Bool"
#define str_is_int(str,str_len) str_get_int ((str), (str_len), NULL)

extern _Bool str_is_in_range (const char *str, uint32_t str_len, char first_c, char last_c);

extern StrText str_pointer (const void *p);

extern const char *type_name (unsigned item, 
                              const char * const *name, // the address in which a pointer to name is found, if item is in range
                              unsigned num_names);

extern void str_print_null_seperated_data (const char *data, unsigned len, _Bool add_newline, _Bool remove_equal_asterisk);

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

#endif
