// ------------------------------------------------------------------
//   strings.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <stdint.h>

#define IS_DIGIT(c)   ((c)>='0' && (c)<='9')
#define IS_CLETTER(c) ((c)>='A' && (c)<='Z')
#define IS_SLETTER(c) ((c)>='a' && (c)<='z')
#define IS_LETTER(c) (IS_CLETTER(c) || IS_SLETTER(c))

extern void str_to_lowercase (char *s);

extern char *str_size (int64_t size, char *str /* out */);
extern char *str_uint_commas (int64_t n, char *str /* out */);
extern char *str_uint (int64_t n, char *str /* out */, unsigned *len);

#define POINTER_STR_LEN 19
extern char *str_pointer (const void *p, char *str /* POINTER_STR_LEN bytes allocated by caller*/);

