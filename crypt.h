// ------------------------------------------------------------------
//   crypt.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef CRYPT_INCLUDED
#define CRYPT_INCLUDED

#include "genozip.h"
#include "sections.h"

// encryption types - these values are part of the genozip file format and cannot be easily changed
#define NUM_ENCRYPTION_TYPES   2
#define ENCRYPTION_TYPE_NONE   0
#define ENCRYPTION_TYPE_AES256 1
#define ENCRYPTION_TYPE_NAMES { "No encryption", "AES 256 bit" }

extern void crypt_set_password (char *new_password);
extern const char *crypt_get_password(void);
extern bool crypt_have_password (void);
extern bool crypt_prompt_for_password(void);
extern unsigned crypt_padded_len (unsigned len);
extern bool crypt_get_encrypted_len (unsigned *data_encrypted_len /* in/out */, unsigned *padding_len /* out */);
extern void crypt_do (VBlockP vb, uint8_t *data, unsigned data_len, uint32_t vb_i, SectionType sec_type, bool is_header);
extern void crypt_continue (VBlockP vb, uint8_t *data, unsigned data_len);
extern void crypt_pad (uint8_t *data, unsigned data_len, unsigned padding_len);
extern unsigned crypt_max_padding_len(void);

extern const char *encryption_name (unsigned encryption_type);

#endif