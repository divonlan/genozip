// ------------------------------------------------------------------
//   crypt.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef CRYPT_INCLUDED
#define CRYPT_INCLUDED

#include "genozip.h"
#include "sections.h"

extern void crypt_set_password (char *new_password);
extern const char *crypt_get_password(void);
extern bool crypt_have_password (void);
extern bool crypt_prompt_for_password(void);
extern unsigned crypt_padded_len (unsigned len);
extern bool crypt_get_encrypted_len (unsigned *data_encrypted_len /* in/out */, unsigned *padding_len /* out */);
extern void crypt_do (VariantBlockP vb, uint8_t *data, unsigned data_len, uint32_t vb_i, SectionType sec_type, bool is_header);
extern void crypt_continue (VariantBlockP vb, uint8_t *data, unsigned data_len);
extern void crypt_pad (uint8_t *data, unsigned data_len, unsigned padding_len);
extern unsigned crypt_max_padding_len(void);

// genozip v1 compatibility
extern void v1_crypt_do (VariantBlockP vb, uint8_t *data, unsigned data_len, uint32_t vb_i, int16_t sec_i); // used to generate an aes key unique to each block

#endif