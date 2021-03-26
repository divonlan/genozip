// ------------------------------------------------------------------
//   crypt.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef CRYPT_INCLUDED
#define CRYPT_INCLUDED

#include "genozip.h"
#include "sections.h"

extern void crypt_set_password (char *new_password);
extern const char *crypt_get_password(void);
extern bool crypt_have_password (void);
extern bool crypt_prompt_for_password(void);
extern uint32_t crypt_padded_len (uint32_t len);
extern bool crypt_get_encrypted_len (uint32_t *data_encrypted_len /* in/out */, uint32_t *padding_len /* out */);
extern void crypt_do (VBlockP vb, uint8_t *data, uint32_t data_len, uint32_t vb_i, SectionType sec_type, bool is_header);
extern void crypt_continue (VBlockP vb, uint8_t *data, uint32_t data_len);
extern void crypt_pad (uint8_t *data, uint32_t data_len, uint32_t padding_len);
extern uint32_t crypt_max_padding_len(void);

extern const char *encryption_name (EncryptionType encryption_type);

#endif