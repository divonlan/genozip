// ------------------------------------------------------------------
//   crypt.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "sections.h"

extern void crypt_set_password (char *new_password);
extern rom crypt_get_password(void);
extern bool crypt_have_password (void);
extern bool crypt_prompt_for_password(void);
extern uint32_t crypt_padded_len (uint32_t len);
extern bool crypt_get_encrypted_len (uint32_t *data_encrypted_len /* in/out */, uint32_t *padding_len /* out */);
extern void crypt_do (VBlockP vb, uint8_t *data, uint32_t data_len, VBIType vb_i, SectionType sec_type, bool is_header);
extern void crypt_continue (VBlockP vb, uint8_t *data, uint32_t data_len);
extern void crypt_pad (uint8_t *data, uint32_t data_len, uint32_t padding_len);
extern uint32_t crypt_max_padding_len(void);

extern rom encryption_name (EncryptionType encryption_type);
