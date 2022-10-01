// ------------------------------------------------------------------
//   encrypt.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "aes.h"
#include "crypt.h"
#include "vblock.h"
#include "md5.h"
#include "strings.h"

static char *password = NULL;

void crypt_set_password (char *new_password)
{
    password = new_password;
}

rom crypt_get_password()
{
    return password;
}

bool crypt_have_password ()
{
    return password != NULL;
}

bool crypt_prompt_for_password()
{
    // to do: consider canceling tty echo while getting password: https://stackoverflow.com/questions/1196418/getting-a-password-in-c-without-using-getpass-3

    // we can only ask for the password if the user hasn't redirected stdin or stdout
    ASSINP0 (isatty (1) && isatty(2), "Error: this file is encrypted, please use --password");

#define MAX_PASSWORD_LEN 100
    password = CALLOC (MAX_PASSWORD_LEN+1); // allocated once, never freed
    printf("\n\nPassword: ");
    
    (void)!fgets (password, MAX_PASSWORD_LEN, stdin); // (void)! to avoid compiler "warning: ignoring return value"
    password[strlen(password)-1] = 0; // truncate the newline

    if (!strlen (password)) { // exit if user clicked Enter only
        printf ("Goodbye!\n"); 
        exit (EXIT_OK); // not an error
    }

    return true;
}

uint32_t crypt_padded_len (uint32_t len)
{
    return len + (password ? (AES_BLOCKLEN - (len % AES_BLOCKLEN)) % AES_BLOCKLEN : 0);
}

// returns true if encrypted
bool crypt_get_encrypted_len (uint32_t *data_encrypted_len /* in/out */, uint32_t *padding_len /* optional out */)
{
    if (!password) return false;

    uint32_t pad_len = (AES_BLOCKLEN - (*data_encrypted_len % AES_BLOCKLEN)) % AES_BLOCKLEN;
    *data_encrypted_len += pad_len; // we are guaranteed there's room for our padding, bc we kept encryption_padding_reserve bytes for it

    if (padding_len) *padding_len = pad_len;
    
    return true;
}

uint32_t crypt_max_padding_len()
{
    return AES_BLOCKLEN-1;
}

// 256 bit AES is a concatenation of 2 MD5 hashes of the password - each one of length 128 bit
// each hash is a hash of the password bound with a constant string
// we add data_len to the hash to give it a near-uniqueness for each section
static void crypt_generate_aes_key (VBlockP vb,                
                                    VBIType vb_i, SectionType sec_type, bool is_header, // used to generate an aes key unique to each block
                                    uint8_t *aes_key /* out */)
{
    rom salt   = "frome";     
    rom pepper = "vaughan";   
    static uint32_t pw_len=0, salt_len=0, pepper_len=0;

    ASSERTNOTNULL (password);

    if (!pw_len) { // first call
        pw_len     = strlen (password);
        salt_len   = strlen (salt);
        pepper_len = strlen (pepper);
    }

    buf_alloc (vb, &vb->spiced_pw, 0, pw_len + sizeof (uint32_t) + sizeof (uint8_t) + sizeof (uint8_t) + salt_len + pepper_len, char, 1, "spiced_pw");
    buf_add (&vb->spiced_pw, password, pw_len);

    uint8_t sec_type_byte  = (uint8_t)sec_type; // convert to a byte, as enum and bool might be represented differently by different compilers
    uint8_t is_header_byte = (uint8_t)is_header;

    // add some salt to the password, mixed with vb_i and sec_i for uniqueness
    buf_add (&vb->spiced_pw, &vb_i, sizeof (uint32_t));
    buf_add (&vb->spiced_pw, &sec_type_byte, sizeof (uint8_t));
    buf_add (&vb->spiced_pw, &is_header_byte, sizeof (uint8_t));
    buf_add (&vb->spiced_pw, salt, salt_len);
    Digest salty_hash = md5_do (STRb(vb->spiced_pw));

    // add some pepper
    buf_add (&vb->spiced_pw, pepper, pepper_len);
    Digest peppered_hash = md5_do (STRb(vb->spiced_pw));

    // get hash
    memcpy (aes_key, salty_hash.bytes, sizeof(Digest)); // first half of key
    memcpy (aes_key + sizeof(Digest), peppered_hash.bytes, sizeof(Digest)); // 2nd half of key

    buf_free (vb->spiced_pw);
}

// we generate a different key for each block by salting the password with vb_i, sec_type and is_header
void crypt_do (VBlockP vb, uint8_t *data, uint32_t data_len, 
               VBIType vb_i, SectionType sec_type, bool is_header)  // used to generate an aes key unique to each block
{
    // generate an AES key just for this one section - combining the pasword with vb_i and sec_i
    uint8_t aes_key[AES_KEYLEN]; 
    crypt_generate_aes_key (vb, vb_i, sec_type, is_header, aes_key);

    aes_initialize (vb, aes_key);

    // encrypt in-place
    aes_xcrypt_buffer (vb, STRa(data));
}

void crypt_continue (VBlockP vb, uint8_t *data, uint32_t data_len)
{
    aes_xcrypt_buffer (vb, STRa(data));
}

// pad data to AES_BLOCKLEN boundary
void crypt_pad (uint8_t *data, uint32_t data_len, uint32_t padding_len)
{
    if (!padding_len) return; // nothing to do

    // use md5 to generate non-trival padding - the hash of the last 100 bytes of data
    uint32_t src_len = MIN_(data_len, 100);
    Digest hash = md5_do (&data[data_len - src_len], src_len);
    
    memcpy (&data[data_len-padding_len], hash.bytes, padding_len); // luckily the length of MD5 hash and AES block are both 16 bytes - so one hash is sufficient for the padding
}

rom encryption_name (EncryptionType encryption_type)
{
    static rom names[NUM_ENCRYPTION_TYPES] = ENC_NAMES;
    return type_name (encryption_type, &names[encryption_type], ARRAY_LEN(names));
}
