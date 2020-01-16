// ------------------------------------------------------------------
//   encrypt.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"

static char *password = NULL;

void crypt_set_password (char *new_password)
{
    password = new_password;
}

bool crypt_have_password ()
{
    return password != NULL;
}

bool crypt_prompt_for_password()
{
    // to do: consider canceling tty echo while getting password: https://stackoverflow.com/questions/1196418/getting-a-password-in-c-without-using-getpass-3

    // we can only ask for the password if the user hasn't redirected stdin or stdout
    ASSERT0 (isatty (1) && isatty(2), "Error: the file is decrypted, please use --password");

#define MAX_PASSWORD_LEN 100
    password = calloc (MAX_PASSWORD_LEN+1, 1); // allocated once, never freed
    printf("\n\nPassword: ");
    
    fgets (password, MAX_PASSWORD_LEN, stdin);
    password[strlen(password)-1] = 0; // truncate the newline

    ASSERT0 (strlen(password), "Goodbye!"); // exit if user clicked Enter only

    return true;
}

unsigned crypt_padded_len (unsigned len)
{
    return len + (password ? (AES_BLOCKLEN - (len % AES_BLOCKLEN)) % AES_BLOCKLEN : 0);
}

// returns true if encrypted
bool crypt_get_encrypted_len (unsigned *data_encrypted_len /* in/out */, unsigned *padding_len /* optional out */)
{
    if (!password) return false;

    unsigned pad_len = (AES_BLOCKLEN - (*data_encrypted_len % AES_BLOCKLEN)) % AES_BLOCKLEN;
    *data_encrypted_len += pad_len; // we are guaranteed there's room for our padding, bc we kept encryption_padding_reserve bytes for it

    if (padding_len) *padding_len = pad_len;
    
    return true;
}

unsigned crypt_max_padding_len()
{
    return AES_BLOCKLEN-1;
}

// 256 bit AES is a concatenation of 2 MD5 hashes of the password - each one of length 128 bit
// each hash is a hash of the password concatenated with a constant string
// we add data_len to the hash to give it a near-uniqueness for each section
static void crypt_generate_aes_key (VariantBlock *vb,                
                                    uint32_t vb_i, int16_t sec_i, // used to generate an aes key unique to each block
                                    uint8_t *aes_key /* out */)
{
    const char *salt   = "frome";    
    const char *pepper = "vaughan";

    unsigned pw_len = strlen(password);

    if (!buf_is_allocated (&vb->flavored_password)) {
        buf_alloc (vb, &vb->flavored_password, pw_len + 4 + 2 + MAX(strlen(salt), strlen(pepper)), 1, "flavored_password", 0);
        memcpy (vb->flavored_password.data, password, pw_len);
    }
    
    // add some salt to the password, mixed with data_len for uniqueness
    char *flavoured_pw = vb->flavored_password.data;
    char *next = flavoured_pw + pw_len;

    memcpy (next, &vb_i, 4);
    next += 4;

    memcpy (next, &sec_i, 2);
    next += 2;

    Md5Hash salty_hash;
    memcpy (next, salt, strlen(salt));
    md5_do (flavoured_pw, next - flavoured_pw + strlen(salt), &salty_hash);

    // now some pepper
    Md5Hash peppered_hash;
    memcpy (next, pepper, strlen(pepper));
    md5_do (flavoured_pw, next - flavoured_pw + strlen(pepper), &peppered_hash);

    // get hash
    memcpy (aes_key, salty_hash.bytes, sizeof(Md5Hash)); // first half of key
    memcpy (aes_key + sizeof(Md5Hash), peppered_hash.bytes, sizeof(Md5Hash)); // 2nd half of key
}

void crypt_do (VariantBlock *vb, uint8_t *data, unsigned data_len, uint32_t vb_i, int16_t sec_i) // used to generate an aes key unique to each block
{
    //fprintf (stderr, "id:%u vb_i=%d sec_i=%d data_len=%u\n", vb->id, vb_i, sec_i, data_len);

    // generate an AES key just for this one section - combining the pasword with vb_i and sec_i
    uint8_t aes_key[AES_KEYLEN]; 
    crypt_generate_aes_key (vb, vb_i, sec_i, aes_key);

    aes_initialize (vb, aes_key);

    // encrypt in-place
    aes_xcrypt_buffer (vb, data, data_len);
}

void crypt_continue (VariantBlock *vb, uint8_t *data, unsigned data_len)
{
    //fprintf (stderr, "continue: data_len=%u\n", data_len);
    aes_xcrypt_buffer (vb, data, data_len);
}

// pad data to AES_BLOCKLEN boundary
void crypt_pad (uint8_t *data, unsigned data_len, unsigned padding_len)
{
    if (!padding_len) return; // nothing to do

    // use md5 to generate non-trival padding - the hash of the last 100 bytes of data
    Md5Hash hash;
    unsigned src_len = MIN (data_len, 100);
    md5_do (&data[data_len-src_len], src_len, &hash);
    
    memcpy (&data[data_len-padding_len], hash.bytes, padding_len); // luckily the length of MD5 hash and AES block are both 16 bytes - so one hash is sufficient for the padding
}
