// ------------------------------------------------------------------
//   encrypt.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "aes.h"
#include "crypt.h"
#include "vb.h"
#include "md5.h"

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
    ASSERT0 (isatty (1) && isatty(2), "Error: this file is encrypted, please use --password");

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
    static unsigned pw_len=0, salt_len=0, pepper_len=0;

    if (!pw_len) { // first call
        pw_len     = strlen (password);
        salt_len   = strlen (salt);
        pepper_len = strlen (pepper);
    }

    buf_alloc (vb, &vb->spiced_pw, pw_len + sizeof (uint32_t) + sizeof (int16_t) + salt_len + pepper_len, 1, "spiced_pw", 0);
    buf_add (&vb->spiced_pw, password, pw_len);

    Md5Hash salty_hash, peppered_hash;
    
    // add some salt to the password, mixed with vb_i and sec_i for uniqueness
    buf_add (&vb->spiced_pw, &vb_i, sizeof (uint32_t));
    buf_add (&vb->spiced_pw, &sec_i, sizeof (int16_t));
    buf_add (&vb->spiced_pw, salt, salt_len);
    md5_do (vb->spiced_pw.data, vb->spiced_pw.len, &salty_hash);

    // add some pepper
    buf_add (&vb->spiced_pw, pepper, pepper_len);
    md5_do (vb->spiced_pw.data, vb->spiced_pw.len, &peppered_hash);

    // get hash
    memcpy (aes_key, salty_hash.bytes, sizeof(Md5Hash)); // first half of key
    memcpy (aes_key + sizeof(Md5Hash), peppered_hash.bytes, sizeof(Md5Hash)); // 2nd half of key

    buf_free (&vb->spiced_pw);
}

// we generate a different key for each block by salting the password with vb_i and sec_i
// for sec_i we use (-1-section_i) for the section header and section_i for the section body
// the VCF header section: vb_i=0 (for all components) and sec_i=0 (i.e: 0 for the body, (-1 - 0)=-1 for header)
// the Variant Data section: vb_i={global consecutive number starting at 1}, sec_i=0 (body=0, header=-1)
// Other sections: vb_i same as Variant Data, sec_i consecutive running starting at 1 
void crypt_do (VariantBlock *vb, uint8_t *data, unsigned data_len, uint32_t vb_i, int16_t sec_i) // used to generate an aes key unique to each block
{
    //fprintf (stderr, "id:%u vb_i=%d sec_i=%d data_len=%u\n", vb->id, vb_i, sec_i, data_len);

    // generate an AES key just for this one section - combining the pasword with vb_i and sec_i
    uint8_t aes_key[AES_KEYLEN]; 
    crypt_generate_aes_key (vb, vb_i, sec_i, aes_key);

    //printf ("vb_i=%u sec_i=%d key= %s\n", vb_i, sec_i, aes_display_key (aes_key)); // DEBUG

    aes_initialize (vb, aes_key);

    // encrypt in-place
    //printf ("BFRE: data len=%u: %s\n", data_len, aes_display_data (data, data_len));
    aes_xcrypt_buffer (vb, data, data_len);
    //printf ("AFTR: data len=%u: %s\n", data_len, aes_display_data (data, data_len));
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
