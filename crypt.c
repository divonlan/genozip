// ------------------------------------------------------------------
//   encrypt.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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

const char *crypt_get_password()
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
static void crypt_generate_aes_key (VBlock *vb,                
                                    uint32_t vb_i, SectionType sec_type, bool is_header, // used to generate an aes key unique to each block
                                    uint8_t *aes_key /* out */)
{
    const char *salt   = "frome";     
    const char *pepper = "vaughan";   
    static unsigned pw_len=0, salt_len=0, pepper_len=0;

    ASSERT0 (password, "Error in crypt_generate_aes_key - password is NULL");

    if (!pw_len) { // first call
        pw_len     = strlen (password);
        salt_len   = strlen (salt);
        pepper_len = strlen (pepper);
    }

    buf_alloc (vb, &vb->spiced_pw, pw_len + sizeof (uint32_t) + sizeof (uint8_t) + sizeof (uint8_t) + salt_len + pepper_len, 1, "spiced_pw", 0);
    buf_add (&vb->spiced_pw, password, pw_len);

    Md5Hash salty_hash, peppered_hash;
    uint8_t sec_type_byte  = (uint8_t)sec_type; // convert to a byte, as enum and bool might be represented differently by different compilers
    uint8_t is_header_byte = (uint8_t)is_header;

    // add some salt to the password, mixed with vb_i and sec_i for uniqueness
    buf_add (&vb->spiced_pw, &vb_i, sizeof (uint32_t));
    buf_add (&vb->spiced_pw, &sec_type_byte, sizeof (uint8_t));
    buf_add (&vb->spiced_pw, &is_header_byte, sizeof (uint8_t));
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

// we generate a different key for each block by salting the password with vb_i, sec_type and is_header
void crypt_do (VBlock *vb, uint8_t *data, unsigned data_len, 
               uint32_t vb_i, SectionType sec_type, bool is_header)  // used to generate an aes key unique to each block
{
    // generate an AES key just for this one section - combining the pasword with vb_i and sec_i
    uint8_t aes_key[AES_KEYLEN]; 
    crypt_generate_aes_key (vb, vb_i, sec_type, is_header, aes_key);

    //fprintf (stderr, "command:%d id:%d vb_i=%d sec_type=%s%s data_len=%u key=%s\n", command, vb->id, vb_i, st_name (sec_type), is_header ? "(Hdr)": "", data_len, aes_display_key (aes_key));

    aes_initialize (vb, aes_key);

    // encrypt in-place
    //fprintf (stderr, "BFRE: data len=%u: %s\n", data_len, aes_display_data (data, data_len));
    aes_xcrypt_buffer (vb, data, data_len);
    //fprintf (stderr, "AFTR: data len=%u: %s\n", data_len, aes_display_data (data, data_len));
}

void crypt_continue (VBlock *vb, uint8_t *data, unsigned data_len)
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

const char *encryption_name (unsigned encryption_type)
{
    static const char *names[NUM_ENCRYPTION_TYPES] = ENCRYPTION_TYPE_NAMES;
    return type_name (encryption_type, &names[encryption_type], sizeof(names)/sizeof(names[0]));
}
