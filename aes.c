// ------------------------------------------------------------------
//   aes.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt
//   Copyright claimed on additions and modifications vs public domain.
//
// This is an implementation of AES based on https://github.com/kokke/tiny-AES-c
// which has been heavily modified. The unmodified terms if the license of tiny-AES-c are as follows:
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
// This is free and unencumbered software released into the public domain.
//
// Anyone is free to copy, modify, publish, use, compile, sell, or
// distribute this software, either in source code form or as a compiled
// binary, for any purpose, commercial or non-commercial, and by any
// means.
//
// In jurisdictions that recognize copyright laws, the author or authors
// of this software dedicate any and all copyright interest in the
// software to the public domain. We make this dedication for the benefit
// of the public at large and to the detriment of our heirs and
// successors. We intend this dedication to be an overt act of
// relinquishment in perpetuity of all present and future rights to this
// software under copyright law.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
//
// For more information, please refer to <http://unlicense.org/>
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#include "genozip.h" 
#include "vblock.h"
#include "aes.h"

#define Nk (AES_KEYLEN/4)
#define Nr 14 // value is for 256 bit key AES

// state - array holding the intermediate results during decryption. can be viewed as a matrix or vector
#define Nb 4 
typedef union {
    uint8_t m[Nb][Nb];
    uint8_t v[Nb*Nb];
} AesStateType;

static const uint8_t sbox[256] = {
    //0     1    2      3     4    5     6     7      8    9     A      B    C     D     E     F
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 };

// The round constant word array, Rcon[i], contains the values given by 
// x to the power (i-1) being powers of x (x is denoted as {02}) in the field GF(2^8)
static const uint8_t Rcon[11] = { 0x8d, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36 };

// This function shifts the 4 bytes in a word to the left once. [a0,a1,a2,a3] becomes [a1,a2,a3,a0]
static void inline aes_rotate_word (uint8_t *w)
{
    const uint8_t u8tmp = w[0];
    w[0] = w[1];
    w[1] = w[2];
    w[2] = w[3];
    w[3] = u8tmp;
}

// this function that takes a four-byte input word and applies the S-box to each of the four bytes to produce an output word.
static void inline aes_substitute_word (uint8_t *w)
{
    w[0] = sbox[w[0]];
    w[1] = sbox[w[1]];
    w[2] = sbox[w[2]];
    w[3] = sbox[w[3]];
}

static void inline aes_assign_word (uint8_t *dst, const uint8_t *src)
{
    *(uint32_t *)dst = *(uint32_t *)src;
}

static void inline aes_xor_word (uint8_t *dst, const uint8_t *w1, const uint8_t *w2)
{
    *(uint32_t *)dst = *(uint32_t *)w1 ^ *(uint32_t *)w2;
}

// This function produces Nb(Nr+1) round keys. The round keys are used in each round to decrypt the states. 
static void inline aes_expand_key(uint8_t* aes_round_key, const uint8_t* Key)
{
    // The first round key is the key itself.
    memcpy (aes_round_key, Key, AES_KEYLEN);

    // All other round keys are found from the previous round keys.
    for (unsigned i = Nk; i < Nb * (Nr + 1); ++i) {
        
        uint8_t tempa[4]; 
  
        aes_assign_word (tempa, &aes_round_key[(i-1)*4]);

        if (i % Nk == 0) {
            aes_rotate_word (tempa);
            aes_substitute_word (tempa);
            tempa[0] = tempa[0] ^ Rcon[i/Nk];
        }

        if (i % Nk == 4) aes_substitute_word (tempa);

        aes_xor_word (&aes_round_key[i*4], &aes_round_key[(i-Nk)*4], tempa);
    }
}

// This function adds the round key to state.
// The round key is added to the state by an XOR function.
static inline void aes_add_round_key(uint8_t round, AesStateType* state, const uint8_t* aes_round_key)
{
    for (uint8_t i=0; i < Nb*Nb; i++) 
        state->v[i] ^= aes_round_key[(round * Nb * 4) + i];
}

// The aes_sub_bytes Function Substitutes the values in the
// state matrix with values in an S-box.
static inline void aes_sub_bytes (AesStateType* state)
{
    for (uint8_t i=0; i < Nb*Nb; i++) 
        state->v[i] = sbox[state->v[i]];
}

// The aes_shift_rows() function shifts the rows in the state to the left.
// Each row is shifted with different offset.
// Offset = Row number. So the first row is not shifted.
static inline void aes_shift_rows(AesStateType* state)
{
    uint8_t temp;

    // Rotate first row 1 columns to left  
    temp           = state->m[0][1];
    state->m[0][1] = state->m[1][1];
    state->m[1][1] = state->m[2][1];
    state->m[2][1] = state->m[3][1];
    state->m[3][1] = temp;

    // Rotate second row 2 columns to left  
    temp           = state->m[0][2];
    state->m[0][2] = state->m[2][2];
    state->m[2][2] = temp;

    temp           = state->m[1][2];
    state->m[1][2] = state->m[3][2];
    state->m[3][2] = temp;

    // Rotate third row 3 columns to left
    temp           = state->m[0][3];
    state->m[0][3] = state->m[3][3];
    state->m[3][3] = state->m[2][3];
    state->m[2][3] = state->m[1][3];
    state->m[1][3] = temp;
}

static inline uint8_t xtime(uint8_t x)
{
    return ((x<<1) ^ (((x>>7) & 1) * 0x1b));
}

// aes_mix_columns function mixes the columns of the state matrix
static inline void aes_mix_columns (AesStateType* state)
{
    uint8_t Tmp, Tm, t;
    for (uint8_t i=0; i < 4; ++i) {  
        t   = state->m[i][0];
        Tmp = state->m[i][0] ^ state->m[i][1] ^ state->m[i][2] ^ state->m[i][3] ;
        Tm  = state->m[i][0] ^ state->m[i][1] ; Tm = xtime(Tm);  state->m[i][0] ^= Tm ^ Tmp ;
        Tm  = state->m[i][1] ^ state->m[i][2] ; Tm = xtime(Tm);  state->m[i][1] ^= Tm ^ Tmp ;
        Tm  = state->m[i][2] ^ state->m[i][3] ; Tm = xtime(Tm);  state->m[i][2] ^= Tm ^ Tmp ;
        Tm  = state->m[i][3] ^ t ;              Tm = xtime(Tm);  state->m[i][3] ^= Tm ^ Tmp ;
    }
}

// aes_cipher is the main function that encrypts the PlainText.
static void aes_cipher(AesStateType* state, const uint8_t* aes_round_key)
{
    // Add the First round key to the state before starting the rounds.
    aes_add_round_key(0, state, aes_round_key); 
    
    // There will be Nr rounds.
    // The first Nr-1 rounds are identical.
    // These Nr-1 rounds are executed in the loop below.
    for (uint8_t round=1; round < Nr; ++round) {
        aes_sub_bytes(state);
        aes_shift_rows(state);
        aes_mix_columns(state);
        aes_add_round_key(round, state, aes_round_key);
    }
    
    // The last round is given below.
    // The aes_mix_columns function is not here in the last round.
    aes_sub_bytes(state);
    aes_shift_rows(state);
    aes_add_round_key(Nr, state, aes_round_key);
}

// Symmetrical operation: same function for encrypting as for decrypting. Note any IV/nonce should never be reused with the same key 
void aes_xcrypt_buffer (VBlock *vb, uint8_t *data, uint32_t length)
{
    AesStateType buffer; 

    for (unsigned i=0; i < length; i++, vb->bi++) {
 
        if (vb->bi == AES_BLOCKLEN) { /* we need to regen xor complement in buffer */
            memcpy (buffer.v, vb->aes_iv, AES_BLOCKLEN);
            aes_cipher (&buffer, vb->aes_round_key);

            /* Increment aes_iv and handle overflow */
            for (vb->bi = (AES_BLOCKLEN - 1); vb->bi >= 0; --vb->bi) {
                if (vb->aes_iv[vb->bi] == 255) {
                    vb->aes_iv[vb->bi] = 0;
                    continue;
                } 
                vb->aes_iv[vb->bi] += 1;
                break;   
            }
            vb->bi = 0;
        }

        data[i] = (data[i] ^ buffer.v[vb->bi]);
    }
}

void aes_initialize (VBlock *vb, const uint8_t* key)
{
    static const uint8_t iv_ad_120[AES_BLOCKLEN] = { 0, 1, 2, 10, 11, 12, 20, 21, 22, 100, 101, 102, 110, 111, 112, 120 };

    aes_expand_key(vb->aes_round_key, key);
    memcpy (vb->aes_iv, iv_ad_120, AES_BLOCKLEN);
    vb->bi = AES_BLOCKLEN;
}

char *aes_display_data (const uint8_t *data, unsigned data_len)
{
    char *str = MALLOC (data_len * 2 + 1);

    for (unsigned i=0; i < data_len; i++) 
        sprintf (&str[i*2], "%2.2x", data[i]);

    str[data_len*2] = 0;

    return str;
}

char *aes_display_key (const uint8_t* key) { return aes_display_data (key, AES_KEYLEN); }
