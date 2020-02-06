// This is an implementation of MD5 based on: https://github.com/WaterJuice/WjCryptLib/blob/master/lib/WjCryptLib_Md5.c
// which has been modified
// 
// The original source, as specified there, is "This is free and unencumbered software released into the public domain - June 2013 waterjuice.org". 
// All modifications are (c) 2020 Divon Lan and are subject to license.

#include <memory.h>
#include "genozip.h"

#define F( x, y, z )            ( (z) ^ ((x) & ((y) ^ (z))) )
#define G( x, y, z )            ( (y) ^ ((z) & ((x) ^ (y))) )
#define H( x, y, z )            ( (x) ^ (y) ^ (z) )
#define I( x, y, z )            ( (y) ^ ((x) | ~(z)) )

#define STEP( f, a, b, c, d, x, t, s )                          \
    (a) += f((b), (c), (d)) + (x) + (t);                        \
    (a) = (((a) << (s)) | (((a) & 0xffffffff) >> (32 - (s))));  \
    (a) += (b);

void md5_display_ctx (const Md5Context *x)
{
    static unsigned iteration=1;

    printf ("%2u: %08x %08x %08x %08x %08x %08x ", iteration, x->hi, x->lo, x->a, x->b, x->c, x->d);
    for (unsigned i=0; i<64; i++) printf ("%2.2x", x->buffer.bytes[i]);
    printf ("\n");

    iteration++;
}


static const void *md5_transform (Md5Context *ctx, const void *data, uintmax_t size)
{
    const uint32_t *ptr = (uint32_t *)data;

    uint32_t a = ctx->a;
    uint32_t b = ctx->b;
    uint32_t c = ctx->c;
    uint32_t d = ctx->d;

    do {
        uint32_t saved_a = a;
        uint32_t saved_b = b;
        uint32_t saved_c = c;
        uint32_t saved_d = d;

#ifdef __LITTLE_ENDIAN__
        // in little endian machines, we can access the data directly - we don't need to copy memory
        const uint32_t *block = ptr;
#else
        // in big endian machines - we need to flip the data to little endian - so we do it in a copy
        uint32_t block[16]; 
        for (unsigned i=0; i < 16; i++) block[i] = LTEN32(ptr[i]);
#endif
        // Round 1
        STEP( F, a, b, c, d, block[0],  0xd76aa478, 7  )
        STEP( F, d, a, b, c, block[1],  0xe8c7b756, 12 )
        STEP( F, c, d, a, b, block[2],  0x242070db, 17 )
        STEP( F, b, c, d, a, block[3],  0xc1bdceee, 22 )
        STEP( F, a, b, c, d, block[4],  0xf57c0faf, 7  )
        STEP( F, d, a, b, c, block[5],  0x4787c62a, 12 )
        STEP( F, c, d, a, b, block[6],  0xa8304613, 17 )
        STEP( F, b, c, d, a, block[7],  0xfd469501, 22 )
        STEP( F, a, b, c, d, block[8],  0x698098d8, 7  )
        STEP( F, d, a, b, c, block[9],  0x8b44f7af, 12 )
        STEP( F, c, d, a, b, block[10], 0xffff5bb1, 17 )
        STEP( F, b, c, d, a, block[11], 0x895cd7be, 22 )
        STEP( F, a, b, c, d, block[12], 0x6b901122, 7  )
        STEP( F, d, a, b, c, block[13], 0xfd987193, 12 )
        STEP( F, c, d, a, b, block[14], 0xa679438e, 17 )
        STEP( F, b, c, d, a, block[15], 0x49b40821, 22 )

        // Round 2
        STEP( G, a, b, c, d, block[1],  0xf61e2562, 5  )
        STEP( G, d, a, b, c, block[6],  0xc040b340, 9  )
        STEP( G, c, d, a, b, block[11], 0x265e5a51, 14 )
        STEP( G, b, c, d, a, block[0],  0xe9b6c7aa, 20 )
        STEP( G, a, b, c, d, block[5],  0xd62f105d, 5  )
        STEP( G, d, a, b, c, block[10], 0x02441453, 9  )
        STEP( G, c, d, a, b, block[15], 0xd8a1e681, 14 )
        STEP( G, b, c, d, a, block[4],  0xe7d3fbc8, 20 )
        STEP( G, a, b, c, d, block[9],  0x21e1cde6, 5  )
        STEP( G, d, a, b, c, block[14], 0xc33707d6, 9  )
        STEP( G, c, d, a, b, block[3],  0xf4d50d87, 14 )
        STEP( G, b, c, d, a, block[8],  0x455a14ed, 20 )
        STEP( G, a, b, c, d, block[13], 0xa9e3e905, 5  )
        STEP( G, d, a, b, c, block[2],  0xfcefa3f8, 9  )
        STEP( G, c, d, a, b, block[7],  0x676f02d9, 14 )
        STEP( G, b, c, d, a, block[12], 0x8d2a4c8a, 20 )

        // Round 3
        STEP( H, a, b, c, d, block[5],  0xfffa3942, 4  )
        STEP( H, d, a, b, c, block[8],  0x8771f681, 11 )
        STEP( H, c, d, a, b, block[11], 0x6d9d6122, 16 )
        STEP( H, b, c, d, a, block[14], 0xfde5380c, 23 )
        STEP( H, a, b, c, d, block[1],  0xa4beea44, 4  )
        STEP( H, d, a, b, c, block[4],  0x4bdecfa9, 11 )
        STEP( H, c, d, a, b, block[7],  0xf6bb4b60, 16 )
        STEP( H, b, c, d, a, block[10], 0xbebfbc70, 23 )
        STEP( H, a, b, c, d, block[13], 0x289b7ec6, 4  )
        STEP( H, d, a, b, c, block[0],  0xeaa127fa, 11 )
        STEP( H, c, d, a, b, block[3],  0xd4ef3085, 16 )
        STEP( H, b, c, d, a, block[6],  0x04881d05, 23 )
        STEP( H, a, b, c, d, block[9],  0xd9d4d039, 4  )
        STEP( H, d, a, b, c, block[12], 0xe6db99e5, 11 )
        STEP( H, c, d, a, b, block[15], 0x1fa27cf8, 16 )
        STEP( H, b, c, d, a, block[2],  0xc4ac5665, 23 )

        // Round 4
        STEP( I, a, b, c, d, block[0],  0xf4292244, 6  )
        STEP( I, d, a, b, c, block[7],  0x432aff97, 10 )
        STEP( I, c, d, a, b, block[14], 0xab9423a7, 15 )
        STEP( I, b, c, d, a, block[5],  0xfc93a039, 21 )
        STEP( I, a, b, c, d, block[12], 0x655b59c3, 6  )
        STEP( I, d, a, b, c, block[3],  0x8f0ccc92, 10 )
        STEP( I, c, d, a, b, block[10], 0xffeff47d, 15 )
        STEP( I, b, c, d, a, block[1],  0x85845dd1, 21 )
        STEP( I, a, b, c, d, block[8],  0x6fa87e4f, 6  )
        STEP( I, d, a, b, c, block[15], 0xfe2ce6e0, 10 )
        STEP( I, c, d, a, b, block[6],  0xa3014314, 15 )
        STEP( I, b, c, d, a, block[13], 0x4e0811a1, 21 )
        STEP( I, a, b, c, d, block[4],  0xf7537e82, 6  )
        STEP( I, d, a, b, c, block[11], 0xbd3af235, 10 )
        STEP( I, c, d, a, b, block[2],  0x2ad7d2bb, 15 )
        STEP( I, b, c, d, a, block[9],  0xeb86d391, 21 )

        a += saved_a;   
        b += saved_b;
        c += saved_c;
        d += saved_d;

        ptr += 16;
    } while (size -= 64);

    ctx->a = a;
    ctx->b = b;
    ctx->c = c;
    ctx->d = d;

    return ptr;
}

void md5_initialize (Md5Context *ctx)
{
    memset (ctx, 0, sizeof(Md5Context));

    ctx->a = 0x67452301;
    ctx->b = 0xefcdab89;
    ctx->c = 0x98badcfe;
    ctx->d = 0x10325476;

    ctx->lo = 0;
    ctx->hi = 0;
}

void md5_update (Md5Context *ctx, const void *data, unsigned len, bool initialize)
{
    if (initialize) md5_initialize (ctx);

    uint32_t    saved_lo;
    uint32_t    used;
    uint32_t    free;

    saved_lo = ctx->lo;
    if ((ctx->lo = (saved_lo + len) & 0x1fffffff) < saved_lo) 
        ctx->hi++;
    
    ctx->hi += (uint32_t)(len >> 29);

    used = saved_lo & 0x3f;

    if (used) {
        free = 64 - used;

        if (len < free) {
            memcpy (&ctx->buffer.bytes[used], data, len);
            goto finish;
        }

        memcpy (&ctx->buffer.bytes[used], data, free);
        data += free;
        len -= free;
        md5_transform (ctx, ctx->buffer.bytes, 64);
    }

    if (len >= 64) {
        data = md5_transform (ctx, data, len & ~(unsigned long)0x3f);
        len &= 0x3f;
    }

    memcpy (ctx->buffer.bytes, data, len);

finish:
    //md5_display_ctx (ctx);
    return;
}

void md5_finalize (Md5Context *ctx, Md5Hash *digest)
{
    uint32_t    used;
    uint32_t    free;

    used = ctx->lo & 0x3f;

    ctx->buffer.bytes[used++] = 0x80;

    free = 64 - used;

    if (free < 8) {
        memset (&ctx->buffer.bytes[used], 0, free);
        md5_transform (ctx, ctx->buffer.bytes, 64);
        used = 0;
        free = 64;
    }

    memset (&ctx->buffer.bytes[used], 0, free - 8);

    ctx->lo <<= 3;
    ctx->buffer.words[14] = LTEN32 (ctx->lo);
    ctx->buffer.words[15] = LTEN32 (ctx->hi);

    md5_transform (ctx, ctx->buffer.bytes, 64);
    digest->words[0] = LTEN32 (ctx->a);
    digest->words[1] = LTEN32 (ctx->b);
    digest->words[2] = LTEN32 (ctx->c);
    digest->words[3] = LTEN32 (ctx->d);
}

void md5_do (const void *data, unsigned len, Md5Hash *digest)
{
    Md5Context ctx;

    md5_update (&ctx, data, len, true);

    md5_finalize (&ctx, digest);
}

const char *md5_display (const Md5Hash *digest, bool prefix_space)
{
    char *str = malloc (34); // we're going to leak this memory - nevermind, it is small and rare

    const uint8_t *b = digest->bytes; 
    
    if (digest->ulls[0] || digest->ulls[1])
        sprintf (str, "%s%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x", prefix_space ? " " : "",
                 b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12], b[13], b[14], b[15]);
    else
        sprintf (str, "%sN/A                             ", prefix_space ? " " : "");
    
    return str;
}