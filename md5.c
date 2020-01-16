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
    for (unsigned i=0; i<64; i++) printf ("%2.2x", x->buffer[i]);
    printf (" ");
    for (unsigned i=0; i<16; i++) printf ("%8.8x", x->block[i]);
    printf ("\n");

    iteration++;
}


static void *md5_transform (Md5Context *ctx, void const *data, uintmax_t size)
{
    uint8_t*     ptr;
    uint32_t     a;
    uint32_t     b;
    uint32_t     c;
    uint32_t     d;
    uint32_t     saved_a;
    uint32_t     saved_b;
    uint32_t     saved_c;
    uint32_t     saved_d;

    #define GET(n) (ctx->block[(n)])
    #define SET(n) (ctx->block[(n)] =             \
            ((uint32_t)ptr[(n)*4 + 0] << 0 )      \
        |   ((uint32_t)ptr[(n)*4 + 1] << 8 )      \
        |   ((uint32_t)ptr[(n)*4 + 2] << 16)      \
        |   ((uint32_t)ptr[(n)*4 + 3] << 24) )

    ptr = (uint8_t*)data;

    a = ctx->a;
    b = ctx->b;
    c = ctx->c;
    d = ctx->d;

    do
    {
        saved_a = a;
        saved_b = b;
        saved_c = c;
        saved_d = d;

        // Round 1
        STEP( F, a, b, c, d, SET(0),  0xd76aa478, 7 )
        STEP( F, d, a, b, c, SET(1),  0xe8c7b756, 12 )
        STEP( F, c, d, a, b, SET(2),  0x242070db, 17 )
        STEP( F, b, c, d, a, SET(3),  0xc1bdceee, 22 )
        STEP( F, a, b, c, d, SET(4),  0xf57c0faf, 7 )
        STEP( F, d, a, b, c, SET(5),  0x4787c62a, 12 )
        STEP( F, c, d, a, b, SET(6),  0xa8304613, 17 )
        STEP( F, b, c, d, a, SET(7),  0xfd469501, 22 )
        STEP( F, a, b, c, d, SET(8 ),  0x698098d8, 7 )
        STEP( F, d, a, b, c, SET(9 ),  0x8b44f7af, 12 )
        STEP( F, c, d, a, b, SET(10 ), 0xffff5bb1, 17 )
        STEP( F, b, c, d, a, SET(11 ), 0x895cd7be, 22 )
        STEP( F, a, b, c, d, SET(12 ), 0x6b901122, 7 )
        STEP( F, d, a, b, c, SET(13 ), 0xfd987193, 12 )
        STEP( F, c, d, a, b, SET(14 ), 0xa679438e, 17 )
        STEP( F, b, c, d, a, SET(15 ), 0x49b40821, 22 )

        // Round 2
        STEP( G, a, b, c, d, GET(1),  0xf61e2562, 5 )
        STEP( G, d, a, b, c, GET(6),  0xc040b340, 9 )
        STEP( G, c, d, a, b, GET(11), 0x265e5a51, 14 )
        STEP( G, b, c, d, a, GET(0),  0xe9b6c7aa, 20 )
        STEP( G, a, b, c, d, GET(5),  0xd62f105d, 5 )
        STEP( G, d, a, b, c, GET(10), 0x02441453, 9 )
        STEP( G, c, d, a, b, GET(15), 0xd8a1e681, 14 )
        STEP( G, b, c, d, a, GET(4),  0xe7d3fbc8, 20 )
        STEP( G, a, b, c, d, GET(9),  0x21e1cde6, 5 )
        STEP( G, d, a, b, c, GET(14), 0xc33707d6, 9 )
        STEP( G, c, d, a, b, GET(3),  0xf4d50d87, 14 )
        STEP( G, b, c, d, a, GET(8),  0x455a14ed, 20 )
        STEP( G, a, b, c, d, GET(13), 0xa9e3e905, 5 )
        STEP( G, d, a, b, c, GET(2),  0xfcefa3f8, 9 )
        STEP( G, c, d, a, b, GET(7),  0x676f02d9, 14 )
        STEP( G, b, c, d, a, GET(12), 0x8d2a4c8a, 20 )

        // Round 3
        STEP( H, a, b, c, d, GET(5),  0xfffa3942, 4 )
        STEP( H, d, a, b, c, GET(8),  0x8771f681, 11 )
        STEP( H, c, d, a, b, GET(11), 0x6d9d6122, 16 )
        STEP( H, b, c, d, a, GET(14), 0xfde5380c, 23 )
        STEP( H, a, b, c, d, GET(1),  0xa4beea44, 4 )
        STEP( H, d, a, b, c, GET(4),  0x4bdecfa9, 11 )
        STEP( H, c, d, a, b, GET(7),  0xf6bb4b60, 16 )
        STEP( H, b, c, d, a, GET(10), 0xbebfbc70, 23 )
        STEP( H, a, b, c, d, GET(13), 0x289b7ec6, 4 )
        STEP( H, d, a, b, c, GET(0),  0xeaa127fa, 11 )
        STEP( H, c, d, a, b, GET(3),  0xd4ef3085, 16 )
        STEP( H, b, c, d, a, GET(6),  0x04881d05, 23 )
        STEP( H, a, b, c, d, GET(9),  0xd9d4d039, 4 )
        STEP( H, d, a, b, c, GET(12), 0xe6db99e5, 11 )
        STEP( H, c, d, a, b, GET(15), 0x1fa27cf8, 16 )
        STEP( H, b, c, d, a, GET(2),  0xc4ac5665, 23 )

        // Round 4
        STEP( I, a, b, c, d, GET(0),  0xf4292244, 6 )
        STEP( I, d, a, b, c, GET(7),  0x432aff97, 10 )
        STEP( I, c, d, a, b, GET(14), 0xab9423a7, 15 )
        STEP( I, b, c, d, a, GET(5),  0xfc93a039, 21 )
        STEP( I, a, b, c, d, GET(12), 0x655b59c3, 6 )
        STEP( I, d, a, b, c, GET(3),  0x8f0ccc92, 10 )
        STEP( I, c, d, a, b, GET(10), 0xffeff47d, 15 )
        STEP( I, b, c, d, a, GET(1),  0x85845dd1, 21 )
        STEP( I, a, b, c, d, GET(8),  0x6fa87e4f, 6 )
        STEP( I, d, a, b, c, GET(15), 0xfe2ce6e0, 10 )
        STEP( I, c, d, a, b, GET(6),  0xa3014314, 15 )
        STEP( I, b, c, d, a, GET(13), 0x4e0811a1, 21 )
        STEP( I, a, b, c, d, GET(4),  0xf7537e82, 6 )
        STEP( I, d, a, b, c, GET(11), 0xbd3af235, 10 )
        STEP( I, c, d, a, b, GET(2),  0x2ad7d2bb, 15 )
        STEP( I, b, c, d, a, GET(9),  0xeb86d391, 21 )

        a += saved_a;
        b += saved_b;
        c += saved_c;
        d += saved_d;

        ptr += 64;
    } while( size -= 64 );

    ctx->a = a;
    ctx->b = b;
    ctx->c = c;
    ctx->d = d;

    #undef GET
    #undef SET

    return ptr;
}

static void md5_initialize (Md5Context *ctx)
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
    printf ("xxx data: len=%u\n", len);
    //for (unsigned i=0; i<len; i++) printf ("%c", ((uint8_t*)data)[i]);
    //printf("\n"); 

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
            memcpy (&ctx->buffer[used], data, len);
            goto finish;
        }

        memcpy (&ctx->buffer[used], data, free);
        data += free;
        len -= free;
        md5_transform (ctx, ctx->buffer, 64);
    }

    if (len >= 64) {
        data = md5_transform (ctx, data, len & ~(unsigned long)0x3f);
        len &= 0x3f;
    }

    memcpy (ctx->buffer, data, len);

finish:
    md5_display_ctx (ctx);
    return;
}

void md5_finalize (Md5Context *ctx, Md5Hash *digest)
{
    uint32_t    used;
    uint32_t    free;

    used = ctx->lo & 0x3f;

    ctx->buffer[used++] = 0x80;

    free = 64 - used;

    if(free < 8)
    {
        memset( &ctx->buffer[used], 0, free );
        md5_transform( ctx, ctx->buffer, 64 );
        used = 0;
        free = 64;
    }

    memset( &ctx->buffer[used], 0, free - 8 );

    ctx->lo <<= 3;
    ctx->buffer[56] = (uint8_t)( ctx->lo );
    ctx->buffer[57] = (uint8_t)( ctx->lo >> 8 );
    ctx->buffer[58] = (uint8_t)( ctx->lo >> 16 );
    ctx->buffer[59] = (uint8_t)( ctx->lo >> 24 );
    ctx->buffer[60] = (uint8_t)( ctx->hi );
    ctx->buffer[61] = (uint8_t)( ctx->hi >> 8 );
    ctx->buffer[62] = (uint8_t)( ctx->hi >> 16 );
    ctx->buffer[63] = (uint8_t)( ctx->hi >> 24 );

    md5_transform( ctx, ctx->buffer, 64 );

    digest->bytes[0]  = (uint8_t)( ctx->a );
    digest->bytes[1]  = (uint8_t)( ctx->a >> 8 );
    digest->bytes[2]  = (uint8_t)( ctx->a >> 16 );
    digest->bytes[3]  = (uint8_t)( ctx->a >> 24 );
    digest->bytes[4]  = (uint8_t)( ctx->b );
    digest->bytes[5]  = (uint8_t)( ctx->b >> 8 );
    digest->bytes[6]  = (uint8_t)( ctx->b >> 16 );
    digest->bytes[7]  = (uint8_t)( ctx->b >> 24 );
    digest->bytes[8]  = (uint8_t)( ctx->c );
    digest->bytes[9]  = (uint8_t)( ctx->c >> 8 );
    digest->bytes[10] = (uint8_t)( ctx->c >> 16 );
    digest->bytes[11] = (uint8_t)( ctx->c >> 24 );
    digest->bytes[12] = (uint8_t)( ctx->d );
    digest->bytes[13] = (uint8_t)( ctx->d >> 8 );
    digest->bytes[14] = (uint8_t)( ctx->d >> 16 );
    digest->bytes[15] = (uint8_t)( ctx->d >> 24 );
}

void md5_do (const void *data, unsigned len, Md5Hash *digest)
{
    Md5Context ctx;

    md5_update (&ctx, data, len, true);

    md5_finalize (&ctx, digest);
}

const char *md5_display (const Md5Hash *digest, bool prefix_space)
{
    static char str[34]; // not thread-safe, but only used by main thread (inc. space and \0)

    const uint8_t *b = digest->bytes; 
    
    sprintf (str, "%s%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x", prefix_space ? " " : "",
             b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12], b[13], b[14], b[15]);

    return str;
}