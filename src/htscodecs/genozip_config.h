// Genozip configuration for htscodecs

#ifdef __x86_64__
    #define HAVE_SSE4_1 1
    #define HAVE_SSSE3  1
    #define HAVE_POPCNT 1

    #ifdef __AVX2__
        #define HAVE_AVX2 1
    #endif
#endif
