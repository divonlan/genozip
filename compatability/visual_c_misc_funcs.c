// ------------------------------------------------------------------
//   visual_c_misc_funcs.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifdef _MSC_VER

#include <windows.h>
#include "../genozip.h"

void usleep(uint64_t usec) 
{ 
    Sleep (usec / 1000); // convert microseconds to milliseconds
}

double log2 (double n)  
{  
    static double loge2 = 0; // really a constant
    if (!loge2) loge2 = log((double)2);

    return log(n) / loge2;  
}  

int my_round (double n)
{
    double floor_n = floor(n);

    return (int)floor_n + (n-floor_n >= 0.5);
}

#endif