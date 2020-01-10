// ------------------------------------------------------------------
//   visual_c_gettime.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifdef _MSC_VER

#include <time.h>

#define CLOCK_REALTIME 0

struct my_timespec { 
    time_t tv_sec;
    long tv_nsec;
};

extern int clock_gettime (int unused, struct my_timespec *ts);

#endif