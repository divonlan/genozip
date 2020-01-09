// ------------------------------------------------------------------
//   visual_c_gettime.h
//   Copyright (C) 2020 Divon Lan <genozip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifdef _WIN32

#define CLOCK_REALTIME 0

extern int clock_gettime (int unused, struct timeval *tv);

#endif