#if defined __APPLE__

#include <time.h>
#include <sys/types.h>
#include <sys/_types/_timespec.h>
#include <mach/mach.h>
#include <mach/clock.h>

#ifndef CLOCK_REALTIME

#ifndef mach_time_h
#define mach_time_h

// only supports a single timer
#define TIMER_ABSTIME -1
#define CLOCK_REALTIME CALENDAR_CLOCK
#define CLOCK_MONOTONIC SYSTEM_CLOCK

typedef int clockid_t;

extern int clock_gettime(clockid_t clk_id, struct timespec *tp);
#endif

#endif

#endif

/*  Copyright (c) 2015-2018 Alf Watt - Open Source - https://opensource.org/licenses/MIT */

