// ------------------------------------------------------------------
//   progress.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <time.h>
#include "genozip.h"
#include "progress.h"
#include "profiler.h" // for TimeSpecType
#if defined __APPLE__ 
#include "compatibility/mac_gettime.h"
#endif

static bool ever_start_time_initialized = false, test_mode, show_progress;
static TimeSpecType ever_start_time, component_start_time;
static unsigned last_len=0; // so progress_update knows how many characters to erase
static double last_percent=0;
static unsigned last_seconds_so_far=0;
static const char *component_name=NULL;

static void progress_human_time (unsigned secs, char *str /* out */)
{
    unsigned hours = secs / 3600;
    unsigned mins  = (secs % 3600) / 60;
             secs  = secs % 60;

    if (hours) 
        sprintf (str, "%u %s %u %s", hours, hours==1 ? "hour" : "hours", mins, mins==1 ? "minute" : "minutes");
    else if (mins)
        sprintf (str, "%u %s %u %s", mins, mins==1 ? "minute" : "minutes", secs, secs==1 ? "second" : "seconds");
    else 
        sprintf (str, "%u %s", secs, secs==1 ? "second" : "seconds");
}

const char *progress_ellapsed_time (bool ever)
{
    TimeSpecType tb; 
    clock_gettime(CLOCK_REALTIME, &tb); 
    
    TimeSpecType start = ever ? ever_start_time : component_start_time;

    int seconds_so_far = ((tb.tv_sec - start.tv_sec)*1000 + (tb.tv_nsec - start.tv_nsec) / 1000000) / 1000; 

    static char time_str[70];
    progress_human_time (seconds_so_far, time_str);

    return time_str;
}

void progress_new_component (const char *new_component_name, 
                             bool txt_file_size_unknown,
                             int new_test_mode) // true, false or -1 for unchanged
{
    // (re) initialize if new component
    if (new_component_name != component_name) { // pointer comparison
        clock_gettime(CLOCK_REALTIME, &component_start_time); 

        if (!ever_start_time_initialized) {
            ever_start_time = component_start_time;
            ever_start_time_initialized = true;
        }

        if (new_test_mode != -1) 
            test_mode = new_test_mode;
        else
            last_len = 2; // unbind mode

        show_progress  = !flag_quiet && !!isatty(2);
        component_name = new_component_name; 
    }

    if (!show_progress) return; 

    const char *progress = txt_file_size_unknown ? "Compressing...\b\b\b\b\b\b\b\b\b\b\b\b\b\b" : "0\%"; // we can't show % when compressing from stdin as we don't know the file size

    if (test_mode) 
        fprintf (stderr, "testing: %s%s --test %s : %s", global_cmd, strstr (global_cmd, "genozip") ? " --decompress" : "", 
                 new_component_name, progress); 
    else
        fprintf (stderr, "%s %s : %s", global_cmd, new_component_name, progress); 
             
    last_len = strlen (progress); // so progress_update knows how many characters to erase
}

void progress_update (uint64_t sofar, uint64_t total, bool done)
{
    if (!show_progress && !flag_debug_progress) return; 

    TimeSpecType tb; 
    clock_gettime(CLOCK_REALTIME, &tb); 
    
    int seconds_so_far = ((tb.tv_sec-component_start_time.tv_sec)*1000 + (tb.tv_nsec-component_start_time.tv_nsec) / 1000000) / 1000; 

    double percent;
    if (total > 10000000) // gentle handling of really big numbers to avoid integer overflow
        percent = MIN (((double)(sofar/100000ULL)*100) / (double)(total/100000ULL), 100.0); // divide by 100000 to avoid number overflows
    else
        percent = MIN (((double)sofar*100) / (double)total, 100.0); // divide by 100000 to avoid number overflows
    
    // need to update progress indicator, max once a second or if 100% is reached
    const char *eraser = "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";

    // case: we've reached 99% prematurely... we under-estimated the time
    if (!done && percent > 99 && (last_seconds_so_far < seconds_so_far)) {
        const char *progress = "Finalizing...";

        // note we have spaces at the end to make sure we erase the previous string, if it is longer than the current one
        fprintf (stderr, "%.*s%s            %.12s", last_len, eraser, progress, eraser);

        last_len = strlen (progress);
    }

    // case: we're making progress... show % and time remaining
    else if (!done && percent && (last_seconds_so_far < seconds_so_far)) { 

        if (!done) { 

            // time remaining
            char time_str[70], progress[100];
            unsigned secs = (100.0 - percent) * ((double)seconds_so_far / (double)percent);
            progress_human_time (secs, time_str);
            sprintf (progress, "%u%% (%s)", (unsigned)percent, time_str);

            // note we have spaces at the end to make sure we erase the previous string, if it is longer than the current one
            if (!flag_debug_progress)
                fprintf (stderr, "%.*s%s            %.12s", last_len, eraser, progress, eraser);
            else
                fprintf (stderr, "%u%% (%s) sofar=%"PRIu64" total=%"PRIu64" seconds_so_far=%d\n", (unsigned)percent, time_str, sofar, total, seconds_so_far);

            last_len = strlen (progress);
        }
    }

    // case: we're done - caller will print the "Done" message after finalizing the genozip header etc
    else if (done && !flag_quiet) 
        fprintf (stderr, "%.*s", last_len, eraser); 

    last_percent = percent;
    last_seconds_so_far = seconds_so_far;
}
