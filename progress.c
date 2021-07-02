// ------------------------------------------------------------------
//   progress.c
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <time.h>
#include "genozip.h"
#include "file.h"
#include "progress.h"
#include "profiler.h" // for TimeSpecType
#include "flags.h"

bool progress_newline_since_update = false; // global: might be not 100% with threads, but that's ok

static bool ever_start_time_initialized = false, test_mode, show_progress;
static TimeSpecType ever_start_time, component_start_time;
static double last_percent=0;
static unsigned last_seconds_so_far=0;
static const char *component_name=NULL;
static unsigned last_len=0; // so we know how many characters to erase on next update

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

static const char *progress_ellapsed_time (bool ever)
{
    TimeSpecType tb; 
    clock_gettime(CLOCK_REALTIME, &tb); 
    
    TimeSpecType start = ever ? ever_start_time : component_start_time;

    int seconds_so_far = ((tb.tv_sec - start.tv_sec)*1000 + (tb.tv_nsec - start.tv_nsec) / 1000000) / 1000; 

    static char time_str[70];
    progress_human_time (seconds_so_far, time_str);

    return time_str;
}

char *progress_new_component (const char *new_component_name, 
                              const char *message, // can be NULL
                              int new_test_mode)   // true, false or -1 for unchanged
{
    char *prefix = NULL;

    // (re) initialize if new component
    if (!component_name || strcmp (new_component_name, component_name)) {
        clock_gettime (CLOCK_REALTIME, &component_start_time); 

        if (!ever_start_time_initialized) {
            ever_start_time = component_start_time;
            ever_start_time_initialized = true;
        }

        if (new_test_mode != -1) 
            test_mode = new_test_mode;

        // if !show_progress - we don't show the advancing %, but we still show the filename, done status, compression ratios etc
        show_progress  = !flag.quiet && !!isatty(2);
        component_name = new_component_name; 

        if (!flag.quiet) {
            prefix = MALLOC (100 + strlen (global_cmd) + strlen (component_name));
            if (test_mode) { 
                char genounzip_str[strlen(global_cmd)+3];
                if (strstr (global_cmd, "genozip")) {
                    char *zip = strstr (global_cmd, "zip");
                    sprintf (genounzip_str, "%.*sun%s", (int)(zip-global_cmd), global_cmd, zip);
                }
                else 
                    strcpy (genounzip_str, global_cmd);

                sprintf (prefix, "testing: %s %s : ", genounzip_str, component_name);
            }
            else
                sprintf (prefix, "%s %s : ", global_cmd, component_name); 
        }
    }

//    if (message) 
    if (!flag.reading_chain) 
        progress_update_status (&prefix, message ? message : "");

    return prefix;
}

void progress_update (char **prefix, uint64_t sofar, uint64_t total, bool done)
{
    char time_str[70], progress[200];
    if (!show_progress && !flag.debug_progress) return; 

    TimeSpecType tb; 
    clock_gettime(CLOCK_REALTIME, &tb); 
    
    int seconds_so_far = ((tb.tv_sec-component_start_time.tv_sec)*1000 + (tb.tv_nsec-component_start_time.tv_nsec) / 1000000) / 1000; 

    double percent;
    if (total > 10000000) // gentle handling of really big numbers to avoid integer overflow
        percent = MIN (((double)(sofar/100000ULL)*100) / (double)(total/100000ULL), 100.0); // divide by 100000 to avoid number overflows
    else
        percent = MIN (((double)sofar*100) / (double)total, 100.0); // divide by 100000 to avoid number overflows

    // need to update progress indicator, max once a second or if 100% is reached

    // case: we've reached 99% prematurely... we under-estimated the time
    if (!done && percent > 99 && (last_seconds_so_far < seconds_so_far)) {
        if (!flag.debug_progress)
            progress_update_status (prefix, "Finalizing...");
        else {
            sprintf (progress, "Finalizing... %u%% sofar=%"PRIu64" total=%"PRIu64, (unsigned)percent, sofar, total);            
            progress_update_status (prefix, progress);
        }
    }
    
    // case: we're making progress... show % and time remaining
    else if (!done && percent && (last_seconds_so_far < seconds_so_far)) { 

        if (!done) { 

            // time remaining
            unsigned secs = (100.0 - percent) * ((double)seconds_so_far / (double)percent);
            progress_human_time (secs, time_str);

            if (!flag.debug_progress)
                sprintf (progress, "%u%% (%s)", (unsigned)percent, time_str);
            else
                sprintf (progress, "%u%% (%s) sofar=%"PRIu64" total=%"PRIu64" seconds_so_far=%d", (unsigned)percent, time_str, sofar, total, seconds_so_far);            

            progress_update_status (prefix, progress);
        }
    }

    // case: we're done - caller will print the "Done" message after finalizing the genozip header etc
    else {}

    last_percent = percent;
    last_seconds_so_far = seconds_so_far;
}

void progress_update_status (char **prefix, const char *status)
{
    if (flag.quiet) return;

    static const char *eraser = "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
    static const char *spaces = "                                                                                ";

    if (prefix && *prefix) {
        iprintf ("%s", *prefix);
        FREE (*prefix);
    }

    if (is_info_stream_terminal && !flag.debug_progress) 
        iprintf ("%.*s%.*s%.*s%s", last_len, eraser, last_len, spaces, last_len, eraser, status);
    else 
        iprintf ("%s\n", status); // if we're outputting to a log file or debugging progress, show every status on its own line

    last_len = strlen (status);

    progress_newline_since_update = false;
}

void progress_finalize_component (const char *status)
{
    if (!flag.quiet) {
        progress_update_status (NULL, status);
        iprint0 ("\n");
    }

    component_name = NULL;
    last_len = 0;
}

#define FINALIZE(format, ...) { \
    char s[500]; \
    sprintf (s, format, __VA_ARGS__);  \
    if (!digest_is_zero (md5) && flag.md5) sprintf (&s[strlen(s)], "\t%s = %s", digest_name(), digest_display (md5).s); \
    progress_finalize_component (s);  \
}

void progress_finalize_component_time (const char *status, Digest md5)
{
    FINALIZE ("%s (%s)", status, progress_ellapsed_time (false));
}

void progress_finalize_component_time_ratio (const char *me, double ratio, Digest md5)
{
    if (component_name)
        FINALIZE ("Done (%s, %s compression ratio: %1.1f)", progress_ellapsed_time (false), me, ratio)
    else
        FINALIZE ("Time: %s, %s compression ratio: %1.1f", progress_ellapsed_time (false), me, ratio);
}

void progress_finalize_component_time_ratio_better (const char *me, double ratio, const char *better_than, double ratio_than, Digest md5)
{
    if (component_name) 
        FINALIZE ("Done (%s, %s compression ratio: %1.1f - better than %s by a factor of %1.1f)", 
                  progress_ellapsed_time (false), me, ratio, better_than, ratio_than)
    else
        FINALIZE ("Time: %s, %s compression ratio: %1.1f - better than %s by a factor of %1.1f", 
                  progress_ellapsed_time (false), me, ratio, better_than, ratio_than)
}

void progress_concatenated_md5 (const char *me, Digest md5)
{
    ASSERT0 (!component_name, "expecting component_name=NULL");

    FINALIZE ("Bound %s", me);
}
