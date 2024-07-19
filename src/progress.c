// ------------------------------------------------------------------
//   progress.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <time.h>
#include "genozip.h"
#include "file.h"
#include "progress.h"
#include "profiler.h" // for TimeSpecType
#include "flags.h"

static bool progress_newline_since_update = false; // note: might be not 100% with threads, but that's ok

static bool test_mode;
static TimeSpecType component_start_time;
static float last_percent=0;
static int last_seconds_so_far=-1;
static rom component_name=NULL;
static unsigned last_len=0; // so we know how many characters to erase on next update
static uint32_t last_secs_remaining=0xffff0000;

static StrText progress_ellapsed_time (void)
{
    TimeSpecType tb; 
    clock_gettime (CLOCK_REALTIME, &tb); 
    
    int seconds_so_far = ((tb.tv_sec - component_start_time.tv_sec)*1000 + ((int64_t)tb.tv_nsec - (int64_t)component_start_time.tv_nsec) / 1000000) / 1000; 

    return str_human_time (seconds_so_far, false);
}

static void progress_update_status (rom prefix, rom status)
{
    if (flag.quiet) return;

    #define eraser "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
    #define spaces "                                                                                                                                                                "

    ASSERT (STRLEN(eraser) == STRLEN(spaces), "eraser.len=%u != spaces.len=%u", STRLEN(eraser), STRLEN(spaces)); // will be optimized out
    
    if (prefix && prefix[0]) 
        iprintf ("%s", prefix);

    iprintf ("%.*s%.*s%.*s%s", MIN_(last_len, STRLEN(eraser)), eraser, MIN_(last_len, STRLEN(spaces)), spaces, last_len, eraser, status);

    last_len = strlen (status);

    progress_newline_since_update = false;
}

// bring the cursor down to a newline, if needed
void progress_newline(void) 
{
    if (!flag.quiet && !progress_newline_since_update) { 
        fputc ('\n', stderr);
        progress_newline_since_update = true;
    }
}

void progress_erase (void)
{
    progress_update_status (NULL, "");
}

void progress_new_component (rom new_component_name, 
                             rom message, // can be NULL
                             bool new_test_mode,
                             TimeSpecType *start_time) // optional: if time already started (e.g. when reading txt_header)
{
    StrTextSuperLong prefix = {};

    // (re) initialize if new component
    if (!component_name || strcmp (new_component_name, component_name)) {
        if (start_time)
            component_start_time = *start_time;
        else
            clock_gettime (CLOCK_REALTIME, &component_start_time); 

        test_mode = new_test_mode;
        component_name = new_component_name; 
        
        if (!flag.quiet) {
            
            if (test_mode) { 
                int global_cmd_len = strlen (global_cmd);
                int component_name_len = strlen (component_name);

                if (64 + global_cmd_len + component_name_len < sizeof(prefix)) { // fits in prefix 
                    char genounzip_str[component_name_len + 32];
                    if (strstr (global_cmd, "genozip")) {
                        char *zip = strstr (global_cmd, "zip");
                        snprintf (genounzip_str, sizeof(genounzip_str), "%.*sun%s", (int)(zip-global_cmd), global_cmd, zip);
                    }
                    else 
                        strcpy (genounzip_str, global_cmd);

                    snprintf (prefix.s, sizeof (prefix.s), "testing: %s %s : ", genounzip_str, component_name);
                }

                else
                    strcpy (prefix.s, "testing: ");

            }
            else if (flag.make_reference)
                snprintf (prefix.s, sizeof(prefix.s), "%saking reference file: %s %s : ",
                         txt_file->is_remote ? "Downloading & m" : "M", global_cmd, component_name);
            
            else
                snprintf (prefix.s, sizeof(prefix.s),"%s %s : ", global_cmd, component_name); 
        }
    }

    last_secs_remaining = 0xffff0000;
    progress_update_status (prefix.s, message ? message : "");
}

void progress_update (rom task, uint64_t sofar, uint64_t total, bool done)
{
    char progress_str[200];
    if (flag.quiet && !flag.debug_progress) return; 

    TimeSpecType tb; 
    clock_gettime(CLOCK_REALTIME, &tb); 
    
    int seconds_so_far = ((tb.tv_sec-component_start_time.tv_sec)*1000 + (tb.tv_nsec-component_start_time.tv_nsec) / 1000000) / 1000; 

    double percent;
    if (total > 10000000) // gentle handling of really big numbers to avoid integer overflow
        percent = MIN_(((double)(sofar/100000ULL)*100) / (double)(total/100000ULL), 100.0); // divide by 100000 to avoid number overflows
    else
        percent = MIN_(((double)sofar*100) / (double)total, 100.0); // divide by 100000 to avoid number overflows

    // need to update progress indicator, max once a second or if 100% is reached

    // case: we've reached 99% prematurely... we under-estimated the time
    if (!done && percent > 99 && (last_seconds_so_far < seconds_so_far)) {
        if (!flag.debug_progress)
            progress_update_status (NULL, "Finalizing...");
        else {
            snprintf (progress_str, sizeof(progress_str), "Finalizing... %u%% task=%s sofar=%.20s total=%.20s", 
                      (unsigned)percent, task, str_int_commas(sofar).s, str_int_commas(total).s);            
            progress_update_status (NULL, progress_str);
        }
    }
    
    // case: we're making progress... show % and time remaining
    else if (!done && percent && (last_seconds_so_far < seconds_so_far)) { 
        // time remaining
        uint32_t secs_remaining = (100.0 - percent) * ((double)seconds_so_far / (double)percent);

        if (!done && (percent != last_percent || secs_remaining <= last_secs_remaining || secs_remaining >= last_secs_remaining+15)) { // timer doesn't go up unless estimate changed by a good fews seconds

            if (!flag.debug_progress)
                snprintf (progress_str, sizeof(progress_str), "%u%% (%s)", (unsigned)percent, str_human_time (secs_remaining, false).s);
            else
                snprintf (progress_str, sizeof(progress_str), "%u%% (%s) task=%s sofar=%.20s total=%.20s seconds_so_far=%d", 
                         (unsigned)percent, str_human_time (secs_remaining, false).s, task, str_int_commas(sofar).s, str_int_commas(total).s, seconds_so_far);            

            progress_update_status (NULL, progress_str);

            last_secs_remaining = secs_remaining;
        }
    }

    // case: we're done - caller will print the "Done" message after finalizing the genozip header etc
    else {}

    last_percent = percent;
    last_seconds_so_far = seconds_so_far;
}

void progress_finalize_component (rom status)
{
    if (!flag.quiet
        && component_name /* not already finalized */) {
        progress_update_status (NULL, status);
        iprint0 ("\n");
    }

    component_name = NULL;
    last_len = 0;
}

bool progress_has_component (void)
{
    return component_name != NULL;
}

#define FINALIZE(format, ...) {                                                               \
    StrTextLong s; int s_len=0;                                                               \
    SNPRINTF (s, format, __VA_ARGS__);                                                        \
    if (IS_ZIP && flag.md5) SNPRINTF (s, "\t%s = %s", digest_name(), digest_display (md5).s); \
    progress_finalize_component (s.s);                                                        \
}

void progress_finalize_component_time (rom status, Digest md5/*ZIP only*/)
{
    FINALIZE ("%s (%s)", status, progress_ellapsed_time().s);
}

void progress_finalize_component_time_ratio (rom me, float ratio, Digest md5)
{
    if (component_name)
        FINALIZE ("Done (%s, %s compression ratio: %1.1f)", progress_ellapsed_time().s, me, ratio)
    else
        FINALIZE ("Time: %s, %s compression ratio: %1.1f", progress_ellapsed_time().s, me, ratio);
}

void progress_finalize_component_time_ratio_better (rom me, float ratio, rom better_than, float ratio_than, Digest md5)
{
    if (component_name) 
        FINALIZE ("Done (%s, %s compression ratio: %1.1f - better than %s by a factor of %1.1f)", 
                  progress_ellapsed_time().s, me, ratio, better_than, ratio_than)
    else
        FINALIZE ("Time: %s, %s compression ratio: %1.1f - better than %s by a factor of %1.1f", 
                  progress_ellapsed_time().s, me, ratio, better_than, ratio_than)
}

void progress_concatenated_md5 (rom me, Digest md5)
{
    ASSERT0 (!component_name, "expecting component_name=NULL");

    FINALIZE (flag.pair ? "Paired %s"
            : flag.deep ? "Deep %s + FASTQ" // SAM or BAM
            :             "Unexpeced"
            , me);
}


