// ------------------------------------------------------------------
//   dispatcher.c
//   Copyright (C) 2020 Divon Lan <genozip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"

typedef struct {
    pthread_t thread_id;
    VariantBlock *vb;
    void (*func)(VariantBlock *);
} Thread;

typedef struct {
    VariantBlockPool *vb_pool;
    VariantBlock *pseudo_vb;
    Buffer compute_threads_buf;
    Thread *compute_threads;
    VariantBlock *next_vb; // next vb to be dispatched
    VariantBlock *processed_vb; // last vb for which caller got the processing results

    bool input_exhausted;

    unsigned next_thread_to_dispatched;
    unsigned next_thread_to_be_joined;

    unsigned num_running_compute_threads;
    unsigned next_vb_i;
    unsigned max_threads;
    File *vcf_file;
    File *z_file;
    bool test_mode;

    // progress indicator stuff
    struct timespec start_time; 
    bool show_progress;
    double last_percent;
    unsigned last_len;
    unsigned last_seconds_so_far;
    const char *filename;
} DispatcherData;

static struct timespec profiler_timer; // wallclock

Dispatcher dispatcher_init (unsigned max_threads, unsigned pool_id, File *vcf_file, File *z_file,
                            bool test_mode, bool show_progress, const char *filename)
{
    clock_gettime(CLOCK_REALTIME, &profiler_timer);

    DispatcherData *dd = calloc (1, sizeof(DispatcherData));
    ASSERT0 (dd, "failed to calloc DispatcherData");

    clock_gettime(CLOCK_REALTIME, &dd->start_time); 

    dd->max_threads   = max_threads;
    dd->vcf_file      = vcf_file;
    dd->z_file        = z_file;
    dd->test_mode     = test_mode;
    dd->show_progress = show_progress && !flag_quiet && !!isatty(2);
    dd->filename      = filename;
    dd->last_len      = 2;

    dd->vb_pool = vb_construct_pool (MAX (2,max_threads)+1 /* +1 for pseudo-vb */, pool_id);

    dd->pseudo_vb = vb_get_vb (dd->vb_pool, vcf_file, z_file, 0);

    buf_alloc (dd->pseudo_vb, &dd->compute_threads_buf, sizeof(Thread) * MAX (1, max_threads-1), 1, "compute_threads_buf", 0);
    dd->compute_threads = (Thread *)dd->compute_threads_buf.data;

    if (dd->show_progress) 
        fprintf (stderr, "%s%s %s: 0%%", dd->test_mode ? "testing " : "", global_cmd, dd->filename);

    return dd;
}

void dispatcher_finish (Dispatcher dispatcher)
{
    DispatcherData *dd = dispatcher;

    COPY_TIMER (dd->pseudo_vb->profile.wallclock);

    if (flag_show_time) 
        profiler_print_report (&dd->pseudo_vb->profile, dd->max_threads, dd->filename, dd->next_vb_i-1);

    // must be before vb_cleanup_memory() 
    if (flag_show_memory) buf_display_memory_usage (false);    

    buf_free (&dd->compute_threads_buf);
    vb_release_vb (&dd->pseudo_vb);

    // free memory allocations that assume subsequent files will have the same number of samples.
    // this is only true if the files are being concatenated
    if (!flag_concat_mode) vb_cleanup_memory(dd->vb_pool); 

    // free dictionary memories that were abandoned when dictionary realloced (they were not freed so that they can remain
    // with their overlaying buffers) 
    buf_free_abandoned_memories();

    free (dd);
}

static void *dispatcher_thread_entry (void *thread_)
{
    Thread *th = (Thread *)thread_;

    th->func (th->vb);
    
    return NULL;
}

VariantBlock *dispatcher_generate_next_vb (Dispatcher dispatcher)
{
    DispatcherData *dd = dispatcher;

    dd->next_vb_i++;

    dd->next_vb = vb_get_vb (dd->vb_pool, dd->vcf_file, dd->z_file, dd->next_vb_i);
    return dd->next_vb;
}

void dispatcher_compute (Dispatcher dispatcher, void (*func)(VariantBlock *))
{
    DispatcherData *dd = dispatcher;
    Thread *th = &dd->compute_threads[dd->next_thread_to_dispatched];

    th->vb = dd->next_vb;
    th->func = func;

    //printf ("dispatching line = %u\n", next_vb->first_line);
    if (dd->max_threads > 1 
#if defined _WIN32 && ! defined _WIN64 // note: _WIN32 is defined for both Windows 32 & 64 bit
        // note: in 32bit Windows, the first variant block is always done on the main thread, so we can measure
        // memory consumption and reduce the number of compute threads and/or num_lines in subsequent vbs
        && dd->next_vb_i > 1
#endif
        ) {

        unsigned err = pthread_create(&th->thread_id, NULL, dispatcher_thread_entry, th);
        ASSERT (!err, "Error: failed to create thread for next_vb_i=%u, err=%u", dd->next_vb->variant_block_i, err);

        dd->next_thread_to_dispatched = (dd->next_thread_to_dispatched + 1) % (dd->max_threads-1);
    }
    else {     // single thread
        func(dd->next_vb);            

#ifdef _WIN32
        if (dd->max_threads > 1) {
            // adjust max_threads and/or num_lines, now that we know how memory a vb consumes
            unsigned vb_memory = buf_vb_memory_consumption (dd->next_vb);
            long long num_threads_that_fit_in_memory = (dd->test_mode ? MAX_32BIT_WINDOWS_MEMORY/2 : MAX_32BIT_WINDOWS_MEMORY) / vb_memory; // in test mode, both zip and piz are running - each side gets half of the memory
            dd->max_threads = MIN (dd->max_threads, num_threads_that_fit_in_memory); // TO DO - play with num_lines to, not just compute threads
#ifdef DEBUG
            char str[30]; 
            printf ("\nvb_memory=%s max_threads=%u\n", buf_human_readable_size (vb_memory, str), dd->max_threads);
#endif
        }
#endif
    }
    dd->next_vb = NULL;
    dd->num_running_compute_threads++;
}

VariantBlock *dispatcher_get_next_processed_vb (Dispatcher dispatcher)
{
    DispatcherData *dd = dispatcher;

    if (dd->max_threads > 1 && !dd->num_running_compute_threads) return NULL; // no running compute threads 

    Thread *th = &dd->compute_threads[dd->next_thread_to_be_joined];

    if (dd->max_threads > 1) 
        // wait for thread to complete (possibly it completed already)
        pthread_join(th->thread_id, NULL);

    dd->processed_vb = th->vb;
    
    memset (th, 0, sizeof(Thread));
    dd->num_running_compute_threads--;
    dd->next_thread_to_be_joined = (dd->next_thread_to_be_joined + 1) % MAX (1, dd->max_threads-1);

    return dd->processed_vb;
}

bool dispatcher_has_free_thread (Dispatcher dispatcher)
{
    DispatcherData *dd = dispatcher;
    return dd->num_running_compute_threads < MAX(1, dd->max_threads-1);
}

VariantBlock *dispatcher_get_next_vb (Dispatcher dispatcher)
{
    DispatcherData *dd = dispatcher;
    return dd->next_vb;
}

VariantBlock *dispatcher_get_pseudo_vb (Dispatcher dispatcher)
{
    DispatcherData *dd = dispatcher;
    return dd->pseudo_vb;
}

static void dispatcher_human_time (unsigned secs, char *str /* out */)
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

static void dispatcher_show_progress (Dispatcher dispatcher, const File *file, long long vcf_data_written_so_far,
                                      uint64_t bytes_compressed /* possibly 0 */)
{
    DispatcherData *dd = dispatcher;

    if (!dd->show_progress) return; 

    struct timespec tb; 
    clock_gettime(CLOCK_REALTIME, &tb); 
    
    int seconds_so_far = ((tb.tv_sec-dd->start_time.tv_sec)*1000 + (tb.tv_nsec-dd->start_time.tv_nsec) / 1000000) / 1000; 

    static double ratio_so_far = 1;

    uint64_t total, sofar;
    if (file->vcf_data_size) { // if we have the VCF data size, go by it 
        total = file->vcf_data_size;
        sofar = vcf_data_written_so_far;
    
    } else if (file->disk_size) {
        total = file->disk_size; // in case of .vcf.gz
        sofar = file->disk_so_far; 

        // sofar can be -1LL bc it seems that gzoffset64 returns that for values over 2GB on Windows -
        // as a work around, use the previous ratio as an estimate
        if (sofar > 9000000000000000000ULL)
            sofar = ratio_so_far * (double)vcf_data_written_so_far * 0.97; // underestimate by bit (97%) - better err on the side of under-estimation
        else
            ratio_so_far = (double)sofar / (double)vcf_data_written_so_far;
    } 
    else // in case of a file over a pipe
        return; // we can't show anything if we don't know the file size

    double percent = MIN (((double)(sofar/100000ULL)*100) / (double)(total/100000ULL), 100.0); // divide by 100000 to avoid number overflows

    // need to update progress indicator, max once a second or if 100% is reached
    if (percent && (dd->last_seconds_so_far < seconds_so_far || percent == 100)) { 

        const char *eraser = "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";

        if (!dispatcher_is_done (dispatcher)) {

            // time remaining
            char time_str[70], progress[100];
            unsigned secs = (100.0 - percent) * ((double)seconds_so_far / (double)percent);
            dispatcher_human_time (secs, time_str);

            sprintf (progress, "%u%% (%s)", (unsigned)percent, time_str);

            // note we have spaces at the end to make sure we erase the previous string, if it is longer than the current one
            fprintf (stderr, "%.*s%s        %.8s", dd->last_len, eraser, progress, eraser);

            dd->last_len = strlen (progress);

        } else if (!dd->test_mode) {
            char time_str[70];
            dispatcher_human_time (seconds_so_far, time_str);

            if (bytes_compressed) {
                if (file->vcf_data_size == file->disk_size) // source file was plain VCF
                    fprintf (stderr, "%.*sDone (%s, compression ratio: %1.1f)           \n", dd->last_len, eraser, time_str, (double)total / (double)bytes_compressed);
                else // source was .vcf.gz
                    fprintf (stderr, "%.*sDone (%s, VCF compression ratio: %1.1f ; ratio vs gzip: %1.1f)\n", 
                             dd->last_len, eraser, time_str, 
                             (double)file->vcf_data_size / (double)bytes_compressed,  // compression vs vcf data size
                             (double)file->disk_size     / (double)bytes_compressed); // compression vs gzipped size
            } else
                fprintf (stderr, "%.*sDone (%s)                         \n", dd->last_len, eraser, time_str);

        } else
            fprintf (stderr, "%.*s", dd->last_len, eraser); // test result comes after
    }

    dd->last_percent = percent;
    dd->last_seconds_so_far = seconds_so_far;
}

void dispatcher_finalize_one_vb (Dispatcher dispatcher, const File *file, long long vcf_data_written_so_far,
                                 uint64_t bytes_compressed /* possibly 0 */)
{
    DispatcherData *dd = dispatcher;

#ifdef DEBUG
    buf_test_overflows(dd->processed_vb);
#endif

    if (flag_show_time) 
        profiler_add (&dd->pseudo_vb->profile, &dd->processed_vb->profile);

    vb_release_vb (&dd->processed_vb); // cleanup vb and get it ready for another usage (without freeing memory)

    dispatcher_show_progress (dispatcher, file, vcf_data_written_so_far, bytes_compressed);
}                           

void dispatcher_input_exhausted (Dispatcher dispatcher)
{
    DispatcherData *dd = dispatcher;
    dd->input_exhausted = true;

    vb_release_vb (&dd->next_vb);
}    

bool dispatcher_is_done (Dispatcher dispatcher)
{
    DispatcherData *dd = dispatcher;

    return dd->input_exhausted && !dd->next_vb && !dd->processed_vb && !dd->num_running_compute_threads;
}

bool dispatcher_is_final_processed_vb (Dispatcher dispatcher)
{
    DispatcherData *dd = dispatcher;

    return dd->input_exhausted && !dd->next_vb && !dd->num_running_compute_threads;
}

bool dispatcher_is_input_exhausted (Dispatcher dispatcher)
{
    DispatcherData *dd = dispatcher;

    return dd->input_exhausted;
}
