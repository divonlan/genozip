// ------------------------------------------------------------------
//   zriter.c
//   Copyright (C) 2023-2025 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <errno.h>
#include <pthread.h>
#include "genozip.h"
#include "buffer.h"
#include "file.h"
#include "vblock.h"
#include "profiler.h"
#include "threads.h"
#include "zriter.h"
#include "arch.h"
#include "tar.h"

static uint32_t zriter_thread_num = 0; // for --show-threads

typedef struct {
    Buffer data; // note: "pointer" points to section_list buffer, and then set to NULL when thread has completed
    Buffer section_list;
    pthread_t thread_id;
    bool completed;
} ZriterThread, *ZriterThreadP;

static int64_t zriter_tell (void)
{
    int64_t offset = ftello64 ((FILE *)z_file->file);
    ASSERT (offset >= 0 , "ftello64 failed for %s (FILE*=%p remote=%s redirected=%s): %s", 
            z_file->name, z_file->file, TF(z_file->is_remote), TF(z_file->redirected), strerror (errno));

    // in a z_file that is being tarred, update the offset to the beginning of the file data in the tar file
    offset -= tar_file_offset(); // 0 if not using tar

    return offset;    
}

void zriter_flush (void)
{
    if (flag.zip_no_z_file) return;

    ASSERT (!fflush ((FILE*)z_file->file), 
            "fflush to %s on a %s filesystem failed: %s", z_name, arch_get_z_filesystem().s, strerror (errno));
}

void zriter_wait_for_bg_writing (void)
{
    if (flag.zip_no_z_file) return;
    
    for_buf (ZriterThreadP, zt_p, z_file->zriter_threads) {
        if (flag_show_threads) iprintf ("zriter: JOINING: thread_id=%"PRIu64"\n", (uint64_t)(*zt_p)->thread_id);
        PTHREAD_JOIN ((*zt_p)->thread_id, "zriter_thread_entry"); // blocking
        if (flag_show_threads) iprintf ("zriter: JOINED: thread_id=%"PRIu64"\n", (uint64_t)(*zt_p)->thread_id);

        buf_destroy ((*zt_p)->data);
        buf_destroy ((*zt_p)->section_list);
        FREE (*zt_p);
    }

    z_file->zriter_threads.len = 0;

    zriter_flush(); // important - eg flush pair-1 to disk before pair-2 re-reads via a different FILE 
}

static void *zriter_thread_entry (void *zt_)
{
    threads_set_zriter_thread();

    mutex_lock (z_file->zriter_mutex); // only one thread can write a time

    START_TIMER;

    ZriterThreadP zt = (ZriterThreadP)zt_;
    ASSERTNOTZERO (zt->section_list.len);

    // note: no fflush() needed between subsequent background threads bc they are writting sequentially.
    int64_t bytes_written = fwrite (zt->data.data, 1, zt->data.len, (FILE *)z_file->file); 
    
    // error if failed to write to file
    ASSERT (bytes_written == zt->data.len, "wrote only %"PRId64" of the expected %"PRId64" bytes to %s on a %s filesystem: %s", 
            bytes_written, zt->data.len, z_file->name, arch_get_z_filesystem().s, strerror(errno));

    sections_list_concat (&zt->section_list); // note: must be before incrementing disk_so_far

    z_file->disk_so_far += zt->data.len; // length of GENOZIP data writen to disk (protected by zriter_mutex)

    // sanity
    uint64_t actual_disk_so_far = zriter_tell();
    ASSERT (actual_disk_so_far == z_file->disk_so_far, "Expecting actual_disk_so_far=%"PRIu64" == z_file->disk_so_far=%"PRIu64,
            actual_disk_so_far, z_file->disk_so_far);

    COPY_TIMER_EVB (write_bg); // "write" profiler resource protected by zriter_mutex

    store_release (zt->completed, true); // signal that thread is ready to be joined

    mutex_unlock (z_file->zriter_mutex);

    threads_unset_zriter_thread();
    return NULL;
}

static void zriter_write_background (BufferP data, BufferP section_list)
{
    ASSERTNOTNULL (section_list); // we must have section_list (even if empty) to write in background

    // join all threads that have completed
    for (int32_t i=0; i < z_file->zriter_threads.len32; i++) { // note: can't use for_buf, bc data.len decreases in the loop
        ZriterThreadP *zt_p = B(ZriterThreadP, z_file->zriter_threads, i);

        if (load_acquire ((*zt_p)->completed)) { // completed
            // free resources

            if (flag_show_threads) iprintf ("zriter: JOINING: thread_id=%"PRIu64"\n", (uint64_t)(*zt_p)->thread_id);
            PTHREAD_JOIN ((*zt_p)->thread_id, "zriter_thread_entry");
            if (flag_show_threads) iprintf ("zriter: JOINED: thread_id=%"PRIu64"\n", (uint64_t)(*zt_p)->thread_id);

            buf_destroy ((*zt_p)->data);
            buf_destroy ((*zt_p)->section_list);
            FREE (*zt_p);

            // remove entry from array
            memmove (zt_p, zt_p+1, (z_file->zriter_threads.len32 - i - 1) * sizeof (ZriterThreadP));
            z_file->zriter_threads.len--;
            i--;
        }
    }

    // zriter_threads is an array of pointers (ZriterThreadP) - so the pointers themselves survive array reallocs
    buf_alloc (evb, &z_file->zriter_threads, 1, 10, ZriterThreadP, 2, "z_file->zriter_threads");
    ZriterThreadP *zt_p = &BNXT(ZriterThreadP, z_file->zriter_threads);
    *zt_p = CALLOC (sizeof(ZriterThread));

    // grab buffers, to save from destruction in vb_release_vb 
    buf_grab (evb, (*zt_p)->data, data->name, *data);
    buf_grab (evb, (*zt_p)->section_list, section_list->name, *section_list);

    // make sure changes *zt_p will be visible to the new background thread (does pthread_create already do this?)
    __atomic_thread_fence (__ATOMIC_RELEASE); 

    // note: we pass "data" to the thread (with a pointer to section_list) instead of zt, bc z_file->zriter_threads might be realloced
    unsigned err = pthread_create (&(*zt_p)->thread_id, NULL, zriter_thread_entry, *zt_p);
    ASSERT (!err, "failed to create thread zriter thread: %s", strerror(err));

    if (flag_show_threads) iprintf ("zriter: CREATE: thread_id=%"PRIu64" zriter_num=%u data_len=%u sections=%u\n", 
                                    (uint64_t)(*zt_p)->thread_id, zriter_thread_num++, (*zt_p)->data.len32, (*zt_p)->section_list.len32);
}

static void zriter_write_foreground (BufferP data, BufferP section_list, int64_t offset_in_z_file)
{
    mutex_lock (z_file->zriter_mutex);

    START_TIMER; // not including mutex wait time

    if (offset_in_z_file != -1) {
        zriter_flush(); // just in case...
        file_seek (z_file, offset_in_z_file, SEEK_SET, WRITE, HARD_FAIL); 
    }
    
    int64_t bytes_written = fwrite (data->data, 1, data->len, (FILE *)z_file->file); // use fwrite - let libc manage write buffers for us

    // error if failed to write to file
    ASSERT (bytes_written == data->len, "wrote only %"PRId64" of the expected %"PRId64" bytes to %s on a %s filesystem: %s", 
            bytes_written, data->len, z_file->name, arch_get_z_filesystem().s, strerror(errno));

    // writing to an offset - return to the end of the file
    if (offset_in_z_file != -1) {
        zriter_flush(); // its not clear why, but without this fflush the bytes immediately after the first header get corrupted (at least on Windows with gcc)
        file_seek (z_file, 0, SEEK_END, WRITE, HARD_FAIL); 
    }

    // appending to end of file
    else {
        ASSERTNOTNULL (section_list);
        sections_list_concat (section_list); // note: must be before incrementing disk_so_far

        z_file->disk_so_far += data->len;   // length of GENOZIP data writen to disk (pr)

        // sanity
        uint64_t actual_disk_so_far = zriter_tell();
        ASSERT (actual_disk_so_far == z_file->disk_so_far, "Expecting actual_disk_so_far=%"PRIu64" == z_file->disk_so_far=%"PRIu64,
                actual_disk_so_far, z_file->disk_so_far);
    }

    COPY_TIMER_EVB (write_fg); // "write" profiler resource protected by zriter_mutex

    mutex_unlock (z_file->zriter_mutex);
}

void zriter_write (BufferP data,
                   BufferP section_list,     // section list to append to z_file->section_list (non-NULL iff offset_in_z_file==-1)
                   int64_t offset_in_z_file, // -1 means append to end of file
                   bool background)          // write in a separate thread to enable concurrent reading (only one thread at a time writes)
{
    START_TIMER;

    ASSERTMAINTHREAD;

    if (!data->len || flag.zip_no_z_file) return; // nothing to do

    ASSERTNOTNULL (z_file);
    ASSERTNOTNULL (z_file->file);
    ASSERTNOTNULL (data->data);

    // cases where we override a caller request for background writing
    if (data->shared || (section_list && section_list->shared) || // shared buffers cannot be grabbed (eg txt_data is shared in SAM PRIM)
        flag.no_zriter) // user requested or as set in flags_update(); 
        background = false; 

    // case: foreground write
    if (!background) 
        zriter_write_foreground (data, section_list, offset_in_z_file);

    // case: background appending
    else {
        ASSERT0 (offset_in_z_file == -1, "cannot write to an offset in the background");
        zriter_write_background (data, section_list);
    }

    z_file->zriter_last_was_fg = !background;

    COPY_TIMER_EVB (zriter_write);
}
