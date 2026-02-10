// ------------------------------------------------------------------
//   zriter.c
//   Copyright (C) 2023-2026 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <errno.h>
#include "genozip.h"
#include "buffer.h"
#include "file.h"
#include "vblock.h"
#include "profiler.h"
#include "zriter.h"
#include "arch.h"
#include "tar.h"
#include "threads.h"

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

static void zriter_write_do (BufferP data, BufferP section_list, int64_t offset_in_z_file)
{
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

    zriter_write_do (data, section_list, offset_in_z_file);

    COPY_TIMER_EVB (zriter_write);
}
