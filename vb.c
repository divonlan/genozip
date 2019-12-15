// ------------------------------------------------------------------
//   vb.c
//   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// vb stands for VariantBlock - i.e. one block of VARIANTS_PER_BLOCK data lines in the vcf file

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "vczip.h"

unsigned vb_num_samples_in_sb(const VariantBlock *vb, unsigned sb_i)
{
    return sb_i < vb->num_sample_blocks-1 ? vb->num_samples_per_block 
                                          : global_num_samples % vb->num_samples_per_block;
} 

unsigned vb_num_sections(VariantBlock *vb) 
{
    return 1 + vb->has_genotype_data + (vb->phase_type == PHASE_MIXED_PHASED) + (vb->num_haplotypes_per_line > 0);
}

// cleanup vb and get it ready for another usage (without freeing memory held in the Buffers)
void vb_release_vb (VariantBlock **vb_p) 
{
    if (! *vb_p) return; // nothing to release

    VariantBlock *vb = *vb_p;
    *vb_p = NULL;

    for (unsigned i=0; i < VARIANTS_PER_BLOCK; i++) {
        DataLine *dl = &vb->data_lines[i];
         
        dl->line_i = dl->genotype_data.len = 0;
        dl->phase_type = PHASE_UNKNOWN;
        dl->has_haplotype_data = dl->has_genotype_data = 0;

        buf_free(&dl->line);
        buf_free(&dl->variant_data);
        buf_free(&dl->genotype_data);
        buf_free(&dl->haplotype_data);
        buf_free(&dl->phase_data);
    }

    for (unsigned i=0; i < vb->num_sample_blocks; i++) {
        if (vb->haplotype_sections_data) buf_free(&vb->haplotype_sections_data[i]);
        if (vb->genotype_sections_data)  buf_free(&vb->genotype_sections_data[i]);
        if (vb->phase_sections_data)     buf_free(&vb->phase_sections_data[i]);
    }
            

    vb->num_lines = vb->first_line = vb->variant_block_i = vb->longest_line_genotype_data = 0;
    vb->vcf_data_size = vb->ploidy = vb->num_sample_blocks = vb->num_haplotypes_per_line = vb->last_pos = 0;
    vb->has_genotype_data = vb->has_haplotype_data = vb->ready_to_dispatch = false;
    vb->phase_type = PHASE_UNKNOWN;
    vb->vcf_file = vb->z_file = NULL;

    memset(vb->add_bytes, 0, sizeof(vb->add_bytes));
    memset(vb->vcf_section_bytes, 0, sizeof(vb->vcf_section_bytes));
    memset(vb->z_section_bytes, 0, sizeof(vb->z_section_bytes));

    memset (&vb->profile, 0, sizeof (vb->profile));

    buf_free(&vb->line_variant_data);
    buf_free(&vb->line_genotype_data);
    buf_free(&vb->line_haplotype_data);
    buf_free(&vb->line_phase_data);
    buf_free(&vb->next_gt_in_sample);
    buf_free(&vb->gt_line_lengths);
    buf_free(&vb->genotype_one_section_data);
    buf_free(&vb->genotype_section_lens_buf);

    buf_free(&vb->variant_data_section_data);
    buf_free(&vb->haplotype_permutation_index);
    buf_free(&vb->genotype_sample_block_line_starts_buf);
    buf_free(&vb->genotype_sample_block_line_lengths_buf);
    buf_free(&vb->helper_index_buf);
    buf_free(&vb->vardata_header_buf);
    buf_free(&vb->compressed);
    buf_free(&vb->z_data);
    buf_free(&vb->z_section_headers);
    buf_free(&vb->ht_columns_data);

    for (unsigned i=0; i < NUM_COMPRESS_BUFS; i++)
        buf_free (&vb->compress_bufs[i]);
        
    vb->z_file = vb->vcf_file = NULL; // we're not freeing them, must setting the point to null

    // vb->buffer_list : Note: we DON'T free this because the buffers listed are still available and going to be re-used

    vb->in_use = false; // released the VB back into the pool - it may now be reused
}

VariantBlockPool *vb_construct_pool (unsigned num_vbs, unsigned vb_id_prefix)
{
    VariantBlockPool *pool = (VariantBlockPool *)calloc (1, sizeof (VariantBlockPool) + num_vbs * sizeof (VariantBlock)); // note we can't use Buffer yet, because we don't have VBs yet...
    ASSERT (pool, "Error: failed to calloc pool", "");

    pool->num_vbs      = num_vbs;
    pool->vb_id_prefix = vb_id_prefix;

    return pool;
}

// allocate an unused vb from the pool. seperate pools for compress and decompress
VariantBlock *vb_get_vb(VariantBlockPool *pool, 
                        File *vcf_file, File *z_file,
                        unsigned variant_block_i)
{
    VariantBlock *vb;

    if (!pool) { // should only be used for unit testing - memory-leaks a VB
        vb = calloc (1, sizeof(VariantBlock));
        ASSERT (vb, "Error: failed to calloc vb", "");
    }
    else {
        // see if there's a VB avaiable for recycling
        unsigned vb_i; for (vb_i=0; vb_i < pool->num_vbs; vb_i++) {
        
            vb = &pool->vb[vb_i];
        
            if (!vb->in_use) {
                vb->id = pool->vb_id_prefix + vb_i;
                break;
            }
        }

        ASSERT (vb_i < pool->num_vbs, "Error: VB pool vb_id_prefix=%u maxed out", pool->vb_id_prefix)
    }

    vb->in_use          = true;
    vb->variant_block_i = variant_block_i;
    vb->vcf_file        = vcf_file;
    vb->z_file          = z_file;

    return vb;
}

static void vb_human_time (unsigned secs, char *str /* out */)
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


void vb_show_progress(double *last_percent, const File *file, long long vcf_data_written_so_far,
                      long long bytes_compressed, // may be 0
                      unsigned seconds_so_far, bool done,
                      bool test_mode)
{
    static double ratio_so_far = 1;

#ifndef PROFILER

    long long total, sofar;
    if (file->vcf_data_size) { // if we have the VCF data size, go by it 
        total = file->vcf_data_size;
        sofar = vcf_data_written_so_far;
    
    } else if (file->disk_size) {
        total = file->disk_size; // in case of .vcf.gz
        
        sofar = file->disk_so_far; 

        // sofar can be -1LL bc it seems that gzoffset64 returns that for values over 2GB on Windows -
        // as a work around, use the previous ratio as an estimate
        if (sofar > 9000000000000000000ULL)
            sofar = ratio_so_far * vcf_data_written_so_far;
        else
            ratio_so_far = (double)sofar / (double)vcf_data_written_so_far;
    } 
    else // in case of a file over a pipe
        return; // we can't show anything if we don't know the file size

    static unsigned last_len;
    static unsigned last_seconds_so_far;
    
    if (! *last_percent) { // we need to re-initialize these because multiple files might be processed in a loop
        last_len = 2;  // initial string is "0%". 
        last_seconds_so_far = 0;
    }

    double percent = MIN (((double)(sofar/100000ULL)*100) / (double)(total/100000ULL), 100.0); // divide by 100000 to avoid number overflows
    
    // need to update progress indicator, max once a second or if 100% is reached
    if (percent && (last_seconds_so_far < seconds_so_far || percent == 100)) { 

        const char *eraser = "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";

        if (!done) {

            // time remaining
            char time_str[70], progress[70];
            unsigned secs = (100.0 - percent) * ((double)seconds_so_far / (double)percent);
            vb_human_time (secs, time_str);

            sprintf (progress, "%u%% (%s)", (unsigned)percent, time_str);

            // note we have spaces at the end to make sure we erase the previous string, if it is longer than the current one
            fprintf (stderr, "%.*s%s        %.8s", last_len, eraser, progress, eraser);

            last_len = strlen (progress);

        } else if (!test_mode) {
            char time_str[70];
            vb_human_time (seconds_so_far, time_str);

            if (bytes_compressed) {
                if (file->vcf_data_size == file->disk_size) // source file was plain VCF
                    fprintf (stderr, "%.*sDone (%s, compression ratio: %1.1f)           \n", last_len, eraser, time_str, (double)total / (double)bytes_compressed);
                else // source was .vcf.gz
                    fprintf (stderr, "%.*sDone (%s, VCF compression ratio: %1.1f ; ratio vs gzip: %1.1f)\n", 
                             last_len, eraser, time_str, 
                             (double)file->vcf_data_size / (double)bytes_compressed,  // compression vs vcf data size
                             (double)file->disk_size     / (double)bytes_compressed); // compression vs gzipped size
            } else
                fprintf (stderr, "%.*sDone (%s)                         \n", last_len, eraser, time_str);

        } else
            fprintf (stderr, "%.*s", last_len, eraser); // test result comes after
    }

    *last_percent = percent;
    last_seconds_so_far = seconds_so_far;
#endif
}

void vb_adjust_max_threads_by_vb_size(const VariantBlock *vb, unsigned *max_threads, bool test_mode)
{
#if __WIN32__
    // adjust max_threads and/or num_lines, now that we know how memory a vb consumes
    unsigned vb_memory = buf_vb_memory_consumption (vb);
    long long num_threads_that_fit_in_memory = (test_mode ? MAX_MEMORY/2 : MAX_MEMORY) / vb_memory; // in test mode, both zip and piz are running - each side gets half of the memory

    *max_threads = MIN (*max_threads, num_threads_that_fit_in_memory); // TO DO - play with num_lines to, not just compute threads
#ifdef DEBUG
    char str[30]; 
    printf ("\nvb_memory=%s max_threads=%u\n", buf_human_readable_size (vb_memory, str), *max_threads);
#endif
#endif
}