// ------------------------------------------------------------------
//   zip.c
//   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <inttypes.h>

#include "vczip.h"

// read entire variant block to memory. this is called from the dispatcher thread
static void zip_read_variant_block (File *vcf_file,
                                    unsigned *line_i,   // in/out next line to be read
                                    Buffer *first_data_line,    // first line might be supplied by caller 
                                    VariantBlock *vb)
{
    START_TIMER;

    unsigned first_line= *line_i;

    unsigned vb_line_i;
    vb->vcf_data_size = 0; // size of variant block as it appears in the source file
    for (vb_line_i=0; vb_line_i < VARIANTS_PER_BLOCK; vb_line_i++) 
    {
        DataLine *dl = &vb->data_lines[vb_line_i];

        if (vb_line_i > 0 || !first_data_line) {// first line might be supplied by caller

            // allocate line Buffer and read line from file 
            bool success = vcffile_get_line (vb, first_line + vb_line_i, &dl->line); 
            if (!success) break; // no more lines - we're done
        }
        else {
            buf_copy (vb, &dl->line, first_data_line, 0, 0);
            dl->line.len = first_data_line->len;
            buf_free (first_data_line); 
        }
        dl->line_i = first_line + vb_line_i;

        (*line_i)++;
        
        // count bytes in the source file, for sanity check after we reconstruct
        vb->vcf_data_size += dl->line.len;
    }

    vb->num_lines = vb_line_i;

    COPY_TIMER (vb->profile.read);
}

// concatenated variant data (up to and including the FORMAT field) into a single string
// with \n separating variants
static void zip_generate_variant_data_section (VariantBlock *vb)
{
    START_TIMER;

    for (unsigned vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) 
        vb->variant_data_section_data.len += vb->data_lines[vb_line_i].variant_data.len;

    buf_alloc (vb, &vb->variant_data_section_data, vb->variant_data_section_data.len + 1 /* for the \0 */, 1.1, "variant_data_section_data", vb->first_line);

    // concatenate variant data of all lines
    unsigned offset, vb_line_i; 
    for (offset=0, vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {

        unsigned len = vb->data_lines[vb_line_i].variant_data.len;

        memcpy (vb->variant_data_section_data.data + offset, vb->data_lines[vb_line_i].variant_data.data, len);

        offset += len;
    }
    vb->variant_data_section_data.data[offset] = '\0';

    COPY_TIMER (vb->profile.zip_generate_variant_data_section)
}

static unsigned zip_get_genotype_vb_start_len(VariantBlock *vb)
{
    buf_alloc (vb, &vb->genotype_section_lens_buf, sizeof(unsigned) * vb->num_sample_blocks, 1, "section_lens_buf", 0);
    unsigned *section_lens = (unsigned *)vb->genotype_section_lens_buf.data;
    memset (section_lens, 0, vb->genotype_section_lens_buf.size); // initialize

    // offsets into genotype data of individual lines
    buf_alloc (vb, &vb->genotype_sample_block_line_starts_buf, 
               vb->num_lines * vb->num_sample_blocks * sizeof(char*), 
               0, "genotype_sample_block_line_starts_buf", vb->first_line);
    char    **genotype_sample_block_line_starts  = (char**)vb->genotype_sample_block_line_starts_buf.data; 
    
    // length of a single line in a sample block
    buf_alloc (vb, &vb->genotype_sample_block_line_lengths_buf, 
               vb->num_lines * vb->num_sample_blocks * sizeof(unsigned), 
               0, "genotype_sample_block_line_lengths_buf", vb->first_line);
    unsigned *genotype_sample_block_line_lengths = (unsigned *)vb->genotype_sample_block_line_lengths_buf.data; 
    
    // calculate offsets and lengths of genotype data of each sample block
    for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {

        char *next = vb->data_lines[line_i].genotype_data.data;
        char *after = vb->data_lines[line_i].genotype_data.data + vb->data_lines[line_i].genotype_data.len;

        for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

            unsigned num_samples_in_sb = vb_num_samples_in_sb (vb, sb_i);

            char *this_start = next;
            genotype_sample_block_line_starts[line_i * vb->num_sample_blocks + sb_i] = next;

            // find the end of this sample block genotype data (which consists of the gt type of all samples in the sample block)
            unsigned sample_i ; for (sample_i = 0; sample_i < num_samples_in_sb && next < after; )
                if (*(next++) == '\t') sample_i++;

            ASSERT ((next == after || sb_i < vb->num_sample_blocks-1) && sample_i == num_samples_in_sb, 
                    "Error: corrupt genotype line data: sb_i=%u line_i(0-based)=%u", sb_i, line_i);

            genotype_sample_block_line_lengths[line_i * vb->num_sample_blocks + sb_i] = next - this_start;
            
            section_lens[sb_i] += next - this_start;
        }
    }

    // find maximum section length and return it
    unsigned max_len = 0;
    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) 
        if (section_lens[sb_i] > max_len)
            max_len = section_lens[sb_i];

    return max_len;
}

// split genotype data to sample groups, within a sample group genotypes are separated by a tab
static void zip_generate_genotype_one_section (VariantBlock *vb, unsigned sb_i)
{
    START_TIMER;

    // build sample block genetype data
    char *dst_next = vb->genotype_one_section_data.data;
    
    const char *first_gt_in_sb = ((char**)vb->genotype_sample_block_line_starts_buf.data)[sb_i];
    unsigned first_gt_in_sb_len = strchr (first_gt_in_sb, '\t') - first_gt_in_sb + 1; // length including the \t

    // move the GT items from the line data to the permuted data - with each 
    // sample block of gt data containing the data in TRANSPOSED order - i.e. first
    // the gt data for all the variants for sample 1, then all of samples 2 etc.
    unsigned num_samples_in_sb = vb_num_samples_in_sb (vb, sb_i);
    for (unsigned sample_i=0; sample_i < num_samples_in_sb; sample_i++) {
        for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {

            char **src_next_in_this_sb_p = &((char**)vb->genotype_sample_block_line_starts_buf.data)[line_i * vb->num_sample_blocks + sb_i];

            // copy string up to and including the \t, and push the index ahead
            unsigned len = strcpy_tab (dst_next, *src_next_in_this_sb_p); 
            *src_next_in_this_sb_p += len;
            dst_next               += len;
            vb->genotype_one_section_data.len += len;
        }
    }

    COPY_TIMER (vb->profile.zip_generate_genotype_sections)
}

// split phase data to sample groups, in each group is a string of / | or -
static void zip_generate_phase_sections (VariantBlock *vb)
{   
    START_TIMER;

    // we allocate memory for the Buffer array only once the first time this VariantBlock
    // is used. Subsequent blocks reusing the memory will have the same number of samples (by VCF spec)
    if (!vb->phase_sections_data) 
        vb->phase_sections_data = calloc (vb->num_sample_blocks, sizeof(Buffer)); // allocate once, never free

    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned num_samples_in_sb = vb_num_samples_in_sb (vb, sb_i); 

        // allocate memory for phase data for each sample block - one character per sample
        buf_alloc (vb, &vb->phase_sections_data[sb_i], vb->num_lines * num_samples_in_sb, 
                   0, "phase_sections_data[sb_i]", vb->first_line);

        // build sample block genetype data
        char *next = vb->phase_sections_data[sb_i].data;
        
        for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {
            
            DataLine *dl = &vb->data_lines[line_i];
            if (dl->phase_type == PHASE_MIXED_PHASED) 
                memcpy (next, &dl->phase_data.data[sb_i * num_samples_in_sb], num_samples_in_sb);
            else
                memset (next, (char)dl->phase_type, num_samples_in_sb);

            next += num_samples_in_sb;
        }
        vb->phase_sections_data[sb_i].len = num_samples_in_sb * vb->num_lines;
    }
          
    // add back the phase data bytes that weren't actually "saved"
    vb->add_bytes[SEC_PHASE_DATA] += global_num_samples * vb->num_lines;

    COPY_TIMER (vb->profile.zip_generate_phase_sections)
}

typedef struct {
    int num_alt_alleles;
    unsigned index_in_original_line;
    unsigned index_in_sorted_line;
} HaploTypeSortHelperIndex;

static int sort_by_alt_allele_comparator(const void *p, const void *q)  
{ 
    int l = ((HaploTypeSortHelperIndex *)p)->num_alt_alleles; 
    int r = ((HaploTypeSortHelperIndex *)q)->num_alt_alleles;  
    return (l - r); 
}

static int sort_by_original_index_comparator(const void *p, const void *q)  
{ 
    int l = ((HaploTypeSortHelperIndex *)p)->index_in_original_line; 
    int r = ((HaploTypeSortHelperIndex *)q)->index_in_original_line;  
    return (l - r); 
}

// sort haplogroups by alt allele count within the variant group, create an index for it, and split
// it to sample groups. for each sample a haplotype is just a string of 1 and 0 etc (could be other alleles too)
static void zip_generate_haplotype_sections (VariantBlock *vb)
{
    START_TIMER;

    // we allocate memory for the Buffer array only once the first time this VariantBlock
    // is used. Subsequent blocks reusing the memory will have the same number of samples (by VCF spec)
    if (!vb->haplotype_sections_data) 
        vb->haplotype_sections_data = calloc (vb->num_sample_blocks, sizeof(Buffer)); // allocate once, never free

    // create a permutation index for the whole variant block, and permuted haplotypes for each sample block        
    buf_alloc (vb, &vb->haplotype_permutation_index, vb->num_haplotypes_per_line * sizeof(unsigned), 
               0, "haplotype_permutation_index", vb->first_line);

    buf_alloc (vb, &vb->helper_index_buf, vb->num_haplotypes_per_line * sizeof(HaploTypeSortHelperIndex), 0,
               "helper_index_buf", vb->first_line);
    memset (vb->helper_index_buf.data, 0, vb->helper_index_buf.size);
    HaploTypeSortHelperIndex *helper_index = (HaploTypeSortHelperIndex *)vb->helper_index_buf.data;

    { 
        START_TIMER; // this loop, testing in 1KGP data, consumes 10% of all compute time
        // build index array 
        for (unsigned ht_i=0; ht_i < vb->num_haplotypes_per_line; ht_i++) {

            helper_index[ht_i].index_in_original_line = ht_i;

            for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {

                // we count as alt alleles : 1 - 99 (ascii 49 to 147)
                //             ref alleles : 0 . (unknown) - (missing) * (ploidy padding)
                char one_ht = vb->data_lines[line_i].haplotype_data.data[ht_i];
                if (one_ht >= '1')
                    helper_index[ht_i].num_alt_alleles++;
            }
        }
        COPY_TIMER (vb->profile.count_alt_alleles);
    }

    // sort the helper index array by number of alt alleles
    qsort (helper_index, vb->num_haplotypes_per_line, sizeof (HaploTypeSortHelperIndex), sort_by_alt_allele_comparator);

    // now build per-sample-block haplotype array, picking haplotypes by the order of the helper index array
    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned num_haplotypes_in_sample_block = 
            vb->ploidy * vb_num_samples_in_sb (vb, sb_i); 

        // allocate memory for haplotype data for each sample block - one character per haplotype
        buf_alloc (vb, &vb->haplotype_sections_data[sb_i], vb->num_lines * num_haplotypes_in_sample_block, 
                   0, "haplotype_sections_data[sb_i]", vb->first_line);

        // build sample block haplptype data - 
        // -- using the helper index to access the haplotypes in sorted order
        // -- transposing the array
        char *next = vb->haplotype_sections_data[sb_i].data;
        
        {   // this loop, tested with 1KGP data, takes up to 1/5 of total compute time, so its highly optimized
            START_TIMER;
            unsigned helper_index_sb_i = sb_i * vb->num_samples_per_block * vb->ploidy;

            for (unsigned ht_i=0; ht_i < num_haplotypes_in_sample_block; ht_i++) {
                
                unsigned haplotype_data_char_i = helper_index[helper_index_sb_i + ht_i].index_in_original_line;
                char **ht_data_ptr = &vb->data_lines[0].haplotype_data.data; // this pointer moves sizeof(DataLine) bytes each iteration - i.e. to the exact same field in the next line

                for (unsigned line_i=0; line_i < vb->num_lines; line_i++, ht_data_ptr += sizeof(DataLine)/sizeof(char*)) 
                    *(next++) = (*ht_data_ptr)[haplotype_data_char_i];
            }
            COPY_TIMER (vb->profile.sample_haplotype_data);
        }
        vb->haplotype_sections_data[sb_i].len = num_haplotypes_in_sample_block * vb->num_lines;
    }

    // final step - build the reverse index that will allow access by the original index to the sorted array
    // this will be included in the dv file
    for (unsigned ht_i=0; ht_i < vb->num_haplotypes_per_line ; ht_i++)
        helper_index[ht_i].index_in_sorted_line = ht_i;

    // sort array back to its original order
    qsort (helper_index, vb->num_haplotypes_per_line, sizeof (HaploTypeSortHelperIndex), sort_by_original_index_comparator);

    // construct file index
    unsigned *hp_index = (unsigned *)vb->haplotype_permutation_index.data;
    for (unsigned ht_i=0; ht_i < vb->num_haplotypes_per_line ; ht_i++)
        hp_index[ht_i] = helper_index[ht_i].index_in_sorted_line;

    buf_free (&vb->helper_index_buf);

    COPY_TIMER (vb->profile.zip_generate_haplotype_sections);
}

// this function receives all lines of a variant block and processes them
// in memory to the compressed format. This thread then terminates the I/O thread writes the output.
static void *zip_compress_variant_block (VariantBlock *vb)
{ 
    START_TIMER;

    // allocate memory for the final compressed data of this vb. allocate 20% of the 
    // vb size on the original file - this is normally enough. if not, we will realloc downstream
    buf_alloc (vb, &vb->z_data, vb->vcf_data_size / 5, 1.2, "z_data", 0);

    // initalize variant block data (everything else is initialzed to 0 via calloc)
    vb->phase_type = PHASE_UNKNOWN;  // phase type of this block
    vb->num_samples_per_block = SAMPLES_PER_BLOCK;
    vb->num_sample_blocks = ceil((float)global_num_samples / (float)vb->num_samples_per_block);

    // if testing, make a copy of the original lines read from the file, to compare to the result of vunblocking later
    Buffer *lines_orig = NULL;
    unsigned max_genotype_section_len;

    // split each lines in this variant block to its components
    segregate_all_data_lines (vb, lines_orig);

    // if block has haplotypes - handle them now
    if (vb->has_haplotype_data)
        zip_generate_haplotype_sections (vb); 

    // if block has genetype data - calculate starts, lengths and allocate memory
    if (vb->has_genotype_data) 
        max_genotype_section_len = zip_get_genotype_vb_start_len (vb);
    
    // if block has phase data - handle it
    if (vb->phase_type == PHASE_MIXED_PHASED) 
        zip_generate_phase_sections (vb);

    // permute & write variant data for the variant block
    zip_generate_variant_data_section (vb);
    zfile_compress_variant_data (vb);

    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {
        
        if (vb->has_genotype_data) {

            buf_alloc (vb, &vb->genotype_one_section_data, max_genotype_section_len, 1.05, "genotype_one_section_data", sb_i);

            // we compress each section at a time to save memory
            zip_generate_genotype_one_section (vb, sb_i); 

            zfile_compress_section_data (vb, SEC_GENOTYPE_DATA, &vb->genotype_one_section_data);

            buf_free (&vb->genotype_one_section_data);
        }

        if (vb->phase_type == PHASE_MIXED_PHASED)
            zfile_compress_section_data (vb, SEC_PHASE_DATA, &vb->phase_sections_data[sb_i]);

        if (vb->has_haplotype_data)
            zfile_compress_section_data (vb, SEC_HAPLOTYPE_DATA, &vb->haplotype_sections_data[sb_i]);
    }

    // update dv file size
    ((SectionHeaderVariantData *)vb->z_data.data)->z_data_bytes = vb->z_data.len;

    // updates stats for --showcontent - reduce by the bytes added by us over the bytes in the original file (might be a negative value if we saved bytes rather than added)
    for (SectionType sec_i=0; sec_i < NUM_SEC_TYPES; sec_i++)
        vb->vcf_section_bytes[sec_i] -= vb->add_bytes[sec_i];

    COPY_TIMER (vb->profile.zip_compress_variant_block);

    return NULL;
}

// this is the compute thread entry 
void *zip_compress_variant_block_thread_entry (void *vb_)
{
    zip_compress_variant_block ((VariantBlock*)vb_);
    
    return NULL;
}

// this is the main dispatcher function. It first processes the VCF header, then proceeds to read 
// a variant block from the input file and send it off to a thread for computation. When the thread
// completes, this function proceeds to write the output to the output file. It can dispatch
// several threads in parallel.
void zip_dispatcher (char *vcf_basename, File *vcf_file, 
                     File *z_file, bool concat_mode, bool test_mode, unsigned max_threads)
{
    START_WALLCLOCK;

    VariantBlockPool *vb_pool = vb_construct_pool (max_threads+1 /* +1 for pseudo-vb */, 100 /* for zip */);

    VariantBlock *pseudo_vb = vb_get_vb (vb_pool, vcf_file, z_file, 0);

    unsigned line_i = 0; // last line read (first line in file = 1, consistent with script line numbers)

    // first compress the VCF header
    Buffer *first_data_line = NULL; // contains a value only for first variant block, otherwise empty. 
    
    // read the vcf header, assign the global variables, and write the compressed header to the VCZ file
    bool success = vcf_header_vcf_to_vcz (pseudo_vb, concat_mode, &line_i, &first_data_line);
    if (!success) goto finish;

    if (!first_data_line) goto finish; // VCF file has only a header or is an empty file - no data - we're done

    typedef struct {
        pthread_t thread_id;
        VariantBlock *vb;
    } Thread;
    
    static Buffer compute_threads_buf = EMPTY_BUFFER;
    buf_alloc (pseudo_vb, &compute_threads_buf, sizeof(Thread) * (max_threads-1), 1, "compute_threads_buf", 0);
    Thread *compute_threads = (Thread *)compute_threads_buf.data;

    // Note: we have one VB for every compute thread, and then one VB who is being populated
    // for dispatching for the next thread - total MAX_COMPUTE_THREADS+1

    VariantBlock *next_vb = NULL; // the next variant block to be sent to a compute thread

    unsigned next_thread_to_dispatched   = 0;
    unsigned next_thread_to_be_joined    = 0;

    unsigned num_running_compute_threads = 0;
    unsigned input_exhausted = false;

    // this is the dispatcher loop. In each iteration, it can do one of 3 things, in this order of priority:
    // 1. In there is a new variant block avaialble, and a compute thread available to take it - dispatch it
    // 2. If there is no new variant block available, but input is not exhausted yet - read one
    // 3. Wait for the first thread (by sequential order) to complete the compute and output the results

    unsigned variant_block_i = 0;
    bool show_progress = vcf_basename && !!isatty(2);

#ifndef PROFILER
    if (show_progress) 
        fprintf (stderr, "%svczip %s: 0%%", test_mode ? "testing " : "", vcf_basename);
#endif

    double last_pc_complete = 0.0;

    do {
        // PRIORITY 1: is there a block available and a compute thread available? in that case dispatch it
       if (next_vb && next_vb->ready_to_dispatch && num_running_compute_threads < max_threads - 1) {

            compute_threads[next_thread_to_dispatched].vb = next_vb;

            // note: the first variant block is always done on the main thread, so we can measure
            // memory consumption and reduce the number of compute threads and/or num_lines in subsequent vbs
            if (max_threads > 1 
#if __WIN32__
                && variant_block_i > 1
#endif
               ) {
                unsigned err = pthread_create(&compute_threads[next_thread_to_dispatched].thread_id, NULL, 
                                              zip_compress_variant_block_thread_entry, next_vb);
                ASSERT (!err, "Error: failed to create thread for processing variant block that started on line %u in the VCF file, err=%u", next_vb->first_line, err);
            }
            else {     // single thread
                zip_compress_variant_block(next_vb);            

                if (max_threads > 1)
                    vb_adjust_max_threads_by_vb_size(next_vb, &max_threads, test_mode);
            }

            next_vb = NULL;
            num_running_compute_threads++;

            next_thread_to_dispatched = (next_thread_to_dispatched + 1) % (max_threads-1);

            continue;
        }

        // PRIORITY 2: If there is no new variant block available, but input is not exhausted yet - read one
        if (!next_vb && !input_exhausted) {

            variant_block_i++;

            next_vb = vb_get_vb (vb_pool, vcf_file, z_file, variant_block_i);
            next_vb->first_line = line_i;

            zip_read_variant_block (vcf_file, &line_i, first_data_line, next_vb);
            first_data_line = NULL;

            // input exhausted, and there are no lines - this can happen if the previous variant block was exactly VARIANTS_PER_BLOCK
            if (!next_vb->num_lines) {
                input_exhausted = true;
                vb_release_vb (&next_vb);
                continue;
            }

            // this is the last block, and input is exhausted. this block contains some lines.
            if (next_vb->num_lines < VARIANTS_PER_BLOCK)
                input_exhausted = true;

            next_vb->ready_to_dispatch = true;

            continue;
        }

        // PRIORITY 3: Wait for the first thread (by sequential order) to complete the compute and output the results
        // dispatch computing the compressed to a separate compute thread
        if (num_running_compute_threads) {

            Thread *th = &compute_threads[next_thread_to_be_joined];

            if (max_threads > 1)
                // wait for thread to complete (possibly it completed already)
                pthread_join(th->thread_id, NULL);
      
            START_TIMER;
            
            fwrite (th->vb->z_data.data, 1, th->vb->z_data.len, z_file->file);

            COPY_TIMER (th->vb->profile.write);

            z_file->disk_so_far     += (long long)th->vb->z_data.len;
            z_file->vcf_data_so_far += (long long)th->vb->vcf_data_size;
            z_file->num_lines       += (long long)th->vb->num_lines;

            // update section stats
            for (SectionType sec_i=1; sec_i < NUM_SEC_TYPES; sec_i++) {
                z_file->section_bytes[sec_i]   += th->vb->z_section_bytes[sec_i];
                vcf_file->section_bytes[sec_i] += th->vb->vcf_section_bytes[sec_i];
            }

#ifdef DEBUG
            buf_test_overflows(th->vb);
#endif
            th->thread_id = 0;
            num_running_compute_threads--;
            next_thread_to_be_joined = (next_thread_to_be_joined + 1) % (max_threads - 1);

            bool done = input_exhausted && !num_running_compute_threads && !next_vb;

#ifdef PROFILER
            profiler_add (&pseudo_vb->profile, &th->vb->profile);

            fprintf (stderr, "COMPLETED     thread #%u block #%u line %u %s\n", 
                     next_thread_to_be_joined, th->vb->variant_block_i, th->vb->first_line, profiler_print_short (&th->vb->profile));
#else   
            if (show_progress) {
                int wallclock_ms=0;
                COPY_WALLCLOCK (wallclock_ms);

                vb_show_progress (&last_pc_complete, vcf_file, z_file->vcf_data_so_far,
                                  z_file->disk_so_far - z_file->disk_at_beginning_of_this_vcf_file,
                                  wallclock_ms/1000, 
                                  done, test_mode);
            }
#endif
            if (done) {
                ASSERT (!z_file->vcf_data_size || 
                        z_file->vcf_data_size /* read from VCF file metadata */ == z_file->vcf_data_so_far, /* actually read */
                        "Error: VCF file uncompressed is of different length than expected - expected: %" PRIu64 " actual: %" PRIu64 "",
                        z_file->vcf_data_size, z_file->vcf_data_so_far);

                z_file->vcf_concat_data_size += z_file->vcf_data_so_far; // we completed one VCF file - add to the count
            }

            vb_release_vb (&th->vb); // cleanup vb and get it ready for another usage (without freeing memory)
            continue;
        }

    } while (!input_exhausted || num_running_compute_threads || next_vb);

#ifdef PROFILER
    int wallclock_ms=0;
    COPY_WALLCLOCK (wallclock_ms);

    fprintf (stderr, "PROFILER:\nwallclock: %u\n%s\n", wallclock_ms, profiler_print_report (&pseudo_vb->profile));
#endif

finish:
    z_file->disk_size = z_file->disk_so_far;
    vcf_file->vcf_data_size = z_file->vcf_data_size = vcf_file->vcf_data_so_far;
    buf_free (&compute_threads_buf);
    vb_release_vb (&pseudo_vb);
}
