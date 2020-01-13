// ------------------------------------------------------------------
//   zip.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"

// read entire variant block to memory. this is called from the dispatcher thread
static void zip_read_variant_block (File *vcf_file,
                                    unsigned *line_i,   // in/out next line to be read
                                    Buffer *first_data_line,    // first line might be supplied by caller 
                                    VariantBlock *vb)
{
    START_TIMER;

    unsigned first_line= *line_i;

    unsigned vb_line_i;
    vb->vb_data_size = 0; // size of variant block as it appears in the source file
    for (vb_line_i=0; vb_line_i < VARIANTS_PER_BLOCK; vb_line_i++) 
    {
        DataLine *dl = &vb->data_lines[vb_line_i];

        if (vb_line_i > 0 || !first_data_line) { // first line might be supplied by caller

            // allocate line Buffer and read line from file 
            bool success = vcffile_get_line (vb, first_line + vb_line_i, &dl->line, "dl->line"); 
            if (!success) break; // no more lines - we're done
        }
        else {
            buf_copy (vb, &dl->line, first_data_line, 1, 0, 0, "dl->line", vb->variant_block_i);
            buf_free (first_data_line); 
        }
        dl->line_i = first_line + vb_line_i;

        (*line_i)++;
        
        // count bytes in the source file, for sanity check after we reconstruct
        vb->vb_data_size += dl->line.len;
    }

    vb->num_lines = vb_line_i;

    COPY_TIMER (vb->profile.zip_read_variant_block);
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

#define SBL(line_i,sb_i) ((line_i) * vb->num_sample_blocks + (sb_i))

static unsigned zip_get_genotype_vb_start_len (VariantBlock *vb)
{
    buf_alloc (vb, &vb->genotype_section_lens_buf, sizeof(unsigned) * vb->num_sample_blocks, 1, "section_lens_buf", 0);
    unsigned section_0_len = 0; // all sections are the same length except the last that might be shorter

    // offsets into genotype data of individual lines
    buf_alloc (vb, &vb->gt_sb_line_starts_buf, 
               vb->num_lines * vb->num_sample_blocks * sizeof(uint32_t*), 
               0, "gt_sb_line_starts_buf", vb->first_line);
    uint32_t **gt_sb_line_starts = (uint32_t**)vb->gt_sb_line_starts_buf.data; 
    
    // each entry is the length of a single line in a sample block
    buf_alloc (vb, &vb->gt_sb_line_lengths_buf, 
               vb->num_lines * vb->num_sample_blocks * sizeof(unsigned), 
               0, "gt_sb_line_lengths_buf", vb->first_line);
    unsigned *gt_sb_line_lengths = (unsigned *)vb->gt_sb_line_lengths_buf.data; 
    
    // calculate offsets and lengths of genotype data of each sample block
    for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {

        uint32_t *gt_data  = (uint32_t*)vb->data_lines[line_i].genotype_data.data;
        unsigned num_subfields_in_cell = vb->data_lines[line_i].num_subfields;

        for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

            unsigned num_samples_in_sb = vb_num_samples_in_sb (vb, sb_i);

            gt_sb_line_starts[SBL(line_i, sb_i)] = 
                &gt_data[SAMPLES_PER_BLOCK * sb_i * num_subfields_in_cell];

            unsigned num_subfields_in_sample_line = num_subfields_in_cell * num_samples_in_sb; // number of uint32_t
            gt_sb_line_lengths[SBL(line_i, sb_i)] = num_subfields_in_sample_line;

            if (!sb_i) section_0_len += num_subfields_in_sample_line;
        }
    }

    return section_0_len; // in subfields (uint32_t)
}

// split genotype data to sample groups, within a sample group genotypes are separated by a tab
static void zip_generate_genotype_one_section (VariantBlock *vb, unsigned sb_i)
{
    START_TIMER;

    // build sample block genetype data
    uint8_t *dst_next = (uint8_t *)vb->genotype_one_section_data.data;
    
    // move the GT items from the line data to the permuted data - with each 
    // sample block of gt data containing the data in TRANSPOSED order - i.e. first
    // the gt data for all the variants for sample 1, then all of samples 2 etc.
    unsigned num_samples_in_sb = vb_num_samples_in_sb (vb, sb_i);
    for (unsigned sample_i=0; sample_i < num_samples_in_sb; sample_i++) {
        for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {

            DataLine *dl = &vb->data_lines[line_i];
            unsigned **sb_lines = (uint32_t**)vb->gt_sb_line_starts_buf.data;
            unsigned *this_line = sb_lines[SBL(line_i, sb_i)];
            
            // lookup word indices in the global dictionary for all the subfields
            const uint8_t *dst_start = dst_next;
            for (unsigned sf=0; sf < dl->num_subfields; sf++) { // iterate on the order as in the line
            
                uint32_t node_index = this_line[dl->num_subfields * sample_i + sf];

                if (node_index <= SEG_MAX_INDEX) { // normal index
                    MtfContext *ctx = &vb->mtf_ctx[dl->sf_i[sf]];
                    ASSERT (node_index < ctx->mtf.len, 
                            "Error: node index out of range: node_index=%u len=%u", node_index, ctx->mtf.len);
                    
                    MtfNode *node = &((MtfNode*)ctx->mtf.data)[node_index];
                    Base250 index = node->word_index;

                    if (index.num_numerals == 1) // shortcut for most common case
                        *(dst_next++) = index.numerals[0];
                    
                    else {
                        memcpy (dst_next, &index.numerals, index.num_numerals);
                        dst_next += index.num_numerals;
                    }
                }
                else if (node_index == SEG_MISSING_SF) 
                    *(dst_next++) = BASE250_MISSING_SF;

                else  // node_index == SEG_EMPTY_SF
                    *(dst_next++) = BASE250_EMPTY_SF;
            }
            
            vb->genotype_one_section_data.len += dst_next - dst_start;
            vb->add_bytes[SEC_GENOTYPE_DATA] += dst_next - dst_start - sizeof(uint32_t) * dl->num_subfields; 
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
        vb->phase_sections_data = (Buffer *)calloc (vb->num_sample_blocks, sizeof(Buffer)); // allocate once, never free

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
        vb->haplotype_sections_data = (Buffer *)calloc (vb->num_sample_blocks, sizeof(Buffer)); // allocate once, never free

    // TO DO : there's a bug in this allocation ^ if a user is compressing multiple VCF, they 
    // may not have the same num_sample_blocks. We should remember how many we allocated and add more
    // if needed.
    
    // create a permutation index for the whole variant block, and permuted haplotypes for each sample block        
    buf_alloc (vb, &vb->haplotype_permutation_index, vb->num_haplotypes_per_line * sizeof(unsigned), 
               0, "haplotype_permutation_index", vb->first_line);

    buf_alloc (vb, &vb->helper_index_buf, vb->num_haplotypes_per_line * sizeof(HaploTypeSortHelperIndex), 0,
               "helper_index_buf", vb->first_line);
    memset (vb->helper_index_buf.data, 0, vb->helper_index_buf.size);
    HaploTypeSortHelperIndex *helper_index = (HaploTypeSortHelperIndex *)vb->helper_index_buf.data;

    { 
        START_TIMER; 
        // build index array 
        for (unsigned ht_i=0; ht_i < vb->num_haplotypes_per_line; ht_i++) 
            helper_index[ht_i].index_in_original_line = ht_i;

        for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {
            for (unsigned ht_i=0; ht_i < vb->num_haplotypes_per_line; ht_i++) {

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
                   0, "haplotype_sections_data", vb->first_line);

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
    // this will be included in the genozip file
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
static void zip_compress_variant_block (VariantBlock *vb)
{ 
    START_TIMER;

    // allocate memory for the final compressed data of this vb. allocate 20% of the 
    // vb size on the original file - this is normally enough. if not, we will realloc downstream
    buf_alloc (vb, &vb->z_data, vb->vb_data_size / 5, 1.2, "z_data", 0);

    // initalize variant block data (everything else is initialzed to 0 via calloc)
    vb->phase_type = PHASE_UNKNOWN;  // phase type of this block
    vb->num_samples_per_block = SAMPLES_PER_BLOCK;
    vb->num_sample_blocks = ceil((float)global_num_samples / (float)vb->num_samples_per_block);

    // if testing, make a copy of the original lines read from the file, to compare to the result of vunblocking later
    Buffer *lines_orig = NULL;
    unsigned max_genotype_section_len=0; // length in subfields

    // clone global dictionaries while granted exclusive access to the global dictionaries
    mtf_clone_ctx (vb);

    // split each lines in this variant block to its components
    seg_all_data_lines (vb, lines_orig);

    // for the first vb only - sort dictionaries so that the most frequent entries get single digit
    // base-250 indices. This can be done only before any dictionary is written to disk, but likely
    // beneficial to all vbs as they are likely to more-or-less have the same frequent entries
    if (vb->variant_block_i == 1)
        mtf_sort_dictionaries_vb_1(vb);

    // if block has haplotypes - handle them now
    if (vb->has_haplotype_data)
        zip_generate_haplotype_sections (vb); 

    // if block has genetype data 
    if (vb->has_genotype_data) 
        // calculate starts, lengths and allocate memory
        max_genotype_section_len = zip_get_genotype_vb_start_len (vb); // length in subfields
    
    // if block has phase data - handle it
    if (vb->phase_type == PHASE_MIXED_PHASED) 
        zip_generate_phase_sections (vb);

    // generate & write variant data for the variant block
    zip_generate_variant_data_section (vb);

    unsigned variant_data_header_pos = vb->z_data.len;
    zfile_compress_variant_data (vb);

    // merge new words added in this vb into the z_file.mtf_ctx, ahead of zip_generate_genotype_one_section()
    // writing indices based on the merged dictionaries. dictionaries are compressed. 
    // all this is done while holding exclusive access to the z_file dictionaries
    unsigned num_dictionary_sections = mtf_merge_in_vb_ctx(vb);

    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {
        
        if (vb->has_genotype_data) {

            // in the worst case scenario, each subfield is represnted by 5 bytes in Base250
            // TO DO: this is a huge waste of memory, as normally each subfield consumes ~ 1.5 B, better allocate less and realloc
            buf_alloc (vb, &vb->genotype_one_section_data, max_genotype_section_len * 5, 1, "genotype_one_section_data", sb_i);

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

    // note: this updates the z_data in memory (not on disk)
    zfile_update_compressed_variant_data_header (vb, variant_data_header_pos, num_dictionary_sections);

    // updates stats for --showcontent - reduce by the bytes added by us over the bytes in the original file (might be a negative value if we saved bytes rather than added)
    for (unsigned sec_i=0; sec_i < NUM_SEC_TYPES; sec_i++) {
        //printf ("sec: %u, bytes before compression: %u, bytes added by alg: %d\n", sec_i, vb->vcf_section_bytes[sec_i], vb->add_bytes[sec_i]);
        vb->vcf_section_bytes[sec_i] -= vb->add_bytes[sec_i];
    }
    
    COPY_TIMER (vb->profile.compute);
}

// this is the main dispatcher function. It first processes the VCF header, then proceeds to read 
// a variant block from the input file and send it off to a thread for computation. When the thread
// completes, this function proceeds to write the output to the output file. It can dispatch
// several threads in parallel.
void zip_dispatcher (const char *vcf_basename, File *vcf_file, 
                     File *z_file, bool test_mode, unsigned max_threads, bool is_last_file)
{
    extern int flag_show_alleles; // defined in main()
    Dispatcher dispatcher = dispatcher_init (max_threads, POOL_ID_ZIP, vcf_file, z_file, test_mode, !flag_show_alleles, vcf_basename);

    unsigned line_i = 0; // last line read (first line in file = 1, consistent with script line numbers)

    // first compress the VCF header
    Buffer *first_data_line = NULL; // contains a value only for first variant block, otherwise empty. 
    
    // read the vcf header, assign the global variables, and write the compressed header to the GENOZIP file
    bool success = vcf_header_vcf_to_genozip (dispatcher_get_pseudo_vb (dispatcher), &line_i, &first_data_line);
    if (!success) goto finish;

    if (!first_data_line) goto finish; // VCF file has only a header or is an empty file - no data - we're done

    mtf_initialize_mutex (z_file);

    // this is the dispatcher loop. In each iteration, it can do one of 3 things, in this order of priority:
    // 1. In there is a new variant block avaialble, and a compute thread available to take it - dispatch it
    // 2. If there is no new variant block available, but input is not exhausted yet - read one
    // 3. Wait for the first thread (by sequential order) to complete the compute and output the results
    VariantBlock *next_vb;
    do {
        next_vb = dispatcher_get_next_vb (dispatcher);
        
        // PRIORITY 1: is there a block available and a compute thread available? in that case dispatch it
        if (next_vb && next_vb->ready_to_dispatch && dispatcher_has_free_thread (dispatcher)) 
            dispatcher_compute (dispatcher, zip_compress_variant_block);
        
        // PRIORITY 2: If there is no new variant block available, but input is not exhausted yet - read one
        else if (!next_vb && !dispatcher_is_input_exhausted (dispatcher)) {

            next_vb = dispatcher_generate_next_vb (dispatcher);

            next_vb->first_line = line_i;

            zip_read_variant_block (vcf_file, &line_i, first_data_line, next_vb);
            first_data_line = NULL;

            if (next_vb->num_lines)  // we found some data 
                next_vb->ready_to_dispatch = true;
            else
                dispatcher_input_exhausted (dispatcher);
        }

        // PRIORITY 3: Wait for the first thread (by sequential order) to complete the compute and output the results
        // dispatch computing the compressed to a separate compute thread
        else {
            VariantBlock *processed_vb = dispatcher_get_next_processed_vb (dispatcher);
            if (!processed_vb) continue;
            
            START_TIMER;
            fwrite (processed_vb->z_data.data, 1, processed_vb->z_data.len, (FILE *)z_file->file);
            COPY_TIMER (processed_vb->profile.write);

            z_file->disk_so_far     += (long long)processed_vb->z_data.len;
            z_file->vcf_data_so_far += (long long)processed_vb->vb_data_size;
            z_file->num_lines       += (long long)processed_vb->num_lines;

            // update section stats
            for (unsigned sec_i=1; sec_i < NUM_SEC_TYPES; sec_i++) {
                z_file->section_bytes[sec_i]   += processed_vb->z_section_bytes[sec_i];
                vcf_file->section_bytes[sec_i] += processed_vb->vcf_section_bytes[sec_i];
            }

            if (dispatcher_is_final_processed_vb (dispatcher)) {

                ASSERT (!z_file->vcf_data_size || 
                        z_file->vcf_data_size /* read from VCF file metadata */ == z_file->vcf_data_so_far, /* actually read */
                        "Error: VCF file length inconsistency - read from VCF file metadata: %" PRIu64 " actually read: %" PRIu64 "",
                        z_file->vcf_data_size, z_file->vcf_data_so_far);

                vcf_file->vcf_data_size = z_file->vcf_data_size = vcf_file->vcf_data_so_far;

                z_file->vcf_concat_data_size += z_file->vcf_data_so_far; // we completed one VCF file - add to the count
            }

            dispatcher_finalize_one_vb (dispatcher, vcf_file, z_file->vcf_data_so_far,
                                        z_file->disk_so_far - z_file->disk_at_beginning_of_this_vcf_file);

        }
    } while (!dispatcher_is_done (dispatcher));

    // go back and update some fields in the vcf header's section header - only if we can go back
    // (i.e. not output redirected). we might need to re-encrypt.
    extern int flag_concat_mode; // set in main()
    if ((is_last_file || !flag_concat_mode) && z_file && z_file->type == GENOZIP) 
        zfile_update_vcf_header_section_header (dispatcher_get_pseudo_vb (dispatcher));

finish:
    z_file->disk_size = z_file->disk_so_far;
    vcf_file->vcf_data_size = z_file->vcf_data_size = vcf_file->vcf_data_so_far; // just in case its not set already
    
    dispatcher_finish (dispatcher);
}
