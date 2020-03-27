// ------------------------------------------------------------------
//   samples.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "buffer.h"
#include "samples.h"

// referring to sample strings from the --samples command line option
static Buffer cmd_samples_buf = EMPTY_BUFFER; // an array of (char *)
static bool cmd_is_negative_samples = false;

// referring to samples in the vcf file
char *vcf_samples_is_included;                // a bytemap indicating for each sample if it is included
static char **vcf_sample_names;               // an array of char * to null-terminated names of samples 
static char *vcf_sample_names_data;           // vcf_sample_names point into here

void samples_add (const char *samples_str)
{
    ASSERT0 (samples_str, "Error: samples_str is NULL");

    bool is_negated = samples_str[0] == '^';

    bool is_conflicting_negation = (cmd_samples_buf.len && (cmd_is_negative_samples != is_negated));
    ASSERT0 (!is_conflicting_negation, "Error: inconsistent negation - all samples listed must either be negated or not");

    cmd_is_negative_samples = is_negated;

    // make a copy of the string and leave the original one for error message. 
    // we don't free this memory as chrom fields in regions will be pointing to it
    char *next_region_token = malloc (strlen (samples_str)+1); 
    strcpy (next_region_token, samples_str + is_negated); // drop the ^ if there is one

    while (1) {
        char *one_sample = strtok_r (next_region_token, ",", &next_region_token);
        if (!one_sample) break;

        bool is_duplicate = false;
        for (unsigned s=0; s < cmd_samples_buf.len; s++)
            if (!strcmp (one_sample, *ENT (char *, &cmd_samples_buf, s))) {
                is_duplicate = true;
                break;
            }
        if (is_duplicate) continue; // skip duplicates "genocat -s sample1,sample2,sample1"

        buf_alloc (evb, &cmd_samples_buf, MAX (cmd_samples_buf.len + 1, 100) * sizeof (char*), 2, "cmd_samples_buf", 0);

        *NEXTENT (char *, &cmd_samples_buf) = one_sample;
    }
}

// accept a sample from the vcf file's samples as consistent with the --samples requested
static void samples_accept (Buffer *vcf_header_buf, const char *sample_str)
{
    buf_add_string (evb, vcf_header_buf, sample_str);
    buf_add_string (evb, vcf_header_buf, "\t");
    global_number_displayed_samples++;
}

// processes the vcf header sample line according to the --samples option, removing samples that are not required
// and building a bytemap PIZ filter the samples
void samples_digest_vcf_header (Buffer *vcf_header_buf)
{
    static const char *standard = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";

    int i=vcf_header_buf->len-2 ; for (; i >= 0; i--) // -2 - skip last newline
        if (vcf_header_buf->data[i] == '\n') { 
            bool header_matches_standard = !memcmp (&vcf_header_buf->data[i+1], standard, MIN (strlen (standard), vcf_header_buf->len-(i+1)));
            if (!header_matches_standard) flag_samples = false;
            RETURNW0 (header_matches_standard,, "Warning: found non-standard VCF sample header line. Ingoring --samples");

            break;
        }

    vcf_samples_is_included = malloc (global_num_samples);
    memset (vcf_samples_is_included, cmd_is_negative_samples, global_num_samples); // 0 if not included unless list says so (positive) and vice versa

    unsigned vcf_names_start_index = i + 1 + strlen(standard);
    unsigned vcf_names_data_len = vcf_header_buf->len - vcf_names_start_index;
    vcf_sample_names_data = malloc (vcf_names_data_len);
    memcpy (vcf_sample_names_data, &vcf_header_buf->data[vcf_names_start_index], vcf_names_data_len);
    vcf_sample_names_data[vcf_names_data_len-1] = '\t'; // change last separator from \n to \t

    vcf_sample_names = malloc (global_num_samples * sizeof (char *));

    vcf_header_buf->len = vcf_names_start_index;
    char *next_token = vcf_sample_names_data;
    
    // go through the vcf file's samples and add those that are consistent with the --samples requested
    global_number_displayed_samples = 0;
    for (unsigned i=0; i < global_num_samples; i++) {
        vcf_sample_names[i] = strtok_r (next_token, "\t", &next_token);

        for (unsigned s=0; s < cmd_samples_buf.len; s++)
            if (!strcmp (vcf_sample_names[i], *ENT(char *, &cmd_samples_buf, s))) { // found
                vcf_samples_is_included[i] = !cmd_is_negative_samples;
                if (!cmd_is_negative_samples) samples_accept (vcf_header_buf, vcf_sample_names[i]);

                // remove this sample from the --samples list as we've found it already
                memcpy (ENT(char *, &cmd_samples_buf, s), ENT(char *, &cmd_samples_buf, s+1), (cmd_samples_buf.len-s-1) * sizeof (char *));
                cmd_samples_buf.len--;
            }
            else if (cmd_is_negative_samples) samples_accept (vcf_header_buf, vcf_sample_names[i]);
    }

    vcf_header_buf->data[vcf_header_buf->len-1] = '\n'; // change last separator from \t

    // warn about any --samples items that were not found in the vcf file (all the ones that still remain in the buffer)
    for (unsigned s=0; s < cmd_samples_buf.len; s++) 
        ASSERTW (false, "Warning: requested sample '%s' is not found in the VCF file, ignoring it", *ENT(char *, &cmd_samples_buf, s));

    //for (i=0; i<global_num_samples; i++) fprintf (stderr, "%u ", vcf_samples_is_included[i]); 
}

// PIZ only: calculates whether a sample block is included, based on --samples. this is called once per sample block
// in zfile_read_one_vb to set vb->is_sb_included, and thereafter vb->is_sb_included is used
bool samples_is_sb_included (uint32_t num_samples_per_block, uint32_t sb_i)
{
    if (!flag_samples) return true; // all sample blocks are included if --samples if not specified

    for (uint32_t sample_i = sb_i * num_samples_per_block; 
         sample_i < MIN ((sb_i+1) * num_samples_per_block, global_num_samples);
         sample_i++)

         if (samples_am_i_included(sample_i)) return true; // at least one sample is included - means the sample block is included

    return false; // no sample in this sample block is included - means the whole sample block is excluded
}
