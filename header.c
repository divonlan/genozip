// ------------------------------------------------------------------
//   header_vcf.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "zfile.h"
#include "txtfile.h"
#include "header.h"
#include "vblock.h"
#include "crypt.h"
#include "version.h"
#include "endianness.h"
#include "file.h"
#include "samples.h"
#include "dispatcher.h"
#include "zip.h"
#include "piz.h"
#include "zfile.h"
#include "txtfile.h"
#include "strings.h"

static bool is_first_txt = true; 

// these names go into the dictionary names on disk. to preserve backward compatibility, they should not be changed.
// (names are not longer than 8=DICT_ID_LEN as the code assumes it)
const char *field_names[NUM_DATATYPES][MAX_NUM_FIELDS_PER_DATA_TYPE] = FIELD_NAMES;

const unsigned datatype_last_field[NUM_DATATYPES]      = DATATYPE_LAST_FIELD;
const unsigned chrom_did_i_by_dt[NUM_DATATYPES]        = CHROM_DID_I_BY_DT; 

// -----------
// VCF stuff
// -----------
uint32_t global_vcf_num_samples           = 0; // number of samples in the file
uint32_t global_vcf_num_displayed_samples = 0; // PIZ only: number of samples to be displayed - might be less that global_vcf_num_samples if --samples is used
Buffer global_vcf_header_line = EMPTY_BUFFER;  // header line of first VCF file read - use to compare to subsequent files to make sure they have the same header during concat

static bool header_vcf_set_globals(const char *filename, Buffer *vcf_header)
{
    static const char *vcf_header_line_filename = NULL; // file from which the header line was taken

    // count tabs in last line which should be the field header line
    unsigned tab_count = 0;
    for (int i=vcf_header->len-1; i >= 0; i--) {
        
        if (vcf_header->data[i] == '\t')
            tab_count++;
        
        // if this is the beginning of field header line 
        else if (vcf_header->data[i] == '#' && (i==0 || vcf_header->data[i-1] == '\n' || vcf_header->data[i-1] == '\r')) {
        
            // if first vcf file - copy the header to the global
            if (!buf_is_allocated (&global_vcf_header_line)) {
                buf_copy (evb, &global_vcf_header_line, vcf_header, 1, i, vcf_header->len - i, "global_vcf_header_line", 0);
                vcf_header_line_filename = filename;
            }

            // ZIP only: subsequent files - if we're in concat mode just compare to make sure the header is the same
            else if (flag_concat && 
                     (vcf_header->len-i != global_vcf_header_line.len || memcmp (global_vcf_header_line.data, &vcf_header->data[i], global_vcf_header_line.len))) {

                fprintf (stderr, "%s: skipping %s: it has a different VCF header line than %s, see below:\n"
                                 "========= %s =========\n"
                                 "%.*s"
                                 "========= %s ==========\n"
                                 "%.*s"
                                 "=======================================\n", 
                         global_cmd, filename, vcf_header_line_filename,
                         vcf_header_line_filename, (uint32_t)global_vcf_header_line.len, global_vcf_header_line.data,
                         filename, (uint32_t)vcf_header->len-i, &vcf_header->data[i]);
                return false;
            }

            //count samples
            global_vcf_num_samples = (tab_count >= 9) ? tab_count-8 : 0; // note: a VCF file without samples would have tab_count==7 (8 fields) and is perfectly legal
            global_vcf_num_displayed_samples = global_vcf_num_samples;

            ASSERT0 (tab_count != 8, "Error: invalid VCF file - field header line contains a FORMAT field but no samples");

            ASSERT (tab_count >= 7, "Error: invalid VCF file - field header line contains only %d fields, expecting at least 8", tab_count+1);

            // if --samples is used, update vcf_header and global_vcf_num_displayed_samples
            if (flag_samples) samples_digest_vcf_header (vcf_header);

            return true; 
        }
    }

    ABORT ("Error: invalid VCF file - it does not contain a field header line; tab_count=%u", tab_count+1);
    return false; // avoid complication warnings
}

// genocat: remove FORMAT and sample names from the vcf header line, in case of --drop-genotypes
void header_vcf_trim_header_line (Buffer *vcf_header_buf)
{
    static const char *standard = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    for (int i=vcf_header_buf->len-2; i >= 0; i--) // -2 - skip last newline
        if (vcf_header_buf->data[i] == '\n') { 
            if (!memcmp (&vcf_header_buf->data[i+1], standard, MIN (strlen (standard), vcf_header_buf->len-(i+1)))) {
                vcf_header_buf->len = i + strlen(standard) + 2; // fixed length of standard VCF header up to INFO
                vcf_header_buf->data[vcf_header_buf->len-1] = '\n';   
            }                  
            return;
        }
    // no newline found - nothing to do
}

// genocat: remove all lines but the last from the vcf header, in case of --header-one
void header_vcf_keep_only_last_line (Buffer *vcf_header_buf)
{
    for (int i=vcf_header_buf->len-2; i >= 0; i--) // -2 - skip last newline
        if (vcf_header_buf->data[i] == '\n') {
            vcf_header_buf->len = vcf_header_buf->len - i - 1;
            memcpy (vcf_header_buf->data, &vcf_header_buf->data[i+1], vcf_header_buf->len);
            break;
        }
}

#define V1_VCF_HEADER // select the vcf_header functions of v1.c
#include "v1_vcf.c"

// ---------------------------
// Methods for all data types
// ---------------------------

// PIZ: called before reading each genozip file
void header_initialize(void)
{
    is_first_txt = true;

    // initialize VCF stuff
    global_vcf_num_samples           = 0;
    global_vcf_num_displayed_samples = 0;
    buf_free (&global_vcf_header_line);
}

// ZIP: reads VCF or SAM header and writes its compressed form to the GENOZIP file
bool header_txt_to_genozip (uint32_t *txt_line_i)
{    
    // data type muliplexors
    static const char first_char     [NUM_DATATYPES] = TXT_HEADER_LINE_FIRST_CHAR;
    static const bool header_allowed [NUM_DATATYPES] = TXT_HEADER_IS_ALLOWED;
    static const bool header_required[NUM_DATATYPES] = TXT_HEADER_IS_REQUIRED;
    
    z_file->disk_at_beginning_of_this_txt_file = z_file->disk_so_far;

    if (header_allowed[txt_file->data_type])
        txtfile_read_header (is_first_txt, header_required[txt_file->data_type], first_char[txt_file->data_type]); // reads into evb->txt_data and evb->lines.len
    
    *txt_line_i += (uint32_t)evb->lines.len;

    // for vcf, we need to check if the samples are the same before approving concatanation.
    // other data types can concatenate without restriction
    bool can_concatenate = (txt_file->data_type == DT_VCF) ? header_vcf_set_globals(txt_file->name, &evb->txt_data)
                                                           : true;
    if (!can_concatenate) { 
        // this is the second+ file in a concatenation list, but its samples are incompatible
        buf_free (&evb->txt_data);
        return false;
    }

    // we always write the txt_header section, even if we don't actually have a header, because the section
    // header contains the data about the file
    if (z_file) zfile_write_txt_header (&evb->txt_data, is_first_txt); // we write all headers in concat mode too, to support --split

    txt_file->section_bytes[SEC_TXT_HEADER] = evb->txt_data.len;
    z_file  ->section_bytes[SEC_TXT_HEADER] = evb->z_section_bytes[SEC_TXT_HEADER]; // comes from comp_compress
    z_file  ->num_sections [SEC_TXT_HEADER]++;
    z_file  ->num_txt_components_so_far++; // when compressing

    buf_free (&evb->txt_data);
    
    is_first_txt = false;

    return true; // everything's good
}

// PIZ: returns offset of header within data, EOF if end of file
bool header_genozip_to_txt (Md5Hash *digest) // NULL if we're just skipped this header (2nd+ header in concatenated file)
{
    z_file->disk_at_beginning_of_this_txt_file = z_file->disk_so_far;
    static Buffer header_section = EMPTY_BUFFER;

    int header_offset = zfile_read_section (evb, 0, NO_SB_I, &header_section, "header_section", 
                                            sizeof(SectionHeaderTxtHeader), SEC_TXT_HEADER, NULL);
    if (header_offset == EOF) {
        buf_free (&header_section);
        return false; // empty file (or in case of split mode - no more components) - not an error
    }

    // handle the GENOZIP header of the txt header section
    SectionHeaderTxtHeader *header = (SectionHeaderTxtHeader *)header_section.data;

    ASSERT (!digest || BGEN32 (header->h.compressed_offset) == crypt_padded_len (sizeof(SectionHeaderTxtHeader)), 
            "Error: invalid txt header's header size: header->h.compressed_offset=%u, expecting=%u", BGEN32 (header->h.compressed_offset), (unsigned)sizeof(SectionHeaderTxtHeader));

    // in split mode - we open the output txt file of the component
    if (flag_split) {
        ASSERT0 (!txt_file, "Error: not expecting txt_file to be open already in split mode");
        txt_file = file_open (header->txt_filename, WRITE, TXT_FILE, z_file->data_type);
        txt_file->txt_data_size_single = BGEN64 (header->txt_data_size);       
    }

    txt_file->max_lines_per_vb = BGEN32 (header->max_lines_per_vb);

    if (is_first_txt || flag_split) 
        z_file->num_lines = BGEN64 (header->num_lines);

    if (flag_split) *digest = header->md5_hash_single; // override md5 from genozip header
        
    // now get the text of the VCF header itself
    static Buffer header_buf = EMPTY_BUFFER;
    zfile_uncompress_section (evb, header, &header_buf, "header_buf", SEC_TXT_HEADER);

    bool is_vcf = (z_file->data_type == DT_VCF);

    bool can_concatenate = is_vcf ? header_vcf_set_globals(z_file->name, &header_buf) : true;
    if (!can_concatenate) {
        buf_free (&header_section);
        buf_free (&header_buf);
        return false;
    }

    if (is_vcf && flag_drop_genotypes) header_vcf_trim_header_line (&header_buf); // drop FORMAT and sample names

    if (is_vcf && flag_header_one) header_vcf_keep_only_last_line (&header_buf);  // drop lines except last (with field and samples name)

    // write vcf header if not in concat mode, or, in concat mode, we write the vcf header, only for the first genozip file
    if ((is_first_txt || flag_split) && !flag_no_header)
        txtfile_write_to_disk (&header_buf);
    
    buf_free (&header_section);
    buf_free (&header_buf);

    z_file->num_txt_components_so_far++;
    is_first_txt = false;

    return true;
}

const char *dt_name (DataType dt)
{
    static const char *names[NUM_DATATYPES] = DATATYPE_NAMES;
    return type_name (dt, &names[dt], sizeof(names)/sizeof(names[0]));
}
