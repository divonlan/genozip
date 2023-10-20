// ------------------------------------------------------------------
//   filename.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.


#ifdef __linux__
#include <linux/limits.h>
#elif defined __APPLE__
#include <sys/syslimits.h>
#endif
#include "genozip.h"
#include "filename.h"
#include "file.h"
#include "codec.h"
#include "url.h"
#include "tar.h"

rom filename_z_by_flag (void)
{
    char *z_filename = (char*)MALLOC(strlen (flag.out_filename) + 1);
    
    strcpy (z_filename, flag.out_filename);
    return z_filename;
}

rom filename_z_normal (rom txt_filename, DataType dt, FileType txt_ft)
{
    if (!txt_filename && (flag.biopsy || flag.biopsy_line.line_i != NO_LINE))
        txt_filename = "dummy"; // we don't have a txt_filename, but that's ok, because we don't need it

    ASSINP0 (txt_filename || !tar_zip_is_tar(), "Piping from stdin is not supported when using --tar");
    ASSINP0 (txt_filename, "Use --output to specify the output filename (with a .genozip extension)");

    unsigned dn_len = flag.out_dirname ? strlen (flag.out_dirname) : 0;
    bool is_url = url_is_url (txt_filename);
    rom basename = (dn_len || is_url) ? filename_base (txt_filename, false, "", 0, 0) : NULL;

    rom local_txt_filename = basename ? basename : txt_filename; 

    unsigned fn_len = strlen (local_txt_filename);

    char *z_filename = (char *)MALLOC (fn_len + dn_len + 30); // add enough the genozip extension e.g. 23andme.genozip

    // if the file has an extension matching its type, replace it with the genozip extension, if not, just add the genozip extension
    rom genozip_ext = file_exts[file_get_z_ft_by_txt_in_ft (dt, txt_ft)];

    if (chain_is_loaded && (TXT_DT(VCF) || TXT_DT(BCF)))
        genozip_ext = DVCF_GENOZIP_;

    if (filename_has_ext (local_txt_filename, file_exts[txt_ft]))
        sprintf (z_filename, "%s%s%.*s%s", 
                 (dn_len ? flag.out_dirname : ""), (dn_len ? "/" : ""), 
                 (int)(fn_len - strlen (file_exts[txt_ft])), local_txt_filename, genozip_ext); 
    else 
        sprintf (z_filename, "%s%s%s%s", 
                 (dn_len ? flag.out_dirname : ""), (dn_len ? "/" : ""), 
                 local_txt_filename, genozip_ext); 

    FREE (basename);

    return z_filename;
}

// rules for pair filename: two filenames need to differ by exactly one character, which is '1' and '2',
// if test_only we validate the rule
// if !test_only, only fn1 needs to be provided, we assume the rule is valid, and we create the output file with "1+2" (caller should free)
char *filename_z_pair (rom fn1, rom fn2, bool test_only)
{
    FileType ft1, ft2;
    rom rn1=0, rn2=0;

    unsigned dn_len = flag.out_dirname ? strlen (flag.out_dirname) : 0;

    // case: new directory or url - take only the basename
    rom bn1 = (dn_len || url_is_url (fn1)) ? filename_base (fn1, 0,0,0,0) : NULL;
    rom bn2 = (dn_len || url_is_url (fn1)) ? filename_base (fn2, 0,0,0,0) : NULL;

    file_get_raw_name_and_type (bn1 ? bn1 : fn1, &rn1, &ft1);
    file_get_raw_name_and_type (bn2 ? bn2 : fn2, &rn2, &ft2);
    
    #define MY_RET(x) ({ FREE (rn1); FREE (rn2); FREE (bn1); FREE (bn2); return (x); })

    unsigned len = strlen (rn1);
    if (len != strlen (rn2)) MY_RET (NULL);

    int df = -1;
    for (unsigned i=0; i < len; i++)
        if (rn1[i] != rn2[i]) {
            if (df >= 0) MY_RET (NULL); // 2nd differing character
            df = i;
        }
    

    if (!((rn1[df] == '1' && rn2[df] == '2') || (rn1[df] == '2' && rn2[df] == '1'))) 
        MY_RET (NULL); // one of them must be '1' and the other '2'

    if (test_only) 
        MY_RET ((char *)1);

    char *pair_fn = MALLOC (len + dn_len + 20);
    sprintf (pair_fn, "%s%s%.*s1+2%s" FASTQ_GENOZIP_, 
             (dn_len ? flag.out_dirname : ""), (dn_len ? "/" : ""), 
             df, rn1, &rn1[df+1]);
    
    MY_RET (pair_fn);
}

char *filename_z_deep (rom sam_name)
{
    unsigned dn_len = flag.out_dirname ? strlen (flag.out_dirname) : 0;

    // case: new directory - take only the basename
    rom bn = (dn_len || url_is_url (sam_name)) ? filename_base (sam_name, 0,0,0,0) : NULL;

    rom raw_name = 0;
    file_get_raw_name_and_type (bn ? bn : sam_name, &raw_name, NULL);
    
    char *deep_fn = MALLOC (strlen (raw_name) + dn_len + 20);
    sprintf (deep_fn, "%s%s%s" DEEP_GENOZIP_, 
             (dn_len ? flag.out_dirname : ""), (dn_len ? "/" : ""), raw_name);
    
    FREE (raw_name);
    FREE (bn);

    return deep_fn;
}

// PIZ: guess original filename from uncompressed txt filename and compression algoritm (allocated memory)
rom filename_guess_original (ConstFileP file)
{
    if (!file) return "(is-NULL)";

    if (file->codec == CODEC_NONE) return file->name;

    unsigned len = strlen (file->name) + 10;
    char *org_name = MALLOC (len);
    strcpy (org_name, file->name);

    rom ext = codec_args[file->codec].ext;

    // remove existing extension if needed (eg when replacing .sam with .bam)
    if (ext[0] == '-') {
        char *last_dot = strrchr (org_name, '.');
        if (last_dot) *last_dot = 0;
    }

    int org_len = strlen (org_name); 
    if (org_len >= 4 && !strcmp (&org_name[org_len-4], ".bam")) {
        // don't add extension to .bam
    }
    else if (!filename_has_ext (org_name, &ext[1])) // add new extension
        strcpy (&org_name[org_len], &ext[1]);

    return org_name;
}

bool filename_has_ext (rom filename, rom extension)
{
    if (!filename) return false;

    unsigned ext_len = strlen (extension);
    unsigned fn_len  = strlen (filename);
    
    return fn_len >= ext_len && !strncmp (&filename[fn_len-ext_len], extension, ext_len);
}

// get basename of a filename - we write our own basename for Visual C and Windows compatibility
rom filename_base (rom filename, bool remove_exe, rom default_basename,
                   char *basename /* optional pre-allocated memory */, unsigned basename_size /* basename bytes */)
{
    if (!filename) filename = default_basename;

    unsigned len = strlen (filename);
    if (remove_exe && filename_has_ext (filename, ".exe")) len -= 4; // for Windows

    // get start of basename
    rom start = filename;
    for (int i=len-1; i >= 0; i--)
        if (filename[i]=='/' || filename[i]=='\\') {
            start = &filename[i+1];
            break;
        }

    len -= (start-filename);

    // remove url parameter, eg https://myserver.com/dir/file.txt?param=1 --> file.txt
    if (url_is_url (filename)) { 
        char *param = strchr (start, '?');
        if (param) len = param - start;
    }

    if (!basename) 
        basename = (char *)MALLOC (len + 1); // +1 for \0
    
    // case: basename too long for pre-allocated memory - use default
    else if (len+1 > basename_size) {
        start = default_basename;
        len = strlen (default_basename);
    }

    sprintf (basename, "%.*s", (int)len, start);

    return basename;
}

// eg "file.fastq.gz" -> "file.fastq"
void filename_remove_codec_ext (char *filename, FileType ft)
{
    rom codec_ext;
    if (!((codec_ext = strchr (&file_exts[ft][1], '.')))) return; // this file type's extension doesn't have a codec ext (it is, eg ".fastq" not ".fastq.gz")
    
    unsigned codec_ext_len = strlen (codec_ext);
    unsigned fn_len = strlen (filename);

    // make sure filename actually has the codec (eg its not eg "(stdin)")
    if (fn_len > codec_ext_len && !strcmp (&filename[fn_len-codec_ext_len], codec_ext))
        filename[fn_len-codec_ext_len] = 0; // shorten string
}

// change eg "c:\dir\file" to "/c/dir/file" and add full path (allocates memory - caller should free)
char *filename_make_unix (char *filename)
{
    char path[PATH_MAX];
    unsigned len = strlen (filename);
    
    if (
#ifdef _WIN32
        (len >= 2 && filename[1] == ':') || // full name - starting with eg C:
#endif
        (len >= 1 && filename[0] == '/') || // full name - starting with /
        !getcwd (path, sizeof (path))) // path too long
        path[0] = 0;  // don't store path

    char *full_fn = MALLOC (strlen (z_name) + strlen (path) + 2);
    sprintf (full_fn, "%s%s%s", path, *path ? "/" : "", z_name);

#ifdef _WIN32 // convert to Unix-style filename
    len = strlen (full_fn);

    if (len >= 2 && full_fn[1] == ':') {
        full_fn[1] = full_fn[0];
        full_fn[0] = '/';
    }

    for (unsigned i=0; i < len; i++)
        if (full_fn[i] == '\\') full_fn[i] = '/';
#endif

    return full_fn;
}

PairType filename_is_fastq_pair (STRp(fn1), STRp(fn2))
{
    if (fn1_len != fn2_len) return false; 

    int mismatches = 0, mm_i=0;
    for (int i=0; i < fn1_len && mismatches < 2; i++)
        if (fn1[i] != fn2[i]) {
            mismatches++;
            mm_i = i;
        }
    
    // its predicted to a pair if filenames are the same, except for '1'⇄'2' switch
    if  (mismatches == 1 && ((fn1[mm_i] == '1' && fn2[mm_i] == '2'))) return PAIR_R1; // fn1 is PAIR_1
    if  (mismatches == 1 && ((fn1[mm_i] == '2' && fn2[mm_i] == '1'))) return PAIR_R2; 
    return NOT_PAIRED;
}