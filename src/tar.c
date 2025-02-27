// ------------------------------------------------------------------
//   tar.c
//   Copyright (C) 2021-2025 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <pwd.h>
#include <grp.h>
#endif

#include "genozip.h"
#include "flags.h"
#include "file.h"
#include "endianness.h"
#include "buffer.h"
#include "arch.h"
#include "tar.h"
#include "license.h"

#define MAX_TAR_UID_GID 07777777 // must fit in 8 characters, inc. \0, printed in octal
#define MAX_TAR_MODE 07777777
#define MAX_TAR_MTIME 077777777777ULL

#define NOBODY  65534            // used for UID / GID in case real values are beyond MAX_TAR_UID_GID. https://wiki.ubuntu.com/nobody

// each file in the tar file is comprised a 512-byte "ustar" header followed by 512-byte blocks. tar file is terminated
// by two zero 512B blocks. see: http://manpages.ubuntu.com/manpages/bionic/man5/tar.5.html
typedef struct __attribute__ ((packed)){
    char name[100];     // nul-terminated
    char mode[8];       // nul-terminated octal in ASCII
    char uid[8];        // nul-terminated octal in ASCII
    char gid[8];        // Unix: nul-terminated octal in ASCII;  Windows: 0
    union {
        char std[12];   // If < 8GB, nul-terminated octal in ASCII (POSIX standard tar)
        struct {        // If >= 8GB, non-standard, supported by GNU tar and others
            char marker;// must be 0x80
            char unused[3];
            char size_bgen[8];
        } gnu;
    } size;
    char mtime[12];     // nul-terminated octal in ASCII
    char checksum[8];   // sum of header bytes (as unsigned char), in six-digit zero-padded octal, followed by nul and space. checksum is calculated while the "checksum" field is all-spaces.
    char typeflag[1];   // always '0' - regular file (genozip doesn't support symlinks, directories etc)
    char linkname[100]; // hard link destination - nul-termianted
    char magic[6];      // POSIX "ustar\0" (different from GNU which is "ustar ")
    char version[2];    // POSIX "00" (different from GNU which is " \0")
    char uname[32];
    char gname[32];
    char devmajor[8];
    char devminor[8];
    char prefix[155];   // nul-termianted. if given, full filename is prefix + '/' + name, otherwise it is just name
    char pad[12];
} HeaderPosixUstar;

static FILE *tar_file = NULL;
static rom tar_name = NULL;
static int64_t t_offset = 0; // after tar_open_file: file start offset within tar (after the file tar header)
static HeaderPosixUstar hdr = {};

#define GENOZIP_TAR_DIR_NAME "genozip-linux-x86_64"

static void tar_add_hard_link (rom fn_on_disk, rom fn_in_tar_src, rom fn_in_tar_dst); // forward

void tar_set_tar_name (rom tar_filename)
{
    tar_name = tar_filename;
}

rom tar_get_tar_name (void)
{
    return tar_name;
}

void tar_initialize (BufferP input_files_buf)
{
    ASSERTNOTNULL (tar_name);

    ASSINP (flag.force || !file_exists (tar_name), "file %s already exists, use --force (or -f) to overwrite", tar_name);

    tar_file = fopen (tar_name, flag.pair ? "wb+" : "wb"); // if --pair, when compressing pair2, we go back and read pair1
    ASSINP (tar_file, "cannot create tar file %s: %s", tar_name, strerror (errno)); 

    if (flag.is_linux && license_allow_distribution()) {
        // verify that top-level no file or directory is named "genozip-linux-x86_64" 
        for_buf (rom, fn_p, *input_files_buf) {
            int fn_len = strlen (*fn_p);
            if (str_issame_(*fn_p, fn_len, _S(GENOZIP_TAR_DIR_NAME)) || // filename is "genozip-linux-x86_64"
                (fn_len > STRLEN(GENOZIP_TAR_DIR_NAME) && !memcmp (*fn_p, GENOZIP_TAR_DIR_NAME "/", STRLEN(GENOZIP_TAR_DIR_NAME "/")))) // filename starts with "genozip-linux-x86_64/"
                return; // name already exist - we can't add our directory
        }

        rom exec = arch_get_executable().s;
        tar_copy_file (exec, GENOZIP_TAR_DIR_NAME "/genozip");
        tar_add_hard_link (exec, GENOZIP_TAR_DIR_NAME "/genozip", GENOZIP_TAR_DIR_NAME "/genocat");
        tar_add_hard_link (exec, GENOZIP_TAR_DIR_NAME "/genozip", GENOZIP_TAR_DIR_NAME "/genounzip");
        tar_add_hard_link (exec, GENOZIP_TAR_DIR_NAME "/genozip", GENOZIP_TAR_DIR_NAME "/genols");
    }
}

static void tar_copy_metadata_from_file (rom fn)
{
    struct stat64 st;    
    int ret = stat64 (fn, &st);
    ASSERT (!ret, "stat64 failed on %s: %s", fn, strerror(errno));

    // change UID and/or GID to NOBODY if they go beyond the maximum. This can happen, for example, if
    // SSSD (https://sssd.io/) uses UIDs from Microsoft Active Directory.
    if ((uint32_t)st.st_uid > MAX_TAR_UID_GID) {
        WARN_ONCE ("UID of %s (and perhaps others) is %u - beyond the maximum allowed by the tar file format; recording uid as %u (nobody). --quiet to suppress this message.",
                   fn, st.st_uid, NOBODY);
        st.st_uid = NOBODY;
    }

    if ((uint32_t)st.st_gid > MAX_TAR_UID_GID) {
        WARN_ONCE ("GID of %s (and perhaps others) is %u - beyond the maximum allowed by the tar file format; recording gid as %u (nogroup). --quiet to suppress this message.",
                   fn, st.st_gid, NOBODY);
        st.st_gid = NOBODY;
    }

    if ((uint32_t)st.st_mode > MAX_TAR_MODE) {
        WARN_ONCE ("mode of %s (and perhaps others) is 0%o - beyond the maximum allowed by the tar file format; recording mode as 0%o. --quiet to suppress this message.",
                   fn, st.st_mode, st.st_mode & 0777);
        st.st_mode &= 0777;
    }

    if ((uint64_t)st.st_mtime > MAX_TAR_MTIME) {
        WARN_ONCE ("mtime of %s (and perhaps others) is %"PRIu64" - beyond the maximum allowed by the tar file format; recording mtime as 0. --quiet to suppress this message.",
                   fn, (uint64_t)st.st_mtime/*note: 32 bit on Darwin*/);
        st.st_mtime = 0;
    }

    // convert to nul-terminated octal in ASCII
    snprintf (hdr.uid,   sizeof (hdr.uid),   "%.*o", (int)sizeof(hdr.uid)-1,  st.st_uid);
    snprintf (hdr.gid,   sizeof (hdr.gid),   "%.*o", (int)sizeof(hdr.gid)-1,  st.st_gid);
    snprintf (hdr.mode,  sizeof (hdr.mode),  "%.*o", (int)sizeof(hdr.mode)-1, st.st_mode);
    snprintf (hdr.mtime, sizeof (hdr.mtime), "%.*"PRIo64, (int)sizeof(hdr.mtime)-1, (uint64_t)st.st_mtime); // mtime is 64b on Windows and Linux, 32b on MacOS

    // case: windows - use current username, no gname
    #ifdef _WIN32        
        DWORD len = sizeof (hdr.uname);
        GetUserName (hdr.uname, &len); // fails if user name is longer than uname - that's ok, as uname is optional
        GetUserName (hdr.gname, &len); // in Windows, we set gname=uname

    // case: Unix - get uname and gname from uid and gid
    #else        
        struct passwd *pwd = getpwuid (st.st_uid); // NULL if failed - that's ok, as uname is optional
        if (pwd && strlen (pwd->pw_name) < sizeof (hdr.uname)) // room for \0
            strcpy (hdr.uname, pwd->pw_name);

        struct group *grp = getgrgid (st.st_gid);  // NULL if failed - that's ok, as gname is optional
        if (grp && strlen (grp->gr_name) < sizeof (hdr.gname)) // room for \0
            strcpy (hdr.gname, grp->gr_name);
    #endif
}

static void tar_fwrite (const void *data, uint32_t size, rom object)
{
    if (!size) return; // nothing to do

    ASSERTNOTNULL (tar_file);
    ASSERTNOTNULL (data);
    
    uint32_t bytes = fwrite (data, 1, size, tar_file); 

    ASSERT (bytes == size, "Error writing %s to %s on filesystem=%s - requested %u bytes but wrote only %u: (%u)%s", 
            object, tar_name, arch_get_filesystem_type (txt_file).s, size, bytes, errno, strerror (errno));
}

// filenames that have a last component longer than 99 characters don't fit in POSIX tar. Instead, we use a GNU-specific extension
// of storing the filename as a pseudo-file in the tarball. See: https://itecnote.com/tecnote/r-what-exactly-is-the-gnu-tar-longlink-trick/
static void tar_write_gnu_long_filename (STRp(fn_in_tar)/* length includes \0 */)
{
    HeaderPosixUstar ll_hdr = {
        .name     = "././@LongLink",
        .mode     = "0000644",
        .uid      = "0000000", 
        .gid      = "0000000", 
        .mtime    = "00000000000",
        .typeflag = { 'L' }, // long file name
        .magic    = "ustar",
        .version  = {'0', '0'},
        .uname    = "root",
        .gname    = "root",
        .checksum = {' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ' }, // required checksum value, while calculating the checksum (8 spaces)
    };
    snprintf (ll_hdr.size.std, sizeof (ll_hdr.size.std), "%o", fn_in_tar_len); // fn_in_tar_len includes \0

    // header checksum
    unsigned checksum = 0;
    for (unsigned i=0; i < 512; i++) checksum += ((bytes)&ll_hdr)[i];
    snprintf (ll_hdr.checksum, sizeof (ll_hdr.checksum), "%06o", checksum);

    tar_fwrite (&ll_hdr, 512, "LongLink header");
    tar_fwrite (STRa(fn_in_tar), "long filename");

    // pad to full block
    if (fn_in_tar_len % 512)
        tar_fwrite ((char[512]){}, ROUNDUP512(fn_in_tar_len) - fn_in_tar_len, "long filename padding");

    t_offset += 512 + ROUNDUP512(fn_in_tar_len);
}

// open z_file within tar, for writing
FILE *tar_open_file (rom fn_on_disk, rom fn_in_tar)
{
    ASSERTNOTNULL (tar_file);

    // initialize - all unspecified fields are initialized to 0. this function initializes all fields except mtime and checksum
    hdr = (HeaderPosixUstar){
        .typeflag = { '0' }, // regular file
        .magic    = "ustar",
        .version  = {'0', '0'},
        .checksum = {' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ' }, // required checksum value, while calculating the checksum (8 spaces)
    };
    
    // remove leading /
    if (fn_in_tar[0] == '/') {
        WARN_ONCE ("FYI: within the tar file, leading '/' are removed from file names%s", "");
        fn_in_tar++;
    }

    // filename: case: name is up to 99 characters 
    unsigned fn_in_tar_len = strlen (fn_in_tar);
    const char *sep;
    if (fn_in_tar_len <= 99) 
        memcpy (hdr.name, fn_in_tar, fn_in_tar_len);
    
    // filename: case: filename is longer than 99, but last component is at most 99 - split between "name" and "prefix" field (separated at a '/', and excluding the seperating '/')
    else if ((sep = strchr (&fn_in_tar[fn_in_tar_len-100], '/'))) {
        memcpy (hdr.prefix, fn_in_tar, sep - fn_in_tar);
        memcpy (hdr.name, sep+1, (&fn_in_tar[fn_in_tar_len] - (sep+1)));
    }

    // filename: case: last component is longer than 99 - use Gnu LongLink extension
    else {
        memcpy (hdr.name, fn_in_tar, 99); // fallback for extracting using non-GNU tar: nul-terminated truncated filename 
        tar_write_gnu_long_filename (fn_in_tar, fn_in_tar_len+1);
    }

    // copy mode, uid, gid, uname, gname, mtime from an existing file
    rom fn = txt_file ? txt_name : fn_on_disk/*copying an existing .genozip file*/;
    tar_copy_metadata_from_file (fn); 

    if (flag.debug_tar) 
        iprintf ("tar_open_file: t_offset=%"PRIu64" ftell=%"PRIu64" data_start=%"PRIu64" %s\n", 
                 t_offset, ftello64 (tar_file), t_offset + 512, fn);

    tar_fwrite (&hdr, 512, fn_in_tar);
    t_offset += 512; // past tar header

    return tar_file;
}

static void tar_add_hard_link (rom fn_on_disk, rom fn_in_tar_src, rom fn_in_tar_dst)
{
    ASSERTNOTNULL (tar_file);

    // initialize - all unspecified fields are initialized to 0. this function initializes all fields except mtime and checksum
    hdr = (HeaderPosixUstar){
        .typeflag = { '1' }, // hard link
        .magic    = "ustar",
        .version  = {'0', '0'},
        .checksum = {' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ' }, // required checksum value, while calculating the checksum (8 spaces)
    };
    
    memcpy (hdr.linkname, fn_in_tar_src, strlen (fn_in_tar_src));
    memcpy (hdr.name,     fn_in_tar_dst, strlen (fn_in_tar_dst));

    // copy mode, uid, gid, uname, gname, mtime from an existing file
    tar_copy_metadata_from_file (fn_on_disk);     // case: we're copying an exiting genozip file - take from that genozip file

    // header checksum
    unsigned checksum = 0;
    for (unsigned i=0; i < 512; i++) checksum += ((bytes)&hdr)[i];
    snprintf (hdr.checksum, sizeof (hdr.checksum), "%06o", checksum);

    if (flag.debug_tar) 
        iprintf ("tar_add_hard_link: t_offset=%"PRIu64" ftell=%"PRIu64" %s\n", t_offset, ftello64 (tar_file), hdr.name);

    tar_fwrite (&hdr, 512, fn_in_tar_dst);
    t_offset += 512; // past tar header
}

bool tar_zip_is_tar (void)
{
    return IS_ZIP && !!tar_name;
}

int64_t tar_file_offset (void) 
{
    return t_offset; // 0 if not using tar
}

void tar_close_file (void **file)
{
    ASSERTNOTNULL (tar_file);

    int64_t tar_size = ftello64 (tar_file);
    ASSERT (tar_size >= 0, "ftello64 failed for %s", tar_name);

    int64_t z_size = tar_size - t_offset;

    // file consist of full 512-byte records. pad it to make it so:
    int64_t padding_len = ROUNDUP512(tar_size) - tar_size;

    if (padding_len) {
        tar_fwrite ((char[512]){}, padding_len, "file padding");        
        tar_size += padding_len; 
    }

    // size - case < 8GB store in nul-terminated octal in ASCII (33 bits = 11 octal numerals x 3 bits each) - standard tar
    if (z_size < (1ULL << 33))
        snprintf (hdr.size.std, sizeof (hdr.size.std), "%.*"PRIo64, (int)sizeof(hdr.size)-1, z_size);
    
    // case >= 8GB: store as binary
    else { 
        hdr.size.gnu.marker = 0x80; 
        uint64_t size_bgen = BGEN64 (z_size);
        memcpy (hdr.size.gnu.size_bgen, &size_bgen, 8);
    }

    // header checksum
    unsigned checksum = 0;
    for (unsigned i=0; i < 512; i++) checksum += ((bytes)&hdr)[i];
    snprintf (hdr.checksum, sizeof (hdr.checksum), "%06o", checksum);
    
    // update header
    ASSERT (!fseeko64 (tar_file, t_offset-512, SEEK_SET), "fseek(%"PRId64") of %s failed (1): %s", t_offset-512, tar_name, strerror (errno));
    tar_fwrite (&hdr, 512, "file header");
    ASSERT (!fseeko64 (tar_file, 0, SEEK_END), "fseek(END) of %s failed (2): %s", tar_name, strerror (errno));

    // flush to finalize z_file within tar file, before deleting txt file, and also before spawning a process to test it
    if (flag.replace || flag.test) 
        fflush (tar_file); 

    t_offset = tar_size; // next file start offset

    if (flag.debug_tar) 
        iprintf ("tar_close_file: t_offset=%"PRIu64" ftell=%"PRIu64"\n", t_offset, ftello64 (tar_file));

    if (file) *file = NULL;
}

void tar_copy_file (rom fn_on_disk, rom fn_in_tar)
{
    ASSERTNOTNULL (tar_file);

    if (flag.debug_tar) 
        iprintf ("tar_copy_file: %s\n", fn_in_tar);

    tar_open_file (fn_on_disk, fn_in_tar);

    FILE *src_file = fopen (fn_on_disk, "rb");

    #define BLOCK_SIZE (1 MB)
    char *data = MALLOC (BLOCK_SIZE);

    int64_t size;
    int64_t bytes_copied = 0;
    while ((size = fread (data, 1, BLOCK_SIZE, src_file))) {
        tar_fwrite (data, size, fn_on_disk);
        bytes_copied += size;
    }

    FREE (data);
    FCLOSE (src_file, fn_on_disk);

    tar_close_file(0);

    ASSERT (bytes_copied == file_get_size (fn_on_disk), "File %s has size %"PRId64" but copied only %"PRId64" bytes to the tar file", 
            fn_on_disk, file_get_size (fn_on_disk), bytes_copied);
}

void tar_finalize (void)
{
    ASSERTNOTNULL (tar_file);

    if (flag.debug_tar) 
        iprintf ("tar_finalize EOF block: t_offset=%"PRIu64" ftell=%"PRIu64"\n", t_offset, ftello64 (tar_file));

    // tar file format: two empty tar blocks as EOF
    tar_fwrite ((char[1024]){}, 1024, "EOF tar blocks");

    FCLOSE (tar_file, "tar_file");
    t_offset = 0;
}
