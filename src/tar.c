// ------------------------------------------------------------------
//   tar.c
//   Copyright (C) 2021-2023 Genozip Limited. Patent pending.
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
    char linkname[100]; // we don't support hard links, so this is always 0
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
static int64_t file_offset = 0; // after tar_open_file: file start offset within tar (after the file tar header)
static HeaderPosixUstar hdr = {};

void tar_set_tar_name (rom tar_filename)
{
    tar_name = tar_filename;
}

void tar_initialize (void)
{
    ASSERTNOTNULL (tar_name);

    ASSINP (flag.force || !file_exists (tar_name), "file %s already exists, use --force (or -f) to overwrite", tar_name);

    tar_file = fopen (tar_name, flag.pair ? "wb+" : "wb"); // if --pair, when compressing pair2, we go back and read pair1
    ASSINP (tar_file, "cannot create tar file %s: %s", tar_name, strerror (errno)); 
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
    sprintf (hdr.uid,   "%.*o", (int)sizeof(hdr.uid)-1,  st.st_uid);
    sprintf (hdr.gid,   "%.*o", (int)sizeof(hdr.gid)-1,  st.st_gid);
    sprintf (hdr.mode,  "%.*o", (int)sizeof(hdr.mode)-1, st.st_mode);
    sprintf (hdr.mtime, "%.*"PRIo64, (int)sizeof(hdr.mtime)-1, (uint64_t)st.st_mtime); // mtime is 64b on Windows and Linux, 32b on MacOS

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

// filenames that have a last component longer than 99 characters don't fit in POSIX tar. Instead, we use a GNU-specific extension
// of storing the filename as a pseudo-file in the tarball. See: https://itecnote.com/tecnote/r-what-exactly-is-the-gnu-tar-longlink-trick/
static void tar_write_gnu_long_filename (const char *z_fn, unsigned z_fn_len/* including \0 */)
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
    sprintf (ll_hdr.size.std, "%o", z_fn_len); // +1 for \0

    // header checksum
    unsigned checksum = 0;
    for (unsigned i=0; i < 512; i++) checksum += ((unsigned char*)&ll_hdr)[i];
    sprintf (ll_hdr.checksum, "%06o", checksum);

    ASSERT (fwrite (&ll_hdr, 512,   1, tar_file) == 1, "failed to write LongLink header of %s to %s", z_fn, tar_name); 
    ASSERT (fwrite (z_fn, z_fn_len, 1, tar_file) == 1, "failed to write long filename of %s to %s", z_fn, tar_name); 

    // pad to full block
    if (z_fn_len % 512) {
        char padding[512] = "";
        ASSERT (fwrite (padding, ROUNDUP512(z_fn_len) - z_fn_len, 1, tar_file) == 1, "failed to write long filename padding of %s to %s", z_fn, tar_name); 
    }

    file_offset += 512 + ROUNDUP512(z_fn_len);
}

// open z_file within tar, for writing
FILE *tar_open_file (rom z_fn)
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
    if (z_fn[0] == '/') {
        WARN_ONCE ("FYI: within the tar file, leading '/' are removed from file names%s", "");
        z_fn++;
    }

    // filename: case: name is up to 99 characters 
    unsigned z_fn_len = strlen (z_fn);
    const char *sep;
    if (z_fn_len <= 99) 
        memcpy (hdr.name, z_fn, z_fn_len);
    
    // filename: case: filename is longer than 99, but last component is at most 99 - split between "name" and "prefix" field (separated at a '/', and excluding the seperating '/')
    else if ((sep = strchr (&z_fn[z_fn_len-100], '/'))) {
        memcpy (hdr.prefix, z_fn, sep - z_fn);
        memcpy (hdr.name, sep+1, (&z_fn[z_fn_len] - (sep+1)));
    }

    // filename: case: last component is longer than 99 - use Gnu LongLink extension
    else {
        memcpy (hdr.name, z_fn, 99); // fallback for extracting using non-GNU tar: nul-terminated truncated filename 
        tar_write_gnu_long_filename (z_fn, z_fn_len+1);
    }

    // copy mode, uid, gid, uname, gname, mtime from an existing file
    if (!txt_file) tar_copy_metadata_from_file (z_fn);     // case: we're copying an exiting genozip file - take from that genozip file
    else           tar_copy_metadata_from_file (txt_name); // case: zipping a txt_file on disk - take from that txt file

    ASSERT (fwrite (&hdr, 512, 1, tar_file) == 1, "failed to write header of %s to %s", z_fn, tar_name); // place holder - we will update this upon close
    file_offset += 512; // past tar header

    return tar_file;
}

bool tar_is_tar (void)
{
    return !!tar_name;
}

int64_t tar_file_offset (void) 
{
    return file_offset; // 0 if not using tar
}

void tar_close_file (void **file)
{
    ASSERTNOTNULL (tar_file);

    if (flag.replace) fflush (tar_file); // flush z_file before deleting txt file

    int64_t tar_size = ftello64 (tar_file);
    ASSERT (tar_size >= 0, "ftello64 failed for %s", tar_name);

    int64_t z_size = tar_size - file_offset;

    // file consist of full 512-byte records. pad it to make it so:
    int64_t padding_len = ROUNDUP512(tar_size) - tar_size;

    if (padding_len) {
        char padding[512] = "";
        ASSERT (fwrite (padding, padding_len, 1, tar_file) == 1, "failed to write file padding to %s", tar_name);
        
        tar_size += padding_len; 
    }

    // size - case < 8GB store in nul-terminated octal in ASCII (33 bits = 11 octal numerals x 3 bits each) - standard tar
    if (z_size < (1ULL << 33))
        sprintf (hdr.size.std, "%.*"PRIo64, (int)sizeof(hdr.size)-1, z_size);
    
    // case >= 8GB: store as binary
    else { 
        hdr.size.gnu.marker = 0x80; 
        uint64_t size_bgen = BGEN64 (z_size);
        memcpy (hdr.size.gnu.size_bgen, &size_bgen, 8);
    }

    // header checksum
    unsigned checksum = 0;
    for (unsigned i=0; i < 512; i++) checksum += ((unsigned char*)&hdr)[i];
    sprintf (hdr.checksum, "%06o", checksum);
    
    // update header
    ASSERT (!fseeko64 (tar_file, file_offset-512, SEEK_SET), "fseek(%"PRId64") of %s failed (1): %s", file_offset-512, tar_name, strerror (errno));
    ASSERT (fwrite (&hdr, 512, 1, tar_file) == 1, "failed to write file header to %s", tar_name);
    ASSERT (!fseeko64 (tar_file, 0, SEEK_END), "fseek(END) of %s failed (2): %s", tar_name, strerror (errno));

    file_offset = tar_size; // next file start offset

    if (file) *file = NULL;
}

void tar_copy_file (rom z_fn)
{
    ASSERTNOTNULL (tar_file);

    tar_open_file (z_fn);

    FILE *src_file = fopen (z_fn, "rb");

    #define BLOCK_SIZE (1 << 20)
    char *data = MALLOC (BLOCK_SIZE);

    int64_t size;
    int64_t bytes_copied = 0;
    while ((size = fread (data, 1, BLOCK_SIZE, src_file))) {
        ASSERT (fwrite (data, 1, size, tar_file) == size, "failed to copy %s to %s - failed to write %"PRId64" bytes", z_fn, tar_name, size);
        bytes_copied += size;
    }

    FREE (data);
    FCLOSE (src_file, z_fn);

    tar_close_file(0);

    ASSERT (bytes_copied == file_get_size (z_fn), "File %s has size %"PRId64" but copied only %"PRId64" bytes to the tar file", 
            z_fn, file_get_size (z_fn), bytes_copied);
}

void tar_finalize (void)
{
    ASSERTNOTNULL (tar_file);

    // tar file format: two empty tar blocks as EOF
    char s[1024] = "";
    ASSERT (fwrite (s, 1024, 1, tar_file) == 1, "failed to EOF tar blocks to %s", tar_name);

    FCLOSE (tar_file, "tar_file");
    file_offset = 0;
}
