#define GENOZIP_CODE_VERSION "15.0.86"
// ⇑ MUST be first line, analyzed by scripts ⇑

#pragma once

typedef struct __attribute__ ((packed)) {
    uint8_t major;
    uint16_t minor; // up to 16383 - limited by GenozipHeader.genozip_minor_ver. Populated for files since 15.0.28
} Version;
#define NO_VERSION ((Version){})

#define VER_GE(a, b) (((a).minor >= (b).minor && (a).major == (b).major) || (a).major > (b).major)
#define VER_GE_(a, major, minor) VER_GE((a), ((Version){(major), (minor)}))

extern bool version_is_devel (void);
extern Version code_version (void);
extern Version file_version (void);
extern Version str_to_version (const char *version_str); // can't use rom due to dependencies
extern StrText STRver (Version v);
static inline StrText STRver_(int major, int minor) { return STRver ((Version){ major, minor }); }
extern StrText version_str (void);
extern void version_background_test_for_newer (void);
extern void version_print_notice_if_has_newer (void);
extern const char *genozip_update_msg (void);
