#define GENOZIP_CODE_VERSION "15.0.32"

extern int code_version_major (void);
extern int code_version_minor (void);
extern StrText version_str (void);
extern void version_background_test_for_newer (void);
extern void version_print_notice_if_has_newer (void);
extern const char *genozip_update_msg (void);
