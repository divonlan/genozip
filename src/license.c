// ------------------------------------------------------------------
//   license.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#else 
#include <sys/ipc.h>
#include <sys/sem.h>
#include <fcntl.h>
#endif
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "genozip.h"
#include "strings.h"
#include "arch.h"
#include "profiler.h" // for TimeSpecType
#include "url.h"
#include "license.h"
#include "buffer.h"
#include "flags.h"
#include "md5.h"
#include "website.h"
#include "version.h"
#include "file.h"

// these field names appear in the license file starting V12.0.7. They are also consumed in register.sh.
#define LIC_FIELD_TYPE         "License type"
#define LIC_FIELD_VERSION      "Genozip license version"
#define LIC_FIELD_INSTITUTION  "License granted to"
#define LIC_FIELD_NAME         "Accepted by (name)"
#define LIC_FIELD_EMAIL        "Accepted by (email)"
#define LIC_FIELD_MACHINE_TIME "Machine time"
#define LIC_FIELD_TIMESTAMP    "Timestamp of acceptance"
#define LIC_FIELD_IP           "IP address of acceptance"
#define LIC_FIELD_NUMBER       "License number"
#define LIC_FIELD_ALLOW_STATS  "Allow_stats"

// IF YOU'RE CONSIDERING MODIFYING THIS CODE TO BYPASS LICENSE LIMITATIONS, DON'T! It would be a violation of the license,
// and might put you personally as well as your organization at legal and financial risk - see "Severly unauthorized use of Genozip"
// section of the license. Rather, please contact sales@genozip.com to discuss which license would be appropriate for your case.
#define EVAL_NUM_FILES 100 // if updating, also update the web page get-genozip 
#define EVAL_NUM_FILES_STR "100"
#define DEEP_or_PAIR_NUM_FILES 10

#include "text_license.h"

// note: value of lic_type is written to license file, so types might be added, but existing types should not change values
typedef enum __attribute__ ((__packed__))/*1 byte*/ { LIC_TYPE_NONE, LIC_TYPE_ACADEMIC, LIC_TYPE_EVAL,       LIC_TYPE_STANDARD_LEGACY/*to v14*/, LIC_TYPE_STANDARD/*since v15*/, LIC_TYPE_ENTERPRISE, LIC_TYPE_FIRST, NUM_LIC_TYPES } LicenseType; 
static rom lic_types[NUM_LIC_TYPES] =               { "",            "Academic",        "30-day evaluation", "Standard(Legacy)",                 "Standard",                     "Enterprise",        "First"                       }; 
// note: these strings ^ are referred to in register.sh

static struct {
    bool initialized;
    LicenseType lic_type; 
    char name[256], institution[1024], email[256], ip[ARCH_IP_LEN], version[20], allow_stats[4];
    StrTime timestamp;
    int64_t machine_time; // timestamp expressed as seconds since epoch
    uint32_t license_num;
} rec = {};

// returns new counter value, or 0 if counter doesn't work
static int counter_increment (rom filename, int inc /*may be negative*/) 
{
    #define MAX_COUNTER 0xfff  // 12 bit counter

#ifdef __linux__
    #define COUNTER_MAGIC 270512 // 20 bit magic
    #define IS_MAGICAL(st) ((st.st_mtim.tv_nsec & 0xfffff) == COUNTER_MAGIC)
    #define COUNTST(st) (st.st_mtim.tv_nsec >> 20)
    #define NSEC(counter) (((counter) << 20) | COUNTER_MAGIC)

    struct stat st;
    if (stat (filename, &st)) return 0; // fail silently 

    int old_count = IS_MAGICAL(st) ? COUNTST(st) : 0;

    int new_count = MIN_(MAX_(old_count + inc, 0), MAX_COUNTER);
    
    int ret = utimensat (0/*ignored*/, filename, (const struct timespec[]){ { .tv_nsec = UTIME_OMIT         }, 
                                                                            { .tv_sec  = st.st_mtim.tv_sec,
                                                                              .tv_nsec = NSEC(new_count)    } }, 0); 
    // verify that it worked, before returning new_count
    return (ret == 0 && !stat (filename, &st) && IS_MAGICAL(st) && COUNTST(st) == new_count) ? new_count : 0;

#elif _WIN32
    #define COUNTER_MAGIC 0x125 // 11 bit magic
    #define IS_MAGICAL(mtime) ((mtime.dwLowDateTime & 0x7ff) == COUNTER_MAGIC)
    #define COUNTST(mtime) ((mtime.dwLowDateTime >> 11) & 0xfff)

    HANDLE h = CreateFileA (filename, GENERIC_READ | FILE_WRITE_ATTRIBUTES, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (h == INVALID_HANDLE_VALUE) { // very unlikely, but might legitimately happen if two genozip processes are in here concurrently
        ASSERT (!flag.debug, "CreateFileA failed: %s", str_win_error());
        return 0;
    }

    FILETIME mtime; // 64 bits, in 100ns granularity. We use the 23 LSb.
    if (!GetFileTime (h, NULL, NULL, &mtime)) {
        ASSERT (!flag.debug, "GetFileTime failed: %s", str_win_error());
        goto fail;
    }

    int old_count = IS_MAGICAL(mtime) ? COUNTST(mtime) : 0;

    int new_count = MIN_(MAX_(old_count + inc, 0), MAX_COUNTER);

    mtime.dwLowDateTime = (mtime.dwLowDateTime & 0xff800000/*9 MSb*/) | (new_count << 11) | COUNTER_MAGIC;
    if (!SetFileTime (h, NULL, NULL, &mtime)) {
        ASSERT (!flag.debug, "SetFileTime failed: %s", str_win_error());
        goto fail;
    }

    // verify that it worked, before returning new_count
    if (GetFileTime (h, NULL, NULL, &mtime) && IS_MAGICAL(mtime) && COUNTST(mtime) == new_count) {
        CloseHandle (h);
        return new_count;
    }

fail:
    CloseHandle (h);
    return 0;

#else
    // TO DO: Support for MacOS (bug 810)
    return 0;
#endif
}

static bool counter_has_reached (rom filename, uint32_t reached_this)
{
    return counter_increment (filename, 0) >= reached_this;
}

static uint32_t license_calc_number (ConstBufferP license_data)
{
    char data_no_ws[license_data->len];
    unsigned data_no_ws_len = str_remove_whitespace (STRb(*license_data), data_no_ws);        

    return md5_do (data_no_ws, data_no_ws_len).words[0];
}

static void license_generate (BufferP license_data)
{
    for (int i=0; i < ARRAY_LEN(license); i++) {
        buf_append_string (evb, license_data, license[i]); // allocs one extra char
        BNXTc (*license_data) = '\n';
    }

    bufprintf (evb, license_data,     // note: the license data includes the Genozip version
               LIC_FIELD_TYPE": %d\n" // added v14
               LIC_FIELD_INSTITUTION": %s\n"
               LIC_FIELD_NAME": %s\n"
               LIC_FIELD_EMAIL": %s\n"
               LIC_FIELD_MACHINE_TIME": %"PRIu64"\n"
               LIC_FIELD_TIMESTAMP": %s\n"
               LIC_FIELD_IP": %s\n"
               LIC_FIELD_ALLOW_STATS": %s\n", 
               rec.lic_type, rec.institution, rec.name, rec.email, rec.machine_time, rec.timestamp.s, rec.ip,
               rec.allow_stats);
    
    rec.initialized = true;
    rec.license_num = license_calc_number (license_data);
    strcpy (rec.version, GENOZIP_CODE_VERSION);

    bufprintf (evb, license_data, LIC_FIELD_NUMBER": %u\n", rec.license_num);
}

void license_set_filename (rom filename)
{
    struct stat sb;
    ASSINP (!stat (filename, &sb), "Failed to access license file %s: %s", filename, strerror (errno));

    flag.license_filename = filename;
}

static rom license_get_filename (bool create_folder_if_needed)
{
    if (flag.license_filename) return flag.license_filename; // non-standard filename set with --licfile

#ifdef _WIN32
    ASSINP0 (getenv ("APPDATA"), "cannot store license, because APPDATA env var is not defined");

    char folder[500];
    sprintf (folder, "%s/genozip", getenv ("APPDATA"));

    if (create_folder_if_needed) {
        int ret = _mkdir (folder); 
        ASSERT (ret >= 0 || errno == EEXIST, "failed to created the folder %s", folder);
    }

#else
    rom folder = getenv ("HOME");
    ASSINP0 (folder, "cannot calculate license file name, because $HOME env var is not defined");
#endif    

    static char filename[1024]; // avoid MALLOC so we don't need to leak memory
    ASSERT (strlen (folder) + STRLEN("/.genozip_license") + 1 < sizeof (filename), "Directory name %s too long for license file", folder);
    sprintf (filename, "%s/.genozip_license", folder);

    return filename;
}

static rom license_load_field (rom field, STRps(line))
{
    unsigned field_len = strlen (field);

    for (int i=n_lines-1; i >= 0; i--)
        if (line_lens[i] > field_len+2 && !memcmp (lines[i], field, field_len) && lines[i][field_len] == ':' && lines[i][field_len+1] == ' ')
            return &lines[i][field_len+2];

    return ""; // not found
}

bool license_is_registered (void)
{
    rom filename = license_get_filename (true);
    return file_exists (filename);
}

bool license_is_eval (void)
{
    return rec.lic_type == LIC_TYPE_EVAL;
}

// IF YOU'RE CONSIDERING MANIPULATING SEMAPHORES TO BYPASS THE REGISTRTION, DON'T! It would be a violation of the license,
// and might put you personally as well as your organization at legal and financial risk - see "Severly unauthorized use of Genozip"
// section of the license. Rather, please contact sales@genozip.com to discuss which license would be appropriate for your case.
static bool license_is_first (void)
{
#ifndef _WIN32
    rom home = getenv ("HOME");
    if (!home) return false;

    key_t key = ftok (home, 20010802);
    int sem = semget (key, 1, IPC_CREAT | IPC_EXCL | 0600); // note: in Linux, a new semaphone is initialized to 0. This is not guaranteed to be true on other platforms.
    if (sem >= 0) {
        rec.lic_type = LIC_TYPE_FIRST;
        rec.initialized = true;
        return true;
    }

    else if (errno == EEXIST) 
        return false;

    else {
        ASSERT (!flag.debug, "Error in semget: %s", strerror (errno));
        return false; 
    }
    
    return true;

#else
    rom exec_path = arch_get_executable (SOFT_FAIL);
    struct stat st;

    if (exec_path && !stat (exec_path, &st) && (time(NULL) - st.st_ctime > 1000000)) 
        return false;

    rom lic_filename = license_get_filename(true);
    
    char first_fn[strlen(lic_filename) + 10];
    sprintf (first_fn, "%s.first", lic_filename);

    if (file_exists (first_fn)) return false;

    FILE *fp = fopen (first_fn, "w");
    fclose (fp);

    return file_exists (first_fn);

#endif
}

// IF YOU'RE CONSIDERING TAMPERING WITH THIS CODE TO BYPASS THE REGISTRTION, DON'T! It would be a violation of the license,
// and might put you personally as well as your organization at legal and financial risk - see "Severly unauthorized use of Genozip"
// section of the license. Rather, please contact sales@genozip.com to discuss which license would be appropriate for your case.
void license_load (void)
{
    if (rec.initialized) return;

    rom filename = license_get_filename (true);
 
    if (!file_exists (filename)) {
        if (IS_ZIP && !license_is_first()) 
            license_register ();

        return;
    }

    file_split_lines (filename, "license", false/*we allow UTF-8*/);
    
    char license_num_str[30] = "", lic_type_str[16]="", machine_time_str[24]="";
    #define COPY_FIELD(var,field) strncpy (var, license_load_field (field, STRas(line)), sizeof (var)-1)

    COPY_FIELD (lic_type_str, LIC_FIELD_TYPE);  // added v14

    // note: licenses prior to v14 don't have this field
    if (!str_get_int_range8 (lic_type_str, strlen (lic_type_str), 1, NUM_LIC_TYPES-1, &rec.lic_type)) 
        rec.lic_type = LIC_TYPE_ACADEMIC; // fallback 

    COPY_FIELD (rec.version,      LIC_FIELD_VERSION);
    COPY_FIELD (rec.institution,  LIC_FIELD_INSTITUTION);
    COPY_FIELD (rec.name,         LIC_FIELD_NAME);
    COPY_FIELD (rec.email,        LIC_FIELD_EMAIL);
    COPY_FIELD (machine_time_str, LIC_FIELD_MACHINE_TIME);
    COPY_FIELD (rec.timestamp.s,  LIC_FIELD_TIMESTAMP);
    COPY_FIELD (rec.ip,           LIC_FIELD_IP);
    COPY_FIELD (rec.allow_stats,  LIC_FIELD_ALLOW_STATS);
    COPY_FIELD (license_num_str,  LIC_FIELD_NUMBER);

    if (!str_get_uint32 (license_num_str, strlen (license_num_str), &rec.license_num)) {
        WARN0 ("old format license");
        goto reregister;
    }

    if (!str_get_int (machine_time_str, strlen (machine_time_str), &rec.machine_time)) {
        WARN0 ("no machine time");
        goto reregister;
    }

    data.len -= line_lens[n_lines-1] + 2;
    if (rec.license_num != license_calc_number (&data)) {
        WARN ("license verification failed (%u, %u)", rec.license_num, license_calc_number (&data));
        goto reregister;
    }

    if (rec.lic_type == LIC_TYPE_EVAL) {
        #define BUY "To purchase a Standard or Enterprise License: " WEBSITE_BUY " or contact " EMAIL_SALES "\n"
        ASSINP0 (time(0) - rec.machine_time < (30*24*60*60),
                 "You reached the end of the 30 evaluation period of Evaluation License.\n" BUY);

        ASSINP (!counter_has_reached (filename, EVAL_NUM_FILES),
                "You reached the maximum number of files (%u) compressible with the Evaluation License.\n" BUY, EVAL_NUM_FILES);
    }

    else if (rec.lic_type == LIC_TYPE_STANDARD_LEGACY) 
        ASSINP (!flag.deep || !counter_has_reached (filename, DEEP_or_PAIR_NUM_FILES),
                "You reached %u --deep compressions, which is the maximum number granted with the %s License.\n"
                "To upgrade to a Enterprise License: " WEBSITE_BUY " or contact " EMAIL_SALES "\n", DEEP_or_PAIR_NUM_FILES, lic_types[rec.lic_type]);

    else if (rec.lic_type == LIC_TYPE_STANDARD) 
        ASSINP ((!flag.deep && !flag.pair) || !counter_has_reached (filename, DEEP_or_PAIR_NUM_FILES),
                "You reached %u --deep / --pair compressions, which is the maximum number granted with the %s License.\n"
                "To upgrade to a Enterprise License: " WEBSITE_BUY " or contact " EMAIL_SALES "\n", DEEP_or_PAIR_NUM_FILES, lic_types[rec.lic_type]);

    rec.initialized = true;

    buf_destroy (data);
    return;

reregister:
    file_remove (filename, true);

    // if stdin or stderr is redirected - we cannot start an interactive registration flow
    ASSINP0 (isatty(0) && isatty(2), "Genozip license terms & conditions have changed, please re-register by running: genozip --register");

    fprintf (stderr, "Genozip license terms & conditions have changed, please re-register:\n\n");
    license_register();
}

static bool license_submit (rom os, unsigned cores, bool is_script, rom user_host, rom dist)
{
    // reference: https://stackoverflow.com/questions/18073971/http-post-to-a-google-form/47444396#47444396

    // FORM_ID is in the url when you preview your form
    #define PREFIX "https://docs.google.com/forms/d/e/1FAIpQLSc6pSBIOBsS5Pu-JNvfnLWV2Z1W7k-4f2pKRo5rTbiU_sCnIw/formResponse"
    
    /* To get entry IDs - in Chrome browser: 1. open form 2. click on eye icon to Preview 2. right-click Inspect 3. go to "console" tab 4. run this code:
    function loop(e){
    if(e.children)
    for(let i=0;i<e.children.length;i++){
        let c = e.children[i], n = c.getAttribute('name');
        if(n) console.log(`${c.getAttribute('aria-label')}: ${n}`);
        loop(e.children[i]);
     }
    }; loop(document.body);
    */

    // note: identical to register.sh
    char *url_format = PREFIX
                       "?entry.344252538=%s"
                       "&entry.926671216=%s"
                       "&entry.1734045469=%s"
                       "&entry.2009586582=%s"
                       "&entry.119966790=%s"  
                       "&entry.81542373=%s"
                       "&entry.1668073218=%u"
                       "&entry.1943454647=%s"
                       "&entry.1763961212=%s"
                       "&entry.1655649315=%u"
                       "&entry.186159495=%s"
                       "&entry.1598028195=%s"
                       "&entry.1384715202=%s";

    char *institutionE = url_esc_all_or_none (rec.institution);
    char *nameE        = url_esc_all_or_none (rec.name);
    char *emailE       = url_esc_all_or_none (rec.email);
    char *lic_typeE    = url_esc_non_valid_chars (lic_types[rec.lic_type]);
    char *osE          = url_esc_non_valid_chars (os);
    char *user_hostE   = url_esc_non_valid_chars (user_host);

    char url[sizeof (rec)*3 + 200];
    sprintf (url, url_format, institutionE, nameE, emailE, lic_typeE, rec.allow_stats, osE, cores, rec.ip, user_hostE, rec.license_num, rec.version, dist, is_script ? "script" : "");

    bool success = url_read_string (url, NULL, 0) >= 0;
    
    FREE (institutionE); FREE (nameE); FREE (emailE); FREE (osE); FREE (user_hostE);
    return success;
}

static bool license_verify_email (STRc(response), rom unused)
{
    rom domain = strchr (response, '@');
    if (!domain) return false;

    domain++;
    int domain_len = strlen (domain);

    if (domain_len < 5 || !strchr (response, '.')) return false;

    static rom disposable_email_domains[] = { // https://gist.github.com/nocturnalgeek/1b8fa44283314544c487
        "mailinator.com", "binkmail.com", "bobmail.info", "chammy.info", "devnullmail.com",
        "letthemeatspam.com", "mailinater.com", "mailinator.net", "mailinator2.com", "notmailinator.com",
        "reallymymail.com", "reconmail.com", "safetymail.info", "sendspamhere.com", "sogetthis.com",
        "spambooger.com", "spamherelots.com", "spamhereplease.com", "spamthisplease.com",
        "streetwisemail.com", "suremail.info", "thisisnotmyrealemail.com", "tradermail.info",
        "veryrealemail.com", "zippymail.info", "keshitv.com", "jollyfree.com", "bbitj.com"
    };

    for (int i=0; i < ARRAY_LEN(disposable_email_domains); i++)
        if (!strcmp (domain, disposable_email_domains[i])) {   
            fprintf (stderr, "\nHey, we worked very hard to build this software. Please respect us by providing your real email address.\n");
            return false;
        }

    return true;
}

static bool license_verify_code (STRc(response), rom unused)
{
    if (response_len != 6 || !str_is_numeric (STRa(response))) {
        fprintf (stderr, "\nExpecting a 6 digit number\n");
        return false;
    }

    Digest digest = md5_do (STRa(response));

    if (digest.bytes[0]!='2' && digest.bytes[1] >= 7) {
        fprintf (stderr, "\nIncorrect verification code\n");
        return false;
    }

    return true;
}

static bool license_is_consumer_email (rom email)
{
    rom domain = strchr (email, '@') + 1;
    int domain_len = strlen (domain);

    if (domain_len > 5 && !memcmp (domain, "live.",    5)) return true; // "live.com", "live.cn", "live.fr", etc.
    if (domain_len > 6 && !memcmp (domain, "yahoo.",   6)) return true; // "yahoo.com", "yahoo.com.br", "yahoo.fr", "yahoo.co.uk", etc.
    if (domain_len > 8 && !memcmp (domain, "hotmail.", 8)) return true; // "hotmail.com", "hotmail.fr", "hotmail.co.jp", etc

    static rom consumer_domains[] = { // except multi-domainers: yahoo, live, hotmail
        "gmail.com", "googlemail.com", "outlook.com", "protonmail.com", "msn.com", 
        "rediffmail.com", "ymail.com", "icloud.com",                    // Global
        "aol.com", "comcast.net",                                       // US
        "qq.com", "163.com",                                            // China 
        "web.de", "gmx.de",                                             // Germany 
        "wanadoo.fr", "orange.fr", "free.fr",                           // France
        "hanmail.net", "naver.com",                                     // Korea
        "list.ru", "mail.ru", "inbox.ru", "yandex.ru", "rambler.ru",    // Russia 
        "yaani.com",                                                    // Turkey
        "libero.it",                                                    // Italy
        "uol.com.br", "bol.com.br",                                     // Brazil
    };

    for (int i=0; i < ARRAY_LEN(consumer_domains); i++)
        if (!strcmp (domain, consumer_domains[i])) return true;

    return false;
}

static bool license_verify_name (STRc(response), rom unused)
{
    if (!strchr (response, ' ')) {
        fprintf (stderr, "Please enter your full name\n");
        return false;
    }
    
    return true;
}

static bool license_verify_license (STRc(response), rom unused)
{
    return strlen (response) == 1 && (*response >= '1' && *response <= '4');
}

static void license_exit_if_not_confirmed (rom query, DefAnswerType def_answer)
{
    if (!str_query_user_yn (query, def_answer)) {
        fprintf (stderr, "\nYou have not registered. You may register at any time in the future.\n\nWishing you a wonderful day from the Genozip team! https://genozip.com\n");
        exit_ok();
    }
}

// UI flow to generate a license registration

// IF YOU'RE CONSIDERING TAMPERING WITH THIS CODE TO BYPASS THE REGISTRTION, DON'T! It would be a violation of the license,
// and might put you personally as well as your organization at legal and financial risk - see "Severly unauthorized use of Genozip"
// section of the license. Rather, please contact sales@genozip.com to discuss which license would be appropriate for your case.
void license_register (void)
{
    char lic_type[100];
    rom os, dist, user_host;
    unsigned cores;
    bool is_script;

    if (!flag.do_register) flag.do_register = "";

    str_split (flag.do_register, strlen (flag.do_register), 11, '|', field, true);
    str_nul_separate (field);

    // if stdin or stderr is redirected - we cannot ask the user an interactive question
    ASSINP0 (isatty(0) && isatty(2), "Use of Genozip is free for academic purposes, but requires registration. Please run: genozip --register.\n"
                                     "If you are unable to register (for example because this is a batch-job machine) please see: " WEBSITE_USING_ON_HPC);

    rom filename = license_get_filename (true);

    if (!n_fields) {

        fprintf (stderr, "Welcome to Genozip!\n\n"
                         "- Genozip is a commercial product, however it is free for academic research use with some limitations.\n"
                         "  To check eligibility see: " WEBSITE_GET_GENOZIP "\n\n"
                         "- For non-academic use, you may use Genozip for free for a 30 days evaluation period.\n\n");

        if (file_exists (filename)) 
            license_exit_if_not_confirmed ("You are already registered. Are you sure you want to re-register again?", QDEF_NONE);
        else 
            license_exit_if_not_confirmed ("Would you like to register now?", QDEF_YES);
    }

    file_remove (filename, true); // remove old license, if one exists
    
    if (n_fields) { // fields correspond to register.sh
        rec.lic_type = atoi (fields[0]);
        strncpy (rec.institution, fields[1], sizeof(rec.institution)-1);
        strncpy (rec.name,        fields[2], sizeof(rec.name)-1);
        strncpy (rec.email,       fields[3], sizeof(rec.email)-1);
        strncpy (rec.ip,          fields[4], sizeof(rec.ip)-1);
        strncpy (rec.allow_stats, fields[5], sizeof(rec.allow_stats)-1);
        os        = fields[6];
        dist      = fields[7]; 
        user_host = fields[9];
        cores     = atoi(fields[10]);
        is_script = true;
    }
    else {
        fprintf (stderr, "\nLicense details -\n");
    
        str_query_user ("\nInstitution / Company name: ", rec.institution, sizeof(rec.institution), false, NULL, NULL);

        str_query_user ("\nYour name: ", rec.name, sizeof(rec.name), false, license_verify_name, NULL);
        
        str_query_user ("\nYour email address (we will send your verification code to this email):\n", rec.email, sizeof(rec.email), false, license_verify_email, NULL);
        int len = strlen (rec.email);
        
        if (license_is_consumer_email (rec.email) && len < sizeof (rec.email)-20) {
            rec.email[len] = ' ';
            fprintf (stderr, "\nHmm... that looks like a personal email address. Please enter your email address at your institution / company\n");
            
            str_query_user ("\nYour email address: ",&rec.email[len+1], sizeof(rec.email)-len-1, false, license_verify_email, NULL);
        }

        str_query_user ("\nWhat type of license do you require? (see: " WEBSITE_PRICING_FAQ ")\n\n"
                        "1. Academic License (free): Free for officially recognized research institutions (excluding data obtained commercially)\n\n"
                        "2. Evaluation License (free): Free use for 30-day (limited to " EVAL_NUM_FILES_STR " files)\n\n"
                        "3. Standard License (paid): I have already paid for a Standard License\n\n"
                        "4. Enterprise License (paid): I have already paid for a Enterprise License\n\n"
                        "Remember your Mom taught you to be honest!\n\n"
                        "Please enter 1, 2, 3 or 4: ",
                        lic_type, sizeof(lic_type), false, license_verify_license, NULL);
    
        int n = lic_type[0] - '0';
        rec.lic_type = n==1?LIC_TYPE_ACADEMIC : n==2?LIC_TYPE_EVAL : n==3?LIC_TYPE_STANDARD : LIC_TYPE_ENTERPRISE;  
    
        if (rec.lic_type == LIC_TYPE_STANDARD || rec.lic_type == LIC_TYPE_ENTERPRISE) {
            bool stats_consent = str_query_user_yn ("\nGenozip optionally logs aggregate statistics and metadata on the Genozip server,\n"
                                                    "helping us provide you with technical support if needed, and also improve our compression algorithms\n"
                                                    "(see "WEBSITE_LOGS"). May we have your permission for this? ", QDEF_NONE);
            strcpy (rec.allow_stats, stats_consent ? "Yes" : "No");
        }
        else
            strcpy (rec.allow_stats, "Yes");            

        fprintf (stderr, "\n\nPlease read the terms and conditions of the license:\n\n"); 
        license_display(); 
        fprintf (stderr, "\n"); 

        license_exit_if_not_confirmed ("Do you accept the terms and conditions of the license?", QDEF_NONE);

        os           = arch_get_os();          
        dist         = arch_get_distribution();
        cores        = arch_get_num_cores();
        user_host    = arch_get_user_host();
        is_script    = false;
        memcpy (rec.ip, arch_get_ip_addr ("Failed to register the license"), ARCH_IP_LEN);  
    }

    rec.timestamp = str_time();    
    rec.machine_time = time (0);   

    static Buffer license_data = EMPTY_BUFFER;
    license_generate (&license_data);
    printf (".\n"); 

    if (!n_fields) {
        // fprintf (stderr, "\nThank you. To complete your license registration, genozip will now submit the following information to the genozip licensing server:\n\n");

        // note: text needs to match scripts/register.sh
        fprintf (stderr, "=====================================================================\n");
        fprintf (stderr, LIC_FIELD_TYPE       ": %s\n", lic_types[rec.lic_type]);
        fprintf (stderr, LIC_FIELD_INSTITUTION": %s\n", rec.institution);
        fprintf (stderr, LIC_FIELD_NAME       ": %s\n", rec.name);
        fprintf (stderr, LIC_FIELD_EMAIL      ": %s\n", rec.email);
        fprintf (stderr, "System info: OS=%s cores=%u IP=%s\n", os, cores, rec.ip);
        fprintf (stderr, "Username: %s\n", user_host);
        fprintf (stderr, "Genozip info: version=%s distribution=%s\n", GENOZIP_CODE_VERSION, dist);
        fprintf (stderr, "Genozip license number: %u\n", rec.license_num);
        if (rec.lic_type == LIC_TYPE_STANDARD || rec.lic_type == LIC_TYPE_ENTERPRISE) fprintf (stderr, "Send statistics: %s\n", rec.allow_stats);
        fprintf (stderr, "I accept the terms and conditions of the Genozip license\n");
        fprintf (stderr, "=====================================================================\n\n");
        
        // license_exit_if_not_confirmed ("Proceed with completing the registration?", QDEF_YES);
    }
        
    bool submitted = license_submit (os, cores, is_script, user_host, dist);

    ASSINP0 (submitted,
             "Failed to register the license, possibly because the Internet is not accessible or the registration server "
             "(which is hosted on a Google server) is not accessible. If this problem persists, you can register manually by "
             "sending an email to register@genozip.com - copy & paste the lines between the \"======\" into the email message.\n");

    if (!is_script) {
        char query[sizeof(rec.email)+256];
        char code[7] = "";
        sprintf (query, "\nA 6-digit verification code was emailed to %s.\n\n"
                        "(If you did not receive it within 2 minutes, please contact register@genozip.com)\n\n"
                        "Please enter it: ", rec.email);
        str_query_user (query, code, sizeof(code), false, license_verify_code, NULL);
    }
    
    ASSINP (file_put_data (filename, STRb(license_data), S_IRUSR | S_IRGRP), 
            "Failed to write license file %s: %s. If this is unexpected, email "EMAIL_SUPPORT" for help.", filename, strerror (errno));

    if (!n_fields) {
        fprintf (stderr, "\nA Genozip %s License has been granted.\n\n"
                         "Getting started: " WEBSITE_QUICK_GUIDE "\n\n"
                         "Support: " EMAIL_SUPPORT "\n\n", lic_types[rec.lic_type]);

        if (rec.lic_type == LIC_TYPE_ACADEMIC)
            fprintf (stderr, "Please take a moment now to make a note to not forget to cite Genozip:\n"
                             "Lan, D., et al. (2022) Genozip 14 - advances in compression of BAM and CRAM files (pre-print), doi: https://doi.org/10.1101/2022.09.12.507582\n"
                             "Lan, D., et al. (2021) Genozip: a universal extensible genomic data compressor, Bioinformatics, 37, 2225-2230\n"
                             "Lan, D., et al. (2020) genozip: a fast and efficient compression tool for VCF files, Bioinformatics, 36, 4091-4092\n\n");

        else if (rec.lic_type == LIC_TYPE_EVAL)
            fprintf (stderr, "We will contact you in 30 days to ask whether you are interested to proceed with purchasing a Standard License.\n\n");
    }

    buf_destroy (license_data);
}

// IF YOU'RE CONSIDERING TAMPERING WITH THIS CODE TO BYPASS THE REGISTRTION OR STATS SUBMISSION, DON'T! It would be a violation of the license,
// and might put you personally as well as your organization at legal and financial risk - see "Severly unauthorized use of Genozip"
// section of the license. Rather, please contact sales@genozip.com to discuss which license would be appropriate for your case.
uint32_t license_get_number (void)
{
    return rec.license_num;
}

bool license_allow_stats (void)
{
    switch (rec.lic_type) {
        case LIC_TYPE_STANDARD:
        case LIC_TYPE_STANDARD_LEGACY:
        case LIC_TYPE_ENTERPRISE: 
            return !strcmp (rec.allow_stats, "Yes");

        case LIC_TYPE_FIRST: 
            return false; // no consent yet

        case LIC_TYPE_ACADEMIC:
        case LIC_TYPE_EVAL:
            return true;
        
        default:
            ABORT_R ("Invalid license_type=%d", rec.lic_type);
    }
}

bool license_allow_tip (void)
{
    return rec.lic_type != LIC_TYPE_FIRST;
}

rom license_get_one_line (void)
{
    static char s[sizeof (rec) + sizeof (rec.name) + 200];

    sprintf (s, "License v%s type: %s granted to: %s for use by: %s accepted by: %s <%s> on %s from IP=%s", 
             rec.version, lic_types[rec.lic_type], rec.institution, rec.name, rec.name, rec.email, rec.timestamp.s, rec.ip);

    return s;
}

void license_display (void)
{
    rom filename = license_get_filename (false);
    static Buffer license_data = {};
    
    if (file_exists (filename) && !flag.force) 
        file_get_file (evb, filename, &license_data, "license_data", 0, false/*bc we allow UTF-8*/, true);

    bool html = (flag.lic_width < 0);
    int width = ABS(flag.lic_width);
    rom newline = html ? "<br><br>\n" : "\n\n";
    rom html_header = html ? "<!DOCTYPE html> <html lang=\"en\"> <body style=\"color:white;\">" : NULL;

    // case: user has already accepted the license and it is new style license - display the license file
    if (license_data.len > 100) {
        str_split (license_data.data, license_data.len, 0, '\n', line, false);
        str_nul_separate (line);
        str_print_text (lines, n_lines-1, "", newline, html_header, width);
    }
    
    // case: license not yet accepted or old style (up to 12.0.6) license - display the current version license
    else
        str_print_text (license, ARRAY_LEN(license), "", newline, html_header, width);  // Makefile sets lic_width to a fixed width for Windows Installer and for Docs
}

StrNotice license_print_default_notice (void)
{
    StrNotice notice = {}; 

    if (IS_PIZ) license_load(); // happens when PIZ is run in --test after ZIP

    if (rec.lic_type == LIC_TYPE_ACADEMIC)
        sprintf (notice.s, "Genozip is licensed to %s for use by %s, for academic research purposes only (see "WEBSITE_PRICING_FAQ"). Other use is prohibited. To get a non-academic license, email " EMAIL_SALES ".",
                 rec.institution, rec.name);

    else if (rec.lic_type == LIC_TYPE_EVAL) {
        int days_left = 30 - (int)(time(0)-rec.machine_time)/24/60/60;

        sprintf (notice.s, "Genozip is licensed to %s for use by %s for a 30-day evaluation. %u day%s left. To get a license, email " EMAIL_SALES ".",
                 rec.institution, rec.name, days_left, days_left!=1 ? "s" : "");
    }
    
    return notice;
}

void license_show_deep_notice (void)
{
    if (flag.deep && (rec.lic_type == LIC_TYPE_STANDARD || rec.lic_type == LIC_TYPE_ACADEMIC))
        iprintf ("Note: using --deep requires a Genozip Enterprise License. It is provided to you on an evaluation basis, limited to compressing %u files\n"
                 "To upgrade to a Enterprise License: " WEBSITE_BUY " or contact " EMAIL_SALES "\n", EVAL_NUM_FILES);
}

void license_one_file_compressed (DataType dt)
{    
    if (rec.lic_type  == LIC_TYPE_EVAL || 
        (rec.lic_type == LIC_TYPE_STANDARD_LEGACY && flag.deep) || 
        (rec.lic_type == LIC_TYPE_STANDARD && (flag.deep || flag.pair)))
        counter_increment (license_get_filename (false), 1);
}

