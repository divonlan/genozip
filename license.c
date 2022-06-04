// ------------------------------------------------------------------
//   license.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifdef _WIN32
#include <direct.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "genozip.h"
#include "strings.h"
#include "arch.h"
#include "profiler.h" // for TimeSpecType
#include "md5.h"
#include "url.h"
#include "license.h"
#include "version.h"
#include "buffer.h"
#include "flags.h"
#include "md5.h"
#include "website.h"
#include "version.h"
#include "file.h"

// these field names appear in the license file starting V12.0.7
#define LIC_FIELD_TYPE         "License type"
#define LIC_FIELD_VERSION      "Genozip license version"
#define LIC_FIELD_INSTITUTION  "License granted to"
#define LIC_FIELD_NAME         "Accepted by (name)"
#define LIC_FIELD_EMAIL        "Accepted by (email)"
#define LIC_FIELD_MACHINE_TIME "Machine time"
#define LIC_FIELD_TIMESTAMP    "Timestamp of acceptance"
#define LIC_FIELD_IP           "IP address of acceptance"
#define LIC_FIELD_NUMBER       "License number"

#include "text_license.h"

static rom license_filename = NULL;  // non-standard filename set with --licfile

static struct {
    bool initialized;
    LicenseType lic_type; 
    char name[256], email[256], institution[1024], ip[ARCH_IP_LEN], version[20];
    StrText timestamp;
    int64_t machine_time; // timestamp expressed as seconds since epoch
    uint32_t license_num;
} rec = {};

static rom lic_types[4] = { "", "Academic", "30-day evaluation", "Paid" }; // these strings are referred to in register.sh

static uint32_t license_calc_number (ConstBufferP license_data)
{
    char data_no_ws[license_data->len];
    unsigned data_no_ws_len = str_remove_whitespace (license_data->data, license_data->len, data_no_ws);        

    return md5_do (data_no_ws, data_no_ws_len).words[0];
}

static void license_generate (BufferP license_data)
{
    for (int i=0; i < ARRAY_LEN(license); i++) {
        buf_add_string (evb, license_data, license[i]); // allocs one extra char
        BNXTc (*license_data) = '\n';
    }

    bufprintf (evb, license_data,     // note: the license data includes the Genozip version
               LIC_FIELD_TYPE": %d\n" // added v14
               LIC_FIELD_INSTITUTION": %s\n"
               LIC_FIELD_NAME": %s\n"
               LIC_FIELD_EMAIL": %s\n"
               LIC_FIELD_MACHINE_TIME": %"PRIu64"\n"
               LIC_FIELD_TIMESTAMP": %s\n"
               LIC_FIELD_IP": %s\n",
               rec.lic_type, rec.institution, rec.name, rec.email, rec.machine_time, rec.timestamp.s, rec.ip);
    
    rec.initialized = true;
    rec.license_num = license_calc_number (license_data);
    strcpy (rec.version, GENOZIP_CODE_VERSION);

    bufprintf (evb, license_data, LIC_FIELD_NUMBER": %u\n", rec.license_num);
}

void license_set_filename (rom filename)
{
    struct stat sb;
    ASSINP (!stat (filename, &sb), "Failed to access license file %s: %s", filename, strerror (errno));

    license_filename = filename;
}

static rom get_license_filename (bool create_folder_if_needed)
{
    if (license_filename) return license_filename; // non-standard filename set with --licfile

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

    char *filename = MALLOC (strlen(folder) + 50);
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

// IF YOU'RE CONSIDERING TAMPERING WITH THIS CODE TO BYPASS THE REGISTRTION, DON'T! It would be a violation of the license,
// and might put you personally as well as your organization at legal and financial risk - see "Severly unauthorized use of Genozip"
// section of the license. Rather, please contact sales@genozip.com to discuss which license would be appropriate for your case.
static void license_load (void)
{
    if (rec.initialized) return;

    rom filename = get_license_filename (true);
    
    if (!file_exists (filename)) {
        flag.do_register = "";
        license_register ();
        return;
    }

    file_split_lines (filename, "license");
    
    char license_num_str[30] = "", lic_type_str[16]="", machine_time_str[24]="";
    #define COPY_FIELD(var,field) strncpy (var, license_load_field (field, n_lines, lines, line_lens), sizeof (var)-1)

    COPY_FIELD (lic_type_str,    LIC_FIELD_TYPE);  // added v14

    // licenses prior to v14 don't have this field, requiring re-registration
    if (!str_get_int_range8 (lic_type_str, strlen (lic_type_str), 1, NUM_LIC_TYPES-1, &rec.lic_type)) goto reregister;

    COPY_FIELD (rec.version,      LIC_FIELD_VERSION);
    COPY_FIELD (rec.institution,  LIC_FIELD_INSTITUTION);
    COPY_FIELD (rec.name,         LIC_FIELD_NAME);
    COPY_FIELD (rec.email,        LIC_FIELD_EMAIL);
    COPY_FIELD (machine_time_str, LIC_FIELD_MACHINE_TIME);
    COPY_FIELD (rec.timestamp.s,  LIC_FIELD_TIMESTAMP);
    COPY_FIELD (rec.ip,           LIC_FIELD_IP);
    COPY_FIELD (license_num_str,  LIC_FIELD_NUMBER);

    if (!str_get_uint32 (license_num_str, strlen (license_num_str), &rec.license_num)) goto reregister;
    if (!str_get_int (machine_time_str, strlen (machine_time_str), &rec.machine_time)) goto reregister;

    data.len -= line_lens[n_lines-1] + 2;
    if (rec.license_num != license_calc_number (&data)) goto reregister;

    ASSINP0 (rec.lic_type != LIC_TYPE_EVAL || time(0)-rec.machine_time < (30*24*60*60),
             "Your 30 evaluation period is over. Please contact sales@genozip.com to purchase a license or to request an extension of the evaluation period");

    rec.initialized = true;

    buf_destroy (data);
    return;

reregister:
    file_remove (filename, true);

    // if stdin or stderr is redirected - we cannot start an interactive registration flow
    ASSINP0 (isatty(0) && isatty(2), "Genozip license terms & conditions have changed, please re-register by running: genozip --register");

    fprintf (stderr, "Genozip license terms & conditions have changed, please re-register:\n\n");
    flag.do_register = "";
    license_register();
}

static bool license_submit (char update, rom os, unsigned cores, rom endianity, rom user_host, rom dist)
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
                       "?entry.344252538=%.100s"
                       "&entry.926671216=%.50s"
                       "&entry.1734045469=%.50s"
                       "&entry.2009586582=%.20s"
                       "&entry.119966790=%c"
                       "&entry.81542373=%s"
                       "&entry.1668073218=%u"
                       "&entry.1943454647=%s"
                       "&entry.1763961212=%s"
                       "&entry.1655649315=%u"
                       "&entry.186159495=%s"
                       "&entry.1598028195=%s"
                       "&entry.1384715202=%s";

    char *institutionE = url_esc_non_valid_chars (rec.institution);
    char *nameE        = url_esc_non_valid_chars (rec.name);
    char *emailE       = url_esc_non_valid_chars (rec.email);
    char *lic_typeE    = url_esc_non_valid_chars (lic_types[rec.lic_type]);
    char *osE          = url_esc_non_valid_chars (os);
    char *user_hostE   = url_esc_non_valid_chars (user_host);

    char url[sizeof (rec) + 200];
    sprintf (url, url_format, institutionE, nameE, emailE, lic_typeE, update, osE, cores, rec.ip, user_hostE, rec.license_num, rec.version, dist, endianity);

    bool success = url_read_string (url, NULL, 0) >= 0;
    
    FREE (institutionE); FREE (nameE); FREE (emailE); FREE (osE); FREE (user_hostE);
    return success;
}

static bool license_verify_email (char *response, unsigned response_size, rom unused)
{
    // sanity check that this is an email address
    return strlen (response) > 3 && strchr (response, '@') && strchr (response, '.');
}

static bool license_verify_name (char *response, unsigned response_size, rom unused)
{
    if (!strchr (response, ' ')) {
        fprintf (stderr, "Please enter your full name\n");
        return false;
    }
    
    return true;
}

static bool license_verify_license (char *response, unsigned response_size, rom unused)
{
    return strlen (response) == 1 && (*response=='1' || *response=='2' || *response=='3');
}

static void license_exit_if_not_confirmed (rom response)
{
    if (response[0] == 'N') {
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
    char confirm[100], lic_type[100], update[100];
    rom os, dist, endianity, user_host;
    unsigned cores;

    str_split (flag.do_register, strlen (flag.do_register), 11, '|', field, true);
    str_nul_separate (field);

    // if stdin or stderr is redirected - we cannot ask the user an interactive question
    ASSINP0 (isatty(0) && isatty(2), "Use of Genozip is free for academic purposes, but requires registration. Please run: genozip --register.\n"
                                     "If you are unable to register (for example because this is a batch-job machine) please see: " WEBSITE_USING_ON_HPC);

    rom filename = get_license_filename (true);

    if (!n_fields) {

        fprintf (stderr, "Welcome to Genozip!\n\n"
                         "- Genozip is FREE for for academic research and some other purposes defined in the license (see "WEBSITE_LICENSE"), but requires registration.\n"
                         "- For clinical, mixed research/clinical, commercial and other cases (examples: "WEBSITE_COMMERCIAL") you may evaluate Genozip for free for 30 days.\n\n");

        if (file_exists (filename)) 
            str_query_user ("You are already registered. Are you sure you want to re-register again? (y or n) ", confirm, sizeof(confirm), str_verify_y_n, NULL);
        else 
            str_query_user ("Would you like to register now? ([y] or n) ", confirm, sizeof(confirm), str_verify_y_n, "Y");

        license_exit_if_not_confirmed (confirm);
    }

    file_remove (filename, true); // remove old license, if one exists
    
    if (n_fields) { // fields correspond to register.sh
        rec.lic_type = atoi (fields[0]);
        strncpy (rec.institution, fields[1], sizeof(rec.institution)-1);
        strncpy (rec.name,        fields[2], sizeof(rec.name)-1);
        strncpy (rec.email,       fields[3], sizeof(rec.email)-1);
        strncpy (update,          fields[4], sizeof(update)-1);
        strncpy (rec.ip,          fields[5], sizeof(rec.ip)-1);
        os           = fields[6];
        dist         = fields[7]; 
        endianity    = fields[8];
        user_host    = fields[9];
        cores        = atoi(fields[10]);
    }
    else {
        fprintf (stderr, "\nLicense details -\n");
    
        str_query_user ("\nInstitution / Company name: ", rec.institution, sizeof(rec.institution), str_verify_not_empty, NULL);

        str_query_user ("\nYour name: ", rec.name, sizeof(rec.name), license_verify_name, NULL);
        
        str_query_user ("\nYour email address: ", rec.email, sizeof(rec.email), license_verify_email, NULL);
        
        str_query_user ("\nWhat type of license do you require?\n\n"
                        "1. Academic license (free, available to recognized research institutions, but excluding mixed clinical/research use)\n\n"
                        "2. Non-academic license - 30 days evaluation (free for 30 days, for clinical or mixed clinical/research settings, biotech/agrotech/SaaS companies and other non-academic uses)\n\n"
                        "3. I have already paid for a non-academic license\n\n"
                        "Remember your Mom taught you to be honest!\n\n"
                        "Please enter 1, 2 or 3: ",
                        lic_type, sizeof(lic_type), license_verify_license, NULL);
    
        rec.lic_type = lic_type[0] - '0';
    
        str_query_user ("\nShall we update you by email when new features are added to genozip? ([y] or n) ", 
                        update, sizeof(update), str_verify_y_n, "Y");

        fprintf (stderr, "\n\nPlease read the terms and conditions of the license:\n\n"); 
        license_display(); 
        fprintf (stderr, "\n"); 

        str_query_user ("Do you accept the terms and conditions of the license? (y or n) ", confirm, sizeof(confirm), str_verify_y_n, NULL);
        license_exit_if_not_confirmed (confirm);

        os           = arch_get_os();
        dist         = arch_get_distribution();
        cores        = arch_get_num_cores();
        endianity    = arch_get_endianity();
        user_host    = arch_get_user_host();
        memcpy (rec.ip, arch_get_ip_addr ("Failed to register the license"), ARCH_IP_LEN);
    }

    rec.timestamp = str_time();
    rec.machine_time = time (0);

    static Buffer license_data = EMPTY_BUFFER;
    license_generate (&license_data);

    if (!n_fields) {
        fprintf (stderr, "\nThank you. To complete your license registration, genozip will now submit the following information to the genozip licensing server:\n\n");

        // note: text needs to match scripts/register.sh
        fprintf (stderr, "=====================================================================\n");
        fprintf (stderr, LIC_FIELD_TYPE       ": %s\n", lic_types[rec.lic_type]);
        fprintf (stderr, LIC_FIELD_INSTITUTION": %s\n", rec.institution);
        fprintf (stderr, LIC_FIELD_NAME       ": %s\n", rec.name);
        fprintf (stderr, LIC_FIELD_EMAIL      ": %s\n", rec.email);
        fprintf (stderr, "Send new feature updates: %s\n", update[0]=='Y' ? "Yes" : "No");
        fprintf (stderr, "System info: OS=%s cores=%u endianity=%s IP=%s\n", os, cores, endianity, rec.ip);
        fprintf (stderr, "Username: %s\n", user_host);
        fprintf (stderr, "Genozip info: version=%s distribution=%s\n", GENOZIP_CODE_VERSION, dist);
        fprintf (stderr, "Genozip license number: %u\n", rec.license_num);
        fprintf (stderr, "I accept the terms and conditions of the Genozip license\n");
        fprintf (stderr, "=====================================================================\n\n");
        
        str_query_user ("Proceed with completing the registration? ([y] or n) ", confirm, sizeof(confirm), str_verify_y_n, "Y");
        license_exit_if_not_confirmed (confirm);
    }
        
    bool submitted = license_submit (update[0], os, cores, endianity, user_host, dist);

    ASSINP0 (submitted,
             "Failed to register the license, possibly because the Internet is not accessible or the registration server "
             "(which is hosted on a Google server) is not accessible. If this problem persists, you can register manually by "
             "sending an email to register@genozip.com - copy & paste the lines between the \"======\" into the email message.\n");

    ASSINP (file_put_data (filename, STRb(license_data), S_IRUSR), 
            "Failed to write license file %s: %s. If this is unexpected, email "EMAIL_SUPPORT" for help.", filename, strerror (errno));

    if (!n_fields) {
        fprintf (stderr, "\nSUCCESS. A Genozip license has been granted:\n"
                         "License type: %s\nLicensee: %s\nFor use by %s\n\n" 
                         "Documentation: " GENOZIP_URL "\n\n"
                         "Support: " EMAIL_SUPPORT "\n\n", lic_types[rec.lic_type], rec.institution, rec.name);

        if (lic_type[0]=='1')
            fprintf (stderr, "Please take a moment now to make a note to not forget to cite Genozip:\n"
                             "Lan, D., et al. (2021) Genozip: a universal extensible genomic data compressor, Bioinformatics, 37, 2225-2230\n"
                             "Lan, D., et al. (2020) genozip: a fast and efficient compression tool for VCF files, Bioinformatics, 36, 4091-4092\n\n");

        else if (lic_type[0]=='2')
            fprintf (stderr, "We will contact you in 30 days to ask whether you are interested to proceed with purchasing a license.\n\n");
    }

    buf_destroy (license_data);
}

// IF YOU'RE CONSIDERING TAMPERING WITH THIS CODE TO BYPASS THE REGISTRTION, DON'T! It would be a violation of the license,
// and might put you personally as well as your organization at legal and financial risk - see "Severly unauthorized use of Genozip"
// section of the license. Rather, please contact sales@genozip.com to discuss which license would be appropriate for your case.
uint32_t license_get_number (void)
{
    license_load();
    return rec.license_num;
}

LicenseType license_get_type (void)
{
    license_load();
    return rec.lic_type;
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
    rom filename = get_license_filename (false);
    static Buffer license_data = {};
    
    if (file_exists (filename) && !flag.force) 
        file_get_file (evb, filename, &license_data, "license_data", true);

    // case: user has already accepted the license and it is new style license - display the license file
    if (license_data.len > 100) {
        str_split (license_data.data, license_data.len, 0, '\n', line, false);
        str_nul_separate (line);
        str_print_text (lines, n_lines-1, "", "\n\n", flag.lic_width);
    }
    
    // case: license not yet accepted or old style (up to 12.0.6) license - display the current version license
    else
        str_print_text (license, sizeof(license) / sizeof(char*), "", "\n\n", flag.lic_width);  // Makefile sets lic_width to a fixed width for Windows Installer and for Docs
}

