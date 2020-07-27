// ------------------------------------------------------------------
//   license.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifdef _WIN32
#include <direct.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "genozip.h"
#include "strings.h"
#include "arch.h"
#include "text_license.h"
#include "profiler.h" // for TimeSpecType
#include "md5.h"
#include "url.h"
#include "license.h"

void license_display (void)
{
    str_print_text (license, sizeof(license) / sizeof(char*), "", "\n\n", 
                    flag_force ? 60 : 0); // --license --force means output license in Windows installer format (used by Makefile) - 60 is width of InstallForge license text field
}

static uint32_t license_generate_num(void)
{
    TimeSpecType timer; 
    clock_gettime(CLOCK_REALTIME, &timer); 

    static Md5Hash md5;
    md5 = md5_do (&timer, sizeof (timer));
    
    if (!md5.words[0]) return license_generate_num(); // chance of 1 in 4 billion that we will need to try again

    return md5.words[0]; // 32 bit
}

static char *get_license_filename (bool create_folder_if_needed)
{
#ifdef _WIN32
    ASSERT0 (getenv ("APPDATA"), "%s: cannot store license, because APPDATA env var is not defined");

    char folder[500];
    sprintf (folder, "%s/genozip", getenv ("APPDATA"));

    if (create_folder_if_needed) {
        int ret = _mkdir (folder); 
        ASSERT (ret >= 0 || errno == EEXIST, "Error: failed to created the folder %s", folder);
    }

#else
    const char *folder = getenv ("HOME");
    ASSERT0 (folder, "%s: cannot store license, because HOME env var is not defined");
#endif    

    char *filename = malloc (strlen(folder) + 50);
    sprintf (filename, "%s/.genozip_license", folder);

    return filename;
}

static void license_store_locally (uint32_t license_num)
{
    char *filename = get_license_filename (true);

    FILE *fp = fopen (filename, "wb");
    ASSERT (fp, "Error: failed to open %s for writing: %s", filename, strerror (errno));

    int res = fprintf (fp, "%u\n", license_num);
    ASSERT (res > 0, "Error: failed to write to %s: %s", filename, strerror (errno));

    fclose (fp);

    free (filename);
}

static uint32_t licence_retrieve_locally (void)
{
    char *filename = get_license_filename (true);

    FILE *fp = fopen (filename, "rb");
    if (!fp) return 0; // no license

    uint32_t license_num=0; // initialize - so that we return 0 if fscanf fails to read a number
    fscanf (fp, "%u", &license_num);
    
    fclose (fp);
    free (filename);

    return license_num;
}

static void license_submit (const char *institution, const char *name, const char *email, 
                            char commerical, char update, 
                            const char *os, unsigned cores, const char *ip, uint32_t license_num)
{
    // reference: https://stackoverflow.com/questions/18073971/http-post-to-a-google-form/47444396#47444396

    #define PREFIX "https://docs.google.com/forms/d/e/1FAIpQLSc6pSBIOBsS5Pu-JNvfnLWV2Z1W7k-4f2pKRo5rTbiU_sCnIw/formResponse"
    char *url_format = PREFIX
                       "?entry.344252538=%.100s"
                       "&entry.926671216=%.50s"
                       "&entry.1734045469=%.50s"
                       "&entry.2009586582=%c"
                       "&entry.119966790=%c"
                       "&entry.81542373=%s"
                       "&entry.1668073218=%u"
                       "&entry.1943454647=%s"
                       "&entry.1655649315=%u";

    char *institutionE = url_esc_non_valid_chars (institution);
    char *nameE        = url_esc_non_valid_chars (name);
    char *emailE       = url_esc_non_valid_chars (email);
    char *osE          = url_esc_non_valid_chars (os);

    char url[600];
    sprintf (url, url_format, institutionE, nameE, emailE, commerical, update, osE, cores, ip, license_num);

    url_read_string (url, NULL, 0);

    free (institutionE); free (nameE); free (emailE); free (osE);
}

static bool license_verify_email (char *response, unsigned response_size, const char *unused)
{
    // sanity check that this is an email address
    return strlen (response) > 3 && strchr (response, '@') && strchr (response, '.');
}

static bool license_verify_name (char *response, unsigned response_size, const char *unused)
{
    if (!strchr (response, ' ')) {
        fprintf (stderr, "Please enter your full name\n");
        return false;
    }
    
    return true;
}

static void license_exit_if_not_confirmed (const char *response)
{
    if (response[0] == 'N') {
        fprintf (stderr, "\nYou have not registered. You may register at any time in the future.\nWishing you a wonderful day from the genozip team!\n");
        exit_ok;
    }
}

// load license if it exists, or register a new one

// IF YOU'RE CONSIDERING EDITING THIS CODE TO BYPASS THE REGISTRTION, DON'T! It would be a violation of the license,
// and might put you personally as well as your organization at legal and financial risk.
// Rather, please contact sales@genozip.com to discuss which license would be appropriate for your case.
uint32_t license_get (void)
{
    uint32_t license_num = licence_retrieve_locally();
    if (license_num && !flag_register) return license_num; // license exists (and user didn't request --register) - we're done!

    // UI flow to generate a new license for the user

    // if stdin or stderr is redirected - we cannot ask the user an interactive question
    ASSERT (isatty(0) && isatty(2), "%s: Use of genozip is free, but requires registration. Please run: genozip --register", global_cmd);

    fprintf (stderr, "Welcome to genozip!\n\n"
                     "The use of genozip for non-commercial purposes (as defined in the license) is FREE, but requires registration\n\n");

    char confirm[1000], institution[1000], name[1000], email[1000], commercial[1000], update[1000];

    str_query_user ("Would you like to register now? ([y] or n)", confirm, sizeof(confirm), str_verify_y_n, "Y");
    license_exit_if_not_confirmed (confirm);

    fprintf (stderr, "\nDetails of the person to be granted the license -\n");
    str_query_user ("\nInstitution / Company name: ", institution, sizeof(institution), str_verify_not_empty, NULL);
    
    str_query_user ("\nYour name: ", name, sizeof(name), license_verify_name, NULL);
    
    str_query_user ("\nYour email address: ", email, sizeof(email), license_verify_email, NULL);
    
    str_query_user ("\nDo you require a commercial license? If yes, we will contact you (this will not stop the registration now) (y or [n]) ", 
                    commercial, sizeof(commercial), str_verify_y_n, "N");

    str_query_user ("\nShall we update you by email when new features are added to genozip? ([y] or n) ", 
                    update, sizeof(update), str_verify_y_n, "Y");

    fprintf (stderr, "\n\nPlease read the terms and conditions of the non-commercial license:\n\n"); 
    license_display(); 
    fprintf (stderr, "\n"); 

    str_query_user ("Do you accept the terms & conditions of the license? (y or n) ", confirm, sizeof(confirm), str_verify_y_n, NULL);
    license_exit_if_not_confirmed (confirm);

    const char *os = arch_get_os();
    const char *ip = arch_get_ip_addr();
    unsigned cores = arch_get_num_cores();
    license_num    = license_generate_num();

    fprintf (stderr, "\nThank you. To complete your license registration, genozip will now submit the following information to the genozip licensing server:\n\n");
    fprintf (stderr, "Licensee institution / company name: %s\n", institution);
    fprintf (stderr, "Licensee name: %s\n", name);
    fprintf (stderr, "Licensee email address: %s\n", email);
    fprintf (stderr, "Commercial: %s\n", commercial[0]=='Y' ? "Yes" : "No");
    fprintf (stderr, "Send new feature updates: %s\n", update[0]=='Y' ? "Yes" : "No");
    fprintf (stderr, "System info: OS=%s cores=%u IP=%s\n", os, cores, ip);
    fprintf (stderr, "genozip license number: %u\n\n", license_num);

    str_query_user ("Proceed with completing the registration? ([y] or n) ", confirm, sizeof(confirm), str_verify_y_n, "Y");
    license_exit_if_not_confirmed (confirm);
    
    license_submit (institution, name, email, commercial[0], update[0], os, cores, ip, license_num);

    license_store_locally (license_num);

    fprintf (stderr, "\nCongratulations! Your genozip non-commerical license has been granted.\n\n");

    return license_num;
}
