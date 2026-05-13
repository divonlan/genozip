// ------------------------------------------------------------------
//   tip.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "flags.h"
#include "arch.h"
#include "license.h"
#include "sam.h"
#include "file.h"
#include "tip.h"

static bool dt_encountered[NUM_DATATYPES] = {};

void tip_dt_encountered (DataType dt)
{
    dt_encountered[dt] = true;
}

// called in genozip, or, if --test and genounzip --test does not return, after the testing is completed 
void tip_print_genozip (void)
{
    #define E(dt) dt_encountered[DT_##dt]

    if (flag.no_tip || // protect from spamming the user with more than one tip
        !license_allow_tip()) return; 

    StrText1K acadmic_notice = license_academic_tip(); // empty unless acadmic license

    if (!is_info_stream_terminal) {
        if (acadmic_notice.s[0]) iprintf ("\n%s\n", acadmic_notice.s);
        return;
    }

    rom valid_tips[256];
    int n=0;

    if (acadmic_notice.s[0])
        for (int i=0; i < 5; i++)
            valid_tips[n++] = acadmic_notice.s; // 5X more likely than other tips

    // italics from: https://en.wikipedia.org/wiki/Mathematical_Alphanumeric_Symbols
    if (!flag.default_make_ref) {
        valid_tips[n++] = _TIP "Want to understand how Genozip works? See the paper: " PAPER2 "\n";
        valid_tips[n++] = _TIP "Use 'genocat --downsample=𝑟𝑎𝑡𝑒' to downsample your data. " WEBSITE_DOWNSAMPLING "\n";
        valid_tips[n++] = _TIP "Use 'genozip --is-exactable' to know if 'genounzip --bgzf=exact' will be able to re-compress the file with identical gz compression. " WEBSITE_GZ "\n";
        valid_tips[n++] = _TIP "Use 'genozip --tar mydata.tar 𝑓𝑖𝑙𝑒𝑠 𝑎𝑛𝑑 𝑑𝑖𝑟𝑒𝑐𝑡𝑜𝑟𝑖𝑒𝑠' to archive entire directories. " WEBSITE_ARCHIVING "\n";
        valid_tips[n++] = _TIP "genozip files are an excellent way to share and publish data – uncompressing genozip files is always free\n";
        valid_tips[n++] = _TIP "genozip can compress a file directly from a URL: 'genozip ftp://mysite/sample.fastq.gz'. " WEBSITE_COMPRESS_URL "\n";
        valid_tips[n++] = _TIP "Use --password=𝑝𝑎𝑠𝑠𝑤𝑜𝑟𝑑 to also encrypt the data. " WEBSITE_ENCRYPTION "\n";
    }

    if (!(E(SAM) || E(BAM) || E(VCF) || E(BCF)) && (license_is_academic() || license_is_eval()))
        valid_tips[n++] = "Please take a moment now to make a note to not forget to cite Genozip:\n" PAPER2_CITATION "\n";

    if (dist_is_github())  
        valid_tips[n++] = _TIP "Do you like Genozip? Please support it by starring it on github: " GITHUB_REPO "\n";

    if (E(SAM) || E(BAM) || E(FASTQ)) 
        valid_tips[n++] = _TIP "Use 'genocat --coverage myfile.genozip' get coverage information for BAM or FASTQ. " WEBSITE_COVERAGE "\n";

    if (E(BCF))
        valid_tips[n++] = _TIP "genozip compresses VCF files 5-10X faster than it compresses BCF files";
    
    if ((E(SAM) || E(BAM)) && (license_is_academic() || license_is_eval()))
        valid_tips[n++] = _TIP "Please take a moment now to make a note to not forget to cite Genozip:\n" PAPER3_CITATION "\n";

    bool is_deep = flag.deep || ((IS_PIZ && z_file && (Z_DT(BAM) || Z_DT(SAM))) && z_file->z_flags.is_deep);
    
    if (is_deep && (license_is_academic() || license_is_eval()))
        valid_tips[n++] = _TIP "Please take a moment now to make a note to not forget to cite Genozip:\n" PAPER4_CITATION "\n";

    if (!is_deep && (E(SAM) || E(BAM) || E(FASTQ)) && sam_get_deep_tip()) 
        for (int i=0; i < 5; i++) // 5X more likely
            valid_tips[n++] = sam_get_deep_tip();

    if (E(FASTQ) && !flag.bam_assist && !flag.deep)
        for (int i=0; i < 5; i++) // 5X more likely
            valid_tips[n++] = _TIP "Use '--bamass "_BAMFILE"' to allow genozip to consult the BAM alignments corresponding to this FASTQ file for faster and better compression. " WEBSITE_BAMASS "\n";

    if ((E(VCF) || E(BCF)) && (license_is_academic() || license_is_eval()))
        valid_tips[n++] = _TIP "Please take a moment now to make a note to not forget to cite Genozip:\n" PAPER1_CITATION "\n";

    if (!flag.optimize && !flag.make_reference) {
        if (E(SAM) || E(BAM) || E(CRAM)) valid_tips[n++] = _TIP "--optimize allows genozip to make minor modifications to the data that usually have no impact on downstream analysis, for better compression. " WEBSITE_OPT_BAM "\n";
        else if (E(VCF) || E(BCF))       valid_tips[n++] = _TIP "--optimize allows genozip to make minor modifications to the data that usually have no impact on downstream analysis, for better compression. " WEBSITE_OPT_VCF "\n";
        else if (E(FASTQ))               valid_tips[n++] = _TIP "--optimize allows genozip to make minor modifications to the data that usually have no impact on downstream analysis, for better compression. " WEBSITE_OPT_FASTQ "\n";
    }

    if (!flag.best && !flag.fast && !flag.low_memory && !flag.make_reference) 
        valid_tips[n++] = _TIP "Use --best for the best compression (slower)\n";

    if (flag.default_make_ref)
        valid_tips[n++] = _TIP "Use --make-reference=[minimal|tiny|small|medium|large] (default: medium) to control the size of the reference file. Larger reference files result is faster and better compression, but using them consumes more RAM. " WEBSITE_REFERENCE "\n";

    if (license_is_eval() && !flag.show_stats && !flag.make_reference)
        valid_tips[n++] = _TIP "Use --stats to see detailed compression statistics\n";

    if (license_is_eval() || license_is_standard() || license_is_enterprise())
        valid_tips[n++] = _TIP "\"genozip --sendto\" lets your clients compress files using 𝑦𝑜𝑢𝑟 Genozip license ahead of sending them to you. " WEBSITE_PREMIUM "\n";

    if ((arch_get_max_resident_set() > 100 GB || flag.is_windows || flag.is_wsl || flag.is_mac) &&
        !flag.make_reference)
        valid_tips[n++] = _TIP "Use --low-memory to reduce RAM consumption at the expense of compression size and time\n";

    iprintf ("\n%s\n", valid_tips[clock() % n]); // "randomly" select one of the valid tips

    flag.no_tip = true;
    sam_destroy_deep_tip();
}

// called when uncompressing with genounzip (not when just testing)
void tip_print_genounzip (void)
{
    if (!license_is_activated()) 
        TIP ("Use genozip to compress FASTQ, BAM, VCF and other genomic files. %s", WEBSITE_QUICK_GUIDE);

    flag.no_tip = true;
}