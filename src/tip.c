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

// called in ZIP, or sometimes in PIZ if ZIP ends to non-returnable --test. 
void tip_print (void)
{
    #define E(dt) dt_encountered[DT_##dt]

    if (flag.no_tip || // protect from spamming the user with more than one tip
        !license_allow_tip()) return; 

    StrTextLong acadmic_notice = license_academic_tip(); // empty unless acadmic license

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
    valid_tips[n++] = _TIP "Want to understand how Genozip works? See the paper: " PAPER2 "\n";
    valid_tips[n++] = _TIP "use 'genocat --downsample=𝑟𝑎𝑡𝑒' to downsample your data, see: " WEBSITE_DOWNSAMPLING "\n";
    valid_tips[n++] = _TIP "use --password=𝑝𝑎𝑠𝑠𝑤𝑜𝑟𝑑 to also encrypt the data. See: " WEBSITE_ENCRYPTION "\n";
    valid_tips[n++] = _TIP "use 'genozip --tar mydata.tar 𝑓𝑖𝑙𝑒𝑠 𝑎𝑛𝑑 𝑑𝑖𝑟𝑒𝑐𝑡𝑜𝑟𝑖𝑒𝑠' to archive entire directories. See: " WEBSITE_ARCHIVING "\n";
    valid_tips[n++] = _TIP "genozip files are an excellent way to share and publish data – uncompressing genozip files is always free\n";
    valid_tips[n++] = _TIP "use, for example, 'genozip ftp://mysite/sample.fastq.gz' to compress a file directly from a URL. See: " WEBSITE_COMPRESS_URL "\n";

    if (!(E(SAM) || E(BAM) || E(VCF) || E(BCF)) && (license_is_academic() || license_is_eval()))
        valid_tips[n++] = "Please take a moment now to make a note to not forget to cite Genozip:\n" PAPER2_CITATION "\n";

    if (dist_is_github())  
        valid_tips[n++] = _TIP "Do you like Genozip? Please support it by starring it on github: " GITHUB_REPO "\n";

    if (E(SAM) || E(BAM) || E(FASTQ)) 
        valid_tips[n++] = _TIP "use 'genocat --coverage myfile.genozip' get coverage information for BAM or FASTQ. See: " WEBSITE_COVERAGE "\n";

    if (E(BCF))
        valid_tips[n++] = _TIP "genozip compresses VCF files 5-10X faster than it compresses BCF files";
    
    if ((E(SAM) || E(BAM)) && (license_is_academic() || license_is_eval()))
        valid_tips[n++] = _TIP "Please take a moment now to make a note to not forget to cite Genozip:\n" PAPER3_CITATION "\n";

    bool is_deep = flag.deep || ((IS_PIZ && z_file && (Z_DT(SAM) || Z_DT(BAM))) && z_file->z_flags.dts2_deep);
    
    if (is_deep && (license_is_academic() || license_is_eval()))
        valid_tips[n++] = _TIP "Please take a moment now to make a note to not forget to cite Genozip:\n" PAPER4_CITATION "\n";

    if (!is_deep && (E(SAM) || E(BAM) || E(FASTQ)) && sam_get_deep_tip()) 
        for (int i=0; i < 5; i++) // 5X more likely
            valid_tips[n++] = sam_get_deep_tip();

    if (E(FASTQ) && !flag.bam_assist && !flag.deep)
        for (int i=0; i < 5; i++) // 5X more likely
            valid_tips[n++] = _TIP "use '--bamass "_BAMFILE"' to allow genozip to consult the BAM alignments corresponding to this FASTQ file for faster and better compression. See: " WEBSITE_BAMASS "\n";

    if ((E(VCF) || E(BCF)) && (license_is_academic() || license_is_eval()))
        valid_tips[n++] = _TIP "Please take a moment now to make a note to not forget to cite Genozip:\n" PAPER1_CITATION "\n";

    if (!flag.optimize) {
        if (E(SAM) || E(BAM) || E(CRAM)) valid_tips[n++] = _TIP "--optimize allows genozip to make minor modifications to the data that usually have no impact on downstream analysis, for better compression, see: " WEBSITE_OPT_BAM "\n";
        else if (E(VCF) || E(BCF))       valid_tips[n++] = _TIP "--optimize allows genozip to make minor modifications to the data that usually have no impact on downstream analysis, for better compression, see: " WEBSITE_OPT_VCF "\n";
        else if (E(FASTQ))               valid_tips[n++] = _TIP "--optimize allows genozip to make minor modifications to the data that usually have no impact on downstream analysis, for better compression, see: " WEBSITE_OPT_FASTQ "\n";
    }

    if (!flag.best && !flag.fast && !flag.low_memory && !flag.make_reference) 
        valid_tips[n++] = _TIP "using --best achieves the best compression (slower)\n";

    if (license_is_eval() && !flag.show_stats)
        valid_tips[n++] = _TIP "add --stats to see detailed compression statistics\n";

    if (license_is_eval() || license_is_standard() || license_is_enterprise())
        valid_tips[n++] = _TIP "\"genozip --sendto\" lets your clients compress files using 𝑦𝑜𝑢𝑟 Genozip license ahead of sending them to you. See: " WEBSITE_PREMIUM "\n";

    if (arch_get_max_resident_set() > 100 GB || flag.is_windows || flag.is_wsl || flag.is_mac)
        valid_tips[n++] = _TIP "using --low-memory reduces RAM consumption at the expense of compression size and time\n";

    iprintf ("\n%s\n", valid_tips[clock() % n]); // "randomly" select one of the valid tips

    flag.no_tip = true;
    sam_destroy_deep_tip();
}
