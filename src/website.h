// ------------------------------------------------------------------
//   website.h
//   Copyright (C) 2021-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#define GENOZIP_URL "https://genozip.com"
#define GENOZIP_WWW_URL "https://www.genozip.com"
#define REPO        "/divonlan/genozip"
#define GITHUB_RAW  "https://raw.githubusercontent.com" REPO "/master"
#define GITHUB_REPO "https://github.com" REPO

#define WEBSITE_QUICK_GUIDE  GENOZIP_URL "/compression"
#define WEBSITE_GENOZIP      GENOZIP_URL "/genozip"
#define WEBSITE_GENOUNZIP    GENOZIP_URL "/genounzip"
#define WEBSITE_GENOCAT      GENOZIP_URL "/genocat"
#define WEBSITE_GENOLS       GENOZIP_URL "/genols"
#define WEBSITE_DIGEST       GENOZIP_URL "/digest"
#define WEBSITE_ARCHIVING    GENOZIP_URL "/archiving"
#define WEBSITE_LICENSE      GENOZIP_URL "/license"
#define WEBSITE_GET_GENOZIP  GENOZIP_URL "/get-genozip"
#define WEBSITE_PUBLICATIONS GENOZIP_URL "/publications"
#define WEBSITE_INSTALLING   GENOZIP_URL "/installing"
#define WEBSITE_USING_ON_HPC GENOZIP_URL "/using-on-hpc"
#define WEBSITE_ATTRIBUTIONS GENOZIP_URL "/attributions"
#define WEBSITE_TELEMETRY    GENOZIP_URL "/telemetry"
#define WEBSITE_COVERAGE     GENOZIP_URL "/coverage"
#define WEBSITE_DOWNSAMPLING GENOZIP_URL "/downsampling"
#define WEBSITE_ENCRYPTION   GENOZIP_URL "/encryption"
#define WEBSITE_INSTITUTIONS GENOZIP_URL "/institutions"
#define WEBSITE_PREMIUM      GENOZIP_URL "/premium"
#define WEBSITE_STUDENT      GENOZIP_URL "/student"
#define WEBSITE_COMPARE      GENOZIP_URL "/compare"

#define EMAIL_SUPPORT  "support@genozip.com"
#define EMAIL_SALES    "sales@genozip.com"
#define EMAIL_REGISTER "register@genozip.com"

#define GITHUB_LATEST_RELEASE    GITHUB_REPO "/releases/latest"
#define GITHUB_INSTALLERS        GITHUB_RAW  "/installers"  
#define GITHUB_WINDOWS_INSTALLER GITHUB_INSTALLERS "/genozip-installer.exe"
#ifdef __linux__
#define TARBALL_DIRNAME          "genozip-linux-x86_64" /* directory contained in tarball - defined in Makefile */
#elif defined(__APPLE__) && defined(__x86_64__)
#define TARBALL_DIRNAME          "genozip-osx-x86"
#elif defined(__APPLE__) && defined(__aarch64__)
#define TARBALL_DIRNAME          "genozip-osx-arm"
#endif

#define TARBALL_NAME             TARBALL_DIRNAME ".tar"
#define GITHUB_GENOZIP_TARBALL   GITHUB_INSTALLERS "/" TARBALL_NAME
// #define WINDOWS_UPDATE_NAME      "genozip.windows"
// #define GITHUB_WINDOWS_UPDATE    GITHUB_INSTALLERS "/" WINDOWS_UPDATE_NAME // .windows and not .exe to avoid antivirus blocking updates
#define GITHUB_LICENSE_TXT       GITHUB_RAW "/LICENSE.txt"

#define PAPER1 "https://www.researchgate.net/publication/341408805_genozip_a_fast_and_efficient_compression_tool_for_VCF_files"
#define PAPER1_CITATION "Lan, D., et al. (2020) genozip: a fast and efficient compression tool for VCF files, Bioinformatics, 36, 4091-4092"

#define PAPER2 "https://www.researchgate.net/publication/349347156_Genozip_-_A_Universal_Extensible_Genomic_Data_Compressor"
#define PAPER2_CITATION "Lan, D., et al. (2021) Genozip: a universal extensible genomic data compressor, Bioinformatics, 37, 2225-2230"

#define PAPER3 "https://www.researchgate.net/publication/363555511_Genozip_14_-_advances_in_compression_of_BAM_and_CRAM_files"
#define PAPER3_CITATION "Lan, D., et al. (2022) Genozip 14 - advances in compression of BAM and CRAM files (pre-print), doi: https://doi.org/10.1101/2022.09.12.507582"

#define PAPER4 "https://www.researchgate.net/publication/372195570_Deep_FASTQ_and_BAM_co-compression_in_Genozip_15"
#define PAPER4_CITATION "Lan, D., et al. (2023) Deep FASTQ and BAM co-compression in Genozip 15 (pre-print), doi: https://doi.org/10.1101/2023.07.07.548069"
