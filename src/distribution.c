// ------------------------------------------------------------------
//   distribution.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "version.h"

rom get_distribution (void)
{
    return version_is_devel() ? "devel" : DISTRIBUTION;
}

bool dist_is_conda        (void) { return !strcmp (get_distribution(), "conda"        ); }
bool dist_is_github       (void) { return !strcmp (get_distribution(), "github"       ); }
bool dist_is_installforge (void) { return !strcmp (get_distribution(), "InstallForge" ); }
