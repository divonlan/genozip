// ------------------------------------------------------------------
//   help-text.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "website.h"

static const char *help_genozip[] = {
    "",
    "Try (for example): genozip myfile.bam",
    "",
    "Please see the genozip manual here: "WEBSITE_GENOZIP
};

static const char *help_genozip_developer[] = {
    "",
    "Options useful mostly for developers of genozip:",
    "",
    "Please see "WEBSITE_ADVANCED
};

static const char *help_genounzip[] = {
    "",
    "Please see the genounzip manual here: "WEBSITE_GENOUNZIP
};

static const char *help_genocat[] = {
    "",
    "Please see the genocat manual here: "WEBSITE_GENOCAT
};

static const char *help_genols[] = {
    "",
    "Please see the genols manual here: "WEBSITE_GENOLS
};

static const char *help_footer[] = {
    "",
    "THIS SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.",
    "",
    "Technical questions, bug reports and feature requests: " EMAIL_SUPPORT,
    "Commercial license inquiries: " EMAIL_SALES,
    "Requests for support for compression of additional public or proprietary genomic file formats: " EMAIL_SALES,
    "",
    "Citing: Do you find Genozip useful? Please cite:",
    "  Lan, D., et al. (2021) Genozip: a universal extensible genomic data compressor. Bioinformatics, 37, 2225-2230, https://doi.org/10.1093/bioinformatics/btab102",
    "  Lan, D., et al. (2020) genozip: a fast and efficient compression tool for VCF files. Bioinformatics, 36, 4091-4092, https://doi.org/10.1093/bioinformatics/btaa290",
    "",
};

