// ------------------------------------------------------------------
//   help-text.h
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "website.h"

static rom help_genozip[] = {
    "",
    "Try (for example): genozip myfile.bam",
    "",
    "Please see the genozip manual here: "WEBSITE_GENOZIP
};

static rom help_genounzip[] = {
    "",
    "Try (for example): genounzip myfile.bam.genozip",
    "or:                genounzip --replace myfile.bam.genozip",
    "",
    "Please see the genounzip manual here: "WEBSITE_GENOUNZIP
};

static rom help_genocat[] = {
    "",
    "Please see the genocat manual here: "WEBSITE_GENOCAT
};

static rom help_genols[] = {
    "",
    "Please see the genols manual here: "WEBSITE_GENOLS
};

static rom help_attributions[] = {
    "",
    "Please see the attributions list here: "WEBSITE_ATTRIBUTIONS
};

static rom help_footer[] = {
    "",
    "Technical questions, bug reports and feature requests: " EMAIL_SUPPORT,
    "License inquiries: " EMAIL_SALES,
    "",
    "THIS SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS, COPYRIGHT HOLDERS OR DISTRIBUTORS OF THIS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
};

