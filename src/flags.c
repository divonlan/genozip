// ------------------------------------------------------------------
//   flags.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <getopt.h>
#ifndef _WIN32
#include <sys/types.h>
#include <signal.h>
#include <errno.h>
#endif
#include "genozip.h"
#include "flags.h"
#include "data_types.h"
#include "filename.h"
#include "file.h"
#include "regions.h"
#include "dict_id.h"
#include "crypt.h"
#include "strings.h"
#include "vcf.h"
#include "fastq.h"
#include "stream.h"
#include "bgzf.h"
#include "bases_filter.h"
#include "website.h"
#include "reference.h"
#include "license.h"
#include "tar.h"
#include "biopsy.h"
#include "endianness.h"
#include "stats.h"
#include "version.h"
#include "arch.h"
#include "user_message.h"

// flags - factory default values (all others are 0)
Flags flag = { 
#ifdef DEBUG
    .debug        = true,
#endif
#ifdef __LITTLE_ENDIAN__
    .is_lten      = true,
#endif
#ifdef _WIN32
    .is_windows   = true,
#elif defined __APPLE__
    .is_mac       = true,
#elif defined __linux__
    .is_linux     = true,
#endif
    .out_dt       = DT_NONE, 
    .bgzf         = BGZF_NOT_INITIALIZED, 
    .lines_first  = NO_LINE, 
    .lines_last   = NO_LINE,
    .biopsy_line  = { .line_i = NO_LINE },
    .show_stats_comp_i = COMP_NONE,
    .one_vb_comp_i     = COMP_NONE,
    .show_time_comp_i  = COMP_NONE,
};

bool option_is_short[256] = { }; // indexed by character of short option.
FILE *info_stream = 0;      // stdout, stderr or log file - where non-error messages should go
bool is_info_stream_terminal = true;

static Buffer command_line    = {},
              debugger_params = {};

static pid_t pipe_in_pid = 0;
static char pipe_in_process_name[100] = "";

#define MAX_LINE ((int64_t)1 << 62)
static void flags_set_lines (rom optarg)
{
    ASSINP0 (flag.lines_first==NO_LINE && !flag.tail, "Only one of these can be used: --lines, --head, --tail");

    flag.lines_first = 0; // defaults
    flag.lines_last  = MAX_LINE;

    rom hyphen = strchr (optarg, '-');
    if (hyphen) {
        if (hyphen > optarg) { // start exists
            ASSINP (str_get_int_range64 (optarg, hyphen-optarg, 1, MAX_LINE, &flag.lines_first),
                    "--lines: start_line=%.*s must be a positive integer", (int)(hyphen-optarg), optarg);
            flag.lines_first--; // 0 based, but command line option is 1-based
        }

        if (*(hyphen+1)) { // end exists (no \0 after hyphen)
            ASSINP (str_get_int_range64 (hyphen+1, 0, 1, MAX_LINE, &flag.lines_last),
                    "--lines: end_line=%s must be a positive integer", hyphen+1);
            flag.lines_last--; // 0 based, but command line option is 1-based
        }

        ASSINP (flag.lines_first <= flag.lines_last, "--lines <start-line>-<end-line>: start_line=%"PRId64" cannot be larger than end_line=%"PRId64, 
                flag.lines_first+1, flag.lines_last+1);
    }
    else { // no hyphen - 10 lines
        ASSINP (str_get_int_range64 (optarg, 0, 1, MAX_LINE, &flag.lines_first),
                "--lines: line=%s must be a positive integer", optarg);
        flag.lines_first--; // 0 based, but command line option is 1-based
        flag.lines_last = flag.lines_first+9;
    }
}
static void flags_set_out_filename (char *optarg)
{
    int len = strlen (optarg);

    // case: a final / means this is a directory (possibly to be created)
    if (optarg[len-1] == '/') {
        flag.out_dirname = CALLOC (len);
        memcpy ((char*)flag.out_dirname, optarg, len-1);
        file_mkdir (flag.out_dirname);
    }

    // case: a existing directory name
    else if (file_is_dir (optarg)) {
        flag.out_dirname = CALLOC (len+1);
        memcpy ((char*)flag.out_dirname, optarg, len);
    }

    // case: not a directory - treat as a filename
    else
        flag.out_filename = optarg;

    // in genocat, cannot specify just a directory
    ASSINP (!is_genocat || flag.out_filename, "Error: %s is a directory", flag.out_dirname);
}

static void flags_set_head (rom optarg)
{
    flag.explicit_head = true;

    if (!optarg) 
        flags_set_lines ("1-10"); // default

    else {        
        ASSINP0 (str_get_int_range64 (optarg, 0, 1, MAX_LINE, NULL),
                 "The optional argument of --head must be a positive integer value");

        char s[strlen(optarg) + 4];
        sprintf (s, "1-%s", optarg);
        flags_set_lines (s);
    }
}

static void flags_set_tail (rom optarg)
{
    ASSINP0 (flag.lines_first==NO_LINE && !flag.tail, "Only one of these can be used: --lines, --head, --tail");

    if (!optarg) 
        flag.tail = 10; // default

    else 
        ASSINP0 (str_get_int_range64 (optarg, 0, 1, MAX_LINE, &flag.tail),
                "The optional argument of --tail must be a positive integer value");
}

static void flags_set_show_time (rom optarg)
{
    int64_t value;

    // unrestricted show time
    if (!optarg) {
        flag.show_time = ""; 
        flag.show_time_comp_i = COMP_ALL;
    }

    // single component
    else if (str_get_int (optarg, strlen (optarg), &value)) {
        flag.show_time = "";
        flag.show_time_comp_i = value;
    }

    // single ProfileRec
    else {
        flag.show_time = optarg; 
        flag.show_time_comp_i = COMP_ALL;
    }
}

static void flags_set_explicit_bgzf (rom level_str)
{
    if (level_str && !strcmp (level_str, "exact"))
        flag.bgzf = BGZF_BY_ZFILE; 
    
    else 
        ASSINP0 (level_str && str_get_int_range32 (level_str, strlen (level_str), 0, 12, &flag.bgzf),
                "--bgzf expects a value between 0 (no compression) and 12 (best, yet slowest, compression), or \"exact\" (without the quotes: try to recover compression level of original file). If you're not sure what value to choose, 6 is a popular option.");

    // note: flag.bgzf is immutable from now onwards
}

static void flags_set_downsample (rom optarg)
{
    char *comma = strchr (optarg, ',');
    unsigned downsample_len = comma ? (comma - optarg) : strlen(optarg);

    ASSINP (str_get_int_range32 (optarg, downsample_len, 2, 1000000000, (int32_t *)&flag.downsample), 
            "--downsample: bad rate \"%.*s\", expecting an integer from 2 to 1000000000", downsample_len, optarg);

    ASSINP (!comma || str_get_int_range32 (comma+1, 0, 0, flag.downsample-1, (int32_t *)&flag.shard), 
            "--downsample: bad shard \"%s\", expecting an integer from 0 to %u", comma+1, (int)flag.downsample-1);
}

static void flags_set_file_from (rom optarg)
{
    ASSINP0 (!flag.files_from, "--files-from (or -T) can be only be used once on the command line");

    flag.files_from = optarg;
}

static void flags_set_grep (rom optarg)
{
    ASSINP0 (!flag.grep, "--grep or --grep-w can be only be used once on the command line");
    
    flag.grep = optarg; 
    flag.grep_len = strlen (flag.grep);
}

static void flag_set_interleaved (rom optarg)
{
    ASSINP0 (!flag.interleaved, "--interleaved can be only be used once on the command line");

    if (!optarg || str_case_compare (optarg, "both", NULL)) flag.interleaved = INTERLEAVE_BOTH; // default
    else if (str_case_compare (optarg, "either", NULL))     flag.interleaved = INTERLEAVE_EITHER;
    else ABORTINP0 ("Invalid value for --interleaved : accepted values are 'either' or 'both' (default: 'both')");
}

// parse argument --biopsy_line=1/1000  <-- first number is vb_i, second is line_i within VB
static void flag_set_biopsy_line (rom optarg)
{
    ASSINP0 (optarg, "--biopsy_line expects an argument, eg: --biopsy_line=1/1000 (meaning: vblock_i=1, line_i=1000 (0-based line within VB))");

    str_split_ints (optarg, optarg ? strlen (optarg) : 0, 2, '/', arg, true);
    ASSINP0 (n_args==2, "--biopsy_line expects an argument, eg: --biopsy_line=1/1000 (meaning: vblock_i=1, line_i=1000 (0-based line within VB))");

    flag.biopsy_line = (struct biopsy_line){ .vb_i = args[0], .line_i = args[1] };
}

static void flag_set_show_deep (rom optarg)
{
    if (optarg) {
        flag.show_deep = 2;

        str_split_ints (optarg, strlen (optarg), 3, ',', hash, true);
        ASSINP0 (n_hashs, "--show-deep takes an optional argument of deep_hash - format is 3 comma-seperator integers: qname_hash,seq_hash,qual_hash");

        flag.debug_deep_hash = (DeepHash){ hashs[0], hashs[1], hashs[2] };
    }
    else
        flag.show_deep = 1;
}

static void flag_set_stats (StatsType stats_type, rom optarg)
{
    // if -W or -w appear multiple times, we print the stats header to stderr to allow use with "| grep"
    flag.show_stats = (ABS(flag.show_stats) == stats_type) ? -stats_type : stats_type;

    if (optarg)
        flag.show_stats_comp_i = atoi (optarg);
}

static void flag_set_show_containers (rom optarg)
{
    if (!optarg)
        flag.show_containers = SHOW_CONTAINERS_ALL_VBs;

    else if (!(flag.show_containers = atoi (optarg))) { // only a specific vblock_i
        flag.show_containers = SHOW_CONTAINERS_ALL_VBs; 
        flag.dict_id_show_containers = dict_id_make (optarg, strlen (optarg), DTYPE_PLAIN); // specific field
    }
}

void flags_init_from_command_line (int argc, char **argv)
{
    // process command line options
    while (1) {

        #define _i  {"input",            required_argument, 0, 'i'                    }
        #define _I  {"input-size",       required_argument, 0, 'I'                    }
        #define _d  {"decompress",       no_argument,       &command,             PIZ }
        #define _f  {"force",            no_argument,       &flag.force,            1 }
        #define _h  {"help",             optional_argument, 0, 140                    }
        #define _l  {"list",             no_argument,       &flag.list,             1 }
        #define _L1 {"license",          optional_argument, 0, 'L'                    } // US spelling
        #define _L2 {"licence",          optional_argument, 0, 'L'                    } // British spelling
        #define _Lf {"licfile",          required_argument, 0, 26                     }
        #define _qq {"no-tip",           no_argument,       &flag.no_tip,           1 } // note: supported in genocat/unzip as well, in case called from ref_fasta_to_ref
        #define _Q  {"noisy",            no_argument,       &flag.noisy,            1 }
        #define _q  {"quiet",            no_argument,       &flag.quiet,            1 }
        #define _DL {"replace",          no_argument,       &flag.replace,          1 }
        #define _nb {"no-bgzf",          no_argument,       &flag.no_bgzf,          1 }
        #define _nz {"no-zriter",        no_argument,       &flag.no_zriter,        1 }
        #define _nc {"no-cache",         no_argument,       &flag.no_cache,         1 }
        #define _nu {"no-upgrade",       no_argument,       &flag.no_upgrade,       1 }
        #define _hc {"hold-cache",       required_argument, 0, 145                    } // undocumented
        #define _V  {"version",          no_argument,       &command,         VERSION }
        #define _z  {"bgzf",             required_argument, 0, 'z'                    }
        #define _zb {"bam",              no_argument,       &flag.out_dt,      DT_BAM }
        #define _zB {"BAM",              no_argument,       &flag.out_dt,      DT_BAM }
        #define _zs {"sam",              no_argument,       &flag.out_dt,      DT_SAM }
        #define _zS {"SAM",              no_argument,       &flag.out_dt,      DT_SAM }
        #define _zq {"fastq",            optional_argument, 0, 130                    }
        #define _zQ {"FASTQ",            optional_argument, 0, 130                    }
        #define _zf {"fq",               optional_argument, 0, 130                    }
        #define _zF {"FQ",               optional_argument, 0, 130                    }
        #define _zc {"bcf",              no_argument,       &flag.out_dt,      DT_BCF }
        #define _zC {"BCF",              no_argument,       &flag.out_dt,      DT_BCF }
        #define _zv {"vcf",              no_argument,       &flag.out_dt,      DT_VCF }
        #define _zV {"VCF",              no_argument,       &flag.out_dt,      DT_VCF }
        #define _m  {"md5",              no_argument,       &flag.md5,              1 }
        #define _t  {"test",             no_argument,       &flag.test,             1 }
        #define _Nt {"no-test",          no_argument,       &flag.no_test,          1 }
        #define _fa {"fast",             no_argument,       &flag.fast,             1 }
        #define _bs {"best",             no_argument,       &flag.best,             1 }
        #define _lm {"low-memory",       no_argument,       &flag.low_memory,       1 }
        #define _9  {"optimize",         no_argument,       &flag.optimize,         1 } // US spelling
        #define _8  {"optimise",         no_argument,       &flag.optimize,         1 } // British spelling
        #define _9s {"optimize-sort",    no_argument,       &flag.optimize_sort,    1 }
        #define _8s {"optimise-sort",    no_argument,       &flag.optimize_sort,    1 }
        #define _9P {"optimize-phred",   no_argument,       &flag.optimize_phred,   1 }
        #define _8P {"optimise-phred",   no_argument,       &flag.optimize_phred,   1 }
        #define _9V {"optimize-VQSLOD",  no_argument,       &flag.optimize_VQSLOD,  1 }
        #define _8V {"optimise-VQSLOD",  no_argument,       &flag.optimize_VQSLOD,  1 }
        #define _9Q {"optimize-QUAL",    no_argument,       &flag.optimize_QUAL,    1 } 
        #define _8Q {"optimise-QUAL",    no_argument,       &flag.optimize_QUAL,    1 } 
        #define _9f {"optimize-Vf",      no_argument,       &flag.optimize_Vf,      1 }
        #define _8f {"optimise-Vf",      no_argument,       &flag.optimize_Vf,      1 }
        #define _9Z {"optimize-ZM",      no_argument,       &flag.optimize_ZM,      1 }
        #define _8Z {"optimise-ZM",      no_argument,       &flag.optimize_ZM,      1 }
        #define _9D {"optimize-DESC",    no_argument,       &flag.optimize_DESC,    1 }
        #define _8D {"optimise-DESC",    no_argument,       &flag.optimize_DESC,    1 }
        #define _9G {"GL-to-PL",         no_argument,       &flag.GL_to_PL,         1 }
        #define _9g {"GP-to-PP",         no_argument,       &flag.GP_to_PP,         1 }
        #define _al {"add-line-numbers", no_argument,       &flag.add_line_numbers, 1 }
        #define _Sd {"secure-DP",        no_argument,       &flag.secure_DP,        1 }
        #define _pe {"pair",             no_argument,       (int*)&flag.pair,  PAIRED } 
        #define _DP {"deep",             no_argument,       &flag.deep,             1 } 
        #define _St {"sendto",           required_argument, 0, 151                    }
        #define _um {"user-message",     required_argument, 0, 152                    }
        #define _th {"threads",          required_argument, 0, '@'                    }
        #define _u  {"prefix",           required_argument, 0, 'u'                    }
        #define _o  {"output",           required_argument, 0, 'o'                    }
        #define _T  {"files-from",       required_argument, 0, 'T'                    }
        #define _TT {"tar",              required_argument, 0, 27                     }
        #define _p  {"password",         required_argument, 0, 'p'                    }
        #define _B  {"vblock",           required_argument, 0, 'B'                    }
        #define _r  {"regions",          required_argument, 0, 'r'                    }
        #define _R  {"regions-file",     required_argument, 0, 'R'                    }
        #define _qf {"qnames-file",      required_argument, 0, 144                    }
        #define _qF {"qname-file",       required_argument, 0, 144                    }
        #define _Qf {"qnames",           required_argument, 0, 149                    }
        #define _QF {"qname",            required_argument, 0, 149                    }
        #define _SF {"seqs-file",        required_argument, 0, 147                    }
        #define _Rg {"gpos",             no_argument,       &flag.gpos,             1 }
        #define _s  {"samples",          required_argument, 0, 's'                    }
        #define _sf {"FLAG",             required_argument, 0, 17                     }
        #define _sq {"MAPQ",             required_argument, 0, 18                     }
        #define _il {"interleaved",      optional_argument, 0, 29                     } // both --interleave and --interleaved will be accepted
        #define _e  {"reference",        required_argument, 0, 'e'                    }
        #define _E  {"REFERENCE",        required_argument, 0, 'E'                    }
        #define _b  {"bytes",            no_argument,       &flag.bytes,            1 }
        #define _LC {"cache",            no_argument,       &flag.ls_cache,         1 }
        #define _me {"make-reference",   no_argument,       &flag.make_reference,   1 }
        #define _x  {"index",            no_argument,       &flag.index_txt,        1 }
        #define _D  {"subdirs"  ,        no_argument,       &flag.subdirs,          1 }
        #define _g  {"grep",             required_argument, 0, 25                     }
        #define _gw {"grep-w",           required_argument, 0, 'g'                    }
        #define _n  {"lines",            required_argument, 0, 'n'                    }
        #define _nh {"head",             optional_argument, 0, 22                     }
        #define _nt {"tail",             optional_argument, 0, 23                     }
        #define _G  {"drop-genotypes",   no_argument,       &flag.drop_genotypes,   1 }
        #define _H1 {"no-header",        no_argument,       &flag.no_header,        1 }
        #define _1  {"header-one",       no_argument,       &flag.header_one,       1 }
        #define _H0 {"header-only",      no_argument,       &flag.header_only,      1 }
        #define _H2 {"seq-only",         no_argument,       &flag.seq_only,         1 }
        #define _H3 {"qual-only",        no_argument,       &flag.qual_only,        1 }
        #define _aa {"allow-ambiguous",  no_argument,       &flag.allow_ambiguous,  1 } 
        #define _IU {"IUPAC",            required_argument, 0, 24                     }
        #define _iu {"bases",            required_argument, 0, 24                     }
        #define _GT {"GT-only",          no_argument,       &flag.gt_only,          1 }
        #define _Gt {"gt-only",          no_argument,       &flag.gt_only,          1 }
        #define _So {"snps-only",        no_argument,       &flag.snps_only,        1 }
        #define _Io {"indels-only",      no_argument,       &flag.indels_only,      1 }
        #define _ds {"downsample",       required_argument, 0, 9                      }
        #define _PG {"no-PG",            no_argument,       &flag.no_pg,            1 }
        #define _pg {"no-pg",            no_argument,       &flag.no_pg,            1 }
        #define _fs {"sequential",       no_argument,       &flag.sequential,       1 }  
        #define _rg {"register",         optional_argument, 0, 28                     }
        #define _sL {"show-lines",       no_argument,       &flag.show_lines,       1 } 
        #define _ss {"stats",            optional_argument, 0, 'w',                   } 
        #define _SS {"STATS",            optional_argument, 0, 'W'                    } 
        #define _lc {"list-chroms",      no_argument,       &flag.list_chroms,      1 } // identical to --show-dict=CHROM 
        #define _lh {"chroms",           no_argument,       &flag.list_chroms,      1 } // identical to --show-dict=CHROM 
        #define _lH {"contigs",          no_argument,       &flag.list_chroms,      1 } 
        #define _s2 {"show-b250",        optional_argument, 0, 2,                     }
        #define _sd {"show-dict",        optional_argument, 0, 3                      }
        #define _s7 {"dump-b250",        required_argument, 0, 5                      }
        #define _S7 {"dump-local",       required_argument, 0, 6                      }
        #define _S0 {"show-singletons",  required_argument, 0, 150                    }
        #define _S8 {"show-counts",      required_argument, 0, 16                     }
        #define _S9 {"dump-section",     required_argument, 0, 7                      }        
        #define _sa {"show-alleles",     no_argument,       &flag.show_alleles,     1 }
        #define _st {"show-time",        optional_argument, 0, 1                      } 
        #define _sm {"show-memory",      no_argument,       &flag.show_memory,      1 } 
        #define _sh {"show-headers",     optional_argument, 0, 10                     } 
        #define _si {"show-index",       no_argument,       &flag.show_index,       1 } 
        #define _sG {"show-ranges",      no_argument,       &flag.show_ranges,      1 } 
        #define _Si {"show-ref-index",   no_argument,       &flag.show_ref_index,   1 } 
        #define _Sh {"show-ref-hash" ,   no_argument,       &flag.show_ref_hash,    1 } 
        #define _sr {"show-gheader",     optional_argument, 0, 135                    }  
        #define _DT {"show-data-type",   no_argument,       &flag.show_data_type,   1 }  
        #define _sT {"show-threads",     no_argument,       &flag.show_threads,     1 }  
        #define _wM {"show-wrong-MD",    no_argument,       &flag.show_wrong_md,    1 }  
        #define _wm {"show-wrong-XG",    no_argument,       &flag.show_wrong_xg,    1 }  
        #define _WM {"show-wrong-XM",    no_argument,       &flag.show_wrong_xm,    1 }  
        #define _WB {"show-wrong-XB",    no_argument,       &flag.show_wrong_xb,    1 }  
        #define _su {"show-uncompress",  no_argument,       &flag.show_uncompress,  1 }  
        #define _sv {"show-vblocks",     optional_argument, 0, 141                    }  
        #define _sn {"show-filename",    no_argument,       &flag.show_filename,    1 }  
        #define _ai {"analyze-insertions", no_argument,     &flag.analyze_ins,      1 }  
        #define _ov {"one-vb",           required_argument, 0, 8                      }  
        #define _R1 {"R1",               no_argument,       &flag.one_component,    1 }  // comp_i+1
        #define _R2 {"R2",               no_argument,       &flag.one_component,    2 }  
        #define _RX {"R",                required_argument, 0, 146                    }  
        #define _sR {"show-reference",   no_argument,       &flag.show_reference,   1 }  
        #define _sC {"show-ref-contigs", no_argument,       &flag.show_ref_contigs, 1 }  
        #define _rI {"show-ref-iupacs",  no_argument,       &flag.show_ref_iupacs,  1 }  
        #define _rA {"show-chrom2ref",   no_argument,       &flag.show_chrom2ref,   1 }  
        #define _rS {"show-ref-seq",     no_argument,       &flag.show_ref_seq,     1 }  
        #define _cn {"show-containers",  optional_argument, 0, 132                    }  
        #define _dd {"debug",            no_argument,       &flag.debug,            1 }
        #define _dg {"debug-seg",        optional_argument, 0, 133                    }  
        #define _hC {"show-txt-contigs", no_argument,       &flag.show_txt_contigs, 1 }
        #define _sI {"show-is-set",      required_argument, 0, '~',                   }  
        #define _sA {"show-aliases",     no_argument,       &flag.show_aliases,     1 }  
        #define _sB {"show-sag",         optional_argument, 0, 136                    }  
        #define _sP {"show-depn",        no_argument,       &flag.show_depn,        1 }  
        #define _sc {"show-codec",       no_argument,       &flag.show_codec,       1 }  
        #define _Sc {"show-cache",       no_argument,       &flag.show_cache,       1 }  
        #define _AL {"show-aligner",     no_argument,       &flag.show_aligner,     1 }  
        #define _sb {"show-bgzf",        no_argument,       &flag.show_bgzf,        1 }
        #define _Sb {"show-gz",          no_argument,       &flag.show_gz,          1 }
        #define _s5 {"show-digest",      no_argument,       &flag.show_digest,      1 }
        #define _S5 {"log-digest",       no_argument,       &flag.log_digest,       1 }
        #define _s6 {"show-plan",        no_argument,       &flag.show_recon_plan,  1 }
        #define _sM {"show-mutex",       optional_argument, 0, 4                      }
        #define _dS {"seg-only",         no_argument,       &flag.seg_only,         1 }  
        #define _bS {"show-bam",         no_argument,       &flag.show_bam,         1 }  
        #define _xt {"xthreads",         no_argument,       &flag.xthreads,         1 }  
        #define _dm {"debug-memory",     optional_argument, 0, 12                     }  
        #define _dp {"debug-progress",   no_argument,       &flag.debug_progress,   1 }  
        #define _dv {"debug-valgrind",   no_argument,       &flag.debug_valgrind,   1 }  
        #define _TR {"debug-tar",        no_argument,       &flag.debug_tar,        1 }  
        #define _dL {"debug-LONG",       no_argument,       &flag.debug_LONG,       1 }  
        #define _dD {"show-qual",        no_argument,       &flag.show_qual,        1 }  
        #define _dq {"debug-qname",      no_argument,       &flag.debug_qname,      1 }  
        #define _dB {"show-buddy",       no_argument,       &flag.show_buddy,       1 }  
        #define _dt {"debug-threads",    no_argument,       &flag.debug_threads,    1 }  
        #define _dw {"debug-stats",      no_argument,       &flag.debug_stats,      1 }  
        #define _dM {"debug-generate",   no_argument,       &flag.debug_generate  , 1 }  
        #define _dr {"debug-recon-size", no_argument,       &flag.debug_recon_size, 1 }  
        #define _dR {"debug-read-ctxs",  no_argument,       &flag.debug_read_ctxs,  1 }  
        #define _dP {"debug-sag",        no_argument,       &flag.debug_sag,        1 }  
        #define _dG {"debug-gencomp",    no_argument,       &flag.debug_gencomp,    1 }  
        #define _dN {"no-gencomp",       no_argument,       &flag.no_gencomp,       1 }  
        #define _dF {"force-gencomp",    no_argument,       &flag.force_gencomp,    1 }  
        #define _DF {"force-deep",       no_argument,       &flag.force_deep,       1 }  
        #define _fP {"force-PLy",        no_argument,       &flag.force_PLy,        1 }  
        #define _dQ {"no-domq",          no_argument,       &flag.no_domqual,       1 }  
        #define _dC {"no-pacb",          no_argument,       &flag.no_pacb,          1 }  
        #define _dH {"no-homp",          no_argument,       &flag.no_homp,          1 }  
        #define _dO {"no-longr",         no_argument,       &flag.no_longr,         1 }  
        #define _dU {"no-smux",          no_argument,       &flag.no_smux,          1 }  
        #define _fQ {"force-domq",       no_argument,       &flag.force_qual_codec, CODEC_DOMQ  }  
        #define _fC {"force-pacb",       no_argument,       &flag.force_qual_codec, CODEC_PACB  }  
        #define _fO {"force-longr",      no_argument,       &flag.force_qual_codec, CODEC_LONGR }  
        #define _fS {"force-smux" ,      no_argument,       &flag.force_qual_codec, CODEC_SMUX  }  
        #define _fH {"force-homp" ,      no_argument,       &flag.force_qual_codec, CODEC_HOMP  }  
        #define _fN {"force-normq",      no_argument,       &flag.force_qual_codec, CODEC_NORMQ }  
        #define _dl {"debug-lines",      no_argument,       &flag.debug_lines,      1 }  
        #define _dc {"verify-codec",     no_argument,       &flag.verify_codec,     1 }          
        #define _oe {"echo",             no_argument,       &flag.echo,             1 }
        #define _dh {"show-hash",        no_argument,       &flag.show_hash,        1 }  
        #define _SH {"show-segconf-has", no_argument,       &flag.show_segconf_has, 1 }  
        #define _ct {"count",            optional_argument, 0, 20                     }  
        #define _SX {"coverage",         optional_argument, 0, 13                     }  
        #define _ix {"idxstats",         no_argument,       &flag.idxstats,         1 }
        #define _vl {"validate",         optional_argument, 0, 19                     }  
        #define _lg {"log",              required_argument, 0, 15                     }  
        #define _bi {"biopsy",           required_argument, 0, 134,                   }
        #define _bl {"biopsy-line",      required_argument, 0, 137,                   }
        #define _sk {"skip-segconf",     no_argument,       &flag.skip_segconf,     1 }
        #define _MR {"match-chrom-to-reference",   no_argument, &flag.match_chrom_to_reference,   1 }
        #define _TL {"truncate-partial-last-line", no_argument, &flag.truncate_partial_last_line, 1 }
        #define _VV {"check-latest",     no_argument,       &flag.check_latest,     1 }
        #define _DV {"debug-latest",     no_argument,       &flag.debug_latest,     1 }
        #define _Dp {"debug-peek",       no_argument,       &flag.debug_peek,       1 }
        #define _Dh {"debug-huffman",    no_argument,       &flag.debug_huffman,    1 }
        #define _sp {"debug-split",      no_argument,       &flag.debug_split,      1 }
        #define _np {"show-snips",       no_argument,       &flag.show_snips,       1 }
        #define _Ds {"submit-stats",     no_argument,       &flag.stats_submit,     1 }
        #define _DS {"debug-submit",     no_argument,       &flag.debug_submit,     1 }
        #define _DD {"debug-debug",      no_argument,       &flag.debug_debug,      1 }
        #define _RC {"recover",          no_argument,       &flag.recover,          1 }
        #define _NE {"no-eval",          no_argument,       &flag.no_eval,          1 }
        #define _Dd {"show-deep",        optional_argument, 0, 139,                   }
        #define _to {"t_offset",         required_argument, 0, 142,                   }
        #define _ts {"t_size",           required_argument, 0, 143,                   }
        #define _lp {"license-prepare",  required_argument, 0, 148,                   }
        #define _00 {0, 0, 0, 0                                                       }

        typedef const struct option Option;
        static Option genozip_lo[]    = { _lg, _i, _I, _d, _f, _h,     _D,    _L1, _L2, _q, _Q, _qq, _t, _Nt, _DL, _nb, _nz, _nc,_nu,  _V, _z,                                                             _m, _th,     _o, _p, _e, _E,                                                                       _H1,                                         _sL, _ss, _SS,      _sd, _sT, _sb, _Sb, _lc, _lh, _lH, _s2, _s7, _S7, _S0, _S8, _S9, _sa, _st, _sm, _sh, _si, _Si, _Sh, _sr, _su, _sv, _sn, _ai,                    _B, _xt, _dm, _dp, _dL, _dD, _dq, _dB, _dt, _dw, _dM, _dr, _dR, _dP, _dG, _dN, _dF, _DF, _dQ, _dH, _dO, _dC, _fQ, _fC, _fO, _fS, _fH, _fN, _dU, _dl, _dc, _dg,      _dh,_dS, _bS, _9, _9s, _9P, _9G, _9g, _9V, _9Q, _9f, _9Z, _9D, _8, _8s, _8P, _8V, _8Q, _8f, _8Z, _8D, _pe, _fa, _bs, _lm,                        _nh, _rg, _sR,      _sC, _hC, _rA, _rI, _rS, _me, _s5, _S5, _sM, _sA, _sB, _sP, _sc, _Sc, _AL, _sI, _cn,                               _s6,          _oe, _aa, _al, _Lf, _dd, _T, _TT, _MR, _TL, _wM, _wm, _WM, _WB, _bi, _bl, _sk, _VV, _DV,           _Ds, _DS, _sp, _DD, _DP, _SH, _Dd,      _to, _ts,      _hc, _dv, _TR, _NE, _lp, _Sd, _St, _um,      _fP, _00 };
        static Option genounzip_lo[]  = { _lg,         _d, _f, _h, _x, _D,    _L1, _L2, _q, _Q, _qq, _t,      _DL,           _nc,      _V, _z,                                                             _m, _th, _u, _o, _p, _e,                                                                                                                        _sL, _ss, _SS, _sG, _sd, _sT, _sb,      _lc, _lh, _lH, _s2, _s7, _S7, _S0, _S8, _S9, _sa, _st, _sm, _sh, _si, _Si, _Sh, _sr, _su, _sv, _sn,      _ov,                   _xt, _dm, _dp,      _dD,      _dB, _dt,                _dR,                                                                                      _dc,                                                                                                                                _lm,                                  _sR,      _sC, _hC, _rA,      _rS, _me, _s5, _S5, _sM, _sA, _sB,           _Sc, _AL, _sI, _cn, _pg,                          _s6,          _oe,                _dd, _T, _TT,                                                        _Dp, _Dh,           _sp, _DD,           _Dd,      _to, _ts, _RC, _hc, _dv, _TR, _NE,                     _np,      _00 };
        static Option genocat_lo[]    = { _lg,         _d, _f, _h,     _D,    _L1, _L2, _q, _Q, _qq,                         _nc,      _V, _z, _zb, _zB, _zs, _zS, _zq, _zQ, _zf, _zF, _zc, _zC, _zv, _zV,     _th,     _o, _p, _e,     _il, _r, _R, _Rg, _qf, _qF, _Qf, _QF, _SF, _s, _sf, _sq, _G, _1, _H0, _H1, _H2, _H3, _Gt, _So, _Io, _IU, _iu, _GT, _sL, _ss, _SS, _sG, _sd, _sT, _sb,      _lc, _lh, _lH, _s2, _s7, _S7, _S0, _S8, _S9, _sa, _st, _sm, _sh, _si, _Si, _Sh, _sr, _su, _sv, _sn,      _ov, _R1, _R2, _RX,    _xt, _dm, _dp,      _dD,      _dB, _dt,                _dR,                                                                                      _dc,      _ds,                                                                                                                      _lm, _fs, _g, _gw, _n, _nt, _nh,      _sR,      _sC, _hC, _rA, _rI, _rS, _me, _s5, _S5, _sM, _sA, _sB,           _Sc, _AL, _sI, _cn, _pg, _PG, _SX, _ix, _ct, _vl, _s6,          _oe,      _al,      _dd, _T,                                                             _Dp, _Dh,           _sp, _DD,           _Dd, _DT,           _RC, _hc, _dv, _TR, _NE,                     _np,      _00 };
        static Option genols_lo[]     = { _lg,             _f, _h,        _l, _L1, _L2, _q,                                            _V,                                                                                  _p,                                                                                                                                                                                                                 _st, _sm,                                                                          _dm,                          _dt,                                                                                                                                                                                                                                                                                                                                 _S5, _sM,                                                                            _b, _LC, _oe,                _dd, _T,                                                                                 _sp, _DD,                                         _dv,      _NE,                               _00 };
        static Option *long_options[NUM_EXE_TYPES] = { genozip_lo, genounzip_lo, genocat_lo, genols_lo }; // same order as ExeType

        // include the option letter here for the short version (eg "-t") to work. ':' indicates an argument.
        static rom short_options[NUM_EXE_TYPES] = { // same order as ExeType
            "z:i:I:cdfhLqQt^Vz:m@:o:p:B:9wWFe:E:C:23K:T:DbX", // genozip (note: includes some genounzip options to be used in combination with -d)
            "cz:fhLqdQt^V@:u:o:p:me:wWxT:D",                  // genounzip
            "z:hLV@:dp:qQ1r:R:s:H1Go:fg:e:E:wWxvk:K:n:yT:D",  // genocat
            "hLVp:qfblT:",                                    // genols
        };

        int option_index = -1;
        int c = getopt_long (argc, argv, short_options[exe_type], long_options[exe_type], &option_index);

        if (c == -1) break; // no more options

        option_is_short[c] = (option_index == -1); // true if user provided a short option - eg -c rather than --stdout

        switch (c) {
            case 140:
                flag.help = optarg;
                goto verify_command;

            case LICENSE :
                flag.lic_param = optarg;
                goto verify_command;

            case PIZ :  
            case VERSION :
verify_command:
                ASSINP (command<0 || command==c, "can't have both -%c and -%c", command, c); 
                command=c; 
                break;

            case '^' : flag.replace       = 1       ; break;
            case '@' : flag.threads_str   = optarg  ; break;
            case '1' : flag.header_one    = 1       ; break;
            case '2' : flag.pair          = PAIRED  ; break;
            case '3' : flag.deep          = 1       ; break;
            case '9' : flag.optimize      = 1       ; break;
            case 'X' : flag.no_test       = 1       ; break; 
            case 'b' : flag.best          = 1       ; break; 
            case 'B' : flag.vblock        = optarg  ; break;
            case 'D' : flag.subdirs       = 1       ; break;  
            case 'e' : ref_set_reference (gref, optarg, REF_EXTERNAL,  true); break;
            case 'E' : ref_set_reference (gref, optarg, REF_EXT_STORE, true); break;
            case 'f' : flag.force         = 1       ; break;
            case 'F' : flag.fast          = 1       ; break;
            case 'g' : flag.grepw = true; // fall-through
            case 25  : flags_set_grep (optarg)      ; break;
            case 'G' : flag.drop_genotypes= 1       ; break;
            case 'H' : flag.no_header     = 1       ; break;
            case 'i' : file_set_input_type (optarg) ; break;
            case 'I' : file_set_input_size (optarg) ; break;
            case 'l' : flag.list          = 1       ; break;
            case 'h' : flag.header_only   = 1       ; break;
            case 'm' : flag.md5           = 1       ; break;
            case 'n' : flags_set_lines (optarg)     ; break;
            case 'o' : flags_set_out_filename(optarg); break;
            case 'p' : crypt_set_password (optarg)  ; break;
            case 'q' : flag.quiet         = 1       ; break;
            case 'Q' : flag.noisy         = 1       ; break;
            case 'r' : flag.regions       = 1       ; regions_add (optarg); break;
            case 'R' : flag.regions       = 1       ; regions_add_by_file (optarg); break;
            case 'T' : flags_set_file_from (optarg) ; break;
            case 's' : flag.samples       = 1       ; vcf_samples_add (optarg); break;
            case 't' : flag.test          = 1       ; break; 
            case 'u' : flag.unbind        = optarg  ; break;
            case 'w' : flag_set_stats (STATS_SHORT, optarg); break;
            case 'W' : flag_set_stats (STATS_LONG,  optarg); break;
            case 'z' : flags_set_explicit_bgzf (optarg)    ; break;
            case 'x' : flag.index_txt     = 1       ; break;
            case '~' : flag.show_is_set   = optarg  ; break;
            case 1   : flags_set_show_time (optarg) ; break; // show_time with or without a specific member of ProfilerRec
            case 2   : if (optarg) flag.dict_id_show_one_b250 = dict_id_make (optarg, strlen (optarg), DTYPE_PLAIN); 
                       else        flag.show_b250 = 1;
                       break;
            case 3   : if (optarg) flag.show_one_dict = optarg;
                       else        flag.show_dict = 1;
                       break;
            case 4   : flag.show_mutex    = optarg ? optarg : (char*)1; break;
            case 5   : flag.dump_one_b250_dict_id   = dict_id_make (optarg, strlen (optarg), DTYPE_PLAIN); break;
            case 6   : flag.dump_one_local_dict_id  = dict_id_make (optarg, strlen (optarg), DTYPE_PLAIN); break;
            case 150 : flag.show_singletons_dict_id = dict_id_make (optarg, strlen (optarg), DTYPE_PLAIN); break;
            case 7   : flag.dump_section  = optarg  ; break;
            case 8   : ASSINP0 (str_get_int_range32 (optarg, 0, 1, 10000000, (int32_t*)&flag.one_vb), 
                                "--one-vb expects a 1-based VBlock number"); break;
            case 9   : flags_set_downsample (optarg); break;
            case 10  : sections_set_show_headers (optarg); break; // +1 so SEC_NONE maps to 0
            case 12  : flag.debug_memory  = optarg ? atoi (optarg) : 1; break;
            case 13  : flag.show_coverage = !optarg                 ? COV_CHROM 
                                          : !strcmp (optarg, "all") ? COV_ALL 
                                          : !strcmp (optarg, "one") ? COV_ONE 
                                          :                           COV_ALL; break;
            case 15  : flag.log_filename  = optarg;  break;
            case 16  : flag.show_one_counts = dict_id_make (optarg, strlen (optarg), DTYPE_PLAIN); break;
            case 131 : flag.show_one_counts = dict_id_make (cSTR("o$TATUS"), DTYPE_PLAIN); break;
            case 17  : sam_set_FLAG_filter (optarg) ; break; // filter by SAM FLAG
            case 18  : sam_set_MAPQ_filter (optarg) ; break; // filter by SAM MAPQ
            case 19  : flag.validate = optarg ? VLD_REPORT_VALID : VLD_REPORT_INVALID ; break;
            case 20  : flag.count = optarg ? CNT_VBs : CNT_TOTAL; break;
            case 22  : flags_set_head (optarg)      ; break;
            case 23  : flags_set_tail (optarg)      ; break;
            case 24  : iupac_set (optarg)           ; break;
            case 26  : license_set_filename (optarg); break;
            case 27  : tar_set_tar_name (optarg)    ; break;
            case 28  : flag.do_register = optarg ? optarg : ""; break;
            case 29  : flag_set_interleaved (optarg); break;
            case 130 : flag.out_dt = DT_FASTQ       ; break; 
            case 132 : flag_set_show_containers (optarg); break; 
            case 133 : flag.debug_seg=1;
                       if (optarg) flag.dict_id_debug_seg = dict_id_make (optarg, strlen (optarg), DTYPE_PLAIN); 
                       break;
            case 134 : biopsy_init (optarg);        ; break;
            case 135 : flag.show_gheader = optarg ? atoi(optarg) : 1; break; // =1 show gheader as in file, =2 show shows section list after possible modiciation by writer_create_plan 
            case 136 : flag.show_sag = optarg ? atoi(optarg)+1 : -1; break;   //-1=show all, >=1 - show grp_i=show_sag-1 
            case 137 : flag_set_biopsy_line (optarg); break;
            case 138 : /* 138 is not currently used - available */; break;
            case 139 : flag_set_show_deep (optarg)  ; break;
            case 141 : flag.show_vblocks = (optarg ? optarg : "") ; break;
            case 142 : flag.t_offset = atoll (optarg); break;
            case 143 : flag.t_size   = atoll (optarg); break;
            case 144 : flag.qnames_file = optarg; break;
            case 149 : flag.qnames_opt = optarg; break;
            case 147 : seq_filter_initialize (optarg); break;
            case 145 : ref_cache_hold (optarg); // doesn't return
            case 146 : flag.one_component = atoi (optarg); break;
            case 148 : license_prepare (optarg); break;
            case 151 : ASSINP (str_get_int_range64 (optarg, strlen (optarg), 1, 0xffffffff, &flag.sendto), "Expecting the value of --sendto=%s to a number", optarg); break;
            case 152 : user_message_init (optarg); break;
            case 0   : break; // a long option that doesn't have short version will land here - already handled so nothing to do
                 
            default  : // unrecognized option 
                fprintf (stderr, "Usage: %s [OPTIONS] filename1 filename2...\nManual page: "GENOZIP_URL"/%s\n", global_cmd, global_cmd);
                exit (EXIT_GENERAL_ERROR);  
        }
    }
    flag.explicit_out_dt = (flag.out_dt >= 0);
    
    flag.to_stdout = is_genocat && !flag.out_filename;

    flag.debug_or_test |= flag.debug; // in case of explicit --debug
    
    command = (command != NO_COMMAND)           ? command // set explictly with -d -V etc
            : is_genols                         ? LIST
            : (is_genocat && flag.show_headers) ? SHOW_HEADERS
            : (is_genounzip || is_genocat)      ? PIZ
            :                                     ZIP;

}

static void flags_warn_if_duplicates (int num_files, rom *filenames)
{
    if (num_files <= 1) return; // nothing to do

    # define BASENAME_LEN 256
    char basenames[num_files * BASENAME_LEN];

    for (unsigned i=0; i < num_files; i++)
        filename_base (filenames[i], false, "", &basenames[i*BASENAME_LEN], BASENAME_LEN);

    qsort (basenames, num_files, BASENAME_LEN, (int (*)(const void *, const void *))strcmp);

    for (unsigned i=1; i < num_files; i++) 
        ASSERTW (strncmp(&basenames[(i-1) * BASENAME_LEN], &basenames[i * BASENAME_LEN], BASENAME_LEN), 
                 "Warning: two files with the same name '%s' - if you later use genounzip, these files will overwrite each other", 
                 &basenames[i * BASENAME_LEN]);
}

// called from main()->flags_update() after flags are assigned, and before analyzing input files
static void flags_test_conflicts (unsigned num_files /* optional */)
{
    #define CONFLICT(flag1, flag2, text1, text2) \
        ASSINP (!(flag1) || !(flag2), "option %s is incompatible with %s", text1, text2);

    #define NEED_DECOMPRESS(flag1, text) \
        ASSINP (!(flag1), "option %s can only be used if --decompress is used too", text);

    CONFLICT (flag.subdirs,     flag.out_filename,   OT("subdirs", "D"),   OT("output", "o")); 
    CONFLICT (flag.to_stdout,   flag.index_txt,      OT("stdout", "c"),    OT("index", "x"));
    CONFLICT (flag.quiet,       flag.noisy,          OT("quiet", "q"),     OT("noisy", "Q"));
    CONFLICT (flag.test,        flag.no_test,        OT("test", "t"),      OT("no-test", "X"));
    CONFLICT (flag.test,        flag.bgzf != BGZF_NOT_INITIALIZED, OT("test", "t"),      OT("bgzf", "z"));
    CONFLICT (flag.dump_one_b250_dict_id.num, flag.dump_one_local_dict_id.num, "--dump-one-b250", "--dump-one-local");

    if (IS_ZIP) {
        #define CONFLICT_NO_DIGEST(f,opt) \
            CONFLICT (flag.test, flag.f, OT("test", "t"), opt); \
            CONFLICT (flag.md5,  flag.f, OT("md5", "m"),  opt);

        CONFLICT_NO_DIGEST (add_line_numbers, "--add-line-numbers");
        CONFLICT_NO_DIGEST (optimize, "optimize");
        CONFLICT_NO_DIGEST (optimize_sort, "optimize-sort");
        CONFLICT_NO_DIGEST (optimize_phred, "optimize-phred");
        CONFLICT_NO_DIGEST (GL_to_PL, "GL-to-PL");
        CONFLICT_NO_DIGEST (GP_to_PP, "GP-to-PP");
        CONFLICT_NO_DIGEST (optimize_VQSLOD, "optimize-VQSLOD");
        CONFLICT_NO_DIGEST (optimize_QUAL, "optimize-QUAL");
        CONFLICT_NO_DIGEST (optimize_Vf, "optimize-Vf");
        CONFLICT_NO_DIGEST (optimize_ZM, "optimize-ZM");
        CONFLICT_NO_DIGEST (optimize_DESC, "--optimize-DESC");
        CONFLICT_NO_DIGEST (match_chrom_to_reference, "--match-chrom-to-reference");

        CONFLICT (flag.test,        flag.biopsy,         OT("test", "t"),      "--biopsy");
        CONFLICT (flag.best,        flag.fast,           OT("best", "b"),      OT("fast", "F"));
        CONFLICT (flag.best,        flag.low_memory,     OT("best", "b"),      "--low-memory");
        CONFLICT (flag.test,        flag.seg_only,       "--seg-only",         OT("test", "t"));
        CONFLICT (flag.test,        flag.make_reference, "--make-reference",   OT("test", "t"));
        CONFLICT (flag.reference,   flag.make_reference, "--make-reference",   OT("reference", "e"));
        CONFLICT (flag.no_gencomp,  flag.force_gencomp,  "--no-gencomp",       "--force-gencomp");
        CONFLICT (flag.biopsy,      flag.out_filename,   "--biopsy",           OT("output", "o"));
        CONFLICT (flag.biopsy_line.line_i>=0, flag.out_filename, "--biopsy-line", OT("output", "o"));
        CONFLICT (flag.deep,        flag.pair,           OT("deep", "3"),      OT("pair", "2"));
        CONFLICT (flag.deep,        flag.subdirs,        OT("deep", "3"),      OT("subdirs", "D"));
        CONFLICT (flag.deep,        tar_zip_is_tar(),    OT("deep", "3"),      "--tar");
        CONFLICT (flag.pair,        tar_zip_is_tar(),    OT("pair", "2"),      "--tar");
        CONFLICT (flag.deep,        flag.optimize_DESC,  OT("deep", "3"),      "--optimize-DESC");

        // see comment in fastq_deep_zip_finalize
        WARN_IF (flag.no_test && flag.deep, "It is not recommended to use %s in combination with %s, as Deep relies on testing to catch rare edge cases", OT("no-test", "X"), OT("deep", "3"));
    }

    if (IS_PIZ) {
        CONFLICT (flag.unbind,      flag.out_filename,   OT("unbind", "u"),    OT("output", "o"));
        CONFLICT (flag.samples,     flag.drop_genotypes, OT("samples", "s"),   OT("drop-genotypes", "G"));
        CONFLICT (flag.one_vb,      flag.interleaved,    "--interleaved",      "--one-vb");
        CONFLICT (flag.one_component, flag.interleaved,  "--interleaved",      "--R1/--R2");
        CONFLICT (flag.snps_only,   flag.indels_only,    "--snps_only",        "--indels-only");
        CONFLICT (flag.header_only, flag.seq_only,       "--seq-only",         "--header-only");
        CONFLICT (flag.header_only, flag.qual_only,      "--qual-only",        "--header-only");
        CONFLICT (flag.seq_only,    flag.qual_only,      "--seq-only",         "--qual-only");
        CONFLICT (flag.header_only, flag.interleaved,    "--interleaved",      "--header-only");
        CONFLICT (flag.regions,     flag.interleaved,    "--interleaved",      "--regions");
        CONFLICT (flag.header_only, flag.no_header==1,   OT("no-header", "H"), "--header-only");
        CONFLICT (flag.no_header,   flag.header_one,     OT("no-header", "H"), OT("header-one", "1"));
        CONFLICT (flag.test,        flag.out_filename,   OT("output", "o"),    OT("test", "t"));
        CONFLICT (flag.test,        flag.replace,        OT("replace", "^"),   OT("test", "t"));
        CONFLICT (flag.show_coverage, flag.idxstats,     "--coverage",         "--idxstats");
        CONFLICT (flag.show_coverage, flag.regions==1,   "--coverage",         OT("regions", "r"));
        CONFLICT (flag.show_coverage, flag.out_filename, "--coverage",         OT("output", "o"));
        CONFLICT (flag.show_coverage, flag.grep,         "--coverage",         OT("grep", "g"));
        CONFLICT (flag.show_coverage, flag.count,        "--coverage",         "--count");
        CONFLICT (flag.show_coverage, flag.test,         "--coverage",         OT("test", "t"));
        CONFLICT (flag.idxstats,    flag.test,           "--idxstats",         OT("test", "t"));
        CONFLICT (flag.idxstats,    flag.count,          "--idxstats",         "--count");
        CONFLICT (flag.idxstats,    flag.show_coverage,  "--idxstats",         "--coverage");
        CONFLICT (flag.count,       flag.out_filename,   "--count",            OT("output", "o"));
        CONFLICT (flag.count,       flag.header_only,    "--count",            OT("header-only", "H"));
        CONFLICT (flag.count,       flag.header_one,     "--count",            OT("header-one", "1"));
        CONFLICT (flag.count,       flag.explicit_head,  "--count",            "--head"); 
        CONFLICT (flag.count,       flag.tail,           "--count",            "--tail"); 
        CONFLICT (flag.count,       flag.lines_first != NO_LINE, "--count",    OT("lines", "n")); 
        CONFLICT (flag.show_stats,  flag.grep,           OT("stats", "w"),     OT("grep", "g"));
        CONFLICT (flag.show_stats,  flag.lines_first != NO_LINE, OT("stats", "w"), OT("lines", "n"));
        CONFLICT (flag.show_stats,  flag.downsample,     OT("stats", "w"),     "--downsample");
        CONFLICT (flag.show_stats,  flag.header_only,    OT("stats", "w"),     "--header-only");
        CONFLICT (flag.show_stats,  flag.MAPQ,           OT("stats", "w"),     "--MAPQ");
        CONFLICT (flag.show_stats,  flag.FLAG,           OT("stats", "w"),     "--FLAG");
        CONFLICT (flag.show_stats,  flag.one_component,  OT("stats", "w"),     "--R1/--R2");
        CONFLICT (flag.show_stats,  flag.one_vb,         OT("stats", "w"),     "--one-vb");
        CONFLICT (flag.show_stats,  flag.sequential,     OT("stats", "w"),     "--sequential");
        CONFLICT (flag.show_stats,  flag.tail,           OT("stats", "w"),     "--tail");
        CONFLICT (flag.show_stats,  flag.drop_genotypes, OT("stats", "w"),     OT("drop-genotypes", "G"));
        CONFLICT (flag.show_stats,  flag.regions,        OT("stats", "w"),     OT("regions", "r"));
        CONFLICT (flag.show_stats,  flag.samples,        OT("stats", "w"),     OT("samples", "s"));
        CONFLICT (flag.qnames_file, flag.qnames_opt,     "--qnames",           "--qnames-file");
        
        // options that require --reference
        if (!flag.reference) {
            ASSINP0 (!flag.show_ref_iupacs, "--show-ref-iupacs requires --reference");
        }
        // options that require absence of --reference 
        else {
            ASSINP (!flag.show_reference, "%s is incompatible with --show-reference. To use --show-reference on a reference file, omit --reference", OT("reference", "e"));
            ASSINP (!flag.show_ref_seq, "%s is incompatible with --show-ref-seq. To use --show-ref-seq on a reference file, omit --reference", OT("reference", "e"));
        }
    }

    // some genozip flags are allowed only in combination with --decompress 
    if (is_genozip && IS_ZIP) {
        char s[20] = {};
        NEED_DECOMPRESS (flag.bgzf != BGZF_NOT_INITIALIZED, OT("bgzf", "z"));
        NEED_DECOMPRESS (!OUT_DT(NONE), str_tolower (dt_name (flag.out_dt), s));
        NEED_DECOMPRESS (flag.unbind, OT("unbind", "u"));
        NEED_DECOMPRESS (flag.show_is_set, "--show_is_set");
        NEED_DECOMPRESS (flag.t_offset, "--t_offset");
        NEED_DECOMPRESS (flag.t_size, "--t_size");
    }

    ASSINP (!IS_REF_EXT_STORE || !is_genocat, "option %s supported only for viewing the reference file itself", OT("REFERENCE", "E"));

    // --output only allowed with a single file, or two FASTQs and --pair, or 3 files and deep
    ASSINP0 (!flag.out_filename || num_files <= 1 || (IS_ZIP && num_files==2 && flag.pair) || (IS_ZIP && num_files>=2 && flag.deep),     
        IS_PIZ ? "--output can only be used with a single input file" // note: flag.out_dt not set yet at this point
               : "--output can only be used with a single input file, or when using --tar (archiving multiple files), --pair (a pair of FASTQs) or --deep (co-compressing related BAM and FASTQs).");

    ASSINP0 (!is_genounzip || !flag.one_vb || flag.test, "--one-vb can only be used with: (1) genocat, or (2) genounzip --test");
}

//--------------------------------
// Data-type specific conditions
//--------------------------------

#define FLAG_ONLY_FOR_2DTs(dt1, dt2, F, Fname)\
    ASSINP0 (!flag.F || (dt == DT_##dt1) || (dt == DT_##dt2) || (DT_##dt1 == DT_SAM && dt == DT_BAM), \
                "--" Fname " is only supported for " #dt1 " and " #dt2 " files");

#define FLAG_ONLY_FOR_DT(dt1, F, Fname)\
    ASSINP0 (!flag.F || (dt == DT_##dt1) || (DT_##dt1 == DT_SAM && dt == DT_BAM), \
                "--" Fname " is only supported for " #dt1 " files");

#define FLAG_NOT_FOR_DT(dt1, F, Fname)\
    ASSINP0 (!flag.F || !((dt == DT_##dt1) || (DT_##dt1 == DT_SAM && dt == DT_BAM)), \
                "--" Fname " is not supported for " #dt1 " files");

void flags_zip_verify_dt_specific (DataType dt)
{
    // SAM
    FLAG_ONLY_FOR_DT(BAM,          show_bam,      "show_bam");
    FLAG_ONLY_FOR_DT(SAM,          optimize_ZM,   "optimize-ZM");
    FLAG_ONLY_FOR_2DTs(SAM, FASTQ, optimize_QUAL, "optimize-QUAL");

    if (flag.show_bam || flag.analyze_ins) 
        flag.seg_only = flag.xthreads = flag.quiet = true; 
    
    // FASTQ
    FLAG_ONLY_FOR_2DTs(SAM, FASTQ, pair, "pair");
    FLAG_ONLY_FOR_DT(FASTQ, optimize_DESC,  "optimize-DESC");
    FLAG_ONLY_FOR_2DTs(FASTQ, FASTA, multiseq, "multiseq");

    // VCF
    FLAG_ONLY_FOR_DT(VCF, add_line_numbers, "add-line-numbers");
    FLAG_ONLY_FOR_DT(VCF, optimize_phred,   "optimize-phred");
    FLAG_ONLY_FOR_DT(VCF, GL_to_PL,         "GL-to-PL");
    FLAG_ONLY_FOR_DT(VCF, GP_to_PP,         "GP-to-PP");
    FLAG_ONLY_FOR_DT(VCF, optimize_VQSLOD,  "optimize-VQSLOD");
    FLAG_ONLY_FOR_2DTs(VCF, GFF, optimize_sort, "optimize-sort");

    // FASTA
    ASSINP0 (!flag.make_reference || dt == DT_REF || dt == DT_GENERIC, 
             "--make-reference is only supported for FASTA files"); // data_type is DT_REF, but error says FASTA

    // GFF
    FLAG_ONLY_FOR_DT(GFF, optimize_Vf,      "optimize-Vf");
    
    // GENERIC
    FLAG_NOT_FOR_DT (GENERIC, debug_lines,  "debug-lines"); // GENERIC doesn't have lines 
}

static void flags_piz_verify_dt_specific (DataType dt)
{
    if (flag.deep && (dt==DT_SAM || dt==DT_BAM) && (flag.sequential || flag.one_component >= 3)) 
        dt = DT_FASTQ;

    // SAM
    FLAG_ONLY_FOR_2DTs(SAM, FASTQ, qname_filter,  "qnames-file");
    FLAG_ONLY_FOR_2DTs(SAM, FASTQ, seq_filter,    "seqs-file");
    FLAG_ONLY_FOR_2DTs(SAM, FASTQ, bases,         "bases");
    FLAG_ONLY_FOR_2DTs(SAM, FASTQ, interleaved,   "interleaved");
    FLAG_ONLY_FOR_2DTs(SAM, FASTQ, show_coverage, "coverage");
    FLAG_ONLY_FOR_2DTs(SAM, FASTQ, idxstats,      "idxstats");
    FLAG_ONLY_FOR_DT(SAM,          FLAG,          "FLAG");
    FLAG_ONLY_FOR_DT(SAM,          MAPQ,          "MAPQ");

    // VCF
    FLAG_ONLY_FOR_DT(VCF,          samples,       "samples");
    FLAG_ONLY_FOR_DT(VCF,          drop_genotypes,"drop-genotypes");
    FLAG_ONLY_FOR_DT(VCF,          gt_only,       "GT-only");
    FLAG_ONLY_FOR_DT(VCF,          snps_only,     "snps-only");
    FLAG_ONLY_FOR_DT(VCF,          indels_only,   "indels-only");

    // FASTQ
    FLAG_ONLY_FOR_DT(FASTQ,          seq_only,    "seq-only");
    FLAG_ONLY_FOR_DT(FASTQ,          qual_only,   "qual-only");

    ASSINP0 (!(flag.deep && flag.header_only), "--header-only is not supported yet for FASTQ components of Deep files (bug 857)");

    // FASTA
    FLAG_ONLY_FOR_DT(FASTA,          sequential,  "sequential");

    // REFERENCE
    FLAG_NOT_FOR_DT (REF,            show_reference, "show_reference");
}

// ZIP: --pair: verify an even number of fastq files, --output, and --reference/--REFERENCE
static void flags_verify_pair_rules (unsigned num_files, rom *filenames)
{
    // verify non-zero even number of files
    ASSINP (num_files && num_files % 2 == 0, "when using %s, expecting an even number of FASTQ input files, each consecutive two being a pair", OT("pair", "2"));
    ASSINP (flag.reference || getenv ("GENOZIP_REFERENCE"), "either --reference or --REFERENCE must be specified when using %s", OT("pair", "2"));

    // verify all are fastq
    for (unsigned i=0; i < num_files; i++)
        ASSINP (txtfile_get_file_dt (filenames[i]) == DT_FASTQ, "when using %s, all input files are expected to be FASTQ files, but %s is not", OT("pair", "2"), filenames[i]);

    // if which --output is missing, we check if every pair of files has a consistent name (except in deep, where default z_name is based on the BAM filename)
    if (!flag.out_filename && !flag.deep) 
        for (unsigned i=0; i < num_files; i += 2) 
            ASSINP (filename_z_pair (filenames[i], filenames[i+1], true),  
                    "to use %s without specifying --output, the naming of the files needs to be consistent and include the numbers 1 and 2 respectively, but these files don't: %s %s", 
                    OT("pair", "2"), filenames[i], filenames[i+1]);
}

// ZIP: --deep: verify conditions
static void flags_zip_verify_deep_rules (unsigned num_files, rom *filenames)
{
    // verify one or two fastq files and one sam/bam/cram file and no stdin (bc we can't confirm its data type yet)
    int n_dt[NUM_DATATYPES] = {};
    int sam_i=-1;
    for (int i=0; i < num_files; i++) {
        rom basename = filename_base (filenames[i], false, 0, 0, 0);
        
        DataType dt = txtfile_get_file_dt (basename);
        if (dt != DT_FASTQ && dt != DT_BAM && dt != DT_SAM) goto validation;
        
        FREE (basename);
        n_dt[dt]++;

        if (dt == DT_SAM || dt == DT_BAM) sam_i = i;
    }

validation:
    ASSINP (n_dt[DT_SAM] + n_dt[DT_BAM] == 1 && 
            n_dt[DT_FASTQ] >= 1 && n_dt[DT_SAM] + n_dt[DT_BAM] + n_dt[DT_FASTQ] == num_files,
            "when using %s, expecting one SAM/BAM/CRAM and one or more FASTQ files", OT("deep", "3"));

    ASSINP (n_dt[DT_FASTQ] <= MAX_NUM_TXT_FILES_IN_ZFILE-1, "Genozip currently supports Deep compression with up to %u FASTQ files, but you specified %u. If you need us to increase this limit, please contact " EMAIL_SUPPORT ".",
            MAX_NUM_TXT_FILES_IN_ZFILE-1, n_dt[DT_FASTQ]);

    ASSINP (flag.reference || getenv ("GENOZIP_REFERENCE"), "either --reference or --REFERENCE must be specified when using %s", OT("deep", "3"));

    // move SAM/BAM to be first
    if (sam_i != 0) SWAP (filenames[sam_i], filenames[0]);

    // case two FASTQs: unless already explicitly set with --pair, implicitly set pair if predicted by filenames    
    if (!flag.pair && n_dt[DT_FASTQ] == 2) { 
        flag.pair = filename_is_fastq_pair (filenames[1], strlen(filenames[1]), filenames[2], strlen(filenames[2]));
        if (flag.pair == PAIR_R2) 
            SWAP (filenames[1], filenames[2]);

        if (flag.pair) flag.pair = PAIRED;
    }

    double bam_size_gb = (double)file_get_size (filenames[0]) / (double)(1 GB);
    double ram_gb = arch_get_physical_mem_size();
    ASSERTW (ram_gb >= bam_size_gb * 2, "FYI: --deep requires a signifcant amount of RAM - usually between 1 and 2 times the size of the BAM file. Your BAM file is %.1f GB while your RAM has only %.1f GB - memory might not be sufficient\n",
             bam_size_gb, ram_gb);

    if (flag.pair)
        flags_verify_pair_rules (num_files - 1, &filenames[1]);
}

static unsigned flags_get_longest_filename (unsigned num_files, rom *filenames)
{
    unsigned len=0;
    for (unsigned i=0; i < num_files; i++) {
        unsigned len_i = strlen (filenames[i]);
        len = MAX_(len, len_i);
    }
    return len;
}

// called from main() after flags are assigned, and before analyzing input files
void flags_update (unsigned num_files, rom *filenames)
{
    // in genocat when showing the file itself, change info_stream to stderr
    if (is_genocat && !(flag.out_filename ||
        // cases when we show some section metadata, and DO NOT show the actual file
        flag.show_dict || flag.show_b250 || flag.show_headers || flag.show_bgzf || flag.dict_id_show_one_b250.num || 
        flag.show_one_dict || flag.show_one_counts.num || flag.show_reference || flag.list_chroms || flag.show_coverage ||
        flag.show_ranges || flag.show_aliases || flag.show_index || flag.show_gheader || flag.show_recon_plan || 
        flag.show_ref_contigs || flag.show_data_type || flag.dump_one_b250_dict_id.num ||
        flag.dump_one_local_dict_id.num || flag.show_singletons_dict_id.num || flag.show_stats || flag.show_ref_index || flag.show_ref_hash || 
        flag.show_chrom2ref || flag.show_ref_seq || flag.show_txt_contigs || flag.count)) {

        info_stream = stderr; // if this is genocat to stdout, send info to stderr instead
        is_info_stream_terminal = isatty (fileno (info_stream)); 
    }

    if (flag.biopsy || flag.biopsy_line.line_i != NO_LINE) {
        ASSERTW (!flag.test,         "FYI: the %s option is ignored when taking a biopsy", OT("test",   "t"));
        ASSERTW (!flag.show_stats,   "FYI: the %s option is ignored when taking a biopsy", OT("stats",  "w"));
        ASSERTW (!flag.out_filename, "FYI: the %s option is ignored when taking a biopsy", OT("output", "o"));
        ASSERTW (!flag.md5,          "FYI: the %s option is ignored when taking a biopsy", OT("md5",    "m"));
        flag.test = false;
        flag.out_filename = NULL;
        flag.md5 = false;
    }
    
    flags_test_conflicts (num_files);

    // verify stuff needed for --pair and --deep
    if (flag.deep && IS_ZIP) flags_zip_verify_deep_rules (num_files, filenames); // --pair and --deep are only available in ZIP
    else if (flag.pair) flags_verify_pair_rules (num_files, filenames); // note: if also deep, called already from flags_zip_verify_deep_rules

    // don't show progress for flags that output throughout the process. no issue with flags that output only in the end
    flag.explicit_quiet = flag.quiet;
    if (flag.show_dict || flag.show_b250 || flag.show_headers || flag.show_threads || flag.show_bgzf || flag.show_mutex || flag.show_containers ||
        flag.dict_id_show_one_b250.num || flag.show_one_dict || flag.show_one_counts.num || flag.show_sag || flag.show_depn || flag.show_singletons_dict_id.num ||
        flag.show_reference || flag.show_digest || flag.list_chroms || flag.show_coverage == COV_ONE || flag.show_ranges || flag.show_snips ||
        flag.show_alleles || flag.show_vblocks || flag.show_codec || flag.show_cache || flag.debug_gencomp || flag.show_qual || flag.show_aligner ||
        flag.show_buddy || flag.debug_peek || flag.show_aliases || (flag.show_index && command==PIZ) || flag.count || flag.biopsy || flag.show_gz)
        flag.quiet=true; // don't show progress or warnings

    // override these ^ if user chose to be --noisy
    if (flag.noisy) flag.quiet=false;

    // note: must be here and not in flags_update_zip_one_file, so its before the z_file creation
    flag.zip_no_z_file = IS_ZIP && 
                         (flag.seg_only || flag.biopsy || flag.biopsy_line.line_i != NO_LINE || flag.show_segconf_has || flag.show_bam || flag.show_gz);

    if (IS_ZIP) {
        // flag.explicit_no_zriter = flag.no_zriter;
        // flag.no_zriter |= global_max_threads==1 || 
        //                   flag.pair || flag.deep; // bug 983
        flag.no_zriter = true;
    }

    flag.multiple_files = (num_files > 1);

    flag.longest_filename = flags_get_longest_filename (num_files, filenames);

    // if using the -o option - check that we don't have duplicate filenames (even in different directory) as they
    // will overwrite each other if extracted genounzip
    if (IS_ZIP && flag.out_filename && !flag.quiet) flags_warn_if_duplicates (num_files, filenames);

    // if reference is not specified, attempt to set it from env var $GENOZIP_REFERENCE
    if (!flag.reference) 
        ref_set_reference (gref, NULL, REF_EXTERNAL, false);

    ASSINP (!flag.make_reference || num_files <= 1, "you can specify only one FASTA file when using --make-reference.\n"
            "To create a reference from multiple FASTAs use something like this:\n"
            "cat *.fa | %s --make-reference --input fasta --output myref.ref.genozip -", global_cmd);

    // where progress, metadata etc messages should go. data always goes to stdout and errors/warning always go to stderr.
    if (flag.log_filename) {
        if (!(info_stream = fopen (flag.log_filename, "a"))) // cannot use ASSERT here bc no info_stream yet...
            fprintf (stderr, "FYI: Failed to open log file %s for writing: %s", flag.log_filename, strerror (errno));
        else
            iprintf ("\n%s\n%s\n", str_time().s, command_line.data);
    }

    if (is_genocat) flag.only_headers = flag.show_headers;

    // cases where we don't need to load the reference file, even if the genozip file normally needs it
    // note: we don't exclude due to collect_coverage here, instead we do it in main_load_reference
    // note: this is here and not in flags_update_piz_one_z_file bc it is consumed by main_load_reference
    // note: this is in flags_update and not flags_update_piz, bc reference file is loaded first 
    flag.genocat_no_ref_file = is_genocat &&
        (flag.show_stats || flag.show_dict || flag.show_b250 || flag.list_chroms || flag.show_one_dict ||
         flag.dump_one_b250_dict_id.num || // all other sections (except CHROM) are blocked from reading in piz_default_skip_section
         flag.show_index || flag.dump_section || flag.show_one_counts.num || flag.show_data_type ||
         flag.show_aliases || flag.show_txt_contigs || flag.show_gheader || flag.show_recon_plan || flag.show_ref_contigs ||
        (flag.count && !flag.bases && !flag.grep) ||
         flag.collect_coverage); // note: this is updated in flags_update_piz_one_z_file

    flag.no_tip |= flag.quiet || (flag.test_i && !flag.noisy) || tar_zip_is_tar() || !isatty(0) || !isatty (1) ||
                   flag.show_bam || flag.biopsy || flag.biopsy_line.line_i != NO_LINE;
}

// ZIP: called ONCE, only for the MAIN component, after opening txt and z files, but before calling zip_one_file 
void flags_update_zip_one_file (void)
{
    DataType dt = z_file->data_type;

    // --make-reference implies --B1 (unless --vblock says otherwise) and not encrypted. 
    // in addition, txtfile_read_vblock() limits each VB to have exactly one contig.
    if (flag.make_reference) {
        ASSINP (!has_password(), "option --make-reference is incompatible with %s", OT("password", "p"));
    }

    z_file->z_flags.txt_is_bin = DTPT (is_binary); // this will go into SectionHeaderGenozipHeader and is determined by the component (eg in Deep it is determined by the SAM/BAM)
    
    DO_ONCE
        if (flag.sendto) WARN0 ("Note: compressing using --sendto means that this file can only be decomrpessed by its intended recipient.\n");

    // if --optimize was selected, all optimizations are turned on
    if (flag.optimize) switch (dt) {
        case DT_BCF   :
        case DT_VCF   : flag.GP_to_PP = flag.GL_to_PL = flag.optimize_phred = flag.optimize_VQSLOD = flag.optimize_sort = true; break;
        case DT_GFF   : flag.optimize_sort = flag.optimize_Vf = true; break;
        case DT_BAM   :
        case DT_SAM   : flag.optimize_QUAL = flag.optimize_ZM = true; break;
        case DT_FASTQ : flag.optimize_QUAL = flag.optimize_DESC = true; break; // note: optimize_DESC is not turned on by --optimize when --deep, bc DT=SAM
        default: break;
    }
    
    // if any optimization flag is on, we turn on flag.optimize
    else flag.optimize = flag.optimize_sort || flag.optimize_phred || flag.GL_to_PL || flag.GP_to_PP || flag.optimize_VQSLOD ||
                         flag.optimize_QUAL || flag.optimize_Vf || flag.optimize_ZM || flag.optimize_DESC;

    // cases where txt data is modified during Seg - digest is not stored, it cannot be tested with --test and other limitations 
    // note: this flag is also set when the file header indicates that it's a Luft file. See vcf_header_get_dual_coords().
    flag.zip_txt_modified = flag.zip_txt_modified  // this is needed, eg, so when compressing a Luft file, the rejects file inherits data_modified 
                      || flag.optimize       // we're modifying data to make it more compressible
                      || flag.match_chrom_to_reference
                      || (flag.add_line_numbers && dt == DT_VCF)
                      || flag.lines_last != NO_LINE  // --head advanced option to compress only a few lines
                      || flag.biopsy_line.line_i != NO_LINE;

    // add v14: --test if user selected --test, OR this condition is true
    flag.explicit_test = flag.test;
    flag.test |= !flag.no_test && !flag.make_reference && !flag.zip_txt_modified && !flag.biopsy && flag.biopsy_line.line_i == NO_LINE &&
                 !flag.show_bam && !flag.seg_only && 
                 !flag.debug && !flag.show_time && !flag.show_memory;

    ASSINP0 (!flag.test_i || flag.test || flag.no_test || flag.debug || flag.make_reference || flag.zip_no_z_file || flag.zip_txt_modified, 
             "When running with GENOZIP_TEST one of: --test, --no-test, --debug, --make-reference must be set");

    ASSINP0 (!flag.match_chrom_to_reference || flag.reference, "--match-chrom-to-reference requires using --reference as well"); 

    flag.bind = flag.deep                         ? BIND_DEEP    // one SAM/BAM (1-3 components) and one or two FASTQs
              : (dt == DT_FASTQ && flag.pair)     ? BIND_FQ_PAIR // FQ_COMP_R1 and FQ_COMP_R2 components
              : (dt == DT_SAM || dt == DT_BAM)    ? BIND_SAM     // SAM_COMP_MAIN component and possibly SAM_COMP_PRIM and/or SAM_COMP_DEPN. If no PRIM/DEPN lines exist, we will cancel it sam_zip_generate_recon_plan
              :                                     BIND_NONE;

    // if biopsy, we seg only for speed. current limitation: in paired files we do the whole thing. TO DO: fix this.
    if (flag.biopsy_line.line_i != NO_LINE && flag.bind != BIND_FQ_PAIR)
        flag.seg_only = true;

    flags_zip_verify_dt_specific (dt);
}

bool flags_is_genocat_global_area_only (void)
{
    return is_genocat &&
        (flag.show_stats || flag.show_dict || flag.list_chroms || flag.show_one_dict ||
         flag.show_index || flag.show_one_counts.num || IS_SHOW_HEADERS ||
         flag.show_reference || flag.show_ref_contigs || flag.show_ranges ||
         flag.show_ref_index || flag.show_ref_hash || flag.show_chrom2ref || 
         flag.show_ref_seq || flag.show_aliases || flag.show_gheader==1 || flag.show_data_type ||
         (Z_DT(REF) && flag.show_ref_iupacs));    
}

static void flags_piz_set_flags_for_deep (void)
{
    if (is_genounzip) {
        flag.deep = true;

        if (!flag.unbind) flag.unbind = "";
    }

    // limitations when genocatting a deep file: user must select between --sam, --bam, --fastq (=interleaved), --interleaved, --R1, --R2 or --one-vb
    // to do: setting the output implicitly with -o to a FASTQ, SAM or BAM file doesn't work yet, this is because
    // flag.explicit_out_dt is set only when opening the TXT file, but we error here before.
    else if (is_genocat && (!flag.no_writer || flag.show_recon_plan)) {
        bool sam   = (flag.explicit_out_dt || flag.one_vb_comp_i <= SAM_COMP_DEPN) && OUT_DT(SAM);
        bool bam   = (flag.explicit_out_dt || flag.one_vb_comp_i <= SAM_COMP_DEPN) && OUT_DT(BAM);
        bool fastq = (flag.explicit_out_dt && OUT_DT(FASTQ)) ||               // --fastq or --output with FASTQ filename
                        flag.one_component >= 1 || // --R1, --R2, --R
                        (flag.one_vb_comp_i != COMP_NONE && flag.one_vb_comp_i >= SAM_COMP_FQ00) || // --one-vb to FASTQ VB
                        flag.interleaved;                                     // --interleaved

        ASSINP0 (sam || bam || fastq || flags_is_genocat_global_area_only(), 
                 "This file was compressed with --deep. Therefore, genocat requires one of these options: --R1, --R2, --R, --interleaved, --fastq, --sam, --bam, --one-vb");

        flag.one_component = (flag.one_vb && 
                              flag.one_vb_comp_i <= SAM_COMP_DEPN) ? (SAM_COMP_MAIN + 1)      // even for a PRIM / DEPN VB
                            : flag.one_vb                          ? (flag.one_vb_comp_i + 1) // +1 bc one_component is 1-based
                            : (fastq && flag.one_component)        ? (flag.one_component + 3) // --R1, --R2, --R --> adjust to SAM_COMP_FQXX
                            : (sam || bam)                         ? (SAM_COMP_MAIN + 1)
                            :                                        0;

        CONFLICT (flag.one_component-1 >= SAM_COMP_FQ00, sam, "--R", "--sam");
        CONFLICT (flag.one_component-1 >= SAM_COMP_FQ00, bam, "--R", "--bam");
        CONFLICT (flag.interleaved, sam, "--interleaved", "--sam"); // not 100% accurate, sam might be set by reasons other than a --sam option. good enough
        CONFLICT (flag.interleaved, bam, "--interleaved", "--bam");
        CONFLICT (flag.one_component-1 >= SAM_COMP_FQ00, flag.one_vb, "--R", "--one-vb");
        CONFLICT (flag.interleaved,                      flag.one_vb, "--interleaved", "--one-vb");
                                
        if (fastq) {
            flag.deep = flag.deep_fq_only = true; // note: in case of SAM/BAM-only reconstruction, we don't need deep.
            flag.out_dt = DT_FASTQ;               // SAM component not written. Used in TRANSLATIONS to specifiy toplevel of TOP2NONE
            
            ASSINP (!flag.interleaved || z_file->num_components == 5, 
                    "--interleaved can't be used because %s does not contain two FASTQs", z_name);
            
            // determined what will be outputted with --fq
            if (!flag.one_component && !flag.one_vb) {
                if (z_file->num_components == 4) // only one FQ file
                    flag.one_component = 4;
                else if (z_file->num_components == 5) { // two FASTQs
                    if (!flag.interleaved) flag.interleaved = INTERLEAVE_BOTH;
                }
                else
                    ABORTINP0 ("--fastq can't be used because file contains more than two FASTQ. Used --R1 or --R2 instead.");
            }

            ASSINP0 (!flag.interleaved || z_file->num_components == 5, // exactly 2 FASTQ components
                     "Cannot interleave FASTQ files, because file does not contain exactly two FASTQ components");

            ASSINP (!flag.interleaved || z_file->comp_num_lines[SAM_COMP_FQ00] == z_file->comp_num_lines[SAM_COMP_FQ01],
                    "Cannot interleave FASTQ files, you must specify --R1 or --R2. This is because R1 has %"PRIu64" lines, but R2 has %"PRIu64,
                    z_file->comp_num_lines[SAM_COMP_FQ00], z_file->comp_num_lines[SAM_COMP_FQ01]);
        }
    }
}

static void flags_piz_set_out_dt (void)
{
    DataType dt_by_filename;

    // case: data type already set by command line option (eg --sam): nothing more to do
    if (!OUT_DT(NONE))
        {}

    // case: set by filename (eg --outfile file.sam) - this is considered explicit setting
    else if (is_genocat && flag.out_filename && 
             ({ dt_by_filename = file_get_data_type (file_get_type (flag.out_filename), true); true; }) &&

             (  (z_file->data_type == dt_by_filename        ) ||
                (Z_DT(SAM)    && dt_by_filename == DT_BAM   ) ||
                (Z_DT(BAM)    && dt_by_filename == DT_SAM   ) ||
                (Z_DT(SAM)    && dt_by_filename == DT_FASTQ ) ||
                (Z_DT(BAM)    && dt_by_filename == DT_FASTQ ) ||
                (Z_DT(VCF)    && dt_by_filename == DT_BCF   ) ||
                (Z_DT(ME23)   && dt_by_filename == DT_VCF   ))) {
        
        flag.out_dt = dt_by_filename;
        flag.explicit_out_dt = true;   // specifiying output data type by filename is considered explicit
    }

    // case: binary file (BAM is recorded in Z_FILE as SAM+txt_is_bin) - if both textual and binary output are possible (SAM/BAM for now):
    // we implicitly translate to textual if genocat without --bgzf is used
    else if (z_file->z_flags.txt_is_bin) {

        // PIZ of a genozip file with is_binary (e.g. BAM) is determined here unless the user overrides with --sam, --bam, --fastq or --bgzf
        if (is_genocat && 
            flag.bgzf == BGZF_NOT_INITIALIZED && // no point outputting as text if BGZF compression was requested in command line 
            DTPZ(bin_type) && DTPZ(txt_type))    // both binary and textual formats exist
            
            flag.out_dt = DTPZ (txt_type); // output in textual format
        else
            flag.out_dt = DTPZ (bin_type); // output in binary format
    }

    // case: data type not set explicitly with command line option (eg --sam or --output) and
    // not chosen for binary file - set implicitly 
    else
        flag.out_dt = z_file->data_type;
}

static bool flags_has_reconstructor_filter (void)
{
    return is_genocat && 
        (((Z_DT(VCF) || Z_DT(BCF)) && (flag.snps_only || flag.indels_only)) || 
         // FASTA specific line droppers
         (Z_DT(FASTA) && (flag.sequential || flag.header_only_fast || flag.header_one || flag.no_header)) || 
         // FASTQ specific line droppers
         (Z_DT(FASTQ) && (flag.bases || flag.qname_filter || flag.qnames_opt || flag.seq_filter)) || 
         // SAM specific line droppers
         (Z_DT(SAM)   && (flag.sam_flag_filter || flag.sam_mapq_filter || flag.bases || flag.qnames_file || flag.qnames_opt || flag.seq_filter || (!flag.deep && OUT_DT(FASTQ)))) || 
         // General filters
         flag.grep || flag.regions);
}

// true if --count happens during writer_create_plan - no reconstruction needed
bool flags_writer_counts (void)
{
    return is_genocat && flag.count==CNT_TOTAL && !flags_has_reconstructor_filter();
}

// called after reading the GENOZIP_HEADER but before reading sections - all we have is the section list
void flags_update_piz_no_ref_file (void)
{
    bool never_need_ref = is_genocat && 
        (flags_writer_counts() || flag.show_stats || flag.show_dict || flag.list_chroms || flag.show_one_dict ||
         flag.show_index || flag.show_data_type ||
         flag.show_aliases || flag.show_gheader || flag.show_recon_plan || flag.show_ref_contigs); 

    // this affects reading an implicit reference specified in the file's SEC_GENOZIP_HEADER. We do it here instead of flags_update because
    // we first need to unset header_only for FASTA and FASTQ
    flag.genocat_no_ref_file |= never_need_ref || (is_genocat && flag.header_only);

    // In SAM with SA:Z: RNAME/POS etc might be reconstructed from SA:Z (saggy); SA:Z requires NM:I ; 
    // NM:I requires MD:Z ; and MD:Z requires SQBITMAP. If SQBITMAP is loaded then SEQ is reconstructed which 
    // requires a reference (see bug 804)
    if (flag.genocat_no_ref_file && !never_need_ref && Z_DT(SAM) && is_there_any_section_with_dict_id (_OPTION_SA_Z))
        flag.genocat_no_ref_file = false;
}

// PIZ: called after opening a z_file and reading the header before opening first txt_file. 
void flags_update_piz_one_z_file (int z_file_i /* -1 if unknown - called form file_open_txt_write */)
{
    if (is_genocat && !flag.out_filename) 
        flag.to_stdout = true; // genocat without -o  

    flags_piz_set_out_dt();

    ASSINP (!flag.test || z_file->z_flags.has_digest, 
            "--test cannot be used with %s, as it was compressed without a digest. See " WEBSITE_DIGEST, z_name);

    flag.collect_coverage = flag.show_coverage || flag.idxstats;
    
    // non-translated mode needed for coverage collection
    if (flag.collect_coverage && Z_DT(SAM)) 
        flag.out_dt = DT_SAM; 
                
    ASSINP0 (!is_genounzip || !flag.to_stdout, "Cannot use --stdout with genounzip, use genocat instead");

    flag.pair = sections_is_paired() ? PAIRED : NOT_PAIRED; // also updates z_file->z_flags in case of backward compatability issues

    ASSINP (!is_genounzip || !flag.out_filename || z_file->num_txt_files == 1, 
            "Cannot use --output because %s contains multiple components. Your options:"
            "\n- Use --prefix to set a prefix for output filenames"
            "\n- Use --output with the target directory name (Tip: end name with / for genozip to mkdir if needed)",
            z_name);
    
    bool is_deep_file = Z_DT(SAM) && VER(15) && z_file->z_flags.dts2_deep; 

    // case: Deep file
    if (is_deep_file) 
        flags_piz_set_flags_for_deep();    

    flags_piz_verify_dt_specific (flag.out_dt); // after deep flags are set

    if (Z_DT(FASTA)) {
        // --downsample in FASTA implies --sequential
        if (flag.downsample)
            flag.sequential = 1;
    }

    else if (Z_DT(FASTQ)) {
        // --R1/--R2 and --interleaved is only possible for on FASTQ data compressed with --pair
        ASSINP (!flag.one_component || flag.pair, 
                "--R%c is not supported for %s because it only works on FASTQ data that was compressed with --pair", '0'+flag.one_component, z_name);

        ASSINP (!flag.interleaved || flag.pair, 
                "--R%c is not supported for %s because it only works on FASTQ data that was compressed with --pair", '0'+flag.one_component, z_name);

        // genocat paired FASTQ: if none of --interleaved, --R1 or --R2 are specified, INTERLEAVE_BOTH is the default
        if (is_genocat && !flag.interleaved && flag.pair && !flag.one_component && !flag.one_vb) 
            flag.interleaved = INTERLEAVE_BOTH; // set before setting flag.maybe_lines_out_of_order

        // genounzip paired FASTQ: always unbind - if user didn't specify prefix, then no prefix
        if (is_genounzip && flag.pair) 
            if (!flag.unbind) flag.unbind = ""; 
    }

    else if (Z_DT(SAM)) {
        ASSINP (!flag.one_component || is_deep_file, "--R%s is supported only for FASTQ files and files compressed with --deep", 
                cond_int (flag.one_component <= 2, "", flag.one_component));

        // --unbind and --component not allowed in SAM/BAM starting v14 (but are allowed starting v15 with --deep)
        if (Z_DT(SAM) && VER(14) && !is_deep_file) {
            ASSERT0 (!flag.unbind, "--unbind not supported for SAM/BAM files");
            ASSERT0 (!flag.one_component, "--component not supported for SAM/BAM files");
        }
        
        // BAM limitations
        if (OUT_DT(BAM) && flag.no_header == 1) {
            ASSINP0 (flag.force, "Using --no-header will result in an invalid BAM file. If this is your intention, use --force.");
            WARN_ONCE0 ("Warning: Outputting a BAM without a header. This will not be a valid BAM file. Use --quiet to suppress this warning.");
        }
    }

    // check validity of --one-vb
    ASSINP (flag.one_vb <= z_file->num_vbs, "%s: --one-vb=%u but file has only %u VBlocks", z_name, flag.one_vb, z_file->num_vbs);

    if (flag.one_vb) 
        flag.one_vb_comp_i = sections_vb_header (flag.one_vb)->comp_i;

    // Check if the reconstructed data type is the same as the source data type
    bool is_binary = z_file->z_flags.txt_is_bin;
    flag.reconstruct_as_src = (OUT_DT(SAM)   && Z_DT(SAM) && !is_binary) || 
                              (OUT_DT(BAM)   && Z_DT(SAM) && is_binary ) ||
                              (OUT_DT(FASTQ) && Z_DT(SAM) && flag.deep ) ||
                              (flag.out_dt == z_file->data_type && !Z_DT(SAM));

    ASSINP (is_genocat || flag.reconstruct_as_src, 
            "genozip file %s is of type %s, but output file is of type %s. Translating between types is not possible in genounzip, use genocat instead",
            z_name, z_dt_name(), dt_name (flag.out_dt));             

    // for FASTA/Q we convert a "header_only" flag to "header_only_fast" because in FASTA/Q
    // this flag affects the reconstruction of the body rather than the txt header (FASTA/Q don't have a txt header)
    if (flag.header_only && (OUT_DT(FASTA) || OUT_DT(FASTQ))) {
        flag.header_only      = false;
        flag.header_only_fast = true;
    }

    // genocat quits after reading global area - it doesn't read TXT_HEADER or VBs
    flag.genocat_global_area_only = flags_is_genocat_global_area_only();

    // if this flag is set, data will be read and uncompressed (unless blocked in piz_is_skip_section), 
    // but not reconstructed or written
    flag.genocat_no_reconstruct = is_genocat && 
        (flag.genocat_global_area_only || flag.dump_section || 
         flag.show_headers || flag.show_txt_contigs || 
         (flag.show_recon_plan && !Z_DT(FASTA) && flag.to_stdout)); // escape to show recon plan with "no_writer": - use -o

    flag.has_reconstructor_filter = flags_has_reconstructor_filter();

    // if this flag is set, no data will be written, although it still could be read and reconstructed 
    // (unless blocked in flag.genocat_no_reconstruct or piz_default_skip_section) 
    flag.no_writer = is_genocat &&
        (flag.genocat_no_reconstruct || flag.collect_coverage || 
         flag.count==CNT_VBs || 
         (flag.count==CNT_TOTAL && flag.has_reconstructor_filter) || // note: if (flag.count==CNT_TOTAL && !flag.has_reconstructor_filter), writer reports during create_plan and no_writer needs to be false
         flag.show_singletons_dict_id.num || flag.dump_one_local_dict_id.num || flag.dump_one_b250_dict_id.num || flag.show_b250 || flag.show_sag);

    flag.no_writer |= is_genounzip && flag.test;

    // usually, if no_writer, we also don't need the writer thread, but there are exceptions to that rule
    flag.no_writer_thread = flag.no_writer && 
                            !(z_sam_gencomp && flag.test); // exception: is SAM with generated components, we need the Writer thread to calculate the digest

    ASSINP (!flag.no_writer || !flag.index_txt, "%s cannot be used because no data is written", OT("index", "x"));

    flag.genocat_no_dicts = is_genocat &&
        (flag.show_stats || flag.show_index || flag.show_one_counts.num || // note: we need dicts for dump_b250 as we need to reconstruct
         flag.show_reference || flag.show_ref_contigs || 
         flag.show_ref_index || flag.show_ref_hash || flag.show_chrom2ref || 
         flag.show_ref_seq || flag.show_aliases || flag.show_txt_contigs || flag.show_gheader || flag.show_data_type ||
         (flag.show_recon_plan && !Z_DT(FASTA))); // we need dicts to generate the FASTA plan (filter for grep etc)

    // don't show progress or warning when outputing to stdout (note: we are "quiet" even if output doesn't go to the terminal
    // because often it will be piped and ultimately go the terminal - a user can override this with --noisy)
    if (flag.to_stdout && !flag.validate && !flag.noisy && !flag.no_writer) flag.quiet=true; 

    // cases we skip displaying the txt header (we still read the section for inspection)
    if (!flag.no_header && is_genocat && 
        flag.no_writer && 
        (flag.genocat_global_area_only || flag.count || flag.show_headers ||
        (z_file_i >= 1 && !flag.no_writer && !flag.header_only))) // when using genocat to concatenate multiple files - don't show the header for the 2nd+ file
        flag.no_header = 2; // 2 = assigned here and not from command line

    bool maybe_txt_header_modified = is_genocat && 
        (flag.no_header || flag.lines_first != NO_LINE || // options that may cause dropping of the txt header
         ((Z_DT(VCF) || Z_DT(BCF)) && (flag.header_one || flag.samples || flag.drop_genotypes))); // VCF specific options that modify the txt header

    bool maybe_vb_dropped_by_writer = is_genocat && // dropped by piz_dispatch_one_vb
        (flag.lines_first != NO_LINE || // decided by writer_create_plan
         flag.tail             || // decided by writer_create_plan
         flag.downsample       || // decided by writer_create_plan
         flag.regions          || // decided by writer_create_plan
         flag.one_vb           || // decided by writer_init_vb_info
         flag.header_only);       // decided by writer_init_vb_info

    bool maybe_vb_dropped_after_read = is_genocat && // dropped by piz_dispatch_one_vb
        ((flag.grep && (Z_DT(FASTQ) || Z_DT(FASTA))) || // decided by piz_read_one_vb
         (flag.regions && Z_DT(FASTA))); // decided by piz_read_one_vb

    flag.maybe_lines_dropped_by_reconstructor = is_genocat && 
        (flag.has_reconstructor_filter ||
         flag.collect_coverage || flag.count); // no-writer, but nevertheless modify the txt_data

    flag.maybe_lines_dropped_by_writer = is_genocat && 
         (flag.downsample || flag.lines_first != NO_LINE || flag.tail);

    flag.maybe_vb_modified_by_reconstructor = is_genocat && 
         // translating to another data
        (!flag.reconstruct_as_src || 
         // lines may be dropped by reconstructor
         flag.maybe_lines_dropped_by_reconstructor || 
         // VCF specific VB modifiers
         ((Z_DT(VCF) || Z_DT(BCF)) && (flag.samples || flag.drop_genotypes || flag.gt_only)) || 
         // FASTA specific modifiers
         (Z_DT(FASTA) && (false)) || 
         // FASTQ specific modifiers
         (Z_DT(FASTQ) && (flag.header_only_fast || flag.seq_only || flag.qual_only)) || // FASTQ "line" is for lines, so these are line modifications, not drops
         // SAM specific modifiers
         (Z_DT(SAM)   && (flag.add_line_numbers)));

    // cases where Writer may re-order lines resulting in different ordering that within the VBs
    flag.maybe_lines_out_of_order = is_genocat && 
         ((Z_DT(FASTQ) && flag.pair && flag.interleaved) || 
          (Z_DT(SAM)   && (z_file->z_flags.has_gencomp || flag.interleaved)));

    // true if the PIZ output txt file will NOT be identical to the source file as recorded in z_file
    flag.piz_txt_modified = maybe_txt_header_modified                  || 
                            maybe_vb_dropped_after_read                ||
                            maybe_vb_dropped_by_writer                 ||
                            flag.maybe_lines_dropped_by_reconstructor  ||
                            flag.maybe_vb_modified_by_reconstructor    || 
                            flag.maybe_lines_out_of_order;

    // calculation depends on flag.piz_txt_modified
    bool pg_line_added_to_header = ((OUT_DT(SAM) || OUT_DT(BAM)) && !flag.reconstruct_as_src)
                                || ((OUT_DT(VCF) || OUT_DT(BCF)) && flag.piz_txt_modified);

    if (pg_line_added_to_header && !flag.no_pg && is_genocat) 
        flag.piz_txt_modified = true;

    ASSINP0 (!is_genounzip || !flag.piz_txt_modified, "Data modification flags are not allowed in genounzip, use genocat instead");

    ASSINP0 (!flag.test || !flag.piz_txt_modified, "--test cannot be used when other flags specify data modification. See " WEBSITE_DIGEST);

    ASSINP0 (!flag.piz_txt_modified || flag.bgzf != BGZF_BY_ZFILE, "Cannot use --bgzf=exact, because other flags might cause output file to differ from original file");

    // cases where we don't read unnecessary contexts, and should just reconstruct them as an empty
    // string (in other cases, it would be an error)
    flag.missing_contexts_allowed = flag.collect_coverage || flag.count || flag.drop_genotypes;

    ASSINP0 (!flag.interleaved || flag.deep_fq_only || flag.pair, 
             "--interleaved is supported only for paired FASTQ files and files compressed with --deep");
    
    ASSINP0 (!flag.unbind || is_deep_file || flag.pair, 
             "--prefix is supported only for paired FASTQ files and files compressed with --deep");

    // downsample not possible for Generic or Chain
    ASSINP (!flag.downsample || (!OUT_DT(CHAIN) && !OUT_DT(GENERIC)), 
            "%s: --downsample is not supported for %s files", z_name, dt_name (flag.out_dt));

    // --sequential only possible on FASTA
    ASSINP (!flag.sequential || OUT_DT(FASTA), 
            "--sequential is not supported for %s because it only works on FASTA data, but this file has %s data",
            z_name, z_dt_name());

    // --coverage is only possible on SAM/BAM and FASTQ
    ASSINP (!flag.show_coverage || OUT_DT(BAM) || OUT_DT(SAM) || OUT_DT(FASTQ), // note: if genozip file has BAM data, it will be translated to SAM bc it is always stdout
            "--coverage is not supported for %s because it only works on SAM, BAM and FASTQ data, but this file has %s data",
            z_name, z_dt_name());

    // --idxstats is only possible on SAM/BAM and FASTQ
    ASSINP (!flag.idxstats || OUT_DT(BAM) || OUT_DT(SAM) || OUT_DT(FASTQ), // note: if genozip file has BAM data, it will be translated to SAM bc it is always stdout
            "--idxstats is not supported for %s because it only works on SAM, BAM and FASTQ data, but this file has %s data",
            z_name, z_dt_name());

    // --add-lines-numbers is only possible on SAM
    ASSINP0 (!flag.add_line_numbers || OUT_DT(SAM), "--add_line_numbers works on SAM/BAM data, when outputting it as SAM");
    
    // // --seq-only and --qual-only only work on FASTQ
    // ASSINP (!flag.seq_only  || OUT_DT(FASTQ), "--seq-only is not supported for %s because it only works on FASTQ data, but this file has %s data", z_name, z_dt_name());
    // ASSINP (!flag.qual_only || OUT_DT(FASTQ), "--qual-only is not supported for %s because it only works on FASTQ data, but this file has %s data", z_name, z_dt_name());

    // --header-one only works on FASTA and VCF
    ASSINP (!flag.header_one || OUT_DT(FASTA) || OUT_DT(VCF) || OUT_DT(BCF), "--header-one is not supported for %s files", z_dt_name());

    // -- grep doesn't work with binary files
    ASSINP (!flag.grep || !out_dt_is_binary, "--grep is not supported when outputting %s data%s", 
            dt_name (flag.out_dt), OUT_DT(BAM) ? ". Tip: add --sam": "");

    // --gpos requires --reference
    ASSINP0 (!flag.gpos || flag.reference, "--gpos requires --reference");
    
    // translator limitations
    ASSINP0 ((!OUT_DT(SAM) && !OUT_DT(BAM)) || Z_DT(SAM) || Z_DT(BAM),
             "--sam and --bam are only allowed for SAM, BAM or CRAM files");

    ASSINP0 (!OUT_DT(VCF) || Z_DT(VCF) || Z_DT(BCF) || Z_DT(ME23),
             "--vcf is only allowed for 23andMe files");

    ASSINP0 ((!OUT_DT(VCF) && !OUT_DT(BCF)) || !Z_DT(ME23) || flag.reference, 
            "--reference must be specified when translating 23andMe to VCF");

    ASSINP0 (!OUT_DT(FASTQ) || Z_DT(FASTQ) || Z_DT(SAM) || Z_DT(BAM),
             "--fastq is only allowed for SAM or BAM files");

    // version limitations

    // sam/bam genozip files generated in v9-11 had a critical bug when translating to fastq
    ASSINP (VER(12) || !(Z_DT(SAM) && OUT_DT(FASTQ)),
            "%s was created with genozip version %u, SAM/BAM to FASTQ translation is supported only for files created with genozip version 12 or later",
            z_name, z_file->genozip_version);

    // num_lines in VbHeader populated since v12 (in v14 moved to SectionEnt)
    ASSINP (VER(12) || (flag.lines_first == NO_LINE && !flag.tail && !flag.downsample),
            "%s was created with genozip version %u, --head, --tail, --lines and --downsample are supported only for files created with genozip version 12 or later",
            z_name, z_file->genozip_version);

    flags_test_conflicts(0); // test again after updating flags
}

static void flags_store_piped_in_details (void)
{
#ifdef __linux__    
    char cmd[200];
    unsigned pid = getpid();
    sprintf (cmd, "lsof -n -P 2> /dev/null | grep $(ls /proc/%u/fd/0 -l | cut -d[ -f2| cut -d] -f1) 2> /dev/null | grep -v %u", pid, pid);
    StreamP strm = stream_create (0, DEFAULT_PIPE_SIZE, 0, 0, 0, 0, 0, "to get counter-process details", 
                                  "bash", "-c", cmd, NULL);
    char str[500];
    unsigned len = fread (str, 1, sizeof (str)-1, stream_from_stream_stdout (strm));
    str[len] = 0;

    char format[20];
    sprintf (format, "%%%us %%u", (unsigned)sizeof (pipe_in_process_name)-1);
    sscanf (str, format, pipe_in_process_name, &pipe_in_pid);
#endif
}

void flags_store_command_line (int argc, char **argv)
{
    unsigned pw_len=0;
    rom pw=0;

    if ((pw = crypt_get_password())) pw_len  = strlen (pw);

    rom slash;
    for (int i=0; i < argc; i++) {

        unsigned arg_len = strlen (argv[i]);

        if (pw && !strcmp (argv[i], pw)) { // "-p 123", "--pass 123" etc
            ASSINP (arg_len < BUFPRINTF_MAX_LEN-10, "argument %u longer than maximum allowed %u characters: %s", i, BUFPRINTF_MAX_LEN-10, argv[i]);
            bufprintf (evb, &command_line, "***%s", (i < argc-1 ? " ": "")); // hide password
        }

        else if (pw && (arg_len >= pw_len + 2) &&  // check for -p123 or eg -fmp123
                 !strcmp (&argv[i][arg_len-pw_len], pw) && // not air-tight test, but good enough (eg "-ofilenamep123" will incorrectly trigger)
                 argv[i][0] == '-' &&
                 argv[i][arg_len-pw_len-1] == 'p') {
            ASSINP (arg_len < BUFPRINTF_MAX_LEN-10, "argument %u longer than maximum allowed %u characters: %s", i, BUFPRINTF_MAX_LEN-10, argv[i]);
            bufprintf (evb, &command_line, "%.*s***%s", arg_len-pw_len, argv[i], (i < argc-1 ? " ": "")); // hide password
        }

        else if (!i && (slash = strrchr (argv[0], '/')))
            bufprintf (evb, &command_line, "%s ", slash+1);

        else if (!i && (slash = strrchr (argv[0], '\\')))
            bufprintf (evb, &command_line, "%s ", slash+1);

        else {
            buf_append_string (evb, &command_line, argv[i]); // can't use bufprintf because argv[i] length is unbound
            if (i < argc-1) BNXTc (command_line) = ' ';   // buf_append_string allocs one char extra
        }

        if (!i) {
            uint32_t i = debugger_params.len32; // IS THIS A BUG? (block-local i override function scope i - loop will not enter)
            bufprintf (evb, &debugger_params, "\"program\": \"%s\",\n", argv[0]);
            for (; i < debugger_params.len32; i++) 
                if (*Bc (debugger_params, i)=='\\')
                    *Bc (debugger_params, i) = '/';
        }
        else {
            bufprint0 (evb, &debugger_params, i==1 ? "\"args\" : [\"" : "\"");
            buf_append_string (evb, &debugger_params, argv[i]); // can't use bufprintf because argv[i] length is unbound
            bufprintf (evb, &debugger_params, "\"%s", i < argc-1 ? ", ": "],");
        }
    }

    if (IS_ZIP && !isatty(0))
        flags_store_piped_in_details();    

    if (flag.echo || flag.debug_or_test) {
        fprintf (stderr, "\n%s: %s\n", str_time().s, command_line.data);
        flag.echo = 2; // done
    }

    if (flag.debug_or_test) 
        fprintf (stderr, "%s\n", debugger_params.data);
}

rom flags_command_line (void)
{
    return command_line.data;
}

rom flags_pipe_in_process_name (void)
{
    return pipe_in_process_name;
}

unsigned flags_pipe_in_pid (void)
{
    return pipe_in_pid;
}

// attempt to check if piped-in process died. 
// note: we won't detect this, if the piped-in process died right away, before completed flags_store_piped_in_details()
bool flags_pipe_in_process_died (void)
{
#ifdef __linux__
    usleep(100000); // wait for dying process to be dead
    return pipe_in_pid && kill (pipe_in_pid, 0);
#else
    return false;
#endif
}

void flags_finalize (void)
{
    buf_destroy (command_line);
    buf_destroy (debugger_params);

    FREE (flag.out_dirname);
}

rom pair_type_name (PairType p)
{
    rom names[] = PAIR_TYPE_NAMES;
    ASSERT (p >= 0 && p <= 3, "Invalid pair_type=%d", p);

    return names[p];
}