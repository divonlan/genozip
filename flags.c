// ------------------------------------------------------------------
//   flags.c
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <getopt.h>
#ifndef _WIN32
#include <sys/types.h>
#include <signal.h>
#include <errno.h>
#endif
#include "genozip.h"
#include "flags.h"
#include "data_types.h"
#include "file.h"
#include "regions.h"
#include "dict_id.h"
#include "crypt.h"
#include "strings.h"
#include "vcf.h"
#include "fastq.h"
#include "stream.h"
#include "bgzf.h"
#include "kraken.h"
#include "bases_filter.h"
#include "website.h"
#include "reference.h"
#include "license.h"
#include "tar.h"
#include "biopsy.h"
#include "endianness.h"
#include "stats.h"

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
    .bgzf         = FLAG_BGZF_BY_ZFILE,
    .kraken_taxid = TAXID_NONE,
    .lines_first  = NO_LINE, 
    .lines_last   = NO_LINE,
    .biopsy_line  = { .line_i = NO_LINE },
    .show_stats_comp_i = COMP_NONE,
};

bool option_is_short[256] = { }; // indexed by character of short option.
FILE *info_stream = 0;      // stdout, stderr or log file - where non-error messages should go
bool is_info_stream_terminal = true;

static Buffer command_line    = EMPTY_BUFFER,
               debugger_params = EMPTY_BUFFER;

static pid_t pipe_in_pid = 0;
static char pipe_in_process_name[100] = "";

// command line options that get assigned to other flags (they are not flags themselves)
static int option_noisy=0;

static void flags_show_flags (void)
{
    #define TF_(f) (flag.f ? "true" : "false") 
    #define S_(f)  (flag.f ? flag.f : "(none)")

    iprintf ("fast=%s\n", TF_(fast));
    iprintf ("best=%s\n", TF_(best));
    iprintf ("make_reference=%s\n", TF_(make_reference));
    iprintf ("multiseq=%s\n", TF_(multiseq));
    iprintf ("md5=%s\n", TF_(md5));
    iprintf ("vblock=%s\n", S_(vblock));
    iprintf ("optimize=%s\n", TF_(optimize));
    iprintf ("optimize_sort=%s\n", TF_(optimize_sort));
    iprintf ("optimize_phred=%s\n", TF_(optimize_phred));
    iprintf ("GL_to_PL=%s\n", TF_(GL_to_PL));
    iprintf ("GP_to_PP=%s\n", TF_(GP_to_PP));
    iprintf ("optimize_VQSLOD=%s\n", TF_(optimize_VQSLOD));
    iprintf ("optimize_QUAL=%s\n", TF_(optimize_QUAL));
    iprintf ("optimize_Vf=%s\n", TF_(optimize_Vf));
    iprintf ("optimize_ZM=%s\n", TF_(optimize_ZM));
    iprintf ("pair=%d\n", flag.pair);
    iprintf ("undocumented_dts_paired=%s\n", TF_(undocumented_dts_paired));
    iprintf ("bgzf=%d\n", flag.bgzf);
    iprintf ("out_dt=%s\n", dt_name (flag.out_dt));
    iprintf ("explicit_out_dt=%s\n", TF_(explicit_out_dt));
    iprintf ("header_one=%s\n", TF_(header_one));
    iprintf ("no_header=%d\n", flag.no_header);
    iprintf ("header_only=%s\n", TF_(header_only));
    iprintf ("header_only_fast=%s\n", TF_(header_only_fast));
    iprintf ("seq_only=%s\n", TF_(seq_only));
    iprintf ("qual_only=%s\n", TF_(qual_only));
    iprintf ("regions=%s\n", TF_(regions));
    iprintf ("gpos=%s\n", TF_(gpos));
    iprintf ("samples=%s\n", TF_(samples));
    iprintf ("drop_genotypes=%s\n", TF_(drop_genotypes));
    iprintf ("gt_only=%s\n", TF_(gt_only));
    iprintf ("snps_only=%s\n", TF_(snps_only));
    iprintf ("indels_only=%s\n", TF_(indels_only));
    iprintf ("sequential=%s\n", TF_(sequential));
    iprintf ("no_pg=%s\n", TF_(no_pg));
    iprintf ("extended_translation=%s\n", TF_(extended_translation));
    iprintf ("interleaved=%s\n", flag.interleaved==INTERLEAVE_NONE ? "none" : flag.interleaved==INTERLEAVE_BOTH ? "both" : flag.interleaved==INTERLEAVE_EITHER ? "either" : "invalid value");
    iprintf ("luft=%s\n", TF_(luft));
    iprintf ("sort=%s\n", TF_(sort));
    iprintf ("unsorted=%s\n", TF_(unsorted));
    iprintf ("kraken_taxid=%d\n", flag.kraken_taxid);
    iprintf ("grep=%s grepw=%s grep_len=%u\n", S_(grep), TF_(grepw), flag.grep_len);
    iprintf ("lines_first=%"PRId64"\n", flag.lines_first);
    iprintf ("lines_last=%"PRId64"\n", flag.lines_last);
    iprintf ("tail=%"PRId64"\n", flag.tail);
    iprintf ("one_vb=%d\n", flag.one_vb);
    iprintf ("one_component=%d\n", flag.one_component);
    iprintf ("downsample=%d\n", flag.downsample);
    iprintf ("shard=%d\n", flag.shard);
    iprintf ("sam_flag_filter=%d\n", flag.sam_flag_filter);
    iprintf ("sam_mapq_filter=%d\n", flag.sam_mapq_filter);
    iprintf ("FLAG=%d\n", flag.FLAG);
    iprintf ("MAPQ=%d\n", flag.MAPQ);
    iupac_show();    
    iprintf ("bytes=%d\n", flag.bytes);
    iprintf ("lic_width=%d\n", flag.lic_width);
    iprintf ("force=%s\n", TF_(force));
    iprintf ("quiet=%s\n", TF_(quiet));
    iprintf ("to_stdout=%s\n", TF_(to_stdout));
    iprintf ("replace=%s\n", TF_(replace));
    iprintf ("do_register=%s\n", TF_(do_register));
    iprintf ("test=%s\n", TF_(test));
    iprintf ("index_txt=%s\n", TF_(index_txt));
    iprintf ("subdirs=%s\n", TF_(subdirs));
    iprintf ("list=%s\n", TF_(list));
    iprintf ("threads_str=%s\n", S_(threads_str));
    iprintf ("out_filename=%s\n", S_(out_filename));
    iprintf ("out_dirname=%s\n", S_(out_dirname));
    iprintf ("reference=%d\n", flag.reference);
    iprintf ("show_stats=%d\n", flag.show_stats);
    iprintf ("show_lift=%d\n", flag.show_lift);
    iprintf ("show_lines=%d\n", flag.show_lines);
    iprintf ("validate=%d\n", flag.validate);
    iprintf ("list_chroms=%s\n", TF_(list_chroms));
    iprintf ("show_sex=%s\n", TF_(show_sex));
    iprintf ("idxstats=%s\n", TF_(idxstats));
    iprintf ("count=%d\n", flag.count);
    iprintf ("show_coverage=%d\n", flag.show_coverage);
    iprintf ("show_dvcf=%s\n", TF_(show_dvcf));
    iprintf ("show_ostatus=%s\n", TF_(show_ostatus));
    iprintf ("show_memory=%s\n", TF_(show_memory));
    iprintf ("show_dict=%s\n", TF_(show_dict));
    iprintf ("show_sa=%d\n", flag.show_sa);
    iprintf ("show_depn=%s\n", TF_(show_depn));
    iprintf ("show_b250=%s\n", TF_(show_b250));
    iprintf ("show_aliases=%s\n", TF_(show_aliases));
    iprintf ("show_digest=%s\n", TF_(show_digest));
    iprintf ("log_digest=%s\n", TF_(log_digest));
    iprintf ("show_recon_plan=%s\n", TF_(show_recon_plan));
    iprintf ("show_index=%s\n", TF_(show_index));
    iprintf ("show_gheader=%s\n", TF_(show_gheader));
    iprintf ("show_ref_contigs=%s\n", TF_(show_ref_contigs));
    iprintf ("show_ref_diff=%s\n", TF_(show_ref_diff));
    iprintf ("show_chain_contigs=%s\n", TF_(show_chain_contigs));
    iprintf ("show_ref_seq=%s\n", TF_(show_ref_seq));
    iprintf ("show_reference=%s\n", TF_(show_reference));
    iprintf ("show_ref_hash=%s\n", TF_(show_ref_hash));
    iprintf ("show_ref_index=%s\n", TF_(show_ref_index));
    iprintf ("show_chrom2ref=%s\n", TF_(show_chrom2ref));
    iprintf ("show_ref_iupacs=%s\n", TF_(show_ref_iupacs));
    iprintf ("show_chain=%s\n", TF_(show_chain));
    iprintf ("show_codec=%s\n", TF_(show_codec));
    iprintf ("show_aligner=%s\n", TF_(show_aligner));
    iprintf ("show_containers=%d\n", flag.show_containers);
    iprintf ("dict_id_show_containers=%s\n", dis_dict_id (flag.dict_id_show_containers).s);
    iprintf ("show_alleles=%s\n", TF_(show_alleles));
    iprintf ("show_bgzf=%s\n", TF_(show_bgzf));
    iprintf ("show_txt_contigs=%s\n", TF_(show_txt_contigs));
    iprintf ("show_vblocks=%s\n", TF_(show_vblocks));
    iprintf ("show_threads=%s\n", TF_(show_threads));
    iprintf ("show_wrong_md=%s\n", TF_(show_wrong_md));
    iprintf ("show_wrong_xg=%s\n", TF_(show_wrong_xg));
    iprintf ("show_wrong_xm=%s\n", TF_(show_wrong_xm));
    iprintf ("show_kraken=%d\n", flag.show_kraken);
    iprintf ("show_uncompress=%s\n", TF_(show_uncompress));
    iprintf ("debug_progress=%s\n", TF_(debug_progress));
    iprintf ("debug_LONG=%s\n", TF_(debug_LONG));
    iprintf ("show_qual=%s\n", TF_(show_qual));
    iprintf ("debug_qname=%s\n", TF_(debug_qname));
    iprintf ("show_hash=%s\n", TF_(show_hash));
    iprintf ("debug_memory=%s\n", TF_(debug_memory));
    iprintf ("debug_threads=%s\n", TF_(debug_threads));
    iprintf ("debug_stats=%s\n", TF_(debug_stats));
    iprintf ("debug_generate=%s\n", TF_(debug_generate)); 
    iprintf ("debug_recon_size=%s\n", TF_(debug_recon_size)); 
    iprintf ("debug_seg=%s\n", TF_(debug_seg)); 
    iprintf ("debug_read_ctxs=%s\n", TF_(debug_read_ctxs));
    iprintf ("debug_sa=%s\n", TF_(debug_sa));
    iprintf ("debug_gencomp=%s\n", TF_(debug_gencomp));
    iprintf ("check_lastest=%s\n", TF_(check_latest));
    iprintf ("debug_lastest=%s\n", TF_(debug_latest));
    iprintf ("no_gencomp=%s\n", TF_(no_gencomp));
    iprintf ("no_domqual=%s\n", TF_(no_domqual));
    iprintf ("debug_lines=%s\n", TF_(debug_lines));
    iprintf ("verify_codec=%s\n", TF_(verify_codec)); 
    iprintf ("seg_only=%s\n", TF_(seg_only));
    iprintf ("show_bam=%s\n", TF_(show_bam));
    iprintf ("biopsy=%s\n", TF_(biopsy));
    iprintf ("biopsy_line=%d,%u\n", flag.biopsy_line.vb_i, flag.biopsy_line.line_i);
    iprintf ("xthreads=%s\n", TF_(xthreads));
    iprintf ("show_flags=%s\n", TF_(show_flags));
    iprintf ("echo=%s\n", TF_(echo));
    iprintf ("show_headers=%d\n", flag.show_headers);
    iprintf ("help=%s\n", S_(help));
    iprintf ("dump_section=%s\n", flag.dump_section ? flag.dump_section : "(dump_section)");
    iprintf ("show_is_set=%s\n", S_(show_is_set));
    iprintf ("show_time=%s\n", S_(show_time));
    iprintf ("show_mutex=%s\n", S_(show_mutex));
    iprintf ("dict_id_debug_seg=%s\n", dis_dict_id (flag.dict_id_debug_seg).s);
    iprintf ("dict_id_show_one_b250=%s\n", dis_dict_id (flag.dict_id_show_one_b250).s);
    iprintf ("show_one_counts=%s\n", dis_dict_id (flag.show_one_counts).s);
    iprintf ("dump_one_b250_dict_id=%s\n", dis_dict_id (flag.dump_one_b250_dict_id).s);
    iprintf ("dump_one_local_dict_id=%s\n", dis_dict_id (flag.dump_one_local_dict_id).s);
    iprintf ("show_one_dict=%s\n", S_(show_one_dict));
    iprintf ("debug=%s\n", TF_(debug));
    iprintf ("debug_top=%s\n", TF_(debug_top));
    iprintf ("windows=%s\n", TF_(is_windows));
    iprintf ("apple=%s\n", TF_(is_mac));
    iprintf ("linux=%s\n", TF_(is_linux));
    iprintf ("aligner_available=%s\n", TF_(aligner_available));
    iprintf ("reading_reference=%s\n", flag.reading_reference==gref ? "GREF" : flag.reading_reference==prim_ref ? "CHAIN_SRC" : "NONE");
    iprintf ("zip_comp_i=%s\n", comp_name(flag.zip_comp_i));
    iprintf ("genocat_no_ref_file=%s\n", TF_(genocat_no_ref_file));
    iprintf ("genocat_no_dicts=%s\n", TF_(genocat_no_dicts));
    iprintf ("genocat_global_area_only=%s\n", TF_(genocat_global_area_only));
    iprintf ("genocat_no_reconstruct=%s\n", TF_(genocat_no_reconstruct));
    iprintf ("no_writer=%s\n", TF_(no_writer));
    iprintf ("no_writer_thread=%s\n", TF_(no_writer_thread));
    iprintf ("multiple_files=%s\n", TF_(multiple_files));
    iprintf ("reconstruct_as_src=%s\n", TF_(reconstruct_as_src));
    iprintf ("maybe_txt_header_modified=%s\n", TF_(maybe_txt_header_modified));
    iprintf ("maybe_vb_dropped_by_writer=%s\n", TF_(maybe_vb_dropped_by_writer));
    iprintf ("maybe_vb_dropped_after_read=%s\n", TF_(maybe_vb_dropped_after_read));
    iprintf ("missing_contexts_allowed=%s\n", TF_(missing_contexts_allowed));
    iprintf ("maybe_vb_modified_by_reconstructor=%s\n", TF_(maybe_vb_modified_by_reconstructor));
    iprintf ("maybe_lines_dropped_by_reconstructor=%s\n", TF_(maybe_lines_dropped_by_reconstructor));
    iprintf ("maybe_lines_dropped_by_writer=%s\n", TF_(maybe_lines_dropped_by_writer));
    iprintf ("maybe_lines_out_of_order=%s\n", TF_(maybe_lines_out_of_order));
    iprintf ("data_modified=%s\n", TF_(data_modified));
    iprintf ("explicit_ref=%s\n", TF_(explicit_ref));
    iprintf ("collect_coverage=%s\n", TF_(collect_coverage));
    iprintf ("dvcf_rename=%s\n", S_(dvcf_rename));
    iprintf ("dvcf_drop=%s\n", S_(dvcf_drop));
    iprintf ("single_coord=%s\n", TF_(single_coord));
    iprintf ("reading_chain=%s\n", S_(reading_chain));
    iprintf ("reading_kraken=%s\n", S_(reading_kraken));
    iprintf ("reading_chain=%s\n", S_(reading_chain));
    iprintf ("preprocessing=%s\n", TF_(preprocessing));
    iprintf ("unbind=%s\n", S_(unbind));
    iprintf ("log_filename=%s\n", S_(log_filename));
    iprintf ("bind=%d\n", flag.bind);
    iprintf ("stdin_size=%"PRIu64"\n", flag.stdin_size);
    iprintf ("longest_filename=%d\n", flag.longest_filename);
    iprintf ("match-chrom-to-reference=%s", TF_(match_chrom_to_reference));

    #undef TF
    #undef S
}

#define MAX_LINE ((int64_t)1 << 62)
static void flag_set_lines (rom optarg)
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
static void flag_set_out_filename (char *optarg)
{
    int len = strlen (optarg);

    // case: a final / means this is a directory (possibly to be created)
    if (optarg[len-1] == '/') {
        flag.out_dirname = CALLOC (len);
        memcpy ((char*)flag.out_dirname, optarg, len-1);
        file_mkdir (flag.out_dirname);
    }

    // case: a existing directory name
    else if (file_is_dir (optarg)) 
        flag.out_dirname = optarg;
    
    // case: not a directory - treat as a filename
    else
        flag.out_filename = optarg;
}

static void flag_set_head (rom optarg)
{
    if (!optarg) 
        flag_set_lines ("1-10"); // default

    else {        
        ASSINP0 (str_get_int_range64 (optarg, 0, 1, MAX_LINE, NULL),
                "The optional argument of --head must be a positive integer value");

        char s[strlen(optarg) + 4];
        sprintf (s, "1-%s", optarg);
        flag_set_lines (s);
    }
}

static void flag_set_tail (rom optarg)
{
    ASSINP0 (flag.lines_first==NO_LINE && !flag.tail, "Only one of these can be used: --lines, --head, --tail");

    if (!optarg) 
        flag.tail = 10; // default

    else 
        ASSINP0 (str_get_int_range64 (optarg, 0, 1, MAX_LINE, &flag.tail),
                "The optional argument of --tail must be a positive integer value");
}

static void flags_set_bgzf (rom level_str)
{
    int64_t level_64=BGZF_COMP_LEVEL_DEFAULT;

    ASSINP0 (level_str && str_get_int (level_str, strlen (level_str), &level_64),
             "--bgzf expects a value between 0 (no compression) and 12 (best, yet slowest, compression). If you're not sure what value to choose, 6 is a popular option.");

    flag.bgzf = (int)level_64;
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

static void flag_set_stats (StatsType stats_type, rom optarg)
{
    flag.show_stats = stats_type;

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
        #define _I  {"stdin-size",       required_argument, 0, 'I'                    }
        #define _d  {"decompress",       no_argument,       &command, PIZ             }
        #define _f  {"force",            no_argument,       &flag.force,            1 }
        #define _h  {"help",             optional_argument, 0, 'h'                    }
        #define _l  {"list",             no_argument,       &flag.list,             1 }
        #define _L1 {"license",          optional_argument, 0, 'L'                    } // US spelling
        #define _L2 {"licence",          optional_argument, 0, 'L'                    } // British spelling
        #define _Lf {"licfile",          required_argument, 0, 26                     }
        #define _q  {"quiet",            no_argument,       &flag.quiet,            1 }
        #define _Q  {"noisy",            no_argument,       &option_noisy,          1 }
        #define _so {"sort",             no_argument,       &flag.sort,             1 }
        #define _SO {"unsorted",         no_argument,       &flag.unsorted,         1 }
        #define _DL {"replace",          no_argument,       &flag.replace,          1 }
        #define _V  {"version",          no_argument,       &command, VERSION         }
        #define _z  {"bgzf",             required_argument, 0, 'z'                    }
        #define _zb {"bam",              no_argument,       &flag.out_dt,      DT_BAM }
        #define _zB {"BAM",              no_argument,       &flag.out_dt,      DT_BAM }
        #define _zs {"sam",              no_argument,       &flag.out_dt,      DT_SAM }
        #define _zS {"SAM",              no_argument,       &flag.out_dt,      DT_SAM }
        #define _zq {"fastq",            optional_argument, 0, 130                    }
        #define _zQ {"FASTQ",            optional_argument, 0, 130                    }
        #define _zf {"fq",               optional_argument, 0, 130                    }
        #define _zF {"FQ",               optional_argument, 0, 130                    }
        #define _za {"fasta",            no_argument,       &flag.out_dt,    DT_FASTA }
        #define _zA {"FASTA",            no_argument,       &flag.out_dt,    DT_FASTA }
        #define _zc {"bcf",              no_argument,       &flag.out_dt,      DT_BCF }
        #define _zC {"BCF",              no_argument,       &flag.out_dt,      DT_BCF }
        #define _zv {"vcf",              no_argument,       &flag.out_dt,      DT_VCF }
        #define _zV {"VCF",              no_argument,       &flag.out_dt,      DT_VCF }
        #define _zy {"phylip",           no_argument,       &flag.out_dt,   DT_PHYLIP }
        #define _zY {"Phylip",           no_argument,       &flag.out_dt,   DT_PHYLIP }
        #define _m  {"md5",              no_argument,       &flag.md5,              1 }
        #define _t  {"test",             no_argument,       &flag.test,             1 }
        #define _fa {"fast",             no_argument,       &flag.fast,             1 }
        #define _bs {"best",             optional_argument, 0, 'b'                    }
        #define _9  {"optimize",         no_argument,       &flag.optimize,         1 } // US spelling
        #define _99 {"optimise",         no_argument,       &flag.optimize,         1 } // British spelling
        #define _9s {"optimize-sort",    no_argument,       &flag.optimize_sort,    1 }
        #define _9P {"optimize-phred",   no_argument,       &flag.optimize_phred,   1 }
        #define _9G {"GL-to-PL",         no_argument,       &flag.GL_to_PL,         1 }
        #define _9g {"GP-to-PP",         no_argument,       &flag.GP_to_PP,         1 }
        #define _9V {"optimize-VQSLOD",  no_argument,       &flag.optimize_VQSLOD,  1 }
        #define _9Q {"optimize-QUAL",    no_argument,       &flag.optimize_QUAL,    1 } 
        #define _9f {"optimize-Vf",      no_argument,       &flag.optimize_Vf,      1 }
        #define _9Z {"optimize-ZM",      no_argument,       &flag.optimize_ZM,      1 }
        #define _9D {"optimize-DESC",    no_argument,       &flag.optimize_DESC,    1 }
        #define _al {"add-line-numbers", no_argument,       &flag.add_line_numbers, 1 }
        #define _pe {"pair",             no_argument,       &flag.pair,   PAIR_READ_1 } 
        #define _pt {"dts_paired",       no_argument,       &flag.undocumented_dts_paired, 1 }  // undocumented flag to uncompress paired files older than 9.0.13 when genozip_header.dts_paired was introduced. A user will get an error message instructing her to use it.
        #define _th {"threads",          required_argument, 0, '@'                    }
        #define _u  {"prefix",           required_argument, 0, 'u'                    }
        #define _o  {"output",           required_argument, 0, 'o'                    }
        #define _T  {"files-from",       required_argument, 0, 'T'                    }
        #define _TT {"tar",              required_argument, 0, 27                     }
        #define _p  {"password",         required_argument, 0, 'p'                    }
        #define _B  {"vblock",           required_argument, 0, 'B'                    }
        #define _r  {"regions",          required_argument, 0, 'r'                    }
        #define _R  {"regions-file",     required_argument, 0, 'R'                    }
        #define _Rg {"gpos",             no_argument,       &flag.gpos,             1 }
        #define _s  {"samples",          required_argument, 0, 's'                    }
        #define _sf {"FLAG",             required_argument, 0, 17                     }
        #define _sq {"MAPQ",             required_argument, 0, 18                     }
        #define _il {"interleaved",      optional_argument, 0, 29                     } // both --interleave and --interleaved will be accepted
        #define _e  {"reference",        required_argument, 0, 'e'                    }
        #define _E  {"REFERENCE",        required_argument, 0, 'E'                    }
        #define _lo {"luft",             no_argument,       &flag.luft,             1 }
        #define _ch {"chain",            required_argument, 0, 'C'                    }
        #define _kr {"kraken",           required_argument, 0, 'K'                    }
        #define _kR {"taxid",            required_argument, 0, 'k'                    }
        #define _b  {"bytes",            no_argument,       &flag.bytes,            1 }
        #define _me {"make-reference",   no_argument,       &flag.make_reference,   1 }
        #define _mf {"multiseq",         no_argument,       &flag.multiseq,         1 }
        #define _mF {"multi-seq",        no_argument,       &flag.multiseq,         1 }
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
        #define _sl {"show-lifts",       no_argument,       &flag.show_lift,        1 } 
        #define _sL {"show-lines",       no_argument,       &flag.show_lines,       1 } 
        #define _ss {"stats",            optional_argument, 0, 138,                   } 
        #define _SS {"STATS",            optional_argument, 0, 139                    } 
        #define _lc {"list-chroms",      no_argument,       &flag.list_chroms,      1 } // identical to --show-dict=CHROM 
        #define _lh {"chroms",           no_argument,       &flag.list_chroms,      1 } // identical to --show-dict=CHROM 
        #define _lH {"contigs",          no_argument,       &flag.list_chroms,      1 } 
        #define _s2 {"show-b250",        optional_argument, 0, 2,                     }
        #define _sd {"show-dict",        optional_argument, 0, 3                      }
        #define _s7 {"dump-b250",        required_argument, 0, 5                      }
        #define _S7 {"dump-local",       required_argument, 0, 6                      }
        #define _S8 {"show-counts",      required_argument, 0, 16                     }
        #define _S9 {"dump-section",     required_argument, 0, 7                      }        
        #define _sa {"show-alleles",     no_argument,       &flag.show_alleles,     1 }
        #define _st {"show-time",        optional_argument, 0, 1                      } 
        #define _sm {"show-memory",      no_argument,       &flag.show_memory,      1 } 
        #define _sS {"show-dvcf",        no_argument,       &flag.show_dvcf,        1 } 
        #define _SC {"single-coord",     no_argument,       &flag.single_coord,     1 } 
        #define _XS {"show-rename-tags", no_argument,       &flag.show_rename_tags, 1 } 
        #define _sY {"show-ostatus",     no_argument,       &flag.show_ostatus,     1 } 
        #define _sy {"show-ostatus-counts", no_argument,    0, 131                    } 
        #define _sh {"show-headers",     optional_argument, 0, 10                     } 
        #define _si {"show-index",       no_argument,       &flag.show_index,       1 } 
        #define _sG {"show-ranges",      no_argument,       &flag.show_ranges,      1 } 
        #define _Si {"show-ref-index",   no_argument,       &flag.show_ref_index,   1 } 
        #define _Sh {"show-ref-hash" ,   no_argument,       &flag.show_ref_hash,    1 } 
        #define _sr {"show-gheader",     optional_argument, 0, 135                    }  
        #define _sT {"show-threads",     no_argument,       &flag.show_threads,     1 }  
        #define _wM {"show-wrong-MD",    no_argument,       &flag.show_wrong_md,    1 }  
        #define _wm {"show-wrong-XG",    no_argument,       &flag.show_wrong_xg,    1 }  
        #define _WM {"show-wrong-XM",    no_argument,       &flag.show_wrong_xm,    1 }  
        #define _sF {"show-flags",       no_argument,       &flag.show_flags,       1 }  
        #define _su {"show-uncompress",  no_argument,       &flag.show_uncompress,  1 }  
        #define _sK {"show-kraken",      optional_argument, 0, 21                     }  
        #define _sv {"show-vblocks",     no_argument,       &flag.show_vblocks,     1 }  
        #define _sH {"show-chain",       no_argument,       &flag.show_chain,       1 }  
        #define _sn {"show-filename",    no_argument,       &flag.show_filename,    1 }  
        #define _ov {"one-vb",           required_argument, 0, 8                      }  
        #define _R1 {"R1",               no_argument,       &flag.one_component,    1 }  // comp_i+1
        #define _R2 {"R2",               no_argument,       &flag.one_component,    2 }  
        #define _sR {"show-reference",   no_argument,       &flag.show_reference,   1 }  
        #define _sC {"show-ref-contigs", no_argument,       &flag.show_ref_contigs, 1 }  
        #define _sD {"show-ref-diff",    no_argument,       &flag.show_ref_diff,    1 }  
        #define _cC {"show-chain-contigs", no_argument,     &flag.show_chain_contigs,1}  
        #define _rI {"show-ref-iupacs",  no_argument,       &flag.show_ref_iupacs,  1 }  
        #define _rA {"show-chrom2ref",   no_argument,       &flag.show_chrom2ref,   1 }  
        #define _rS {"show-ref-seq",     no_argument,       &flag.show_ref_seq,     1 }  
        #define _cn {"show-containers",  optional_argument, 0, 132                    }  
        #define _dg {"debug-seg",        optional_argument, 0, 133                    }  
        #define _hC {"show-txt-contigs", no_argument,       &flag.show_txt_contigs, 1 }
        #define _sI {"show-is-set",      required_argument, 0, '~',                   }  
        #define _sA {"show-aliases",     no_argument,       &flag.show_aliases,     1 }  
        #define _sB {"show-sa",          optional_argument, 0, 136                    }  
        #define _sP {"show-depn",        no_argument,       &flag.show_depn,        1 }  
        #define _sc {"show-codec",       no_argument,       &flag.show_codec,       1 }  
        #define _AL {"show-aligner",     no_argument,       &flag.show_aligner,     1 }  
        #define _sb {"show-bgzf",        no_argument,       &flag.show_bgzf,        1 }
        #define _s5 {"show-digest",      no_argument,       &flag.show_digest,      1 }
        #define _S5 {"log-digest",       no_argument,       &flag.log_digest,       1 }
        #define _s6 {"show-plan",        no_argument,       &flag.show_recon_plan,  1 }
        #define _sM {"show-mutex",       optional_argument, 0, 4                      }
        #define _dS {"seg-only",         no_argument,       &flag.seg_only,         1 }  
        #define _bS {"show-bam",         no_argument,       &flag.show_bam,         1 }  
        #define _xt {"xthreads",         no_argument,       &flag.xthreads,         1 }  
        #define _dm {"debug-memory",     optional_argument, 0, 12                     }  
        #define _dp {"debug-progress",   no_argument,       &flag.debug_progress,   1 }  
        #define _dL {"debug-LONG",       no_argument,       &flag.debug_LONG,       1 }  
        #define _dD {"show-qual",        no_argument,       &flag.show_qual,        1 }  
        #define _dq {"debug-qname",      no_argument,       &flag.debug_qname,      1 }  
        #define _dt {"debug-threads",    no_argument,       &flag.debug_threads,    1 }  
        #define _dw {"debug-stats",      no_argument,       &flag.debug_stats,      1 }  
        #define _dM {"debug-generate",   no_argument,       &flag.debug_generate  , 1 }  
        #define _dr {"debug-recon-size", no_argument,       &flag.debug_recon_size, 1 }  
        #define _dR {"debug-read-ctxs",  no_argument,       &flag.debug_read_ctxs,  1 }  
        #define _dP {"debug-sa",         no_argument,       &flag.debug_sa,         1 }  
        #define _dG {"debug-gencomp",    no_argument,       &flag.debug_gencomp,    1 }  
        #define _dN {"no-gencomp",       no_argument,       &flag.no_gencomp,       1 }  
        #define _dQ {"no-domqual",       no_argument,       &flag.no_domqual,       1 }  
        #define _dl {"debug-lines",      no_argument,       &flag.debug_lines,      1 }  
        #define _dc {"verify-codec",     no_argument,       &flag.verify_codec,     1 }          
        #define _oe {"echo",             no_argument,       &flag.echo,             1 }
        #define _dh {"show-hash",        no_argument,       &flag.show_hash,        1 }  
        #define _sx {"sex",              no_argument,       &flag.show_sex,         1 }  
        #define _ct {"count",            optional_argument, 0, 20                     }  
        #define _SX {"coverage",         optional_argument, 0, 13                     }  
        #define _ix {"idxstats",         no_argument,       &flag.idxstats,         1 }
        #define _vl {"validate",         optional_argument, 0, 19                     }  
        #define _lg {"log",              required_argument, 0, 15                     }  
        #define _Xr {"dvcf-rename",      required_argument, 0, 128,                   }
        #define _Xd {"dvcf-drop",        required_argument, 0, 129,                   }
        #define _bi {"biopsy",           required_argument, 0, 134,                   }
        #define _bl {"biopsy-line",      required_argument, 0, 137,                   }
        #define _MR {"match-chrom-to-reference", no_argument, &flag.match_chrom_to_reference, 1 }
        #define _VV {"check-latest",     no_argument,       &flag.check_latest,     1 }
        #define _DV {"debug-latest",     no_argument,       &flag.debug_latest,     1 }
        #define _00 {0, 0, 0, 0                                                       }

        typedef const struct option Option;
        static Option genozip_lo[]    = { _lg, _i, _I, _d, _f, _h,     _D,    _L1, _L2, _q, _Q, _t, _DL, _V, _z, _zb, _zB, _zs, _zS, _zq, _zQ, _za, _zA, _zf, _zF, _zc, _zC, _zv, _zV, _zy, _zY, _m, _th,     _o, _p, _e, _E, _ch,                                                                                                _sl, _sL, _ss, _SS,      _sd, _sT,  _sF, _sK, _sb, _lc, _lh, _lH, _s2, _s7, _S7, _S8, _S9, _sa, _st, _sm, _sh, _si, _Si, _Sh, _sr, _su, _sv,      _sn,               _B, _xt, _dm, _dp, _dL, _dD, _dq, _dt, _dw, _dM, _dr, _dR, _dP, _dG, _dN, _dQ, _dl, _dc, _dg,                _sy,    _dh,_dS, _bS, _9, _99, _9s, _9P, _9G, _9g, _9V, _9Q, _9f, _9Z, _9D, _pe,      _fa, _bs,                              _nh, _rg, _sR,      _sC,           _hC, _rA, _rI, _rS, _me, _mf, _mF,     _s5, _S5, _sM, _sA, _sB, _sP, _sc, _AL, _sI, _cn,                                    _so, _SO, _s6, _kr,         _oe, _aa, _al, _Lf, _T, _TT, _Xr, _Xd, _XS, _MR, _wM, _wm, _WM, _bi, _bl, _VV, _DV, _00 };
        static Option genounzip_lo[]  = { _lg,         _d, _f, _h, _x, _D,    _L1, _L2, _q, _Q, _t, _DL, _V, _z, _zb, _zB, _zs, _zS, _zq, _zQ, _za, _zA, _zf, _zF, _zc, _zC, _zv, _zV, _zy, _zY, _m, _th, _u, _o, _p, _e,                                                                                                                   _ss, _SS, _sG, _sd, _sT,  _sF,      _sb, _lc, _lh, _lH, _s2, _s7, _S7, _S8, _S9, _sa, _st, _sm, _sh, _si, _Si, _Sh, _sr, _su, _sv,      _sn, _ov,              _xt, _dm, _dp,      _dD,      _dt,                _dR,                          _dc,                     _sy,                                                                             _pt,                                                  _sR,      _sC,           _hC, _rA,      _rS,                    _s5, _S5, _sM, _sA, _sB,           _AL, _sI, _cn, _pg, _PG, _sx, _SX, _ix,                     _s6,              _oe,                _T,                                                             _00 };
        static Option genocat_lo[]    = { _lg,         _d, _f, _h, _x, _D,    _L1, _L2, _q, _Q,          _V, _z, _zb, _zB, _zs, _zS, _zq, _zQ, _za, _zA, _zf, _zF, _zc, _zC, _zv, _zV, _zy, _zY,     _th,     _o, _p, _e,           _lo, _il, _r, _R, _Rg, _s, _sf, _sq, _G, _1, _H0, _H1, _H2, _H3, _Gt, _So, _Io, _IU, _iu, _GT,          _ss, _SS, _sG, _sd, _sT,  _sF, _sK, _sb, _lc, _lh, _lH, _s2, _s7, _S7, _S8, _S9, _sa, _st, _sm, _sh, _si, _Si, _Sh, _sr, _su, _sv, _sH, _sn, _ov, _R1, _R2,    _xt, _dm, _dp,      _dD,      _dt,                _dR,                          _dc, _ds, _sS, _SC, _sY, _sy,                                                                             _pt,                 _fs, _g, _gw, _n, _nt, _nh,      _sR,      _sC, _sD, _cC, _hC, _rA, _rI, _rS,                    _s5, _S5, _sM, _sA, _sB,           _AL, _sI, _cn, _pg, _PG, _sx, _SX, _ix, _ct, _vl,      _SO, _s6, _kr, _kR,    _oe,      _al,      _T,                                                             _00 };
        static Option genols_lo[]     = { _lg,             _f, _h,        _l, _L1, _L2, _q,              _V,                                                                                                      _p,                                                                                                                                                 _sF,                                                        _st, _sm,                                                                     _dm,                     _dt,                                                                                                                                                                                                                                                                                   _S5, _sM,                                                                                               _b, _oe,                _T,                                                             _00 };
        static Option *long_options[] = { genozip_lo, genounzip_lo, genols_lo, genocat_lo }; // same order as ExeType

        // include the option letter here for the short version (eg "-t") to work. ':' indicates an argument.
        static rom short_options[NUM_EXE_TYPES] = { // same order as ExeType
            "z:i:I:cdfhLqQt^Vzm@:o:p:B:9wWFe:E:C:2K:T:Db",   // genozip (note: includes some genounzip options to be used in combination with -d)
            "cz:fhLqdQt^V@:u:o:p:me:wWxT:D",                 // genounzip
            "hLVp:qfblT:",                                   // genols
            "z:hLV@:dp:qQ1r:R:s:H1Go:fg:e:E:wWxvk:K:n:yT:D"  // genocat
        };

        int option_index = -1;
        int c = getopt_long (argc, argv, short_options[exe_type], long_options[exe_type], &option_index);

        if (c == -1) break; // no more options

        option_is_short[c] = (option_index == -1); // true if user provided a short option - eg -c rather than --stdout

        switch (c) {
            case HELP:
                flag.help = optarg;
                goto verify_command;

            case LICENSE :
                flag.lic_width = optarg ? atoi (optarg) : 0;
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
            case '2' : flag.pair    = PAIR_READ_1   ; break;
            case '9' : flag.optimize      = 1       ; break;
            case 'b' : flag.bytes         = 1       ;        // genols
                       flag.best = optarg ? 2 : 1   ; break; // genozp: --best=NO_REF (or any other string) overrides reference requirement for SAM/BAM/FASTQ
            case 'B' : flag.vblock        = optarg  ; break;
            case 'C' : flag.reading_chain = optarg  ; break;  
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
            case 'k' : kraken_set_taxid (optarg)    ; break; // a 1-7 digit NCBI taxonomy ID
            case 'K' : flag.reading_kraken = optarg ; break;  
            case 'l' : flag.list          = 1       ; break;
            case 'm' : flag.md5           = 1       ; break;
            case 'n' : flag_set_lines (optarg)      ; break;
            case 'o' : flag_set_out_filename(optarg); break;
            case 'p' : crypt_set_password (optarg)  ; break;
            case 'q' : flag.quiet         = 1       ; break;
            case 'Q' : option_noisy       = 1       ; break;
            case 'r' : flag.regions       = 1       ; regions_add (optarg); break;
            case 'R' : flag.regions       = 1       ; regions_add_by_file (optarg); break;
            case 'T' : flags_set_file_from (optarg) ; break;
            case 's' : flag.samples       = 1       ; vcf_samples_add (optarg); break;
            case 't' : flag.test          = 1       ; break; 
            case 'u' : flag.unbind        = optarg  ; break;
            case 'v' : flag.luft          = 1       ; break;
            case 'y' : flag.show_dvcf     = 1       ; break;
            case 'w' : flag.show_stats    = STATS_SHORT ; break;
            case 'W' : flag.show_stats    = STATS_LONG  ; break;
            case 'x' : flag.index_txt     = 1       ; break;
            case 'z' : flags_set_bgzf (optarg)      ; break;
            case '~' : flag.show_is_set   = optarg  ; break;
            case 1   : flag.show_time = optarg ? optarg : "" ; break; // show_time with or without a specific member of ProfilerRec
            case 2   : if (optarg) flag.dict_id_show_one_b250 = dict_id_make (optarg, strlen (optarg), DTYPE_PLAIN); 
                       else        flag.show_b250 = 1;
                       break;
            case 3   : if (optarg) flag.show_one_dict = optarg;
                       else        flag.show_dict = 1;
                       break;
            case 4   : flag.show_mutex    = optarg ? optarg : (char*)1; break;
            case 5   : flag.dump_one_b250_dict_id  = dict_id_make (optarg, strlen (optarg), DTYPE_PLAIN); break;
            case 6   : flag.dump_one_local_dict_id = dict_id_make (optarg, strlen (optarg), DTYPE_PLAIN); break;
            case 7   : flag.dump_section  = optarg  ; break;
            case 8   : ASSINP0 (str_get_int_range32 (optarg, 0, 1, 10000000, (int32_t*)&flag.one_vb), 
                                "--one-vb expects a 1-based VBlock number"); break;
            case 9   : flags_set_downsample (optarg); break;
            case 10  : flag.show_headers = 1 + sections_st_by_name (optarg); break; // +1 so SEC_NONE maps to 0
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
            case 20  : flag.count = optarg ? COUNT_VBs : CNT_TOTAL; break;
            case 21  : kraken_set_show_kraken (optarg); break;
            case 22  : flag_set_head (optarg)       ; break;
            case 23  : flag_set_tail (optarg)       ; break;
            case 24  : iupac_set (optarg)           ; break;
            case 26  : license_set_filename (optarg); break;
            case 27  : tar_set_tar_name (optarg)    ; break;
            case 28  : flag.do_register = optarg ? optarg : ""; break;
            case 29  : flag_set_interleaved (optarg); break;
            case 128 : ASSERTRUNONCE ("--dvcf-rename option can only appear once, see:  " WEBSITE_DVCF);
                       flag.dvcf_rename = optarg    ; break;
            case 129 : ASSERTRUNONCE ("--dvcf-drop option can only appear once, see:  " WEBSITE_DVCF);
                       flag.dvcf_drop = optarg      ; break;
            case 130 : flag.out_dt = DT_FASTQ; flag.extended_translation = !!optarg; break; 
            case 132 : flag_set_show_containers (optarg); break; 
            case 133 : flag.debug_seg=1;
                       if (optarg) flag.dict_id_debug_seg = dict_id_make (optarg, strlen (optarg), DTYPE_PLAIN); 
                       break;
            case 134 : biopsy_init (optarg);        ; break;
            case 135 : flag.show_gheader = optarg ? atoi(optarg) : 1; break; // =1 show gheader as in file, =2 show shows section list after possible modiciation by writer_create_plan 
            case 136 : flag.show_sa = optarg ? atoi(optarg)+1 : -1; break;   //-1=show all, >=1 - show grp_i=show_sa-1 
            case 137 : flag_set_biopsy_line (optarg); break;
            case 138 : flag_set_stats (STATS_SHORT, optarg); break;
            case 139 : flag_set_stats (STATS_LONG,  optarg); break;
            case 0   : break; // a long option that doesn't have short version will land here - already handled so nothing to do
                 
            default  : // unrecognized option 
                fprintf (stderr, "Usage: %s [OPTIONS] filename1 filename2...\nManual page: "GENOZIP_URL"/%s.html\n", global_cmd, global_cmd);
                exit (EXIT_GENERAL_ERROR);  
        }
    }
    flag.explicit_out_dt = (flag.out_dt >= 0);
}

static void flags_warn_if_duplicates (int num_files, rom *filenames)
{
    if (num_files <= 1) return; // nothing to do

    # define BASENAME_LEN 256
    char basenames[num_files * BASENAME_LEN];

    for (unsigned i=0; i < num_files; i++)
        file_basename (filenames[i], false, "", &basenames[i*BASENAME_LEN], BASENAME_LEN);

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

    char s[20]; 

    CONFLICT (flag.unbind,      flag.out_filename,   OT("unbind", "u"),    OT("output", "o"));
    CONFLICT (flag.subdirs,     flag.out_filename,   OT("subdirs", "D"),   OT("output", "o")); 
    CONFLICT (flag.to_stdout,   flag.index_txt,      OT("stdout", "c"),    OT("index", "x"));
    CONFLICT (flag.one_vb,      flag.interleaved,     "--interleaved",       "--one-vb");
    CONFLICT (flag.one_component, flag.interleaved,   "--interleaved",       "--R1/--R2");
    CONFLICT (flag.one_component, flag.luft,         OT("luft", "v"),      "--R1/--R2");
    CONFLICT (flag.xthreads,    flag.interleaved,     "--interleaved",       "--xthreads");
    CONFLICT (flag.snps_only,   flag.indels_only,    "--snps_only",        "--indels-only");
    CONFLICT (flag.header_only, flag.seq_only,       "--seq-only",         "--header-only");
    CONFLICT (flag.header_only, flag.qual_only,      "--qual-only",        "--header-only");
    CONFLICT (flag.seq_only,    flag.qual_only,      "--seq-only",         "--qual-only");
    CONFLICT (flag.header_only, flag.interleaved,     "--interleaved",       "--header-only");
    CONFLICT (flag.regions,     flag.interleaved,     "--interleaved",       "--regions");
    CONFLICT (flag.header_only, flag.no_header==1,   OT("no-header", "H"), "--header-only");
    CONFLICT (flag.no_header,   flag.header_one,     OT("no-header", "H"), OT("header-one", "1"));
    CONFLICT (flag.quiet,       option_noisy,        OT("quiet", "q"),     OT("noisy", "Q"));
    CONFLICT (flag.test,        flag.biopsy,         OT("test", "t"),      "--biopsy");
    CONFLICT (flag.test,        flag.bgzf != FLAG_BGZF_BY_ZFILE, OT("test", "t"),      OT("bgzf", "z"));
    CONFLICT (flag.test,        flag.add_line_numbers, OT("test", "t"),    "--add-line-numbers");
    CONFLICT (flag.md5,         flag.add_line_numbers, OT("md5", "m"),     "--add-line-numbers");
    CONFLICT (flag.test,        flag.optimize,       OT("test", "t"),      OT("optimize", "9"));
    CONFLICT (flag.md5,         flag.optimize,       OT("md5", "m"),       OT("optimize", "9"));
    CONFLICT (flag.test,        flag.match_chrom_to_reference, OT("test", "t"), "--match-chrom-to-reference");
    CONFLICT (flag.md5,         flag.match_chrom_to_reference, OT("md5", "m"),  "--match-chrom-to-reference");
    CONFLICT (flag.test,        flag.reading_chain,  OT("test", "t"),      OT("chain", "C"));
    CONFLICT (flag.md5,         flag.reading_chain,  OT("md5", "m"),       OT("chain", "C"));
    CONFLICT (flag.make_reference, flag.reading_chain, "--make-reference", OT("chain", "C"));
    CONFLICT (flag.samples,     flag.drop_genotypes, OT("samples", "s"),   OT("drop-genotypes", "G"));
    CONFLICT (flag.best,        flag.fast,           OT("best", "b"),      OT("fast", "F"));
    CONFLICT (flag.show_sex,    flag.regions==1,     "--sex",              OT("regions", "r"));
    CONFLICT (flag.show_sex,    flag.out_filename,   "--sex",              OT("output", "o"));
    CONFLICT (flag.show_sex,    flag.grep,           "--sex",              OT("grep", "g"));
    CONFLICT (flag.show_sex,    flag.idxstats,       "--sex",              "--idxstats");
    CONFLICT (flag.show_sex,    flag.count,          "--sex",              "--count");
    CONFLICT (flag.show_sex,    flag.test,           "--show_sex",         OT("test", "t"));
    CONFLICT (flag.show_coverage, flag.idxstats,     "--coverage",         "--idxstats");
    CONFLICT (flag.show_coverage, flag.regions==1,   "--coverage",         OT("regions", "r"));
    CONFLICT (flag.show_coverage, flag.out_filename, "--coverage",         OT("output", "o"));
    CONFLICT (flag.show_coverage, flag.grep,         "--coverage",         OT("grep", "g"));
    CONFLICT (flag.show_coverage, flag.count,        "--coverage",         "--count");
    CONFLICT (flag.show_coverage, flag.test,         "--coverage",         OT("test", "t"));
    CONFLICT (flag.idxstats,    flag.test,           "--idxstats",         OT("test", "t"));
    CONFLICT (flag.idxstats,    flag.count,          "--idxstats",         "--count");
    CONFLICT (flag.idxstats,    flag.show_sex,       "--idxstats",         "--sex");
    CONFLICT (flag.idxstats,    flag.show_coverage,  "--idxstats",         "--coverage");
    CONFLICT (flag.count,       flag.out_filename,   "--count",            OT("output", "o"));
    CONFLICT (flag.count,       flag.header_only,    "--count",            OT("header-only", "H"));
    CONFLICT (flag.count,       flag.header_one,     "--count",            OT("header-one", "1"));
    CONFLICT (flag.test,        flag.make_reference, "--make-reference",   OT("test", "t"));
    CONFLICT (flag.sort,        flag.unsorted,       "--sort",             "--unsorted");
    CONFLICT (flag.show_ostatus,flag.show_dvcf,      "--show-ostatus",     "--show-dvcf");
    CONFLICT (flag.multiseq,    flag.make_reference, "--make-reference",   "--multiseq");
    CONFLICT (flag.reference,   flag.make_reference, "--make-reference",   OT("reference", "e"));
    CONFLICT (flag.dump_one_b250_dict_id.num, flag.dump_one_local_dict_id.num, "--dump-one-b250", "--dump-one-local");
    CONFLICT (flag.show_stats,  flag.grep,           OT("stats", "w"),     OT("grep", "g"));
    CONFLICT (flag.show_stats,  flag.lines_first != NO_LINE, OT("stats", "w"),     OT("lines", "n"));
    CONFLICT (flag.show_stats,  flag.downsample,     OT("stats", "w"),     "--downsample");
    CONFLICT (flag.show_stats,  flag.header_only,    OT("stats", "w"),     "--header-only");
    CONFLICT (flag.show_stats,  flag.MAPQ,           OT("stats", "w"),     "--MAPQ");
    CONFLICT (flag.show_stats,  flag.FLAG,           OT("stats", "w"),     "--FLAG");
    CONFLICT (flag.show_stats,  flag.one_component,  OT("stats", "w"),     "--R1/--R2");
    CONFLICT (flag.show_stats,  flag.one_vb,         OT("stats", "w"),     "--one-vb");
    CONFLICT (flag.show_stats,  flag.kraken_taxid != TAXID_NONE, OT("stats", "w"), "--taxid");
    CONFLICT (flag.show_stats,  flag.sequential,     OT("stats", "w"),     "--sequential");
    CONFLICT (flag.show_stats,  flag.tail,           OT("stats", "w"),     "--tail");
    CONFLICT (flag.show_stats,  flag.drop_genotypes, OT("stats", "w"),     OT("drop-genotypes", "G"));
    CONFLICT (flag.show_stats,  flag.regions,        OT("stats", "w"),     OT("regions", "r"));
    CONFLICT (flag.show_stats,  flag.samples,        OT("stats", "w"),     OT("samples", "s"));
    CONFLICT (flag.biopsy,      flag.out_filename,   "--biopsy",           OT("output", "o"));
    CONFLICT (flag.biopsy_line.line_i>=0, flag.out_filename, "--biopsy-line", OT("output", "o"));
    
    if (command == PIZ) {
        CONFLICT (flag.test, flag.out_filename,      OT("output", "o"),    OT("test", "t"));
        CONFLICT (flag.test, flag.replace,           OT("replace", "^"),   OT("test", "t"));

        // options that require --reference
        if (!flag.reference) {
            ASSINP0 (!flag.show_ref_iupacs, "--show-ref-iupacs requires --reference");
            ASSINP0 (!flag.show_ref_diff, "--show-ref-diff requires two --reference arguments");
        }
        // options that require absence of --reference 
        else {
            ASSINP (!flag.show_reference, "%s is incompatible with --show-reference. To use --show-reference on a reference file, omit --reference", OT("reference", "e"));
            ASSINP (!flag.show_ref_seq, "%s is incompatible with --show-ref-seq. To use --show-ref-seq on a reference file, omit --reference", OT("reference", "e"));
        }
    }

    // some genozip flags are allowed only in combination with --decompress 
    if (exe_type == EXE_GENOZIP && command == ZIP) {
        NEED_DECOMPRESS (flag.bgzf >= 0, OT("bgzf", "z"));
        NEED_DECOMPRESS (flag.out_dt != DT_NONE, str_tolower (dt_name (flag.out_dt), s));
        NEED_DECOMPRESS (flag.unbind, OT("unbind", "u"));
        NEED_DECOMPRESS (flag.show_aliases, "--show_aliases");
        NEED_DECOMPRESS (flag.show_is_set, "--show_is_set");
    }

    ASSINP (flag.reference != REF_EXT_STORE || exe_type != EXE_GENOCAT, "option %s supported only for viewing the reference file itself", OT("REFERENCE", "E"));

    // --output only allowed with a single file, or two FASTQs and --pair
    ASSINP0 (!flag.out_filename || num_files <= 1 || (num_files==2 && flag.pair), 
        flag.out_dt == DT_FASTQ ? "--output can only be used with a single input file, or using --pair and two input files. Use --tar to archive multiple files. See " WEBSITE_ARCHIVING
                                : "--output can only be used with a single input file. Use --tar to archive multiple files. See " WEBSITE_ARCHIVING);

    ASSINP0 (flag.reading_chain || !flag.dvcf_rename, "--dvcf-rename can only be used in combination with --chain. See: " WEBSITE_DVCF);
    ASSINP0 (flag.reading_chain || !flag.dvcf_drop, "--dvcf-drop can only be used in combination with --chain. See: " WEBSITE_DVCF);

    ASSINP0 (exe_type!=EXE_GENOUNZIP || !flag.one_vb || flag.test, "--one-vb can only be used with: (1) genocat, or (2) genounzip --test");
}

// ZIP: --pair: verify an even number of fastq files, --output, and --reference/--REFERENCE
static void flags_verify_pair_rules (unsigned num_files, rom *filenames)
{
    // verify even number of files
    ASSINP (num_files % 2 == 0, "when using %s, expecting an even number of FASTQ input files, each consecutive two being a pair", OT("pair", "2"));
    ASSINP (flag.reference, "either --reference or --REFERENCE must be specified when using %s", OT("pair", "2"));

    // verify all are fastq
    for (unsigned i=0; i < num_files; i++)
        ASSINP (txtfile_get_file_dt (filenames[i]) == DT_FASTQ, "when using %s, all input files are expected to be FASTQ files, but %s is not", OT("pair", "2"), filenames[i]);

    // if which --output is missing, we check if every pair of files has a consistent name
    if (!flag.out_filename) 
        for (unsigned i=0; i < num_files; i += 2) 
            ASSINP (file_get_fastq_pair_filename (filenames[i], filenames[i+1], true),  
                    "to use %s without specifying --output, the naming of the files needs to be consistent and include the numbers 1 and 2 respectively, but these files don't: %s %s", 
                    OT("pair", "2"), filenames[i], filenames[i+1]);
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
    if (flag.biopsy || flag.biopsy_line.line_i != NO_LINE) {
        ASSERTW (!flag.test,         "FYI: the %s option is ignored when taking a biopsy", OT("test", "t"));
        ASSERTW (!flag.out_filename, "FYI: the %s option is ignored when taking a biopsy", OT("output", "o"));
        flag.test = false;
        flag.out_filename = NULL;
    }
    
    flags_test_conflicts (num_files);

    flag.debug_top = flag.echo || getenv ("GENOZIP_TEST");

    // verify stuff needed for --pair
    if (flag.pair) flags_verify_pair_rules (num_files, filenames); // --pair is only available in ZIP
        
    // don't show progress for flags that output throughout the process. no issue with flags that output only in the end
    if (flag.show_dict || flag.show_b250 || flag.show_headers || flag.show_threads || flag.show_bgzf || flag.show_mutex ||
        flag.dict_id_show_one_b250.num || flag.show_one_dict || flag.show_one_counts.num || flag.show_sa || flag.show_depn || 
        flag.show_reference || flag.show_digest || flag.list_chroms || flag.show_coverage == COV_ONE || flag.show_ranges ||
        flag.show_alleles || flag.show_vblocks || flag.show_codec || flag.debug_gencomp || flag.show_qual || flag.show_aligner ||
        (flag.show_index && command==PIZ))
        flag.quiet=true; // don't show progress or warnings

    // override these ^ if user chose to be --noisy
    if (option_noisy) flag.quiet=false;

    // if genocat --sex, we only need the X,Y,1 chromosomes
    if (flag.show_sex && !flag.show_coverage && exe_type == EXE_GENOCAT) {
        flag.regions = 2;    // set, but not to 1, to avoid CONFLICT error
        regions_add ("X,chrX,ChrX,Y,chrY,ChrY,1,chr1,Chr1");
    }

    flag.multiple_files = (num_files > 1);

    flag.longest_filename = flags_get_longest_filename (num_files, filenames);

    // if using the -o option - check that we don't have duplicate filenames (even in different directory) as they
    // will overwrite each other if extracted genounzip
    if (command == ZIP && flag.out_filename && !flag.quiet) flags_warn_if_duplicates (num_files, filenames);

    // if reference is not specified, attempt to set it from env var $GENOZIP_REFERENCE
    if (!flag.reference) {
        ref_set_reference (gref, NULL, REF_EXTERNAL, false);
        if (flag.reading_chain) ref_set_reference (prim_ref, NULL, REF_EXTERNAL, false);
    }

    // genocat --show-chain <chain-file> - actually load the chain
    if ((flag.show_chain || flag.show_chain_contigs) && exe_type == EXE_GENOCAT && num_files==1)
        flag.reading_chain = filenames[0];

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

    if (exe_type == EXE_GENOCAT) flag.only_headers = flag.show_headers;

    // cases where we don't need to load the reference file, even if the genozip file normally needs it
    // note: we don't exclude due to collect_coverage here, instead we do it in main_load_reference
    // note: this is here and not in flags_update_piz_one_file bc it is consumed by main_load_reference
    // note: this is in flags_update and not flags_update_piz, bc reference file is loaded first 
    flag.genocat_no_ref_file = exe_type == EXE_GENOCAT &&
        (flag.show_stats || flag.show_dict || flag.show_b250 || flag.list_chroms || flag.show_one_dict ||
         flag.dump_one_b250_dict_id.num || // all other sections (except CHROM) are blocked from reading in piz_default_skip_section
         flag.show_index || flag.dump_section || flag.show_one_counts.num || flag.show_flags ||
         flag.show_aliases || flag.show_txt_contigs || flag.show_gheader || flag.show_recon_plan || flag.show_ref_contigs ||
         (flag.count && !flag.bases) ||
         flag.collect_coverage); // note: this is updated in flags_update_piz_one_file
}

// ZIP: called for each file, after opening txt and z files, but before calling zip_one_file 
void flags_update_zip_one_file (void)
{
    DataType dt = z_file->data_type;

    // --make-reference implies --md5 --B1 (unless --vblock says otherwise), and not encrypted. 
    // in addition, txtfile_read_vblock() limits each VB to have exactly one contig.
    if (flag.make_reference) {
        ASSINP (!crypt_have_password(), "option --make-reference is incompatible with %s", OT("password", "p"));
        flag.md5 = true;
    }

    ASSINP0 (chain_is_loaded || dt == DT_CHAIN || !ref_get_filename (prim_ref), 
             "Two --reference arguments were specified. This is not expected, unless --chain is specified as well, or compressing a chain file");

    ASSINP0 (!flag.reading_chain || !ref_get_filename (gref) || ref_get_filename (prim_ref),
             "When using --chain, you either specify two --reference arguments or none at all. If specified, the first is the reference file in primary "
             "coordinates (i.e. those of the original file), and the second is the reference file in luft coordinates (i.e. the coordinates aftering lifting). See "WEBSITE_DVCF);
    
    // if --optimize was selected, all optimizations are turned on
    if (flag.optimize) switch (dt) {
        case DT_BCF   :
        case DT_VCF   : flag.GP_to_PP = flag.GL_to_PL = flag.optimize_phred = flag.optimize_VQSLOD = flag.optimize_sort = true; break;
        case DT_GFF3  : flag.optimize_sort = flag.optimize_Vf = true; break;
        case DT_BAM   :
        case DT_SAM   : flag.optimize_QUAL = flag.optimize_ZM = true; break;
        case DT_FASTQ : flag.optimize_QUAL = flag.optimize_DESC = true; break;
        default: break;
    }
    
    // if any optimization flag is on, we turn on flag.optimize
    else flag.optimize = flag.optimize_sort || flag.optimize_phred || flag.GL_to_PL || flag.GP_to_PP || flag.optimize_VQSLOD ||
                         flag.optimize_QUAL || flag.optimize_Vf || flag.optimize_ZM || flag.optimize_DESC;

    // cases where txt data is modified during Seg - digest is not stored, it cannot be tested with --test and other limitations 
    // note: this flag is also set when the file header indicates that it's a Luft file. See vcf_header_get_dual_coords().
    // note: we don't set this for comp_i > 1, bc DVCF: when compressing a primary DC file it is reconstructed without modification
    flag.data_modified = flag.data_modified  // this is needed, eg, so when compressinga Luft file, the rejects file inherits data_modified 
                      || flag.optimize       // we're modifying data to make it more compressible
                      || flag.match_chrom_to_reference
                      || (chain_is_loaded && dt == DT_VCF) // converting a standard VCF to a dual-coordinates VCF
                      || (flag.add_line_numbers && dt == DT_VCF)
                      || (kraken_is_loaded && (dt == DT_SAM || dt == DT_BAM)) // adding a tx:i optional field
                      || flag.lines_last != NO_LINE  // --head advanced option to compress only a few lines
                      || flag.biopsy_line.line_i != NO_LINE;

    if (chain_is_loaded && dt == DT_VCF && !flag.show_one_counts.num && !flag.quiet)
        flag.show_one_counts = dict_id_typeless ((DictId)_VCF_oSTATUS);

    info_stream = stdout; // always stdout in zip
    is_info_stream_terminal = isatty (fileno (info_stream)); 

    ASSINP0 (!flag.match_chrom_to_reference || flag.reference, "--match-chrom-to-reference requires using --reference as well"); 

    flag.bind = (dt == DT_FASTQ && flag.pair)     ? BIND_FQ_PAIR 
              : (dt == DT_SAM || dt == DT_BAM)    ? BIND_SAM
              : (dt == DT_VCF && chain_is_loaded) ? BIND_DVCF
              :                                     BIND_NONE;

    // if biopsy, we seg only for speed. current limitation: in paired files we do the whole thing. TO DO: fix this.
    if (flag.biopsy_line.line_i != NO_LINE && flag.pair != BIND_FQ_PAIR)
        flag.seg_only = true;

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

    // SAM
    FLAG_ONLY_FOR_DT(SAM,          show_bam,      "show_bam");
    FLAG_ONLY_FOR_DT(SAM,          optimize_ZM,   "optimize-ZM");
    FLAG_ONLY_FOR_2DTs(SAM, FASTQ, optimize_QUAL, "optimize-QUAL");
    if (dt == DT_SAM || dt == DT_BAM) {
        flag.bind = BIND_SAM; // initialized here, if no PRIM/DEPN lines exist, we will cancel it sam_zip_generate_recon_plan
        if (flag.show_bam) 
            flag.seg_only = flag.xthreads = flag.quiet = true; 
    }
    
    // FASTQ
    FLAG_ONLY_FOR_DT(FASTQ, pair,           "pair");
    FLAG_ONLY_FOR_DT(FASTQ, optimize_DESC,  "optimize-DESC");
    FLAG_ONLY_FOR_2DTs(FASTQ, FASTA, multiseq, "multiseq");

    // VCF
    FLAG_ONLY_FOR_DT(VCF, sort,             "sort");
    FLAG_ONLY_FOR_DT(VCF, dvcf_rename,      "dvcf-rename");
    FLAG_ONLY_FOR_DT(VCF, show_lift,        "show-lifts");
    FLAG_ONLY_FOR_DT(VCF, show_chain,       "show-chain");
    FLAG_ONLY_FOR_DT(VCF, show_rename_tags, "show-rename-tags");
    FLAG_ONLY_FOR_DT(VCF, unsorted,         "unsorted");
    FLAG_ONLY_FOR_DT(VCF, add_line_numbers, "add-line-numbers");
    FLAG_ONLY_FOR_DT(VCF, optimize_phred,   "optimize-phred");
    FLAG_ONLY_FOR_DT(VCF, GL_to_PL,         "GL-to-PL");
    FLAG_ONLY_FOR_DT(VCF, GP_to_PP,         "GP-to-PP");
    FLAG_ONLY_FOR_DT(VCF, optimize_VQSLOD,  "optimize-VQSLOD");
    FLAG_ONLY_FOR_2DTs(VCF, GFF3, optimize_sort, "optimize-sort");

    // FASTA
    ASSINP0 (!flag.make_reference || dt == DT_REF, "--make-reference is only supported for FASTA files"); // data_type is DT_REF

    // GFF3
    FLAG_ONLY_FOR_DT(GFF3, optimize_Vf,     "optimize-Vf");
    
    // GENERIC
    FLAG_NOT_FOR_DT (GENERIC, debug_lines,   "debug-lines"); // GENERIC doesn't have lines 

    if (flag.show_flags) flags_show_flags();
}

bool flags_is_genocat_global_area_only (void)
{
    return exe_type == EXE_GENOCAT &&
        (flag.show_stats || flag.show_dict || flag.list_chroms || flag.show_one_dict ||
         flag.show_index || flag.show_one_counts.num || command == SHOW_HEADERS ||
         flag.show_reference || flag.show_ref_contigs || flag.show_ranges ||
         flag.show_ref_index || flag.show_ref_hash || flag.show_chrom2ref || 
         flag.show_ref_seq || flag.show_aliases || flag.show_gheader==1 ||
         (z_file->data_type == DT_REF && flag.show_ref_iupacs));    
}

// PIZ: called after opening z_file and reading the header before opening txt_file
void flags_update_piz_one_file (int z_file_i /* -1 if unknown */)
{
    DataType dt = z_file->data_type;

    if (flag.out_dt == DT_NONE) {

        // handle native binary formats (BAM). note on BCF and CRAM: we used bcftools/samtools as an external 
        // compressor, so that genozip sees the text, not binary, data of these files - the same as if the file were compressed with eg bz2
        if (z_file->z_flags.txt_is_bin) {
            
            // PIZ of a genozip file with is_binary (e.g. BAM) is determined here unless the user overrides with --sam, --bam, --fastq or --bgzf
            if (exe_type == EXE_GENOCAT && !flag.out_filename && isatty(1) && flag.bgzf == FLAG_BGZF_BY_ZFILE) 
                flag.out_dt = dt; // output in textual format
            else
                flag.out_dt = DTPZ (bin_type);   // output in binary format
        }
        else
            flag.out_dt = dt;
    }

    ASSINP (!ref_get_filename (prim_ref) || !flag.explicit_ref, "More than one --reference argument was specified when reading %s", z_name);

    flag.collect_coverage = flag.show_sex || flag.show_coverage || flag.idxstats;
    
    // non-translated mode needed for coverage collection
    if (flag.collect_coverage && dt == DT_SAM) 
        flag.out_dt = DT_SAM; 

    // --luft is only possible on dual-coordinates files
    ASSINP (!flag.luft || z_is_dvcf, "--luft is not possible for %s because it is not a dual-coordinates file", z_name);

    // genounzip not possible on dual coord files (bc it unbinds and we need to concatenate)
    ASSINP (!z_is_dvcf || exe_type == EXE_GENOCAT, "Cannot access dual-coordinates file %s with genounzip, use genocat instead", z_name);
                
    ASSINP0 (exe_type != EXE_GENOUNZIP || !flag.to_stdout, "Cannot use --stdout with genounzip, use genocat instead");

    ASSINP (exe_type != EXE_GENOUNZIP || !flag.out_filename || 
            z_file->num_components <= (1 + z_is_dvcf + 2 * z_sam_gencomp), // allow single-component, DVCF, and SAM with genocomp
            "Cannot use --output because %s is a bound file containing multiple components. Use --prefix to set a prefix for output filenames, use genocat to output as a single concatenated file, or use genols to see the components' metadata",
            z_name);
             
    // --phylip implies sequential and header_one
    if (dt == DT_FASTA && flag.out_dt == DT_PHYLIP) {
        flag.sequential = 1;
        flag.header_one = flag.no_header = flag.header_only = 0;
    }

    // --downsample in FASTA implies --sequential
    if (dt == DT_FASTA && flag.downsample)
        flag.sequential = 1;

    if (!z_is_dvcf && flag.luft) {
        WARN ("%s: ignoring the --luft option, because file was not compressed with --chain", z_name);
        flag.luft = 0;
    }

    // check validity of --one-vb
    ASSINP (flag.one_vb <= z_file->num_vbs, "%s: --one-vb=%u but file has only %u VBlocks", z_name, flag.one_vb, z_file->num_vbs);

    // Check if the reconstructed data type is the same as the source data type
    bool is_binary = z_file->z_flags.txt_is_bin;
    flag.reconstruct_as_src = (flag.out_dt == DT_SAM && dt==DT_SAM && !is_binary) || 
                              (flag.out_dt == DT_BAM && dt==DT_SAM && is_binary ) ||
                              (flag.out_dt == dt     && dt!=DT_SAM);

    ASSINP (exe_type == EXE_GENOCAT || flag.reconstruct_as_src, 
            "genozip file %s is of type %s, but output file is of type %s. Translating between types is not possible in genounzip, use genocat instead",
            z_name, (dt==DT_SAM && is_binary) ? "BAM" : dt_name (dt), dt_name (flag.out_dt));             

    // for FASTA/Q we convert a "header_only" flag to "header_only_fast" because in FASTA/Q
    // this flag affects the reconstruction of the body rather than the txt header (FASTA/Q don't have a txt header)
    if (flag.header_only && (flag.out_dt == DT_FASTA || flag.out_dt == DT_FASTQ)) {
        flag.header_only      = false;
        flag.header_only_fast = true;
    }

    // this affects reading an implicit reference specified in the file's SEC_GENOZIP_HEADER. We do it here instead of flags_update because
    // we first need to unset header_only for FASTA and FASTQ
    flag.genocat_no_ref_file |= (exe_type == EXE_GENOCAT && flag.header_only && !flag.luft);

    // genocat quits after reading global area - it doesn't read TXT_HEADER or VBs
    flag.genocat_global_area_only = flags_is_genocat_global_area_only();

    // if this flag is set, data will be read and uncompressed (unless blocked in piz_is_skip_section), 
    // but not reconstructed or written
    flag.genocat_no_reconstruct = exe_type == EXE_GENOCAT && 
        (flag.genocat_global_area_only || flag.dump_section || 
         flag.show_headers || flag.show_txt_contigs || (flag.show_recon_plan && dt != DT_FASTA));

    // if this flag is set, no data will be written, although it still could be read and reconstructed 
    // (unless blocked in flag.genocat_no_reconstruct or piz_default_skip_section) 
    flag.no_writer = exe_type == EXE_GENOCAT &&
        (flag.genocat_no_reconstruct || flag.collect_coverage || flag.show_kraken || flag.count ||
         flag.dump_one_local_dict_id.num || flag.dump_one_b250_dict_id.num || flag.show_b250 || flag.show_sa);

    flag.no_writer |= exe_type == EXE_GENOUNZIP && flag.test;

    // usually, if no_writer, we also don't need the writer thread, but there are exceptions to that rule
    flag.no_writer_thread = flag.no_writer && 
                            !(z_sam_gencomp && flag.test); // exception: is SAM with generated components, we need the Writer thread to calculate the digest

    ASSINP0 (!flag.no_writer || !flag.index_txt, "--index cannot be used with this command line option combination");

    flag.genocat_no_dicts = exe_type == EXE_GENOCAT &&
        (flag.show_stats || flag.show_index || flag.show_one_counts.num || // note: we need dicts for dump_b250 as we need to reconstruct
         flag.show_reference || flag.show_ref_contigs || 
         flag.show_ref_index || flag.show_ref_hash || flag.show_chrom2ref || 
         flag.show_ref_seq || flag.show_aliases || flag.show_txt_contigs || flag.show_gheader || 
         (flag.show_recon_plan && dt != DT_FASTA)); // we need dicts to generate the FASTA plan (filter for grep etc)

    // don't show progress or warning when outputing to stdout (note: we are "quiet" even if output doesn't go to the terminal
    // because often it will be piped and ultimately go the terminal - a user can override this with --noisy)
    if (flag.to_stdout && !flag.validate && !option_noisy && !flag.no_writer) flag.quiet=true; 

    // genocat of dual coords implies sorting unless overridden with --unsorted
    if (exe_type == EXE_GENOCAT && z_is_dvcf && !flag.unsorted && (!flag.no_writer || flag.show_recon_plan)) 
        flag.sort = true;

    // cases we skip displaying the txt header (we still read the section for inspection)
    if (!flag.no_header && exe_type == EXE_GENOCAT && 
        (flag.out_dt != DT_BAM || flag.no_writer) && // it is not possible to have BAM without a header
        (flag.genocat_global_area_only || flag.count || flag.show_headers ||
        (z_file_i >= 1 && !flag.no_writer && !flag.header_only))) // when using genocat to concatenate multiple files - don't show the header for the 2nd+ file
        flag.no_header = 2; // 2 = assigned here and not from command line

    flag.maybe_txt_header_modified = exe_type == EXE_GENOCAT && 
        (flag.no_header || flag.lines_first != NO_LINE || // options that may cause dropping of the txt header
         flag.luft ||                               // --luft modifies the txt header
         (dt == DT_VCF && (flag.header_one || flag.samples || flag.drop_genotypes)) || // VCF specific options that modify the txt header
         (z_file->num_components > (1 + z_is_dvcf))); // txtheaders are dropped if concatenating

    flag.maybe_vb_dropped_by_writer = exe_type == EXE_GENOCAT && // dropped by piz_dispatch_one_vb
        (flag.lines_first != NO_LINE || // decided by writer_create_plan
         flag.tail             || // decided by writer_create_plan
         flag.downsample       || // decided by writer_create_plan
         flag.regions          || // decided by writer_create_plan
         flag.one_component    || // decided by writer_init_comp_info 
         flag.one_vb           || // decided by writer_init_vb_info
         flag.header_only);       // decided by writer_init_vb_info

    flag.maybe_vb_dropped_after_read = exe_type == EXE_GENOCAT && // dropped by piz_dispatch_one_vb
        ((flag.grep && (dt == DT_FASTQ || dt == DT_FASTA)) || // decided by piz_read_one_vb
         (flag.regions && dt == DT_FASTA)); // decided by piz_read_one_vb

    flag.maybe_lines_dropped_by_reconstructor = exe_type == EXE_GENOCAT && 
         ((dt == DT_VCF   && (flag.snps_only || flag.indels_only || z_is_dvcf)) || // vcf_lo_piz_TOPLEVEL_cb_filter_line will drop lines of the wrong coordinate
         // FASTA specific line droppers
         (dt == DT_FASTA && (flag.sequential || flag.header_only_fast || flag.header_one || flag.no_header)) || 
         // FASTQ specific line droppers
         (dt == DT_FASTQ && flag.bases) || 
         // SAM specific line droppers
         (dt == DT_SAM   && (flag.sam_flag_filter || flag.sam_mapq_filter || flag.bases || flag.out_dt == DT_FASTQ)) || 
         // General filters
         flag.kraken_taxid != TAXID_NONE || flag.grep || flag.regions || 
         // no-writer, but nevertheless modify the txt_data
         flag.collect_coverage || flag.count);

    flag.maybe_lines_dropped_by_writer = exe_type == EXE_GENOCAT && 
         (flag.downsample || flag.lines_first != NO_LINE || flag.tail);

    flag.maybe_vb_modified_by_reconstructor = exe_type == EXE_GENOCAT && 
         // translating to another data
        (!flag.reconstruct_as_src || 
         // lines may be dropped by reconstructor
         flag.maybe_lines_dropped_by_reconstructor || 
         // VCF specific VB modifiers
         (dt == DT_VCF   && (flag.samples || flag.drop_genotypes || flag.gt_only || flag.show_dvcf || flag.show_ostatus)) || 
         // FASTA specific modifiers
         (dt == DT_FASTA && (false)) || 
         // FASTQ specific modifiers
         (dt == DT_FASTQ && (flag.header_only_fast || flag.seq_only || flag.qual_only)) || // FASTQ "line" is for lines, so these are line modifications, not drops
         // SAM specific modifiers
         (dt == DT_SAM   && (flag.add_line_numbers)));

    // cases where Writer may re-order lines resulting in different ordering that within the VBs
    flag.maybe_lines_out_of_order = exe_type == EXE_GENOCAT && 
        ((dt == DT_VCF   && (z_is_dvcf || sections_get_comp_recon_plan_sec(0,0))) || // DVCF or genozip --sort
         (dt == DT_FASTQ && z_file->z_flags.dts_paired && flag.interleaved) || 
         (dt == DT_SAM   && z_file->z_flags.has_gencomp));

    // true if the PIZ output txt file will NOT be identical to the source file as recorded in z_file
    flag.data_modified = flag.maybe_txt_header_modified             || 
                         flag.maybe_vb_dropped_after_read           ||
                         flag.maybe_vb_dropped_by_writer            ||
                         flag.maybe_lines_dropped_by_reconstructor  ||
                         flag.maybe_vb_modified_by_reconstructor    || 
                         flag.maybe_lines_out_of_order;

    // calculation depends on flag.data_modified
    bool pg_line_added_to_header = ((flag.out_dt == DT_SAM || flag.out_dt == DT_BAM) && !flag.reconstruct_as_src)
                                || (flag.out_dt == DT_VCF && (flag.data_modified || z_is_dvcf));

    if (pg_line_added_to_header && !flag.no_pg && exe_type == EXE_GENOCAT) 
        flag.maybe_txt_header_modified = flag.data_modified = true;

    ASSINP0 (exe_type != EXE_GENOUNZIP || !flag.data_modified, "Data modification flags are not allowed in genounzip, use genocat instead");

    ASSINP0 (!flag.test || !flag.data_modified, "--test cannot be used when other flags specify data modification. See " WEBSITE_DIGEST);

    // cases where we don't read unnecessary contexts, and should just reconstruct them as an empty
    // string (in other cases, it would be an error)
    flag.missing_contexts_allowed = flag.collect_coverage || flag.count || flag.drop_genotypes;

    if (Z_DT(DT_FASTQ)) {
        bool is_paired_fastq = fastq_piz_is_paired(); // also updates z_file->z_flags in case of backward compatability issues

        // --R1/--R2 and --interleaved is only possible for on FASTQ data compressed with --pair
        ASSINP (!flag.one_component || is_paired_fastq, 
                "--R%c is not supported for %s because it only works on FASTQ data that was compressed with --pair", '0'+flag.one_component, z_name);

        ASSINP (!flag.interleaved || is_paired_fastq, 
                "--R%c is not supported for %s because it only works on FASTQ data that was compressed with --pair", '0'+flag.one_component, z_name);

        // genocat paired FASTQ: if none of --interleaved, --R1 or --R2 are specified, INTERLEAVE_BOTH is the default
        if (exe_type == EXE_GENOCAT && !flag.interleaved && is_paired_fastq && !flag.one_component)
            flag.interleaved = INTERLEAVE_BOTH;

        // genounzip paired FASTQ: always unbind - if user didn't specify prefix, then no prefix
        if (exe_type == EXE_GENOUNZIP && is_paired_fastq && !flag.unbind)
            flag.unbind = ""; 
    }
    else {
        ASSINP  (!flag.one_component, "--R%c is supported only for FASTQ files", '0'+flag.one_component);
        ASSINP0 (!flag.interleaved, "--interleaved is supported only for FASTQ files");
        ASSINP0 (!flag.unbind, "--prefix is supported only for FASTQ files");
    }

    // downsample not possible for Generic or Chain
    ASSINP (!flag.downsample || (flag.out_dt != DT_CHAIN && flag.out_dt != DT_GENERIC), 
            "%s: --downsample is not supported for %s files", z_name, dt_name (flag.out_dt));

    // --sequential only possible on FASTA (including --phylip)
    ASSINP (!flag.sequential || flag.out_dt == DT_FASTA || flag.out_dt == DT_PHYLIP, 
            "--sequential is not supported for %s because it only works on FASTA data, but this file has %s data",
            z_name, dt_name (dt));

    // --sex is only possible on SAM/BAM and FASTQ
    ASSINP (!flag.show_sex || flag.out_dt == DT_BAM || flag.out_dt == DT_SAM || flag.out_dt == DT_FASTQ, // note: if genozip file has BAM data, it will be translated to SAM bc it is always stdout
            "--sex is not supported for %s because it only works on SAM, BAM and FASTQ data, but this file has %s data",
            z_name, dt_name (dt));

    // --coverage is only possible on SAM/BAM and FASTQ
    ASSINP (!flag.show_coverage || flag.out_dt == DT_BAM || flag.out_dt == DT_SAM || flag.out_dt == DT_FASTQ, // note: if genozip file has BAM data, it will be translated to SAM bc it is always stdout
            "--coverage is not supported for %s because it only works on SAM, BAM and FASTQ data, but this file has %s data",
            z_name, dt_name (dt));

    // --idxstats is only possible on SAM/BAM and FASTQ
    ASSINP (!flag.idxstats || flag.out_dt == DT_BAM || flag.out_dt == DT_SAM || flag.out_dt == DT_FASTQ, // note: if genozip file has BAM data, it will be translated to SAM bc it is always stdout
            "--idxstats is not supported for %s because it only works on SAM, BAM and FASTQ data, but this file has %s data",
            z_name, dt_name (dt));

    // --add-lines-numbers is only possible on SAM
    ASSINP0 (!flag.add_line_numbers || flag.out_dt == DT_SAM, "--add_line_numbers works on SAM/BAM data, when outputting it as SAM");

    // --show-chain on CHAIN
    ASSINP (!flag.show_chain || flag.out_dt == DT_CHAIN, "--show-chain is not supported for %s because it only works on Chain files, but this file has %s data", z_name, dt_name (dt));
    
    // --seq-only and --qual-only only work on FASTQ
    ASSINP (!flag.seq_only  || flag.out_dt == DT_FASTQ, "--seq-only is not supported for %s because it only works on FASTQ data, but this file has %s data", z_name, dt_name (dt));
    ASSINP (!flag.qual_only || flag.out_dt == DT_FASTQ, "--qual-only is not supported for %s because it only works on FASTQ data, but this file has %s data", z_name, dt_name (dt));

    // --unbind and --component not allowed in SAM/BAM starting v14
    if (Z_DT(DT_SAM) && z_file->genozip_version >= 14) {
        ASSERT0 (!flag.unbind, "--unbind not supported for SAM/BAM files");
        ASSERT0 (!flag.one_component, "--component not supported for SAM/BAM files");
    }
    
    // BAM limitations
    if (flag.out_dt == DT_BAM && flag.no_header == 1) {
        if (!flag.force) {
            flag.no_header = 0;
            WARN_ONCE0 ("ignoring --no-header: cannot output a BAM without a header. Try using --sam. Use --force to override.");
        }
        else
            WARN_ONCE0 ("Warning: Outputting a BAM without a header. This will not be a valid BAM file. Did you intend to combine with --sam? Use --quiet to suppress this warning.");
    }

    ASSINP0 (flag.out_dt != DT_BAM || z_file_i <= 0 || flag.test, "Cannot concatenate multiple BAM files. Try using --sam."); // TODO: bug 349

    // -- grep doesn't work with binary files
    ASSINP (!flag.grep || !out_dt_is_binary, "--grep is not supported when outputting %s data%s", 
            dt_name (flag.out_dt), flag.out_dt == DT_BAM ? ". Tip: add --sam": "");

    // --gpos requires --reference
    ASSINP0 (!flag.gpos || flag.reference, "--gpos requires --reference");
    
    // if its not a dual coordintes file, --single-coord is ignored
    if (!z_is_dvcf) flag.single_coord = false;

    // translator limitations
    ASSINP0 ((flag.out_dt != DT_SAM && flag.out_dt != DT_BAM) || dt == DT_SAM || dt == DT_BAM,
             "--sam and --bam are only allowed for SAM, BAM or CRAM files");

    ASSINP0 (flag.out_dt != DT_VCF || dt == DT_VCF || dt == DT_ME23,
             "--vcf is only allowed for 23andMe files");

    ASSINP0 (flag.out_dt != DT_PHYLIP || dt == DT_PHYLIP || dt == DT_FASTA,
             "--phylip is only allowed for multi-FASTA files");

    ASSINP0 (flag.out_dt != DT_FASTA || dt == DT_PHYLIP || dt == DT_FASTA,
             "--fasta is only allowed for PHYLIP files");

    ASSINP0 (flag.out_dt != DT_FASTQ || dt == DT_FASTQ || dt == DT_SAM || dt == DT_BAM,
             "--fastq is only allowed for SAM or BAM files");

    // version limitations

    // sam/bam genozip files generated in v9-11 had a critical bug when translating to fastq
    ASSINP (z_file->genozip_version >= 12 || !(dt == DT_SAM && flag.out_dt == DT_FASTQ),
            "%s was created with genozip version %u, SAM/BAM to FASTQ translation is supported only for files created with genozip version 12 or later",
            z_name, z_file->genozip_version);

    // version 12 broke backward compatability of being able to translate old multifasta to phylip
    ASSINP (z_file->genozip_version >= 12 || !(dt == DT_FASTA && flag.out_dt == DT_PHYLIP),
            "%s was created with genozip version %u, MULTIFASTA to PHYLIP translation is supported only for files created with genozip version 12 or later",
            z_name, z_file->genozip_version);

    // num_lines in VbHeader populated since v12 (in v14 moved to SectionEnt)
    ASSINP (z_file->genozip_version >= 12 || (flag.lines_first == NO_LINE && !flag.tail && !flag.downsample),
            "%s was created with genozip version %u, --head, --tail, --lines and --downsample are supported only for files created with genozip version 12 or later",
            z_name, z_file->genozip_version);

    flags_test_conflicts(0); // test again after updating flags

    info_stream = (!flag.to_stdout || flag.no_writer) ? stdout : stderr;
    is_info_stream_terminal = isatty (fileno (info_stream)); 

    if (flag.show_flags) flags_show_flags();
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
            buf_add_string (evb, &command_line, argv[i]);       // can't use bufprintf because argv[i] length is unbound
            if (i < argc-1) BNXTc (command_line) = ' '; // buf_add_string allocs one char extra
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
            buf_add_string (evb, &debugger_params, argv[i]); // can't use bufprintf because argv[i] length is unbound
            bufprintf (evb, &debugger_params, "\"%s", i < argc-1 ? ", ": "],");
        }
    }

    if (command == ZIP && !isatty(0))
        flags_store_piped_in_details();    

    if (flag.echo || getenv ("GENOZIP_TEST") || flag.debug) {
        fprintf (stderr, "\n%s: %s\n", str_time().s, command_line.data);
        flag.echo = 2; // done
    }

    if (getenv ("GENOZIP_TEST") || flag.debug) 
        flags_display_debugger_params();
}

rom flags_command_line (void)
{
    return command_line.data;
}

void flags_display_debugger_params (void)
{
    fprintf (stderr, "%s\n", debugger_params.data);
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

