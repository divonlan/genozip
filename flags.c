// ------------------------------------------------------------------
//   flags.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <getopt.h>
#include "genozip.h"
#include "flags.h"
#include "data_types.h"
#include "file.h"
#include "regions.h"
#include "dict_id.h"
#include "crypt.h"
#include "strings.h"
#include "fastq.h"

Flags flag = { .out_dt = DT_NONE };
bool option_is_short[256] = { }; // indexed by character of short option.
FILE *info_stream;               // either stdout or stderr - where non-error messages should go

static Buffer command_line = EMPTY_BUFFER;

// command line options that get assigned to flags
static int option_noisy=0, option_best=0;

void flags_init_from_command_line (int argc, char **argv)
{
    // process command line options
    while (1) {

        #define _i  {"input",         required_argument, 0, 'i'                    }
        #define _I  {"stdin-size",    required_argument, 0, 'I'                    }
        #define _c  {"stdout",        no_argument,       &flag.to_stdout,        1 }
        #define _d  {"decompress",    no_argument,       &command, PIZ             }
        #define _f  {"force",         no_argument,       &flag.force,            1 }
        #define _h  {"help",          optional_argument, 0, 'h'                    }
        #define _l  {"list",          no_argument,       &command, LIST            }
        #define _L1 {"license",       no_argument,       &command, LICENSE         } // US spelling
        #define _L2 {"licence",       no_argument,       &command, LICENSE         } // British spelling
        #define _q  {"quiet",         no_argument,       &flag.quiet,            1 }
        #define _Q  {"noisy",         no_argument,       &option_noisy,          1 }
        #define _DL {"replace",       no_argument,       &flag.replace,          1 }
        #define _V  {"version",       no_argument,       &command, VERSION         }
        #define _z  {"bgzf",          no_argument,       &flag.bgzf,             1 }
        #define _z0 {"plain",         no_argument,       &flag.plain,            1 }
        #define _zb {"bam",           no_argument,       &flag.out_dt,           DT_BAM }
        #define _zB {"BAM",           no_argument,       &flag.out_dt,           DT_BAM }
        #define _zs {"sam",           no_argument,       &flag.out_dt,           DT_SAM }
        #define _zS {"SAM",           no_argument,       &flag.out_dt,           DT_SAM }
        #define _zq {"fastq",         no_argument,       &flag.out_dt,           DT_FASTQ }
        #define _zQ {"FASTQ",         no_argument,       &flag.out_dt,           DT_FASTQ }
        #define _zf {"fq",            no_argument,       &flag.out_dt,           DT_FASTQ }
        #define _zF {"FQ",            no_argument,       &flag.out_dt,           DT_FASTQ }
        #define _za {"fasta",         no_argument,       &flag.out_dt,           DT_FASTA }
        #define _zA {"FASTA",         no_argument,       &flag.out_dt,           DT_FASTA }
        #define _zc {"bcf",           no_argument,       &flag.out_dt,           DT_BCF }
        #define _zC {"BCF",           no_argument,       &flag.out_dt,           DT_BCF }
        #define _zv {"vcf",           no_argument,       &flag.out_dt,           DT_VCF }
        #define _zV {"VCF",           no_argument,       &flag.out_dt,           DT_VCF }
        #define _zy {"phylip",        no_argument,       &flag.out_dt,           DT_PHYLIP }
        #define _zY {"Phylip",        no_argument,       &flag.out_dt,           DT_PHYLIP }
        #define _m  {"md5",           no_argument,       &flag.md5,              1 }
        #define _t  {"test",          no_argument,       &flag.test,             1 }
        #define _fa {"fast",          no_argument,       &flag.fast,             1 }
        #define _bs {"best",          no_argument,       &option_best,           1 }
        #define _9  {"optimize",      no_argument,       &flag.optimize,         1 } // US spelling
        #define _99 {"optimise",      no_argument,       &flag.optimize,         1 } // British spelling
        #define _9s {"optimize-sort", no_argument,       &flag.optimize_sort,    1 }
        #define _9P {"optimize-PL",   no_argument,       &flag.optimize_PL,      1 }
        #define _9G {"optimize-GL",   no_argument,       &flag.optimize_GL,      1 }
        #define _9g {"optimize-GP",   no_argument,       &flag.optimize_GP,      1 }
        #define _9V {"optimize-VQSLOD", no_argument,     &flag.optimize_VQSLOD,  1 }
        #define _9Q {"optimize-QUAL", no_argument,       &flag.optimize_QUAL,    1 } 
        #define _9f {"optimize-Vf",   no_argument,       &flag.optimize_Vf,      1 }
        #define _9Z {"optimize-ZM",   no_argument,       &flag.optimize_ZM,      1 }
        #define _9D {"optimize-DESC", no_argument,       &flag.optimize_DESC,    1 }
        #define _gt {"gtshark",       no_argument,       &flag.gtshark,          1 } 
        #define _pe {"pair",          no_argument,       &flag.pair,   PAIR_READ_1 } 
        #define _th {"threads",       required_argument, 0, '@'                    }
        #define _u  {"unbind",        optional_argument, 0, 'u'                    }
        #define _o  {"output",        required_argument, 0, 'o'                    }
        #define _p  {"password",      required_argument, 0, 'p'                    }
        #define _B  {"vblock",        required_argument, 0, 'B'                    }
        #define _r  {"regions",       required_argument, 0, 'r'                    }
        #define _s  {"samples",       required_argument, 0, 's'                    }
        #define _il {"interleave",    no_argument,       &flag.interleave,       1 }
        #define _e  {"reference",     required_argument, 0, 'e'                    }
        #define _E  {"REFERENCE",     required_argument, 0, 'E'                    }
        #define _b  {"bytes",         no_argument,       &flag.bytes,            1 }
        #define _me {"make-reference",no_argument,       &flag.make_reference,   1 }
        #define _x  {"index",         no_argument,       &flag.index_txt,        1 }
        #define _g  {"grep",          required_argument, 0, 'g'                    }
        #define _G  {"drop-genotypes",no_argument,       &flag.drop_genotypes,   1 }
        #define _H1 {"no-header",     no_argument,       &flag.no_header,        1 }
        #define _H0 {"header-only",   no_argument,       &flag.header_only,      1 }
        #define _1  {"header-one",    no_argument,       &flag.header_one,       1 }
        #define _GT {"GT-only",       no_argument,       &flag.gt_only,          1 }
        #define _Gt {"gt-only",       no_argument,       &flag.gt_only,          1 }
        #define _ds {"downsample",    required_argument, 0, 9                      }
        #define _PG {"no-PG",         no_argument,       &flag.no_pg,            1 }
        #define _pg {"no-pg",         no_argument,       &flag.no_pg,            1 }
        #define _fs {"sequential",    no_argument,       &flag.sequential,       1 }  
        #define _rg {"register",      no_argument,       &flag.do_register,      1 }
        #define _ss {"show-stats",    no_argument,       &flag.show_stats,       1 } 
        #define _SS {"SHOW-STATS",    no_argument,       &flag.show_stats,       2 } 
        #define _lc {"list-chroms",   no_argument,       &flag.list_chroms,      1 } // identical to --show-dict=CHROM 
        #define _lC {"list-contigs",  no_argument,       &flag.list_chroms,      1 } 
        #define _s2 {"show-b250",     optional_argument, 0, 2,                     }
        #define _sd {"show-dict",     optional_argument, 0, 3                      }
        #define _s7 {"dump-b250",     required_argument, 0, 5                      }
        #define _S7 {"dump-local",    required_argument, 0, 6                      }
        #define _S9 {"dump-section",  required_argument, 0, 7                      }        
        #define _sa {"show-alleles",  no_argument,       &flag.show_alleles,     1 }
        #define _st {"show-time",     optional_argument, 0, 1                      } 
        #define _sm {"show-memory",   no_argument,       &flag.show_memory ,     1 } 
        #define _sh {"show-headers",  optional_argument, 0, 10                     } 
        #define _si {"show-index",    no_argument,       &flag.show_index  ,     1 } 
        #define _Si {"show-ref-index",no_argument,       &flag.show_ref_index,   1 } 
        #define _Sh {"show-ref-hash" ,no_argument,       &flag.show_ref_hash,    1 } 
        #define _sr {"show-gheader",  no_argument,       &flag.show_gheader,     1 }  
        #define _sT {"show-threads",  no_argument,       &flag.show_threads,     1 }  
        #define _sv {"show-vblocks",  no_argument,       &flag.show_vblocks,     1 }  
        #define _ov {"one-vb",        required_argument, 0, 8                      }  
        #define _sR {"show-reference",no_argument,       &flag.show_reference,   1 }  
        #define _sC {"show-ref-contigs", no_argument,    &flag.show_ref_contigs, 1 }  
        #define _rA {"show-ref-alts", no_argument,       &flag.show_ref_alts,    1 }  
        #define _rS {"show-ref-seq",  no_argument,       &flag.show_ref_seq,     1 }  
        #define _cn {"show-containers", no_argument,     &flag.show_containers,  1 }  
        #define _hC {"show-txt-contigs", no_argument,    &flag.show_txt_contigs, 1 }
        #define _sI {"show-is-set",   required_argument, 0, '~',                   }  
        #define _sA {"show-aliases",  no_argument,       &flag.show_aliases,     1 }  
        #define _sc {"show-codec",    no_argument,       &flag.show_codec,       1 }  
        #define _sb {"show-bgzf",     no_argument,       &flag.show_bgzf,        1 }
        #define _s5 {"show-digest",   no_argument,       &flag.show_digest,      1 }
        #define _sM {"show-mutex",    optional_argument, 0, 4                      }
        #define _dS {"test-seg",      no_argument,       &flag.test_seg,         1 }  
        #define _xt {"xthreads",      no_argument,       &flag.xthreads,         1 }  
        #define _dm {"debug-memory",  no_argument,       &flag.debug_memory,     1 }  
        #define _dp {"debug-progress",no_argument,       &flag.debug_progress,   1 }  
        #define _dh {"show-hash",     no_argument,       &flag.show_hash,        1 }  
        #define _00 {0, 0, 0, 0                                                    }

        typedef const struct option Option;
        static Option genozip_lo[]    = { _i, _I, _c, _d, _f, _h,    _l, _L1, _L2, _q, _Q, _t, _DL, _V, _z, _z0, _zb, _zB, _zs, _zS, _zq, _zQ, _za, _zA, _zf, _zF, _zc, _zC, _zv, _zV, _zy, _zY, _m, _th, _u, _o, _p, _e, _E,                                          _ss, _SS, _sd, _sT, _sb, _lc, _lC, _s2, _s7, _S7, _S9, _sa, _st, _sm, _sh, _si, _Si, _Sh, _sr, _sv,     _B, _xt, _dm, _dp,      _dh,_dS, _9, _99, _9s, _9P, _9G, _9g, _9V, _9Q, _9f, _9Z, _9D, _pe, _fa, _bs,              _rg, _sR,      _sC, _hC, _rA, _rS, _me,      _s5, _sM, _sA, _sc, _sI, _gt,                    _00 };
        static Option genounzip_lo[]  = {         _c,     _f, _h, _x,    _L1, _L2, _q, _Q, _t, _DL, _V, _z, _z0, _zb, _zB, _zs, _zS, _zq, _zQ, _za, _zA, _zf, _zF, _zc, _zC, _zv, _zV, _zy, _zY, _m, _th, _u, _o, _p, _e,                                              _ss, _SS, _sd, _sT, _sb, _lc, _lC, _s2, _s7, _S7, _S9, _sa, _st, _sm, _sh, _si, _Si, _Sh, _sr, _sv,         _xt, _dm, _dp,                                                                                                      _sR, _sC, _hC, _rA, _rS,           _s5, _sM, _sA,      _sI,      _cn, _pg, _PG,     _00 };
        static Option genocat_lo[]    = {         _c,     _f, _h, _x,    _L1, _L2, _q, _Q,          _V,     _z0, _zb, _zB, _zs, _zS, _zq, _zQ, _za, _zA, _zf, _zF, _zc, _zC, _zv, _zV, _zy, _zY,     _th,     _o, _p,         _il, _r, _s, _G, _1, _H0, _H1, _Gt, _GT, _ss, _SS, _sd, _sT, _sb, _lc, _lC, _s2, _s7, _S7, _S9, _sa, _st, _sm, _sh, _si, _Si, _Sh, _sr, _sv, _ov,    _xt, _dm, _dp, _ds,                                                                                   _fs, _g,      _sR, _sC, _hC, _rA, _rS,           _s5, _sM, _sA,      _sI,      _cn, _pg, _PG,     _00 };
        static Option genols_lo[]     = {                 _f, _h,        _L1, _L2, _q,              _V,                                                                                                   _u,     _p, _e,                                                                                                          _st, _sm,                                       _dm,                                                                                                                                                                  _sM,                                    _b, _00 };
        static Option *long_options[] = { genozip_lo, genounzip_lo, genols_lo, genocat_lo }; // same order as ExeType

        // include the option letter here for the short version (eg "-t") to work. ':' indicates an argument.
        static const char *short_options[NUM_EXE_TYPES] = { // same order as ExeType
            "i:I:cdfhlLqQt^Vzm@:o:p:B:9wWFe:E:2zu", // genozip (note: includes some genounzip options to be used in combination with -d)
            "czfhLqQt^V@:uo:p:me:wWx",              // genounzip
            "hLVp:qfub",                            // genols
            "hLV@:p:qQ1r:s:H1Go:fg:e:E:wWx"         // genocat
        };

        int option_index = -1;
        int c = getopt_long (argc, argv, short_options[exe_type], long_options[exe_type], &option_index);

        if (c == -1) break; // no more options

        option_is_short[c] = (option_index == -1); // true if user provided a short option - eg -c rather than --stdout

        switch (c) {
            case HELP:
                flag.help = optarg; // fall through

            case PIZ : case LIST : case LICENSE : 
            case VERSION :
                ASSINP (command<0 || command==c, "%s: can't have both -%c and -%c", global_cmd, command, c); 
                command=c; 
                break;

            case 'i' : file_set_input_type (optarg) ; break;
            case 'I' : file_set_input_size (optarg) ; break;
            case 'c' : flag.to_stdout     = 1       ; break;
            case 'F' : flag.fast          = 1       ; break;
            case 'z' : flag.bgzf          = 1       ; break;
            case 'f' : flag.force         = 1       ; break;
            case '^' : flag.replace       = 1       ; break;
            case 'q' : flag.quiet         = 1       ; break;
            case 'Q' : option_noisy       = 1       ; break;
            case '9' : flag.optimize      = 1       ; break;
            case 'w' : flag.show_stats    = 1       ; break;
            case 'W' : flag.show_stats    = 2       ; break;
            case '2' : flag.pair    = PAIR_READ_1   ; break;
            case 't' : flag.test          = 1       ; break; 
            case 'x' : flag.index_txt     = 1       ; break;
            case 'b' : flag.bytes         = 1       ; break;         
            case 'r' : flag.regions       = 1       ; regions_add     (optarg); break;
            case 's' : flag.samples       = 1       ; vcf_samples_add (optarg); break;
            case 'e' : flag.reference     = REF_EXTERNAL  ; ref_set_reference (optarg); break;
            case 'E' : flag.reference     = REF_EXT_STORE ; ref_set_reference (optarg); break;
            case 'm' : flag.md5           = 1       ; break;
            case 'u' : flag.unbind    = optarg ? optarg : "" ; break; // unbind with or without a prefix
            case '\1': flag.show_time = optarg ? optarg : "" ; break; // show_time with or without a specific member of ProfilerRec
            case 'G' : flag.drop_genotypes= 1       ; break;
            case 'H' : flag.no_header     = 1       ; break;
            case '1' : flag.header_one    = 1       ; break;
            case '@' : flag.threads_str   = optarg  ; break;
            case 'o' : flag.out_filename  = optarg  ; break;
            case 'g' : flag.grep          = optarg  ; break;
            case '~' : flag.show_is_set   = optarg  ; break;
            case 7   : flag.dump_section  = optarg  ; break;
            case 4   : flag.show_mutex    = optarg ? optarg : (char*)1; break;
            case 2   : if (optarg) flag.dict_id_show_one_b250  = dict_id_make (optarg, strlen (optarg)); 
                       else        flag.show_b250 = 1;
                       break;
            case 3   : if (optarg) flag.dict_id_show_one_dict  = dict_id_make (optarg, strlen (optarg)); 
                       else        flag.show_dict = 1;
                       break;
            case 5   : flag.dump_one_b250_dict_id  = dict_id_make (optarg, strlen (optarg)); break;
            case 6   : flag.dump_one_local_dict_id = dict_id_make (optarg, strlen (optarg)); break;
            case 8   : flag.one_vb = atoi (optarg);  break;
            case 9   : flag.downsample = atoi (optarg); break;
            case 10  : flag.show_headers = 1 + sections_st_by_name (optarg); break; // +1 so SEC_NONE maps to 0
            case 'B' : vb_set_global_max_memory_per_vb (optarg); 
                       flag.vblock = true;
                       break;
            case 'p' : crypt_set_password (optarg) ; break;

            case 0   : // a long option - already handled; except for 'o' and '@'

                if (long_options[exe_type][option_index].val == 'o') 
                    flag.out_filename = optarg;

                if (long_options[exe_type][option_index].val == 'p') 
                    crypt_set_password (optarg);

                else if (long_options[exe_type][option_index].val == '@') 
                    flag.threads_str = optarg;

                else 
                    ASSINP (long_options[exe_type][option_index].flag != &command || 
                            long_options[exe_type][option_index].val  == command ||
                            command < 0, 
                            "%s: can't have both --%s and -%c", global_cmd, long_options[exe_type][option_index].name, command);
                break; 

            case '?' : // unrecognized option - error message already displayed by libc
            default  :
                fprintf (stderr, "Usage: %s [OPTIONS] filename1 filename2...\nTry %s --help for more information.\n", global_cmd, global_cmd);
                exit(1);  
        }
    }
}

static void flags_warn_if_duplicates (int num_txt_files, const char **filenames)
{
    if (num_txt_files <= 1) return; // nothing to do

    # define BASENAME_LEN 256
    char basenames[num_txt_files * BASENAME_LEN];

    for (unsigned i=0; i < num_txt_files; i++)
        file_basename (filenames[i], false, "", &basenames[i*BASENAME_LEN], BASENAME_LEN);

    qsort (basenames, num_txt_files, BASENAME_LEN, (int (*)(const void *, const void *))strcmp);

    for (unsigned i=1; i < num_txt_files; i++) 
        ASSERTW (strncmp(&basenames[(i-1) * BASENAME_LEN], &basenames[i * BASENAME_LEN], BASENAME_LEN), 
                 "Warning: two files with the same name '%s' - if you later unbind with 'genounzip --unbind %s', these files will overwrite each other", 
                 &basenames[i * BASENAME_LEN], flag.out_filename);
}

static void flags_test_conflicts (void)
{
    ASSINP (!flag.to_stdout   || !flag.out_filename,                "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("output", "o"));
    ASSINP (!flag.to_stdout   || !flag.unbind,                      "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("unbind", "u"));
    ASSINP (!flag.unbind      || !flag.out_filename,                "%s: option %s is incompatable with %s", global_cmd, OT("unbind",  "u"), OT("output", "o"));
    ASSINP (!flag.to_stdout   || !flag.replace,                     "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("replace", "^"));

    // if flag.md5 we need to seek back and update the md5 in the txt header section - this is not possible with flag.to_stdout
    ASSINP (!flag.to_stdout   || !flag.md5,                         "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("md5", "m"));
    ASSINP (!flag.to_stdout   || !flag.index_txt,                   "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("index", "x"));
    ASSINP (!flag.test        || !flag.out_filename || command != PIZ, "%s: option %s is incompatable with %s", global_cmd, OT("output", "o"),  OT("test", "t"));
    ASSINP (!flag.test        || !flag.replace || command != PIZ,   "%s: option %s is incompatable with %s", global_cmd, OT("replace", "^"), OT("test", "t"));
    ASSINP (!flag.to_stdout   || command != ZIP,                    "%s: option %s only works for decompressing files, not compressing", global_cmd, OT("stdout", "c"));
    ASSINP (!flag.one_vb      || !flag.interleave,                  "%s: option %s is incompatable with %s", global_cmd, "interleave", "one-vb");
    ASSINP (!flag.xthreads    || !flag.interleave,                  "%s: option %s is incompatable with %s", global_cmd, "interleave", "xthreads");
    ASSINP (!flag.header_only || !flag.interleave,                  "%s: option %s is incompatable with %s", global_cmd, "interleave", "header-only");
    ASSINP (!flag.header_only || !flag.no_header,                   "%s: option %s is incompatable with %s", global_cmd, OT("no-header", "H"), "header-only");
    ASSINP (!flag.no_header   || !flag.header_one,                  "%s: option %s is incompatable with %s", global_cmd, OT("no-header", "H"), OT("header-one", "1"));
    ASSINP (!flag.quiet       || !option_noisy,                     "%s: option %s is incompatable with %s", global_cmd, OT("quiet", "q"), OT("noisy", "Q"));
    ASSINP (!flag.test        || !flag.optimize,                    "%s: option %s is incompatable with %s", global_cmd, OT("test", "t"), OT("optimize", "9"));
    ASSINP (!flag.md5         || !flag.optimize,                    "%s: option %s is incompatable with %s", global_cmd, OT("md5", "m"), OT("optimize", "9"));
    ASSINP (!flag.samples     || !flag.drop_genotypes,              "%s: option %s is incompatable with %s", global_cmd, OT("samples", "s"), OT("drop-genotypes", "G"));
    ASSINP (!option_best      || !flag.fast,                        "%s: option %s is incompatable with --best", global_cmd, OT("fast", "t"));
    ASSINP (!flag.plain       || !flag.bgzf,                        "%s: option %s is incompatable with --plain", global_cmd, OT("bgzf", "z"));
    ASSINP (!flag.test        || !flag.make_reference,              "%s: option %s is incompatable with --make-reference", global_cmd, OT("test", "t"));
    ASSINP (flag.reference != REF_EXTERNAL  || !flag.make_reference,"%s: option %s is incompatable with --make-reference", global_cmd, OT("reference", "e"));
    ASSINP (flag.reference != REF_EXT_STORE || !flag.make_reference,"%s: option %s is incompatable with --make-reference", global_cmd, OT("REFERENCE", "E"));
    ASSINP (flag.reference != REF_EXT_STORE || exe_type != EXE_GENOCAT, "%s: option %s supported only for viewing the reference file itself", global_cmd, OT("REFERENCE", "E"));
    ASSINP (flag.reference != REF_EXTERNAL  || !flag.show_ref_seq,  "%s: option %s is incompatable with --show-ref-seq: use genocat --show-ref-seq on the reference file itself instead", global_cmd, OT("reference", "e"));
    ASSINP (!flag.dump_one_b250_dict_id.num || !flag.dump_one_local_dict_id.num, "%s: option --dump-one-b250 is incompatable with --dump-one-local", global_cmd);

    // some genozip flags are allowed only in combination with --decompress 
    if (exe_type == EXE_GENOZIP && command == ZIP) {
        char s[20]; 
        ASSINP (!flag.bgzf,        "%s: option %s can only be used if --decompress is used too", global_cmd, OT("bgzf", "z"));
        ASSINP (flag.out_dt == DT_NONE, "%s: option --%s can only be used if --decompress is used too", global_cmd, str_tolower (dt_name (flag.out_dt), s));
        ASSINP (!flag.unbind,       "%s: option %s can only be used if --decompress is used too", global_cmd, OT("unbind", "u"));
        ASSINP (!flag.show_aliases, "%s: option --show_aliases can only be used if --decompress is used too", global_cmd);
        ASSINP (!flag.show_is_set,  "%s: option --show_is_set can only be used if --decompress is used too", global_cmd);
    }
}

// --pair: verify an even number of fastq files, --output, and --reference/--REFERENCE
static void flags_verify_pair_rules (unsigned num_txt_files, const char **filenames)
{
    // ZIP only
    if (command != ZIP) {
        flag.pair = false;
        return;
    }

    // verify even number of files
    ASSINP (num_txt_files % 2 == 0, "%s: when using %s, expecting an even number of FASTQ input files, each consecutive two being a pair", global_cmd, OT("pair", "2"));
    ASSINP (flag.reference,     "%s: either --reference or --REFERENCE must be specified when using %s", global_cmd, OT("pair", "2"));

    // verify all are fastq
    for (unsigned i=0; i < num_txt_files; i++)
        ASSERT (txtfile_get_file_dt (filenames[i]) == DT_FASTQ, "%s: when using %s, all input files are expected to be FASTQ files, but %s is not", global_cmd, OT("pair", "2"), filenames[i]);

    // if which --output is missing, we check if every pair of files has a consistent name
    if (!flag.out_filename) 
        for (unsigned i=0; i < num_txt_files; i += 2) 
            ASSINP (file_get_fastq_pair_filename (filenames[i], filenames[i+1], true),  
                    "%s: to use %s without specifying --output, the naming of the files needs to be consistent and include the numbers 1 and 2 respectively, but these files don't: %s %s", 
                    global_cmd, OT("pair", "2"), filenames[i], filenames[i+1]);
}

void flags_update (unsigned num_txt_files, const char **filenames)
{
    flags_test_conflicts();

    // verify stuff needed for --pair
    if (flag.pair) flags_verify_pair_rules (num_txt_files, filenames);
    
    // don't show progress or warning when outputing to stdout (note: we are "quiet" even if output doesn't go to the terminal
    // because often it will be piped and ultimately go the terminal)
    if (flag.to_stdout) flag.quiet=true; 
    
    // don't show progress for flags that output throughout the process. no issue with flags that output only in the end
    if (flag.show_dict || flag.show_b250 || flag.show_headers || flag.show_threads || flag.show_bgzf || flag.show_mutex ||
        flag.dict_id_show_one_b250.num || flag.dict_id_show_one_dict.num || flag.show_reference || flag.show_digest ||
        flag.show_alleles || flag.show_vblocks || flag.show_codec || (flag.show_index && command==PIZ))
        flag.quiet=true; // don't show progress

    // override these ^ if user chose to be --noisy
    if (option_noisy) flag.quiet=false;

    if (command == ZIP && flag.test) flag.md5=true; // test implies md5

    // default values, if not overridden by the user
    if (!flag.vblock) vb_set_global_max_memory_per_vb (flag.fast ? TXT_DATA_PER_VB_FAST : TXT_DATA_PER_VB_DEFAULT); 

    // --make-reference implies --md5 --B1 (unless --vblock says otherwise), and not encrypted. 
    // in addition, txtfile_read_vblock() limits each VB to have exactly one contig.
    if (flag.make_reference) {
        ASSINP (!crypt_have_password(), "%s: option --make-reference is incompatable with %s", global_cmd, OT("password", "p"));

        flag.md5 = true;
        if (!flag.vblock) vb_set_global_max_memory_per_vb("1");

        ASSINP (num_txt_files <= 1, "%s: you can specify only one FASTA file when using --make-reference.\n"
                "To create a reference from multiple FASTAs use something like this:\n"
                "cat *.fa | %s --make-reference --input fasta --output myref.ref.genozip -", global_cmd, global_cmd);
    }

    // if --optimize was selected, all optimizations are turned on
    if (flag.optimize)
        flag.optimize_sort = flag.optimize_PL = flag.optimize_GL = flag.optimize_GP   = flag.optimize_VQSLOD = 
        flag.optimize_QUAL = flag.optimize_Vf = flag.optimize_ZM = flag.optimize_DESC = true;
    
    // if any optimization flag is on, we turn on flag.optimize
    if (flag.optimize_sort || flag.optimize_PL || flag.optimize_GL || flag.optimize_GP   || flag.optimize_VQSLOD ||
        flag.optimize_QUAL || flag.optimize_Vf || flag.optimize_ZM || flag.optimize_DESC)
        flag.optimize = true;

    // if using the -o option - check that we don't have duplicate filenames (even in different directory) as they
    // will overwrite each other if extracted with --unbind
    if (command == ZIP && flag.out_filename && !flag.quiet) flags_warn_if_duplicates (num_txt_files, filenames);

    flag.multiple_files = (num_txt_files > 1);

    if (command == ZIP && num_txt_files > 1) 
        flag.bind = flag.out_filename ? BIND_ALL   : 
                    flag.pair         ? BIND_PAIRS : 
                                        BIND_NONE  ;

    // cases where genocat is used to view some information, but not the file contents
    flag.genocat_info_only = exe_type == EXE_GENOCAT &&
        (flag.show_stats || flag.show_dict || flag.show_b250 || flag.list_chroms || flag.dict_id_show_one_dict.num ||
        flag.show_index || flag.dump_one_local_dict_id.num || flag.dump_one_b250_dict_id.num || flag.dump_section || flag.show_headers ||
        flag.show_reference || flag.show_ref_contigs || flag.show_ref_index || flag.show_ref_hash || flag.show_ref_alts || 
        flag.show_ref_seq || flag.show_aliases || flag.show_txt_contigs);

    // where progress, metadata etc messages should go. data always goes to stdout and errors/warning always go to stderr.
    info_stream = (!flag.to_stdout || flag.genocat_info_only) ? stdout : stderr;
}

void flags_update_zip_one_file (void)
{
}

// PIZ: called after opening z_file and reading the header before opening txt_file
void flags_update_piz_one_file (void)
{
    if (flag.out_dt == DT_NONE) {

        // handle native binary formats (BAM). note on BCF and CRAM: we used bcftools/samtools as an external 
        // compressor, so that genozip sees the text, not binary, data of these files - the same as if the file were compressed with eg bz2
        if (z_file->z_flags.txt_is_bin) {
            
            // PIZ of a genozip file with is_binary (e.g. BAM) is determined here unless the user overrides with --sam or --fastq
            if (exe_type == EXE_GENOCAT) {
                flag.out_dt = z_file->data_type; // output in textual format
                flag.no_pg  = true; // if a user is genocatting a BAM file, we don't add @PG (only when he genounzips it) (flag is ignored if not BAM, so no harm)
            }
            else
                flag.out_dt = DTPZ (bin_type);   // output in binary format
        }
        else
            flag.out_dt = z_file->data_type;
    }

    // .bcf will be bgzipped by bcftools, ignore --bgzf flag as we don't need an additional bgzf step
    if (flag.out_dt == DT_BCF) flag.bgzf=0;

    // phylip implied sequential and header_one
    if (z_file->data_type == DT_FASTA && flag.out_dt == DT_PHYLIP) {
        flag.sequential = 1;
        flag.header_one = flag.no_header = flag.header_only = 0;
    }

    // for FASTA and FASTQ we convert a "header_only" flag to "header_only_fast" as flag.header_only has some additional logic
    // that doesn't work for FASTA / FASTQ
    if (flag.header_only && (flag.out_dt == DT_FASTA || flag.out_dt == DT_FASTQ)) {
        flag.header_only = false;
        flag.header_only_fast  = true;
    }

    // Case were we set flag.bgzf even if the user didn't explicitly ask for it with --bgzf
    // note: another case not cover here is when the output file type is GZ or BGZ or BAM - handled in file_open_txt_write
    if (!flag.plain && // user didn't not explicitly tell us to refrain from BGZF
        (z_file->z_flags.bgzf && exe_type != EXE_GENOCAT) && // source file had BGZF 
        (!flag.to_stdout || (flag.show_headers && exe_type == EXE_GENOCAT))) // we are outputing to a disk file OR user asked to see the headers (so we want to show her BGZF too)
        flag.bgzf = true;   

    // if interleaving bgzf is always false
    if (flag.interleave)
        flag.bgzf = false;

    // Note on BAM/SAM: BAM is stored as binary SAM, so trans_containers=true for BAM->BAM , but false for BAM->SAM
    flag.trans_containers = dt_get_translation().trans_containers; 

    // Check if the reconstructed data type is the same as the source data type
    bool is_binary = z_file->z_flags.txt_is_bin;
    flag.reconstruct_as_src = (flag.out_dt == DT_SAM            && z_file->data_type==DT_SAM && !is_binary) || 
                              (flag.out_dt == DT_BAM            && z_file->data_type==DT_SAM && is_binary ) ||
                              (flag.out_dt == z_file->data_type && z_file->data_type!=DT_SAM);

    // true if the output file of genounzip or genocat will NOT be identical to the source file as recorded in z_file
    // note: this does not account for changes to the data done at the compression stage with --optimize
    flag.data_modified = !flag.reconstruct_as_src || // translating to another data
                         flag.header_one || flag.no_header || flag.header_only || flag.header_only_fast || flag.grep || 
                         flag.regions || flag.samples || flag.drop_genotypes || flag.gt_only || flag.sequential || 
                         flag.one_vb || flag.downsample || flag.interleave;

    // interleaving is only possible for on FASTQ data compressed with --pair
    ASSINP (!flag.interleave || (flag.out_dt == DT_FASTQ && fastq_piz_is_paired()), // untranslated FASTQ that is paired
            "%s: Error: --interleave is not supported for %s because it only works on FASTQ data that was compressed with --pair", 
            global_cmd, z_name);

    flags_test_conflicts(); // test again after updating flags
}

void flags_store_command_line (int argc, char **argv)
{
    unsigned len=0, pw_len=0;
    const char *pw=0;

    for (int i=0; i < argc; i++)
        len += strlen (argv[i]) + 1; // +1 for seperator (' ' or '\0')

    if ((pw = crypt_get_password())) pw_len  = strlen (pw);

    for (int i=0; i < argc; i++) {

        unsigned arg_len = strlen (argv[i]);
        ASSINP (arg_len < BUFPRINTF_MAX_LEN-10, "%s: argument %u longer than maximum allowed %u characters: %s", 
                global_cmd, i, BUFPRINTF_MAX_LEN-10, argv[i]);

        if (pw && !strcmp (argv[i], pw)) // "-p 123", "--pass 123" etc
            bufprintf (evb, &command_line, "***%s", (i < argc-1 ? " ": "")) // hide password

        else if (pw && (arg_len >= pw_len + 2) &&  // check for -p123 or eg -fmp123
                    !strcmp (&argv[i][arg_len-pw_len], pw) && // not air-tight test, but good enough (eg "-ofilenamep123" will incorrectly trigger)
                    argv[i][0] == '-' &&
                    argv[i][arg_len-pw_len-1] == 'p')
            bufprintf (evb, &command_line, "%.*s***%s", arg_len-pw_len, argv[i], (i < argc-1 ? " ": "")) // hide password

        else
            bufprintf (evb, &command_line, "%s%s", argv[i], (i < argc-1 ? " ": ""))
    }
}

const BufferP flags_command_line (void)
{
    return &command_line;
}
