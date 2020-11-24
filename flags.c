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

Flags flag = { .out_dt = DT_NONE };

// command line options that get assigned to flags
static int option_noisy=0, option_best=0;

void flags_init_from_command_line (int argc, char **argv, bool *is_short)
{
    // process command line options
    while (1) {

        #define _i  {"input-type",    required_argument, 0, 'i'                    }
        #define _I  {"stdin-size",    required_argument, 0, 'I'                    }
        #define _c  {"stdout",        no_argument,       &flag.to_stdout,           1 }
        #define _d  {"decompress",    no_argument,       &command, PIZ             }
        #define _f  {"force",         no_argument,       &flag.force,            1 }
        #define _h  {"help",          no_argument,       &command, HELP            }
        #define _l  {"list",          no_argument,       &command, LIST            }
        #define _L1 {"license",       no_argument,       &command, LICENSE         } // US spelling
        #define _L2 {"licence",       no_argument,       &command, LICENSE         } // British spelling
        #define _q  {"quiet",         no_argument,       &flag.quiet,            1 }
        #define _Q  {"noisy",         no_argument,       &option_noisy,          1 }
        #define _DL {"replace",       no_argument,       &flag.replace,          1 }
        #define _V  {"version",       no_argument,       &command, VERSION         }
        #define _z  {"bgzip",         no_argument,       &flag.bgzf,             1 }
        #define _z0 {"plain",         no_argument,       &flag.plain,            1 }
        #define _zb {"bam",           no_argument,       &flag.out_dt,           DT_BAM }
        #define _zB {"BAM",           no_argument,       &flag.out_dt,           DT_BAM }
        #define _zs {"sam",           no_argument,       &flag.out_dt,           DT_SAM }
        #define _zS {"SAM",           no_argument,       &flag.out_dt,           DT_SAM }
        #define _zq {"fastq",         no_argument,       &flag.out_dt,           DT_FASTQ }
        #define _zQ {"FASTQ",         no_argument,       &flag.out_dt,           DT_FASTQ }
        #define _zf {"fq",            no_argument,       &flag.out_dt,           DT_FASTQ }
        #define _zF {"FQ",            no_argument,       &flag.out_dt,           DT_FASTQ }
        #define _zc {"bcf",           no_argument,       &flag.out_dt,           DT_BCF }
        #define _zC {"BCF",           no_argument,       &flag.out_dt,           DT_BCF }
        #define _zv {"vcf",           no_argument,       &flag.out_dt,           DT_VCF }
        #define _zV {"VCF",           no_argument,       &flag.out_dt,           DT_VCF }
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
        #define _e  {"reference",     required_argument, 0, 'e'                    }
        #define _E  {"REFERENCE",     required_argument, 0, 'E'                    }
        #define _b  {"bytes",         no_argument,       &flag.bytes,            1 }
        #define _me {"make-reference",no_argument,       &flag.make_reference,   1 }
        #define _g  {"grep",          required_argument, 0, 'g'                    }
        #define _G  {"drop-genotypes",no_argument,       &flag.drop_genotypes,   1 }
        #define _H1 {"no-header",     no_argument,       &flag.no_header,        1 }
        #define _H0 {"header-only",   no_argument,       &flag.header_only,      1 }
        #define _1  {"header-one",    no_argument,       &flag.header_one,       1 }
        #define _GT {"GT-only",       no_argument,       &flag.gt_only,          1 }
        #define _Gt {"gt-only",       no_argument,       &flag.gt_only,          1 }
        #define _fs {"sequential",    no_argument,       &flag.sequential,       1 }  
        #define _rg {"register",      no_argument,       &flag.do_register,      1 }
        #define _ss {"show-stats",    no_argument,       &flag.show_stats,       1 } 
        #define _SS {"SHOW-STATS",    no_argument,       &flag.show_stats,       2 } 
        #define _lc {"list-chroms",   no_argument,       &flag.list_chroms,      1 } // identical to --show-dict=CHROM 
        #define _lC {"list-contigs",  no_argument,       &flag.list_chroms,      1 } 
        #define _s2 {"show-b250",     optional_argument, 0, '\2'                   }
        #define _sd {"show-dict",     optional_argument, 0, '\3'                   }
        #define _s7 {"dump-b250",     required_argument, 0, '\5'                   }
        #define _S7 {"dump-local",    required_argument, 0, '\6'                   }
        #define _S9 {"dump-section",  required_argument, 0, '\7'                   }        
        #define _sa {"show-alleles",  no_argument,       &flag.show_alleles,     1 }
        #define _st {"show-time",     optional_argument, 0, '\1'                   } 
        #define _sm {"show-memory",   no_argument,       &flag.show_memory ,     1 } 
        #define _sh {"show-headers",  no_argument,       &flag.show_headers,     1 } 
        #define _si {"show-index",    no_argument,       &flag.show_index  ,     1 } 
        #define _Si {"show-ref-index",no_argument,       &flag.show_ref_index,   1 } 
        #define _Sh {"show-ref-hash" ,no_argument,       &flag.show_ref_hash,    1 } 
        #define _sr {"show-gheader",  no_argument,       &flag.show_gheader,     1 }  
        #define _sT {"show-threads",  no_argument,       &flag.show_threads,     1 }  
        #define _sv {"show-vblocks",  no_argument,       &flag.show_vblocks,     1 }  
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
        #define _s5 {"show-md5",      no_argument,       &flag.show_md5,         1 }
        #define _dS {"test-seg",      no_argument,       &flag.test_seg,         1 }  
        #define _dm {"debug-memory",  no_argument,       &flag.debug_memory,     1 }  
        #define _dp {"debug-progress",no_argument,       &flag.debug_progress,   1 }  
        #define _dh {"show-hash",     no_argument,       &flag.show_hash,        1 }  
        #define _00 {0, 0, 0, 0                                                    }

        typedef const struct option Option;
        static Option genozip_lo[]    = { _i, _I, _c, _d, _f, _h, _l, _L1, _L2, _q, _Q, _t, _DL, _V, _z, _z0, _zb, _zB, _zs, _zS, _zq, _zQ, _zf, _zF, _zc, _zC, _zv, _zV, _m, _th, _u, _o, _p, _e, _E,                                     _ss, _SS, _sd, _sT, _sb, _lc, _lC, _s2, _s7, _S7, _S9, _sa, _st, _sm, _sh, _si, _Si, _Sh, _sr, _sv, _B, _dm, _dp, _dh,_dS, _9, _99, _9s, _9P, _9G, _9g, _9V, _9Q, _9f, _9Z, _9D, _pe, _fa, _bs,         _rg, _sR, _sC, _hC, _rA, _rS, _me,      _s5, _sA, _sc, _sI, _gt,          _00 };
        static Option genounzip_lo[]  = {         _c,     _f, _h,     _L1, _L2, _q, _Q, _t, _DL, _V, _z, _z0, _zb, _zB, _zs, _zS, _zq, _zQ, _zf, _zF, _zc, _zC, _zv, _zV, _m, _th, _u, _o, _p, _e,                                         _ss, _SS, _sd, _sT, _sb, _lc, _lC, _s2, _s7, _S7, _S9, _sa, _st, _sm, _sh, _si, _Si, _Sh, _sr, _sv,     _dm, _dp,                                                                                            _sR, _sC, _hC, _rA, _rS,           _s5, _sA,      _sI,      _cn,     _00 };
        static Option genocat_lo[]    = {                 _f, _h,     _L1, _L2, _q, _Q,          _V,     _z0, _zb, _zB, _zs, _zS, _zq, _zQ, _zf, _zF, _zc, _zC, _zv, _zV,     _th,     _o, _p,         _r, _s, _G, _1, _H0, _H1, _Gt, _GT, _ss, _SS, _sd, _sT, _sb, _lc, _lC, _s2, _s7, _S7, _S9, _sa, _st, _sm, _sh, _si, _Si, _Sh, _sr, _sv,     _dm, _dp,                                                                                   _fs, _g, _sR, _sC, _hC, _rA, _rS,           _s5, _sA,      _sI,      _cn,     _00 };
        static Option genols_lo[]     = {                 _f, _h,     _L1, _L2, _q,              _V,                                                                          _u,     _p, _e,                                                                                                                              _st, _sm,                                   _dm,                                                                                                                                              _b, _00 };
        static Option *long_options[] = { genozip_lo, genounzip_lo, genols_lo, genocat_lo }; // same order as ExeType

        // include the option letter here for the short version (eg "-t") to work. ':' indicates an argument.
        static const char *short_options[NUM_EXE_TYPES] = { // same order as ExeType
            "i:I:cdfhlLqQt^Vzm@:o:p:B:9wWFe:E:2zu", // genozip (note: includes some genounzip options to be used in combination with -d)
            "czfhLqQt^V@:uo:p:me:wW",               // genounzip
            "hLVp:qfub",                            // genols
            "hLV@:p:qQ1r:s:H1Go:fg:e:E:wW"          // genocat
        };

        int option_index = -1;
        int c = getopt_long (argc, argv, short_options[exe_type], long_options[exe_type], &option_index);

        if (c == -1) break; // no more options

        is_short[c] = (option_index == -1); // true if user provided a short option - eg -c rather than --stdout

        switch (c) {
            case PIZ : case LIST : case LICENSE : 
            case VERSION  : case HELP :
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
            case '\7': flag.dump_section  = optarg  ; break;
            case '\2': if (optarg) flag.dict_id_show_one_b250  = dict_id_make (optarg, strlen (optarg)); 
                       else        flag.show_b250 = 1;
                       break;
            case '\3': if (optarg) flag.dict_id_show_one_dict  = dict_id_make (optarg, strlen (optarg)); 
                       else        flag.show_dict = 1;
                       break;
            case '\5': flag.dump_one_b250_dict_id  = dict_id_make (optarg, strlen (optarg)); break;
            case '\6': flag.dump_one_local_dict_id = dict_id_make (optarg, strlen (optarg)); break;
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
                fprintf(stderr, "Usage: %s [OPTIONS] filename1 filename2...\nTry %s --help for more information.\n", global_cmd, global_cmd);
                exit(1);  
        }
    }
}

// if two filenames differ by one character only, which is '1' and '2', creates a combined filename with "1+2"
static char *flags_get_fastq_pair_filename (const char *fn1, const char *fn2)
{
    FileType ft1, ft2;
    char *rn1, *rn2;

    file_get_raw_name_and_type ((char *)fn1, &rn1, &ft1);
    file_get_raw_name_and_type ((char *)fn2, &rn2, &ft2);
    
    unsigned len = strlen (rn1);
    if (len != strlen (rn2)) return NULL;

    int df = -1;
    for (unsigned i=0; i < len; i++)
        if (rn1[i] != rn2[i]) {
            if (df >= 0) return NULL; // 2nd differing character
            df = i;
        }

    if (!((rn1[df] == '1' && rn2[df] == '2') || (rn1[df] == '2' && rn2[df] == '1'))) return NULL; // one of them must be '1' and the other '2'

    char *pair_fn = MALLOC (len+20);
    sprintf (pair_fn, "%.*s1+2%s" FASTQ_GENOZIP_, df, rn1, &rn1[df+1]);
    
    return pair_fn;
}

static void flags_warn_if_duplicates (int num_files, char **filenames)
{
    if (num_files <= 1) return; // nothing to do

    # define BASENAME_LEN 256
    char basenames[num_files * BASENAME_LEN];

    for (unsigned i=0; i < num_files; i++)
        file_basename (filenames[i], false, "", &basenames[i*BASENAME_LEN], BASENAME_LEN);

    qsort (basenames, num_files, BASENAME_LEN, (int (*)(const void *, const void *))strcmp);

    for (unsigned i=1; i < num_files; i++) 
        ASSERTW (strncmp(&basenames[(i-1) * BASENAME_LEN], &basenames[i * BASENAME_LEN], BASENAME_LEN), 
                 "Warning: two files with the same name '%s' - if you later unbind with 'genounzip --unbind %s', these files will overwrite each other", 
                 &basenames[i * BASENAME_LEN], flag.out_filename);
}

void flags_update (unsigned num_files, char **filenames, const bool *is_short)
{
    ASSINP (!flag.to_stdout   || !flag.out_filename,                "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("output", "o"));
    ASSINP (!flag.to_stdout   || !flag.unbind,                      "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("unbind", "u"));
    ASSINP (!flag.unbind      || !flag.out_filename,                "%s: option %s is incompatable with %s", global_cmd, OT("unbind",  "u"), OT("output", "o"));
    ASSINP (!flag.to_stdout   || !flag.replace,                     "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("replace", "^"));

    // if flag.md5 we need to seek back and update the md5 in the txt header section - this is not possible with flag.to_stdout
    ASSINP (!flag.to_stdout   || !flag.md5,                         "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("md5", "m"));
    ASSINP (!flag.test        || !flag.out_filename || command != PIZ, "%s: option %s is incompatable with %s", global_cmd, OT("output", "o"),  OT("test", "t"));
    ASSINP (!flag.test        || !flag.replace || command != PIZ,   "%s: option %s is incompatable with %s", global_cmd, OT("replace", "^"), OT("test", "t"));
    ASSINP (!flag.test        || !flag.to_stdout  || command != ZIP,   "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("test", "t"));
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
    ASSINP (flag.out_dt != DT_VCF || flag.reference,                "%s: --reference must be specified when using --vcf", global_cmd);

    // some genozip flags are allowed only in combination with --decompress 
    if (exe_type == EXE_GENOZIP && command == ZIP) {
        char s[20]; 
        ASSINP (!flag.bgzf,        "%s: option %s can only be used if --decompress is used too", global_cmd, OT("bgzip", "z"));
        ASSINP (flag.out_dt == DT_NONE, "%s: option --%s can only be used if --decompress is used too", global_cmd, str_tolower (dt_name (flag.out_dt), s));
        ASSINP (!flag.unbind,       "%s: option %s can only be used if --decompress is used too", global_cmd, OT("unbind", "u"));
        ASSINP (!flag.show_aliases, "%s: option --show_aliases can only be used if --decompress is used too", global_cmd);
        ASSINP (!flag.show_is_set,  "%s: option --show_is_set can only be used if --decompress is used too", global_cmd);
    }

    // --pair (ZIP only): verify an even number of fastq files, --output, and --reference/--REFERENCE
    if (flag.pair) {

        for (unsigned i=0; i < num_files; i++)
            ASSERT (txtfile_get_file_dt (filenames[i]) == DT_FASTQ, "%s: when using %s, all input files are expected to be FASTQ files, but %s is not", global_cmd, OT("pair", "2"), filenames[i]);

        // in case of a flag.pair with 2 files, in which --output is missing, we attempt to figure it out if possible
        if (!flag.out_filename && num_files == 2) 
            flag.out_filename = flags_get_fastq_pair_filename (filenames[0], filenames[1]);

        ASSINP (flag.out_filename,  "%s: --output must be specified when using %s", global_cmd, OT("pair", "2"));
        ASSINP (flag.reference,     "%s: either --reference or --REFERENCE must be specified when using %s", global_cmd, OT("pair", "2"));
        ASSINP (num_files % 2 == 0, "%s: when using %s, expecting an even number of FASTQ input files, each consecutive two being a pair", global_cmd, OT("pair", "2"));
    }

    // deal with final setting of flag.quiet that suppresses warnings 
    
    // don't show progress or warning when outputing to stdout (note: we are "quiet" even if output doesn't go to the terminal
    // because often it will be piped and ultimately go the terminal)
    if (flag.to_stdout) flag.quiet=true; 
    
    // don't show progress for flags that output throughout the process. no issue with flags that output only in the end
    if (flag.show_dict || flag.show_b250 || flag.show_headers || flag.show_threads || flag.show_bgzf ||
        flag.dict_id_show_one_b250.num || flag.dict_id_show_one_dict.num || flag.show_reference ||
        flag.show_alleles || flag.show_vblocks || flag.show_codec || (flag.show_index && command==PIZ))
        flag.quiet=true; // don't show progress

    // override these ^ if user chose to be --noisy
    if (option_noisy) flag.quiet=false;

    if (flag.test) flag.md5=true; // test implies md5

    // default values, if not overridden by the user
    if (!flag.vblock) vb_set_global_max_memory_per_vb (flag.fast ? TXT_DATA_PER_VB_FAST : TXT_DATA_PER_VB_DEFAULT); 

    // --make-reference implies --md5 --B1 (unless --vblock says otherwise), and not encrypted. 
    // in addition, txtfile_read_vblock() limits each VB to have exactly one contig.
    if (flag.make_reference) {
        ASSINP (!crypt_have_password(), "%s: option --make-reference is incompatable with %s", global_cmd, OT("password", "p"));

        flag.md5 = true;
        if (!flag.vblock) vb_set_global_max_memory_per_vb("1");

        ASSINP (num_files <= 1, "%s: you can specify only one FASTA file when using --make-reference.\n"
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
    if (command == ZIP && flag.out_filename && !flag.quiet) flags_warn_if_duplicates (num_files, filenames);

    flag.multiple_files = (num_files > 1);

    flag.bind = (command == ZIP) && (flag.out_filename != NULL) && (num_files > 1);

    // cases where genocat is used to view some information, but not the file contents
    flag.genocat_info_only = exe_type == EXE_GENOCAT &&
        (flag.show_stats || flag.show_dict || flag.show_b250 || flag.list_chroms || flag.dict_id_show_one_dict.num ||
        flag.show_index || flag.dump_one_local_dict_id.num || flag.dump_one_b250_dict_id.num || flag.dump_section || flag.show_headers ||
        flag.show_reference || flag.show_ref_contigs || flag.show_ref_index || flag.show_ref_hash || flag.show_ref_alts || 
        flag.show_ref_seq || flag.show_aliases || flag.show_txt_contigs);
}

void flags_update_zip_one_file (void)
{
}

// PIZ: called after opening z_file and reading the header before opening txt_file
void flags_update_piz_one_file (void)
{
    // handle native binary formats (BAM). note on BCF and CRAM: we used bcftools/samtools as an external 
    // compressor, so that genozip sees the text, not binary, data of these files - the same as if the file were compressed with eg bz2
    if (command == PIZ && flag.out_dt == DT_NONE && (z_file->flags & SEC_GENOZIP_HEADER_FL_TXT_IS_BIN)) {
        if (z_file->data_type == DT_SAM) 
            // genounzip of a SAM genozip file with is_binary outputs BAM unless the user overrides with --sam or --fastq
            flag.out_dt = DT_BAM;
        // future binary data types here
    }

    // in case translating from SAM.genozip to BAM
    // is this correct ??????
    if (flag.out_dt == DT_BAM) z_file->flags |= SEC_GENOZIP_HEADER_FL_TXT_IS_BIN; // reconstructed file is in binary form
    
    if (flag.out_dt == DT_NONE) 
        flag.out_dt = z_file->data_type;

    // .bcf will be bgzipped by bcftools, ignore --bgzip flag as we don't need an additional bgzip step
    if (flag.out_dt == DT_BCF) flag.bgzf=0;

    // BAM or SEC_GENOZIP_HEADER_FL_BGZF imply bgzf, unless user specifically asked for plain or we're outputting to stdout
    if ((flag.out_dt == DT_BAM || (z_file->flags & SEC_GENOZIP_HEADER_FL_BGZF)) && (!flag.plain && !flag.to_stdout)) 
        flag.bgzf=true;   

    // Note: BAM is stored as binary SAM, so do_translate=true for BAM->BAM , but false for BAM->SAM
    flag.do_translate = dt_get_translation().is_alt_toplevel; 

    // Check if the reconstructed data type is the same as the source data type
    bool is_binary = (z_file->flags & SEC_GENOZIP_HEADER_FL_TXT_IS_BIN);
    flag.reconstruct_as_src = (flag.out_dt == DT_SAM            && z_file->data_type==DT_SAM && !is_binary) || 
                              (flag.out_dt == DT_BAM            && z_file->data_type==DT_SAM && is_binary ) ||
                              (flag.out_dt == z_file->data_type && z_file->data_type!=DT_SAM);
}