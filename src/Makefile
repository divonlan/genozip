# ------------------------------------------------------------------
#   Makefile
#   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
#   Please see terms and conditions in the file LICENSE.txt
#
#   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
#   and subject to penalties specified in the license.

# Note for Windows: to run this make, you need mingw (for the gcc compiler) and cygwin (for Unix-like tools):
# Mingw: http://mingw-w64.org/doku.php 
# Cygwin: https://www.cygwin.com/

OBJDIR         := objdir         # fallback if not win, linux, mac
OBJDIR_LINUX   := objdir.linux
OBJDIR_WINDOWS := objdir.windows
OBJDIR_MAC     := objdir.mac


ifeq ($(OS),Windows_NT)
	Windows := 1
	OBJDIR=$(OBJDIR_WINDOWS)
	WSL=wsl
else
	uname := $(shell uname -s)

	ifeq ($(uname),Linux)
		Linux := 1
		OBJDIR=$(OBJDIR_LINUX)
	endif

    ifeq ($(uname),Darwin) # Mac
		Mac := 1
		OBJDIR=$(OBJDIR_MAC)
    endif
endif

# note: modern nasm/yasm have AS_FEATURE_LEVEL=10. HAVE_AS_KNOWS_AVX512 is true if level 10.
ASFLAGS = -DAS_FEATURE_LEVEL=10 -DHAVE_AS_KNOWS_AVX512

CPU = $(shell uname -m | sed -e 's/amd/x86_/')

ifeq ($(CPU),x86_64) # or amd64
	AS = nasm

	ifdef Windows
		x86_64 := 1
		ASFLAGS += -f win64
	endif

	ifdef Linux
		x86_64 := 1
		ASFLAGS += -f elf64 -DINTEL_CET_ENABLED
	endif

	ifdef Mac
		# Mac is not tested yet with .asm - fix this after testing
		# noarch := 1
		x86_64 := 1   
		ASFLAGS += -f macho64 --prefix=_
	endif	
else 
	ifeq ($(CPU),aarch64) # ARM 64
		aarch64 := 1
		AS = $(CC) -D__ASSEMBLY__ -c
		DEBUGFLAGS_AS = -g 
	else 
		noarch := 1
	endif
endif

ifndef_any_of = $(filter undefined,$(foreach v,$(1),$(origin $(v))))
ifdef_any_of  = $(filter-out undefined,$(foreach v,$(1),$(origin $(v))))

ifdef BUILD_PREFIX
IS_CONDA=1
endif

LDFLAGS    += -lpthread -lm 
CFLAGS     += -Wall -D_LARGEFILE64_SOURCE=1 -D_7ZIP_ST -D_GNU_SOURCE -std=gnu11 -fms-extensions -msse4.1 -msse4.2
CXXFLAGS   += -Wall -D_LARGEFILE64_SOURCE=1 -D_GNU_SOURCE -Ibsc -O3 -DLIBBSC_CPU_FEATURE=LIBBSC_CPU_FEATURE_SSE41 -msse4.1
OPTFLAGS   += -O3 
DEBUGFLAGS += -DDEBUG -ggdb3 -O0  

ifdef x86_64
	CFLAGS += -DGENOZIP_ALLOW_UNALIGNED_ACCESS 
	CXXFLAGS += -DLIBBSC_ALLOW_UNALIGNED_ACCESS
endif

ifdef IS_CONDA 
	DISTOBJ = $(OBJDIR)/distribution.conda.o # note: matches the name in version.c 
	CFLAGS += -O3 

	ifdef Windows
		CC=gcc # in Windows, override conda's default Visual C with gcc 
		CC=g++
		LDFLAGS += -L$(PREFIX)/Library/lib
		CFLAGS  += -I$(PREFIX)/Library/include 
	endif

	ifeq ($(target_architecture),"osx-arm64")
		$(error Genozip cannot yet be built for osx-arm64)
	endif

else
	ifdef DISTRIBUTION # custom destribution, run eg: DISTRIBUTION=special make genozip
		DISTOBJ = $(OBJDIR)/distribution.$(DISTRIBUTION).o 
	else
		DISTOBJ = $(OBJDIR)/distribution.github.o # note: matches the name in version.c 
	endif

	ifdef GCC 
		CC=$(GCC) 
	else
		CC=gcc
	endif

	ifdef GXX 
		CXX=$(GXX) 
	else
		CXX=g++
	endif

	# set -march. note: if conda, we leave it to conda to set
	ifdef Windows
		CFLAGS += -march=haswell 
	endif

	ifdef Linux
		CFLAGS += -march=haswell # see: https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html  
	endif

	ifdef Mac
		CFLAGS += -march=native # conda used "core2". We use the highest available optimizations, because it is typically compiled for use on a specific Mac
				# CFLAGS += -mmacosx-version-min=10.13 # for utimensat (in conda we use utimes instead)
	endif

endif 

WSL_GCC_8_5 = /home_local_fs/divon/miniconda3/envs/distribution/bin/x86_64-conda-linux-gnu-gcc # path to gcc 8.5 on WSL (generate with: conda create -n distribution gcc_linux-64=8.5.0)
WSL_GCC_DEFAULT = /home_local_fs/divon/miniconda3/bin/gcc

DEV_DIR  = ../../genozip
TEST_DIR = $(DEV_DIR)/private/test
SCRIPTS  = $(DEV_DIR)/private/scripts
UTILS    = $(DEV_DIR)/private/utils
RELEASES = $(DEV_DIR)/private/releases
DATA     = $(DEV_DIR)/data
LATEST_SRC = ../../genozip-latest/src

SRC_DIRS = secure zlib bzlib lzma bsc libdeflate_1.7 libdeflate_1.19 libdeflate_1.19/x86 libdeflate_1.19/arm 	\
		  htscodecs igzip igzip/aarch64 igzip/x86_64 igzip/noarch

MY_SRCS = genozip.c genols.c context.c container.c strings.c stats.c arch.c tip.c seg_id.c zip_dyn_int.c     	\
		  data_types.c bits.c progress.c writer.c zriter.c tar.c chrom.c qname.c tokenizer.c mutex.c threads.c	\
          zip.c piz.c reconstruct.c recon_history.c recon_peek.c seg.c zfile.c aligner.c flags.c digest.c    	\
		  reference.c contigs.c ref_lock.c refhash.c ref_make.c ref_contigs.c ref_iupacs.c ref_cache.c     	 	\
		  vcf_piz.c vcf_seg.c vcf_vblock.c vcf_header.c vcf_info.c vcf_samples.c vcf_hgvs.c      				\
		  vcf_format_GT.c vcf_format_PS_PID.c vcf_dbsnp.c vcf_giab.c vcf_vep.c      							\
		  vcf_refalt.c vcf_format.c vcf_illum_gtyping.c vcf_gwas.c vcf_vagrent.c vcf_svaba.c vcf_pbsv.c			\
		  vcf_icgc.c vcf_snpeff.c vcf_cosmic.c vcf_mastermind.c vcf_isaac.c	vcf_manta.c vcf_pos.c vcf_ultima.c	\
		  vcf_platypus.c vcf_info_AC_AF_AN.c vcf_format_GQ.c vcf_gatk.c vcf_sv.c								\
		  sam_seg.c sam_piz.c sam_shared.c sam_header.c sam_md.c sam_nm.c sam_tlen.c sam_cigar.c sam_fields.c  	\
		  sam_sa.c bam_seg.c bam_seq.c bam_show.c sam_pacbio.c sam_ultima.c sam_xcons.c cram.c agilent.c   		\
		  sam_seq.c sam_qual.c sam_sag_zip.c sam_sag_piz.c sam_sag_load.c sam_sag_ingest.c sam_sag_scan.c    	\
		  sam_bwa.c sam_bowtie2.c sam_bsseeker2.c sam_bsbolt.c sam_bismark.c sam_gem3.c sam_tmap.c sam_hisat2.c \
		  sam_blasr.c sam_dragen.c sam_minimap2.c sam_cellranger.c sam_biobambam.c sam_pos.c sam_deep.c 		\
		  sam_star.c sam_abra2.c																				\
		  fastq.c fastq_desc.c fastq_seq.c fastq_qual.c fastq_deep.c fastq_saux.c deep.c						\
		  fasta.c gff.c bed.c me23.c locs.c generic.c lookback.c compressor.c									\
		  buffer.c buf_struct.c buf_list.c random_access.c sections.c base64.c bgzf.c coverage.c txtheader.c  	\
		  codec.c codec_bz2.c codec_lzma.c codec_acgt.c codec_domq.c codec_bsc.c codec_pacb.c					\
		  codec_pbwt.c codec_none.c codec_htscodecs.c codec_longr.c codec_normq.c codec_homp.c codec_t0.c		\
		  codec_smux.c																							\
	      txtfile.c profiler.c file.c filename.c dispatcher.c crypt.c aes.c md5.c segconf.c biopsy.c gencomp.c	\
		  vblock.c regions.c optimize.c dict_id.c aliases.c hash.c stream.c url.c bases_filter.c dict_io.c		\
		  recon_plan_io.c version.c huffman.c user_message.c b250.c qname_filter.c
		  
ZLIB_SRCS  = zlib/gzlib.c zlib/zutil.c zlib/deflate.c zlib/trees.c

BZLIB_SRCS = bzlib/blocksort.c bzlib/bzlib.c bzlib/compress.c bzlib/crctable.c bzlib/decompress.c bzlib/huffman.c bzlib/randtable.c

LZMA_SRCS  = lzma/LzmaEnc.c lzma/LzmaDec.c lzma/LzFind.c

BSC_C_SRCS = bsc/libsais.c

BSC_CXX_SRCS   = bsc/bwt.cpp bsc/coder.cpp bsc/detectors.cpp bsc/libbsc.cpp bsc/lzp.cpp 						\
		  bsc/platform.cpp bsc/preprocessing.cpp bsc/qlfc.cpp bsc/qlfc_model.cpp bsc/adler32.cpp

HTSCODECS_SRC = htscodecs/rANS_static4x16pr.c htscodecs/rle.c htscodecs/pack.c htscodecs/arith_dynamic.c

LIBDEFLATE_1_7_SRCS = libdeflate_1.7/deflate_compress.c libdeflate_1.7/utils.c  

LIBDEFLATE_1_19_SRCS = libdeflate_1.19/deflate_compress.c libdeflate_1.19/deflate_decompress.c 					\
		  libdeflate_1.19/utils.c libdeflate_1.19/x86/cpu_features.c libdeflate_1.19/arm/cpu_features.c 		\
		  libdeflate_1.19/crc32.c libdeflate_1.19/adler32.c

CONDA_DEVS = Makefile ../.gitignore 

CONDA_DOCS = ../LICENSE.txt ../AUTHORS ../README.md

CONDA_INCS += dict_id_gen.h aes.h dispatcher.h optimize.h profiler.h dict_id.h aliases.h txtfile.h zip.h bits.h progress.h website.h 		\
            endianness.h md5.h sections.h text_help.h strings.h hash.h stream.h url.h flags.h segconf.h biopsy.h huffman.h 					\
            buffer.h buf_struct.h buf_list.h file.h context.h context_struct.h container.h seg.h text_license.h version.h compressor.h 		\
            crypt.h genozip.h piz.h vblock.h zfile.h random_access.h regions.h reconstruct.h tar.h qname.h qname_flavors.h codec.h  		\
		 	lookback.h tokenizer.h codec_longr_alg.c gencomp.h dict_io.h recon_plan_io.h tip.h deep.h filename.h stats.h multiplexer.h 		\
		 	reference.h ref_private.h refhash.h ref_iupacs.h aligner.h mutex.h bgzf.h coverage.h threads.h local_type.h 					\
			arch.h license.h file_types.h data_types.h base64.h txtheader.h writer.h zriter.h bases_filter.h genols.h contigs.h chrom.h 	\
			vcf.h vcf_private.h sam.h sam_private.h me23.h fasta.h fasta_private.h gff.h bed.h locs.h generic.h 							\
			fastq.h fastq_private.h user_message.h mac_compat.h b250.h zip_dyn_int.h qname_filter.h 										\
			\
			zlib/gzguts.h zlib/zconf.h zlib/deflate.h zlib/trees.h zlib/zlib.h zlib/zutil.h													\
			\
			lzma/7zTypes.h lzma/LzFind.h lzma/LzHash.h lzma/LzmaDec.h lzma/LzmaEnc.h 														\
			\
			bzlib/bzlib.h bzlib/bzlib_private.h																								\
			\
			htscodecs/rANS_static4x16.h htscodecs/rle.h htscodecs/pack.h htscodecs/arith_dynamic.h htscodecs/c_simple_model.h				\
			htscodecs/rANS_word.h htscodecs/htscodecs_endian.h htscodecs/rANS_word.h htscodecs/utils.h htscodecs/varint.h 					\
			htscodecs/varint2.h htscodecs/utils.h 																							\
			\
			bsc/adler32.h bsc/bwt.h bsc/coder.h bsc/coder_tables.h bsc/filter_tables.h bsc/filters.h bsc/libbsc.h bsc/libsais.h bsc/lzp.h 	\
			bsc/platform.h bsc/predictor.h bsc/qlfc.h bsc/qlfc_model.h bsc/rangecoder.h														\
 			\
			libdeflate_1.7/unaligned.h libdeflate_1.7/deflate_compress.h libdeflate_1.7/arm_matchfinder_impl.h 								\
			libdeflate_1.7/deflate_constants.h libdeflate_1.7/bt_matchfinder.h libdeflate_1.7/hc_matchfinder.h libdeflate_1.7/common_defs.h	\
			libdeflate_1.7/lib_common.h libdeflate_1.7/x86_matchfinder_impl.h libdeflate_1.7/compiler_gcc.h libdeflate_1.7/libdeflate.h 	\
			libdeflate_1.7/matchfinder_common.h 																							\
 			\
			libdeflate_1.19/adler32_vec_template.h libdeflate_1.19/crc32_tables.h 									 						\
			libdeflate_1.19/crc32_multipliers.h libdeflate_1.19/deflate_constants.h libdeflate_1.19/lib_common.h 							\
			libdeflate_1.19/bt_matchfinder.h libdeflate_1.19/hc_matchfinder.h libdeflate_1.19/common_defs.h	libdeflate_1.19/libdeflate.h	\
			libdeflate_1.19/gzip_constants.h libdeflate_1.19/decompress_template.h libdeflate_1.19/cpu_features_common.h					\
			libdeflate_1.19/matchfinder_common.h libdeflate_1.19/deflate_compress.h libdeflate_1.19/ht_matchfinder.h 						\
			libdeflate_1.19/x86/adler32_impl.h libdeflate_1.19/x86/crc32_impl.h libdeflate_1.19/x86/matchfinder_impl.h						\
			libdeflate_1.19/x86/cpu_features.h libdeflate_1.19/x86/crc32_pclmul_template.h libdeflate_1.19/x86/decompress_impl.h			\
			libdeflate_1.19/arm/adler32_impl.h libdeflate_1.19/arm/cpu_features.h libdeflate_1.19/arm/crc32_impl.h                          \
			libdeflate_1.19/arm/matchfinder_impl.h libdeflate_1.19/arm/crc32_pmull_helpers.h libdeflate_1.19/arm/crc32_pmull_wide.h			\
			\
			igzip/bitbuf2.h igzip/crc.h igzip/encode_df.h igzip/flatten_ll.h igzip/huff_codes.h igzip/huffman.h igzip/igzip_checksums.h 	\
		  	igzip/igzip_level_buf_structs.h igzip/igzip_lib.h igzip/igzip_wrapper.h igzip/repeated_char_result.h igzip/static_inflate.h 	\
		  	igzip/unaligned.h																												\
			\
			igzip/x86_64/bitbuf2.asm igzip/x86_64/data_struct2.asm igzip/x86_64/heap_macros.asm igzip/x86_64/huffman.asm 					\
			igzip/x86_64/igzip_compare_types.asm igzip/x86_64/igzip_decode_block_stateless.asm igzip/x86_64/igzip_update_histogram.asm 		\
			igzip/x86_64/inflate_data_structs.asm igzip/x86_64/lz0a_const.asm igzip/x86_64/multibinary.asm igzip/x86_64/options.asm 		\
			igzip/x86_64/reg_sizes.asm igzip/x86_64/stdmac.asm																				\
			\
			igzip/aarch64/aarch64_label.h igzip/aarch64/aarch64_multibinary.h igzip/aarch64/bitbuf2_aarch64.h 								\
			igzip/aarch64/crc32_aarch64_common.h igzip/aarch64/crc32_gzip_refl_pmull.h igzip/aarch64/crc32_norm_common_pmull.h 				\
			igzip/aarch64/crc32_refl_common_pmull.h igzip/aarch64/crc_common_pmull.h igzip/aarch64/data_struct_aarch64.h 					\
			igzip/aarch64/huffman_aarch64.h igzip/aarch64/lz0a_const_aarch64.h igzip/aarch64/options_aarch64.h 								\
			igzip/aarch64/stdmac_aarch64.h igzip/aarch64/crc32_common_crc_ext_cortex_a72.S igzip/aarch64/crc32_common_mix_neoverse_n1.S 	\
			igzip/aarch64/crc32_gzip_refl_3crc_fold.S igzip/aarch64/crc32_gzip_refl_crc_ext.S igzip/aarch64/crc32_gzip_refl_pmull.S			\
			igzip/aarch64/crc32_mix_default.S igzip/aarch64/crc32_mix_default_common.S igzip/aarch64/crc32_mix_neoverse_n1.S 				\
			igzip/aarch64/crc32c_mix_default.S igzip/aarch64/crc32c_mix_neoverse_n1.S igzip/aarch64/crc_multibinary_arm.S 

INSTALLERS := ../installers

IGZIP_SRCS = igzip/igzip.c igzip/hufftables_c.c igzip/igzip_base.c igzip/igzip_icf_base.c igzip/adler32_base.c igzip/flatten_ll.c 			\
		  igzip/encode_df.c igzip/igzip_icf_body.c igzip/crc_base.c igzip/huff_codes.c igzip/igzip_inflate.c 

IGZIP_X86_ASM_SRCS := igzip/x86_64/igzip_body.asm igzip/x86_64/igzip_finish.asm igzip/x86_64/igzip_icf_body_h1_gr_bt.asm 					\
		  igzip/x86_64/igzip_icf_finish.asm igzip/x86_64/rfc1951_lookup.asm igzip/x86_64/adler32_sse.asm 									\
		  igzip/x86_64/adler32_avx2_4.asm igzip/x86_64/igzip_multibinary.asm igzip/x86_64/igzip_update_histogram_01.asm 					\
		  igzip/x86_64/igzip_update_histogram_04.asm igzip/x86_64/igzip_decode_block_stateless_01.asm 										\
		  igzip/x86_64/igzip_decode_block_stateless_04.asm igzip/x86_64/igzip_inflate_multibinary.asm 										\
		  igzip/x86_64/encode_df_04.asm igzip/x86_64/encode_df_06.asm igzip/x86_64/proc_heap.asm igzip/x86_64/igzip_deflate_hash.asm 		\
		  igzip/x86_64/igzip_gen_icf_map_lh1_06.asm igzip/x86_64/igzip_gen_icf_map_lh1_04.asm igzip/x86_64/igzip_set_long_icf_fg_04.asm 	\
		  igzip/x86_64/igzip_set_long_icf_fg_06.asm igzip/x86_64/crc_multibinary.asm igzip/x86_64/crc32_gzip_refl_by8.asm 					\
		  igzip/x86_64/crc32_gzip_refl_by8_02.asm igzip/x86_64/crc32_gzip_refl_by16_10.asm

IGZIP_ARM_ASM_SRCS := igzip/aarch64/igzip_inflate_multibinary_arm64.S igzip/aarch64/igzip_multibinary_arm64.S 								\
		  igzip/aarch64/igzip_isal_adler32_neon.S igzip/aarch64/igzip_deflate_body_aarch64.S igzip/aarch64/igzip_deflate_finish_aarch64.S	\
		  igzip/aarch64/isal_deflate_icf_body_hash_hist.S igzip/aarch64/isal_deflate_icf_finish_hash_hist.S									\
		  igzip/aarch64/igzip_set_long_icf_fg.S igzip/aarch64/encode_df.S igzip/aarch64/isal_update_histogram.S igzip/aarch64/gen_icf_map.S \
		  igzip/aarch64/igzip_deflate_hash_aarch64.S igzip/aarch64/igzip_decode_huffman_code_block_aarch64.S 

ifdef x86_64
IGZIP_ARCH_DIR := igzip/x86_64
endif
ifdef aarch64
IGZIP_ARCH_DIR := igzip/aarch64
IGZIP_SRCS     += igzip/aarch64/crc_aarch64_dispatcher.c igzip/aarch64/igzip_multibinary_aarch64_dispatcher.c igzip/proc_heap_base.c
endif
ifdef noarch
IGZIP_ARCH_DIR := igzip/noarch
IGZIP_SRCS     += igzip/noarch/igzip_base_aliases.c igzip/noarch/proc_heap_base.c igzip/noarch/crc_base_aliases.c
endif

ifeq ($(CC),cl) # Microsoft Visual C
	$(error Only the gcc and clang compilers are currently supported)
endif

ifneq ($(wildcard secure/.*),) # "secure" directory exists
	MY_SRCS += secure/license.c
	CONDA_LIC_OBJS := $(OBJDIR_LINUX)/secure/license.o $(OBJDIR_MAC)/secure/license.o
else
	LDFLAGS += $(OBJDIR)/secure/license.o
endif

ifdef Windows
	EXE = .exe
	LDFLAGS += -static-libgcc -lIPHLPAPI -lUrlmon -lDbgHelp
endif

ifdef Linux
	LDFLAGS += -lrt -fuse-ld=gold -rdynamic # rt is required by pthreads ; rdynamic for backtrace_symbols_fd to show function names
endif

ifdef Mac # Mac
	MY_SRCS += mac_compat.c
	CFLAGS += -mmacosx-version-min=10.9 -DNO_UTIMENSAT # Conda requires support from 10.9 but utimensat() was only introduced in MacOS 10.13
	ifeq ($(CC),gcc) 
		LDFLAGS += -rdynamic
	endif
endif

C_SRCS = $(MY_SRCS) $(ZLIB_SRCS) $(BZLIB_SRCS) $(LZMA_SRCS) $(LIBDEFLATE_1_7_SRCS) $(LIBDEFLATE_1_19_SRCS) $(BSC_C_SRCS) \
		 $(HTSCODECS_SRC) $(IGZIP_SRCS)	

CXX_SRCS = $(BSC_CXX_SRCS) 

X86_SRCS = $(IGZIP_X86_ASM_SRCS) # note: these are the sources which pass to the compiler - they may include other .asm files. These other files are only included in CONDA_INCS.

ARM_SRCS = $(IGZIP_ARM_ASM_SRCS)

ASFLAGS += -I$(IGZIP_ARCH_DIR)

# ifdef Linux
# 	C_SRCS += $(ONION_SRCS)
# endif

ifdef x86_64
IGZIP_DEPS   = $(addprefix $(OBJDIR)/,$(X86_SRCS:.asm=.d)) 
IGZIP_OBJS   = $(addprefix $(OBJDIR)/,$(X86_SRCS:.asm=.o))
# BSC_OPT_OBJS = $(addprefix $(OBJDIR)/,bsc/qlfc_AVX2.o bsc/qlfc_AVX.o bsc/qlfc_SSE41.o bsc/adler32_AVX.o bsc/adler32_SSSE3.o)
endif

ifdef aarch64
IGZIP_DEPS   = $(addprefix $(OBJDIR)/,$(ARM_SRCS:.S=.d)) 
IGZIP_OBJS   = $(addprefix $(OBJDIR)/,$(ARM_SRCS:.S=.o))
# BSC_OPT_OBJS = $(addprefix $(OBJDIR)/,bsc/qlfc_A64.o bsc/adler32_A64.o)
endif

OBJS       := $(addprefix $(OBJDIR)/,$(C_SRCS:.c=.o))       $(addprefix $(OBJDIR)/,$(CXX_SRCS:.cpp=.o))       $(IGZIP_OBJS) $(BSC_OPT_OBJS)
DEBUG_OBJS := $(addprefix $(OBJDIR)/,$(C_SRCS:.c=.debug-o)) $(addprefix $(OBJDIR)/,$(CXX_SRCS:.cpp=.debug-o)) $(IGZIP_OBJS) 
OPT_OBJS   := $(addprefix $(OBJDIR)/,$(C_SRCS:.c=.opt-o))   $(addprefix $(OBJDIR)/,$(CXX_SRCS:.cpp=.opt-o))   $(IGZIP_OBJS) $(BSC_OPT_OBJS)  # optimized but with debug info, for debugging issues that only manifest with compiler optimization
DEPS       := $(addprefix $(OBJDIR)/,$(C_SRCS:.c=.d))       $(addprefix $(OBJDIR)/,$(CXX_SRCS:.cpp=.d))       $(IGZIP_DEPS)

EXECUTABLES       = genozip$(EXE)       genounzip$(EXE)       genocat$(EXE)       genols$(EXE)
DEBUG_EXECUTABLES = genozip-debug$(EXE) genounzip-debug$(EXE) genocat-debug$(EXE) genols-debug$(EXE)
OPT_EXECUTABLES   = genozip-opt$(EXE)   genounzip-opt$(EXE)   genocat-opt$(EXE)   genols-opt$(EXE)

all   : CFLAGS += $(OPTFLAGS) 
all   : CXXFLAGS += $(OPTFLAGS) 
all   : $(OBJDIR) $(EXECUTABLES) 
	@chmod +x test.sh

# ifdef Windows
# debug opt : LDFLAGS += -lbacktrace
# endif

debug : CFLAGS += $(DEBUGFLAGS) 
debug : CXXFLAGS += $(DEBUGFLAGS) 
debug : $(OBJDIR) $(DEBUG_EXECUTABLES)

opt   : CFLAGS += -ggdb3 $(OPTFLAGS) -DDEBUG 
opt   : $(OBJDIR) $(OPT_EXECUTABLES)

ifdef Windows
objdir.mac/secure/license.o : secure/license.c secure/license_buf.c secure/license_stats.c secure/license_sys.c secure/license_zfile.c
	@echo "getting $>"
	@wsl mac-get-license-o.sh
endif

ifdef Mac
objdir.mac/secure/license.o : CFLAGS := $(patsubst -march=% ,-march=x86-64 ,$(CFLAGS)) # plain vanilla x86-64, so it runs on Rosetta x86-64 simulator on M1/M2/M3
endif

-include $(DEPS)

$(OBJDIR):
	@echo Making directory $@ and its subdirs
	@mkdir -p $@ $(addprefix $@/,$(SRC_DIRS)) 

$(OBJDIR)/%.d: %.c | $(OBJDIR) # directory is an "order only prerequesite": https://www.gnu.org/savannah-checkouts/gnu/make/manual/html_node/Prerequisite-Types.html#Prerequisite-Types
	@echo Calculating dependencies $<
	@$(CC) $(CFLAGS) -MM -MT $@ $< -MF $(@:%.o=%.d)

$(OBJDIR)/%.d: %.cpp | $(OBJDIR) # directory is an "order only prerequesite": https://www.gnu.org/savannah-checkouts/gnu/make/manual/html_node/Prerequisite-Types.html#Prerequisite-Types
	@echo Calculating dependencies $<
	@$(CXX) $(CXXFLAGS) -MM -MT $@ $< -MF $(@:%.o=%.d)

$(OBJDIR)/distribution.%.o : distribution.c 
	@echo Compiling $< to $@
	@$(CC) -c -o $@ $< $(CFLAGS) -DDISTRIBUTION=\"$*\"

# $(OBJDIR)/bsc/qlfc_%.o : bsc/qlfc.cpp 
# 	@echo Compiling $< to $@
# 	$(CXX) -c -o $@ $< $(CXXFLAGS) -DLIBBSC_DYNAMIC_CPU_DISPATCH -DLIBBSC_CPU_FEATURE=LIBBSC_CPU_FEATURE_$*

# $(OBJDIR)/bsc/adler32_%.o : bsc/adler32.cpp 
# 	@echo Compiling $< to $@
# 	$(CXX) -c -o $@ $< $(CXXFLAGS) -DLIBBSC_DYNAMIC_CPU_DISPATCH -DLIBBSC_CPU_FEATURE=LIBBSC_CPU_FEATURE_$*

$(OBJDIR)/%.o: %.c $(OBJDIR)/%.d
	@echo Compiling $<
	@$(CC) -c -o $@ $< $(CFLAGS) 

$(OBJDIR)/%.debug-o: %.c $(OBJDIR)/%.d
	@echo "Compiling $< (debug)"
	@$(CC) -c -o $@ $< $(CFLAGS)

$(OBJDIR)/%.opt-o: %.c $(OBJDIR)/%.d
	@echo "Compiling $< (opt)"
	@$(CC) -c -o $@ $< $(CFLAGS)

$(OBJDIR)/%.o: %.cpp $(OBJDIR)/%.d
	@echo Compiling $<
	@$(CXX) -c -o $@ $< $(CXXFLAGS) 

$(OBJDIR)/%.debug-o: %.cpp $(OBJDIR)/%.d
	@echo "Compiling $< (debug)"
	@$(CXX) -c -o $@ $< $(CXXFLAGS)

$(OBJDIR)/%.opt-o: %.cpp $(OBJDIR)/%.d
	@echo "Compiling $< (opt)"
	@$(CXX) -c -o $@ $< $(CXXFLAGS)

$(OBJDIR)/%.d: %.asm
	@echo Calculating dependencies $<
	@$(AS) $(ASFLAGS) -M $< > $@

$(OBJDIR)/%.d: %.S
	@echo Calculating dependencies $<
	@$(AS) $(ASFLAGS) -M $< > $@

$(OBJDIR)/%.o: %.asm # nasm or yasm
	@echo Compiling $<
	@$(AS) $(ASFLAGS) -o $@ $< 

$(OBJDIR)/%.o: %.S   # GNU Assembler 
	@echo Compiling $<
	@$(AS) $(ASFLAGS) -o $@ $< 

%.S: %.c $(OBJDIR)/%.d
	@echo "Generating $@"
	@$(CC) -S -o $@ $< $(CFLAGS)

%.E: %.c $(OBJDIR)/%.d
	@echo "Generating $@"
	@$(CC) -E -o $@ $< $(CFLAGS)

GENDICT_OBJS := $(addprefix $(OBJDIR)/,$(GENDICT_SRCS:.c=.o))

# dict_id_gen.h generation:
# Step 1: dict_id_gen.sh generates dict_id_gen.c, including all the GENDICT definitions from the data type include files (eg vcf.h)
# Step 2: dict_id_gen.sh compiles dict_id_gen.c and generate dict_id_gen[.exe]
# Step 3. dict_id_gen.sh generates dict_id_gen.h: it uses dict_id_gen[.exe]to generate the field constant, and then adds the fields enum and mapping
ifdef Windows
dict_id_gen.h : $(shell grep -w "pragma GENDICT" *.h | cut -d: -f1 | uniq) dict_id_gen.sh
	@echo Generating $@
	@./dict_id_gen.sh $(CC)

endif # ugly hack to avoid conda failure due to bash issues in dict_id_gen.sh - pre-generate on Windows and check in to github

genozip$(EXE): dict_id_gen.h $(OBJS) $(DISTOBJ)
	@echo Linking $@
	@$(CC) -o $@ $(OBJS) $(DISTOBJ) $(CFLAGS) $(LDFLAGS) 
	@if ! ln -f $@ ../$@; then ln -fs $@ ../$@; fi  # hard linking not possible on some filesystems
 
genozip-debug$(EXE): dict_id_gen.h $(DEBUG_OBJS) $(DISTOBJ)
	@echo Linking $@
	@$(CC) -o $@ $(DEBUG_OBJS) $(CFLAGS) $(LDFLAGS) $(DISTOBJ)
	@if ! ln -f $@ ../$@; then ln -fs $@ ../$@; fi

genozip-opt$(EXE): dict_id_gen.h $(OPT_OBJS) $(DISTOBJ)
	@echo Linking $@
	@$(CC) -o $@ $(OPT_OBJS) $(CFLAGS) $(LDFLAGS) $(DISTOBJ)
	@if ! ln -f $@ ../$@; then ln -fs $@ ../$@; fi

genounzip$(EXE) genocat$(EXE) genols$(EXE): genozip$(EXE)
	@echo Hard linking $@
	@rm -f $@ 
	@ln -f $^ $@
	@if ! ln -f $^ $@; then ln -fs $^ $@; fi
	@if ! ln -f $@ ../$@; then ln -fs $@ ../$@; fi

genounzip-debug$(EXE) genocat-debug$(EXE) genols-debug$(EXE): genozip-debug$(EXE)
	@echo Hard linking $@
	@rm -f $@ 
	@if ! ln -f $^ $@; then ln -fs $^ $@; fi
	@if ! ln -f $@ ../$@; then ln -fs $@ ../$@; fi

genounzip-opt$(EXE) genocat-opt$(EXE) genols-opt$(EXE): genozip-opt$(EXE)
	@echo Hard linking $@
	@rm -f $@ 
	@if ! ln -f $^ $@; then ln -fs $^ $@; fi
	@if ! ln -f $@ ../$@; then ln -fs $@ ../$@; fi

../LICENSE.txt: text_license.h version.h # not dependent on genozip.exe, so we don't generate it every compilation
	@$(MAKE) -j genozip$(EXE) # recursive call to make genozip.exe with the latest version
	@echo Generating $@
	@./genozip$(EXE) --license=txt,100 > $@
	
# this is used by build.sh to install on conda for Linux and Mac. Installation for Windows in in bld.bat
install: genozip$(EXE)
	@echo Installing in $(PREFIX)/bin
	@if ( test ! -d $(PREFIX)/bin ) ; then mkdir -p $(PREFIX)/bin ; fi
	@cp -f genozip$(EXE) $(PREFIX)/bin/genozip$(EXE)
ifndef Windows
	@chmod a+x $(PREFIX)/bin/genozip$(EXE)
endif
	@cp -f $(PREFIX)/bin/genozip$(EXE) $(PREFIX)/bin/genounzip$(EXE)
	@cp -f $(PREFIX)/bin/genozip$(EXE) $(PREFIX)/bin/genocat$(EXE)
	@cp -f $(PREFIX)/bin/genozip$(EXE) $(PREFIX)/bin/genols$(EXE)

version = $(shell head -n1 version.h |cut -d\" -f2)

SH_VERIFY_ALL_COMMITTED = (( `git status |grep 'modified\|Untracked files'|grep -v .gitkeep |wc -l` == 0 )) || 				\
                          (echo ERROR: there are some uncommitted changes: ; echo ; git status ; exit 1) ;     				\
						  (( `cd secure; git status |grep 'modified\|Untracked files'|grep -v .gitkeep |wc -l` == 0 )) || 	\
                          (echo ERROR: there are some uncommitted changes: ; echo ; cd secure ; git status ; exit 1) 

test:
	@cat test.sh | tr -d "\r" | bash -

clean-installers:
	@rm -fR $(INSTALLERS)/*

clean-debug:
	@echo Cleaning up debug
	@rm -f $(DEBUG_OBJS) $(DEBUG_EXECUTABLES) $(OBJDIR)/*.debug-o
	@rm -f $(OPT_OBJS) $(OPT_EXECUTABLES) $(OBJDIR)/*.opt-o 

clean:  # careful not to delete license.o
	@echo Cleaning up optimized
	@rm -f $(EXECUTABLES) $(OBJS) $(DEPS) $(OBJDIR)/*.o $(OBJDIR)/*.d

clean-distribution: # clean all platforms
	@bash --norc --noprofile secure/prepare.sh # before clean removes executables
	@(cd secure ; git stage secure.data ; exit 0) > /dev/null
	@(cd secure ; git commit -m $(version) ; exit 0) > /dev/null
	@(cd secure ; git push > /dev/null)
	@$(SH_VERIFY_ALL_COMMITTED)
	@echo Cleaning up all platforms for distribution
	@rm -f $(OBJDIR_LINUX)/*.o $(OBJDIR_LINUX)/*.d $(OBJDIR_WINDOWS)/*.o $(OBJDIR_WINDOWS)/*.d $(OBJDIR_MAC)/*.o $(OBJDIR_MAC)/*.d
	
clean-opt:
	@echo Cleaning up opt
	@rm -f $(OPT_OBJS) $(EXECUTABLES) $(OPT_EXECUTABLES)/*.opt-o

clean-test:
	@rm -f $(TEST_DIR)/*.genozip $(TEST_DIR)/*.good $(TEST_DIR)/*.bad $(TEST_DIR)/*.local $(TEST_DIR)/*.b250 $(TEST_DIR)/tmp/* $(TEST_DIR)/*.DEPN $(TEST_DIR)/back-compat/*/*.bad

clean-dev: clean-installers clean-test # use only in dev - deletes license.o
	@echo Cleaning up
	@rm -f $(DEPS) $(filter-out LICENSE.txt,$(WINDOWS_INSTALLER_OBJS)) *.d .archive.tar.gz *.stackdump $(EXECUTABLES) $(OPT_EXECUTABLES) $(DEBUG_EXECUTABLES) 
	@rm -f *.S *.good *.bad $(DATA)/*.good $(DATA)/*.bad *.local genozip.threads-log.* *.b250 vgcore.*
	@rm -R $(OBJDIR)
	@mkdir $(OBJDIR) $(addprefix $(OBJDIR)/,$(SRC_DIRS))
	@touch $(OBJDIR)/.gitkeep

.PHONY: clean clean-debug clean-optimized clean-installers clean-distribution \
        git-pull macos mac/.remote_mac_timestamp delete-arch testfiles test-backup \
		dict_id_gen$(EXE) \
		objdir.linux objdir.windows objdir.mac \
		push-build increment-version $(INSTALLERS)/LICENSE.html  

# builds latest for local OS
genozip-latest$(EXE): 
	@echo "building latest"
	@(cd $(LATEST_SRC); rm -Rf $(OBJDIR) ; if [[ `pwd` = *genozip-latest* ]]; then git reset --hard ; git pull ; $(MAKE) clean ; touch dict_id_gen.h ; $(MAKE) -j; cp $@ ../../genozip/private/releases/genozip-$(version)$(EXE) ; fi )

# currently, I build for conda from my Windows machine so I don't bother supporting other platforms
ifdef Windows

genozip-latest:
	@run-on-wsl.sh GCC=$(WSL_GCC_DEFAULT) make genozip-latest  # make -j doesn't work well on WSL - filesystem clock issues (try: hwclock -s)

# increments minor version, eg. 1.0.1 -> 1.0.2. 
# To increment a major version, manually edit version.h and set minor version to -1 e.g. 1.1.-1 (careful! no newlines or spaces)
# and re-compile so that genozip --version gets updated
# IMPORTANT: the first number in the version indicates the genozip file format version and goes into
# the genozip file header SectionHeaderTxtHeader.genozip_version
#increment-version: $(C_SRCS) $(CONDA_DEVS) $(CONDA_DOCS) $(CONDA_INCS) # note: target name is not "version.h" so this is not invoked during "make all" or "make debug"
increment-version: # note: target name is not "version.h" so this is not invoked during "make all" or "make debug"
	@echo "Incrementing version.h"
	@bash $(SCRIPTS)/increment-version.sh

.archive.tar.gz : $(C_SRCS) $(X86_SRCS) $(ARM_SRCS) $(CONDA_LIC_OBJS) $(CONDA_DEVS) $(CONDA_DOCS) $(CONDA_INCS) ../LICENSE.txt
	@echo Creating github tag genozip-$(version) and archive
	@$(SH_VERIFY_ALL_COMMITTED)
	@git push > /dev/null
	@(git tag genozip-$(version) > /dev/null ; exit 0) # ignore errors if already exists (due to re-run of make)
	@(cd secure; git tag genozip-$(version) > /dev/null ; exit 0)
	@git push origin genozip-$(version) > /dev/null
	@(cd secure; git push origin genozip-$(version) > /dev/null)
	@curl https://github.com/divonlan/genozip/archive/genozip-$(version).tar.gz --ssl-no-revoke --silent --location -o $@ > /dev/null

conda/meta.yaml: conda/meta.template.yaml .archive.tar.gz ../README.md
	@echo "Generating conda/meta.yaml"
	@bash conda/generate_meta.sh > $@

conda/README.md: conda/README.template.md html-to-md.sed ../README.md
	@echo "Generating conda/README.md"
	@bash conda/generate_README.sh > $@

#CONDA_RECIPE_DIR = ../staged-recipes/recipes/genozip # initial stage-recipes step, keeping here for future reference
CONDA_FEEDSTOCK  = ../../genozip-feedstock
CONDA_RECIPE_DIR = $(CONDA_FEEDSTOCK)/recipe

# publish to conda-forge 
conda/.conda-timestamp: conda/meta.yaml conda/README.md conda/build.sh conda/bld.bat 
	@echo "Publishing to conda-forge"
	@$(SH_VERIFY_ALL_COMMITTED)
	@echo " "
	@echo "Copying $^ to conda feedstock"
	@(cd $(CONDA_FEEDSTOCK); git reset --hard; git pull)
	@cp conda/README.md $(CONDA_FEEDSTOCK)
	@cp conda/meta.yaml conda/build.sh conda/bld.bat $(CONDA_RECIPE_DIR)
	@echo "Committing my files to branch genozip on my fork"
	@(cd $(CONDA_FEEDSTOCK); git pull; git commit -m "update" recipe/meta.yaml README.md recipe/build.sh recipe/bld.bat; git push) > /dev/null
	@echo " "
#	@echo "Submitting pull request to conda-forge"
#	@(cd $(CONDA_RECIPE_DIR); git request-pull master https://github.com://conda-forge/genozip-feedstock master)
#	@(cd $(CONDA_RECIPE_DIR); git request-pull master https://github.com/divonlan/genozip-feedstock master)
	@touch $@
	@echo "CONDA: Using a browser:"
# 'into conda-forge:master from divonlan:genozip'"
	@echo "  (1) Go to https://github.com/divonlan/genozip-feedstock/tree/genozip "
	@echo "  (2) Click 'Contribute' -> 'Open pull request' OR 'Create pull request'"
	@echo "  (3) Click 'Resolve conflicts' -> select the top option in each of the 2 conflicts, then commit"

windows/genozip.exe: CFLAGS += $(OPTFLAGS)
windows/genozip.exe: $(OBJS) $(OBJDIR)/distribution.InstallForge.o 
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

windows/genounzip.exe windows/genocat.exe windows/genols.exe: windows/genozip.exe 
	@echo Generating $@
	@cp $< $@

windows/LICENSE.for-installer.txt: genozip$(EXE) version.h
	@echo Generating $@
	@./genozip$(EXE) --license=txt,60 > $@

windows/genozip-installer.ifp: windows/genozip-installer.template.ifp windows/generate-ifp.sh windows/genozip.exe ../LICENSE.txt
	@echo "Generating $@"
	@(cd windows ; bash generate-ifp.sh)

WINDOWS_INSTALLER_OBJS = windows/genozip.exe windows/genounzip.exe windows/genocat.exe windows/genols.exe windows/LICENSE.for-installer.txt ../LICENSE.txt windows/genozip.ico windows/genozip-installer.ifp

# this must be run ONLY has part of "make distribution" or else versions will be out of sync
$(INSTALLERS)/genozip-installer.exe: $(WINDOWS_INSTALLER_OBJS) 
	@echo 'Creating Windows installer'
	@echo 'WINDOWS: Using the UI:'
	@echo '  (1) Open windows/genozip-installer.ifp'
	@echo '  (2) Click Build'
	@echo '  (3) Exit the UI (close the window)'
	@rm -f ../installers/genozip-installer.exe
	@("../../../../Program Files (x86)/solicus/InstallForge/bin/ifbuilderenvx86.exe" ; exit 0)

$(INSTALLERS)/genozip-linux-x86_64.tar.build: genozip-linux-x86_64/LICENSE.txt 
	@(mkdir genozip-linux-x86_64 >& /dev/null ; exit 0)
	@run-on-wsl.sh GCC=$(WSL_GCC_8_5) make $(INSTALLERS)/genozip-linux-x86_64.tar # make -j doesn't work well on WSL - filesystem clock issues

mac/.remote_mac_timestamp: # to be run from Windows to build on a remote mac
	@echo "Creating Mac installer"
	@echo "Logging in to remote mac" 
	@# Get IP address - check if the previous IP address still works or ask for a new one. Assuming a LAN on an Android hotspot.
	@ip=`cat mac/.mac_ip_address` ; a=`echo $$ip|cut -d. -f4`; (( `ping  -n 1 $$ip | grep "round trip times" | wc -l` > 0 )) || read -p "IP Address: 192.168.43." a ; ip=192.168.43.$$a ; echo $$ip > mac/.mac_ip_address
	@[ -f mac/.mac_username ] || ( echo ERROR: file mac/.mac_username missing && exit 1 )
	@ssh `cat mac/.mac_ip_address` -l `cat mac/.mac_username`  "cd genozip ; echo "Pulling from git" ; git pull >& /dev/null ; make -j mac/.from_remote_timestamp" # pull before make as Makefile might have to be pulled
	@touch $@

$(INSTALLERS)/LICENSE.html: genozip$(EXE)
	@./genozip$(EXE) --license=html,5000 > $@

BUILD_FILES = version.h Makefile \
	objdir.linux/secure/license.o objdir.windows/secure/license.o objdir.mac/secure/license.o \
	../RELEASE_NOTES.txt ../LICENSE.txt \
	$(INSTALLERS)/genozip-installer.exe $(INSTALLERS)/genozip-linux-x86_64.tar $(INSTALLERS)/LICENSE.html

push-build: 
	@git stage $(BUILD_FILES)
	@git commit -m $(version) 
	@git push > /dev/null

distribution-maintenance: clean-distribution increment-version testfiles ../LICENSE.txt $(INSTALLERS)/LICENSE.html \
                          $(INSTALLERS)/genozip-linux-x86_64.tar.build $(INSTALLERS)/genozip-installer.exe objdir.mac/secure/license.o \
                          push-build conda/.conda-timestamp genozip-latest genozip-latest.exe
	@(printf "\nClick Enter after committing the feedstock conflict resolution: <Click>")
	@(read x)
	@(cd $(CONDA_FEEDSTOCK) ; git pull)
	@echo "  Next:"
	@echo "  (1) Wait for all 4 tests to pass (~8 minutes)"
	@echo "  (2) Click 'Merge pull request' and the 'Confirm merge' (DONT CLICK 'Delete branch')"
	@echo "  (3) Open: https://github.com/divonlan/genozip/releases/new"
	@echo "  (4) Set 'Tag version' and 'Release title' are both: genozip-$(version)"
	@echo "  (5) Copy the notes for the version from RELEASE NOTES -> Publish release"

test-backup: genozip.exe
	@echo "Compressing test/ files for in preparation for backup (except cram and bcf)"
	@rm -f test/*.genozip
	@(cd test; genozip.exe -f `ls -1d *|grep -v / |grep -v cram | grep -v bcf`)

# license copied on Windows, not Linux due to file mode issues on NTFS causing git to think LICENSE.txt has changed
genozip-linux-x86_64/LICENSE.txt: ../LICENSE.txt
	@echo Generating $@
	@(mkdir genozip-linux-x86_64 >& /dev/null ; exit 0)
	@cp -f $< $@

endif # Windows

ifdef Linux

# note: getpwuid and getgrgid will cause dymanically loading of the locally installed glibc in the --tar option, or segfault. that's normally fine.
genozip-linux-x86_64/genozip: CFLAGS += $(OPTFLAGS)
genozip-linux-x86_64/genozip: $(OBJS) $(OBJDIR)/distribution.linux-x86_64.o
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

genozip-linux-x86_64/genounzip genozip-linux-x86_64/genocat genozip-linux-x86_64/genols: genozip-linux-x86_64/genozip 
	@echo Generating $@
	@ln -f $< $@

genozip-linux-x86_64/LICENSE.txt: ../LICENSE.txt
	@cp -f $< $@

LINUX_TARGZ_OBJS = genozip-linux-x86_64/genozip genozip-linux-x86_64/genounzip genozip-linux-x86_64/genocat genozip-linux-x86_64/genols genozip-linux-x86_64/LICENSE.txt

# this must be run ONLY as part of "make distribution" or else versions will be out of sync
# this is run on Linux by run-on-wsl.sh
$(INSTALLERS)/genozip-linux-x86_64.tar: $(LINUX_TARGZ_OBJS) 
	@echo "Creating $@"
	@tar cf $@ -z genozip-linux-x86_64

endif # Linux

ifdef Mac

MACDWNDIR = mac/darwinpkg
MACLIBDIR = $(MACDWNDIR)/Library/genozip
MACSCTDIR = mac/scripts
MACRSSDIR = mac/Resources

$(MACRSSDIR)/welcome.html: mac/welcome.template.html
	@sed -e "s/__VERSION__/$(version)/g" $< > $@ 

$(MACSCTDIR)/postinstall: mac/postinstall.template.sh
	@sed -e "s/__FILES__/$(EXECUTABLES)/g" $< > $@

$(MACLIBDIR)/uninstall.sh: mac/uninstall.template.sh
	@sed -e "s/__VERSION__/$(version)/g" $< | sed -e "s/__FILES__/$(EXECUTABLES)/g" > $@ 

$(MACRSSDIR)/README.html: ../README.md
	@cp -f $< $@

$(MACRSSDIR)/LICENSE.txt: ../LICENSE.txt
	@cp -f $< $@

$(MACLIBDIR)/%: %
	@cp -f $< $@

pkg_identifier  := genozip-$(version)

#app_specific_pw := $(shell cat .altool_app_specifc_password)

apple_id        := $(shell /usr/libexec/PlistBuddy -c "print :Accounts:0:AccountID" ~/Library/Preferences/MobileMeAccounts.plist)

signer_name     := ${shell security find-identity -v|grep "3rd Party Mac Developer Installer" | cut -d ":" -f2 | cut -d"(" -f1 }

mac/genozip.pkg: $(MACLIBDIR)/genozip $(MACLIBDIR)/genounzip $(MACLIBDIR)/genocat $(MACLIBDIR)/genols $(MACLIBDIR)/uninstall.sh \
                 $(MACSCTDIR)/postinstall
	@echo "Building Mac package $@"
	@$(SH_VERIFY_ALL_COMMITTED)
	@chmod -R 755 $(MACLIBDIR) $(MACSCTDIR)				 
	@pkgbuild --identifier $(pkg_identifier) --version $(version) --scripts $(MACSCTDIR) --root $(MACDWNDIR) mac/genozip.pkg > /dev/null

mac/genozip_installer.unsigned.pkg: mac/genozip.pkg mac/Distribution \
                                    $(MACRSSDIR)/welcome.html $(MACRSSDIR)/README.html $(MACRSSDIR)/LICENSE.txt
	@echo "Building Mac product $@"
	@productbuild --distribution mac/Distribution --resources $(MACRSSDIR) --package-path mac $@ > /dev/null

mac/genozip_installer.pkg: mac/genozip_installer.unsigned.pkg
	@bash -c "echo -n Unlocking the keychain. Password:"
	@# note: keychain requires unlocking if logged in remotely (through SSH)
	@read pw ; security -v unlock-keychain -p $$pw `security list-keychains|grep login|cut -d\" -f2`
	@echo "Signing Mac product $@"
	@# note: productsign needs a "3rd party mac developer" certificate, and the Apple developer CA certificate, installed in the keychain. see: https://developer.apple.com/developer-id/. I keep them on Drive for backup.
	@productsign --sign "$(signer_name)" $< $@ > /dev/null
	@echo "Verifying the signature"
	@(( `pkgutil --check-signature $@ | grep "signed by a developer certificate issued by Apple (Development)" | wc -l ` > 0 )) || (echo Error: signature verification failed ; exit 1)
	@#@echo 'Committing Mac installer and pushing to repo'
	@#@(git stage $@ ; exit 0)
	@#(git commit -m mac_installer_for_version_$(version) $@ ; exit 0)
	@#git push

mac/.from_remote_timestamp: mac/genozip_installer.pkg
	@echo "Notarizing Mac app"
	@xcrun altool --notarize-app --primary-bundle-id $(pkg_identifier) --username $(apple_id) --password $(app_specific_pw) --file $< >& .notarize.out ; exit 0
	@grep ERROR\: .notarize.out ; exit 0
	@(( `grep ERROR\: .notarize.out | wc -l` == 0 )) || (echo "See .notarize.out" ; exit 1)
	@(( `grep "No errors uploading" .notarize.out | wc -l` == 0 )) || (grep "RequestUUID" .notarize.out ; exit 0)
	@touch $@

endif # Darwin
