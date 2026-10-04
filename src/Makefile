# ------------------------------------------------------------------
#   Makefile
#   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
#   Please see terms and conditions in the file LICENSE.txt
#
#   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
#   and subject to penalties specified in the license.

# Note for Windows: to run this make, you need mingw: https://code.visualstudio.com/docs/cpp/config-mingw

OBJDIR_LINUX   := objdir.linux
OBJDIR_WINDOWS := objdir.windows
OBJDIR_OSX_ARM := objdir.osx-arm

# macros
define verify_dir
	$(if $(strip $($(1))),,$(error $(1) is not defined))
	$(if $(wildcard $($(1))/.),,$(error directory $(1)='$($(1))' does not exist))
endef

define set_path
	$(eval $(1) := $(subst \,/,$(strip $(2))))
	$(call verify_dir,$(1))
	$(info ℹ️ $(1)='$($(1))')
endef

#------------------------------------------------------------------------------------------
# Set architecture, conda, one of: { win_x86, linux_x86, osx_arm64 } and: { arm64, x86_64 }
#------------------------------------------------------------------------------------------
ifeq ($(OS),Windows_NT)
	architecture := Windows $(shell uname -m)
else
	architecture := $(shell uname -sm | sed -e 's/amd/x86_/')
	architecture := $(subst aarch64,arm64,$(architecture))
endif

ifeq ($(architecture),Windows x86_64)	
	win_x86   := 1
	x86_64    := 1

else ifeq ($(architecture),Linux x86_64)
	linux_x86 := 1
	x86_64    := 1

else ifeq ($(architecture),Darwin arm64)
	osx_arm64 := 1		
	arm64     := 1

else
	fail := $(error Architecture unsupported or could not be identified: architecture="$(architecture)" OS="$(OS)") 
endif

ifdef BUILD_PREFIX # conda-forge
	conda     := 1
endif

#--------------------------------------------------
# Settings that is only on conda-forge 
#--------------------------------------------------
ifdef conda 
	ifdef win_x86
		CC  ?= gcc 
		CXX ?= g++

		LDFLAGS += -L$(PREFIX)/Library/lib
		CFLAGS  += -I$(PREFIX)/Library/include 
	endif
	
	DISTOBJ := distribution.conda.o # note: matches the name in version.c 

#--------------------------------------------------
# Settings if non-conda (local build)
#--------------------------------------------------
else 	
	ifdef DISTRIBUTION # custom destribution, run eg: DISTRIBUTION=special make genozip
		DISTOBJ := distribution.$(DISTRIBUTION).o 
	else ifeq ($(USER),divon)
		DISTOBJ := distribution.devel.o
		vendor  := 1 
	else ifeq ($(USERNAME),divon) # check USER and USERNAME separately - in some Windows environments both are set
		DISTOBJ := distribution.devel.o
		vendor  := 1 
	else
		DISTOBJ := distribution.github.o # note: matches the name in version.c 
	endif


	ifdef vendor
        ifndef GENOZIP_HOME
            GENOZIP_HOME := $(HOME)/genozip
		endif
        $(call set_path,GENOZIP_HOME,$(GENOZIP_HOME))
        $(call set_path,LATEST_SRC,$(GENOZIP_HOME)/../genozip-latest/src)
        $(call set_path,TEST_DIR,$(GENOZIP_HOME)/private/test)
        $(call set_path,RELEASES,$(GENOZIP_HOME)/private/releases)
        $(call set_path,INSTALLERS,../installers) # must be relative path
		
		ifdef linux_x86 # Linux but not conda
			# fast, multi-threaded linker, comment out if needed
			LDFLAGS += -fuse-ld=lld

			# need old gcc linked to an old sysroot(=glibc) for older Linux machines
			# up to 15.0.79 used distribution-gcc.8.5 (generate with: conda create -n distribution gcc_linux-64=8.5.0)
			# starting 15.0.80 distribution-gcc.12.4-sysroot.2.17 
			# starting 15.0.88 distribution-gcc.16.1-sysroot.2.17 (see "development environment.txt")
			# To see the glibc version genozip requires: 'objdump -T genozip | grep GLIBC_' (the highest version listed)
			GCC_BIN := /home_local_fs/divon/miniforge3/envs/distribution_gcc-16.1_sysroot-2.17/bin
            $(call verify_dir,GCC_BIN)

			CC  := $(GCC_BIN)/x86_64-conda-linux-gnu-gcc
			CXX := $(GCC_BIN)/x86_64-conda-linux-gnu-g++

		else ifdef osx_arm64 # Mac but not conda
			CFLAGS += -DDEFAULT_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX13.sdk 
		endif
	endif

	# not vendor, or not linux
	CC  ?= gcc
	CXX ?= g++
	
endif 

#--------------------------------------------------------
# Flags that are on all platforms and regardless of conda
#--------------------------------------------------------
LDFLAGS    += -lpthread -lm 
CFLAGS     += -Wall -Wstrict-aliasing=2 -D_LARGEFILE64_SOURCE=1 -D_GNU_SOURCE -std=gnu11 -D_7ZIP_ST -fms-extensions
CXXFLAGS   += -Wall -D_LARGEFILE64_SOURCE=1 -D_GNU_SOURCE -std=gnu++11 -Ibsc 
OPTFLAGS   += -O3 
DEBUGFLAGS += -DDEBUG -ggdb3 -O0 -fstack-check -fstack-protector-strong

$(info ℹ️ Compiler version=$(shell $(CC) --version | head -1))
ifdef linux_x86
    $(info ℹ️ Target sysroot/glibc=$(shell $(CC) -dM -E -include features.h - < /dev/null | egrep -w "__GLIBC__|__GLIBC_MINOR__" | grep -v __GLIBC_PREREQ | awk '{printf "%s%s", sep, $$3; sep="."}' ))
endif
$(info ℹ️ make version=$(shell $(MAKE) --version | head -1))
$(info ℹ️ tar version=$(shell tar --version | head -1))

#--------------------------------
# OS-specific stuff
#--------------------------------
ifdef win_x86
	OBJDIR  := $(OBJDIR_WINDOWS)
	EXE     := .exe
	LDFLAGS := -static -static-libgcc -static-libstdc++ $(LDFLAGS) -lIPHLPAPI -lUrlmon -lDbgHelp # link MingW DLLs staticly (order is important!)
	CFLAGS += -D_FILE_OFFSET_BITS=64

	WSL     := wsl
	NASM    := nasm
	ASMC    := run-on-wsl.sh asmc
	MASM    := $(ASMC) -c -10 -win64

	NASM_FLAGS := -f win64 -DAS_FEATURE_LEVEL=10 -DHAVE_AS_KNOWS_AVX512 # note: modern nasm/yasm have AS_FEATURE_LEVEL=10. HAVE_AS_KNOWS_AVX512 is true if level 10.

	# Use GCC, not Windows, bitfields: allows packing together of members of different types into a bitfield
	CFLAGS  += -mno-ms-bitfields 

	# note: -lasan needed for -fsanitize=address is not supported by MinGW: https://stackoverflow.com/questions/55018627/cannot-find-lasan-using-address-sanitizer-in-mingw-in-windows-mingw

	ifeq ($(CC),cl) # Microsoft Visual C
		fail := $(error Only the gcc and clang compilers are currently supported)
	endif

    $(info ℹ️ NASM version=$(shell $(NASM) --version | head -1))

else ifdef linux_x86
	OBJDIR  := $(OBJDIR_LINUX)
	NASM    := nasm
	ASMC    := asmc
	MASM    := $(ASMC) -c -10 -elf64 -DABI_LINUX 
	
	NASM_FLAGS := -f elf64 -DINTEL_CET_ENABLED -DAS_FEATURE_LEVEL=10 -DHAVE_AS_KNOWS_AVX512
	
	LDFLAGS += -lrt -rdynamic # rt is required by pthreads ; rdynamic for backtrace_symbols_fd to show function names

	# DEBUGFLAGS += -fsanitize=undefined 
	# DEBUGFLAGS += -fsanitize=thread -Dsanitize_thread=1
	
	# Notes on Address Sanitizer:
	# 1. not compatable with valgrind - comment out if using valgrind 
	# 2. Run with: setarch x86_64 -R genozip-debug ... (see internal-docs/address-sanitizer)
	# DEBUGFLAGS += -fsanitize=address -Dsanitize_address=1 

    $(info ℹ️ NASM version=$(shell $(NASM) --version | head -1))

else ifdef osx_arm64
	OBJDIR  := $(OBJDIR_OSX_ARM)

	CFLAGS     += -march=armv8-a+crc
	DEBUGFLAGS += -fsanitize=address -Dsanitize_address -fsanitize=undefined
endif

# check if asmc is available (eg it's not available on conda-forge): we'll use pre-compiled objects it not
ifdef MASM
ifneq ($(shell $(ASMC) -h >/dev/null 2>&1 && echo yes),yes)
    undefine MASM
    $(info ℹ️ MASM version=N/A)
else
    $(info ℹ️ MASM version=$(shell $(MASM) -h | head -1))
endif
endif

SRC_DIRS := \
	secure zlib bzlib lzma bsc libdeflate_1.7 libdeflate_1.19 libdeflate_1.19/x86 libdeflate_1.19/arm 	\
	htscodecs igzip igzip/aarch64 igzip/x86_64 igzip/noarch

C_SRCS := \
	genozip.c genols.c context.c context_stats.c container.c strings.c arch.c tip.c seg_id.c dyn_int.c	\
	data_types.c bits.c progress.c writer.c zriter.c tar.c chrom.c qname.c tokenizer.c mutex.c 			\
	threads.c error.c context_validate.c crc64.c stats.c 												\
	zip.c piz.c reconstruct.c recon_history.c recon_peek.c seg.c zfile.c aligner.c flags.c specials.c   \
	reference.c contigs.c ref_lock.c refhash_make.c refhash_load.c ref_make.c ref_contigs.c 			\
	ref_iupacs.c ref_cache.c digest.c 																	\
	vcf_piz.c vcf_seg.c vcf_vblock.c vcf_header.c vcf_info.c vcf_samples.c vcf_hgvs.c vcf_modify.c     	\
	vcf_format_GT.c vcf_format_PS_PID.c vcf_dbsnp.c vcf_giab.c vcf_vep.c vcf_qual.c vcf_1000G.c 		\
	vcf_me.c vcf_refalt.c vcf_format.c vcf_illum_gtyping.c vcf_gwas.c vcf_vagrent.c vcf_svaba.c 		\
	vcf_pbsv.c vcf_icgc.c vcf_snpeff.c vcf_cosmic.c vcf_mastermind.c vcf_isaac.c vcf_manta.c vcf_pos.c	\
	vcf_ultima.c vcf_platypus.c vcf_info_AC_AF_AN.c vcf_format_GQ.c vcf_gatk.c vcf_sv.c vcf_gnomad.c 	\
	vcf_freebayes.c vcf_local_alleles.c vcf_copy_sample.c vcf_info_DP.c									\
	sam_seg.c sam_piz.c sam_shared.c sam_header.c sam_md.c sam_nm.c sam_tlen.c sam_cigar.c sam_fields.c	\
	sam_sa.c bam_seg.c bam_show.c sam_pacbio.c sam_ultima.c sam_xcons.c cram.c agilent.c sam_hisat2.c	\
	sam_seq.c sam_qual.c sam_sag_zip.c sam_sag_piz.c sam_sag_load.c sam_sag_ingest.c sam_sag_scan.c    	\
	sam_bwa.c sam_bowtie2.c sam_bsseeker2.c sam_bsbolt.c sam_bismark.c sam_gem3.c sam_tmap.c 			\
	sam_blasr.c sam_dragen.c sam_10xGenomics.c sam_biobambam.c sam_pos.c sam_deep.c sam_cpu.c bai.c		\
	sam_star.c sam_abra2.c sam_minimap2.c sam_modify.c sam_recon_plan_piz.c sam_parse.c					\
	gencomp_zip.c gencomp_piz.c																			\
	fastq.c fastq_desc.c fastq_seq.c fastq_qual.c fastq_deep.c fastq_saux.c fastq_bamass.c fastq_aux.c  \
	fastq_parse.c fastq_10xGen.c deep.c																	\
	fasta.c gff.c bed.c me23.c locs.c generic.c lookback.c compressor.c base64.c						\
	buffer.c buf_struct.c buf_list.c random_access.c sections.c sections_show.c 						\
	codec.c codec_bz2.c codec_lzma.c codec_acgt.c codec_domq.c codec_bsc.c codec_pacb.c codec_tmpl.c	\
	codec_pbwt.c codec_none.c codec_htscodecs.c codec_longr.c codec_normq.c codec_homp.c codec_t0.c		\
	codec_smux.c codec_oq.c mgzip.c mgzip_sections.c													\
	txtfile.c profiler.c file.c filename.c dispatcher.c crypt.c aes.c md5.c segconf.c biopsy.c 			\
	vblock.c regions.c dict_id.c aliases.c hash.c stream.c url.c bases_filter.c dict_io.c				\
	version.c huffman.c user_message.c b250.c qname_filter.c txtheader.c 								\
	\
	zlib/zutil.c zlib/deflate.c zlib/trees.c 															\
	\
	bzlib/blocksort.c bzlib/bzlib.c bzlib/compress.c bzlib/crctable.c bzlib/decompress.c 				\
	bzlib/huffman.c bzlib/randtable.c 																	\
	\
	lzma/LzmaEnc.c lzma/LzmaDec.c lzma/LzFind.c lzma/CpuArch.c 											\
	\
	bsc/libsais.c 																						\
	\
	htscodecs/rANS_static4x16pr.c htscodecs/rle.c htscodecs/pack.c htscodecs/arith_dynamic.c 			\
	\
	libdeflate_1.7/deflate_compress.c libdeflate_1.7/utils.c  											\
	\
	libdeflate_1.19/deflate_compress.c libdeflate_1.19/deflate_decompress.c 							\
	libdeflate_1.19/gzip_decompress.c libdeflate_1.19/utils.c libdeflate_1.19/x86/cpu_features.c 		\
	libdeflate_1.19/arm/cpu_features.c libdeflate_1.19/crc32.c libdeflate_1.19/adler32.c 				\
	\
 	xxhash/xxhash.c 																					\
	\
	igzip/igzip.c igzip/hufftables_c.c igzip/igzip_base.c igzip/igzip_icf_base.c igzip/adler32_base.c 	\
	igzip/flatten_ll.c igzip/encode_df.c igzip/igzip_icf_body.c igzip/crc_base.c igzip/huff_codes.c 	\
	igzip/igzip_inflate.c 

CXX_SRCS := \
	bsc/bwt.cpp bsc/coder.cpp bsc/detectors.cpp bsc/libbsc.cpp bsc/lzp.cpp bsc/platform.cpp 			\
    bsc/preprocessing.cpp bsc/qlfc.cpp bsc/qlfc_model.cpp bsc/adler32.cpp 

# others sources that are not listed in OBJS need for remote builds (conda and mac)
OTHER_SRCS := \
	dict_id_gen.h aes.h dispatcher.h profiler.h dict_id.h aliases.h txtfile.h zip.h bits.h progress.h 	\
	website.h endianness.h md5.h sections.h strings.h hash.h stream.h url.h flags.h segconf.h biopsy.h 	\
	huffman.h sorter.h buffer.h buf_struct.h buf_list.h file.h context.h context_struct.h container.h 	\
	seg.h text_license.h version.h compressor.h crypt.h genozip.h piz.h vblock.h zfile.h 				\
	random_access.h regions.h reconstruct.h tar.h qname.h qname_flavors.h codec.h lookback.h 			\
	tokenizer.h codec_longr_alg.c gencomp.h dict_io.h tip.h deep.h filename.h stats.h multiplexer.h 	\
	local_type.h reference.h ref_private.h refhash.h refhash_friend.h ref_iupacs.h aligner.h mutex.h 	\
	mgzip.h threads.h error.h arch.h license.h file_types.h data_types.h base64.h txtheader.h writer.h 	\
	writer_private.h zriter.h bases_filter.h genols.h contigs.h chrom.h vcf.h vcf_private.h sam.h 		\
	sam_private.h sam_friend.h me23.h fasta.h fasta_friend.h gff.h bed.h locs.h generic.h fastq.h 		\
	fastq_private.h user_message.h mac_compat.h b250.h dyn_int.h qname_filter.h bai.h pointer_rel.h			\
	flags_assign.h aligner_layered.h 																	\
	\
	zlib/gzguts.h zlib/zconf.h zlib/deflate.h zlib/trees.h zlib/zlib.h zlib/zutil.h						\
	\
	lzma/7zTypes.h lzma/LzFind.h lzma/LzHash.h lzma/LzmaDec.h lzma/LzmaEnc.h lzma/Precomp.h 			\
	lzma/Compiler.h lzma/CpuArch.h lzma/7zWindows.h lzma/7zAsm.masm lzma/7zAsm.S						\
	\
	bzlib/bzlib.h bzlib/bzlib_private.h																	\
	\
	htscodecs/rANS_static4x16.h htscodecs/rle.h htscodecs/pack.h htscodecs/arith_dynamic.h 				\
	htscodecs/c_simple_model.h htscodecs/rANS_word.h htscodecs/htscodecs_endian.h htscodecs/rANS_word.h \
	htscodecs/utils.h htscodecs/varint.h htscodecs/varint2.h htscodecs/utils.h							\
	\
	bsc/adler32.h bsc/bwt.h bsc/coder.h bsc/coder_tables.h bsc/filter_tables.h bsc/filters.h 			\
	bsc/libbsc.h bsc/libsais.h bsc/lzp.h bsc/platform.h bsc/predictor.h bsc/qlfc.h bsc/qlfc_model.h 	\
	bsc/rangecoder.h																					\
	\
	libdeflate_1.7/unaligned.h libdeflate_1.7/deflate_compress.h libdeflate_1.7/arm_matchfinder_impl.h 	\
	libdeflate_1.7/deflate_constants.h libdeflate_1.7/bt_matchfinder.h libdeflate_1.7/hc_matchfinder.h 	\
	libdeflate_1.7/common_defs.h libdeflate_1.7/lib_common.h libdeflate_1.7/x86_matchfinder_impl.h 		\
	libdeflate_1.7/compiler_gcc.h libdeflate_1.7/libdeflate.h libdeflate_1.7/matchfinder_common.h		\
	\
	libdeflate_1.19/adler32_vec_template.h libdeflate_1.19/crc32_tables.h 								\
	libdeflate_1.19/crc32_multipliers.h libdeflate_1.19/deflate_constants.h 							\
	libdeflate_1.19/lib_common.h libdeflate_1.19/bt_matchfinder.h libdeflate_1.19/hc_matchfinder.h 		\
	libdeflate_1.19/common_defs.h	libdeflate_1.19/libdeflate.h libdeflate_1.19/gzip_constants.h 		\
	libdeflate_1.19/decompress_template.h libdeflate_1.19/cpu_features_common.h 						\
	libdeflate_1.19/matchfinder_common.h libdeflate_1.19/deflate_compress.h 							\
	libdeflate_1.19/ht_matchfinder.h libdeflate_1.19/x86/adler32_impl.h 								\
	libdeflate_1.19/x86/crc32_impl.h libdeflate_1.19/x86/matchfinder_impl.h 							\
	libdeflate_1.19/x86/cpu_features.h libdeflate_1.19/x86/crc32_pclmul_template.h 						\
	libdeflate_1.19/x86/decompress_impl.h libdeflate_1.19/arm/adler32_impl.h 							\
	libdeflate_1.19/arm/cpu_features.h libdeflate_1.19/arm/crc32_impl.h 								\
	libdeflate_1.19/arm/matchfinder_impl.h libdeflate_1.19/arm/crc32_pmull_helpers.h 					\
	libdeflate_1.19/arm/crc32_pmull_wide.h																\
	\
	igzip/bitbuf2.h igzip/crc.h igzip/encode_df.h igzip/flatten_ll.h igzip/huff_codes.h igzip/huffman.h \
	igzip/igzip_checksums.h igzip/igzip_level_buf_structs.h igzip/igzip_lib.h igzip/igzip_wrapper.h 	\
	igzip/repeated_char_result.h igzip/static_inflate.h igzip/unaligned.h								\
	\
	igzip/x86_64/bitbuf2.asm igzip/x86_64/data_struct2.asm igzip/x86_64/heap_macros.asm 				\
	igzip/x86_64/huffman.asm igzip/x86_64/igzip_compare_types.asm 										\
	igzip/x86_64/igzip_decode_block_stateless.asm igzip/x86_64/igzip_update_histogram.asm 				\
	igzip/x86_64/inflate_data_structs.asm igzip/x86_64/lz0a_const.asm igzip/x86_64/multibinary.asm 		\
	igzip/x86_64/options.asm igzip/x86_64/reg_sizes.asm igzip/x86_64/stdmac.asm							\
	\
	igzip/aarch64/aarch64_label.h igzip/aarch64/aarch64_multibinary.h igzip/aarch64/bitbuf2_aarch64.h 	\
	igzip/aarch64/crc32_aarch64_common.h igzip/aarch64/crc32_gzip_refl_pmull.h 							\
	igzip/aarch64/crc32_norm_common_pmull.h igzip/aarch64/crc32_refl_common_pmull.h 					\
	igzip/aarch64/crc_common_pmull.h igzip/aarch64/data_struct_aarch64.h 								\
	igzip/aarch64/huffman_aarch64.h igzip/aarch64/lz0a_const_aarch64.h igzip/aarch64/options_aarch64.h	\
	igzip/aarch64/stdmac_aarch64.h igzip/aarch64/crc32_common_crc_ext_cortex_a72.S 						\
	igzip/aarch64/crc32_common_mix_neoverse_n1.S														\
	\
	xxhash/xxhash.h xxhash/xxh_x86dispatch.h

#--------------------------------------------
# CPU-specific source files compilation flags
#--------------------------------------------
ifdef x86_64
	C_SRCS += \
		xxhash/xxh_x86dispatch.c
	
	NASM_SRCS := \
		igzip/x86_64/igzip_body.asm igzip/x86_64/igzip_finish.asm igzip/x86_64/igzip_icf_body_h1_gr_bt.asm 					\
		igzip/x86_64/igzip_icf_finish.asm igzip/x86_64/rfc1951_lookup.asm igzip/x86_64/adler32_sse.asm 									\
		igzip/x86_64/adler32_avx2_4.asm igzip/x86_64/igzip_multibinary.asm igzip/x86_64/igzip_update_histogram_01.asm 					\
		igzip/x86_64/igzip_update_histogram_04.asm igzip/x86_64/igzip_decode_block_stateless_01.asm 										\
		igzip/x86_64/igzip_decode_block_stateless_04.asm igzip/x86_64/igzip_inflate_multibinary.asm 										\
		igzip/x86_64/encode_df_04.asm igzip/x86_64/encode_df_06.asm igzip/x86_64/proc_heap.asm igzip/x86_64/igzip_deflate_hash.asm 		\
		igzip/x86_64/igzip_gen_icf_map_lh1_06.asm igzip/x86_64/igzip_gen_icf_map_lh1_04.asm igzip/x86_64/igzip_set_long_icf_fg_04.asm 	\
		igzip/x86_64/igzip_set_long_icf_fg_06.asm igzip/x86_64/crc_multibinary.asm igzip/x86_64/crc32_gzip_refl_by8.asm 					\
		igzip/x86_64/crc32_gzip_refl_by8_02.asm igzip/x86_64/crc32_gzip_refl_by16_10.asm

	ifdef MASM # asmc available
		MASM_SRCS := lzma/LzmaDecOpt.masm 
		CFLAGS += -DZ7_LZMA_DEC_OPT

	else # asmc not available (e.g. on conda-forge): use pre-assembled objects
		LDFLAGS += $(OBJDIR)/lzma/LzmaDecOpt.o
	endif

	CFLAGS     += -march=haswell 

	CXXFLAGS   += -DLIBBSC_ALLOW_UNALIGNED_ACCESS -march=haswell -DLIBBSC_CPU_FEATURE=LIBBSC_CPU_FEATURE_AVX2 

	NASM_FLAGS += -Iigzip/x86_64

else ifdef osx_arm64
	C_SRCS += 			\
		mac_compat.c	\
		\
		igzip/aarch64/crc_aarch64_dispatcher.c igzip/aarch64/igzip_multibinary_aarch64_dispatcher.c igzip/proc_heap_base.c

ifdef conda # we can't yet compile igzip .S on conda - appears to be some Assembler syntax difference between the versions of clang
	C_SRCS += \
		igzip/noarch/igzip_base_aliases.c igzip/proc_heap_base.c igzip/noarch/crc_base_aliases.c

	CFLAGS += -Iigzip/noarch
else
	S_SRCS := \
		igzip/aarch64/igzip_inflate_multibinary_arm64.S igzip/aarch64/igzip_multibinary_arm64.S 								\
		igzip/aarch64/igzip_isal_adler32_neon.S igzip/aarch64/igzip_deflate_body_aarch64.S igzip/aarch64/igzip_deflate_finish_aarch64.S	\
		igzip/aarch64/isal_deflate_icf_body_hash_hist.S igzip/aarch64/isal_deflate_icf_finish_hash_hist.S									\
		igzip/aarch64/igzip_set_long_icf_fg.S igzip/aarch64/encode_df.S igzip/aarch64/isal_update_histogram.S igzip/aarch64/gen_icf_map.S \
		igzip/aarch64/igzip_deflate_hash_aarch64.S igzip/aarch64/igzip_decode_huffman_code_block_aarch64.S 								\
		igzip/aarch64/crc32_gzip_refl_3crc_fold.S igzip/aarch64/crc32_gzip_refl_crc_ext.S igzip/aarch64/crc32_gzip_refl_pmull.S			\
		igzip/aarch64/crc32_mix_default.S igzip/aarch64/crc32_mix_default_common.S igzip/aarch64/crc32_mix_neoverse_n1.S 					\
		igzip/aarch64/crc32c_mix_default.S igzip/aarch64/crc32c_mix_neoverse_n1.S igzip/aarch64/crc_multibinary_arm.S 
endif

	# Luckicly, we can always compile LZMA .S
	S_SRCS   += lzma/LzmaDecOpt.S 
	CFLAGS   += -DZ7_LZMA_DEC_OPT

	CXXFLAGS += -DLIBBSC_ALLOW_UNALIGNED_ACCESS  
	
	S_FLAGS  := -Iigzip/aarch64

else # another CPU
	C_SRCS += \
		igzip/noarch/igzip_base_aliases.c igzip/proc_heap_base.c igzip/noarch/crc_base_aliases.c

	CFLAGS += -Iigzip/noarch
endif

#-------------------------------------------------------------
# secure sources
#-------------------------------------------------------------
ifneq ($(wildcard secure/.*),) # case: "secure" directory exists
	C_SRCS += secure/license.c
else # case: "secure" directory doesn't exists: conda, or user is compiling from github
	LDFLAGS += $(OBJDIR)/secure/license.o
endif

#-------------------------------------------------------------
# list of all source object, dependencies and executable files
#-------------------------------------------------------------
# object files that are shared between optimized/debug/opt builds
SHA_OBJS          := $(addprefix $(OBJDIR)/,$(MASM_SRCS:.masm=.o)) $(addprefix $(OBJDIR)/,$(NASM_SRCS:.asm=.o))  $(addprefix $(OBJDIR)/,$(S_SRCS:.S=.o))
DEPS              := $(addprefix $(OBJDIR)/,$(C_SRCS:.c=.d))       $(addprefix $(OBJDIR)/,$(CXX_SRCS:.cpp=.d))  $(addprefix $(OBJDIR)/,$(NASM_SRCS:.asm=.d)) $(addprefix $(OBJDIR)/,$(S_SRCS:.S=.d)) 
# object files for installers, distribution.c will be added later
INSTALLER_OBJS    := $(addprefix $(OBJDIR)/,$(C_SRCS:.c=.o))       $(addprefix $(OBJDIR)/,$(CXX_SRCS:.cpp=.o))       $(SHA_OBJS)
OBJS              := $(INSTALLER_OBJS) $(OBJDIR)/$(DISTOBJ)
DEBUG_OBJS        := $(addprefix $(OBJDIR)/,$(C_SRCS:.c=.debug-o)) $(addprefix $(OBJDIR)/,$(CXX_SRCS:.cpp=.debug-o)) $(SHA_OBJS) $(OBJDIR)/$(DISTOBJ)
# optimized but with debug info, for debugging issues that only manifest with compiler optimization
OPT_OBJS          := $(addprefix $(OBJDIR)/,$(C_SRCS:.c=.opt-o))   $(addprefix $(OBJDIR)/,$(CXX_SRCS:.cpp=.opt-o))   $(SHA_OBJS) $(OBJDIR)/$(DISTOBJ)

EXECUTABLES       := genozip$(EXE)       genounzip$(EXE)       genocat$(EXE)       genols$(EXE)       
DEBUG_EXECUTABLES := genozip-debug$(EXE) genounzip-debug$(EXE) genocat-debug$(EXE) genols-debug$(EXE) 
OPT_EXECUTABLES   := genozip-opt$(EXE)   genounzip-opt$(EXE)   genocat-opt$(EXE)   genols-opt$(EXE)   

all   : CFLAGS += $(OPTFLAGS)
all   : CXXFLAGS += $(OPTFLAGS) 
all   : $(OBJDIR) $(EXECUTABLES) 
ifndef win_x86
	@chmod +x test.sh
endif

# ifdef win_x86
# debug opt : LDFLAGS += -lbacktrace
# endif

debug : CFLAGS += $(DEBUGFLAGS) 
debug : CXXFLAGS += $(DEBUGFLAGS) 
debug : $(OBJDIR) $(DEBUG_EXECUTABLES)

opt   : CFLAGS += -ggdb3 $(OPTFLAGS) -DDEBUG
opt   : CXXFLAGS += $(OPTFLAGS) 
opt   : $(OBJDIR) $(OPT_EXECUTABLES)

# no calculating dependencies when make command includes one of these as a target
NO_DEPS_TARGETS := clean clean-debug clean-opt clean-installers clean-dev clean-test \
                   increment-version dict_id_gen.h $(OSX).tar genozip-latest
	
ifeq ($(filter $(NO_DEPS_TARGETS),$(MAKECMDGOALS)),)
-include $(DEPS)
endif

# --------------------------------------------------------------------
# Generated sources: dict_id_gen.h, flags_assign.h, context_validate.c
# --------------------------------------------------------------------
ifndef conda # conda doesn't need to generate - files available on github (+ bash issues generating)
#  generation:
# Step 1: dict_id_gen.sh generates dict_id_gen.c, including all the GENDICT definitions from the data type include files (eg vcf.h)
# Step 2: dict_id_gen.sh compiles dict_id_gen.c and generate dict_id_gen[.exe]
# Step 3. dict_id_gen.sh generates dict_id_gen.h: it uses dict_id_gen[.exe]to generate the field constant, and then adds the fields enum and mapping
GENERATED_SRCS := dict_id_gen.h flags_assign.h context_validate.c

dict_id_gen.h : bed.h fasta.h fastq.h generic.h genozip.h gff.h locs.h me23.h reference.h sam.h vcf.h dict_id_gen.sh
	@echo Generating $@								;\
	bash --noprofile --norc ./dict_id_gen.sh $(CC)

flags_assign.h : flags.h flags_assign.sh
	@echo Generating $@								;\
	bash --noprofile --norc ./flags_assign.sh

context_validate.c: context_validate.sh context_struct.h
	@echo Generating $@								;\
	bash --noprofile --norc ./context_validate.sh > $@
endif 

# -----------------------------------------------------------------------
# Compilation rules: { dependecies, C, C++ } x { optimized, debug, opt }
# -----------------------------------------------------------------------
$(OBJDIR):
	@echo Making directory $@ and its subdirs
	@mkdir -p $@ $(addprefix $@/,$(SRC_DIRS)) 

$(OBJDIR)/%.o: %.c $(OBJDIR)/%.d
	@echo "Compiling ($(architecture)) $<" 			;\
	$(CC) -c -o $@ $< $(CFLAGS) 

$(OBJDIR)/%.debug-o: %.c $(OBJDIR)/%.d
	@echo "Compiling ($(architecture)) $< (debug)" 	;\
	$(CC) -c -o $@ $< $(CFLAGS)

$(OBJDIR)/%.opt-o: %.c $(OBJDIR)/%.d
	@echo "Compiling ($(architecture)) $< (opt)"	;\
	$(CC) -c -o $@ $< $(CFLAGS)

$(OBJDIR)/%.d: %.c $(GENERATED_SRCS) | $(OBJDIR) # directory is an "order only prerequesite": https://www.gnu.org/savannah-checkouts/gnu/make/manual/html_node/Prerequisite-Types.html#Prerequisite-Types
	@echo "Calculating dependencies ($(architecture)) $< ";\
	$(CC) $(CFLAGS) -MM -MT $@ $< -MF $(@:%.o=%.d)

$(OBJDIR)/%.d: %.cpp $(GENERATED_SRCS) | $(OBJDIR) # directory is an "order only prerequesite": https://www.gnu.org/savannah-checkouts/gnu/make/manual/html_node/Prerequisite-Types.html#Prerequisite-Types
	@echo "Calculating dependencies ($(architecture))$< ";\
	$(CXX) $(CXXFLAGS) -MM -MT $@ $< -MF $(@:%.o=%.d)

$(OBJDIR)/distribution.%.o : distribution.c 
	@echo "Compiling ($(architecture)) $< to $@"	;\
	$(CC) -c -o $@ $< $(CFLAGS) -DDISTRIBUTION=\"$*\"

$(OBJDIR)/%.o: %.cpp $(OBJDIR)/%.d
	@echo "Compiling ($(architecture)) $<"			;\
	$(CXX) -c -o $@ $< $(CXXFLAGS) 

$(OBJDIR)/%.debug-o: %.cpp $(OBJDIR)/%.d
	@echo "Compiling ($(architecture)) $< (debug)"	;\
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(OBJDIR)/%.opt-o: %.cpp $(OBJDIR)/%.d
	@echo "Compiling ($(architecture)) $< (opt)"	;\
	$(CXX) -c -o $@ $< $(CXXFLAGS)

# ---------------------------------------------------------------
# Assember rules: three differnet formats: Intel, Microsoft, GNU
# ---------------------------------------------------------------
ifdef x86_64
$(OBJDIR)/%.d: %.asm  # Dependencies: nasm (Intel) format
	@echo "Calculating dependencies ($(architecture)) $< " ;\
	$(NASM) $(NASM_FLAGS) -M $< > $@

$(OBJDIR)/%.o: %.asm $(OBJDIR)/%.d # Assembling: nasm (Intel) format
	@echo "Assembling ($(architecture)) $< " 	;\
	$(NASM) $(NASM_FLAGS) -o $@ $< 

ifdef MASM # asmc available
$(OBJDIR)/%.o: %.masm # Assembling: masm (Microsoft) format (note: asmc doesn't calculate dependencies)
	@echo "Assembling ($(architecture)) $< " 	;\
	$(MASM) -Fo $@ $< > /dev/null 2>&1 			;\
	status=$$?									;\
	src=$<; err=$$(basename $${src%.masm}.err)	;\
	if [ -f $$err ]; then cat $$err > /dev/stderr ; rm -f $$err; fi ;\
	exit $$status
endif

else ifdef osx_arm64
$(OBJDIR)/%.d: %.S    # Dependencies: GNU format
	@echo "Calculating dependencies ($(architecture)) $< ";\
	$(CC) -D__ASSEMBLY__ -M $(S_FLAGS) $< > $@

$(OBJDIR)/%.o: %.S $(OBJDIR)/%.d  # Assembling: GNU format
	@echo "Assembling ($(architecture)) $<" 		;\
	$(CC) -D__ASSEMBLY__ -c $(S_FLAGS) -o $@ $< 

else 
	fail := $(error Assembler flavor not recognized) 
endif

# ---------------------------------------------------------------
# Preprocssing-only (for debugging). Usage, e.g.: make buffer.E
# ---------------------------------------------------------------
%.E: %.c $(OBJDIR)/%.d
	@echo "Generating $@"							;\
	$(CC) -E -o $@ $< $(CFLAGS)

# --------------------------------------------------------------------
# Executables
# --------------------------------------------------------------------
genozip$(EXE): $(GENERATED_SRCS) $(OBJS)
	@echo Linking $@ 								;\
	$(CC) -o $@ $(OBJS) $(CFLAGS) $(LDFLAGS) 		;\
	if ! ln -f $@ ../$@; then ln -fs $@ ../$@; fi  # hard linking not possible on some filesystems
 
genozip-debug$(EXE): $(GENERATED_SRCS) $(DEBUG_OBJS)
	@echo Linking $@ 								;\
	$(CC) -o $@ $(DEBUG_OBJS) $(CFLAGS) $(LDFLAGS)	;\
	if ! ln -f $@ ../$@; then ln -fs $@ ../$@; fi

genozip-opt$(EXE): $(GENERATED_SRCS) $(OPT_OBJS) 
	@echo Linking $@ 								;\
	$(CC) -o $@ $(OPT_OBJS) $(CFLAGS) $(LDFLAGS) 	;\
	if ! ln -f $@ ../$@; then ln -fs $@ ../$@; fi

genounzip$(EXE) genocat$(EXE) genols$(EXE): genozip$(EXE)
	@echo Hard linking $@ 							;\
	rm -f $@ 										;\
	ln -f $^ $@ 									;\
	if ! ln -f $^ $@; then ln -fs $^ $@; fi 		;\
	if ! ln -f $@ ../$@; then ln -fs $@ ../$@; fi

genounzip-debug$(EXE) genocat-debug$(EXE) genols-debug$(EXE): genozip-debug$(EXE)
	@echo Hard linking $@ 							;\
	rm -f $@ 										;\
	if ! ln -f $^ $@; then ln -fs $^ $@; fi 		;\
	if ! ln -f $@ ../$@; then ln -fs $@ ../$@; fi

genounzip-opt$(EXE) genocat-opt$(EXE) genols-opt$(EXE): genozip-opt$(EXE)
	@echo Hard linking $@ 							;\
	rm -f $@ 										;\
	if ! ln -f $^ $@; then ln -fs $^ $@; fi 		;\
	if ! ln -f $@ ../$@; then ln -fs $@ ../$@; fi

test:
	@cat test.sh | tr -d "\r" | bash -

# -----------------------------------------------------------------------------------------------
# Rules for building a new Genozip release: Dependencies are here, recipes are in distribution.mk
# -----------------------------------------------------------------------------------------------
ifdef vendor

# distribution platform names
LINUX := genozip-linux-x86_64
OSX   := genozip-osx-arm
WIN   := windows
CONDA := conda

# macro that is evaluated each time when called
DIST_MAKE = @\
	case "$$(pwd)" in \
        */genozip/src) ;; \
        *) echo "❌ make distribution must run from $(GENOZIP_HOME)/src, but pwd=$$(pwd)"; exit 1 ;; \
    esac; \
	$(MAKE) -f secure/distribution.mk --no-print-directory $@

# pass to distribution.mk :
export CC CXX CFLAGS CXXFLAGS LDFLAGS OPTFLAGS EXE SRC_DIRS C_SRCS CXX_SRCS NASM_SRCS S_SRCS OTHER_SRCS
export OBJDIR OBJDIR_LINUX OBJDIR_WINDOWS OBJDIR_OSX_ARM OBJS DEPS DISTOBJ SCRIPTS LATEST_SRC 
export LINUX OSX WIN CONDA RELEASES
export EXECUTABLES OPT_EXECUTABLES DEBUG_EXECUTABLES INSTALLERS CONDA_LIC_OBJS TEST_DIR
export architecture win_x86 linux_x86 osx_arm64 x86_64 arm64 vendor

# ------------------------------------
# Build rules that run only on Windows
# ------------------------------------
ifdef win_x86

# Windows distribution.mk targets with no dependencies
initialize-build-distribution push-build finalize-distribution genozip-latest.exe genozip-latest: 
	+$(DIST_MAKE) 

#CONDA_RECIPE_DIR := ../staged-recipes/recipes/genozip # initial stage-recipes step, keeping here for future reference
CONDA_FEEDSTOCK  := ../../genozip-feedstock
CONDA_RECIPE_DIR := $(CONDA_FEEDSTOCK)/recipe

OTHER_FILES_FOR_CONDA := \
	Makefile ../.gitignore ../LICENSE.txt ../AUTHORS ../README.md \
	$(OBJDIR_LINUX)/secure/license.o  $(OBJDIR_WINDOWS)/secure/license.o $(OBJDIR_OSX_ARM)/secure/license.o \
	$(OBJDIR_LINUX)/lzma/LzmaDecOpt.o $(OBJDIR_WINDOWS)/lzma/LzmaDecOpt.o

../LICENSE.txt: text_license.h version.h $(WIN)/genozip.exe
	@echo Generating $@ ;\
	$(WIN)/genozip.exe --license=txt,100 > $@ 

increment-version: # note: target name is not "version.h" so this is not invoked during "make all" or "make debug"
	+$(DIST_MAKE) 

.archive.tar.gz : $(C_SRCS) $(NASM_SRCS) $(S_SRCS) $(OTHER_SRCS) $(OTHER_FILES_FOR_CONDA)
	+$(DIST_MAKE) 

$(CONDA)/meta.yaml: $(CONDA)/meta.template.yaml .archive.tar.gz ../README.md
	+$(DIST_MAKE) 

$(CONDA)/README.md: $(CONDA)/README.template.md html-to-md.sed ../README.md
	+$(DIST_MAKE) 

$(CONDA)/.conda-timestamp: $(CONDA)/meta.yaml $(CONDA)/README.md $(CONDA)/build.sh $(CONDA)/bld.bat 
	+$(DIST_MAKE) 

WINDOWS_INSTALLER_OBJS := \
	$(WIN)/genozip.exe $(WIN)/genounzip.exe $(WIN)/genocat.exe $(WIN)/genols.exe \
	$(WIN)/LICENSE.for-installer.txt ../LICENSE.txt $(WIN)/genozip.ico $(WIN)/genozip-installer.ifp

$(WIN)/genozip.exe: CFLAGS += $(OPTFLAGS)
#$(WIN)/genozip.exe: $(INSTALLER_OBJS) $(OBJDIR)/distribution.MicrosoftStore.o 
$(WIN)/genozip.exe: $(INSTALLER_OBJS) $(OBJDIR)/distribution.InstallForge.o 
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

$(WIN)/genounzip.exe $(WIN)/genocat.exe $(WIN)/genols.exe: $(WIN)/genozip.exe 
	+$(DIST_MAKE) 

$(WIN)/LICENSE.for-installer.txt: $(WIN)/genozip.exe version.h
	+$(DIST_MAKE) 

$(WIN)/genozip-installer.ifp: $(WIN)/genozip-installer.template.ifp $(WIN)/generate-ifp.sh $(WIN)/genozip.exe ../LICENSE.txt
	+$(DIST_MAKE) 

$(info ℹ️ WIN='$(WIN)')
$(info ℹ️ INSTALLERS='$(INSTALLERS)')

$(INSTALLERS)/genozip-installer.exe: $(WINDOWS_INSTALLER_OBJS) 
	+$(DIST_MAKE) 

# Windows making LICENSE.txt for Linux
$(LINUX)/LICENSE.txt: ../LICENSE.txt
	+$(DIST_MAKE) 

# Windows rule for making Linux build 
$(INSTALLERS)/$(LINUX).tar.build: $(LINUX)/LICENSE.txt 
	+$(DIST_MAKE) 

$(INSTALLERS)/$(OSX).tar : $(C_SRCS) $(S_SRCS) $(OTHER_SRCS) $(wildcard secure/*.c) $(wildcard secure/*.h) ../LICENSE.txt Makefile
	+$(DIST_MAKE) 

$(INSTALLERS)/LICENSE.html: $(WIN)/genozip.exe
	+$(DIST_MAKE) 

.NOTPARALLEL: distribution 
distribution: 							\
	initialize-build-distribution 		\
	increment-version 					\
	../LICENSE.txt 						\
	$(INSTALLERS)/LICENSE.html 			\
	$(INSTALLERS)/$(OSX).tar 			\
	$(INSTALLERS)/$(LINUX).tar.build	\
	$(INSTALLERS)/genozip-installer.exe	\
	push-build							\
	$(CONDA)/.conda-timestamp 			\
	genozip-latest.exe 					\
	genozip-latest 						\
	finalize-distribution				\

# ------------------------------------
# Build rules that run only on Linux
# ------------------------------------
else ifdef linux_x86

distribution initialize-build-distribution finalize-distribution push-build increment-version $(CONDA)/.conda-timestamp .archive.tar.gz $(CONDA)/meta.yaml:
	@echo "❌ \"make $@\" must be run from Windows"

# Linux distribution.mk targets with no dependencies
genozip-latest: 
	+$(DIST_MAKE) 

# note: getpwuid and getgrgid will cause dymanically loading of the locally installed glibc in the --tar option, or segfault. that's normally fine.
$(LINUX)/genozip: CFLAGS += $(OPTFLAGS)
$(LINUX)/genozip: $(INSTALLER_OBJS) $(OBJDIR)/distribution.linux-x86_64.o
	@echo Linking $@ ;\
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

$(LINUX)/genounzip $(LINUX)/genocat $(LINUX)/genols: $(LINUX)/genozip 
	+$(DIST_MAKE) 

# this must be run ONLY as part of "make distribution" or else versions will be out of sync
# this is run on Linux by run-on-wsl.sh
$(INSTALLERS)/$(LINUX).tar: $(LINUX)/genozip $(LINUX)/genounzip $(LINUX)/genocat $(LINUX)/genols $(LINUX)/LICENSE.txt
	+$(DIST_MAKE) 

# ------------------------------------
# Build rules that run only on Mac
# ------------------------------------
else ifdef osx_arm64

distribution initialize-build-distribution finalize-distribution push-build increment-version $(CONDA)/.conda-timestamp .archive.tar.gz $(CONDA)/meta.yaml:
	@echo "❌ \"make $@\" must be run from Windows"

$(OSX)/genozip: CFLAGS += $(OPTFLAGS)
$(OSX)/genozip: $(INSTALLER_OBJS) $(OBJDIR)/distribution.osx-%.o
	@echo Linking $@ ;\
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

$(OSX)/genounzip $(OSX)/genocat $(OSX)/genols: $(OSX)/genozip 
	+$(DIST_MAKE) 

$(OSX)/LICENSE.txt: ../LICENSE.txt
	+$(DIST_MAKE) 

# Mac distribution.mk targets with no dependencies
genozip-latest: 
	+$(DIST_MAKE) 

$(INSTALLERS)/$(OSX).tar: $(OSX)/genozip $(OSX)/genounzip $(OSX)/genocat $(OSX)/genols $(OSX)/LICENSE.txt
	+$(DIST_MAKE) 

endif # osx-arm64

else # not vendor

distribution: 
	@echo "Target $@ intended for vendor use only. USER=$(USER) USERNAME=$(USERNAME)"

endif # vendor

#-----------------------
# cleaners
#-----------------------

# macro rm of files and showing removed files in one line
define rm
	@/bin/echo -n "$(1) : " ; rm -fv $(2) | cut -d\' -f2 | xargs -r -n1 basename | tr "\n" " "; echo
endef

clean-debug:
	$(call rm,Cleaning up debug,	\
		$(DEBUG_OBJS) $(DEBUG_EXECUTABLES) $(OBJDIR)/*.debug-o)

clean-opt:
	$(call rm,Cleaning up opt,		\
		$(OPT_OBJS) $(OPT_EXECUTABLES) $(OBJDIR)/*.opt-o)

clean:  # careful not to delete license.o and LzmaDecOpt.o (these are not in $(OBJS) if not used)
	$(call rm,Cleaning up optimized，dependencies and generated sources, \
		$(OBJS) $(EXECUTABLES) $(DEPS) $(GENERATED_SRCS) $(OBJDIR)/*.o $(OBJDIR)/*.d)
	
clean-test:
	$(call rm,Cleaning up Test, 	\
		$(TEST_DIR)/*.genozip $(TEST_DIR)/*.bad $(TEST_DIR)/*.local $(TEST_DIR)/*.b250 $(TEST_DIR)/outdir.*/* \
	    $(TEST_DIR)/*.DEPN $(TEST_DIR)/back-compat/*/*.bad $(TEST_DIR)/*.bad-offset-* $(TEST_DIR)/*/*.bad-offset-*
			 
clean-dev: clean-installers clean-test # use only in dev - deletes license.o and LzmaDecOpt.o
	$(call rm,Cleaning dev stuff,													\
		$(DEPS) $(filter-out LICENSE.txt,$(WINDOWS_INSTALLER_OBJS)) .archive.tar.gz \
		$(EXECUTABLES) $(OPT_EXECUTABLES) $(DEBUG_EXECUTABLES) $(OBJDIR)/* $(OBJDIR)/*/* $(OBJDIR)/*/*/*			
	mkdir $(OBJDIR) $(addprefix $(OBJDIR)/,$(SRC_DIRS))	;\
	touch $(OBJDIR)/.gitkeep

clean-installers:
	$(call rm,Cleaning installers, \
		$(INSTALLERS)/*)

.PHONY: clean clean-debug clean-opt dict_id_gen$(EXE) $(OBJDIR_LINUX) $(OBJDIR_WINDOWS) $(OBJDIR_OSX_ARM) \
        distribution genozip-latest genozip-latest.exe

