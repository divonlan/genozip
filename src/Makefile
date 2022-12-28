# ------------------------------------------------------------------
#   Makefile
#   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
#   Please see terms and conditions in the file LICENSE.txt
#
#   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
#   and subject to penalties specified in the license.

# Note for Windows: to run this make, you need mingw (for the gcc compiler) and cygwin (for Unix-like tools):
# Mingw: http://mingw-w64.org/doku.php 
# Cygwin: https://www.cygwin.com/

uname := $(shell uname -s)

ifdef BUILD_PREFIX
IS_CONDA=1
endif

LDFLAGS    += -lpthread -lm 
CFLAGS     += -Wall -D_LARGEFILE64_SOURCE=1 -D_7ZIP_ST -std=gnu11 -fms-extensions 
OPTFLAGS   += -Ofast 
DEBUGFLAGS += -DDEBUG -g -O0 

ifdef IS_CONDA 
	CFLAGS += -DDISTRIBUTION=\"conda\" # note: matches the name in version.c -DIS_CONDA
	CFLAGS += -Ofast 

	ifeq ($(OS),Windows_NT)
		CC=gcc # in Windows, override conda's default Visual C with gcc 
		LDFLAGS += -L$(PREFIX)/Library/lib
		CFLAGS  += -I$(PREFIX)/Library/include 
	endif

else
	CC=gcc

	# set -march. note: if conda, we leave it to conda to set
	ifeq ($(OS),Windows_NT)
		CFLAGS += -march=native 
	else
		ifeq ($(uname),Linux)
			CFLAGS += -march=nocona # similar to conda - relatively low optimization, since different hosts in a cluster might have different capabilities
		endif

		ifeq ($(uname),Darwin) # Mac
			CFLAGS += -march=native # conda used "core2". We use the highest available optimizations, because it is typically compiled for use on a specific Mac
		endif
	endif

endif 

DEV_DIR  = ../../genozip
DEV_SRC  = $(DEV_DIR)/src
TEST_DIR = $(DEV_DIR)/private/test
SCRIPTS  = $(DEV_DIR)/private/scripts
UTILS    = $(DEV_DIR)/private/utils
RELEASES = $(DEV_DIR)/private/releases
DATA     = $(DEV_DIR)/data
PROD_SRC = ../../genozip-prod/src

SRC_DIRS = zlib bzlib lzma bsc libdeflate htscodecs compatibility

MY_SRCS = genozip.c genols.c context.c container.c strings.c stats.c arch.c license.c tip.c      				\
		  data_types.c bits.c progress.c writer.c tar.c chrom.c qname.c tokenizer.c 							\
          zip.c piz.c reconstruct.c recon_history.c recon_peek.c seg.c zfile.c aligner.c flags.c digest.c    	\
		  reference.c contigs.c ref_lock.c refhash.c ref_make.c ref_contigs.c ref_iupacs.c mutex.c threads.c	\
		  vcf_piz.c vcf_seg.c vcf_vblock.c vcf_header.c vcf_info.c vcf_samples.c vcf_liftover.c 		\
		  vcf_format_GT.c vcf_format_PS_PID.c vcf_info_QD.c vcf_info_SF.c                                                                           \
		  vcf_refalt.c vcf_tags.c vcf_linesort.c vcf_format.c vcf_illum_gtyping.c           		\
		  sam_seg.c sam_piz.c sam_shared.c sam_header.c sam_md.c sam_tlen.c sam_cigar.c sam_fields.c  			\
		  sam_sa.c bam_seg.c bam_seq.c bam_show.c                                             					\
		  sam_seq.c sam_qual.c sam_sag_zip.c sam_sag_piz.c sam_sag_load.c sam_sag_ingest.c sam_sag_scan.c    	\
		  sam_bwa.c sam_bowtie2.c sam_bsseeker2.c sam_bsbolt.c sam_bismark.c sam_gem3.c sam_tmap.c sam_hisat2.c \
		  sam_blasr.c sam_minimap2.c sam_cellranger.c sam_biobambam.c sam_pos.c									\
		  fasta.c fastq.c gff.c me23.c phylip.c chain.c kraken.c locs.c generic.c 								\
		  buffer.c random_access.c sections.c base64.c bgzf.c coverage.c txtheader.c lookback.c 				\
		  compressor.c codec.c codec_bz2.c codec_lzma.c codec_acgt.c codec_domq.c codec_hapmat.c codec_bsc.c	\
		  codec_pbwt.c codec_none.c codec_htscodecs.c codec_longr.c codec_normq.c    							\
	      txtfile.c profiler.c file.c dispatcher.c crypt.c aes.c md5.c segconf.c biopsy.c gencomp.c 			\
		  vblock.c regions.c  optimize.c dict_id.c hash.c stream.c url.c bases_filter.c dict_io.c 				\
		  recon_plan_io.c version.c

CONDA_COMPATIBILITY_SRCS =  compatibility/mac_gettime.c

ZLIB_SRCS  = zlib/gzlib.c zlib/gzread.c zlib/inflate.c zlib/inffast.c zlib/zutil.c zlib/inftrees.c zlib/deflate.c zlib/trees.c

BZLIB_SRCS = bzlib/blocksort.c bzlib/bzlib.c bzlib/compress.c bzlib/crctable.c bzlib/decompress.c bzlib/huffman.c bzlib/randtable.c

LZMA_SRCS  = lzma/LzmaEnc.c lzma/LzmaDec.c lzma/LzFind.c

BSC_SRCS   = bsc/divsufsort.c bsc/bwt.c bsc/coder.c bsc/libbsc.c bsc/lzp.c bsc/qlfc_model.c bsc/qlfc.c

HTSCODECS_SRC = htscodecs/rANS_static4x16pr.c htscodecs/rle.c htscodecs/pack.c htscodecs/arith_dynamic.c

DEFLATE_SRCS = libdeflate/deflate_compress.c libdeflate/deflate_decompress.c libdeflate/utils.c libdeflate/x86_cpu_features.c \
             libdeflate/arm_cpu_features.c libdeflate/crc32.c libdeflate/adler32.c

CONDA_DEVS = Makefile ../.gitignore 

CONDA_DOCS = ../LICENSE.txt ../AUTHORS ../README.md

CONDA_INCS = dict_id_gen.h aes.h dispatcher.h optimize.h profiler.h dict_id.h txtfile.h zip.h bits.h progress.h website.h \
             endianness.h md5.h sections.h text_help.h strings.h hash.h stream.h url.h flags.h segconf.h biopsy.h \
             buffer.h file.h context.h context_struct.h container.h seg.h text_license.h version.h compressor.h codec.h stats.h \
             crypt.h genozip.h piz.h vblock.h zfile.h random_access.h regions.h reconstruct.h tar.h qname.h qname_flavors.h \
			 lookback.h tokenizer.h codec_longr_alg.c gencomp.h dict_io.h recon_plan_io.h tip.h \
			 reference.h ref_private.h refhash.h ref_iupacs.h aligner.h mutex.h bgzf.h coverage.h threads.h \
			 arch.h license.h file_types.h data_types.h base64.h txtheader.h writer.h bases_filter.h genols.h contigs.h chrom.h \
			 vcf.h vcf_private.h sam.h sam_private.h me23.h fasta.h fasta_private.h fastq.h gff.h phylip.h chain.h kraken.h locs.h generic.h \
             compatibility/mac_gettime.h  \
			 zlib/gzguts.h zlib/inffast.h zlib/inffixed.h zlib/inflate.h zlib/inftrees.h zlib/zconf.h \
			 zlib/deflate.h zlib/trees.h \
			 zlib/zlib.h zlib/zutil.h \
			 lzma/7zTypes.h lzma/LzFind.h lzma/LzHash.h lzma/LzmaDec.h lzma/LzmaEnc.h \
			 bzlib/bzlib.h bzlib/bzlib_private.h \
			 htscodecs/rANS_static4x16.h htscodecs/rle.h htscodecs/pack.h htscodecs/arith_dynamic.h htscodecs/c_simple_model.h\
			 htscodecs/rANS_word.h htscodecs/htscodecs_endian.h htscodecs/rANS_word.h htscodecs/utils.h htscodecs/varint.h \
			 htscodecs/varint2.h htscodecs/utils.h \
			 bsc/bwt.h bsc/coder.h bsc/divsufsort.h bsc/libbsc.h bsc/lzp.h bsc/platform.h \
			 bsc/qlfc_model.h bsc/qlfc.h bsc/rangecoder.h bsc/tables.h \
 			 libdeflate/adler32_vec_template.h  libdeflate/crc32_table.h          libdeflate/unaligned.h \
 			 libdeflate/arm_adler32_impl.h      libdeflate/crc32_vec_template.h   libdeflate/x86_adler32_impl.h \
 			 libdeflate/arm_cpu_features.h      libdeflate/decompress_template.h  libdeflate/x86_cpu_features.h \
			 libdeflate/arm_crc32_impl.h        libdeflate/deflate_compress.h     libdeflate/x86_crc32_impl.h \
			 libdeflate/arm_matchfinder_impl.h  libdeflate/deflate_constants.h    libdeflate/x86_crc32_pclmul_template.h \
			 libdeflate/bt_matchfinder.h        libdeflate/hc_matchfinder.h       libdeflate/x86_decompress_impl.h \
			 libdeflate/common_defs.h           libdeflate/lib_common.h           libdeflate/x86_matchfinder_impl.h \
			 libdeflate/compiler_gcc.h          libdeflate/libdeflate.h \
			 libdeflate/cpu_features_common.h   libdeflate/matchfinder_common.h

LINUXDIR = genozip-linux-x86_64 # directory for creating the Linux binaries distribution

INSTALLERS = ../installers

ifeq ($(CC),cl) # Microsoft Visual C
	$(error Only the gcc compiler is currently supported)
endif

OBJDIR=objdir # fallback if not win, linux, mac
OBJDIR_LINUX=objdir.linux
OBJDIR_WINDOWS=objdir.windows
OBJDIR_MAC=objdir.mac

ifeq ($(OS),Windows_NT)
EXE = .exe
	LDFLAGS += -static -static-libgcc
	OBJDIR=$(OBJDIR_WINDOWS)
	WSL=wsl
else
    ifeq ($(uname),Linux)
        LDFLAGS += -lrt # required by pthreads
    	OBJDIR=$(OBJDIR_LINUX)
	endif

    ifeq ($(uname),Darwin) # Mac
		MY_SRCS += compatibility/mac_gettime.c
    	OBJDIR=$(OBJDIR_MAC)
    endif
endif

C_SRCS = $(MY_SRCS) $(ZLIB_SRCS) $(BZLIB_SRCS) $(BSC_SRCS) $(LZMA_SRCS) $(DEFLATE_SRCS) $(HTSCODECS_SRC)

OBJS       := $(addprefix $(OBJDIR)/, $(C_SRCS:.c=.o))
DEBUG_OBJS := $(addprefix $(OBJDIR)/, $(C_SRCS:.c=.debug-o)) 
OPT_OBJS   := $(addprefix $(OBJDIR)/, $(C_SRCS:.c=.opt-o))   # optimized but with debug info, for debugging issues that only manifest with compiler optimization
DEPS       := $(addprefix $(OBJDIR)/, $(C_SRCS:.c=.d)) 

EXECUTABLES       = genozip$(EXE)       genounzip$(EXE)       genocat$(EXE)       genols$(EXE)
DEBUG_EXECUTABLES = genozip-debug$(EXE) genounzip-debug$(EXE) genocat-debug$(EXE) genols-debug$(EXE)
OPT_EXECUTABLES   = genozip-opt$(EXE)   genounzip-opt$(EXE)   genocat-opt$(EXE)   genols-opt$(EXE)

all   : CFLAGS += $(OPTFLAGS) -DDISTRIBUTION=\"$(DISTRIBUTION)\"  
all   : $(OBJDIR) $(EXECUTABLES) 
	@chmod +x test.sh

debug : CFLAGS += $(DEBUGFLAGS) -DDISTRIBUTION=\"debug\"
debug : $(OBJDIR) $(DEBUG_EXECUTABLES)

opt   : CFLAGS += -g $(OPTFLAGS) -DDEBUG -DDISTRIBUTION=\"opt\"
opt   : $(OBJDIR) $(OPT_EXECUTABLES)

docker : CFLAGS += $(OPTFLAGS) -DDISTRIBUTION=\"Docker\"
docker : $(OBJDIR) $(EXECUTABLES) ../LICENSE.txt

-include $(DEPS)

$(OBJDIR): 
	@echo Making directory $@
	@mkdir $@ $(addprefix $@/, $(SRC_DIRS))

$(OBJDIR)/%.d: %.c | $(OBJDIR) # directory is an "order only prerequesite": https://www.gnu.org/savannah-checkouts/gnu/make/manual/html_node/Prerequisite-Types.html#Prerequisite-Types
	@echo Calculating dependencies $<
	@$(CC) $(CFLAGS) -MM -MT $@ $< -MF $(@:%.o=%.d)

$(OBJDIR)/%.o: %.c $(OBJDIR)/%.d
	@echo Compiling $<
	@$(CC) -c -o $@ $< $(CFLAGS)

$(OBJDIR)/%.debug-o: %.c $(OBJDIR)/%.d
	@echo "Compiling $< (debug)"
	@$(CC) -c -o $@ $< $(CFLAGS)

$(OBJDIR)/%.opt-o: %.c $(OBJDIR)/%.d
	@echo "Compiling $< (opt)"
	@$(CC) -c -o $@ $< $(CFLAGS)

%.S: %.c $(OBJDIR)/%.d
	@echo "Generating $@"
	@$(CC) -S -O3-o $@ $< $(CFLAGS)

%.E: %.c $(OBJDIR)/%.d
	@echo "Generating $@"
	@$(CC) -E -o $@ $< $(CFLAGS)

GENDICT_OBJS := $(addprefix $(OBJDIR)/, $(GENDICT_SRCS:.c=.o))

# dict_id_gen.h generation:
# Step 1: dict_id_gen.sh generates dict_id_gen.c, including all the GENDICT definitions from the data type include files (eg vcf.h)
# Step 2: dict_id_gen.sh compiles dict_id_gen.c and generate dict_id_gen[.exe]
# Step 3. dict_id_gen.sh generates dict_id_gen.h: it uses dict_id_gen[.exe]to generate the field constant, and then adds the fields enum and mapping
ifeq ($(OS),Windows_NT)

dict_id_gen.h : $(shell grep -w "pragma GENDICT" *.h | cut -d: -f1 | uniq) dict_id_gen.sh
	@echo Generating $@
	@./dict_id_gen.sh $(CC)

endif # ugly hack to avoid conda failure due to bash issues in dict_id_gen.sh - pre-generate on Windows and check in to github

genozip$(EXE): dict_id_gen.h $(OBJS)
	@echo Linking $@
	@$(CC) -o $@ $(OBJS) $(CFLAGS) $(LDFLAGS)
	@ln -f $@ ../$@
 
genozip-debug$(EXE): dict_id_gen.h $(DEBUG_OBJS)
	@echo Linking $@
	@$(CC) -o $@ $(DEBUG_OBJS) $(CFLAGS) $(LDFLAGS) 
	@ln -f $@ ../$@

genozip-opt$(EXE): dict_id_gen.h $(OPT_OBJS)
	@echo Linking $@
	@$(CC) -o $@ $(OPT_OBJS) $(CFLAGS) $(LDFLAGS)
	@ln -f $@ ../$@

genounzip$(EXE) genocat$(EXE) genols$(EXE): genozip$(EXE)
	@echo Hard linking $@
	@rm -f $@ 
	@ln -f $^ $@
	@ln -f $@ ../$@

genounzip-debug$(EXE) genocat-debug$(EXE) genols-debug$(EXE): genozip-debug$(EXE)
	@echo Hard linking $@
	@rm -f $@ 
	@ln -f $^ $@
	@ln -f $@ ../$@

genounzip-opt$(EXE) genocat-opt$(EXE) genols-opt$(EXE): genozip-opt$(EXE)
	@echo Hard linking $@
	@rm -f $@ 
	@ln -f $^ $@
	@ln -f $@ ../$@

../LICENSE.txt: text_license.h version.h # not dependent on genozip.exe, so we don't generate it every compilation
	@make -j genozip$(EXE) # recursive call to make genozip.exe with the latest version
	@echo Generating $@
	@./genozip$(EXE) --license=100 --force > $@
	@rm -f $(OBJDIR)/arch.o
	
# this is used by build.sh to install on conda for Linux and Mac. Installation for Windows in in bld.bat
install: genozip$(EXE)
	@echo Installing in $(PREFIX)/bin
	@if ( test ! -d $(PREFIX)/bin ) ; then mkdir -p $(PREFIX)/bin ; fi
	@cp -f genozip$(EXE) $(PREFIX)/bin/genozip$(EXE)
ifneq ($(OS),Windows_NT)
	@chmod a+x $(PREFIX)/bin/genozip$(EXE)
endif
	@cp -f $(PREFIX)/bin/genozip$(EXE) $(PREFIX)/bin/genounzip$(EXE)
	@cp -f $(PREFIX)/bin/genozip$(EXE) $(PREFIX)/bin/genocat$(EXE)
	@cp -f $(PREFIX)/bin/genozip$(EXE) $(PREFIX)/bin/genols$(EXE)

version = $(shell head -n1 version.h |cut -d\" -f2)

SH_VERIFY_ALL_COMMITTED = (( `git status |grep 'modified\|Untracked files'|grep -v .gitkeep |wc -l` == 0 )) || \
                          (echo ERROR: there are some uncommitted changes: ; echo ; git status ; exit 1)

test:
	@cat test.sh | tr -d "\r" | bash -

clean-installers:
	@rm -fR $(INSTALLERS)/*

clean-debug:
	@echo Cleaning up debug
	@rm -f $(DEBUG_OBJS) $(DEBUG_EXECUTABLES) $(OBJDIR)/*.debug-o
	@rm -f $(OPT_OBJS) $(OPT_EXECUTABLES) $(OBJDIR)/*.opt-o

clean-optimized:
	@echo Cleaning up optimized
	@rm -f $(OBJS) $(EXECUTABLES) $(OBJDIR)/*.o

clean-distribution: # clean all platforms
	@echo Cleaning up all platforms for distribution
	@rm -f $(OBJDIR_LINUX)/*.o $(OBJDIR_LINUX)/*.d $(OBJDIR_WINDOWS)/*.o $(OBJDIR_WINDOWS)/*.d $(OBJDIR_MAC)/*.o $(OBJDIR_MAC)/*.d 

clean-opt:
	@echo Cleaning up opt
	@rm -f $(OPT_OBJS) $(EXECUTABLES) $(OPT_EXECUTABLES)/*.opt-o

clean-test:
	@rm -f $(TEST_DIR)/*.good $(TEST_DIR)/*.bad $(TEST_DIR)/*.local $(TEST_DIR)/*.b250 $(TEST_DIR)/tmp/* $(TEST_DIR)/*.DEPN $(TEST_DIR)/back-compat/*/*.bad

clean: clean-installers clean-test
	@echo Cleaning up
	@rm -f $(DEPS) $(filter-out LICENSE.txt,$(WINDOWS_INSTALLER_OBJS)) *.d .archive.tar.gz *.stackdump $(EXECUTABLES) $(OPT_EXECUTABLES) $(DEBUG_EXECUTABLES) 
	@rm -f *.S *.good *.bad $(DATA)/*.good $(DATA)/*.bad *.local genozip.threads-log.* *.b250 
	@rm -R $(OBJDIR)
	@mkdir $(OBJDIR) $(addprefix $(OBJDIR)/, $(SRC_DIRS))

.PHONY: clean clean-debug clean-optimized clean-installers clean-distribution git-pull macos mac/.remote_mac_timestamp delete-arch testfiles test-backup genozip-linux-x86_64/clean genozip-prod genozip-prod.exe dict_id_gen$(EXE) push-build increment-version $(INSTALLERS)/LICENSE.html rm_arch_o

# builds prod for local OS
genozip-prod$(EXE): 
	@echo "building prod"
	@(cd $(PROD_SRC) ; git stash ; git pull ; rm -Rf $(OBJDIR) ; make -j clean ; touch dict_id_gen.h ; make -j)
	@cp $(PROD_SRC)/genozip$(EXE) $(DEV_DIR)/genozip-prod$(EXE)
	@cp $(PROD_SRC)/genozip$(EXE) $(RELEASES)/genozip-$(version)$(EXE)
	@cp $(PROD_SRC)/genounzip$(EXE) $(DEV_DIR)/genounzip-prod$(EXE)
	@cp $(PROD_SRC)/genocat$(EXE) $(DEV_DIR)/genocat-prod$(EXE)

# currently, I build for conda from my Windows machine so I don't bother supporting other platforms
ifeq ($(OS),Windows_NT)

# When running on Windows, builds prod for Linux
genozip-prod:
	@run-on-wsl.sh make -j genozip-prod

# increments minor version, eg. 1.0.1 -> 1.0.2. 
# To increment a major version, manually edit version.h and set minor version to -1 e.g. 1.1.-1 (careful! no newlines or spaces)
# and re-compile so that genozip --version gets updated
# IMPORTANT: the first number in the version indicates the genozip file format version and goes into
# the genozip file header SectionHeaderTxtHeader.genozip_version
#increment-version: $(C_SRCS) $(CONDA_COMPATIBILITY_SRCS) $(CONDA_DEVS) $(CONDA_DOCS) $(CONDA_INCS) # note: target name is not "version.h" so this is not invoked during "make all" or "make debug"
increment-version: # note: target name is not "version.h" so this is not invoked during "make all" or "make debug"
	@echo "Incrementing version.h"
	@bash $(SCRIPTS)/increment-version.sh

decrement-version:
	@echo "Do manually:"
	@echo "Remove tag: git push --delete origin genozip-a.b.c"
	@echo "Change version.h to the last version that still has a tag"

.archive.tar.gz : $(C_SRCS) $(CONDA_COMPATIBILITY_SRCS) $(CONDA_DEVS) $(CONDA_DOCS) $(CONDA_INCS) ../LICENSE.txt
	@echo Creating github tag genozip-$(version) and archive
	@$(SH_VERIFY_ALL_COMMITTED)
	@git push > /dev/null
	@(git tag genozip-$(version) > /dev/null ; exit 0) # ignore errors if already exists (due to re-run of make)
	@git push origin genozip-$(version) > /dev/null
	@curl https://github.com/divonlan/genozip/archive/genozip-$(version).tar.gz --ssl-no-revoke --silent --location -o $@ > /dev/null
	@echo GITHUB: go to here: https://github.com/divonlan/genozip/releases/new
	@echo "1. Set 'Tag version' and 'Release title' are both: genozip-$(version)"
	@echo "2. Copy the notes for the version from RELEASE NOTES"
	@echo "3. Update 'latest' release to new tag"

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
	@cp conda/README.md $(CONDA_FEEDSTOCK)
	@cp conda/meta.yaml conda/build.sh conda/bld.bat $(CONDA_RECIPE_DIR)
	@echo "Committing my files to branch genozip on my fork"
	@(cd $(CONDA_FEEDSTOCK); git pull; git commit -m "update" recipe/meta.yaml README.md recipe/build.sh recipe/bld.bat; git push) > /dev/null
	@echo " "
	@echo "Submitting pull request to conda-forge"
#	@(cd $(CONDA_RECIPE_DIR); git request-pull master https://github.com://conda-forge/genozip-feedstock master)
#	@(cd $(CONDA_RECIPE_DIR); git request-pull master https://github.com/divonlan/genozip-feedstock master)
	@touch $@
	@echo "CONDA: Using a browser:"
	@echo "  (1) Go to https://github.com/divonlan/genozip-feedstock/tree/genozip "
	@echo "  (2) Click 'Contribute' -> 'Open pull request'"
	@echo "  (3) Go to https://dev.azure.com/conda-forge/feedstock-builds/_build "
	@echo "  (4) Click on genozip and wait (~8 min) for the test to complete. Fix any issues."
	@echo "  (5) Go back to the tab in (2) and click 'Merge pull request' and the 'Confirm merge' (DONT CLICK 'Delete branch')"
	@echo "  (6) Go to https://dev.azure.com/conda-forge/feedstock-builds/_build and watch the build - it should be fine"
	@echo "  (7) In ~30 minutes users will be able to 'conda update genozip'"

# Building Windows InstallForge with distribution flag: we delete arch.o to force it to re-compile with DISTRIBUTION=InstallForge.
windows/%.exe: CFLAGS += $(OPTFLAGS) -DDISTRIBUTION=\"InstallForge\" # note: matches the name in version.c
windows/%.exe: $(OBJS) %.exe
	@echo Linking $@
	@(mkdir windows >& /dev/null ; exit 0)
	@$(CC) -o $@ $(OBJS) $(CFLAGS) $(LDFLAGS)

windows/LICENSE.for-installer.txt: genozip$(EXE) version.h
	@echo Generating $@
	@(mkdir windows >& /dev/null ; exit 0)
	@./genozip$(EXE) --license=60 --force > $@

rm_arch_o : 
	@(rm -f $(OBJDIR)/arch.o >& /dev/null ; exit 0) # remove this arch.o which contains DISTRIBUTION

WINDOWS_INSTALLER_OBJS = windows/genozip.exe windows/genounzip.exe windows/genocat.exe windows/genols.exe windows/LICENSE.for-installer.txt ../LICENSE.txt

$(INSTALLERS)/LICENSE.html: genozip$(EXE)
	@./genozip$(EXE) --license=-5000 --force > $@

# this must be run ONLY has part of "make distribution" or else versions will be out of sync
$(INSTALLERS)/genozip-installer.exe: rm_arch_o $(WINDOWS_INSTALLER_OBJS) 
	@(mkdir windows >& /dev/null ; exit 0)
	@echo 'Creating Windows installer'
	@echo 'WINDOWS: Using the UI:'
	@echo '  (1) Open genozip-installer.ifp'
	@echo '  (2) Set General-Program version to $(version)'
	@echo '  (3) In Dialogs->License upload license from windows/LICENSE-for-installer.txt'
	@echo '  (3) Click Save, then click Build'
	@echo '  (4) Optionally: Click Yes, and copy the resulting files to releases/* and also c:\bin'	
	@echo '  (5) Exit the UI (close the window)'
	@if [ `cd ..; basename ${PWD}` != genozip ] ; then cp $(WINDOWS_INSTALLER_OBJS) $(DEV_SRC)/windows ; cp $(DEV_SRC)/genozip-installer.ifp $(DEV_SRC)/genozip-installer.ifp.save ; fi # so this works for genozip-prod too - because InstallForge uses absolute paths 
	@($(UTILS)/InstallForge/InstallForge.exe ; exit 0)
	@rm -f $(OBJDIR)/arch.o # remove this arch.o which contains DISTRIBUTION
	@if [ `cd ..; basename ${PWD}` != genozip ] ; then mv $(DEV_SRC)/genozip-installer.ifp . ; mv $(DEV_SRC)/genozip-installer.ifp.save $(DEV_SRC)/genozip-installer.ifp ; fi # we always edit the version in the genozip dir
#	@(C:\\\\Program\\ Files\\ \\(x86\\)\\\\solicus\\\\InstallForge\\\\bin\\\\ifbuilderenvx86.exe ; exit 0)

$(INSTALLERS)/genozip-linux-x86_64.tar.build: genozip-linux-x86_64/LICENSE.txt 
	@(mkdir genozip-linux-x86_64 >& /dev/null ; exit 0)
	@run-on-wsl.sh make $(INSTALLERS)/genozip-linux-x86_64.tar # make -j doesn't work well on WSL - filesystem clock issues

mac/.remote_mac_timestamp: # to be run from Windows to build on a remote mac
	@echo "Creating Mac installer"
	@echo "Logging in to remote mac" 
	@# Get IP address - check if the previous IP address still works or ask for a new one. Assuming a LAN on an Android hotspot.
	@ip=`cat mac/.mac_ip_address` ; a=`echo $$ip|cut -d. -f4`; (( `ping  -n 1 $$ip | grep "round trip times" | wc -l` > 0 )) || read -p "IP Address: 192.168.43." a ; ip=192.168.43.$$a ; echo $$ip > mac/.mac_ip_address
	@[ -f mac/.mac_username ] || ( echo ERROR: file mac/.mac_username missing && exit 1 )
	@ssh `cat mac/.mac_ip_address` -l `cat mac/.mac_username`  "cd genozip ; echo "Pulling from git" ; git pull >& /dev/null ; make -j mac/.from_remote_timestamp" # pull before make as Makefile might have to be pulled
	@touch $@

BUILD_FILES = version.h genozip-installer.ifp ../LICENSE.txt Makefile

BUILD_FILES_DOCS = genozip-installer.exe genozip-linux-x86_64.tar LICENSE.html

push-build: 
	@(git stage $(BUILD_FILES) ; exit 0) > /dev/null
	@(git commit -m $(version) ; exit 0) > /dev/null
	@git push > /dev/null
	@(cd $(INSTALLERS); git stage $(BUILD_FILES_DOCS) ; exit 0) > /dev/null
	@(git commit -m $(version) ; exit 0) > /dev/null
	@git push > /dev/null

distribution: clean-distribution increment-version testfiles ../LICENSE.txt $(INSTALLERS)/LICENSE.html $(INSTALLERS)/genozip-linux-x86_64.tar.build $(INSTALLERS)/genozip-installer.exe push-build conda/.conda-timestamp genozip-prod.exe genozip-prod
	@(cd $(CONDA_FEEDSTOCK) ; git pull)

distribution-maintenance: clean-distribution increment-version ../LICENSE.txt testfiles $(INSTALLERS)/genozip-linux-x86_64.tar.build $(INSTALLERS)/genozip-installer.exe \
                          push-build conda/.conda-timestamp genozip-prod.exe genozip-prod
	@(cd $(CONDA_FEEDSTOCK) ; git pull)

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

# Building Linux pre-compiled binaries: we delete arch.o to force it to re-compile with DISTRIBUTION=linux-x86_64.
ifeq ($(uname),Linux)

genozip-linux-x86_64/clean:
	@rm -fR genozip-linux-x86_64
	@mkdir genozip-linux-x86_64

# note: getpwuid and getgrgid will cause dymanically loading of the locally installed glibc in the --tar option, or segfault. that's normally fine.
genozip-linux-x86_64/genozip: CFLAGS += $(OPTFLAGS) -DDISTRIBUTION=\"linux-x86_64\" # note: matches the name in version.c
genozip-linux-x86_64/genozip: clean-optimized $(OBJS)  # clean first, as we will compile without march=native
	@echo Linking $@
	@$(CC) -static -o $@ $(OBJS) $(CFLAGS) $(LDFLAGS)

genozip-linux-x86_64/genounzip genozip-linux-x86_64/genocat genozip-linux-x86_64/genols: genozip-linux-x86_64/genozip 
	@echo Generating $@
	@ln -f $< $@

LINUX_TARGZ_OBJS = genozip-linux-x86_64/genozip genozip-linux-x86_64/genounzip genozip-linux-x86_64/genocat genozip-linux-x86_64/genols 

# this must be run ONLY as part of "make distribution" or else versions will be out of sync
$(INSTALLERS)/genozip-linux-x86_64.tar: version.h genozip-linux-x86_64/clean $(LINUX_TARGZ_OBJS) # run on Linux by make-linux-from-windows.sh
	@echo "Creating $@"
	@tar cf $@ -z genozip-linux-x86_64
	@rm -f $(OBJDIR)/arch.o # remove this arch.o which contains DISTRIBUTION

endif # Linux

ifeq ($(uname),Darwin)

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

app_specific_pw := $(shell cat .altool_app_specifc_password)

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