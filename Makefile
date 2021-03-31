# ------------------------------------------------------------------
#   Makefile
#   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
#   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

# Note for Windows: to run this make, you need mingw (for the gcc compiler) and cygwin (for Unix-like tools):
# Mingw: http://mingw-w64.org/doku.php 
# Cygwin: https://www.cygwin.com/

ifdef BUILD_PREFIX
IS_CONDA=1
endif

LDFLAGS     += -lpthread -lm 

ifdef IS_CONDA 
	CFLAGS  += -Wall -I. -D_LARGEFILE64_SOURCE=1 -DDISTRIBUTION=\"conda\"
	LDFLAGS += -lbz2 # conda - dynamic linking with bz2

	ifeq ($(OS),Windows_NT)
		CC=gcc # in Windows, override conda's default Visual C with gcc 
		LDFLAGS += -L$(PREFIX)/Library/lib
		CFLAGS  += -I$(PREFIX)/Library/include 
	endif
else
	CC=gcc
	CFLAGS = -Wall -I. -Izlib -Ibzlib -Ilibdeflate -D_LARGEFILE64_SOURCE=1 -march=native 
endif 

SRC_DIRS = zlib bzlib lzma bsc libdeflate compatibility

MY_SRCS = genozip.c base250.c context.c container.c strings.c stats.c arch.c license.c data_types.c bit_array.c progress.c \
          zip.c piz.c reconstruct.c seg.c zfile.c aligner.c flags.c digest.c mutex.c liftover.c sorter.c threads.c \
		  txtheader.c reference.c ref_lock.c refhash.c ref_make.c ref_contigs.c ref_alt_chroms.c  \
		  vcf_piz.c vcf_seg.c vcf_shared.c vcf_header.c \
          sam_seg.c sam_piz.c sam_seg_bam.c sam_shared.c sam_header.c \
		  fasta.c fastq.c gff3_seg.c me23.c phylip.c chain.c generic.c \
		  buffer.c random_access.c sections.c base64.c bgzf.c coverage.c \
		  compressor.c codec.c codec_bz2.c codec_lzma.c codec_acgt.c codec_domq.c codec_hapmat.c codec_bsc.c\
		  codec_gtshark.c codec_pbwt.c codec_none.c \
	      txtfile.c profiler.c file.c dispatcher.c crypt.c aes.c md5.c \
		  vblock.c regions.c  optimize.c dict_id.c hash.c stream.c url.c

CONDA_COMPATIBILITY_SRCS =  compatibility/mac_gettime.c

ZLIB_SRCS  = zlib/gzlib.c zlib/gzread.c zlib/inflate.c zlib/inffast.c zlib/zutil.c zlib/inftrees.c zlib/deflate.c zlib/trees.c

BZLIB_SRCS = bzlib/blocksort.c bzlib/bzlib.c bzlib/compress.c bzlib/crctable.c bzlib/decompress.c bzlib/huffman.c bzlib/randtable.c

LZMA_SRCS  = lzma/LzmaEnc.c lzma/LzmaDec.c lzma/LzFind.c

BSC_SRCS   = bsc/divsufsort.c bsc/bwt.c bsc/coder.c bsc/libbsc.c bsc/lzp.c bsc/qlfc_model.c bsc/qlfc.c

DEFLATE_SRCS = libdeflate/deflate_compress.c libdeflate/deflate_decompress.c libdeflate/utils.c libdeflate/x86_cpu_features.c \
             libdeflate/arm_cpu_features.c libdeflate/crc32.c libdeflate/adler32.c

CONDA_DEVS = Makefile .gitignore 

CONDA_DOCS = LICENSE.non-commercial.txt LICENSE.commercial.txt AUTHORS README.md

CONDA_INCS = aes.h dispatcher.h optimize.h profiler.h dict_id.h txtfile.h zip.h bit_array.h progress.h \
             base250.h endianness.h md5.h sections.h text_help.h strings.h hash.h stream.h url.h flags.h \
             buffer.h file.h context.h container.h seg.h text_license.h version.h compressor.h codec.h stats.h \
             crypt.h genozip.h piz.h vblock.h zfile.h random_access.h regions.h reconstruct.h liftover.h  \
			 reference.h ref_private.h refhash.h aligner.h mutex.h bgzf.h coverage.h sorter.h threads.h \
			 arch.h license.h data_types.h base64.h txtheader.h \
			 vcf.h vcf_private.h sam.h sam_private.h me23.h fasta.h fastq.h gff3.h phylip.h chain.h generic.h \
             compatibility/mac_gettime.h  \
			 zlib/gzguts.h zlib/inffast.h zlib/inffixed.h zlib/inflate.h zlib/inftrees.h zlib/zconf.h \
			 zlib/deflate.h zlib/trees.h	 \
			 zlib/zlib.h zlib/zutil.h \
			 lzma/7zTypes.h lzma/Compiler.h lzma/LzFind.h lzma/LzFindMt.h lzma/LzHash.h lzma/LzmaDec.h lzma/LzmaEnc.h \
			 lzma/Precomp.h lzma/Threads.h \
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

BAM_FILES = test/basic.bam test/minimal.bam

ifeq ($(CC),cl) # Microsoft Visual C
	$(error Only the gcc compiler is currently supported)
endif

OBJDIR=objdir # fallback if not win, linux, mac

ifeq ($(OS),Windows_NT)
# Windows
	EXE = .exe
	LDFLAGS += -static -static-libgcc
	LZMA_SRCS += lzma/Threads.c lzma/LzFindMt.c
	OBJDIR=objdir.windows
else
	CFLAGS += -D_7ZIP_ST
    uname := $(shell uname -s)
    ifeq ($(uname),Linux)
# Linux
        LDFLAGS += -lrt
    	OBJDIR=objdir.linux
	endif
    ifeq ($(uname),Darwin)
# Mac
		MY_SRCS += compatibility/mac_gettime.c
    	OBJDIR=objdir.mac
    endif
endif

ifndef IS_CONDA 
	# local - static link everything
	C_SRCS = $(MY_SRCS) $(ZLIB_SRCS) $(BZLIB_SRCS) $(BSC_SRCS) $(LZMA_SRCS) $(DEFLATE_SRCS)

#	ifneq ($(shell uname -a | grep ppc64),)
#		CFLAGS += -mcpu=native 
#	endif

else  # conda
	# use packages for bzip2
	C_SRCS = $(MY_SRCS) $(ZLIB_SRCS) $(LZMA_SRCS) $(BSC_SRCS) $(DEFLATE_SRCS)
endif

OBJS       := $(addprefix $(OBJDIR)/, $(C_SRCS:.c=.o))
DEBUG_OBJS := $(addprefix $(OBJDIR)/, $(C_SRCS:.c=.debug-o)) 
OPT_OBJS   := $(addprefix $(OBJDIR)/, $(C_SRCS:.c=.opt-o))   # optimized but with debug info, for debugging issues that only manifest with compiler optimization
DEPS       := $(addprefix $(OBJDIR)/, $(C_SRCS:.c=.d)) 

EXECUTABLES       = genozip$(EXE)       genounzip$(EXE)       genocat$(EXE)       genols$(EXE)
DEBUG_EXECUTABLES = genozip-debug$(EXE) genounzip-debug$(EXE) genocat-debug$(EXE) genols-debug$(EXE)
OPT_EXECUTABLES   = genozip-opt$(EXE)   genounzip-opt$(EXE)   genocat-opt$(EXE)   genols-opt$(EXE)

ifeq ($(CC),gcc)
	OPTFLAGS += -Ofast -std=gnu99
	DEBUGFLAGS += -std=gnu99 -DDEBUG -g -O0
else
	OPTFLAGS += -O2
	DEBUGFLAGS += -DDEBUG -g -O0
endif

all   : CFLAGS += $(OPTFLAGS) 
all   : $(OBJDIR) $(EXECUTABLES) LICENSE.non-commercial.txt
	@chmod +x test.sh

debug : CFLAGS += $(DEBUGFLAGS)
debug : $(OBJDIR) $(DEBUG_EXECUTABLES) LICENSE.non-commercial.txt

opt   : CFLAGS += -g $(OPTFLAGS)
opt   : $(OBJDIR) $(OPT_EXECUTABLES) LICENSE.non-commercial.txt

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

test/%.bam : test/%.sam
	@echo "Generating $@ from $<"
	@wsl bash -c "exec /home/divon/miniconda3/bin/samtools  view $< -OBAM -o $@"

genozip$(EXE): $(OBJS)
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)
 
genozip-debug$(EXE): $(DEBUG_OBJS)
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) 

genozip-opt$(EXE): $(OPT_OBJS)
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

genounzip$(EXE) genocat$(EXE) genols$(EXE): genozip$(EXE)
	@echo Hard linking $@
	@rm -f $@ 
	@ln $^ $@

genounzip-debug$(EXE) genocat-debug$(EXE) genols-debug$(EXE): genozip-debug$(EXE)
	@echo Hard linking $@
	@rm -f $@ 
	@ln $^ $@

genounzip-opt$(EXE) genocat-opt$(EXE) genols-opt$(EXE): genozip-opt$(EXE)
	@echo Hard linking $@
	@rm -f $@ 
	@ln $^ $@

LICENSE.non-commercial.txt: text_license.h # not dependent on genozip.exe, so we don't generate it every compilation
	@make ./genozip$(EXE) # recursive call 
	@echo Generating $@
	@./genozip$(EXE) --license=100 > $@

SPHINX = /home/divon/miniconda3/bin/sphinx-build
DOCS = docs/genozip.rst docs/genounzip.rst docs/genocat.rst docs/genols.rst docs/developer.rst docs/index.rst docs/license.rst \
       docs/publications.rst docs/installing.rst docs/contact.rst docs/examples.rst docs/source.rst docs/logo.png \
	   docs/opt-help.rst docs/opt-piz.rst docs/opt-quiet.rst docs/opt-stats.rst docs/opt-threads.rst docs/opt-translation.rst \
	   docs/manual.rst docs/sex-assignment.rst docs/sex-assignment-alg-sam.rst docs/sex-assignment-alg-fastq.rst \
	   docs/fastq-to-bam-pipeline.rst docs/coverage.rst docs/algorithms.rst docs/losslessness.rst docs/idxstats.rst \
	   docs/downsampling.rst docs/dual-coordinates-vcf.rst docs/pipelines.rst docs/capabilities.rst

docs/conf.py: docs/conf.template.py version.h
	@sed -e "s/__VERSION__/$(version)/g" $< |sed -e "s/__YEAR__/`date +'%Y'`/g" > $@ 

docs/LICENSE.for-docs.txt: genozip$(EXE)
	@./genozip$(EXE) --license=74 > $@

docs/_build/html/.buildinfo: docs/LICENSE.for-docs.txt docs/conf.py $(DOCS)
	@echo Building HTML docs
	@wsl $(SPHINX) -M html docs docs/_build -q -a 

docs: docs/_build/html/.buildinfo 

docs-debug: docs/_build/html/.buildinfo
	@(C:\\\\Program\\ Files\\ \\(x86\\)\\\\Google\\\\Chrome\\\\Application\\\\chrome.exe file:///C:/Users/USER/projects/genozip/docs/_build/html/index.html; exit 0)
	
testfiles : $(BAM_FILES)

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

SH_VERIFY_ALL_COMMITTED = (( `git status|grep 'modified\|Untracked files'|grep -v .gitkeep |wc -l ` == 0 )) || \
                          (echo ERROR: there are some uncommitted changes: ; echo ; git status ; exit 1)

# currently, I build for conda from my Windows machine so I don't bother supporting other platforms
ifeq ($(OS),Windows_NT)

# increments minor version, eg. 1.0.1 -> 1.0.2. 
# To increment a major version, manually edit version.h and set minor version to -1 e.g. 1.1.-1 (careful! no newlines or spaces)
# and re-compile so that genozip --version gets updated
# IMPORTANT: the first number in the version indicates the genozip file format version and goes into
# the genozip file header SectionHeaderTxtHeader.genozip_version
increment-version: $(C_SRCS) $(CONDA_COMPATIBILITY_SRCS) $(CONDA_DEVS) $(CONDA_DOCS) $(CONDA_INCS) # note: target name is not "version.h" so this is not invoked during "make all" or "make debug"
	@echo "Incrementing version.h"
	@$(SH_VERIFY_ALL_COMMITTED)
	@bash increment-version.sh
	@git commit -m "increment version" version.h 

decrement-version:
	@echo "Do manually:"
	@echo "Remove tag: git push --delete origin genozip-a.b.c"
	@echo "Change version.h to the last version that still has a tag"

.archive.tar.gz : increment-version $(C_SRCS) $(CONDA_COMPATIBILITY_SRCS) $(CONDA_DEVS) $(CONDA_DOCS) $(CONDA_INCS) 
	@echo Creating github tag genozip-$(version) and archive
	@$(SH_VERIFY_ALL_COMMITTED)
	@git push > /dev/null
	@git tag genozip-$(version) > /dev/null
	@git push origin genozip-$(version) > /dev/null
	@curl https://github.com/divonlan/genozip/archive/genozip-$(version).tar.gz --silent --location -o $@ > /dev/null
	@echo GITHUB: go to here: https://github.com/divonlan/genozip/releases/new
	@echo "1. Set 'Tag version' and 'Release title' are both: genozip-$(version)"
	@echo "2. Copy the notes for the version from RELEASE NOTES"

conda/meta.yaml: conda/meta.template.yaml .archive.tar.gz README.md
	@echo "Generating conda/meta.yaml"
	@bash conda/generate_meta.sh > $@

conda/README.md: conda/README.template.md html-to-md.sed README.md
	@echo "Generating conda/README.md"
	@bash conda/generate_README.sh > $@

#CONDA_RECIPE_DIR = ../staged-recipes/recipes/genozip # initial stage-recipes step, keeping here for future reference
CONDA_FEEDSTOCK  = ../genozip-feedstock
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
	@echo "  (1) Go to https://github.com/conda-forge/genozip-feedstock/pulls"
	@echo "  (2) Click 'Compare and pull request' then 'Create pull request' and wait 30 seconds for the test to start"
	@echo "      Fallback: if you can't see 'Compare & pull', manually created a pull request 'into conda-forge:master from divonlan:genozip'"
	@echo "  (3) Go to https://dev.azure.com/conda-forge/feedstock-builds/_build"
	@echo "  (4) Click on genozip and wait (~5 min) for the test to complete. Fix any issues."
	@echo "  (5) Go back to the tab in (2) and click 'Merge pull request' and the 'Confirm merge' (DONT CLICK 'Delete branch')"
	@echo "  (6) Go to https://dev.azure.com/conda-forge/feedstock-builds/_build and watch the build - it should be fine"
	@echo "  (7) In ~30 minutes users will be able to 'conda update genozip'"

# Building Windows InstallForge with distribution flag: we delete arch.o to force it to re-compile with DISTRIBUTION=InstallForge.
delete-arch: 
	@rm -f arch.o

windows/%.exe: CFLAGS += -DDISTRIBUTION=\"InstallForge\"
windows/%.exe: delete-arch $(OBJS) %.exe
	@echo Linking $@
	@$(CC) -o $@ $(OBJS) $(CFLAGS) $(LDFLAGS)

windows/readme.txt: $(EXECUTABLES)
	@echo Generating $@
	@./genozip$(EXE)   --help  > $@
	@printf '%.s-' {1..120}   >> $@
	@./genounzip$(EXE) --help >> $@
	@printf '%.s-' {1..120}   >> $@
	@./genols$(EXE)    --help >> $@
	@printf '%.s-' {1..120}   >> $@
	@./genocat$(EXE)   --help >> $@

windows/LICENSE.for-installer.txt: genozip$(EXE)
	@echo Generating $@
	@./genozip$(EXE) --license=60 > $@

WINDOWS_INSTALLER_OBJS = windows/genozip.exe windows/genounzip.exe windows/genocat.exe windows/genols.exe \
                         windows/LICENSE.for-installer.txt windows/readme.txt

# this must be run ONLY has part of "make distribution" or else versions will be out of sync
docs/genozip-installer.exe: $(WINDOWS_INSTALLER_OBJS) LICENSE.commercial.txt LICENSE.non-commercial.txt 
	@echo 'Creating Windows installer'
	@$(SH_VERIFY_ALL_COMMITTED)
	@echo 'WINDOWS: Using the UI:'
	@echo '  (1) Open genozip-installer.ifp'
	@echo '  (2) Set General-Program version to $(version)'
	@echo '  (3) If changed, update license in Dialogs->License from windows/LICENSE-for-installer.txt'
	@echo '  (3) Click Save, then click Build'
	@echo '  (4) Optionally: Click Yes, and copy the resulting files to releases/* and also c:\bin'	
	@echo '  (5) Exit the UI (close the window)'
	@cp $(WINDOWS_INSTALLER_OBJS) ../genozip/windows # so this works for genozip-prod too - because InstallForge uses absolute paths
	@(C:\\\\Program\\ Files\\ \\(x86\\)\\\\solicus\\\\InstallForge\\\\InstallForge.exe ; exit 0)
	@echo 'Committing Windows installer and pushing to repo'
	@mv ../genozip/windows/genozip-installer.exe docs  # so this works for genozip-prod too - because InstallForge uses absolute paths
	@(git stage genozip-installer.ifp $@ ; exit 0) > /dev/null
	@(git commit -m windows_files_for_version_$(version) genozip-installer.ifp $@ ; exit 0) > /dev/null
	@git push > /dev/null
	@rm -f arch.o # remove this arch.o which contains DISTRIBUTION

mac/.remote_mac_timestamp: # to be run from Windows to build on a remote mac
	@echo "Creating Mac installer"
	@$(SH_VERIFY_ALL_COMMITTED)
	@echo "Pushing all committed changes to github"
	@(( `git push |& grep "Everything up-to-date"  | wc -l` > 0 )) || (echo "Pushed some stuff... waiting 5 seconds for it to settle in" ; sleep 5)
	@echo "Logging in to remote mac" 
	@# Get IP address - check if the previous IP address still works or ask for a new one. Assuming a LAN on an Android hotspot.
	@ip=`cat mac/.mac_ip_address` ; a=`echo $$ip|cut -d. -f4`; (( `ping  -n 1 $$ip | grep "round trip times" | wc -l` > 0 )) || read -p "IP Address: 192.168.43." a ; ip=192.168.43.$$a ; echo $$ip > mac/.mac_ip_address
	@[ -f mac/.mac_username ] || ( echo ERROR: file mac/.mac_username missing && exit 1 )
	@ssh `cat mac/.mac_ip_address` -l `cat mac/.mac_username`  "cd genozip ; echo "Pulling from git" ; git pull >& /dev/null ; make mac/.from_remote_timestamp" # pull before make as Makefile might have to be pulled
	@touch $@

distribution: CFLAGS := $(filter-out -march=native,$(CFLAGS))
distribution: docs testfiles conda/.conda-timestamp docs/genozip-installer.exe # mac/.remote_mac_timestamp
	
endif # Windows

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

$(MACRSSDIR)/README.html: README.md
	@cp -f $< $@

$(MACRSSDIR)/LICENSE.non-commercial.txt: LICENSE.non-commercial.txt
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
                                    $(MACRSSDIR)/welcome.html $(MACRSSDIR)/README.html $(MACRSSDIR)/LICENSE.non-commercial.txt
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

test:
	@cat test.sh | tr -d "\r" | bash -

clean-test.sh-files: 
	@rm -f test/tmp/* test/*.rejects.*

clean-debug:
	@echo Cleaning up debug
	@rm -f $(DEBUG_OBJS) $(DEBUG_EXECUTABLES) $(OBJDIR)/*.debug-o
	@rm -f $(OPT_OBJS) $(OPT_EXECUTABLES) $(OBJDIR)/*.opt-o

clean-optimized:
	@echo Cleaning up optimized
	@rm -f $(OBJS) $(EXECUTABLES) $(OBJDIR)/*.o

clean: clean-test.sh-files
	@echo Cleaning up
	@rm -f $(DEPS) $(WINDOWS_INSTALLER_OBJS) *.d .archive.tar.gz *.stackdump $(EXECUTABLES) $(OPT_EXECUTABLES) $(DEBUG_EXECUTABLES)
	@rm -f *.good *.bad *.local *.b250 test/*.good test/*.bad test/*.local test/*.b250
	@rm -Rf $(OBJDIR)

.PHONY: clean clean-debug clean-optimized clean-test.sh-files git-pull macos mac/.remote_mac_timestamp delete-arch docs testfiles

