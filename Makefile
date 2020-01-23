# ------------------------------------------------------------------
#   Makefile
#   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
#   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

# Note for Windows: to run this make, you need mingw (for the gcc compiler) and cygwin (for Unix-like tools):
# Mingw: http://mingw-w64.org/doku.php  (Windows 32 bit version also works)
# Cygwin: https://www.cygwin.com/

ifdef BUILD_PREFIX
IS_CONDA=1
endif

CFLAGS       += -D_LARGEFILE64_SOURCE=1 -Wall -I.
LDFLAGS      += -lpthread -lm

ifdef IS_CONDA 
	LDFLAGS += -lbz2 -lz  # conda - dynamic linking with bz2 and zlib

	ifeq ($(OS),Windows_NT)
		CC=gcc # in Windows, override conda's default Visual C with gcc 
        LDFLAGS += -L$(PREFIX)/Library/lib
		CFLAGS  += -I$(PREFIX)/Library/include 
	endif
else
	CC=gcc
	CFLAGS += -Izlib -Ibzlib
endif 

ifeq ($(CC),gcc)
	CFLAGS += -Ofast
else
	CFLAGS += -O2
endif

MY_SRCS = genozip.c base250.c move_to_front.c vcf_header.c zip.c piz.c gloptimize.c buffer.c \
	      vcffile.c squeeze.c zfile.c segregate.c profiler.c file.c vb.c dispatcher.c crypt.c aes.c md5.c bzlib_mod.c

CONDA_COMPATIBILITY_SRCS = compatability/win32_pthread.c compatability/visual_c_gettime.c compatability/visual_c_misc_funcs.c compatability/mac_gettime.c

EXT_SRCS = bzlib/blocksort.c bzlib/bzlib.c bzlib/compress.c bzlib/crctable.c bzlib/decompress.c bzlib/huffman.c bzlib/randtable.c \
           zlib/gzlib.c zlib/gzread.c zlib/inflate.c zlib/inffast.c zlib/zutil.c zlib/inftrees.c zlib/crc32.c zlib/adler32.c 
                
CONDA_DEVS = Makefile .gitignore test-file.vcf 

CONDA_DOCS = LICENSE.non-commercial.txt LICENSE.commercial.txt AUTHORS README.md

CONDA_INCS = genozip.h lic-text.h  \
             compatability/visual_c_getopt.h compatability/visual_c_stdbool.h compatability/visual_c_unistd.h \
	         compatability/visual_c_gettime.h compatability/visual_c_stdint.h compatability/visual_c_misc_funcs.h \
	         compatability/win32_pthread.h \
      	     compatability/mac_gettime.h  # doesn't include version.h bc it would create a circular dependency with .version

ifeq ($(CC),cl)
	MY_SRCS += compatability/visual_c_gettime.c compatability/visual_c_misc_funcs.c 
endif

ifeq ($(OS),Windows_NT)
# Windows
	EXE = .exe
	MY_SRCS += compatability/win32_pthread.c
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
# Linux
        LDFLAGS += -lrt -s
    endif
    ifeq ($(UNAME_S),Darwin)
# Mac
		MY_SRCS += compatability/mac_gettime.c
    endif
endif

ifndef IS_CONDA 
# local - static link everything
SRCS = $(MY_SRCS) $(EXT_SRCS)
else 
# conda - use packages for zlib and bzip2
SRCS = $(MY_SRCS) 
endif

OBJS       := $(SRCS:.c=.o)
DEBUG_OBJS := $(SRCS:.c=.debug-o)

DEPS       := $(SRCS:.c=.d)

all: genozip$(EXE) genounzip$(EXE) genocat$(EXE) genols$(EXE)

debug : CFLAGS  += -DDEBUG -g
debug : genozip-debug$(EXE)

-include $(DEPS)

%.d: %.c
	@echo Calculating dependencies $<
	@$(CC) $(CFLAGS) -MM -MT $@ $< -MF $(@:%.o=%.d)

%.o: %.c %.d
	@echo Compiling $<
	@$(CC) -c -o $@ $< $(CFLAGS)

%.debug-o: %.c %.d
	@echo "Compiling $< (debug)"
	@$(CC) -c -o $@ $< $(CFLAGS)

genozip$(EXE): $(OBJS)
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)
 
genozip-debug$(EXE): $(DEBUG_OBJS)
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

genounzip$(EXE) genocat$(EXE) genols$(EXE): genozip$(EXE)
	@echo Hard linking $@
	@rm -f $@ 
	@ln $^ $@

LICENSE.non-commercial.txt: lic-text.h
	@echo Generating $@
	@./genozip$(EXE) --license > $@

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

# currently, I build for conda from my Windows machine so I don't bother supporting other platforms
ifeq ($(OS),Windows_NT)

# increments minor version, eg. 1.0.1 -> 1.0.2. 
# To increment a major version, manually edit .version and set minor version to -1 e.g. 1.1.-1 (careful! no newlines or spaces)
# and re-compile so that genozip --version gets update
# IMPORTANT: the first number in the version indicates the genozip file format version and goes into
# the genozip file header SectionHeaderVCFHeader.genozip_version
increment-version: $(MY_SRCS) $(EXT_SRCS) $(CONDA_COMPATIBILITY_SRCS) $(CONDA_DEVS) $(CONDA_DOCS) $(CONDA_INCS)
	@if (( `git status|grep 'Changes not staged for commit\|Untracked files'|wc -l` > 0 )) ; then echo "Making $@: ERROR: Please 'git commit' everything first" ; exit 1 ; fi
	@echo $(shell cut -d. -f1-2 .version).$(shell expr 1 + `cut -d. -f3 .version`) > .version 
	@git commit -m "increment version" .version
	
version.h : .version
	@echo \#define GENOZIP_CODE_VERSION \"$(shell cat .version)\"             > $@   # override previous
	@echo \#define GENOZIP_FILE_FORMAT_VERSION $(shell cut -d. -f1 .version) >> $@
	@git commit -m "increment version" version.h

.archive.tar.gz : increment-version version.h
	@if (( `git status|grep 'Changes not staged for commit\|Untracked files'|wc -l` > 0 )); then echo "Making $@: ERROR: Please 'git commit' everything first" ; exit 1 ; fi
	@echo Creating github tag genozip-$(shell cat .version) and archive
	@git push 
	@git tag genozip-$(shell cat .version)
	@git push origin genozip-$(shell cat .version)
	@curl https://github.com/divonlan/genozip/archive/genozip-$(shell cat .version).tar.gz --silent --location -o $@

conda/meta.yaml: conda/meta.template.yaml .archive.tar.gz 
	@echo "Generating meta.yaml (for conda)"
	@cat conda/meta.template.yaml | \
		sed s/SHA256/$(shell openssl sha256 .archive.tar.gz | cut -d= -f2 | cut -c2-)/ | \
		sed s/VERSION/$(shell cat .version)/g | \
		grep -v "^#" \
		> $@
 
#CONDA_RECIPE_DIR = ../staged-recipes/recipes/genozip # initial stage-recipes step, keeping here for future reference
CONDA_RECIPE_DIR = ../genozip-feedstock/recipe

# publish to conda-forge 
conda/.conda-timestamp: conda/meta.yaml conda/build.sh conda/bld.bat
	@if (( `git status|grep 'Changes not staged for commit\|Untracked files'|wc -l` > 0 )); then echo "Making $@: ERROR: Please 'git commit' everything first" ; exit 1 ; fi
	@echo " "
#	@echo Rebasing my staged-recipes fork, and pushing changes to genozip branch of the fork
#	@(cd ../staged-recipes/; git checkout master; git pull --rebase upstream master ; git push origin master --force ; git checkout genozip)  # needed for initial stage-recipes step, keeping here for future reference
	@echo " "
	@echo "Copying meta.yaml build.sh bld.bat to conda-forge"
	@cp conda/meta.yaml conda/build.sh conda/bld.bat $(CONDA_RECIPE_DIR)
	@echo "Committing my files to the branch and pushing them"
	@(cd $(CONDA_RECIPE_DIR); git commit -m "update" meta.yaml build.sh bld.bat; git push)
	@echo " "
	@echo "Submitting pull request to conda-forge"
#	@(cd $(CONDA_RECIPE_DIR); git request-pull master https://github.com/divonlan/staged-recipes genozip)
	@(cd $(CONDA_RECIPE_DIR); git request-pull master https://github.com/divonlan/genozip-feedstock master)
	@touch $@
	@echo " "
	@echo "Check status on: https://dev.azure.com/conda-forge/feedstock-builds/_build"
#	@echo "and: https://github.com/conda-forge/staged-recipes/pull/10617"
	@echo "and: https://github.com/conda-forge/genozip-feedstock/pulls"

WINDOWS_INSTALL_FILES = windows/genozip.exe windows/genounzip.exe windows/genocat.exe windows/genols.exe LICENSE.commercial.txt LICENSE.non-commercial.txt windows/readme.txt test-file.vcf

windows/%.exe: %.exe
	@echo Copying $<
	@cp -f $< $@ 

windows/readme.txt: genozip$(EXE)
	@echo Generating $@
	@./genozip$(EXE) --help > $@
	
windows/LICENSE.for-installer.txt: lic-text.h
	@echo Generating $@
	@./genozip$(EXE) --license --force > $@

# this must be run ONLY has part of "make distribution" or else versions will be out of sync
windows/genozip-installer.exe: $(WINDOWS_INSTALL_FILES) windows/LICENSE.for-installer.txt
	@echo 'Committing Windows files and pushing all changes to repo'
	@git stage $(WINDOWS_INSTALL_FILES) $@
	@(git commit -m windows_files_for_version_$(shell cat .version) $(WINDOWS_INSTALL_FILES) $@ ; exit 0)
	@echo Verifying that all files are committed to the repo
	@(exit `git status|grep 'Changes not staged for commit\|Untracked files'|wc -l`)
	@echo 'Using the UI:'
	@echo '  (1) Open windows/genozip.ifp'
	@echo '  (2) Set General-Program version to $(shell cat .version)'
	@echo '  (3) Verify the files Setup-Files, and the license from LICENSE.for-installer.txt are up to date'
	@echo '  (4) Click Save, then click Build, then click No to the popup question'
	@echo '  (5) Exit the UI (close the window)'
	@(C:\\\\Program\\ Files\\ \\(x86\\)\\\\solicus\\\\InstallForge\\\\InstallForge.exe ; exit 0)
	@(git stage windows/genozip.ifp $@ ; exit 0)
	@(git commit -m windows_files_for_version_$(shell cat .version) $(WINDOWS_INSTALL_FILES) windows/genozip.ifp $@ ; exit 0)
	@git push


endif

.PHONY: clean clean-debug clean-all increment-version

distribution: conda/.conda-timestamp windows/genozip-installer.exe

clean:
	@echo Cleaning up
	@rm -f $(DEPS) $(OBJS) genozip$(EXE) genounzip$(EXE) genocat$(EXE) genols$(EXE) 

clean-debug:
	@echo Cleaning up debug
	@rm -f $(DEPS) $(DEBUG_OBJS) genozip-debug$(EXE) 
