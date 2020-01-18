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

# use gcc, unless its conda - let it define its own compiler (but for Windows, we override conda's default Visual C with gcc)
ifndef IS_CONDA 
CC=gcc
else ifeq ($(OS),Windows_NT)
CC=gcc
endif 

VERSION_1 := $(shell cat .version | cut -d. -f1)  # eg 1      -- this goes into the file header as the version of the file format
VERSION_2 := \"$(shell cat .version)\"            # eg "1.0"  -- this is reported in genozip --version

CFLAGS       += -D_LARGEFILE64_SOURCE=1 -Wall -I. -DVERSION_1=$(VERSION_1) -DVERSION_2=$(VERSION_2)
DEBUG_CLAGS  += -D_LARGEFILE64_SOURCE=1 -Wall -I. -DVERSION_1=$(VERSION_1) -DVERSION_2=$(VERSION_2) -Izlib -Ibzlib -DDEBUG -g
LDFLAGS      += -lpthread -lm

ifdef IS_CONDA 
# conda - dynamic linking with bz2 and zlib
    LDFLAGS += -lbz2 -lz 
    ifdef LIB
        LDFLAGS += -L$(PREFIX)/Library/lib     # this appears (at least) in the Windows build in conda-forge
    endif
	ifdef INCLUDE
		CFLAGS  += -I$(PREFIX)/Library/include  # this appears (at least) in the Windows build in conda-forge
	endif
else
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
                
CONDA_DEVS = Makefile .gitignore test-file.vcf .version

CONDA_DOCS = LICENSE.non-commercial.txt LICENSE.commercial.txt AUTHORS README.md

CONDA_INCS = genozip.h lic-text.h \
       compatability/visual_c_getopt.h compatability/visual_c_stdbool.h compatability/visual_c_unistd.h\
	   compatability/mac_gettime.h \
	   compatability/win32_pthread.h compatability/visual_c_gettime.h \
	   compatability/visual_c_stdint.h compatability/visual_c_misc_funcs.h

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

debug: genozip-debug$(EXE)

-include $(DEPS)

%.d: %.c
	@echo Calculating dependencies $<
	@$(CC) $(CFLAGS) -MM -MT $@ $< -MF $(@:%.o=%.d)

%.o: %.c %.d
	@echo Compiling $<
	@$(CC) -c -o $@ $< $(CFLAGS)

%.debug-o: %.c %.d
	@echo "Compiling $< (debug)"
	@$(CC) -c -o $@ $< $(DEBUG_CLAGS)

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
.version: 
# double check that everything is committed (we check several times)
	@echo Verifying that all files are committed to the repo
	@(exit `git status|grep 'Changes not staged for commit\|Untracked files'|wc -l`)
	@echo Verifying that something has changed since version $(shell cat .version)
## THIS DOESN'T WORK
	@([ `git log $(shell cat .version)..HEAD --oneline | wc -l` ] ; exit $$?)
	@echo $(shell cut -d. -f1-2 $@).$(shell expr 1 + `cut -d. -f3 $@`) > $@
	@git commit -m "increment version" .version

.archive.tar.gz: .version
# double check that everything is committed (we check several times)
	@echo Verifying that all files are committed to the repo
	@(exit `git status|grep 'Changes not staged for commit\|Untracked files'|wc -l`)
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
 
OLD_C_COMPILE_AS_C := $(shell echo $(OLD_C_SRCS) | tr ' ' '\n' | sed 's/^/-Tc /g'|tr '\n' ' ')
C99_COMPILE_AS_CPP := $(shell echo $(C99_SRCS)   | tr ' ' '\n' | sed 's/^/-Tp /g'|tr '\n' ' ')

OLD_C_WIN_SRCS  := $(shell echo $(OLD_C_COMPILE_AS_C) | sed 's/\\//\\\\\\\\\\\\\\\\/g' ) # crazy! we need 16 blackslashes to end up with a single one in the bld.bat file
C99_WIN_SRCS    := $(shell echo $(C99_COMPILE_AS_CPP) | sed 's/\\//\\\\\\\\\\\\\\\\/g' ) 
 
# publish to conda-forge 
conda/.conda-timestamp: conda/meta.yaml conda/build.sh conda/bld.bat $(MY_SRCS) $(CONDA_INCS) $(CONDA_DOCS) $(CONDA_DEVS) $(CONDA_COMPATIBILITY_SRCS)
# double check that everything is committed (we check several times)
	@echo Verifying that all files are committed to the repo
	@(exit `git status|grep 'Changes not staged for commit\|Untracked files'|wc -l`)
	@echo " "
	@echo Rebasing my staged-recipes fork, and pushing changes to genozip branch of the fork
	@(cd ../staged-recipes/; git checkout master; git pull --rebase upstream master ; git push origin master --force ; git checkout genozip)
	@echo " "
	@echo "Copying meta.yaml build.sh bld.bat to staged-recipes"
	@cp conda/meta.yaml conda/build.sh conda/bld.bat ../staged-recipes/recipes/genozip/
	@echo "Committing my files to the branch and pushing them"
	@(cd ../staged-recipes/recipes/genozip; git commit -m "update" meta.yaml build.sh bld.bat; git push)
	@echo " "
	@echo "Submitting pull request to conda-forge"
	@(cd ../staged-recipes/recipes/genozip; git request-pull master https://github.com/divonlan/staged-recipes genozip)
	@touch $@
	@echo " "
	@echo "Check status on: https://dev.azure.com/conda-forge/feedstock-builds/_build"
	@echo "and: https://github.com/conda-forge/staged-recipes/pull/10617"
	@echo "(if you don't see it there, try https://github.com/divonlan/staged-recipes - select Branch: genozip-branch + New pull request)"

WINDOWS_INSTALL_FILES = windows/genozip.exe windows/genounzip.exe windows/genocat.exe windows/genols.exe LICENSE.commercial.txt LICENSE.non-commercial.txt windows/readme.txt test-file.vcf

windows/%.exe: %.exe
	@echo Copying $<
	@cp -f $< $@ 

# this must be run AFTER conda and BEFORE any other changes - or else the version will not yet be updated
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
	@echo '  (4) Click Save, then click Build'
	@echo '  (5) Exit the UI (close the window)'
	@(C:\\\\Program\\ Files\\ \\(x86\\)\\\\solicus\\\\InstallForge\\\\InstallForge.exe ; exit 0)
	@(git stage windows/genozip.ifp $@ ; exit 0)
	@(git commit -m windows_files_for_version_$(shell cat .version) $(WINDOWS_INSTALL_FILES) $@ ; exit 0)
	@git push


endif

windows/LICENSE.for-installer.txt: lic-text.h
	@echo Generating $@
	@./genozip$(EXE) --license --force > $@

LICENSE.non-commercial.txt: lic-text.h
	@echo Generating $@
	@./genozip$(EXE) --license > $@

windows/readme.txt: genozip$(EXE)
	@echo Generating $@
	@./genozip$(EXE) --help > $@

.PHONY: clean clean-debug clean-all .version

clean:
	@echo Cleaning up
	@rm -f $(DEPS) $(OBJS) genozip$(EXE) genounzip$(EXE) genocat$(EXE) genols$(EXE) 

clean-debug:
	@echo Cleaning up debug
	@rm -f $(DEPS) $(DEBUG_OBJS) genozip-debug$(EXE) 
