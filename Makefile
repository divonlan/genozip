# ------------------------------------------------------------------
#   Makefile
#   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
#   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

# Note for Windows: to run this make, you need mingw (for the gcc compiler) and cygwin (for Unix-like tools):
# Mingw: http://mingw-w64.org/doku.php  (Windows 32 bit version also works)
# Cygwin: https://www.cygwin.com/

VERSION = 1.0.0

# use gcc, unless its conda - let it define its own compiler
ifndef CONDA_BUILD_SYSROOT 
CC=gcc
endif

CFLAGS       = -D_LARGEFILE64_SOURCE=1 -Wall 
CFLAGS_DEBUG = -D_LARGEFILE64_SOURCE=1 -Wall -DDEBUG -g
LDFLAGS      = -lpthread -lm

ifeq ($(CC),gcc)
	CFLAGS += -Ofast
else
	CFLAGS += -O2
endif

DEVS = Makefile .gitignore genozip.code-workspace \
       .vscode/c_cpp_properties.json .vscode/launch.json .vscode/settings.json .vscode/tasks.json \
       test-file.vcf 

DOCS = LICENSE.non-commercial.txt LICENSE.commercial.txt AUTHORS README.md \
       bzlib/LICENSE.bzlib bzlib/README.bzlib \
	   zlib/README.zlib

INCS = genozip.h lic-text.h \
       bzlib/bzlib.h bzlib/bzlib_private.h \
	   zlib/crc32.h zlib/gzguts.h zlib/inffast.h zlib/inffixed.h zlib/inflate.h zlib/inftrees.h zlib/zconf.h zlib/zlib.h zlib/zutil.h \
	   crypto/aes.h crypto/WjCryptLib_Md5.h \
       compatability/visual_c_getopt.h compatability/visual_c_stdbool.h compatability/visual_c_unistd.h\
	   compatability/mac_gettime.h \
	   compatability/win32_pthread.h compatability/visual_c_gettime.h \
	   compatability/visual_c_stdint.h compatability/visual_c_misc_funcs.h

OLD_C_SRCS = bzlib/blocksort.c bzlib/bzlib.c bzlib/compress.c bzlib/crctable.c bzlib/decompress.c bzlib/huffman.c bzlib/randtable.c \
       zlib/gzlib.c zlib/gzread.c zlib/inflate.c zlib/inffast.c zlib/zutil.c zlib/inftrees.c zlib/crc32.c zlib/adler32.c \
	   crypto/aes.c crypto/WjCryptLib_Md5.c
       
C99_SRCS = genozip.c base250.c move_to_front.c vcf_header.c zip.c piz.c gloptimize.c buffer.c \
	   vcffile.c squeeze.c zfile.c segregate.c profiler.c file.c vb.c dispatcher.c crypt.c \
       compatability/mac_gettime.c compatability/win32_pthread.c compatability/visual_c_gettime.c \
	   compatability/visual_c_misc_funcs.c

SRCS = $(C99_SRCS) $(OLD_C_SRCS)

ifeq ($(OS),Windows_NT)
# Windows
	EXE = .exe
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
# Linux
        LDFLAGS += -lrt -s
    endif
    ifeq ($(UNAME_S),Darwin)
# Mac
    endif
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
	@$(CC) -c -o $@ $< $(CFLAGS_DEBUG)

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

# this is used by build.sh and bld.bat to build on conda
install: genozip$(EXE)
	@echo "Installing in $(PREFIX)/bin"
	@if ( test ! -d $(PREFIX)/bin ) ; then mkdir -p $(PREFIX)/bin ; fi
	@cp -f genozip$(EXE) $(PREFIX)/bin/genozip$(EXE)
ifneq ($(OS),Windows_NT)
	@chmod a+x $(PREFIX)/bin/genozip$(EXE)
endif
	@cp -f $(PREFIX)/bin/genozip$(EXE) $(PREFIX)/bin/genounzip$(EXE)
	@cp -f $(PREFIX)/bin/genozip$(EXE) $(PREFIX)/bin/genocat$(EXE)
	@cp -f $(PREFIX)/bin/genozip$(EXE) $(PREFIX)/bin/genols$(EXE)

TARBALL := conda/genozip-$(VERSION).tar.gz

$(TARBALL): $(SRCS) $(INCS) $(DOCS) $(DEVS) 
	@echo "Archiving to $@"
	@tar --create --gzip --file $(TARBALL) $^
	@echo "Committing $(TARBALL) & pushing changes to genozip/master"
	@(git commit -m "update archive" $(TARBALL) ; git push)

# currently, I build for conda from my Windows machine so I don't bother supporting other platforms
ifeq ($(OS),Windows_NT)
 
conda/meta.yaml: conda/meta.template.yaml $(TARBALL)
	@echo "Generating meta.yaml (for conda)"
	@cat conda/meta.template.yaml | \
		sed s/'{{ sha256 }}'/$$(openssl sha256 $(TARBALL)|cut -d= -f2|cut -c2-)/ | \
		sed s/'{{ version }}'/$(VERSION)/g | \
		grep -v "^#" \
		> $@

OLD_C_COMPILE_AS_C := $(shell echo $(OLD_C_SRCS) | tr ' ' '\n' | sed 's/^/-Tc /g'|tr '\n' ' ')
C99_COMPILE_AS_CPP := $(shell echo $(C99_SRCS)   | tr ' ' '\n' | sed 's/^/-Tp /g'|tr '\n' ' ')

OLD_C_WIN_SRCS  := $(shell echo $(OLD_C_COMPILE_AS_C) | sed 's/\\//\\\\\\\\\\\\\\\\/g' ) # crazy! we need 16 blackslashes to end up with a single one in the bld.bat file
C99_WIN_SRCS    := $(shell echo $(C99_COMPILE_AS_CPP) | sed 's/\\//\\\\\\\\\\\\\\\\/g' ) 

conda/bld.bat: conda/bld.template.bat Makefile
	@echo "Building $@ (for conda)"
	@sed s/'{{ src }}'/'$(C99_WIN_SRCS) $(OLD_C_WIN_SRCS)'/ conda/bld.template.bat > $@

# publish to conda-forge
conda: $(TARBALL) conda/meta.yaml conda/build.sh conda/bld.bat
	@echo "Copying meta.yaml build.sh bld.bat to staged-recipes"
	@cp conda/meta.yaml conda/build.sh conda/bld.bat ../staged-recipes/recipes/genozip/
	@echo "Committing files & pushing changes to my forked staged-recipies/genozip-branch"
	@(cd ../staged-recipes/recipes/genozip; git commit -m "update" meta.yaml build.sh bld.bat; git push)
	@echo "Submitting pull request to conda-forge"
	@(cd ../staged-recipes/recipes/genozip; git request-pull master https://github.com/divonlan/staged-recipes genozip-branch)
	@echo " "
	@echo "Check status on: https://dev.azure.com/conda-forge/feedstock-builds/_build"
	@echo "and: https://github.com/conda-forge/staged-recipes/pull/10543"
	@echo "(if you don't see it there, try https://github.com/divonlan/staged-recipes - select Branch: genozip-branch + New pull request)"

WINDOWS_INSTALL_FILES = genozip.exe genounzip.exe genocat.exe genols.exe LICENSE.commercial.txt LICENSE.non-commercial.txt windows/readme.txt test-file.vcf

windows-installer: $(WINDOWS_INSTALL_FILES) windows/LICENSE.for-installer.txt
	@echo 'Copying files to windows'
	@cp genozip.exe genounzip.exe genocat.exe genols.exe windows
	@echo 'Using the UI:'
	@echo '  (1) Open windows/genozip.ifp'
	@echo '  (2) Make sure the version is set to $(VERSION)'
	@echo '  (3) Make sure the files list, and the license from LICENSE.for-installer.txt are up to date'
	@echo '  (4) Click Build'
	@C:\\Program\ Files\ \(x86\)\\solicus\\InstallForge\\InstallForge.exe
	@echo


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

.PHONY: clean clean-debug clean-all 

clean:
	@echo Cleaning up
	@rm -f $(DEPS) $(OBJS) genozip$(EXE) genounzip$(EXE) genocat$(EXE) genols$(EXE) 

clean-debug:
	@echo Cleaning up debug
	@rm -f $(DEPS) $(DEBUG_OBJS) genozip-debug$(EXE) 
