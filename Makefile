# ------------------------------------------------------------------
#   Makefile
#   Copyright (C) 2019-2020 Divon Lan <genozip@blackpawventures.com>
#   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

# Note for Windows: to run this make, you need mingw (for the gcc compiler) and cygwin (for Unix-like tools):
# Mingw: http://mingw-w64.org/doku.php  (Windows 32 bit version also works)
# Cygwin: https://www.cygwin.com/

VERSION = 1.0.0

# use gcc, unless its conda - let it define its own compiler
ifndef CONDA_BUILD_SYSROOT 
CC=gcc
endif

CFLAGS       = -Ibzlib -Izlib -D_LARGEFILE64_SOURCE=1 -Wall -Ofast 
CFLAGS_DEBUG = -Ibzlib -Izlib -D_LARGEFILE64_SOURCE=1 -Wall -DDEBUG -g
LDFLAGS = -lpthread -lm

DEVS = Makefile .gitignore genozip.code-workspace \
       .vscode/c_cpp_properties.json .vscode/launch.json .vscode/settings.json .vscode/tasks.json \
       test-file.vcf 

DOCS = LICENSE.non-commercial.txt LICENSE.commercial.txt AUTHORS README.md \
       bzlib/LICENSE.bzlib bzlib/README.bzlib \
	   zlib/README.zlib

INCS = genozip.h lic-text.h \
       bzlib/bzlib.h bzlib/bzlib_private.h \
	   zlib/crc32.h zlib/gzguts.h zlib/inffast.h zlib/inffixed.h zlib/inflate.h zlib/inftrees.h zlib/zconf.h zlib/zlib.h zlib/zutil.h \
	   mac/mach_gettime.h

SRCS = base250.c move_to_front.c vcf_header.c zip.c piz.c gloptimize.c buffer.c main.c \
	   vcffile.c squeeze.c zfile.c segregate.c profiler.c file.c vb.c dispatcher.c \
       bzlib/blocksort.c bzlib/bzlib.c bzlib/compress.c bzlib/crctable.c bzlib/decompress.c bzlib/huffman.c bzlib/randtable.c \
       zlib/gzlib.c zlib/gzread.c zlib/inflate.c zlib/inffast.c zlib/zutil.c zlib/inftrees.c zlib/crc32.c zlib/adler32.c \
       mac/mach_gettime.c

OBJS       := $(SRCS:.c=.o)
DEBUG_OBJS := $(SRCS:.c=.debug-o)

DEPS       := $(SRCS:.c=.d)

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

all: genozip$(EXE) genounzip$(EXE) genocat$(EXE)

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

all: genozip$(EXE) genounzip$(EXE) genocat$(EXE)

genozip$(EXE): $(OBJS)
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

genozip-debug$(EXE): $(DEBUG_OBJS)
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

genounzip$(EXE) genocat$(EXE): genozip$(EXE)
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

TARBALL := archive/genozip-$(VERSION).tar.gz

$(TARBALL): $(SRCS) $(INCS) $(DOCS) $(DEVS)
	@echo "Archiving to $@"
	@tar --create --gzip --file $(TARBALL) $^
	@echo "Committing $(TARBALL) & pushing changes to genozip/master"
	@(git commit -m "update archive" $(TARBALL) ; git push)

# currently, I build for conda from my Windows machine so I don't bother supporting other platforms
ifeq ($(OS),Windows_NT)

conda/meta.yaml: conda/meta.yaml.template $(TARBALL)
	@echo "Generating meta.yaml (for conda)"
	@cat conda/meta.yaml.template | \
		sed s/'{{ sha256 }}'/$$(openssl sha256 $(TARBALL)|cut -d= -f2|cut -c2-)/ | \
		sed s/'{{ version }}'/$(VERSION)/g | \
		grep -v "^#" \
		> $@

# we make the build scripts dependents on the archives, so if file list changes, we need to re-generate build
#UNIX_SRCS := $(shell echo $(SRCS) | sed 's/\\//\\\\\\//g' ) # a list of files that look like: zlib\/inflate.c
#build.sh: conda/build.sh.template 
#	@echo "Generating $@ (for conda)"
#	@sed 's/%BUILD/\\$$CC $(CFLAGS) $(UNIX_SRCS) -o genozip/' conda/$@.template > $@
	
#WIN_SRCS  := $(shell echo $(SRCS) | sed 's/\\//\\\\\\\\\\\\\\\\/g' ) # crazy! we need 16 blackslashes to end up with a single one in the bld.bat file
#bld.bat: conda/bld.bat.template
#	@echo "Generating $@ (for conda)"
#	@sed 's/%BUILD/%GCC% $(CFLAGS) $(WIN_SRCS) -o genozip.exe/' conda/$@.template > $@

# publish to conda-forge
conda: $(TARBALL) meta.yaml build.sh bld.bat
	@echo "Copying meta.yaml build.sh bld.bat to staged-recipes"
	@cp conda/meta.yaml conda/build.sh conda/bld.bat ../staged-recipes/recipes/genozip/
	@echo "Committing files & pushing changes to my forked staged-recipies/genozip-branch"
	@(cd ../staged-recipes/recipes/genozip; git commit -m "update" meta.yaml build.sh bld.bat; git push)
	@echo "Submitting pull request to conda-forge"
	@(cd ../staged-recipes/recipes/genozip; git request-pull master https://github.com/divonlan/staged-recipes genozip-branch)
	@echo " "
	@echo "Check status on: https://dev.azure.com/conda-forge/feedstock-builds/_build"
	@echo "(if you don't see it there, try https://github.com/divonlan/staged-recipes - select Branch: genozip-branch + New pull request)"

endif

LICENSE.non-commercial.txt: genozip$(EXE)
	@echo Generating $@
	@./genozip$(EXE) -L > $@

.PHONY: clean clean-debug clean-all build.sh bld.bat meta.yaml

clean:
	@echo Cleaning up
	@rm -f $(DEPS) $(OBJS) genozip$(EXE) genounzip$(EXE) genocat$(EXE) 

clean-debug:
	@echo Cleaning up debug
	@rm -f $(DEPS) $(DEBUG_OBJS) genozip-debug$(EXE) 
