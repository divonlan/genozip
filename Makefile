# ------------------------------------------------------------------
#   Makefile
#   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
#   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

# Note for Windows: to run this make, you need mingw (for the gcc compiler) and cygwin (for Unix-like tools):
# Mingw: http://mingw-w64.org/doku.php  (Windows 32 bit version also works)
# Cygwin: https://www.cygwin.com/

VERSION = 1.0.0

EXE =
ifeq ($(OS),Windows_NT)
	EXE = .exe
endif

CC=gcc
CFLAGS       = -Ibzlib -Izlib -D_LARGEFILE64_SOURCE=1 -Wall -Ofast
CFLAGS_DEBUG = -Ibzlib -Izlib -D_LARGEFILE64_SOURCE=1 -Wall -DDEBUG -g
LIBS = -lpthread -lm

DEVS = Makefile .gitignore vczip.code-workspace \
       .vscode/c_cpp_properties.json .vscode/launch.json .vscode/settings.json .vscode/tasks.json \
       data/test-file.vcf \
       conda/build.sh.template conda/bld.bat.template conda/meta.yaml.template

DOCS = LICENSE.non-commercial.txt LICENSE.commercial.txt AUTHORS README.md \
       bzlib/LICENSE.bzlib bzlib/README.bzlib \
	   zlib/README.zlib

INCS = vczip.h license.h \
       bzlib/bzlib.h bzlib/bzlib_private.h \
	   zlib/crc32.h zlib/gzguts.h zlib/inffast.h zlib/inffixed.h zlib/inflate.h zlib/inftrees.h zlib/zconf.h zlib/zlib.h zlib/zutil.h \
	   mac/mach_gettime.h

SRCS = vcf_header.c zip.c piz.c gloptimize.c buffer.c main.c vcffile.c squeeze.c zfile.c segregate.c profiler.c file.c vb.c \
       bzlib/blocksort.c bzlib/bzlib.c bzlib/compress.c bzlib/crctable.c bzlib/decompress.c bzlib/huffman.c bzlib/randtable.c \
       zlib/gzlib.c zlib/gzread.c zlib/inflate.c zlib/inffast.c zlib/zutil.c zlib/inftrees.c zlib/crc32.c zlib/adler32.c \
       mac/mach_gettime.c

OBJS       := $(SRCS:.c=.o)
DEBUG_OBJS := $(SRCS:.c=.debug-o)

DEPS       := $(SRCS:.c=.d)

all: vczip$(EXE) vcpiz$(EXE) vccat$(EXE)

debug: vczip-debug$(EXE)

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

all: vczip$(EXE) vcpiz$(EXE) vccat$(EXE)

vczip$(EXE): $(OBJS)
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

vczip-debug$(EXE): $(DEBUG_OBJS)
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

vcpiz$(EXE) vccat$(EXE): vczip$(EXE)
	@echo Hard linking $@
	@rm -f $@ 
	@ln $^ $@

TARBALL := archive/vczip-$(VERSION).tar.gz

$(TARBALL): $(SRCS) $(INCS) $(DOCS) $(DEVS)
	@echo "Archiving to $@ - WARNING: Make sure you have no un-pushed changes locally (Makefile doesn't verify this)"
	@tar --create --gzip --file $(TARBALL) $^

meta.yaml: conda/meta.yaml.template $(TARBALL)
	@echo "Generating meta.yaml (for conda)"
	@cat conda/meta.yaml.template | \
		sed s/%SHA256/$$(openssl sha256 $(TARBALL)|cut -d= -f2|cut -c2-)/ | \
		sed s/%VERSION/$(VERSION)/g | \
		grep -v "^#" \
		> $@

# we make the build scripts dependents on the archives, so if file list changes, we need to re-generate build
UNIX_SRCS := $(shell echo $(SRCS) | sed 's/\\//\\\\\\//g' ) # a list of files that look like: zlib\/inflate.c
build.sh: $(TARBALL) conda/build.sh.template 
	@echo "Generating $@ (for conda)"
	@sed 's/%BUILD/$(CC) $(CFLAGS) $(LIBS) $(UNIX_SRCS) -o vczip/' conda/$@.template > $@
	
WIN_SRCS  := $(shell echo $(SRCS) | sed 's/\\//\\\\\\\\\\\\\\\\/g' ) # crazy! we need 16 blackslashes to end up with a single one in the bld.bat file
bld.bat: $(TARBALL) conda/bld.bat.template
	@echo "Generating $@ (for conda)"
	@sed 's/%BUILD/$(CC).exe $(CFLAGS) $(LIBS) $(WIN_SRCS) -o vczip.exe/' conda/$@.template > $@

conda: $(TARBALL) meta.yaml build.sh bld.bat

LICENSE.non-commercial.txt: vczip$(EXE)
	@echo Generating $@
	@./vczip$(EXE) -L > $@

.PHONY: clean clean-debug clean-all

clean:
	@echo Cleaning up
	@rm -f $(OBJS) vczip$(EXE) vcpiz$(EXE) vccat$(EXE) 

clean-debug:
	@echo Cleaning up debug
	@rm -f $(DEBUG_OBJS) vczip-debug$(EXE) 

clean-all: clean clean-debug
	@echo Cleaning up dependencies
	@rm -f $(DEPS) 
