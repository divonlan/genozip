# ------------------------------------------------------------------
#   Makefile
#   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
#   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

CC=gcc
CFLAGS       = -Ibzlib -Izlib -D_LARGEFILE64_SOURCE=1 -Wall -Ofast
CFLAGS_DEBUG = -Ibzlib -Izlib -D_LARGEFILE64_SOURCE=1 -Wall -DDEBUG -g
LIBS = -lpthread -lm

SLASH :=
ifeq ($(OS),Windows_NT)
	SLASH := \\
else
	SLASH := /
endif

RM :=
ifeq ($(OS),Windows_NT)
	RM := del /q /f 
else
	RM := rm -f
endif

AFTER_RM :=
ifeq ($(OS),Windows_NT)
	AFTER_RM := >nul 2>&1
endif

MAC   = mac$(SLASH)
ZLIB  = zlib$(SLASH)
BZLIB = bzlib$(SLASH)

INCS = vczip.h licenese.h \
       $(BZLIB)bzlib.h $(BZLIB)bzlib_private.h \
	   $(ZLIB)crc32.h $(ZLIB)gzguts.h $(ZLIB)inffast.h $(ZLIB)inffixed.h $(ZLIB)inflate.h $(ZLIB)inftrees.h $(ZLIB)zconf.h $(ZLIB)zlib.h $(ZLIB)zutil.h \
	   $(MAC)mach_gettime.h

SRCS = vcf_header.c zip.c piz.c gloptimize.c buffer.c main.c vcffile.c squeeze.c zfile.c segregate.c profiler.c file.c vb.c \
       $(BZLIB)blocksort.c $(BZLIB)bzlib.c $(BZLIB)compress.c $(BZLIB)crctable.c $(BZLIB)decompress.c $(BZLIB)huffman.c $(BZLIB)randtable.c \
       $(ZLIB)gzlib.c $(ZLIB)gzread.c $(ZLIB)inflate.c $(ZLIB)inffast.c $(ZLIB)zutil.c $(ZLIB)inftrees.c $(ZLIB)crc32.c $(ZLIB)adler32.c \
       $(MAC)mach_gettime.c

OBJS  := $(SRCS:.c=.o)
DEBUG_OBJS := $(SRCS:.c=.debug-o)

DEPS  := $(SRCS:.c=.d)

EXE :=
ifeq ($(OS),Windows_NT)
	EXE := .exe
endif

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
	@echo Compiling debug $<
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
	@$(RM) $@ $(AFTER_RM)
	@ln $^ $@

archive: 
vczip.tar.gz: 

.PHONY: clean

LICENSE.non-commercial.txt: vczip$(EXE)
	@echo Generating $@
	@vczip$(EXE) -L > $@

clean:
	@echo Cleaning up
	@$(RM) $(OBJS) vczip$(EXE) vcpiz$(EXE) vccat$(EXE) $(AFTER_RM)

clean-debug:
	@echo Cleaning up debug
	@$(RM) $(DEBUG_OBJS) vczip-debug$(EXE) $(AFTER_RM)

clean-all: clean clean-debug
	@echo Cleaning up dependencies
	@$(RM) $(DEPS) $(AFTER_RM)
