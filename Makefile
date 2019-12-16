# ------------------------------------------------------------------
#   Makefile
#   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
#   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

CC=gcc
CFLAGS       = -Ibzlib -Izlib -D_LARGEFILE64_SOURCE=1 -Ofast
CFLAGS_DEBUG = -Ibzlib -Izlib -D_LARGEFILE64_SOURCE=1 -DDEBUG -g
LIBS = -pthread -lm

SLASH :=
ifeq ($(OS),Windows_NT)
	SLASH := \\
else
	SLASH := /
endif

ZLIB = zlib$(SLASH)
BZLIB = bzlib$(SLASH)

DEPS = vczip.h
SRC = vcf_header.c zip.c piz.c gloptimize.c buffer.c main.c vcffile.c squeeze.c zfile.c segregate.c profiler.c file.c vb.c \
      $(BZLIB)blocksort.c $(BZLIB)bzlib.c $(BZLIB)compress.c $(BZLIB)crctable.c $(BZLIB)decompress.c $(BZLIB)huffman.c $(BZLIB)randtable.c \
      $(ZLIB)gzlib.c $(ZLIB)gzread.c $(ZLIB)inflate.c $(ZLIB)inffast.c $(ZLIB)zutil.c $(ZLIB)inftrees.c $(ZLIB)crc32.c $(ZLIB)adler32.c 

OBJ  := $(SRC:.c=.o)

DEBUG_OBJ := $(SRC:.c=.debug-o)

EXE :=
ifeq ($(OS),Windows_NT)
	EXE := .exe
endif

all: vczip$(EXE) vcpiz$(EXE) vccat$(EXE)

debug: vczip-debug$(EXE)

%.o: %.c vczip.h 
	@echo Compiling $<
	@$(CC) -c -o $@ $< $(CFLAGS)

%.debug-o: %.c vczip.h 
	@echo Compiling debug $<
	@$(CC) -c -o $@ $< $(CFLAGS_DEBUG)

all: vczip$(EXE) vcpiz$(EXE) vccat$(EXE)

vczip$(EXE): $(OBJ)
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

vczip-debug$(EXE): $(DEBUG_OBJ)
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

vcpiz$(EXE) vccat$(EXE): vczip$(EXE)
	@echo Hard linking $@
ifeq ($(OS),Windows_NT)
	@if exist $@ (del $@)
else
	@rm -f $@
endif
	@ln $^ $@

.PHONY: clean license

LICENSE.non-commercial.txt: vczip$(EXE)
	@echo Generating $@
	@vczip$(EXE) -L > $@

clean:
	@echo Cleaning up
ifeq ($(OS),Windows_NT)
	@del /q /f $(OBJ) vczip$(EXE) vcpiz$(EXE) vccat$(EXE) >nul 2>&1
else
	@rm -f $(OBJ) vczip$(EXE) vcpiz$(EXE) vccat$(EXE)
endif

clean-debug:
	@echo Cleaning up debug
ifeq ($(OS),Windows_NT)
	@del /q /f $(DEBUG_OBJ) vczip-debug$(EXE) >nul 2>&1
else
	@rm -f $(DEBUG_OBJ) vczip-debug$(EXE)
endif
