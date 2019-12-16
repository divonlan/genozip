#/bin/bash --norc

# ------------------------------------------------------------------
#   make.sh
#   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
#   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

CC=gcc
CFLAGS=-Ibzlib -Izlib -D_LARGEFILE64_SOURCE=1 -Ofast
LIBS = -pthread -lm

DEPS = vczip.h
OBJ = vcf_header.o zip.o piz.o gloptimize.o buffer.o main.o vcffile.o squeeze.o zfile.o segregate.o profiler.o file.o vb.o\
      bzlib/blocksort.o bzlib/bzlib.o bzlib/compress.o bzlib/crctable.o bzlib/decompress.o bzlib/huffman.o bzlib/randtable.o \
      zlib/gzlib.o zlib/gzread.o zlib/inflate.o zlib/inffast.o zlib/zutil.o zlib/inftrees.o zlib/crc32.o zlib/adler32.o 

EXE :=
ifeq ($(OS),Windows_NT)
	EXE := .exe
endif

all: vczip$(EXE) vcpiz$(EXE) vccat$(EXE)

.c.o:
	@echo Compiling $<
	@$(CC) -c -o $@ $< $(CFLAGS)

all: vczip$(EXE) vcpiz$(EXE) vccat$(EXE) LICENSE.non-commercial.txt

vczip$(EXE): $(OBJ)
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

vcpiz$(EXE) vccat$(exe): vczip$(EXE)
	@echo Hard linking $@
ifeq ($(OS),Windows_NT)
	if exist $@ (del $@)
	mklink /h $@ $^ 
else
	@rm -f $@
	@ln $^ $@
endif

.PHONY: clean license

LICENSE.non-commercial.txt: vczip$(EXE)
	@echo Generating $@
	@vczip$(EXE) -L > $@

clean:
	@echo Cleaning up
	@rm -f $(OBJ) vczip$(EXE) vcpiz$(EXE) vccat$(EXE)

