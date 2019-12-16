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

all: vczip vcpiz vccat

.c.o:
	@echo Compiling $<
	@$(CC) -c -o $@ $< $(CFLAGS)

all: vczip vcpiz vccat LICENSE.non-commercial.txt

vczip: $(OBJ)
	@echo Linking $@
	@$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

vcpiz vccat: vczip
	@echo Hard linking $@
	@rm -f $@
	@ln $^ $@

.PHONY: clean license

LICENSE.non-commercial.txt: vczip
	@echo Generating $@
	@vczip -L > $@

clean:
	@echo Cleaning up
	@rm -f $(OBJ) vczip vcpiz vccat

