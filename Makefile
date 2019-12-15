#/bin/bash --norc

# ------------------------------------------------------------------
#   make.sh
#   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
#   Please see terms and conditions in the file LICENSE.txt

CC=gcc
CFLAGS=-Ibzlib -Izlib -D_LARGEFILE64_SOURCE=1 -Ofast
LIBS = -pthread -lm

DEPS = vczip.h
OBJ = vcf_header.o zip.o piz.o gloptimize.c buffer.o main.o vcffile.o squeeze.o zfile.o segregate.o profiler.o file.o vb.o\
      bzlib/blocksort.o bzlib/bzlib.o bzlib/compress.o bzlib/crctable.o bzlib/decompress.o bzlib/huffman.o bzlib/randtable.o \
      zlib/gzlib.o zlib/gzread.o zlib/inflate.o zlib/inffast.o zlib/zutil.o zlib/inftrees.o zlib/crc32.o zlib/adler32.o 

ODIR=obj

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

vczip: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

vcpiz: vczip
	ln $^ $o

vccat: vczip
	ln $^ $o

.PHONY: clean

all: vczip vcpiz vccat

clean:
	rm -f *.o zlib/*.o bzlib/*.o *~ core $(INCDIR)/*~ vczip


