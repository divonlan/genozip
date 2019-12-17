gcc.exe -Ibzlib -Izlib -D_LARGEFILE64_SOURCE=1 -Wall -Ofast -lpthread -lm vcf_header.c zip.c piz.c gloptimize.c buffer.c main.c vcffile.c squeeze.c zfile.c segregate.c profiler.c file.c vb.c bzlib\blocksort.c bzlib\bzlib.c bzlib\compress.c bzlib\crctable.c bzlib\decompress.c bzlib\huffman.c bzlib\randtable.c zlib\gzlib.c zlib\gzread.c zlib\inflate.c zlib\inffast.c zlib\zutil.c zlib\inftrees.c zlib\crc32.c zlib\adler32.c mac\mach_gettime.c  -o vczip.exe

md %PREFIX%\bin >nul 2>&1
copy %SRC_DIR%\vczip.exe %PREFIX%\bin\vczip.exe
copy %SRC_DIR%\vczip.exe %PREFIX%\bin\vcpiz.exe
copy %SRC_DIR%\vczip.exe %PREFIX%\bin\vccat.exe
copy %SRC_DIR%\LICENSE.non-commerical.txt %PREFIX%
copy %SRC_DIR%\LICENSE.commerical.txt %PREFIX%

