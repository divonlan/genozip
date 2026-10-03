# ------------------------------------------------------------------
#   Makefile
#   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
#   Please see terms and conditions in the file LICENSE.txt
#
#   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
#   and subject to penalties specified in the license.

ifeq ($(OS),Windows_NT)
EXE = .exe
endif

all debug opt genozip$(EXE) genozip-debug$(EXE) genounzip$(EXE) genocat$(EXE) genols$(EXE) genozip-latest genozip-latest.exe install distribution-maintenance distribution clean clean-debug clean-opt clean-dev clean-test clean-installers clean-dumps :
	@cd src ; $(MAKE) --no-print-directory $@ # note: cd because make -C doesn't work well on Mac

LICENSE.txt: 
	@cd src ; $(MAKE) --no-print-directory ../LICENSE.txt 

.PHONY: all debug opt clean clean-debug clean-opt clean-dev clean-test clean-installers clean-dumps genozip$(EXE) genozip-debug$(EXE) genounzip$(EXE) genocat$(EXE) genols$(EXE) genozip-latest genozip-latest.exe install distribution-maintenance distribution LICENSE.txt
