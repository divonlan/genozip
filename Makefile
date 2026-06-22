ifeq ($(OS),Windows_NT)
EXE = .exe
endif

all debug opt clean clean-debug clean-optimized genozip$(EXE) genounzip$(EXE) genocat$(EXE) genols$(EXE) genozip-latest genozip-latest.exe install distribution-maintenance distribution:
	@(cd src ; $(MAKE) $(MAKECMDGOALS)) # note: cd because make -C doesn't work well on Mac

LICENSE.txt: 
	@(cd src ; $(MAKE) ../LICENSE.txt)

.PHONY: all debug opt clean clean-debug clean-optimized genozip$(EXE) genounzip$(EXE) genocat$(EXE) genols$(EXE) genozip-latest genozip-latest.exe install LICENSE.txt
