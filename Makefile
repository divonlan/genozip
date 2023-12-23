ifeq ($(OS),Windows_NT)
EXE = .exe
endif

all debug opt clean clean-debug clean-optimized git-pull genozip$(EXE) genounzip$(EXE) genocat$(EXE) genols$(EXE) genozip-prod genozip-prod.exe install distribution-maintenance distribution:
	@(cd src ; $(MAKE) $(MAKECMDGOALS)) # note: cd because make -C doesn't work well on Mac

.PHONY: all debug opt clean clean-debug clean-optimized git-pull genozip$(EXE) genounzip$(EXE) genocat$(EXE) genols$(EXE) genozip-prod genozip-prod.exe install
