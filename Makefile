ifeq ($(OS),Windows_NT)
EXE = .exe
endif

all debug opt clean clean-debug clean-optimized git-pull genozip$(EXE) genounzip$(EXE) genocat$(EXE) genols$(EXE) genozip-prod genozip-prod.exe :
	$(MAKE) -C src $(MAKECMDGOALS)

.PHONY: all debug opt clean clean-debug clean-optimized git-pull genozip$(EXE) genounzip$(EXE) genocat$(EXE) genols$(EXE) genozip-prod genozip-prod.exe
