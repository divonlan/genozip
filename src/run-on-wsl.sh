# Used to run stuff via wsl from Makefile on Windows as it was too complicated to correctly escape this command within Makefile
MSYS2_ARG_CONV_EXCL="*" wsl "PATH=\$PATH:/home/divon/miniconda3/bin/" $*
