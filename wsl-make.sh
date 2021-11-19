# this script is run as the "Build Task" menu in Visual Studio Code, for building Linux, as defined in .vscode/tasks.json
make -j -C ~/genozip $1
