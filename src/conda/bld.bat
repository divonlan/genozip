make OS=Windows_NT all

if not exist "%LIBRARY_BIN%" mkdir "%LIBRARY_BIN%"

copy src\genozip.exe   "%LIBRARY_BIN%"
copy src\genounzip.exe "%LIBRARY_BIN%"
copy src\genols.exe    "%LIBRARY_BIN%"
copy src\genocat.exe   "%LIBRARY_BIN%"

if errorlevel 1 exit /b 1
exit /b 0
