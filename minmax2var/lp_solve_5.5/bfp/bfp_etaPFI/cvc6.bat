@echo off

set c=cl

REM determine platform (win32/win64)
echo main(){printf("SET PLATFORM=win%%d\n", (int) (sizeof(void *)*8));}>platform.c
%c% /nologo platform.c /Feplatform.exe
del platform.c
platform.exe >platform.bat
del platform.exe
call platform.bat
del platform.bat

if "%PLATFORM%" == "win32" goto ok1
echo.
echo This batch file is intended for 32 bit compilation with MS Visual C 6
echo For newer versions use cvc8*.bat
goto done
:ok1

if not exist bin\%PLATFORM%\*.* md bin\%PLATFORM%

set src=../../colamd/colamd.c ../../shared/commonlib.c bfp_etaPFI.c lp_etaPFI.c ../../lp_utils.c

%c% -I.. -I../.. -I../../colamd -I../../shared /LD /MD /O2 /Zp8 /Gz -D_WINDLL -D_USRDLL -DWIN32 -DRoleIsExternalInvEngine -DINVERSE_ACTIVE=INVERSE_ETAPFI %src% ../lp_BFP.def -o bin\%PLATFORM%\bfp_etaPFI.dll

if exist *.obj del *.obj
:done
set PLATFORM=
