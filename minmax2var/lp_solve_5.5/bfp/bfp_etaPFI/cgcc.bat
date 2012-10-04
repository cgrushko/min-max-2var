@echo off

set c=gcc

REM determine platform (win32/win64)
echo main(){printf("SET PLATFORM=win%%d\n", (int) (sizeof(void *)*8));}>platform.c
%c% platform.c -o platform.exe
del platform.c
platform.exe >platform.bat
del platform.exe
call platform.bat
del platform.bat

if not exist bin\%PLATFORM%\*.* md bin\%PLATFORM%

set src=../../colamd/colamd.c ../../shared/commonlib.c bfp_etaPFI.c lp_etaPFI.c ../../lp_utils.c

%c% -DINLINE=static -I.. -I../.. -I../../colamd -I../../shared -s -O3 -shared -mno-cygwin -enable-stdcall-fixup -D_WINDLL -D_USRDLL -DWIN32 -DRoleIsExternalInvEngine -DINVERSE_ACTIVE=INVERSE_ETAPFI %src% ../lp_BFP.def -o bin\%PLATFORM%\bfp_etaPFI.dll

%c% -DINLINE=static -I.. -I../.. -I../../colamd -I../../shared -s -O3 -shared -D_WINDLL -D_USRDLL -DWIN32 -DRoleIsExternalInvEngine -DINVERSE_ACTIVE=INVERSE_ETAPFI %src% -o bin\%PLATFORM%\libbfp_etaPFI.so

%c% -DINLINE=static -I.. -I../.. -I../../colamd -I../../shared -s -O3 -c -DRoleIsExternalInvEngine -DINVERSE_ACTIVE=INVERSE_ETAPFI %src%
ar rv bin\%PLATFORM%\libbfp_etaPFI.a *.o

if exist *.o del *.o

set PLATFORM=
