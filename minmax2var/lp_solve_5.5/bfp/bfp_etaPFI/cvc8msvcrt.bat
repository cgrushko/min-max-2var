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

if not exist bin\%PLATFORM%\*.* md bin\%PLATFORM%

set src=../../colamd/colamd.c ../../shared/commonlib.c bfp_etaPFI.c lp_etaPFI.c ../../lp_utils.c

if not exist bin\%PLATFORM%\*.* md bin\%PLATFORM%

%c% -I.. -I../.. -I../../colamd -I../../shared /LD /MD /O2 /Zp8 /Gz -D_WINDLL -D_USRDLL -DWIN32 -DRoleIsExternalInvEngine -DINVERSE_ACTIVE=INVERSE_ETAPFI -D_CRT_SECURE_NO_DEPRECATE %src% /Febin\%PLATFORM%\bfp_etaPFI.dll
editbin /LARGEADDRESSAWARE bin\%PLATFORM%\bfp_etaPFI.dll

rem http://msdn2.microsoft.com/en-us/library/ms235229.aspx
rem for vs2005 need to embed manifest in dll with manifest tool -  #2 on the next line does this.
mt /outputresource:"bin\%PLATFORM%\bfp_etaPFI.dll;#2" /manifest "bin\%PLATFORM%\bfp_etaPFI.dll.manifest"

if exist *.obj del *.obj
set PLATFORM=
