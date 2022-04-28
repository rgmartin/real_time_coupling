@echo off
SET THEFILE=E:\Documentos\IB\Maestría\REPO\PUMITA\PumitaTester_CAREM_nuevo\PumitaTester\reactor.dll
echo Linking %THEFILE%
C:\lazarus\fpc\3.0.4\bin\i386-win32\ld.exe -b pei-i386 -m i386pe  --gc-sections  -s --dll  --entry _DLLMainCRTStartup   --base-file base.$$$ -o "E:\Documentos\IB\Maestría\REPO\PUMITA\PumitaTester_CAREM_nuevo\PumitaTester\reactor.dll" "E:\Documentos\IB\Maestría\REPO\PUMITA\PumitaTester_CAREM_nuevo\PumitaTester\link.res"
if errorlevel 1 goto linkend
C:\lazarus\fpc\3.0.4\bin\i386-win32\dlltool.exe -S C:\lazarus\fpc\3.0.4\bin\i386-win32\as.exe -D "E:\Documentos\IB\Maestría\REPO\PUMITA\PumitaTester_CAREM_nuevo\PumitaTester\reactor.dll" -e exp.$$$ --base-file base.$$$ 
if errorlevel 1 goto linkend
C:\lazarus\fpc\3.0.4\bin\i386-win32\ld.exe -b pei-i386 -m i386pe  -s --dll  --entry _DLLMainCRTStartup   -o "E:\Documentos\IB\Maestría\REPO\PUMITA\PumitaTester_CAREM_nuevo\PumitaTester\reactor.dll" "E:\Documentos\IB\Maestría\REPO\PUMITA\PumitaTester_CAREM_nuevo\PumitaTester\link.res" exp.$$$
if errorlevel 1 goto linkend
C:\lazarus\fpc\3.0.4\bin\i386-win32\postw32.exe --subsystem console --input "E:\Documentos\IB\Maestría\REPO\PUMITA\PumitaTester_CAREM_nuevo\PumitaTester\reactor.dll" --stack 16777216
if errorlevel 1 goto linkend
goto end
:asmend
echo An error occurred while assembling %THEFILE%
goto end
:linkend
echo An error occurred while linking %THEFILE%
:end
