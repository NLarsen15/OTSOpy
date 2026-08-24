```bat
@echo off
setlocal EnableDelayedExpansion

rem ============================================================
rem CONFIGURATION
rem ============================================================

set "FORTRAN_DIR=..\fortran"
set "LIBS_DIR=..\..\OTSO\_core\libs"
set "MODULE_NAME=MiddleMan"
set "EXTENSION_NAME=_MiddleMan"

rem Python versions
set "PYTHON_VERSIONS=3.9 3.10 3.11 3.12 3.13 3.14"

echo.
echo ============================================================
echo        MiddleMan Multi-Version Build
echo ============================================================
echo.

rem ============================================================
rem BUILD ALL OR SPECIFIC PYTHON VERSION
rem ============================================================

if "%~1"=="" (

    echo No Python version specified.
    echo Building all versions:
    echo %PYTHON_VERSIONS%
    echo.

    for %%V in (%PYTHON_VERSIONS%) do (
        call :build_version %%V
    )

) else (

    echo Building Python %~1 only...
    echo.

    call :build_version %~1
)

echo.
echo ============================================================
echo        BUILD PROCESS COMPLETE
echo ============================================================
echo.

echo Built libraries are in:
echo %LIBS_DIR%
echo.

if exist "%LIBS_DIR%" (
    echo Listing built libraries:
    dir /b "%LIBS_DIR%\%MODULE_NAME%.*" 2>nul
) else (
    echo No library directory exists.
)

echo.
echo Done.

endlocal
goto :eof


rem ============================================================
rem BUILD ONE PYTHON VERSION
rem ============================================================

:build_version

set "PYVER=%~1"
set "ENV_NAME=middleman_py%PYVER%"

echo.
echo.
echo ============================================================
echo        BUILDING PYTHON %PYVER%
echo ============================================================
echo.


rem ============================================================
rem CHECK CONDA
rem ============================================================

echo Checking conda...

where conda

if errorlevel 1 (
    echo ERROR: conda was not found in PATH.
    goto :skip_version
)

echo.


rem ============================================================
rem CREATE OR UPDATE CONDA ENVIRONMENT
rem ============================================================

echo Checking if conda environment %ENV_NAME% exists...

call conda env list | findstr /R /C:"%ENV_NAME% " >nul 2>&1

if errorlevel 1 (

    echo Environment %ENV_NAME% does not exist.
    echo Creating environment...

    call conda create -n "%ENV_NAME%" ^
        python=%PYVER% ^
        numpy ^
        pandas ^
        ninja ^
        meson ^
        pip ^
        -y

    if errorlevel 1 (
        echo ERROR: Failed to create conda environment.
        goto :skip_version
    )

) else (

    echo Environment %ENV_NAME% already exists.
    echo Updating packages...

    call conda install -n "%ENV_NAME%" ^
        numpy ^
        pandas ^
        ninja ^
        meson ^
        pip ^
        -y

    if errorlevel 1 (
        echo ERROR: Failed to update conda environment.
        goto :skip_version
    )
)

echo.


rem ============================================================
rem ACTIVATE CONDA ENVIRONMENT
rem ============================================================

echo Activating %ENV_NAME%...

call conda activate "%ENV_NAME%"

if errorlevel 1 (
    echo ERROR: Failed to activate %ENV_NAME%.
    goto :skip_version
)

echo Environment activated.
echo.


rem ============================================================
rem INSTALL F90WRAP
rem ============================================================

echo Installing/checking f90wrap...

python -m pip install --upgrade f90wrap

if errorlevel 1 (
    echo ERROR: Failed to install f90wrap.
    goto :skip_version
)

echo.
echo f90wrap installed.
echo.

where f90wrap

if errorlevel 1 (
    echo ERROR: f90wrap executable cannot be found.
    goto :skip_version
)

echo.


rem ============================================================
rem PYTHON DEBUG INFORMATION
rem ============================================================

echo ============================================================
echo DEBUG: PYTHON INFORMATION
echo ============================================================
echo.

echo Python executable:
where python

echo.
echo Python version:
python --version

echo.
echo Python details:
python -c "import sys; import sysconfig; import platform; print('Executable:',sys.executable); print('Version:',sys.version); print('Platform:',platform.platform()); print('Machine:',platform.machine()); print('EXT_SUFFIX:',sysconfig.get_config_var('EXT_SUFFIX'))"

echo.
echo ============================================================
echo.


rem ============================================================
rem COMPILER INFORMATION
rem ============================================================

echo ============================================================
echo DEBUG: COMPILER INFORMATION
echo ============================================================
echo.

echo GCC:
where gcc

if errorlevel 1 (
    echo WARNING: gcc not found.
) else (
    gcc --version
)

echo.
echo GFortran:
where gfortran

if errorlevel 1 (
    echo WARNING: gfortran not found.
) else (
    gfortran --version
)

echo.
echo ============================================================
echo.


rem ============================================================
rem CHANGE TO FORTRAN DIRECTORY
rem ============================================================

echo Changing to Fortran directory:

cd /d "%~dp0%FORTRAN_DIR%"

if errorlevel 1 (
    echo ERROR: Cannot change to:
    echo %~dp0%FORTRAN_DIR%
    goto :skip_version
)

echo.
echo Current directory:
cd
echo.


rem ============================================================
rem CHECK fortran_files.txt
rem ============================================================

if not exist "fortran_files.txt" (

    echo.
    echo ============================================================
    echo ERROR: fortran_files.txt NOT FOUND
    echo ============================================================
    echo.

    echo Current directory:
    cd

    echo.
    echo Directory contents:
    dir /a

    goto :skip_version
)


rem ============================================================
rem BUILD FILE LIST
rem ============================================================

set "FILES="

for /f "usebackq delims=" %%F in ("fortran_files.txt") do (
    if not "%%F"=="" (
        set "FILES=!FILES! %%F"
    )
)


rem ============================================================
rem DEBUG: FILE LIST
rem ============================================================

echo.
echo ============================================================
echo DEBUG: FORTRAN SOURCE FILE LIST
echo ============================================================
echo.

echo Contents of fortran_files.txt:
echo ------------------------------------------------------------
type "fortran_files.txt"
echo ------------------------------------------------------------

echo.
echo FILES variable:
echo ------------------------------------------------------------
echo !FILES!
echo ------------------------------------------------------------

echo.


rem ============================================================
rem DEBUG: VERIFY SOURCE FILES
rem ============================================================

echo ============================================================
echo DEBUG: VERIFYING SOURCE FILES
echo ============================================================
echo.

for %%F in (!FILES!) do (

    if exist "%%F" (
        echo FOUND:
        echo   %%~fF
    ) else (
        echo MISSING:
        echo   %%F
    )

)

echo.
echo ============================================================
echo.


rem ============================================================
rem DEBUG: DIRECTORY BEFORE F90WRAP
rem ============================================================

echo.
echo ============================================================
echo DEBUG: BEFORE F90WRAP
echo ============================================================
echo.

echo Current directory:
cd

echo.
echo ------------------------------------------------------------
echo ALL DIRECTORY CONTENTS
echo ------------------------------------------------------------
dir /a
echo ------------------------------------------------------------

echo.
echo ------------------------------------------------------------
echo ALL FILE/DIRECTORY NAMES
echo ------------------------------------------------------------
dir /a /b
echo ------------------------------------------------------------

echo.
echo ------------------------------------------------------------
echo FORTRAN FILES
echo ------------------------------------------------------------
dir /a *.f *.f90 *.f95 *.for *.f77 *.F *.F90 *.F95 2>nul
echo ------------------------------------------------------------

echo.
echo ------------------------------------------------------------
echo fortran_files.txt
echo ------------------------------------------------------------
if exist "fortran_files.txt" (
    type "fortran_files.txt"
) else (
    echo fortran_files.txt DOES NOT EXIST
)
echo ------------------------------------------------------------

echo.
echo ------------------------------------------------------------
echo FILES VARIABLE
echo ------------------------------------------------------------
echo !FILES!
echo ------------------------------------------------------------

echo.
echo ------------------------------------------------------------
echo F90WRAP EXECUTABLE
echo ------------------------------------------------------------
where f90wrap
echo ------------------------------------------------------------

echo.
echo ------------------------------------------------------------
echo F90WRAP VERSION
echo ------------------------------------------------------------
f90wrap --version 2>&1

if errorlevel 1 (
    echo f90wrap --version returned an error.
    echo Trying f90wrap --help:
    f90wrap --help 2>&1
)

echo ------------------------------------------------------------

echo.
echo ============================================================
echo END DEBUG: BEFORE F90WRAP
echo ============================================================
echo.


rem ============================================================
rem CLEAN OLD BUILD FILES
rem ============================================================

echo Cleaning previous generated files...

del /q "%MODULE_NAME%.*.pyd" 2>nul
del /q "%MODULE_NAME%.*.so" 2>nul
del /q "%EXTENSION_NAME%.*.pyd" 2>nul
del /q "%EXTENSION_NAME%.*.so" 2>nul

if exist "f90wrap" (
    echo Removing previous f90wrap directory...
    rmdir /s /q "f90wrap"
)

echo Cleanup complete.
echo.


rem ============================================================
rem RUN F90WRAP
rem ============================================================

echo.
echo ============================================================
echo DEBUG: RUNNING F90WRAP
echo ============================================================
echo.

echo Expected source:
echo MiddleMan.f95
echo.

if exist "MiddleMan.f95" (
    echo FOUND: MiddleMan.f95
) else (
    echo ERROR: MiddleMan.f95 NOT FOUND
    goto :skip_version
)

echo.
echo Exact command:
echo ------------------------------------------------------------
echo f90wrap -m "%MODULE_NAME%" MiddleMan.f95
echo ------------------------------------------------------------
echo.

echo Starting f90wrap...
echo.

f90wrap -m "%MODULE_NAME%" MiddleMan.f95

set "F90WRAP_ERROR=%ERRORLEVEL%"

echo.
echo f90wrap exit code:
echo %F90WRAP_ERROR%
echo.


rem ============================================================
rem DEBUG: DIRECTORY AFTER F90WRAP
rem ============================================================

echo.
echo ============================================================
echo DEBUG: AFTER F90WRAP
echo ============================================================
echo.

echo Current directory:
cd

echo.
echo ------------------------------------------------------------
echo ALL DIRECTORY CONTENTS
echo ------------------------------------------------------------
dir /a
echo ------------------------------------------------------------

echo.
echo ------------------------------------------------------------
echo ALL FILE/DIRECTORY NAMES
echo ------------------------------------------------------------
dir /a /b
echo ------------------------------------------------------------

echo.
echo ------------------------------------------------------------
echo RECURSIVE FILE LIST
echo ------------------------------------------------------------
dir /s /a /b
echo ------------------------------------------------------------

echo.
echo ------------------------------------------------------------
echo ALL FORTRAN FILES
echo ------------------------------------------------------------
dir /s /a /b *.f *.f90 *.f95 *.for *.f77 *.F *.F90 *.F95 2>nul
echo ------------------------------------------------------------

echo.
echo ------------------------------------------------------------
echo FILES CONTAINING "WRAP"
echo ------------------------------------------------------------
dir /s /a /b *wrap* 2>nul
echo ------------------------------------------------------------

echo.
echo ------------------------------------------------------------
echo FILES CONTAINING "%MODULE_NAME%"
echo ------------------------------------------------------------
dir /s /a /b *%MODULE_NAME%* 2>nul
echo ------------------------------------------------------------


rem ============================================================
rem INSPECT F90WRAP DIRECTORY
rem ============================================================

echo.
echo ============================================================
echo DEBUG: F90WRAP DIRECTORY
echo ============================================================
echo.

if exist "f90wrap" (

    echo f90wrap directory EXISTS.

    echo.
    echo Contents:
    dir /a "f90wrap"

    echo.
    echo Recursive contents:
    dir /s /a /b "f90wrap"

    echo.
    echo F90 files:
    dir /s /a /b "f90wrap\*.f90" 2>nul

    echo.
    echo F95 files:
    dir /s /a /b "f90wrap\*.f95" 2>nul

    echo.
    echo Python files:
    dir /s /a /b "f90wrap\*.py" 2>nul

) else (

    echo f90wrap directory DOES NOT EXIST.

)

echo.
echo ============================================================
echo END DEBUG: AFTER F90WRAP
echo ============================================================
echo.


rem ============================================================
rem CHECK F90WRAP ERROR
rem ============================================================

if not "%F90WRAP_ERROR%"=="0" (

    echo.
    echo ============================================================
    echo F90WRAP FAILED
    echo ============================================================
    echo.

    echo Exit code:
    echo %F90WRAP_ERROR%

    echo.
    echo No f2py build will be attempted.

    goto :skip_version
)


rem ============================================================
rem FIND WRAPPER FILES
rem ============================================================

echo.
echo ============================================================
echo DEBUG: WRAPPER DISCOVERY
echo ============================================================
echo.

set "WRAPPERS="

for %%F in (f90wrap_*.f90) do (
    if exist "%%F" (
        echo FOUND WRAPPER:
        echo   %%~fF
        set "WRAPPERS=!WRAPPERS! "%%F""
    )
)

echo.
echo Wrapper variable:
echo ------------------------------------------------------------
echo !WRAPPERS!
echo ------------------------------------------------------------
echo.


if not defined WRAPPERS (

    echo.
    echo ============================================================
    echo ERROR: NO F90WRAP .F90 FILES FOUND
    echo ============================================================
    echo.

    echo Looking for ANY generated files:
    dir /s /a /b f90wrap\* 2>nul

    goto :skip_version
)


rem ============================================================
rem RUN F2PY
rem ============================================================

echo.
echo ============================================================
echo RUNNING F2PY
echo ============================================================
echo.

set "F90FLAGS=-O3 -fopenmp -funroll-loops -fautomatic -frecursive"
set "F77FLAGS=-O3 -fopenmp -funroll-loops -fautomatic -frecursive"
set "OPTFLAGS=-O3 -fopenmp -funroll-loops -fautomatic -frecursive"

echo F90 flags:
echo %F90FLAGS%

echo.
echo F77 flags:
echo %F77FLAGS%

echo.
echo OPT flags:
echo %OPTFLAGS%

echo.
echo Wrapper files:
echo !WRAPPERS!

echo.
echo Source files:
echo !FILES!

echo.

python -m numpy.f2py -c ^
    --backend=meson ^
    --f90flags="%F90FLAGS%" ^
    --f77flags="%F77FLAGS%" ^
    --opt="%OPTFLAGS%" ^
    !FILES! ^
    !WRAPPERS! ^
    -m "%EXTENSION_NAME%" ^
    -lgomp

set "F2PY_ERROR=%ERRORLEVEL%"

echo.
echo f2py exit code:
echo %F2PY_ERROR%
echo.


if not "%F2PY_ERROR%"=="0" (

    echo.
    echo ============================================================
    echo F2PY FAILED
    echo ============================================================
    echo.

    echo Directory contents:
    dir /a /b

    goto :skip_version
)


rem ============================================================
rem FIND GENERATED Pyd
rem ============================================================

echo.
echo ============================================================
echo SEARCHING FOR GENERATED EXTENSION
echo ============================================================
echo.

set "OUTPUT_FILE="

echo Looking for:
echo %EXTENSION_NAME%.*.pyd
echo.

for %%F in ("%EXTENSION_NAME%.*.pyd") do (
    if exist "%%~F" (
        set "OUTPUT_FILE=%%~nxF"
    )
)

if not defined OUTPUT_FILE (
    echo No %EXTENSION_NAME% .pyd found.
    echo.
    echo Looking for:
    echo %MODULE_NAME%.*.pyd

    for %%F in ("%MODULE_NAME%.*.pyd") do (
        if exist "%%~F" (
            set "OUTPUT_FILE=%%~nxF"
        )
    )
)

echo.

if not defined OUTPUT_FILE (

    echo ============================================================
    echo ERROR: NO WINDOWS Pyd FILE FOUND
    echo ============================================================
    echo.

    echo All files generated:
    dir /s /a /b

    goto :skip_version
)

echo FOUND:
echo !OUTPUT_FILE!
echo.


rem ============================================================
rem CREATE LIBRARY DIRECTORY
rem ============================================================

cd /d "%~dp0"

if not exist "%LIBS_DIR%" (
    mkdir "%LIBS_DIR%"
)

if not exist "%LIBS_DIR%" (
    echo ERROR: Could not create:
    echo %LIBS_DIR%
    goto :skip_version
)


rem ============================================================
rem COPY LIBRARY
rem ============================================================

echo.
echo Copying:
echo %FORTRAN_DIR%\!OUTPUT_FILE!

echo.
echo To:
echo %LIBS_DIR%

copy /Y "%FORTRAN_DIR%\!OUTPUT_FILE!" "%LIBS_DIR%\" >nul

if errorlevel 1 (
    echo ERROR: Failed to copy library.
    goto :skip_version
)

echo.
echo ============================================================
echo SUCCESS
echo ============================================================
echo Python:
echo %PYVER%
echo.
echo Extension:
echo !OUTPUT_FILE!
echo.
echo Destination:
echo %LIBS_DIR%
echo ============================================================
echo.

call conda deactivate >nul 2>&1

goto :eof


rem ============================================================
rem SKIP VERSION
rem ============================================================

:skip_version

echo.
echo ============================================================
echo SKIPPING PYTHON %PYVER% DUE TO ERRORS
echo ============================================================
echo.

call conda deactivate >nul 2>&1

goto :eof
```
