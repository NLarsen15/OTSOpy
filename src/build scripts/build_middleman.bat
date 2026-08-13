@echo off
setlocal enabledelayedexpansion
set FORTRAN_DIR=..\fortran
set LIBS_DIR=..\..\OTSO\_core\libs
set MODULE_NAME=MiddleMan

rem Check if a specific Python version was passed as parameter
if "%1"=="" (
    rem No parameter - build for all versions (local use)
    set PYTHON_VERSIONS=3.9 3.10 3.11 3.12 3.13 3.14
    echo Starting multi-version build for %MODULE_NAME%...
    echo Building for Python versions: !PYTHON_VERSIONS!
    echo.
    rem Loop through each Python version
    for %%v in (!PYTHON_VERSIONS!) do call :build_version %%v
) else (
    rem Parameter provided - build for specific version (GitHub Actions use)
    echo Building %MODULE_NAME% for Python %1 only...
    call :build_version %1
)

echo ========================================
echo Build complete!
echo ========================================
echo Built libraries are in: %LIBS_DIR%
echo.
echo Listing all built libraries:
dir /b "%LIBS_DIR%\%MODULE_NAME%.*" 2>nul

echo.
echo Done.
endlocal
goto :eof

:build_version
set PYVER=%1
echo ========================================
echo Building for Python %PYVER%
echo ========================================

rem Create or update conda environment
echo Checking if conda environment middleman_py%PYVER% exists...
call conda info --envs | findstr "middleman_py%PYVER%" >nul
if errorlevel 1 (
    echo Creating new conda environment for Python %PYVER%...
    call conda create -n middleman_py%PYVER% python=%PYVER% numpy pandas ninja meson -y
    if errorlevel 1 goto :skip_version
) else (
    echo Environment middleman_py%PYVER% already exists, ensuring packages are installed...
    call conda install -n middleman_py%PYVER% numpy pandas ninja meson -y
)

rem Activate the environment
echo Activating environment middleman_py%PYVER%...
call conda activate middleman_py%PYVER%
if errorlevel 1 goto :skip_version

rem Verify Python version
python -c "import sys; print('Using Python', sys.version)"

rem Change to fortran directory for compilation
cd /d "%FORTRAN_DIR%"
echo Running f2py to build %MODULE_NAME% for Python %PYVER%...

rem Build the file list from fortran_files.txt
set FILES=
for /f "delims=" %%i in (fortran_files.txt) do set FILES=!FILES! %%i

rem Clean any existing build files
del %MODULE_NAME%.*.pyd 2>nul
del %MODULE_NAME%.*.so 2>nul

python -m numpy.f2py -c !FILES! -m %MODULE_NAME% --backend=meson

rem Find the generated .pyd or .so file
set OUTPUT_FILE=
for %%f in (%MODULE_NAME%.*.pyd) do set OUTPUT_FILE=%%f
if not defined OUTPUT_FILE for %%f in (%MODULE_NAME%.*.so) do set OUTPUT_FILE=%%f

if not defined OUTPUT_FILE goto :skip_version

rem Create libs directory if it doesn't exist
cd /d "%~dp0"
if not exist "%LIBS_DIR%" mkdir "%LIBS_DIR%"

echo Moving "%FORTRAN_DIR%\%OUTPUT_FILE%" to "%LIBS_DIR%"...
move "%FORTRAN_DIR%\%OUTPUT_FILE%" "%LIBS_DIR%" >nul
if errorlevel 1 goto :skip_version

echo Successfully built and moved %MODULE_NAME% for Python %PYVER%
call conda deactivate
echo.
goto :eof

:skip_version
echo Skipping Python %PYVER% due to errors.
call conda deactivate 2>nul
echo.
goto :eof