#!/bin/bash
# Script to build MiddleMan with f2py for multiple Python versions using venv

set -e

FORTRAN_DIR="../fortran"
LIBS_DIR="../../OTSO/_core/libs"
MODULE_NAME="MiddleMan"

# Define Python versions to build for
PYTHON_VERSIONS="3.9 3.10 3.11 3.12 3.13 3.14"

# Function to build for a specific Python version
build_version() {
    local PYVER=$1
    echo "========================================"
    echo "Building for Python $PYVER"
    echo "========================================"
    
    # Use pyenv to install and manage Python versions
    if ! command -v pyenv >/dev/null 2>&1; then
        echo "pyenv is not installed. Please install pyenv and try again."
        return 1
    fi

    # Install Python version with pyenv if not already installed
    if ! pyenv versions --bare | grep -qx "$PYVER"; then
        echo "Installing Python $PYVER with pyenv..."
        pyenv install "$PYVER" || { echo "Failed to install Python $PYVER with pyenv."; return 1; }
    else
        echo "Python $PYVER already installed with pyenv."
    fi

    # Get the full path to the pyenv Python executable

    # Find the full installed version directory (e.g., 3.9.25 for 3.9)
    FULL_PYVER=$(pyenv versions --bare | grep -E "^$PYVER(\\.|$)" | head -n 1)
    if [ -z "$FULL_PYVER" ]; then
        echo "Could not find a pyenv version directory for $PYVER."
        pyenv versions
        return 1
    fi

    PYENV_BIN_DIR="$(pyenv root)/versions/$FULL_PYVER/bin"
    PYTHON_EXE="$PYENV_BIN_DIR/python$PYVER"
    if [ ! -x "$PYTHON_EXE" ]; then
        # Fallback for some pyenv layouts
        PYTHON_EXE="$PYENV_BIN_DIR/python3"
    fi
    if [ ! -x "$PYTHON_EXE" ]; then
        echo "Could not find Python executable for $PYVER (full: $FULL_PYVER) in pyenv."
        echo "Contents of $PYENV_BIN_DIR for debugging:"
        ls -l "$PYENV_BIN_DIR"
        return 1
    fi

    # Create virtual environment
    local VENV_NAME="middleman_py${PYVER}"
    local VENV_DIR="../venvs/$VENV_NAME"

    echo "Creating virtual environment $VENV_NAME..."
    mkdir -p "../venvs"

    if [ ! -d "$VENV_DIR" ]; then
        "$PYTHON_EXE" -m venv "$VENV_DIR"
        if [ $? -ne 0 ]; then
            echo "Failed to create virtual environment for Python $PYVER, skipping..."
            return 1
        fi
    else
        echo "Virtual environment $VENV_NAME already exists."
    fi
    
    # Activate virtual environment
    echo "Activating virtual environment..."
    source "$VENV_DIR/bin/activate" || source "$VENV_DIR/Scripts/activate"
    
    # Install required packages
    echo "Installing/updating required packages..."
    PYTHON_VERSION_MAJOR=$(python -c "import sys; print(sys.version_info.major)")
    PYTHON_VERSION_MINOR=$(python -c "import sys; print(sys.version_info.minor)")
    if [ "$PYTHON_VERSION_MAJOR" -lt 3 ] || { [ "$PYTHON_VERSION_MAJOR" -eq 3 ] && [ "$PYTHON_VERSION_MINOR" -lt 12 ]; }; then
        python -m pip install --upgrade pip
        python -m pip install --force-reinstall --no-cache-dir setuptools==68.0.0
    else
        python -m pip install --upgrade pip setuptools
    fi
    python -m pip install --upgrade numpy pandas ninja meson
    
    # Verify Python version
    python --version
    
    # Change to fortran directory for compilation
    cd "$FORTRAN_DIR"
    echo "Running f2py to build $MODULE_NAME for Python $PYVER..."
    
    # Read file list from fortran_files.txt like the bat file does
    if [ -f "fortran_files.txt" ]; then
        FILES=$(cat fortran_files.txt | tr '\n' ' ')
        echo "Building with files from fortran_files.txt: $FILES"
        python -m numpy.f2py -c $FILES -m "$MODULE_NAME" --backend=meson
    else
        echo "fortran_files.txt not found, using all Fortran files..."
        python -m numpy.f2py -c *.f95 *.f *.for -m "$MODULE_NAME" --backend=meson
    fi
    
    # Find the generated shared object file
    OUTPUT_FILE=$(ls ${MODULE_NAME}.*.so 2>/dev/null || ls ${MODULE_NAME}.*.pyd 2>/dev/null || echo "")
    
    if [ -z "$OUTPUT_FILE" ]; then
        echo "Error: Could not find generated $MODULE_NAME shared object file for Python $PYVER."
        cd - > /dev/null
        deactivate
        return 1
    fi
    
    cd - > /dev/null
    
    # Create libs directory if it doesn't exist
    mkdir -p "$LIBS_DIR"
    
    echo "Moving $FORTRAN_DIR/$OUTPUT_FILE to $LIBS_DIR/"
    mv "$FORTRAN_DIR/$OUTPUT_FILE" "$LIBS_DIR/"
    
    if [ $? -eq 0 ]; then
        echo "Successfully built and moved $MODULE_NAME for Python $PYVER"
    else
        echo "Error: Failed to move $OUTPUT_FILE for Python $PYVER"
        deactivate
        return 1
    fi
    
    # Deactivate virtual environment
    deactivate
    echo
}

# Check if a specific Python version was passed as parameter
if [ -z "$1" ]; then
    # No parameter - build for all versions (local use)
    echo "Starting multi-version build for $MODULE_NAME..."
    echo "Building for Python versions: $PYTHON_VERSIONS"
    echo
    
    # Loop through each Python version
    for PYVER in $PYTHON_VERSIONS; do
        build_version "$PYVER" || echo "Skipping Python $PYVER due to errors."
    done
else
    # Parameter provided - build for specific version (CI use)
    echo "Building $MODULE_NAME for Python $1 only..."
    build_version "$1"
fi

echo "========================================"
echo "Build complete!"
echo "========================================"
echo "Built libraries are in: $LIBS_DIR"
echo
echo "Listing all built libraries:"
ls -la "$LIBS_DIR"/${MODULE_NAME}.* 2>/dev/null || echo "No libraries found"

echo
echo "Done."