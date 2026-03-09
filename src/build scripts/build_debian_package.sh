#!/bin/bash
# Script to build OTSO Debian package

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DEBIAN_DIR="$PROJECT_ROOT/debian"
BUILD_DIR="$PROJECT_ROOT/build"
DIST_DIR="$PROJECT_ROOT/dist"

# Package information (will be read from setup.py)
PACKAGE_NAME="python3-otso"
PACKAGE_VERSION=""
PACKAGE_DESCRIPTION=""
PACKAGE_AUTHOR=""
PACKAGE_EMAIL=""
PACKAGE_URL=""

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}       OTSO Debian Package Builder      ${NC}"
echo -e "${BLUE}========================================${NC}"

# Function to print colored output
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to extract information from setup.py
extract_package_info() {
    log_info "Extracting package information from setup.py..."
    
    cd "$PROJECT_ROOT"
    
    # Extract version
    PACKAGE_VERSION=$(python3 -c "
import sys
sys.path.insert(0, '.')
from setup import setup
import setuptools
# Monkey patch setup to capture arguments
original_setup = setuptools.setup
captured_args = {}
def capture_setup(**kwargs):
    captured_args.update(kwargs)
    return original_setup(**kwargs)
setuptools.setup = capture_setup
# Import setup.py to capture arguments
exec(open('setup.py').read())
print(captured_args.get('version', '1.1.2'))
" 2>/dev/null || echo "1.1.2")
    
    # Extract other info using more reliable method
    PACKAGE_DESCRIPTION=$(python3 -c "
exec(open('setup.py').read())
" 2>/dev/null | grep -oP "description='\K[^']+" || echo "Geomagnetic Cutoff Computation Tool")
    
    # Set defaults if extraction failed
    PACKAGE_AUTHOR="Nicholas Larsen"
    PACKAGE_EMAIL="nlarsen1505@gmail.com"
    PACKAGE_URL="https://github.com/NLarsen15/OTSOpy"
    
    log_info "Package: $PACKAGE_NAME"
    log_info "Version: $PACKAGE_VERSION"
    log_info "Description: $PACKAGE_DESCRIPTION"
}

# Function to detect distribution
detect_distro() {
    if [ -f /etc/os-release ]; then
        . /etc/os-release
        echo $ID
    elif [ -f /etc/redhat-release ]; then
        echo "rhel"
    elif [ -f /etc/debian_version ]; then
        echo "debian"
    else
        echo "unknown"
    fi
}

# Function to get package manager and install command
get_package_manager() {
    local distro=$(detect_distro)
    case $distro in
        fedora|rhel|centos|rocky|almalinux)
            echo "dnf"
            ;;
        ubuntu|debian|pop|mint)
            echo "apt"
            ;;
        opensuse*|sles)
            echo "zypper"
            ;;
        arch|manjaro)
            echo "pacman"
            ;;
        *)
            echo "unknown"
            ;;
    esac
}

# Function to check required tools
check_dependencies() {
    log_info "Checking required dependencies..."
    
    local missing_deps=()
    local missing_packages=()
    
    # Check for required commands
    for cmd in python3 dpkg-deb fakeroot; do
        if ! command -v $cmd >/dev/null 2>&1; then
            missing_deps+=($cmd)
            case $cmd in
                dpkg-deb)
                    missing_packages+=(dpkg)
                    ;;
                *)
                    missing_packages+=($cmd)
                    ;;
            esac
        fi
    done
    
    if [ ${#missing_deps[@]} -ne 0 ]; then
        log_error "Missing required dependencies: ${missing_deps[*]}"
        
        local pm=$(get_package_manager)
        local distro=$(detect_distro)
        
        case $pm in
            dnf)
                log_info "Install them with: sudo dnf install ${missing_packages[*]}"
                ;;
            apt)
                log_info "Install them with: sudo apt-get install ${missing_packages[*]}"
                ;;
            zypper)
                log_info "Install them with: sudo zypper install ${missing_packages[*]}"
                ;;
            pacman)
                log_info "Install them with: sudo pacman -S ${missing_packages[*]}"
                ;;
            *)
                log_info "Please install the following packages using your system package manager: ${missing_packages[*]}"
                ;;
        esac
        
        exit 1
    fi
    
    log_info "All dependencies satisfied"
}

# Function to clean previous builds
clean_build() {
    log_info "Cleaning previous builds..."
    
    rm -rf "$DEBIAN_DIR"
    rm -rf "$BUILD_DIR"
    rm -rf "$DIST_DIR"
    rm -rf "$PROJECT_ROOT"/*.egg-info
    
    # Clean Python cache
    find "$PROJECT_ROOT" -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
    find "$PROJECT_ROOT" -name "*.pyc" -delete 2>/dev/null || true
    
    log_info "Cleanup completed"
}

# Function to create debian directory structure
create_debian_structure() {
    log_info "Creating Debian package structure..."
    
    # Create main debian directory
    mkdir -p "$DEBIAN_DIR"
    
    # Create DEBIAN control directory
    mkdir -p "$DEBIAN_DIR/DEBIAN"
    
    # Create package directory structure
    local pkg_dir="$DEBIAN_DIR/usr/lib/python3/dist-packages"
    local bin_dir="$DEBIAN_DIR/usr/bin"
    local doc_dir="$DEBIAN_DIR/usr/share/doc/$PACKAGE_NAME"
    
    mkdir -p "$pkg_dir"
    mkdir -p "$bin_dir"
    mkdir -p "$doc_dir"
    
    log_info "Directory structure created"
}

# Function to create control file
create_control_file() {
    log_info "Creating DEBIAN/control file..."
    
    # Calculate installed size (rough estimate)
    local installed_size=$(du -sk "$PROJECT_ROOT/OTSO" | cut -f1)
    
    cat > "$DEBIAN_DIR/DEBIAN/control" << EOF
Package: $PACKAGE_NAME
Version: $PACKAGE_VERSION-1
Section: python
Priority: optional
Architecture: all
Depends: python3 (>= 3.12), python3-numpy (>= 2.2.0), python3-pandas (>= 2.2.0), python3-requests, python3-psutil, python3-tqdm
Maintainer: $PACKAGE_AUTHOR <$PACKAGE_EMAIL>
Description: $PACKAGE_DESCRIPTION
 OTSO (Open-source Trajectory Simulation and Optimization) is a comprehensive
 geomagnetic cutoff computation tool designed for space physics applications.
 It provides functionality for particle trajectory tracing, geomagnetic field
 modeling, and cutoff rigidity calculations.
 .
 Key features:
  - Geomagnetic field modeling using IGRF
  - Particle trajectory simulation
  - Coordinate system transformations
  - Station management for ground-based observations
  - Command-line interface for automated processing
Homepage: $PACKAGE_URL
Installed-Size: $installed_size
EOF
    
    log_info "Control file created"
}

# Function to create postinst script
create_postinst_script() {
    log_info "Creating post-installation script..."
    
    cat > "$DEBIAN_DIR/DEBIAN/postinst" << 'EOF'
#!/bin/bash
set -e

case "$1" in
    configure)
        # Update Python package cache
        if command -v python3 >/dev/null 2>&1; then
            python3 -c "import sys; sys.path_importer_cache.clear()" 2>/dev/null || true
        fi
        
        # Create symbolic links for console scripts if they don't exist
        for script in OTSO.clean OTSO.addstation OTSO.removestation OTSO.liststations OTSO.IGRFupdate OTSO.serverdownload; do
            if [ ! -f "/usr/bin/$script" ] && [ -f "/usr/lib/python3/dist-packages/OTSO/otso_cli.py" ]; then
                # Note: The actual console scripts will be handled by setuptools entry points
                echo "Console script $script available via Python module"
            fi
        done
        ;;
    *)
        ;;
esac

exit 0
EOF
    
    chmod 755 "$DEBIAN_DIR/DEBIAN/postinst"
    log_info "Post-installation script created"
}

# Function to create prerm script
create_prerm_script() {
    log_info "Creating pre-removal script..."
    
    cat > "$DEBIAN_DIR/DEBIAN/prerm" << 'EOF'
#!/bin/bash
set -e

case "$1" in
    remove|upgrade|deconfigure)
        # Clean up Python cache for this package
        find /usr/lib/python3/dist-packages/OTSO -name "*.pyc" -delete 2>/dev/null || true
        find /usr/lib/python3/dist-packages/OTSO -name "__pycache__" -type d -exec rmdir {} + 2>/dev/null || true
        ;;
    *)
        ;;
esac

exit 0
EOF
    
    chmod 755 "$DEBIAN_DIR/DEBIAN/prerm"
    log_info "Pre-removal script created"
}

# Function to copy package files
copy_package_files() {
    log_info "Copying package files..."
    
    # Copy OTSO package
    cp -r "$PROJECT_ROOT/OTSO" "$DEBIAN_DIR/usr/lib/python3/dist-packages/"
    
    # Remove non-Linux libraries from the package
    log_info "Removing non-Linux libraries from package..."
    local libs_dir="$DEBIAN_DIR/usr/lib/python3/dist-packages/OTSO/_core/libs"
    
    if [ -d "$libs_dir" ]; then
        # Remove Windows libraries (.pyd files)
        find "$libs_dir" -name "*.pyd" -delete 2>/dev/null || true
        
        # Remove macOS libraries (darwin.so files)
        find "$libs_dir" -name "*darwin.so" -delete 2>/dev/null || true
        
        # Count remaining libraries
        local remaining_libs=$(find "$libs_dir" -name "MiddleMan.*" | wc -l)
        log_info "Kept $remaining_libs Linux library files"
        
        # List remaining libraries for verification
        if [ $remaining_libs -gt 0 ]; then
            log_info "Remaining libraries:"
            find "$libs_dir" -name "MiddleMan.*" -exec basename {} \; | sed 's/^/  - /'
        fi
    fi
    
    # Copy documentation
    cp "$PROJECT_ROOT/README.md" "$DEBIAN_DIR/usr/share/doc/$PACKAGE_NAME/"
    cp "$PROJECT_ROOT/LICENSE" "$DEBIAN_DIR/usr/share/doc/$PACKAGE_NAME/"
    
    # Create changelog
    cat > "$DEBIAN_DIR/usr/share/doc/$PACKAGE_NAME/changelog" << EOF
$PACKAGE_NAME ($PACKAGE_VERSION-1) unstable; urgency=low

  * Version $PACKAGE_VERSION release
  * Built from source on $(date)

 -- $PACKAGE_AUTHOR <$PACKAGE_EMAIL>  $(date -R)
EOF
    
    # Compress changelog
    gzip -9 "$DEBIAN_DIR/usr/share/doc/$PACKAGE_NAME/changelog"
    
    # Set proper permissions
    find "$DEBIAN_DIR" -type f -exec chmod 644 {} \;
    find "$DEBIAN_DIR" -type d -exec chmod 755 {} \;
    chmod 755 "$DEBIAN_DIR/DEBIAN/postinst"
    chmod 755 "$DEBIAN_DIR/DEBIAN/prerm"
    
    log_info "Files copied and permissions set"
}

# Function to build the debian package
build_package() {
    log_info "Building Debian package..."
    
    cd "$PROJECT_ROOT"
    
    # Create output directory
    mkdir -p "$PROJECT_ROOT/src/debian"
    
    local deb_filename="${PACKAGE_NAME}_${PACKAGE_VERSION}-1_all.deb"
    local output_path="$PROJECT_ROOT/src/debian/$deb_filename"
    
    # Build the package
    fakeroot dpkg-deb --build "$DEBIAN_DIR" "$output_path"
    
    if [ -f "$output_path" ]; then
        log_info "Package built successfully: $output_path"
        
        # Display package information
        echo -e "\n${BLUE}Package Information:${NC}"
        dpkg --info "$output_path"
        
        # Display file size
        local file_size=$(du -h "$output_path" | cut -f1)
        log_info "Package size: $file_size"
        
    else
        log_error "Package build failed"
        exit 1
    fi
}

# Function to verify the package
verify_package() {
    log_info "Verifying package contents..."
    
    local deb_filename="${PACKAGE_NAME}_${PACKAGE_VERSION}-1_all.deb"
    local package_path="$PROJECT_ROOT/src/debian/$deb_filename"
    
    if [ -f "$package_path" ]; then
        echo -e "\n${BLUE}Package Contents:${NC}"
        dpkg --contents "$package_path"
        
        # Check for common issues
        echo -e "\n${BLUE}Package Validation:${NC}"
        
        # Check if main module exists
        if dpkg --contents "$package_path" | grep -q "usr/lib/python3/dist-packages/OTSO/__init__.py"; then
            log_info "✓ Main OTSO module found"
        else
            log_warn "✗ Main OTSO module not found"
        fi
        
        # Check if documentation exists
        if dpkg --contents "$package_path" | grep -q "usr/share/doc/$PACKAGE_NAME/README.md"; then
            log_info "✓ Documentation found"
        else
            log_warn "✗ Documentation not found"
        fi
        
        log_info "Package verification completed"
    else
        log_error "Package file not found for verification"
    fi
}

# Main execution
main() {
    echo -e "${GREEN}Starting Debian package build process...${NC}\n"
    
    extract_package_info
    check_dependencies
    clean_build
    create_debian_structure
    create_control_file
    create_postinst_script
    create_prerm_script
    copy_package_files
    build_package
    verify_package
    
    echo -e "\n${GREEN}========================================${NC}"
    echo -e "${GREEN}     Debian Package Build Complete!    ${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo -e "Package location: ${YELLOW}$PROJECT_ROOT/src/debian/${PACKAGE_NAME}_${PACKAGE_VERSION}-1_all.deb${NC}"
    echo -e "\nTo install the package:"
    echo -e "${BLUE}sudo dpkg -i $PROJECT_ROOT/src/debian/${PACKAGE_NAME}_${PACKAGE_VERSION}-1_all.deb${NC}"
    echo -e "${BLUE}sudo apt-get install -f${NC}  # Fix dependencies if needed"
}

# Handle command line arguments
case "${1:-build}" in
    "build"|"")
        main
        ;;
    "clean")
        log_info "Cleaning build artifacts..."
        clean_build
        log_info "Clean completed"
        ;;
    "help"|"-h"|"--help")
        echo "OTSO Debian Package Builder"
        echo ""
        echo "Usage: $0 [command]"
        echo ""
        echo "Commands:"
        echo "  build    Build the Debian package (default)"
        echo "  clean    Clean build artifacts"
        echo "  help     Show this help message"
        echo ""
        ;;
    *)
        log_error "Unknown command: $1"
        echo "Use '$0 help' for usage information"
        exit 1
        ;;
esac