import sys, os
import subprocess
import platform

def setup():
    #print("Running setup tasks...")
    DeleteLibs()
    #DarwinLib()
    Generatefile()
    #print("Setup complete")

def DarwinLib():
    """Detect if user is on macOS and automatically add Rpath to .so and .pyd files in OTSO package"""
    if platform.system() == "darwin":  # Check if running on macOS
        so_files = []
        script_dir = os.path.dirname(os.path.abspath(__file__))
        package_path = os.path.join(script_dir, 'Parameters', 'functions')

        # Walk through the package directory to find shared object files
        for root, _, files in os.walk(package_path):
            for file in files:
                if file.endswith(('.so', '.pyd')):  # Collect relevant files
                    so_files.append(os.path.join(root, file))

        # Apply Rpath modification to each detected file
        for lib_path in so_files:
            try:
                subprocess.run(['install_name_tool', '-add_rpath', '@loader_path', package_path], check=True)
                print(f"Updated Rpath for {package_path}")
            except subprocess.CalledProcessError as e:
                print(f"Error updating Rpath for {package_path}: {e}")


def DeleteLibs():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    package_path = os.path.join(script_dir, 'Parameters', 'functions')
    system_type = platform.system().lower()  # Get OS type (windows, linux, darwin for macOS)
    python_version = f"{sys.version_info.major}{sys.version_info.minor}"  # Extract Python version (e.g., cp310 for Python 3.10)

    if system_type == "windows":
        system_type = "win"  # Normalize to "win" for Windows systems

    for filename in os.listdir(package_path):
        file_path = os.path.join(package_path, filename)
        
        # Ensure it's a file before processing
        if os.path.isfile(file_path):
            if filename.endswith(('.so', '.pyd')):
                # Check if file matches system type and Python version
                if system_type not in filename or python_version not in filename:
                    #print(f"Would delete: {file_path}")
                    os.remove(file_path)

def Generatefile():
    """Creates a file in the package directory to signal that setup has been completed"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    marker_path = os.path.join(script_dir, 'setup_complete.txt')

    try:
        with open(marker_path, 'w') as f:
            f.write('Setup has been completed.\n')
        #print(f"Setup marker file created at: {marker_path}")
    except Exception as e:
        print(f"Error creating setup marker file: {e}")

