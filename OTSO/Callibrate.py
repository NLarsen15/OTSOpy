import sys, os
import subprocess
import platform

def setup():
    #print("Running setup tasks...")
    a = DeleteLibs()
    #DarwinLib()
    Generatefile(a)
    #print("Setup complete")

def DarwinLib():
    if platform.system() == "darwin": 
        so_files = []
        script_dir = os.path.dirname(os.path.abspath(__file__))
        package_path = os.path.join(script_dir, 'Parameters', 'functions')

        for root, _, files in os.walk(package_path):
            for file in files:
                if file.endswith(('.so', '.pyd')): 
                    so_files.append(os.path.join(root, file))

        for lib_path in so_files:
            try:
                subprocess.run(['install_name_tool', '-add_rpath', '@loader_path', package_path], check=True)
                print(f"Updated Rpath for {package_path}")
            except subprocess.CalledProcessError as e:
                print(f"Error updating Rpath for {package_path}: {e}")


def DeleteLibs():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    package_path = os.path.join(script_dir, 'Parameters', 'functions')
    system_type = platform.system().lower()
    python_version = f"{sys.version_info.major}{sys.version_info.minor}"

    if system_type == "windows":
        system_type = "win"

    error_messages = []

    for filename in os.listdir(package_path):
        file_path = os.path.join(package_path, filename)
        if os.path.isfile(file_path) and filename.endswith(('.so', '.pyd')):
            if system_type not in filename or python_version not in filename:
                try:
                    os.remove(file_path)
                except Exception as e:
                    error_messages.append(f"Could not delete {file_path}: {e}")

    return "\n".join(error_messages) if error_messages else ""


def Generatefile(a):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    marker_path = os.path.join(script_dir, 'setup_complete.txt')

    try:
        with open(marker_path, 'w') as f:
            f.write('Setup has been completed.\n')
            f.write(a)
        #print(f"Setup marker file created at: {marker_path}")
    except Exception as e:
        print(f"Error creating setup marker file: {e}")

