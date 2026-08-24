import glob
import os

from .extract import LIVEDATA_DIR

def remove_files() -> None:
    file_list = ["Kp_data.txt", "Magnetic_data.csv", "Magnetic_data.json", "space_data.csv", "space_data.json"]
    script_dir = LIVEDATA_DIR

    file_paths = [os.path.join(script_dir, file_name) for file_name in file_list]
    file_paths += glob.glob(os.path.join(script_dir, "Dst_data_*.txt"))
    for file_path in file_paths:
        try:
            os.remove(file_path)
        except FileNotFoundError:
            print(f"File not found: {file_path}")
        except PermissionError:
            print(f"Permission denied: {file_path}")
        except Exception as e:
            print(f"Error deleting {file_path}: {e}")