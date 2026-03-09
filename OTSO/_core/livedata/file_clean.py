import os

def remove_files() -> None:
    file_list = ["Dst_data.txt", "Kp_data.txt", "Magnetic_data.csv", "Magnetic_data.json", "space_data.csv", "space_data.json"]
    script_dir = os.path.dirname(__file__)
    for file_name in file_list:
        file_path = os.path.join(script_dir, file_name)
        try:
            os.remove(file_path)
        except FileNotFoundError:
            print(f"File not found: {file_path}")
        except PermissionError:
            print(f"Permission denied: {file_path}")
        except Exception as e:
            print(f"Error deleting {file_path}: {e}")