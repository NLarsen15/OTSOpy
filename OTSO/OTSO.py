import sys, os
import subprocess
import platform
import shutil

def clean():
    print("Cleaning OTSO...")
    Delete()
    print("OTSO cleaned")


def Delete():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    setupfile = os.path.join(script_dir, 'setup_complete.txt')

    if os.path.exists(file_path):
        try:
            os.remove(file_path)
        except Exception as e:

    server_data_folder_path = os.path.join(script_dir, 'Parameters', 'functions', 'ServerData')
    if os.path.exists(server_data_folder_path):
        try:
            shutil.rmtree(server_data_folder_path)
        except Exception as e:
    
