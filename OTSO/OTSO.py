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

    if os.path.exists(setupfile):
        os.remove(setupfile)

    server_data_folder_path = os.path.join(script_dir, 'Parameters', 'functions', 'ServerData')
    if os.path.exists(server_data_folder_path):
        shutil.rmtree(server_data_folder_path)
    
