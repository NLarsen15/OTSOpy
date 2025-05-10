import sys, os, csv
import subprocess
import platform
import shutil
import pandas as pd

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

def AddStation(Name, Latitude, Longitude):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    stationfile = os.path.join(script_dir, 'Parameters', 'functions', 'StationList.csv')

    df = pd.read_csv(stationfile)

    existing_row = df[df["Name"] == Name]

    if not existing_row.empty:
        print(f"Station '{Name}' already exists:")
        print(existing_row)

        user_input = input("Do you want to overwrite this entry? (y/n): ").strip().lower()
        
        if user_input != "y":
            print("No changes were made.")
            return

        df = df[df["Name"] != Name]

    new_station = pd.DataFrame([[Name, Latitude, Longitude]], columns=["Name", "Latitude", "Longitude"])
    df = pd.concat([df, new_station], ignore_index=True)

    df = df.sort_values(by="Name")

    df.to_csv(stationfile, index=False)

    print(f"Station '{Name}' has been added/updated successfully.")


def RemoveStation(Name):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    stationfile = os.path.join(script_dir, 'Parameters', 'functions', 'StationList.csv')

    df = pd.read_csv(stationfile)

    existing_row = df[df["Name"] == Name]

    if existing_row.empty:
        print(f"Station '{Name}' not found. No changes were made.")
        return

    print(f"Station '{Name}' found:")
    print(existing_row)

    user_input = input("Do you want to delete this entry? (y/n): ").strip().lower()
    
    if user_input != "y":
        print("No changes were made.")
        return

    df = df[df["Name"] != Name]

    df.to_csv(stationfile, index=False)

    print(f"Station '{Name}' has been removed successfully.")


def ListStations():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    stationfile = os.path.join(script_dir, 'Parameters', 'functions', 'StationList.csv')

    df = pd.read_csv(stationfile)

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)

    print(df.to_string(index=False))




