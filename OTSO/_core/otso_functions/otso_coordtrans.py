import time
import multiprocessing as mp
import pandas as pd
import sys
import queue
import numpy as np
from tqdm import tqdm

from ..utils.igrf_utils import compute_gauss_coefficients
from ..custom_classes import date
from ..fortran_calls import coordtrans
from ..readme_generators import coordtrans_readme

def OTSO_coordtrans(Locations,Dates,CoordIN,CoordOUT,corenum,Verbose):
    
    if CoordIN not in ["GDZ","GEO","GSM","GSE","SM","GEI","MAG","SPH","RLL"]:
         print("Please select a valid CoordIN: ""GDZ"", ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""SPH"", ""RLL""")
         exit()
    if CoordOUT not in ["GDZ","GEO","GSM","GSE","SM","GEI","MAG","SPH","RLL"]:
         print("Please select a valid CoordOUT: ""GDZ"", ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""SPH"", ""RLL""")
         exit()

    ChildProcesses = []
    results = []
    DateArrayList = []
    hlist = []
    glist = []

    N_COEFFS = (13 + 1) * (13 + 2) // 2   # maxdegree 13, as passed to the Fortran below
    gauss_cache = {}   # IGRF coefficients once per distinct date, not once per location
    for x in Dates:
          DateCreate = date.Date(x)
          DateArray = DateCreate.GetDate()
          DateArrayList.append(DateArray)
          key = tuple(np.asarray(DateArray).ravel().tolist())
          if key not in gauss_cache:
              gauss_cache[key] = compute_gauss_coefficients(DateArray)
          gausscoefs = gauss_cache[key]
          # Trim to degree 13 (105 coefficients) like every other OTSO function: GEOPACK's RECALC_08
          # allocates its G/H arrays once, at the first call's length, so passing the untrimmed
          # 136 after another function had run (105) overflowed them and corrupted the heap.
          hlist.append(gausscoefs['H_coefficients'][:N_COEFFS])
          glist.append(gausscoefs['G_coefficients'][:N_COEFFS])

    IndexList = np.array_split(np.arange(len(Locations)), corenum)
    LocationsList = np.array_split(Locations, corenum)
    DateArrayList = np.array_split(DateArrayList, corenum)
    glist = np.array_split(glist, corenum)
    hlist = np.array_split(hlist, corenum)

    maxdegree = 13


    start = time.time()

    if Verbose:
        print("OTSO Coordtrans Computation Started")

    total_stations = len(Locations)
    results = []

    if corenum == 1:
        num_batches = len(LocationsList)
        
        progress_bar = None
        if Verbose and tqdm is not None:
            progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" transformation")
        elif Verbose:
            print(f"Processing {total_stations} coordinates...")

        class SimpleQueue:
            def __init__(self):
                self.items = []
            def put(self, item):
                self.items.append(item)
            def get_all(self):
                return self.items

        simple_queue = SimpleQueue()
        processed = 0
        
        for Data, Date, g, h, Index in zip(LocationsList, DateArrayList, glist, hlist, IndexList):
            coordtrans.FortranCoordtrans(Data, Date, CoordIN, CoordOUT, simple_queue, maxdegree, g, h, Index)
            
            processed += len(Data)
            if Verbose:
                if progress_bar is not None:
                    progress_bar.update(len(Data))
                else:
                    percent_complete = (processed / total_stations) * 100
                    sys.stdout.write(f"\r{percent_complete:.2f}% complete")
                    sys.stdout.flush()
        
        results = simple_queue.get_all()
        
        if progress_bar is not None:
            progress_bar.close()

    else:
        try:
            if not mp.get_start_method(allow_none=True):
                mp.set_start_method('spawn')
        except RuntimeError:
            pass
        
        ProcessQueue = mp.Manager().Queue()
        for Data,Date,g,h,Index in zip(LocationsList,DateArrayList,glist,hlist,IndexList):
            Child = mp.Process(target=coordtrans.FortranCoordtrans,  args=(Data, Date, CoordIN, CoordOUT, ProcessQueue, maxdegree, g, h, Index))
            ChildProcesses.append(Child)

        for a in ChildProcesses:
            a.start()

        processed = 0

        progress_bar = None
        if Verbose and tqdm is not None:
            progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" transformation")
        elif Verbose:
            print(f"Processing {total_stations} coordinates...")

        while processed < total_stations:
          try:
            result_df = ProcessQueue.get(timeout=0.001)
            results.append(result_df)
            processed += len(result_df)          # one DataFrame per worker chunk

            if Verbose:
                if progress_bar is not None:
                    progress_bar.update(len(result_df))
                else:
                    percent_complete = (processed / total_stations) * 100
                    sys.stdout.write(f"\r{percent_complete:.2f}% complete")
                    sys.stdout.flush()

          except queue.Empty:
            if all(not Child.is_alive() for Child in ChildProcesses) and ProcessQueue.empty():
                print("ERROR: An OTSO worker process stopped before finishing. Please check the error message above.")
                exit()
          
          time.sleep(0.0001)

        if progress_bar is not None:
            progress_bar.close()

        for b in ChildProcesses:
            b.join()

    combined_df = pd.concat(results)
    combined_df[combined_df.columns[1:]] = combined_df[combined_df.columns[1:]].astype(float)
    sorted_df = combined_df.sort_values(by=combined_df.columns[:4].tolist())

    stop = time.time()
    Printtime = round((stop-start),3)
    if Verbose:
        print("\nOTSO Coordtrans Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")
    
    README = coordtrans_readme.READMECoordtrans(CoordIN,CoordOUT,Printtime)
    
    return [sorted_df, README]