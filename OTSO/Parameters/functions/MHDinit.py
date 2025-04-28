from . import MiddleMan as OTSOLib
import time
import os
import shutil
from datetime import datetime
import pandas as pd
import glob
import csv
import numpy as np
import multiprocessing as mp
from . import date
import gc
import sys
import psutil
import tracemalloc

def MHDinitialise(MHDfile):
      data = pd.read_csv(MHDfile)
      x1 = data["X"].values.astype(np.float64)
      y1 = data["Y"].values.astype(np.float64)
      z1 = data["Z"].values.astype(np.float64)
      bx = data["Bx"].values.astype(np.float64)
      by = data["By"].values.astype(np.float64)
      bz = data["Bz"].values.astype(np.float64)
  
      XU = np.unique(x1)
      YU = np.unique(y1)
      ZU = np.unique(z1)
      
      XUlen, YUlen, ZUlen = len(XU), len(YU), len(ZU)
  
      ix = np.searchsorted(XU, x1) - 1
      iy = np.searchsorted(YU, y1) - 1
      iz = np.searchsorted(ZU, z1) - 1
  
      MHDposition = np.zeros((XUlen, YUlen, ZUlen, 3), dtype=np.float64)
      MHDB = np.zeros((XUlen, YUlen, ZUlen, 3), dtype=np.float64)
  
      for i in range(len(x1)):
        
        #if x1[i]**2 + y1[i]**2 + z1[i]**2 < 1.0:
        #  bx[i], by[i], bz[i] = 0.0, 0.0, 0.0
        MHDB[ix[i], iy[i], iz[i], :] = [bx[i], by[i], bz[i]]
        MHDposition[ix[i], iy[i], iz[i], :] = [x1[i], y1[i], z1[i]]
  
      # Chunking
      min_chunk = 10
      n_x_split = max(1, XUlen // min_chunk)
      n_y_split = max(1, YUlen // min_chunk)
      n_z_split = max(1, ZUlen // min_chunk)
  
      chunk_x = XUlen // n_x_split
      chunk_y = YUlen // n_y_split
      chunk_z = ZUlen // n_z_split
  
      num_regions = n_x_split * n_y_split * n_z_split
      region_info = []
  
      for r in range(num_regions):
          rx = r % n_x_split
          ry = (r // n_x_split) % n_y_split
          rz = r // (n_x_split * n_y_split)
  
          sx = rx * chunk_x
          sy = ry * chunk_y
          sz = rz * chunk_z
  
          ex = (rx + 1) * chunk_x if rx < n_x_split - 1 else XUlen
          ey = (ry + 1) * chunk_y if ry < n_y_split - 1 else YUlen
          ez = (rz + 1) * chunk_z if rz < n_z_split - 1 else ZUlen
  
          cx = np.mean(XU[sx:ex])
          cy = np.mean(YU[sy:ey])
          cz = np.mean(ZU[sz:ez])
          dist = np.sqrt(cx**2 + cy**2 + cz**2)
  
          region_info.append((dist, r+1, sx+1, ex+1, sy+1, ey+1, sz+1, ez+1))
  
      region_info.sort()
      region_order = [int(r[1]) for r in region_info]
      start_x = [int(r[2]) for r in region_info]
      end_x   = [int(r[3]) for r in region_info]
      start_y = [int(r[4]) for r in region_info]
      end_y   = [int(r[5]) for r in region_info]
      start_z = [int(r[6]) for r in region_info]
      end_z   = [int(r[7]) for r in region_info]
  
      # Convert all arrays to lists
      XU_list = XU.tolist()
      YU_list = YU.tolist()
      ZU_list = ZU.tolist()
  
      MHDposition_list = MHDposition.tolist()
      MHDB_list = MHDB.tolist()
      
      minX_global = np.min(MHDposition[:, :, :, 0])  # Min value in X axis for the whole grid
      maxX_global = np.max(MHDposition[:, :, :, 0])  # Max value in X axis for the whole grid
      minY_global = np.min(MHDposition[:, :, :, 1])  # Min value in Y axis for the whole grid
      maxY_global = np.max(MHDposition[:, :, :, 1])  # Max value in Y axis for the whole grid
      minZ_global = np.min(MHDposition[:, :, :, 2])  # Min value in Z axis for the whole grid
      maxZ_global = np.max(MHDposition[:, :, :, 2])  # Max value in Z axis for the whole grid

      #print("MHD file Loaded")

      del data
      
  
      # Call Fortran
      OTSOLib.mhdstartupsorted(
          XU_list, YU_list, ZU_list,
          MHDposition_list, MHDB_list, n_x_split, n_y_split, n_z_split,
          minX_global,maxX_global,minY_global,maxY_global,minZ_global,maxZ_global,
          region_order, start_x, end_x,
          start_y, end_y, start_z, end_z,
          num_regions,XUlen, YUlen, ZUlen)