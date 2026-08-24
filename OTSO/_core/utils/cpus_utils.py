import os
import sys
import ctypes
from ctypes import wintypes


def set_process_affinity(cpus):

    # ---------------------------------------------------------
    # Windows
    # ---------------------------------------------------------
    if sys.platform == "win32":

        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)

        kernel32.GetCurrentProcess.restype = wintypes.HANDLE

        kernel32.SetProcessAffinityMask.argtypes = [
            wintypes.HANDLE,
            ctypes.c_size_t,
        ]
        kernel32.SetProcessAffinityMask.restype = wintypes.BOOL

        mask = 0

        for cpu in cpus:
            if cpu < 0:
                raise ValueError(f"Invalid CPU: {cpu}")

            mask |= 1 << cpu

        handle = kernel32.GetCurrentProcess()

        if not kernel32.SetProcessAffinityMask(handle, mask):
            raise ctypes.WinError(ctypes.get_last_error())

        return


    # ---------------------------------------------------------
    # Linux
    # ---------------------------------------------------------
    if sys.platform.startswith("linux"):

        os.sched_setaffinity(0, cpus)

        return


    # ---------------------------------------------------------
    # macOS
    # ---------------------------------------------------------
    if sys.platform == "darwin":

        cpu_list = sorted(cpus)

        os.environ["OMP_NUM_THREADS"] = str(len(cpu_list))
        os.environ["OMP_PROC_BIND"] = "TRUE"
        os.environ["OMP_PLACES"] = "{" + ",".join(
            str(cpu) for cpu in cpu_list
        ) + "}"

        return


    raise RuntimeError(
        f"CPU affinity is not implemented for {sys.platform}"
    )