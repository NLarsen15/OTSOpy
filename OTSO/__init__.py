import os
import sys


def _ensure_setup_complete():
    from .Callibrate import setup as setup_func
    package_dir = os.path.dirname(os.path.abspath(__file__))
    setup_marker = os.path.join(package_dir, 'setup_complete.txt')
    if not os.path.exists(setup_marker):
        setup_func()
    return

def cutoff(*args, **kwargs):
    _ensure_setup_complete()
    from .Cutoff import cutoff as cutoff_func
    return cutoff_func(*args, **kwargs)

def cone(*args, **kwargs):
    _ensure_setup_complete()
    from .Cone import cone as cone_func
    return cone_func(*args, **kwargs)

def planet(*args, **kwargs):
    _ensure_setup_complete()
    from .Planet import planet as planet_func
    return planet_func(*args, **kwargs)

def trajectory(*args, **kwargs):
    _ensure_setup_complete()
    from .Trajectory import trajectory as trajectory_func
    return trajectory_func(*args, **kwargs)

def flight(*args, **kwargs):
    _ensure_setup_complete()
    from .Flight import flight as flight_func
    return flight_func(*args, **kwargs)

def magfield(*args, **kwargs):
    _ensure_setup_complete()
    from .Magfield import magfield as magfield_func
    return magfield_func(*args, **kwargs)

def coordtrans(*args, **kwargs):
    _ensure_setup_complete()
    from .Coordtrans import coordtrans as coordtrans_func
    return coordtrans_func(*args, **kwargs)

def trace(*args, **kwargs):
    _ensure_setup_complete()
    from .Trace import trace as trace_func
    return trace_func(*args, **kwargs)

def clean(*args, **kwargs):
    from .OTSO import clean as clean_func
    return clean_func(*args, **kwargs)