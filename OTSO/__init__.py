def cutoff(*args, **kwargs):
    from .Cutoff import cutoff as cutoff_func
    return cutoff_func(*args, **kwargs)

def cone(*args, **kwargs):
    from .Cone import cone as cone_func
    return cone_func(*args, **kwargs)

def planet(*args, **kwargs):
    from .Planet import planet as planet_func
    return planet_func(*args, **kwargs)

def trajectory(*args, **kwargs):
    from .Trajectory import trajectory as trajectory_func
    return trajectory_func(*args, **kwargs)

def flight(*args, **kwargs):
    from .Flight import flight as flight_func
    return flight_func(*args, **kwargs)

def magfield(*args, **kwargs):
    from .Magfield import magfield as magfield_func
    return magfield_func(*args, **kwargs)

def coordtrans(*args, **kwargs):
    from .Coordtrans import coordtrans as coordtrans_func
    return coordtrans_func(*args, **kwargs)

def trace(*args, **kwargs):
    from .Trace import trace as trace_func
    return trace_func(*args, **kwargs)

def setup(*args, **kwargs):
    from .OTSO import setup as setup_func
    return setup_func(*args, **kwargs)