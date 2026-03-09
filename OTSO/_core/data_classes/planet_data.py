
import numpy as np
from typing import Optional, Sequence

from dataclasses import dataclass, field


@dataclass
class PlanetData:
    startaltitude: float
    cutoff_comp: str
    minaltitude: float
    maxdistance: float
    maxtime: float
    serverdata: str
    livedata: str
    vx: float
    vy: float
    vz: float
    bx: float
    by: float
    bz: float
    density: float
    pdyn: float
    Dst: float
    G1: float
    G2: float
    G3: float
    W1: float
    W2: float
    W3: float
    W4: float
    W5: float
    W6: float
    kp: float
    by_avg: float
    bz_avg: float
    n_index: float
    b_index: float
    sym_h_corrected: float
    anti: str
    year: int   
    month: int
    day: int
    hour: int
    minute: int
    second: int
    internalmag: str
    externalmag: str
    boberg: bool
    bobergtype: str
    intmodel: str
    startrigidity: float
    endrigidity: float
    rigiditystep: float
    rigidityscan: str
    gyropercent: float
    mpause: str
    corenum: int
    azimuth: float
    zenith: float
    asymptotic: str
    asymlevels: Sequence[float]
    unit: str
    latstep: float
    longstep: float
    maxlat: float
    minlat: float
    maxlong: float
    minlong: float
    g: Sequence[float]
    h: Sequence[float]
    MHDfile: str
    MHDcoordsys: str
    spheresize: float
    inputcoord: str
    Verbose: bool
    AdaptiveExternalModel: bool
    array_of_lats_and_longs: Optional[Sequence[float]]
    grid_params_user_set: Optional[bool]
    mindist: float
    delim: str
    adapt: bool
    Berr: float
    totalbetacheck: bool
    maxsteps: int

    rigidityarray: Optional[Sequence[float]] = field(default=None)
    datearray: Optional[Sequence[int]] = field(default=None)
    model: Optional[Sequence[int]] = field(default=None)
    integrationmodel: Optional[int] = field(default=None)
    particlearray: Optional[Sequence[int]] = field(default=None)
    IOPT: Optional[int] = field(default=None)
    windarray: Optional[Sequence[float]] = field(default=None)
    magnetopause: Optional[int] = field(default=None)
    coordinatesystem: Optional[str] = field(default=None)
    maxsteppercent: Optional[float] = field(default=None)
    endparams: Optional[Sequence[float]] = field(default=None)
    station_array: Optional[Sequence[np.ndarray]] = field(default=None)
    Rcomp: Optional[int] = field(default=None)
    Rscan: Optional[int] = field(default=None)
    Kp: Optional[float] = field(default=None)
    LatitudeList_meta: Optional[Sequence[float]] = field(default=None)
    LongitudeList_meta: Optional[Sequence[float]] = field(default=None)
    coordinate_pairs: Optional[Sequence[Sequence[float]]] = field(default=None)
    custom_coords_provided: bool = field(default=False)