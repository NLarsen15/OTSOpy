
import numpy as np
from typing import Optional, Sequence
from datetime import date

from dataclasses import dataclass, field


@dataclass
class FlightData:
    latitudes: Sequence[float]
    longitudes: Sequence[float]
    dates: Sequence[date]
    altitudes: Sequence[float]
    cutoff_comp: str
    minaltitude: float
    maxdistance: float
    maxtime: float
    serverdata: str
    livedata: str
    vx: Sequence[float]
    vy: Sequence[float]
    vz: Sequence[float]
    bx: Sequence[float]
    by: Sequence[float]
    bz: Sequence[float]
    density: Sequence[float]
    pdyn: Sequence[float]
    Dst: Sequence[float]
    G1: Sequence[float]
    G2: Sequence[float]
    G3: Sequence[float]
    W1: Sequence[float]
    W2: Sequence[float]
    W3: Sequence[float]
    W4: Sequence[float]
    W5: Sequence[float]
    W6: Sequence[float]
    kp: Sequence[float]
    by_avg: Sequence[float]
    bz_avg: Sequence[float]
    n_index: Sequence[float]
    b_index: Sequence[float]
    sym_h_corrected: Sequence[float]
    Anum: int
    anti: str
    internalmag: str
    externalmag: str
    boberg: bool
    bobergtype: str
    intmodel: str
    startrigidity: float
    endrigidity: float
    rigiditystep: float
    rigidityscan: str
    coordsystem: str
    gyropercent: float
    mpause: str
    corenum: int
    azimuth: float
    zenith: float
    g: Sequence[float]
    h: Sequence[float]
    asymptotic: str
    asymlevels: Sequence[float]
    unit: str
    MHDfile: str
    MHDcoordsys: str
    spheresize: float
    inputcoord: str
    Verbose: bool
    AdaptiveExternalModel: bool
    mindist: float
    delim: str
    adapt: bool
    Berr: float
    totalbetacheck: bool
    maxsteps: int

    rigidityarray: Optional[Sequence[float]] = field(default=None)
    datearraylist: Optional[Sequence[Sequence[int]]] = field(default=None)
    model: Optional[Sequence[int]] = field(default=None)
    integrationmodel: Optional[int] = field(default=None)
    particlearray: Optional[Sequence[int]] = field(default=None)
    IOPTlist: Optional[Sequence[int]] = field(default=None)
    windarraylist: Optional[Sequence[Sequence[float]]] = field(default=None)
    magnetopause: Optional[int] = field(default=None)
    coordinatesystem: Optional[str] = field(default=None)
    maxsteppercent: Optional[float] = field(default=None)
    endparams: Optional[Sequence[float]] = field(default=None)
    station_array: Optional[Sequence[np.ndarray]] = field(default=None)
    Rcomp: Optional[int] = field(default=None)
    Rscan: Optional[int] = field(default=None)
    Kplist: Optional[Sequence[float]] = field(default=None)
    glist: Optional[Sequence[float]] = field(default=None)
    hlist: Optional[Sequence[float]] = field(default=None)
