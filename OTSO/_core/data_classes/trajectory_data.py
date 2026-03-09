
import numpy as np
from typing import Union, Optional, Sequence

from dataclasses import dataclass, field


@dataclass
class TrajectoryData:
    Stations: Union[str, Sequence[str]]
    customlocations: Sequence[list]
    rigidity: float
    startaltitude: float
    minaltitude: float
    maxdistance: float
    maxtime: float
    zenith: float
    azimuth: float
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
    Anum: int
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
    coordsystem: str
    gyropercent: float
    mpause: str
    corenum: int
    g: Sequence[float]
    h: Sequence[float]
    MHDfile: str
    MHDcoordsys: str
    spheresize: float
    inputcoord: str
    Verbose: bool
    AdaptiveExternalModel: bool
    mindist: float
    adapt: bool
    Berr: float
    totalbetacheck: bool
    maxsteps: int

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
    Kp: Optional[float] = field(default=None)
