
from typing import Union, Optional, Sequence

from dataclasses import dataclass, field


@dataclass
class MagfieldData:
    locations: Union[str, Sequence[str]]
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
    inputcoord: str
    coordout: str
    g: Sequence[float]
    h: Sequence[float]
    corenum: int
    MHDfile: str
    MHDcoordsys: str
    Verbose: bool

    datearray: Optional[Sequence[int]] = field(default=None)
    model: Optional[Sequence[int]] = field(default=None)
    IOPT: Optional[int] = field(default=None)
    windarray: Optional[Sequence[float]] = field(default=None)
    Kp: Optional[float] = field(default=None)
