"""Equal-area global grids for planet (grid_params["grid_type"]).
"""
import numpy as np


def healpix(nside):
    """HEALPix pixel centres (RING scheme): 12 * nside**2 equal-area points on 4*nside - 1
    iso-latitude rings. Returns an (N, 2) array of [lat, lon] in degrees, lon in 0-360."""
    out = []
    for i in range(1, 4 * nside):
        if i < nside:                         # north polar cap
            z, n, shift = 1 - i**2 / (3 * nside**2), 4 * i, 0.5
        elif i <= 3 * nside:                  # equatorial belt
            z, n = 4 / 3 - 2 * i / (3 * nside), 4 * nside
            shift = 0.5 if (i - nside) % 2 == 0 else 0.0
        else:                                 # south polar cap
            k = 4 * nside - i
            z, n, shift = -(1 - k**2 / (3 * nside**2)), 4 * k, 0.5
        lat = np.degrees(np.arcsin(z))
        out += [(lat, (j + shift) * 360 / n) for j in range(n)]
    return np.array(out)


def fibonacci(n):
    """Golden-angle (Fibonacci) spiral: n near-equal-area points for any n, with no rings or
    hierarchy. Returns an (N, 2) array of [lat, lon] in degrees, lon in 0-360."""
    i = np.arange(n) + 0.5
    lat = np.degrees(np.arcsin(1 - 2 * i / n))
    lon = np.degrees(np.pi * (1 + 5**0.5) * i) % 360
    return np.column_stack([lat, lon])


def dipole_pole(g, h):
    """Geographic latitude and longitude (degrees) of the north centred-dipole pole from the
    degree-1 Gauss coefficients (g[1] = g10, g[2] = g11, h[2] = h11, OTSO's ordering)."""
    g10, g11, h11 = g[1], g[2], h[2]
    b0 = np.sqrt(g10**2 + g11**2 + h11**2)
    return 90 - np.degrees(np.arccos(-g10 / b0)), np.degrees(np.arctan2(-h11, -g11))


def to_geomagnetic_frame(points, pole):
    """Treat `points` as centred-dipole [lat, lon] and rotate them to geographic [lat, lon]."""
    la, lo = np.radians(points[:, 0]), np.radians(points[:, 1])
    X = np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    t, p = np.radians(90 - pole[0]), np.radians(pole[1])
    Ry = np.array([[np.cos(t), 0, np.sin(t)], [0, 1, 0], [-np.sin(t), 0, np.cos(t)]])
    Rz = np.array([[np.cos(p), -np.sin(p), 0], [np.sin(p), np.cos(p), 0], [0, 0, 1]])
    G = X @ (Rz @ Ry).T
    return np.column_stack([np.degrees(np.arcsin(np.clip(G[:, 2], -1, 1))),
                            np.degrees(np.arctan2(G[:, 1], G[:, 0])) % 360])


WGS84_E2 = 0.00669437999014  # first eccentricity squared


def geocentric_to_geodetic(points):
    """Convert geocentric [lat, lon] to WGS84 geodetic [lat, lon] at the surface."""
    lat = np.degrees(np.arctan(np.tan(np.radians(points[:, 0])) / (1 - WGS84_E2)))
    return np.column_stack([lat, points[:, 1]])


def healpix_geomag(nside, g, h):
    """HEALPix rings laid out in centred-dipole coordinates for the date of g, h."""
    return to_geomagnetic_frame(healpix(nside), dipole_pole(g, h))


GRID_TYPES = ("regular", "healpix", "healpix_geomag", "fibonacci")


def generate_grid(grid_type, nside, g, h, geodetic=False):
    """[lat, lon] pairs for an equal-area grid_type (not "regular"). Every type gives
    12 * nside**2 points, so fibonacci maps match the HEALPix sizes. Points are built as
    geocentric directions; geodetic=True converts their latitudes for GDZ input."""
    if grid_type not in GRID_TYPES:
        raise ValueError(f"grid_type must be one of {GRID_TYPES}, got {grid_type!r}.")
    if isinstance(nside, bool) or not isinstance(nside, (int, np.integer)) or nside < 1:
        raise ValueError(f"nside must be a positive integer, got {nside!r}.")
    if grid_type == "healpix":
        points = healpix(nside)
    elif grid_type == "fibonacci":
        points = fibonacci(12 * nside**2)
    else:
        points = healpix_geomag(nside, g, h)
    return geocentric_to_geodetic(points) if geodetic else points
