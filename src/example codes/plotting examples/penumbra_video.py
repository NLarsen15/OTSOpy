import os
import shutil
import subprocess

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox, TransformedBbox

from OTSO import cone, trajectory


# ============================================================
# Location and video settings
# ============================================================

station = "OULU"
custom_location = None

output_file = "penumbra_video.mp4"
preview_time = None

shown_above = 3
shown_penumbra = 14
shown_below = 3
seconds_per_trajectory = 2.0
fps = 30
show_otso_branding = True

logo_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "images", "OTSO_logo.png")


# ============================================================
# OTSO input parameters
# ============================================================

date_inputs = {
    "year": 2024,
    "month": 1,
    "day": 1,
    "hour": 12,
    "minute": 0,
    "second": 0
}

magfield_inputs = {
    "internalmag": "IGRF",
    "externalmag": "TSY89c",
    "boberg": False,
    "bobergtype": "EXTENSION",
    "magnetopause": "Kobel",
    "spheresize": 10,
    "AdaptiveExternalModel": False
}

rigidity_inputs = {
    "startrigidity": 20,
    "endrigidity": 0,
    "rigiditystep": 0.005
}

solar_wind_inputs = {
    "vx": -500,
    "vy": 0,
    "vz": 0,
    "bx": 0,
    "by": 5,
    "bz": 5,
    "by_avg": 0,
    "bz_avg": 0,
    "density": 1,
    "pdyn": 0
}

geomagnetic_inputs = {
    "Dst": 0,
    "kp": 0,
    "n_index": 0,
    "b_index": 0,
    "sym_h_corrected": 0
}

tsyganenko_inputs = {
    "G1": 0,
    "G2": 0,
    "G3": 0,
    "W1": 0,
    "W2": 0,
    "W3": 0,
    "W4": 0,
    "W5": 0,
    "W6": 0
}

integration_inputs = {
    "intmodel": "Boris-Buneman",
    "gyropercent": 10,
    "minaltitude": 20,
    "maxdistance": 100,
    "maxtime": 0,
    "mintrapdist": 6.6,
    "startaltitude": 20,
    "betaerror": 0.001,
    "totalbetacheck": True,
    "adaptivestep": False,
    "maxsteps": 1000000
}

particle_inputs = {
    "Anum": 1,
    "anti": "YES",
    "zenith": 0,
    "azimuth": 0
}

computation_inputs = {
    "corenum": 1,
    "Verbose": False
}

coord_inputs = {
    "coordsystem": "GEO",
    "inputcoord": "GDZ"
}

data_retrieval_inputs = {
    "serverdata": "OFF",
    "livedata": "OFF"
}

common_inputs = dict(
    datetime_params=date_inputs,
    magfield_params=magfield_inputs,
    solar_wind_params=solar_wind_inputs,
    geomagnetic_params=geomagnetic_inputs,
    tsyganenko_params=tsyganenko_inputs,
    integration_params=integration_inputs,
    computation_params=computation_inputs,
    coordinate_params=coord_inputs,
    data_retrieval_params=data_retrieval_inputs
)


# ============================================================
# OTSO runs
# ============================================================

def location_args():
    if custom_location is not None:
        return {"Stations": [], "customlocations": [custom_location]}, custom_location[0]
    return {"Stations": [station]}, station


def location_coords():
    if custom_location is not None:
        return float(custom_location[1]), float(custom_location[2])
    import OTSO
    stations = pd.read_csv(os.path.join(os.path.dirname(OTSO.__file__), "_core", "custom_classes", "StationList.csv"))
    row = stations[stations["Name"] == station].iloc[0]
    return float(row["Latitude"]), float(row["Longitude"])


def run_cone():
    args, name = location_args()
    results = cone(**args, rigidity_params=rigidity_inputs, particle_params=particle_inputs, **common_inputs)
    cone_df, cutoff_df = results[0], results[1]
    rigidity = pd.to_numeric(cone_df["R [GV]"], errors="coerce").to_numpy()
    column = next(c for c in cone_df.columns if str(c).strip().upper() == name.upper())
    filt = cone_df[column].astype(str).str.split(";").str[0].astype(float).astype(int).to_numpy()
    keep = np.isfinite(rigidity)
    order = np.argsort(-rigidity[keep])
    rows = {str(i).strip().upper(): i for i in cutoff_df.index}
    cutoff_column = next(c for c in cutoff_df.columns if str(c).strip().upper() == name.upper())
    Ru, Rc, Rl = (float(cutoff_df.loc[rows[k], cutoff_column]) for k in ("RU", "RC", "RL"))
    return rigidity[keep][order], filt[keep][order] == 1, Ru, Rc, Rl


def choose_rigidities(rigidity, allowed, Ru, Rl):
    step = rigidity_inputs["rigiditystep"]
    inside = (rigidity <= Ru + step / 2) & (rigidity >= Rl - step / 2)
    r_in, a_in = rigidity[inside], allowed[inside]
    breaks = np.nonzero(np.diff(a_in.astype(int)))[0] + 1
    bands = np.split(np.arange(len(r_in)), breaks)
    if len(bands) > shown_penumbra:
        bands = [bands[i] for i in np.unique(np.linspace(0, len(bands) - 1, shown_penumbra).round().astype(int))]
    penumbra = [r_in[b[len(b) // 2]] for b in bands]
    width = max(Ru - Rl, 0.1 * Ru, 5 * step)
    top = min(Ru + width, rigidity_inputs["startrigidity"])
    bottom = max(Rl - width, step)
    above = np.linspace(top, Ru + width / (shown_above + 1), shown_above)
    below = np.linspace(Rl - width / (shown_below + 1), bottom, shown_below)
    chosen = np.concatenate([above, penumbra, below])
    return np.round(chosen / step) * step, top, bottom


def run_trajectories(rigidities):
    args, _ = location_args()
    paths, allowed = [], []
    for i, R in enumerate(rigidities):
        params = dict(particle_inputs, rigidity=float(R))
        result = trajectory(**args, particle_params=params, **common_inputs)[0][0]
        df = result["trajectory"]
        p = df[["GEO_X [Re]", "GEO_Y [Re]", "GEO_Z [Re]"]].to_numpy(float)
        p = p[np.isfinite(p).all(axis=1)]
        paths.append(p[:: max(1, len(p) // 3000)])
        allowed.append(int(result["Filter"]) == 1)
        print(f"trajectory {i + 1}/{len(rigidities)}: R = {R:.3f} GV, Filter {int(result['Filter'])}", flush=True)
    return paths, np.array(allowed)


# ============================================================
# Drawing utilities
# ============================================================

def camera_dir(elev, azim):
    e, a = np.radians(elev), np.radians(azim)
    return np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])


def hide(points, view):
    along = points @ view
    perp = np.linalg.norm(points - along[:, None] * view, axis=1)
    out = points.copy()
    out[(along < 0) & (perp < 1.0)] = np.nan
    return out


def on_sphere(lat, lon, r=1.003):
    la, lo = np.radians(lat), np.radians(lon)
    return r * np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def surface_lines():
    try:
        import cartopy.feature as cf
        segments = []
        for geom in cf.COASTLINE.with_scale("110m").geometries():
            for line in getattr(geom, "geoms", [geom]):
                lon, lat = np.array(line.coords).T
                segments += [on_sphere(lat, lon), np.full((1, 3), np.nan)]
    except ImportError:
        segments = []
        for lat in range(-60, 61, 30):
            segments += [on_sphere(lat, np.linspace(-180, 180, 181)), np.full((1, 3), np.nan)]
        for lon in range(-180, 180, 30):
            segments += [on_sphere(np.linspace(-90, 90, 91), lon), np.full((1, 3), np.nan)]
    return np.vstack(segments)


def add_logo(fig, x=0.065, y=0.04, h=0.15):
    if not os.path.exists(logo_file):
        fig.text(x, y + 0.045, "Generated with", color="0.85", fontsize=15, ha="center", va="bottom")
        fig.text(x, y, "OTSO", color="white", fontsize=26, weight="bold", ha="center", va="bottom")
        return
    img = mpimg.imread(logo_file)
    rows, cols = np.nonzero(img[..., 3] > 0.05)
    img = img[rows.min():rows.max() + 1, cols.min():cols.max() + 1]
    fig_w, fig_h = fig.get_size_inches()
    w = h * img.shape[1] / img.shape[0] * fig_h / fig_w
    ax = fig.add_axes([x - w / 2, y + 0.035, w, h])
    ax.imshow(img, interpolation="lanczos")
    ax.set_axis_off()
    fig.text(x, y, "OTSO", color="white", fontsize=26, weight="bold", ha="center", va="bottom")


def video_writer():
    encoders = subprocess.run(["ffmpeg", "-hide_banner", "-encoders"], capture_output=True, text=True).stdout
    if "libx264" in encoders:
        return FFMpegWriter(fps=fps, codec="libx264", extra_args=["-crf", "16", "-pix_fmt", "yuv420p"])
    return FFMpegWriter(fps=fps, codec="libopenh264", extra_args=["-g", "15", "-b:v", "40M", "-rc_mode", "quality",
                                                                  "-allow_skip_frames", "0", "-pix_fmt", "yuv420p"])


# ============================================================
# Video
# ============================================================

def render(name, scan_R, scan_allowed, Ru, Rc, Rl, top, bottom, rigidities, paths, allowed):
    draw_s, hold_s, fade_s, end_hold = 2.0, 0.5, 1.2, 5.0
    r_view = 11.0
    colours = {True: "#3ec1ff", False: "#ff8c42"}

    paths = [np.where(np.linalg.norm(p, axis=1)[:, None] > r_view, np.nan, p) for p in paths]
    site_lat, site_lon = location_coords()
    site = on_sphere(site_lat, site_lon, 1.004)
    t0 = 1.0 + seconds_per_trajectory * np.arange(len(rigidities))
    life = draw_s + hold_s + fade_s
    T = t0[-1] + life + end_hold
    t_end = t0[-1] + draw_s + hold_s
    elev = float(np.clip(0.4 * site_lat, -30, 30))
    azim_a, azim_b = site_lon - 70, site_lon + 30
    lines_xyz = surface_lines()

    in_bar = (scan_R <= top + 1e-9) & (scan_R >= bottom - 1e-9)
    bar_R, bar_allowed = scan_R[in_bar], scan_allowed[in_bar]
    reveal = np.interp(-bar_R, -rigidities, t0 + draw_s)

    fig = plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor="black")
    ax = fig.add_axes([0, 0.12, 1, 0.88], projection="3d", computed_zorder=False)
    ax.set_facecolor("black")
    ax.set_proj_type("ortho")
    ax.set_xlim(-4, 4); ax.set_ylim(-4, 4); ax.set_zlim(-4, 4)
    ax.set_box_aspect((1, 1, 1), zoom=1.3)
    ax.set_axis_off()
    u, v = np.mgrid[0:2 * np.pi:96j, 0:np.pi:48j]
    ax.plot_surface(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v), color="#163e66", shade=True,
                    linewidth=0, antialiased=False, zorder=1)
    surface, = ax.plot([], [], [], color="#cfe3f5", lw=0.8, alpha=0.9, zorder=2)
    site_marker, = ax.plot([], [], [], "o", ms=9, mfc="#ffe14d", mec="black", zorder=6)
    clip = TransformedBbox(Bbox([[0, 0.2], [1, 1]]), fig.transFigure)
    lines, heads = [], []
    for _ in range(int(np.ceil(life / seconds_per_trajectory)) + 1):
        ln, = ax.plot([], [], [], lw=1.6, zorder=4)
        hd, = ax.plot([], [], [], "o", ms=5, zorder=5)
        ln.set_clip_box(clip); hd.set_clip_box(clip)
        lines.append(ln); heads.append(hd)

    box = dict(facecolor="black", alpha=0.65, edgecolor="none", pad=3)
    fig.text(0.03, 0.935, f"{name} ({abs(site_lat):.2f}°{'N' if site_lat >= 0 else 'S'}, "
             f"{abs(site_lon):.2f}°{'E' if site_lon >= 0 else 'W'})", color="white", fontsize=24, weight="bold", bbox=box)
    fig.text(0.03, 0.9, f"vertical protons backtraced at decreasing rigidity  ·  {magfield_inputs['internalmag']} + "
             f"{magfield_inputs['externalmag']}  ·  {date_inputs['year']}-{date_inputs['month']:02d}-{date_inputs['day']:02d}",
             color="0.8", fontsize=15, bbox=box)
    fig.text(0.03, 0.855, "●  allowed (escapes the magnetosphere)", color=colours[True], fontsize=15, bbox=box)
    fig.text(0.03, 0.825, "●  forbidden (returns to Earth / trapped)", color=colours[False], fontsize=15, bbox=box)
    fig.text(0.03, 0.795, f"●  {name}", color="#ffe14d", fontsize=15, bbox=box)
    if show_otso_branding:
        add_logo(fig)
    counter = fig.text(0.97, 0.92, "", color="white", fontsize=40, ha="right", family="monospace", weight="bold", bbox=box)
    regime = fig.text(0.97, 0.875, "", color="0.85", fontsize=18, ha="right", bbox=box)

    bx = fig.add_axes([0.14, 0.075, 0.72, 0.045])
    bx.set_facecolor("#111111")
    bx.set_ylim(0, 1)
    bx.set_yticks([])
    for s in bx.spines.values():
        s.set_color("0.5")
    bx.tick_params(colors="0.85", labelsize=13)
    bx.set_xlabel("rigidity (GV)", color="0.85", fontsize=14)
    step = rigidity_inputs["rigiditystep"]
    bx.set_xlim(bar_R[0] + step / 2, bar_R[-1] - step / 2)
    rects = []
    for R, a in zip(bar_R, bar_allowed):
        rect = Rectangle((R - step / 2, 0), step, 1, color=colours[bool(a)], visible=False)
        bx.add_patch(rect)
        rects.append(rect)
    marker, = bx.plot([], [], "v", color="white", ms=12, clip_on=False, transform=bx.get_xaxis_transform())
    cut_artists = []
    for R, label, ha in ((Ru, "R$_U$", "right"), (Rc, "R$_C$", "center"), (Rl, "R$_L$", "left")):
        cut_artists.append(bx.axvline(R, color="white", lw=2, ls="--" if label == "R$_C$" else "-", visible=False))
        cut_artists.append(bx.text(R, 1.35, f"{label} = {R:.3f} GV", color="white", fontsize=14, ha=ha,
                                   transform=bx.get_xaxis_transform(), visible=False))
    summary = fig.text(0.5, 0.165, "", color="white", fontsize=17, ha="center")

    def update(frame):
        t = frame / fps
        azim = azim_a + (azim_b - azim_a) * t / T
        ax.view_init(elev=elev, azim=azim)
        view = camera_dir(elev, azim)
        c = hide(lines_xyz, view)
        c[(lines_xyz @ view) < 0] = np.nan
        surface.set_data_3d(c[:, 0], c[:, 1], c[:, 2])
        s = site[None, :] if site @ view > 0 else np.full((1, 3), np.nan)
        site_marker.set_data_3d(s[:, 0], s[:, 1], s[:, 2])
        active = sorted((t - t0[k], k) for k in range(len(rigidities)) if 0 <= t - t0[k] < life)
        for slot, (ln, hd) in enumerate(zip(lines, heads)):
            if slot >= len(active):
                ln.set_data_3d([], [], []); hd.set_data_3d([], [], [])
                continue
            age, k = active[slot]
            p = paths[k]
            frac = min(1.0, age / draw_s)
            m = max(2, int(frac * len(p)))
            alpha = 1.0 if age < draw_s + hold_s else max(0.0, 1 - (age - draw_s - hold_s) / fade_s)
            q = hide(p[:m], view)
            ln.set_data_3d(q[:, 0], q[:, 1], q[:, 2])
            ln.set_color(colours[bool(allowed[k])]); ln.set_alpha(alpha)
            h = hide(p[m - 1:m], view) if frac < 1 else np.full((1, 3), np.nan)
            hd.set_data_3d(h[:, 0], h[:, 1], h[:, 2])
            hd.set_color("white"); hd.set_alpha(alpha)
        for rect, r in zip(rects, reveal):
            rect.set_visible(t >= r)
        launched = np.nonzero(t0 <= t)[0]
        if len(launched):
            k = launched[-1]
            R = rigidities[k]
            marker.set_data([R], [1.25])
            counter.set_text(f"R = {R:.3f} GV")
            counter.set_color(colours[bool(allowed[k])] if t - t0[k] >= draw_s else "white")
            regime.set_text("above the penumbra: all allowed" if R > Ru else
                            "inside the penumbra: allowed and forbidden bands" if R >= Rl else
                            "below the penumbra: all forbidden")
        else:
            counter.set_text(""); regime.set_text("")
        finished = t >= t_end
        for a in cut_artists:
            a.set_visible(finished)
        if finished:
            marker.set_data([], [])
            summary.set_text(f"penumbra {Rl:.3f}–{Ru:.3f} GV  ·  effective cutoff R$_C$ = {Rc:.3f} GV")
        else:
            summary.set_text("")
        return []

    if preview_time is not None:
        update(int(preview_time * fps))
        fig.savefig(os.path.splitext(output_file)[0] + f"_{preview_time:g}s.png", facecolor="black")
        plt.close(fig)
        return
    anim = FuncAnimation(fig, update, frames=int(T * fps), blit=False)
    anim.save(output_file, writer=video_writer(), savefig_kwargs={"facecolor": "black"})
    plt.close(fig)


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":

    if shutil.which("ffmpeg") is None:
        print("ffmpeg is required to write the video")
        exit()

    _, name = location_args()

    scan_R, scan_allowed, Ru, Rc, Rl = run_cone()
    print(f"{name}: Ru = {Ru:.3f} GV, Rc = {Rc:.3f} GV, Rl = {Rl:.3f} GV")

    rigidities, top, bottom = choose_rigidities(scan_R, scan_allowed, Ru, Rl)
    paths, allowed = run_trajectories(rigidities)

    render(name, scan_R, scan_allowed, Ru, Rc, Rl, top, bottom, rigidities, paths, allowed)
    print("written", output_file if preview_time is None else "preview frame")
