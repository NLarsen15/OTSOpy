"""Cross-platform consistency check used by the Test Suite workflow.

    python tests/ci_reference.py run OUT.json      compute a fixed set of outputs on this machine
    python tests/ci_reference.py compare DIR       compare every *.json in DIR (searched recursively)

The prebuilt Fortran libraries are compiled separately on every OS, architecture and Python
version, with different gfortran versions. They should all give the same numbers; this catches
a library that was built wrongly, or a compiler or flag that changes results.
"""
import datetime
import json
import os
import platform
import sys
from pathlib import Path

# Tolerances. Field and coordinate values should agree to rounding error. Rc can move by a
# rigidity step or two if a trajectory near the penumbra goes the other way. Ru and Rl depend on
# single allowed/forbidden trajectories inside the penumbra, so they are reported but not checked.
TOL = {"coordtrans": 1e-7, "magfield_rel": 1e-6, "cutoff": 0.05}

STATIONS = ["OULU", "ROME", "MXCO"]


def run(out_path):
    from OTSO import coordtrans, cutoff, magfield

    results = {"meta": {"label": os.environ.get("REF_LABEL", platform.platform()),
                        "python": platform.python_version(), "machine": platform.machine(),
                        "system": platform.system()}}

    dates = [datetime.datetime(1850, 6, 21, 12), datetime.datetime(2000, 10, 12, 8), datetime.datetime(2024, 5, 10, 18)]
    ct = {}
    for d in dates:
        df = coordtrans(Locations=[[1.0, 0.0, 0.0], [3.0, -2.0, 1.5]], dates=[d, d],
                        CoordIN="GEO", CoordOUT="GSM", corenum=1, Verbose=False)[0]
        ct[d.isoformat()] = df.iloc[:, -3:].to_numpy(dtype=float).tolist()
    results["coordtrans"] = ct

    mf = {}
    for ext in ["TSY89c", "TSY96", "TSY04", "TSY89_refit"]:
        df = magfield(
            Locations=[[6.0, 0.0, 0.0], [-8.0, 2.0, 1.0], [0.0, 5.0, 3.0]],
            datetime_params={"year": 2015, "month": 3, "day": 17, "hour": 12},
            magfield_params={"internalmag": "IGRF", "externalmag": ext, "boberg": ext == "TSY89_refit",
                             "bobergtype": "DST_MIDPOINT"},
            solar_wind_params={"vx": -500, "vy": 0, "vz": 0, "bx": 0, "by": 2, "bz": -5,
                               "by_avg": 2, "bz_avg": -5, "density": 5, "pdyn": 0},
            geomagnetic_params={"Dst": -50, "kp": 4.3, "n_index": 0, "b_index": 0, "sym_h_corrected": -50},
            tsyganenko_params={"G1": 5, "G2": 5, "G3": 5, "W1": 0.5, "W2": 0.5, "W3": 0.5,
                               "W4": 0.5, "W5": 0.5, "W6": 0.5},
            coordinate_params={"inputcoord": "GSM", "coordsystem": "GSM", "coordout": "GSM"},
            computation_params={"corenum": 1, "Verbose": False},
        )[0]
        mf[ext] = df.iloc[:, -3:].to_numpy(dtype=float).tolist()
    results["magfield"] = mf

    co = {}
    for ext in ["TSY89c", "TSY89_refit"]:
        df = cutoff(
            Stations=STATIONS, cutoff_comp="Vertical",
            datetime_params={"year": 2024, "month": 5, "day": 1, "hour": 0},
            magfield_params={"internalmag": "IGRF", "externalmag": ext, "boberg": ext == "TSY89_refit",
                             "bobergtype": "DST_MIDPOINT", "magnetopause": "Kobel"},
            rigidity_params={"startrigidity": 20, "endrigidity": 0, "rigiditystep": 0.01, "rigidityscan": "ON"},
            geomagnetic_params={"Dst": -20, "kp": 2.3, "n_index": 0, "b_index": 0, "sym_h_corrected": 0},
            computation_params={"corenum": 2, "Verbose": False},
            data_retrieval_params={"serverdata": "OFF", "livedata": "OFF"},
        )[0]
        co[ext] = {st: {r: float(df.loc[r, st]) for r in ("Ru", "Rc", "Rl")} for st in df.columns}
    results["cutoff"] = co

    Path(out_path).write_text(json.dumps(results, indent=1))
    print(f"wrote {out_path}")


def _flat(res):
    """Yield (kind, key, value) for every number in a results file."""
    for d, rows in res["coordtrans"].items():
        for i, row in enumerate(rows):
            for j, v in enumerate(row):
                yield "coordtrans", f"{d} p{i} c{j}", v
    for ext, rows in res["magfield"].items():
        for i, row in enumerate(rows):
            for j, v in enumerate(row):
                yield "magfield", f"{ext} p{i} B{'xyz'[j]}", v
    for ext, stations in res["cutoff"].items():
        for st, vals in stations.items():
            for r, v in vals.items():
                yield ("cutoff" if r == "Rc" else "penumbra"), f"{ext} {st} {r}", v


def compare(directory):
    files = sorted(Path(directory).rglob("*.json"))
    if not files:
        print("No result files found")
        return 1
    runs = {f: json.loads(f.read_text()) for f in files}
    # Baseline: Linux x64 if present (the configuration the developers check by hand), else the first file.
    base_file = next((f for f, r in runs.items() if r["meta"]["system"] == "Linux"
                      and r["meta"]["machine"] in ("x86_64", "AMD64")), files[0])
    base = dict(((k, n), v) for k, n, v in _flat(runs[base_file]))
    scale = {}
    for (k, n), v in base.items():
        if k == "magfield":
            p = n.rsplit(" ", 1)[0]
            scale[p] = max(scale.get(p, 0.0), abs(v))

    lines = [f"Baseline: **{runs[base_file]['meta']['label']}**", "",
             "| run | coordtrans max abs diff | magfield max rel diff | Rc max diff (GV) | worst Rc | Ru/Rl max diff (GV, not checked) | status |",
             "|---|---|---|---|---|---|---|"]
    failed = []
    for f, res in runs.items():
        worst = {"coordtrans": 0.0, "magfield": 0.0, "cutoff": 0.0, "penumbra": 0.0}
        worst_cut = ""
        for k, n, v in _flat(res):
            b = base.get((k, n))
            if b is None:
                continue
            d = abs(v - b)
            if k == "magfield":
                d /= max(scale[n.rsplit(" ", 1)[0]], 1e-12)
            if d > worst[k]:
                worst[k] = d
                if k == "cutoff":
                    worst_cut = f"{n}: {v:.2f} vs {b:.2f}"
        ok = (worst["coordtrans"] <= TOL["coordtrans"] and worst["magfield"] <= TOL["magfield_rel"]
              and worst["cutoff"] <= TOL["cutoff"])
        if not ok:
            failed.append(res["meta"]["label"])
        lines.append(f"| {res['meta']['label']} | {worst['coordtrans']:.1e} | {worst['magfield']:.1e} | "
                     f"{worst['cutoff']:.2f} | {worst_cut} | {worst['penumbra']:.2f} | {'ok' if ok else '**FAIL**'} |")

    lines += ["", f"Tolerances: coordtrans {TOL['coordtrans']}, magfield relative {TOL['magfield_rel']}, "
                  f"Rc {TOL['cutoff']} GV."]
    report = "\n".join(lines)
    print(report)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a") as fh:
            fh.write("## Cross-platform consistency\n\n" + report + "\n")
    if failed:
        print("\nOutside tolerance: " + ", ".join(failed))
        return 1
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in ("run", "compare"):
        print(__doc__)
        sys.exit(2)
    if sys.argv[1] == "run":
        run(sys.argv[2])
    else:
        sys.exit(compare(sys.argv[2]))
