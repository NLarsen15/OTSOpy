"""Release checks that need no compiled code, run by the Test Suite workflow (and runnable locally).

    python tests/ci_checks.py

1. Every supported Python version has a prebuilt library for every platform in OTSO/_core/libs,
   built for the right CPU architecture (macOS: Apple Silicon only).
2. The version in setup.py and CITATION.cff agree.

Exits non-zero if anything is missing or wrong.
"""
import os
import re
import struct
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIBS = ROOT / "OTSO" / "_core" / "libs"

ELF_MACHINE = {0x3E: "x86_64", 0xB7: "aarch64"}
MACHO_CPU = {0x01000007: "x86_64", 0x0100000C: "arm64"}
PE_MACHINE = {0x8664: "x86_64", 0xAA64: "arm64"}


def supported_versions():
    text = (ROOT / "setup.py").read_text(encoding="utf-8-sig")
    m = re.search(r"python_requires\s*=\s*['\"]>=\s*3\.(\d+)\s*,\s*<\s*3\.(\d+)['\"]", text)
    lo, hi = int(m.group(1)), int(m.group(2))
    return [f"3{minor}" for minor in range(lo, hi)]


def arches(path):
    data = path.read_bytes()
    if data[:4] == b"\x7fELF":
        return {ELF_MACHINE.get(struct.unpack_from("<H", data, 18)[0], "?")}
    if data[:4] == b"\xca\xfe\xba\xbe":  # Mach-O universal (fat) binary, big-endian header
        n = struct.unpack_from(">I", data, 4)[0]
        return {MACHO_CPU.get(struct.unpack_from(">I", data, 8 + 20 * i)[0], "?") for i in range(n)}
    if data[:4] == b"\xcf\xfa\xed\xfe":  # thin 64-bit Mach-O, little-endian
        return {MACHO_CPU.get(struct.unpack_from("<I", data, 4)[0], "?")}
    if data[:2] == b"MZ":
        pe = struct.unpack_from("<I", data, 0x3C)[0]
        return {PE_MACHINE.get(struct.unpack_from("<H", data, pe + 4)[0], "?")}
    return {"?"}


def check_libraries(versions):
    expected = {  # file name pattern -> required architectures
        "_MiddleMan.cpython-{v}-x86_64-linux-gnu.so": {"x86_64"},
        "_MiddleMan.cpython-{v}-aarch64-linux-gnu.so": {"aarch64"},
        "_MiddleMan.cpython-{v}-darwin.so": {"arm64"},
        "_MiddleMan.cp{v}-win_amd64.pyd": {"x86_64"},
    }
    problems, rows = [], []
    for v in versions:
        row = [f"3.{v[1:]}"]
        for pattern, want in expected.items():
            path = LIBS / pattern.format(v=v)
            if not path.exists():
                row.append("missing")
                problems.append(f"missing {path.name}")
                continue
            have = arches(path)
            if not want <= have:
                row.append(f"{'+'.join(sorted(have))} (need {'+'.join(sorted(want))})")
                problems.append(f"{path.name} is {'+'.join(sorted(have))}, needs {'+'.join(sorted(want))}")
            else:
                row.append("ok " + "+".join(sorted(have)))
        rows.append(row)

    header = ["Python", "Linux x86_64", "Linux aarch64", "macOS arm64", "Windows x86_64"]
    table = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    table += ["| " + " | ".join(r) + " |" for r in rows]
    return problems, "\n".join(table)


def check_version():
    setup_v = re.search(r"version\s*=\s*['\"]([^'\"]+)['\"]",
                        (ROOT / "setup.py").read_text(encoding="utf-8-sig")).group(1)
    cff_v = re.search(r"^version:\s*['\"]?([^'\"\s]+)", (ROOT / "CITATION.cff").read_text(encoding="utf-8-sig"),
                      re.M).group(1)
    ok = setup_v == cff_v
    return ([] if ok else [f"setup.py version {setup_v} != CITATION.cff version {cff_v}"],
            f"setup.py {setup_v}, CITATION.cff {cff_v}")


def main():
    versions = supported_versions()
    lib_problems, table = check_libraries(versions)
    ver_problems, ver_line = check_version()
    report = f"### Prebuilt libraries\n\n{table}\n\n### Version\n\n{ver_line}\n"
    print(report)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a") as fh:
            fh.write("## Release checks\n\n" + report)
    problems = lib_problems + ver_problems
    for p in problems:
        print("ERROR:", p)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
