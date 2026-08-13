import subprocess
import os

try:
    subprocess.run(["python", "-m", "numpy.f2py", "-c", "hello.f90", "-m", "hello", "--backend=meson"], check=True, capture_output=True)
    for root, dirs, files in os.walk("."):
        if "compile_commands.json" in files or "build.ninja" in files:
            for f in files:
                if f == "build.ninja":
                    with open(os.path.join(root, f)) as f_in:
                        for line in f_in.readlines():
                            if "-O" in line:
                                print(line.strip())
                                break
except Exception as e:
    print(e)
