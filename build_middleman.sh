#!/usr/bin/env bash

set -euo pipefail

FORTRAN_DIR="src/fortran"
LIBS_DIR="OTSO/_core/libs"

MODULE_NAME="MiddleMan"
EXTENSION_NAME="_MiddleMan"

echo "=== Entering ${FORTRAN_DIR} ==="
cd "$FORTRAN_DIR"

echo "=== Cleaning previous build ==="

rm -f MiddleMan.py
rm -f f90wrap_*.f90
rm -f ${EXTENSION_NAME}*.so
rm -f ${EXTENSION_NAME}*.pyd
rm -rf build
rm -rf .f2py
rm -rf __pycache__
rm -rf *.egg-info

echo
echo "=== Running f90wrap ==="

f90wrap -m "${MODULE_NAME}" MiddleMan.f95

if [ ! -f "${MODULE_NAME}.py" ]; then
    echo "ERROR: ${MODULE_NAME}.py was not generated."
    exit 1
fi

echo
echo "=== Patching Python wrapper to use relative imports ==="

python3 <<'EOF'
from pathlib import Path

p = Path("MiddleMan.py")
text = p.read_text()

text = text.replace(
    "import _MiddleMan",
    "from . import _MiddleMan"
)

p.write_text(text)

print("Patched MiddleMan.py")
EOF

echo
echo "=== Building extension ==="

WRAPPERS=$(ls f90wrap*.f90 2>/dev/null)

python3 -m numpy.f2py \
-c \
--backend=meson \
--opt="-O3" \
--f90flags="-O3 -fopenmp -fautomatic -frecursive -fno-fast-math -ffp-contract=off" \
--f77flags="-O3 -fopenmp -fautomatic -frecursive -fno-fast-math -ffp-contract=off" \
MiddleMan.f95 \
SolarWindModule.f95 \
SharedParameters.f95 \
GEOPACK_DP.f \
GeopackModule.f95 \
CustomGaussModule.f95 \
old_commonsModule.f95 \
Particle.f95 \
Conversion.f95 \
Velocity.f95 \
CoordinateTransforms.f95 \
CoordTrans.for \
onera_desp_lib.f \
init_nouveau.f \
heliospheric_transformation.f \
date_util.f \
igrf_coef.f \
MagneticFieldFunctions.f95 \
Tsyganenko87l.f \
Tsyganenko87s.f \
TSY89_refit.f95 \
TSY89a.f95 \
TSY89c.f95 \
MagneticField.f95 \
MagnetopauseFunctions.f95 \
MagnetopauseModule.f95 \
IntegrationFunctions.f95 \
VectorCalc.f95 \
integration_utils.f95 \
LorentzRelativity.f95 \
Asymptotic.f95 \
Termination_checks.f95 \
Acceleration.f95 \
TSYmodules.f95 \
T96.f \
Tsyg_01.for \
t01_s.f \
Tsyganenko04.f \
TSY15_N.f \
TSY15_B.f \
TA16_RBF.f \
MHDModule.f95 \
MHDBfield.f95 \
betachecking.f95 \
igrf14.f \
4RK.f95 \
5RK.f95 \
6RK.f95 \
Boris_Buneman.f95 \
Vay.f95 \
HC.f95 \
field_trace_boris.f95 \
${WRAPPERS} \
-m "${EXTENSION_NAME}" \
-lgomp

echo
echo "=== Finding generated extension ==="

SO_FILE=$(find . -maxdepth 1 \( -name "${EXTENSION_NAME}*.so" -o -name "${EXTENSION_NAME}*.pyd" \) | head -n 1)

if [ -z "$SO_FILE" ]; then
    echo "ERROR: Extension module was not created."
    exit 1
fi

cd - >/dev/null

mkdir -p "$LIBS_DIR"

echo
echo "=== Installing ==="

mv "${FORTRAN_DIR}/MiddleMan.py" "$LIBS_DIR/"
mv "${FORTRAN_DIR}/${SO_FILE#./}" "$LIBS_DIR/"

echo
echo "Installed:"
ls -l "$LIBS_DIR"

echo
echo "Done."