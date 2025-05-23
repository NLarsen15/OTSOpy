cmake_minimum_required(VERSION 3.10)
project(FortranToPython LANGUAGES Fortran)

set(CMAKE_Fortran_COMPILER "C:/rtools40/mingw64/bin/gfortran.exe")
set(CMAKE_C_COMPILER "C:/rtools40/mingw64/bin/gcc.exe")
set(CMAKE_CXX_COMPILER "C:/rtools40/mingw64/bin/g++.exe")

# Find Python and NumPy
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

# Set the Fortran source files
set(FORTRAN_SOURCES
    "MiddleMan.f95"
    "Acceleration.f95"
    "Asymptotic.f95"
    "Conversion.f95"
    "CoordTrans.for"
    "CoordinateTransforms.f95"
    "CustomGaussModule.f95"
    "Fill_IMF_gaps.for"
    "Fill_SW_gaps.for"
    "GEOPACK_DP.f"
    "GeopackModule.f95"
    "Integration.f95"
    "IntegrationFunctions.f95"
    "LorentzRelativity.f95"
    "MagneticField.f95"
    "MagneticFieldFunctions.f95"
    "MagnetopauseEscape.f95"
    "MagnetopauseFunctions.f95"
    "MagnetopauseModule.f95"
    "Particle.f95"
    "ParticleModule.f95"
    "Prepare_input_4.for"
    "Prepare_intervals_1.for"
    "Rigidity.f95"
    "SolarWindModule.f95"
    "T96.f"
    "TSY89DP.f"
    "Tsyg_01.for"
    "Tsyganenko04.f"
    "Tsyganenko87l.f"
    "Tsyganenko87s.f"
    "VectorCalc.f95"
    "Velocity.f95"
    "date_util.f"
    "dipole.f"
    "heliospheric_transformation.f"
    "igrf_coef.f"
    "init_nouveau.f"
    "onera_desp_lib.f"
    "t01_s.f"
    "TSY89_BobergExtension.f95"
    "MHDModule.f95"
    "MHDBfield.f95"
    "GridModule.f95"
)

# Generate the .pyf file using F2PY
set(PYF_FILE "${CMAKE_CURRENT_BINARY_DIR}/MiddleMan.pyf")

add_custom_command(
    OUTPUT "${PYF_FILE}"
    COMMAND ${Python_EXECUTABLE} -m numpy.f2py -m MiddleMan -h "${PYF_FILE}" --overwrite-signature ${FORTRAN_SOURCES}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    DEPENDS ${FORTRAN_SOURCES}
)

add_custom_target(GeneratePYF ALL DEPENDS "${PYF_FILE}")

# Build the Fortran shared library
add_library(MiddleMan SHARED ${FORTRAN_SOURCES})

# Ensure FortranLib is compiled before the Python module
add_dependencies(MiddleMan GeneratePYF)

# Compile Fortran into a Python extension module using F2PY
set(PYTHON_MODULE "${CMAKE_CURRENT_BINARY_DIR}/MiddleMan.pyd")

add_custom_command(
    OUTPUT "${PYTHON_MODULE}"
    COMMAND ${Python_EXECUTABLE} -m numpy.f2py --backend=meson -c "${PYF_FILE}" ${FORTRAN_SOURCES}
            -L"${CMAKE_CURRENT_BINARY_DIR}"  MiddleMan.f95 -m MiddleMan
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    DEPENDS "${PYF_FILE}" MiddleMan
)

add_custom_target(GeneratePythonModule ALL DEPENDS "${PYTHON_MODULE}")
