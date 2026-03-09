project('MiddleMan', 'fortran', version : '1.0')

# Simple meson script that uses f2py directly
# This creates a run target to build the extension with f2py

sources = [
  'MiddleMan.f95',
  'TSY89_refit',
  'Acceleration.f95',
  'Asymptotic.f95',
  'Conversion.f95',
  'CoordTrans.for',
  'CoordinateTransforms.f95',
  'CustomGaussModule.f95',
  'Fill_IMF_gaps.for',
  'Fill_SW_gaps.for',
  'GEOPACK_DP.f',
  'GeopackModule.f95',
  'Integration.f95',
  'IntegrationFunctions.f95',
  'LorentzRelativity.f95',
  'MagneticField.f95',
  'MagneticFieldFunctions.f95',
  'MagnetopauseEscape.f95',
  'MagnetopauseFunctions.f95',
  'MagnetopauseModule.f95',
  'Particle.f95',
  'ParticleModule.f95',
  'Prepare_input_4.for',
  'Prepare_intervals_1.for',
  'Rigidity.f95',
  'SolarWindModule.f95',
  'T96.f',
  'Tsyg_01.for',
  'Tsyganenko04.f',
  'Tsyganenko87l.f',
  'Tsyganenko87s.f',
  'VectorCalc.f95',
  'Velocity.f95',
  'date_util.f',
  'dipole.f',
  'heliospheric_transformation.f',
  'igrf_coef.f',
  'init_nouveau.f',
  'onera_desp_lib.f',
  't01_s.f',
  'TSY89a.f95',
  'TSY89c.f95',
  'MHDModule.f95',
  'MHDBfield.f95',
  'GridModule.f95',
  'TSY15_B.f',
  'TSY15_N.f',
  'TA16_RBF.f',
]

# Create a run target that builds the extension with f2py
f2py_build = run_target('build_f2py',
  command: [
    'python', '-m', 'numpy.f2py', '-c', 
    '-m', 'MiddleMan'
  ] + sources,
  depends: [],
  env: environment()
)
