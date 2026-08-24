module SharedParameters
    implicit none

    ! Physical constants
    real(8), parameter      :: q_0 = 1.602176634e-19        ! Elementary charge (C)
    real(8), parameter      :: mp_0 = 1.67262192595e-27     ! Proton mass (kg)
    real(8), parameter      :: mm_0 = 1.883531627e-28       ! muon mass (kg)
    real(8), parameter      :: me_0 = 9.1093837139e-31      ! electron mass (kg)
    real(8), parameter      :: c = 299792458.0              ! Speed of light (m/s)
    real(8), parameter      :: Re_km = 6371.2               ! Earth radius (km)
    real(8), parameter      :: Re_m = 6371200               ! Earth radius (m)
    real(8), parameter      :: Lasth = 1E-6                 ! Last step size (s)
    real(8), parameter      :: Joule2MeV = 6.241509074461E9 ! Conversion factor for Joules to MeV
    integer(4), parameter   :: backsavelim = 3              ! Number of steps back wards

    integer(4) :: coordtranscounter

    integer(4) :: model(4)
    integer(4) :: year
    integer(4) :: day
    integer(4) :: hour
    integer(4) :: minute
    integer(4) :: secondINT

    character(3) :: CoordINMHD
    character(3) :: CoordOUTMHD

    real(8) :: mintrapdist
    real(8) :: spheresize
    real(8) :: q_0A

    logical :: adaptivestep
    logical :: totalbetacheck
    logical :: trapdistcheck

    real(8) :: CurrentGyro

    !$omp threadprivate(CurrentGyro)

    save
end module SharedParameters