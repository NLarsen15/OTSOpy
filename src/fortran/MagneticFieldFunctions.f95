! ************************************************************************************************************************************
! MagneticFieldFunctions.f95 - Module file containing pointers to different magnetospheric magnetic field models. This module
! assigns and external and internal model to appropriate pointers to be repeatedly used throughout the computations. This avoids 
! reapeating slower IF ELSE statements over the thousands of itterations needed for the compuations.
!
! For information on each of the magnetic field models used in this file please refer to the respective files within the
! libary folder.
!
! ************************************************************************************************************************************
! subroutine MagneticField:
! Subroutine that computes the external and internal magnetic field strengths and combines them to obtain a total magnetic field
! strength at any given point within the magnetosphere.
!
! INPUT:
! X1 - Position of the CR [GDZ coordinates]
!
! OUTPUT:
! BfieldFinal - Magnetic field strength [T]
!
! ************************************************************************************************************************************
module MagneticFieldFunctions
USE Particle
USE GEOPACK1
USE GEOPACK2
USE SolarWind
implicit none

procedure (funcInternal), pointer :: InternalMagPointer => null ()
procedure (funcExternal), pointer :: ExternalMagPointer => null ()

abstract interface
function funcInternal(DUMMY1)
   real(8) :: funcInternal(3)
   real(8), intent (in) :: DUMMY1(3)
end function funcInternal
function funcExternal(DUMMY2)
    real(8) :: funcExternal(3)
    real(8), intent (in) :: DUMMY2(3)
 end function funcExternal
end interface

 
 contains

  function functionNoInt(x) ! Internal IGRF model
    real(8) :: functionNoInt(3), INTERNALGSM(3)
    real(8), intent (in) :: x(3)

    INTERNALGSM(1) = 0
    INTERNALGSM(2) = 0
    INTERNALGSM(3) = 0

    functionNoInt = INTERNALGSM
    return
  end function functionNoInt

  function functionIGRF(x) ! Internal IGRF model
    real(8) :: functionIGRF(3), INTERNALGSW(3), INTERNALGSM(3)
    real(8), intent (in) :: x(3)
  
    call IGRF_GSW_08(x(1), x(2), x(3), INTERNALGSW(1), INTERNALGSW(2), INTERNALGSW(3))
    call GSWGSM_08(INTERNALGSM(1), INTERNALGSM(2), INTERNALGSM(3), INTERNALGSW(1), INTERNALGSW(2), INTERNALGSW(3), 1)
    functionIGRF = INTERNALGSM
    return
  end function functionIGRF

  function functionDIP(x) ! Internal dipole model
    real(8) :: functionDIP(3), INTERNALGSW(3), INTERNALGSM(3)
    real(8), intent (in) :: x(3)
  
    call DIP_08(x(1), x(2), x(3), INTERNALGSW(1), INTERNALGSW(2), INTERNALGSW(3))
    call GSWGSM_08(INTERNALGSM(1), INTERNALGSM(2), INTERNALGSM(3), INTERNALGSW(1), INTERNALGSW(2), INTERNALGSW(3), 1)

    functionDIP = INTERNALGSM
    return
  end function functionDIP

  function functionCustom(x) ! Custom Spherical Haormincs model
    real(8) :: functionCustom(3), INTERNALSPH(3), INTERNALGSW(3), INTERNALGEO(3) 
    real(8) :: INTERNALGSM(3), GEOPosition(3), GEOSPHPosition(3)
    real(8) :: R, Theta, Phi
    real(8), intent (in) :: x(3)

    call CoordinateTransform("GSM", "GEO", year, day, secondTotal, x, GEOPosition)

    call GAUSSCUSTOM(GEOPosition(1), GEOPosition(2), GEOPosition(3), INTERNALGEO(1), INTERNALGEO(2), INTERNALGEO(3))

    call GEOGSW_08(INTERNALGEO(1),INTERNALGEO(2),INTERNALGEO(3),INTERNALGSW(1),INTERNALGSW(2),INTERNALGSW(3),1)
    call GSWGSM_08(INTERNALGSM(1), INTERNALGSM(2), INTERNALGSM(3), INTERNALGSW(1), INTERNALGSW(2), INTERNALGSW(3), 1)

    functionCustom = INTERNALGSM
    return
  end function functionCustom

  function functionCustomNonStandard(x) ! Custom Spherical Haormincs model
    real(8) :: functionCustomNonStandard(3), INTERNALSPH(3), INTERNALGSW(3), INTERNALGEO(3) 
    real(8) :: INTERNALGSM(3), GEOPosition(3), GEOSPHPosition(3)
    real(8) :: R, Theta, Phi
    real(8), intent (in) :: x(3)
  
    call GAUSSCUSTOM(x(1), x(2), x(3), INTERNALGEO(1), INTERNALGEO(2), INTERNALGEO(3))

    functionCustomNonStandard = INTERNALGEO
    return
  end function functionCustomNonStandard

  function functionNoEx(x) !No external field
    real(8) :: functionNoEx(3), TSYGSM(3), TSYGSM1(3)
    real(8), intent (in) :: x(3)
  
    TSYGSM(1) = 0
    TSYGSM(2) = 0
    TSYGSM(3) = 0

    functionNoEx = TSYGSM
  
    return
  end function functionNoEx

 function function87S(x) ! Tsyganenko 1987 short
   real(8) :: function87S(3), TSYGSM(3), TSYGSM1(3)
   real(8), intent (in) :: x(3)
 
   call TSY87S(IOPT, x(1), x(2), x(3), TSYGSM1(1), TSYGSM1(2), TSYGSM1(3))
   if (model(1) == 4) THEN
   call CoordinateTransformVec("GSM", "GEO", year, day, secondTotal, TSYGSM1, TSYGSM)
   else
   TSYGSM=TSYGSM1
   end if
   function87S = TSYGSM
 
   return
 end function function87S
 
  function function87L(x) ! Tsyganenko 1987 long
    real(8) :: function87L(3), TSYGSM(3), TSYGSM1(3)
    real(8), intent (in) :: x(3)
  
    call TSY87L(IOPT, x(1), x(2), x(3), TSYGSM1(1), TSYGSM1(2), TSYGSM1(3))
    if (model(1) == 4) THEN
    call CoordinateTransformVec("GSM", "GEO", year, day, secondTotal, TSYGSM1, TSYGSM)
    else
    TSYGSM=TSYGSM1
    end if
    function87L = TSYGSM
  
    return
  end function function87L

  function function89(x) ! Tsyganenko 1989
    real(8) :: function89(3), TSYGSM(3), TSYGSM1(3)
    real(8), intent (in) :: x(3)
  
    call T89D_DP(IOPT, parmod, PSI, x(1), x(2), x(3), TSYGSM1(1), TSYGSM1(2), TSYGSM1(3))
    if (model(1) == 4) THEN
    call CoordinateTransformVec("GSM", "GEO", year, day, secondTotal, TSYGSM1, TSYGSM)
    else
    TSYGSM=TSYGSM1
    end if
    function89 = TSYGSM
  
    return
  end function function89

  function function89Boberg(x) ! Tsyganenko 1989 Boberg Extension
    real(8) :: function89Boberg(3), TSYGSM(3), TSYGSM1(3)
    real(8), intent (in) :: x(3)
  
    call T89DBoberg(IOPT, parmod, PSI, DSTBob, KpIndex, x(1), x(2), x(3), TSYGSM1(1), TSYGSM1(2), TSYGSM1(3))

    if (model(1) == 4) THEN
    call CoordinateTransformVec("GSM", "GEO", year, day, secondTotal, TSYGSM1, TSYGSM)
    else
    TSYGSM=TSYGSM1
    end if

    function89Boberg = TSYGSM
  
    return
  end function function89Boberg

  function function96(x) ! Tsyganenko 1996
    real(8) :: function96(3), TSYGSM(3), TSYGSM1(3)
    real(8), intent (in) :: x(3)
    real(8), dimension(10) :: parmod2
    real(8) :: GSMx(3), PSItemp, TSYfield(3)

    parmod2 = real(parmod,4)
    GSMx(1) = real(x(1),4) 
    GSMx(2) = real(x(2),4)
    GSMx(3) = real(x(3),4)
    PSItemp = real(PSI,8)
  
    call T96_01(IOPT, parmod2, PSItemp, GSMx(1), GSMx(2), GSMx(3), TSYfield(1), TSYfield(2), TSYfield(3))
    TSYGSM1(1) = TSYfield(1)
    TSYGSM1(2) = TSYfield(2)
    TSYGSM1(3) = TSYfield(3)
    if (model(1) == 4) THEN
    call CoordinateTransformVec("GSM", "GEO", year, day, secondTotal, TSYGSM1, TSYGSM)
    else
    TSYGSM=TSYGSM1
    end if
    function96 = TSYGSM
  
    return
  end function function96

  function function01(x) ! Tsyganenko 2001
    real(8) :: function01(3), TSYGSM(3), TSYGSM1(3)
    real(8), intent (in) :: x(3)
    real, dimension(10) :: parmod2
    real :: GSMx(3), PSItemp, TSYfield(3)

    parmod2 = real(parmod,4)
    GSMx(1) = real(x(1),4) 
    GSMx(2) = real(x(2),4)
    GSMx(3) = real(x(3),4)
    PSItemp = real(PSI,8)

    call T01_01(IOPT, parmod2, GSMx(1), GSMx(2), GSMx(3), TSYfield(1), TSYfield(2), TSYfield(3))
    TSYGSM1(1) = TSYfield(1)
    TSYGSM1(2) = TSYfield(2)
    TSYGSM1(3) = TSYfield(3)
    if (model(1) == 4) THEN
    call CoordinateTransformVec("GSM", "GEO", year, day, secondTotal, TSYGSM1, TSYGSM)
    else
    TSYGSM=TSYGSM1
    end if
    function01 = TSYGSM
    
    return
  end function function01

  function function01S(x) ! Tsyganenko 2001 storm-time variation
    real(8) :: function01S(3), TSYGSM(3), TSYGSM1(3)
    real(8), intent (in) :: x(3)
  
    call T01_S(parmod, x(1), x(2), x(3), TSYGSM1(1), TSYGSM1(2), TSYGSM1(3))
    if (model(1) == 4) THEN
    call CoordinateTransformVec("GSM", "GEO", year, day, secondTotal, TSYGSM1, TSYGSM)
    else
    TSYGSM=TSYGSM1
    end if
    function01S = TSYGSM
  
    return
  end function function01S

  function function04(x) ! Tsyganenko 2004 
    real(8) :: function04(3), TSYGSM(3), TSYGSM1(3)
    real(8), intent (in) :: x(3)
  
    call T04_S(parmod, x(1), x(2), x(3), TSYGSM1(1), TSYGSM1(2), TSYGSM1(3))
    if (model(1) == 4) THEN
    call CoordinateTransformVec("GSM", "GEO", year, day, secondTotal, TSYGSM1, TSYGSM)
    else
    TSYGSM=TSYGSM1
    end if
    function04 = TSYGSM
  
    return
  end function function04


  function function15N(x) ! Tsyganenko 2015 N-index 
    real(8) :: function15N(3), TSYGSM(3), TSYGSM1(3)
    real(8), intent (in) :: x(3)
  
    call TA_2015_N(0, parmod, PSI, x(1), x(2), x(3), TSYGSM1(1), TSYGSM1(2), TSYGSM1(3))
    
    if (model(1) == 4) THEN
    call CoordinateTransformVec("GSM", "GEO", year, day, secondTotal, TSYGSM1, TSYGSM)
    else
    TSYGSM=TSYGSM1
    end if
    function15N = TSYGSM
  
    return
  end function function15N


  function function15B(x) ! Tsyganenko 2015 B-index 
    real(8) :: function15B(3), TSYGSM(3), TSYGSM1(3)
    real(8), intent (in) :: x(3)
  
    call TA_2015_B(0, parmod, PSI, x(1), x(2), x(3), TSYGSM1(1), TSYGSM1(2), TSYGSM1(3))
    
    if (model(1) == 4) THEN
    call CoordinateTransformVec("GSM", "GEO", year, day, secondTotal, TSYGSM1, TSYGSM)
    else
    TSYGSM=TSYGSM1
    end if
    function15B = TSYGSM
  
    return
  end function function15B

  function function16(x) ! Tsyganenko 2015 B-index 
    real(8) :: function16(3), TSYGSM(3), TSYGSM1(3)
    real(8), intent (in) :: x(3)
  
    call RBF_MODEL_2016(0, parmod, PSI, x(1), x(2), x(3), TSYGSM1(1), TSYGSM1(2), TSYGSM1(3))
    
    if (model(1) == 4) THEN
    call CoordinateTransformVec("GSM", "GEO", year, day, secondTotal, TSYGSM1, TSYGSM)
    else
    TSYGSM=TSYGSM1
    end if
    function16 = TSYGSM
  
    return
  end function function16

  function functionMHD(x) ! MHD 
    real(8) :: functionMHD(3), MHDexternal(3), PositionMAG
    real(8), intent (in) :: x(3)
  
    call MHDField(x, MHDexternal)
    
    functionMHD = MHDexternal
  
    PositionMAG = ((x(1)**2)+(x(2)**2)+(x(3)**2))*(1/2)
  
    IF (ISNAN(functionMHD(1))) THEN
      functionMHD(1) = 0.0
    END IF
    IF (ISNAN(functionMHD(2))) THEN
      functionMHD(2) = 0.0
    END IF
    IF (ISNAN(functionMHD(3))) THEN
      functionMHD(3) = 0.0
    END IF
    
    return
  end function functionMHD

! ************************************************************************************************************************************
! subroutine MagneticFieldAssign:
! Subroutine that assigns the functions for specific magnetic field models to an internal and external pointer. To be used within
! the MagneticField subroutine (MagneticField.f95).
!
! INPUT:
! mode - integer array of length 2 containg information on the models to be used. (e.g. [1,3] = IGRF and Tsyganenko1989)
!
! OUTPUT:
! InternalMagPointer and ExternalMagPointer are assigned appropriate magnetic field models to be used.
!
! ************************************************************************************************************************************
  subroutine MagneticFieldAssign(mode)
  implicit none
  integer(8) :: mode(2)

  call RECALC_08(year, day, hour, minute, secondINT, SW(1), SW(2), SW(3))
  
  IF (mode(1) == 0) THEN
    InternalMagPointer => functionNoInt ! NO INTERNAL FIELD
  ELSE IF (mode(1) == 1) THEN
    InternalMagPointer => functionIGRF ! IGRF
  ELSE IF (mode(1) == 2) THEN
    InternalMagPointer => functionDIP  ! DIPOLE
  ELSE IF (mode(1) == 3) THEN
    InternalMagPointer => functionCustom ! Custom Gauss
  ELSE IF (mode(1) == 4) THEN
    InternalMagPointer => functionCustomNonStandard ! Custom Gauss Non-Standard
  ELSE
    print *, "Please enter valid internal magnetic field model"
  END IF

  IF (mode(2) == 0) THEN
    ExternalMagPointer => functionNoEx ! NO EXTERNAL FIELD
  ELSE IF (mode(2) == 1) THEN
    ExternalMagPointer => function87S  ! TSYGANENKO 87 SHORT
  ELSE IF (mode(2) == 2) THEN
    ExternalMagPointer => function87L  ! TSYGANENKO 87 LONG
  ELSE IF (mode(2) == 3) THEN
    ExternalMagPointer => function89   ! TSYGANENKO 89
  ELSE IF (mode(2) == 4) THEN
    ExternalMagPointer => function96   ! TSYGANENKO 96
  ELSE IF (mode(2) == 5) THEN
    ExternalMagPointer => function01   ! TSYGANENKO 01
  ELSE IF (mode(2) == 6) THEN
    ExternalMagPointer => function01S  ! TSYGANENKO 01 STORM
  ELSE IF (mode(2) == 7) THEN
    ExternalMagPointer => function04  ! TSYGANENKO 04
  ELSE IF (mode(2) == 8) THEN
    ExternalMagPointer => function89Boberg ! TSYGANENKO 89 BOBERG EXTENSION
  ELSE IF (mode(2) == 9) THEN
    ExternalMagPointer => function15N  ! TSYGANENKO 15 N-index
  ELSE IF (mode(2) == 10) THEN
    ExternalMagPointer => function15B  ! TSYGANENKO 15 B-index
  ELSE IF (mode(2) == 11) THEN
    ExternalMagPointer => function16  ! TSYGANENKO 16 RBF
  ELSE IF (mode(2) == 99) THEN
    ExternalMagPointer => functionMHD  ! MHD
  ELSE
    print *, "Please enter valid external magnetic field model"
  END IF


    
  end subroutine MagneticFieldAssign
 
end module MagneticFieldFunctions