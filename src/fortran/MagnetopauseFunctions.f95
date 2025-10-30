! MagnetopauseFunctions.f95 - Module file containing pointers to different magnetopauuse models. This module
! assigns a model to appropriate pointers to be repeatedly used throughout the computations. This avoids 
! reapeating slower IF ELSE statements over the thousands of itterations needed for the compuations.
!
! For information on each of the magnetopause models used in this file please refer to the respective 
! papers listed below:
!
! ************************************************************************************************************************************
module MagnetopauseFunctions
USE Particle
USE GEOPACK1
USE GEOPACK2
USE SolarWind
USE MagnetoPause
implicit none

procedure (funcPause), pointer :: PausePointer => null ()

abstract interface
function funcPause()
   integer(4) :: funcPause
end function funcPause
end interface

 
 contains

  function functionSphere() ! Sphere
    integer(4) :: functionSphere
    real(8) :: GSEPosition(3), x1, y1, z1, TestResult

    call CoordinateTransform("GDZ", "GSE", year, day, secondTotal, Position, GSEPosition)

    if (model(1) == 4) then
      call CoordinateTransform("GDZ", "GEO", year, day, secondTotal, Position, GSEPosition)
    end if
    
    TestResult = -1
    Result = 0
    x1 = GSEPosition(1)
    y1 = GSEPosition(2)
    z1 = GSEPosition(3)

    TestResult = (z1**2 + y1**2 + x1**2)**(0.5) - (spheresize)

    IF (TestResult < 0) THEN
      functionSphere = 0
    ELSE IF (TestResult >= 0) THEN
      functionSphere = 1
      IF (FinalStep == 0) THEN
        FinalStep = 1
        TestResult = -1
      END IF
    END IF
    
    return
  end function functionSphere

  function functionDisabled() ! No Magnetopause
    integer(4) :: functionDisabled

    functionDisabled = 0
    
    return
  end function functionDisabled

  function functionAberratedFormisano() ! Aberrated Formisano Model
    integer(4) :: functionAberratedFormisano
    real(8) :: a11, a22, a33, a13, a23, a34, a14, a12, a24, a44
    real(8) :: GSEPosition(3), x1, y1, z1, TestResult

    call CoordinateTransform("GDZ", "GSE", year, day, secondTotal, Position, GSEPosition)

    TestResult = -1
    Result = 0
    x1 = GSEPosition(1)
    y1 = GSEPosition(2)
    z1 = GSEPosition(3)

    a11 = 0.65
    a22 = 1
    a33 = 1.16
    a13 = -0.28
    a23 = -0.11
    a34 = -0.36
    a14 = 21.41
    a12 = 0.03
    a24 = 0.46
    a44 = -221

    TestResult = a11*(x1**2) + a22*(y1**2) + a33*(z1**2) + a12*(x1*y1) + a13*x1*z1 + a23*y1*z1 + a14*x1 + a24*y1 + a34*z1 + (a44)
    IF (x1 < -60) THEN
        TestResult = 1
        IF (FinalStep == 0) THEN
          FinalStep = 1
          TestResult = -1
        END IF
    END IF

    IF (TestResult < 0) THEN
      functionAberratedFormisano = 0
    ELSE IF (TestResult >= 0) THEN
      functionAberratedFormisano = 1
    END IF

    return
  end function functionAberratedFormisano

  function functionSibeck() ! Sibeck Model
    integer(4) :: functionSibeck
    real(8) :: a11, a22, a33, a14
    real(8) :: GSMPosition(3), x1, y1, z1, TestResult, p, Rvalue

    call CoordinateTransform("GDZ", "GSM", year, day, secondTotal, Position, GSMPosition)
    
    TestResult = -1
    Result = 0
    x1 = GSMPosition(1)
    y1 = GSMPosition(2)
    z1 = GSMPosition(3)

    p = PDYN * 10.0**(9.0)
    IF (p >= 0.54 .AND. p < 0.87) THEN
        a11 = 0.19
        a22 = 19.3
        a33 = -272.4
        a14 =  0.705
    ELSE IF (p >= 0.87 .AND. p < 1.47) THEN
        a11 = 0.19
        a22 = 18.7
        a33 = -243.9
        a14 = 1.17
    ELSE IF (p >= 1.47 .AND. p < 2.60) THEN
        a11 = 0.14
        a22 = 18.2
        a33 = -217.2
        a14 = 2.04
    ELSE IF (p >= 2.60 .AND. p < 4.90) THEN
        a11 = 0.15
        a22 = 17.3
        a33 = -187.4
        a14 = 3.75
    ELSE IF (p >= 4.90 .AND. p < 9.90) THEN
        a11 = 0.18
        a22 = 14.2
        a33 = -139.2
        a14 = 7.4
    ELSE IF (p < 0.54 .OR. p > 9.90)  THEN
        a11 = 0.14
        a22 = 18.2
        a33 = -217.2
        a14 = 2.04
        IF (PDYN == 0) THEN
            p = a14
        END IF
    END IF

    Rvalue = (z1**2 + y1**2)**(0.5)
    TestResult = Rvalue**2 + (a11)*(x1**2) + a22*x1*((a14/p)**(1.0/6.0)) + (a33)*((a14/p)**(1.0/3.0))

    IF (x1 < -50) THEN
        TestResult = 1
        IF (FinalStep == 0) THEN
          FinalStep = 1
          TestResult = -1
        END IF
    END IF

    IF (TestResult < 0) THEN
      functionSibeck = 0
    ELSE IF (TestResult >= 0) THEN
      functionSibeck = 1
    END IF

    return
  end function functionSibeck

  function functionKobel() ! Kobel Model
    integer(4) :: functionKobel
    real(8) :: Ak, Bk(7), Fk(7), rho2, kpar, sink, cosk
    real(8) :: x1rot, y1rot, z1rot, rhorot
    real(8) :: GSMPosition(3), x1, y1, z1, TestResult, DIP

    TestResult = -1
    Result = 0
    dip = PSI

    Ak = -0.0545
    Bk(1) = 11.7
    Bk(2) = 11.1
    Bk(3) = 10.8
    Bk(4) = 10.4
    Bk(5) = 10.4
    Bk(6) = 10.2
    Bk(7) = 10.2

    Fk(1) = 20.0
    Fk(2) = 15.0
    Fk(3) = 6.67
    Fk(4) = 10.0
    Fk(5) = 5.0
    Fk(6) = 6.0
    Fk(7) = 6.0

    call CoordinateTransform("GDZ", "GSM", year, day, secondTotal, Position, GSMPosition)
    x1 = GSMPosition(1)
    y1 = GSMPosition(2)
    z1 = GSMPosition(3)

    rho2 = y1*y1 + z1*z1

    if (rho2 > 900) THEN
        TestResult = 1
        IF (FinalStep == 0) THEN
          FinalStep = 1
          TestResult = -1
        END IF
    end if

    if (x1 < -60) THEN
        TestResult = 1
        IF (FinalStep == 0) THEN
          FinalStep = 1
          TestResult = -1
        END IF
    end if

    kpar = DIP/Fk(IOPT)
    sink = SIN(kpar)
    cosk = COS(kpar)

    x1rot = x1*cosk + z1*sink
    y1rot = y1
    z1rot = (-x1)*sink + z1*cosk

    rhorot = y1rot*y1rot + z1rot*z1rot

    

    IF (x1rot > Ak*rhorot + Bk(IOPT)) THEN
        TestResult = 1
        IF (FinalStep == 0) THEN
          FinalStep = 1
          TestResult = -1
        END IF
    END IF

    IF (TestResult < 0) THEN
      functionKobel = 0
    ELSE IF (TestResult >= 0) THEN
      functionKobel = 1
    END IF

    return
  end function functionKobel

  function functionLin() !Lin et al 2010 Model
  IMPLICIT NONE

  integer(4) :: functionLin
  real(8) :: GSEPosition(3)

  DOUBLE PRECISION :: X, Y, Z, r_pos, theta, phi, r_mp, Bz
  DOUBLE PRECISION :: r0, Psum, tilt
  DOUBLE PRECISION :: b0, b1, b2, b3, b, inner, f_theta_phi
  DOUBLE PRECISION :: cn, cs, dn, ds, en, es, theta_n, theta_s
  DOUBLE PRECISION :: phi_n, phi_s, cos_yn, cos_ys, y_n, y_s, Qn, Qs, Q
  DOUBLE PRECISION :: pi, TestResult, Result
  
  ! Lin et al. (2010) coefficients
  DOUBLE PRECISION, PARAMETER :: a0=12.544, a1=-0.194, a2=0.305, a3=0.0573, a4=2.178
  DOUBLE PRECISION, PARAMETER :: a5=0.0571, a6=-0.999, a7=16.473, a8=0.00152, a9=0.382
  DOUBLE PRECISION, PARAMETER :: a10=0.0431, a11=-0.00763, a12=-0.210, a13=0.0405
  DOUBLE PRECISION, PARAMETER :: a14=-4.430, a15=-0.636, a16=-2.600, a17=0.832
  DOUBLE PRECISION, PARAMETER :: a18=-5.328, a19=1.103, a20=-0.907, a21=1.450

  pi = 4.0d0*ATAN(1.0d0)

  TestResult = -1
  Result = 0

  call CoordinateTransform("GDZ", "GSM", year, day, secondTotal, Position, GSEPosition)

  X = GSEPosition(1)
  Y = GSEPosition(2)
  Z = GSEPosition(3)

  Bz = IMF(3)
  r_pos = SQRT(X**2 + Y**2 + Z**2)
  IF (r_pos == 0.0d0) THEN
     !PRINT *, 'Error: Position cannot be the origin.'
     STOP
  END IF

  theta = ACOS(X / r_pos)            ! Polar angle from +X axis
  phi = ATAN2(Y, Z)
  IF (phi .LT. 0.0d0) phi = phi + 2.0d0*pi
  tilt  = PSI                        ! Convert degrees to radians

  Psum  = (Pdyn*10**9) + Magpressure


  ! --- Subsolar stand-off distance r0 ---
  r0 = a0 * Psum**a1 * (1.0d0 + a2*(EXP(a3*Bz)-1.0d0)/(EXP(a4*Bz)+1.0d0))

  ! --- b parameter ---
  b0 = a6 + a7*(EXP(a8*Bz)-1.0d0)/(EXP(a9*Bz)+1.0d0)
  b1 = a10
  b2 = a11 + a12*tilt
  b3 = a13
  b  = b0 + b1*COS(phi) + b2*SIN(phi) + b3*SIN(phi)**2


  inner = COS(theta/2.0d0) + a5*SIN(2.0d0*theta)*(1.0d0-EXP(-theta))
  IF (inner < 1.0d-1) inner = 1.0d-1   ! clamp to 0.1 (avoid huge tail values)
  f_theta_phi = inner**b


  ! --- Q terms for north/south ---
  cn = a14*Psum**a15
  cs = cn
  dn = a16 + a17*tilt + a18*tilt**2
  ds = a16 - a17*tilt + a18*tilt**2
  en = a21
  es = a21
  theta_n = a19 + a20*tilt
  theta_s = a19 - a20*tilt
  phi_n = pi/2.0d0
  phi_s = 3.0d0*pi/2.0d0

  cos_yn = COS(theta)*COS(theta_n) + SIN(theta)*SIN(theta_n)*COS(phi - phi_n)
  cos_ys = COS(theta)*COS(theta_s) + SIN(theta)*SIN(theta_s)*COS(phi - phi_s)
  y_n = ACOS(MIN(MAX(cos_yn,-1.0d0),1.0d0))
  y_s = ACOS(MIN(MAX(cos_ys,-1.0d0),1.0d0))

  Qn = cn*EXP(dn*(y_n**en))
  Qs = cs*EXP(ds*(y_s**es))
  Q  = Qn + Qs

  ! --- Final magnetopause radius ---
  r_mp = r0*f_theta_phi + Q

  ! --- Debug prints ---
  !PRINT *, 'Lin Model Debug:'
  !PRINT *, '  Position (GSE): X =', X, ' Y =', Y, ' Z =', Z
  !PRINT *, '  Distance from origin:', r_pos
  !PRINT *, '  Magnetopause radius:', r_mp

  ! --- Check if position is outside ---
  if (r_pos > r_mp) THEN
     TestResult = 1
     !PRINT *, '  Status: OUTSIDE magnetosphere (r_pos > r_mp)'
     IF (FinalStep == 0) THEN
         FinalStep = 1
         TestResult = -1
     END IF
  ELSE
     TestResult = -1
     !PRINT *, '  Status: INSIDE magnetosphere (r_pos <= r_mp)'
  END IF

  ! --- Additional check for particles beyond -60 Re in X ---
  IF (X < -60.0d0) THEN
     TestResult = 1
     !PRINT *, '  Status: OUTSIDE magnetosphere (X < -60 Re)'
     IF (FinalStep == 0) THEN
         FinalStep = 1
         TestResult = -1
     END IF
  END IF

  IF (TestResult < 0) THEN
      functionLin = 0
      !PRINT *, '  Final result: INSIDE (returning 0)'
  ELSE IF (TestResult >= 0) THEN
      functionLin = 1
      !PRINT *, '  Final result: OUTSIDE (returning 1)'
  END IF

  return
  end function functionLin
 
  function functionTSY() ! Magnetopause models used within the Tsyganenko models
    integer(4) :: functionTSY
    real(8) :: GSMPosition(3), x1, y1, z1, TestResult

    call CoordinateTransform("GDZ", "GSM", year, day, secondTotal, Position, GSMPosition)
    TestResult = -1
    Result = 0
    x1 = GSMPosition(1)
    y1 = GSMPosition(2)
    z1 = GSMPosition(3)
  
    IF (SubResult == 1) THEN
      TestResult = 1
      IF (FinalStep == 0) THEN
        FinalStep = 1
        TestResult = -1
      END IF
    END IF

    IF (x1 < -50) THEN
      TestResult = 1
      IF (FinalStep == 0) THEN
        FinalStep = 1
        TestResult = -1
      END IF
    END IF

    IF (TestResult < 0) THEN
      functionTSY = 0
    ELSE IF (TestResult >= 0) THEN
      functionTSY = 1
    END IF
  
    return
  end function functionTSY

! ************************************************************************************************************************************
! subroutine MagnetopuaseAssign:
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
  subroutine MagnetopauseAssign(Pause)
  USE Particle
  implicit none
  integer(4) :: Pause

  IF (model(2) == 4) THEN
    Pause = 4
  ELSE IF (model(2) == 5) THEN
    Pause = 5
  ELSE IF (model(2) == 6) THEN
    Pause = 6
  ELSE IF (model(2) == 7) THEN
    Pause = 7
  ELSE IF (model(2) == 9) THEN
    Pause = 8
  ELSE IF (model(2) == 10) THEN
    Pause = 9
  ELSE IF (model(2) == 11) THEN
    Pause = 10
  END IF

  IF (Pause == 0) THEN
    PausePointer => functionSphere ! 25Re Sphere
  ELSE IF (Pause == 1) THEN
    PausePointer => functionAberratedFormisano  ! Aberrated Formisano Model
  ELSE IF (Pause == 2) THEN
    PausePointer => functionSibeck  ! Sibeck Model
  ELSE IF (Pause == 3) THEN
    PausePointer => functionKobel   ! Kobel Model
  ELSE IF (Pause == 4) THEN
    PausePointer => functionTSY  ! TSYGANENKO models
  ELSE IF (Pause == 5) THEN
    PausePointer => functionTSY  ! TSYGANENKO models
  ELSE IF (Pause == 6) THEN
    PausePointer => functionTSY  ! TSYGANENKO models
  ELSE IF (Pause == 7) THEN
    PausePointer => functionTSY  ! TSYGANENKO models
  ELSE IF (Pause == 8) THEN
    PausePointer => functionLin  ! Lin et al 2010 model
  ELSE IF (Pause == 9) THEN
    PausePointer => functionLin  ! Lin et al 2010 model
  ELSE IF (Pause == 10) THEN
    PausePointer => functionLin  ! Lin et al 2010 model
  END IF

  IF (Pause == 99) THEN
    PausePointer => functionDisabled  ! No Magnetopause
  END IF


    
  end subroutine MagnetopauseAssign
 
 end module MagnetopauseFunctions