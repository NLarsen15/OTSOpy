! ************************************************************************************************************************************
! Integration.f95 - File responsible for preforming the numerical integration for the equations of motion for the CR to compute 
! the trajectory. Subroutines to compute a rough time-step value are included also.
! ************************************************************************************************************************************

! ************************************************************************************************************************************
! subroutine TimeStep:
! Subroutine that calculates the time-step size for the integration. Value is taken to be 5% of the gyration time for the CR in
! the magnetic field present at that point in space, magnetic field is assumed unform.
! Similar method is used in Smart and Shea (1981) except the time step is taken to be 1% of gyration time in their work.
!
! INPUT:
! B - Magnetic field stength [T]
!
! OUTPUT:
! Time-step value (h) is updated within the particle module
!
! REFERENCES:
! (1) Smart, D.F., Shea, M.A., 1981. Optimum Step Length Control for Cosmic-Ray Trajectory Calculations,
!     in: International Cosmic Ray Conference, p. 208.
!
! ************************************************************************************************************************************
subroutine TimeStep(Velocity, B, MaxGyroPercent, R, hOLD, h)
use SharedParameters
implicit none

real(8), intent(in) :: B(3), Velocity(3)
real(8), intent(in) :: MaxGyroPercent, R
real(8), intent(inout) :: h, hOLD

real(8), parameter :: PI = 3.1415926535897932384626433832795d0

real(8) :: Max
real(8) :: GyroRadius, Tgyro, Vmag, Bmag

Bmag = sqrt(B(1)**2 + B(2)**2 + B(3)**2) * 1.0d4
Vmag = sqrt(Velocity(1)**2 + Velocity(2)**2 + Velocity(3)**2) / 1.0d3

GyroRadius = (1.0d10/c) * (R/Bmag)

Tgyro = 2.0d0*PI*GyroRadius/Vmag

Max = MaxGyroPercent * Tgyro

IF (adaptivestep .eqv. .TRUE.) THEN
    h = 0.01d0 * Tgyro
ELSE
    h = Max
END IF

IF (HOLD == 0.0) THEN
    HOLD = h
ELSE IF (h > HOLD * 1.1) THEN
    h = HOLD*1.1
END IF

IF (h > Max) THEN
    h = Max
END IF

HOLD = h

end subroutine TimeStep

! ************************************************************************************************************************************
! subroutine FirstTimeStep:
! Subroutine that calculates the initial time-step size for the begining of the integration by calling the TimeStep subroutine.
!
! If FixedStep > 0 and adaptivestep is off, that literal value (in seconds) is used directly as h, bypassing the
! gyropercent-derived TimeStep calculation entirely. This lets a fixed step be pinned to an exact value (e.g. to
! match a reference implementation's dt) rather than only being reachable indirectly via gyropercent.
!
! INPUT:
! Values from particle module
! FixedStep - fixed step size in seconds; <= 0 means "not set", fall back to gyropercent-derived step
!
! OUTPUT:
! Time-step value (h) is updated within the particle module
!
! ************************************************************************************************************************************
subroutine FirstTimeStep(PositionArray, VelocityArray, &
 secondTotal, MaxGyroPercent, &
 R, h, hOLD, firsth, FixedStep)
use SharedParameters
implicit none

real(8), intent(in) :: VelocityArray(2,3), secondTotal
real(8), intent(in) :: PositionArray(3,3)
real(8), intent(inout) :: h, hOLD, firsth
real(8), intent(in) :: MaxGyroPercent, R
real(8), intent(in) :: FixedStep
real(8) :: Bfield(3), Bmag
real(8) :: GEOPosition(3), GSMPosition(3)

if ((.not. adaptivestep) .and. (FixedStep > 0.0d0)) then
    h     = FixedStep
    hOLD  = FixedStep
    firsth = FixedStep
    return
end if

GEOPosition = PositionArray(2,:)
GSMPosition = PositionArray(3,:)

if (model(1) == 4) then
call MagneticField(GEOPosition, secondTotal, Bfield)
else
call MagneticField(GSMPosition, secondTotal, Bfield)
end if

hOLD = 0.0

call TimeStep(VelocityArray(1,:), Bfield, MaxGyroPercent, R, hOLD, h)

Bmag = ((Bfield(1)**2.0 + Bfield(2)**2.0 + Bfield(3)**2.0))**(0.5)
if (Bmag == 0) then
    h = 10**(-6)
end if

firsth = h

end subroutine FirstTimeStep

! ************************************************************************************************************************************
! subroutine NewMax:
! Subroutine that calls the functions in order to calculate the maximum time-step value.
!
! INPUT:
! Values from particle module
!
! OUTPUT:
! Max - Maximum time-step value [s]
!
! ************************************************************************************************************************************
subroutine NewMax(VelocityArray, PositionArray, &
MaxGyroPercent, secondTotal, R, Max, BfieldOut)
USE SharedParameters
implicit none

real(8), intent(in) :: R, secondTotal, MaxGyroPercent
real(8), intent(in) :: VelocityArray(2,3), PositionArray(3,3)
real(8), intent(out) :: Max
real(8), intent(out) :: BfieldOut(3)
real(8) :: Bfield(3), Bmag
real(8) :: GEOPosition(3), GSMPosition(3), Velocity(3)

GEOPosition = PositionArray(2,:)
GSMPosition = PositionArray(3,:)

Velocity = VelocityArray(1,:)

if (model(1) == 4 .or. model(1) == 1 .or. model(1) == 5) then
call MagneticField(GEOPosition, secondTotal, Bfield)
else
call MagneticField(GSMPosition, secondTotal, Bfield)
end if

BfieldOut = Bfield

call TimeStepMax(Bfield, Velocity, R, MaxGyroPercent, Max)

Bmag = ((Bfield(1)**2.0 + Bfield(2)**2.0 + Bfield(3)**2.0))**(0.5)
if (Bmag == 0) then
Max = 10**(-4)
end if

if (Max < 10E-7) then
Max = 10E-7
end if

end subroutine NewMax

! ************************************************************************************************************************************
! subroutine TimeStepMax:
! Subroutine that calculates the maximum time-step size for the integration for a given magnetic field strength.
! The max value is taken to be 5% of the gyration time for the CR in the magnetic field present at that point in space, 
! magnetic field is assumed unform. This allows the time-step to grow as the CR enters weaker areas of magnetic field strength.
! Smart and Shea (1981) method is used except the max time step is taken to be 1.5% of gyration time in their work.
! Changing the value of the MaxMulti variable will alter the size of the maximum time-step and will impact the accuracy of the 
! computation accordingly.
!
! INPUT:
! B - Magnetic field stength [T]
!
! OUTPUT:
! Time-step value (h) is updated within the particle module
!
! REFERENCES:
! (1) Smart, D.F., Shea, M.A., 1981. Optimum Step Length Control for Cosmic-Ray Trajectory Calculations,
!     in: International Cosmic Ray Conference, p. 208.
!
! ************************************************************************************************************************************
subroutine TimeStepMax(B, Velocity, R, MaxGyroPercent, Max)
use SharedParameters
implicit none
real(8), intent(in) :: B(3), Velocity(3), R
real(8), intent(in) :: MaxGyroPercent
real(8), intent(out) :: Max

real(8), parameter :: PI = 3.1415926535897932384626433832795d0

real(8) :: GyroRadius, Tgyro, Vmag, Bmag

! See TimeStep for the derivation of these formulas.

Bmag = sqrt(B(1)**2 + B(2)**2 + B(3)**2) * 1.0d4
Vmag = sqrt(Velocity(1)**2 + Velocity(2)**2 + Velocity(3)**2) / 1.0d3

GyroRadius = (1.0d10/c) * (R/Bmag)
Tgyro = 2.0d0*PI*GyroRadius/Vmag

Max = MaxGyroPercent * Tgyro

end subroutine TimeStepMax

subroutine mindistexamine(X, XGDZ, mindistcheck)
    USE SharedParameters
    implicit none
    real(8), intent(in) :: X(3), XGDZ(3)
    logical, intent(inout) :: mindistcheck
    real(8) :: radial

    radial = (X(1)**2.0 + X(2)**2.0 + X(3)**2.0)**(0.5)

    if (radial >= mintrapdist) then
        mindistcheck = .true.
    !    MDP(1) = XGDZ(1)
    !    MDP(2) = XGDZ(2)
    !    MDP(3) = XGDZ(3)
    end if

end subroutine mindistexamine

! ************************************************************************************************************************************
! subroutine FieldJumpCheck:
! Subroutine to check for sudden increases in magnetic field strenght, currently disabled due to its impact on computation time.
!
! INPUT:
! Bnorm0     - |B| at the start of the step [T]
! PositionRe - Candidate new position [Re, GSM/GEO as appropriate for the caller]
! t          - Time at the candidate new position [s]
!
! OUTPUT:
! h        - Timestep, scaled down in place if the step is rejected
! rejected - .true. if the field grew by more than FieldJumpLimit and the step should be retried with the new h
!
! ************************************************************************************************************************************
subroutine FieldJumpCheck(Bnorm0, PositionRe, t, h, rejected)
    implicit none

    real(8), intent(in)    :: Bnorm0
    real(8), intent(in)    :: PositionRe(3)
    real(8), intent(in)    :: t
    real(8), intent(inout) :: h
    logical, intent(out)   :: rejected

    real(8), parameter :: FieldJumpLimit = 1.10d0   ! reject if |B| grows by more than 10%

    real(8) :: Bfield(3), BnormEnd
    real(8) :: h_before

    rejected = .false.
    h_before = h

    GOTO 10

    if (Bnorm0 <= tiny(1.0d0)) then
        !print *, "FieldJumpCheck: Bnorm0 ~ 0, skipping check"
        return
    end if

    call MagneticField(PositionRe, t, Bfield)
    BnormEnd = sqrt(dot_product(Bfield,Bfield))

    !print *, "FieldJumpCheck: Bnorm0 =", Bnorm0, " BnormEnd =", BnormEnd, &
    !         " ratio =", BnormEnd/Bnorm0

    if (BnormEnd > FieldJumpLimit*Bnorm0) then
        h = h * (Bnorm0/BnormEnd)
        rejected = .true.
        !print *, "FieldJumpCheck: REJECTED -- h ", h_before, " -> ", h
    end if

    10 return

end subroutine FieldJumpCheck

subroutine OldVariables(VelocityArray, PositionArray, secondTotal, &
OLDPositionArray, OLDVelocityArray, OLDsecondTotal)
USE SharedParameters
implicit none
real(8), intent(in) :: VelocityArray(2,3), PositionArray(3,3), secondTotal
real(8), intent(inout) :: OLDPositionArray(3,3), OLDVelocityArray(2,3)
real(8), intent(inout) :: OLDsecondTotal

OLDPositionArray = PositionArray

OLDVelocityArray = VelocityArray

OLDsecondTotal = secondTotal

end subroutine OldVariables

! ************************************************************************************************************************************
! subroutine BorisRotate:
! The exact-tan(theta/2) relativistic Buneman-Boris velocity rotation.
! ************************************************************************************************************************************
subroutine BorisRotate(v0, Bfield, Bnorm, M, Q, qhalf, gamma0, Vnew)

USE SharedParameters
implicit none

real(8), intent(in)  :: v0(3), Bfield(3), Bnorm, M, Q, qhalf, gamma0
real(8), intent(out) :: Vnew(3)

real(8) :: u_minus(3)
real(8) :: TT, Tvec(3)
real(8) :: Uvec(3), UU, YY, S
real(8) :: tt_new, Tvec2(3), s_scalar
real(8) :: crossed1(3), crossed2(3)
real(8) :: dotUT, gammaNew

u_minus = gamma0 * v0

TT   = gamma0 * tan(qhalf * Bnorm / gamma0)
Tvec = TT * (Bfield / Bnorm)

call VecCross(v0, Tvec, crossed1)
Uvec = crossed1 + u_minus

dotUT = dot_product(Uvec, Tvec)
UU    = dotUT*dotUT / (c*c)
YY    = sqrt(1.0d0 + dot_product(Uvec,Uvec)/(c*c))

S = YY*YY - TT*TT

gammaNew = sqrt( 0.5d0*( S + sqrt(S*S + 4.0d0*(TT*TT + UU)) ) )

tt_new = tan(qhalf * Bnorm / gammaNew)
Tvec2  = tt_new * (Bfield / Bnorm)
s_scalar = 1.0d0 / (1.0d0 + tt_new*tt_new)

call VecCross(Uvec, Tvec2, crossed2)

Vnew = (s_scalar/gammaNew) * ( Uvec + Tvec2*dot_product(Uvec,Tvec2) + crossed2 )

end subroutine BorisRotate
