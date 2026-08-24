subroutine ParticleMass(R, Atomic, &
    M, Q, Z, A, E_0, lambda)
    use SharedParameters
    implicit none
    real(8), intent(in) :: R
    integer(8), intent(in) :: Atomic
    real(8), intent(out) :: M, Q, Z, A, E_0, lambda

IF (Atomic == -1) THEN ! Muon
    A = 1.0
    Z = 1.0
    M = mm_0
    Q = -1.0 * q_0A
    E_0 = (M * (c**2)) * Joule2MeV
    lambda = (((R*Z/(E_0 * A))**2) + 1)**(0.5)

ELSE IF (Atomic == 0) THEN ! Electron
    A = 1.0
    Z = 1.0
    M = me_0
    Q = -1.0 * q_0A
    E_0 = (M * (c**2)) * Joule2MeV
    lambda = (((R*Z/(E_0 * A))**2) + 1)**(0.5)

 ELSE IF (Atomic == 1) THEN ! Hydrogen
    A = 1.0
    Z = 1.0
    M = mp_0 * A
    Q = q_0A * Z
    E_0 = (M * (c**2)) * Joule2MeV
    lambda = (((R*Z/(E_0 * A))**2) + 1)**(0.5)

 ELSE IF (Atomic == 2) THEN ! Helium
    A = 4.0
    Z = 2.0
    M = mp_0 * A
    Q = q_0A * Z
    E_0 = (M * (c**2)) * Joule2MeV
    lambda = (((R*Z/(E_0 * A))**2) + 1)**(0.5)

 ELSE IF (Atomic == 3) THEN ! Lithium
    A = 7.0
    Z = 3.0
    M = mp_0 * A
    Q = q_0A * Z
    E_0 = (M * (c**2)) * Joule2MeV
    lambda = (((R*Z/(E_0 * A))**2) + 1)**(0.5)

 ELSE IF (Atomic == 4) THEN ! Beryllium
    A = 9.0
    Z = 4.0
    M = mp_0 * A
    Q = q_0A * Z
    E_0 = (M * (c**2)) * Joule2MeV
    lambda = (((R*Z/(E_0 * A))**2) + 1)**(0.5)

 ELSE IF (Atomic > 4) THEN
    print *, "Values above Z=4 not supported yet."
    stop
 END IF

 end subroutine ParticleMass

subroutine ParticleVelocities(StartPosition, &
    Rigidity, Date, PositionArray, lambda, secondTotal, &
    inputcoord, VelocityArray)

use solarwind
use SharedParameters
implicit none

! ============================================================
! INPUT / OUTPUT
! ============================================================

real(8), intent(in)    :: StartPosition(5)
real(8), intent(in)    :: Date(6)
real(8), intent(in)    :: Rigidity
real(8), intent(in)    :: lambda

real(8), intent(inout) :: secondTotal

character(len=3), intent(in) :: inputcoord

real(8), intent(inout) :: PositionArray(3,3)
real(8), intent(inout) :: VelocityArray(2,3)

! ============================================================
! LOCAL VARIABLES
! ============================================================

real(8) :: V

real(8) :: tempposition(3)
real(8) :: temppositionGDZ(3)

real(8) :: Norm(3)
real(8) :: NormUnit(3)

real(8) :: GEOSPHposition(3)

real(8) :: GSMPosition(3)
real(8) :: GEOPosition(3)
real(8) :: Position(3)
real(8) :: TrimStartPosition(3)

real(8) :: RadialUnit(3)

real(8) :: radialMag
real(8) :: normMag

real(8) :: dotNormRadial
real(8) :: angleNormRadial

real(8) :: velocityMagNormal
real(8) :: velocityMagGEO
real(8) :: velocityMagGSM

real(8) :: velocityDifference(3)

real(8) :: VelocityNormal(3)
real(8) :: VelocityGEO(3)
real(8) :: VelocityGSM(3)

real(8) :: w

real(8), parameter :: PI = 3.1415926535897932384626433832795d0

secondTotal = real(Date(6), kind=8)

VelocityNormal = 0.0d0
VelocityGEO    = 0.0d0
VelocityGSM    = 0.0d0

Norm        = 0.0d0
NormUnit    = 0.0d0
RadialUnit  = 0.0d0

Position       = 0.0d0
GEOPosition    = 0.0d0
GSMPosition    = 0.0d0
GEOSPHposition = 0.0d0

TrimStartPosition = real(dnint(StartPosition(:3) * 1.0d8), kind=8) / 1.0d8

call CoordinateTransform(inputcoord, "GDZ", year, day, secondTotal, &
                        TrimStartPosition, Position)

call CoordinateTransform(inputcoord, "GSM", year, day, secondTotal, &
                        TrimStartPosition, GSMPosition)

call CoordinateTransform(inputcoord, "GEO", year, day, secondTotal, &
                        TrimStartPosition, GEOPosition)

PositionArray(1,:) = Position
PositionArray(2,:) = GEOPosition
PositionArray(3,:) = GSMPosition

PositionArray(1,:) = real(dnint(PositionArray(1,:) * 1.0d8), kind=8) / 1.0d8
PositionArray(2,:) = real(dnint(PositionArray(2,:) * 1.0d8), kind=8) / 1.0d8
PositionArray(3,:) = real(dnint(PositionArray(3,:) * 1.0d8), kind=8) / 1.0d8

call Rigidity2velocity(lambda, V)

tempposition(1) = StartPosition(1)
tempposition(2) = StartPosition(2)
tempposition(3) = StartPosition(3)

call NormalVector(tempposition, inputcoord, Norm, model, &
                  year, day, secondTotal)

call VelocityComponents(V, Norm, VelocityNormal)

w = sqrt( &
    VelocityNormal(1)**2 + &
    VelocityNormal(2)**2 + &
    VelocityNormal(3)**2 )

call CoordinateTransform(inputcoord, "GEO", year, day, secondTotal, &
                        tempposition, temppositionGDZ)

call CoordinateTransform("GEO", "SPH", year, day, secondTotal, &
                         temppositionGDZ, GEOSPHposition)

call AzimuthZenith2GEO( &
    w, &
    GEOSPHposition(2), &
    GEOSPHposition(3), &
    StartPosition(4), &
    StartPosition(5), &
    VelocityGEO)

if (model(1) /= 4 .and. model(1) /= 1 .and. model(1) /= 5) then

    call CoordinateTransformVec( &
        "GEO", "GSM", year, day, secondTotal, &
        VelocityGEO, VelocityGSM)

else

    VelocityGSM = VelocityGEO

end if

VelocityArray(1,:) = VelocityGSM
VelocityArray(2,:) = VelocityGEO

end subroutine ParticleVelocities

 subroutine AntiAssignCharge(Anti)
    use SharedParameters
    implicit none
    integer(8), intent(in) :: Anti

    if (Anti == 1) then
        q_0A = -1.0 * q_0
    else
        q_0A = q_0
    end if

end subroutine AntiAssignCharge