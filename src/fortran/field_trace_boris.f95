! ************************************************************************************************************************************
! subroutine Boris_FieldTrace_Advanced:
! Advanced Boris integration scheme for field line tracing with improved accuracy and stability.
! Uses the Boris leapfrog method adapted for field line following.
!
! INPUT:
! Bsign - Direction sign for field line tracing (+1 or -1)
!
! OUTPUT:
! Bfield - Magnetic field vector at the new position [nT]
! Updates Position in the Particle module
!
! ************************************************************************************************************************************
subroutine Boris_FieldTrace_Advanced(PositionArray, &
    secondTotal, DistanceTraveled, Bsign, Bfield)

    use SharedParameters
    implicit none

    real(8), intent(inout) :: PositionArray(3,3)
    real(8), intent(inout) :: secondTotal
    real(8), intent(inout) :: DistanceTraveled

    real(8), intent(in) :: Bsign
    real(8), intent(out) :: Bfield(3)

    real(8) :: x0(3), xhalf(3), x1(3)
    real(8) :: Xnew(3), xGSM(3)
    real(8) :: x0GDZ(3), xhalfGDZ(3), x1GDZ(3), XnewGDZ(3)
    real(8) :: x0GSM(3), xhalfGSM(3), x1GSM(3), XnewGSM(3)
    real(8) :: Bfield0(3), Bfield0Mag, Bfield0unit(3)
    real(8) :: BfieldHalf(3), BfieldHalfMag, BfieldHalfunit(3)
    real(8) :: Bfield1(3), Bfield1Mag, Bfield1unit(3)
    real(8) :: v0(3), vhalf(3), v1(3)  ! Direction vectors (not actual velocities)
    real(8) :: dt, dt_half
    real(8) :: speed_norm ! Normalized speed for field line following

    !print *, "========================================="
    !print *, "Boris_FieldTrace_Advanced: Starting advanced field line trace"
    !print *, "========================================="

    xGSM = PositionArray(3,:)

    !print *, "Initial position (GSM):", xGSM(1), xGSM(2), xGSM(3)
    !print *, "Field direction sign (Bsign):", Bsign

    ! Step size and time step for Boris method
    dt = 10000.0  ! meters (same as other methods for consistency)
    dt_half = dt / 2.0
    speed_norm = 1.0  ! Normalized speed for field line following

    !print *, "Time step dt:", dt, "meters"

    if (model(1) == 4) then
        xGSM = PositionArray(2,:)
    end if

    call MagneticField(xGSM, secondTotal, Bfield0)

    !print *, "Initial B-field (GSM):", Bfield0(1), Bfield0(2), Bfield0(3)

    ! Convert position to meters for calculations
    xGSM(1) = xGSM(1)*Re_m
    xGSM(2) = xGSM(2)*Re_m
    xGSM(3) = xGSM(3)*Re_m

    !print *, "Position in GSM (meters):", xGSM(1), xGSM(2), xGSM(3)

    ! Calculate initial field direction
    Bfield0Mag = ((Bfield0(1)**2.0 + Bfield0(2)**2.0 + Bfield0(3)**2.0))**(0.5)
    Bfield0unit(1) = Bfield0(1)/Bfield0Mag
    Bfield0unit(2) = Bfield0(2)/Bfield0Mag
    Bfield0unit(3) = Bfield0(3)/Bfield0Mag

    !print *, "B-field magnitude:", Bfield0Mag
    !print *, "B-field unit vector:", Bfield0unit(1), Bfield0unit(2), Bfield0unit(3)

    ! Initialize velocity (direction vector) from field
    v0(1) = Bsign * Bfield0unit(1) * speed_norm
    v0(2) = Bsign * Bfield0unit(2) * speed_norm
    v0(3) = Bsign * Bfield0unit(3) * speed_norm

    ! Boris Step 1: Advance position by half time step
    !print *, "-----------------------------------------"
    !print *, "Boris Advanced Step 1: Half position advance"
    xhalf(1) = xGSM(1) + v0(1) * dt_half
    xhalf(2) = xGSM(2) + v0(2) * dt_half
    xhalf(3) = xGSM(3) + v0(3) * dt_half

    !print *, "Half position (GSM meters):", xhalf(1), xhalf(2), xhalf(3)

    ! Transform to get field at half position
    xhalfGSM(1) = xhalf(1)/Re_m
    xhalfGSM(2) = xhalf(2)/Re_m
    xhalfGSM(3) = xhalf(3)/Re_m

    call MagneticField(xhalfGSM, secondTotal, BfieldHalf)

    BfieldHalfMag = ((BfieldHalf(1)**2.0 + BfieldHalf(2)**2.0 + BfieldHalf(3)**2.0))**(0.5)
    BfieldHalfunit(1) = BfieldHalf(1)/BfieldHalfMag
    BfieldHalfunit(2) = BfieldHalf(2)/BfieldHalfMag
    BfieldHalfunit(3) = BfieldHalf(3)/BfieldHalfMag

    !print *, "B-field at half position:", BfieldHalf(1), BfieldHalf(2), BfieldHalf(3)
    !print *, "B-field magnitude at half position:", BfieldHalfMag

    ! Boris Step 2: Update velocity using field at half position
    !print *, "-----------------------------------------"
    !print *, "Boris Advanced Step 2: Velocity update with midpoint field"
    vhalf(1) = Bsign * BfieldHalfunit(1) * speed_norm
    vhalf(2) = Bsign * BfieldHalfunit(2) * speed_norm
    vhalf(3) = Bsign * BfieldHalfunit(3) * speed_norm

    ! Boris Step 3: Complete the step with updated velocity
    !print *, "-----------------------------------------"
    !print *, "Boris Advanced Step 3: Complete position advance"
    Xnew(1) = xhalf(1) + vhalf(1) * dt_half
    Xnew(2) = xhalf(2) + vhalf(2) * dt_half
    Xnew(3) = xhalf(3) + vhalf(3) * dt_half

    !print *, "New position (GSM meters):", Xnew(1), Xnew(2), Xnew(3)

    ! Convert back to coordinate system units
    XnewGSM(1) = Xnew(1)/Re_m
    XnewGSM(2) = Xnew(2)/Re_m
    XnewGSM(3) = Xnew(3)/Re_m

    !print *, "New position (GSM):", XnewGSM(1), XnewGSM(2), XnewGSM(3)

    call CoordinateTransform("GSM", "GDZ", year, day, secondTotal, XnewGSM, XnewGDZ)

    if (model(1) == 4) then
        call CoordinateTransform("GEO", "GDZ", year, day, secondTotal, XnewGSM, XnewGDZ)
    end if

    !print *, "Final position (GDZ):", XnewGDZ(1), XnewGDZ(2), XnewGDZ(3)

    ! Get magnetic field at final position
    call MagneticField(XnewGSM, secondTotal, Bfield)

    !print *, "Final B-field at new position:", Bfield(1), Bfield(2), Bfield(3)

    ! Update position in particle
    PositionArray(1,:) = XnewGDZ
    PositionArray(2,:) = XnewGSM
    PositionArray(3,:) = XnewGSM

    DistanceTraveled = DistanceTraveled + dt

    !print *, "Updated distance traveled:", DistanceTraveled, "meters"
    !print *, "========================================="
    !print *, "Boris_FieldTrace_Advanced: Complete"
    !print *, "========================================="

    end subroutine Boris_FieldTrace_Advanced
