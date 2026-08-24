Subroutine Termination_checks(End, PositionArray, Escaped, DistanceTraveled, Steps, &
 TimeElapsed, termtype)
USE SolarWind
USE Magnetopause
USE SharedParameters
implicit none

logical :: Particle_Finished
logical, intent(in) :: Escaped
real(8), intent(in) :: DistanceTraveled
integer(8), intent(in) :: Steps
real(8), intent(in) :: End(4), PositionArray(3,3), TimeElapsed
integer(4), intent(inout) :: termtype
Real(8) :: Position(3)

    Particle_Finished = .false.

    Position = PositionArray(1,:)

    ! check for particle escape
    IF (Escaped) THEN
        CALL outside_magnetopause(Particle_Finished, termtype)
        if (Particle_Finished) return
    END IF

    ! Check for earth encounter
    CALL EarthEncounter(End, Position, Particle_Finished, termtype)
    if (Particle_Finished) return

    ! check for max distance traveled
    IF (End(2) .ne. 0) THEN
        CALL MaxDistance(End, DistanceTraveled, Particle_Finished, termtype)
        if (Particle_Finished) return
    END IF

    ! check for max steps
    IF (End(4) .ne. 0) THEN
        CALL MaxSteps(End, Steps, Particle_Finished, termtype)
        if (Particle_Finished) return
    END IF
    
    IF (End(3) .ne. 0) THEN
        CALL MaxTime(End, TimeElapsed, Particle_Finished, termtype)
        if (Particle_Finished) return
    END IF

end subroutine Termination_checks

subroutine EarthEncounter(End, Position, Particle_Finished, termtype)
    USE SharedParameters
    logical, intent(inout) :: Particle_Finished
    real(8), intent(in) :: Position(3)
    real(8), intent(in) :: End(4)
    integer(4), intent(inout) :: termtype

    IF (Position(1) < End(1)) THEN
        termtype = -1
        Particle_Finished = .true.
    END IF
end subroutine EarthEncounter


subroutine MaxDistance(End, DistanceTraveled, Particle_Finished, termtype)
    USE SharedParameters
    logical, intent(inout) :: Particle_Finished
    real(8), intent(in) :: DistanceTraveled
    real(8), intent(in) :: End(4)
    integer(4), intent(inout) :: termtype

    IF (DistanceTraveled/Re_m > End(2)) THEN
            exceededistance = DistanceTraveled/Re_m - End(2)
            termtype = -2
            Particle_Finished = .true.
    END IF

end subroutine MaxDistance

subroutine MaxSteps(End, Steps, Particle_Finished, termtype)
    USE SharedParameters
    logical, intent(inout) :: Particle_Finished
    integer(8), intent(in) :: Steps
    real(8), intent(in) :: End(4)
    integer(4), intent(inout) :: termtype

    IF (REAL(Steps) >= End(4)) THEN
        termtype = -3
        Particle_Finished = .true.
    END IF

end subroutine MaxSteps

subroutine MaxTime(End, TimeElapsed, Particle_Finished, termtype)
    USE SharedParameters
    logical, intent(inout) :: Particle_Finished
    real(8), intent(in) :: TimeElapsed
    real(8), intent(in) :: End(4)
    integer(4), intent(inout) :: termtype

    IF (TimeElapsed > End(3)) THEN
        termtype = -4
        Particle_Finished = .true.
    END IF

end subroutine MaxTime

subroutine outside_magnetopause(Particle_Finished, termtype)
    USE SharedParameters
    logical, intent(inout) :: Particle_Finished
    integer(4), intent(inout) :: termtype

    termtype = 1
    Particle_Finished = .true.

end subroutine outside_magnetopause

subroutine FinalStepCheck(Finalstep, PositionArray, VelocityArray, &
OldPositionArray, OldVelocityArray, termtype, secondTotal, OldsecondTotal)

implicit none

real(8), intent(inout) :: PositionArray(3,3), OldPositionArray(3,3)
real(8), intent(inout) :: VelocityArray(2,3), OldVelocityArray(2,3)
real(8), intent(inout) :: secondTotal, OldsecondTotal
logical, intent(inout) :: Finalstep
integer(4), intent(inout) :: termtype

IF (FinalStep .neqv. .true.) THEN
    !print *, "Final Step Triggered"
    !print *, Result
    !print *, FinalStep
    FinalStep = .true.
    
    IF (termtype .ne. 1) THEN
    PositionArray(1,:) = OldPositionArray(1,:)
    PositionArray(2,:) = OldPositionArray(2,:)
    PositionArray(3,:) = OldPositionArray(3,:)

    VelocityArray(1,:) = OldVelocityArray(1,:)

    !print *, "RESET"
    !print *, "Before Crossing", OLDGSMPosition(1), OLDGSMPosition(2), OLDGSMPosition(3)
    !print *, "Distance Before Crossing", (OLDGSMPosition(1)**2 + OLDGSMPosition(2)**2 + OLDGSMPosition(3)**2)**0.5
    !print *, "After Crossing", GSMPosition(1), GSMPosition(2), GSMPosition(3)
    !print *, "Distance After Crossing", (GSMPosition(1)**2 + GSMPosition(2)**2 + GSMPosition(3)**2)**0.5

    secondTotal = OLDsecondTotal
    termtype = 0

    END IF
END IF

end subroutine FinalStepCheck