! **********************************************************************************************************************
! Subroutine Cutoff:
!            subroutine that calculates the trajectory of a cosmic ray across a range of rigidities
!            within different input magnetic field models and determines the cutoff rigidity.
!            Will create a csv file in which the asymptotic cone data is stored.
!            Output can be the vertical cutoff or apparent cutoff
!            Apparent cutoff computation is roughly 9 times as long as vertical
!
! **********************************************************************************************************************
subroutine cutoff(PositionIN, StartRigidity, EndRigidity, RigidityStep, Date, mode, IntMode, & 
    AtomicNumber, Anti, I, Wind, Pause, FileName, CoordSystem, GyroPercent, End, Rcomputation, scanchoice, &
    gOTSO, hOTSO, Rigidities)
    USE Particle
    USE SolarWind
    USE MagneticFieldFunctions
    USE MagnetopauseFunctions
    USE IntegrationFunctions
    USE GEOPACK1
    USE GEOPACK2
    USE Magnetopause
    USE CUSTOMGAUSS
    implicit none
    
    real(8) :: PositionIN(5), StartRigidity, EndRigidity, RigidityScan, RigidityStep, Date(6), End(3)
    real(8) :: Wind(17), Re, Lat, Long, GyroPercent, EndLoop
    real(8) :: Geofile(3), RuMemory(9), RlMemory(9), RefMemory(9), Rigidity(3)
    real(8) :: Zenith(9), Azimuth(9), sumrl, sumru, sumref
    integer(8) :: mode(2), IntMode, Anti, AtomicNumber
    integer(4) :: I, Limit, bool, Pause, stepNum, loop, Rcomputation, scanchoice, scan, LastCheck
    character(len=50) :: FileName
    character(len=3) :: CoordSystem
    real(8) :: gOTSO(105), hOTSO(105)

    real(8), intent(out) :: Rigidities(3)
    
    R = real(StartRigidity, kind = selected_real_kind(15,307))
    Re = 6371.2
    Limit = 0
    Acount = 0
    Result = 0
    stepNum = 0
    loop = 1
    forbiddencount = 0
    NeverFail = 0
    Step = RigidityStep
    SubResult = 0
    MaxGyroPercent = GyroPercent
    sumrl = 0
    sumref = 0
    sumru = 0
    LastCheck = 0
    FailCheck = 0

    if (mode(1) == 4) then
    Ginput = gOTSO
    Hinput = hOTSO
    end if

    Rigidity(1) = StartRigidity
    Rigidity(2) = EndRigidity
    Rigidity(3) = RigidityStep

    RigidityScan = 0.50
    RigidityStep = 0.50

    Zenith(1) = 0
    Zenith(2) = 30
    Zenith(3) = 30
    Zenith(4) = 30
    Zenith(5) = 30
    Zenith(6) = 30
    Zenith(7) = 30
    Zenith(8) = 30
    Zenith(9) = 30

    Azimuth(1) = 0
    Azimuth(2) = 0
    Azimuth(3) = 45
    Azimuth(4) = 90
    Azimuth(5) = 135
    Azimuth(6) = 180
    Azimuth(7) = 225
    Azimuth(8) = 270
    Azimuth(9) = 315

    IF (scanchoice == 1) THEN
        scan = 0
        RigidityStep = RigidityScan
    ELSE
        scan = 1
        RigidityStep = Rigidity(3)
    END IF


    IF (Rcomputation == 0) THEN
        EndLoop = 1.0
    ELSE
        EndLoop = 9.0
    END IF

    RigidityStep = real(RigidityStep, kind = selected_real_kind(10,307))

    !open(unit=10,file=FileName,status='replace')
    !write(10,"(a)")"Zenith,Azimuth,Ru,Rc,Rl"

    do while (loop <= EndLoop)
   
    100 do while (R > EndRigidity)

    PositionIN(4) = Zenith(loop)
    PositionIN(5) = Azimuth(loop)
    
    IF (R < Rigidity(3)) THEN
        R = EndRigidity
        GOTO 50
    END IF

    R = real(R, kind = selected_real_kind(10,307))
    RigidityStep = real(RigidityStep, kind = selected_real_kind(10,307))
    
    call CreateParticle(PositionIN, R, Date, AtomicNumber, Anti, mode)
    
    call initializeWind(Wind, I, mode)
    call initializeCustomGauss(mode)

    call MagneticFieldAssign(mode)
    call MagnetopauseAssign(Pause)
    call IntegrationAssign(IntMode)

    call FirstTimeStep()
    test = 0

    do while (Result == 0)
    
    call IntegrationPointer()

    call EscapeCheck()

    !call FinalStepCheck()
    
    IF (Position(1) < End(1) ) THEN
        bool = -1
        Limit = 1
        forbiddencount = forbiddencount + 1
        NeverFail = 1
        FailCheck = 1
        call AsymptoticDirection(Lat, Long)
        call CoordinateTransform("GDZ", CoordSystem, year, day, secondTotal, Position, GEOfile)
        !print *, R, " ", "Returned to Earth"
        EXIT
    END IF

    IF (End(2) == 0) THEN
        
    ELSE IF (DistanceTraveled/1000.0 > End(2) * Re) THEN
        bool = 0
        Limit = 1
        forbiddencount = forbiddencount + 1
        NeverFail = 1
        FailCheck = 1
        call AsymptoticDirection(Lat, Long)
        call CoordinateTransform("GDZ", CoordSystem, year, day, secondTotal, Position, GEOfile)
        !print *, R, " ", "Trapped"
        EXIT
    END IF

    IF (End(3) == 0) THEN
        
    ELSE IF (TimeElapsed > End(3)) THEN
        bool = 0
        Limit = 1
        forbiddencount = forbiddencount + 1
        NeverFail = 1
        FailCheck = 1
        call AsymptoticDirection(Lat, Long)
        call CoordinateTransform("GDZ", CoordSystem, year, day, secondTotal, Position, GEOfile)
        !print *, R, " ", "Time Elapsed"
        EXIT
    END IF
    
    IF (Result == 1) THEN
        call AsymptoticDirection(Lat, Long)
        call CoordinateTransform("GDZ", CoordSystem, year, day, secondTotal, Position, GEOfile)
        forbiddencount = 0
        bool = 1
        RL = R
        !print *, R, " ", "Escaped"
        IF (Limit == 0) THEN
            RU = R
        ELSE IF (Limit == 1) THEN
            Acount = Acount + 1
        END IF
        EXIT
    END IF
    
    end do

    stepNum = stepNum + 1

    R = (StartRigidity - (stepNum*RigidityStep))
    IF(R < EndRigidity) THEN
        R = EndRigidity
    ELSE IF (R < RigidityStep) THEN
        IF(NeverFail == 0) THEN
            IF(LastCheck == 1) THEN
                R = EndRigidity
            END IF
            IF(LastCheck == 0) THEN
                R = Rigidity(3)
                LastCheck = 1
                stepNum = stepNum - 1
            END IF
        END IF
    END IF
    50 Result = 0

    PositionIN(4) = Zenith(loop)
    PositionIN(5) = Azimuth(loop)

    end do

    call EffectiveRigidity(RigidityStep)

    IF (scan == 0) THEN
        scan = 1
        StartRigidity = RU + 1.5*RigidityScan
        EndRigidity = RL - 1.5*RigidityScan
        IF (EndRigidity < 0) THEN
            EndRigidity = 0
        END IF
        R = StartRigidity
        RigidityStep = Rigidity(3)
        Limit = 0
        Acount = 0
        Result = 0
        forbiddencount = 0
        SubResult = 0
        stepNum = 0
        IF(NeverFail == 1) THEN
            NeverFail = 0
            GOTO 100
        ELSE IF(NeverFail == 0) THEN
            RU = 0
            RL = 0
            Ref = 0
        END IF
    END IF

    RlMemory(loop) = Rl
    RefMemory(loop) = Ref
    RuMemory(loop) = Ru

    IF (scanchoice == 1) THEN
        scan = 0
        RigidityStep = RigidityScan
        StartRigidity = Rigidity(1)
        EndRigidity = Rigidity(2)
    ELSE
        scan = 1
        RigidityStep = Rigidity(3)
    END IF

    !write(10,'(*(G0.6,:""))')PositionIN(4), ",", PositionIN(5), ",", Ru, ",", Ref, ",", Rl 

    loop = loop + 1
    PositionIN(4) = Zenith(loop)
    PositionIN(5) = Azimuth(loop)
    R = real(StartRigidity, kind = selected_real_kind(15,307))
    Re = 6371.2
    Limit = 0
    Acount = 0
    Result = 0
    stepNum = 0
    forbiddencount = 0
    NeverFail = 0
    SubResult = 0
    LastCheck = 0
    FailCheck = 0

    end do

    IF(Rcomputation /= 0) THEN
        sumrl = RLMemory(1)/2.0
        do i = 2, 9
            sumrl = sumrl + RLMemory(i)/16.0
        end do
        sumru = RUMemory(1)/2.0
        do i = 2, 9
            sumru = sumru + RUMemory(i)/16.0
        end do
        sumref = RefMemory(1)/2.0
        do i = 2, 9
            sumref = sumref + RefMemory(i)/16.0
        end do

        RU = sumru
        RL = sumrl
        Ref = sumref
    END IF

    !print *, "Ru:", RU
    !print *, "Rl:", RL
    !print *, "Rc:", Ref

    Rigidities(1) = RU
    Rigidities(2) = Ref
    Rigidities(3) = RL

    !write(10,'(*(G0.6,:""))')"Ru:", RU, ",  Rc:", Ref, ",  Rl:", RL 

    Close(10, STATUS='KEEP') 

end subroutine cutoff

! **********************************************************************************************************************
! Subroutine Cone:
!            subroutine that calculates the trajectory of a cosmic ray across a range of rigidities
!            within different input magnetic field models and determines the Asymptotic Latitude and
!            Longitude. Will create a csv file in which the asymptotic cone data is stored.
!            Accepted Condition: cosmic ray encounters the model magnetopause
!            Forbidden Conditions: - cosmic ray encounters the Earth (20km above Earth's surface)
!                                  - cosmic ray travels over 100Re without escaping or encountering Earth
!                                  - cosmic ray is simulated for a given period of time
!
! **********************************************************************************************************************
subroutine cone(PositionIN, StartRigidity, EndRigidity, RigidityStep, Date, mode, IntMode, & 
    AtomicNumber, Anti, I, Wind, Pause, FileName, CoordSystem, GyroPercent, End, &
    length, gOTSO, hOTSO, ConeArray, Rigidities)
    USE Particle
    USE SolarWind
    USE MagneticFieldFunctions
    USE MagnetopauseFunctions
    USE IntegrationFunctions
    USE GEOPACK1
    USE GEOPACK2
    USE Magnetopause
    USE CUSTOMGAUSS
    implicit none
    
    real(8) :: PositionIN(5), StartRigidity, EndRigidity, RigidityStep, Date(6), End(3)
    real(8) :: Wind(17), Re, Lat, Long, GyroPercent, Rigidities(3)
    real(8) :: Geofile(3)
    integer(8) :: mode(2), IntMode, Anti, AtomicNumber
    integer(4) :: I, Limit, bool, Pause, stepNum
    character(len=50) :: FileName
    character(len=3) :: CoordSystem
    character(len=100) :: ConeArray(1, length)
    character(len=100) :: temp_array(1, length)
    character(len=100) :: temp_string
    integer(4) :: length, old_size
    real(8) :: gOTSO(105), hOTSO(105)

    intent(out) :: ConeArray, Rigidities
    
    R = real(StartRigidity, kind = selected_real_kind(15,307))
    Re = 6371.2
    old_size = 0
    Limit = 0
    Acount = 0
    Result = 0
    stepNum = 0
    forbiddencount = 0
    NeverFail = 0
    Step = RigidityStep
    SubResult = 0
    MaxGyroPercent = GyroPercent

    if (mode(1) == 4) then
    Ginput = gOTSO
    Hinput = hOTSO
    end if

    RigidityStep = real(RigidityStep, kind = selected_real_kind(10,307))

    IF (PositionIN(4) > 90.0) THEN
        print *, "ERROR: Please enter a zenith angle between 0 and 90 degrees"
        stop
    END IF
    
    IF (PositionIN(5) < 0) THEN
        print *, "ERROR: Please enter an azimuth angle between 0 and 360 degrees"
        print *, "N = 0, E = 90, S = 180, and W = 270 (degrees)"
        stop
    ELSE IF (PositionIN(5) > 360) THEN
        print *, "ERROR: Please enter an azimuth angle between 0 and 360 degrees"
        print *, "N = 0, E = 90, S = 180, and W = 270 (degrees)"
        stop
    END IF

    do while (R > EndRigidity)

    R = real(R, kind = selected_real_kind(10,307))
    RigidityStep = real(RigidityStep, kind = selected_real_kind(10,307))

    call CreateParticle(PositionIN, R, Date, AtomicNumber, Anti, mode)
    call initializeWind(Wind, I, mode)
    call initializeCustomGauss(mode)
    call MagneticFieldAssign(mode)
    call MagnetopauseAssign(Pause)
    call IntegrationAssign(IntMode)

    call FirstTimeStep()
    
    do while (Result == 0) 
    
    call IntegrationPointer()

    call EscapeCheck()

    call FinalStepCheck()
    
    IF (Position(1) < End(1) ) THEN
        bool = -1
        Limit = 1
        forbiddencount = forbiddencount + 1
        NeverFail = 1
        FailCheck = 1
        call AsymptoticDirection(Lat, Long)
        call CoordinateTransform("GDZ", CoordSystem, year, day, secondTotal, Position, GEOfile)
        !print *, R, lat, Long, "Forbidden", "      Encountered Earth" !(Prints the outputs to the command module while running (Can lead to delays with multi-core proccessing))
        EXIT
    END IF

    IF (End(2) == 0) THEN
        
    ELSE IF (DistanceTraveled/1000.0 > End(2) * Re) THEN
        bool = 0
        Limit = 1
        forbiddencount = forbiddencount + 1
        NeverFail = 1
        FailCheck = 1
        call AsymptoticDirection(Lat, Long)
        call CoordinateTransform("GDZ", CoordSystem, year, day, secondTotal, Position, GEOfile)
        !print *, R, lat, Long, "Forbidden",  "      Exceeded Travel Distance Without Escape" !(Prints the outputs to the command module while running (Can lead to delays with multi-core proccessing))
        EXIT
    END IF

    IF (End(3) == 0) THEN
        
    ELSE IF (TimeElapsed > End(3)) THEN
        bool = 0
        Limit = 1
        forbiddencount = forbiddencount + 1
        NeverFail = 1
        FailCheck = 1
        call AsymptoticDirection(Lat, Long)
        call CoordinateTransform("GDZ", CoordSystem, year, day, secondTotal, Position, GEOfile)
        !print *, R, lat, Long, "Forbidden",  "      Maximum Time Exceeded" !(Prints the outputs to the command module while running (Can lead to delays with multi-core proccessing))
        EXIT
    END IF
    
    IF (Result == 1) THEN
        call AsymptoticDirection(Lat, Long)
        call CoordinateTransform("GDZ", CoordSystem, year, day, secondTotal, Position, GEOfile)
        forbiddencount = 0
        bool = 1
        !print *, R, lat, Long !(Prints the outputs to the command module while running (Can lead to delays with multi-core proccessing))
        RL = R
        FailCheck = 0
        IF (Limit == 0) THEN
            RU = R
        ELSE IF (Limit == 1) THEN
            Acount = Acount + 1
        END IF
        EXIT
    END IF

    end do
    
    write(temp_string, '(F7.3, 1X, I5, 1X, F7.3, 1X, F7.3)') R, bool, Lat, Long

    ConeArray(1, old_size + 1) = temp_string
    old_size = old_size + 1
    
    stepNum = stepNum + 1

    R = (StartRigidity - (stepNum*RigidityStep))
    IF (EndRigidity < R .AND. R < RigidityStep) THEN
        R = RigidityStep
    ELSEIF (R < EndRigidity) THEN
        R = EndRigidity
    END IF
    Result = 0

    end do

    call EffectiveRigidity(RigidityStep)
    Rigidities(1) = RU
    Rigidities(2) = Ref
    Rigidities(3) = RL

end subroutine cone


! **********************************************************************************************************************
! Subroutine Trajectory:
!            subroutine that calculates the trajectory of a cosmic ray within different input
!            magnetic field models.
!            The output data is in the user given coordinate system.
!            Will also state if the cosmic ray has an allowed or forbidden trajectory.
!            Accepted Condition: cosmic ray encounters the magnetopause
!            Forbidden Conditions: - cosmic ray encounters the Earth (20km above Earth's surface)
!                                  - cosmic ray travels over 100Re without escaping or encountering Earth
!                                  - cosmic ray is simulated for a given period of time
!
! **********************************************************************************************************************
subroutine trajectory(PositionIN, Rigidity, Date, mode, IntMode, & 
    AtomicNumber, Anti, I, Wind, Pause, FileName, CoordSystem, GyroPercent, End, &
    gOTSO, hOTSO, bool, Lat, Long)
USE Particle
USE GEOPACK1
USE GEOPACK2
USE SolarWind
USE MagneticFieldFunctions
USE MagnetopauseFunctions
USE IntegrationFunctions
USE Magnetopause
USE CUSTOMGAUSS
implicit none

real(8) :: PositionIN(5), Rigidity, Date(6), End(3)
real(8) :: Wind(17), Re, GyroPercent
real(8) :: Xnew(3), XnewConverted(3)
integer(8) :: mode(2), IntMode, Anti, AtomicNumber
integer(4) :: I, Limit, Pause
character(len=50) :: FileName
character(len=3) :: CoordSystem

real(8), intent(out) :: Lat, Long
integer(4), intent(out) :: bool
real(8) :: gOTSO(105), hOTSO(105)


Re = 6371.2
Limit = 0
Acount = 0
Result = 0
SubResult = 0
MaxGyroPercent = GyroPercent

if (mode(1) == 4) then
Ginput = gOTSO
Hinput = hOTSO
end if

IF (PositionIN(4) > 90.0) THEN
    print *, "ERROR: Please enter a zenith angle between 0 and 90 degrees"
    stop
END IF

IF (PositionIN(5) < 0) THEN
    print *, "ERROR: Please enter an azimuth angle between 0 and 360 degrees"
    print *, "N = 0, E = 90, S = 180, and W = 270 (degrees)"
    stop
ELSE IF (PositionIN(5) > 360) THEN
    print *, "ERROR: Please enter an azimuth angle between 0 and 360 degrees"
    print *, "N = 0, E = 90, S = 180, and W = 270 (degrees)"
    stop
END IF

call CreateParticle(PositionIN, Rigidity, Date, AtomicNumber, Anti, mode)

call initializeWind(Wind, I, mode)
call initializeCustomGauss(mode)

call MagneticFieldAssign(mode)
call MagnetopauseAssign(Pause)
call IntegrationAssign(IntMode)

call FirstTimeStep()

!open(unit=10,file=FileName,status='replace')
!write(10,"(a)")"X,Y,Z"

IF (Rigidity == 0) THEN
  Result = 1
END IF

do while (Result == 0)
call IntegrationPointer
call EscapeCheck()
call FinalStepCheck()
Xnew(1) = XnewTemp(1)/1000
Xnew(2) = XnewTemp(2)/1000
Xnew(3) = XnewTemp(3)/1000

call CoordinateTransform("GSM", CoordSystem, year, day, secondTotal, Xnew, XnewConverted)

if (model(1) == 4) then
    if (CoordSystem == "GEO") then
        XnewConverted = Xnew
    else
        call CoordinateTransform("GEO", CoordSystem, year, day, secondTotal, Xnew, XnewConverted)
    end if
end if

!write(10,'(*(G0.6,:,","))') XnewConverted

IF (Position(1) < End(1) ) THEN
    !print *, "This is Forbidden", "      Encountered Earth"
    call AsymptoticDirection(Lat, Long)
    bool = -1
    !print *, "Final Position (Latitude, Longitude)"
    !print *, Position
    !print *, "Asymptotic Directions (Latitude, Longitude)"
    !print *, Lat, Long
    Limit = 1
    EXIT
END IF

IF (End(2) == 0) THEN
    
ELSE IF ( DistanceTraveled/1000.0 > End(2)*Re) THEN
    !print *, "This is Forbidden", "      Exceeded Travel Distance Without Escape"
    !print *, DistanceTraveled
    call AsymptoticDirection(Lat, Long)
    bool = 0
    !print *, "Final Position (Latitude, Longitude)"
    !print *, Position
    !print *, "Asymptotic Directions (Latitude, Longitude)"
    !print *, Lat, Long
    Limit = 1
    EXIT
END IF

IF (End(3) == 0) THEN

ELSE IF ( TimeElapsed > End(3)) THEN
    !print *, "This is Forbidden", "      Exceeded Maximum Time"
    !print *, TimeElapsed
    call AsymptoticDirection(Lat, Long)
    bool = 0
    !print *, "Final Position (Latitude, Longitude)"
    !print *, Position
    !print *, "Asymptotic Directions (Latitude, Longitude)"
    !print *, Lat, Long
    Limit = 1
    EXIT
END IF


IF (Result == 1)  THEN
    !print *, "This is Allowed", "      Successfully Escaped"
    call AsymptoticDirection(Lat, Long)
    bool = 1
    !print *, "Escape Position (Altitude [km], Latitude, Longitude)"
    !print *, Position
    !print *, "Asymptotic Directions (Latitude, Longitude)"
    !print *, Lat, Long
    EXIT
END IF

end do

!Close(10, STATUS='KEEP') 

end subroutine trajectory


! **********************************************************************************************************************
! Subroutine Planet:
!            subroutine that calculates the trajectory of a cosmic ray across a range of rigidities
!            within different input magnetic field models and determines the effective cutoff rigidity for
!            a range of latitude and longitudes. Typically done over the entire planet.
!            Will create a csv file with the calculated rigidities for the locations.
!            This code works with the Planet.py tool to assign large amounts of locations across multiple cores.
!            A planet.csv file will be produced storing the data for latitude, longitude, and cutoff
!
! **********************************************************************************************************************
subroutine planet(PositionIN, Rigidity, Date, mode, IntMode, AtomicNumber, Anti, I, Wind, Pause, &
     FileName, GyroPercent, End, Rcomputation, scanchoice, gOTSO, hOTSO, Rigidities)
    USE Particle
    USE SolarWind
    USE MagneticFieldFunctions
    USE MagnetopauseFunctions
    USE IntegrationFunctions
    USE GEOPACK1
    USE GEOPACK2
    USE Magnetopause
    USE CUSTOMGAUSS
    implicit none
    
    real(8) :: PositionIN(5), StartRigidity, EndRigidity, RigidityScan, RigidityStep, Date(6), End(3)
    real(8) :: Wind(17), Re, Lat, Long, GyroPercent, EndLoop
    real(8) :: Geofile(3), RuMemory(9), RlMemory(9), RefMemory(9), Rigidity(3)
    real(8) :: Zenith(9), Azimuth(9), sumrl, sumru, sumref
    integer(8) :: mode(2), IntMode, Anti, AtomicNumber
    integer(4) :: I, Limit, bool, Pause, stepNum, loop, Rcomputation, scanchoice, scan, LastCheck
    character(len=50) :: FileName
    character(len=3) :: CoordSystem
    real(8) :: gOTSO(105), hOTSO(105)

    real(8), intent(out) :: Rigidities(3)
    
    R = real(StartRigidity, kind = selected_real_kind(15,307))
    Re = 6371.2
    Limit = 0
    Acount = 0
    Result = 0
    stepNum = 0
    loop = 1
    forbiddencount = 0
    NeverFail = 0
    Step = RigidityStep
    SubResult = 0
    MaxGyroPercent = GyroPercent
    sumrl = 0
    sumref = 0
    sumru = 0
    LastCheck = 0
    FailCheck = 0

    if (mode(1) == 4) then
    Ginput = gOTSO
    Hinput = hOTSO
    end if

    Rigidity(1) = StartRigidity
    Rigidity(2) = EndRigidity
    Rigidity(3) = RigidityStep

    RigidityScan = 0.50
    RigidityStep = 0.50

    Zenith(1) = 0
    Zenith(2) = 30
    Zenith(3) = 30
    Zenith(4) = 30
    Zenith(5) = 30
    Zenith(6) = 30
    Zenith(7) = 30
    Zenith(8) = 30
    Zenith(9) = 30

    Azimuth(1) = 0
    Azimuth(2) = 0
    Azimuth(3) = 45
    Azimuth(4) = 90
    Azimuth(5) = 135
    Azimuth(6) = 180
    Azimuth(7) = 225
    Azimuth(8) = 270
    Azimuth(9) = 315

    IF (scanchoice == 1) THEN
        scan = 0
        RigidityStep = RigidityScan
    ELSE
        scan = 1
        RigidityStep = Rigidity(3)
    END IF


    IF (Rcomputation == 0) THEN
        EndLoop = 1.0
    ELSE
        EndLoop = 9.0
    END IF

    RigidityStep = real(RigidityStep, kind = selected_real_kind(10,307))

    do while (loop <= EndLoop)
   
    100 do while (R > EndRigidity)

    PositionIN(4) = Zenith(loop)
    PositionIN(5) = Azimuth(loop)
    
    IF (R < Rigidity(3)) THEN
        R = EndRigidity
        GOTO 50
    END IF

    R = real(R, kind = selected_real_kind(10,307))
    RigidityStep = real(RigidityStep, kind = selected_real_kind(10,307))
    
    call CreateParticle(PositionIN, R, Date, AtomicNumber, Anti, mode)
    
    call initializeWind(Wind, I, mode)
    call initializeCustomGauss(mode)

    call MagneticFieldAssign(mode)
    call MagnetopauseAssign(Pause)
    call IntegrationAssign(IntMode)

    call FirstTimeStep()
    test = 0

    do while (Result == 0)
    
    call IntegrationPointer()

    call EscapeCheck()

    !call FinalStepCheck()
    
    IF (Position(1) < End(1) ) THEN
        bool = -1
        Limit = 1
        forbiddencount = forbiddencount + 1
        NeverFail = 1
        FailCheck = 1
        call AsymptoticDirection(Lat, Long)
        call CoordinateTransform("GDZ", CoordSystem, year, day, secondTotal, Position, GEOfile)
        !print *, R, " ", "Returned to Earth"
        EXIT
    END IF

    IF (End(2) == 0) THEN
        
    ELSE IF (DistanceTraveled/1000.0 > End(2) * Re) THEN
        bool = 0
        Limit = 1
        forbiddencount = forbiddencount + 1
        NeverFail = 1
        FailCheck = 1
        call AsymptoticDirection(Lat, Long)
        call CoordinateTransform("GDZ", CoordSystem, year, day, secondTotal, Position, GEOfile)
        !print *, R, " ", "Trapped"
        EXIT
    END IF

    IF (End(3) == 0) THEN
        
    ELSE IF (TimeElapsed > End(3)) THEN
        bool = 0
        Limit = 1
        forbiddencount = forbiddencount + 1
        NeverFail = 1
        FailCheck = 1
        call AsymptoticDirection(Lat, Long)
        call CoordinateTransform("GDZ", CoordSystem, year, day, secondTotal, Position, GEOfile)
        !print *, R, " ", "Time Elapsed"
        EXIT
    END IF
    
    IF (Result == 1) THEN
        call AsymptoticDirection(Lat, Long)
        call CoordinateTransform("GDZ", CoordSystem, year, day, secondTotal, Position, GEOfile)
        forbiddencount = 0
        bool = 1
        RL = R
        !print *, R, " ", "Escaped"
        IF (Limit == 0) THEN
            RU = R
        ELSE IF (Limit == 1) THEN
            Acount = Acount + 1
        END IF
        EXIT
    END IF
    
    end do

    stepNum = stepNum + 1

    R = (StartRigidity - (stepNum*RigidityStep))
    IF(R < EndRigidity) THEN
        R = EndRigidity
    ELSE IF (R < RigidityStep) THEN
        IF(NeverFail == 0) THEN
            IF(LastCheck == 1) THEN
                R = EndRigidity
            END IF
            IF(LastCheck == 0) THEN
                R = Rigidity(3)
                LastCheck = 1
                stepNum = stepNum - 1
            END IF
        END IF
    END IF
    50 Result = 0

    PositionIN(4) = Zenith(loop)
    PositionIN(5) = Azimuth(loop)

    end do

    call EffectiveRigidity(RigidityStep)

    IF (scan == 0) THEN
        scan = 1
        StartRigidity = RU + 1.5*RigidityScan
        EndRigidity = RL - 1.5*RigidityScan
        IF (EndRigidity < 0) THEN
            EndRigidity = 0
        END IF
        R = StartRigidity
        RigidityStep = Rigidity(3)
        Limit = 0
        Acount = 0
        Result = 0
        forbiddencount = 0
        SubResult = 0
        stepNum = 0
        IF(NeverFail == 1) THEN
            NeverFail = 0
            GOTO 100
        ELSE IF(NeverFail == 0) THEN
            RU = 0
            RL = 0
            Ref = 0
        END IF
    END IF

    RlMemory(loop) = Rl
    RefMemory(loop) = Ref
    RuMemory(loop) = Ru

    IF (scanchoice == 1) THEN
        scan = 0
        RigidityStep = RigidityScan
        StartRigidity = Rigidity(1)
        EndRigidity = Rigidity(2)
    ELSE
        scan = 1
        RigidityStep = Rigidity(3)
    END IF

    !write(10,'(*(G0.6,:""))')PositionIN(4), ",", PositionIN(5), ",", Ru, ",", Ref, ",", Rl 

    loop = loop + 1
    PositionIN(4) = Zenith(loop)
    PositionIN(5) = Azimuth(loop)
    R = real(StartRigidity, kind = selected_real_kind(15,307))
    Re = 6371.2
    Limit = 0
    Acount = 0
    Result = 0
    stepNum = 0
    forbiddencount = 0
    NeverFail = 0
    SubResult = 0
    LastCheck = 0
    FailCheck = 0

    end do

    IF(Rcomputation /= 0) THEN
        sumrl = RLMemory(1)/2.0
        do i = 2, 9
            sumrl = sumrl + RLMemory(i)/16.0
        end do
        sumru = RUMemory(1)/2.0
        do i = 2, 9
            sumru = sumru + RUMemory(i)/16.0
        end do
        sumref = RefMemory(1)/2.0
        do i = 2, 9
            sumref = sumref + RefMemory(i)/16.0
        end do

        RU = sumru
        RL = sumrl
        Ref = sumref
    END IF

    !print *, "Ru:", RU
    !print *, "Rl:", RL
    !print *, "Rc:", Ref

    Rigidities(1) = RU
    Rigidities(2) = Ref
    Rigidities(3) = RL
        
end subroutine planet

! **********************************************************************************************************************
! Subroutine Trajectory:
!            subroutine that calculates the trajectory of a cosmic ray within different input
!            magnetic field models.
!            The output data is in the user given coordinate system.
!            Will also state if the cosmic ray has an allowed or forbidden trajectory.
!            Accepted Condition: cosmic ray encounters the magnetopause
!            Forbidden Conditions: - cosmic ray encounters the Earth (20km above Earth's surface)
!                                  - cosmic ray travels over 100Re without escaping or encountering Earth
!                                  - cosmic ray is simulated for a given period of time
!
! **********************************************************************************************************************
subroutine trajectory_full(PositionIN, Rigidity, Date, mode, IntMode, & 
    AtomicNumber, Anti, I, Wind, Pause, FileName, CoordSystem, GyroPercent, &
    End, gOTSO, hOTSO)
USE Particle
USE GEOPACK1
USE GEOPACK2
USE SolarWind
USE MagneticFieldFunctions
USE MagnetopauseFunctions
USE IntegrationFunctions
USE Magnetopause
USE CUSTOMGAUSS
implicit none

real(8) :: PositionIN(5), Rigidity, Date(6), End(3)
real(8) :: Wind(17), Re, Lat, Long, GyroPercent
real(8) :: Xnew(3), XnewConverted(3)
integer(8) :: mode(2), IntMode, Anti, AtomicNumber
integer(4) :: I, Limit, Pause
character(len=50) :: FileName
character(len=3) :: CoordSystem
real(8) :: gOTSO(105), hOTSO(105)

Re = 6371.2
Limit = 0
Acount = 0
Result = 0
SubResult = 0
MaxGyroPercent = GyroPercent

if (mode(1) == 4) then
Ginput = gOTSO
Hinput = hOTSO
end if

IF (PositionIN(4) > 90.0) THEN
    print *, "ERROR: Please enter a zenith angle between 0 and 90 degrees"
    stop
END IF

IF (PositionIN(5) < 0) THEN
    print *, "ERROR: Please enter an azimuth angle between 0 and 360 degrees"
    print *, "N = 0, E = 90, S = 180, and W = 270 (degrees)"
    stop
ELSE IF (PositionIN(5) > 360) THEN
    print *, "ERROR: Please enter an azimuth angle between 0 and 360 degrees"
    print *, "N = 0, E = 90, S = 180, and W = 270 (degrees)"
    stop
END IF

call CreateParticle(PositionIN, Rigidity, Date, AtomicNumber, Anti, mode)

call initializeWind(Wind, I, mode)
call initializeCustomGauss(mode)

call MagneticFieldAssign(mode)
call MagnetopauseAssign(Pause)
call IntegrationAssign(IntMode)

call FirstTimeStep()

open(unit=10,file=FileName,status='replace')
write(10,"(a)")"X,Y,Z"

do while (Result == 0) 
call IntegrationPointer
call EscapeCheck()
call FinalStepCheck()
Xnew(1) = XnewTemp(1)/1000
Xnew(2) = XnewTemp(2)/1000
Xnew(3) = XnewTemp(3)/1000

call CoordinateTransform("GSM", CoordSystem, year, day, secondTotal, Xnew, XnewConverted)

if (model(1) == 4) then
    if (CoordSystem == "GEO") then
        XnewConverted = Xnew
    else
        call CoordinateTransform("GEO", CoordSystem, year, day, secondTotal, Xnew, XnewConverted)
    end if
end if

write(10,'(*(G0.6,:,","))') XnewConverted/Re

IF (Position(1) < End(1) ) THEN
    !print *, "This is Forbidden", "      Encountered Earth"
    call AsymptoticDirection(Lat, Long)
    !print *, "Final Position (Latitude, Longitude)"
    !print *, Position
    !print *, "Asymptotic Directions (Latitude, Longitude)"
    !print *, Lat, Long
    Limit = 1
    EXIT
END IF

IF (End(2) == 0) THEN
    
ELSE IF ( DistanceTraveled/1000.0 > End(2)*Re) THEN
    !print *, "This is Forbidden", "      Exceeded Travel Distance Without Escape"
    !print *, DistanceTraveled
    call AsymptoticDirection(Lat, Long)
    !print *, "Final Position (Latitude, Longitude)"
    !print *, Position
    !print *, "Asymptotic Directions (Latitude, Longitude)"
    !print *, Lat, Long
    Limit = 1
    EXIT
END IF

IF (End(3) == 0) THEN

ELSE IF ( TimeElapsed > End(3)) THEN
    !print *, "This is Forbidden", "      Exceeded Maximum Time"
    !print *, TimeElapsed
    call AsymptoticDirection(Lat, Long)
    !print *, "Final Position (Latitude, Longitude)"
    !print *, Position
    !print *, "Asymptotic Directions (Latitude, Longitude)"
    !print *, Lat, Long
    Limit = 1
    EXIT
END IF


IF (Result == 1)  THEN
    !print *, "This is Allowed", "      Successfully Escaped"
    call AsymptoticDirection(Lat, Long)
    !print *, "Escape Position (Altitude [km], Latitude, Longitude)"
    !print *, Position
    !print *, "Asymptotic Directions (Latitude, Longitude)"
    !print *, Lat, Long
    EXIT
END IF

end do

Close(10, STATUS='KEEP') 

end subroutine trajectory_full



subroutine GETTSY04DATAWINDOWS(OMNIYEAR, length)
    USE GEOPACK1
    implicit none
    
    integer(4) :: OMNIYEAR,inputyear,length
    character(len=length) :: DIRECTORY2
    

    inputyear = OMNIYEAR
    
    CALL FILLIMFGAPS(inputyear, length, DIRECTORY2)
    CALL FILLSWGAPS(inputyear, DIRECTORY2, length)
    CALL PREPAREINTERVALS1(inputyear, DIRECTORY2, length)
    CALL PREPAREINPUT4(inputyear, DIRECTORY2, length)
    
    
    end subroutine GETTSY04DATAWINDOWS


    subroutine GETTSY04DATALINUX(OMNIYEAR, DIRECTORY, length)
    USE GEOPACK1
    implicit none
    
    integer(4) :: OMNIYEAR,inputyear,length
    character(len=length) :: DIRECTORY
    

    inputyear = OMNIYEAR
    
    CALL FILLIMFGAPSLINUX(inputyear, length, DIRECTORY)
    CALL FILLSWGAPS(inputyear, DIRECTORY, length)
    CALL PREPAREINTERVALS1(inputyear, DIRECTORY, length)
    CALL PREPAREINPUT4(inputyear, DIRECTORY, length)
    
    
    end subroutine GETTSY04DATALINUX


! **********************************************************************************************************************
! Subroutine MagStrength:
!            subroutine that will tell you the strength of the magnetic field at any given point within the 
!            magnetosphere. Output is in GSM coordinates.
!
! **********************************************************************************************************************
subroutine MagStrength(Pin, Date, mode, I, Wind, CoordIN, gOTSO, hOTSO, Bfield)
    USE Particle
    USE SolarWind
    USE MagneticFieldFunctions
    USE MagnetopauseFunctions
    USE GEOPACK1
    USE GEOPACK2
    USE CUSTOMGAUSS
    implicit none
    
    real(8) :: Pin(3), Pout(3), Wind(17), Date(6)
    character(len = 3) :: CoordIN
    integer(4) :: I
    integer(8) :: mode(2)
    real(8) :: gOTSO(105), hOTSO(105)

    real(8), intent(out) :: Bfield(3) 

    if (mode(1) == 4) then
       Ginput = gOTSO
       Hinput = hOTSO
    end if

    year = INT(Date(1))
    day = INT(Date(2))
    hour = INT(Date(3))
    minute = INT(Date(4))
    secondINT = INT(Date(5))
    secondTotal = real(Date(6))

    call initializeWind(Wind, I, mode)
    call initializeCustomGauss(mode)

    call MagneticFieldAssign(mode)
    
    call CoordinateTransform(CoordIN, "GSM", year, day, secondTotal, Pin, Pout)

    call MagFieldCheck(Pout, Bfield)
    
    end subroutine MagStrength

! **********************************************************************************************************************
! Subroutine CoordTrans:
!            subroutine that uses the IBREM database of coordinate transforms to convert coordinates into
!            a new coordinate system.
!
! **********************************************************************************************************************
subroutine CoordTrans(Pin, year, day, hour, minute, secondINT, secondTotal, CoordIN, CoordOUT, Pout)
    implicit none
    
    real(8) :: sec, Pin(3), secondTotal
    character(len = 3) :: CoordIN, CoordOUT 
    integer(8) :: year, day, hour, minute, secondINT

    real(8), intent(out) :: Pout(3)
    
    year = INT(year)
    day = INT(day)
    hour = INT(hour)
    minute = INT(minute)
    secondINT = INT(secondINT)
    secondTotal = real(secondTotal)

    call RECALC_08(year, day, hour, minute, secondINT, -500, 0, 0)
    
    call CoordinateTransform(CoordIN, CoordOUT, year, day, secondTotal, Pin, Pout)
    
    end subroutine CoordTrans

! **********************************************************************************************************************
! Subroutine FieldTrace:
!            subroutine that traces the magnetic field lines within different inputted
!            magnetic field models. The field lines are output in csv files named within a zip file.
! **********************************************************************************************************************
subroutine FieldTrace(PositionIN, Rigidity, Date, mode, IntMode, & 
    AtomicNumber, Anti, I, Wind, Pause, CoordSystem, GyroPercent, &
    End, FileName, gOTSO,hOTSO)
USE Particle
USE GEOPACK1
USE GEOPACK2
USE SolarWind
USE MagneticFieldFunctions
USE MagnetopauseFunctions
USE IntegrationFunctions
USE Magnetopause
USE CUSTOMGAUSS
implicit none

real(8) :: PositionIN(5), Rigidity, Date(6), End(3)
real(8) :: Wind(17), Re, GyroPercent, Pin(3), Pout(3)
real(8) :: Xnew(3), XnewConverted(3), Bfield(3)
integer(8) :: mode(2), IntMode, Anti, AtomicNumber
integer(4) :: I, Limit, Pause
character(len=3) :: CoordSystem
character(len=50) :: FileName
real(8) :: gOTSO(105), hOTSO(105)

Re = 6371.2
Limit = 0
Acount = 0
Result = 0
SubResult = 0
MaxGyroPercent = GyroPercent

if (mode(1) == 4) then
   Ginput = gOTSO
   Hinput = hOTSO
end if

call CreateParticle(PositionIN, Rigidity, Date, AtomicNumber, Anti, mode)

call initializeWind(Wind, I, mode)
call initializeCustomGauss(mode)

call MagneticFieldAssign(mode)
call MagnetopauseAssign(Pause)
call IntegrationAssign(IntMode)

open(unit=10,file=FileName,status='replace')
write(10,"(a)")"X,Y,Z,Bx,By,Bz"

call CoordinateTransform("GDZ", "GSM", year, day, secondTotal, Pin, Pout)

call MagFieldCheck(Pout, Bfield)

call CoordinateTransform("GSM", CoordSystem, year, day, secondTotal, Pout, XnewConverted)

if (model(1) == 4) then
    if (CoordSystem == "GEO") then
        XnewConverted = Xnew
    else
        call CoordinateTransform("GEO", CoordSystem, year, day, secondTotal, Xnew, XnewConverted)
    end if
end if

write(10,'(*(G0.6,:,","))') XnewConverted, Bfield



do while (Result == 0) 
call RK4_FieldTrace(Bfield)
call EscapeCheck()
Xnew(1) = XnewTemp(1)/1000
Xnew(2) = XnewTemp(2)/1000
Xnew(3) = XnewTemp(3)/1000

call CoordinateTransform("GSM", CoordSystem, year, day, secondTotal, Xnew, XnewConverted)

if (model(1) == 4) then
    if (CoordSystem == "GEO") then
        XnewConverted = Xnew
    else
        call CoordinateTransform("GEO", CoordSystem, year, day, secondTotal, Xnew, XnewConverted)
    end if
end if

IF ( DistanceTraveled/1000.0 > End(2)*Re) THEN
    Limit = 1
    EXIT
END IF

write(10,'(*(G0.6,:,","))') XnewConverted, Bfield

IF (Position(1) < End(1) ) THEN
    EXIT
END IF

IF (Result == 1)  THEN
    EXIT
END IF
end do
Close(10, STATUS='KEEP') 

end subroutine FieldTrace