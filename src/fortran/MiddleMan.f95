module MiddleMan
    use iso_fortran_env, only: real32, real64
    implicit none

type :: FortranData
sequence
    real(8) :: StartRigidity
    real(8) :: EndRigidity
    real(8) :: RigidityStep
    real(8) :: GyroPercent
    real(8) :: FixedStepSize
    real(8) :: sphere
    real(8) :: trapdist
    real(8) :: Berr
    real(8) :: transmissionRres

    integer(4) :: gaussianlength

    real :: PositionIN(5)
    real :: Date(6)
    real :: Wind(25)
    real :: end(4)

    integer(8) :: IntMode
    integer(8) :: AtomicNumber
    integer(8) :: Anti
    integer(4) :: IOPT
    integer(4) :: mode(4)
    integer(4) :: Pause
    integer(4) :: Rcomputation
    integer(4) :: scanchoice
    integer(4) :: FortranThreads
    integer(4) :: n
    integer(4) :: transmissionsamples
    integer(4) :: maxdegree

    character(len=3) :: CoordSystem
    character(len=3) :: MHDCoordSys
    character(len=3) :: inputcoord
    
    logical :: adapt
    logical :: totalbetacheck
    logical :: trapdistcheck

end type FortranData

type :: ParticleData

    ! Position and velocity
    real(8) :: Velocity(3) = 0.0d0
    real(8) :: GEOVelocity(3) = 0.0d0

    real(8) :: PositionArray(3,3) = 0.0d0 
    ! PositionArray(1,:) = Position
    ! PositionArray(2,:) = GEOPosition
    ! PositionArray(3,:) = GSMPosition

    real(8) :: OldPositionArray(3,3) = 0.0d0
    ! OLDPositionArray(1,:) = OLDPosition
    ! OLDPositionArray(2,:) = OLDGEOPosition
    ! OLDPositionArray(3,:) = OLDGSMPosition

    real(8) :: VelocityArray(2,3) = 0.0d0
    ! VelocityArray(1,:) = Velocity
    ! VelocityArray(2,:) = GEOVelocity

    real(8) :: OldVelocityArray(2,3) = 0.0d0
    ! OldVelocityArray(1,:) = OLDVelocity
    ! OldVelocityArray(2,:) = OLDGEOVelocity

    
    ! Particle properties
    real(8) :: M = 0.0d0 !mass
    real(8) :: Q = 0.0d0 !charge
    real(8) :: Z = 0.0d0 !charge number
    real(8) :: A = 0.0d0 !mass number

    real(8) :: Lat = 0.0d0 !asymptotic latitude
    real(8) :: Long = 0.0d0 !asymptotic longitude

    real(8) :: E_0 = 0.0d0 !rest energy

    real(8) :: R = 0.0d0 !rigidity
    real(8) :: lambda = 0.0d0 !relativistic factor

    ! Time
    real(8) :: secondTotal = 0.0d0
    real(8) :: OLDsecondTotal = 0.0d0
    real(8) :: TimeElapsed = 0.0d0

    ! Integration
    real(8) :: h = 0.0d0
    real(8) :: hOLD = 0.0d0
    real(8) :: Lasth = 1E-4
    real(8) :: firsth = 1E-6
    real(8) :: MaxGyroPercent = 0.0d0

    ! Cached field from the previous step's NewMax evaluation, reused as this
    ! step's own starting-position field call when valid (same position/time).
    real(8) :: CachedBfield(3) = 0.0d0
    logical :: CachedBfieldValid = .false.

    real(8) :: DistanceTraveled = 0.0d0

    ! Minimum distance position
    real(8) :: MDP(3) = 0.0d0

    real(8) :: BetaError = 0.0d0
    real(8) :: OriginalBeta = 0.0d0
    real(8) :: CurrentBeta = 0.0d0

    ! Counters
    logical :: FinalStep = .false.
    logical :: mindistcheck = .false.
    logical :: Escaped = .false.
    logical :: TotalBetaCheckTrigger =.false.


    integer(8) :: steps = 0
    integer(4) :: counter = 0
    integer(4) :: Termtype = 0


end type ParticleData


    contains

! **********************************************************************************************************************
! Subroutine Cutoff:
!            subroutine that calculates the trajectory of a cosmic ray across a range of rigidities
!            within different input magnetic field models and determines the cutoff rigidity.
!            Will create a csv file in which the asymptotic cone data is stored.
!            Output can be the vertical cutoff or apparent cutoff
!            Apparent cutoff computation is roughly 9 times as long as vertical
!
! **********************************************************************************************************************
subroutine cutoff(Data, g8, h8, Rigidities, Allowed)

    USE omp_lib
    USE solarwind
    USE SharedParameters
    USE GEOPACK1
    USE CUSTOMGAUSS
    USE MagneticFieldFunctions
    USE MagnetopauseFunctions
    USE IntegrationFunctions
    implicit none

    type(ParticleData) :: Particle

    real(8), intent(out) :: Rigidities(:)
    integer(4), intent(out) :: Allowed(:)
    type(FortranData), intent(in) :: Data

    integer     :: loop, start_time, current_time, rate
    integer(4)  :: n
    real(8)     :: r, delay
    real(8)     :: Ri
    real(8)     :: Wind8(25), Date8(6), PositionIN8(5), End8(4)
    real(8)     :: g8(Data%gaussianlength), h8(Data%gaussianlength)
    integer     :: thread_id
    real(8)     :: test_value
    integer     :: i, l
    real(8)     :: x
    real(8)     :: BfieldFinal(3)
    real(8)     :: StartVelocity(3)

    !------------------------------------------------------------------
    ! Start timing
    !------------------------------------------------------------------
    Wind8            = real(Data%Wind,       kind=8)
    Date8            = real(Data%Date,       kind=8)
    PositionIN8      = real(Data%PositionIN, kind=8)
    End8             = real(Data%End,        kind=8)
    model            = Data%mode
    mintrapdist      = Data%trapdist
    totalbetacheck   = Data%totalbetacheck
    trapdistcheck    = Data%trapdistcheck
    adaptivestep     = Data%adapt
    year             = int(Data%Date(1))
    day              = int(Data%Date(2))
    hour             = int(Data%Date(3))
    minute           = int(Data%Date(4))
    secondINT        = int(Data%Date(5))
    Ginput           = g8
    Hinput           = h8
    Gaussianlen      = Data%gaussianlength
    degreemax        = Data%maxdegree
    spheresize       = Data%sphere
    n                = Data%n

    if (model(2) == 99) then
        CoordINMHD = Data%MHDCoordSys
        CoordOUTMHD = "GSM"
    end if

    call omp_set_num_threads(Data%FortranThreads)

    call initializeWind(Wind8, Data%IOPT, Data%mode)

    call initializeCustomGauss(model)

    call MagneticFieldAssign(model)

    call MagnetopauseAssign(Data%Pause)

    call IntegrationAssign(Data%IntMode)

    call AntiAssignCharge(Data%Anti)

    !$omp parallel do schedule(dynamic,1) &
    !$omp& private(loop, Particle, thread_id)
    do loop = 1, n

        CurrentGyro      = Data%GyroPercent

        thread_id = omp_get_thread_num()

10      Particle = ParticleData()  ! Initialize the ParticleData type

        Particle%R = real(loop, kind=8)*Data%RigidityStep + Data%EndRigidity

        Particle%MaxGyroPercent = CurrentGyro
        Particle%steps = 0
        Particle%DistanceTraveled = 0.0
        Particle%BetaError = Data%Berr
        Particle%FinalStep = .false.
        Particle%Termtype = 0
        Particle%Escaped = .false.

        call ParticleMass(Particle%R, Data%AtomicNumber, &
            Particle%M, Particle%Q, Particle%Z, &
            Particle%A, Particle%E_0, Particle%lambda)

        call ParticleVelocities(PositionIN8, Particle%R, Date8, Particle%PositionArray, &
            Particle%lambda, Particle%secondTotal, &
            Data%inputcoord, Particle%VelocityArray)

        call BetaComp(Particle%VelocityArray, Particle%OriginalBeta)
        
        call FirstTimeStep(Particle%PositionArray,Particle%VelocityArray, &
            Particle%secondTotal, Particle%MaxGyroPercent, Particle%R, &
            Particle%h, Particle%hOLD, Particle%firsth, Data%FixedStepSize, &
            Data%IntMode, Particle%M, Particle%Q)

        DO WHILE (Particle%Termtype == 0)

50          call IntegrationPointer(Particle%VelocityArray, Particle%PositionArray, &
                Particle%h, Particle%BetaError, Particle%FinalStep, &
                Particle%M, Particle%Q, Particle%secondTotal, Particle%mindistcheck, &
                Particle%DistanceTraveled, Particle%steps, Particle%TimeElapsed, Particle%counter, &
                Particle%OLDPositionArray, Particle%OLDVelocityArray, &
                Particle%OLDsecondTotal, Particle%MDP, Particle%MaxGyroPercent, Particle%R, Particle%firsth, &
                Particle%CachedBfield, Particle%CachedBfieldValid)

            if (totalbetacheck) then

                call betacheck(Particle%OriginalBeta, Particle%VelocityArray, Particle%Betaerror, &
                  Particle%TotalBetaCheckTrigger)

                if (Particle%TotalBetaCheckTrigger) then

                    Particle%TotalBetaCheckTrigger = .false.

                    GOTO 10
                end if

            end if

            call PausePointer(Particle%PositionArray, Particle%secondTotal, Particle%FinalStep, Particle%Escaped)

            call Termination_checks(End8, Particle%PositionArray, Particle%Escaped, &
                Particle%DistanceTraveled, Particle%steps, Particle%TimeElapsed, Particle%Termtype)

            if (Particle%Termtype .ne. 0) then
                 !if (Particle%Finalstep .neqv. .true.) then
                 !call FinalStepCheck(Particle%Finalstep, Particle%PositionArray, &
                 !                    Particle%VelocityArray, Particle%OldPositionArray, &
                 !                    Particle%OldVelocityArray, &
                 !                    Particle%termtype, Particle%secondTotal, &
                 !                    Particle%OldsecondTotal)
                 !GOTO 50
                 !end if
                 GOTO 100
            end if

100         if (Particle%Termtype .ne. 1) then
                    Allowed(loop) = Particle%Termtype !forbidden
                    rigidities(loop) = Particle%R

            else
                    Allowed(loop) = Particle%Termtype !allowed
                    rigidities(loop) = Particle%R

            endif

        end do

    end do
    !$omp end parallel do

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
subroutine cone(Data, g8, h8, Rigidities, Allowed, Asymlat, Asymlong)
    USE omp_lib
    USE solarwind
    USE SharedParameters
    USE GEOPACK1
    USE CUSTOMGAUSS
    USE MagneticFieldFunctions
    USE MagnetopauseFunctions
    USE IntegrationFunctions
    implicit none

    type(ParticleData) :: Particle

    real(8), intent(out) :: Rigidities(:)
    integer(4), intent(out) :: Allowed(:)
    real(8), intent(out) :: Asymlat(:)
    real(8), intent(out) :: Asymlong(:)
    type(FortranData), intent(in) :: Data

    integer     :: loop, start_time, current_time, rate
    integer(4)  :: n
    real(8)     :: r, delay
    real(8)     :: Ri
    real(8)     :: Wind8(25), Date8(6), PositionIN8(5), End8(4)
    real(8)     :: g8(Data%gaussianlength), h8(Data%gaussianlength)
    integer     :: thread_id
    real(8)     :: test_value
    integer     :: i, l
    real(8)     :: x
    real(8)     :: BfieldFinal(3)
    real(8)     :: StartVelocity(3)

    !------------------------------------------------------------------
    ! Start timing
    !------------------------------------------------------------------
    Wind8            = real(Data%Wind,       kind=8)
    Date8            = real(Data%Date,       kind=8)
    PositionIN8      = real(Data%PositionIN, kind=8)
    End8             = real(Data%End,        kind=8)
    model            = Data%mode
    mintrapdist      = Data%trapdist
    totalbetacheck   = Data%totalbetacheck
    trapdistcheck    = Data%trapdistcheck
    adaptivestep     = Data%adapt
    year             = int(Data%Date(1))
    day              = int(Data%Date(2))
    hour             = int(Data%Date(3))
    minute           = int(Data%Date(4))
    secondINT        = int(Data%Date(5))
    Ginput           = g8
    Hinput           = h8
    Gaussianlen      = Data%gaussianlength
    degreemax        = Data%maxdegree
    spheresize       = Data%sphere
    n                = Data%n

    if (model(2) == 99) then
        CoordINMHD = Data%MHDCoordSys
        CoordOUTMHD = "GSM"
    end if

    call omp_set_num_threads(Data%FortranThreads)

    call initializeWind(Wind8, Data%IOPT, Data%mode)
    call initializeCustomGauss(model)
    call MagneticFieldAssign(model)
    call MagnetopauseAssign(Data%Pause)
    call IntegrationAssign(Data%IntMode)

    call AntiAssignCharge(Data%Anti)

    !$omp parallel do schedule(dynamic,1) &
    !$omp& private(loop, Particle, thread_id)
    do loop = 1, n

        CurrentGyro      = Data%GyroPercent

        thread_id = omp_get_thread_num()

10      Particle = ParticleData()  ! Initialize the ParticleData type

        Particle%R = real(loop, kind=8)*Data%RigidityStep + Data%EndRigidity

        Particle%MaxGyroPercent = CurrentGyro
        Particle%steps = 0
        Particle%DistanceTraveled = 0.0
        Particle%BetaError = Data%Berr
        Particle%FinalStep = .false.
        Particle%Termtype = 0
        Particle%Escaped = .false.

        call ParticleMass(Particle%R, Data%AtomicNumber, &
            Particle%M, Particle%Q, Particle%Z, &
            Particle%A, Particle%E_0, Particle%lambda)

        call ParticleVelocities(PositionIN8, Particle%R, Date8, Particle%PositionArray, &
            Particle%lambda, Particle%secondTotal, &
            Data%inputcoord, Particle%VelocityArray)
        
        call BetaComp(Particle%VelocityArray, Particle%OriginalBeta)
        
        call FirstTimeStep(Particle%PositionArray,Particle%VelocityArray, &
            Particle%secondTotal, Particle%MaxGyroPercent, Particle%R, &
            Particle%h, Particle%hOLD, Particle%firsth, Data%FixedStepSize, &
            Data%IntMode, Particle%M, Particle%Q)

        DO WHILE (Particle%Termtype == 0)

50          call IntegrationPointer(Particle%VelocityArray, Particle%PositionArray, &
                Particle%h, Particle%BetaError, Particle%FinalStep, &
                Particle%M, Particle%Q, Particle%secondTotal, Particle%mindistcheck, &
                Particle%DistanceTraveled, Particle%steps, Particle%TimeElapsed, Particle%counter, &
                Particle%OLDPositionArray, Particle%OLDVelocityArray, &
                Particle%OLDsecondTotal, Particle%MDP, Particle%MaxGyroPercent, Particle%R, Particle%firsth, &
                Particle%CachedBfield, Particle%CachedBfieldValid)

            if (totalbetacheck) then

                call betacheck(Particle%OriginalBeta, Particle%VelocityArray, Particle%Betaerror, &
                  Particle%TotalBetaCheckTrigger)

                if (Particle%TotalBetaCheckTrigger) then

                    Particle%TotalBetaCheckTrigger = .false.

                    GOTO 10
                end if

            end if

            call PausePointer(Particle%PositionArray, Particle%secondTotal, Particle%FinalStep, Particle%Escaped)

            call Termination_checks(End8, Particle%PositionArray, Particle%Escaped, &
                Particle%DistanceTraveled, Particle%steps, Particle%TimeElapsed, Particle%Termtype)

            if (Particle%Termtype .ne. 0) then
                 if (Particle%Finalstep .neqv. .true.) then
                 call FinalStepCheck(Particle%Finalstep, Particle%PositionArray, &
                                     Particle%VelocityArray, Particle%OldPositionArray, &
                                     Particle%OldVelocityArray, &
                                     Particle%termtype, Particle%secondTotal, &
                                     Particle%OldsecondTotal)
                 GOTO 50
                 end if
                 GOTO 100
            end if

100         if (Particle%Termtype .ne. 1) then
                    Allowed(loop) = Particle%Termtype !forbidden
                    rigidities(loop) = Particle%R

                    call AsymptoticDirection(Particle%PositionArray, &
                    Particle%VelocityArray, &
                    Particle%secondTotal, Data%CoordSystem, &
                    Particle%Lat, Particle%Long)
                    
                    Asymlat(loop) = Particle%Lat
                    Asymlong(loop) = Particle%Long
            else
                    Allowed(loop) = Particle%Termtype !allowed
                    rigidities(loop) = Particle%R

                    call AsymptoticDirection(Particle%PositionArray, &
                    Particle%VelocityArray, &
                    Particle%secondTotal, Data%CoordSystem, &
                    Particle%Lat, Particle%Long)

                    Asymlat(loop) = Particle%Lat
                    Asymlong(loop) = Particle%Long
            endif

        end do

    end do
    !$omp end parallel do

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
subroutine trajectory_full(Data, g8, h8, Rigidity, TrajectoryFile, &
    TrajectoryFilelen, Filter, Alat, Along)

    USE omp_lib
    USE solarwind
    USE SharedParameters
    USE GEOPACK1
    USE CUSTOMGAUSS
    USE MagneticFieldFunctions
    USE MagnetopauseFunctions
    USE IntegrationFunctions
    implicit none

    type(ParticleData) :: Particle
    integer, intent(in) :: TrajectoryFilelen
    character(len=TrajectoryFilelen), intent(in) :: TrajectoryFile
    real(8), intent(in) :: Rigidity
    type(FortranData), intent(in) :: Data
    integer, intent(out) :: Filter
    real(8), intent(out) :: Alat, Along

    real(8)     :: Ri
    real(8)     :: Wind8(25), Date8(6), PositionIN8(5), End8(4)
    real(8)     :: g8(Data%gaussianlength), h8(Data%gaussianlength)
    integer     :: thread_id
    real(8)     :: test_value
    integer     :: i, l
    real(8)     :: x
    real(8)     :: BfieldFinal(3)
    real(8)     :: StartVelocity(3)
    integer     :: unit
    integer     :: ios
    real(8)     :: Xnew(3), XnewConverted(3)
    real(8)     :: Vnew(3), VnewConverted(3)

    character(len=256) :: iomsg

    unit = 10

    !------------------------------------------------------------------
    ! Start timing
    !------------------------------------------------------------------
    Wind8            = real(Data%Wind,       kind=8)
    Date8            = real(Data%Date,       kind=8)
    PositionIN8      = real(Data%PositionIN, kind=8)
    End8             = real(Data%End,        kind=8)
    model            = Data%mode
    mintrapdist      = Data%trapdist
    totalbetacheck   = Data%totalbetacheck
    trapdistcheck    = Data%trapdistcheck
    adaptivestep     = Data%adapt
    year             = int(Data%Date(1))
    day              = int(Data%Date(2))
    hour             = int(Data%Date(3))
    minute           = int(Data%Date(4))
    secondINT        = int(Data%Date(5))
    Ginput           = g8
    Hinput           = h8
    Gaussianlen      = Data%gaussianlength
    degreemax        = Data%maxdegree
    spheresize       = Data%sphere
    CurrentGyro      = Data%GyroPercent

10  continue

    open(newunit=unit, file=trim(TrajectoryFile), &
        status='old', position='append', action='write', &
        iostat=ios, iomsg=iomsg)

    if (ios /= 0) then
        print *, "Error opening file: ", trim(TrajectoryFile)
        print *, "IOSTAT =", ios
        print *, "IOMSG  =", trim(iomsg)
        stop
    end if

    if (model(2) == 99) then
        CoordINMHD = Data%MHDCoordSys
        CoordOUTMHD = "GSM"
    end if

    call omp_set_num_threads(Data%FortranThreads)

    call initializeWind(Wind8, Data%IOPT, Data%mode)

    call initializeCustomGauss(model)

    call MagneticFieldAssign(model)

    call MagnetopauseAssign(Data%Pause)

    call IntegrationAssign(Data%IntMode)

    call AntiAssignCharge(Data%Anti)

    Particle = ParticleData()  ! Initialize the ParticleData type

    Particle%R = real(dnint(Rigidity * 1.0d8), kind=8) / 1.0d8

    Particle%MaxGyroPercent = CurrentGyro
    Particle%steps = 0
    Particle%DistanceTraveled = 0.0
    Particle%BetaError = Data%Berr
    Particle%FinalStep = .false.
    Particle%Termtype = 0
    Particle%Escaped = .false.

    call ParticleMass(Particle%R, Data%AtomicNumber, &
        Particle%M, Particle%Q, Particle%Z, &
        Particle%A, Particle%E_0, Particle%lambda)

    call ParticleVelocities(PositionIN8, Particle%R, Date8, Particle%PositionArray, &
        Particle%lambda, Particle%secondTotal, &
        Data%inputcoord, Particle%VelocityArray)

    call BetaComp(Particle%VelocityArray, Particle%OriginalBeta)
    
    call FirstTimeStep(Particle%PositionArray,Particle%VelocityArray, &
        Particle%secondTotal, Particle%MaxGyroPercent, Particle%R, &
        Particle%h, Particle%hOLD, Particle%firsth, Data%FixedStepSize, &
        Data%IntMode, Particle%M, Particle%Q)

    DO WHILE (Particle%Termtype == 0)

50          call IntegrationPointer(Particle%VelocityArray, Particle%PositionArray, &
                 Particle%h, Particle%BetaError, Particle%FinalStep, &
                 Particle%M, Particle%Q, Particle%secondTotal, Particle%mindistcheck, &
                 Particle%DistanceTraveled, Particle%steps, Particle%TimeElapsed, Particle%counter, &
                 Particle%OLDPositionArray, Particle%OLDVelocityArray, &
                 Particle%OLDsecondTotal, Particle%MDP, Particle%MaxGyroPercent, Particle%R, Particle%firsth, &
                Particle%CachedBfield, Particle%CachedBfieldValid)

                 if (model(1) == 4 .or. model(1) == 1 .or. model(1) == 5) then
                     Xnew = Particle%PositionArray(2,:)
                     Vnew = Particle%VelocityArray(1,:)/1000   ! slot 1, not slot 2 -- slot 2 is never updated past init
                     call CoordinateTransform("GEO", Data%CoordSystem, year, day, Particle%secondTotal, Xnew, XnewConverted)
                     call CoordinateTransformVec("GEO", Data%CoordSystem, year, day, Particle%secondTotal, Vnew, VnewConverted)
                 else
                     Xnew = Particle%PositionArray(1,:)
                     Vnew = Particle%VelocityArray(1,:)/1000
                     call CoordinateTransform("GDZ", Data%CoordSystem, year, day, Particle%secondTotal, Xnew, XnewConverted)
                     call CoordinateTransformVec("GSM", Data%CoordSystem, year, day, Particle%secondTotal, Vnew, VnewConverted)
                
                 end if

                 if (totalbetacheck) then

                    call betacheck(Particle%OriginalBeta, Particle%VelocityArray, Particle%Betaerror, &
                    Particle%TotalBetaCheckTrigger)

                    if (Particle%TotalBetaCheckTrigger) then
                    close(unit)

                    open(unit, file=trim(TrajectoryFile), status='replace', action='write')

                    Particle%TotalBetaCheckTrigger = .false.

                    GOTO 10
                    end if

                end if


                write(unit,'(ES25.16,",",ES25.16,",",ES25.16,",",ES25.16,",",ES25.16,",",ES25.16)') &
                     XnewConverted(1), &
                     XnewConverted(2), &
                     XnewConverted(3), &
                     VnewConverted(1), &
                     VnewConverted(2), &
                     VnewConverted(3)

             call PausePointer(Particle%PositionArray, Particle%secondTotal, Particle%FinalStep, Particle%Escaped)

             call Termination_checks(End8, Particle%PositionArray, Particle%Escaped, &
                 Particle%DistanceTraveled, Particle%steps, Particle%TimeElapsed, Particle%Termtype)

    end do

    Filter = Particle%Termtype
    call AsymptoticDirection(Particle%PositionArray, &
        Particle%VelocityArray, &
        Particle%secondTotal, Data%CoordSystem, &
        Particle%Lat, Particle%Long)
    Alat = Particle%Lat
    Along = Particle%Long

end subroutine trajectory_full
!**********************************************************************************************************************
subroutine trajectory(Data, g8, h8, Rigidities, RigiditiesLen, &
    Allowed, Asymlat, Asymlong)
    USE omp_lib
    USE solarwind
    USE SharedParameters
    USE GEOPACK1
    USE CUSTOMGAUSS
    USE MagneticFieldFunctions
    USE MagnetopauseFunctions
    USE IntegrationFunctions
    implicit none

    type(ParticleData) :: Particle

    integer(4), intent(in) :: RigiditiesLen

    real(8), intent(out) :: Rigidities(RigiditiesLen)
    integer(4), intent(out) :: Allowed(RigiditiesLen)
    real(8), intent(out) :: Asymlat(RigiditiesLen)
    real(8), intent(out) :: Asymlong(RigiditiesLen)
    type(FortranData), intent(in) :: Data

    integer     :: loop, start_time, current_time, rate
    integer(4)  :: n
    real(8)     :: r, delay
    real(8)     :: Ri
    real(8)     :: Wind8(25), Date8(6), PositionIN8(5), End8(4)
    real(8)     :: g8(Data%gaussianlength), h8(Data%gaussianlength)
    integer     :: thread_id
    real(8)     :: test_value
    integer     :: i, l
    real(8)     :: x
    real(8)     :: BfieldFinal(3)
    real(8)     :: StartVelocity(3)

    !------------------------------------------------------------------
    ! Start timing
    !------------------------------------------------------------------
    Wind8            = real(Data%Wind,       kind=8)
    Date8            = real(Data%Date,       kind=8)
    PositionIN8      = real(Data%PositionIN, kind=8)
    End8             = real(Data%End,        kind=8)
    model            = Data%mode
    mintrapdist      = Data%trapdist
    totalbetacheck   = Data%totalbetacheck
    trapdistcheck    = Data%trapdistcheck
    adaptivestep     = Data%adapt
    year             = int(Data%Date(1))
    day              = int(Data%Date(2))
    hour             = int(Data%Date(3))
    minute           = int(Data%Date(4))
    secondINT        = int(Data%Date(5))
    Ginput           = g8
    Hinput           = h8
    Gaussianlen      = Data%gaussianlength
    degreemax        = Data%maxdegree
    spheresize       = Data%sphere
    n                = Data%n

    if (model(2) == 99) then
        CoordINMHD = Data%MHDCoordSys
        CoordOUTMHD = "GSM"
    end if

    call omp_set_num_threads(Data%FortranThreads)

    call initializeWind(Wind8, Data%IOPT, Data%mode)
    call initializeCustomGauss(model)
    call MagneticFieldAssign(model)
    call MagnetopauseAssign(Data%Pause)
    call IntegrationAssign(Data%IntMode)

    call AntiAssignCharge(Data%Anti)

    !$omp parallel do schedule(dynamic,1) &
    !$omp& private(loop, Particle, thread_id)
    do loop = 1, RigiditiesLen

        CurrentGyro      = Data%GyroPercent

        thread_id = omp_get_thread_num()

10      Particle = ParticleData()  ! Initialize the ParticleData type

        Particle%R = Rigidities(loop)

        Particle%MaxGyroPercent = CurrentGyro
        Particle%steps = 0
        Particle%DistanceTraveled = 0.0
        Particle%BetaError = Data%Berr
        Particle%FinalStep = .false.
        Particle%Termtype = 0
        Particle%Escaped = .false.

        call ParticleMass(Particle%R, Data%AtomicNumber, &
            Particle%M, Particle%Q, Particle%Z, &
            Particle%A, Particle%E_0, Particle%lambda)

        call ParticleVelocities(PositionIN8, Particle%R, Date8, Particle%PositionArray, &
            Particle%lambda, Particle%secondTotal, &
            Data%inputcoord, Particle%VelocityArray)

        call BetaComp(Particle%VelocityArray, Particle%OriginalBeta)
        
        call FirstTimeStep(Particle%PositionArray,Particle%VelocityArray, &
            Particle%secondTotal, Particle%MaxGyroPercent, Particle%R, &
            Particle%h, Particle%hOLD, Particle%firsth, Data%FixedStepSize, &
            Data%IntMode, Particle%M, Particle%Q)

        DO WHILE (Particle%Termtype == 0)

50          call IntegrationPointer(Particle%VelocityArray, Particle%PositionArray, &
                Particle%h, Particle%BetaError, Particle%FinalStep, &
                Particle%M, Particle%Q, Particle%secondTotal, Particle%mindistcheck, &
                Particle%DistanceTraveled, Particle%steps, Particle%TimeElapsed, Particle%counter, &
                Particle%OLDPositionArray, Particle%OLDVelocityArray, &
                Particle%OLDsecondTotal, Particle%MDP, Particle%MaxGyroPercent, Particle%R, Particle%firsth, &
                Particle%CachedBfield, Particle%CachedBfieldValid)

            if (totalbetacheck) then

                call betacheck(Particle%OriginalBeta, Particle%VelocityArray, Particle%Betaerror, &
                  Particle%TotalBetaCheckTrigger)

                if (Particle%TotalBetaCheckTrigger) then

                    Particle%TotalBetaCheckTrigger = .false.

                    GOTO 10
                end if
            end if

            call PausePointer(Particle%PositionArray, Particle%secondTotal, Particle%FinalStep, Particle%Escaped)

            call Termination_checks(End8, Particle%PositionArray, Particle%Escaped, &
                Particle%DistanceTraveled, Particle%steps, Particle%TimeElapsed, Particle%Termtype)

            if (Particle%Termtype .ne. 0) then
                 if (Particle%Finalstep .neqv. .true.) then
                 call FinalStepCheck(Particle%Finalstep, Particle%PositionArray, &
                                     Particle%VelocityArray, Particle%OldPositionArray, &
                                     Particle%OldVelocityArray, &
                                     Particle%termtype, Particle%secondTotal, &
                                     Particle%OldsecondTotal)
                 GOTO 50
                 end if
                 GOTO 100
            end if

100         if (Particle%Termtype .ne. 1) then
                    Allowed(loop) = Particle%Termtype !forbidden

                    call AsymptoticDirection(Particle%PositionArray, &
                    Particle%VelocityArray, &
                    Particle%secondTotal, Data%CoordSystem, &
                    Particle%Lat, Particle%Long)
                    
                    Asymlat(loop) = Particle%Lat
                    Asymlong(loop) = Particle%Long
            else
                    Allowed(loop) = Particle%Termtype !allowed

                    call AsymptoticDirection(Particle%PositionArray, &
                    Particle%VelocityArray, &
                    Particle%secondTotal, Data%CoordSystem, &
                    Particle%Lat, Particle%Long)

                    Asymlat(loop) = Particle%Lat
                    Asymlong(loop) = Particle%Long
            endif

        end do

    end do
    !$omp end parallel do

end subroutine trajectory
! **********************************************************************************************************************
! Subroutine Transmission:
!            subroutine that calculates the trajectory of a cosmic ray across a range of rigidities
!            within different input magnetic field models and determines the cutoff rigidity.
!            Will create a csv file in which the asymptotic cone data is stored.
!            Output can be the vertical cutoff or apparent cutoff
!            Apparent cutoff computation is roughly 9 times as long as vertical
!
! **********************************************************************************************************************
subroutine transmission(Data, g8, h8, Rigidities, Transmissions)

    USE omp_lib
    USE solarwind
    USE SharedParameters
    USE GEOPACK1
    USE CUSTOMGAUSS
    USE MagneticFieldFunctions
    USE MagnetopauseFunctions
    USE IntegrationFunctions
    implicit none

    type(ParticleData) :: Particle

    real(8), intent(out) :: Rigidities(:)
    real(8), intent(out) :: Transmissions(:)
    type(FortranData), intent(in) :: Data

    integer     :: loop, start_time, current_time, rate
    integer(4)  :: n, num
    real(8)     :: r, delay
    real(8)     :: Ri
    real(8)     :: Wind8(25), Date8(6), PositionIN8(5), End8(4)
    real(8)     :: g8(Data%gaussianlength), h8(Data%gaussianlength)
    integer     :: thread_id
    real(8)     :: test_value
    integer     :: i, l
    real(8)     :: x
    real(8)     :: BfieldFinal(3)
    real(8)     :: StartVelocity(3)
    real(8)     :: centralR
    real(8)     :: Rres
    real(8)     :: lowr
    real(8)     :: allowedcount

    !------------------------------------------------------------------
    ! Start timing
    !------------------------------------------------------------------
    Wind8            = real(Data%Wind,       kind=8)
    Date8            = real(Data%Date,       kind=8)
    PositionIN8      = real(Data%PositionIN, kind=8)
    End8             = real(Data%End,        kind=8)
    model            = Data%mode
    mintrapdist      = Data%trapdist
    totalbetacheck   = Data%totalbetacheck
    trapdistcheck    = Data%trapdistcheck
    adaptivestep     = Data%adapt
    year             = int(Data%Date(1))
    day              = int(Data%Date(2))
    hour             = int(Data%Date(3))
    minute           = int(Data%Date(4))
    secondINT        = int(Data%Date(5))
    Ginput           = g8
    Hinput           = h8
    Gaussianlen      = Data%gaussianlength
    degreemax        = Data%maxdegree
    spheresize       = Data%sphere
    n                = Data%n

    if (model(2) == 99) then
        CoordINMHD = Data%MHDCoordSys
        CoordOUTMHD = "GSM"
    end if

    Rres = (2.0*Data%transmissionRres) / real(Data%transmissionsamples, kind=8)

    call omp_set_num_threads(Data%FortranThreads)

    call initializeWind(Wind8, Data%IOPT, Data%mode)

    call initializeCustomGauss(model)

    call MagneticFieldAssign(model)

    call MagnetopauseAssign(Data%Pause)

    call IntegrationAssign(Data%IntMode)

    call AntiAssignCharge(Data%Anti)

    !$omp parallel do schedule(dynamic,1) &
    !$omp& private(loop, Particle, thread_id, centralR, lowr, i, allowedcount, num)
    do loop = 1, n

        thread_id = omp_get_thread_num()

        allowedcount = 0

        centralR = real(loop, kind=8)*Data%RigidityStep + Data%EndRigidity
        lowr = centralR - Data%transmissionRres
        if (lowr < 0.0) lowr = 0.0

        do i = 1, Data%transmissionsamples

            CurrentGyro      = Data%GyroPercent

10          Particle = ParticleData()  ! Initialize the ParticleData type

            Particle%R = lowr + (i-1)*Rres

            !print *, "Thread ", thread_id, " is processing transmission rigidity ",&
            ! centralR, " sample ", i
            !print *, "Particle rigidity: ", Particle%R
            !print *, "Rres: ", Rres
            !print *, "max rigidity: ", centralR + Data%transmissionRres
            !print *, "Transmission res: ", Data%transmissionRres

            if (Particle%R > centralR + Data%transmissionRres + 1E-10) then
            GOTO 200
            end if

            Particle%MaxGyroPercent = CurrentGyro
            Particle%steps = 0
            Particle%DistanceTraveled = 0.0
            Particle%BetaError = Data%Berr
            Particle%FinalStep = .false.
            Particle%Termtype = 0
            Particle%Escaped = .false.
    
            call ParticleMass(Particle%R, Data%AtomicNumber, &
                Particle%M, Particle%Q, Particle%Z, &
                Particle%A, Particle%E_0, Particle%lambda)
    
            call ParticleVelocities(PositionIN8, Particle%R, Date8, Particle%PositionArray, &
                Particle%lambda, Particle%secondTotal, &
                Data%inputcoord, Particle%VelocityArray)

            call BetaComp(Particle%VelocityArray, Particle%OriginalBeta)
            
            call FirstTimeStep(Particle%PositionArray,Particle%VelocityArray, &
                Particle%secondTotal, Particle%MaxGyroPercent, Particle%R, &
                Particle%h, Particle%hOLD, Particle%firsth, Data%FixedStepSize, &
                Data%IntMode, Particle%M, Particle%Q)

            DO WHILE (Particle%Termtype == 0)
    
    50          call IntegrationPointer(Particle%VelocityArray, Particle%PositionArray, &
                    Particle%h, Particle%BetaError, Particle%FinalStep, &
                    Particle%M, Particle%Q, Particle%secondTotal, Particle%mindistcheck, &
                    Particle%DistanceTraveled, Particle%steps, Particle%TimeElapsed, Particle%counter, &
                    Particle%OLDPositionArray, Particle%OLDVelocityArray, &
                    Particle%OLDsecondTotal, Particle%MDP, Particle%MaxGyroPercent, Particle%R, Particle%firsth, &
                Particle%CachedBfield, Particle%CachedBfieldValid)

                if (totalbetacheck) then

                call betacheck(Particle%OriginalBeta, Particle%VelocityArray, Particle%Betaerror, &
                  Particle%TotalBetaCheckTrigger)

                if (Particle%TotalBetaCheckTrigger) then

                    Particle%TotalBetaCheckTrigger = .false.

                    GOTO 10
                end if

                end if
    
                call PausePointer(Particle%PositionArray, Particle%secondTotal, Particle%FinalStep, Particle%Escaped)
    
                call Termination_checks(End8, Particle%PositionArray, Particle%Escaped, &
                    Particle%DistanceTraveled, Particle%steps, Particle%TimeElapsed, Particle%Termtype)
    
                if (Particle%Termtype .ne. 0) then
                     if (Particle%Finalstep .neqv. .true.) then
                     call FinalStepCheck(Particle%Finalstep, Particle%PositionArray, &
                                         Particle%VelocityArray, Particle%OldPositionArray, &
                                         Particle%OldVelocityArray, &
                                         Particle%termtype, Particle%secondTotal, &
                                         Particle%OldsecondTotal)
                     GOTO 50
                     end if
                     GOTO 100
                end if
    
    100         if (Particle%Termtype == 1) then
                        allowedcount = allowedcount + 1
                endif
    
            end do

        end do

    200 Rigidities(loop) = centralR
    
    if (i > Data%transmissionsamples) then
        num = Data%transmissionsamples
    end if
    Transmissions(loop) = allowedcount/real(num, kind=8)

    end do
    !$omp end parallel do

end subroutine transmission
! **********************************************************************************************************************
! **********************************************************************************************************************
! Subroutine MagStrength:
!            subroutine that will tell you the strength of the magnetic field at any given point within the 
!            magnetosphere. Output is in GSM coordinates.
!
! **********************************************************************************************************************
subroutine MagStrength(Pin, Data, CoordIN, CoordOUT, g8, h8, Bfield)
    USE SharedParameters
    USE SolarWind
    USE MagneticFieldFunctions
    USE GEOPACK1
    USE GEOPACK2
    USE CUSTOMGAUSS
    USE Interpolation
    implicit none

    type(FortranData), intent(in) :: Data
    
    real(8) :: Pin(3), Pout(3), Bfieldtemp(3)
    character(len = 3) :: CoordIN, CoordOUT

    real(8)     :: Wind8(25), Date8(6)
    real(8)     :: g8(Data%gaussianlength), h8(Data%gaussianlength)
    real(8)     :: BfieldFinal(3)
    real(8)     :: secondTotal

    real(8), intent(out) :: Bfield(3) 

    !------------------------------------------------------------------
    ! Start timing
    !------------------------------------------------------------------
    Wind8            = real(Data%Wind,       kind=8)
    Date8            = real(Data%Date,       kind=8)
    model            = Data%mode
    year             = int(Data%Date(1))
    day              = int(Data%Date(2))
    hour             = int(Data%Date(3))
    minute           = int(Data%Date(4))
    secondINT        = int(Data%Date(5))
    secondTotal      = real(Date8(6))
    Ginput           = g8
    Hinput           = h8
    Gaussianlen      = Data%gaussianlength
    degreemax        = Data%maxdegree

    if (model(2) == 99) then
        CoordINMHD = Data%MHDCoordSys
        CoordOUTMHD = "GSM"
    end if

    year = INT(Data%Date(1))
    day = INT(Data%Date(2))
    hour = INT(Data%Date(3))
    minute = INT(Data%Date(4))
    secondINT = INT(Data%Date(5))
    secondTotal = real(Date8(6))

    call initializeWind(Wind8, Data%IOPT, model)
    call initializeCustomGauss(model)

    call MagneticFieldAssign(model)

    if (model(1) == 4 .or. model(1) == 1 .or. model(1) == 5) then
        call CoordinateTransform(CoordIN, "GEO", year, day, secondTotal, Pin, Pout)
    else
        call CoordinateTransform(CoordIN, "GSM", year, day, secondTotal, Pin, Pout)
    end if

    call MagFieldCheck(Pout, secondTotal, Bfieldtemp)

    if (model(1) == 4 .or. model(1) == 1 .or. model(1) == 5) then
        call CoordinateTransformVec("GEO", CoordOUT, year, day, secondTotal, Bfieldtemp, Bfield)
    else
        call CoordinateTransformVec("GSM", CoordOUT, year, day, secondTotal, Bfieldtemp, Bfield)
    end if
    
    end subroutine MagStrength
! **********************************************************************************************************************
! **********************************************************************************************************************
! Subroutine CoordTrans:
!            subroutine that uses the IBREM database of coordinate transforms to convert coordinates into
!            a new coordinate system.
!
! **********************************************************************************************************************
subroutine CoordTrans(Pin, Data, CoordIN, CoordOUT, g8, h8, Pout)
    USE SharedParameters
    USE CUSTOMGAUSS
    USE GEOPACK1
    USE GEOPACK2
    implicit none
    
    type(FortranData), intent(in) :: Data
    real(8) :: Pin(3)
    real(8) :: secondTotal
    character(len = 3) :: CoordIN, CoordOUT
    real(8)     :: g8(Data%gaussianlength), h8(Data%gaussianlength)
    real(8), intent(out) :: Pout(3)

    call initializeCustomGauss(model)

    Ginput = g8
    Hinput = h8

    year = INT(Data%Date(1))
    day = INT(Data%Date(2))
    hour = INT(Data%Date(3))
    minute = INT(Data%Date(4))
    secondINT = INT(Data%Date(5))
    secondTotal = real(Data%Date(6),kind=8)
    Gaussianlen      = Data%gaussianlength
    degreemax        = Data%maxdegree

    call RECALC_08(INT(year,kind=8), INT(day,kind=8), INT(hour,kind=8), INT(minute,kind=8), &
     INT(secondINT,kind=8), real(-500.0,kind=8), real(0.0,kind=8), real(0.0,kind=8))
    
    call CoordinateTransform(CoordIN, CoordOUT, year, day, secondTotal, Pin, Pout)
    
    end subroutine CoordTrans
! **********************************************************************************************************************
! **********************************************************************************************************************  
! **********************************************************************************************************************
! Subroutine FieldTrace:
!            subroutine that traces the magnetic field lines within different inputted
!            magnetic field models. The field lines are output in csv files named within a zip file.
! **********************************************************************************************************************
subroutine FieldTrace(Data, FileName, FileNamelen, g8, h8)
USE omp_lib
USE SharedParameters
USE GEOPACK1
USE GEOPACK2
USE CUSTOMGAUSS
USE SolarWind
USE MagneticFieldFunctions
USE MagnetopauseFunctions
USE IntegrationFunctions
implicit none

! -----------------------------
! Inputs
! -----------------------------
integer(4) :: FileNamelen
character(len=FileNamelen) :: FileName
type(ParticleData) :: Particle
type(FortranData), intent(in) :: Data

! -----------------------------
! Locals
! -----------------------------
integer(4) :: idir, nPlus, nMinus, inte, result
real(8) :: Bsign
real(8) :: Wind8(25), Date8(6), PositionIN8(5), End8(4)
real(8) :: g8(Data%gaussianlength), h8(Data%gaussianlength)

real(8), allocatable :: LinePlus(:,:), LineMinus(:,:)
real(8) :: Xnew(3), XnewConverted(3), Bfield(3)
integer :: io_unit

integer(4) :: MaxSteps
MaxSteps = 1000000   ! safety cap

allocate(LinePlus(MaxSteps,6))
allocate(LineMinus(MaxSteps,6))

!------------------------------------------------------------------
! Start timing
!------------------------------------------------------------------
Wind8            = real(Data%Wind,       kind=8)
Date8            = real(Data%Date,       kind=8)
PositionIN8      = real(Data%PositionIN, kind=8)
End8             = real(Data%End,        kind=8)
model            = Data%mode
mintrapdist      = Data%trapdist
totalbetacheck   = Data%totalbetacheck
trapdistcheck    = Data%trapdistcheck
adaptivestep     = Data%adapt
year             = int(Data%Date(1))
day              = int(Data%Date(2))
hour             = int(Data%Date(3))
minute           = int(Data%Date(4))
secondINT        = int(Data%Date(5))
Ginput           = g8
Hinput           = h8
Gaussianlen      = Data%gaussianlength
degreemax        = Data%maxdegree
spheresize       = Data%sphere

nPlus  = 0
nMinus = 0

! -----------------------------
! Initialization (unchanged)
! -----------------------------

call initializeWind(Wind8, Data%IOPT, Data%mode)
call initializeCustomGauss(model)
call MagneticFieldAssign(model)
call MagnetopauseAssign(Data%Pause)
call IntegrationAssign(Data%IntMode)

open(newunit=io_unit,file=FileName,status='old',position='append')

! -----------------------------
! Trace BOTH directions
! -----------------------------
do idir = 1, 2

    if (idir == 1) then
        Bsign = -1.0   ! toward south footprint
        nMinus = 0
    else
        Bsign =  1.0   ! toward north footprint
        nPlus = 0
    end if

    Particle = ParticleData()

    Particle%R = 1

    call ParticleVelocities(PositionIN8, Particle%R, Date8, Particle%PositionArray, &
            Particle%lambda, Particle%secondTotal, &
            Data%inputcoord, Particle%VelocityArray)
    Result = 0
    Particle%DistanceTraveled = 0.0d0

    do while (Result == 0)

        call Boris_FieldTrace_Advanced(Particle%PositionArray, Particle%secondTotal, &
        Particle%DistanceTraveled, Bsign, Bfield)
        
        call Termination_checks(End8, Particle%PositionArray, Particle%Escaped, &
                Particle%DistanceTraveled, Particle%steps, Particle%TimeElapsed, Particle%Termtype)


        Xnew = Particle%PositionArray(2,:)


        call CoordinateTransform("GSM", Data%CoordSystem, year, day, &
                                 Particle%secondTotal, &
                                 Xnew, XnewConverted)

        if (idir == 1) then
            nMinus = nMinus + 1
            if (nMinus <= MaxSteps) then
                LineMinus(nMinus,:) = (/ XnewConverted, Bfield*1E9 /)
            else
                exit
            end if
        else
            nPlus = nPlus + 1
            if (nPlus <= MaxSteps) then
                LinePlus(nPlus,:) = (/ XnewConverted, Bfield*1E9 /)
            else
                exit
            end if
        end if

        if (Particle%DistanceTraveled/1000.d0 > End8(2)*Re_km) exit
        if (Result == 1) exit
        if (Particle%PositionArray(1,1) < End8(1) ) exit

    end do
    Particle = ParticleData()
end do


do inte = min(nMinus,MaxSteps), 1, -1
    write(io_unit,'(*(G0.6,:,","))') LineMinus(inte,:)
    !print *, "nMinus = ", nMinus, "inte = ", inte
    !print *, LineMinus(inte,:)
end do

do inte = 2, min(nPlus,MaxSteps)
    write(io_unit,'(*(G0.6,:,","))') LinePlus(inte,:)
    !print *, "nPlus = ", nPlus, "inte = ", inte
    !print *, LinePlus(inte,:)
end do

close(io_unit)

deallocate(LinePlus, LineMinus)

end subroutine FieldTrace
! **********************************************************************************************************************








































subroutine MHDstartupSorted(XU, YU, ZU, MHDposition_in, MHDB_in, nx_split, ny_split, nz_split, &
                            mix,max,miy,may,miz,maz, &
                            region_order_in, start_x, end_x, start_y, end_y, start_z, end_z, &
                            num_regions,XUlen, YUlen, ZUlen)

  use Interpolation
  implicit none

  ! Inputs
  integer(4) :: XUlen, YUlen, ZUlen, num_regions
  real(8) :: XU(XUlen), YU(YUlen), ZU(ZUlen)
  real(8) :: MHDposition_in(XUlen, YUlen, ZUlen, 3)
  real(8) :: MHDB_in(XUlen, YUlen, ZUlen, 3)
  integer(4) :: region_order_in(num_regions)
  integer(4) :: start_x(num_regions), end_x(num_regions)
  integer(4) :: start_y(num_regions), end_y(num_regions)
  integer(4) :: start_z(num_regions), end_z(num_regions)
  integer :: nx_split,ny_split,nz_split,i,j,dx,dy,dz,idx,temp
  !integer :: search_range
  real :: mix,max,miy,may,miz,maz

  ! Save grid sizes
  n_x = XUlen
  n_y = YUlen
  n_z = ZUlen
  regions = num_regions

  MinX = mix
  MaxX = max
  MinY = miy
  MaxY = may
  MinZ = miz
  MaxZ = maz

  ! Allocate and copy fields
  allocate(MHDposition(n_x, n_y, n_z, 3))
  allocate(MHDB(n_x, n_y, n_z, 3))
  MHDposition = MHDposition_in
  MHDB = MHDB_in

  ! Store region bounds
  allocate(start_idx_x_region(regions)); start_idx_x_region = start_x
  allocate(end_idx_x_region(regions));   end_idx_x_region   = end_x
  allocate(start_idx_y_region(regions)); start_idx_y_region = start_y
  allocate(end_idx_y_region(regions));   end_idx_y_region   = end_y
  allocate(start_idx_z_region(regions)); start_idx_z_region = start_z
  allocate(end_idx_z_region(regions));   end_idx_z_region   = end_z

  ! Store region processing order
  allocate(region_order(regions)); region_order = region_order_in

  ! You could optionally compute resolution here too:
  x_res = abs(XU(2) - XU(1))
  y_res = abs(YU(2) - YU(1))
  z_res = abs(ZU(2) - ZU(1))

  n_x_split = nx_split
  n_y_split = ny_split
  n_z_split = nz_split

end subroutine MHDstartupSorted

subroutine gse2gswTSY15(date, position_gse, Wind, gOTSO, &
    hOTSO, glen, position_gsw)
    USE GEOPACK1
    USE GEOPACK2
    USE CUSTOMGAUSS
    implicit none
    
    real(8) :: sec, Pin(3), secondTotal
    real(8) :: Wind(3)
    real(8) :: XGSW, YGSW, ZGSW, XGSE, YGSE, ZGSE
    real(8) :: date(6), position_gse(3)
    real(8) :: gOTSO(glen), hOTSO(glen)
    integer(8) :: year, day, hour, minute, secondINT
    real(8), intent(inout) :: position_gsw(3)
    integer(4) :: maxdegree, glen

    degreemax = 13

    year = INT(Date(1))
    day = INT(Date(2))
    hour = INT(Date(3))
    minute = INT(Date(4))
    secondINT = INT(Date(5))
    secondTotal = real(Date(6))

    Ginput = gOTSO
    Hinput = hOTSO
    Gaussianlen = glen

    XGSE = position_gse(1)
    YGSE = position_gse(2)
    ZGSE = position_gse(3)

    call RECALC_08(year, day, hour, minute, secondINT, Wind(1), Wind(2), Wind(3))
    
    call GSWGSE_08(XGSW,YGSW,ZGSW,XGSE,YGSE,ZGSE,-1)

    position_gsw(1) = XGSW
    position_gsw(2) = YGSW      
    position_gsw(3) = ZGSW


end subroutine gse2gswTSY15


end module MiddleMan