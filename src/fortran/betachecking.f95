    subroutine BetaComp(Velocity, Beta)
    use SharedParameters
    implicit none
    real(8), intent(in) :: Velocity(2,3) 
    real(8), intent(out) :: Beta

    Beta = sqrt(Velocity(1,1)**2 + Velocity(1,2)**2 + Velocity(1,3)**2) / c

    end subroutine BetaComp


    subroutine BetaCheck(OriginalBeta, Velocity, BetaError, BetaCheckResult)
    use SharedParameters
    implicit none
    real(8), intent(in) :: OriginalBeta, Velocity(2,3)
    real(8), intent(in) :: BetaError
    logical, intent(inout) :: BetaCheckResult

    real(8) :: CurrentBeta, CurrentError


    CurrentBeta = sqrt(Velocity(1,1)**2 + Velocity(1,2)**2 + Velocity(1,3)**2) / c
    CurrentError = abs((CurrentBeta - OriginalBeta)/OriginalBeta)*100

    IF (CurrentError > BetaError*100) THEN
        !print *, "Beta changed by more than ", BetaError*100, " percent from original value. Original Beta: ", OriginalBeta, " Current Beta: ", CurrentBeta
         BetaCheckResult = .True.
        IF (adaptivestep .eqv. .TRUE.) THEN
            CurrentGyro = CurrentGyro / 1.25
            !print *, "Reducing gyrosize"
        END IF

        IF (adaptivestep .eqv. .FALSE.) THEN
            print *, "ERROR: Beta change exceeded allowed error in fixed step mode."
            print *, "We recommend enabling adaptive step size or stricter error controls."
            stop
        END IF
    END IF
    end subroutine BetaCheck