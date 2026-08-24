module IntegrationFunctions
USE GEOPACK1
USE GEOPACK2
USE SolarWind
implicit none

procedure (funcIntegration), pointer :: IntegrationPointer => null ()

abstract interface
subroutine funcIntegration(DUMMY1, DUMMY2, &
 DUMMY5, DUMMY6, DUMMY7, DUMMY8, DUMMY9, &
 DUMMY10, DUMMY11, DUMMY12, DUMMY13, DUMMY14, DUMMY15, &
 DUMMY16, DUMMY17, DUMMY20, DUMMY21, &
 DUMMY22, DUMMY23, DUMMY24, DUMMY25, DUMMY26)
implicit none
real(8), intent(inout) :: DUMMY1(2,3) !VelocityArray
real(8), intent(inout) :: DUMMY2(3,3) !PositionArray
real(8), intent(inout) :: DUMMY5 !h (stepsize)
real(8), intent(in)    :: DUMMY6 !BetaError

logical, intent(in)    :: DUMMY7 !FinalStep logical

real(8), intent(in)    :: DUMMY8 !M (mass)
real(8), intent(in)    :: DUMMY9 !Q (charge)
real(8), intent(inout) :: DUMMY10 !totalSeconds
logical, intent(inout) :: DUMMY11 !mindistcheck

real(8), intent(inout) :: DUMMY12 !DistanceTraveled
integer(8), intent(inout) :: DUMMY13 !steps
real(8), intent(inout) :: DUMMY14 !TimeElapsed
integer, intent(inout) :: DUMMY15 !counter
real(8), intent(inout) :: DUMMY16(3,3) !OLDPositionArray
real(8), intent(inout) :: DUMMY17(2,3) !OLDVelocityArray
real(8), intent(inout) :: DUMMY20 !OLDsecondTotal
real(8), intent(inout) :: DUMMY21(3) !MDP
real(8), intent(in) :: DUMMY22 !MaxGyroPercent
real(8), intent(in) :: DUMMY23 !R
real(8), intent(in) :: DUMMY24 !firsth
real(8), intent(inout) :: DUMMY25(3) !CachedBfield
logical, intent(inout) :: DUMMY26 !CachedBfieldValid


end subroutine funcIntegration
end interface

contains

  subroutine function4RK(VelocityArray, PositionArray, &
    h, BetaError, FinalStep, M, Q, secondTotal, &
    mindistcheck, DistanceTraveled, steps, TimeElapsed, counter, &
    OLDPositionArray, OLDVelocityArray, &
    OLDsecondTotal, MDP, MaxGyroPercent, R, firsth, &
    CachedBfield, CachedBfieldValid) ! 4th Order Runge-Kutta
    real(8), intent(inout) :: VelocityArray(2,3), PositionArray(3,3)
    real(8), intent(inout) :: h, secondTotal
    logical, intent(inout) :: mindistcheck
    real(8), intent(in)    :: BetaError, M, Q
    logical, intent(in)    :: FinalStep
    real(8), intent(inout) :: DistanceTraveled, TimeElapsed
    integer(8), intent(inout) :: steps
    integer, intent(inout) :: counter
    real(8), intent(inout) :: OLDPositionArray(3,3)
    real(8), intent(inout) :: OLDVelocityArray(2,3)
    real(8), intent(inout) :: OLDsecondTotal
    real(8), intent(inout) :: MDP(3)
    real(8), intent(in) :: MaxGyroPercent
    real(8), intent(in) :: R
    real(8), intent(in) :: firsth
    real(8), intent(inout) :: CachedBfield(3)
    logical, intent(inout) :: CachedBfieldValid

    call RK4(VelocityArray, PositionArray, &
    h, BetaError, FinalStep, M, Q, secondTotal, &
    mindistcheck, DistanceTraveled, steps, TimeElapsed, counter, &
    OLDPositionArray, OLDVelocityArray, &
    OLDsecondTotal, MDP, MaxGyroPercent, R, firsth, &
    CachedBfield, CachedBfieldValid)
    return
  end subroutine function4RK

  subroutine function6RK(VelocityArray, PositionArray, &
    h, BetaError, FinalStep, M, Q, secondTotal, &
    mindistcheck, DistanceTraveled, steps, TimeElapsed, counter, &
    OLDPositionArray, OLDVelocityArray, &
    OLDsecondTotal, MDP, MaxGyroPercent, R, firsth, &
    CachedBfield, CachedBfieldValid) ! 6th Order Runge-Kutta
    real(8), intent(inout) :: VelocityArray(2,3), PositionArray(3,3)
    real(8), intent(inout) :: h, secondTotal
    logical, intent(inout) :: mindistcheck
    real(8), intent(in)    :: BetaError, M, Q
    logical, intent(in)    :: FinalStep
    real(8), intent(inout) :: DistanceTraveled, TimeElapsed
    integer(8), intent(inout) :: steps
    integer, intent(inout) :: counter
    real(8), intent(inout) :: OLDPositionArray(3,3)
    real(8), intent(inout) :: OLDVelocityArray(2,3)
    real(8), intent(inout) :: OLDsecondTotal
    real(8), intent(inout) :: MDP(3)
    real(8), intent(in) :: MaxGyroPercent
    real(8), intent(in) :: R
    real(8), intent(in) :: firsth
    real(8), intent(inout) :: CachedBfield(3)
    logical, intent(inout) :: CachedBfieldValid

    call RK6(VelocityArray, PositionArray, &
    h, BetaError, FinalStep, M, Q, secondTotal, &
    mindistcheck, DistanceTraveled, steps, TimeElapsed, counter, &
    OLDPositionArray, OLDVelocityArray, &
    OLDsecondTotal, MDP, MaxGyroPercent, R, firsth, &
    CachedBfield, CachedBfieldValid)
    return
  end subroutine function6RK

  subroutine function5RK(VelocityArray, PositionArray, &
    h, BetaError, FinalStep, M, Q, secondTotal, &
    mindistcheck, DistanceTraveled, steps, TimeElapsed, counter, &
    OLDPositionArray, OLDVelocityArray, &
    OLDsecondTotal, MDP, MaxGyroPercent, R, firsth, &
    CachedBfield, CachedBfieldValid) ! 5th Order Runge-Kutta (frozen-field, 6-stage)
    real(8), intent(inout) :: VelocityArray(2,3), PositionArray(3,3)
    real(8), intent(inout) :: h, secondTotal
    logical, intent(inout) :: mindistcheck
    real(8), intent(in)    :: BetaError, M, Q
    logical, intent(in)    :: FinalStep
    real(8), intent(inout) :: DistanceTraveled, TimeElapsed
    integer(8), intent(inout) :: steps
    integer, intent(inout) :: counter
    real(8), intent(inout) :: OLDPositionArray(3,3)
    real(8), intent(inout) :: OLDVelocityArray(2,3)
    real(8), intent(inout) :: OLDsecondTotal
    real(8), intent(inout) :: MDP(3)
    real(8), intent(in) :: MaxGyroPercent
    real(8), intent(in) :: R
    real(8), intent(in) :: firsth
    real(8), intent(inout) :: CachedBfield(3)
    logical, intent(inout) :: CachedBfieldValid

    call RK5(VelocityArray, PositionArray, &
    h, BetaError, FinalStep, M, Q, secondTotal, &
    mindistcheck, DistanceTraveled, steps, TimeElapsed, counter, &
    OLDPositionArray, OLDVelocityArray, &
    OLDsecondTotal, MDP, MaxGyroPercent, R, firsth, &
    CachedBfield, CachedBfieldValid)
    return
  end subroutine function5RK

  subroutine functionBorisBuneman(VelocityArray, PositionArray, &
    h, BetaError, FinalStep, M, Q, secondTotal, &
    mindistcheck, DistanceTraveled, steps, TimeElapsed, counter, &
    OLDPositionArray, OLDVelocityArray, &
    OLDsecondTotal, MDP, MaxGyroPercent, R, firsth, &
    CachedBfield, CachedBfieldValid) ! Boris-Buneman Method
    real(8), intent(inout) :: VelocityArray(2,3), PositionArray(3,3)
    real(8), intent(inout) :: h, secondTotal
    logical, intent(inout) :: mindistcheck
    real(8), intent(in)    :: BetaError, M, Q
    logical, intent(in)    :: FinalStep
    real(8), intent(inout) :: DistanceTraveled, TimeElapsed
    integer(8), intent(inout) :: steps
    integer, intent(inout) :: counter
    real(8), intent(inout) :: OLDPositionArray(3,3)
    real(8), intent(inout) :: OLDVelocityArray(2,3)
    real(8), intent(inout) :: OLDsecondTotal
    real(8), intent(inout) :: MDP(3)
    real(8), intent(in) :: MaxGyroPercent
    real(8), intent(in) :: R
    real(8), intent(in) :: firsth
    real(8), intent(inout) :: CachedBfield(3)
    logical, intent(inout) :: CachedBfieldValid

    call BorisBuneman(VelocityArray, PositionArray, &
    h, BetaError, FinalStep, M, Q, secondTotal, &
    mindistcheck, DistanceTraveled, steps, TimeElapsed, counter, &
    OLDPositionArray, OLDVelocityArray, &
    OLDsecondTotal, MDP, MaxGyroPercent, R, firsth, &
    CachedBfield, CachedBfieldValid)

    return
  end subroutine functionBorisBuneman

  subroutine functionVay(VelocityArray, PositionArray, &
    h, BetaError, FinalStep, M, Q, secondTotal, &
    mindistcheck, DistanceTraveled, steps, TimeElapsed, counter, &
    OLDPositionArray, OLDVelocityArray, &
    OLDsecondTotal, MDP, MaxGyroPercent, R, firsth, &
    CachedBfield, CachedBfieldValid) ! Vay Method
    real(8), intent(inout) :: VelocityArray(2,3), PositionArray(3,3)
    real(8), intent(inout) :: h, secondTotal
    logical, intent(inout) :: mindistcheck
    real(8), intent(in)    :: BetaError, M, Q
    logical, intent(in)    :: FinalStep
    real(8), intent(inout) :: DistanceTraveled, TimeElapsed
    integer(8), intent(inout) :: steps
    integer, intent(inout) :: counter
    real(8), intent(inout) :: OLDPositionArray(3,3)
    real(8), intent(inout) :: OLDVelocityArray(2,3)
    real(8), intent(inout) :: OLDsecondTotal
    real(8), intent(inout) :: MDP(3)
    real(8), intent(in) :: MaxGyroPercent
    real(8), intent(in) :: R
    real(8), intent(in) :: firsth
    real(8), intent(inout) :: CachedBfield(3)
    logical, intent(inout) :: CachedBfieldValid

    call Vay(VelocityArray, PositionArray, &
    h, BetaError, FinalStep, M, Q, secondTotal, &
    mindistcheck, DistanceTraveled, steps, TimeElapsed, counter, &
    OLDPositionArray, OLDVelocityArray, &
    OLDsecondTotal, MDP, MaxGyroPercent, R, firsth, &
    CachedBfield, CachedBfieldValid)
    return
  end subroutine functionVay

  subroutine functionHC(VelocityArray, PositionArray, &
    h, BetaError, FinalStep, M, Q, secondTotal, &
    mindistcheck, DistanceTraveled, steps, TimeElapsed, counter, &
    OLDPositionArray, OLDVelocityArray, &
    OLDsecondTotal, MDP, MaxGyroPercent, R, firsth, &
    CachedBfield, CachedBfieldValid) ! Higuera-Cary Method
    real(8), intent(inout) :: VelocityArray(2,3), PositionArray(3,3)
    real(8), intent(inout) :: h, secondTotal
    logical, intent(inout) :: mindistcheck
    real(8), intent(in)    :: BetaError, M, Q
    logical, intent(in)    :: FinalStep
    real(8), intent(inout) :: DistanceTraveled, TimeElapsed
    integer(8), intent(inout) :: steps
    integer, intent(inout) :: counter
    real(8), intent(inout) :: OLDPositionArray(3,3)
    real(8), intent(inout) :: OLDVelocityArray(2,3)
    real(8), intent(inout) :: OLDsecondTotal
    real(8), intent(inout) :: MDP(3)
    real(8), intent(in) :: MaxGyroPercent
    real(8), intent(in) :: R
    real(8), intent(in) :: firsth
    real(8), intent(inout) :: CachedBfield(3)
    logical, intent(inout) :: CachedBfieldValid

    call HC(VelocityArray, PositionArray, &
    h, BetaError, FinalStep, M, Q, secondTotal, &
    mindistcheck, DistanceTraveled, steps, TimeElapsed, counter, &
    OLDPositionArray, OLDVelocityArray, &
    OLDsecondTotal, MDP, MaxGyroPercent, R, firsth, &
    CachedBfield, CachedBfieldValid)
    return
  end subroutine functionHC

  subroutine IntegrationAssign(Intmode)
  implicit none
  integer(8), intent(in) :: IntMode
    
  IF (IntMode == 1) THEN
    IntegrationPointer => function4RK ! 4th Order Runge-Kutta
  ELSE IF (IntMode == 3) THEN
    IntegrationPointer => functionVay  ! Vay Method
  ELSE IF (IntMode == 4) THEN
    IntegrationPointer => functionHC  ! Higuera-Cary Method
  ELSE IF (IntMode == 5) THEN
    IntegrationPointer => function6RK  ! 6th Order Runge-Kutta
  ELSE IF (IntMode == 6) THEN
    IntegrationPointer => function5RK  ! 5th Order Runge-Kutta (frozen-field, 6-stage)
  ELSE IF (IntMode == 7) THEN
    IntegrationPointer => functionBorisBuneman  ! Boris-Buneman Method
  ELSE
    print *, "Please enter valid integration method"
  END IF
  
  
      
    end subroutine IntegrationAssign
   
 end module IntegrationFunctions