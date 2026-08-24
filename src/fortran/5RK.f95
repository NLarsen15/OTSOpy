! ************************************************************************************************************************************
! subroutine RK5:
! 6-stage, 5th order explicit Runge-Kutta method 
!
! INPUT:
! Variables used from the particle module
!
! OUTPUT:
! Xnew - New position of CR [GSM coordinates]
! Vnew - New velocity of CR [m/s]
! XnewGDZ - New position of CR [GDZ coordinates]
!
! ************************************************************************************************************************************
subroutine RK5(VelocityArray, PositionArray, h, BetaError, &
               FinalStep, M, Q, secondTotal, mindistcheck, DistanceTraveled, &
               steps, TimeElapsed, counter, OLDPositionArray, OLDVelocityArray, &
               OLDsecondTotal, MDP, MaxGyroPercent, R, firsth)

use SharedParameters
implicit none

real(8), intent(inout) :: VelocityArray(2,3)
real(8), intent(inout) :: PositionArray(3,3)
real(8), intent(inout) :: h
real(8), intent(inout) :: secondTotal

logical, intent(inout) :: mindistcheck
real(8), intent(inout) :: DistanceTraveled
real(8), intent(inout) :: TimeElapsed
integer(8), intent(inout) :: steps
integer, intent(inout) :: counter

real(8), intent(in) :: BetaError
real(8), intent(in) :: M, Q
logical, intent(in) :: FinalStep

real(8), intent(inout) :: OLDPositionArray(3,3)
real(8), intent(inout) :: OLDVelocityArray(2,3)
real(8), intent(inout) :: OLDsecondTotal

real(8), intent(inout) :: MDP(3)

real(8), intent(in) :: MaxGyroPercent
real(8), intent(in) :: R
real(8), intent(in) :: firsth

real(8) :: d
real(8) :: B0(3), crossed(3)

real(8) :: Vs1(3), Vs2(3), Vs3(3), Vs4(3), Vs5(3), Vs6(3)
real(8) :: K1(3), K2(3), K3(3), K4(3), K5(3), K6(3)

real(8) :: x0(3), x0_Re(3), v0(3)

real(8) :: Xnew(3), Vnew(3)
real(8) :: Xnew_Re(3)
real(8) :: XnewGDZ(3)

real(8) :: h_try
real(8) :: t0
real(8) :: Ym

real(8) :: v0MAG, Vend
real(8) :: Vabs, Verr, LOWVerr, SpeedChange
real(8) :: Max

real(8) :: Bnorm0
logical :: FieldRejected

logical :: StepAccepted

Verr    = BetaError / 100.0d0
LOWVerr = BetaError / 10000.0d0

if (h <= 0.0d0) then
    h = 1.0d-6
end if

if (FinalStep) then
    h_try = firsth
else
    h_try = h
end if

if (h_try <= 0.0d0) then
    h_try = 1.0d-6
end if

h = h_try

if (model(1) == 4 .or. model(1) == 1 .or. model(1) == 5) then
    x0_Re = PositionArray(2,:)
else
    x0_Re = PositionArray(3,:)
end if

v0 = VelocityArray(1,:)

x0 = x0_Re * Re_m

v0MAG = sqrt(dot_product(v0,v0))

StepAccepted = .false.

do while (.not. StepAccepted)

    if (h <= tiny(1.0d0)) then
        print *, "WARNING: RK5 timestep became extremely small."
        print *, "h = ", h
        return
    end if

    t0 = secondTotal

    call MagneticField(x0_Re, t0, B0)

    Bnorm0 = sqrt(dot_product(B0,B0))

    Ym = 1.0d0 / sqrt(1.0d0 - (v0MAG/c)**2)
    d  = Q / (Ym*M)

    !-------------------------------------------------------------
    ! Stage 1  (c1 = 0)
    !-------------------------------------------------------------
    Vs1 = v0
    call VecCross(Vs1, B0, crossed)
    K1 = h*d*crossed

    !-------------------------------------------------------------
    ! Stage 2  (c2 = 1/3)
    !-------------------------------------------------------------
    Vs2 = v0 + K1/3.0d0
    call VecCross(Vs2, B0, crossed)
    K2 = h*d*crossed

    !-------------------------------------------------------------
    ! Stage 3  (c3 = 2/5)
    !-------------------------------------------------------------
    Vs3 = v0 + (4.0d0*K1 + 6.0d0*K2)/25.0d0
    call VecCross(Vs3, B0, crossed)
    K3 = h*d*crossed

    !-------------------------------------------------------------
    ! Stage 4  (c4 = 1)
    !-------------------------------------------------------------
    Vs4 = v0 + (K1 - 12.0d0*K2 + 15.0d0*K3)/4.0d0
    call VecCross(Vs4, B0, crossed)
    K4 = h*d*crossed

    !-------------------------------------------------------------
    ! Stage 5  (c5 = 2/3)
    !-------------------------------------------------------------
    Vs5 = v0 + (6.0d0*K1 + 90.0d0*K2 - 50.0d0*K3 + 8.0d0*K4)/81.0d0
    call VecCross(Vs5, B0, crossed)
    K5 = h*d*crossed

    !-------------------------------------------------------------
    ! Stage 6  (c6 = 4/5)
    !-------------------------------------------------------------
    Vs6 = v0 + (6.0d0*K1 + 36.0d0*K2 + 10.0d0*K3 + 8.0d0*K4)/75.0d0
    call VecCross(Vs6, B0, crossed)
    K6 = h*d*crossed

    !-------------------------------------------------------------
    ! FINAL VELOCITY
    !
    ! b = [23/192, 0, 125/192, 0, -81/192, 125/192]
    !-------------------------------------------------------------

    Vnew = v0 + ( 23.0d0*K1 + 125.0d0*K3 - 81.0d0*K5 + 125.0d0*K6 ) / 192.0d0

    !-------------------------------------------------------------
    ! FINAL POSITION
    !
    ! X_new = X0 + h * sum(b_i * V_stage_i)
    !-------------------------------------------------------------

    Xnew = x0 + h * ( 23.0d0*Vs1 + 125.0d0*Vs3 - 81.0d0*Vs5 + 125.0d0*Vs6 ) / 192.0d0

    Vend = sqrt(dot_product(Vnew,Vnew))

    if (adaptivestep) then

        call FieldJumpCheck(Bnorm0, Xnew/Re_m, secondTotal+h, h, FieldRejected)

        if (FieldRejected) then
            cycle
        end if

    end if

    if (adaptivestep .and. v0MAG > tiny(1.0d0)) then

        SpeedChange = abs(Vend - v0MAG) / v0MAG

        if (SpeedChange > Verr) then
            h = 0.5d0*h
            cycle
        end if

    end if

    if (adaptivestep .and. Vend >= c) then
        h = 0.5d0*h
        cycle
    end if

    StepAccepted = .true.

end do

Xnew_Re = Xnew / Re_m

if (model(1) == 4 .or. model(1) == 1 .or. model(1) == 5) then
    call CoordinateTransform("GEO", "GDZ", &
                             year, day, secondTotal, &
                             Xnew_Re, XnewGDZ)
else
    call CoordinateTransform("GSM", "GDZ", &
                             year, day, secondTotal, &
                             Xnew_Re, XnewGDZ)
end if

if (trapdistcheck .and. .not. mindistcheck) then
    call mindistexamine(Xnew_Re, XnewGDZ, mindistcheck)
end if

Vabs = v0MAG

call TimeCheck(Vabs, h, TimeElapsed)

DistanceTraveled = DistanceTraveled + h * ( &
                    23.0d0*sqrt(dot_product(Vs1,Vs1)) + &
                    125.0d0*sqrt(dot_product(Vs3,Vs3)) - &
                    81.0d0*sqrt(dot_product(Vs5,Vs5)) + &
                    125.0d0*sqrt(dot_product(Vs6,Vs6)) ) / 192.0d0

if (counter == backsavelim) then

    call OldVariables(VelocityArray, PositionArray, secondTotal, &
                      OLDPositionArray, OLDVelocityArray, &
                      OLDsecondTotal)

    counter = 0

end if

secondTotal = secondTotal + h

VelocityArray(1,:) = Vnew

PositionArray(1,:) = XnewGDZ

if (model(1) == 4 .or. model(1) == 1 .or. model(1) == 5) then
    PositionArray(2,:) = Xnew_Re
else
    PositionArray(3,:) = Xnew_Re
end if

if (adaptivestep) then

    if (v0MAG > tiny(1.0d0)) then

        SpeedChange = abs(Vend - v0MAG) / v0MAG

        if (SpeedChange < LOWVerr) then
            h = 1.1d0*h
        end if

    end if

end if

if (.not. adaptivestep) then
    h = firsth
end if

if (adaptivestep) then

    call NewMax(VelocityArray, PositionArray, &
                MaxGyroPercent, secondTotal, R, Max)

    if (h > Max) then
        h = Max
    end if

end if

counter = counter + 1

return

end subroutine RK5
