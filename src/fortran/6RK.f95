! ************************************************************************************************************************************
! subroutine RK6:
! 7-stage, order-6 explicit Runge-Kutta
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
subroutine RK6(VelocityArray, PositionArray, h, BetaError, &
               FinalStep, M, Q, secondTotal, mindistcheck, DistanceTraveled, &
               steps, TimeElapsed, counter, OLDPositionArray, OLDVelocityArray, &
               OLDsecondTotal, MDP, MaxGyroPercent, R, firsth, &
               CachedBfield, CachedBfieldValid)

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

real(8), intent(inout) :: CachedBfield(3)
logical, intent(inout) :: CachedBfieldValid

real(8) :: P0(3), P1(3), P2(3), P3(3), P4(3), P5(3), P6(3), P7(3)
real(8) :: F1(3), F2(3), F3(3), F4(3), F5(3), F6(3), F7(3)
real(8) :: x1(3), x2(3), x3(3), x4(3), x5(3), x6(3), x7(3)
real(8) :: v1(3), v2(3), v3(3), v4(3), v5(3), v6(3), v7(3)
real(8) :: Bfield(3)

real(8) :: x1_Re(3), x2_Re(3), x3_Re(3), x4_Re(3)
real(8) :: x5_Re(3), x6_Re(3), x7_Re(3)

real(8) :: x0(3), x0_Re(3), v0(3)

real(8) :: Xnew(3), Vnew(3), Pnew(3)
real(8) :: Xnew_Re(3)
real(8) :: XnewGDZ(3)

real(8) :: h_try
real(8) :: t0, t1, t2, t3, t4, t5, t6, t7

real(8) :: Ym, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Ynew
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

Ym = 1.0d0 / sqrt(1.0d0 - (v0MAG/c)**2)

P0 = Ym * M * v0

StepAccepted = .false.

do while (.not. StepAccepted)

    if (h <= tiny(1.0d0)) then
        print *, "WARNING: RK6 timestep became extremely small."
        print *, "h = ", h
        return
    end if

    t0 = secondTotal

    !-------------------------------------------------------------
    ! Stage 1  (c1 = 0)
    !-------------------------------------------------------------
    x1    = x0
    P1    = P0
    x1_Re = x0_Re

    Y1 = sqrt(1.0d0 + dot_product(P1,P1)/(M*M*c*c))
    v1 = P1 / (Y1*M)

    if (CachedBfieldValid) then
        Bfield = CachedBfield
    else
        call MagneticField(x1_Re, t0, Bfield)
    end if
    call LorentzForce(v1, Bfield, Q, F1)

    Bnorm0 = sqrt(dot_product(Bfield,Bfield))

    !-------------------------------------------------------------
    ! Stage 2  (c2 = 1/3)
    !-------------------------------------------------------------
    t1 = t0 + h/3.0d0

    x2 = x0 + h*( v1 ) / 3.0d0
    P2 = P0 + h*( F1 ) / 3.0d0

    Y2 = sqrt(1.0d0 + dot_product(P2,P2)/(M*M*c*c))
    v2 = P2 / (Y2*M)

    x2_Re = x2 / Re_m
    call MagneticField(x2_Re, t1, Bfield)
    call LorentzForce(v2, Bfield, Q, F2)

    !-------------------------------------------------------------
    ! Stage 3  (c3 = 2/3)
    !-------------------------------------------------------------
    t2 = t0 + 2.0d0*h/3.0d0

    x3 = x0 + h*( (2.0d0/3.0d0)*v2 )
    P3 = P0 + h*( (2.0d0/3.0d0)*F2 )

    Y3 = sqrt(1.0d0 + dot_product(P3,P3)/(M*M*c*c))
    v3 = P3 / (Y3*M)

    x3_Re = x3 / Re_m
    call MagneticField(x3_Re, t2, Bfield)
    call LorentzForce(v3, Bfield, Q, F3)

    !-------------------------------------------------------------
    ! Stage 4  (c4 = 1/3)
    !-------------------------------------------------------------
    t3 = t0 + h/3.0d0

    x4 = x0 + h*( (1.0d0/12.0d0)*v1 + (1.0d0/3.0d0)*v2 - (1.0d0/12.0d0)*v3 )
    P4 = P0 + h*( (1.0d0/12.0d0)*F1 + (1.0d0/3.0d0)*F2 - (1.0d0/12.0d0)*F3 )

    Y4 = sqrt(1.0d0 + dot_product(P4,P4)/(M*M*c*c))
    v4 = P4 / (Y4*M)

    x4_Re = x4 / Re_m
    call MagneticField(x4_Re, t3, Bfield)
    call LorentzForce(v4, Bfield, Q, F4)

    !-------------------------------------------------------------
    ! Stage 5  (c5 = 1/2)
    !-------------------------------------------------------------
    t4 = t0 + 0.5d0*h

    x5 = x0 + h*( -(1.0d0/16.0d0)*v1 + (9.0d0/8.0d0)*v2 &
               -   (3.0d0/16.0d0)*v3 - (3.0d0/8.0d0)*v4 )
    P5 = P0 + h*( -(1.0d0/16.0d0)*F1 + (9.0d0/8.0d0)*F2 &
               -   (3.0d0/16.0d0)*F3 - (3.0d0/8.0d0)*F4 )

    Y5 = sqrt(1.0d0 + dot_product(P5,P5)/(M*M*c*c))
    v5 = P5 / (Y5*M)

    x5_Re = x5 / Re_m
    call MagneticField(x5_Re, t4, Bfield)
    call LorentzForce(v5, Bfield, Q, F5)

    !-------------------------------------------------------------
    ! Stage 6  (c6 = 1/2)
    !-------------------------------------------------------------
    t5 = t0 + 0.5d0*h

    x6 = x0 + h*( (9.0d0/8.0d0)*v2 - (3.0d0/8.0d0)*v3 &
               -  (3.0d0/4.0d0)*v4 + (1.0d0/2.0d0)*v5 )
    P6 = P0 + h*( (9.0d0/8.0d0)*F2 - (3.0d0/8.0d0)*F3 &
               -  (3.0d0/4.0d0)*F4 + (1.0d0/2.0d0)*F5 )

    Y6 = sqrt(1.0d0 + dot_product(P6,P6)/(M*M*c*c))
    v6 = P6 / (Y6*M)

    x6_Re = x6 / Re_m
    call MagneticField(x6_Re, t5, Bfield)
    call LorentzForce(v6, Bfield, Q, F6)

    !-------------------------------------------------------------
    ! Stage 7  (c7 = 1)
    !-------------------------------------------------------------
    t6 = t0 + h

    x7 = x0 + h*(  (9.0d0/44.0d0)*v1  - (9.0d0/11.0d0)*v2 &
               +  (63.0d0/44.0d0)*v3  + (18.0d0/11.0d0)*v4 &
               -  (16.0d0/11.0d0)*v6 )
    P7 = P0 + h*(  (9.0d0/44.0d0)*F1  - (9.0d0/11.0d0)*F2 &
               +  (63.0d0/44.0d0)*F3  + (18.0d0/11.0d0)*F4 &
               -  (16.0d0/11.0d0)*F6 )

    Y7 = sqrt(1.0d0 + dot_product(P7,P7)/(M*M*c*c))
    v7 = P7 / (Y7*M)

    x7_Re = x7 / Re_m
    t7 = t0 + h
    call MagneticField(x7_Re, t7, Bfield)
    call LorentzForce(v7, Bfield, Q, F7)

    !-------------------------------------------------------------
    ! FINAL RK6 POSITION, MOMENTUM AND VELOCITY
    !
    ! b = [11/120, 0, 27/40, 27/40, -4/15, -4/15, 11/120]
    !-------------------------------------------------------------

    Xnew = x0 + h * ( (11.0d0/120.0d0)*v1 + (27.0d0/40.0d0)*v3 &
                    + (27.0d0/40.0d0)*v4  - (4.0d0/15.0d0)*v5 &
                    - (4.0d0/15.0d0)*v6   + (11.0d0/120.0d0)*v7 )

    Pnew = P0 + h * ( (11.0d0/120.0d0)*F1 + (27.0d0/40.0d0)*F3 &
                    + (27.0d0/40.0d0)*F4  - (4.0d0/15.0d0)*F5 &
                    - (4.0d0/15.0d0)*F6   + (11.0d0/120.0d0)*F7 )

    Ynew = sqrt(1.0d0 + dot_product(Pnew,Pnew)/(M*M*c*c))
    Vnew = Pnew / (Ynew*M)

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
                    (11.0d0/120.0d0)*sqrt(dot_product(v1,v1)) + &
                    (27.0d0/40.0d0)*sqrt(dot_product(v3,v3)) + &
                    (27.0d0/40.0d0)*sqrt(dot_product(v4,v4)) - &
                    (4.0d0/15.0d0)*sqrt(dot_product(v5,v5)) - &
                    (4.0d0/15.0d0)*sqrt(dot_product(v6,v6)) + &
                    (11.0d0/120.0d0)*Vend )

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
                MaxGyroPercent, secondTotal, R, Max, CachedBfield)

    CachedBfieldValid = .true.

    if (h > Max) then
        h = Max
    end if

end if

counter = counter + 1

return

end subroutine RK6
