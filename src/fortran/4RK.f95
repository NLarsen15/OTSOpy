! ************************************************************************************************************************************
! subroutine RK4:
! 4th order Runge-Kutta method for solving the differential equations of motion for the CR's trajectory.
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
subroutine RK4(VelocityArray, PositionArray, h, BetaError, &
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

real(8) :: P0(3), P1(3), P2(3), P3(3), P4(3)
real(8) :: F1(3), F2(3), F3(3), F4(3)
real(8) :: x0(3), x1(3), x2(3), x3(3), x4(3)
real(8) :: v0(3), v1(3), v2(3), v3(3), v4(3)
real(8) :: Bfield(3)

real(8) :: x0_Re(3)
real(8) :: x1_Re(3)
real(8) :: x2_Re(3)
real(8) :: x3_Re(3)
real(8) :: x4_Re(3)

real(8) :: Xnew(3)
real(8) :: Vnew(3)
real(8) :: Pnew(3)

real(8) :: Xnew_Re(3)
real(8) :: XnewGDZ(3)

real(8) :: h_try

real(8) :: t0
real(8) :: t1
real(8) :: t2
real(8) :: t3

real(8) :: Ym, Y1, Y2, Y3, Y4, Ynew

real(8) :: v0MAG
real(8) :: Vend

real(8) :: Vabs
real(8) :: Verr
real(8) :: LOWVerr
real(8) :: SpeedChange

real(8) :: Bnorm0
logical :: FieldRejected

real(8) :: Max

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
        print *, "WARNING: RK4 timestep became extremely small."
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
    ! Stage 2  (c2 = 1/2)
    !-------------------------------------------------------------
    t1 = t0 + 0.5d0*h

    x2 = x0 + 0.5d0*h*v1
    P2 = P0 + 0.5d0*h*F1

    Y2 = sqrt(1.0d0 + dot_product(P2,P2)/(M*M*c*c))
    v2 = P2 / (Y2*M)

    x2_Re = x2 / Re_m
    call MagneticField(x2_Re, t1, Bfield)
    call LorentzForce(v2, Bfield, Q, F2)

    !-------------------------------------------------------------
    ! Stage 3  (c3 = 1/2)
    !-------------------------------------------------------------
    t2 = t0 + 0.5d0*h

    x3 = x0 + 0.5d0*h*v2
    P3 = P0 + 0.5d0*h*F2

    Y3 = sqrt(1.0d0 + dot_product(P3,P3)/(M*M*c*c))
    v3 = P3 / (Y3*M)

    x3_Re = x3 / Re_m
    call MagneticField(x3_Re, t2, Bfield)
    call LorentzForce(v3, Bfield, Q, F3)

    !-------------------------------------------------------------
    ! Stage 4  (c4 = 1)
    !-------------------------------------------------------------
    t3 = t0 + h

    x4 = x0 + h*v3
    P4 = P0 + h*F3

    Y4 = sqrt(1.0d0 + dot_product(P4,P4)/(M*M*c*c))
    v4 = P4 / (Y4*M)

    x4_Re = x4 / Re_m
    call MagneticField(x4_Re, t3, Bfield)
    call LorentzForce(v4, Bfield, Q, F4)

    !-------------------------------------------------------------
    ! FINAL RK4 POSITION, MOMENTUM AND VELOCITY
    !
    ! b = [1/6, 1/3, 1/3, 1/6]
    !-------------------------------------------------------------

    Xnew = x0 + h * ( v1 + 2.0d0*v2 + 2.0d0*v3 + v4 ) / 6.0d0

    Pnew = P0 + h * ( F1 + 2.0d0*F2 + 2.0d0*F3 + F4 ) / 6.0d0

    Ynew = sqrt(1.0d0 + dot_product(Pnew,Pnew)/(M*M*c*c))
    Vnew = Pnew / (Ynew*M)

    Vend = sqrt(dot_product(Vnew,Vnew))

    if (adaptivestep) then

        call FieldJumpCheck(Bnorm0, Xnew/Re_m, secondTotal+h, h, FieldRejected)

        if (FieldRejected) then
            cycle
        end if

    end if

    if (adaptivestep) then

        if (v0MAG > tiny(1.0d0)) then

            SpeedChange = abs(Vend - v0MAG) / v0MAG

            if (SpeedChange > Verr) then
                h = 0.5d0*h
                cycle
            end if

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
                    sqrt(dot_product(v1,v1)) + 2.0d0*sqrt(dot_product(v2,v2)) + &
                    2.0d0*sqrt(dot_product(v3,v3)) + Vend ) / 6.0d0

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

end subroutine RK4
