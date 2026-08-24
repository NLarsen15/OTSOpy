! ************************************************************************************************************************************
! subroutine HC:
! Higuera-Cary method for solving the differential equations of motion for the CR's trajectory.
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
subroutine HC(VelocityArray, PositionArray, &
    h, BetaError, FinalStep, M, Q, secondTotal, &
    mindistcheck, DistanceTraveled, steps, TimeElapsed, counter, &
    OLDPositionArray, OLDVelocityArray, &
    OLDsecondTotal, MDP, MaxGyroPercent, R, firsth)

USE SharedParameters
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

real(8) :: Bfield(3)

real(8) :: xGSM(3)
real(8) :: x0(3)
real(8) :: Xnew(3)

real(8) :: Xnew_Re(3)
real(8) :: XnewGDZ(3)

real(8) :: v0(3)
real(8) :: Vnew(3)

real(8) :: Uminus(3)
real(8) :: Uplus(3)
real(8) :: UnNew(3)

real(8) :: tau(3)
real(8) :: tnew(3)
real(8) :: finalcalc1(3)

real(8) :: Fcalc2(3)
real(8) :: Fcalc3(3)
real(8) :: Fcalc4(3)
real(8) :: Fcalc5(3)

real(8) :: Vabs1
real(8) :: Vabs2

real(8) :: UminusDot
real(8) :: UplusDot
real(8) :: UnNewDot

real(8) :: Ustar
real(8) :: tauDot
real(8) :: tnewDot

real(8) :: lamminus
real(8) :: lamplus
real(8) :: lamnew

real(8) :: Sigma
real(8) :: Sfinal

real(8) :: scaler
real(8) :: h_used

real(8) :: Fcalc1

real(8) :: Verr
real(8) :: LOWVerr

real(8) :: Max

real(8) :: Bnorm0
logical :: FieldRejected

Verr = BetaError / 100.0d0
LOWVerr = BetaError / 10000.0d0

if (h <= 0.0d0) then
    h = 1.0d-6
end if

h_used = h

if (FinalStep) then
    h_used = Lasth
end if

if (model(1) == 4 .or. model(1) == 1 .or. model(1) == 5) then
    xGSM = PositionArray(2,:)
else
    xGSM = PositionArray(3,:)
end if

xGSM = xGSM * Re_m

x0 = xGSM

v0 = VelocityArray(1,:)

Vabs1 = sqrt(dot_product(v0,v0))

lamminus = 1.0d0 / sqrt(1.0d0 - (Vabs1/c)**2)

Uminus = lamminus * v0

10 continue

call MagneticField(x0 / Re_m, secondTotal, Bfield)

Bnorm0 = sqrt(dot_product(Bfield,Bfield))

scaler = (Q * h_used) / (2.0d0 * M)

call VecScale(scaler, Bfield, tau)

call VecDot(Uminus, tau, Ustar)

Ustar = Ustar / c

call VecDot(Uminus, Uminus, UminusDot)

lamminus = sqrt(1.0d0 + UminusDot/(c*c))

call VecDot(tau, tau, tauDot)

Sigma = lamminus*lamminus - tauDot

lamplus = sqrt( &
    (Sigma + &
    sqrt(Sigma*Sigma + &
    4.0d0*(tauDot + Ustar*Ustar))) / 2.0d0 )

call VecScale(1.0d0/lamplus, tau, tnew)

call VecDot(tnew, tnew, tnewDot)

Sfinal = 1.0d0 / (1.0d0 + tnewDot)

call VecDot(Uminus, tnew, Fcalc1)

call VecScale(Fcalc1, tnew, Fcalc2)

call VecCross(Uminus, tnew, Fcalc3)

call VecAddition(Fcalc2, Fcalc3, Fcalc4)

call VecAddition(Fcalc4, Uminus, Fcalc5)

call VecScale(Sfinal, Fcalc5, Uplus)

call VecCross(Uplus, tnew, finalcalc1)

call VecAddition(Uplus, finalcalc1, UnNew)

call VecDot(UnNew, UnNew, UnNewDot)

lamnew = sqrt(1.0d0 + UnNewDot/(c*c))

Vnew = UnNew / lamnew

Vabs2 = sqrt(dot_product(Vnew,Vnew))

Xnew = x0 + h_used*Vnew

Xnew_Re = Xnew / Re_m

call CoordinateTransform("GSM", "GDZ", &
                         year, day, &
                         secondTotal + h_used, &
                         Xnew_Re, XnewGDZ)

if (model(1) == 4 .or. model(1) == 1 .or. model(1) == 5) then

    call CoordinateTransform("GEO", "GDZ", &
                             year, day, &
                             secondTotal + h_used, &
                             Xnew_Re, XnewGDZ)

end if

if (trapdistcheck .and. .not. mindistcheck) then

    call mindistexamine(Xnew_Re, XnewGDZ, mindistcheck)

end if

if (adaptivestep .eqv. .TRUE.) then

    call FieldJumpCheck(Bnorm0, Xnew_Re, secondTotal + h_used, h_used, FieldRejected)

    if (FieldRejected) then

        if (.not. FinalStep) then
            h = h_used
        end if

        goto 10

    end if

end if

if (adaptivestep .eqv. .TRUE.) then

    if (Vabs2 >= c) then

        h_used = 0.5d0*h_used

        if (.not. FinalStep) then
            h = h_used
        end if

        goto 10

    end if

    if (((Vabs2 - Vabs1) / Vabs1) > Verr) then

        h_used = 0.5d0*h_used

        if (.not. FinalStep) then
            h = h_used
        end if

        goto 10

    end if

end if

secondTotal = secondTotal + h_used

call TimeCheck(Vabs1, h_used, TimeElapsed)

DistanceTraveled = DistanceTraveled + h_used*Vabs1

if (counter == backsavelim) then

    call OldVariables(VelocityArray, PositionArray, &
                      secondTotal, &
                      OLDPositionArray, &
                      OLDVelocityArray, &
                      OLDsecondTotal)

    counter = 0

end if

VelocityArray(1,:) = Vnew

PositionArray(1,:) = XnewGDZ

if (model(1) == 4 .or. model(1) == 1 .or. model(1) == 5) then
    PositionArray(2,:) = Xnew_Re
else
    PositionArray(3,:) = Xnew_Re
end if

if (adaptivestep .eqv. .TRUE.) then

    if (((Vabs2 - Vabs1) / Vabs1) < LOWVerr) then

        h = 1.1d0*h_used

    else

        h = h_used

    end if

end if

if (adaptivestep .eqv. .FALSE.) then

    h = firsth

end if

if (adaptivestep .eqv. .TRUE.) then

    call NewMax(VelocityArray, PositionArray, &
                MaxGyroPercent, &
                secondTotal, &
                R, &
                Max)

    if (h > Max) then
        h = Max
    end if

end if

if (FinalStep) then
    h = Lasth
end if

counter = counter + 1

end subroutine HC
