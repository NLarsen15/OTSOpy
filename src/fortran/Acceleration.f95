! ************************************************************************************************************************************
! Acceleration.f95 - File containing subroutines used in the calculation of the CR's acceleration
!
! ************************************************************************************************************************************
! subroutine AccelerationCalc:
! Subroutine that calculates the acceleration of the CR.
!
! INPUT:
! xnew - Position of the CR [GSM coordinates]
! vnew - Velocity of the CR [m/s]
!
! OUTPUT:
! acceleration - Acceleration experienced by the CR [m/s^2]
!
! ************************************************************************************************************************************
subroutine AccelerationCalc(xnew, vnew, secondTotal, Q, M, acceleration)
implicit none

real(8) :: xnew(3), vnew(3), Bfield(3)
real(8) :: secondTotal, Q, M
real(8) :: acceleration(3)

call MagneticField(xnew, secondTotal, Bfield)

call Lorentz(vnew, Bfield, Q, M, acceleration)

end subroutine AccelerationCalc