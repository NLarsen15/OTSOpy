! ************************************************************************************************************************************
! LorentzRelativity.f95 - File containing the subroutine responsible for using the lorentz equation to calculate the force and
! acceleration that a charged particle experiences.
!
! ************************************************************************************************************************************
! subroutine Lorentz:
! Subroutine that calulates the acceleration that a charged particle moving at relativistic speeds 
! experiences within a magnetic field.
!
! INPUT: vnew - Velocity [m/s]
!        b - Magnetic field strength [nT]
!
! OUTPUT: Accel - Acceleration expereinced by CR [m/s^2]
! 
! ************************************************************************************************************************************
subroutine Lorentz(vnew, b, Q, M, Accel)
use SharedParameters
implicit none

real(8) :: vnew(3), vabs, b(3), Q, M
real(8) :: Accel(3), lam
   
vabs = (vnew(1)*vnew(1) + vnew(2)*vnew(2) + vnew(3)*vnew(3))**(0.5)

IF (vabs > c) THEN
    print *, "ERROR: Particle has exceeded light speed. Please check timestep value"
    stop
END IF

lam = (1 - ((vabs/c)**2))**(-0.5)

Accel(1) = (Q*(vnew(2)*b(3) - b(2)*vnew(3)))/(lam*M)
Accel(2) = (Q*(vnew(3)*b(1) - b(3)*vnew(1)))/(lam*M)
Accel(3) = (Q*(vnew(1)*b(2) - b(1)*vnew(2)))/(lam*M)
    
    
end subroutine Lorentz

! ************************************************************************************************************************************
! subroutine LorentzForce:
! Raw Lorentz force F = Q*(v x B) on a charged particle, i.e. dP/dt for the relativistic momentum P.
! Unlike Lorentz (above), this does NOT divide by gamma*M, since in a momentum-space integrator (e.g. RK6)
! the state derivative of P is the force itself, not an acceleration.
!
! INPUT:  v - Velocity [m/s]
!         B - Magnetic field strength [T]
!         Q - Charge [C]
!
! OUTPUT: F - Force [N]
!
! ************************************************************************************************************************************
subroutine LorentzForce(v, B, Q, F)
implicit none

real(8) :: v(3), B(3), Q
real(8) :: F(3)

F(1) = Q*(v(2)*B(3) - B(2)*v(3))
F(2) = Q*(v(3)*B(1) - B(3)*v(1))
F(3) = Q*(v(1)*B(2) - B(1)*v(2))

end subroutine LorentzForce


subroutine TimeCheck(Vabs, h, TimeElapsed)
USE SharedParameters
implicit none
real(8), intent(in) :: Vabs, h
real(8), intent(inout) :: TimeElapsed
real(8) :: lam

lam = (1 - ((Vabs/c)**2))**(-0.5)
    
TimeElapsed = TimeElapsed + h*lam        
        
end subroutine TimeCheck