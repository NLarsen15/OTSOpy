! ************************************************************************************************************************************
! Velocity.f95 - File responsible for determining the initial velocity of the CR by determining its speed and finding the velocity
! components along the x, y, and z axis using the normal vector at the origin point on Earth.
!
! ************************************************************************************************************************************
! subroutine Rigidity2Velocity:
! Subroutine used to calculate the velocity of the CR based on its rigidity value.
!
! INPUT:
! Uses variables within the particle module
!
! OUTPUT:
! V - Speed of the CR [m/s]
!
! ************************************************************************************************************************************
subroutine Rigidity2velocity(V)
USE particle
implicit none

real(8) :: V

V = c * (1.0-(1.0/(lambda**2.0)))**(0.5)

end subroutine Rigidity2velocity


! ************************************************************************************************************************************
! subroutine VelocityComponents:
! Subroutine used to calculate the x, y, and z axis components of the velocity.
!
! INPUT:
! Vabs - Speed of the CR [m/s]
! NormalVector - Normal vector at the given point on Earth
!
! OUTPUT:
! x, y, and z components of the velocty are updated in the particle module.
!
! ************************************************************************************************************************************
subroutine VelocityComponents(Vabs, NormalVector)
USE particle
implicit none

real(8) :: Vabs, NormalVector(3)
real(8) :: Normabs

Normabs = (NormalVector(1)**2 + NormalVector(2)**2 + NormalVector(3)**2)**(0.5)

Velocity(1) = Vabs*(NormalVector(1)/Normabs)
Velocity(2) = Vabs*(NormalVector(2)/Normabs)
Velocity(3) = Vabs*(NormalVector(3)/Normabs)

end subroutine VelocityComponents


! ************************************************************************************************************************************
! subroutine NormalVector:
! Subroutine used to calculate the normal vector at the given origin point for the CR.
!
! INPUT:
! StartPosition - Origin point on Earth for the CR [Geodetic coordinates (GDZ)]
!
! OUTPUT:
! NormOUT - Normal vector at the given point on Earth
!
! ************************************************************************************************************************************
subroutine NormalVector(EnteredPosition, inputcoord, NormOUT)
USE particle
implicit none
    
real(8) :: Earth, EnteredPosition(3)
real(8) :: NormOUT(3), xDT(3), xDTConvert(3), xINConvert(3)
real(8) :: rMag
character(len=3) :: inputcoord
    
!f2py intent(in) xIN, year, day, sec
!f2py intent(out) NormOUT
    
Earth = 6371200.0

if (inputcoord == "GDZ") then
    ! Geodetic: step +0.00001 Re along geodetic altitude
    xDT(1) = EnteredPosition(1) + 1
    xDT(2) = EnteredPosition(2)
    xDT(3) = EnteredPosition(3)
    if (model(1) == 4) then
        call CoordinateTransform("GDZ", "GEO", year, day, secondTotal, EnteredPosition, xINConvert)
        call CoordinateTransform("GDZ", "GEO", year, day, secondTotal, xDT, xDTConvert)
    else
        call CoordinateTransform("GDZ", "GSM", year, day, secondTotal, EnteredPosition, xINConvert)
        call CoordinateTransform("GDZ", "GSM", year, day, secondTotal, xDT, xDTConvert)
    end if

else if (inputcoord == "SPH") then
    ! Spherical: step +0.00001 Re along geocentric radius directly
    xDT(1) = EnteredPosition(1) + 0.00001
    xDT(2) = EnteredPosition(2)
    xDT(3) = EnteredPosition(3)
    if (model(1) == 4) then
        call CoordinateTransform("SPH", "GEO", year, day, secondTotal, EnteredPosition, xINConvert)
        call CoordinateTransform("SPH", "GEO", year, day, secondTotal, xDT, xDTConvert)
    else
        call CoordinateTransform("SPH", "GSM", year, day, secondTotal, EnteredPosition, xINConvert)
        call CoordinateTransform("SPH", "GSM", year, day, secondTotal, xDT, xDTConvert)
    end if

else
    ! Cartesian (GEO, GSM, SM, MAG, etc.): step +0.00001 Re along geocentric radial
    rMag = sqrt(EnteredPosition(1)**2 + EnteredPosition(2)**2 + EnteredPosition(3)**2)
    xDT(1) = EnteredPosition(1) * (1.0d0 + 0.00001d0 / rMag)
    xDT(2) = EnteredPosition(2) * (1.0d0 + 0.00001d0 / rMag)
    xDT(3) = EnteredPosition(3) * (1.0d0 + 0.00001d0 / rMag)
    if (model(1) == 4) then
        call CoordinateTransform(inputcoord, "GEO", year, day, secondTotal, EnteredPosition, xINConvert)
        call CoordinateTransform(inputcoord, "GEO", year, day, secondTotal, xDT, xDTConvert)
    else
        call CoordinateTransform(inputcoord, "GSM", year, day, secondTotal, EnteredPosition, xINConvert)
        call CoordinateTransform(inputcoord, "GSM", year, day, secondTotal, xDT, xDTConvert)
    end if

end if

NormOUT(1) = (xDTConvert(1) - xINConvert(1)) * Earth
NormOUT(2) = (xDTConvert(2) - xINConvert(2)) * Earth
NormOUT(3) = (xDTConvert(3) - xINConvert(3)) * Earth
    
end subroutine NormalVector