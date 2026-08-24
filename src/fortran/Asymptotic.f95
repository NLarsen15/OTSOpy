! ************************************************************************************************************************************
! Asymptotic.f95 - File responsible for calculating the asymptotic positions at a given point in the magnetosphere.
!
! ************************************************************************************************************************************
! subroutine AsymptoticDirection:
! Subroutine used to calculate the asymptotic latitude and longitude based on the coordinates of the CR.
!
! INPUT:
! Values from the particle module are used
!
! OUTPUT:
! lat - asymptotic latitude
! long - asymptotic longitude
!
! ************************************************************************************************************************************
subroutine AsymptoticDirection(PositionArray, VelocityArray, secondTotal, &
CoordSystem, Lat, Long)
use SharedParameters
implicit none
real(8) :: E, top, bottom, GSE(3)
real(8) :: THETAtemp, PHItemp, Alttemp, XGSEsph(3)
real(8) :: Xtemp, Ytemp, Ztemp, tempSPH(3), tempSPH2(3)
real(8) :: GSW(3), GEO(3), XGEO(3), theta, XGSE(3)
real(8) :: GEOsph(3), NewLat, tr, tTHETA, tPHI, tempposition(3), tempposition2(3)
real(8) :: Alttemp2, latnew, longnew
real(8), parameter :: pi  = 4 * atan(1.0_8)

character(len=3), intent(in) :: CoordSystem
real(8), intent(in) :: PositionArray(3,3), VelocityArray(2,3), secondTotal
real(8), intent(out) :: Lat, Long


call CoordinateTransform("GDZ", "GEO", year, day, secondTotal, PositionArray(1,:), XGEO)

call GSWGSM_08 (VelocityArray(1,1),VelocityArray(1,2),VelocityArray(1,3),GSW(1),GSW(2),GSW(3),-1)

call GEOGSW_08(GEO(1),GEO(2),GEO(3),GSW(1),GSW(2),GSW(3),-1)

call SPHCAR_08(GEOsph(1),GEOsph(2),GEOsph(3),XGEO(1), XGEO(2), XGEO(3), -1)

if (model(1) == 4 .or. model(1) == 1 .or. model(1) == 5) then
    GEO = VelocityArray(1,:)
end if

call BCARSP_08(XGEO(1),XGEO(2),XGEO(3),GEO(1),GEO(2),GEO(3),tr,tTHETA,tPHI)

theta = GEOsph(2)

top = -tTHETA*sin(theta) + tr*cos(theta)
bottom = (tPHI**2 + (tTHETA*cos(theta) + tr*sin(theta))**2)**(0.5)

Lat = atan2(top,bottom)

E = atan2(tPHI,(tTHETA*cos(theta) + tr*sin(theta)))

Long = GEOsph(3) + E

Long = Long * (180 / pi)
Lat = Lat * (180 / pi)

IF(Long > 360) THEN
    Long = Long - 360
ELSE IF(Long < 0) THEN
    Long = 360 + Long
END IF

IF (CoordSystem .NE. "GEO") THEN
THETAtemp = (pi/2) - Lat * (pi/180)
PHItemp = Long * (pi/180)
Alttemp = 1
tempSPH(1) = Alttemp
tempSPH(2) = THETAtemp
tempSPH(3) = PHItemp

CALL SPHCAR_08(Alttemp, THETAtemp, PHItemp, Xtemp, Ytemp, Ztemp, 1)
tempposition2(1) = Xtemp
tempposition2(2) = Ytemp
tempposition2(3) = Ztemp

call CoordinateTransform("GEO", CoordSystem, year, day, secondTotal, tempposition2, XGSE)

CALL SPHCAR_08(Alttemp, latnew, longnew, XGSE(1), XGSE(2), XGSE(3), -1)

Long = longnew * (180 / pi)
Lat  = 90 - (latnew  * (180 / pi))

IF(Long > 360) THEN
    Long = Long - 360
ELSE IF(Long < 0) THEN
    Long = 360 + Long
END IF

END IF

end subroutine AsymptoticDirection