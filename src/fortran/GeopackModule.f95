! ************************************************************************************************************************************
! GeopackModule.f95 - Module file containing the GEOPACK1 and GEOPACK2 modules. This file is used to update some of the older
! fortran code within the magnetic field models provided by the IRBEM library by replacing the geopack COMMON blocks within them with these modules. 
! These modules make the magnetic field models compatable with the newer OTSO code.
!
! Values within these modules are asigned when calling the RECALC_08 subroutine from the GEOPACK_DP.f file.
!
! ************************************************************************************************************************************
module GEOPACK1
    implicit none
    save

    real(8) :: ST0, CT0, SL0, CL0, CTCL, STCL, CTSL, STSL, SFI, CFI, &
               SPS, CPS, DS3, CGST, SGST, PSI

    real(8) :: A11, A21, A31, A12, A22, A32, A13, A23, A33
    real(8) :: E11, E21, E31, E12, E22, E32, E13, E23, E33
    real(8) :: EE11, EE12, EE13, EE21, EE22, EE23, EE31, EE32, EE33
    real(8) :: L11, L12, L13, L21, L22, L23, L31, L32, L33

end module GEOPACK1

module GEOPACK2
    implicit none
    save

    real(8), dimension(105) :: G, H, REC

end module GEOPACK2