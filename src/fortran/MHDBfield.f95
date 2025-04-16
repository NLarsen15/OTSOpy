subroutine MHDField(InputPosition,outputBfield)
    USE Interpolation
    USE particle
    implicit none
    real(8) :: x_target, y_target, z_target, InputPosition(3)
    real(8) :: Bx_target, By_target, Bz_target, outputBfieldTemp(3)
    real(8) :: outputBfield(3),targetposition(3)
    integer, allocatable :: x_values(:), y_values(:), z_values(:)
    character(len=3) :: CoordIN

    if (model(1) == 4) then
    CoordIN = "GEO"
    else
    CoordIN = "GSM"
    end if

    call CoordinateTransform(CoordIN, CoordINMHD, year, day, secondINT, InputPosition, targetposition)

    x_target = targetposition(1)
    y_target = targetposition(2)
    z_target = targetposition(3)

    call Interpolate(x_target, y_target, z_target, MHDposition, MHDB, n_x, n_y, n_z, Bx_target,By_target,Bz_target)
    
    outputBfieldTemp(1) = Bx_target
    outputBfieldTemp(2) = By_target
    outputBfieldTemp(3) = Bz_target

    
    call CoordinateTransformVec(CoordINMHD, CoordOUTMHD, year, day, secondINT, outputBfieldTemp, outputBfield)

end subroutine MHDField
  
  