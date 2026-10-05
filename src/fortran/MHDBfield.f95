subroutine MHDField(InputPosition,secondTotal,outputBfield)
    USE Interpolation
    USE SharedParameters
    implicit none
    real(8) :: x_target, y_target, z_target, InputPosition(3),InputPositionTEMP(3)
    real(8) :: Bx_target, By_target, Bz_target, outputBfieldTemp(3)
    real(8) :: outputBfield(3),targetposition(3)
    real(8) :: secondTotal
    integer, allocatable :: x_values(:), y_values(:), z_values(:)
    character(len=3) :: CoordIN

    if (model(1) == 4 .or. model(1) == 1 .or. model(1) == 5) then
    CoordIN = "GEO"
    else
    CoordIN = "GSM"
    end if

    InputPositionTEMP(1) = InputPosition(1)
    InputPositionTEMP(2) = InputPosition(2)
    InputPositionTEMP(3) = InputPosition(3)

    call CoordinateTransform(CoordIN, CoordINMHD, year, day, secondTotal, InputPositionTEMP, targetposition)

    x_target = targetposition(1)
    y_target = targetposition(2)
    z_target = targetposition(3)

    call Interpolate(x_target, y_target, z_target, n_x, n_y, n_z, Bx_target,By_target,Bz_target)
    
    outputBfieldTemp(1) = Bx_target
    outputBfieldTemp(2) = By_target
    outputBfieldTemp(3) = Bz_target

    ! Return B in the same frame the position came in (GEO for GEO-native internal fields, else GSM), so it can be
    ! summed with the internal field. (Previously always GSM via CoordOUTMHD, which was wrong with IGRF.)
    call CoordinateTransformVec(CoordINMHD, CoordIN, year, day, secondTotal, outputBfieldTemp, outputBfield)

end subroutine MHDField
  
  