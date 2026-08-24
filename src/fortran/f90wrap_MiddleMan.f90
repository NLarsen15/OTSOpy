! Module middleman defined in file MiddleMan.f95

subroutine f90wrap_fortrandata__get__startrigidity(this, f90wrap_StartRigidity)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_StartRigidity
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_StartRigidity = this_ptr%p%StartRigidity
end subroutine f90wrap_fortrandata__get__startrigidity

subroutine f90wrap_fortrandata__set__startrigidity(this, f90wrap_StartRigidity)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_StartRigidity
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%StartRigidity = f90wrap_StartRigidity
end subroutine f90wrap_fortrandata__set__startrigidity

subroutine f90wrap_fortrandata__get__endrigidity(this, f90wrap_EndRigidity)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_EndRigidity
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_EndRigidity = this_ptr%p%EndRigidity
end subroutine f90wrap_fortrandata__get__endrigidity

subroutine f90wrap_fortrandata__set__endrigidity(this, f90wrap_EndRigidity)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_EndRigidity
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%EndRigidity = f90wrap_EndRigidity
end subroutine f90wrap_fortrandata__set__endrigidity

subroutine f90wrap_fortrandata__get__rigiditystep(this, f90wrap_RigidityStep)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_RigidityStep
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_RigidityStep = this_ptr%p%RigidityStep
end subroutine f90wrap_fortrandata__get__rigiditystep

subroutine f90wrap_fortrandata__set__rigiditystep(this, f90wrap_RigidityStep)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_RigidityStep
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%RigidityStep = f90wrap_RigidityStep
end subroutine f90wrap_fortrandata__set__rigiditystep

subroutine f90wrap_fortrandata__get__gyropercent(this, f90wrap_GyroPercent)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_GyroPercent
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_GyroPercent = this_ptr%p%GyroPercent
end subroutine f90wrap_fortrandata__get__gyropercent

subroutine f90wrap_fortrandata__set__gyropercent(this, f90wrap_GyroPercent)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_GyroPercent
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%GyroPercent = f90wrap_GyroPercent
end subroutine f90wrap_fortrandata__set__gyropercent

subroutine f90wrap_fortrandata__get__fixedstepsize(this, f90wrap_FixedStepSize)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_FixedStepSize
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_FixedStepSize = this_ptr%p%FixedStepSize
end subroutine f90wrap_fortrandata__get__fixedstepsize

subroutine f90wrap_fortrandata__set__fixedstepsize(this, f90wrap_FixedStepSize)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_FixedStepSize
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%FixedStepSize = f90wrap_FixedStepSize
end subroutine f90wrap_fortrandata__set__fixedstepsize

subroutine f90wrap_fortrandata__get__sphere(this, f90wrap_sphere)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_sphere
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_sphere = this_ptr%p%sphere
end subroutine f90wrap_fortrandata__get__sphere

subroutine f90wrap_fortrandata__set__sphere(this, f90wrap_sphere)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_sphere
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%sphere = f90wrap_sphere
end subroutine f90wrap_fortrandata__set__sphere

subroutine f90wrap_fortrandata__get__trapdist(this, f90wrap_trapdist)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_trapdist
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_trapdist = this_ptr%p%trapdist
end subroutine f90wrap_fortrandata__get__trapdist

subroutine f90wrap_fortrandata__set__trapdist(this, f90wrap_trapdist)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_trapdist
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%trapdist = f90wrap_trapdist
end subroutine f90wrap_fortrandata__set__trapdist

subroutine f90wrap_fortrandata__get__berr(this, f90wrap_Berr)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_Berr
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Berr = this_ptr%p%Berr
end subroutine f90wrap_fortrandata__get__berr

subroutine f90wrap_fortrandata__set__berr(this, f90wrap_Berr)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_Berr
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Berr = f90wrap_Berr
end subroutine f90wrap_fortrandata__set__berr

subroutine f90wrap_fortrandata__get__transmissionrres(this, f90wrap_transmissionRres)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_transmissionRres
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_transmissionRres = this_ptr%p%transmissionRres
end subroutine f90wrap_fortrandata__get__transmissionrres

subroutine f90wrap_fortrandata__set__transmissionrres(this, f90wrap_transmissionRres)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_transmissionRres
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%transmissionRres = f90wrap_transmissionRres
end subroutine f90wrap_fortrandata__set__transmissionrres

subroutine f90wrap_fortrandata__get__gaussianlength(this, f90wrap_gaussianlength)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_gaussianlength
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_gaussianlength = this_ptr%p%gaussianlength
end subroutine f90wrap_fortrandata__get__gaussianlength

subroutine f90wrap_fortrandata__set__gaussianlength(this, f90wrap_gaussianlength)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_gaussianlength
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%gaussianlength = f90wrap_gaussianlength
end subroutine f90wrap_fortrandata__set__gaussianlength

subroutine f90wrap_fortrandata__array__positionin(this, nd, dtype, dshape, dloc)
    use middleman
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in) :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%PositionIN)
    dloc = loc(this_ptr%p%PositionIN)
end subroutine f90wrap_fortrandata__array__positionin

subroutine f90wrap_fortrandata__array__date(this, nd, dtype, dshape, dloc)
    use middleman
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in) :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%Date)
    dloc = loc(this_ptr%p%Date)
end subroutine f90wrap_fortrandata__array__date

subroutine f90wrap_fortrandata__array__wind(this, nd, dtype, dshape, dloc)
    use middleman
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in) :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%Wind)
    dloc = loc(this_ptr%p%Wind)
end subroutine f90wrap_fortrandata__array__wind

subroutine f90wrap_fortrandata__array__end(this, nd, dtype, dshape, dloc)
    use middleman
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in) :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%end)
    dloc = loc(this_ptr%p%end)
end subroutine f90wrap_fortrandata__array__end

subroutine f90wrap_fortrandata__get__intmode(this, f90wrap_IntMode)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(8), intent(out) :: f90wrap_IntMode
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_IntMode = this_ptr%p%IntMode
end subroutine f90wrap_fortrandata__get__intmode

subroutine f90wrap_fortrandata__set__intmode(this, f90wrap_IntMode)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(8), intent(in) :: f90wrap_IntMode
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%IntMode = f90wrap_IntMode
end subroutine f90wrap_fortrandata__set__intmode

subroutine f90wrap_fortrandata__get__atomicnumber(this, f90wrap_AtomicNumber)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(8), intent(out) :: f90wrap_AtomicNumber
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_AtomicNumber = this_ptr%p%AtomicNumber
end subroutine f90wrap_fortrandata__get__atomicnumber

subroutine f90wrap_fortrandata__set__atomicnumber(this, f90wrap_AtomicNumber)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(8), intent(in) :: f90wrap_AtomicNumber
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%AtomicNumber = f90wrap_AtomicNumber
end subroutine f90wrap_fortrandata__set__atomicnumber

subroutine f90wrap_fortrandata__get__anti(this, f90wrap_Anti)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(8), intent(out) :: f90wrap_Anti
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Anti = this_ptr%p%Anti
end subroutine f90wrap_fortrandata__get__anti

subroutine f90wrap_fortrandata__set__anti(this, f90wrap_Anti)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(8), intent(in) :: f90wrap_Anti
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Anti = f90wrap_Anti
end subroutine f90wrap_fortrandata__set__anti

subroutine f90wrap_fortrandata__get__iopt(this, f90wrap_IOPT)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_IOPT
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_IOPT = this_ptr%p%IOPT
end subroutine f90wrap_fortrandata__get__iopt

subroutine f90wrap_fortrandata__set__iopt(this, f90wrap_IOPT)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_IOPT
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%IOPT = f90wrap_IOPT
end subroutine f90wrap_fortrandata__set__iopt

subroutine f90wrap_fortrandata__array__mode(this, nd, dtype, dshape, dloc)
    use middleman
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in) :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%mode)
    dloc = loc(this_ptr%p%mode)
end subroutine f90wrap_fortrandata__array__mode

subroutine f90wrap_fortrandata__get__pause(this, f90wrap_Pause)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_Pause
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Pause = this_ptr%p%Pause
end subroutine f90wrap_fortrandata__get__pause

subroutine f90wrap_fortrandata__set__pause(this, f90wrap_Pause)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_Pause
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Pause = f90wrap_Pause
end subroutine f90wrap_fortrandata__set__pause

subroutine f90wrap_fortrandata__get__rcomputation(this, f90wrap_Rcomputation)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_Rcomputation
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Rcomputation = this_ptr%p%Rcomputation
end subroutine f90wrap_fortrandata__get__rcomputation

subroutine f90wrap_fortrandata__set__rcomputation(this, f90wrap_Rcomputation)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_Rcomputation
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Rcomputation = f90wrap_Rcomputation
end subroutine f90wrap_fortrandata__set__rcomputation

subroutine f90wrap_fortrandata__get__scanchoice(this, f90wrap_scanchoice)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_scanchoice
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_scanchoice = this_ptr%p%scanchoice
end subroutine f90wrap_fortrandata__get__scanchoice

subroutine f90wrap_fortrandata__set__scanchoice(this, f90wrap_scanchoice)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_scanchoice
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%scanchoice = f90wrap_scanchoice
end subroutine f90wrap_fortrandata__set__scanchoice

subroutine f90wrap_fortrandata__get__fortranthreads(this, f90wrap_FortranThreads)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_FortranThreads
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_FortranThreads = this_ptr%p%FortranThreads
end subroutine f90wrap_fortrandata__get__fortranthreads

subroutine f90wrap_fortrandata__set__fortranthreads(this, f90wrap_FortranThreads)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_FortranThreads
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%FortranThreads = f90wrap_FortranThreads
end subroutine f90wrap_fortrandata__set__fortranthreads

subroutine f90wrap_fortrandata__get__n(this, f90wrap_n)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_n
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_n = this_ptr%p%n
end subroutine f90wrap_fortrandata__get__n

subroutine f90wrap_fortrandata__set__n(this, f90wrap_n)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_n
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%n = f90wrap_n
end subroutine f90wrap_fortrandata__set__n

subroutine f90wrap_fortrandata__get__transmissionsamples(this, f90wrap_transmissionsamples)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_transmissionsamples
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_transmissionsamples = this_ptr%p%transmissionsamples
end subroutine f90wrap_fortrandata__get__transmissionsamples

subroutine f90wrap_fortrandata__set__transmissionsamples(this, f90wrap_transmissionsamples)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_transmissionsamples
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%transmissionsamples = f90wrap_transmissionsamples
end subroutine f90wrap_fortrandata__set__transmissionsamples

subroutine f90wrap_fortrandata__get__maxdegree(this, f90wrap_maxdegree)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_maxdegree
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_maxdegree = this_ptr%p%maxdegree
end subroutine f90wrap_fortrandata__get__maxdegree

subroutine f90wrap_fortrandata__set__maxdegree(this, f90wrap_maxdegree)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_maxdegree
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%maxdegree = f90wrap_maxdegree
end subroutine f90wrap_fortrandata__set__maxdegree

subroutine f90wrap_fortrandata__get__coordsystem(this, f90wrap_CoordSystem)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    character*(3), intent(out) :: f90wrap_CoordSystem
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_CoordSystem = this_ptr%p%CoordSystem
end subroutine f90wrap_fortrandata__get__coordsystem

subroutine f90wrap_fortrandata__set__coordsystem(this, f90wrap_CoordSystem)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    character*(3), intent(in) :: f90wrap_CoordSystem
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%CoordSystem = f90wrap_CoordSystem
end subroutine f90wrap_fortrandata__set__coordsystem

subroutine f90wrap_fortrandata__get__mhdcoordsys(this, f90wrap_MHDCoordSys)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    character*(3), intent(out) :: f90wrap_MHDCoordSys
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_MHDCoordSys = this_ptr%p%MHDCoordSys
end subroutine f90wrap_fortrandata__get__mhdcoordsys

subroutine f90wrap_fortrandata__set__mhdcoordsys(this, f90wrap_MHDCoordSys)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    character*(3), intent(in) :: f90wrap_MHDCoordSys
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%MHDCoordSys = f90wrap_MHDCoordSys
end subroutine f90wrap_fortrandata__set__mhdcoordsys

subroutine f90wrap_fortrandata__get__inputcoord(this, f90wrap_inputcoord)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    character*(3), intent(out) :: f90wrap_inputcoord
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_inputcoord = this_ptr%p%inputcoord
end subroutine f90wrap_fortrandata__get__inputcoord

subroutine f90wrap_fortrandata__set__inputcoord(this, f90wrap_inputcoord)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    character*(3), intent(in) :: f90wrap_inputcoord
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%inputcoord = f90wrap_inputcoord
end subroutine f90wrap_fortrandata__set__inputcoord

subroutine f90wrap_fortrandata__get__adapt(this, f90wrap_adapt)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_adapt
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_adapt = this_ptr%p%adapt
end subroutine f90wrap_fortrandata__get__adapt

subroutine f90wrap_fortrandata__set__adapt(this, f90wrap_adapt)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_adapt
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%adapt = f90wrap_adapt
end subroutine f90wrap_fortrandata__set__adapt

subroutine f90wrap_fortrandata__get__totalbetacheck(this, f90wrap_totalbetacheck)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_totalbetacheck
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_totalbetacheck = this_ptr%p%totalbetacheck
end subroutine f90wrap_fortrandata__get__totalbetacheck

subroutine f90wrap_fortrandata__set__totalbetacheck(this, f90wrap_totalbetacheck)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_totalbetacheck
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%totalbetacheck = f90wrap_totalbetacheck
end subroutine f90wrap_fortrandata__set__totalbetacheck

subroutine f90wrap_fortrandata__get__trapdistcheck(this, f90wrap_trapdistcheck)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_trapdistcheck
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_trapdistcheck = this_ptr%p%trapdistcheck
end subroutine f90wrap_fortrandata__get__trapdistcheck

subroutine f90wrap_fortrandata__set__trapdistcheck(this, f90wrap_trapdistcheck)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(FortranData_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_trapdistcheck
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%trapdistcheck = f90wrap_trapdistcheck
end subroutine f90wrap_fortrandata__set__trapdistcheck

subroutine f90wrap_middleman__fortrandata_initialise(this)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    type(fortrandata_ptr_type) :: this_ptr
    integer(c_int), intent(out), dimension(4) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_middleman__fortrandata_initialise

subroutine f90wrap_middleman__fortrandata_finalise(this)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    type(fortrandata_ptr_type) :: this_ptr
    integer(c_int), intent(in), dimension(4) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_middleman__fortrandata_finalise

subroutine f90wrap_particledata__array__velocity(this, nd, dtype, dshape, dloc)
    use middleman
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in) :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%Velocity)
    dloc = loc(this_ptr%p%Velocity)
end subroutine f90wrap_particledata__array__velocity

subroutine f90wrap_particledata__array__geovelocity(this, nd, dtype, dshape, dloc)
    use middleman
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in) :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%GEOVelocity)
    dloc = loc(this_ptr%p%GEOVelocity)
end subroutine f90wrap_particledata__array__geovelocity

subroutine f90wrap_particledata__array__positionarray(this, nd, dtype, dshape, dloc)
    use middleman
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in) :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%PositionArray)
    dloc = loc(this_ptr%p%PositionArray)
end subroutine f90wrap_particledata__array__positionarray

subroutine f90wrap_particledata__array__oldpositionarray(this, nd, dtype, dshape, dloc)
    use middleman
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in) :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%OldPositionArray)
    dloc = loc(this_ptr%p%OldPositionArray)
end subroutine f90wrap_particledata__array__oldpositionarray

subroutine f90wrap_particledata__array__velocityarray(this, nd, dtype, dshape, dloc)
    use middleman
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in) :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%VelocityArray)
    dloc = loc(this_ptr%p%VelocityArray)
end subroutine f90wrap_particledata__array__velocityarray

subroutine f90wrap_particledata__array__oldvelocityarray(this, nd, dtype, dshape, dloc)
    use middleman
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in) :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%OldVelocityArray)
    dloc = loc(this_ptr%p%OldVelocityArray)
end subroutine f90wrap_particledata__array__oldvelocityarray

subroutine f90wrap_particledata__get__m(this, f90wrap_M)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_M
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_M = this_ptr%p%M
end subroutine f90wrap_particledata__get__m

subroutine f90wrap_particledata__set__m(this, f90wrap_M)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_M
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%M = f90wrap_M
end subroutine f90wrap_particledata__set__m

subroutine f90wrap_particledata__get__q(this, f90wrap_Q)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_Q
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Q = this_ptr%p%Q
end subroutine f90wrap_particledata__get__q

subroutine f90wrap_particledata__set__q(this, f90wrap_Q)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_Q
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Q = f90wrap_Q
end subroutine f90wrap_particledata__set__q

subroutine f90wrap_particledata__get__z(this, f90wrap_Z)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_Z
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Z = this_ptr%p%Z
end subroutine f90wrap_particledata__get__z

subroutine f90wrap_particledata__set__z(this, f90wrap_Z)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_Z
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Z = f90wrap_Z
end subroutine f90wrap_particledata__set__z

subroutine f90wrap_particledata__get__a(this, f90wrap_A)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_A
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_A = this_ptr%p%A
end subroutine f90wrap_particledata__get__a

subroutine f90wrap_particledata__set__a(this, f90wrap_A)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_A
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%A = f90wrap_A
end subroutine f90wrap_particledata__set__a

subroutine f90wrap_particledata__get__lat(this, f90wrap_Lat)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_Lat
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Lat = this_ptr%p%Lat
end subroutine f90wrap_particledata__get__lat

subroutine f90wrap_particledata__set__lat(this, f90wrap_Lat)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_Lat
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Lat = f90wrap_Lat
end subroutine f90wrap_particledata__set__lat

subroutine f90wrap_particledata__get__long_bn(this, f90wrap_Long)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_Long
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Long = this_ptr%p%Long
end subroutine f90wrap_particledata__get__long_bn

subroutine f90wrap_particledata__set__long_bn(this, f90wrap_Long)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_Long
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Long = f90wrap_Long
end subroutine f90wrap_particledata__set__long_bn

subroutine f90wrap_particledata__get__e_0(this, f90wrap_E_0)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_E_0
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_E_0 = this_ptr%p%E_0
end subroutine f90wrap_particledata__get__e_0

subroutine f90wrap_particledata__set__e_0(this, f90wrap_E_0)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_E_0
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%E_0 = f90wrap_E_0
end subroutine f90wrap_particledata__set__e_0

subroutine f90wrap_particledata__get__r(this, f90wrap_R)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_R
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_R = this_ptr%p%R
end subroutine f90wrap_particledata__get__r

subroutine f90wrap_particledata__set__r(this, f90wrap_R)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_R
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%R = f90wrap_R
end subroutine f90wrap_particledata__set__r

subroutine f90wrap_particledata__get__lambda_(this, f90wrap_lambda)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_lambda
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_lambda = this_ptr%p%lambda
end subroutine f90wrap_particledata__get__lambda_

subroutine f90wrap_particledata__set__lambda_(this, f90wrap_lambda)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_lambda
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%lambda = f90wrap_lambda
end subroutine f90wrap_particledata__set__lambda_

subroutine f90wrap_particledata__get__secondtotal(this, f90wrap_secondTotal)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_secondTotal
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_secondTotal = this_ptr%p%secondTotal
end subroutine f90wrap_particledata__get__secondtotal

subroutine f90wrap_particledata__set__secondtotal(this, f90wrap_secondTotal)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_secondTotal
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%secondTotal = f90wrap_secondTotal
end subroutine f90wrap_particledata__set__secondtotal

subroutine f90wrap_particledata__get__oldsecondtotal(this, f90wrap_OLDsecondTotal)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_OLDsecondTotal
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_OLDsecondTotal = this_ptr%p%OLDsecondTotal
end subroutine f90wrap_particledata__get__oldsecondtotal

subroutine f90wrap_particledata__set__oldsecondtotal(this, f90wrap_OLDsecondTotal)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_OLDsecondTotal
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%OLDsecondTotal = f90wrap_OLDsecondTotal
end subroutine f90wrap_particledata__set__oldsecondtotal

subroutine f90wrap_particledata__get__timeelapsed(this, f90wrap_TimeElapsed)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_TimeElapsed
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_TimeElapsed = this_ptr%p%TimeElapsed
end subroutine f90wrap_particledata__get__timeelapsed

subroutine f90wrap_particledata__set__timeelapsed(this, f90wrap_TimeElapsed)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_TimeElapsed
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%TimeElapsed = f90wrap_TimeElapsed
end subroutine f90wrap_particledata__set__timeelapsed

subroutine f90wrap_particledata__get__h(this, f90wrap_h)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_h
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_h = this_ptr%p%h
end subroutine f90wrap_particledata__get__h

subroutine f90wrap_particledata__set__h(this, f90wrap_h)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_h
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%h = f90wrap_h
end subroutine f90wrap_particledata__set__h

subroutine f90wrap_particledata__get__hold(this, f90wrap_hOLD)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_hOLD
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_hOLD = this_ptr%p%hOLD
end subroutine f90wrap_particledata__get__hold

subroutine f90wrap_particledata__set__hold(this, f90wrap_hOLD)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_hOLD
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%hOLD = f90wrap_hOLD
end subroutine f90wrap_particledata__set__hold

subroutine f90wrap_particledata__get__lasth(this, f90wrap_Lasth)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_Lasth
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Lasth = this_ptr%p%Lasth
end subroutine f90wrap_particledata__get__lasth

subroutine f90wrap_particledata__set__lasth(this, f90wrap_Lasth)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_Lasth
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Lasth = f90wrap_Lasth
end subroutine f90wrap_particledata__set__lasth

subroutine f90wrap_particledata__get__firsth(this, f90wrap_firsth)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_firsth
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_firsth = this_ptr%p%firsth
end subroutine f90wrap_particledata__get__firsth

subroutine f90wrap_particledata__set__firsth(this, f90wrap_firsth)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_firsth
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%firsth = f90wrap_firsth
end subroutine f90wrap_particledata__set__firsth

subroutine f90wrap_particledata__get__maxgyropercent(this, f90wrap_MaxGyroPercent)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_MaxGyroPercent
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_MaxGyroPercent = this_ptr%p%MaxGyroPercent
end subroutine f90wrap_particledata__get__maxgyropercent

subroutine f90wrap_particledata__set__maxgyropercent(this, f90wrap_MaxGyroPercent)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_MaxGyroPercent
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%MaxGyroPercent = f90wrap_MaxGyroPercent
end subroutine f90wrap_particledata__set__maxgyropercent

subroutine f90wrap_particledata__get__distancetraveled(this, f90wrap_DistanceTraveled)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_DistanceTraveled
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_DistanceTraveled = this_ptr%p%DistanceTraveled
end subroutine f90wrap_particledata__get__distancetraveled

subroutine f90wrap_particledata__set__distancetraveled(this, f90wrap_DistanceTraveled)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_DistanceTraveled
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%DistanceTraveled = f90wrap_DistanceTraveled
end subroutine f90wrap_particledata__set__distancetraveled

subroutine f90wrap_particledata__array__mdp(this, nd, dtype, dshape, dloc)
    use middleman
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in) :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%MDP)
    dloc = loc(this_ptr%p%MDP)
end subroutine f90wrap_particledata__array__mdp

subroutine f90wrap_particledata__get__betaerror(this, f90wrap_BetaError)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_BetaError
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_BetaError = this_ptr%p%BetaError
end subroutine f90wrap_particledata__get__betaerror

subroutine f90wrap_particledata__set__betaerror(this, f90wrap_BetaError)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_BetaError
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%BetaError = f90wrap_BetaError
end subroutine f90wrap_particledata__set__betaerror

subroutine f90wrap_particledata__get__originalbeta(this, f90wrap_OriginalBeta)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_OriginalBeta
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_OriginalBeta = this_ptr%p%OriginalBeta
end subroutine f90wrap_particledata__get__originalbeta

subroutine f90wrap_particledata__set__originalbeta(this, f90wrap_OriginalBeta)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_OriginalBeta
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%OriginalBeta = f90wrap_OriginalBeta
end subroutine f90wrap_particledata__set__originalbeta

subroutine f90wrap_particledata__get__currentbeta(this, f90wrap_CurrentBeta)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_CurrentBeta
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_CurrentBeta = this_ptr%p%CurrentBeta
end subroutine f90wrap_particledata__get__currentbeta

subroutine f90wrap_particledata__set__currentbeta(this, f90wrap_CurrentBeta)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_CurrentBeta
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%CurrentBeta = f90wrap_CurrentBeta
end subroutine f90wrap_particledata__set__currentbeta

subroutine f90wrap_particledata__get__finalstep(this, f90wrap_FinalStep)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_FinalStep
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_FinalStep = this_ptr%p%FinalStep
end subroutine f90wrap_particledata__get__finalstep

subroutine f90wrap_particledata__set__finalstep(this, f90wrap_FinalStep)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_FinalStep
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%FinalStep = f90wrap_FinalStep
end subroutine f90wrap_particledata__set__finalstep

subroutine f90wrap_particledata__get__mindistcheck(this, f90wrap_mindistcheck)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_mindistcheck
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_mindistcheck = this_ptr%p%mindistcheck
end subroutine f90wrap_particledata__get__mindistcheck

subroutine f90wrap_particledata__set__mindistcheck(this, f90wrap_mindistcheck)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_mindistcheck
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%mindistcheck = f90wrap_mindistcheck
end subroutine f90wrap_particledata__set__mindistcheck

subroutine f90wrap_particledata__get__escaped(this, f90wrap_Escaped)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_Escaped
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Escaped = this_ptr%p%Escaped
end subroutine f90wrap_particledata__get__escaped

subroutine f90wrap_particledata__set__escaped(this, f90wrap_Escaped)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_Escaped
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Escaped = f90wrap_Escaped
end subroutine f90wrap_particledata__set__escaped

subroutine f90wrap_particledata__get__totalbetachecktrigger(this, f90wrap_TotalBetaCheckTrigger)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_TotalBetaCheckTrigger
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_TotalBetaCheckTrigger = this_ptr%p%TotalBetaCheckTrigger
end subroutine f90wrap_particledata__get__totalbetachecktrigger

subroutine f90wrap_particledata__set__totalbetachecktrigger(this, f90wrap_TotalBetaCheckTrigger)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_TotalBetaCheckTrigger
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%TotalBetaCheckTrigger = f90wrap_TotalBetaCheckTrigger
end subroutine f90wrap_particledata__set__totalbetachecktrigger

subroutine f90wrap_particledata__get__steps(this, f90wrap_steps)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    integer(8), intent(out) :: f90wrap_steps
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_steps = this_ptr%p%steps
end subroutine f90wrap_particledata__get__steps

subroutine f90wrap_particledata__set__steps(this, f90wrap_steps)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    integer(8), intent(in) :: f90wrap_steps
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%steps = f90wrap_steps
end subroutine f90wrap_particledata__set__steps

subroutine f90wrap_particledata__get__counter(this, f90wrap_counter)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_counter
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_counter = this_ptr%p%counter
end subroutine f90wrap_particledata__get__counter

subroutine f90wrap_particledata__set__counter(this, f90wrap_counter)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_counter
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%counter = f90wrap_counter
end subroutine f90wrap_particledata__set__counter

subroutine f90wrap_particledata__get__termtype(this, f90wrap_Termtype)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_Termtype
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Termtype = this_ptr%p%Termtype
end subroutine f90wrap_particledata__get__termtype

subroutine f90wrap_particledata__set__termtype(this, f90wrap_Termtype)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    integer(c_int), intent(in)   :: this(4)
    type(ParticleData_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_Termtype
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Termtype = f90wrap_Termtype
end subroutine f90wrap_particledata__set__termtype

subroutine f90wrap_middleman__particledata_initialise(this)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    type(particledata_ptr_type) :: this_ptr
    integer(c_int), intent(out), dimension(4) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_middleman__particledata_initialise

subroutine f90wrap_middleman__particledata_finalise(this)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    
    type particledata_ptr_type
        type(particledata), pointer :: p => NULL()
    end type particledata_ptr_type
    type(particledata_ptr_type) :: this_ptr
    integer(c_int), intent(in), dimension(4) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_middleman__particledata_finalise

subroutine f90wrap_middleman__cutoff(f90wrap_n0, f90wrap_n1, f90wrap_n2, f90wrap_n3, data, g8, h8, rigidities, allowed)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer :: f90wrap_n0
    !f2py intent(hide), depend(g8) :: f90wrap_n0 = shape(g8,0)
    integer :: f90wrap_n1
    !f2py intent(hide), depend(h8) :: f90wrap_n1 = shape(h8,0)
    integer :: f90wrap_n2
    !f2py intent(hide), depend(rigidities) :: f90wrap_n2 = shape(rigidities,0)
    integer :: f90wrap_n3
    !f2py intent(hide), depend(allowed) :: f90wrap_n3 = shape(allowed,0)
    type(fortrandata_ptr_type) :: data_ptr
    integer(c_int), intent(in), dimension(4) :: data
    real(8), dimension(f90wrap_n0) :: g8
    real(8), dimension(f90wrap_n1) :: h8
    real(8), intent(inout), dimension(f90wrap_n2) :: rigidities
    integer(4), intent(inout), dimension(f90wrap_n3) :: allowed
    data_ptr = transfer(data, data_ptr)
    call cutoff(Data=data_ptr%p, g8=g8, h8=h8, Rigidities=rigidities, Allowed=allowed)
end subroutine f90wrap_middleman__cutoff

subroutine f90wrap_middleman__cone(f90wrap_n0, f90wrap_n1, f90wrap_n2, f90wrap_n3, f90wrap_n4, f90wrap_n5, data, g8, h8, &
    rigidities, allowed, asymlat, asymlong)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer :: f90wrap_n0
    !f2py intent(hide), depend(g8) :: f90wrap_n0 = shape(g8,0)
    integer :: f90wrap_n1
    !f2py intent(hide), depend(h8) :: f90wrap_n1 = shape(h8,0)
    integer :: f90wrap_n2
    !f2py intent(hide), depend(rigidities) :: f90wrap_n2 = shape(rigidities,0)
    integer :: f90wrap_n3
    !f2py intent(hide), depend(allowed) :: f90wrap_n3 = shape(allowed,0)
    integer :: f90wrap_n4
    !f2py intent(hide), depend(asymlat) :: f90wrap_n4 = shape(asymlat,0)
    integer :: f90wrap_n5
    !f2py intent(hide), depend(asymlong) :: f90wrap_n5 = shape(asymlong,0)
    type(fortrandata_ptr_type) :: data_ptr
    integer(c_int), intent(in), dimension(4) :: data
    real(8), dimension(f90wrap_n0) :: g8
    real(8), dimension(f90wrap_n1) :: h8
    real(8), intent(inout), dimension(f90wrap_n2) :: rigidities
    integer(4), intent(inout), dimension(f90wrap_n3) :: allowed
    real(8), intent(inout), dimension(f90wrap_n4) :: asymlat
    real(8), intent(inout), dimension(f90wrap_n5) :: asymlong
    data_ptr = transfer(data, data_ptr)
    call cone(Data=data_ptr%p, g8=g8, h8=h8, Rigidities=rigidities, Allowed=allowed, Asymlat=asymlat, Asymlong=asymlong)
end subroutine f90wrap_middleman__cone

subroutine f90wrap_middleman__trajectory_full(f90wrap_n0, f90wrap_n1, data, g8, h8, rigidity, trajectoryfile, &
    trajectoryfilelen, filter, alat, along)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer :: f90wrap_n0
    !f2py intent(hide), depend(g8) :: f90wrap_n0 = shape(g8,0)
    integer :: f90wrap_n1
    !f2py intent(hide), depend(h8) :: f90wrap_n1 = shape(h8,0)
    type(fortrandata_ptr_type) :: data_ptr
    integer(c_int), intent(in), dimension(4) :: data
    real(8), dimension(f90wrap_n0) :: g8
    real(8), dimension(f90wrap_n1) :: h8
    real(8), intent(in) :: rigidity
    character*(1024), intent(in) :: trajectoryfile
    integer, intent(in) :: trajectoryfilelen
    integer, intent(out) :: filter
    real(8), intent(out) :: alat
    real(8), intent(out) :: along
    data_ptr = transfer(data, data_ptr)
    call trajectory_full(Data=data_ptr%p, g8=g8, h8=h8, Rigidity=rigidity, TrajectoryFile=trajectoryfile, &
        TrajectoryFilelen=trajectoryfilelen, Filter=filter, Alat=alat, Along=along)
end subroutine f90wrap_middleman__trajectory_full

subroutine f90wrap_middleman__trajectory(f90wrap_n0, f90wrap_n1, f90wrap_n2, f90wrap_n3, f90wrap_n4, f90wrap_n5, data, &
    g8, h8, rigidities, rigiditieslen, allowed, asymlat, asymlong)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer :: f90wrap_n0
    !f2py intent(hide), depend(g8) :: f90wrap_n0 = shape(g8,0)
    integer :: f90wrap_n1
    !f2py intent(hide), depend(h8) :: f90wrap_n1 = shape(h8,0)
    integer :: f90wrap_n2
    !f2py intent(hide), depend(rigidities) :: f90wrap_n2 = shape(rigidities,0)
    integer :: f90wrap_n3
    !f2py intent(hide), depend(allowed) :: f90wrap_n3 = shape(allowed,0)
    integer :: f90wrap_n4
    !f2py intent(hide), depend(asymlat) :: f90wrap_n4 = shape(asymlat,0)
    integer :: f90wrap_n5
    !f2py intent(hide), depend(asymlong) :: f90wrap_n5 = shape(asymlong,0)
    type(fortrandata_ptr_type) :: data_ptr
    integer(c_int), intent(in), dimension(4) :: data
    real(8), dimension(f90wrap_n0) :: g8
    real(8), dimension(f90wrap_n1) :: h8
    real(8), intent(inout), dimension(f90wrap_n2) :: rigidities
    integer(4), intent(in) :: rigiditieslen
    integer(4), intent(inout), dimension(f90wrap_n3) :: allowed
    real(8), intent(inout), dimension(f90wrap_n4) :: asymlat
    real(8), intent(inout), dimension(f90wrap_n5) :: asymlong
    data_ptr = transfer(data, data_ptr)
    call trajectory(Data=data_ptr%p, g8=g8, h8=h8, Rigidities=rigidities, RigiditiesLen=rigiditieslen, Allowed=allowed, &
        Asymlat=asymlat, Asymlong=asymlong)
end subroutine f90wrap_middleman__trajectory

subroutine f90wrap_middleman__transmission(f90wrap_n0, f90wrap_n1, f90wrap_n2, f90wrap_n3, data, g8, h8, rigidities, &
    transmissions)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer :: f90wrap_n0
    !f2py intent(hide), depend(g8) :: f90wrap_n0 = shape(g8,0)
    integer :: f90wrap_n1
    !f2py intent(hide), depend(h8) :: f90wrap_n1 = shape(h8,0)
    integer :: f90wrap_n2
    !f2py intent(hide), depend(rigidities) :: f90wrap_n2 = shape(rigidities,0)
    integer :: f90wrap_n3
    !f2py intent(hide), depend(transmissions) :: f90wrap_n3 = shape(transmissions,0)
    type(fortrandata_ptr_type) :: data_ptr
    integer(c_int), intent(in), dimension(4) :: data
    real(8), dimension(f90wrap_n0) :: g8
    real(8), dimension(f90wrap_n1) :: h8
    real(8), intent(inout), dimension(f90wrap_n2) :: rigidities
    real(8), intent(inout), dimension(f90wrap_n3) :: transmissions
    data_ptr = transfer(data, data_ptr)
    call transmission(Data=data_ptr%p, g8=g8, h8=h8, Rigidities=rigidities, Transmissions=transmissions)
end subroutine f90wrap_middleman__transmission

subroutine f90wrap_middleman__magstrength(f90wrap_n0, f90wrap_n1, pin, data, coordin, coordout, g8, h8, bfield)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer :: f90wrap_n0
    !f2py intent(hide), depend(g8) :: f90wrap_n0 = shape(g8,0)
    integer :: f90wrap_n1
    !f2py intent(hide), depend(h8) :: f90wrap_n1 = shape(h8,0)
    real(8), dimension(3) :: pin
    type(fortrandata_ptr_type) :: data_ptr
    integer(c_int), intent(in), dimension(4) :: data
    character*(1024) :: coordin
    character*(1024) :: coordout
    real(8), dimension(f90wrap_n0) :: g8
    real(8), dimension(f90wrap_n1) :: h8
    real(8), dimension(3), intent(inout) :: bfield
    data_ptr = transfer(data, data_ptr)
    call MagStrength(Pin=pin, Data=data_ptr%p, CoordIN=coordin, CoordOUT=coordout, g8=g8, h8=h8, Bfield=bfield)
end subroutine f90wrap_middleman__magstrength

subroutine f90wrap_middleman__coordtrans(f90wrap_n0, f90wrap_n1, pin, data, coordin, coordout, g8, h8, pout)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer :: f90wrap_n0
    !f2py intent(hide), depend(g8) :: f90wrap_n0 = shape(g8,0)
    integer :: f90wrap_n1
    !f2py intent(hide), depend(h8) :: f90wrap_n1 = shape(h8,0)
    real(8), dimension(3) :: pin
    type(fortrandata_ptr_type) :: data_ptr
    integer(c_int), intent(in), dimension(4) :: data
    character*(1024) :: coordin
    character*(1024) :: coordout
    real(8), dimension(f90wrap_n0) :: g8
    real(8), dimension(f90wrap_n1) :: h8
    real(8), dimension(3), intent(inout) :: pout
    data_ptr = transfer(data, data_ptr)
    call CoordTrans(Pin=pin, Data=data_ptr%p, CoordIN=coordin, CoordOUT=coordout, g8=g8, h8=h8, Pout=pout)
end subroutine f90wrap_middleman__coordtrans

subroutine f90wrap_middleman__fieldtrace(f90wrap_n0, f90wrap_n1, data, filename, filenamelen, g8, h8)
    use middleman
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    
    type fortrandata_ptr_type
        type(fortrandata), pointer :: p => NULL()
    end type fortrandata_ptr_type
    integer :: f90wrap_n0
    !f2py intent(hide), depend(g8) :: f90wrap_n0 = shape(g8,0)
    integer :: f90wrap_n1
    !f2py intent(hide), depend(h8) :: f90wrap_n1 = shape(h8,0)
    type(fortrandata_ptr_type) :: data_ptr
    integer(c_int), intent(in), dimension(4) :: data
    character*(1024) :: filename
    integer(4) :: filenamelen
    real(8), dimension(f90wrap_n0) :: g8
    real(8), dimension(f90wrap_n1) :: h8
    data_ptr = transfer(data, data_ptr)
    call FieldTrace(Data=data_ptr%p, FileName=filename, FileNamelen=filenamelen, g8=g8, h8=h8)
end subroutine f90wrap_middleman__fieldtrace

subroutine f90wrap_middleman__mhdstartupsorted(f90wrap_n0, f90wrap_n1, f90wrap_n2, f90wrap_n3, f90wrap_n4, f90wrap_n5, &
    f90wrap_n6, f90wrap_n7, f90wrap_n8, f90wrap_n9, f90wrap_n10, f90wrap_n11, f90wrap_n12, f90wrap_n13, f90wrap_n14, &
    f90wrap_n15, xu, yu, zu, mhdposition_in, mhdb_in, nx_split, ny_split, nz_split, mix, max_bn, miy, may, miz, maz, &
    region_order_in, start_x, end_x, start_y, end_y, start_z, end_z, num_regions, xulen, yulen, zulen)
    use middleman
    implicit none
    
    integer :: f90wrap_n0
    !f2py intent(hide), depend(xu) :: f90wrap_n0 = shape(xu,0)
    integer :: f90wrap_n1
    !f2py intent(hide), depend(yu) :: f90wrap_n1 = shape(yu,0)
    integer :: f90wrap_n2
    !f2py intent(hide), depend(zu) :: f90wrap_n2 = shape(zu,0)
    integer :: f90wrap_n3
    !f2py intent(hide), depend(mhdposition_in) :: f90wrap_n3 = shape(mhdposition_in,0)
    integer :: f90wrap_n4
    !f2py intent(hide), depend(mhdposition_in) :: f90wrap_n4 = shape(mhdposition_in,1)
    integer :: f90wrap_n5
    !f2py intent(hide), depend(mhdposition_in) :: f90wrap_n5 = shape(mhdposition_in,2)
    integer :: f90wrap_n6
    !f2py intent(hide), depend(mhdb_in) :: f90wrap_n6 = shape(mhdb_in,0)
    integer :: f90wrap_n7
    !f2py intent(hide), depend(mhdb_in) :: f90wrap_n7 = shape(mhdb_in,1)
    integer :: f90wrap_n8
    !f2py intent(hide), depend(mhdb_in) :: f90wrap_n8 = shape(mhdb_in,2)
    integer :: f90wrap_n9
    !f2py intent(hide), depend(region_order_in) :: f90wrap_n9 = shape(region_order_in,0)
    integer :: f90wrap_n10
    !f2py intent(hide), depend(start_x) :: f90wrap_n10 = shape(start_x,0)
    integer :: f90wrap_n11
    !f2py intent(hide), depend(end_x) :: f90wrap_n11 = shape(end_x,0)
    integer :: f90wrap_n12
    !f2py intent(hide), depend(start_y) :: f90wrap_n12 = shape(start_y,0)
    integer :: f90wrap_n13
    !f2py intent(hide), depend(end_y) :: f90wrap_n13 = shape(end_y,0)
    integer :: f90wrap_n14
    !f2py intent(hide), depend(start_z) :: f90wrap_n14 = shape(start_z,0)
    integer :: f90wrap_n15
    !f2py intent(hide), depend(end_z) :: f90wrap_n15 = shape(end_z,0)
    real(8), dimension(f90wrap_n0) :: xu
    real(8), dimension(f90wrap_n1) :: yu
    real(8), dimension(f90wrap_n2) :: zu
    real(8), dimension(f90wrap_n3,f90wrap_n4,f90wrap_n5, 3) :: mhdposition_in
    real(8), dimension(f90wrap_n6,f90wrap_n7,f90wrap_n8, 3) :: mhdb_in
    integer :: nx_split
    integer :: ny_split
    integer :: nz_split
    real :: mix
    real :: max_bn
    real :: miy
    real :: may
    real :: miz
    real :: maz
    integer(4), dimension(f90wrap_n9) :: region_order_in
    integer(4), dimension(f90wrap_n10) :: start_x
    integer(4), dimension(f90wrap_n11) :: end_x
    integer(4), dimension(f90wrap_n12) :: start_y
    integer(4), dimension(f90wrap_n13) :: end_y
    integer(4), dimension(f90wrap_n14) :: start_z
    integer(4), dimension(f90wrap_n15) :: end_z
    integer(4) :: num_regions
    integer(4) :: xulen
    integer(4) :: yulen
    integer(4) :: zulen
    call MHDstartupSorted(XU=xu, YU=yu, ZU=zu, MHDposition_in=mhdposition_in, MHDB_in=mhdb_in, nx_split=nx_split, &
        ny_split=ny_split, nz_split=nz_split, mix=mix, max=max_bn, miy=miy, may=may, miz=miz, maz=maz, &
        region_order_in=region_order_in, start_x=start_x, end_x=end_x, start_y=start_y, end_y=end_y, start_z=start_z, &
        end_z=end_z, num_regions=num_regions, XUlen=xulen, YUlen=yulen, ZUlen=zulen)
end subroutine f90wrap_middleman__mhdstartupsorted

subroutine f90wrap_middleman__gse2gswtsy15(f90wrap_n0, f90wrap_n1, date, position_gse, wind, gotso, hotso, glen, &
    position_gsw)
    use middleman
    implicit none
    
    integer :: f90wrap_n0
    !f2py intent(hide), depend(gotso) :: f90wrap_n0 = shape(gotso,0)
    integer :: f90wrap_n1
    !f2py intent(hide), depend(hotso) :: f90wrap_n1 = shape(hotso,0)
    real(8), dimension(6) :: date
    real(8), dimension(3) :: position_gse
    real(8), dimension(3) :: wind
    real(8), dimension(f90wrap_n0) :: gotso
    real(8), dimension(f90wrap_n1) :: hotso
    integer(4) :: glen
    real(8), intent(inout), dimension(3) :: position_gsw
    call gse2gswTSY15(date=date, position_gse=position_gse, Wind=wind, gOTSO=gotso, hOTSO=hotso, glen=glen, &
        position_gsw=position_gsw)
end subroutine f90wrap_middleman__gse2gswtsy15

! End of module middleman defined in file MiddleMan.f95

