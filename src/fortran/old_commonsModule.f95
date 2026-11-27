module Soleil
  implicit none
  real*8 :: xS,yS,zS,cgsta,sgsta
!$omp threadprivate(xS,yS,zS,cgsta,sgsta)
end module Soleil

module sundip
  implicit none
  real*8 :: xD,yD,zD
  real*8 :: xSD,ySD,zSD
  real*8 :: xSSD,ySSD,zSSD
  real*8 :: xSDD,ySDD,zSDD
!$omp threadprivate(xD,yD,zD,xSD,ySD,zSD, &
!$omp              xSSD,ySSD,zSSD,xSDD,ySDD,zSDD)
end module sundip

module dipigrf
  implicit none

  real*8 :: Bo
  real*8 :: xc, yc, zc
  real*8 :: ct, st, cp, sp

!$omp threadprivate(Bo, xc, yc, zc, ct, st, cp, sp)

end module dipigrf

module rconst
  implicit none

  real*8 :: rad, pi

!$omp threadprivate(rad, pi)

end module rconst

module MODEL
  implicit none

  real*8 :: GH1(144)

!$omp threadprivate(GH1)

end module MODEL

module dgrf
  implicit none

  real*8 :: Ga(66), Ha(66)

!$omp threadprivate(Ga, Ha)

end module dgrf

module GENER
  implicit none

  real*8 :: ERA, AQUAD, BQUAD

!$omp threadprivate(ERA, AQUAD, BQUAD)

end module GENER

module sunMAH
  implicit none

  real*8 :: T0, Lamda0
  real*8 :: st, ct, sl, cl, so, co

!$omp threadprivate(T0, Lamda0, st, ct, sl, cl, so, co)

end module sunMAH

module irbem_cache
  implicit none

  logical :: dtd_valid = .false.
  integer(8) :: dtd_ghver = -1
  integer(8) :: dtd_gen = 0

  logical :: gsm_valid = .false.
  integer(8) :: gsm_gen = -1
  integer(4) :: gsm_iyr = 0, gsm_iday = 0
  integer(8) :: gsm_secs = 0
  real(8) :: gsm_psi = 0.0d0

!$omp threadprivate(dtd_valid, dtd_ghver, dtd_gen, gsm_valid, gsm_gen, &
!$omp              gsm_iyr, gsm_iday, gsm_secs, gsm_psi)

end module irbem_cache
