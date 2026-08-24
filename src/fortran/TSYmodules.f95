module TSY87module

   implicit none

   real(8) :: TILT

   integer(4) :: IP = 100

   real(8) :: FC = 0.3183099031D0
   real(8) :: RT = 30.0D0

   real(8) :: PA(32)

   real(8) :: C1
   real(8) :: RRC2
   real(8) :: DSTR
   real(8) :: XN
   real(8) :: RH
   real(8) :: X1 = 4.0D0
   real(8) :: X2 = 5.0D0
   real(8) :: DY
   real(8) :: B0
   real(8) :: B1
   real(8) :: XN21
   real(8) :: XN2
   real(8) :: XNR
   real(8) :: XN22
   real(8) :: ADLN

!$omp threadprivate( &
!$omp& TILT,IP,FC,RT,PA,C1,RRC2,DSTR,XN,RH,X1,X2,DY,B0,B1, &
!$omp& XN21,XN2,XNR,XN22,ADLN)

end module TSY87module

module TSY89module
    implicit none

    integer(4) :: IOP = 1000

    real(8) :: A(30)

    real(8) :: AK1, AK2, AK3, AK4, AK5
    real(8) :: AK6, AK7, AK8, AK9
    real(8) :: AK10, AK11, AK12, AK13
    real(8) :: AK14, AK15, AK16, AK17

    real(8) :: DYC, DYC2
    real(8) :: DX
    real(8) :: HA02
    real(8) :: RDX2M, RDX2
    real(8) :: RDYC2
    real(8) :: HLWC2M
    real(8) :: DRDYC2, DRDYC3
    real(8) :: HXLW2M

    real(8) :: ADR, D0, DD, RC, G, AT, DT
    real(8) :: DEL, P, Q, SX, GAM
    real(8) :: HXLD2M

    real(8) :: W1, W2, W3, W4, W5, W6
    real(8) :: DBLDEL

    real(8) :: AK610, AK711, AK812, AK913

    real(8) :: RDXL, HRDXL

    real(8) :: A6H, A9T, YNP, YND

    real(8) :: ADSL, XGHS, H, HS, GAMH

    real(8) :: SXA, SYA, SZA

!$omp threadprivate( &
!$omp& IOP, A, &
!$omp& AK1, AK2, AK3, AK4, AK5, AK6, AK7, AK8, AK9, &
!$omp& AK10, AK11, AK12, AK13, AK14, AK15, AK16, AK17, &
!$omp& DYC, DYC2, DX, HA02, RDX2M, RDX2, RDYC2, HLWC2M, &
!$omp& DRDYC2, DRDYC3, HXLW2M, ADR, D0, DD, RC, G, AT, DT, &
!$omp& DEL, P, Q, SX, GAM, HXLD2M, W1, W2, W3, W4, W5, W6, &
!$omp& DBLDEL, AK610, AK711, AK812, AK913, RDXL, HRDXL, &
!$omp& A6H, A9T, YNP, YND, ADSL, XGHS, H, HS, GAMH, &
!$omp& SXA, SYA, SZA)

end module TSY89module


module TSY96_Warp
   implicit none

   real(8) :: CPSS
   real(8) :: SPSS
   real(8) :: DPSRR
   real(8) :: RPS
   real(8) :: WARP
   real(8) :: D
   real(8) :: XS
   real(8) :: ZS
   real(8) :: DXSX
   real(8) :: DXSY
   real(8) :: DXSZ
   real(8) :: DZSX
   real(8) :: DZSY
   real(8) :: DZSZ
   real(8) :: DZETAS
   real(8) :: DDZETADX
   real(8) :: DDZETADY
   real(8) :: DDZETADZ
   real(8) :: ZSWW

!$omp threadprivate( &
!$omp& CPSS,SPSS,DPSRR,RPS,WARP,D,XS,ZS,DXSX,DXSY,DXSZ, &
!$omp& DZSX,DZSY,DZSZ,DZETAS,DDZETADX,DDZETADY,DDZETADZ,ZSWW)

end module TSY96_Warp

module TSY96_Coord11

   implicit none

   real(8) :: XX1(12) = (/ &
      -11.D0, -7.D0, -7.D0, &
       -3.D0, -3.D0, &
        1.D0,  1.D0,  1.D0, &
        5.D0,  5.D0, &
        9.D0,  9.D0 /)

   real(8) :: YY1(12) = (/ &
       2.D0, 0.D0, 4.D0, 2.D0, 6.D0, 0.D0, &
       4.D0, 8.D0, 2.D0, 6.D0, 0.D0, 4.D0 /)

!$omp threadprivate(XX1,YY1)

end module TSY96_Coord11

module TSY96_RHDR
   implicit none

   real(8) :: RH = 9.0D0
   real(8) :: DR = 4.0D0

!$omp threadprivate(RH,DR)

end module TSY96_RHDR

module TSY96_LoopDip1

   implicit none

   real(8) :: TILT = 1.00891D0

   real(8) :: XCENTRE(2) = (/ &
      2.28397D0, -5.60831D0 /)

   real(8) :: RADIUS(2) = (/ &
      1.86106D0, 7.83281D0 /)

   real(8) :: DIPX = 1.12541D0

   real(8) :: DIPY = 0.945719D0

!$omp threadprivate(TILT,XCENTRE,RADIUS,DIPX,DIPY)

end module TSY96_LoopDip1

module TSY96_Coord21

   implicit none

   real(8) :: XX2(14) = (/ &
      -10.D0, -7.D0, -4.D0, -4.D0, &
       0.D0,  4.D0,  4.D0,  7.D0, 10.D0, &
       0.D0,  0.D0,  0.D0,  0.D0,  0.D0 /)

   real(8) :: YY2(14) = (/ &
       3.D0, 6.D0, 3.D0, 9.D0, 6.D0, &
       3.D0, 9.D0, 6.D0, 3.D0, &
       0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /)

   real(8) :: ZZ2(14) = (/ &
      20.D0, 20.D0, 4.D0, 20.D0, 4.D0, 4.D0, &
      20.D0, 20.D0, 20.D0, 2.D0, 3.D0, &
      4.5D0, 7.D0, 10.D0 /)

!$omp threadprivate(XX2,YY2,ZZ2)

end module TSY96_Coord21

module TSY96_DX1
   implicit none

   real(8) :: DX = -0.16D0
   real(8) :: SCALEIN = 0.08D0
   real(8) :: SCALEOUT = 0.4D0

!$omp threadprivate(DX,SCALEIN,SCALEOUT)

end module TSY96_DX1

module TSY96_intercon
   implicit none

   integer(4) :: M = 0

   real(8) :: RP(3), RR(3)
   real(8) :: P(3), R(3)

!$omp threadprivate(M,RP,RR,P,R)

end module TSY96_intercon

module TSY96_TKSI
   implicit none

   integer(4) :: M = 0
   real(8) :: TDZ3

!$omp threadprivate(M,TDZ3)

end module TSY96_TKSI

module TSY96_Dipole
   implicit none

   integer(4) :: M = 0
   real(8) :: PSI, SPS, CPSS

   !$omp threadprivate(M,PSI,SPS,CPSS)

end module TSY96_Dipole

module TSY96_BIRK
    implicit none
    
    real(8) :: PSI = 10.D0
    real(8) :: SPS, CPS
    
    !$omp threadprivate(PSI,SPS,CPS)
    
    end module TSY96_BIRK

module TSY01_Tail
   implicit none

   real(8) :: DXSHIFT1, DXSHIFT2, D, DELTADY

   !$omp threadprivate(DXSHIFT1,DXSHIFT2,D,DELTADY)

end module TSY01_Tail

module TSY01_BIRKPAR
   implicit none

   real(8) :: XKAPPA1, XKAPPA2

   !$omp threadprivate(XKAPPA1, XKAPPA2)

end module TSY01_BIRKPAR

module TSY01_RCPAR
   implicit none

   real(8) :: SC_SY, SC_AS, PHI

   !$omp threadprivate(SC_SY, SC_AS, PHI)

end module TSY01_RCPAR

module TSY01_G
   implicit none

   real(8) :: GA

   !$omp threadprivate(GA)

end module TSY01_G

module TSY01_RH0
   implicit none

   real(8) :: RH0 = 8.0D0

   !$omp threadprivate(RH0)

end module TSY01_RH0

module TSY01_BDIP
   implicit none

   real(8) :: BOLD

   !$omp threadprivate(BOLD)

end module TSY01_BDIP

module TSY01_DPHI_B_RHO0
   implicit none

   real(8) :: DPHI,Ba,RHO_0,XKAPPA

   !$omp threadprivate(DPHI,Ba,RHO_0,XKAPPA)

end module TSY01_DPHI_B_RHO0

module TSY01_MODENUM
   implicit none

   integer(4) :: MO

   !$omp threadprivate(MO)

end module TSY01_MODENUM

module TSY01_DTHETA
   implicit none

   real(8) :: DTHETA

   !$omp threadprivate(DTHETA)

end module TSY01_DTHETA

module TSY01_WHERE_IN_MAGNETOPAUSE2001
   implicit none

   integer(4) :: LOCATION

   !$omp threadprivate(LOCATION)

end module TSY01_WHERE_IN_MAGNETOPAUSE2001

module TSY01_DIP
      implicit none

      integer(4) :: M = 0
      real(8) :: PSI = 5.D0

      save

      !$omp threadprivate(M, PSI)

end module TSY01_DIP

module TSY01_INIT
      implicit none

      real(8) :: PDYN, DST_AST, BYIMF, BZIMF, G1, G2, G3, PSS, B0DIP

      save PDYN, DST_AST, BYIMF, BZIMF, G1, G2, G3, PSS, B0DIP

      !$omp threadprivate(PDYN, DST_AST, BYIMF, BZIMF, G1, G2, G3, PSS, B0DIP)
end module TSY01_INIT

module TSY04_Tail
   implicit none

   real(8) :: DXSHIFT1, DXSHIFT2, D, DELTADY

   !$omp threadprivate(DXSHIFT1,DXSHIFT2,D,DELTADY)

end module TSY04_Tail

module TSY04_BIRKPAR
   implicit none

   real(8) :: XKAPPA1, XKAPPA2

   !$omp threadprivate(XKAPPA1, XKAPPA2)

end module TSY04_BIRKPAR

module TSY04_RCPAR
   implicit none

   real(8) :: SC_SY, SC_AS, PHI

   !$omp threadprivate(SC_SY, SC_AS, PHI)

end module TSY04_RCPAR

module TSY04_G
   implicit none

   real(8) :: GA

   !$omp threadprivate(GA)

end module TSY04_G

module TSY04_RH0
   implicit none

   real(8) :: RH0 = 7.5D0

   !$omp threadprivate(RH0)

end module TSY04_RH0

module TSY04_INIT
      implicit none

      real(8) :: PDYN, DST_AST, BYIMF, BZIMF, PSS
      real(8) :: W1, W2, W3, W4, W5, W6

      save PDYN, DST_AST, BYIMF, BZIMF, PSS, W1, W2, W3, W4, W5, W6

      !$omp threadprivate(PDYN, DST_AST, BYIMF, BZIMF, PSS, W1, W2, W3, W4, W5, W6)
end module TSY04_INIT

module TSY04_DPHI_B_RHO0
   implicit none

   real(8) :: DPHI,Ba,RHO_0,XKAPPA

   !$omp threadprivate(DPHI,Ba,RHO_0,XKAPPA)

end module TSY04_DPHI_B_RHO0

module TSY04_MODENUM
   implicit none

   integer(4) :: MO

   !$omp threadprivate(MO)

end module TSY04_MODENUM

module TSY04_DTHETA
   implicit none

   real(8) :: DTHETA

   !$omp threadprivate(DTHETA)

end module TSY04_DTHETA

module TSY15N_XYZD
   implicit none

   real(8) :: CURDPHI,DD(15),SP(25),CP(25),XXN(15,25),YYN(15,25)
   real(8) :: ZZN(15,25),XXS(15,25),YYS(15,25),ZZS(15,25)
   real(8) :: T0 = 10.D0, DT0 = 10.D0, PS0 = 10.D0

   save

   !$omp threadprivate(CURDPHI,DD,SP,CP,XXN,YYN,ZZN,XXS,YYS,ZZS,T0,DT0,PS0)

end module TSY15N_XYZD

module TSY15N_DIP
   implicit none

   real(8) :: SPS, CPS
   real(8) :: PSI = 5.0D0
   integer(4)  :: M = 0

   save

   !$omp threadprivate(SPS,CPS,PSI,M)

end module TSY15N_DIP

module TSY15N_KEEP
   implicit none

   integer(4) :: KEEP = 0

   save

   !$omp threadprivate(KEEP)

end module TSY15N_KEEP

module TSY15N_PRC_COEFFS

   implicit none

   real(8) :: &
      ABR01 = 18.83023014D0, &
      ABR02 = 468934.6659D0, &
      ABR03 = -468945.9056D0, &
      ABR04 = 0.9310123492D-03, &
      ABR05 = 341.5640336D0, &
      ABR06 = -341.5724899D0, &
      ABR07 = 0.2834275159D-01, &
      ABR08 = 2815.457002D0, &
      ABR09 = -2815.394238D0, &
      ABR10 = 4273.405633D0, &
      ABR11 = 142912593.7D0, &
      ABR12 = -142912703.7D0, &
      ABR13 = -.3144179212D-08, &
      ABR14 = -.9652596615D-05, &
      ABR15 = 0.9661074043D-05, &
      ABR16 = 0.2089306536D-02, &
      ABR17 = -14080.28783D0, &
      ABR18 = 14080.43837D0, &
      ABR19 = 0.1113334263D0, &
      ABR20 = -0.1062192861D0, &
      ABR21 = 0.3149232975D0, &
      ABR22 = 11.36515559D0, &
      ABR23 = 4.723154529D0, &
      ABR24 = 3.984559924D0, &
      ABR25 = 4.585739484D0, &
      ABR26 = 4.032706815D0, &
      ABR27 = 3.741667937D0, &
      ABR28 = 0.5798779413D-01, &
      ABR29 = 2.417509610D0, &
      ABR30 = 8.492134117D0, &
      ABR31 = 7.761662895D0, &
      ABR32 = 7.389416256D0, &
      ABR33 = 7.460193567D0, &
      ABR34 = 1.110098139D0, &
      ABR35 = 1.110074621D0, &
      ABR36 = 2.821936969D0, &
      ABR37 = 9.037995133D0, &
      ABR38 = 5.282888651D0, &
      ABR39 = 6.262172568D0, &
      ABR40 = 0.8249294900D0, &
      ABR41 = 1.195457593D0, &
      ABR42 = 2.469064281D0, &
      ABR43 = 5.071438197D0, &
      ABR44 = 4.822507104D0, &
      ABR45 = 0.1568808225D0, &
      ABR46 = 0.1079224544D0, &
      ABR47 = 0.1616535726D0, &
      ABR48 = 2.078602607D0, &
      ABR49 = 0.3735552854D0, &
      ABR50 = 4.279143730D0, &
      ABR51 = 1.863086547D0, &
      ABR52 = 1.080625720D0, &
      ABR53 = 1.104468013D0, &
      ABR54 = 6.456298650D0, &
      ABR55 = 1.970461884D0, &
      ABR56 = 1.974588436D0

   real(8) :: &
      ABT01 = 3.356791597D0, &
      ABT02 = -5.883419687D0, &
      ABT03 = -2.058543379D0, &
      ABT04 = .1100303924D-03, &
      ABT05 = -.1943189958D-02, &
      ABT06 = .2793918244D-02, &
      ABT07 = 0.2049870425D-03, &
      ABT08 = 0.5220714000D-03, &
      ABT09 = 0.1019940452D-01, &
      ABT10 = -75.87370859D0, &
      ABT11 = -47.96110111D0, &
      ABT12 = 24.36085698D0, &
      ABT13 = -.5362635518D-05, &
      ABT14 = 0.1925145540D-04, &
      ABT15 = -0.3548180534D-05, &
      ABT16 = -.3037668516D-03, &
      ABT17 = .7799787312D-03, &
      ABT18 = -.1411180236D-01, &
      ABT19 = .6526098636D-01, &
      ABT20 = -0.6165246724D-01, &
      ABT21 = 0.7321218781D0, &
      ABT22 = 11.63639868D0, &
      ABT23 = 4.570394366D0, &
      ABT24 = 5.991114966D0, &
      ABT25 = 7.610650545D0, &
      ABT26 = 4.318030252D0, &
      ABT27 = 4.599681445D0, &
      ABT28 = 0.3488975773D0, &
      ABT29 = 1.700478186D0, &
      ABT30 = 5.759626972D0, &
      ABT31 = 1.666470246D0, &
      ABT32 = 7.067086748D0, &
      ABT33 = 7.451917780D0, &
      ABT34 = 0.4817635507D0, &
      ABT35 = 2.681223071D0, &
      ABT36 = 4.434829894D0, &
      ABT37 = 5.173486046D0, &
      ABT38 = 0.4858312077D0, &
      ABT39 = 5.763669434D0, &
      ABT40 = 3.765430495D0, &
      ABT41 = 1.565144414D0, &
      ABT42 = 2.389620597D0, &
      ABT43 = 6.040651474D0, &
      ABT44 = 5.200117858D0, &
      ABT45 = 0.1443594185D0, &
      ABT46 = 0.1078706697D0, &
      ABT47 = 0.1655199312D0, &
      ABT48 = 1.108410791D0, &
      ABT49 = 0.2563543278D0, &
      ABT50 = 1.647886458D0, &
      ABT51 = 2.002342423D0, &
      ABT52 = 1.537733423D0, &
      ABT53 = 1.010006109D0, &
      ABT54 = 2.360742095D0, &
      ABT55 = 1.907047647D0, &
      ABT56 = 3.557905953D0

   save

!$omp threadprivate( &
!$omp ABR01,ABR02,ABR03,ABR04,ABR05,ABR06,ABR07,ABR08,ABR09,ABR10, &
!$omp ABR11,ABR12,ABR13,ABR14,ABR15,ABR16,ABR17,ABR18,ABR19,ABR20, &
!$omp ABR21,ABR22,ABR23,ABR24,ABR25,ABR26,ABR27,ABR28,ABR29,ABR30, &
!$omp ABR31,ABR32,ABR33,ABR34,ABR35,ABR36,ABR37,ABR38,ABR39,ABR40, &
!$omp ABR41,ABR42,ABR43,ABR44,ABR45,ABR46,ABR47,ABR48,ABR49,ABR50, &
!$omp ABR51,ABR52,ABR53,ABR54,ABR55,ABR56, &
!$omp ABT01,ABT02,ABT03,ABT04,ABT05,ABT06,ABT07,ABT08,ABT09,ABT10, &
!$omp ABT11,ABT12,ABT13,ABT14,ABT15,ABT16,ABT17,ABT18,ABT19,ABT20, &
!$omp ABT21,ABT22,ABT23,ABT24,ABT25,ABT26,ABT27,ABT28,ABT29,ABT30, &
!$omp ABT31,ABT32,ABT33,ABT34,ABT35,ABT36,ABT37,ABT38,ABT39,ABT40, &
!$omp ABT41,ABT42,ABT43,ABT44,ABT45,ABT46,ABT47,ABT48,ABT49,ABT50, &
!$omp ABT51,ABT52,ABT53,ABT54,ABT55,ABT56)

end module TSY15N_PRC_COEFFS

module TSY15B_XYZD
   implicit none

   real(8) :: CURDPHI,DD(15),SP(25),CP(25),XXN(15,25),YYN(15,25)
   real(8) :: ZZN(15,25),XXS(15,25),YYS(15,25),ZZS(15,25)
   real(8) :: T0 = 10.D0, DT0 = 10.D0, PS0 = 10.D0

   save

   !$omp threadprivate(CURDPHI,DD,SP,CP,XXN,YYN,ZZN,XXS,YYS,ZZS,T0,DT0,PS0)

end module TSY15B_XYZD

module TSY15B_DIP
   implicit none

   real(8) :: SPS, CPS
   real(8) :: PSI = 5.0D0
   integer(4)  :: M = 0

   save

   !$omp threadprivate(SPS,CPS,PSI,M)

end module TSY15B_DIP

module TSY15B_KEEP
   implicit none

   integer(4) :: KEEP = 0

   save

   !$omp threadprivate(KEEP)

end module TSY15B_KEEP

module TSY15B_PRC_COEFFS

   implicit none

   real(8) :: &
      ABR01 = 18.83023014D0, &
      ABR02 = 468934.6659D0, &
      ABR03 = -468945.9056D0, &
      ABR04 = 0.9310123492D-03, &
      ABR05 = 341.5640336D0, &
      ABR06 = -341.5724899D0, &
      ABR07 = 0.2834275159D-01, &
      ABR08 = 2815.457002D0, &
      ABR09 = -2815.394238D0, &
      ABR10 = 4273.405633D0, &
      ABR11 = 142912593.7D0, &
      ABR12 = -142912703.7D0, &
      ABR13 = -.3144179212D-08, &
      ABR14 = -.9652596615D-05, &
      ABR15 = 0.9661074043D-05, &
      ABR16 = 0.2089306536D-02, &
      ABR17 = -14080.28783D0, &
      ABR18 = 14080.43837D0, &
      ABR19 = 0.1113334263D0, &
      ABR20 = -0.1062192861D0, &
      ABR21 = 0.3149232975D0, &
      ABR22 = 11.36515559D0, &
      ABR23 = 4.723154529D0, &
      ABR24 = 3.984559924D0, &
      ABR25 = 4.585739484D0, &
      ABR26 = 4.032706815D0, &
      ABR27 = 3.741667937D0, &
      ABR28 = 0.5798779413D-01, &
      ABR29 = 2.417509610D0, &
      ABR30 = 8.492134117D0, &
      ABR31 = 7.761662895D0, &
      ABR32 = 7.389416256D0, &
      ABR33 = 7.460193567D0, &
      ABR34 = 1.110098139D0, &
      ABR35 = 1.110074621D0, &
      ABR36 = 2.821936969D0, &
      ABR37 = 9.037995133D0, &
      ABR38 = 5.282888651D0, &
      ABR39 = 6.262172568D0, &
      ABR40 = 0.8249294900D0, &
      ABR41 = 1.195457593D0, &
      ABR42 = 2.469064281D0, &
      ABR43 = 5.071438197D0, &
      ABR44 = 4.822507104D0, &
      ABR45 = 0.1568808225D0, &
      ABR46 = 0.1079224544D0, &
      ABR47 = 0.1616535726D0, &
      ABR48 = 2.078602607D0, &
      ABR49 = 0.3735552854D0, &
      ABR50 = 4.279143730D0, &
      ABR51 = 1.863086547D0, &
      ABR52 = 1.080625720D0, &
      ABR53 = 1.104468013D0, &
      ABR54 = 6.456298650D0, &
      ABR55 = 1.970461884D0, &
      ABR56 = 1.974588436D0

   real(8) :: &
      ABT01 = 3.356791597D0, &
      ABT02 = -5.883419687D0, &
      ABT03 = -2.058543379D0, &
      ABT04 = .1100303924D-03, &
      ABT05 = -.1943189958D-02, &
      ABT06 = .2793918244D-02, &
      ABT07 = 0.2049870425D-03, &
      ABT08 = 0.5220714000D-03, &
      ABT09 = 0.1019940452D-01, &
      ABT10 = -75.87370859D0, &
      ABT11 = -47.96110111D0, &
      ABT12 = 24.36085698D0, &
      ABT13 = -.5362635518D-05, &
      ABT14 = 0.1925145540D-04, &
      ABT15 = -0.3548180534D-05, &
      ABT16 = -.3037668516D-03, &
      ABT17 = .7799787312D-03, &
      ABT18 = -.1411180236D-01, &
      ABT19 = .6526098636D-01, &
      ABT20 = -0.6165246724D-01, &
      ABT21 = 0.7321218781D0, &
      ABT22 = 11.63639868D0, &
      ABT23 = 4.570394366D0, &
      ABT24 = 5.991114966D0, &
      ABT25 = 7.610650545D0, &
      ABT26 = 4.318030252D0, &
      ABT27 = 4.599681445D0, &
      ABT28 = 0.3488975773D0, &
      ABT29 = 1.700478186D0, &
      ABT30 = 5.759626972D0, &
      ABT31 = 1.666470246D0, &
      ABT32 = 7.067086748D0, &
      ABT33 = 7.451917780D0, &
      ABT34 = 0.4817635507D0, &
      ABT35 = 2.681223071D0, &
      ABT36 = 4.434829894D0, &
      ABT37 = 5.173486046D0, &
      ABT38 = 0.4858312077D0, &
      ABT39 = 5.763669434D0, &
      ABT40 = 3.765430495D0, &
      ABT41 = 1.565144414D0, &
      ABT42 = 2.389620597D0, &
      ABT43 = 6.040651474D0, &
      ABT44 = 5.200117858D0, &
      ABT45 = 0.1443594185D0, &
      ABT46 = 0.1078706697D0, &
      ABT47 = 0.1655199312D0, &
      ABT48 = 1.108410791D0, &
      ABT49 = 0.2563543278D0, &
      ABT50 = 1.647886458D0, &
      ABT51 = 2.002342423D0, &
      ABT52 = 1.537733423D0, &
      ABT53 = 1.010006109D0, &
      ABT54 = 2.360742095D0, &
      ABT55 = 1.907047647D0, &
      ABT56 = 3.557905953D0

   save

!$omp threadprivate( &
!$omp ABR01,ABR02,ABR03,ABR04,ABR05,ABR06,ABR07,ABR08,ABR09,ABR10, &
!$omp ABR11,ABR12,ABR13,ABR14,ABR15,ABR16,ABR17,ABR18,ABR19,ABR20, &
!$omp ABR21,ABR22,ABR23,ABR24,ABR25,ABR26,ABR27,ABR28,ABR29,ABR30, &
!$omp ABR31,ABR32,ABR33,ABR34,ABR35,ABR36,ABR37,ABR38,ABR39,ABR40, &
!$omp ABR41,ABR42,ABR43,ABR44,ABR45,ABR46,ABR47,ABR48,ABR49,ABR50, &
!$omp ABR51,ABR52,ABR53,ABR54,ABR55,ABR56, &
!$omp ABT01,ABT02,ABT03,ABT04,ABT05,ABT06,ABT07,ABT08,ABT09,ABT10, &
!$omp ABT11,ABT12,ABT13,ABT14,ABT15,ABT16,ABT17,ABT18,ABT19,ABT20, &
!$omp ABT21,ABT22,ABT23,ABT24,ABT25,ABT26,ABT27,ABT28,ABT29,ABT30, &
!$omp ABT31,ABT32,ABT33,ABT34,ABT35,ABT36,ABT37,ABT38,ABT39,ABT40, &
!$omp ABT41,ABT42,ABT43,ABT44,ABT45,ABT46,ABT47,ABT48,ABT49,ABT50, &
!$omp ABT51,ABT52,ABT53,ABT54,ABT55,ABT56)

end module TSY15B_PRC_COEFFS

module TSY16module
      implicit none

      real(8) :: XX(1296),YY(1296),ZZ(1296),ST(1296)
      real(8) :: RHO(1296),ZSP(1296),ZCP(1296),RHBR(1296)
      real(8) :: D = 4.D0
      real(8) :: PI = 3.14159265359D0

      integer(4) :: IOP = 100001

      real(8) :: A0  = 12.544d0,  &
         A1  = -0.194d0,  &
         A2  = 0.305d0,   &
         A3  = 0.0573d0,  &
         A4  = 2.178d0,   &
         A5  = 0.0571d0,  &
         A6  = -0.999d0,  &
         A7  = 16.473d0,  &
         A8  = 0.00152d0, &
         A9  = 0.382d0,   &
         A10 = 0.0431d0,  &
         A11 = -0.00763d0,&
         A12 = -0.210d0,  &
         A13 = 0.0405d0,  &
         A14 = -4.430d0,  &
         A15 = -0.636d0,  &
         A16 = -2.600d0,  &
         A17 = 0.832d0,   &
         A18 = -5.328d0,  &
         A19 = 1.103d0,   &
         A20 = -0.907d0,  &
         A21 = 1.450d0

!$omp threadprivate( &
!$omp& XX, YY, ZZ, ST, RHO, ZSP, ZCP, RHBR, &
!$omp& D, PI, IOP, &
!$omp& A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, &
!$omp& A10, A11, A12, A13, A14, A15, A16, A17, A18, A19, A20, A21 )

end module TSY16module
      