module ihamocc_shared
   use fabm_types, only: rk
   real(rk), parameter :: Xconvxa = 6.97e-07_rk !oxygen.f90 !carbon.f90 !cfc.f90        ! Wanninkhof's a=0.251 converted from [cm hr-1]/[m s-1]^2 to [ms-1]/[m s-1]^2       NIC: from carchm.f90
   real(rk), parameter :: atm2pa = 101325.0_rk !oxygen.f90 !carbon.f90                      ! conversion factor from atmospheres to pascal
   real(rk), parameter :: tzero = 273.15_rk !oxygen.f90 !carbon.f90 !bromo.f90 !cfc.f90                     ! absolute min temperature (*C) 
   
   ! mo_control_bgc parameters
   !---------------------------------------------------------------------------------
   real(rk), parameter :: dtbgc = 86400.0_rk              !  time step length [sec].
   real(rk), parameter :: dtb = 86400.0_rk              !  time step length [sec].

   
   ! mo_chemcon parameters
   !---------------------------------------------------------------------------------
   real(rk), parameter :: BOR1=0.000232_rk !carbon.f90                          !BORON CONCENTRATION IN SEA WATER IN G/KG PER O/OO CL (RILEY AND SKIRROW, 1965, P.250)
   real(rk), parameter :: BOR2=1./10.811_rk !carbon.f90                         !INVERSE OF ATOMIC WEIGHT OF BORON [G**-1] (USED TO CONVERT SPECIFIC TOTAL BORAT INTO CONCENTRATIONS)
   real(rk), parameter :: SALCHL=1./1.80655_rk !carbon.f90                      !CONVERSION FACTOR SALINITY -> CHLORINITY (AFTER WOOSTER ET AL., 1969)
   real(rk), parameter :: rrrcl=salchl*1.025_rk*bor1*bor2
   real(rk), parameter :: tzero=273.15_rk !oxygen.f90 !carbon.f90                           !ZERO DEG CENTIGRADE AT KELVIN SCALE
   real(rk), parameter :: CALCON=0.01028_rk !carbon.f90                                    !SET MEAN TOTAL [CA++] IN SEAWATER (MOLES/KG) (SEE BROECKER A. PENG, 1982, P. 26; [CA++](MOLES/KG)=1.028E-2*(S/35.); Value taken from Sarmiento and Gruber, 2006, p. 365
   real(rk), parameter :: OXYCO=1./22414.4_rk !oxygen.f90                       !INVERS OF NORMAL MOLAL VOLUME OF AN IDEAL GAS [mol/ml] at 0C
   real(rk), parameter :: OX0=-173.4292_rk !oxygen.f90                          !VOLUMETRIC SOLUBILITY CONSTANTS FOR O2 IN ML/L from moist air at one atm total pressure. Table 2 in WEISS, R. F. (1970) THE SOLUBILITY OF NITROGEN OXYGEN AND ARGON IN WATER AND SEAWATER. DEEP-SEA RESEARCH, VOL. 17, 721-735.
   real(rk), parameter :: OX1=249.6339_rk  !oxygen.f90 
   real(rk), parameter :: OX2=143.3483_rk  !oxygen.f90 
   real(rk), parameter :: OX3=-21.8492_rk  !oxygen.f90 
   real(rk), parameter :: OX4=-0.033096_rk !oxygen.f90 
   real(rk), parameter :: OX5=0.014259_rk  !oxygen.f90 
   real(rk), parameter :: OX6=-0.0017_rk   !oxygen.f90 
   real(rk), parameter :: AN0=-172.4965_rk                                      !VOLUMETRIC SOLUBILITY CONSTANTS FOR N2 IN ML/L from moist air at one atm total pressure. Table 2 in WEISS, R. F. (1970) THE SOLUBILITY OF NITROGEN OXYGEN AND ARGON IN WATER AND SEAWATER. DEEP-SEA RESEARCH, VOL. 17, 721-735.
   real(rk), parameter :: AN1=248.4262_rk
   real(rk), parameter :: AN2=143.0738_rk
   real(rk), parameter :: AN3=-21.7120_rk
   real(rk), parameter :: AN4=-0.049781_rk
   real(rk), parameter :: AN5=0.025018_rk
   real(rk), parameter :: AN6=-0.0034861_rk
   real(rk), parameter :: ac1= -162.8301_rk !carbon.f90                         !Constants for CO2 solubility in mol/kg/atm from moist air at one atm total pressure. Table 6 in WEISS, R.F., NITROUS OXIDE SOLUBILITY IN WATER AND SEAWATER, Marine Chemistry, 8, 347-359, 1980
   real(rk), parameter :: ac2= 218.2968_rk !carbon.f90
   real(rk), parameter :: ac3= 90.9241_rk !carbon.f90
   real(rk), parameter :: ac4= -1.47696_rk !carbon.f90
   real(rk), parameter :: bc1= 0.025695_rk !carbon.f90
   real(rk), parameter :: bc2= -0.025225_rk !carbon.f90
   real(rk), parameter :: bc3= 0.0049867_rk !carbon.f90
   real(rk), parameter :: ad1= -60.2409_rk !carbon.f90                          !Constants for CO2 solubility in mol/kg/atm for dry air at one atm total pressure. Table 1 in WEISS, R.F., CARBON DIOXIDE IN WATER AND SEAWATER: THE SOLUBILITY OF A NON - IDEAL GAS, Marine Chemistry, 2, 203-215, 1974
   real(rk), parameter :: ad2= 93.4517_rk !carbon.f90
   real(rk), parameter :: ad3= 23.3585_rk !carbon.f90
   real(rk), parameter :: bd1= 0.023517_rk !carbon.f90
   real(rk), parameter :: bd2= -0.023656_rk !carbon.f90
   real(rk), parameter :: bd3= 0.0047036_rk !carbon.f90
   real(rk), parameter :: al1= -165.8806_rk                                     !Constants for laughing gas solubility in mol/l/atm from moist air at one atm total pressure. Table 2 in WEISS, R.F., NITROUS OXIDE SOLUBILITY IN WATER AND SEAWATER, Marine Chemistry, 8, 347-359, 1980
   real(rk), parameter :: al2= 222.8743_rk
   real(rk), parameter :: al3= 92.0792_rk
   real(rk), parameter :: al4= -1.48425_rk
   real(rk), parameter :: bl1= -0.056235_rk
   real(rk), parameter :: bl2= 0.031619_rk
   real(rk), parameter :: bl3= -0.0048472_rk
   real(rk), parameter :: atn2o=3.e-7_rk                                        !Atmospheric mixing ratio of N2O around 1980 300 ppb
   REAL(rk), DIMENSION(11) :: a0, a1, a2, b0, b1, b2 !carbon.f90                !Constants needed for pressure correction of equilibrium constants. F. Millero, Thermodynamics of the carbon dioxide system in the oceans, Geochimica et Cosmochimica Acta, Vol. 59, No. 4, pp. 661-677, 1995
   DATA a0 /-25.5_rk, -15.82_rk, -29.48_rk, -25.60_rk, -18.03_rk, -9.78_rk, -48.76_rk, &
            -46._rk, -14.51_rk, -23.12_rk, -26.57_rk/
   DATA a1 /0.1271_rk, -0.0219_rk, 0.1622_rk, 0.2324_rk, 0.0466_rk, -0.0090_rk,     &
            0.5304_rk, 0.5304_rk, 0.1211_rk, 0.1758_rk, 0.2020_rk/
   DATA a2 /0.0_rk, 0.0_rk, 2.608e-3_rk, -3.6246e-3_rk, 0.316e-3_rk,             &
           -0.942e-3_rk, 0.0_rk, 0.0_rk, -0.321e-3_rk, -2.647e-3_rk, -3.042e-3_rk/
   DATA b0 /-3.08e-3_rk, 1.13e-3_rk, -2.84e-3_rk, -5.13e-3_rk, -4.53e-3_rk,      &
            -3.91e-3_rk, -11.76e-3_rk, -11.76e-3_rk, -2.67e-3_rk, -5.15e-3_rk,   & 
            -4.08e-3_rk/
   DATA b1 /0.0877e-3_rk, -0.1475e-3_rk, 0.0_rk, 0.0794e-3_rk, 0.09e-3_rk,       &
            0.054e-3_rk, 0.3692e-3_rk, 0.3692e-3_rk, 0.0427e-3_rk,            &
            0.09e-3_rk, 0.0714e-3_rk/
   DATA b2 /0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk/
   real(rk), parameter :: rgas = 83.131_rk !carbon.f90                           !Gas constant, value as used by Millero (1995)

   
end module