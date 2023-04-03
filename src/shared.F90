module ihamocc_shared
   use fabm_types, only: rk
   real(rk), parameter :: tzero = 273.15_rk !ZERO DEG CENTIGRADE AT KELVIN SCALE
   real(rk), parameter :: OXYCO=1._rk/22414.4_rk  !INVERS OF NORMAL MOLAL VOLUME OF AN IDEAL GAS [mol/ml] at 0C
   real(rk), parameter :: Xconvxa = 6.97e-07_rk   ! Wanninkhof's a=0.251 converted from [cm hr-1]/[m s-1]^2 to [ms-1]/[m s-1]^2 
   real(rk), parameter :: dtbgc = 86400.0_rk              !  time step length [sec].
   real(rk), parameter :: atm2pa = 101325.0_rk ! conversion factor from atmospheres to pascal
end module