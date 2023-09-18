#include "fabm_driver.h"

module ihamocc_cisonew

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_cisonew
      type (type_dependency_id) :: id_psao, id_ptho, id_prho, id_prb, id_silica, id_hi_in, id_pddpo, id_phytomi, id_phosy, id_phymor, id_phyrem, id_pommor, id_dimmor, id_phosy
      type (type_dependency_id) :: id_exud, id_graton, id_grawa, id_gratpoc, id_domex, id_pocrem, id_docrem, id_delcar_part, id_remin2o
      type (type_surface_dependency_id) :: id_atmc13, id_atmc14, id_kwco2sol
      type (type_state_variable_id) :: id_sco212, id_calc13, id_calc13, id_calc, id_sco213, id_sco214, id_phy, id_doc
      type (type_diagnostic_variable_id) :: id_co213fxd, id_co213fxu, id_co214fxd, id_co214fxu
      
      real(rk) :: 
      
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_surface
      procedure :: do
      procedure :: get_vertical_movement
   end type type_ihamocc_cisonew

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_cisonew), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      ! Register parameters
      call self%get_parameter(self%c14_t_half, 'c14_t_half', 'd','c14_t_half', default=5730._rk*365._rk) ! 

      
      ! Register environmental dependencies
      call self%register_dependency(self%id_atmc13, 'atmc13', '-', 'surface air carbon dioxide 13 mixing ratio') ! atmospheric co2 13 mixing ratio (i.e. partial presure = mixing ratio*SLP/P_0 [atm]) 
      call self%register_dependency(self%id_atmc14, 'atmc14', '-', 'surface air carbon dioxide 14 mixing ratio') ! atmospheric co2 14 mixing ratio (i.e. partial presure = mixing ratio*SLP/P_0 [atm]) 
      call self%register_dependency(self%id_sco212, 'sco212', 'kmol/m^3', 'Dissolved co2')
      call self%register_dependency(self%id_hi, 'hi', 'mol/kg', 'Hydrogen ion concentration')
      call self%register_dependency(self%id_kwco2sol, 'kwco2sol', 'm/s mol/kg/microatm', 'kwco2sol')
      call self%register_dependency(self%id_psao, standard_variables%practical_salinity)
      call self%register_dependency(self%id_ptho, standard_variables%temperature)
      call self%register_dependency(self%id_prho, standard_variables%density)
      call self%register_dependency(self%id_prb, standard_variables%pressure)
      call self%register_dependency(self%id_pddpo, standard_variables%cell_thickness)
      call self%register_dependency(self%id_pfu10, standard_variables%wind_speed)
      call self%register_dependency(self%id_psicomo, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_ppao, standard_variables%surface_air_pressure) ! surface air pressure in pascal
      call self%register_dependency(self%id_phytomi,'phytomi', 'kmol P/m3', 'minimum concentration of phytoplankton')
      call self%register_dependency(self%id_phosy, 'phosy', 'kmol/m3/d', 'photosynthetic rate')
      
      call self%register_dependency(self%id_phymor, 'phymor', 'kmol/m3/d', 'photosynthetic mortality rate')
      call self%register_dependency(self%id_phyrem, 'phyrem', 'kmol/m3/d', 'photosynthetic remineralization rate')
      call self%register_dependency(self%id_pommor, 'pommor', 'kmol/m3/d', 'zooplankton particulate export from mortality')
      call self%register_dependency(self%id_dimmor, 'dimmor', 'kmol/m3/d', 'zooplankton dissolved inorganic export from mortality')
      call self%register_dependency(self%id_phosy, 'phosy', 'kmol/m3/d', 'photosynthetic rate')
      call self%register_dependency(self%id_exud, 'exud', 'kmol/m3/d', 'phytoplankton exudation rate')
      call self%register_dependency(self%id_graton, 'graton', 'kmol/m3/d', 'zooplankton sloppy feeding inorganic release rate')
      call self%register_dependency(self%id_grawa, 'grawa', 'kmol/m3/d', 'zooplankton assimilation rate')
      call self%register_dependency(self%id_gratpoc, 'gratpoc', 'kmol/m3/d', 'zooplankton sloppy feeding particulate release rate')
      call self%register_dependency(self%id_domex, 'domex', 'kmol/m3/d', 'zooplankton dissolved organic export')
      call self%register_dependency(self%id_pocrem, 'pocrem', 'kmol/m^3 d-1', 'deep remineralization of POC') 
      call self%register_dependency(self%id_docrem, 'docrem', 'kmol/m^3 d-1', 'deep remineralization of DOC')
      call self%register_dependency(self%id_remin2o, 'remin2o', 'kmol/m^3 d-1', 'remin2o')
      call self%register_dependency(self%id_delcar_part, 'delcar_part', '-', 'calcium fraction of export') 

      
      ! Register state variables
      call self%register_state_variable(self%id_sco213, 'sco213', 'kmol/m^3', 'Dissolved co2 13')
      call self%register_state_variable(self%id_sco214, 'sco214', 'kmol/m^3', 'Dissolved co2 14')
      call self%register_state_variable(self%id_calc13, 'calc13', 'kmol/m^3', 'Calcium carbonate 13')
      call self%register_state_variable(self%id_calc14, 'calc14', 'kmol/m^3', 'Calcium carbonate 14')
      call self%register_state_variable(self%id_det14, 'det14', '', 'detritus 14')
      call self%register_state_variable(self%id_det13, 'det13', '', 'detritus 13')
      call self%register_state_variable(self%id_doc14, 'doc14', '', 'dissolved organic carbon 14')
      call self%register_state_variable(self%id_doc13, 'doc13', '', 'dissolved organic carbon 13')
      call self%register_state_variable(self%id_phy14, 'phy14', '', 'phytoplankton 14')
      call self%register_state_variable(self%id_phy13, 'phy13', '', 'phytoplankton 13')
      call self%register_state_variable(self%id_zoo14, 'zoo14', '', 'zooplankton 14')
      call self%register_state_variable(self%id_zoo13, 'zoo13', '', 'zooplankton 13')

      ! Register environmental dependencies
      call self%register_state_dependency(self%id_calc, 'calc',   'kmol/m^3', 'Calcium carbonate)
      call self%register_state_dependency(self%id_phy, 'phy', 'kmol/m^3', 'phytoplankton')
      call self%register_state_dependency(self%id_doc, 'doc','kmol/m^3', 'dissolvecd organic carbon')
      call self%register_state_dependency(self%id_det, 'det','kmol/m^3', 'detritus')



      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_co213fxd, 'co213fxd', 'kmol/m2/s', 'Downwards co2 13 surface flux')
      call self%register_diagnostic_variable(self%id_co213fxu, 'co213fxu', 'kmol/m2/s', 'Downwards co2 13 surface flux')
      call self%register_diagnostic_variable(self%id_co214fxd, 'co214fxd', 'kmol/m2/s', 'Downwards co2 14 surface flux')
      call self%register_diagnostic_variable(self%id_co214fxu, 'co214fxu', 'kmol/m2/s', 'Downwards co2 14 surface flux')


   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_cisonew), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: t, t2, t3, t4, tk, tk100, s, psao, ptho, prho, prb, sco212, alkali, silica, phosph, hi, cu, ac, K1, K2, pco2, scco2, ppao, pfu10, kwco2, rpp0, fluxu, fluxd, ta
      
      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_kwco2sol, kwco2sol)
         _GET_(self%id_id_ptho, ptho)
         _GET_(self%id_psao, psao)
         !_GET_(self%id_pddpo, pddpo)
         _GET_(self%id_prho, prho)
         _GET_(self%id_prb, prb)
         _GET_(self%id_sco212, sco212)
         _GET_(self%id_sco213, sco213)
         _GET_(self%id_sco214, sco214)
         !_GET_(self%id_alkali, alkali)
         _GET_(self%id_hi, hi)
         !_GET_(self%id_silica, silica)
         !_GET_(self%id_phosph, phosph)
         _GET_SURFACE_(self%id_atmc14, atco214)
         _GET_SURFACE_(self%id_atmc13, atco213)
         _GET_SURFACE_(self%id_ppao, ppao)
         _GET_SURFACE_(self%id_pfu10, pfu10)
         _GET_SURFACE_(self%id_psicomo, psicomo)

         ! Carbon chemistry: Calculate equilibrium constants and solve for [H+] and
         ! carbonate alkalinity (ac)
         t    = min(40._rk,max(-3._rk,ptho))
         t2   = t**2
         t3   = t**3
         t4   = t**4
         tk   = t + tzero
         tk100= tk/100.0_rk
         s    = min(40._rk,max( 25._rk,psao))
         rrho = prho/1000.0_rk                ! seawater density [kg/m3]->[g/cm3]
         prb  = prb*10._rk  !convert from dbar to bar. ORIGINAL: ptiestu(i,j,k)*98060*1.027e-6_rk ! pressure in unit bars, 98060 = onem
   
         tc   = sco212 / rrho  ! convert to mol/kg
         ah1  = hi
   
         CALL CARCHM_KEQUI(t,s,prb,Kh,Khd,K1,K2,Kb,Kw,Ks1,Kf,Ksi,             &
                           K1p,K2p,K3p,Kspc,Kspa)
   
         CALL CARCHM_SOLVE(s,tc,ta,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p, &
                           ah1,ac,niter)
   
         ! Determine CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
         cu = ( 2._rk * tc - ac ) / ( 2._rk + K1 / ah1 )
         cb = K1 * cu / ah1
         cc = K2 * cb / ah1
  
         pco2 = cu * 1.e6_rk / Kh ! Determine CO2 pressure and fugacity (in micoatm)   NOTE: equation below for pCO2 needs requires CO2 in mol/kg

         ! Ratio P/P_0, where P is the local SLP and P_0 is standard pressure (1 atm). This is
         ! used in all surface flux calculations where atmospheric concentration is given as a
         ! mixing ratio (i.e. partial presure = mixing ratio*SLP/P_0 [atm])
         rpp0 = ppao/atm2pa
         
         ! Ocean-Atmosphere fluxes for carbon isotopes
         rco213=sco213/(sco212+safediv) ! Fraction DIC13 over total DIC
         rco214=sco214/(sco212+safediv) ! Fraction DIC14 over total DIC

         pco213 = pco2 * rco213 ! Determine water CO213 pressure and fugacity (microatm)
         pco214 = pco2 * rco214 ! Determine water CO214 pressure and fugacity (microatm)
         
         ! fractionation factors for 13C during air-sea gas exchange (Zhang et al. 1995, Orr et al. 2017)
         frac_k    = 0.99912_rk                                 !Constant kinetic fractionation
         frac_aqg  = (0.0049_rk*t - 1.31_rk)/1000._rk + 1._rk  !Gas dissolution fractionation
         frac_dicg = (0.0144_rk*t*(cc/(cc+cu+cb)) - 0.107_rk*t + 10.53_rk)/1000._rk + 1._rk !DIC to CO2 frac
         flux13d=atco213*rpp0*dtbgc*kwco2sol*rrho*frac_aqg*frac_k         
         flux13u=pco213      *dtbgc*kwco2sol*rrho*frac_aqg*frac_k/frac_dicg   
         flux14d=atco214*rpp0*dtbgc*kwco2sol*rrho*(frac_aqg**2._rk)*(frac_k**2._rk)           
         flux14u=pco214      *dtbgc*kwco2sol*rrho*(frac_aqg**2._rk)*(frac_k**2._rk)/(frac_dicg**2._rk)  
       
         _ADD_SURFACE_FLUX_(self%id_sco213, (flux13d-flux13u)/dtbgc)  !Nic: divided by the time step to get instantaneous rate of change
         _ADD_SURFACE_FLUX_(self%id_sco214, (flux14d-flux14u)/dtbgc)  !Nic: divided by the time step to get instantaneous rate of change
         
         _SET_SURFACE_DIAGNOSTIC_(self%id_co213fxd, flux13d) ! Save up- and downward components of carbon fluxes for output
         _SET_SURFACE_DIAGNOSTIC_(self%id_co213fxu, flux13u)
         _SET_SURFACE_DIAGNOSTIC_(self%id_co214fxd, flux14d) ! Save up- and downward components of carbon fluxes for output
         _SET_SURFACE_DIAGNOSTIC_(self%id_co214fxu, flux14u)
      _SURFACE_LOOP_END_
   end subroutine do_surface   
   
   subroutine do(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_cisonew), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: 
      
      _LOOP_BEGIN_
         _GET_(self%id_dissol, dissol)
         _GET_(self%id_calc13, calc13)
         _GET_(self%id_calc14, calc14)
         _GET_(self%id_calc, calc)
         _GET_(self%id_depth, depth)
         _GET_(self%id_det14, det14)
         _GET_(self%id_det13, det13)
         _GET_(self%id_doc14, doc14)
         _GET_(self%id_doc13, doc13)
         _GET_(self%id_phy14, phy14)
         _GET_(self%id_phy13, phy13)
         _GET_(self%id_zoo14, zoo14)
         _GET_(self%id_zoo13, zoo13)
         _GET_(self%id_sco213,sco213)
         _GET_(self%id_sco214,sco214)
         _GET_(self%id_phytomi,phytomi)
         _GET_(self%id_phy,phy)
         _GET_(self%id_phosy,phosy)
         _GET_(self%id_phymor,phymor)
         _GET_(self%id_phyrem,phyrem)
         _GET_(self%id_exud,exud)
         _GET_(self%id_gratpoc,gratpoc)
         _GET_(self%id_graton,graton)
         _GET_(self%id_grawa,grawa)
         _GET_(self%id_pommor,pommor)
         _GET_(self%id_dimmor,dimmor)
         _GET_(self%id_domex,domex)
         _GET_(self%id_delcar_part,delcar_part)
         _GET_(self%id_docrem,docrem)
         _GET_(self%id_pocrem,pocrem)
         _GET_(self%id_remin2o,remin2o)
         _GET_(self%id_co2star,co2star)
         _GET_(self%id_sco212,sco212)
         _GET_(self%id_phy,   phy)
         _GET_(self%id_zoo,   zoo)
         _GET_(self%id_doc,   doc)
         _GET_(self%id_det,   det)
         
         !bgc sources/sinks
         rco213 = sco213/sco212+safediv) ! calculation of 13C and 14C equivalent of biology
         rco214 = sco214/sco212+safediv)
         rphy13 = phy13 /phy+safediv)
         rphy14 = phy14 /phy+safediv)
         rzoo13 = zoo13 /zoo+safediv)
         rzoo14 = zoo14 /zoo+safediv)
         rdoc13 = doc13 /doc+safediv)
         rdoc14 = doc14 /doc+safediv)
         rdet13 = det13 /det+safediv)
         rdet14 = det14 /det+safediv)
         if (phy < phytomi) then
             bifr13 = 1._rk
         else
             phygrowth   = (phy+phosy)/phy ! Growth rate phytoplankton [1/d]
             growth_co2  = phygrowth/(co2star*1.e6_rk+safediv)             ! CO2* in [mol/kg]
             bifr13_perm = (6.03_rk + 5.5_rk*growth_co2)/(0.225_rk + growth_co2)        ! Permil (~20)
             bifr13_perm = max(5._rk,min(26._rk,bifr13_perm))                        ! Limit the range to [5,26]
             bifr13      = (1000._rk - bifr13_perm) / 1000._rk                       ! Fractionation factor 13c (~0.98)
         endif
         bifr14 = bifr13**2
         phosy13 = phosy*bifr13*rco213
         phosy14 = phosy*bifr14*rco214
         
         if (depth<=100_rk) then
             sterzo = 0.0_rk
         else
             sterzo = dimmor+pommor
             dimmor = 0.0_rk
             pommor = 0.0_rk
         endif
         
         export13 = pommor*rzoo13 + phymor*rphy13 + gratpoc*rphy13
         export14 = pommor*rzoo14 + phymor*rphy14 + gratpoc*rphy14
         
         !rocs
         det13_roc  = phymor*rphy13 + gratpoc*rphy13 + pommor*rzoo13 - (pocrem+remin2o)*rdet13 + sterzo*rzoo13
         det14_roc  = phymor*rphy14 + gratpoc*rphy14 + pommor*rzoo14 - (pocrem+remin2o)*rdet14 + sterzo*rzoo14
         doc13_roc  = -docrem*rdoc13 + domex*rzoo13 + exud*rphy13
         doc14_roc  = -docrem*rdoc14 + domex*rzoo14 + exud*rphy14
         phy13_roc  = phosy13 - (graton+grawa+gratpoc)*rphy13 - phymor*rphy13 - exud*rphy13 - phyrem*rphy13
         phy14_roc  = phosy14 - (graton+grawa+gratpoc)*rphy14 - phymor*rphy14 - exud*rphy14 - phyrem*rphy14
         zoo13_roc  = grawa*rphy13 - domex*rzoo13 - (pommor+dimmor+sterzo)*rzoo13
         zoo14_roc  = grawa*rphy14 - domex*rzoo14 - (pommor+dimmor+sterzo)*rzoo14
         sco213_roc = - delcar_part*export13 + rcar*(docrem*rdoc13-phosy13 + graton*rphy13 + dimmor*rzoo13 + (pocrem+remin2o)*rdet13 + phyrem*rphy13)
         sco214_roc = - delcar_part*export14 + rcar*(docrem*rdoc14-phosy14 + graton*rphy14 + dimmor*rzoo14 + (pocrem+remin2o)*rdet14 + phyrem*rphy14)
         calc13_roc = delcar_part*export13
         calc14_roc = delcar_part*export14
               
         ! carbon 14 decay
         c14dec=1._rk-(log(2._rk)/self%c14_t_half)
         
         !dissolution of calcite
         dissol13=dissol*calc13/(calc+safediv)
         dissol14=dissol*calc14/(calc+safediv)
         
         sco214_roc = sco214_roc - sco214*(1._rk-c14dec) + dissol14
         det14_roc  = det14_roc  - det14 *(1._rk-c14dec) 
         calc14_roc = calc14_roc - calc14*(1._rk-c14dec) - dissol14
         doc14_roc  = doc14_roc  - doc14 *(1._rk-c14dec)
         phy14_roc  = phy14_roc  - phy14 *(1._rk-c14dec)
         zoo14_roc  = zoo14_roc  - zoo14 *(1._rk-c14dec)
         
         !update rates of change
         _ADD_SOURCE_(self%id_det13,  det13_roc  /dtbgc)
         _ADD_SOURCE_(self%id_det14,  det14_roc  /dtbgc)
         _ADD_SOURCE_(self%id_doc13,  doc13_roc  /dtbgc)
         _ADD_SOURCE_(self%id_doc14,  doc14_roc  /dtbgc)
         _ADD_SOURCE_(self%id_phy13,  phy13_roc  /dtbgc)
         _ADD_SOURCE_(self%id_phy14,  phy14_roc  /dtbgc)
         _ADD_SOURCE_(self%id_zoo13,  zoo13_roc  /dtbgc)
         _ADD_SOURCE_(self%id_zoo14,  zoo14_roc  /dtbgc)
         _ADD_SOURCE_(self%id_sco213, sco213_roc /dtbgc) 
         _ADD_SOURCE_(self%id_sco214, sco214_roc /dtbgc) 
         _ADD_SOURCE_(self%id_calc13, calc13_roc /dtbgc) 
         _ADD_SOURCE_(self%id_calc14, calc14_roc /dtbgc)          
      _LOOP_END_
   end subroutine do

   subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_ihamocc_cisonew), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
      
      integer :: 
      real(rk) :: wdust, wpoc, wopal, wcal, wnos
      
      _LOOP_BEGIN_
         _GET_(id_wdust,wdust)
         _GET_(id_wpoc,wpoc)
         _GET_(id_wopal,wopal)
         _GET_(id_wcal,wcal)
         
         _ADD_VERTICAL_VELOCITY_(self%id_det13, -wpoc)
         _ADD_VERTICAL_VELOCITY_(self%id_det14, -wpoc)

         _ADD_VERTICAL_VELOCITY_(self%id_calc13, -wcal)
         _ADD_VERTICAL_VELOCITY_(self%id_calc14, -wcal)
      _LOOP_END_
   end subroutine get_vertical_movement
   
end module ihamocc_cisonew
