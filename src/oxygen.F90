#include "fabm_driver.h"

module ihamocc_oxygen

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_oxygen
      type (type_dependency_id) :: id_psao, id_ptho
      type (type_surface_dependency_id) :: id_psicomo, id_pfu10, id_ato2, id_ppao
      type (type_state_variable_id) :: id_oxygen
      type (type_surface_diagnostic_variable_id) ::  id_bromoflx
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_surface
   end type type_ihamocc_oxygen

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_oxygen), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      ! Register state variables
      call self%register_state_variable(self%id_oxygen, 'oxygen', 'kmol/m^3', 'Dissolved oxygen')

      ! Register environmental dependencies
      call self%register_dependency(self%id_psao, standard_variables%practical_salinity)
      call self%register_dependency(self%id_ptho, standard_variables%temperature)
      call self%register_dependency(self%id_psicomo, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_pfu10, standard_variables%wind_speed)
      call self%register_dependency(self%id_ato2, standard_variables%surface_air_oxygen_concentration) ! atmospheric oxygen mixing ratio (i.e. partial presure = mixing ratio*SLP/P_0 [atm])  NOTE: variable does not exist/non-standard!
      call self%register_dependency(self%id_ppao, standard_variables%surface_air_pressure) ! surface air pressure in pascal
      
      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_oxflux, 'oxflux', 'kmol/m2/s', 'oxygen surface flux')

   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_oxygen), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: oxy, t, t2, t3, t4, tk, tk100, s, psao, ptho, oxflux, sco2, satoxy, ato2, rpp0, ppao
            
      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_oxygen, oxygen)
         _GET_(self%id_id_ptho, ptho)
         _GET_(self%id_psao, psao)
         _GET_SURFACE_(self%id_psicomo, psicomo)
         _GET_SURFACE_(self%id_ato2, ato2)
         _GET_SURFACE_(self%id_ppao, ppao)
         _GET_SURFACE_(self%id_pfu10, pfu10)
         
         
         t = min(40._rk,max(-3._rk,ptho))
         t2   = t**2
         t3   = t**3
         t4   = t**4
         tk = t + tzero
         tk100 = tk/100.0_rk
         s = min(40._rk,max( 25._rk,psao))
         sco2  = 1920.4_rk - 135.6_rk *t + 5.2122_rk*t2 - 0.10939_rk *t3 + 0.00093777_rk*t4 ! Schmidt numbers according to Wanninkhof (2014), Table 1
         
         ! solubility of O2 (Weiss, R.F. 1970, Deep-Sea Res., 17, 721-735) for moist air at 1 atm; multiplication with oxyco converts to kmol/m^3/atm
         oxy = ox0+ox1/tk100+ox2*alog(tk100)+ox3*tk100+s*(ox4+ox5*tk100+ox6*tk100**2)
		 satoxy = exp(oxy)*oxyco
         
         kwo2  = (1._rk-psicomo) * Xconvxa * pfu10**2*(660._rk/sco2)**0.5_rk 
         rpp0 = ppao/atm2pa
         
         ! Surface flux of oxygen
		 oxflux=kwo2*(oxygen-satoxy*(ato2/196800)*rpp0) ! originally multiplied by dtbgc (ts in s) to get absolute change. Removed as FABM rates-of-change has units s-1
         _SET_SURFACE_DIAGNOSTIC_(self%id_oxflux, oxflux)
         _ADD_SURFACE_FLUX_(self%id_oxygen, -oxflux) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
      _SURFACE_LOOP_END_
   end subroutine do_surface   