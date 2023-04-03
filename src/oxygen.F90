#include "fabm_driver.h"

module ihamocc_oxygen

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_oxygen
      type (type_dependency_id) :: id_
      type (type_surface_dependency_id) :: id_
      type (type_state_variable_id) :: id_psao, id_ptho
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_surface
   end type type_ihamocc_oxygen

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_oxygen), intent(inout), target :: self
      integer,                  intent(in)            :: configunit

      call self%get_parameter(self%OX0, 'OX0', '','Oxygen vol. solubility constant 0', default=-173.4292_rk) !VOLUMETRIC SOLUBILITY CONSTANTS FOR O2 IN ML/L from moist air at one atm total pressure.  DEEP-SEA RESEARCH, VOL. 17, 721-735.
      call self%get_parameter(self%OX1, 'OX1', '','Oxygen vol. solubility constant 1', default=249.6339_rk) !Table 2 in WEISS, R. F. (1970) THE SOLUBILITY OF NITROGEN OXYGEN AND ARGON IN WATER AND SEAWATER.
      call self%get_parameter(self%OX2, 'OX2', '','Oxygen vol. solubility constant 2', default=143.3483_rk) !DEEP-SEA RESEARCH, VOL. 17, 721-735.
      call self%get_parameter(self%OX3, 'OX3', '','Oxygen vol. solubility constant 3', default=-21.8492_rk)
      call self%get_parameter(self%OX4, 'OX4', '','Oxygen vol. solubility constant 4', default=-0.033096_rk)
      call self%get_parameter(self%OX5, 'OX5', '','Oxygen vol. solubility constant 5', default=0.014259_rk)
      call self%get_parameter(self%OX6, 'OX6', '','Oxygen vol. solubility constant 6', default=-0.0017_rk)
      
      ! Register diagnostic variables
      call self%register_state_variable(self%id_oxygen, 'oxygen', 'mmol O2/m3', 'concentration')

      ! Register environmental dependencies
      call self%register_dependency(self%id_psao, standard_variables%practical_salinity)
      call self%register_dependency(self%id_ptho, standard_variables%temperature)
      call self%register_dependency(self%id_psicomo, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_pfu10, standard_variables%wind_speed)
      call self%register_dependency(self%id_ato2, standard_variables%surface_air_oxygen_concentration) ! atmospheric oxygen mixing ratio (i.e. partial presure = mixing ratio*SLP/P_0 [atm])  NOTE: variable does not exist/non-standard!
      call self%register_dependency(self%id_ppao, standard_variables%surface_air_pressure) ! surface air pressure in pascal
      
      

      
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
         tk100 = tk/100.0
         s = min(40.,max( 25.,psao))
         sco2  = 1920.4 - 135.6 *t + 5.2122*t2 - 0.10939 *t3 + 0.00093777*t4 ! Schmidt numbers according to Wanninkhof (2014), Table 1
         
         ! solubility of O2 (Weiss, R.F. 1970, Deep-Sea Res., 17, 721-735) for moist air at 1 atm; multiplication with oxyco converts to kmol/m^3/atm
         oxy = self%ox0+self%ox1/tk100+self%ox2*alog(tk100)+self%ox3*tk100+s*(self%ox4+self%ox5*tk100+self%ox6*tk100**2)
		 satoxy = exp(oxy)*oxyco
         
         kwo2  = (1.-psicomo) * Xconvxa * pfu10**2*(660./sco2)**0.5 
         rpp0 = ppao/atm2pa
         
         ! Surface flux of oxygen
		 oxflux=kwo2*dtbgc*(oxygen-satoxy*(ato2/196800)*rpp0)
         
         _ADD_SURFACE_FLUX_(self%id_oxy, -oxflux) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
      _SURFACE_LOOP_END_
   end subroutine do_surface   
   
   ! Determine flux, inc. correction for local atmos surface pressure
      o2ex = Kwexch*(atmosp*O2sat-soxy)
