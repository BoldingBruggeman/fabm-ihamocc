#include "fabm_driver.h"

module ihamocc_dms

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_dms
      type (type_dependency_id) :: id_ptho
      type (type_surface_dependency_id) :: id_psicomo, id_pfu10
      type (type_state_variable_id) :: id_dms
      type (type_surface_diagnostic_variable_id) ::  id_atmdms

   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_surface
   end type type_ihamocc_dms

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_dms), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      ! Register state variables
      call self%register_state_variable(self%id_dms, 'dms', 'kmol/m^3', 'dimethyl sulfide concentration')
      ! Register environmental dependencies
      call self%register_dependency(self%id_ptho, standard_variables%temperature)
      call self%register_dependency(self%id_psicomo, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_pfu10, standard_variables%wind_speed)
      
      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_atmdms, 'atmdms', 'kmol/m2/s', 'dimethyl sulfide surface flux')

   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_dms), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: t, t2, t3, t4, ptho, psicomo, pfu10, schdms, kwdms, dmsflux, dms
            
      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_dms, dms)

         _GET_(self%id_id_ptho, ptho)
         _GET_SURFACE_(self%id_psicomo, psicomo)
         _GET_SURFACE_(self%id_pfu10, pfu10)
         
         t = min(40._rk,max(-3._rk,ptho))
         t2   = t**2
         t3   = t**3
         t4   = t**4

         scdms = 2855.7_rk - 177.63_rk*t + 6.0438_rk*t2 - 0.11645_rk *t3 + 0.00094743_rk*t4 
         kwdms = (1._rk-psicomo) * Xconvxa * pfu10**2*(660._rk/scdms)**0.5_rk 
         
         dmsflux = kwdms*dms ! Surface flux of dms
         
         _ADD_SURFACE_FLUX_(self%id_dms, -dmsflux) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
         _SET_SURFACE_DIAGNOSTIC_(self%id_atmdms, dmsflux)
      _SURFACE_LOOP_END_
   end subroutine do_surface
end module ihamocc_dms
