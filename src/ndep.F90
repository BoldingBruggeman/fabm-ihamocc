#include "fabm_driver.h"

module ihamocc_ndep

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_ndep
      type (type_surface_dependency_id) :: id_ndep_in
      type (type_state_variable_id) :: id_alkali, id_ano3
      type (type_surface_diagnostic_variable_id) ::  id_ndep
      
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_surface
   end type type_ihamocc_ndep

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_ndep), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      ! Register state dependencies
      call self%register_state_dependency(self%id_ano3, 'ano3', 'kmol/m^3', 'dissolved nitrate')
      call self%register_state_dependency(self%id_alkali, 'alkali', 'kmol/m^3', 'Alkalinity')

      ! Register environmental dependencies
      call self%register_dependency(self%id_ndep_in, , 'ndep_in', 'kmol m-2 s-1', 'nitrogen deposition flux')
      
      ! Register diagnostics
      call self%register_diagnostic_variable(self%id_ndep, 'ndep', 'kmol/m2/s', 'nitrogen deposition flux')

   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_ndep), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: ndep, 
      
      _SURFACE_LOOP_BEGIN_
         _GET_SURFACE_(self%id_ndep_in, ndep)
         
         _ADD_SURFACE_FLUX_(self%id_ano3, ndep) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
         _ADD_SURFACE_FLUX_(self%id_alkali, -ndep) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
         _SET_SURFACE_DIAGNOSTIC_(self%id_ndep, ndep)
      _SURFACE_LOOP_END_
   end subroutine do_surface
end module ihamocc_ndep
