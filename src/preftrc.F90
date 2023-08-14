#include "fabm_driver.h"

module ihamocc_preftrc

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_preftrc
      type (type_state_variable_id) :: id_oxygen, id_phosph, id_alkali, id_sco212, id_prefo2, id_prefpo4, id_prefalk, id_prefdic
      
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
   end type type_ihamocc_preftrc

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_preftrc), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      ! Register state dependencies
      call self%register_state_dependency(self%id_oxygen, 'oxygen', 'kmol/m^3', 'Dissolved oxygen')
      call self%register_state_dependency(self%id_alkali, 'alkali', 'kmol/m^3', 'Alkalinity')
      call self%register_state_dependency(self%id_sco212, 'sco212', 'kmol/m^3', 'Dissolved co2')
      call self%register_state_dependency(self%id_phosph, 'phosph', 'kmol/m^3', 'dissolved phosphate')

   end subroutine
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_ihamocc_preftrc), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: oxygen, alkali, sco212, phosph
      
      _LOOP_BEGIN_
         _GET_(self%id_oxygen, oxygen)
         _GET_(self%id_alkali, alkali)
         _GET_(self%id_sco212, sco212)
         _GET_(self%id_phosph, phosph)
         
         _SET_(self%id_prefo2,oxygen)
         _SET_(self%id_prefalk,alkali)
         _SET_(self%id_prefsic,sco212)
         _SET_(self%id_prefpo4,phosph)
      _LOOP_END_
   end subroutine do

   
end module ihamocc_preftrc
