module ihamocc_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

!KB   use ihammac_sediment
   use ihamocc_oxygen
   use ihamocc_carbon
   use ihamocc_bromo
   use ihamocc_cfc
   use ihamocc_dms
   use ihamocc_alkalinization
   use ihamocc_iron
   use ihamocc_ndep
   use ihamocc_preftrc
   use ihamocc_phytoplankton
   use ihamocc_zooplankton
   use ihamocc_detritus

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: ihamocc_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
!KB         case ('sediment');            allocate(type_ihammac_sediment::model)
         case ('oxygen');            allocate(type_ihamocc_oxygen::model)
         case ('carbon');            allocate(type_ihamocc_carbon::model)
         case ('bromo');             allocate(type_ihamocc_bromo::model)
         case ('cfc');               allocate(type_ihamocc_cfc::model)
         case ('dms');               allocate(type_ihamocc_dms::model)
         case ('alkalinization');    allocate(type_ihamocc_alkalinization::model)    
         case ('iron');              allocate(type_ihamocc_iron::model)               
         case ('ndep');              allocate(type_ihamocc_ndep::model) 
         case ('preftrc');           allocate(type_ihamocc_preftrc::model)
         case ('phytoplankton');     allocate(type_ihamocc_phytoplankton::model)
         case ('zooplankton');       allocate(type_ihamocc_zooplankton::model)
         case ('detritus');          allocate(type_ihamocc_detritus::model)    
         ! Add new models here
         case default
            call self%type_base_model_factory%create(name, model)
      end select
   end subroutine create

end module
