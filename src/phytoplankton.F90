#include "fabm_driver.h"

module ihamocc_phytoplankton

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_phytoplankton
      type (type_dependency_id) :: id_ptho, id_strahl, id_depth
      type (type_surface_dependency_id) :: id_
      type (type_state_variable_id) :: id_phy, id_silica, id_sco212, id_phosph, id_det, id_oxygen, id_doc
      type (type_diagnostic_variable_id) :: id_phosy, id_exud, id_phymor
      
      real(rk) :: pi_alpha, phytomi, bkphy, dyphy, gammap, ropal, rcalc
      
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
   end type type_ihamocc_phytoplankton

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_phytoplankton), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      ! Register environmental dependencies
      call self%register_dependency(self%id_ptho, standard_variables%temperature)
      call self%register_dependency(self%id_strahl, standard_variables%downwelling_shortwave_flux)
      
      ! Register parameters
      call self%get_parameter(self%phytomi, 'phytomi', 'kmol P/m3','minimum concentration of phytoplankton', default=1.0e-11_rk) ! default value 0.02*0.4
      call self%get_parameter(self%pi_alpha, 'pi_alpha', '-','initial slope of production vs irradiance curve (alpha)', default=0.008_rk) ! default value 0.02*0.4
      call self%get_parameter(self%bkphy, 'bkphy', 'kmol P/m3','phytoplankton half sat. constant', default=4.e-8_rk) !i.e. 0.04 mmol P/m3
      call self%get_parameter(self%dyphy, 'dyphy', '1/d','phytoplankton mortality', default=0.004_rk) ! default value 0.02*0.4dyphy=0.004
      call self%get_parameter(self%gammap, 'gammap', '1/d','phytoplankton exudation rate', default=0.04_rk)
      call self%get_parameter(self%ropal, 'ropal', '-','opal to organic phosphorous production ratio', default=30._rk)
      call self%get_parameter(self%rcalc, 'rcalc', '-','calcium carbonate to organic phosphorous production ratio', default=40._rk)

      ! Register state variables
      call self%register_state_variable(self%id_phy, 'phy', 'kmol/m^3', 'phytoplankton', minimum=self%phytomi)
      !call self%register_state_variable(self%id_fdust, 'fdust', , 'kg/m^3', 'non-aggregated dust deposition')

      ! Register environmental dependencies
      call self%register_state_dependency(self%id_silica, 'silica', 'kmol/m^3', 'Silicid acid (Si(OH)4)')
      call self%register_state_dependency(self%id_sco212, 'sco212', 'kmol/m^3', 'Dissolved co2')
      call self%register_state_dependency(self%id_phosph, 'phosph', 'kmol/m^3', 'Dissolved hosphate')
      call self%register_state_dependency(self%id_ano3, 'ano3', 'kmol/m^3', 'Dissolved nitrate')
      call self%register_state_dependency(self%id_iron, 'iron', 'kmol/m^3', 'Dissolved iron')
      call self%register_state_dependency(self%id_det, 'det', 'kmol/m^3', 'detritus')
      call self%register_state_dependency(self%id_oxygen, 'oxygen', 'kmol/m^3', 'Dissolved oxygen')
      call self%register_state_dependency(self%id_doc, 'doc', 'kmol/m^3', 'Dissolved organic carbon')
      call self%register_dependency(self%id_depth, standard_variables%depth)

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_phytomi, 'phytomi', 'kmol P/m3', 'minimum concentration of phytoplankton')
      call self%register_diagnostic_variable(self%id_phosy, 'phosy', 'kmol/m3/d', 'photosynthetic rate')
      call self%register_diagnostic_variable(self%id_phymor, 'phymor', 'kmol/m3/d', 'photosynthetic mortality rate')
      call self%register_diagnostic_variable(self%id_phyrem, 'phyrem', 'kmol/m3/d', 'photosynthetic remineralization rate')
      call self%register_diagnostic_variable(self%id_exud, 'exud', 'kmol/m3/d', 'phytoplankton exudation rate')

   end subroutine
   
   subroutine do(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_phytoplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: depth
      
      _LOOP_BEGIN_
         _GET_(self%id_ptho, ptho)
         _GET_(self%id_strahl, strahl)  ! NIC: Ask Jorn about most appropriate way to handle light attenuation and -input. Here we assume strahl is the current layer swr flux in W m-2
         _GET_(self%id_phy, avphy)
         _GET_(self%id_silica, avsil)
         _GET_(self%id_sco212, avdic)
         _GET_(self%id_phosph, phosph)
         _GET_(self%id_ano3, ano3)
         _GET_(self%id_iron, iron)
         _GET_(self%id_depth, depth)
         _GET_(self%id_oxygen, oxygen)
         
         temp = min(40._rk,max(-3._rk,ptho))
         temfa = 0.6_rk * 1.066_rk**temp
         phythresh = MAX(0._rk,avphy-2._rk*self%phytomi)
         
         exud = 0.0_rk
         phosy = 0.0_rk
         sterph = 0.0_rk
         phymor = 0.0_rk
         phyrem = 0.0_rk
         if (depth<=100_rk) then
             phofa = self%pi_alpha * strahl
             
             pho = phofa * temfa / sqrt(phofa**2._rk + temfa**2._rk) 
             avanut = MIN(phosph,rnoi*ano3)
             avanfe = MIN(avanut,iron/riron)
             xa = avanfe
             xn = xa/(1._rk+pho*avphy/(xa+self%bkphy))
             phosy = MAX(0._rk,xa-xn)
             phosy = MERGE(avdic/rcar, phosy, avdic <= rcar*phosy)     ! limit phosy by available DIC
             
             phymor = self%dyphy*phythresh
             exud = self%gammap*phythresh
         else
             sterph = 0.5_rk*self%dyphy*phythresh                                ! phytoplankton to detritus
             if (oxygen > 5.0e-8_rk) then
                 phyrem = MIN(0.5_rk*self%dyphy*phythresh,       0.33_rk*oxygen/ro2ut)
             endif
         endif
         _ADD_SOURCE_(self%id_phy, (phosy-phymor-exud-sterph-phyrem)/dtbgc)
         _ADD_SOURCE_(self%id_doc, exud/dtbgc)
         
         _SET_DIAGNOSTIC_(self%id_phosy, phosy)
         _SET_DIAGNOSTIC_(self%id_phymor, phymor+sterph)
         _SET_DIAGNOSTIC_(self%id_exud, exud)
         _SET_DIAGNOSTIC_(self%id_phyrem, phyrem)
         _SET_DIAGNOSTIC_(self%id_phytomi, phytomi)
      _LOOP_END_
   end subroutine do


   
end module ihamocc_phytoplankton
