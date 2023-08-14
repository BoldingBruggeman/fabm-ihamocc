#include "fabm_driver.h"

module ihamocc_zooplankton

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_zooplankton
      type (type_dependency_id) :: id_ptho, id_strahl
      type (type_surface_dependency_id) :: id_
      type (type_state_variable_id) :: id_zoo, id_phy, id_silica, id_sco212, id_phosph, id_det, id_doc
      type (type_diagnostic_variable_id) :: id_phosy, id_exud, id_phymor, id_gratpoc, id_pommor, id_dimmor, id_graton, id_grawa
      
      real(rk) :: grami, grazra, bkzoo, epsher
      
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
   end type type_ihamocc_zooplankton

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_zooplankton), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      ! Register environmental dependencies
      call self%register_dependency(self%id_ptho, standard_variables%temperature)
      call self%register_dependency(self%id_phytomi, 'phytomi', 'kmol P/m3', 'minimum concentration of phytoplankton')
      
      ! Register parameters
      call self%get_parameter(self%grami,  'grami',  'kmol P/m3','minimum concentration of zooplankton', default=1.0e-10_rk) 
      call self%get_parameter(self%grazra, 'grazra', '1/d','zooplankton grazing rate', default=1.2_rk) 
      call self%get_parameter(self%bkzoo,  'bkzoo',  'kmol P/m3','zooplankton half sat. constant', default=8.e-8_rk) 
      call self%get_parameter(self%epsher, 'epsher', '-','ingestion fraction of grazing', default=0.8_rk) 
      call self%get_parameter(self%zinges, 'zinges', '-','assimilation efficiency', default=0.6_rk)
      call self%get_parameter(self%spemor, 'spemor', '1/d','quadratic mortality constant', default=3.e6_rk)
      call self%get_parameter(self%gammaz, 'gammaz', '1/d','excretion rate', default=0.06_rk)
      call self%get_parameter(self%ecan,   'ecan',   'fraction of mortality as PO_4', default=0.95_rk) 

      ! Register state variables
      call self%register_state_variable(self%id_zoo, 'zoo', 'kmol/m^3', 'zooplankton', minimum=self%grami)
      !call self%register_state_variable(self%id_fdust, 'fdust', , 'kg/m^3', 'non-aggregated dust deposition')

      ! Register environmental dependencies
      call self%register_state_dependency(self%id_phy, 'phy', 'kmol/m^3', 'phytoplankton')
      call self%register_state_dependency(self%id_doc, 'doc', 'kmol/m^3', 'Dissolved organic carbon')
      call self%register_dependency(self%id_phosy, 'phosy', 'kmol/m3/d', 'photosynthetic rate')

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_pommor, 'pommor', 'kmol/m3/d', 'zooplankton particulate export from mortality')
      call self%register_diagnostic_variable(self%id_dimmor, 'dimmor', 'kmol/m3/d', 'zooplankton dissolved inorganic export from mortality')
      call self%register_diagnostic_variable(self%id_domex, 'domex', 'kmol/m3/d', 'zooplankton dissolved organic export')
      call self%register_diagnostic_variable(self%id_graton, 'graton', 'kmol/m3/d', 'zooplankton sloppy feeding inorganic release rate')
      call self%register_diagnostic_variable(self%id_grawa, 'grawa', 'kmol/m3/d', 'zooplankton assimilation rate')
      call self%register_diagnostic_variable(self%id_gratpoc, 'gratpoc', 'kmol/m3/d', 'zooplankton sloppy feeding particulate release rate')
   end subroutine
   
   subroutine do(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_zooplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: 
      
      _LOOP_BEGIN_
         _GET_(self%id_ptho, ptho)
         !_GET_(self%id_strahl, strahl)  ! NIC: Ask Jorn about most appropriate way to handle light attenuation and -input. Here we assume strahl is the current layer swr flux in W m-2
         _GET_(self%id_phy, avphy)
         _GET_(self%id_phytomi, phytomi)
         _GET_(self%id_zoo, avgra)
         !_GET_(self%id_silica, avsil)
         !_GET_(self%id_sco212, avdic)
         !_GET_(self%id_phosph, phosph)
         !_GET_(self%id_ano3, ano3)
         !_GET_(self%id_iron, iron)
      
         !temp = min(40._rk,max(-3._rk,ptho))
         !temfa = 0.6_rk * 1.066_rk**temp
         
         grazing = MAX(0.0_rk,avgra*self%grazra*(avphy-phytomi)/(avphy+self%bkzoo)) ! NIC: Changed from BLOM-iHAMOCC, now identical to formulation in Six and Maier-Reimer (1996)
         graton = self%epsher*(1._rk-self%zinges)*grazing
         gratpoc = (1._rk-self%epsher)*grazing
         grawa = epsher*zinges*grazing
         
         zoothresh = MAX(0._rk,avgra-2._rk*self%grami))
         zoomor = self%spemor*zoothresh*zoothresh           ! *10 compared to linear in tropics (tinka)
         excdoc = self%gammaz*zoothresh                     ! excretion of doc by zooplankton
         pommor = zoomor*(1._rk-self%ecan)! + gratpoc
         dimmor = zoomor*self%ecan! + graton
         domex = exdoc
         
         _ADD_SOURCE_(self%id_zoo, (grawa-excod-zoomor)/dtbgc)
         _ADD_SOURCE_(self%id_doc, exdoc/dtbgc)
         
         _SET_DIAGNOSTIC_(self%id_pommor, pommor)
         _SET_DIAGNOSTIC_(self%id_gratpoc, gratpoc)
         _SET_DIAGNOSTIC_(self%id_dimmor, dimmor)
         _SET_DIAGNOSTIC_(self%id_domex, domex)
         _SET_DIAGNOSTIC_(self%id_graton, graton)
         _SET_DIAGNOSTIC_(self%id_grawa, grawa)
      _LOOP_END_
   end subroutine do


   
end module ihamocc_zooplankton
