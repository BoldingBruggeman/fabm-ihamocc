#include "fabm_driver.h"

module ihamocc_cfc

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_cfc
      !type (type_dependency_id) :: id_psao, id_ptho, id_hi, id_Kw
      !type (type_surface_dependency_id) :: id_psicomo, id_pfu10, id_atmbromo, id_ppao
      !type (type_state_variable_id) :: id_oxygen
      !type (type_surface_diagnostic_variable_id) ::  id_bromoflx

   contains
      ! Model procedures
      procedure :: initialize
      !procedure :: do_surface
      !procedure :: do
   end type type_ihamocc_cfc

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_cfc), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      ! Register state variables
      !call self%register_state_variable(self%id_bromo, 'bromoform', 'kmol/m^3', 'Dissolved bromoform')

      ! Register environmental dependencies
      call self%register_dependency(self%id_psao, standard_variables%practical_salinity)
      call self%register_dependency(self%id_ptho, standard_variables%temperature)
      call self%register_dependency(self%id_psicomo, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_pfu10, standard_variables%wind_speed)
      !call self%register_dependency(self%id_atmbromo, standard_variables%surface_air_bromoform_concentration) 
      !call self%register_dependency(self%id_ppao, standard_variables%surface_air_pressure) ! surface air pressure in pascal
      !call self%register_dependency(self%id_hi, 'hi', 'mol/kg', 'Hydrogen ion concentration')
      !call self%register_dependency(self%id_Kw, 'kW', 'mol/kg', 'Water dissociation product')
      
      ! Register diagnostic variables
      !call self%register_diagnostic_variable(self%id_bromoflx, 'bromoflx', 'kmol/m2/s', 'Bromoform surface flux')

   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_cfc), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: t, t2, t3, tk, ptho, ppao, atmbrf, bromo, psicomo, pfu10, flx_bromo, kw_bromo, a_bromo, sch_bromo
            
      _SURFACE_LOOP_BEGIN_
         !_GET_(self%id_bromo, bromo)
         _GET_(self%id_id_ptho, ptho)
         _GET_(self%id_psao, psao)
         _GET_SURFACE_(self%id_psicomo, psicomo)
         !_GET_SURFACE_(self%id_atmbromo, atbrf)
         !_GET_SURFACE_(self%id_ppao, ppao)
         _GET_SURFACE_(self%id_pfu10, pfu10)
         
         
         t = min(40._rk,max(-3._rk,ptho))
         t2   = t**2
         t3   = t**3
         t4   = t**4
         tk = t + tzero
         tk100 = tk/100.0_rk
         s = min(40._rk,max( 25._rk,psao))

         
         sch_11= 3579.2_rk - 222.63_rk*t + 7.5749_rk*t2 - 0.14595_rk *t3 + 0.0011874_rk *t4
         sch_12= 3828.1_rk - 249.86_rk*t + 8.7603_rk*t2 - 0.1716_rk  *t3 + 0.001408_rk  *t4
         sch_sf= 3177.5_rk - 200.57_rk*t + 6.8865_rk*t2 - 0.13335_rk *t3 + 0.0010877_rk *t4
         
         
         a_11 = exp(-229.9261_rk + 319.6552_rk*(100._rk/tk) + 119.4471_rk*log(tk100)  & ! solubility of cfc11,12 (mol/(l*atm)) (Warner and Weiss 1985) and sf6 from eq. 6 of Bullister et al. (2002) These are the alpha in (1b) of the ocmpic2 howto
         &         -1.39165_rk*(tk100)**2 + s*(-0.142382_rk + 0.091459_rk*(tk100)  &
         &         -0.0157274_rk*(tk100)**2)) 
         a_12 = exp(-218.0971_rk + 298.9702_rk*(100._rk/tk) + 113.8049_rk*log(tk100)  &
         &         -1.39165_rk*(tk100)**2 + s*(-0.143566_rk + 0.091015_rk*(tk100)  &
         &         -0.0153924_rk*(tk100)**2)) 
         a_sf = exp(-80.0343_rk  + 117.232_rk *(100._rk/tk) +  29.5817_rk*log(tk100)  &
         &         +s*(0.033518_rk-0.0373942_rk*(tk100)+0.00774862_rk*(tk100)**2)) 

         a_11 = 1e-12_rk * a_11 ! conversion from mol/(l * atm) to kmol/(m3 * pptv) 
         a_12 = 1e-12_rk * a_12
         a_sf = 1e-12_rk * a_sf
      
         kw_11 = (1._rk-psicomo) * Xconvxa * pfu10**2*(660._rk/sch_11)**0.5_rk
         kw_12 = (1._rk-psicomo) * Xconvxa * pfu10**2*(660._rk/sch_12)**0.5_rk
         kw_sf = (1._rk-psicomo) * Xconvxa * pfu10**2*(660._rk/sch_sf)**0.5_rk
       
         ! Surface fluxes for CFC: eqn. (1a) in ocmip2 howto doc(hyc)
         !     flux of CFC: downward direction (mol/m**2/s)
         !      flx11=kw_11*(a_11*cfc11_atm(i,j)*ppair/p0-trc(i,j,1,1))
         !      flx12=kw_12*(a_12*cfc12_atm(i,j)*ppair/p0-trc(i,j,1,2))
         !      unit should be in [kmol cfc m-2]
         !      unit of [cfc11_atm(i,j)*ppair/p0] should be in [pptv]
         !      unit of [flx11-12] is in [kmol / m2]

      IF (pglat(i,j).GE.10) THEN
       atm_cfc11=atm_cfc11_nh
       atm_cfc12=atm_cfc12_nh
       atm_sf6=atm_sf6_nh
      ELSE IF (pglat(i,j).LE.-10) THEN
       atm_cfc11=atm_cfc11_sh
       atm_cfc12=atm_cfc12_sh
       atm_sf6=atm_sf6_sh
      ELSE
       fact=(pglat(i,j)-(-10))/20.
       atm_cfc11=fact*atm_cfc11_nh+(1-fact)*atm_cfc11_sh
       atm_cfc12=fact*atm_cfc12_nh+(1-fact)*atm_cfc12_sh
       atm_sf6=fact*atm_sf6_nh+(1-fact)*atm_sf6_sh
      ENDIF

! Use conversion of 9.86923e-6 [std atm / Pascal]
! Surface flux of cfc11
      flx11=kw_11*dtbgc*                                               &
     & (a_11*atm_cfc11*ppao(i,j)*9.86923*1e-6-ocetra(i,j,1,icfc11))
      ocetra(i,j,1,icfc11)=ocetra(i,j,1,icfc11)+flx11/pddpo(i,j,1)
! Surface flux of cfc12
      flx12=kw_12*dtbgc*                                               &
     & (a_12*atm_cfc12*ppao(i,j)*9.86923*1e-6-ocetra(i,j,1,icfc12))
      ocetra(i,j,1,icfc12)=ocetra(i,j,1,icfc12)+flx12/pddpo(i,j,1)
! Surface flux of sf6
      flxsf=kw_sf*dtbgc*                                               &
     & (a_sf*atm_sf6*ppao(i,j)*9.86923*1e-6-ocetra(i,j,1,isf6))
      ocetra(i,j,1,isf6)=ocetra(i,j,1,isf6)+flxsf/pddpo(i,j,1)
      
      atmflx(i,j,iatmf11)=flx11
       atmflx(i,j,iatmf12)=flx12
       atmflx(i,j,iatmsf6)=flxsf
         
         !_ADD_SURFACE_FLUX_(self%id_bromo, -flx_bromo) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
         !_SET_SURFACE_DIAGNOSTIC_(self%id_bromoflx, -flx_bromo)
      _SURFACE_LOOP_END_
   end subroutine do_surface
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_ihamocc_cfc), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: t, tk, hi, ah1, Kb1, rocbromo, Kb1, lsub, ah1, bromo, Kw
      
      _LOOP_BEGIN_
         !_GET_(self%id_bromo, bromo)
         !_GET_(self%id_hi, hi)
         !_GET_(self%id_Kw,Kw)

         !_ADD_SOURCE_(self%id_bromo,rocbromo)
      _LOOP_END_
   end subroutine do
