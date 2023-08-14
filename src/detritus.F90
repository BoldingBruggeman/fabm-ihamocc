#include "fabm_driver.h"

module ihamocc_detritus

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_detritus
      type (type_dependency_id) :: id_exud, id_depth, id_ptho, id_dimmor, id_pommor, id_phymor, id_phosy, id_exud, id_graton, id_grawa, id_gratpoc, id_hi
      type (type_surface_dependency_id) :: id_
      type (type_state_variable_id) :: id_phy, id_silica, id_nos, id_oxygen, id_dms, id_sco212, id_phosph, id_ano3, id_alkali, id_calc, id_opal, id_iron, id_an2o
      type (type_diagnostic_variable_id) :: id_
      
      logical  :: AGG, WLIN, with_dmsph
      real(rk) :: remido, bkopal, ropal, rcalc, calmax, drempoc, dremopal, dremn2o, dremsul, dms_gamma, dmsp1, dmsp2, dmsp3, dmsp4, dmsp5, dmsp6, relaxfe
      
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
   end type type_ihamocc_detritus

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_detritus), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      
      ! Register parameters
      call self%get_parameter(self%dms_gamma, 'dms_gamma', '-','dms_ph scaling factor', default=0.87_rk)
      call self%get_parameter(self%dmsp6, 'dmsp6', '-','0 half saturation microbial', default=1.e-8_rk) !Parameter are a result from kettle optimisation 02.03.04
      call self%get_parameter(self%dmsp5, 'dmsp5', '-','production with delsil', default=0.025_rk) !Parameter are a result from kettle optimisation 02.03.04. Following Kloster et al., 06 Table 1, but increased by a factor of ~2
      call self%get_parameter(self%dmsp4, 'dmsp4', '-','production with delcar', default=0.125_rk) !Parameter are a result from kettle optimisation 02.03.04. Following Kloster et al., 06 Table 1, but reduced by ~7%
      call self%get_parameter(self%dmsp3, 'dmsp3', '-','dms parameter 3', default=0.0864_rk) !Following Kloster et al., 06 Table 1 with 50% reduction to reduce bacterial removal and increase dms emissions
      call self%get_parameter(self%dmsp2, 'dmsp2', '-','dms parameter 2', default=0.0011_rk) !Following Kloster et al., 06 Table 1
      call self%get_parameter(self%dmsp1, 'dmsp1', '-','dms parameter 1', default=10._rk) !2*5. production with temp
      
      call self%get_parameter(self%rdnit0, 'rdnit0', '-','moles nitrate lost for remineralisation of 1 mole P', default=0.8_rk*ro2ut) !Paulmier et al. 2009, Table 1 and equation 18. Note that their R_0=ro2ut-2*rnit.
      call self%get_parameter(self%rdnit1, 'rdnit1', '-','moles nitrate net  for remineralisation of 1 mole P', default=0.8_rk*ro2ut-rnit) !Paulmier et al. 2009, Table 1 and equation 18. Note that their R_0=ro2ut-2*rnit.
      call self%get_parameter(self%rdnit2, 'rdnit2', '-','moles N2 released  for remineralisation of 1 mole P', default=0.4_rk*ro2ut) !Paulmier et al. 2009, Table 1 and equation 18. Note that their R_0=ro2ut-2*rnit.      

      call self%get_parameter(self%rdn2o1, 'rdn2o1', '-','moles N2O used for remineralisation of 1 mole P', default=2._rk*ro2ut-2.5*rnit) !Paulmier et al. 2009, Table 1 and equation 18. Note that their R_0=ro2ut-2*rnit.      
      call self%get_parameter(self%rdn2o2, 'rdn2o2', '-','moles N2 released  for remineralisation of 1 mole P', default=2._rk*ro2ut-2._rk*rnit) !Paulmier et al. 2009, Table 1 and equation 18. Note that their R_0=ro2ut-2*rnit.      

      call self%get_parameter(self%remido, 'remido', '1/d','DOM remineralization rate', default=0.004_rk)
      call self%get_parameter(self%drempoc, 'drempoc,', '1/d','deep sea poc remineralisation rate', default=0.025_rk)
      call self%get_parameter(self%dremopal,'dremopal', '1/d','deep sea opal remineralisation rate', default=0.003_rk)
      call self%get_parameter(self%dremn2o, 'dremn2o', '1/d','deep sea n2o remineralisation rate', default=0.01_rk)
      call self%get_parameter(self%dremsul, 'dremsul', '1/d','deep sea sulphate remineralisation rate', default=0.005_rk)
      call self%get_parameter(self%bkopal, 'bkopal', 'kmol Si/m3','half sat. constant for opal', default=5.e-6_rk) !i.e. 0.04 mmol P/m3
      call self%get_parameter(self%ropal, 'ropal', '-','opal to organic phosphorous production ratio', default=10.5_rk)
      call self%get_parameter(self%rcalc, 'rcalc', '-','calcium carbonate to organic phosphorous production ratio', default=14._rk)
      call self%get_parameter(self%with_dmsph, 'with_dmsph', '-','turn on ph dependence of dms processes', default=.FALSE.)
      call self%get_parameter(self%relaxfe, 'relaxfe', '-','relaxfe', default=1.3699e-4_rk)
      
      call self%get_parameter(self%AGG, 'AGG', '-','turn on aggregations', default=.FALSE.)
      if self%AGG then
        call self%get_parameter(self%calmax, 'calmax', '-','', default=0.2_rk)
        call self%get_parameter(self%FractDim, 'FractDim', '-','-', default=1.62_rk)
        call self%get_parameter(self%cellmass, 'cellmass', 'nmol P','cellmass', default=0.012_rk/rnit)
        
        call self%register_state_variable(self%id_nos, 'nos','1/g', 'marine snow aggregates per g sea water')
      endif
      call self%get_parameter(self%WLIN, 'WLIN', '-','second sinking scheme', default=.FALSE.)
      
      

      ! Register state variables
      call self%register_state_variable(self%id_det, 'det','kmol/m^3', 'detritus')
      call self%register_state_variable(self%id_doc, 'doc','kmol/m^3', 'dissolvecd organic carbon')

      ! Register environmental dependencies
      call self%register_state_dependency(self%id_phy, 'phy', 'kmol/m^3', 'phytoplankton')
      call self%register_state_dependency(self%id_silica, 'silica', 'kmol/m^3', 'Silicid acid (Si(OH)4)')
      call self%register_state_dependency(self%id_oxygen, 'oxygen', 'kmol/m3', 'Dissolved oxygen')
      call self%register_state_dependency(self%id_id_dms, 'dms',    'kmol/m^3', 'dimethyl sulfide concentration')
      call self%register_state_dependency(self%id_sco212, 'sco212', 'kmol/m^3', 'Dissolved co2')
      call self%register_state_dependency(self%id_phosph, 'phosph', 'kmol/m^3', 'Dissolved hosphate')
      call self%register_state_dependency(self%id_ano3,   'ano3',   'kmol/m^3', 'Dissolved nitrate')
      call self%register_state_dependency(self%id_an2o,   'an2o',   'kmol/m^3', 'laughing gas')
      call self%register_state_dependency(self%id_alkali, 'alkali', 'kmol/m^3', 'Alkalinity')
      call self%register_state_dependency(self%id_calc,   'calc',   'kmol/m^3', 'Calcium carbonate')
      call self%register_state_dependency(self%id_opal,   'opal',   'kmol/m^3', 'Biogenic silica')
      call self%register_state_dependency(self%id_iron,   'iron',   'kmol/m^3', 'dissolved iron')
      


      
      call self%register_dependency(self%id_pi_ph, 'pi_pi', 'mol/kg', 'PI pH') ! NIC: this appears to be just the ph value of seawater(?)
      call self%register_dependency(self%id_hi, 'hi', 'mol/kg', 'Hydrogen ion concentration')
      call self%register_dependency(self%id_phymor, 'phymor', 'kmol/m3/d', 'photosynthetic mortality rate')
      call self%register_dependency(self%id_phyrem, 'phyrem', 'kmol/m3/d', 'photosynthetic remineralization rate')
      call self%register_dependency(self%id_pommor, 'pommor', 'kmol/m3/d', 'zooplankton particulate export from mortality')
      call self%register_dependency(self%id_dimmor, 'dimmor', 'kmol/m3/d', 'zooplankton dissolved inorganic export from mortality')
      call self%register_dependency(self%id_phosy, 'phosy', 'kmol/m3/d', 'photosynthetic rate')
      call self%register_dependency(self%id_exud, 'exud', 'kmol/m3/d', 'phytoplankton exudation rate')
      call self%register_dependency(self%id_graton, 'graton', 'kmol/m3/d', 'zooplankton sloppy feeding inorganic release rate')
      call self%register_dependency(self%id_grawa, 'grawa', 'kmol/m3/d', 'zooplankton assimilation rate')
      call self%register_dependency(self%id_gratpoc, 'gratpoc', 'kmol/m3/d', 'zooplankton sloppy feeding particulate release rate')
      call self%register_dependency(self%id_satoxy, 'satoxy', 'kmol/m^3', 'oxygen solubility')
      
      call self%register_dependency(self%id_depth, standard_variables%depth)
      call self%register_dependency(self%id_ptho, standard_variables%temperature)
      call self%register_dependency(self%id_strahl, standard_variables%downwelling_shortwave_flux)


      !call self%register_dependency(self%id_depth, standard_variables%depth)

      ! Register diagnostic variables
      !call self%register_diagnostic_variable(self%id_phytomi, 'phytomi', 'kmol P/m3', 'minimum concentration of phytoplankton')
      
   end subroutine
   
   subroutine do(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: phy, oxygen, det, doc, phymor, pomex, phosy, avsil, silica, bacfra, export, avsil, avmass, delsil, delcar, avnos
      real(rk) :: anosloss, nos, exud, graton, grawa, nos_roc, zmornos, pocrem, docrem, phyrem, pommor, dimmor, hi, dms_ph, dmsprod
      real(rk) :: dms_bac, dms_uv, strahl, dms, dissopal, opal, iron, pocrem, phyrem, doc_roc, phosph_roc, ano3_roc, det_roc, dms_roc
      real(rk) :: sco212_roc, alkali_roc, oxygen_roc, calc_roc, silica_roc, opal_roc, iron_roc  
      
      _LOOP_BEGIN_                               
         _GET_(self%id_ptho, ptho)               
         _GET_(self%id_depth, depth)             
         _GET_(self%id_phy, phy)                 
         _GET_(self%id_oxygen, oxygen)           
         _GET_(self%id_det, det)                 
         _GET_(self%id_doc, doc)                 
         _GET_(self%id_phymor, phymor)           
         _GET_(self%id_phyrem, phyrem)           
         _GET_(self%id_pommor, pommor)           
         _GET_(self%id_dimmor, dimmor)
         _GET_(self%id_phosy, phosy)
         _GET_(self%id_silica, silica)
         _GET_(self%id_opal, opal)
         _GET_(self%id_gratpoc, gratpoc)
         _GET_(self%id_iron, iron)
         _GET_(self%id_satoxy, satoxy)
         _GET_(self%id_dms, dms)
         _GET_(self%id_ano3, ano3)
         _GET_(self%id_an2o, an2o)
         
         !_GET_(self%id_ptho, ptho)
         bacfra = self%remido*doc
         export = pommor + gratpoc + phymor
         avsil = max(0.0_rk,silica)
         temp = min(40._rk,max(-3._rk,ptho))

         doc_roc     = 0.0_rk
         phosph_roc  = 0.0_rk
         ano3_roc    = 0.0_rk
         det_roc     = 0.0_rk
         dms_roc     = 0.0_rk
         sco212_roc  = 0.0_rk
         alkali_roc  = 0.0_rk
         oxygen_roc  = 0.0_rk
         calc_roc    = 0.0_rk
         silica_roc  = 0.0_rk
         opal_roc    = 0.0_rk
         iron_roc    = 0.0_rk
         gasnit_roc  = 0.0_rk
         an2o_roc    = 0.0_rk
         if (depth<=100_rk) then ! in photic zone
             _GET_(self%id_graton, graton)

             if self%AGG ! if aggregations are turned on
                 delsil = MIN(self%ropal*phosy*avsil/(avsil+self%bkopal),0.5_rk*avsil)
                 delcar = self%rcalc*MIN(self%calmax*phosy,(phosy-delsil/self%ropal))
             else
                 delsil = MIN(self%ropal*export*avsil/(avsil+self%bkopal),0.5_rk*avsil)
                 delcar = self%rcalc * export * self%bkopal/(avsil+self%bkopal)
             endif
             
             ! DMS sources/sinks
             _GET_(self%id_dms, dms)
             _GET_(self%id_strahl, strahl)  ! NIC: Ask Jorn about most appropriate way to handle light attenuation and -input. Here we assume strahl is the current layer swr flux in W m-2
             if (with_dmsph) then
                 _GET_(self%id_hi, hi)
                 _GET_(self%id_pi_ph, pi_ph)
                 dms_ph  = 1._rk + (-log10(hi) - pi_ph)*self%dms_gamma
             else
                 dms_ph  = 1._rk
             endif
             dmsprod = (self%dmsp5*delsil+self%dmsp4*delcar)*(1._rk+1._rk/(temp+self%dmsp1)**2._rk)*dms_ph
             dms_bac = self%dmsp3*abs(temp+3._rk)*dms*(dms/(self%dmsp6+dms))
             dms_uv  = self%dmsp2*strahl*dms
             
             ! net remineralization rate
             dtr = bacfra-phosy+graton+dimmor
             
             dissopal = self%dremopal*opal
             
             if self%AGG then! if aggregations are turned on
                 _GET_(self%id_nos, nos)
                 _GET_(self%id_exud, exud)
                 _GET_(self%id_graton, graton)
                 _GET_(self%id_grawa, grawa)
                 avmass = det + phy
                 
                 nos_roc  = 0.0_rk
                 if (avmass > 0._rk) then
                     avnos = nos
                     anosloss = (phosy-exud-graton-grawa)*avnos/avmass
                     nos_roc = anosloss
                 endif
                 zdis = 0.01_rk / ((self%FractDim + 0.01_rk)*self%cellmass)
                 
                 zmornos = pommor * zdis * 1.e+6_rk
                 nos_roc = nosroc + zmornos
             endif
             
             !
             doc_roc    = doc_roc    - bacfra
             phosph_roc = phosph_roc + dtr
             ano3_roc   = ano3_roc   + dtr*rnit
             det_roc    = det_roc    + export
             dms_roc    = dms_roc    + dmsprod - dms_bac - dms_uv
             sco212_roc = sco212_roc + rcar*dtr - delcar
             alkali_roc = alkali_roc - 2._rk*delcar - (rnit + 1._rk)*dtr
             oxygen_roc = oxygen_roc - dtr*ro2ut
             calc_roc   = calc_roc   + delcar
             silica_roc = silica_roc + dissopal - delsil
             opal_roc   = opal_roc   + delsil - dissopal
             iron_roc   = iron_roc   + dtr*riron
         else ! below mixed layer/photic zone
             if (oxygen > 5.0e-8_rk) then
                 pocrem = MIN(self%drempoc*det,0.33_rk*oxygen/ro2ut)
                 docrem = MIN(self%remido*doc,0.33_rk*oxygen/ro2ut)
             else
                 pocrem = 0._rk
                 docrem = 0._rk
             endif
    
             sterzo = pommor + dimmor
             remin = pocrem + docrem + phyrem
             
             opalrem = self%dremopal*0.1_rk*(temp+3._rk)*opal
             
             aou = satoxy-oxygen
             refra = 1._rk+3._rk*(0.5_rk+sign(0.5_rk,aou-1.97e-4_rk))
             dms_bac = self%dmsp3 * abs(temp+3._rk) * dms * dms / (self%dmsp6+dms)
             
             if self%AGG then! if aggregations are turned on
                 _GET_(self%id_nos, nos)
                 avmass = det + phy
                 if (avmass > 0._rk) then
                     avnos = nos
                     nos_roc = -remin*avnos/avmass
                 endif
                 sterzo = pommor + dimmor
                 zmornos = sterzo * zdis * 1.e+6_rk
                 nos_roc = nos_roc + zmornos
             endif

             !
             dms_roc    = dms_roc    - dms_bac
             gasnit_roc = gasnit_roc - remin*1.e-4_rk*ro2ut*refra
             an2o_roc   = an2o_roc   + remin*1.e-4_rk*ro2ut*refra
             det_roc    = det_roc    - pocrem + phymor + sterzo
             doc_roc    = doc_roc    - docrem
             phosph_roc = phosph_roc + remin
             ano3_roc   = ano3_roc   + remin*rnit
             sco212_roc = sco212_roc + rcar*remin
             alkali_roc = alkali_roc - (rnit + 1._rk)*remin
             oxygen_roc = oxygen_roc - ro2ut*remin - remin*1.e-4_rk*ro2ut*refra*0.5_rk
             iron_roc   = iron_roc   + remin*riron
             opal_roc   = opal_roc   - opalrem
             silica_roc = silica_roc + opalrem
             
             rdnit0 = 0.8_rk*ro2ut
             rdnit1 = 0.8_rk*ro2ut-rnit
             rdnit2 = 0.4_rk*ro2ut
             rdn2o1 = 2._rk*ro2ut-2.5*rnit
             rdn2o2 = 2._rk*ro2ut-2._rk*rnit
             if (oxygen < 5.0e-7_rk) then
                 remin   = 0.05_rk * self%drempoc * MIN(det,0.5_rk * ano3 / rdnit1)
                 remin2o = self%dremn2o * MIN(det,0.003_rk * an2o / rdn2o1)
                 
                 if self%AGG then! if aggregations are turned on
                     _GET_(self%id_nos, nos)
                     avmass = det + phy
                     if (avmass > 0.) then
                         avnos = nos
                         nos_roc = nos_roc - (remin + remin2o)*avnos/avmass
                     endif
                 endif
                 !
                 alkali_roc = alkali_roc + (rdnit1-1._rk)*remin-remin2o
                 sco212_roc = sco212_roc + rcar*(remin+remin2o)
                 det_roc    = det_roc    - (remin+remin2o)
                 phosph_roc = phosph_roc + (remin+remin2o)
                 ano3_roc   = ano3_roc   - rdnit1*remin
                 gasnit_roc = gasnit_roc + rdnit2*remin+rdn2o2*remin2o
                 an2o_roc   = an2o_roc   - rdn2o1*remin2o
                 iron_roc   = iron_roc   + riron*(remin+remin2o)
                 if (ano3 < 3.0e-6_rk) then
                     remin = self%dremsul*det
                     if self%AGG then! if aggregations are turned on
                         _GET_(self%id_nos, nos)
                         avmass = det + phy
                         if (avmass > 0.) then
                             avnos = nos
                             nos_roc = nos_roc - remin*avnos/avmass
                         endif
                     endif
                     
                     !
                     det_roc    = det_roc    - remin
                     alkali_roc = alkali_roc - (rnit+1._rk)*remin
                     sco212_roc = sco212_roc + rcar*remin
                     phosph_roc = phosph_roc + remin
                     ano3_roc   = ano3_roc   + rnit*remin
                     iron_roc   = iron_roc   + riron*remin
                 endif
             endif
         endif
    
         if self%AGG then! if aggregations are turned on
             _GET_(self%id_nos, nos)
             avmass = det + phy
             snow = avmass*1.e6_rk
             if (avmass > 0.) then
                 avnos = nos
                 nos_roc = nos_roc - (remin + remin2o)*avnos/avmass
             endif
         endif

    
                 !!!!!!! REMEMBER TO ADD BIOLOGICAL BROMO SOURCES/SINKS !!!!!!! but where? 

    
             ! update state variables
             _ADD_SOURCE_(self%id_doc, -bacfra/dtbgc)
             _ADD_SOURCE_(self%id_phosph, dtr/dtbgc)
             _ADD_SOURCE_(self%id_ano3, dtr*rnit/dtbgc)
             _ADD_SOURCE_(self%id_det, export/dtbgc)
             _ADD_SOURCE_(self%id_dms, (dmsprod-dms_bac-dms_uv)/dtbgc)
             _ADD_SOURCE_(self%id_sco212, (-delcar+rcar*dtr)/dtbgc)
             _ADD_SOURCE_(self%id_alkali, (-2._rk*delcar-(rnit+1._rk)*dtr)/dtbgc)
             _ADD_SOURCE_(self%id_oxygen, -dtr*ro2ut/dtbgc)
             _ADD_SOURCE_(self%id_calc, delcar/dtbgc)
             
             _ADD_SOURCE_(self%id_silica, (dissopal-delsil)/dtbgc)
             _ADD_SOURCE_(self%id_opal, (delsil-dissopal)/dtbgc)
             _ADD_SOURCE_(self%id_iron, rociron/dtbgc)
             
             
             
    
        
             
             

         
       
      
             
             
             
             
             
             
      
      _LOOP_END_
   end subroutine do


   
end module ihamocc_detritus
