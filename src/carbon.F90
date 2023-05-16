#include "fabm_driver.h"

module ihamocc_carbon

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_carbon
      type (type_dependency_id) :: id_psao, id_ptho, id_prho, id_prb, id_silica, id_hi_in, id_pddpo, id_hi_in
      type (type_surface_dependency_id) :: id_atco2, id_pfu10, id_psicomo
      type (type_state_variable_id) :: id_sco212, id_alkali, id_calc
      type (type_diagnostic_variable_id) :: id_hi, id_co2star, id_co3, id_omegaA, id_omegaC, id_Kw
      type (type_surface_diagnostic_variable_id) ::  id_dicsat, id_co2fxd, id_co2fxu, id_pco2d, id_pco2m, id_kwco2sol, id_kwco2d, co2sold, co2solm
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_surface
      procedure :: do
   end type type_ihamocc_carbon

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_carbon), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      ! Register state variables
      call self%register_state_variable(self%id_sco212, 'sco212', 'kmol/m^3', 'Dissolved co2')
      call self%register_state_variable(self%id_alkali, 'alkali', 'kmol/m^3', 'Alkalinity')
      call self%register_state_variable(self%id_calc,   'calc',   'kmol/m^3', 'Calcium carbonate')
      
      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_Kw, 'Kw', 'mol/kg', 'Water dissociation product')
      call self%register_diagnostic_variable(self%id_hi, 'hi', 'mol/kg', 'Hydrogen ion concentration')
      call self%register_diagnostic_variable(self%id_co2star, 'co2star', 'mol/kg', 'Dissolved CO2 (CO2*)')
      call self%register_diagnostic_variable(self%id_co3, 'co3', 'kmol/m3', 'Dissolved carbonate (CO3)')
      call self%register_diagnostic_variable(self%id_dicsat, 'co3', 'kmol/m3', 'Saturated dic')
      call self%register_diagnostic_variable(self%id_omegaA, 'omegaA', '-', 'omegaA')
      call self%register_diagnostic_variable(self%id_omegaC, 'omegaC', '-', 'omegaC')
      call self%register_diagnostic_variable(self%id_co2fxd, 'co2fxd', 'kmol/m2/s', 'Downwards co2 surface flux')
      call self%register_diagnostic_variable(self%id_co2fxu, 'co2fxu', 'kmol/m2/s', 'Downwards co2 surface flux')
      call self%register_diagnostic_variable(self%id_pco2d,  'pco2d', 'microatm', 'Dry air co2 pressure')
      call self%register_diagnostic_variable(self%id_pco2m,  'pco2m', 'microatm', 'Moist air co2 pressure')
      call self%register_diagnostic_variable(self%id_kwco2sol, 'kwco2sol', 'm/s mol/kg/microatm', 'kwco2sol')
      call self%register_diagnostic_variable(self%id_kwco2d, 'kwco2d', 'm/s', 'kwco2d')
      call self%register_diagnostic_variable(self%id_co2sold, 'co2sold', 'mol/kg/atm', 'co2sold')
      call self%register_diagnostic_variable(self%id_co2solm, 'co2solm', 'mol/kg/atm', 'co2solm')
      
      ! Register environmental dependencies
      call self%register_dependency(self%id_hi_in, 'hi', 'mol/kg', 'Hydrogen ion concentration')
      call self%register_dependency(self%id_psao, standard_variables%practical_salinity)
      call self%register_dependency(self%id_ptho, standard_variables%temperature)
      call self%register_dependency(self%id_prho, standard_variables%density)
      call self%register_dependency(self%id_prb, standard_variables%pressure)
      call self%register_dependency(self%id_pddpo, standard_variables%cell_thickness)
      call self%register_dependency(self%id_pfu10, standard_variables%wind_speed)
      call self%register_dependency(self%id_psicomo, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_atco2, standard_variables%surface_air_carbon_dioxide_concentration) ! atmospheric oxygen mixing ratio (i.e. partial presure = mixing ratio*SLP/P_0 [atm])  NOTE: variable does not exist/non-standard!
      call self%register_dependency(self%id_ppao, standard_variables%surface_air_pressure) ! surface air pressure in pascal
      call self%register_dependency(self%id_silica, 'silica', 'kmol/m^3', 'Silicid acid (Si(OH)4)')
      call self%register_dependency(self%id_phosph, 'phosph', 'kmol/m^3', 'Dissolved hosphate')
      
      ! Register diagnostic variables
      
            
   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_carbon), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: t, t2, t3, t4, tk, tk100, s, psao, ptho, prho, prb, sco212, alkali, silica, phosph, hi, cu, ac, K1, K2, pco2, scco2, ppao, pfu10, kwco2, rpp0, fluxu, fluxd, ta
      
      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_id_ptho, ptho)
         _GET_(self%id_psao, psao)
         _GET_(self%id_pddpo, pddpo)
         _GET_(self%id_prho, prho)
         _GET_(self%id_prb, prb)
         _GET_(self%id_sco212, sco212)
         _GET_(self%id_alkali, alkali)
         _GET_(self%id_hi_in, hi)
         _GET_(self%id_silica, silica)
         _GET_(self%id_phosph, phosph)
         _GET_SURFACE_(self%id_atco2, atco2)
         _GET_SURFACE_(self%id_ppao, ppao)
         _GET_SURFACE_(self%id_pfu10, pfu10)
         _GET_SURFACE_(self%id_psicomo, psicomo)

         ! Carbon chemistry: Calculate equilibrium constants and solve for [H+] and
         ! carbonate alkalinity (ac)
         t    = min(40._rk,max(-3._rk,ptho))
         t2   = t**2
         t3   = t**3
         t4   = t**4
         tk   = t + tzero
         tk100= tk/100.0_rk
         s    = min(40._rk,max( 25._rk,psao))
         rrho = prho/1000.0_rk                ! seawater density [kg/m3]->[g/cm3]
         prb  = prb*10._rk  !convert from dbar to bar. ORIGINAL: ptiestu(i,j,k)*98060*1.027e-6_rk ! pressure in unit bars, 98060 = onem
   
         tc   = sco212 / rrho  ! convert to mol/kg
         ta   = alkali / rrho
         sit  = silica / rrho
         pt   = phosph / rrho
         ah1  = hi
   
         CALL CARCHM_KEQUI(t,s,prb,Kh,Khd,K1,K2,Kb,Kw,Ks1,Kf,Ksi,             &
                           K1p,K2p,K3p,Kspc,Kspa)
   
         CALL CARCHM_SOLVE(s,tc,ta,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p, &
                           ah1,ac,niter)
   
         ! Determine CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
         cu = ( 2._rk * tc - ac ) / ( 2._rk + K1 / ah1 )
  
         pco2 = cu * 1.e6_rk / Kh ! Determine CO2 pressure and fugacity (in micoatm)   NOTE: equation below for pCO2 needs requires CO2 in mol/kg

         scco2 = 2116.8_rk - 136.25_rk*t + 4.7353_rk*t2 - 0.092307_rk*t3 + 0.0007555_rk *t4 ! Schmidt numbers according to Wanninkhof (2014), Table 1

         kwco2 = (1._rk-psicomo) * Xconvxa * pfu10**2._rk*(660._rk/scco2)**0.5_rk    ! Transfer (piston) velocity kw according to Wanninkhof (2014), in units of ms-1 

         ! Ratio P/P_0, where P is the local SLP and P_0 is standard pressure (1 atm). This is
         ! used in all surface flux calculations where atmospheric concentration is given as a
         ! mixing ratio (i.e. partial presure = mixing ratio*SLP/P_0 [atm])
         rpp0 = ppao/atm2pa

         fluxd=atco2*rpp0*kwco2*dtbgc*Kh*1.e-6_rk*rrho ! Kh is in mol/kg/atm. Multiply by rrho (g/cm^3) to get fluxes in kmol/m^2   NOTE: originally multiplied by dtbgc (86400s/d). Removed as FABM rates-of-change has units s-1
         fluxu=pco2      *kwco2*dtbgc*Kh*1.e-6_rk*rrho
         fluxu=min(fluxu,fluxd-(1.e-5_rk - sco212)*pddpo) !JT set limit for CO2 outgassing to avoid negative DIC concentration, set minimum DIC concentration to 1e-5 kmol/m3 

         ! Calculate saturation DIC concentration in mixed layer
         ta = alkali / rrho
         CALL carchm_solve_DICsat(s,atco2*rpp0,ta,sit,pt,Kh,K1,K2,Kb,Kw,Ks1,Kf, &
                                 Ksi,K1p,K2p,K3p,tc_sat,niter)
         dicsat = tc_sat * rrho ! convert mol/kg to kmlo/m^3         
         
         _ADD_SURFACE_FLUX_(self%id_sco212, (fluxd-fluxu)/dtbgc)  !Nic: divided by the time step to get instantaneous rate of change
         
         _SET_SURFACE_DIAGNOSTIC_(self%id_dicsat, dicsat) !NOTE: Nic: Implemented as surface diagnostic. If required further down the water column, a subroutine with a vertical loop will be implemented.
         _SET_SURFACE_DIAGNOSTIC_(self%id_co2fxd, fluxd) ! Save up- and downward components of carbon fluxes for output
         _SET_SURFACE_DIAGNOSTIC_(self%id_co2fxu, fluxu)
         _SET_SURFACE_DIAGNOSTIC_(self%id_pco2d, cu * 1.e6 / Khd) ! Save pco2 w.r.t. dry air for output
         _SET_SURFACE_DIAGNOSTIC_(self%id_pco2m, pco2) !pCO2 wrt moist air
         _SET_SURFACE_DIAGNOSTIC_(self%id_kwco2sol, kwco2*Kh*1e-6) ! Save product of piston velocity and solubility for output
         _SET_SURFACE_DIAGNOSTIC_(self%id_kwco2d, kwco2)
         _SET_SURFACE_DIAGNOSTIC_(self%id_co2sold, Khd)
         _SET_SURFACE_DIAGNOSTIC_(self%id_co2solm, Kh)
      _SURFACE_LOOP_END_
   end subroutine do_surface   
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_ihamocc_carbon), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: t, t2, t3, t4, tk, tk100, s, psao, ptho, prho, prb, sco212, silica, phosph, alkali, hi, cu, ac, K1, K2, cb, cc, co2star, c03, omega, OmegaA, OmegaC, calc
      
      _LOOP_BEGIN_
         _GET_(self%id_id_ptho, ptho)
         _GET_(self%id_psao, psao)
         _GET_(self%id_prho, prho)
         _GET_(self%id_prb, prb)
         _GET_(self%id_sco212, sco212)
         _GET_(self%id_alkali, alkali)
         _GET_(self%id_hi_in, hi)
         _GET_(self%id_silica, silica)
         _GET_(self%id_phosph, phosph)
         _GET_(self%id_calc, calc)

         ! Carbon chemistry: Calculate equilibrium constants and solve for [H+] and
         ! carbonate alkalinity (ac)
         t    = min(40._rk,max(-3._rk,ptho))
         t2   = t**2
         t3   = t**3
         t4   = t**4
         tk   = t + tzero
         tk100= tk/100.0_rk
         s    = min(40._rk,max( 25._rk,psao))
         rrho = prho/1000.0_rk                ! seawater density [kg/m3]->[g/cm3]
         prb  = prb*10._rk  !convert from dbar to bar. ORIGINAL: ptiestu(i,j,k)*98060*1.027e-6_rk ! pressure in unit bars, 98060 = onem
   
         tc   = sco212 / rrho  ! convert to mol/kg
         ta   = alkali / rrho
         sit  = silica / rrho
         pt   = phosph / rrho
         ah1  = hi
   
         CALL CARCHM_KEQUI(t,s,prb,Kh,Khd,K1,K2,Kb,Kw,Ks1,Kf,Ksi,             &
                           K1p,K2p,K3p,Kspc,Kspa)
   
         CALL CARCHM_SOLVE(s,tc,ta,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p, &
                           ah1,ac,niter)
   
         ! Determine CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
         cu = ( 2._rk * tc - ac ) / ( 2._rk + K1 / ah1 )
         cb = K1 * cu / ah1
         cc = K2 * cb / ah1
         co2star=cu
   
         ! Carbonate ion concentration, convert from mol/kg to kmol/m^3 
         co3  = cc * rrho 
   
         ! -----------------------------------------------------------------
         ! Deep ocean processes
         omega = ( calcon * s / 35._rk ) * cc           ! Determine Omega Calcite/Aragonite and dissolution of caco3 based on OmegaC:
         OmegaA = omega / Kspa                          !   omegaC=([CO3]*[Ca])/([CO3]sat*[Ca]sat)
         OmegaC = omega / Kspc                          !   Following Sarmiento and Gruber book, assumed that [Ca]=[Ca]sat
         supsat=co3-co3/OmegaC     !   Thus, [CO3]sat=[CO3]/OmegaC. 
         undsa=MAX(0._rk,-supsat)
         dissol=MIN(undsa,0.05_rk*calc)


!      ! Save bottom level dissociation konstants for use in sediment module    NIC: These can be calculated in the carbon.f90 do_bottom subroutine
!      if( k==kbo(i,j) ) then
!        keqb( 1,i,j)  = K1
!        keqb( 2,i,j)  = K2
!        keqb( 3,i,j)  = Kb
!        keqb( 4,i,j)  = Kw
!        keqb( 5,i,j)  = Ks1
!        keqb( 6,i,j)  = Kf
!        keqb( 7,i,j)  = Ksi
!        keqb( 8,i,j)  = K1p
!        keqb( 9,i,j)  = K2p
!        keqb(10,i,j)  = K3p
!        keqb(11,i,j)  = Kspc
!      end if
         if(ah1.gt.0.) then
            _SET_DIAGNOSTIC_(self%id_hi, max(1.e-20_rk,ah1))
         else
            _SET_DIAGNOSTIC_(self%id_hi, hi)
         endif
         _SET_DIAGNOSTIC_(self%id_co2star, co2star)
         _SET_DIAGNOSTIC_(self%id_co3, co3)
         _SET_DIAGNOSTIC_(self%id_omegaA, OmegaA)
         _SET_DIAGNOSTIC_(self%id_omegaC, OmegaC)
         _SET_DIAGNOSTIC_(self%id_Kw, Kw)
         _ADD_SOURCE_(self%id_calc,-dissol/dtbgc)
         _ADD_SOURCE_(self%id_alkali,(2._rk*dissol)/dtbgc)
         _ADD_SOURCE_(self%id_sco212,dissol/dtbgc)
      _LOOP_END_
   end subroutine do

   subroutine carchm_kequi(temp,saln,prb,Kh,Khd,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,Kspc,Kspa) !Calculate equilibrium constant for the carbonate system
!     *REAL*    *temp*    - potential temperature [degr C].
!     *REAL*    *saln*    - salinity [psu].
!     *REAL*    *prb*     - pressure [bar].
!     *REAL*    *Kh*      - equilibrium constant Kh  =  [CO2]/pCO2, moist air.
!     *REAL*    *Khd*     - equilibrium constant Kh  =  [CO2]/pCO2, dry air.
!     *REAL*    *K1*      - equilibrium constant K1  = [H][HCO3]/[H2CO3].
!     *REAL*    *K2*      - equilibrium constant K2  = [H][CO3]/[HCO3].
!     *REAL*    *Kb*      - equilibrium constant Kb  = [H][BO2]/[HBO2].
!     *REAL*    *Kw*      - equilibrium constant Kw  = [H][OH].
!     *REAL*    *Ks1*     - equilibrium constant Ks1 = [H][SO4]/[HSO4].
!     *REAL*    *Kf*      - equilibrium constant Kf  = [H][F]/[HF].
!     *REAL*    *Ksi*     - equilibrium constant Ksi = [H][SiO(OH)3]/[Si(OH)4].
!     *REAL*    *K1p*     - equilibrium constant K1p = [H][H2PO4]/[H3PO4].
!     *REAL*    *K2p*     - equilibrium constant K2p = [H][HPO4]/[H2PO4].
!     *REAL*    *K3p*     - equilibrium constant K3p = [H][PO4]/[HPO4].
!     *REAL*    *Kspc*    - equilibrium constant Kspc= [Ca2+]T [CO3]T.
!     *REAL*    *Kspa*    - equilibrium constant Kspa= [Ca2+]T [CO3]T.
      IMPLICIT NONE
      REAL(rk),    INTENT(IN)    :: temp,saln,prb
      REAL(rk),    INTENT(OUT)   :: Kh,Khd,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,Kspc,Kspa

      ! Local varibles
      INTEGER                    :: js
      REAL(rk)                   :: tk,tk100,invtk,dlogtk
      REAL(rk)                   :: s,is,is2,sqrtis,s15,s2,sqrts,scl
      REAL(rk)                   :: nKhwe74,deltav,deltak,zprb,zprb2
      REAL(rk)                   :: lnkpok0(11)

      s = MAX(25._rk,saln)
      tk = temp + tzero
      tk100 = tk/100.0_rk
      invtk = 1.0_rk / tk
      dlogtk = log(tk)
      is = 19.924_rk * s / ( 1000._rk - 1.005_rk * s )
      is2 = is * is
      sqrtis = SQRT(is)
      s15    = s**1.5_rk
      s2     = s * s
      sqrts  = SQRT(s)
      scl    = s * salchl
      
      ! Kh = [CO2]/ p CO2
      ! Weiss (1974), refitted for moist air Weiss and Price (1980) [mol/kg/atm]
      nKhwe74 = ac1+ac2/tk100+ac3*log(tk100)+ac4*tk100**2+s*(bc1+bc2*tk100+bc3*tk100**2)
      Kh      = exp( nKhwe74 )
      ! Khd = [CO2]/ p CO2
      ! Weiss (1974) for dry air [mol/kg/atm]
      nKhwe74 = ad1+ad2/tk100+ad3*log(tk100)+s*(bd1+bd2*tk100+bd3*tk100**2)
      Khd     = exp( nKhwe74 )
      ! K1 = [H][HCO3]/[H2CO3]   ; K2 = [H][CO3]/[HCO3]
      ! Millero p.664 (1995) using Mehrbach et al. data on seawater scale
      K1 = 10**( -1.0_rk * ( 3670.7_rk * invtk - 62.008_rk + 9.7944_rk * dlogtk - 0.0118_rk * s + 0.000116_rk * s2 ) )
      K2 = 10**( -1.0_rk * ( 1394.7_rk * invtk + 4.777_rk - 0.0184_rk * s + 0.000118_rk * s2 ) )
      ! Kb = [H][BO2]/[HBO2] !
      ! Millero p.669 (1995) using DATA from Dickson (1990)
      Kb = exp( ( -8966.90_rk - 2890.53_rk  * sqrts - 77.942_rk  * s + 1.728_rk * s15 - 0.0996_rk * s2 ) * invtk +    &
                ( 148.0248_rk + 137.1942_rk * sqrts + 1.62142_rk * s ) +                                        &
                ( -24.4344_rk - 25.085_rk   * sqrts - 0.2474_rk  * s ) * dlogtk + 0.053105_rk * sqrts * tk )
      ! K1p = [H][H2PO4]/[H3PO4] ; K2p = [H][HPO4]/[H2PO4] ; K3p = [H][PO4]/[HPO4]
      ! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
      K1p = exp( -4576.752_rk * invtk + 115.525_rk - 18.453_rk * dlogtk + ( -106.736_rk * invtk + 0.69171_rk ) *      &
                 sqrts + ( -0.65643_rk * invtk - 0.01844_rk ) * s )
      K2p = exp( -8814.715_rk * invtk + 172.0883_rk - 27.927_rk * dlogtk + ( -160.340_rk * invtk + 1.3566_rk ) *      &
                 sqrts + ( 0.37335_rk * invtk - 0.05778_rk ) *s );
      K3p = exp( -3070.75_rk * invtk - 18.141_rk + ( 17.27039_rk * invtk + 2.81197_rk ) * sqrts + ( -44.99486_rk *    &
                 invtk - 0.09984_rk ) * s );
      ! Ksi = [H][SiO(OH)3]/[Si(OH)4]
      ! Millero p.671 (1995) using data from Yao and Millero (1995)
      Ksi = exp( -8904.2_rk * invtk + 117.385_rk - 19.334_rk * dlogtk + ( -458.79_rk * invtk + 3.5913_rk ) * sqrtis   & 
             + ( 188.74_rk * invtk - 1.5998_rk) * is + ( -12.1652_rk * invtk + 0.07871_rk) * is2 +                 &
                 log(1.0-0.001005_rk*s))
      ! Kw = [H][OH] 
      ! Millero p.670 (1995) using composite data
      Kw = exp( -13847.26_rk * invtk + 148.9652_rk - 23.6521_rk * dlogtk + ( 118.67_rk * invtk - 5.977_rk + 1.0495_rk *  &
                dlogtk ) * sqrts - 0.01615_rk * s)
      ! Ks = [H][SO4]/[HSO4]
      ! Dickson (1990, J. chem. Thermodynamics 22, 113)
      Ks1 = exp( -4276.1_rk * invtk + 141.328_rk - 23.093_rk * dlogtk + ( -13856._rk * invtk + 324.57_rk - 47.986_rk *   &
                 dlogtk ) * sqrtis + ( 35474._rk * invtk - 771.54_rk + 114.723_rk * dlogtk ) * is - 2698._rk *     &
                 invtk * is**1.5_rk + 1776._rk * invtk * is2 + log(1.0_rk - 0.001005_rk * s ) )
      ! Kf = [H][F]/[HF]
      ! Dickson and Riley (1979) -- change pH scale to total
      Kf = exp( 1590.2_rk * invtk - 12.641_rk + 1.525_rk * sqrtis + log( 1.0_rk - 0.001005_rk * s ) + log( 1.0_rk + (    &
                0.1400_rk / 96.062_rk ) * scl / Ks1 ) )
      ! Kspc (calcite)
      ! apparent solubility product of calcite : Kspc = [Ca2+]T [CO32-]T
      ! where $[]_T$ refers to the equilibrium total (free + complexed) ion concentration.
      !          Mucci 1983 mol/kg-soln
      Kspc = 10**( -171.9065_rk - 0.077993_rk * tk + 2839.319_rk / tk + 71.595_rk * log10( tk ) + ( - 0.77712_rk +    &
                   0.0028426_rk * tk + 178.34_rk / tk ) * sqrts - 0.07711_rk * s + 0.0041249_rk * s15 );
      ! Kspa (aragonite)
      ! apparent solubility product of aragonite : Kspa = [Ca2+]T [CO32-]T
      ! where $[]_T$ refers to the equilibrium total (free + complexed) ion concentration.
      !          Mucci 1983 mol/kg-soln
      Kspa = 10**( -171.945_rk - 0.077993_rk * tk + 2903.293_rk / tk  + 71.595_rk * log10( tk ) + ( -0.068393_rk +    &
                   0.0017276_rk * tk + 88.135_rk / tk ) * sqrts - 0.10018_rk * s + 0.0059415_rk * s15 );
      
      
      !---------------------- Pressure effect on Ks (Millero, 95) --------------------
      ! index: K1 1, K2 2, Kb 3, Kw 4, Ks 5, Kf 6, Kspc 7, Kspa 8, K1p 9, K2p 10, K3p 11
      DO js = 1,11
         deltav      = a0(js) + a1(js) * temp + a2(js) * temp * temp
         deltak      = b0(js) + b1(js) * temp + b2(js) * temp * temp
         zprb        = prb / ( rgas * tk )
         zprb2       = prb * zprb
         lnkpok0(js) = - ( deltav * zprb + 0.5_rk * deltak * zprb2 )
      ENDDO
      
      K1   = K1   * exp( lnkpok0(1)  )
      K2   = K2   * exp( lnkpok0(2)  )
      Kb   = Kb   * exp( lnkpok0(3)  )
      Kw   = Kw   * exp( lnkpok0(4)  )
      Ks1  = Ks1  * exp( lnkpok0(5)  )
      Kf   = Kf   * exp( lnkpok0(6)  )
      Kspc = Kspc * exp( lnkpok0(7)  )
      Kspa = Kspa * exp( lnkpok0(8)  )
      K1p  = K1p  * exp( lnkpok0(9)  )
      K2p  = K2p  * exp( lnkpok0(10) )
      K3p  = K3p  * exp( lnkpok0(11) )
   
   end subroutine carchm_kequi

   subroutine carchm_solve(saln,tc,ta,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,ah1,ac,niter) !Solve carbon chemistry.
!     *REAL*    *saln*    - salinity [psu].
!     *REAL*    *tc*      - total DIC concentraion [mol/kg].
!     *REAL*    *ta*      - total alkalinity [eq/kg].
!     *REAL*    *sit*     - silicate concentration [mol/kg].
!     *REAL*    *pt*      - phosphate concentration [mol/kg].
!     *REAL*    *K1*      - equilibrium constant K1  = [H][HCO3]/[H2CO3].
!     *REAL*    *K2*      - equilibrium constant K2  = [H][CO3]/[HCO3].
!     *REAL*    *Kb*      - equilibrium constant Kb  = [H][BO2]/[HBO2].
!     *REAL*    *Kw*      - equilibrium constant Kw  = [H][OH].
!     *REAL*    *Ks1*     - equilibrium constant Ks1 = [H][SO4]/[HSO4].
!     *REAL*    *Kf*      - equilibrium constant Kf  = [H][F]/[HF].
!     *REAL*    *Ksi*     - equilibrium constant Ksi = [H][SiO(OH)3]/[Si(OH)4].
!     *REAL*    *K1p*     - equilibrium constant K1p = [H][H2PO4]/[H3PO4].
!     *REAL*    *K2p*     - equilibrium constant K2p = [H][HPO4]/[H2PO4].
!     *REAL*    *K3p*     - equilibrium constant K3p = [H][PO4]/[HPO4].
!     *REAL*    *ah1*     - hydrogen ion concentration.
!     *REAL*    *ac*      - carbonate alkalinity.
!     *INTEGER* *niter*   - maximum number of iteration
      IMPLICIT NONE
      REAL(rk),    INTENT(IN)    :: saln,tc,ta,sit,pt
      REAL(rk),    INTENT(IN)    :: K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p
      REAL(rk),    INTENT(INOUT) :: ah1
      REAL(rk),    INTENT(OUT)   :: ac
      INTEGER, INTENT(IN)        :: niter
        
      ! Parameters to set accuracy of iteration 
      REAL(rk),    PARAMETER     :: eps=5.e-5_rk
      
      ! Local varibles
      INTEGER                    :: jit
      REAL(rk)                   :: s,scl,borat,sti,ft
      REAL(rk)                   :: hso4,hf,hsi,hpo4,ab,aw,ah2o,ah2,erel



      ! Calculate concentrations for borate, sulfate, and fluoride; see Dickson, A.G.,
      ! Sabine, C.L. and Christian, J.R. (Eds.) 2007. Guide to best practices 
      ! for ocean CO2 measurements. PICES Special Publication 3, chapter 5 p. 10
      s = MAX(25._rk,saln)
      scl = s * salchl
      borat = bor1 * scl * bor2           ! Uppstrom (1974)
      sti = 0.14_rk * scl / 96.062_rk     ! Morris & Riley (1966)
      ft = 0.000067_rk * scl / 18.9984_rk ! Riley (1965)


      iflag: DO jit = 1,niter
         hso4 = sti / ( 1._rk + Ks1 / ( ah1 / ( 1._rk + sti / Ks1 ) ) )
         hf   = 1._rk / ( 1._rk + Kf / ah1 )
         hsi  = 1._rk/ ( 1._rk + ah1 / Ksi )
         hpo4 = ( K1p * K2p * ( ah1 + 2._rk * K3p ) - ah1**3 ) /    & 
                ( ah1**3 + K1p * ah1**2 + K1p * K2p * ah1 + K1p * K2p * K3p )
         ab   = borat / ( 1._rk + ah1 / Kb )
         aw   = Kw / ah1 - ah1 / ( 1._rk + sti / Ks1 )
         ac   = ta + hso4 - sit * hsi - ab - aw + ft * hf - pt * hpo4
         ah2o = SQRT( ( tc - ac )**2 + 4._rk * ( ac * K2 / K1 ) * ( 2._rk * tc - ac ) )
         ah2  = 0.5_rk * K1 / ac *( ( tc - ac ) + ah2o )
         erel = ( ah2 - ah1 ) / ah2
         if (abs( erel ).ge.eps) then
            ah1 = ah2
         else
            exit iflag
         endif
      ENDDO iflag

   end subroutine carchm_solve

   subroutine carchm_solve_DICsat(saln,pco2,ta,sit,pt,Kh,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,tc_sat,niter) !Solve DICsat from TALK and pCO2.
!     *REAL*    *saln*    - salinity [psu].
!     *REAL*    *pco2*    - partial pressure of CO2 [ppm].
!     *REAL*    *ta*      - total alkalinity [eq/kg].
!     *REAL*    *sit*     - silicate concentration [mol/kg].
!     *REAL*    *pt*      - phosphate concentration [mol/kg].
!     *REAL*    *Kh*      - equilibrium constant K0  = [H2CO3]/pCO2.
!     *REAL*    *K1*      - equilibrium constant K1  = [H][HCO3]/[H2CO3].
!     *REAL*    *K2*      - equilibrium constant K2  = [H][CO3]/[HCO3].
!     *REAL*    *Kb*      - equilibrium constant Kb  = [H][BO2]/[HBO2].
!     *REAL*    *Kw*      - equilibrium constant Kw  = [H][OH].
!     *REAL*    *Ks1*     - equilibrium constant Ks1 = [H][SO4]/[HSO4].
!     *REAL*    *Kf*      - equilibrium constant Kf  = [H][F]/[HF].
!     *REAL*    *Ksi*     - equilibrium constant Ksi = [H][SiO(OH)3]/[Si(OH)4].
!     *REAL*    *K1p*     - equilibrium constant K1p = [H][H2PO4]/[H3PO4].
!     *REAL*    *K2p*     - equilibrium constant K2p = [H][HPO4]/[H2PO4].
!     *REAL*    *K3p*     - equilibrium constant K3p = [H][PO4]/[HPO4].
!     *REAL*    *tc_sat*  - saturated total DIC concentration [mol/kg].
!     *INTEGER* *niter*   - maximum number of iteration
      IMPLICIT NONE
      REAL(rk),    INTENT(IN)    :: saln,pco2,ta,sit,pt
      REAL(rk),    INTENT(IN)    :: Kh,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p
      REAL(rk),    INTENT(OUT)   :: tc_sat
      INTEGER, INTENT(IN)        :: niter
     
      ! Parameters to set accuracy of iteration 
      REAL(rk),    PARAMETER     :: eps=5.e-5_rk
   
      ! Local varibles
      INTEGER                    :: jit
      REAL(rk)                   :: s,scl,borat,sti,ft
      REAL(rk)                   :: hso4,hf,hsi,hpo4,ab,aw,ah2o,ah2,erel
      REAL(rk)                   :: dic_h2co3,dic_hco3,dic_co3,ah1,ac
   
      ! Calculate concentrations for borate, sulfate, and fluoride; see Dickson, A.G.,
      ! Sabine, C.L. and Christian, J.R. (Eds.) 2007. Guide to best practices 
      ! for ocean CO2 measurements. PICES Special Publication 3, chapter 5 p. 10
      s = MAX(25._rk,saln)
      scl = s * salchl
      borat = bor1 * scl * bor2            ! Uppstrom (1974)
      sti = 0.14_rk * scl / 96.062_rk      ! Morris & Riley (1966)
      ft = 0.000067_rk * scl / 18.9984_rk  ! Riley (1965)
      ah1=1.e-8_rk
      dic_h2co3 = Kh * pco2 * 1.e-6_rk 
      
      iflag: DO jit = 1,niter
         hso4 = sti / ( 1._rk + Ks1 / ( ah1 / ( 1._rk + sti / Ks1 ) ) )
         hf   = 1._rk / ( 1._rk + Kf / ah1 )
         hsi  = 1._rk/ ( 1._rk + ah1 / Ksi )
         hpo4 = ( K1p * K2p * ( ah1 + 2._rk * K3p ) - ah1**3 ) /    & 
                ( ah1**3 + K1p * ah1**2 + K1p * K2p * ah1 + K1p * K2p * K3p )
         ab   = borat / ( 1._rk + ah1 / Kb )
         aw   = Kw / ah1 - ah1 / ( 1._rk + sti / Ks1 )
         ac   = ta + hso4 - sit * hsi - ab - aw + ft * hf - pt * hpo4
         ah2o = SQRT((K1*dic_h2co3)**2 + 4._rk*ac*2._rk*K1*k2*dic_h2co3) 
         ah2  = (K1*dic_h2co3 + ah2o)/(2._rk*ac)
         erel = ( ah2 - ah1 ) / ah2
         if (abs( erel ).ge.eps) then
            ah1 = ah2
         else
            exit iflag
         endif
      ENDDO iflag
      
      dic_hco3  = Kh * K1 *      pco2 * 1.e-6_rk / ah1
      dic_co3   = Kh * K1 * K2 * pco2 * 1.e-6_rk / ah1**2
      tc_sat    = dic_h2co3 + dic_hco3 + dic_co3 
      
   end subroutine carchm_solve_DICsat
end module ihamocc_carbon
