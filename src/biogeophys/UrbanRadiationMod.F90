module UrbanRadiationMod

#include "shr_assert.h"

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Calculate solar and longwave radiation, and turbulent fluxes for urban landunit
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_sys_mod       , only : shr_sys_flush 
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type, subgrid_level_landunit
  use clm_varpar        , only : numrad
  use clm_varctl        , only : iulog
  use abortutils        , only : endrun  
  use UrbanParamsType   , only : urbanparams_type
  use atm2lndType       , only : atm2lnd_type
  use WaterDiagnosticBulkType    , only : waterdiagnosticbulk_type
  use TemperatureType   , only : temperature_type
  use SolarAbsorbedType , only : solarabs_type 
  use SurfaceAlbedoType , only : surfalb_type
  use EnergyFluxType    , only : energyflux_type
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  use PatchType         , only : patch                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanRadiation    ! Urban physics - radiative fluxes
  !
  ! PRIVATE MEMBER FUNCTIONS
  private :: net_longwave     ! Net longwave radiation for road and both walls in urban canyon 

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine UrbanRadiation (bounds                                   , &
       num_nourbanl, filter_nourbanl                                  , &
       num_urbanl, filter_urbanl                                      , &
       num_urbanc, filter_urbanc                                      , &
       num_urbanp, filter_urbanp                                      , &
       atm2lnd_inst, waterdiagnosticbulk_inst, temperature_inst, urbanparams_inst, &
       solarabs_inst, surfalb_inst, energyflux_inst)
    !
    ! !DESCRIPTION: 
    ! Solar fluxes absorbed and reflected by roof and canyon (walls, road).
    ! Also net and upward longwave fluxes.

    ! !USES:
    use clm_varcon          , only : spval, sb, tfrz
    use column_varcon       , only : icol_road_perv, icol_road_imperv, icol_road_tree
    use column_varcon       , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_time_manager    , only : get_step_size_real
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds    
    integer                , intent(in)    :: num_nourbanl       ! number of non-urban landunits in clump
    integer                , intent(in)    :: filter_nourbanl(:) ! non-urban landunit filter
    integer                , intent(in)    :: num_urbanl         ! number of urban landunits in clump
    integer                , intent(in)    :: filter_urbanl(:)   ! urban landunit filter
    integer                , intent(in)    :: num_urbanc         ! number of urban columns in clump
    integer                , intent(in)    :: filter_urbanc(:)   ! urban column filter
    integer                , intent(in)    :: num_urbanp         ! number of urban patches in clump
    integer                , intent(in)    :: filter_urbanp(:)   ! urban pft filter
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(waterdiagnosticbulk_type)  , intent(in)    :: waterdiagnosticbulk_inst
    
!-------------------[kz.8]Ray tracing test------------------------- 
    type(temperature_type) , intent(inout)    :: temperature_inst
!-------------------[kz.8]Ray tracing test-------------------------     
    
    type(urbanparams_type) , intent(in)    :: urbanparams_inst
    type(solarabs_type)    , intent(inout) :: solarabs_inst
    type(surfalb_type)     , intent(in)    :: surfalb_inst
    type(energyflux_type)  , intent(inout) :: energyflux_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp,fl,p,c,l,g              ! indices
    real(r8) :: dtime                      ! land model time step (sec)


    real(r8), parameter :: mpe    = 1.e-06_r8 ! prevents overflow for division by zero
    real(r8), parameter :: snoem  = 0.97_r8   ! snow emissivity (should use value from Biogeophysics1)

    real(r8) :: t_roof(bounds%begl:bounds%endl)         ! roof temperature (K)
    real(r8) :: t_improad(bounds%begl:bounds%endl)      ! imppervious road temperature (K)
    real(r8) :: t_perroad(bounds%begl:bounds%endl)      ! pervious road temperature (K)
    real(r8) :: t_sunwall(bounds%begl:bounds%endl)      ! sunlit wall temperature (K)
    real(r8) :: t_shadewall(bounds%begl:bounds%endl)    ! shaded wall temperature (K)
    real(r8) :: t_tree_perroad(bounds%begl:bounds%endl)      ! pervious road temperature (K)
    real(r8) :: t_br_tree(bounds%begl:bounds%endl)         ! below-roof tree temperature (K)
    real(r8) :: t_ar_tree(bounds%begl:bounds%endl)         ! above-roof tree temperature (K)

    real(r8) :: lwdown(bounds%begl:bounds%endl)         ! atmospheric downward longwave radiation (W/m**2)
    real(r8) :: em_roof_s(bounds%begl:bounds%endl)      ! roof emissivity with snow effects
    real(r8) :: em_improad_s(bounds%begl:bounds%endl)   ! impervious road emissivity with snow effects
    real(r8) :: em_perroad_s(bounds%begl:bounds%endl)   ! pervious road emissivity with snow effects
    real(r8) :: em_tree_urb_s(bounds%begl:bounds%endl)   ! road tree emissivity with snow effects
    real(r8) :: em_br_tree(bounds%begl:bounds%endl)   ! below-roof tree emissivity 
    real(r8) :: em_ar_tree(bounds%begl:bounds%endl)   ! above-roof tree emissivity 
    
    real(r8) :: em_br_tree_s(bounds%begl:bounds%endl)   ! below-roof tree emissivity with snow effects
    real(r8) :: em_ar_tree_s(bounds%begl:bounds%endl)   ! above-roof tree emissivity with snow effects
    !-----------------------------------------------------------------------

    associate(                                                                 & 
         ctype              =>    col%itype                                  , & ! Input:  [integer (:)    ]  column type                                        
         coli               =>    lun%coli                                   , & ! Input:  [integer (:)    ]  beginning column index for landunit                
         colf               =>    lun%colf                                   , & ! Input:  [integer (:)    ]  ending column index for landunit                   
         canyon_hwr         =>    lun%canyon_hwr                             , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width          
         ht_roof            => lun%ht_roof                       , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width          
         wtlunit_roof         => lun%wtlunit_roof                       , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width          
         wtroad_perv        =>    lun%wtroad_perv                            , & ! Input:  [real(r8) (:)   ]  weight of pervious road wrt total road            
         wtroad_tree        =>    lun%wtroad_tree                            , & ! Input:  [real(r8) (:)   ]  weight of road tree wrt total road            
         A_v1            => lun%A_v1                       , & ! Input:  [real(r8) (:)   ]  Leaf area for urban tree canopy below roof   
         A_v2         => lun%A_v2                       , & ! Input:  [real(r8) (:)   ]  Leaf area for urban tree canopy above roof   

         forc_solad         =>    atm2lnd_inst%forc_solad_not_downscaled_grc          , & ! Input:  [real(r8) (:,:) ]  direct beam radiation  (vis=forc_sols , nir=forc_soll ) (W/m**2)
         forc_solai         =>    atm2lnd_inst%forc_solai_grc                , & ! Input:  [real(r8) (:,:) ]  diffuse beam radiation (vis=forc_sols , nir=forc_soll ) (W/m**2)
         forc_solar         =>    atm2lnd_inst%forc_solar_not_downscaled_grc          , & ! Input:  [real(r8) (:)   ]  incident solar radiation (W/m**2)                 
         forc_lwrad         =>    atm2lnd_inst%forc_lwrad_not_downscaled_grc , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)   

         frac_sno           =>    waterdiagnosticbulk_inst%frac_sno_col               , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)       

         t_ref2m            =>    temperature_inst%t_ref2m_patch             , & ! Input:  [real(r8) (:)   ]  2 m height surface air temperature (K)            
         t_grnd             =>    temperature_inst%t_grnd_col                , & ! Input:  [real(r8) (:)   ]  ground temperature (K)                            

         em_roof            =>    urbanparams_inst%em_roof                   , & ! Input:  [real(r8) (:)   ]  roof emissivity                                   
         em_improad         =>    urbanparams_inst%em_improad                , & ! Input:  [real(r8) (:)   ]  impervious road emissivity                        
         em_perroad         =>    urbanparams_inst%em_perroad                , & ! Input:  [real(r8) (:)   ]  pervious road emissivity                          
         em_br_tree         =>    urbanparams_inst%em_tree_urb                , & ! Input:  [real(r8) (:)   ]  below-roof tree emissivity                          
         em_ar_tree         =>    urbanparams_inst%em_tree_urb                , & ! Input:  [real(r8) (:)   ]  above-roof tree emissivity                          
         em_wall            =>    urbanparams_inst%em_wall                   , & ! Input:  [real(r8) (:)   ]  wall emissivity                                   
         tree_elaisai_per_uroad       => urbanparams_inst%tree_elaisai_per_uroad        , & ! Input: [real(r8) (:,:) ]
         albd               =>    surfalb_inst%albd_patch                    , & ! Input:  [real(r8) (:,:) ] pft surface albedo (direct)                         
         albi               =>    surfalb_inst%albi_patch                    , & ! Input:  [real(r8) (:,:) ] pft surface albedo (diffuse)                        

         sabg               =>    solarabs_inst%sabg_patch                   , & ! Output: [real(r8) (:)   ]  solar radiation absorbed by ground (W/m**2)        
         sabv               =>    solarabs_inst%sabv_patch                   , & ! Output: [real(r8) (:)   ]  solar radiation absorbed by vegetation (W/m**2)   
         fsa                =>    solarabs_inst%fsa_patch                    , & ! Output: [real(r8) (:)   ]  solar radiation absorbed (total) (W/m**2)         
         fsa_u              =>    solarabs_inst%fsa_u_patch                  , & ! Output: [real(r8) (:)   ]  urban solar radiation absorbed (total) (W/m**2)   

         eflx_lwrad_out     =>    energyflux_inst%eflx_lwrad_out_patch       , & ! Output: [real(r8) (:)   ]  emitted infrared (longwave) radiation (W/m**2)    
         eflx_lwrad_net     =>    energyflux_inst%eflx_lwrad_net_patch       , & ! Output: [real(r8) (:)   ]  net infrared (longwave) rad (W/m**2) [+ = to atm] 
         eflx_lwrad_net_u   =>    energyflux_inst%eflx_lwrad_net_u_patch     , & ! Output: [real(r8) (:)   ]  urban net infrared (longwave) rad (W/m**2) [+ = to atm]

         sabs_roof_dir      =>    solarabs_inst%sabs_roof_dir_lun      , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by roof per unit ground area per unit incident flux
         sabs_roof_dif      =>    solarabs_inst%sabs_roof_dif_lun      , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by roof per unit ground area per unit incident flux
         sabs_sunwall_dir   =>    solarabs_inst%sabs_sunwall_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_sunwall_dif   =>    solarabs_inst%sabs_sunwall_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_shadewall_dir =>    solarabs_inst%sabs_shadewall_dir_lun , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_shadewall_dif =>    solarabs_inst%sabs_shadewall_dif_lun , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_improad_dir   =>    solarabs_inst%sabs_improad_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_improad_dif   =>    solarabs_inst%sabs_improad_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_perroad_dir   =>    solarabs_inst%sabs_perroad_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
         sabs_perroad_dif   =>    solarabs_inst%sabs_perroad_dif_lun    ,& ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux
         sabs_br_tree_dir      =>    solarabs_inst%sabs_br_tree_dir_lun      , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by below-roof treeetation per unit treeetation area per unit incident flux
         sabs_br_tree_dif      =>    solarabs_inst%sabs_br_tree_dif_lun      , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by below-roof treeetation per unit treeetation area per unit incident flux
         sabs_ar_tree_dir      =>    solarabs_inst%sabs_ar_tree_dir_lun      , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by above-roof treeetation per unit treeetation area per unit incident flux
         sabs_ar_tree_dif      =>    solarabs_inst%sabs_ar_tree_dif_lun      , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by above-roof treeetation per unit treeetation area per unit incident flux
         sabs_tree_dir   =>    solarabs_inst%sabs_tree_dir_lun         , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by road tree per unit ground area per unit incident flux
         sabs_tree_dif   =>    solarabs_inst%sabs_tree_dif_lun         , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by road tree per unit ground area per unit incident flux
         
         lwnet_roof   => solarabs_inst%lwnet_roof_lun   , & ! Output: [real(r8) (:) ]  net (outgoing-incoming) longwave radiation (per unit ground area), roof (W/m**2)
         lwnet_improad     => solarabs_inst%lwnet_improad_lun     , & ! Output: [real(r8) (:) ]  net (outgoing-incoming) longwave radiation (per unit ground area), impervious road (W/m**2)
         lwnet_perroad     => solarabs_inst%lwnet_perroad_lun     , & ! Output: [real(r8) (:) ]  net (outgoing-incoming) longwave radiation (per unit ground area), pervious road (W/m**2)
         lwnet_tree_perroad     => solarabs_inst%lwnet_tree_perroad_lun     , & ! Output: [real(r8) (:) ]  net (outgoing-incoming) longwave radiation (per unit ground area), pervious road beneath the tree canopy (W/m**2)
         lwnet_sunwall     => solarabs_inst%lwnet_sunwall_lun     , & ! Output: [real(r8) (:) ]  net (outgoing-incoming) longwave radiation (per unit wall area), sunlit wall (W/m**2)
         lwnet_shadewall   => solarabs_inst%lwnet_shadewall_lun   , & ! Output: [real(r8) (:) ]  net (outgoing-incoming) longwave radiation (per unit wall area), shaded wall (W/m**2)
         lwnet_br_tree     => solarabs_inst%lwnet_br_tree_lun     , & ! Output: [real(r8) (:) ]  net (outgoing-incoming) longwave radiation (per unit leaf area), below-roof tree (W/m**2)
         lwnet_ar_tree     => solarabs_inst%lwnet_ar_tree_lun     , & ! Output: [real(r8) (:) ]  net (outgoing-incoming) longwave radiation (per unit leaf area), above-roof tree (W/m**2)
         lwnet_tree     => solarabs_inst%lwnet_tree_lun     , & ! Output: [real(r8) (:) ]  net (outgoing-incoming) longwave radiation (per unit leaf area), tree (W/m**2)
         lwnet_canyon      => solarabs_inst%lwnet_canyon_lun      , & ! Output: [real(r8) (:) ]  net (outgoing-incoming) longwave radiation for canyon, per unit ground area (W/m**2)
         lwup_roof    => solarabs_inst%lwup_roof_lun    , & ! Output: [real(r8) (:) ]  upward longwave radiation (per unit ground area), roof (W/m**2)
         lwup_improad      => solarabs_inst%lwup_improad_lun      , & ! Output: [real(r8) (:) ]  upward longwave radiation (per unit ground area), impervious road (W/m**2)
         lwup_perroad      => solarabs_inst%lwup_perroad_lun      , & ! Output: [real(r8) (:) ]  upward longwave radiation (per unit ground area), pervious road (W/m**2)
         lwup_tree_perroad      => solarabs_inst%lwup_tree_perroad_lun      , & ! Output: [real(r8) (:) ]  upward longwave radiation (per unit ground area), pervious road beneath the tree canopy (W/m**2)
         lwup_sunwall      => solarabs_inst%lwup_sunwall_lun      , & ! Output: [real(r8) (:) ]  upward longwave radiation (per unit wall area), sunlit wall (W/m**2)
         lwup_shadewall    => solarabs_inst%lwup_shadewall_lun    , & ! Output: [real(r8) (:) ]  upward longwave radiation (per unit wall area), shaded wall (W/m**2)
         lwup_br_tree      => solarabs_inst%lwup_br_tree_lun      , & ! Output: [real(r8) (:) ]  upward longwave radiation (per unit leaf area), below-roof tree (W/m**2)
         lwup_ar_tree      => solarabs_inst%lwup_ar_tree_lun      , & ! Output: [real(r8) (:) ]  upward longwave radiation (per unit leaf area), above-roof tree (W/m**2)
         lwup_tree      => solarabs_inst%lwup_tree_lun      , & ! Output: [real(r8) (:) ]  upward longwave flux at above-roof tree
         lwup_canyon       => solarabs_inst%lwup_canyon_lun       , & ! Output: [real(r8) (:) ]  upward longwave radiation for canyon, per unit ground area (W/m**2)
         t_veg                  => temperature_inst%t_veg_patch                 , & ! Input: [real(r8) (:)   ]  vegetation temperature (Kelvin)  
         begl               =>    bounds%begl                                , &
         endl               =>    bounds%endl                                , &
!-------------------[kz.9]Ray tracing test------------------------- 
         fww1d_u             =>    urbanparams_inst%fww1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         fww1d_t1             =>    temperature_inst%fww1d_out1                , & ! Output:  [real(r8) (:,:)   ]
         fww1d_t2             =>    temperature_inst%fww1d_out2                , & ! Output:  [real(r8) (:,:)   ]  
         
         fvv1d_u             =>    urbanparams_inst%fvv1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         fvv1d_t1             =>    temperature_inst%fvv1d_out1                , & ! Output:  [real(r8) (:,:)   ] 
         fvv1d_t2             =>    temperature_inst%fvv1d_out2                , & ! Output:  [real(r8) (:,:)   ] 
         
         fwv1d_u             =>    urbanparams_inst%fwv1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         fwv1d_t1             =>    temperature_inst%fwv1d_out1                , & ! Output:  [real(r8) (:,:)   ]  
         fwv1d_t2             =>    temperature_inst%fwv1d_out2                , & ! Output:  [real(r8) (:,:)   ]  
         fvw1d_u             =>    urbanparams_inst%fvw1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         fvw1d_t1             =>    temperature_inst%fvw1d_out1                , & ! Output:  [real(r8) (:,:)   ]  
         fvw1d_t2             =>    temperature_inst%fvw1d_out2                , & ! Output:  [real(r8) (:,:)   ]  
         
         fwr1d_u             =>    urbanparams_inst%fwr1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         fwr1d_t1             =>    temperature_inst%fwr1d_out1                , & ! Output:  [real(r8) (:,:)   ]  
         fwr1d_t2             =>    temperature_inst%fwr1d_out2                , & ! Output:  [real(r8) (:,:)   ]  

         frw1d_u             =>    urbanparams_inst%frw1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         frw1d_t1             =>    temperature_inst%frw1d_out1                , & ! Output:  [real(r8) (:,:)   ]  
         frw1d_t2             =>    temperature_inst%frw1d_out2                , & ! Output:  [real(r8) (:,:)   ]  

         fvr1d_u             =>    urbanparams_inst%fvr1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         fvr1d_t1             =>    temperature_inst%fvr1d_out1                , & ! Output:  [real(r8) (:,:)   ]  
         fvr1d_t2             =>    temperature_inst%fvr1d_out2                , & ! Output:  [real(r8) (:,:)   ]  

         frv1d_u             =>    urbanparams_inst%frv1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         frv1d_t1             =>    temperature_inst%frv1d_out1                , & ! Output:  [real(r8) (:,:)   ]  
         frv1d_t2             =>    temperature_inst%frv1d_out2               , & ! Output:  [real(r8) (:,:)   ]  

         fwg1d_u             =>    urbanparams_inst%fwg1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         fwg1d_t             =>    temperature_inst%fwg1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         fgw1d_u             =>    urbanparams_inst%fgw1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         fgw1d_t             =>    temperature_inst%fgw1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         fgv1d_u             =>    urbanparams_inst%fgv1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         fgv1d_t             =>    temperature_inst%fgv1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         fsw1d_u             =>    urbanparams_inst%fsw1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         fsw1d_t             =>    temperature_inst%fsw1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         fvg1d_u             =>    urbanparams_inst%fvg1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         fvg1d_t             =>    temperature_inst%fvg1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         fsr1d_u             =>    urbanparams_inst%fsr1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         fsr1d_t             =>    temperature_inst%fsr1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         fsv1d_u             =>    urbanparams_inst%fsv1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         fsv1d_t             =>    temperature_inst%fsv1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         fws1d_u             =>    urbanparams_inst%fws1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         fws1d_t             =>    temperature_inst%fws1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         fvs1d_u             =>    urbanparams_inst%fvs1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         fvs1d_t             =>    temperature_inst%fvs1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         frs1d_u             =>    urbanparams_inst%frs1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         frs1d_t             =>    temperature_inst%frs1d_out                , & ! Output:  [real(r8) (:,:)   ]  

         fsg1d_u             =>    urbanparams_inst%fsg1d_out                , & ! Input:  [real(r8) (:)   ]  
         fsg1d_t             =>    temperature_inst%fsg1d_out                , & ! Output:  [real(r8) (:)   ]  
         fts1d_u             =>    urbanparams_inst%fts1d_out                , & ! Input:  [real(r8) (:)   ]  
         fts1d_t             =>    temperature_inst%fts1d_out                , & ! Output:  [real(r8) (:)   ]  

         kww1d_u             =>    urbanparams_inst%kww1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         kww1d_t1             =>    temperature_inst%kww1d_out1                , & ! Output:  [real(r8) (:,:)   ]
         kww1d_t2             =>    temperature_inst%kww1d_out2                , & ! Output:  [real(r8) (:,:)   ]  
         
         kvv1d_u             =>    urbanparams_inst%kvv1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         kvv1d_t1             =>    temperature_inst%kvv1d_out1                , & ! Output:  [real(r8) (:,:)   ] 
         kvv1d_t2             =>    temperature_inst%kvv1d_out2                , & ! Output:  [real(r8) (:,:)   ] 
         
         kwv1d_u             =>    urbanparams_inst%kwv1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         kwv1d_t1             =>    temperature_inst%kwv1d_out1                , & ! Output:  [real(r8) (:,:)   ]  
         kwv1d_t2             =>    temperature_inst%kwv1d_out2                , & ! Output:  [real(r8) (:,:)   ]  
         kvw1d_u             =>    urbanparams_inst%kvw1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         kvw1d_t1             =>    temperature_inst%kvw1d_out1                , & ! Output:  [real(r8) (:,:)   ]  
         kvw1d_t2             =>    temperature_inst%kvw1d_out2                , & ! Output:  [real(r8) (:,:)   ]  
         
         kwr1d_u             =>    urbanparams_inst%kwr1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         kwr1d_t1             =>    temperature_inst%kwr1d_out1                , & ! Output:  [real(r8) (:,:)   ]  
         kwr1d_t2             =>    temperature_inst%kwr1d_out2                , & ! Output:  [real(r8) (:,:)   ]  

         krw1d_u             =>    urbanparams_inst%krw1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         krw1d_t1             =>    temperature_inst%krw1d_out1                , & ! Output:  [real(r8) (:,:)   ]  
         krw1d_t2             =>    temperature_inst%krw1d_out2                , & ! Output:  [real(r8) (:,:)   ]  

         kvr1d_u             =>    urbanparams_inst%kvr1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         kvr1d_t1             =>    temperature_inst%kvr1d_out1                , & ! Output:  [real(r8) (:,:)   ]  
         kvr1d_t2             =>    temperature_inst%kvr1d_out2                , & ! Output:  [real(r8) (:,:)   ]  

         krv1d_u             =>    urbanparams_inst%krv1d_out                , & ! Input:  [real(r8) (:,:,:)   ]  
         krv1d_t1             =>    temperature_inst%krv1d_out1                , & ! Output:  [real(r8) (:,:)   ]  
         krv1d_t2             =>    temperature_inst%krv1d_out2               , & ! Output:  [real(r8) (:,:)   ]  
 
         kwg1d_u             =>    urbanparams_inst%kwg1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         kwg1d_t             =>    temperature_inst%kwg1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         kgw1d_u             =>    urbanparams_inst%kgw1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         kgw1d_t             =>    temperature_inst%kgw1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         kgv1d_u             =>    urbanparams_inst%kgv1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         kgv1d_t             =>    temperature_inst%kgv1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         ksw1d_u             =>    urbanparams_inst%ksw1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         ksw1d_t             =>    temperature_inst%ksw1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         kvg1d_u             =>    urbanparams_inst%kvg1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         kvg1d_t             =>    temperature_inst%kvg1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         ksr1d_u             =>    urbanparams_inst%ksr1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         ksr1d_t             =>    temperature_inst%ksr1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         ksv1d_u             =>    urbanparams_inst%ksv1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         ksv1d_t             =>    temperature_inst%ksv1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         kws1d_u             =>    urbanparams_inst%kws1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         kws1d_t             =>    temperature_inst%kws1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         kvs1d_u             =>    urbanparams_inst%kvs1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         kvs1d_t             =>    temperature_inst%kvs1d_out                , & ! Output:  [real(r8) (:,:)   ]  
         krs1d_u             =>    urbanparams_inst%krs1d_out                , & ! Input:  [real(r8) (:,:)   ]  
         krs1d_t             =>    temperature_inst%krs1d_out                , & ! Output:  [real(r8) (:,:)   ]  

         ksg1d_u             =>    urbanparams_inst%ksg1d_out                , & ! Input:  [real(r8) (:)   ]  
         ksg1d_t             =>    temperature_inst%ksg1d_out                , & ! Output:  [real(r8) (:)   ]  
         kts1d_u             =>    urbanparams_inst%kts1d_out                , & ! Input:  [real(r8) (:)   ]  
         kts1d_t             =>    temperature_inst%kts1d_out                , & ! Output:  [real(r8) (:)   ]    
         
         vfww_k_u                  =>    urbanparams_inst%vfww_k_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfww_k_t1                 =>    temperature_inst%vfww_k_out1               , & ! Output: [real(r8) (:,:)     ]
         vfww_k_t2                 =>    temperature_inst%vfww_k_out2               , & ! Output: [real(r8) (:,:)     ]
         vfvv_k_u                  =>    urbanparams_inst%vfvv_k_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfvv_k_t1                 =>    temperature_inst%vfvv_k_out1               , & ! Output: [real(r8) (:,:)     ]
         vfvv_k_t2                 =>    temperature_inst%vfvv_k_out2               , & ! Output: [real(r8) (:,:)     ]
         vfwv_k_u                  =>    urbanparams_inst%vfwv_k_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfwv_k_t1                 =>    temperature_inst%vfwv_k_out1               , & ! Output: [real(r8) (:,:)     ]
         vfwv_k_t2                 =>    temperature_inst%vfwv_k_out2               , & ! Output: [real(r8) (:,:)     ]
         vfvw_k_u                  =>    urbanparams_inst%vfvw_k_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfvw_k_t1                 =>    temperature_inst%vfvw_k_out1               , & ! Output: [real(r8) (:,:)     ]
         vfvw_k_t2                 =>    temperature_inst%vfvw_k_out2               , & ! Output: [real(r8) (:,:)     ]
         vfwr_k_u                  =>    urbanparams_inst%vfwr_k_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfwr_k_t1                 =>    temperature_inst%vfwr_k_out1               , & ! Output: [real(r8) (:,:)     ]
         vfwr_k_t2                 =>    temperature_inst%vfwr_k_out2               , & ! Output: [real(r8) (:,:)     ]
         vfrw_k_u                  =>    urbanparams_inst%vfrw_k_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfrw_k_t1                 =>    temperature_inst%vfrw_k_out1               , & ! Output: [real(r8) (:,:)     ]
         vfrw_k_t2                 =>    temperature_inst%vfrw_k_out2               , & ! Output: [real(r8) (:,:)     ]
         vfvr_k_u                  =>    urbanparams_inst%vfvr_k_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfvr_k_t1                 =>    temperature_inst%vfvr_k_out1               , & ! Output: [real(r8) (:,:)     ]
         vfvr_k_t2                 =>    temperature_inst%vfvr_k_out2               , & ! Output: [real(r8) (:,:)     ]
         vfrv_k_u                  =>    urbanparams_inst%vfrv_k_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfrv_k_t1                 =>    temperature_inst%vfrv_k_out1               , & ! Output: [real(r8) (:,:)     ]
         vfrv_k_t2                 =>    temperature_inst%vfrv_k_out2               , & ! Output: [real(r8) (:,:)     ]
         vfww_f_u                  =>    urbanparams_inst%vfww_f_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfww_f_t1                 =>    temperature_inst%vfww_f_out1               , & ! Output: [real(r8) (:,:)     ]
         vfww_f_t2                 =>    temperature_inst%vfww_f_out2               , & ! Output: [real(r8) (:,:)     ]
         vfvv_f_u                  =>    urbanparams_inst%vfvv_f_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfvv_f_t1                 =>    temperature_inst%vfvv_f_out1               , & ! Output: [real(r8) (:,:)     ]
         vfvv_f_t2                 =>    temperature_inst%vfvv_f_out2               , & ! Output: [real(r8) (:,:)     ]
         vfwv_f_u                  =>    urbanparams_inst%vfwv_f_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfwv_f_t1                 =>    temperature_inst%vfwv_f_out1               , & ! Output: [real(r8) (:,:)     ]
         vfwv_f_t2                 =>    temperature_inst%vfwv_f_out2               , & ! Output: [real(r8) (:,:)     ]
         vfvw_f_u                  =>    urbanparams_inst%vfvw_f_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfvw_f_t1                 =>    temperature_inst%vfvw_f_out1               , & ! Output: [real(r8) (:,:)     ]
         vfvw_f_t2                 =>    temperature_inst%vfvw_f_out2               , & ! Output: [real(r8) (:,:)     ]
         vfwr_f_u                  =>    urbanparams_inst%vfwr_f_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfwr_f_t1                 =>    temperature_inst%vfwr_f_out1               , & ! Output: [real(r8) (:,:)     ]
         vfwr_f_t2                 =>    temperature_inst%vfwr_f_out2               , & ! Output: [real(r8) (:,:)     ]
         vfrw_f_u                  =>    urbanparams_inst%vfrw_f_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfrw_f_t1                 =>    temperature_inst%vfrw_f_out1               , & ! Output: [real(r8) (:,:)     ]
         vfrw_f_t2                 =>    temperature_inst%vfrw_f_out2               , & ! Output: [real(r8) (:,:)     ]
         vfvr_f_u                  =>    urbanparams_inst%vfvr_f_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfvr_f_t1                 =>    temperature_inst%vfvr_f_out1               , & ! Output: [real(r8) (:,:)     ]
         vfvr_f_t2                 =>    temperature_inst%vfvr_f_out2               , & ! Output: [real(r8) (:,:)     ]
         vfrv_f_u                  =>    urbanparams_inst%vfrv_f_out                , & ! Input:  [real(r8) (:,:,:)   ]
         vfrv_f_t1                 =>    temperature_inst%vfrv_f_out1               , & ! Output: [real(r8) (:,:)     ]
         vfrv_f_t2                 =>    temperature_inst%vfrv_f_out2               , & ! Output: [real(r8) (:,:)     ]

         vfwt_k_u                  =>    urbanparams_inst%vfwt_k_out                , & ! Input:  [real(r8) (:,:)     ]
         vfwt_k_t                  =>    temperature_inst%vfwt_k_out                , & ! Output: [real(r8) (:,:)     ]
         vftw_k_u                  =>    urbanparams_inst%vftw_k_out                , & ! Input:  [real(r8) (:,:)     ]
         vftw_k_t                  =>    temperature_inst%vftw_k_out                , & ! Output: [real(r8) (:,:)     ]
         vftv_k_u                  =>    urbanparams_inst%vftv_k_out                , & ! Input:  [real(r8) (:,:)     ]
         vftv_k_t                  =>    temperature_inst%vftv_k_out                , & ! Output: [real(r8) (:,:)     ]
         vfvt_k_u                  =>    urbanparams_inst%vfvt_k_out                , & ! Input:  [real(r8) (:,:)     ]
         vfvt_k_t                  =>    temperature_inst%vfvt_k_out                , & ! Output: [real(r8) (:,:)     ]
         vfsw_k_u                  =>    urbanparams_inst%vfsw_k_out                , & ! Input:  [real(r8) (:,:)     ]
         vfsw_k_t                  =>    temperature_inst%vfsw_k_out                , & ! Output: [real(r8) (:,:)     ]
         vfsr_k_u                  =>    urbanparams_inst%vfsr_k_out                , & ! Input:  [real(r8) (:,:)     ]
         vfsr_k_t                  =>    temperature_inst%vfsr_k_out                , & ! Output: [real(r8) (:,:)     ]
         vfsv_k_u                  =>    urbanparams_inst%vfsv_k_out                , & ! Input:  [real(r8) (:,:)     ]
         vfsv_k_t                  =>    temperature_inst%vfsv_k_out                , & ! Output: [real(r8) (:,:)     ]
         svfw_k_u                  =>    urbanparams_inst%svfw_k_out                , & ! Input:  [real(r8) (:,:)     ]
         svfw_k_t                  =>    temperature_inst%svfw_k_out                , & ! Output: [real(r8) (:,:)     ]
         svfv_k_u                  =>    urbanparams_inst%svfv_k_out                , & ! Input:  [real(r8) (:,:)     ]
         svfv_k_t                  =>    temperature_inst%svfv_k_out                , & ! Output: [real(r8) (:,:)     ]

         vfst_k_u                  =>    urbanparams_inst%vfst_k_out                , & ! Input:  [real(r8) (:)       ]
         vfst_k_t                  =>    temperature_inst%vfst_k_out                , & ! Output: [real(r8) (:)       ]
         svft_k_u                  =>    urbanparams_inst%svft_k_out                , & ! Input:  [real(r8) (:)       ]
         svft_k_t                  =>    temperature_inst%svft_k_out                , & ! Output: [real(r8) (:)       ]

         svfr_k_u                  =>    urbanparams_inst%svfr_k_out                , & ! Input:  [real(r8) (:,:)     ]
         svfr_k_t                  =>    temperature_inst%svfr_k_out                , & ! Output: [real(r8) (:,:)     ]
  
         vfwt_f_u                  =>    urbanparams_inst%vfwt_f_out                , & ! Input:  [real(r8) (:,:)     ]
         vfwt_f_t                  =>    temperature_inst%vfwt_f_out                , & ! Output: [real(r8) (:,:)     ]
         vftw_f_u                  =>    urbanparams_inst%vftw_f_out                , & ! Input:  [real(r8) (:,:)     ]
         vftw_f_t                  =>    temperature_inst%vftw_f_out                , & ! Output: [real(r8) (:,:)     ]
         vftv_f_u                  =>    urbanparams_inst%vftv_f_out                , & ! Input:  [real(r8) (:,:)     ]
         vftv_f_t                  =>    temperature_inst%vftv_f_out                , & ! Output: [real(r8) (:,:)     ]
         vfvt_f_u                  =>    urbanparams_inst%vfvt_f_out                , & ! Input:  [real(r8) (:,:)     ]
         vfvt_f_t                  =>    temperature_inst%vfvt_f_out                , & ! Output: [real(r8) (:,:)     ]
         vfsw_f_u                  =>    urbanparams_inst%vfsw_f_out                , & ! Input:  [real(r8) (:,:)     ]
         vfsw_f_t                  =>    temperature_inst%vfsw_f_out                , & ! Output: [real(r8) (:,:)     ]
         vfsr_f_u                  =>    urbanparams_inst%vfsr_f_out                , & ! Input:  [real(r8) (:,:)     ]
         vfsr_f_t                  =>    temperature_inst%vfsr_f_out                , & ! Output: [real(r8) (:,:)     ]
         vfsv_f_u                  =>    urbanparams_inst%vfsv_f_out                , & ! Input:  [real(r8) (:,:)     ]
         vfsv_f_t                  =>    temperature_inst%vfsv_f_out                , & ! Output: [real(r8) (:,:)     ]
         svfw_f_u                  =>    urbanparams_inst%svfw_f_out                , & ! Input:  [real(r8) (:,:)     ]
         svfw_f_t                  =>    temperature_inst%svfw_f_out                , & ! Output: [real(r8) (:,:)     ]
         svfv_f_u                  =>    urbanparams_inst%svfv_f_out                , & ! Input:  [real(r8) (:,:)     ]
         svfv_f_t                  =>    temperature_inst%svfv_f_out                , & ! Output: [real(r8) (:,:)     ]
         svfr_f_u                  =>    urbanparams_inst%svfr_f_out                , & ! Input:  [real(r8) (:,:)     ]
         svfr_f_t                  =>    temperature_inst%svfr_f_out                , & ! Output: [real(r8) (:,:)     ]

         svft_f_u                  =>    urbanparams_inst%svft_f_out                , & ! Input:  [real(r8) (:)       ]
         svft_f_t                  =>    temperature_inst%svft_f_out                , & ! Output: [real(r8) (:)       ]
         vfst_f_u                  =>    urbanparams_inst%vfst_f_out                , & ! Input:  [real(r8) (:)       ]
         vfst_f_t                  =>    temperature_inst%vfst_f_out                , & ! Output: [real(r8) (:)       ]

         frontal_ai_u                  =>    urbanparams_inst%frontal_ai_out                , & ! Input:  [real(r8) (:)       ]
         frontal_ai_t                  =>    temperature_inst%frontal_ai_out                , & ! Output: [real(r8) (:)       ]
         plan_ai_u                  =>    urbanparams_inst%plan_ai_out             , & ! Input:  [real(r8) (:)       ]
         plan_ai_t                  =>    temperature_inst%plan_ai_out                 ,& ! Output: [real(r8) (:)       ]
         z_d_town_u                  =>    urbanparams_inst%z_d_town_out                , & ! Input:  [real(r8) (:)       ]
         z_d_town_t                  =>    temperature_inst%z_d_town_out                 ,& ! Output: [real(r8) (:)       ]
         z_0_town_u                  =>    urbanparams_inst%z_0_town_out                , & ! Input:  [real(r8) (:)       ]
         z_0_town_t                  =>    temperature_inst%z_0_town_out                 & ! Output: [real(r8) (:)       ]

!-------------------[kz.9]Ray tracing test------------------------- 
         )

      ! Define fields that appear on the restart file for non-urban landunits 
      do fl = 1,num_nourbanl
         l = filter_nourbanl(fl)

         sabs_roof_dir(l,:)      = spval
         sabs_roof_dif(l,:)      = spval
         sabs_sunwall_dir(l,:)   = spval
         sabs_sunwall_dif(l,:)   = spval
         sabs_shadewall_dir(l,:) = spval
         sabs_shadewall_dif(l,:) = spval
         sabs_improad_dir(l,:)   = spval
         sabs_improad_dif(l,:)   = spval
         sabs_perroad_dir(l,:)   = spval
         sabs_perroad_dif(l,:)   = spval    
         sabs_br_tree_dir(l,:)   = spval   
         sabs_br_tree_dif(l,:)   = spval 
         sabs_ar_tree_dir(l,:)   = spval   
         sabs_ar_tree_dif(l,:)   = spval   
         sabs_tree_dir(l,:)   = spval
         sabs_tree_dif(l,:)   = spval      
      end do

      ! Set input forcing fields
      do fl = 1,num_urbanl
         l = filter_urbanl(fl)
         g = lun%gridcell(l) 

         ! Need to set the following temperatures to some defined value even if it
         ! does not appear in the urban landunit for the net_longwave computation
         ! revise
         t_roof(l)      = 19._r8 + tfrz
         t_sunwall(l)   = 19._r8 + tfrz
         t_shadewall(l) = 19._r8 + tfrz
         t_improad(l)   = 19._r8 + tfrz
         t_perroad(l)   = 19._r8 + tfrz
         t_tree_perroad(l)   = 19._r8 + tfrz
         t_ar_tree(l)   = 19._r8 + tfrz
         t_br_tree(l)   = 19._r8 + tfrz
         
         ! Initial assignment of emissivity
         em_roof_s(l)    = em_roof(l)
         em_improad_s(l) = em_improad(l)
         em_perroad_s(l) = em_perroad(l)
         em_br_tree_s(l) = em_br_tree(l)
         em_ar_tree_s(l) = em_ar_tree(l)
!-------------------[kz.10]Ray tracing test------------------------- 
         fww1d_t1(l,:)=fww1d_u(l,:,1)
         fww1d_t2(l,:)=fww1d_u(l,:,2)
         fvv1d_t1(l,:)=fvv1d_u(l,:,1)
         fvv1d_t2(l,:)=fvv1d_u(l,:,2)
         fwv1d_t1(l,:)=fwv1d_u(l,:,1)
         fwv1d_t2(l,:)=fwv1d_u(l,:,2)
         fvw1d_t1(l,:)=fvw1d_u(l,:,1)
         fvw1d_t2(l,:)=fvw1d_u(l,:,2)
         fwr1d_t1(l,:)=fwr1d_u(l,:,1)
         fwr1d_t2(l,:)=fwr1d_u(l,:,2)
         frw1d_t1(l,:)=frw1d_u(l,:,1)
         frw1d_t2(l,:)=frw1d_u(l,:,2)
         fvr1d_t1(l,:)=fvr1d_u(l,:,1)
         fvr1d_t2(l,:)=fvr1d_u(l,:,2)
         frv1d_t1(l,:)=frv1d_u(l,:,1)    
         frv1d_t2(l,:)=frv1d_u(l,:,2)    
         
         fwg1d_t(l,:)=fwg1d_u(l,:)           
         fgw1d_t(l,:)=fgw1d_u(l,:)           
         fgv1d_t(l,:)=fgv1d_u(l,:)           
         fsw1d_t(l,:)=fsw1d_u(l,:)           
         fvg1d_t(l,:)=fvg1d_u(l,:)           
         fsr1d_t(l,:)=fsr1d_u(l,:)           
         fsv1d_t(l,:)=fsv1d_u(l,:)           
         fws1d_t(l,:)=fws1d_u(l,:)           
         fvs1d_t(l,:)=fvs1d_u(l,:)           
         frs1d_t(l,:)=frs1d_u(l,:)           
         fts1d_t(l)=fts1d_u(l)           
         fsg1d_t(l)=fsg1d_u(l)           

         kww1d_t1(l,:)=kww1d_u(l,:,1)
         kww1d_t2(l,:)=kww1d_u(l,:,2)
         kvv1d_t1(l,:)=kvv1d_u(l,:,1)
         kvv1d_t2(l,:)=kvv1d_u(l,:,2)
         kwv1d_t1(l,:)=kwv1d_u(l,:,1)
         kwv1d_t2(l,:)=kwv1d_u(l,:,2)
         kvw1d_t1(l,:)=kvw1d_u(l,:,1)
         kvw1d_t2(l,:)=kvw1d_u(l,:,2)
         kwr1d_t1(l,:)=kwr1d_u(l,:,1)
         kwr1d_t2(l,:)=kwr1d_u(l,:,2)
         krw1d_t1(l,:)=krw1d_u(l,:,1)
         krw1d_t2(l,:)=krw1d_u(l,:,2)
         kvr1d_t1(l,:)=kvr1d_u(l,:,1)
         kvr1d_t2(l,:)=kvr1d_u(l,:,2)
         krv1d_t1(l,:)=krv1d_u(l,:,1)    
         krv1d_t2(l,:)=krv1d_u(l,:,2)    
         
         kwg1d_t(l,:)=kwg1d_u(l,:)           
         kgw1d_t(l,:)=kgw1d_u(l,:)           
         kgv1d_t(l,:)=kgv1d_u(l,:)           
         ksw1d_t(l,:)=ksw1d_u(l,:)           
         kvg1d_t(l,:)=kvg1d_u(l,:)           
         ksr1d_t(l,:)=ksr1d_u(l,:)           
         ksv1d_t(l,:)=ksv1d_u(l,:)           
         kws1d_t(l,:)=kws1d_u(l,:)           
         kvs1d_t(l,:)=kvs1d_u(l,:)           
         krs1d_t(l,:)=krs1d_u(l,:)           
         kts1d_t(l)=kts1d_u(l)           
         ksg1d_t(l)=ksg1d_u(l)   
         
         vfww_k_t1(l,:) = vfww_k_u(l,:,1)
         vfww_k_t2(l,:) = vfww_k_u(l,:,2)
         vfvv_k_t1(l,:) = vfvv_k_u(l,:,1)
         vfvv_k_t2(l,:) = vfvv_k_u(l,:,2)
         vfwv_k_t1(l,:) = vfwv_k_u(l,:,1)
         vfwv_k_t2(l,:) = vfwv_k_u(l,:,2)
         vfvw_k_t1(l,:) = vfvw_k_u(l,:,1)
         vfvw_k_t2(l,:) = vfvw_k_u(l,:,2)
         vfwr_k_t1(l,:) = vfwr_k_u(l,:,1)
         vfwr_k_t2(l,:) = vfwr_k_u(l,:,2)
         vfrw_k_t1(l,:) = vfrw_k_u(l,:,1)
         vfrw_k_t2(l,:) = vfrw_k_u(l,:,2)
         vfvr_k_t1(l,:) = vfvr_k_u(l,:,1)
         vfvr_k_t2(l,:) = vfvr_k_u(l,:,2)
         vfrv_k_t1(l,:) = vfrv_k_u(l,:,1)
         vfrv_k_t2(l,:) = vfrv_k_u(l,:,2)
         vfww_f_t1(l,:) = vfww_f_u(l,:,1)
         vfww_f_t2(l,:) = vfww_f_u(l,:,2)
         vfvv_f_t1(l,:) = vfvv_f_u(l,:,1)
         vfvv_f_t2(l,:) = vfvv_f_u(l,:,2)
         vfwv_f_t1(l,:) = vfwv_f_u(l,:,1)
         vfwv_f_t2(l,:) = vfwv_f_u(l,:,2)
         vfvw_f_t1(l,:) = vfvw_f_u(l,:,1)
         vfvw_f_t2(l,:) = vfvw_f_u(l,:,2)
         vfwr_f_t1(l,:) = vfwr_f_u(l,:,1)
         vfwr_f_t2(l,:) = vfwr_f_u(l,:,2)
         vfrw_f_t1(l,:) = vfrw_f_u(l,:,1)
         vfrw_f_t2(l,:) = vfrw_f_u(l,:,2)
         vfvr_f_t1(l,:) = vfvr_f_u(l,:,1)
         vfvr_f_t2(l,:) = vfvr_f_u(l,:,2)
         vfrv_f_t1(l,:) = vfrv_f_u(l,:,1)
         vfrv_f_t2(l,:) = vfrv_f_u(l,:,2)

         vfwt_k_t(l,:) = vfwt_k_u(l,:)
         vftw_k_t(l,:) = vftw_k_u(l,:)
         vftv_k_t(l,:) = vftv_k_u(l,:)
         vfvt_k_t(l,:) = vfvt_k_u(l,:)
         vfsw_k_t(l,:) = vfsw_k_u(l,:)
         vfsr_k_t(l,:) = vfsr_k_u(l,:)
         vfsv_k_t(l,:) = vfsv_k_u(l,:)
         svfw_k_t(l,:) = svfw_k_u(l,:)
         svfv_k_t(l,:) = svfv_k_u(l,:)
         vfst_k_t(l)   = vfst_k_u(l)         
         svft_k_t(l)   = svft_k_u(l)
         svfr_k_t(l,:) = svfr_k_u(l,:)
         
         vfwt_f_t(l,:) = vfwt_f_u(l,:)
         vftw_f_t(l,:) = vftw_f_u(l,:)
         vftv_f_t(l,:) = vftv_f_u(l,:)
         vfvt_f_t(l,:) = vfvt_f_u(l,:)
         vfsw_f_t(l,:) = vfsw_f_u(l,:)
         vfsr_f_t(l,:) = vfsr_f_u(l,:)
         vfsv_f_t(l,:) = vfsv_f_u(l,:)
         svfw_f_t(l,:) = svfw_f_u(l,:)
         svfv_f_t(l,:) = svfv_f_u(l,:)
         svfr_f_t(l,:) = svfr_f_u(l,:)
         vfst_f_t(l)   = vfst_f_u(l)
         svft_f_t(l)   = svft_f_u(l)
         frontal_ai_t(l)   = frontal_ai_u(l)
         plan_ai_t(l)   = plan_ai_u(l)
         z_d_town_t(l)   = z_d_town_u(l)
         z_0_town_t(l)   = z_0_town_u(l)

        do fp = 1,num_urbanp
          p = filter_urbanp(fp)
          c = patch%column(p)
          l = patch%landunit(p)

          if (ctype(c) == icol_road_tree) then       
               t_br_tree(l)      = t_veg(p)
               t_ar_tree(l)      = t_veg(p)
          end if
        end do

!-------------------[kz.10]Ray tracing test------------------------- 
         ! Set urban temperatures and emissivity including snow effects.
         ! frac_sno=0 for walls because qflx_through_snow and qflx_through_liq are 0 for them
         ! frac_sno are non-zero for other urban surfaces, computed at
          ! clm_drv
          ! └─ HandleNewSnow                         (HydrologyNoDrainageMod)
          !   └─ UpdateQuantitiesForNewSnow         (SnowHydrologyMod:518)
          !       └─ scf_method%UpdateSnowDepthAndFrac   ← computes frac_sno here

         do c = coli(l),colf(l)
            if (ctype(c) == icol_roof)  then
               t_roof(l)      = t_grnd(c)
               em_roof_s(l) = em_roof(l)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
            else if (ctype(c) == icol_road_imperv) then 
               t_improad(l)   = t_grnd(c)
               em_improad_s(l) = em_improad(l)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
            else if (ctype(c) == icol_road_perv) then
               t_perroad(l)   = t_grnd(c)
               em_perroad_s(l) = em_perroad(l)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
            else if (ctype(c) == icol_road_tree) then
               t_tree_perroad(l) = t_grnd(c)
               em_br_tree_s(l) = em_br_tree(l)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)    
               em_ar_tree_s(l) = em_ar_tree(l)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)                                                         
            else if (ctype(c) == icol_sunwall) then
               t_sunwall(l)   = t_grnd(c)
            else if (ctype(c) == icol_shadewall) then
               t_shadewall(l) = t_grnd(c)
            end if
         end do
         lwdown(l) = forc_lwrad(g)
      end do

      ! Net longwave radiation for road and both walls in urban canyon allowing for multiple re-emission

      if (num_urbanl > 0) then
          call net_longwave (bounds,       &
               num_urbanl, filter_urbanl,  &
               canyon_hwr(begl:endl),      &
               ht_roof(begl:endl),      &
               A_v1(begl:endl), &
               A_v2(begl:endl), &               
               wtlunit_roof(begl:endl),      &
               wtroad_perv(begl:endl),     &
               wtroad_tree(begl:endl),     &
               tree_elaisai_per_uroad(begl:endl),&
               lwdown(begl:endl),          &
               em_roof_s(begl:endl),       &
               em_improad_s(begl:endl),    &
               em_perroad_s(begl:endl),    &
               em_wall(begl:endl),         &
               em_br_tree_s(begl:endl),         &
               em_ar_tree_s(begl:endl),         &
               t_roof(begl:endl),          &
               t_improad(begl:endl),       &
               t_perroad(begl:endl),       &
               t_sunwall(begl:endl),       &
               t_shadewall(begl:endl),     &
               t_tree_perroad(begl:endl),&
               t_br_tree(begl:endl),     &
               t_ar_tree(begl:endl),     &             
               urbanparams_inst, &
               solarabs_inst)          
      end if


      dtime = get_step_size_real()

      ! Determine variables needed for history output and communication with atm
      ! Loop over urban patches in clump

      do fp = 1,num_urbanp
         p = filter_urbanp(fp)
         c = patch%column(p)
         l = patch%landunit(p)
         g = patch%gridcell(p)

         ! Solar absorbed and longwave out and net
         ! per unit ground area (roof, road) and per unit wall area (sunwall, shadewall)
         ! Each urban pft has its own column - this is used in the logic below

         if (ctype(c) == icol_roof) then   
            eflx_lwrad_out(p) = lwup_roof(l)
            eflx_lwrad_net(p) = lwnet_roof(l)
            eflx_lwrad_net_u(p) = lwnet_roof(l)
            sabv(p)   = 0._r8
            sabg(p) = sabs_roof_dir(l,1)*forc_solad(g,1) + &
                 sabs_roof_dif(l,1)*forc_solai(g,1) + &
                 sabs_roof_dir(l,2)*forc_solad(g,2) + &
                 sabs_roof_dif(l,2)*forc_solai(g,2) 

         else if (ctype(c) == icol_sunwall) then   
            eflx_lwrad_out(p)   = lwup_sunwall(l)
            eflx_lwrad_net(p)   = lwnet_sunwall(l)
            eflx_lwrad_net_u(p) = lwnet_sunwall(l)
            sabv(p)   = 0._r8
            sabg(p) = sabs_sunwall_dir(l,1)*forc_solad(g,1) + &
                 sabs_sunwall_dif(l,1)*forc_solai(g,1) + &
                 sabs_sunwall_dir(l,2)*forc_solad(g,2) + &
                 sabs_sunwall_dif(l,2)*forc_solai(g,2) 

         else if (ctype(c) == icol_shadewall) then   
            eflx_lwrad_out(p)   = lwup_shadewall(l)
            eflx_lwrad_net(p)   = lwnet_shadewall(l)
            eflx_lwrad_net_u(p) = lwnet_shadewall(l)
            sabv(p)   = 0._r8
            sabg(p) = sabs_shadewall_dir(l,1)*forc_solad(g,1) + &
                 sabs_shadewall_dif(l,1)*forc_solai(g,1) + &
                 sabs_shadewall_dir(l,2)*forc_solad(g,2) + &
                 sabs_shadewall_dif(l,2)*forc_solai(g,2) 

         else if (ctype(c) == icol_road_perv) then       
            eflx_lwrad_out(p)   = lwup_perroad(l)
            eflx_lwrad_net(p)   = lwnet_perroad(l)
            eflx_lwrad_net_u(p) = lwnet_perroad(l)
            sabv(p)   = 0._r8
            sabg(p) = sabs_perroad_dir(l,1)*forc_solad(g,1) + &
                 sabs_perroad_dif(l,1)*forc_solai(g,1) + &
                 sabs_perroad_dir(l,2)*forc_solad(g,2) + &
                 sabs_perroad_dif(l,2)*forc_solai(g,2) 

         else if (ctype(c) == icol_road_tree) then       
            eflx_lwrad_out(p)   = lwup_tree(l) + lwup_tree_perroad(l)
            eflx_lwrad_net(p)   = lwnet_tree(l) + lwnet_tree_perroad(l)
            eflx_lwrad_net_u(p) = lwnet_tree(l) + lwnet_tree_perroad(l)
            ! A_g = ht_roof(l)/canyon_hwr(l)
            ! convert the original absorbed solar (flux per tree leave area)
            ! to flux per unit road tree area
            ! Ag = ht_roof(l)/canyon_hwr(l)
            if (tree_elaisai_per_uroad(l) > 1e-6_r8) then
              sabv(p) = (sabs_tree_dir(l,1)*forc_solad(g,1) + &
                  sabs_tree_dif(l,1)*forc_solai(g,1) + &
                  sabs_tree_dir(l,2)*forc_solad(g,2) + &
                  sabs_tree_dif(l,2)*forc_solai(g,2))*(A_v1(l)+A_v2(l))/(wtroad_tree(l)*ht_roof(l)/canyon_hwr(l))
            else
              sabv(p) = 0.0_r8
            end if
            ! pervious road absorbed solar flux -- for solar radiation
            ! we assumpe tree leaves are evenly distributed and shade all three road columns equivalently
            sabg(p) = sabs_perroad_dir(l,1)*forc_solad(g,1) + &
                 sabs_perroad_dif(l,1)*forc_solai(g,1) + &
                 sabs_perroad_dir(l,2)*forc_solad(g,2) + &
                 sabs_perroad_dif(l,2)*forc_solai(g,2) 

         else if (ctype(c) == icol_road_imperv) then       
            eflx_lwrad_out(p)   = lwup_improad(l)
            eflx_lwrad_net(p)   = lwnet_improad(l)
            eflx_lwrad_net_u(p) = lwnet_improad(l)
            sabv(p)   = 0._r8
            sabg(p) = sabs_improad_dir(l,1)*forc_solad(g,1) + &
                 sabs_improad_dif(l,1)*forc_solai(g,1) + &
                 sabs_improad_dir(l,2)*forc_solad(g,2) + &
                 sabs_improad_dif(l,2)*forc_solai(g,2) 
         end if

         
         fsa(p)    = sabv(p) + sabg(p)
         fsa_u(p)  = fsa(p)

      end do ! end loop over urban patches

    end associate

  end subroutine UrbanRadiation
  
  !-----------------------------------------------------------------------
    
  subroutine net_longwave (bounds                                                        , &
       num_urbanl, filter_urbanl, canyon_hwr,ht_roof,A_v1,A_v2,wtlunit_roof, wtroad_perv, wtroad_tree,tree_elaisai_per_uroad, &
       lwdown, em_roof, em_improad, em_perroad, em_wall,em_br_tree,em_ar_tree                 , &
       t_roof,  t_improad, t_perroad, t_sunwall, t_shadewall,t_tree_perroad,t_br_tree,t_ar_tree              , &
       urbanparams_inst,solarabs_inst)
    !
    ! !DESCRIPTION: 
    ! Net longwave radiation for road and both walls in urban canyon allowing for 
    ! multiple reflection. Also net longwave radiation for urban roof. 
    !
    ! !USES:
    use clm_varcon , only : sb
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds                  
    integer , intent(in)  :: num_urbanl                      ! number of urban landunits
    integer , intent(in)  :: filter_urbanl(:)                ! urban landunit filter
    real(r8), intent(in)  :: canyon_hwr( bounds%begl: )      ! ratio of building height to street width [landunit]
    real(r8), intent(in)  :: ht_roof( bounds%begl: )      ! building height[landunit]
    real(r8), intent(in)  :: wtlunit_roof( bounds%begl: )      ! weight of roof with respect to urban landunit (-)
    real(r8), intent(in)  :: wtroad_perv( bounds%begl: )     ! weight of pervious road wrt total road [landunit]
    real(r8), intent(in)  :: wtroad_tree( bounds%begl: )     ! weight of road tree wrt total road [landunit]
    real(r8), intent(in)  :: tree_elaisai_per_uroad( bounds%begl: ) ! exposed tree elai+esai per unit total road area [landunit]
    real(r8), intent(in)  :: A_v1( bounds%begl: )      ! Leaf area for urban tree canopy below roof  
    real(r8), intent(in)  :: A_v2( bounds%begl: )     ! Leaf area for urban tree canopy above roof  

    real(r8), intent(in)  :: lwdown( bounds%begl: )          ! atmospheric longwave radiation (W/m**2) [landunit]
    real(r8), intent(in)  :: em_roof( bounds%begl: )         ! roof emissivity [landunit]
    real(r8), intent(in)  :: em_improad( bounds%begl: )      ! impervious road emissivity [landunit]
    real(r8), intent(in)  :: em_perroad( bounds%begl: )      ! pervious road emissivity [landunit]
    real(r8), intent(in)  :: em_wall( bounds%begl: )         ! wall emissivity [landunit]
    real(r8), intent(in)  :: em_br_tree( bounds%begl: )         ! tree emissivity [landunit]
    real(r8), intent(in)  :: em_ar_tree( bounds%begl: )         ! tree emissivity [landunit]
    real(r8), intent(in)  :: t_roof( bounds%begl: )          ! roof temperature (K) [landunit]
    real(r8), intent(in)  :: t_improad( bounds%begl: )       ! impervious road temperature (K) [landunit]
    real(r8), intent(in)  :: t_perroad( bounds%begl: )       ! ervious road temperature (K) [landunit]
    real(r8), intent(in)  :: t_tree_perroad( bounds%begl: )       ! 
    real(r8), intent(in)  :: t_sunwall( bounds%begl: )       ! sunlit wall temperature (K) [landunit]
    real(r8), intent(in)  :: t_shadewall( bounds%begl: )     ! shaded wall temperature (K) [landunit]
    real(r8), intent(in)  :: t_br_tree( bounds%begl: )     ! below-roof tree temperature (K) [landunit]
    real(r8), intent(in)  :: t_ar_tree( bounds%begl: )     ! above-roof tree temperature (K) [landunit]

    type(urbanparams_type) , intent(in) :: urbanparams_inst
    type(solarabs_type)   , intent(inout) :: solarabs_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: lwdown_road(bounds%begl:bounds%endl)         ! atmospheric longwave radiation for total road (W/m**2)
    real(r8) :: lwdown_sunwall(bounds%begl:bounds%endl)      ! atmospheric longwave radiation (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: lwdown_shadewall(bounds%begl:bounds%endl)    ! atmospheric longwave radiation (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: lwdown_br_tree(bounds%begl:bounds%endl)    ! atmospheric longwave radiation (per unit below-roof tree area) for below-roof tree (W/m**2)
    real(r8) :: lwdown_ar_tree(bounds%begl:bounds%endl)    ! atmospheric longwave radiation (per unit above-roof tree area) for above-roof tree (W/m**2)
    real(r8) :: lwdown_roof(bounds%begl:bounds%endl)    ! atmospheric longwave radiation (per unit roof area) for roof (W/m**2)
    real(r8) :: lwtot(bounds%begl:bounds%endl)               ! incoming longwave radiation (W/m**2)

    real(r8) :: improad_a(bounds%begl:bounds%endl)           ! absorbed longwave for improad (W/m**2)
    real(r8) :: improad_r(bounds%begl:bounds%endl)           ! reflected longwave for improad (W/m**2)
    real(r8) :: improad_r_sky(bounds%begl:bounds%endl)       ! improad_r to sky (W/m**2)
    real(r8) :: improad_r_sunwall(bounds%begl:bounds%endl)   ! improad_r to sunlit wall (W/m**2)
    real(r8) :: improad_r_shadewall(bounds%begl:bounds%endl) ! improad_r to shaded wall (W/m**2)
    !real(r8) :: improad_e(bounds%begl:bounds%endl)           ! emitted longwave for improad (W/m**2)
    real(r8) :: improad_e_sky(bounds%begl:bounds%endl)       ! improad_e to sky (W/m**2)
    real(r8) :: improad_e_sunwall(bounds%begl:bounds%endl)   ! improad_e to sunlit wall (W/m**2)
    real(r8) :: improad_e_shadewall(bounds%begl:bounds%endl) ! improad_e to shaded wall (W/m**2)
    real(r8) :: improad_r_br_tree(bounds%begl:bounds%endl) ! improad_r to below-roof tree (W/m**2)
    real(r8) :: improad_r_ar_tree(bounds%begl:bounds%endl) ! improad_r to above-roof tree (W/m**2)
    real(r8) :: improad_e_br_tree(bounds%begl:bounds%endl) ! improad_e to below-roof tree (W/m**2)
    real(r8) :: improad_e_ar_tree(bounds%begl:bounds%endl) ! improad_e to above-roof tree (W/m**2)

    real(r8) :: perroad_a(bounds%begl:bounds%endl)           ! absorbed longwave for perroad (W/m**2)
    real(r8) :: perroad_r(bounds%begl:bounds%endl)           ! reflected longwave for perroad (W/m**2)
    real(r8) :: perroad_r_sky(bounds%begl:bounds%endl)       ! perroad_r to sky (W/m**2)
    real(r8) :: perroad_r_sunwall(bounds%begl:bounds%endl)   ! perroad_r to sunlit wall (W/m**2)
    real(r8) :: perroad_r_shadewall(bounds%begl:bounds%endl) ! perroad_r to shaded wall (W/m**2)
    real(r8) :: perroad_e(bounds%begl:bounds%endl)           ! emitted longwave for perroad (W/m**2)
    real(r8) :: perroad_e_sky(bounds%begl:bounds%endl)       ! perroad_e to sky (W/m**2)
    real(r8) :: perroad_e_sunwall(bounds%begl:bounds%endl)   ! perroad_e to sunlit wall (W/m**2)
    real(r8) :: perroad_e_shadewall(bounds%begl:bounds%endl) ! perroad_e to shaded wall (W/m**2)
    real(r8) :: perroad_r_br_tree(bounds%begl:bounds%endl) ! perroad_r to below-roof tree (W/m**2)
    real(r8) :: perroad_r_ar_tree(bounds%begl:bounds%endl) ! perroad_r to above-roof tree (W/m**2)
    real(r8) :: perroad_e_br_tree(bounds%begl:bounds%endl) ! perroad_e to below-roof tree (W/m**2)
    real(r8) :: perroad_e_ar_tree(bounds%begl:bounds%endl) ! perroad_e to above-roof tree (W/m**2)

    real(r8) :: tree_perroad_a(bounds%begl:bounds%endl)           ! absorbed longwave for tree_perroad (W/m**2)
    real(r8) :: tree_perroad_r(bounds%begl:bounds%endl)           ! reflected longwave for tree_perroad (W/m**2)
    real(r8) :: tree_perroad_r_sky(bounds%begl:bounds%endl)       ! tree_perroad_r to sky (W/m**2)
    real(r8) :: tree_perroad_r_sunwall(bounds%begl:bounds%endl)   ! tree_perroad_r to sunlit wall (W/m**2)
    real(r8) :: tree_perroad_r_shadewall(bounds%begl:bounds%endl) ! tree_perroad_r to shaded wall (W/m**2)
    real(r8) :: tree_perroad_e(bounds%begl:bounds%endl)           ! emitted longwave for tree_perroad (W/m**2)
    real(r8) :: tree_perroad_e_sky(bounds%begl:bounds%endl)       ! tree_perroad_e to sky (W/m**2)
    real(r8) :: tree_perroad_e_sunwall(bounds%begl:bounds%endl)   ! tree_perroad_e to sunlit wall (W/m**2)
    real(r8) :: tree_perroad_e_shadewall(bounds%begl:bounds%endl) ! tree_perroad_e to shaded wall (W/m**2)
    real(r8) :: tree_perroad_r_br_tree(bounds%begl:bounds%endl) ! tree_perroad_r to below-roof tree (W/m**2)
    real(r8) :: tree_perroad_r_ar_tree(bounds%begl:bounds%endl) ! tree_perroad_r to above-roof tree (W/m**2)
    real(r8) :: tree_perroad_e_br_tree(bounds%begl:bounds%endl) ! tree_perroad_e to below-roof tree (W/m**2)
    real(r8) :: tree_perroad_e_ar_tree(bounds%begl:bounds%endl) ! tree_perroad_e to above-roof tree (W/m**2)

    real(r8) :: road_a(bounds%begl:bounds%endl)              ! absorbed longwave for total road (W/m**2)
    real(r8) :: road_r(bounds%begl:bounds%endl)              ! reflected longwave for total road (W/m**2)
    real(r8) :: road_r_sky(bounds%begl:bounds%endl)          ! total road_r to sky (W/m**2)
    real(r8) :: road_r_sunwall(bounds%begl:bounds%endl)      ! total road_r to sunlit wall (W/m**2)
    real(r8) :: road_r_shadewall(bounds%begl:bounds%endl)    ! total road_r to shaded wall (W/m**2)
    real(r8) :: road_e(bounds%begl:bounds%endl)              ! emitted longwave for total road (W/m**2)
    real(r8) :: road_e_sky(bounds%begl:bounds%endl)          ! total road_e to sky (W/m**2)
    real(r8) :: road_e_sunwall(bounds%begl:bounds%endl)      ! total road_e to sunlit wall (W/m**2)
    real(r8) :: road_e_shadewall(bounds%begl:bounds%endl)    ! total road_e to shaded wall (W/m**2)
    real(r8) :: road_r_br_tree(bounds%begl:bounds%endl) ! road_r to below-roof tree (W/m**2)
    real(r8) :: road_r_ar_tree(bounds%begl:bounds%endl) ! road_r to above-roof tree (W/m**2)
    real(r8) :: road_e_br_tree(bounds%begl:bounds%endl) ! road_e to below-roof tree (W/m**2)
    real(r8) :: road_e_ar_tree(bounds%begl:bounds%endl) ! road_e to above-roof tree (W/m**2)

    real(r8) :: sunwall_a(bounds%begl:bounds%endl)           ! absorbed longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_r(bounds%begl:bounds%endl)           ! reflected longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_r_sky(bounds%begl:bounds%endl)       ! sunwall_r to sky (W/m**2)
    real(r8) :: sunwall_r_road(bounds%begl:bounds%endl)      ! sunwall_r to road (W/m**2)
    real(r8) :: sunwall_r_shadewall(bounds%begl:bounds%endl) ! sunwall_r to opposing (shaded) wall (W/m**2)
    real(r8) :: sunwall_e(bounds%begl:bounds%endl)           ! emitted longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_e_sky(bounds%begl:bounds%endl)       ! sunwall_e to sky (W/m**2)
    real(r8) :: sunwall_e_road(bounds%begl:bounds%endl)      ! sunwall_e to road (W/m**2)
    real(r8) :: sunwall_e_shadewall(bounds%begl:bounds%endl) ! sunwall_e to opposing (shaded) wall (W/m**2)
    real(r8) :: sunwall_r_br_tree(bounds%begl:bounds%endl) ! sunwall_r to below-roof tree (W/m**2)
    real(r8) :: sunwall_r_ar_tree(bounds%begl:bounds%endl) ! sunwall_r to above-roof tree (W/m**2)
    real(r8) :: sunwall_e_br_tree(bounds%begl:bounds%endl) ! sunwall_e to below-roof tree (W/m**2)
    real(r8) :: sunwall_e_ar_tree(bounds%begl:bounds%endl) ! sunwall_e to above-roof tree (W/m**2)

    real(r8) :: shadewall_a(bounds%begl:bounds%endl)         ! absorbed longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_r(bounds%begl:bounds%endl)         ! reflected longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_r_sky(bounds%begl:bounds%endl)     ! shadewall_r to sky (W/m**2)
    real(r8) :: shadewall_r_road(bounds%begl:bounds%endl)    ! shadewall_r to road (W/m**2)
    real(r8) :: shadewall_r_sunwall(bounds%begl:bounds%endl) ! shadewall_r to opposing (sunlit) wall (W/m**2)
    real(r8) :: shadewall_e(bounds%begl:bounds%endl)         ! emitted longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_e_sky(bounds%begl:bounds%endl)     ! shadewall_e to sky (W/m**2)
    real(r8) :: shadewall_e_road(bounds%begl:bounds%endl)    ! shadewall_e to road (W/m**2)
    real(r8) :: shadewall_e_sunwall(bounds%begl:bounds%endl) ! shadewall_e to opposing (sunlit) wall (W/m**2)
    real(r8) :: shadewall_r_br_tree(bounds%begl:bounds%endl) ! shadewall_r to below-roof tree (W/m**2)
    real(r8) :: shadewall_r_ar_tree(bounds%begl:bounds%endl) ! shadewall_r to above-roof tree (W/m**2)
    real(r8) :: shadewall_e_br_tree(bounds%begl:bounds%endl) ! shadewall_e to below-roof tree (W/m**2)
    real(r8) :: shadewall_e_ar_tree(bounds%begl:bounds%endl) ! shadewall_e to above-roof tree (W/m**2)

    real(r8) :: br_tree_a(bounds%begl:bounds%endl)         ! absorbed longwave (per unit wall area) for below-roof tree (W/m**2)
    real(r8) :: br_tree_r(bounds%begl:bounds%endl)         ! reflected longwave (per unit wall area) for below-roof tree (W/m**2)
    real(r8) :: br_tree_r_sky(bounds%begl:bounds%endl)     ! br_tree_r to sky (W/m**2)
    real(r8) :: br_tree_r_road(bounds%begl:bounds%endl)    ! br_tree_r to road (W/m**2)
    real(r8) :: br_tree_r_sunwall(bounds%begl:bounds%endl) ! br_tree_r to opposing (sunlit) wall (W/m**2)
    real(r8) :: br_tree_r_shadewall(bounds%begl:bounds%endl) ! br_tree_r to opposing (sunlit) wall (W/m**2)
    real(r8) :: br_tree_e(bounds%begl:bounds%endl)         ! emitted longwave (per unit wall area) for below-roof tree (W/m**2)
    real(r8) :: br_tree_e_sky(bounds%begl:bounds%endl)     ! br_tree_e to sky (W/m**2)
    real(r8) :: br_tree_e_road(bounds%begl:bounds%endl)    ! br_tree_e to road (W/m**2)
    real(r8) :: br_tree_e_sunwall(bounds%begl:bounds%endl) ! br_tree_e to opposing (sunlit) wall (W/m**2)
    real(r8) :: br_tree_e_shadewall(bounds%begl:bounds%endl) ! br_tree_e to opposing (sunlit) wall (W/m**2)
    real(r8) :: br_tree_r_br_tree(bounds%begl:bounds%endl) ! br_tree_r to below-roof tree (W/m**2)
    real(r8) :: br_tree_r_ar_tree(bounds%begl:bounds%endl) ! br_tree_r to above-roof tree (W/m**2)
    real(r8) :: br_tree_e_br_tree(bounds%begl:bounds%endl) ! br_tree_e to below-roof tree (W/m**2)
    real(r8) :: br_tree_e_ar_tree(bounds%begl:bounds%endl) ! br_tree_e to above-roof tree (W/m**2)  

    real(r8) :: ar_tree_a(bounds%begl:bounds%endl)         ! absorbed longwave (per unit wall area) for above-roof tree (W/m**2)
    real(r8) :: ar_tree_r(bounds%begl:bounds%endl)         ! reflected longwave (per unit wall area) for above-roof tree (W/m**2)
    real(r8) :: ar_tree_r_sky(bounds%begl:bounds%endl)     ! ar_tree_r to sky (W/m**2)
    real(r8) :: ar_tree_r_road(bounds%begl:bounds%endl)    ! ar_tree_r to road (W/m**2)
    real(r8) :: ar_tree_r_sunwall(bounds%begl:bounds%endl) ! ar_tree_r to opposing (sunlit) wall (W/m**2)
    real(r8) :: ar_tree_r_shadewall(bounds%begl:bounds%endl) ! ar_tree_r to opposing (sunlit) wall (W/m**2)
    real(r8) :: ar_tree_r_roof(bounds%begl:bounds%endl) ! ar_tree_r to opposing (sunlit) roof (W/m**2) 
    real(r8) :: ar_tree_e(bounds%begl:bounds%endl)         ! emitted longwave (per unit wall area) for above-roof tree (W/m**2)
    real(r8) :: ar_tree_e_sky(bounds%begl:bounds%endl)     ! ar_tree_e to sky (W/m**2)
    real(r8) :: ar_tree_e_road(bounds%begl:bounds%endl)    ! ar_tree_e to road (W/m**2)
    real(r8) :: ar_tree_e_sunwall(bounds%begl:bounds%endl) ! ar_tree_e to opposing (sunlit) wall (W/m**2)
    real(r8) :: ar_tree_e_shadewall(bounds%begl:bounds%endl) ! ar_tree_e to opposing (sunlit) wall (W/m**2)
    real(r8) :: ar_tree_e_roof(bounds%begl:bounds%endl) ! ar_tree_e to opposing (sunlit) wall (W/m**2)
    real(r8) :: ar_tree_r_br_tree(bounds%begl:bounds%endl) ! ar_tree_r to below-roof tree (W/m**2)
    real(r8) :: ar_tree_r_ar_tree(bounds%begl:bounds%endl) ! ar_tree_r to above-roof tree (W/m**2)
    real(r8) :: ar_tree_e_br_tree(bounds%begl:bounds%endl) ! ar_tree_e to below-roof tree (W/m**2)
    real(r8) :: ar_tree_e_ar_tree(bounds%begl:bounds%endl) ! ar_tree_e to above-roof tree (W/m**2)  
 
    real(r8) :: roof_a(bounds%begl:bounds%endl)         ! absorbed longwave (per unit roof area) for roof (W/m**2)
    real(r8) :: roof_r(bounds%begl:bounds%endl)         ! reflected longwave (per unit roof area) for roof (W/m**2)
    real(r8) :: roof_r_sky(bounds%begl:bounds%endl)     ! roof_r to sky (W/m**2)
    real(r8) :: roof_e(bounds%begl:bounds%endl)         ! emitted longwave (per unit roof area) for shaded roof (W/m**2)
    real(r8) :: roof_e_sky(bounds%begl:bounds%endl)     ! roof_e to sky (W/m**2)
    real(r8) :: roof_r_ar_tree(bounds%begl:bounds%endl) ! roof_r to above-roof tree (W/m**2)
    real(r8) :: roof_e_ar_tree(bounds%begl:bounds%endl) ! roof_e to above-roof tree (W/m**2)  

    real(r8)  :: A_s(bounds%begl:bounds%endl)                           ! Area of the street canyon (normalized)
    real(r8)  :: A_g(bounds%begl:bounds%endl)                           ! Aarea of the ground (normalized)
    real(r8)  :: A_r(bounds%begl:bounds%endl)                           ! Area of the roof
    real(r8)  :: A_w(bounds%begl:bounds%endl)                           ! Area of the wall (normalized)
    
    integer  :: l,fl,iter                    ! indices
    integer, parameter  :: n = 50            ! number of interations
    real(r8) :: crit                         ! convergence criterion (W/m**2)
    real(r8) :: err                          ! energy conservation error (W/m**2)
    real(r8) :: err_r                          ! energy conservation error (W/m**2)
    real(r8) :: err_e                          ! energy conservation error (W/m**2)
    real(r8) :: wtroad_imperv(bounds%begl:bounds%endl)       ! weight of impervious road wrt total road
    logical  :: debug_write = .false.                  ! true => write out many intermediate variables for debugging 
    !-----------------------------------------------------------------------
    
    ! Enforce expected array sizes

    SHR_ASSERT_ALL_FL((ubound(canyon_hwr)      == (/bounds%endl/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(wtroad_perv)     == (/bounds%endl/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(lwdown)          == (/bounds%endl/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(em_roof)         == (/bounds%endl/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(em_improad)      == (/bounds%endl/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(em_perroad)      == (/bounds%endl/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(em_wall)         == (/bounds%endl/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(t_roof)          == (/bounds%endl/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(t_improad)       == (/bounds%endl/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(t_perroad)       == (/bounds%endl/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(t_sunwall)       == (/bounds%endl/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(t_shadewall)     == (/bounds%endl/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(lwnet_roof)      == (/bounds%endl/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(lwnet_improad)   == (/bounds%endl/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(lwnet_perroad)   == (/bounds%endl/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(lwnet_sunwall)   == (/bounds%endl/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(lwnet_shadewall) == (/bounds%endl/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(lwnet_canyon)    == (/bounds%endl/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(lwup_roof)       == (/bounds%endl/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(lwup_improad)    == (/bounds%endl/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(lwup_perroad)    == (/bounds%endl/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(lwup_sunwall)    == (/bounds%endl/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(lwup_shadewall)  == (/bounds%endl/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(lwup_canyon)     == (/bounds%endl/)), sourcefile, __LINE__)

    associate(                             & 

      fgw1d =>    urbanparams_inst%fgw1d_out  , & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor from ground to  one wall[landunit, nzcanm]    
      fgv1d =>    urbanparams_inst%fgv1d_out  , & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor from ground to vegetation [landunit, nzcanm]  
      fts1d =>    urbanparams_inst%fts1d_out ,  & ! Input:  [real(r8) (:) ]  Monte carlo view factor from ground to sky[landunit, nzcanm]    

      fsg1d =>    urbanparams_inst%fsg1d_out   ,& ! Input:  [real(r8) (:) ]    Monte carlo view factor from sky to ground[landunit]               
      fsw1d =>    urbanparams_inst%fsw1d_out , & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor from sky to wall[landunit, nzcanm]   
      fsr1d =>    urbanparams_inst%fsr1d_out , & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor from sky to roof[landunit, nzcanm]                      
      fsv1d =>    urbanparams_inst%fsv1d_out , & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor from sky to vegetation[landunit, nzcanm]      
                      
      fwg1d =>    urbanparams_inst%fwg1d_out   ,& ! Input:  [real(r8) (:,:) ]      Monte carlo view factor of from one wall to ground[landunit]  
      fww1d =>    urbanparams_inst%fww1d_out  , & ! Input:  [real(r8) (:,:,:) ]  Monte carlo view factor from on wall to opposing wall [landunit]  
      fwv1d =>    urbanparams_inst%fwv1d_out  , & ! Input:  [real(r8) (:,:,:) ]  Monte carlo view factor from on wall to vegetation[landunit]  
      fws1d =>    urbanparams_inst%fws1d_out  , & ! Input:  [real(r8) (:,:) ]    Monte carlo view factor from on wall to sky[landunit]  

      fvg1d =>    urbanparams_inst%fvg1d_out  , & ! Input:  [real(r8) (:,:) ]      Monte carlo view factor of from vegetation to ground[landunit]  
      fvw1d =>    urbanparams_inst%fvw1d_out  , & ! Input:  [real(r8) (:,:,:) ]  Monte carlo view factor from vegetation to one wall [landunit]  
      fvv1d =>    urbanparams_inst%fvv1d_out ,  & ! Input:  [real(r8) (:,:,:) ]  Monte carlo view factor from vegetation to vegetation[landunit]  
      fvs1d =>    urbanparams_inst%fvs1d_out ,  & ! Input:  [real(r8) (:,:) ]    Monte carlo view factor from vegetation to sky[landunit]  
      fvr1d =>    urbanparams_inst%fvr1d_out  , & ! Input:  [real(r8) (:,:,:) ]    Monte carlo view factor from ovegetation to roof[landunit]  
      
      frs1d =>    urbanparams_inst%frs1d_out , & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor from roof to sky [landunit, nzcanm]                 
      frv1d =>    urbanparams_inst%frv1d_out , & ! Input:  [real(r8) (:,:,:) ]  Monte carlo view factor  from roof to vegetation [landunit, nzcanm]              
      lwnet_roof   => solarabs_inst%lwnet_roof_lun   , & ! Output: [real(r8) (:) ]  net longwave flux at shaded roof
      lwnet_improad     => solarabs_inst%lwnet_improad_lun     , & ! Output: [real(r8) (:) ]  net longwave flux at impervious road
      ! required for validation at London & SLC
      improad_e     => solarabs_inst%lw_emi_improad_lun     , & ! Output: [real(r8) (:) ]  emitted longwave flux at impervious road
      lwnet_perroad     => solarabs_inst%lwnet_perroad_lun     , & ! Output: [real(r8) (:) ]  net longwave flux at pervious road
      lwnet_tree_perroad     => solarabs_inst%lwnet_tree_perroad_lun     , & ! Output: [real(r8) (:) ]  net longwave flux at pervious road beneath the tree canopy

      lwnet_sunwall     => solarabs_inst%lwnet_sunwall_lun     , & ! Output: [real(r8) (:) ]  net longwave flux at sunlit wall
      lwnet_shadewall   => solarabs_inst%lwnet_shadewall_lun   , & ! Output: [real(r8) (:) ]  net longwave flux at shaded wall
      lwnet_br_tree     => solarabs_inst%lwnet_br_tree_lun     , & ! Output: [real(r8) (:) ]  net longwave flux at below-roof tree
      lwnet_ar_tree     => solarabs_inst%lwnet_ar_tree_lun     , & ! Output: [real(r8) (:) ]  net longwave flux at above-roof tree
      lwnet_tree     => solarabs_inst%lwnet_tree_lun     , & ! Output: [real(r8) (:) ]  net longwave flux at above-roof tree
      lwnet_canyon      => solarabs_inst%lwnet_canyon_lun      , & ! Output: [real(r8) (:) ]  net longwave flux at canyon center
      lwdown_road_out      => solarabs_inst%lwdown_road_lun      , & ! Output: [real(r8) (:) ] atmospheric longwave radiation for total road (W/m**2)
      lwdown_sunwall_out      => solarabs_inst%lwdown_sunwall_lun      , & ! Output: [real(r8) (:) ] atmospheric longwave  radiation for sunlit wall (W/m**2)
      lwdown_shadewall_out      => solarabs_inst%lwdown_shadewall_lun      , & ! Output: [real(r8) (:) ] atmospheric longwave  radiation for shaded wall(W/m**2)
      lwdown_roof_out      => solarabs_inst%lwdown_roof_lun      , & ! Output: [real(r8) (:) ] atmospheric longwave  radiation for roof (W/m**2)
      lwdown_br_tree_out      => solarabs_inst%lwdown_br_tree_lun      , & ! Output: [real(r8) (:) ] atmospheric longwave  radiation for below-roof tree (W/m**2)
      lwdown_ar_tree_out      => solarabs_inst%lwdown_ar_tree_lun      , & ! Output: [real(r8) (:) ] atmospheric longwave  radiation for above-roof tree (W/m**2)
      lwup_roof    => solarabs_inst%lwup_roof_lun    , & ! Output: [real(r8) (:) ]  upward longwave flux at shaded roof
      lwup_improad      => solarabs_inst%lwup_improad_lun      , & ! Output: [real(r8) (:) ]  upward longwave flux at impervious road
      lwup_perroad      => solarabs_inst%lwup_perroad_lun      , & ! Output: [real(r8) (:) ]  upward longwave flux at pervious road
      lwup_tree_perroad      => solarabs_inst%lwup_tree_perroad_lun      , & ! Output: [real(r8) (:) ]  upward longwave flux at pervious road beneath the tree canopy
      lwup_sunwall      => solarabs_inst%lwup_sunwall_lun      , & ! Output: [real(r8) (:) ]  upward longwave flux at sunlit wall
      lwup_shadewall    => solarabs_inst%lwup_shadewall_lun    , & ! Output: [real(r8) (:) ]  upward longwave flux at shaded wall
      lwup_br_tree      => solarabs_inst%lwup_br_tree_lun      , & ! Output: [real(r8) (:) ]  upward longwave flux at below-roof tree
      lwup_ar_tree      => solarabs_inst%lwup_ar_tree_lun      , & ! Output: [real(r8) (:) ]  upward longwave flux at above-roof tree
      lwup_tree      => solarabs_inst%lwup_tree_lun      , & ! Output: [real(r8) (:) ]  upward longwave flux at above-roof tree
      lwup_canyon       => solarabs_inst%lwup_canyon_lun       & ! Output: [real(r8) (:) ]  upward longwave flux at canyon center
         )
      
      debug_write=.false.       
     ! Calculate impervious road

     do fl = 1,num_urbanl 
        l = filter_urbanl(fl)
        wtroad_imperv(l) = 1._r8 - wtroad_perv(l) - wtroad_tree(l)
        A_w(l)=ht_roof(l)
        A_r(l)=ht_roof(l)/(canyon_hwr(l)*(1._r8-wtlunit_roof(l))/wtlunit_roof(l))
        A_g(l)=ht_roof(l)/canyon_hwr(l)
        A_s(l)=A_r(l)+A_g(l)        
     end do
     
     do fl = 1,num_urbanl
        l = filter_urbanl(fl)
        ! atmospheric longwave radiation incident on walls and road in urban canyon.
        ! check for conservation (need to convert wall fluxes to ground area).
        ! lwdown (from atmosphere) = lwdown_road + (lwdown_sunwall + lwdown_shadewall)*canyon_hwr

        lwdown_road(l)      = lwdown(l) * fsg1d(l) 
        lwdown_sunwall(l)   = lwdown(l) * fsw1d(l,1)
        lwdown_shadewall(l) = lwdown(l) * fsw1d(l,1) 
        lwdown_roof(l) = lwdown(l) * fsr1d(l,2) 
        lwdown_br_tree(l) = lwdown(l) * fsv1d(l,1) 
        lwdown_ar_tree(l) = lwdown(l) * fsv1d(l,2)

        lwdown_road_out(l)      = lwdown(l) * fsg1d(l) 
        lwdown_sunwall_out(l)   = lwdown(l) * fsw1d(l,1)
        lwdown_shadewall_out(l) = lwdown(l) * fsw1d(l,1) 
        lwdown_roof_out(l) = lwdown(l) * fsr1d(l,2) 
        lwdown_br_tree_out(l) = lwdown(l) * fsv1d(l,1) 
        lwdown_ar_tree_out(l) = lwdown(l) * fsv1d(l,2)
        
   

         ! write(6,*) '----------------incoming longwave radiation------------ '
         ! write(6,*) 'lwdown(l) = ', lwdown(l)
         ! write(6,*) 'fsg1d(l),fsw1d(l,1),fsr1d(l,2),fsv1d(l,1),fsv1d(l,2)', fsg1d(l),fsw1d(l,1),fsr1d(l,2),fsv1d(l,1),fsv1d(l,2)
         ! write(6,*) 'A_s(l),A_g(l), A_w(l), A_v1(l), A_v2(l), A_r(l) = ', A_s(l),A_g(l), A_w(l), A_v1(l), A_v2(l), A_r(l)
      
        
        err = lwdown(l) - (lwdown_road(l)*A_g(l)/A_s(l) + (lwdown_shadewall(l) + lwdown_sunwall(l))* A_w(l)/A_s(l) &
                 +lwdown_br_tree(l)*A_v1(l)/A_s(l)+ lwdown_ar_tree(l)*A_v2(l)/A_s(l) + lwdown_roof(l)*A_r(l)/A_s(l))
         
        if (abs(err) > 0.10_r8 ) then
           write(iulog,*) 'urban incident atmospheric longwave radiation balance error',err
           write(iulog,*) 'l          = ',l
           write(iulog,*) 'lwdown     = ',lwdown(l)
           write(iulog,*) 'fsg1d(l)      = ',fsg1d(l)
           write(iulog,*) 'fsw1d(l,1)      = ',fsw1d(l,1)
           write(iulog,*) 'canyon_hwr = ',canyon_hwr(l)
           write(iulog,*) 'clm model is stopping'
           call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
        endif
     end do

     do fl = 1,num_urbanl
        l = filter_urbanl(fl)
        
        ! initial absorption, reflection, and emission for road and both walls. 
        ! distribute reflected and emitted radiation to sky, road, and walls according 
        ! to appropriate view factor. radiation reflected to road and walls will
        ! undergo multiple reflections within the canyon.

        road_a(l)              = 0.0_r8
        road_r(l)              = 0.0_r8
        road_e(l)              = 0.0_r8
        improad_a(l)           =     em_improad(l)  * lwdown_road(l) 
        improad_r(l)           = (1._r8-em_improad(l)) * lwdown_road(l) 
        improad_r_sky(l)       = improad_r(l) * fts1d(l) ! unweighted view factor; ! this flux towards sky is in respect to road
        improad_r_sunwall(l)   = improad_r(l) * fgw1d(l,1) ! weighted view factor; ! this flux is in respect to wall area
        improad_r_shadewall(l) = improad_r(l) * fgw1d(l,1)
        improad_e(l)           = em_improad(l) * sb * (t_improad(l)**4) 
        improad_e_sky(l)       = improad_e(l) * fts1d(l)
        improad_e_sunwall(l)   = improad_e(l) * fgw1d(l,1)
        improad_e_shadewall(l) = improad_e(l) * fgw1d(l,1)
        
        improad_r_br_tree(l) = improad_r(l) * fgv1d(l,1)
        improad_r_ar_tree(l) = improad_r(l) * fgv1d(l,2)
        improad_e_br_tree(l) = improad_e(l) * fgv1d(l,1)
        improad_e_ar_tree(l) = improad_e(l) * fgv1d(l,2)
        
      !   err_e=improad_e(l) -improad_e_br_tree(l)*A_v1(l)/A_g(l) -improad_e_ar_tree(l)*A_v2(l)/A_g(l) -improad_e_sky(l)-improad_e_sunwall(l)*A_w(l)/A_g(l)-improad_e_shadewall(l)*A_w(l)/A_g(l)
      !   err_r=improad_r(l) -improad_r_br_tree(l)*A_v1(l)/A_g(l) -improad_r_ar_tree(l)*A_v2(l)/A_g(l) -improad_r_sky(l)-improad_r_sunwall(l)*A_w(l)/A_g(l)-improad_r_shadewall(l)*A_w(l)/A_g(l)
      !   if (debug_write) then
      !      write(6,*) '-----improad--------l=---', l
      !      write(6,*) 'err_e,err_r', err_e,err_r
      !   end if   

        road_a(l)              = road_a(l) + improad_a(l)*wtroad_imperv(l)
        road_r(l)              = road_r(l) + improad_r(l)*wtroad_imperv(l)
        road_e(l)              = road_e(l) + improad_e(l)*wtroad_imperv(l)

        perroad_a(l)           =     em_perroad(l)  * lwdown_road(l)
        perroad_r(l)           = (1._r8-em_perroad(l)) * lwdown_road(l)
        perroad_r_sky(l)       = perroad_r(l) * fts1d(l)! unweighted view factor; ! this flux towards sky is in respect to road
        perroad_r_sunwall(l)   = perroad_r(l) * fgw1d(l,1)
        perroad_r_shadewall(l) = perroad_r(l) * fgw1d(l,1)
        perroad_e(l)           = em_perroad(l) * sb * (t_perroad(l)**4) 
        perroad_e_sky(l)       = perroad_e(l) * fts1d(l)
        perroad_e_sunwall(l)   = perroad_e(l) * fgw1d(l,1)
        perroad_e_shadewall(l) = perroad_e(l) * fgw1d(l,1)
        
        perroad_r_br_tree(l) = perroad_r(l) * fgv1d(l,1)
        perroad_r_ar_tree(l) = perroad_r(l) * fgv1d(l,2)
        perroad_e_br_tree(l) = perroad_e(l) * fgv1d(l,1)
        perroad_e_ar_tree(l) = perroad_e(l) * fgv1d(l,2) 

        tree_perroad_a(l)           =     em_perroad(l)  * lwdown_road(l)
        tree_perroad_r(l)           = (1._r8-em_perroad(l)) * lwdown_road(l)
        tree_perroad_r_sky(l)       = tree_perroad_r(l) * fts1d(l)! unweighted view factor; ! this flux towards sky is in respect to road
        tree_perroad_r_sunwall(l)   = tree_perroad_r(l) * fgw1d(l,1)
        tree_perroad_r_shadewall(l) = tree_perroad_r(l) * fgw1d(l,1)
        tree_perroad_e(l)           = em_perroad(l) * sb * (t_tree_perroad(l)**4) 
        tree_perroad_e_sky(l)       = tree_perroad_e(l) * fts1d(l)
        tree_perroad_e_sunwall(l)   = tree_perroad_e(l) * fgw1d(l,1)
        tree_perroad_e_shadewall(l) = tree_perroad_e(l) * fgw1d(l,1)
        
        tree_perroad_r_br_tree(l) = tree_perroad_r(l) * fgv1d(l,1)
        tree_perroad_r_ar_tree(l) = tree_perroad_r(l) * fgv1d(l,2)
        tree_perroad_e_br_tree(l) = tree_perroad_e(l) * fgv1d(l,1)
        tree_perroad_e_ar_tree(l) = tree_perroad_e(l) * fgv1d(l,2) 

      !   err_e=perroad_e(l) -perroad_e_br_tree(l)*A_v1(l)/A_g(l) -perroad_e_ar_tree(l)*A_v2(l)/A_g(l) -perroad_e_sky(l)-perroad_e_sunwall(l)*A_w(l)/A_g(l)-perroad_e_shadewall(l)*A_w(l)/A_g(l)
      !   err_r=perroad_r(l) -perroad_r_br_tree(l)*A_v1(l)/A_g(l) -perroad_r_ar_tree(l)*A_v2(l)/A_g(l) -perroad_r_sky(l)-perroad_r_sunwall(l)*A_w(l)/A_g(l)-perroad_r_shadewall(l)*A_w(l)/A_g (l)
      !   if (debug_write) then
      !      write(6,*) '-----perroad--------l=---', l
      !      write(6,*) 'err_e,err_r', err_e,err_r
      !   end if   
        
        road_a(l)              = road_a(l) + perroad_a(l)*wtroad_perv(l)+ tree_perroad_a(l)*wtroad_tree(l)
        road_r(l)              = road_r(l) + perroad_r(l)*wtroad_perv(l)+tree_perroad_r(l)*wtroad_tree(l)
        road_e(l)              = road_e(l) + perroad_e(l)*wtroad_perv(l)+tree_perroad_e(l)*wtroad_tree(l)

        road_r_sky(l)          = road_r(l) * fts1d(l)! unweighted view factor; ! this flux towards sky is in respect to road
        road_r_sunwall(l)      = road_r(l) * fgw1d(l,1)
        road_r_shadewall(l)    = road_r(l) * fgw1d(l,1)
        road_e_sky(l)          = road_e(l) * fts1d(l)
        road_e_sunwall(l)      = road_e(l) * fgw1d(l,1)
        road_e_shadewall(l)    = road_e(l) * fgw1d(l,1)
        road_r_br_tree(l) = road_r(l) * fgv1d(l,1)
        road_r_ar_tree(l) = road_r(l) * fgv1d(l,2)
        road_e_br_tree(l) = road_e(l) * fgv1d(l,1)
        road_e_ar_tree(l) = road_e(l) * fgv1d(l,2)

      !   err_e=road_e(l) -road_e_br_tree(l)*A_v1(l)/A_g(l) -road_e_ar_tree(l)*A_v2(l)/A_g(l) -road_e_sky(l)-road_e_sunwall(l)*A_w(l)/A_g(l)-road_e_shadewall(l)*A_w(l)/A_g(l)
      !   err_r=road_r(l) -road_r_br_tree(l)*A_v1(l)/A_g(l) -road_r_ar_tree(l)*A_v2(l)/A_g(l) -road_r_sky(l)-road_r_sunwall(l)*A_w(l)/A_g(l)-road_r_shadewall(l)*A_w(l)/A_g(l)
      !   if (debug_write) then
      !      write(6,*) '-----road--------l=---', l
      !      write(6,*) 'err_e,err_r', err_e,err_r
      !   end if   

      !   if (debug_write) then
      !      write(6,*) '-----------------l=', l
      !      write(6,*) 'em_perroad(l),em_improad(l),em_wall(l),em_br_tree(l),em_ar_tree(l),em_roof(l) ', em_perroad(l),em_improad(l),em_wall(l),em_br_tree(l),em_ar_tree(l),em_roof(l)
      !      write(6,*) 'wtroad_perv(l),wtroad_imperv(l)', wtroad_perv(l),wtroad_imperv(l)
      !      write(6,*) 't_perroad(l),t_improad(l),t_sunwall(l),t_shadewall(l),t_br_tree(l),t_ar_tree(l),t_roof(l) ', t_perroad(l),t_improad(l),t_sunwall(l),t_shadewall(l),t_br_tree(l),t_ar_tree(l),t_roof(l)
      !      write(6,*) 'lwdown_road(l),lwdown_sunwall(l),lwdown_shadewall(l), lwdown_br_tree(l), lwdown_ar_tree(l),lwdown_roof(l)', lwdown_road(l),lwdown_sunwall(l),lwdown_shadewall(l),  lwdown_br_tree(l), lwdown_ar_tree(l),lwdown_roof(l)
      !      write(6,*) 'fts1d(l),fgw1d(l,1),fts1d(l), fgw1d(l,1),fgv1d(l,1),fgv1d(l,2)',fts1d(l),fgw1d(l,1),fts1d(l), fgw1d(l,1),fgv1d(l,1),fgv1d(l,2)
      !      write(6,*) 'fws1d(l,1),fwg1d(l,1),fww1d(l,1,1),fws1d(l,1),fwg1d(l,1),fww1d(l,1,1), fwv1d(l,1,1),fwv1d(l,1,2)',fws1d(l,1),fwg1d(l,1),fww1d(l,1,1),fws1d(l,1),fwg1d(l,1),fww1d(l,1,1), fwv1d(l,1,1),fwv1d(l,1,2)
      !   end if   

        sunwall_a(l)           = em_wall(l) * lwdown_sunwall(l)
        sunwall_r(l)           = (1._r8-em_wall(l)) * lwdown_sunwall(l)
        sunwall_e(l)           = em_wall(l) * sb * (t_sunwall(l)**4) 

        sunwall_r_sky(l)       = sunwall_r(l) * fws1d(l,1) ! this flux towards sky is in respect to wall
        sunwall_r_road(l)      = sunwall_r(l) * fwg1d(l,1)
        sunwall_r_shadewall(l) = sunwall_r(l) * fww1d(l,1,1)
        sunwall_e_sky(l)       = sunwall_e(l) * fws1d(l,1)
        sunwall_e_road(l)      = sunwall_e(l) * fwg1d(l,1)
        sunwall_e_shadewall(l) = sunwall_e(l) * fww1d(l,1,1)
        
        sunwall_r_br_tree(l) = sunwall_r(l) * fwv1d(l,1,1)
        sunwall_r_ar_tree(l) = sunwall_r(l) * fwv1d(l,1,2)
        sunwall_e_br_tree(l) = sunwall_e(l) * fwv1d(l,1,1)
        sunwall_e_ar_tree(l) = sunwall_e(l) * fwv1d(l,1,2)
        
      !   err_e=sunwall_e(l) -sunwall_e_br_tree(l)*A_v1(l)/A_w(l) -sunwall_e_ar_tree(l)*A_v2(l)/A_w(l) -sunwall_e_sky(l)-sunwall_e_road(l)*A_g(l)/A_w(l)-sunwall_e_shadewall(l)*A_w(l)/A_w(l)
      !   err_r=sunwall_r(l) -sunwall_r_br_tree(l)*A_v1(l)/A_w(l) -sunwall_r_ar_tree(l)*A_v2(l)/A_w(l) -sunwall_r_sky(l)-sunwall_r_road(l)*A_g(l)/A_w(l)-sunwall_r_shadewall(l)*A_w(l)/A_w(l)
      !   if (debug_write) then
      !      write(6,*) '------sunwall--------l=---', l
      !      write(6,*) 'err_e,err_r', err_e,err_r
      !   end if     
        
        shadewall_a(l)         = em_wall(l) * lwdown_shadewall(l)
        shadewall_r(l)         = (1._r8-em_wall(l)) * lwdown_shadewall(l)
        shadewall_r_sky(l)     = shadewall_r(l) * fws1d(l,1)! this flux towards sky is in respect to wall
        shadewall_r_road(l)    = shadewall_r(l) * fwg1d(l,1)
        shadewall_r_sunwall(l) = shadewall_r(l) * fww1d(l,1,1)
        shadewall_e(l)         = em_wall(l) * sb * (t_shadewall(l)**4) 
        shadewall_e_sky(l)     = shadewall_e(l) * fws1d(l,1)
        shadewall_e_road(l)    = shadewall_e(l) * fwg1d(l,1)
        shadewall_e_sunwall(l) = shadewall_e(l) * fww1d(l,1,1)
        
        shadewall_r_br_tree(l) = shadewall_r(l) * fwv1d(l,1,1)
        shadewall_r_ar_tree(l) = shadewall_r(l) * fwv1d(l,1,2)
        shadewall_e_br_tree(l) = shadewall_e(l) * fwv1d(l,1,1)
        shadewall_e_ar_tree(l) = shadewall_e(l) * fwv1d(l,1,2)

      !   err_e=shadewall_e(l) -shadewall_e_br_tree(l)*A_v1(l)/A_w(l) -shadewall_e_ar_tree(l)*A_v2(l)/A_w(l) -shadewall_e_sky(l)-shadewall_e_road(l)*A_g(l)/A_w(l)-shadewall_e_sunwall(l)*A_w(l)/A_w(l)
      !   err_r=shadewall_r(l) -shadewall_r_br_tree(l)*A_v1(l)/A_w(l) -shadewall_r_ar_tree(l)*A_v2(l)/A_w(l) -shadewall_r_sky(l)-shadewall_r_road(l)*A_g(l)/A_w(l)-shadewall_r_sunwall(l)*A_w(l)/A_w(l)
      !   if (debug_write) then
      !      write(6,*) '-----shadewall--------l=---', l
      !      write(6,*) 'err_e,err_r', err_e,err_r
      !   end if   
        
        br_tree_a(l)         = em_br_tree(l) * lwdown_br_tree(l)
        br_tree_r(l)         = (1._r8-em_br_tree(l)) * lwdown_br_tree(l)
        br_tree_e(l)         = em_br_tree(l) * sb * (t_br_tree(l)**4) 

        br_tree_r_sky(l)     = br_tree_r(l) * fvs1d(l,1)! this flux towards sky is in respect to veg
        br_tree_r_road(l)    = br_tree_r(l) * fvg1d(l,1)
        br_tree_r_sunwall(l) = br_tree_r(l) * fvw1d(l,1,1)
        br_tree_r_shadewall(l) = br_tree_r(l) * fvw1d(l,1,1)
        
        br_tree_e_sky(l)     = br_tree_e(l) * fvs1d(l,1)
        br_tree_e_road(l)    = br_tree_e(l) * fvg1d(l,1)
        br_tree_e_sunwall(l) = br_tree_e(l) * fvw1d(l,1,1)
        br_tree_e_shadewall(l) = br_tree_e(l) * fvw1d(l,1,1)

        br_tree_r_br_tree(l) = br_tree_r(l) * fvv1d(l,1,1)
        br_tree_r_ar_tree(l) = br_tree_r(l) * fvv1d(l,1,2)
        br_tree_e_br_tree(l) = br_tree_e(l) * fvv1d(l,1,1)
        br_tree_e_ar_tree(l) = br_tree_e(l) * fvv1d(l,1,2)

      !   err_e=br_tree_e(l) -br_tree_e_br_tree(l) -br_tree_e_ar_tree(l)*A_v2(l)/A_v1(l) -br_tree_e_sky(l)-br_tree_e_road(l)*A_g(l)/A_v1(l)-br_tree_e_sunwall(l)*A_w(l)/A_v1(l)-br_tree_e_shadewall(l)*A_w(l)/A_v1(l)
      !   err_r=br_tree_r(l) -br_tree_r_br_tree(l) -br_tree_r_ar_tree(l)*A_v2(l)/A_v1(l) -br_tree_r_sky(l)-br_tree_r_road(l)*A_g(l)/A_v1(l)-br_tree_r_sunwall(l)*A_w(l)/A_v1(l)-br_tree_r_shadewall(l)*A_w(l)/A_v1(l)
        
      !    if (debug_write) then
      !       write(6,*) '-----br_tree inputs------- l = ', l
      !       write(6,*) 'em_br_tree(l)      = ', em_br_tree(l)
      !       write(6,*) 'lwdown_br_tree(l)  = ', lwdown_br_tree(l)
      !       write(6,*) 't_br_tree(l)       = ', t_br_tree(l)
      !       write(6,*) 'sb                 = ', sb

      !       write(6,*) 'fvs1d(l,1)         = ', fvs1d(l,1)
      !       write(6,*) 'fvg1d(l,1)         = ', fvg1d(l,1)
      !       write(6,*) 'fvw1d(l,1,1)       = ', fvw1d(l,1,1)
      !       write(6,*) 'fvv1d(l,1,1)       = ', fvv1d(l,1,1)
      !       write(6,*) 'fvv1d(l,1,2)       = ', fvv1d(l,1,2)

      !       write(6,*) 'A_v1(l)            = ', A_v1(l)
      !       write(6,*) 'A_v2(l)            = ', A_v2(l)
      !       write(6,*) 'A_g(l)             = ', A_g(l)
      !       write(6,*) 'A_w(l)             = ', A_w(l)

      !       write(6,*) 'br_tree_a(l)       = ', br_tree_a(l)
      !       write(6,*) 'br_tree_r(l)       = ', br_tree_r(l)
      !       write(6,*) 'br_tree_e(l)       = ', br_tree_e(l)

      !       write(6,*) 'br_tree_r_sky(l)   = ', br_tree_r_sky(l)
      !       write(6,*) 'br_tree_r_road(l)  = ', br_tree_r_road(l)
      !       write(6,*) 'br_tree_r_sunwall(l)= ', br_tree_r_sunwall(l)
      !       write(6,*) 'br_tree_r_shadewall(l)= ', br_tree_r_shadewall(l)
      !       write(6,*) 'br_tree_r_br_tree(l)= ', br_tree_r_br_tree(l)
      !       write(6,*) 'br_tree_r_ar_tree(l)= ', br_tree_r_ar_tree(l)

      !       write(6,*) 'br_tree_e_sky(l)   = ', br_tree_e_sky(l)
      !       write(6,*) 'br_tree_e_road(l)  = ', br_tree_e_road(l)
      !       write(6,*) 'br_tree_e_sunwall(l)= ', br_tree_e_sunwall(l)
      !       write(6,*) 'br_tree_e_shadewall(l)= ', br_tree_e_shadewall(l)
      !       write(6,*) 'br_tree_e_br_tree(l)= ', br_tree_e_br_tree(l)
      !       write(6,*) 'br_tree_e_ar_tree(l)= ', br_tree_e_ar_tree(l)

      !       write(6,*) 'A_v2/A_v1          = ', A_v2(l)/A_v1(l)
      !       write(6,*) 'A_g/A_v1           = ', A_g(l)/A_v1(l)
      !       write(6,*) 'A_w/A_v1           = ', A_w(l)/A_v1(l)
      !    end if
      !   if (debug_write) then
      !      write(6,*) '-----br_tree-------l=---', l
      !      write(6,*) 'err_e,err_r', err_e,err_r
      !   end if    

        ar_tree_a(l)         = em_ar_tree(l) * lwdown_ar_tree(l)
        ar_tree_r(l)         = (1._r8-em_ar_tree(l)) * lwdown_ar_tree(l)
        ar_tree_e(l)         = em_ar_tree(l) * sb * (t_ar_tree(l)**4) 

        ar_tree_r_sky(l)     = ar_tree_r(l) * fvs1d(l,2)! this flux towards sky is in respect to veg
        ar_tree_r_road(l)    = ar_tree_r(l) * fvg1d(l,2)
        ar_tree_r_sunwall(l) = ar_tree_r(l) * fvw1d(l,2,1)
        ar_tree_r_shadewall(l) = ar_tree_r(l) * fvw1d(l,2,1)
        ar_tree_r_roof(l) = ar_tree_r(l) * fvr1d(l,2,2)
        ar_tree_r_br_tree(l) = ar_tree_r(l) * fvv1d(l,2,1)
        ar_tree_r_ar_tree(l) = ar_tree_r(l) * fvv1d(l,2,2)

        ar_tree_e_sky(l)     = ar_tree_e(l) * fvs1d(l,2)
        ar_tree_e_road(l)    = ar_tree_e(l) * fvg1d(l,2)
        ar_tree_e_sunwall(l) = ar_tree_e(l) * fvw1d(l,2,1)
        ar_tree_e_shadewall(l) = ar_tree_e(l) * fvw1d(l,2,1)
        ar_tree_e_roof(l) = ar_tree_e(l) * fvr1d(l,2,2)
        ar_tree_e_br_tree(l) = ar_tree_e(l) * fvv1d(l,2,1)
        ar_tree_e_ar_tree(l) = ar_tree_e(l) * fvv1d(l,2,2)

        !err_e=ar_tree_e(l) -ar_tree_e_ar_tree(l) -ar_tree_e_br_tree(l)*A_v1(l)/A_v2(l) -ar_tree_e_sky(l)-ar_tree_e_road(l)*A_g(l)/A_v2(l)-ar_tree_e_sunwall(l)*A_w(l)/A_v2(l)-ar_tree_e_shadewall(l)*A_w(l)/A_v2(l)-ar_tree_e_roof(l)*A_r(l)/A_v2(l)
        !err_r=ar_tree_r(l) -ar_tree_r_ar_tree(l) -ar_tree_r_br_tree(l)*A_v1(l)/A_v2(l) -ar_tree_r_sky(l)-ar_tree_r_road(l)*A_g(l)/A_v2(l)-ar_tree_r_sunwall(l)*A_w(l)/A_v2(l)-ar_tree_r_shadewall(l)*A_w(l)/A_v2(l)-ar_tree_r_roof(l)*A_r(l)/A_v2(l)
        !if (debug_write) then
        !   write(6,*) '-----ar_tree------l=---', l
        !   write(6,*) 'err_e,err_r', err_e,err_r
        !end if     

        roof_a(l)         = em_roof(l) * lwdown_roof(l) 
        roof_r(l)         = (1._r8-em_roof(l)) * lwdown_roof(l)
        roof_r_sky(l)     = roof_r(l) * frs1d(l,2) ! this flux towards sky is in respect to roof
        
        roof_e(l)         = em_roof(l) * sb * (t_roof(l)**4) 
        roof_e_sky(l)     = roof_e(l) * frs1d(l,2) ! this flux towards sky is in respect to roof
        
        roof_r_ar_tree(l) = roof_r(l) * frv1d(l,2,2)
        roof_e_ar_tree(l) = roof_e(l) * frv1d(l,2,2)  
        
      !   err_e=roof_e(l) -roof_e_ar_tree(l)*A_v2(l)/A_r(l) -roof_e_sky(l)
      !   err_r=roof_r(l) -roof_r_ar_tree(l)*A_v2(l)/A_r(l) -roof_r_sky(l)
      !   if (debug_write) then
      !      write(6,*) '-----roof------l=---', l
      !      write(6,*) 'err_e,err_r', err_e,err_r
      !   end if   
               
      !   if (debug_write) then
      !      write(6,*) '-----------------l=----------', l
      !      write(6,*) 'fvs1d(l,1),fvg1d(l,1),fvw1d(l,1,1),fvv1d(l,1,1),fvv1d(l,1,2)',fvs1d(l,1),fvg1d(l,1),fvw1d(l,1,1),fvv1d(l,1,1),fvv1d(l,1,2)
      !      write(6,*) 'frs1d(l,2),frv1d(l,2,2)',frs1d(l,2),frv1d(l,2,2)
      !   end if   
        ! initialize sum of net and upward longwave radiation for road and both walls

        lwnet_improad(l)   = improad_e(l)   - improad_a(l)
        lwnet_perroad(l)   = perroad_e(l)   - perroad_a(l)
        lwnet_tree_perroad(l)   = tree_perroad_e(l)   - tree_perroad_a(l)
        lwnet_sunwall(l)   = sunwall_e(l)   - sunwall_a(l)
        lwnet_shadewall(l) = shadewall_e(l) - shadewall_a(l)
        lwnet_br_tree(l) = br_tree_e(l) - br_tree_a(l)
        lwnet_ar_tree(l) = ar_tree_e(l) - ar_tree_a(l)
        lwnet_roof(l) = roof_e(l) - roof_a(l)
        
        lwup_improad(l)   = improad_r_sky(l)   + improad_e_sky(l)
        lwup_perroad(l)   = perroad_r_sky(l)   + perroad_e_sky(l)
        lwup_tree_perroad(l)   = tree_perroad_r_sky(l)   + tree_perroad_e_sky(l)
        lwup_sunwall(l)   = sunwall_r_sky(l)   + sunwall_e_sky(l)
        lwup_shadewall(l) = shadewall_r_sky(l) + shadewall_e_sky(l)
        lwup_br_tree(l) = br_tree_r_sky(l) + br_tree_e_sky(l)
        lwup_ar_tree(l) = ar_tree_r_sky(l) + ar_tree_e_sky(l)
        lwup_roof(l) = roof_r_sky(l) + roof_e_sky(l)   
            
      !   if (debug_write) then
      !      write(6,*) '--------- Net longwave (lwnet_*) initial absorption -------'
      !      write(6,*) 'improad_e(l), improad_a(l), lwnet_improad(l) = ', improad_e(l), improad_a(l), lwnet_improad(l)
      !      write(6,*) 'perroad_e(l), perroad_a(l), lwnet_perroad(l) = ', perroad_e(l), perroad_a(l), lwnet_perroad(l)
      !      write(6,*) 'sunwall_e(l), sunwall_a(l), lwnet_sunwall(l) = ', sunwall_e(l), sunwall_a(l), lwnet_sunwall(l)
      !      write(6,*) 'shadewall_e(l), shadewall_a(l), lwnet_shadewall(l) = ', shadewall_e(l), shadewall_a(l), lwnet_shadewall(l)
      !      write(6,*) 'br_tree_e(l), br_tree_a(l), lwnet_br_tree(l) = ', br_tree_e(l), br_tree_a(l), lwnet_br_tree(l)
      !      write(6,*) 'ar_tree_e(l), ar_tree_a(l), lwnet_ar_tree(l) = ', ar_tree_e(l), ar_tree_a(l), lwnet_ar_tree(l)
      !      write(6,*) 'roof_e(l), roof_a(l), lwnet_roof(l) = ', roof_e(l), roof_a(l), lwnet_roof(l)

      !      write(6,*) '--------- Upward longwave (lwup_*) initial absorption --------- '
      !      write(6,*) 'improad_r_sky(l), improad_e_sky(l), lwup_improad(l) = ', improad_r_sky(l), improad_e_sky(l), lwup_improad(l)
      !      write(6,*) 'perroad_r_sky(l), perroad_e_sky(l), lwup_perroad(l) = ', perroad_r_sky(l), perroad_e_sky(l), lwup_perroad(l)
      !      write(6,*) 'sunwall_r_sky(l), sunwall_e_sky(l), lwup_sunwall(l) = ', sunwall_r_sky(l), sunwall_e_sky(l), lwup_sunwall(l)
      !      write(6,*) 'shadewall_r_sky(l), shadewall_e_sky(l), lwup_shadewall(l) = ', shadewall_r_sky(l), shadewall_e_sky(l), lwup_shadewall(l)
      !      write(6,*) 'br_tree_r_sky(l), br_tree_e_sky(l), lwup_br_tree(l) = ', br_tree_r_sky(l), br_tree_e_sky(l), lwup_br_tree(l)
      !      write(6,*) 'ar_tree_r_sky(l), ar_tree_e_sky(l), lwup_ar_tree(l) = ', ar_tree_r_sky(l), ar_tree_e_sky(l), lwup_ar_tree(l)
      !      write(6,*) 'roof_r_sky(l), roof_e_sky(l), lwup_roof(l) = ', roof_r_sky(l), roof_e_sky(l), lwup_roof(l)
      !   end if
        
     end do

     ! now account for absorption and reflection within canyon of fluxes from road and walls 
     ! allowing for multiple reflections
     !
     ! (1) absorption and reflection. note: emission from road and walls absorbed by walls and roads
     !     only occurs in first iteration. zero out for later iterations.
     !
     !     road: fluxes from walls need to be projected to ground area
     !     wall: fluxes from road need to be projected to wall area
     !
     ! (2) add net longwave for ith reflection to total net longwave
     !
     ! (3) distribute reflected radiation to sky, road, and walls according to view factors
     !
     ! (4) add upward longwave radiation to sky from road and walls for ith reflection to total
     !
     ! (5) stop iteration when absorption for ith reflection is less than some nominal amount. 
     !     small convergence criteria is required to ensure radiation is conserved

     do fl = 1,num_urbanl
        l = filter_urbanl(fl)

        
        do iter = 1, n
        !do iter = 1, 1
           ! step (1)

           lwtot(l) =  (sunwall_r_road(l) + sunwall_e_road(l)  &
                + shadewall_r_road(l) + shadewall_e_road(l)    &
                + br_tree_r_road(l) + br_tree_e_road(l)        &
                + ar_tree_r_road(l) + ar_tree_e_road(l))
           !err=0.0_r8
           !err=err+lwtot(l)*A_w(l)/A_s(l)
         !  if (debug_write) then
         !     write(6,*) '-------------- Reflected Longwave Total Radiation Terms for l =------- ', l
         !     write(6,*) 'sunwall_r_road(l) = ', sunwall_r_road(l)
         !     write(6,*) 'sunwall_e_road(l) = ', sunwall_e_road(l)
         !     write(6,*) 'shadewall_r_road(l) = ', shadewall_r_road(l)
         !     write(6,*) 'shadewall_e_road(l) = ', shadewall_e_road(l)
         !     write(6,*) 'br_tree_r_road(l) = ', br_tree_r_road(l)
         !     write(6,*) 'br_tree_e_road(l) = ', br_tree_e_road(l)
         !     write(6,*) 'ar_tree_r_road(l) = ', ar_tree_r_road(l)
         !     write(6,*) 'ar_tree_e_road(l) = ', ar_tree_e_road(l)
         !     write(6,*) 'lwtot(l) = ', lwtot(l) 
         !   end if                
           road_a(l)    = 0.0_r8
           road_r(l)    = 0.0_r8
           improad_r(l) = (1._r8-em_improad(l)) * lwtot(l)
           improad_a(l) =     em_improad(l)  * lwtot(l)
           road_a(l)    = road_a(l) + improad_a(l)*wtroad_imperv(l)
           road_r(l)    = road_r(l) + improad_r(l)*wtroad_imperv(l)
           perroad_r(l) = (1._r8-em_perroad(l)) * lwtot(l) 
           perroad_a(l) =     em_perroad(l)  * lwtot(l) 
           road_a(l)    = road_a(l) + perroad_a(l)*wtroad_perv(l)
           road_r(l)    = road_r(l) + perroad_r(l)*wtroad_perv(l)

           tree_perroad_r(l) = (1._r8-em_perroad(l)) * lwtot(l) 
           tree_perroad_a(l) =     em_perroad(l)  * lwtot(l) 
           road_a(l)    = road_a(l) + tree_perroad_a(l)*wtroad_tree(l)
           road_r(l)    = road_r(l) + tree_perroad_r(l)*wtroad_tree(l)

           lwtot(l) = (road_r_sunwall(l) + road_e_sunwall(l)) &
                + (shadewall_r_sunwall(l) + shadewall_e_sunwall(l)) &
                + (br_tree_r_sunwall(l) + br_tree_e_sunwall(l)) &
                + (ar_tree_r_sunwall(l) + ar_tree_e_sunwall(l)) 
           sunwall_a(l) =     em_wall(l)  * lwtot(l)
           sunwall_r(l) = (1._r8-em_wall(l)) * lwtot(l)

           lwtot(l) = (road_r_shadewall(l) + road_e_shadewall(l)) &
                + (sunwall_r_shadewall(l) + sunwall_e_shadewall(l)) &
                + (br_tree_r_shadewall(l) + br_tree_e_shadewall(l)) &
                + (ar_tree_r_shadewall(l) + ar_tree_e_shadewall(l))
           shadewall_a(l) =     em_wall(l)  * lwtot(l)
           shadewall_r(l) = (1._r8-em_wall(l)) * lwtot(l)

           lwtot(l) = (sunwall_r_br_tree(l) + sunwall_e_br_tree(l)) &
                + (shadewall_r_br_tree(l) + shadewall_e_br_tree(l)) &
                + (br_tree_r_br_tree(l) + br_tree_e_br_tree(l)) &
                + (ar_tree_r_br_tree(l) + ar_tree_e_br_tree(l)) &
                + (road_r_br_tree(l) + road_e_br_tree(l))

           br_tree_a(l) =     em_br_tree(l)  * lwtot(l)
           br_tree_r(l) = (1._r8-em_br_tree(l)) * lwtot(l)
           
         !   if (.NOT. debug_write) then
              !write(6,*) 'to road: sunwall_r_br_tree(l),sunwall_e_br_tree(l),shadewall_r_br_tree(l)',sunwall_r_br_tree(l),sunwall_e_br_tree(l),shadewall_r_br_tree(l)
              !write(6,*) 'to road: shadewall_e_br_tree(l),br_tree_r_br_tree(l),br_tree_e_br_tree(l)',shadewall_e_br_tree(l),br_tree_r_br_tree(l),br_tree_e_br_tree(l)
              !write(6,*) 'to road: ar_tree_r_br_tree(l),ar_tree_e_br_tree(l),road_r_br_tree(l),road_e_br_tree(l)',ar_tree_r_br_tree(l),ar_tree_e_br_tree(l),road_r_br_tree(l),road_e_br_tree(l)
              !write(6,*) 'to road: lwtot(l),br_tree_a(l),br_tree_r(l)',lwtot(l),br_tree_a(l),br_tree_r(l)
         !   end if 

           lwtot(l) = (sunwall_r_ar_tree(l) + sunwall_e_ar_tree(l)) &
                + (shadewall_r_ar_tree(l) + shadewall_e_ar_tree(l)) &
                + (br_tree_r_ar_tree(l) + br_tree_e_ar_tree(l)) &
                + (ar_tree_r_ar_tree(l) + ar_tree_e_ar_tree(l)) &
                + (road_r_ar_tree(l) + road_e_ar_tree(l))&
                + (roof_r_ar_tree(l) + roof_e_ar_tree(l))

           ar_tree_a(l) =     em_ar_tree(l)  * lwtot(l)
           ar_tree_r(l) = (1._r8-em_ar_tree(l)) * lwtot(l)
           
           lwtot(l) = (ar_tree_r_roof(l) + ar_tree_e_roof(l))
           roof_a(l) =     em_roof(l)  * lwtot(l)
           roof_r(l) = (1._r8-em_roof(l)) * lwtot(l) 

           sunwall_e_road(l)      = 0._r8
           sunwall_e_br_tree(l)      = 0._r8
           sunwall_e_ar_tree(l)      = 0._r8
           sunwall_e_shadewall(l) = 0._r8
           
           shadewall_e_road(l)    = 0._r8
           shadewall_e_br_tree(l)    = 0._r8
           shadewall_e_ar_tree(l)    = 0._r8
           shadewall_e_sunwall(l)    = 0._r8
           
           road_e_sunwall(l)      = 0._r8
           road_e_shadewall(l)    = 0._r8
           road_e_br_tree(l)    = 0._r8
           road_e_ar_tree(l)    = 0._r8
           
           roof_e_ar_tree(l) = 0._r8

           br_tree_e_br_tree(l) = 0._r8
           br_tree_e_ar_tree(l) = 0._r8
           br_tree_e_sunwall(l) = 0._r8
           br_tree_e_shadewall(l) = 0._r8
           br_tree_e_road(l) = 0._r8
           
           ar_tree_e_br_tree(l) = 0._r8
           ar_tree_e_ar_tree(l) = 0._r8
           ar_tree_e_sunwall(l) = 0._r8
           ar_tree_e_shadewall(l) = 0._r8
           ar_tree_e_road(l) = 0._r8
           ar_tree_e_roof(l) = 0._r8
        
           ! step (2)

           lwnet_improad(l)   = lwnet_improad(l)   - improad_a(l)
           lwnet_perroad(l)   = lwnet_perroad(l)   - perroad_a(l)
           lwnet_tree_perroad(l)   = lwnet_tree_perroad(l)   - tree_perroad_a(l)
           lwnet_sunwall(l)   = lwnet_sunwall(l)   - sunwall_a(l)
           lwnet_shadewall(l) = lwnet_shadewall(l) - shadewall_a(l)
           lwnet_br_tree(l) = lwnet_br_tree(l) - br_tree_a(l)
           lwnet_ar_tree(l) = lwnet_ar_tree(l) - ar_tree_a(l)
           lwnet_roof(l) = lwnet_roof(l) - roof_a(l)
           
           ! step (3)

           improad_r_sky(l)       = improad_r(l) * fts1d(l)
           improad_r_sunwall(l)   = improad_r(l) * fgw1d(l,1)
           improad_r_shadewall(l) = improad_r(l) * fgw1d(l,1)
           improad_r_br_tree(l) = improad_r(l) * fgv1d(l,1)
           improad_r_ar_tree(l) = improad_r(l) * fgv1d(l,2)
  
         !   err_r=improad_r(l) -improad_r_br_tree(l)*A_v1(l)/A_g(l) -improad_r_ar_tree(l)*A_v2(l)/A_g(l) -improad_r_sky(l)-improad_r_sunwall(l)*A_w(l)/A_g(l)-improad_r_shadewall(l)*A_w(l)/A_g(l)
         !   if (debug_write) then
         !      write(6,*) '--------------improad---------------l=---', l
         !      write(6,*) 'err_e,err_r', err_e,err_r
         !   end if  
           perroad_r_sky(l)       = perroad_r(l) * fts1d(l)
           perroad_r_sunwall(l)   = perroad_r(l) * fgw1d(l,1)
           perroad_r_shadewall(l) = perroad_r(l) * fgw1d(l,1)
           perroad_r_br_tree(l) = perroad_r(l) * fgv1d(l,1)
           perroad_r_ar_tree(l) = perroad_r(l) * fgv1d(l,2)

           tree_perroad_r_sky(l)       = tree_perroad_r(l) * fts1d(l)
           tree_perroad_r_sunwall(l)   = tree_perroad_r(l) * fgw1d(l,1)
           tree_perroad_r_shadewall(l) = tree_perroad_r(l) * fgw1d(l,1)
           tree_perroad_r_br_tree(l) = tree_perroad_r(l) * fgv1d(l,1)
           tree_perroad_r_ar_tree(l) = tree_perroad_r(l) * fgv1d(l,2)

         !   err_r=perroad_r(l) -perroad_r_br_tree(l)*A_v1(l)/A_g(l) -perroad_r_ar_tree(l)*A_v2(l)/A_g(l) -perroad_r_sky(l)-perroad_r_sunwall(l)*A_w(l)/A_g(l)-perroad_r_shadewall(l)*A_w(l)/A_g (l)
         !   if (debug_write) then
         !      write(6,*) '--------------perroad---------------l=---', l
         !      write(6,*) 'err_e,err_r', err_e,err_r
         !   end if   

           road_r_sky(l)          = road_r(l) * fts1d(l)
           road_r_sunwall(l)      = road_r(l) * fgw1d(l,1)
           road_r_shadewall(l)    = road_r(l) * fgw1d(l,1)
           road_r_br_tree(l) = road_r(l) * fgv1d(l,1)
           road_r_ar_tree(l) = road_r(l) * fgv1d(l,2)
  
         !   err_r=road_r(l) -road_r_br_tree(l)*A_v1(l)/A_g(l) -road_r_ar_tree(l)*A_v2(l)/A_g(l) -road_r_sky(l)-road_r_sunwall(l)*A_w(l)/A_g(l)-road_r_shadewall(l)*A_w(l)/A_g(l)
         !   if (debug_write) then
         !      write(6,*) '--------------road---------------l=---', l
         !      write(6,*) 'err_e,err_r', err_e,err_r
         !   end if   

           sunwall_r_sky(l)       = sunwall_r(l) * fws1d(l,1)
           sunwall_r_road(l)      = sunwall_r(l) * fwg1d(l,1)
           sunwall_r_shadewall(l) = sunwall_r(l) * fww1d(l,1,1)
           sunwall_r_br_tree(l) = sunwall_r(l) * fwv1d(l,1,1)
           sunwall_r_ar_tree(l) = sunwall_r(l) * fwv1d(l,1,2)
  
         !   err_r=sunwall_r(l) -sunwall_r_br_tree(l)*A_v1(l)/A_w(l) -sunwall_r_ar_tree(l)*A_v2(l)/A_w(l) -sunwall_r_sky(l)-sunwall_r_road(l)*A_g(l)/A_w(l)-sunwall_r_shadewall(l)*A_w(l)/A_w(l)
         !   if (debug_write) then
         !      write(6,*) '--------------sunwall---------------l=---', l
         !      write(6,*) 'err_e,err_r', err_e,err_r
         !   end if     

           shadewall_r_sky(l)     = shadewall_r(l) * fws1d(l,1)
           shadewall_r_road(l)    = shadewall_r(l) * fwg1d(l,1)
           shadewall_r_sunwall(l) = shadewall_r(l) * fww1d(l,1,1)
           shadewall_r_br_tree(l) = shadewall_r(l) * fwv1d(l,1,1)
           shadewall_r_ar_tree(l) = shadewall_r(l) * fwv1d(l,1,2)

         !   err_r=shadewall_r(l) -shadewall_r_br_tree(l)*A_v1(l)/A_w(l) -shadewall_r_ar_tree(l)*A_v2(l)/A_w(l) -shadewall_r_sky(l)-shadewall_r_road(l)*A_g(l)/A_w(l)-shadewall_r_sunwall(l)*A_w(l)/A_w(l)
         !   if (debug_write) then
         !      write(6,*) '--------------shadewall---------------l=---', l
         !      write(6,*) 'err_e,err_r', err_e,err_r
         !   end if   

           br_tree_r_sky(l)     = br_tree_r(l) * fvs1d(l,1)
           br_tree_r_road(l)    = br_tree_r(l) * fvg1d(l,1)
           br_tree_r_sunwall(l) = br_tree_r(l) * fvw1d(l,1,1)
           br_tree_r_shadewall(l) = br_tree_r(l) * fvw1d(l,1,1)
           br_tree_r_br_tree(l) = br_tree_r(l) * fvv1d(l,1,1)
           br_tree_r_ar_tree(l) = br_tree_r(l) * fvv1d(l,1,2)

         !   err_r=br_tree_r(l) -br_tree_r_br_tree(l) -br_tree_r_ar_tree(l)*A_v2(l)/A_v1(l) -br_tree_r_sky(l)-br_tree_r_road(l)*A_g(l)/A_v1(l)-br_tree_r_sunwall(l)*A_w(l)/A_v1(l)-br_tree_r_shadewall(l)*A_w(l)/A_v1(l)
         !   if (debug_write) then
         !      write(6,*) '--------------br_tree--------------l=---', l
         !      write(6,*) 'err_e,err_r', err_e,err_r
         !   end if    

           ar_tree_r_sky(l)     = ar_tree_r(l) * fvs1d(l,2)
           ar_tree_r_road(l)    = ar_tree_r(l) * fvg1d(l,2)
           ar_tree_r_sunwall(l) = ar_tree_r(l) * fvw1d(l,2,1)
           ar_tree_r_shadewall(l) = ar_tree_r(l) * fvw1d(l,2,1)
           ar_tree_r_roof(l) = ar_tree_r(l) * fvr1d(l,2,2)
           ar_tree_r_br_tree(l) = ar_tree_r(l) * fvv1d(l,2,1)
           ar_tree_r_ar_tree(l) = ar_tree_r(l) * fvv1d(l,2,2)
           
           roof_r_sky(l)     = roof_r(l) * frs1d(l,2)
           roof_r_ar_tree(l) = roof_r(l) * frv1d(l,2,2)

         !   err_r=roof_r(l) -roof_r_ar_tree(l)*A_v2(l)/A_r(l) -roof_r_sky(l)
         !   if (debug_write) then
         !      write(6,*) '--------------roof---------------l=---', l
         !      write(6,*) 'err_e,err_r', err_e,err_r
         !   end if   

           ! step (4)

           lwup_improad(l)   = lwup_improad(l)   + improad_r_sky(l)
           lwup_perroad(l)   = lwup_perroad(l)   + perroad_r_sky(l)
           lwup_tree_perroad(l)   = lwup_tree_perroad(l)   + tree_perroad_r_sky(l)
           lwup_sunwall(l)   = lwup_sunwall(l)   + sunwall_r_sky(l)
           lwup_shadewall(l) = lwup_shadewall(l) + shadewall_r_sky(l)
           lwup_br_tree(l) = lwup_br_tree(l) + br_tree_r_sky(l) 
           lwup_ar_tree(l) = lwup_ar_tree(l) + ar_tree_r_sky(l) 
           lwup_roof(l) = lwup_roof(l) + roof_r_sky(l)

         !   if (debug_write) then
         !      write(6,*) '--- lwnet_* AFTER subtracting absorbed flux ---'
         !      write(6,*) 'lwnet_improad(l) = ', lwnet_improad(l), ', improad_a(l) = ', improad_a(l)
         !      write(6,*) 'lwnet_perroad(l) = ', lwnet_perroad(l), ', perroad_a(l) = ', perroad_a(l)
         !      write(6,*) 'lwnet_sunwall(l) = ', lwnet_sunwall(l), ', sunwall_a(l) = ', sunwall_a(l)
         !      write(6,*) 'lwnet_shadewall(l) = ', lwnet_shadewall(l), ', shadewall_a(l) = ', shadewall_a(l)
         !      write(6,*) 'lwnet_br_tree(l) = ', lwnet_br_tree(l), ', br_tree_a(l) = ', br_tree_a(l)
         !      write(6,*) 'lwnet_ar_tree(l) = ', lwnet_ar_tree(l), ', ar_tree_a(l) = ', ar_tree_a(l)
         !      write(6,*) 'lwnet_roof(l) = ', lwnet_roof(l), ', roof_a(l) = ', roof_a(l)

         !      write(6,*) '--- lwup_canyon(l) components ---'
         !      write(6,*) 'lwup_improad(l) = ', lwup_improad(l)
         !      write(6,*) 'lwup_perroad(l) = ', lwup_perroad(l)
         !      write(6,*) 'lwup_sunwall(l) = ', lwup_sunwall(l)
         !      write(6,*) 'lwup_shadewall(l) = ', lwup_shadewall(l)
         !      write(6,*) 'lwup_br_tree(l) = ', lwup_br_tree(l)
         !      write(6,*) 'lwup_ar_tree(l) = ', lwup_ar_tree(l)
         !      write(6,*) 'lwup_roof(l) = ', lwup_roof(l)
         !   end if
           
           ! step (5)
           
           crit = max(road_a(l), sunwall_a(l), shadewall_a(l), roof_a(l), br_tree_a(l), ar_tree_a(l))
           if (crit < .001_r8) exit
        end do
        if (iter >= n) then
           write (iulog,*) 'urban net longwave radiation error: no convergence'
           write (iulog,*) 'clm model is stopping'
           call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
        endif
        
        ! tree road is in road land unit; later in subgridAveMod.F90, it was aggregated using c2l_scale = 3
        ! this should be treated carefully
        if (tree_elaisai_per_uroad(l) > 1e-6_r8) then
           ! also convert the longwave radiation into flux per unit road_tree column area
            lwnet_tree(l)=(lwnet_br_tree(l)*A_v1(l)+ lwnet_ar_tree(l)*A_v2(l))/(wtroad_tree(l)*A_g(l))
            lwup_tree(l)=(lwup_br_tree(l)*A_v1(l)+ lwup_ar_tree(l)*A_v2(l))/(wtroad_tree(l)*A_g(l))
        else
            lwnet_tree(l)=0.0_r8
            lwup_tree(l)=0.0_r8
        end if

        ! lwup_tree_per_col_area = (lwup_br_tree*A_v1 + lwup_ar_tree*A_v2) / (wtroad_tree*A_g)????

        ! total net longwave radiation for canyon. project wall fluxes to horizontal surface

        lwnet_canyon(l) = 0.0_r8        
        lwnet_canyon(l)=(lwnet_improad(l)*wtroad_imperv(l)+lwnet_perroad(l)*wtroad_perv(l)+lwnet_tree_perroad(l)*wtroad_tree(l))*A_g(l)/A_s(l)&
                       +(lwnet_sunwall(l) + lwnet_shadewall(l))* A_w(l)/A_s(l) &
                        +lwnet_roof(l)*A_r(l)/A_s(l) + lwnet_br_tree(l)*A_v1(l)/A_s(l)+ lwnet_ar_tree(l)*A_v2(l)/A_s(l) 

        ! total emitted longwave for canyon. project wall fluxes to horizontal

        lwup_canyon(l) = 0.0_r8
        lwup_canyon(l) = lwup_canyon(l) + (lwup_improad(l)*wtroad_imperv(l)+ lwup_perroad(l)*wtroad_perv(l)+lwup_tree_perroad(l)*wtroad_tree(l))*A_g(l)/A_s(l)&
                            + (lwup_sunwall(l) + lwup_shadewall(l))* A_w(l)/A_s(l)&
                            +lwup_br_tree(l)*A_v1(l)/A_s(l) + lwup_ar_tree(l)*A_v2(l)/A_s(l) + lwup_roof(l)*A_r(l)/A_s(l)

      !   if (debug_write) then
      !      write(6,*) '--- lwnet_canyon(l) components ---'
      !      write(6,*) 'lwnet_improad(l) = ', lwnet_improad(l)
      !      write(6,*) 'lwnet_perroad(l) = ', lwnet_perroad(l)
      !      write(6,*) 'lwnet_sunwall(l) = ', lwnet_sunwall(l)
      !      write(6,*) 'lwnet_shadewall(l) = ', lwnet_shadewall(l)
      !      write(6,*) 'lwnet_roof(l) = ', lwnet_roof(l)
      !      write(6,*) 'lwnet_br_tree(l) = ', lwnet_br_tree(l)
      !      write(6,*) 'lwnet_ar_tree(l) = ', lwnet_ar_tree(l)
      !      write(6,*) 'lwup_canyon(l) = ', lwnet_ar_tree(l)

      !      write(6,*) '--- lwup_canyon(l) components ---'
      !      write(6,*) 'lwup_improad(l) = ', lwup_improad(l)
      !      write(6,*) 'lwup_perroad(l) = ', lwup_perroad(l)
      !      write(6,*) 'lwup_sunwall(l) = ', lwup_sunwall(l)
      !      write(6,*) 'lwup_shadewall(l) = ', lwup_shadewall(l)
      !      write(6,*) 'lwup_br_tree(l) = ', lwup_br_tree(l)
      !      write(6,*) 'lwup_ar_tree(l) = ', lwup_ar_tree(l)
      !      write(6,*) 'lwup_roof(l) = ', lwup_roof(l)
      !      write(6,*) 'lwup_canyon(l) = ', lwup_roof(l)
      !   end if

        ! conservation check. note: previous conservation check confirms partioning of incident
        ! atmospheric longwave radiation to road and walls is conserved as
        ! lwdown (from atmosphere) = lwdown_improad + lwdown_perroad + (lwdown_sunwall + lwdown_shadewall)*canyon_hwr
        err = lwnet_canyon(l) - (lwup_canyon(l) - lwdown(l))
        !err=0.01_r8
        !write(6,*) 'lwup_canyon(l), lwdown(l), lwnet_canyon(l) ',lwup_canyon(l), lwdown(l), lwnet_canyon(l)
        !write(6,*) 'err',err
        
      !   if (debug_write) then
      !      write(6,*) 'lwup_canyon(l), lwdown(l), lwnet_canyon(l),err ',lwup_canyon(l), lwdown(l), lwnet_canyon(l),err
      !   end if 

        if (abs(err) > .10_r8 ) then
          write (iulog,*) 'urban net longwave radiation balance error',err
          write (iulog,*) 'clm model is stopping'
          call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
        end if


     end do

   end associate

 end subroutine net_longwave

end module UrbanRadiationMod
