module UrbanAlbedoMod

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
  use clm_varcon        , only : isecspday, degpsec
  use clm_varctl        , only : iulog
  use abortutils        , only : endrun  
  use UrbanParamsType   , only : urbanparams_type
  use WaterStateBulkType    , only : waterstatebulk_type
  use WaterDiagnosticBulkType    , only : waterdiagnosticbulk_type
  use SolarAbsorbedType , only : solarabs_type
  use SurfaceAlbedoType , only : surfalb_type
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  use PatchType         , only : patch  
  use GridcellType      , only : grc    
  use atm2lndType       , only : atm2lnd_type          
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanAlbedo       ! Urban physics - albedos  
  !
  ! PRIVATE MEMBER FUNCTIONS
  private :: SnowAlbedo       ! Snow albedos
  private :: incident_direct  ! Direct beam solar rad incident on walls and road in urban canyon
  private :: incident_diffuse ! Diffuse solar rad incident on walls and road in urban canyon 
  private :: gaussian_quadrature   ! Direct beam solar rad incident on walls and road in urban canyon 
  private :: net_solar        ! Solar radiation absorbed by road and both walls in urban canyon 
  
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine UrbanAlbedo (bounds, num_urbanl, filter_urbanl, &
       num_urbanc, filter_urbanc, num_urbanp, filter_urbanp, &
       waterstatebulk_inst, waterdiagnosticbulk_inst, urbanparams_inst, solarabs_inst, surfalb_inst,atm2lnd_inst) 
    !
    ! !DESCRIPTION: 
    ! Determine urban landunit component albedos
    !
    ! Note that this is called with the "inactive_and_active" version of the filters, because
    ! the variables computed here are needed over inactive points that might later become
    ! active (due to landuse change). Thus, this routine cannot depend on variables that are
    ! only computed over active points.
    !
    ! !USES:
    use shr_orb_mod   , only : shr_orb_decl, shr_orb_cosz
    use clm_varcon    , only : sb
    use column_varcon , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon , only : icol_road_perv, icol_road_imperv, icol_road_tree
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds  
    integer                , intent(in)    :: num_urbanl       ! number of urban landunits in clump
    integer                , intent(in)    :: filter_urbanl(:) ! urban landunit filter
    integer                , intent(in)    :: num_urbanc       ! number of urban columns in clump
    integer                , intent(in)    :: filter_urbanc(:) ! urban column filter
    integer                , intent(in)    :: num_urbanp       ! number of urban patches in clump
    integer                , intent(in)    :: filter_urbanp(:) ! urban pft filter
    type(waterstatebulk_type)  , intent(in)    :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(in)    :: waterdiagnosticbulk_inst
    type(urbanparams_type) , intent(inout) :: urbanparams_inst
    type(solarabs_type)    , intent(inout) :: solarabs_inst
    type(surfalb_type)     , intent(inout) :: surfalb_inst
    type(atm2lnd_type), intent(in)        :: atm2lnd_inst
    
    !
    ! !LOCAL VARIABLES:
    integer  :: fl,fp,fc,g,l,p,c,ib                                  ! indices
    integer  :: ic                                                   ! 0=unit incoming direct; 1=unit incoming diffuse
    integer  :: num_solar                                            ! counter
    real(r8) :: coszen             (bounds%begl:bounds%endl)         ! cosine solar zenith angle for next time step (landunit)
    real(r8) :: zen                (bounds%begl:bounds%endl)         ! solar zenith angle (radians)
    real(r8) :: sdir               (bounds%begl:bounds%endl, numrad) ! direct beam solar radiation on horizontal surface
    real(r8) :: sdif               (bounds%begl:bounds%endl, numrad) ! diffuse solar radiation on horizontal surface

    real(r8) :: sdif_road          (bounds%begl:bounds%endl, numrad) ! direct beam solar radiation incident on road
    real(r8) :: sdif_sunwall       (bounds%begl:bounds%endl, numrad) ! direct beam solar radiation (per unit wall area) incident on sunlit wall per unit incident flux
    real(r8) :: sdif_shadewall     (bounds%begl:bounds%endl, numrad) ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux
    real(r8) :: sdif_roof     (bounds%begl:bounds%endl, numrad) ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux
    real(r8) :: sdif_ar_tree      (bounds%begl:bounds%endl, numrad) ! direct beam solar radiation absorbed by above-roof canopy before reaching to roof per unit direct flux
    real(r8) :: sdif_br_tree       (bounds%begl:bounds%endl, numrad) ! direct beam solar radiation absorbed by below-roof canopy before reaching to roof per unit direct flux

    real(r8) :: albsnd_roof        (bounds%begl:bounds%endl, numrad) ! snow albedo for roof (direct)
    real(r8) :: albsni_roof        (bounds%begl:bounds%endl, numrad) ! snow albedo for roof (diffuse)
    real(r8) :: albsnd_improad     (bounds%begl:bounds%endl, numrad) ! snow albedo for impervious road (direct)
    real(r8) :: albsni_improad     (bounds%begl:bounds%endl, numrad) ! snow albedo for impervious road (diffuse)
    real(r8) :: albsnd_perroad     (bounds%begl:bounds%endl, numrad) ! snow albedo for pervious road (direct)
    real(r8) :: albsni_perroad     (bounds%begl:bounds%endl, numrad) ! snow albedo for pervious road (diffuse)
    real(r8) :: albsnd_tree     (bounds%begl:bounds%endl, numrad) ! snow albedo for road tree (direct)
    real(r8) :: albsni_tree     (bounds%begl:bounds%endl, numrad) ! snow albedo for road tree (diffuse)
    real(r8) :: albsnd_br_tree        (bounds%begl:bounds%endl, numrad) ! snow albedo for below-roof tree (direct)
    real(r8) :: albsni_br_tree        (bounds%begl:bounds%endl, numrad) ! snow albedo for below-roof tree (diffuse)
    real(r8) :: albsnd_ar_tree        (bounds%begl:bounds%endl, numrad) ! snow albedo for above-roof tree (direct)
    real(r8) :: albsni_ar_tree        (bounds%begl:bounds%endl, numrad) ! snow albedo for above-roof tree (diffuse)

    real(r8) :: alb_roof_dir_s     (bounds%begl:bounds%endl, numrad) ! direct roof albedo with snow effects
    real(r8) :: alb_roof_dif_s     (bounds%begl:bounds%endl, numrad) ! diffuse roof albedo with snow effects
    real(r8) :: alb_improad_dir_s  (bounds%begl:bounds%endl, numrad) ! direct impervious road albedo with snow effects
    real(r8) :: alb_perroad_dir_s  (bounds%begl:bounds%endl, numrad) ! direct pervious road albedo with snow effects
    real(r8) :: alb_improad_dif_s  (bounds%begl:bounds%endl, numrad) ! diffuse impervious road albedo with snow effects
    real(r8) :: alb_perroad_dif_s  (bounds%begl:bounds%endl, numrad) ! diffuse pervious road albedo with snow effects
    real(r8) :: alb_br_tree_dir_s  (bounds%begl:bounds%endl, numrad) ! diffuse below-roof tree albedo with snow effects
    real(r8) :: alb_ar_tree_dir_s  (bounds%begl:bounds%endl, numrad) ! diffuse above-roof tree albedo with snow effects
    real(r8) :: alb_br_tree_dif_s  (bounds%begl:bounds%endl, numrad) ! diffuse below-roof tree albedo with snow effects
    real(r8) :: alb_ar_tree_dif_s  (bounds%begl:bounds%endl, numrad) ! diffuse above-roof tree albedo with snow effects

    real(r8) :: forc_solad             (bounds%begl:bounds%endl, numrad)         ! forced solar
    ! For now, the albedo of tree is hard coded; I need to read them from surface data 
    real(r8) :: alb_br_tree_dir     (bounds%begl:bounds%endl, numrad) ! direct br tree albedo with snow effects
    real(r8) :: alb_ar_tree_dir     (bounds%begl:bounds%endl, numrad) ! direct br tree albedo with snow effects
    real(r8) :: alb_br_tree_dif     (bounds%begl:bounds%endl, numrad) ! diffuse br tree albedowith snow effects  
    real(r8) :: alb_ar_tree_dif     (bounds%begl:bounds%endl, numrad) ! diffuse br tree albedowith snow effects
    
    integer              :: start_time, end_time, clock_rate    ! Timekeeping variables
    real(r8)            :: elapsed_time                        ! Elapsed time     
    !-----------------------------------------------------------------------

    associate(                                                        &
         ctype              => col%itype                            , & ! Input:  [integer (:)    ]  column type                                        
         coli               => lun%coli                             , & ! Input:  [integer (:)    ]  beginning column index for landunit                
         canyon_hwr         => lun%canyon_hwr                       , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width          
         ht_roof            => lun%ht_roof                       , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width          
         wtlunit_roof         => lun%wtlunit_roof                       , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width          

         wtroad_perv        => lun%wtroad_perv                      , & ! Input:  [real(r8) (:)   ]  weight of pervious road wrt total road    
         wtroad_tree        => lun%wtroad_tree                      , & ! Input:  [real(r8) (:)   ]  weight of pervious road wrt total road    

         frac_sno           => waterdiagnosticbulk_inst%frac_sno_col         , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)       
         h1            => lun%h1                      , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width          
         h2            => lun%h2                       , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width          
         A_v1            => lun%A_v1                       , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width          
         A_v2         => lun%A_v2                       , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width          
         
         alb_roof_dir       => urbanparams_inst%alb_roof_dir        , & ! Output: [real(r8) (:,:) ]  direct roof albedo                              
         alb_roof_dif       => urbanparams_inst%alb_roof_dif        , & ! Output: [real(r8) (:,:) ]  diffuse roof albedo                             
         alb_improad_dir    => urbanparams_inst%alb_improad_dir     , & ! Output: [real(r8) (:,:) ]  direct impervious road albedo                   
         alb_improad_dif    => urbanparams_inst%alb_improad_dif     , & ! Output: [real(r8) (:,:) ]  diffuse imprevious road albedo                  
         alb_perroad_dir    => urbanparams_inst%alb_perroad_dir     , & ! Output: [real(r8) (:,:) ]  direct pervious road albedo                     
         alb_perroad_dif    => urbanparams_inst%alb_perroad_dif     , & ! Output: [real(r8) (:,:) ]  diffuse pervious road albedo                    
         alb_wall_dir       => urbanparams_inst%alb_wall_dir        , & ! Output: [real(r8) (:,:) ]  direct wall albedo                              
         alb_wall_dif       => urbanparams_inst%alb_wall_dif        , & ! Output: [real(r8) (:,:) ]  diffuse wall albedo                             
        

         sabs_roof_dir      =>    solarabs_inst%sabs_roof_dir_lun      , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by roof per unit ground area per unit incident flux
         sabs_sunwall_dir   =>    solarabs_inst%sabs_sunwall_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_shadewall_dir =>    solarabs_inst%sabs_shadewall_dir_lun , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_improad_dir   =>    solarabs_inst%sabs_improad_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_perroad_dir   =>    solarabs_inst%sabs_perroad_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
         sabs_br_tree_dir      =>    solarabs_inst%sabs_br_tree_dir_lun      , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by below-roof treeetation per unit treeetation area per unit incident flux
         sabs_ar_tree_dir      =>    solarabs_inst%sabs_ar_tree_dir_lun      , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by above-roof treeetation per unit treeetation area per unit incident flux
         sabs_roof_dif      =>    solarabs_inst%sabs_roof_dif_lun      , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by roof per unit ground area per unit incident flux
         sabs_sunwall_dif   =>    solarabs_inst%sabs_sunwall_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_shadewall_dif =>    solarabs_inst%sabs_shadewall_dif_lun , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_improad_dif   =>    solarabs_inst%sabs_improad_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_perroad_dif   =>    solarabs_inst%sabs_perroad_dif_lun   ,& ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux
         sabs_br_tree_dif      =>    solarabs_inst%sabs_br_tree_dif_lun      , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by below-roof treeetation per unit treeetation area per unit incident flux
         sabs_ar_tree_dif      =>    solarabs_inst%sabs_ar_tree_dif_lun      , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by above-roof treeetation per unit treeetation area per unit incident flux
         sabs_tree_dir   => solarabs_inst%sabs_tree_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
         sabs_tree_dif   => solarabs_inst%sabs_tree_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux

         sref_roof_dir   => solarabs_inst%sref_roof_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar reflected by roof per unit ground area per unit incident flux
         sref_sunwall_dir     => solarabs_inst%sref_sunwall_dir_lun     , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by sunwall per unit wall area per unit incident flux
         sref_shadewall_dir   => solarabs_inst%sref_shadewall_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar reflected by shadewall per unit wall area per unit incident flux
         sref_improad_dir     => solarabs_inst%sref_improad_dir_lun     , & ! Output: [real(r8) (:,:) ]  direct  solar reflected by impervious road per unit ground area per unit incident flux
         sref_perroad_dir     => solarabs_inst%sref_perroad_dir_lun     , & ! Output: [real(r8) (:,:) ]  direct  solar reflected by pervious road per unit ground area per unit incident flux
         sref_br_tree_dir     => solarabs_inst%sref_br_tree_dir_lun     , & ! Output: [real(r8) (:,:) ]  direct  solar reflected by below-roof tree per unit ground area per unit incident flux
         sref_ar_tree_dir     => solarabs_inst%sref_ar_tree_dir_lun     , & ! Output: [real(r8) (:,:) ]  direct  solar reflected by above-roof tree per unit ground area per unit incident flux
         sref_roof_dif   => solarabs_inst%sref_roof_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by roof per unit ground area per unit incident flux
         sref_sunwall_dif     => solarabs_inst%sref_sunwall_dif_lun     , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by sunwall per unit wall area per unit incident flux
         sref_shadewall_dif   => solarabs_inst%sref_shadewall_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by shadewall per unit wall area per unit incident flux
         sref_improad_dif     => solarabs_inst%sref_improad_dif_lun     , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by impervious road per unit ground area per unit incident flux
         sref_perroad_dif     => solarabs_inst%sref_perroad_dif_lun     , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by pervious road per unit ground area per unit incident flux
         sref_br_tree_dif     => solarabs_inst%sref_br_tree_dif_lun     , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by below-roof tree per unit ground area per unit incident flux
         sref_ar_tree_dif     => solarabs_inst%sref_ar_tree_dif_lun     , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by above-roof tree per unit ground area per unit incident flux
         sref_tree_dir   => solarabs_inst%sref_tree_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar reflected  by pervious road per unit ground area per unit incident flux
         sref_tree_dif   => solarabs_inst%sref_tree_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected  by pervious road per unit ground area per unit incident flux

         sdir_road               => solarabs_inst%sdir_road_lun              , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit direct flux
         sdir_sunwall               => solarabs_inst%sdir_sunwall_lun            , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit direct flux
         sdir_shadewall               => solarabs_inst%sdir_shadewall_lun            , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit direct flux
         sdir_roof                => solarabs_inst%sdir_roof_lun               , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit direct flux
         sdir_ar_tree                => solarabs_inst%sdir_ar_tree_lun               , & ! Output:  [real(r8) (:,:) ]  flux absorbed by above-roof canopy before reaching to roof per unit direct flux
         sdir_br_tree                => solarabs_inst%sdir_br_tree_lun               , & ! Output:  [real(r8) (:,:) ]  flux absorbed by below-roof canopy before reaching to roof per unit direct flux
         sdir_tree                => solarabs_inst%sdir_tree_lun               , & ! Output:  [real(r8) (:,:) ]  flux absorbed by below-roof canopy before reaching to roof per unit direct flux
         sdir_force                => solarabs_inst%sdir_force_lun               , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit direct flux

         fabd               => surfalb_inst%fabd_patch              , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit direct flux
         fabd_sun           => surfalb_inst%fabd_sun_patch          , & ! Output:  [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit direct flux
         fabd_sha           => surfalb_inst%fabd_sha_patch          , & ! Output:  [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit direct flux
         fabi               => surfalb_inst%fabi_patch              , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit diffuse flux
         fabi_sun           => surfalb_inst%fabi_sun_patch          , & ! Output:  [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit diffuse flux
         fabi_sha           => surfalb_inst%fabi_sha_patch          , & ! Output:  [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit diffuse flux
         ftdd               => surfalb_inst%ftdd_patch              , & ! Output:  [real(r8) (:,:) ]  down direct flux below canopy per unit direct flux
         ftid               => surfalb_inst%ftid_patch              , & ! Output:  [real(r8) (:,:) ]  down diffuse flux below canopy per unit direct flux
         ftii               => surfalb_inst%ftii_patch              , & ! Output:  [real(r8) (:,:) ]  down diffuse flux below canopy per unit diffuse flux
         albgrd             => surfalb_inst%albgrd_col              , & ! Output: [real(r8) (:,:) ]  urban col ground albedo (direct) 
         albgri             => surfalb_inst%albgri_col              , & ! Output: [real(r8) (:,:) ]  urban col ground albedo (diffuse)
         albd               => surfalb_inst%albd_patch              , & ! Output  [real(r8) (:,:) ]  urban pft surface albedo (direct)                         
         albi               => surfalb_inst%albi_patch              , & ! Output: [real(r8) (:,:) ]  urban pft surface albedo (diffuse)                        
! add new snicar albedo output for history files
         albd_hst           => surfalb_inst%albd_hst_patch          , & ! Output:  [real(r8) (:,:) ]  surface albedo (direct) for history files
         albi_hst           => surfalb_inst%albi_hst_patch          , & ! Output:  [real(r8) (:,:) ]  surface albedo (diffuse) for history files
         albgrd_hst         => surfalb_inst%albgrd_hst_col          , & ! Output:  [real(r8) (:,:) ]  ground albedo (direct) for history files              
         albgri_hst         => surfalb_inst%albgri_hst_col          , & ! Output:  [real(r8) (:,:) ]  ground albedo (diffuse) for history files
! end add new snicar
         begc               => bounds%begc                          , &
         endc               => bounds%endc                          , &
         
         begl               => bounds%begl                          , &
         krs1d              =>    urbanparams_inst%ksv1d_out , & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor from roof to sky [landunit, nzcanm]                 
         kws1d              =>    urbanparams_inst%kws1d_out  , & ! Input:  [real(r8) (:,:) ]    Monte carlo view factor from wall to sky[landunit]  
         kvs1d              =>    urbanparams_inst%kvs1d_out  , & ! Input:  [real(r8) (:,:) ]    Monte carlo view factor from vegetation to sky[landunit]  
         endl               => bounds%endl                            &
         )

      ! ----------------------------------------------------------------------------
      ! Solar declination and cosine solar zenith angle and zenith angle for 
      ! next time step
      ! ----------------------------------------------------------------------------
      ! Get the clock rate (ticks per second)
      call system_clock(count_rate=clock_rate)
      call system_clock(start_time)
      !write(6,*)'starting timing...'
      
      alb_br_tree_dir(:,1)=0.11_r8
      alb_br_tree_dir(:,2)=0.35_r8
      alb_br_tree_dif(:,1)=0.11_r8
      alb_br_tree_dif(:,2)=0.35_r8 

      alb_ar_tree_dir(:,1)=0.11_r8
      alb_ar_tree_dir(:,2)=0.35_r8
      alb_ar_tree_dif(:,1)=0.11_r8
      alb_ar_tree_dif(:,2)=0.35_r8      
        
      do fl = 1,num_urbanl
         l = filter_urbanl(fl)
         g = lun%gridcell(l)
         coszen(l) = surfalb_inst%coszen_col(coli(l))  ! Assumes coszen for each column are the same
         do ib = 1, numrad
            forc_solad(l,ib) = atm2lnd_inst%forc_solad_downscaled_col(coli(l),ib)  ! Assumes forced solar for each column are the same
         end do
         zen(l)    = acos(coszen(l))
      end do

       ! Initialize output because solar radiation only done if coszen > 0

      do ib = 1, numrad
         do fc = 1,num_urbanc
            c = filter_urbanc(fc)
            albgrd(c,ib) = 0._r8
            albgri(c,ib) = 0._r8
         end do

         do fp = 1,num_urbanp  
            p = filter_urbanp(fp)
            l = patch%landunit(p)
            ! Setting albedos to wall and road view factors ensures that urban
            ! albedo will scale up to 1.0
            c = patch%column(p)
            ! revise
            if (col%itype(c) == icol_sunwall) then
               albd(p,ib) = kws1d(l,1)!vf_sw(l)
               albi(p,ib) = kws1d(l,1)!vf_sw(l)
            else if (col%itype(c) == icol_shadewall) then
               albd(p,ib) = kws1d(l,1)!vf_sw(l)
               albi(p,ib) = kws1d(l,1)!vf_sw(l)
            else if (col%itype(c) == icol_road_perv .or. col%itype(c) == icol_road_imperv) then
               albd(p,ib) = krs1d(l,2)!vf_sr(l)
               albi(p,ib) = krs1d(l,2)!vf_sr(l)
            else if (col%itype(c) == icol_roof) then
               albd(p,ib) = 1._r8
               albi(p,ib) = 1._r8
            else if (col%itype(c) == icol_road_tree) then
               albd(p,ib) = 1._r8
               albi(p,ib) = 1._r8               
            endif
            fabd(p,ib)     = 0._r8
            fabd_sun(p,ib) = 0._r8
            fabd_sha(p,ib) = 0._r8
            fabi(p,ib)     = 0._r8
            fabi_sun(p,ib) = 0._r8
            fabi_sha(p,ib) = 0._r8
            if (coszen(l) > 0._r8) then
               ftdd(p,ib)  = 1._r8
            else
               ftdd(p,ib)  = 0._r8
            end if
            ftid(p,ib)     = 0._r8
            if (coszen(l) > 0._r8) then
               ftii(p,ib)  = 1._r8
            else
               ftii(p,ib)  = 0._r8
            end if
         end do
      end do

      ! ----------------------------------------------------------------------------
      ! Urban Code
      ! ----------------------------------------------------------------------------

      num_solar = 0
      do fl = 1,num_urbanl
         l = filter_urbanl(fl)
         if (coszen(l) > 0._r8) num_solar = num_solar + 1
      end do
      ! revise
      do ib = 1,numrad
         do fl = 1,num_urbanl
            l = filter_urbanl(fl)

            ! Setting sref to wall and road view factors ensures that urban
            ! albedo will scale up to 1.0
            
            sdir_road(l,ib)      = 0._r8
            sdir_sunwall(l,ib)      = 0._r8
            sdir_shadewall(l,ib)      = 0._r8
            sdir_roof(l,ib)   = 0._r8
            sdir_ar_tree(l,ib)   = 0._r8
            sdir_br_tree(l,ib)   = 0._r8
            sdir_tree(l,ib)   = 0._r8

            sdir_force(l,ib)   = 0._r8
            
            sabs_roof_dir(l,ib)      = 0._r8
            sabs_roof_dif(l,ib)      = 0._r8
            sabs_sunwall_dir(l,ib)   = 0._r8
            sabs_sunwall_dif(l,ib)   = 0._r8
            sabs_shadewall_dir(l,ib) = 0._r8
            sabs_shadewall_dif(l,ib) = 0._r8
            sabs_improad_dir(l,ib)   = 0._r8
            sabs_improad_dif(l,ib)   = 0._r8
            sabs_perroad_dir(l,ib)   = 0._r8
            sabs_perroad_dif(l,ib)   = 0._r8
            sabs_br_tree_dir(l,ib)   = 0._r8
            sabs_br_tree_dif(l,ib)   = 0._r8
            sabs_ar_tree_dir(l,ib)   = 0._r8
            sabs_ar_tree_dif(l,ib)   = 0._r8
            sabs_tree_dir(l,ib)   = 0._r8
            sabs_tree_dif(l,ib)   = 0._r8

            sref_roof_dir(l,ib)      = 1._r8
            sref_roof_dif(l,ib)      = 1._r8
            
            ! Setting sref to wall and road view factors ensures that urban
            ! albedo will scale up to 1.0
            sref_sunwall_dir(l,ib)   = kws1d(l,1)
            sref_sunwall_dif(l,ib)   = kws1d(l,1)
            sref_shadewall_dir(l,ib) = kws1d(l,1)
            sref_shadewall_dif(l,ib) = kws1d(l,1)
            sref_improad_dir(l,ib)   = krs1d(l,2)
            sref_improad_dif(l,ib)   = krs1d(l,2)
            sref_perroad_dir(l,ib)   = krs1d(l,2)
            sref_perroad_dif(l,ib)   = krs1d(l,2)
            sref_br_tree_dir(l,ib)   = kvs1d(l,1)
            sref_br_tree_dif(l,ib)   = kvs1d(l,1)
            sref_ar_tree_dir(l,ib)   = kvs1d(l,2)
            sref_ar_tree_dif(l,ib)   = kvs1d(l,2)
            sref_tree_dir(l,ib)   = kvs1d(l,2)
            sref_tree_dif(l,ib)   = kvs1d(l,2)

         end do
      end do

      ! ----------------------------------------------------------------------------
      ! Only do the rest if all coszen are positive 
      ! ----------------------------------------------------------------------------

      if (num_solar > 0)then

         ! Set constants - solar fluxes are per unit incoming flux

         do ib = 1,numrad
            do fl = 1,num_urbanl
               l = filter_urbanl(fl)
               sdir(l,ib) = 1._r8
               sdif(l,ib) = 1._r8
            end do
         end do

         ! Incident direct beam radiation for 
         ! (a) roof and (b) road and both walls in urban canyon

         if (num_urbanl > 0) then
           call incident_direct (bounds, &
               num_urbanl, filter_urbanl, &
               canyon_hwr(begl:endl), &
               wtlunit_roof(begl:endl), &
               ht_roof(begl:endl), &
               coszen(begl:endl), &
               zen(begl:endl), &
               sdir(begl:endl, :), &
               forc_solad(begl:endl, :), &
               sdir_road(begl:endl, :), &
               sdir_sunwall(begl:endl, :), &
               sdir_shadewall(begl:endl, :), &
               sdir_roof(begl:endl, :),&         
               sdir_ar_tree(begl:endl, :),&
               sdir_br_tree(begl:endl, :),&
               sdir_force(begl:endl, :),&
               h1(begl:endl),&
               h2(begl:endl),&
               A_v2(begl:endl),&
               A_v1(begl:endl))   
                       
         end if

         ! Incident diffuse radiation for 
         ! (a) roof and (b) road and both walls in urban canyon.
         if (num_urbanl > 0) then
            call incident_diffuse (bounds, &
                 num_urbanl, filter_urbanl, &
                 canyon_hwr(begl:endl), &
                 ht_roof(begl:endl), &
                 A_v1(begl:endl), &
                 A_v2(begl:endl), &               
                 wtlunit_roof(begl:endl), &
                 sdif(begl:endl, :), &
                 sdif_road(begl:endl, :), &
                 sdif_sunwall(begl:endl, :), &
                 sdif_shadewall(begl:endl, :), &
                 sdif_roof(begl:endl, :), &
                 sdif_ar_tree(begl:endl, :), &
                 sdif_br_tree(begl:endl, :), &
                 urbanparams_inst,solarabs_inst)
         end if
         
         ! Get snow albedos for roof and impervious and pervious road
         if (num_urbanl > 0) then
            ic = 0
            call SnowAlbedo(bounds, &
                 num_urbanc, filter_urbanc, &
                 coszen(begl:endl), &
                 ic, &
                 albsnd_roof(begl:endl, :), &
                 albsnd_improad(begl:endl, :), &
                 albsnd_perroad(begl:endl, :), &
                 albsnd_tree(begl:endl, :), &
                 waterstatebulk_inst)

            ic = 1
            call SnowAlbedo(bounds, &
                 num_urbanc, filter_urbanc, &
                 coszen(begl:endl), &
                 ic, &
                 albsni_roof(begl:endl, :), &
                 albsni_improad(begl:endl, :), &
                 albsni_perroad(begl:endl, :), &
                 albsni_tree(begl:endl, :), &
                 waterstatebulk_inst)
         end if

         ! Combine snow-free and snow albedos
         ! revise
         do ib = 1,numrad
            do fc = 1,num_urbanc
               c = filter_urbanc(fc)
               l = col%landunit(c)
               if (ctype(c) == icol_roof) then    
                  alb_roof_dir_s(l,ib) = alb_roof_dir(l,ib)*(1._r8-frac_sno(c))  &
                       + albsnd_roof(l,ib)*frac_sno(c)
                  alb_roof_dif_s(l,ib) = alb_roof_dif(l,ib)*(1._r8-frac_sno(c))  &
                       + albsni_roof(l,ib)*frac_sno(c)
               else if (ctype(c) == icol_road_imperv) then    
                  alb_improad_dir_s(l,ib) = alb_improad_dir(l,ib)*(1._r8-frac_sno(c))  &
                       + albsnd_improad(l,ib)*frac_sno(c)
                  alb_improad_dif_s(l,ib) = alb_improad_dif(l,ib)*(1._r8-frac_sno(c))  &
                       + albsni_improad(l,ib)*frac_sno(c)
               else if (ctype(c) == icol_road_perv) then    
                  alb_perroad_dir_s(l,ib) = alb_perroad_dir(l,ib)*(1._r8-frac_sno(c))  &
                       + albsnd_perroad(l,ib)*frac_sno(c)
                  alb_perroad_dif_s(l,ib) = alb_perroad_dif(l,ib)*(1._r8-frac_sno(c))  &
                       + albsni_perroad(l,ib)*frac_sno(c)   
               else if (ctype(c) == icol_road_tree) then    
                  alb_br_tree_dir_s(l,ib) = alb_br_tree_dir(l,ib)*(1._r8-frac_sno(c))  &
                      + albsnd_roof(l,ib)*frac_sno(c)
                  alb_br_tree_dif_s(l,ib) = alb_br_tree_dif(l,ib)*(1._r8-frac_sno(c))  &
                      + albsni_roof(l,ib)*frac_sno(c)
                  alb_ar_tree_dir_s(l,ib) = alb_ar_tree_dir(l,ib)*(1._r8-frac_sno(c))  &
                     + albsnd_roof(l,ib)*frac_sno(c)
                  alb_ar_tree_dif_s(l,ib) = alb_ar_tree_dif(l,ib)*(1._r8-frac_sno(c))  &
                     + albsni_roof(l,ib)*frac_sno(c)                        
                                  
               end if
            end do
         end do

                  
         ! Reflected and absorbed solar radiation per unit incident radiation 
         ! for road and both walls in urban canyon allowing for multiple reflection
         ! Reflected and absorbed solar radiation per unit incident radiation for roof

         if (num_urbanl > 0) then
             call net_solar(bounds, &
                  num_urbanl, filter_urbanl, &
                  coszen             (begl:endl), &
                  canyon_hwr         (begl:endl), &
                  ht_roof             (begl:endl), &
                  A_v1(begl:endl), &
                  A_v2(begl:endl), &               
                  wtlunit_roof       (begl:endl), &
                  wtroad_perv        (begl:endl), &
                  wtroad_tree        (begl:endl), &
                  sdir               (begl:endl, :), &
                  sdif               (begl:endl, :), &
                  alb_improad_dir_s  (begl:endl, :), &
                  alb_perroad_dir_s  (begl:endl, :), &
                  alb_wall_dir       (begl:endl, :), &
                  alb_roof_dir_s     (begl:endl, :), &
                  alb_br_tree_dir_s  (begl:endl, :), &
                  alb_ar_tree_dir_s  (begl:endl, :), &
                  alb_improad_dif_s  (begl:endl, :), &
                  alb_perroad_dif_s  (begl:endl, :), &
                  alb_wall_dif       (begl:endl, :), &
                  alb_roof_dif_s     (begl:endl, :), &
                  alb_br_tree_dif_s  (begl:endl, :), &
                  alb_ar_tree_dif_s  (begl:endl, :), &                  
                  sdir_road          (begl:endl, :), &
                  sdir_sunwall       (begl:endl, :), &
                  sdir_shadewall     (begl:endl, :),  &
                  sdir_roof          (begl:endl, :), &
                  sdif_road          (begl:endl, :), &
                  sdif_sunwall       (begl:endl, :), &
                  sdif_shadewall     (begl:endl, :),  &
                  sdif_roof     (begl:endl, :),  &
                  sdir_br_tree     (begl:endl, :),  &
                  sdir_ar_tree     (begl:endl, :),  &
                  sdif_br_tree     (begl:endl, :),  &
                  sdif_ar_tree     (begl:endl, :),  &                  
                  urbanparams_inst, solarabs_inst) 

         end if

         ! ----------------------------------------------------------------------------
         ! Map urban output to surfalb_inst components 
         ! ----------------------------------------------------------------------------

         !  Set albgrd and albgri (ground albedos) and albd and albi (surface albedos)
         ! revise
         do ib = 1,numrad
            do fc = 1,num_urbanc
               c = filter_urbanc(fc)
               l = col%landunit(c)
               if (ctype(c) == icol_roof) then    
                  albgrd(c,ib) = sref_roof_dir(l,ib) 
                  albgri(c,ib) = sref_roof_dif(l,ib) 
               else if (ctype(c) == icol_sunwall) then   
                  albgrd(c,ib) = sref_sunwall_dir(l,ib)
                  albgri(c,ib) = sref_sunwall_dif(l,ib)
               else if (ctype(c) == icol_shadewall) then 
                  albgrd(c,ib) = sref_shadewall_dir(l,ib)
                  albgri(c,ib) = sref_shadewall_dif(l,ib)
               else if (ctype(c) == icol_road_perv) then
                  albgrd(c,ib) = sref_perroad_dir(l,ib)
                  albgri(c,ib) = sref_perroad_dif(l,ib)
               else if (ctype(c) == icol_road_tree) then
                  albgrd(c,ib) = sref_tree_dir(l,ib)
                  albgri(c,ib) = sref_tree_dif(l,ib)
               else if (ctype(c) == icol_road_imperv) then
                  albgrd(c,ib) = sref_improad_dir(l,ib)
                  albgri(c,ib) = sref_improad_dif(l,ib)
               endif
! add new snicar albedo variables for history fields
               if (coszen(l) > 0._r8) then
                  albgrd_hst(c,ib) = albgrd(c,ib)
                  albgri_hst(c,ib) = albgri(c,ib)
               end if
! end add new snicar
            end do
            do fp = 1,num_urbanp
               p = filter_urbanp(fp)
               c = patch%column(p)
               l = patch%landunit(p)
               albd(p,ib) = albgrd(c,ib)
               albi(p,ib) = albgri(c,ib)
! add new snicar albedo variables for history fields
               if (coszen(l) > 0._r8) then
                  albd_hst(p,ib) = albd(p,ib)
                  albi_hst(p,ib) = albi(p,ib)
               end if
! end add new snicar
            end do
         end do
      end if

    end associate
    call system_clock(end_time)        
    ! Calculate elapsed time in seconds
    ! elapsed_time = real(end_time - start_time,r8) / real(clock_rate,r8)
    
    ! write the elapsed time for each iteration  
    ! write(*, '(A, I0, A, F6.3)') 'Finished urban albedo calculation: ', elapsed_time, ' seconds'

  end subroutine UrbanAlbedo

  !-----------------------------------------------------------------------
  subroutine SnowAlbedo (bounds          , &
       num_urbanc, filter_urbanc, coszen, ind , &
       albsn_roof, albsn_improad, albsn_perroad, albsn_tree , &
       waterstatebulk_inst)
    !
    ! !DESCRIPTION:
    ! Determine urban snow albedos
    !
    ! !USES:
    use column_varcon, only : icol_roof, icol_road_perv, icol_road_imperv, icol_road_tree
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds                     
    integer , intent(in) :: num_urbanc                          ! number of urban columns in clump
    integer , intent(in) :: filter_urbanc(:)                    ! urban column filter
    integer , intent(in) :: ind                                 ! 0=direct beam, 1=diffuse radiation
    real(r8), intent(in) :: coszen        ( bounds%begl: )      ! cosine solar zenith angle [landunit]
    real(r8), intent(out):: albsn_roof    ( bounds%begl: , 1: ) ! roof snow albedo by waveband [landunit, numrad]
    real(r8), intent(out):: albsn_improad ( bounds%begl: , 1: ) ! impervious road snow albedo by waveband [landunit, numrad]
    real(r8), intent(out):: albsn_perroad ( bounds%begl: , 1: ) ! pervious road snow albedo by waveband [landunit, numrad]
    real(r8), intent(out):: albsn_tree ( bounds%begl: , 1: ) ! road tree snow albedo by waveband [landunit, numrad]
    type(waterstatebulk_type), intent(in) :: waterstatebulk_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fc,c,l              ! indices
    real(r8) :: h2osno_total(bounds%begc:bounds%endc)  ! total snow water (mm H2O)
    !
    ! These values are derived from Marshall (1989) assuming soot content of 1.5e-5 
    ! (three times what LSM uses globally). Note that snow age effects are ignored here.
    real(r8), parameter :: snal0 = 0.66_r8 ! vis albedo of urban snow
    real(r8), parameter :: snal1 = 0.56_r8 ! nir albedo of urban snow
    !-----------------------------------------------------------------------

    ! this code assumes that numrad = 2 , with the following
    ! index values: 1 = visible, 2 = NIR
    SHR_ASSERT_ALL_FL(numrad == 2, sourcefile, __LINE__)

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(coszen)        == (/bounds%endl/)),         sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(albsn_roof)    == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(albsn_improad) == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(albsn_perroad) == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(albsn_tree) == (/bounds%endl, numrad/)), sourcefile, __LINE__)

    call waterstatebulk_inst%CalculateTotalH2osno(bounds, num_urbanc, filter_urbanc, &
         caller = 'UrbanAlbedoMod:SnowAlbedo', &
         h2osno_total = h2osno_total(bounds%begc:bounds%endc))

    do fc = 1,num_urbanc
       c = filter_urbanc(fc)
       l = col%landunit(c)
       if (coszen(l) > 0._r8 .and. h2osno_total(c) > 0._r8) then
          if (col%itype(c) == icol_roof) then
             albsn_roof(l,1) = snal0
             albsn_roof(l,2) = snal1
          else if (col%itype(c) == icol_road_imperv) then
             albsn_improad(l,1) = snal0
             albsn_improad(l,2) = snal1
          else if (col%itype(c) == icol_road_perv) then
             albsn_perroad(l,1) = snal0
             albsn_perroad(l,2) = snal1
          else if (col%itype(c) == icol_road_tree) then
             albsn_tree(l,1) = snal0
             albsn_tree(l,2) = snal1
          end if
       else
          if (col%itype(c) == icol_roof) then
             albsn_roof(l,1) = 0._r8
             albsn_roof(l,2) = 0._r8
          else if (col%itype(c) == icol_road_imperv) then
             albsn_improad(l,1) = 0._r8
             albsn_improad(l,2) = 0._r8
          else if (col%itype(c) == icol_road_perv) then
             albsn_perroad(l,1) = 0._r8
             albsn_perroad(l,2) = 0._r8
          else if (col%itype(c) == icol_road_tree) then
             albsn_tree(l,1) = 0._r8
             albsn_tree(l,2) = 0._r8       
          end if
       end if
    end do

  end subroutine SnowAlbedo

!-----------------------------------------------------------------------  
  function gaussian_quadrature(f, a, b, x, w,l) result(integral)
    ! calculate integration using gaussian quadrature nodes
    implicit none
    real(r8), intent(in) :: a, b        ! integration boundaries
    real(r8), intent(in) :: x(:), w(:)  ! nodes and weights
    real(r8) :: integral, sum
    integer, intent(in) :: l
    integer :: i, n
    
    interface
        function f(theta, l) result(value)
            implicit none
            integer, intent(in) :: l 
            !real(r8), intent(in) :: theta
            !real(r8) :: value
            real(selected_real_kind(12)), intent(in) :: theta
            real(selected_real_kind(12)) :: value            
        end function f
    end interface

    n = size(x)
    sum = 0.0_r8
    do i = 1, n
        sum = sum + w(i) * f(0.5_r8 * (b - a) * x(i) + 0.5_r8 * (b + a), l)
    end do
    
    integral = 0.5_8 * (b - a) * sum
  end function gaussian_quadrature
  

!---------------------------------------------------------------------------------
  subroutine incident_direct (bounds,num_urbanl, filter_urbanl, &
       canyon_hwr, wtlunit_roof,ht_roof,coszen, zen,sdir, forc_solad,sdir_road, sdir_sunwall,&
        sdir_shadewall, sdir_roof, sdir_ar_tree, sdir_br_tree,&
         sdir_force,h1,h2,A_v2,A_v1)

    !
    ! !DESCRIPTION:
    ! Direct beam solar radiation incident on urban canyon surfaces (road, walls, roof, tree)
    ! 
    !
    !                           Sun
    !                  ........  /
    !             roof ........ /
    !            ------..tree. /---            -
    !                 |...... /.|              |
    !    sunlit wall  |..... /..| shaded wall  h
    !                 |     /   |              |
    !                 |    /    |              |
    !                 ----------               -
    !                    road
    !                 <--- w --->
    !
    ! The incident direct beam solar radiation expressions are derived based on the analytical solution in Masson, 2000 
    ! with additional tree attenuation terms.
    ! 1. Masson, V. (2000) A physically-based scheme for the urban energy budget in 
    ! atmospheric models. Boundary-Layer Meteorology 94:357-397
    ! 2. ***CLMU tree citation***
    !  
    ! 
    ! !USES:
    use clm_varcon, only : rpi

    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds                      
    integer , intent(in)  :: num_urbanl                          ! number of urban landunits
    integer , intent(in)  :: filter_urbanl(:)                    ! urban landunit filter
    real(r8), intent(in)  :: canyon_hwr( bounds%begl: )          ! ratio of building height to street width [landunit]
    real(r8), intent(in)  :: wtlunit_roof( bounds%begl: )        ! weight of roof with respect to landunit  [landunit]
    real(r8), intent(in)  :: ht_roof( bounds%begl: )             ! height of urban roof (m) [landunit]
    real(r8), intent(in) :: h1( bounds%begl: )      ! direct beam solar radiation incident attenuated by above-roof tree before reaching to roof,wall, and road, per unit incident flux [landunit, numrad]
    real(r8), intent(in) :: h2( bounds%begl: )      ! direct beam solar radiation incident attenuated by below-roof tree before reaching to wall and road, per unit incident flux [landunit, numrad]
    real(r8), intent(in) :: A_v1( bounds%begl: )      ! direct beam solar radiation incident attenuated by above-roof tree before reaching to roof,wall, and road, per unit incident flux [landunit, numrad]
    real(r8), intent(in) :: A_v2( bounds%begl: )      ! direct beam solar radiation incident attenuated by below-roof tree before reaching to wall and road, per unit incident flux [landunit, numrad]
    real(r8), intent(in)  :: coszen( bounds%begl: )              ! cosine solar zenith angle [landunit]
    real(r8), intent(in)  :: zen( bounds%begl: )                 ! solar zenith angle (radians) [landunit]
    real(r8), intent(in)  :: sdir( bounds%begl: , 1: )             ! direct beam solar radiation incident on horizontal surface [landunit, numrad]
    real(r8), intent(out) :: sdir_road( bounds%begl: , 1: )      ! direct beam solar radiation incident on road per unit incident flux with tree attenuation [landunit, numrad]
    real(r8), intent(out) :: sdir_sunwall( bounds%begl: , 1: )   ! direct beam solar radiation (per unit wall area) incident on sunlit wall per unit incident flux with tree attenuation [landunit, numrad]
    real(r8), intent(out) :: sdir_shadewall( bounds%begl: , 1: ) ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux with tree attenuation [landunit, numrad]
    real(r8), intent(out) :: sdir_roof( bounds%begl: , 1: )      ! direct beam solar radiation (per unit wall area) incident on roof per unit incident flux with tree attenuation [landunit, numrad] 
    real(r8), intent(out) :: sdir_ar_tree( bounds%begl: , 1: )      ! direct beam solar radiation incident attenuated by above-roof tree before reaching to roof,wall, and road, per unit incident flux [landunit, numrad]
    real(r8), intent(out) :: sdir_br_tree( bounds%begl: , 1: )      ! direct beam solar radiation incident attenuated by below-roof tree before reaching to wall and road, per unit incident flux [landunit, numrad]
    ! temporary incoming incident direct solar flux for diagnosis purpose only   
    real(r8), intent(out) :: sdir_force( bounds%begl: , 1: )       ! direct beam radiation  (vis=forc_sols , nir=forc_soll ) for diagnosis   
    real(r8), intent(in) :: forc_solad(bounds%begl:, :)            ! direct beam radiation  (vis=forc_sols , nir=forc_soll ) for diagnosis     

    ! !LOCAL VARIABLES:
    integer  :: fl,l,g,i,ib                     ! indices
    ! In this scheme, we can perform numerical check of gaussian quadrature integration 
    !logical  :: numchk = .false.                ! true => perform numerical check of analytical solution
    real(r8) :: tanzen(bounds%begl:bounds%endl) ! tan(zenith angle)
    real(r8) :: sinzen(bounds%begl:bounds%endl) ! sin(zenith angle)  
    real(r8)  ::  theta_max                      ! pi/2
    ! check if these minimum angles are really necessary  
    real(r8) :: tanzen_min                      ! minimum tan(zenith angle) (>0)
    real(r8) :: coszen_min                      ! minimum cos(zenith angle) (>0)
    real(r8) :: sinzen_min                      ! minimum sin(zenith angle) (>0)
    real(r8) :: min_zen                         ! minimum threshold for tan(zen),cos(zen), and sin(zen)
    ! I haven't added conservation check yet
    real(r8) :: swall_projected_t               ! direct beam solar radiation (per unit ground area) incident on wall
    real(r8) :: err1(bounds%begl:bounds%endl)   ! energy conservation error
    real(r8) :: sumr                            ! sum of sroad for each orientation (0 <= theta <= pi/2)
    real(r8) :: sumw                            ! sum of swall for each orientation (0 <= theta <= pi/2)
    real(r8) :: num                             ! number of orientations
    
    real(r8) :: theta                           ! canyon orientation relative to sun (0 <= theta <= pi/2)
    real(r8) :: wbui(bounds%begl:bounds%endl)   ! building width
    real(r8) :: wcan(bounds%begl:bounds%endl)   ! street width
    real(r8) :: omega(bounds%begl:bounds%endl)  ! leaf clumping index (Eq. 15 in Krayenhoff et al. 2020)     
    real(r8) :: LAD(bounds%begl:bounds%endl)    ! Leaf area density in the canyon column [m-1]
    real(r8) :: Tree_at(bounds%begl:bounds%endl)       ! tree attenuation term (Kbs * Omega * LAD)
    real(r8) :: Tree_abr(bounds%begl:bounds%endl)      ! tree height above roof (h1+h2-ht_roof)
    integer  :: ca_order(bounds%begl:bounds%endl)      ! order of critical zenith angles [0,1,2,3,4]
    !real(r8) ::  sdir_roof_shaded(bounds%begl:bounds%endl)! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux on the shaded side 
    real(r8) ::  sdir_tree_roof(bounds%begl:bounds%endl,numrad)! direct beam solar radiation (per unit wall area) attenuated by above roof tree before reaching to roof
    real(r8) ::  sdir_tree_road(bounds%begl:bounds%endl,numrad)! direct beam solar radiation (per unit wall area) attenuated by above roof tree before reaching to roof
    real(r8) ::  sdir_tree_sunwall(bounds%begl:bounds%endl,numrad)! direct beam solar radiation (per unit wall area) attenuated by above roof tree before reaching to roof
    real(r8) ::  sdir_ar_tree_wr(bounds%begl:bounds%endl,numrad)! direct beam solar radiation (per unit wall area) attenuated by above roof tree before reaching to roof, wall or road 
    real(r8) ::  sdir_br_tree_wr(bounds%begl:bounds%endl,numrad)! direct beam solar radiation (per unit wall area) attenuated by below roof tree before reaching to wall or road 
    real(r8) ::  sdir_after_ar(bounds%begl:bounds%endl,numrad)! direct beam solar radiation (per unit wall area) attenuated by above roof tree before reaching to wall or road 

    real(r8)  :: Kbs                                    ! Extinction coefficient for vegetation foliage
    ! We could add another function to calculate the gaussian quadrature weights and nodes for any numbers of nodes
    real(r8)  :: Gauss_nodes(4), Gauss_weights(4)       ! Nodes and weights used in gaussian quadrature integration
    
    ! theta : angle between the sun direction and the along-canyon axis; zen: solar zenith angle
    ! below are critical theta and zenith angles for incident direct deam solar calculation on roof, wall, and road
    real(r8) :: theta_ar_wr5(bounds%begl:bounds%endl) ! The 5th critical theta angle for calculating wall & road fluxes in tree-above-roof scenario
    real(r8) :: theta_ar_wr4(bounds%begl:bounds%endl) ! The 4th critical theta angle for calculating wall & road fluxes in tree-above-roof scenario
    real(r8) :: theta_ar_wr3(bounds%begl:bounds%endl) ! The 3th critical theta angle for calculating wall & road fluxes in tree-above-roof scenario
    real(r8) :: theta_ar_wr2(bounds%begl:bounds%endl) ! The 2nd critical theta angle for calculating wall & road fluxes in tree-above-roof scenario
    real(r8) :: theta_ar_wr1(bounds%begl:bounds%endl) ! The 1st critical theta angle for calculating wall & road fluxes in tree-above-roof scenario
    real(r8) :: zen_ar_wr5(bounds%begl:bounds%endl)   ! The 5th critical zenith angle for calculating wall & road fluxes in tree-above-roof scenario
    real(r8) :: zen_ar_wr4(bounds%begl:bounds%endl)   ! The 4th critical zenith angle for calculating wall & road fluxes in tree-above-roof scenario
    real(r8) :: zen_ar_wr3(bounds%begl:bounds%endl)   ! The 3th critical zenith angle for calculating wall & road fluxes in tree-above-roof scenario
    real(r8) :: zen_ar_wr2(bounds%begl:bounds%endl)   ! The 2nd critical zenith angle for calculating wall & road fluxes in tree-above-roof scenario
    real(r8) :: zen_ar_wr1(bounds%begl:bounds%endl)   ! The 1st critical zenith angle for calculating wall & road fluxes in tree-above-roof scenario
    
    real(r8) :: theta_br_wr3(bounds%begl:bounds%endl) ! The 3th critical theta angle for calculating wall & road fluxes in tree-below-roof scenario
    real(r8) :: theta_br_wr2(bounds%begl:bounds%endl) ! The 2nd critical theta angle for calculating wall & road fluxes in tree-below-roof scenario
    real(r8) :: theta_br_wr1(bounds%begl:bounds%endl) ! The 1st critical theta angle for calculating wall & road fluxes in tree-below-roof scenario
    real(r8) :: zen_br_wr3(bounds%begl:bounds%endl)   ! The 3th critical zenith angle for calculating wall & road fluxes in tree-below-roof scenario
    real(r8) :: zen_br_wr2(bounds%begl:bounds%endl)   ! The 2nd critical zenith angle for calculating wall & road fluxes in tree-below-roof scenario
    real(r8) :: zen_br_wr1(bounds%begl:bounds%endl)   ! The 1st critical zenith angle for calculating wall & road fluxes in tree-below-roof scenario
        
    real(r8) :: theta_ar_roof3(bounds%begl:bounds%endl) ! The 3th critical theta angle for calculating roof fluxes in tree-above-roof scenario
    real(r8) :: theta_ar_roof2(bounds%begl:bounds%endl) ! The 2nd critical theta angle for calculating roof fluxes in tree-above-roof scenario
    real(r8) :: theta_ar_roof1(bounds%begl:bounds%endl) ! The 1st critical theta angle for calculating roof fluxes in tree-above-roof scenario
    real(r8) :: zen_ar_roof3(bounds%begl:bounds%endl)   ! The 3th critical zenith angle for calculating roof fluxes in tree-above-roof scenario
    real(r8) :: zen_ar_roof2(bounds%begl:bounds%endl)   ! The 2nd critical zenith angle for calculating roof fluxes in tree-above-roof scenario
    real(r8) :: zen_ar_roof1(bounds%begl:bounds%endl)   ! The 1st critical zenith angle for calculating roof fluxes in tree-above-roof scenario

    
    
    ! these variables are used to calculate Masson2000 solution without tree. They are included to test the consistency between the new and old scheme when LAI=0
    real(r8) :: theta0(bounds%begl:bounds%endl)             ! critical canyon orientation for which road is no longer illuminated
    real(r8) :: sdir_road_o(bounds%begl:bounds%endl , numrad)     ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux on the shaded side 
    real(r8) :: sdir_sunwall_o(bounds%begl:bounds%endl , numrad)  ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux on the shaded side 
    real(r8) :: sdir_shadewall_o(bounds%begl:bounds%endl , numrad) ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux on the shaded side 

    real(r8) :: latdeg                                 ! latdeg
    real(r8) :: londeg                                 ! londeg    
    logical  :: debug_write_dir = .false.                  ! true => write out many intermediate variables for debugging
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(canyon_hwr)     == (/bounds%endl/)),         sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(coszen)         == (/bounds%endl/)),         sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(zen)            == (/bounds%endl/)),         sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sdir)           == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sdir_road)      == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sdir_sunwall)   == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sdir_shadewall) == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(wtlunit_roof)     == (/bounds%endl/)),         sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(ht_roof)     == (/bounds%endl/)),         sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sdir_roof) == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    
    theta_max = rpi / 2.0_r8
    debug_write_dir = .false.!.true.!.false.
    !--------------hard coded values--------------!  
    !ca_order=4
    ! B<W
    !wbui=7.0_r8
    !wcan=10.0_r8
    !h1=4.0_r8
    !h2=7.0_r8
    !ht_roof2=10.0_r8
	  !ca_order=4

    Omega=0.6_r8
    !LAD=0.2_r8    
    LAD=0.4_r8 
    Kbs=0.5_r8
    min_zen = 0.000001_r8
    
    if (debug_write_dir) then
       write (6,'(A,2F10.3)') 'LAD ',LAD
    end if
    
    ! Gaussian nodes and weights when n = 4
    Gauss_nodes = (/ -0.86113631_r8, -0.33998104_r8,  0.33998104_r8,  0.86113631_r8/)
    Gauss_weights = (/ 0.34785485_r8, 0.65214515_r8, 0.65214515_r8, 0.34785485_r8/)
    !-------------------------------------------------!      
    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       
       if (coszen(l) > 0._r8) then
          tanzen(l) = tan(zen(l))
          sinzen(l) = sin(zen(l))
          
          ! Create a new variable solarabs_inst%sdir_force_lun for diagnosis, which is rasier to output
          do ib = 1,numrad
              sdir_force(l,ib)=forc_solad(l,ib)
          end do
          
          ! calculate building height and street width
          wbui(l)=ht_roof(l)/(canyon_hwr(l)*(1._r8-wtlunit_roof(l))/wtlunit_roof(l))
          wcan(l)=ht_roof(l)/canyon_hwr(l) 
          
          tanzen_min= max(tanzen(l),min_zen) !tan(max(zen(l),min_zen))   
          if (debug_write_dir) then
              !write (6,'(A,2F10.3)') 'wbui(l)',wbui(l)            
              !write (6,'(A,2F10.3)') 'tanzen_min',tanzen_min
          end if                    

          if (debug_write_dir) then
              !write (6,'(A,2F10.3)') 'h1(l),h2(l)',h1(l),h2(l)   
          end if    
               
          ! two intermediate variables constantly used by incident direct solar expressions
          Tree_abr(l)=h1 (l)+ h2(l) - ht_roof(l)
          Tree_at(l)=max(Kbs * Omega(l) * LAD (l),min_zen)
          
          if (debug_write_dir) then
              !write (6,'(A,2F10.3)') 'Tree_abr(l), Tree_at(l)',Tree_abr(l), Tree_at(l) 
          end if    
          
          ! Determin the ca_order, which influences the flux expressions in the tree-above-roof scenarios
          if (h2(l)+h1(l) > 2.0_r8*ht_roof(l)) then
              write (iulog,*) 'the maximum tree height should be lower than two times of building height. h1 + h2, ht_roof =',h2(l)+h1(l),ht_roof(l)
              write (iulog,*) 'clm model is stopping'
              call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
          else if (h2(l)+h1(l) > ht_roof(l)) then  ! tree above roof
             if (h2(l) > ht_roof(l)) then
                 if (ht_roof(l) < h1(l) + 0.5_r8 * h2(l)) then
                     ca_order(l) = 1 !tree above roof; zen_ar_wr1<3<2<5<4
                 else  ! H >= h1 + 0.5 * h2
                     ca_order(l) = 2 !tree above roof; zen_ar_wr1<3<2<4<5
                 end if
             else  ! h2 <= H
                 if (ht_roof(l) < h1(l) + 0.5_r8 * h2(l)) then
                     ca_order(l) = 3 !tree above roof; zen_ar_wr1<2<3<5<4
                 else  ! H >= h1 + 0.5_r8 * h2
                     ca_order(l) = 4 !tree above roof; zen_ar_wr1<2<3<4<5
                 end if
             end if
          else if (h2(l)+h1(l) <= ht_roof(l)) then                 
               ca_order(l) = 0 ! tree below roof scenario
          end if
          
          if (h1(l) > ht_roof(l)) then
              write (iulog,*) 'the maximum tree crown bottom should be lower than the building height. h1, ht_roof = ',h1(l),ht_roof(l)
              write (iulog,*) 'clm model is stopping'
              call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
          end if 
          
          if (debug_write_dir) then
              !write (6,'(A,2I10.3)') 'ca_order(l) ',ca_order(l)
          end if             
          ! Compute critical theta values for wall and road calculations under tree-bwlow-roof scenario  
          if (ca_order(l) == 0) then !tree below roof 
             theta_br_wr1(l) = asin(min(wcan(l) / (ht_roof(l)* tanzen_min), 1._r8))  
             theta_br_wr2(l) = asin(min(wcan(l) / (max(ht_roof(l)-h1(l),min_zen) * tanzen_min), 1._r8))
             theta_br_wr3(l) = asin(min(wcan(l) / (max(ht_roof(l)-h1(l)-h2(l),min_zen) * tanzen_min), 1._r8))
            
             zen_br_wr1(l) = atan(wcan(l) / (max(ht_roof(l)-h1(l)-h2(l),min_zen)))
             zen_br_wr2(l) = atan(wcan(l) / (max(ht_roof(l)-h1(l),min_zen)))
             zen_br_wr3(l) = atan(wcan(l) / ((ht_roof(l))))
          else                 !tree above roof           
             ! Compute critical theta values for wall and road calculations under tree-above-roof scenario       
             theta_ar_wr5(l) = asin(min(wcan(l) / (max(Tree_abr(l),min_zen)* tanzen_min), 1._r8))
             theta_ar_wr4(l) = asin(min(wcan(l) / (max(ht_roof(l) - h1(l),min_zen) * tanzen_min), 1._r8))
             theta_ar_wr3(l) = asin(min(wcan(l) / (h2(l) * tanzen_min), 1._r8))
             theta_ar_wr2(l) = asin(min(wcan(l) / (ht_roof(l) * tanzen_min), 1._r8))
             theta_ar_wr1(l) = asin(min(wcan(l) / ((h1(l) + h2(l)) * tanzen_min), 1._r8))     

             ! Compute critical zenith values for wall and road calculations
             zen_ar_wr1(l) = atan(wcan(l) / (h1(l) + h2(l)))
             zen_ar_wr2(l) = atan(wcan(l) / ht_roof(l))
             zen_ar_wr3(l) = atan(wcan(l) / h2(l))
             zen_ar_wr4(l) = atan(wcan(l) / max(ht_roof(l) - h1(l),min_zen))
             zen_ar_wr5(l) = atan(wcan(l) / max(Tree_abr(l),min_zen))
          
             ! Compute critical theta values for roof calculations
             theta_ar_roof3(l) = asin(min((wbui(l) + wcan(l)) / (max(Tree_abr(l),min_zen) * tanzen_min), 1._r8))
             theta_ar_roof2(l) = asin(min(wbui(l) / (max(Tree_abr(l),min_zen) * tanzen_min), 1._r8))
             theta_ar_roof1(l) = asin(min(wcan(l) / (max(Tree_abr(l),min_zen) * tanzen_min), 1._r8))
                      
             ! Compute critical zenith values for roof calculations
             zen_ar_roof1(l) = atan(wcan(l) / max(Tree_abr(l),min_zen))
             zen_ar_roof2(l) = atan(wbui(l) / max(Tree_abr(l),min_zen))
             zen_ar_roof3(l) = atan((wbui(l) + wcan(l)) / max(Tree_abr(l),min_zen))   
          end if
      end if
    end do   
    
    ! The Masson 2000 calculation
    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       if (coszen(l) > 0._r8) then
          theta0(l) = asin(min( (1._r8/((ht_roof(l)/wcan(l))*tan(max(zen(l),0.000001_r8)))), 1._r8 ))
       end if
    end do
        
    do ib = 1,numrad
       do fl = 1,num_urbanl
          l = filter_urbanl(fl)
          g = lun%gridcell(l)
          ! (zen(l)-theta_max)>0.001_r8) seems important because some numerical errors occur when zen is too close to pi/2
          !if ((coszen(l) > 0._r8) .and. ((theta_max - zen(l))>0.001_r8)) then   
          if (coszen(l) > 0._r8) then   
             latdeg=grc%latdeg(g)
             londeg=grc%londeg(g)
             ! Initialization
             sdir_shadewall(l,ib) = 0._r8  
             sdir_road(l,ib) = 0._r8  
             sdir_sunwall(l,ib) = 0._r8  
             sdir_roof(l,ib) = 1.0_r8
             sdir_ar_tree(l,ib) = 0._r8 
             sdir_br_tree(l,ib) = 0._r8 
             
             coszen_min=max(coszen(l),min_zen)
             sinzen_min=max(sinzen(l),min_zen)
                          
             ! printputs for diagnosis
             if (debug_write_dir) then              
                 write (6,'(A,2F10.3)') '-------------------(l)------------------- ',l 
                 write (6,'(A,2F10.3)') 'zen(l) ',zen(l) 
                 write (6,'(A,2F10.3)') 'latdeg,londeg ',latdeg,londeg
                 write (6,'(A,2I10.3)') 'ca_order(l) ',ca_order(l) 
                 write (6,'(A,2F10.3)') 'wbui(l)',wbui(l)
                 write (6,'(A,2F10.3)') 'wcan(l),ht_roof(l)',wcan(l),ht_roof(l)
                 write (6,'(A,2F10.3)') 'h1(l),h2(l)',h1(l),h2(l)
                 write (6,'(A,2F10.3)') ' sinzen_min, coszen_min', sinzen_min,coszen_min
                 if (ca_order(l)==0) then
                    write (6,'(A,2F10.3)') 'theta_br_wr1(l): ',theta_br_wr1(l)
                    write (6,'(A,2F10.3)') 'theta_br_wr2(l): ',theta_br_wr2(l)
                    write (6,'(A,2F10.3)') 'theta_br_wr3(l): ',theta_br_wr3(l)
                    
                    write (6,'(A,2F10.3)') 'zen_br_wr1(l): ',zen_br_wr1(l)
                    write (6,'(A,2F10.3)') 'zen_br_wr2(l): ',zen_br_wr2(l)
                    write (6,'(A,2F10.3)') 'zen_br_wr3(l): ',zen_br_wr3(l)  
                 else
                    write (6,'(A,2F10.3)') 'theta_ar_wr1(l): ',theta_ar_wr1(l) 
                    write (6,'(A,2F10.3)') 'theta_ar_wr2(l): ',theta_ar_wr2(l) 
                    write (6,'(A,2F10.3)') 'theta_ar_wr3(l): ',theta_ar_wr3(l) 
                    write (6,'(A,2F10.3)') 'theta_ar_wr4(l): ',theta_ar_wr4(l) 
                    write (6,'(A,2F10.3)') 'theta_ar_wr5(l): ',theta_ar_wr5(l)  
                   
                    write (6,'(A,2F10.3)') 'theta_ar_roof1(l): ',theta_ar_roof1(l)
                    write (6,'(A,2F10.3)') 'theta_ar_roof2(l): ',theta_ar_roof2(l)
                    write (6,'(A,2F10.3)') 'theta_ar_roof3(l): ',theta_ar_roof3(l)
                 end if
             end if  
                                
             if (ca_order(l)==0) then !tree below roof
                 ! roof is not shaded
                 ! calculate sdir_sunwall(l,ib) and sdir_road(l,ib)
                 if (zen(l) > zen_br_wr3(l)) then
                     sdir_sunwall(l,ib)=(integration_br_wr_01(theta_br_wr1(l),l) &
                            +integration_br_wr_12(theta_br_wr1(l), theta_br_wr2(l),l) &
                            +gaussian_quadrature(integrand_br_wr_23, theta_br_wr2(l), theta_br_wr3(l), Gauss_nodes, Gauss_weights,l)+integration_br_wr_23(theta_br_wr2(l), theta_br_wr3(l),l) &
                            +(theta_max- theta_br_wr3(l)) * wcan(l) / ht_roof(l))/theta_max
                     sdir_road(l,ib)=integration_br_road(theta_br_wr1(l),l)/theta_max
                 else if (zen(l)>zen_br_wr2(l)) then
                     sdir_sunwall(l,ib)=(integration_br_wr_01(theta_br_wr1(l),l) &
                            +integration_br_wr_12(theta_br_wr1(l), theta_br_wr2(l),l) &
                            +gaussian_quadrature(integrand_br_wr_23, theta_br_wr2(l), theta_max, Gauss_nodes, Gauss_weights,l)+integration_br_wr_23(theta_br_wr2(l), theta_max,l))/theta_max 
                     sdir_road(l,ib)=integration_br_road(theta_br_wr1(l),l)/theta_max       
                 else if (zen(l)>zen_br_wr1(l)) then
                     sdir_sunwall(l,ib)=(integration_br_wr_01(theta_br_wr1(l),l) &
                            +integration_br_wr_12(theta_br_wr1(l), theta_max,l))/theta_max
                     sdir_road(l,ib)=integration_br_road(theta_br_wr1(l),l)/theta_max                  
                 else
                     sdir_sunwall(l,ib)=(integration_br_wr_01(theta_max,l))/theta_max
                     sdir_road(l,ib)=integration_br_road(theta_max,l)/theta_max
                 end if                         
             else                     ! tree above roof scenario
                if (zen(l)>zen_ar_wr5(l)) then
                    sdir_after_ar(l,ib) = (tree_partition_05(theta_ar_wr5(l),l)&
                               + gaussian_quadrature(integrand_ar_BW_12, theta_ar_wr5(l), theta_max, Gauss_nodes, Gauss_weights,l)) / theta_max
                else      
                    sdir_after_ar(l,ib) = (tree_partition_05(theta_max,l)) / theta_max            
                end if  
                ! calcluate incident direct solar flux of roof
                if  (wbui(l) > wcan(l)) then        !zen_ar_roof1<2<3
                    if (zen(l)>zen_ar_roof3(l)) then
                        sdir_roof(l,ib)=(integration_ar_BW_WB_01_02(theta_ar_roof1(l),l)&
                                     +gaussian_quadrature(integrand_ar_BW_12, theta_ar_roof1(l), theta_ar_roof2(l),Gauss_nodes, Gauss_weights, l)+integration_ar_BW_12(theta_ar_roof1(l), theta_ar_roof2(l),l)&
                                     +gaussian_quadrature(integrand_ar_BW_WB_13_23, theta_ar_roof2(l), theta_ar_roof3(l), Gauss_nodes, Gauss_weights,l)+&
                                     gaussian_quadrature(integrand_ar_WB_BW_3pi, theta_ar_roof3(l),theta_max,Gauss_nodes, Gauss_weights,l))/theta_max
                    else if (zen(l)>zen_ar_roof2(l)) then
                        sdir_roof(l,ib)=(integration_ar_BW_WB_01_02(theta_ar_roof1(l),l)&
                                      +gaussian_quadrature(integrand_ar_BW_12, theta_ar_roof1(l), theta_ar_roof2(l),Gauss_nodes, Gauss_weights, l)+integration_ar_BW_12(theta_ar_roof1(l), theta_ar_roof2(l),l)&
                                      +gaussian_quadrature(integrand_ar_BW_WB_13_23, theta_ar_roof2(l),theta_max, Gauss_nodes, Gauss_weights,l))/theta_max                                  
                    else if (zen(l)>zen_ar_roof1(l)) then
                        sdir_roof(l,ib)=(integration_ar_BW_WB_01_02(theta_ar_roof1(l),l)&
                                      +gaussian_quadrature(integrand_ar_BW_12, theta_ar_roof1(l), theta_max,Gauss_nodes, Gauss_weights, l)+integration_ar_BW_12(theta_ar_roof1(l), theta_max,l))/theta_max 
                    else
                        sdir_roof(l,ib)=(integration_ar_BW_WB_01_02(theta_max,l))/theta_max
                    end if 
                else   !wbui < wcan and zen_ar_roof2<1<3
                    if (zen(l)>zen_ar_roof3(l)) then
                        sdir_roof(l,ib)=(integration_ar_BW_WB_01_02(theta_ar_roof2(l),l)&
                                    +gaussian_quadrature(integrand_ar_WB_21, theta_ar_roof2(l), theta_ar_roof1(l),Gauss_nodes, Gauss_weights, l)+&
                                   +gaussian_quadrature(integrand_ar_BW_WB_13_23, theta_ar_roof1(l), theta_ar_roof3(l), Gauss_nodes, Gauss_weights,l)+&
                                   gaussian_quadrature(integrand_ar_WB_BW_3pi, theta_ar_roof3(l),theta_max,Gauss_nodes, Gauss_weights,l))/theta_max
                    else if (zen(l)>zen_ar_roof1(l)) then
                        sdir_roof(l,ib)=(integration_ar_BW_WB_01_02(theta_ar_roof2(l),l)&
                                   +gaussian_quadrature(integrand_ar_WB_21, theta_ar_roof2(l), theta_ar_roof1(l),Gauss_nodes, Gauss_weights, l)+&
                                   gaussian_quadrature(integrand_ar_BW_WB_13_23, theta_ar_roof1(l), theta_max, Gauss_nodes, Gauss_weights,l))/theta_max
                    else if (zen(l)>zen_ar_roof2(l)) then
                        sdir_roof(l,ib)=(integration_ar_BW_WB_01_02(theta_ar_roof2(l),l)&
                                  +gaussian_quadrature(integrand_ar_WB_21, theta_ar_roof2(l), theta_max, Gauss_nodes, Gauss_weights, l))/theta_max
                    else
                        sdir_roof(l,ib)=(integration_ar_BW_WB_01_02(theta_max,l))/theta_max  
                    end if            
               end if 
                              
               if (ca_order(l)==4) then !zen_ar_wr1<2<3<4<5
                  if (zen(l)>zen_ar_wr5(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l)&
                                 + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l)+ integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)&
                                 + integration_ar_wr_23_ca43(theta_ar_wr2(l), theta_ar_wr3(l),l) &
                                 + gaussian_quadrature(integrand_ar_wr_34_35_24_25_ca4321, theta_ar_wr3(l), theta_ar_wr4(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_34_35_24_25_ca4321(theta_ar_wr3(l), theta_ar_wr4(l),l)&
                                 + gaussian_quadrature(integrand_ar_wr_45_ca42, theta_ar_wr4(l), theta_ar_wr5(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_45_ca42(theta_ar_wr4(l), theta_ar_wr5(l),l) &
                                 + gaussian_quadrature(integrand_ar_wr_45pi_ca1234, theta_ar_wr5(l), theta_max, Gauss_nodes, Gauss_weights,l)) / theta_max
                                  
                     sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                                + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)) / theta_max                          
                                  
                  else if (zen(l)>zen_ar_wr4(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l)&
                                  + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l)+ integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)&
                                  + integration_ar_wr_23_ca43(theta_ar_wr2(l), theta_ar_wr3(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_34_35_24_25_ca4321, theta_ar_wr3(l), theta_ar_wr4(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_34_35_24_25_ca4321(theta_ar_wr3(l), theta_ar_wr4(l),l)&
                                  + gaussian_quadrature(integrand_ar_wr_45_ca42, theta_ar_wr4(l), theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_45_ca42(theta_ar_wr4(l), theta_max,l)) / theta_max
                                   
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                                 + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)) / theta_max                          
                                
                  else if (zen(l)>zen_ar_wr3(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l)&
                                  + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l)+ integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)&
                                  + integration_ar_wr_23_ca43(theta_ar_wr2(l), theta_ar_wr3(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_34_35_24_25_ca4321, theta_ar_wr3(l), theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_34_35_24_25_ca4321(theta_ar_wr3(l), theta_max,l)) / theta_max
                                   
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                                 + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)) / theta_max                          
                                       
                  else if (zen(l)>zen_ar_wr2(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l)&
                                  + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l)+ integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)&
                                  + integration_ar_wr_23_ca43(theta_ar_wr2(l), theta_max,l) ) / theta_max
                                   
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                                 + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)) / theta_max                          
                           
                  else if (zen(l)>zen_ar_wr1(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l)&
                                  + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_max, Gauss_nodes, Gauss_weights,l)+ integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_max,l)) / theta_max
                                   
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                                 + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_max,l)) / theta_max                          
                  else 
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_max,l)) / theta_max
                                   
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_max,l)) / theta_max                          
                  end if         
               else if (ca_order(l)==3) then       !zen_ar_wr1<2<3<5<4              
                  if (zen(l)>zen_ar_wr4(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)&
                                  + integration_ar_wr_23_ca43(theta_ar_wr2(l), theta_ar_wr3(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_34_35_24_25_ca4321, theta_ar_wr3(l), theta_ar_wr5(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_34_35_24_25_ca4321(theta_ar_wr3(l), theta_ar_wr5(l),l)&
                                  + gaussian_quadrature(integrand_ar_wr_54_ca31, theta_ar_wr5(l), theta_ar_wr4(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_54_ca31(theta_ar_wr5(l), theta_ar_wr4(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_45pi_ca1234, theta_ar_wr4(l), theta_max, Gauss_nodes, Gauss_weights,l)) / theta_max
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                              + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)) / theta_max                          
                  else if (zen(l)>zen_ar_wr5(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)&
                                  + integration_ar_wr_23_ca43(theta_ar_wr2(l), theta_ar_wr3(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_34_35_24_25_ca4321, theta_ar_wr3(l), theta_ar_wr5(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_34_35_24_25_ca4321(theta_ar_wr3(l), theta_ar_wr5(l),l)&
                                  + gaussian_quadrature(integrand_ar_wr_54_ca31, theta_ar_wr5(l),theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_54_ca31(theta_ar_wr5(l),theta_max,l) ) / theta_max
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                              + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)) / theta_max                                            
                  else if (zen(l)>zen_ar_wr3(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)&
                                  + integration_ar_wr_23_ca43(theta_ar_wr2(l), theta_ar_wr3(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_34_35_24_25_ca4321, theta_ar_wr3(l), theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_34_35_24_25_ca4321(theta_ar_wr3(l), theta_max,l)) / theta_max
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                              + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)) / theta_max                                                                              
                  else if (zen(l)>zen_ar_wr2(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)&
                                  + integration_ar_wr_23_ca43(theta_ar_wr2(l), theta_max,l) ) / theta_max
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                              + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)) / theta_max                                                                                   
                  else if (zen(l)>zen_ar_wr1(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l),theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_max,l)) / theta_max
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                                  + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_max,l)) / theta_max                                                                                                                                      
                  else 
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_max,l) ) / theta_max
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_max,l)) / theta_max                                                                                                                                      
                  
                  end if 
               else if (ca_order(l)==2) then  !zen_ar_wr1<3<2<4<5
                  if (zen(l)>zen_ar_wr5(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr3(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr3(l),l)&
                                  + gaussian_quadrature(integrand_ar_wr_32_ca21, theta_ar_wr3(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l)+integration_ar_wr_32_ca21(theta_ar_wr3(l), theta_ar_wr2(l),l)&
                                  + gaussian_quadrature(integrand_ar_wr_34_35_24_25_ca4321, theta_ar_wr2(l), theta_ar_wr4(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_34_35_24_25_ca4321(theta_ar_wr2(l), theta_ar_wr4(l),l)&
                                  + gaussian_quadrature(integrand_ar_wr_45_ca42, theta_ar_wr4(l), theta_ar_wr5(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_45_ca42(theta_ar_wr4(l), theta_ar_wr5(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_45pi_ca1234, theta_ar_wr5(l), theta_max, Gauss_nodes, Gauss_weights,l)) / theta_max
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                            + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)) / theta_max                               
                  else if (zen(l)>zen_ar_wr4(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr3(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr3(l),l)&
                                   + gaussian_quadrature(integrand_ar_wr_32_ca21, theta_ar_wr3(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l)+integration_ar_wr_32_ca21(theta_ar_wr3(l), theta_ar_wr2(l),l)&
                                   + gaussian_quadrature(integrand_ar_wr_34_35_24_25_ca4321, theta_ar_wr2(l), theta_ar_wr4(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_34_35_24_25_ca4321(theta_ar_wr2(l), theta_ar_wr4(l),l)&
                                   + gaussian_quadrature(integrand_ar_wr_45_ca42, theta_ar_wr4(l), theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_45_ca42(theta_ar_wr4(l), theta_max,l)) / theta_max
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                               + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)) / theta_max                                  
                  else if (zen(l)>zen_ar_wr2(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr3(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr3(l),l)&
                                   + gaussian_quadrature(integrand_ar_wr_32_ca21, theta_ar_wr3(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l)+integration_ar_wr_32_ca21(theta_ar_wr3(l), theta_ar_wr2(l),l)&
                                   + gaussian_quadrature(integrand_ar_wr_34_35_24_25_ca4321, theta_ar_wr2(l), theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_34_35_24_25_ca4321(theta_ar_wr2(l), theta_max,l)) / theta_max
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                               + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)) / theta_max                                                                            
                  else if (zen(l)>zen_ar_wr3(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr3(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr3(l),l)&
                                   + gaussian_quadrature(integrand_ar_wr_32_ca21, theta_ar_wr3(l), theta_max, Gauss_nodes, Gauss_weights,l)+integration_ar_wr_32_ca21(theta_ar_wr3(l), theta_max,l)) / theta_max
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                               + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_max,l)) / theta_max                                                          
                  else if (zen(l)>zen_ar_wr1(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_max,l)) / theta_max
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l)&
                               + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_max,l)) / theta_max                                                                             
                  else 
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_max,l)) / theta_max    
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_max,l)) / theta_max                   
                  end if  
               else if (ca_order(l)==1) then !zen_ar_wr1<3<2<5<4
                  if (zen(l)>zen_ar_wr4(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr3(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr3(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_32_ca21, theta_ar_wr3(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l)+ integration_ar_wr_32_ca21(theta_ar_wr3(l), theta_ar_wr2(l),l) + &
                                  + gaussian_quadrature(integrand_ar_wr_34_35_24_25_ca4321, theta_ar_wr2(l), theta_ar_wr5(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_34_35_24_25_ca4321(theta_ar_wr2(l), theta_ar_wr5(l),l)&
                                  + gaussian_quadrature(integrand_ar_wr_54_ca31, theta_ar_wr5(l), theta_ar_wr4(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_54_ca31(theta_ar_wr5(l), theta_ar_wr4(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_45pi_ca1234, theta_ar_wr4(l), theta_max, Gauss_nodes, Gauss_weights,l)) / theta_max
                                  
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l) &
                                  + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)) / theta_max                                  
                                  
                  else if (zen(l)>zen_ar_wr5(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr3(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr3(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_32_ca21, theta_ar_wr3(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l)+ integration_ar_wr_32_ca21(theta_ar_wr3(l), theta_ar_wr2(l),l) + &
                                   + gaussian_quadrature(integrand_ar_wr_34_35_24_25_ca4321, theta_ar_wr2(l), theta_ar_wr5(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_34_35_24_25_ca4321(theta_ar_wr2(l), theta_ar_wr5(l),l)&
                                   + gaussian_quadrature(integrand_ar_wr_54_ca31, theta_ar_wr5(l), theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_54_ca31(theta_ar_wr5(l), theta_max,l)) / theta_max        
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)) / theta_max                                  
                                   
                  else if (zen(l)>zen_ar_wr2(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr3(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr3(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_32_ca21, theta_ar_wr3(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l)+ integration_ar_wr_32_ca21(theta_ar_wr3(l), theta_ar_wr2(l),l) + &
                                   + gaussian_quadrature(integrand_ar_wr_34_35_24_25_ca4321, theta_ar_wr2(l), theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_34_35_24_25_ca4321(theta_ar_wr2(l), theta_max,l)) / theta_max        
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_ar_wr2(l),l)) / theta_max                                                                          
                  else if (zen(l)>zen_ar_wr3(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr3(l), Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_ar_wr3(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_32_ca21, theta_ar_wr3(l),theta_max , Gauss_nodes, Gauss_weights,l)+ integration_ar_wr_32_ca21(theta_ar_wr3(l),theta_max ,l) ) / theta_max        
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_max , Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_max ,l)) / theta_max                                                                                                          
                  else if (zen(l)>zen_ar_wr1(l)) then
                      sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l),theta_max, Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_13_ca4321(theta_ar_wr1(l), theta_max,l) ) / theta_max        
                      sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_ar_wr1(l),l) &
                                   + gaussian_quadrature(integrand_ar_wr_12_road_ca4321, theta_ar_wr1(l), theta_max , Gauss_nodes, Gauss_weights,l) + integration_ar_wr_12_road_ca4321(theta_ar_wr1(l), theta_max ,l)) / theta_max                                                                                                                                                   
                  else 
                     sdir_sunwall(l,ib) = (integration_ar_wr_01_ca4321(theta_max,l)) / theta_max    
                     sdir_road(l,ib)= (integration_ar_wr_road_01_ca4321(theta_max,l)) / theta_max                   
                  end if  
               end if
            
            end if 
            
            ! calculate the incident direct solar flux without tree using analytical solution
            sdir_shadewall_o(l,ib) = 0._r8
            ! incident solar radiation on wall and road integrated over all canyon orientations (0 <= theta <= pi/2)
            sdir_road_o(l,ib) = 1.0_r8 *                                    &
                 (2._r8*theta0(l)/rpi - 2./rpi*(ht_roof(l)/wcan(l))*tanzen(l)*(1._r8-cos(theta0(l))))              
            sdir_sunwall_o(l,ib) = 2._r8 * 1.0_r8  * ((1._r8/(ht_roof(l)/wcan(l)))* &
                 (0.5_r8-theta0(l)/rpi) + (1._r8/rpi)*tanzen(l)*(1._r8-cos(theta0(l))))
            
            ! calculate the tree attenuation terms for each surface
            if (ca_order(l)==0) then !tree below roofs
               sdir_ar_tree(l,ib)=0.0_r8
               sdir_br_tree(l,ib)=(sdir(l,ib)-sdir_road(l,ib)-(sdir_shadewall(l,ib) + sdir_sunwall(l,ib))* canyon_hwr(l))*wcan(l)/max(A_v1(l),min_zen)
            else 
               sdir_tree_roof(l,ib) = (1.0_r8 - sdir_roof(l,ib))*wbui(l)/max(A_v2(l),min_zen)
               sdir_ar_tree(l,ib)=(sdir(l,ib)-sdir_after_ar(l,ib))*wcan(l)/max(A_v2(l),min_zen)+ sdir_tree_roof(l,ib)
               sdir_br_tree(l,ib)=(sdir_after_ar(l,ib)-sdir_sunwall(l,ib)*canyon_hwr(l)-sdir_shadewall(l,ib)*canyon_hwr(l)-sdir_road(l,ib))*wcan(l)/max(A_v1(l),min_zen)
            end if
            
            !write(6,*) 'sdir_ar_tree(l,ib),sdir_br_tree(l,ib)', sdir_ar_tree(l,ib),sdir_br_tree(l,ib)
            
            ! conservation check for road and wall. need to use wall fluxes converted to ground area
            ! calculate the incident direct solar beam by Masson 2000 for consistency check
            ! Need to add tree attenuation term later
            ! swall_projected_t = (sdir_shadewall(l,ib) + sdir_sunwall(l,ib)) * canyon_hwr(l)

            ! conservation check for road and wall. need to use wall fluxes converted to ground area
            ! Need to add tree attenuation term later
            sdir_tree_sunwall(l,ib) = (sdir_sunwall_o(l,ib) - sdir_sunwall(l,ib))*ht_roof(l)/max(A_v1(l)+A_v2(l),min_zen)
            sdir_tree_road(l,ib) = (sdir_road_o(l,ib) - sdir_road(l,ib))*wcan(l)/max(A_v1(l)+A_v2(l),min_zen)   
            swall_projected_t = (sdir_shadewall(l,ib) + sdir_sunwall(l,ib)) * canyon_hwr(l) + sdir_tree_sunwall(l,ib)*(A_v1(l)+A_v2(l))/wcan(l)
            err1(l) = sdir(l,ib) - (sdir_road(l,ib) + swall_projected_t + sdir_tree_road(l,ib)*(A_v1(l)+A_v2(l))/wcan(l))   
            if (abs(err1(l)) > 0.001_r8) then
                write (6,'(A,2F10.3)') 'swall_projected_t, sdir_road(l,ib)',swall_projected_t, sdir_road(l,ib)
                write (6,'(A,2F10.3)') 'sdir_tree_road(l,ib)*Av1/wcan,err1(l)',sdir_tree_road(l,ib)*(A_v1(l)+A_v2(l))/wcan(l), err1(l)
            endif                                                     
          else
            sdir_road(l,ib) = 0._r8
            sdir_sunwall(l,ib) = 0._r8
            sdir_shadewall(l,ib) = 0._r8
            sdir_roof(l,ib)=0._r8
            sdir_ar_tree(l,ib)=0._r8
            sdir_br_tree(l,ib)=0._r8
            sdir_road_o(l,ib) = 0._r8
            sdir_sunwall_o(l,ib) = 0._r8
            sdir_shadewall_o(l,ib) = 0._r8            
         end if
         

         
         !write (6,'(A,2F10.3)') 'A_v1, A_v2',A_v1(l), A_v2(l)
         !write (6,'(A,2F10.3)') 'sdir_shadewall(l,ib),sdir_sunwall(l,ib)',sdir_shadewall(l,ib),sdir_sunwall(l,ib)
         !write (6,'(A,2F10.3)') 'canyon_hwr(l),sdir_tree_sunwall(l,ib)',canyon_hwr(l),sdir_tree_sunwall(l,ib)
         !write (6,'(A,2F10.3)') 'wcan(l)',wcan(l)
         !write (6,'(A,2F10.3)') '(A_v1(l)+A_v2(l))',(A_v1(l)+A_v2(l))
         

                  
         ! printouts
         if (debug_write_dir) then          
             write (6,'(A,2F10.3)') 'sdir_road(l,ib), sdir_road_o(l,ib) ',sdir_road(l,ib),sdir_road_o(l,ib)
             write (6,'(A,2F10.3)') 'sdir_road(l,ib) - sdir_road_o(l,ib) ',sdir_road(l,ib) - sdir_road_o(l,ib)
             write (6,'(A,2F10.3)') 'sdir_sunwall(l,ib), sdir_sunwall_o(l,ib) ',sdir_sunwall(l,ib), sdir_sunwall_o(l,ib)
             write (6,'(A,2F10.3)') 'sdir_sunwall(l,ib) - sdir_sunwall_o(l,ib)) ',sdir_sunwall(l,ib) - sdir_sunwall_o(l,ib)
             !write (6,'(A,2F10.3)') 'sdir_roof(l,ib) ',sdir_roof(l,ib)
             !write (6,'(A,2F10.3)') 'sdir_roof_shaded(l) ',sdir_roof_shaded(l)
             write (6,'(A,2F10.3)') 'err1(l) ',err1(l)

             !write (6,'(A,2F10.3)') 'integration_ar_wr_01_ca4321 ',integration_ar_wr_01_ca4321(theta_ar_wr1(l),l) 
             !write (6,'(A,2F10.3)') 'integrand_ar_wr_12_13_ca4321',gaussian_quadrature(integrand_ar_wr_12_13_ca4321, theta_ar_wr1(l), theta_ar_wr3(l), Gauss_nodes, Gauss_weights,l)
             !write (6,'(A,2F10.3)') 'integration_ar_wr_12_13__ca4321',integration_ar_wr_12_13__ca4321(theta_ar_wr1(l), theta_ar_wr3(l),l)
             !write (6,'(A,2F10.3)') 'integration_ar_wr_32_ca21',integration_ar_wr_32_ca21(theta_ar_wr3(l), theta_ar_wr2(l),l)
             !write (6,'(A,2F10.3)') 'integrand_ar_wr_32_ca21',gaussian_quadrature(integrand_ar_wr_32_ca21, theta_ar_wr3(l), theta_ar_wr2(l), Gauss_nodes, Gauss_weights,l)
             !write (6,'(A,2F10.3)') 'integrand_ar_wr_34_35_24_25_ca4321',gaussian_quadrature(integrand_ar_wr_34_35_24_25_ca4321, theta_ar_wr2(l), theta_ar_wr5(l), Gauss_nodes, Gauss_weights,l) 
             !write (6,'(A,2F10.3)') 'integration_ar_wr_34_35_24_25__ca4321',integration_ar_wr_34_35_24_25__ca4321(theta_ar_wr2(l), theta_ar_wr5(l),l)
             !write (6,'(A,2F10.3)') 'integrand_ar_wr_54_ca31',gaussian_quadrature(integrand_ar_wr_54_ca31, theta_ar_wr5(l), theta_ar_wr4(l), Gauss_nodes, Gauss_weights,l)
             !write (6,'(A,2F10.3)') 'integration_ar_wr_54_ca31',integration_ar_wr_54_ca31(theta_ar_wr5(l), theta_ar_wr4(l),l) 
             !write (6,'(A,2F10.3)') 'integrand_ar_wr_45pi_ca1234',gaussian_quadrature(integrand_ar_wr_45pi_ca1234, theta_ar_wr4(l), theta_max, Gauss_nodes, Gauss_weights,l) 
             !write (6,'(A,2F10.3)') 'sdir_road(l,ib),sdir_road_o(l,ib) ',sdir_road(l,ib),sdir_road_o(l,ib) 
             !write (6,'(A,2F10.3)') 'sdir_sunwall(l,ib),sdir_sunwall_o(l,ib) ',sdir_sunwall(l,ib),sdir_sunwall_o(l,ib)
             !write (6,'(A,2F10.3)') 'sdir_roof(l,ib)', sdir_roof(l,ib)
         end if
       end do
    end do
    
    contains    
      
    ! The integration functions 
    !------------------------------tree below roof wall functions-----------------------------
    function integration_br_wr_01(theta_in1, l) result(value)
      implicit none
      integer, intent(in) :: l                      ! urban landunit index
      real(r8), intent(in) :: theta_in1
      real(r8)            :: value
      real(r8)            :: term1, term2, term3
      real(r8)            :: exp_m_h2

      exp_m_h2 = exp(-Tree_at(l) * (h2(l) / coszen_min))
      term1 = h1(l) * tanzen(l) * exp_m_h2
      term2 = sinzen(l)/Tree_at(l) * (1.0_r8 - exp_m_h2)
      term3 = -Tree_abr(l) * tanzen(l)
      value = (1.0_r8 / ht_roof(l)) * (term1 + term2 + term3) * (1.0_r8 - cos(theta_in1))
    end function integration_br_wr_01

    function integration_br_wr_12(theta_in1, theta_in2, l) result(value)
        implicit none
        integer, intent(in) :: l                      ! urban landunit index
        real(r8), intent(in) :: theta_in1, theta_in2
        real(r8)            :: value
        real(r8)            :: term1, term2
        real(r8)            :: exp_m_h2

        exp_m_h2 = exp(-Tree_at(l) * (h2(l) / coszen_min))
        term1 = exp_m_h2 * wcan(l) * (theta_in2 - theta_in1)
        term2 = sinzen(l)/Tree_at(l)  * (1.0_r8 - exp_m_h2) &
              - exp_m_h2 * (ht_roof(l) - h1(l)) * tanzen(l) &
              - Tree_abr(l) * tanzen(l)
        value = (1.0_r8 / ht_roof(l)) * (term1 + term2* (cos(theta_in1) - cos(theta_in2)))
    end function integration_br_wr_12

    function integrand_br_wr_23(theta, l) result(value)
        implicit none
        integer, intent(in) :: l                      ! urban landunit index
        real(r8), intent(in) :: theta
        real(r8)            :: value
        real(r8)            :: sin_theta, exp_m_W_h1_h2_H

        sin_theta = sin(theta)
        ! KZ exp_m_W_h1_h2_H = exp(-Tree_at(l) * (wcan(l) / (max(sin(theta),min_zen) * sinzen_min) + Tree_abr(l) / coszen_min))
        exp_m_W_h1_h2_H = exp(-Tree_at(l) * (wcan(l) / (sin(theta) * sinzen(l)) + Tree_abr(l) / coszen(l)))
        value = sin_theta * sinzen(l) / (ht_roof(l) * Tree_at(l)) * (1.0_r8 - exp_m_W_h1_h2_H)
    end function integrand_br_wr_23

    function integration_br_wr_23(theta_in2, theta_in3, l) result(value)
        implicit none
        integer, intent(in)  :: l                      ! urban landunit index
        real(r8), intent(in) :: theta_in2, theta_in3
        real(r8)             :: value

        value = -Tree_abr(l) * tanzen(l) * (cos(theta_in2) - cos(theta_in3)) / ht_roof(l)
    end function integration_br_wr_23
        
    !---------------------------------tree aboveroof wall functions---------------------------------
    function integration_ar_wr_01_ca4321(theta_in1,l) result(value) !checked
       implicit none
       integer, intent(in) :: l                      ! urban landunit index
       real(r8), intent(in) :: theta_in1
       real(r8)            :: value 
       real(r8)            :: exp_m_h2, exp_H_h1_h2, term1, term2
       
       exp_m_h2 = exp(-Tree_at(l) * h2(l) / coszen_min)
       exp_H_h1_h2 = exp(Tree_at(l) * (-Tree_abr(l) / coszen_min))
       term1 = h1(l) * tanzen(l) * exp_m_h2
       term2 = sinzen(l) / Tree_at(l) * (exp_H_h1_h2 - exp_m_h2)
       value = 1.0_r8 / ht_roof(l) * (term1 + term2) * (1-cos(theta_in1))       
    end function integration_ar_wr_01_ca4321
    
    function integrand_ar_wr_12_13_ca4321(theta,l) result(value) !checked
       implicit none
       integer, intent(in) :: l                      ! urban landunit index
       real(r8), intent(in) :: theta
       real(r8)            :: value 
       real(r8)            :: exp_m_w_h1

       ! KZ exp_m_w_h1 = exp(-Tree_at(l) * (wcan(l) / (max(sin(theta),min_zen) * sinzen_min) - h1(l) / coszen_min))
       exp_m_w_h1 = exp(-Tree_at(l) * (wcan(l) / (sin(theta) * sinzen(l)) - h1(l) / coszen(l)))
       value = sin(theta) / ht_roof(l)*sinzen(l)/Tree_at(l)* exp_m_w_h1 
    end function integrand_ar_wr_12_13_ca4321

    function integration_ar_wr_12_13_ca4321(theta_in_s, theta_in_l,l) result(value) !checked
       implicit none
       integer, intent(in) :: l                      ! urban landunit index
       real(r8), intent(in) :: theta_in_s, theta_in_l
       real(r8)            :: value 
       real(r8)             :: exp_m_h2, exp_H_h1_h2, term1, term2, term3
       
       exp_m_h2 = exp(-Tree_at(l) * h2(l) /coszen_min)
       exp_H_h1_h2 = exp(Tree_at(l) * (-Tree_abr(l)/ coszen_min))
       term1 = wcan(l) * exp_m_h2 * (theta_in_l-theta_in_s)
       term2 = (h2(l) * tanzen(l)  + 2.0_r8 * sinzen(l)  / Tree_at(l)) * exp_m_h2 * (cos(theta_in_s) - cos(theta_in_l))
       term3 = sinzen(l)  / Tree_at(l) * exp_H_h1_h2 * (cos(theta_in_s) - cos(theta_in_l))
       value =1.0_r8 / ht_roof(l) * (term1 - term2 + term3)
    end function integration_ar_wr_12_13_ca4321

    function integration_ar_wr_23_ca43(theta_in2, theta_in3,l) result(value)
      implicit none
      integer, intent(in) :: l                      ! urban landunit index
      real(r8), intent(in) :: theta_in2, theta_in3
      real(r8) :: value, exp_m_H_h1, exp_H_h1_h2, exp_m_h2, exp_m_2h2, term1, term2
      
      exp_m_H_h1 = exp(-Tree_at(l) * (ht_roof(l) - h1(l)) / coszen_min)
      exp_H_h1_h2 = exp(Tree_at(l) * (-Tree_abr(l)) /coszen_min)
      exp_m_h2 = exp(-Tree_at(l) * h2(l) / coszen_min )
      term1 = sinzen(l)  / Tree_at(l) * (exp_m_H_h1 + exp_H_h1_h2 - 2 * exp_m_h2) - h2(l) * tanzen(l)  * exp_m_h2
      term2 = wcan(l) * exp_m_h2  * (theta_in3- theta_in2)
      value = 1.0_r8 / ht_roof(l) * (term1 * (cos(theta_in2) - cos(theta_in3)) + term2)
    end function integration_ar_wr_23_ca43

    function integrand_ar_wr_34_35_24_25_ca4321(theta,l) result(value) !checked
      implicit none
      integer, intent(in) :: l                      ! urban landunit index
      real(r8), intent(in) :: theta
      real(r8) :: value, sin_theta, exp_m_W
      
      sin_theta = max(sin(theta),min_zen)
      exp_m_W = exp(-Tree_at(l)* wcan(l) / (sinzen_min   * sin_theta))
      value = 1.0_r8 / ht_roof(l) * (h2(l) * tanzen(l)  - wcan(l) / sin_theta - 2 * sinzen (l) / Tree_at(l)) * exp_m_W * sin_theta
    end function integrand_ar_wr_34_35_24_25_ca4321

    function integration_ar_wr_34_35_24_25_ca4321(theta_in_s, theta_in_l,l) result(value) !checked
      implicit none
      integer, intent(in) :: l                      ! urban landunit index
      real(r8), intent(in) :: theta_in_s, theta_in_l
      real(r8) :: value, exp_m_h2_h1_H, exp_m_H_h1
      
      exp_m_h2_h1_H = exp(-Tree_at(l) * Tree_abr(l) / coszen_min)
      exp_m_H_h1 = exp(-Tree_at(l) * (ht_roof(l) - h1(l)) / coszen_min)
      value = 1.0_r8 / ht_roof(l) * sinzen(l) / Tree_at(l) * (exp_m_h2_h1_H + exp_m_H_h1) *(cos(theta_in_s) - cos(theta_in_l))
    end function integration_ar_wr_34_35_24_25_ca4321
    
    function integrand_ar_wr_45_ca42(theta,l) result(value)
      implicit none
      integer, intent(in) :: l                      ! urban landunit index
      real(r8), intent(in) :: theta
      real(r8) :: value, sin_theta, exp_m_W
      
      sin_theta = max(sin(theta),min_zen)
      exp_m_W = exp(-Tree_at(l) * wcan(l) / (sin_theta * sinzen_min ))
      value = 1.0_r8 / ht_roof(l) * (Tree_abr(l) * tanzen(l) - sinzen(l) / Tree_at(l)) * exp_m_W * sin_theta
    end function integrand_ar_wr_45_ca42

    function integration_ar_wr_45_ca42(theta_in4, theta_in5,l) result(value)
      implicit none
      integer, intent(in) :: l                      ! urban landunit index
      real(r8), intent(in) :: theta_in4, theta_in5
      real(r8) :: value, exp_m_h2_h1_H
      
      exp_m_h2_h1_H = exp(-Tree_at(l) * Tree_abr(l) / coszen_min)
      value = 1.0_r8 / ht_roof(l) * sinzen(l) / Tree_at(l) * exp_m_h2_h1_H * (cos(theta_in4) - cos(theta_in5))
    end function integration_ar_wr_45_ca42

    function integrand_ar_wr_45pi_ca1234(theta,l) result(value) !checked
      implicit none
      integer, intent(in) :: l                      ! urban landunit index
      real(r8), intent(in) :: theta
      real(r8) :: value, exp_m_W
      exp_m_W = exp(-Tree_at(l) * wcan(l) / (sinzen_min *max(sin(theta),min_zen)))
      value = wcan(l)/ ht_roof(l) * exp_m_W
    end function integrand_ar_wr_45pi_ca1234

    function integrand_ar_wr_54_ca31(theta,l) result(value) !checked
      implicit none
      integer, intent(in) :: l                      ! urban landunit index
      real(r8), intent(in) :: theta
      real(r8) :: value, sin_theta, exp_m_W
      
      sin_theta = max(sin(theta),min_zen)
      exp_m_W = exp(-Tree_at(l) * wcan(l) / (sin_theta *sinzen_min ))
      value = sin_theta / ht_roof(l) * exp_m_W *((ht_roof(l)-h1(l))*tanzen(l)-sinzen(l)/Tree_at(l))
    end function integrand_ar_wr_54_ca31
 
    function integration_ar_wr_54_ca31(theta_in5, theta_in4,l) result(value) !checked
      implicit none
      integer, intent(in) :: l                      ! urban landunit index
      real(r8), intent(in) :: theta_in4, theta_in5
      real(r8) :: value, exp_m_H_h1
      
      exp_m_H_h1 = exp(-Tree_at(l) * (ht_roof(l) - h1(l)) / coszen_min)
      value =  sinzen(l)/ht_roof(l)/Tree_at(l)* exp_m_H_h1* (cos(theta_in5) - cos(theta_in4))
    end function integration_ar_wr_54_ca31
    
    function integrand_ar_wr_32_ca21(theta,l) result(value) !checked
      implicit none
      integer, intent(in) :: l                      ! urban landunit index
      real(r8), intent(in) :: theta
      real(r8) :: value, sin_theta, exp_m_W, exp_m_W_h1,term1
      
      sin_theta = max(sin(theta),min_zen)
      exp_m_W = exp(-Tree_at(l)* wcan(l) / (sinzen_min   * sin_theta))
      ! KZ exp_m_W_h1 = exp(-Tree_at(l) * (wcan(l) / (sin_theta * sinzen_min ) - h1(l) /coszen_min))
      exp_m_W_h1 = exp(-Tree_at(l) * (wcan(l) / (sin(theta) * sinzen(l)) - h1(l) /coszen(l)))
      term1=h2(l)*tanzen(l)-wcan(l)/sin_theta-2*sinzen(l)/Tree_at(l)
      value = sin_theta / ht_roof(l) * (term1*exp_m_W+sinzen(l)/Tree_at(l)*exp_m_W_h1)
    end function integrand_ar_wr_32_ca21
            
    function integration_ar_wr_32_ca21(theta_in_s, theta_in_l,l) result(value) !checked
       implicit none
       integer, intent(in) :: l                      ! urban landunit index
       real(r8), intent(in) :: theta_in_s, theta_in_l
       real(r8)            :: value 
       real(r8)             :: exp_m_H_h1_h2
       
       exp_m_H_h1_h2 = exp(-Tree_at(l) * Tree_abr(l) / coszen_min )
       value = sinzen(l)/ht_roof(l)/ Tree_at(l) * exp_m_H_h1_h2 * (cos(theta_in_s) - cos(theta_in_l))
    end function integration_ar_wr_32_ca21

    !----------------------------------road-------------------------------------  
    function integration_br_road(theta_in1,l) result(value)
       implicit none
       integer, intent(in) :: l                      ! urban landunit index
       real(r8), intent(in) :: theta_in1
       real(r8) :: value, exp_m_h2

       exp_m_h2 = exp(-Tree_at(l) * h2(l) / coszen_min)
       value = 1.0_r8 / wcan(l) * (wcan(l) * theta_in1 - ht_roof(l) * tanzen(l) * (1.0_r8 - cos(theta_in1))) * exp_m_h2
    end function integration_br_road
        
    function integration_ar_wr_road_01_ca4321(theta_in1,l) result(value) !checked
       implicit none
       integer, intent(in) :: l                      ! urban landunit index
       real(r8), intent(in) :: theta_in1
       real(r8) :: value, exp_m_H_h1, exp_m_h1_h2_H, exp_m_h2, term1, term2
       
       exp_m_H_h1 = exp(-Tree_at(l) * (ht_roof(l) - h1(l)) / coszen_min)
       exp_m_h1_h2_H = exp(-Tree_at(l) * Tree_abr(l) / coszen_min)
       exp_m_h2 = exp(-Tree_at(l) * h2(l) / coszen_min)
       term1 = exp_m_H_h1 * sinzen(l) / Tree_at(l) * (1.0_r8 - exp_m_h1_h2_H) * (1.0_r8 - cos(theta_in1))
       term2 = (wcan(l) * theta_in1 - (h1(l) + h2(l)) * tanzen(l)  * (1.0_r8 - cos(theta_in1)))* exp_m_h2
       value = (1.0_r8 / wcan(l)) * (term1 + term2)
    end function integration_ar_wr_road_01_ca4321

    function integrand_ar_wr_12_road_ca4321(theta,l) result(value) !check
       implicit none
       integer, intent(in) :: l                      ! urban landunit index
       real(r8), intent(in) :: theta
       real(r8) :: value, sin_theta, exp_m_W_h1
      
       sin_theta = max(sin(theta),min_zen)
       ! KZ exp_m_W_h1 = exp(-Tree_at(l) * (wcan(l) / (sin_theta*sinzen_min) - h1(l) / coszen_min))
       exp_m_W_h1 = exp(-Tree_at(l) * (wcan(l) / (sin(theta)*sinzen(l)) - h1(l) / coszen(l)))
       value = 1.0_r8/ wcan(l) * (-sinzen(l)  / Tree_at(l)) * exp_m_W_h1 * sin_theta
    end function integrand_ar_wr_12_road_ca4321

    function integration_ar_wr_12_road_ca4321(theta_in1, theta_in2,l) result(value) !check
       implicit none
       integer, intent(in) :: l                      ! urban landunit index
       real(r8), intent(in) :: theta_in1, theta_in2
       real(r8) :: value, exp_m_H_h1

       exp_m_H_h1 = exp(-Tree_at(l) * (ht_roof(l) - h1(l)) / coszen_min)
       value = 1.0_r8/ wcan(l)*sinzen(l)  / Tree_at(l) * exp_m_H_h1 * (cos(theta_in1) - cos(theta_in2))
    end function integration_ar_wr_12_road_ca4321

    !---------------------------------- roofs----------------------------------   
    function integration_ar_BW_WB_01_02(theta_in,l) result(value)
        implicit none
        real(r8), intent(in) :: theta_in
        integer, intent(in) :: l
        real(r8) :: value, term1, exp_m_h1_h2_H, term2

        exp_m_h1_h2_H = exp(-Tree_at(l) * Tree_abr(l) /coszen_min)
        term1 = sinzen(l) / Tree_at(l)*(1.0_r8 - exp_m_h1_h2_H) * (1.0_r8 - cos(theta_in))
        term2 = wbui(l) * theta_in - Tree_abr(l) * tanzen(l) * (1.0_r8 - cos(theta_in))
        value = 1.0_r8 / wbui(l) * (term1  + term2)
    end function integration_ar_BW_WB_01_02
    
    function integrand_ar_BW_12(theta, l) result(value)
        implicit none
        integer, intent(in) :: l                      ! urban landunit index
        real(r8), intent(in) :: theta
        real(r8) :: value, sin_theta, term1, exp_m_W
        sin_theta = max(sin(theta),min_zen)

        term1 = (Tree_abr(l) * tanzen(l) - wcan(l) / sin_theta - sinzen(l)  / Tree_at(l))
        exp_m_W = exp(-Tree_at(l) * wcan(l) / (sin_theta * sinzen_min))
        value = sin_theta / wbui(l) * term1 * exp_m_W
    end function integrand_ar_BW_12
    
    function integration_ar_BW_12(theta_in1, theta_in2, l) result(value)
        implicit none
        real(r8), intent(in) :: theta_in1, theta_in2
        integer, intent(in) :: l
        real(r8) :: value, term1, term2, term3

        term1 = sinzen(l) / Tree_at(l) * (cos(theta_in1) - cos(theta_in2))
        term2 = wbui(l) * (theta_in2 - theta_in1)
        term3 = Tree_abr(l) * tanzen(l) * (cos(theta_in1) - cos(theta_in2))
        value = 1.0_r8 / wbui(l) * (term1 + term2 - term3)
    end function integration_ar_BW_12
    
    function integrand_ar_BW_WB_13_23(theta, l) result(value)
        implicit none
        integer, intent(in) :: l                      ! urban landunit index
        real(r8), intent(in) :: theta
        real(r8) :: value, sin_theta, term1, term2, exp_B_h1_h2_H , exp_m_W
        
        sin_theta =max(sin(theta),min_zen)
        ! exp_B_h1_h2_H = exp(Tree_at(l) * (wbui(l) / (sinzen_min  * sin_theta) - Tree_abr(l) / coszen_min))
        exp_B_h1_h2_H = exp(Tree_at(l) * (wbui(l) / (sinzen(l)  * sin(theta)) - Tree_abr(l) / coszen(l)))
        exp_m_W = exp(-Tree_at(l) * wcan(l) / (sinzen_min  * sin_theta))
        term1 = (Tree_abr(l) * tanzen(l) - wcan(l) / sin_theta) * exp_m_W
        term2 = sinzen(l) / Tree_at(l) * (exp_B_h1_h2_H - exp_m_W)
        value = sin_theta / wbui(l) * (term1+term2)
    end function integrand_ar_BW_WB_13_23

    function integrand_ar_WB_BW_3pi(theta, l) result(value)
        implicit none
        integer, intent(in) :: l                      ! urban landunit index
        real(r8), intent(in) :: theta
        real(r8) :: value !exp_m_W
        value  = exp(-Tree_at(l) * wcan(l) / (sinzen_min  * max(sin(theta),min_zen)))
    end function integrand_ar_WB_BW_3pi
        
    function integrand_ar_WB_21(theta, l) result(value)
        implicit none
        integer, intent(in) :: l                      ! urban landunit index
        real(r8), intent(in) :: theta
        real(r8) :: value, sin_theta, exp_B_h1_h2_H,exp_m_h1_h2_H
        
        !sin_theta = max(sin(theta),min_zen)
        sin_theta = sin(theta)
        if (debug_write_dir) then 
            write (6,'(A,2F10.3)') 'Tree_at(l),Tree_abr(l)',Tree_at(l),Tree_abr(l)
            write (6,'(A,2F10.3)') 'sinzen(l)',sinzen(l)
            write (6,'(A,2F10.3)') 'sin(theta)* sinzen(l)',sin(theta)* sinzen(l)
            write (6,'(A,2F10.3)') 'Tree_abr(l)/coszen_min',Tree_abr(l)/coszen_min
        end if
        exp_m_h1_h2_H = exp(-Tree_at(l) * Tree_abr(l)/coszen_min)
        if (debug_write_dir) then 
            write (6,'(A,2F10.3)') 'exp_m_h1_h2_H',Tree_abr(l)/coszen_min
        end if        
        exp_B_h1_h2_H = exp(Tree_at(l) * (wbui(l) / (sin(theta) * sinzen(l)/coszen(l)*coszen_min)-Tree_abr(l)/coszen_min))
        !exp_B_h1_h2_H = exp(Tree_at(l) * (wbui(l) / (sin_theta * sinzen(l))-Tree_abr(l)/coszen(l)))
        !exp_B_h1_h2_H = exp(Tree_at(l) * (wbui(l) / (sin_theta * sinzen_min)-Tree_abr(l)/coszen_min))
        !exp_B_h1_h2_H = exp(Tree_at(l) * (wbui(l) / max(sin(theta)* sinzen(l),min_zen)-Tree_abr(l)/coszen_min))
        value= sin_theta / wbui(l) *sinzen(l)/Tree_at(l)*(exp_B_h1_h2_H-exp_m_h1_h2_H)
    end function integrand_ar_WB_21    
    !---------------------------------partition incident solar---------------------------------
    function tree_partition_05(theta_in, l) result(value)
        implicit none
        real(r8), intent(in) :: theta_in
        integer, intent(in) :: l
        real(r8) :: value, term1, term2, term3, exp_m_h1_h2_H 
        exp_m_h1_h2_H = exp(-Tree_at(l) * Tree_abr(l)/coszen_min)
        
        term1 = sinzen(l)/Tree_at(l)*(1.0_r8 - exp_m_h1_h2_H) * (1.0_r8 - cos(theta_in))
        term2 = Tree_abr(l)*tanzen(l)*exp_m_h1_h2_H * (1.0_r8 - cos(theta_in))
        term3 = wcan(l)*exp_m_h1_h2_H*theta_in
        value = 1.0_r8 / wcan(l) * (term1 - term2 + term3)
    end function tree_partition_05  

    function integrand_tree_partition_5pi(theta, l) result(value)
        implicit none
        integer, intent(in) :: l                      ! urban landunit index
        real(r8), intent(in) :: theta
        real(r8) :: value,  exp_m_W!,sin_theta
        
        !sin_theta = max(sin(theta),min_zen)
        exp_m_W = exp(-Tree_at(l) * wcan(l) / (sin(theta) * sinzen(l)))
        value = sin(theta)/ wcan(l) *  sinzen(l)  / Tree_at(l) * (1.0_r8 - exp_m_W)
    end function integrand_tree_partition_5pi
   
  end subroutine incident_direct 
  !-----------------------------------------------------------------------
  subroutine incident_diffuse (bounds, &
       num_urbanl, filter_urbanl, canyon_hwr,ht_roof,A_v1,A_v2,wtlunit_roof, &
       sdif, sdif_road, sdif_sunwall, sdif_shadewall,sdif_roof, sdif_ar_tree, sdif_br_tree, &
       urbanparams_inst,solarabs_inst)
    !
    ! !DESCRIPTION: 
    ! Diffuse solar radiation incident on walls and road in urban canyon with tree 
    ! Conservation check: Total incoming diffuse 
    ! (sdif) = sdif_road + (sdif_shadewall + sdif_sunwall)*canyon_hwr
    ! Multiplication by canyon_hwr scales wall fluxes (per unit wall area) to per unit ground area
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)  :: bounds                      
    integer               , intent(in)  :: num_urbanl                           ! number of urban landunits
    integer               , intent(in)  :: filter_urbanl(:)                     ! urban landunit filter
    real(r8)              , intent(in)  :: canyon_hwr     ( bounds%begl: )      ! ratio of building height to street width [landunit]
    real(r8)              , intent(in)  :: sdif           ( bounds%begl: , 1: ) ! diffuse solar radiation incident on horizontal surface [landunit, numrad]
    real(r8)              , intent(out) :: sdif_road      ( bounds%begl: , 1: ) ! diffuse solar radiation incident on road [landunit, numrad]
    real(r8)              , intent(out) :: sdif_sunwall   ( bounds%begl: , 1: ) ! diffuse solar radiation (per unit wall area) incident on sunlit wall [landunit, numrad]
    real(r8)              , intent(out) :: sdif_shadewall ( bounds%begl: , 1: ) ! diffuse solar radiation (per unit wall area) incident on shaded wall [landunit, numrad]
    real(r8)              , intent(out) :: sdif_roof      ( bounds%begl: , 1: ) ! diffuse solar radiation incident on roof [landunit, numrad]
    real(r8)              , intent(out) :: sdif_ar_tree      ( bounds%begl: , 1: ) ! diffuse solar radiation incident on tree above roof [landunit, numrad]
    real(r8)              , intent(out) :: sdif_br_tree      ( bounds%begl: , 1: ) ! diffuse solar radiation incident on tree below roof [landunit, numrad]
    real(r8)              , intent(in)  :: ht_roof     ( bounds%begl: )      ! ratio of building height to street width [landunit]
    real(r8)              , intent(in)  :: wtlunit_roof     ( bounds%begl: )      ! ratio of building height to street width [landunit]
    real(r8)              ,intent(in)    :: A_v1               ( bounds%begl: )             ! height of urban roof (m) [landunit]
    real(r8)              ,intent(in)    :: A_v2               ( bounds%begl: )             ! height of urban roof (m) [landunit]

    type(urbanparams_type), intent(in)  :: urbanparams_inst
    type(solarabs_type)   , intent(in) :: solarabs_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: l, fl, ib       ! indices      
    real(r8) :: err(bounds%begl:bounds%endl)    ! energy conservation error (W/m**2)
    real(r8) :: scan_projected ! projected diffuse solar radiation (per unit ground area) for the whole stret canyon (W/m**2)
    real(r8)            :: A_s                           ! Area of the street canyon (normalized)
    real(r8)            :: A_g                           ! Aarea of the ground (normalized)
    real(r8)            :: A_r                           ! Area of the roof
    real(r8)            :: A_w                           ! Area of the wall (normalized)
    logical  :: debug_write = .false.                  ! true => write out many intermediate variables for debugging

    !-----------------------------------------------------------------------
    
    
    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(canyon_hwr)     == (/bounds%endl/)),         sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sdif)           == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sdif_road)      == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sdif_sunwall)   == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sdif_shadewall) == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sdif_roof)      == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sdif_ar_tree)   == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sdif_br_tree) == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(ht_roof)     == (/bounds%endl/)),         sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(wtlunit_roof)     == (/bounds%endl/)),         sourcefile, __LINE__)

    associate(                            & 
         ksg1d =>    urbanparams_inst%ksg1d_out ,  & ! Input:  [real(r8) (:) ]  Monte carlo view factor of sky for ground[landunit]               
         ksw1d =>    urbanparams_inst%ksw1d_out , & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor of sky for wall[landunit, nzcanm]                      
         ksr1d =>    urbanparams_inst%ksr1d_out , & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor of sky for roof [landunit, nzcanm]                   
         ksv1d =>    urbanparams_inst%ksv1d_out  & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor of sky for vegetation [landunit, nzcanm]                 
         )

      do ib = 1, numrad
         
         ! diffuse solar and conservation check. need to convert wall fluxes to ground area
         
         do fl = 1,num_urbanl
            l = filter_urbanl(fl)
            
            if (l==3) then   
                debug_write=.false. 
            else 
                debug_write=.false. 
            end if
                        
            A_w=ht_roof(l)
            A_r=ht_roof(l)/(canyon_hwr(l)*(1._r8-wtlunit_roof(l))/wtlunit_roof(l))
            A_g=ht_roof(l)/canyon_hwr(l)
            A_s=A_r+A_g
   
            sdif_road(l,ib)      = sdif(l,ib) * ksg1d(l) !ksg1d=vfst*A_s/A_g slightly different from CLM; flux in respect to road area
            sdif_sunwall(l,ib)   = sdif(l,ib) * ksw1d(l,1) !ksw1d(izcan)=vfsw(izcan)*A_s/A_w_max(izcan)
            sdif_shadewall(l,ib) = sdif(l,ib) * ksw1d(l,1) 
            ! ksr1d(izcan)=vfsr(izcan)*A_s/A_r_max(izcan)
            ! ksv1d(izcan)=vfsv(izcan)*A_s/A_vs_max(izcan)
            sdif_roof(l,ib)      = sdif(l,ib) * ksr1d(l,2)  ! the roof is on the 2nd layer
            sdif_ar_tree(l,ib)   = sdif(l,ib) * ksv1d(l,2) 
            sdif_br_tree(l,ib)   = sdif(l,ib) * ksv1d(l,1) 
            scan_projected=sdif_road(l,ib)*A_g/A_s+(sdif_shadewall(l,ib) + sdif_sunwall(l,ib))* A_w/A_s &
                                 +sdif_roof(l,ib)*A_r/A_s + sdif_br_tree(l,ib)*A_v1(l)/A_s+ sdif_ar_tree(l,ib)*A_v2(l)/A_s

            err(l) = sdif(l,ib) - scan_projected
            if (debug_write) then
               write(6,*) '----------------incoming longwave radiation------------ '
               write(6,*) 'sdif(l,ib)= ', sdif(l,ib)
               write(6,*) 'ksg1d(l),ksw1d(l,1),ksr1d(l,2),ksv1d(l,1),ksv1d(l,2)', ksg1d(l),ksw1d(l,1),ksr1d(l,2),ksv1d(l,1),ksv1d(l,2)
               write(6,*) 'A_s(l),A_g(l), A_w(l), A_v1(l), A_v2(l), A_r(l) = ', A_s,A_g, A_w, A_v1(l), A_v2(l), A_r
            end if                    
         end do

         ! error check

         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            if (abs(err(l)) > 0.001_r8) then
               write (iulog,*) 'urban diffuse solar radiation balance error',err(l) 
               write (iulog,*) 'clm model is stopping'
               call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
            endif
         end do

      end do

    end associate

  end subroutine incident_diffuse
    
  !-----------------------------------------------------------------------
  subroutine net_solar (bounds                                                                 , &
       num_urbanl, filter_urbanl, coszen, canyon_hwr, ht_roof, A_v1,A_v2,wtlunit_roof, wtroad_perv,wtroad_tree, sdir, sdif                  , &
       alb_improad_dir, alb_perroad_dir, alb_wall_dir, alb_roof_dir,alb_br_tree_dir,alb_ar_tree_dir  , &
       alb_improad_dif, alb_perroad_dif, alb_wall_dif, alb_roof_dif,alb_br_tree_dif,alb_ar_tree_dif  , &
       sdir_road, sdir_sunwall, sdir_shadewall,sdir_roof,                                                  &
       sdif_road, sdif_sunwall, sdif_shadewall,sdif_roof,                                                    &
       sdir_br_tree, sdir_ar_tree,                                                  &
       sdif_br_tree, sdif_ar_tree,                                                  &         
       urbanparams_inst, solarabs_inst) 

    !
    ! !DESCRIPTION: 
    ! Solar radiation absorbed by road and both walls in urban canyon allowing 
    ! for multiple reflection.
    !
    ! !ARGUMENTS:
    type (bounds_type), intent(in) :: bounds                            
    integer , intent(in)    :: num_urbanl                               ! number of urban landunits
    integer , intent(in)    :: filter_urbanl(:)                         ! urban landunit filter
    real(r8), intent(in)    :: coszen             ( bounds%begl: )      ! cosine solar zenith angle [landunit]
    real(r8), intent(in)    :: canyon_hwr         ( bounds%begl: )      ! ratio of building height to street width [landunit]
    real(r8), intent(in)    :: wtlunit_roof         ( bounds%begl: )        ! weight of roof with respect to landunit  [landunit]
    real(r8), intent(in)    :: ht_roof               ( bounds%begl: )             ! height of urban roof (m) [landunit]
    real(r8), intent(in)    :: A_v1               ( bounds%begl: )             ! height of urban roof (m) [landunit]
    real(r8), intent(in)    :: A_v2               ( bounds%begl: )             ! height of urban roof (m) [landunit]
    real(r8), intent(in)    :: wtroad_perv        ( bounds%begl: )      ! weight of pervious road wrt total road [landunit]
    real(r8), intent(in)    :: wtroad_tree        ( bounds%begl: )      ! weight of road tree wrt total road [landunit]
    real(r8), intent(in)    :: sdir               ( bounds%begl: , 1: ) ! direct beam solar radiation incident on horizontal surface [landunit, numrad]
    real(r8), intent(in)    :: sdif               ( bounds%begl: , 1: ) ! diffuse solar radiation on horizontal surface [landunit, numrad]
    real(r8), intent(in)    :: alb_improad_dir    ( bounds%begl: , 1: ) ! direct impervious road albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_perroad_dir    ( bounds%begl: , 1: ) ! direct pervious road albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_wall_dir       ( bounds%begl: , 1: ) ! direct  wall albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_roof_dir       ( bounds%begl: , 1: ) ! direct  roof albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_improad_dif    ( bounds%begl: , 1: ) ! diffuse impervious road albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_perroad_dif    ( bounds%begl: , 1: ) ! diffuse pervious road albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_wall_dif       ( bounds%begl: , 1: ) ! diffuse wall albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_roof_dif       ( bounds%begl: , 1: ) ! diffuse roof albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_br_tree_dir    ( bounds%begl: , 1: ) ! diffuse impervious road albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_br_tree_dif    ( bounds%begl: , 1: ) ! diffuse pervious road albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_ar_tree_dir    ( bounds%begl: , 1: ) ! diffuse impervious road albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_ar_tree_dif    ( bounds%begl: , 1: ) ! diffuse pervious road albedo [landunit, numrad]
    real(r8), intent(in)    :: sdir_road          ( bounds%begl: , 1: ) ! direct beam solar radiation incident on road per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdir_sunwall       ( bounds%begl: , 1: ) ! direct beam solar radiation (per unit wall area) incident on sunlit wall per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdir_shadewall     ( bounds%begl: , 1: ) ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdif_road          ( bounds%begl: , 1: ) ! diffuse solar radiation incident on road per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdif_sunwall       ( bounds%begl: , 1: ) ! diffuse solar radiation (per unit wall area) incident on sunlit wall per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdif_shadewall     ( bounds%begl: , 1: ) ! diffuse solar radiation (per unit wall area) incident on shaded wall per unit incident flux [landunit, numrad]

    real(r8), intent(in)    :: sdir_br_tree       ( bounds%begl: , 1: ) ! direct solar radiation (per unit wall area) incident on below-roof tree per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdif_br_tree       ( bounds%begl: , 1: ) ! diffuse solar radiation (per unit wall area) incident on below-roof tree per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdir_ar_tree       ( bounds%begl: , 1: ) ! direct solar radiation (per unit wall area) incident on above-roof tree per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdif_ar_tree      ( bounds%begl: , 1: ) ! diffuse solar radiation (per unit wall area) incident on above-roof tree per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdir_roof       ( bounds%begl: , 1: ) ! direct solar radiation (per unit wall area) incident on roof per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdif_roof       ( bounds%begl: , 1: ) ! diffuse solar radiation (per unit wall area) incident on roof per unit incident flux [landunit, numrad] 

    type(urbanparams_type), intent(in)    :: urbanparams_inst
    type(solarabs_type)   , intent(inout) :: solarabs_inst
    !
    ! !LOCAL VARIABLES
    real(r8) :: wtroad_imperv(bounds%begl:bounds%endl)           ! weight of impervious road wrt total road
    
    real(r8) :: improad_a_dir(bounds%begl:bounds%endl)           ! absorbed direct solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_a_dif(bounds%begl:bounds%endl)           ! absorbed diffuse solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_r_dir(bounds%begl:bounds%endl)           ! reflected direct solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_r_dif(bounds%begl:bounds%endl)           ! reflected diffuse solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_r_sky_dir(bounds%begl:bounds%endl)       ! improad_r_dir to sky per unit incident flux
    real(r8) :: improad_r_sunwall_dir(bounds%begl:bounds%endl)   ! improad_r_dir to sunlit wall per unit incident flux
    real(r8) :: improad_r_shadewall_dir(bounds%begl:bounds%endl) ! improad_r_dir to shaded wall per unit incident flux
    real(r8) :: improad_r_br_tree_dir(bounds%begl:bounds%endl) ! improad_r_dir to below-roof tree per unit incident flux
    real(r8) :: improad_r_ar_tree_dir(bounds%begl:bounds%endl) ! improad_r_dir to above-roof tree per unit incident flux
    real(r8) :: improad_r_sky_dif(bounds%begl:bounds%endl)       ! improad_r_dif to sky per unit incident flux
    real(r8) :: improad_r_sunwall_dif(bounds%begl:bounds%endl)   ! improad_r_dif to sunlit wall per unit incident flux
    real(r8) :: improad_r_shadewall_dif(bounds%begl:bounds%endl) ! improad_r_dif to shaded wall per unit incident flux
    real(r8) :: improad_r_br_tree_dif(bounds%begl:bounds%endl) ! improad_r_dif to below-roof tree per unit incident flux
    real(r8) :: improad_r_ar_tree_dif(bounds%begl:bounds%endl) ! improad_r_dif to above-roof tree per unit incident flux

    real(r8) :: perroad_a_dir(bounds%begl:bounds%endl)           ! absorbed direct solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_a_dif(bounds%begl:bounds%endl)           ! absorbed diffuse solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_r_dir(bounds%begl:bounds%endl)           ! reflected direct solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_r_dif(bounds%begl:bounds%endl)           ! reflected diffuse solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_r_sky_dir(bounds%begl:bounds%endl)       ! perroad_r_dir to sky per unit incident flux 
    real(r8) :: perroad_r_sunwall_dir(bounds%begl:bounds%endl)   ! perroad_r_dir to sunlit wall per unit incident flux
    real(r8) :: perroad_r_shadewall_dir(bounds%begl:bounds%endl) ! perroad_r_dir to shaded wall per unit incident flux
    real(r8) :: perroad_r_br_tree_dir(bounds%begl:bounds%endl) ! perroad_r_dir to tree per unit incident flux
    real(r8) :: perroad_r_ar_tree_dir(bounds%begl:bounds%endl) ! perroad_r_dir to above-roof tree per unit incident flux
    real(r8) :: perroad_r_sky_dif(bounds%begl:bounds%endl)       ! perroad_r_dif to sky per unit incident flux
    real(r8) :: perroad_r_sunwall_dif(bounds%begl:bounds%endl)   ! perroad_r_dif to sunlit wall per unit incident flux
    real(r8) :: perroad_r_shadewall_dif(bounds%begl:bounds%endl) ! perroad_r_dif to shaded wall per unit incident flux
    real(r8) :: perroad_r_br_tree_dif(bounds%begl:bounds%endl) ! perroad_r_dif to below-roof tree per unit incident flux
    real(r8) :: perroad_r_ar_tree_dif(bounds%begl:bounds%endl) ! perroad_r_dif to above-roof tree per unit incident flux

    real(r8) :: road_a_dir(bounds%begl:bounds%endl)              ! absorbed direct solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_a_dif(bounds%begl:bounds%endl)              ! absorbed diffuse solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_r_dir(bounds%begl:bounds%endl)              ! reflected direct solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_r_dif(bounds%begl:bounds%endl)              ! reflected diffuse solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_r_sky_dir(bounds%begl:bounds%endl)          ! road_r_dir to sky per unit incident flux
    real(r8) :: road_r_sunwall_dir(bounds%begl:bounds%endl)      ! road_r_dir to sunlit wall per unit incident flux
    real(r8) :: road_r_shadewall_dir(bounds%begl:bounds%endl)    ! road_r_dir to shaded wall per unit incident flux
    real(r8) :: road_r_br_tree_dir(bounds%begl:bounds%endl) ! road_r_dir to below-roof tree per unit incident flux
    real(r8) :: road_r_ar_tree_dir(bounds%begl:bounds%endl) ! road_r_dir to above-roof tree per unit incident flux
    real(r8) :: road_r_sky_dif(bounds%begl:bounds%endl)          ! road_r_dif to sky per unit incident flux
    real(r8) :: road_r_sunwall_dif(bounds%begl:bounds%endl)      ! road_r_dif to sunlit wall per unit incident flux
    real(r8) :: road_r_shadewall_dif(bounds%begl:bounds%endl)    ! road_r_dif to shaded wall per unit incident flux
    real(r8) :: road_r_br_tree_dif(bounds%begl:bounds%endl) ! road_r_dif to below-roof tree per unit incident flux
    real(r8) :: road_r_ar_tree_dif(bounds%begl:bounds%endl) ! road_r_dif to above-roof tree per unit incident flux

    real(r8) :: sunwall_a_dir(bounds%begl:bounds%endl)           ! absorbed direct solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_a_dif(bounds%begl:bounds%endl)           ! absorbed diffuse solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_r_dir(bounds%begl:bounds%endl)           ! reflected direct solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_r_dif(bounds%begl:bounds%endl)           ! reflected diffuse solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_r_sky_dir(bounds%begl:bounds%endl)       ! sunwall_r_dir to sky per unit incident flux
    real(r8) :: sunwall_r_road_dir(bounds%begl:bounds%endl)      ! sunwall_r_dir to road per unit incident flux
    real(r8) :: sunwall_r_shadewall_dir(bounds%begl:bounds%endl) ! sunwall_r_dir to opposing (shaded) wall per unit incident flux
    real(r8) :: sunwall_r_br_tree_dir(bounds%begl:bounds%endl) ! sunwall_r_dir to below-roof tree per unit incident flux
    real(r8) :: sunwall_r_ar_tree_dir(bounds%begl:bounds%endl) ! sunwall_r_dir to above-roof tree per unit incident flux
    real(r8) :: sunwall_r_sky_dif(bounds%begl:bounds%endl)       ! sunwall_r_dif to sky per unit incident flux
    real(r8) :: sunwall_r_road_dif(bounds%begl:bounds%endl)      ! sunwall_r_dif to road per unit incident flux
    real(r8) :: sunwall_r_shadewall_dif(bounds%begl:bounds%endl) ! sunwall_r_dif to opposing (shaded) wall per unit incident flux
    real(r8) :: sunwall_r_br_tree_dif(bounds%begl:bounds%endl) ! sunwall_r_dif to below-roof tree per unit incident flux
    real(r8) :: sunwall_r_ar_tree_dif(bounds%begl:bounds%endl) ! sunwall_r_dif to above-roof tree per unit incident flux

    real(r8) :: shadewall_a_dir(bounds%begl:bounds%endl)         ! absorbed direct solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_a_dif(bounds%begl:bounds%endl)         ! absorbed diffuse solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_r_dir(bounds%begl:bounds%endl)         ! reflected direct solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_r_dif(bounds%begl:bounds%endl)         ! reflected diffuse solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_r_sky_dir(bounds%begl:bounds%endl)     ! shadewall_r_dir to sky per unit incident flux
    real(r8) :: shadewall_r_road_dir(bounds%begl:bounds%endl)    ! shadewall_r_dir to road per unit incident flux
    real(r8) :: shadewall_r_sunwall_dir(bounds%begl:bounds%endl) ! shadewall_r_dir to opposing (sunlit) wall per unit incident flux
    real(r8) :: shadewall_r_br_tree_dir(bounds%begl:bounds%endl) ! shadewall_r_dir to below-roof tree per unit incident flux
    real(r8) :: shadewall_r_ar_tree_dir(bounds%begl:bounds%endl) ! shadewall_r_dir to above-roof tree per unit incident flux
    real(r8) :: shadewall_r_sky_dif(bounds%begl:bounds%endl)     ! shadewall_r_dif to sky per unit incident flux
    real(r8) :: shadewall_r_road_dif(bounds%begl:bounds%endl)    ! shadewall_r_dif to road per unit incident flux
    real(r8) :: shadewall_r_sunwall_dif(bounds%begl:bounds%endl) ! shadewall_r_dif to opposing (sunlit) wall per unit incident flux
    real(r8) :: shadewall_r_br_tree_dif(bounds%begl:bounds%endl) ! shadewall_r_dif to below-roof tree per unit incident flux
    real(r8) :: shadewall_r_ar_tree_dif(bounds%begl:bounds%endl) ! shadewall_r_dif to above-roof tree per unit incident flux

    real(r8) :: ar_tree_a_dir(bounds%begl:bounds%endl)         ! absorbed direct solar for ar_tree (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: ar_tree_a_dif(bounds%begl:bounds%endl)         ! absorbed diffuse solar for ar_tree (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: ar_tree_r_dir(bounds%begl:bounds%endl)         ! reflected direct solar for ar_tree (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: ar_tree_r_dif(bounds%begl:bounds%endl)         ! reflected diffuse solar for ar_tree (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: ar_tree_r_sky_dir(bounds%begl:bounds%endl)     ! ar_tree_r_dir to sky per unit incident flux
    real(r8) :: ar_tree_r_road_dir(bounds%begl:bounds%endl)    ! ar_tree_r_dir to road per unit incident flux
    real(r8) :: ar_tree_r_sunwall_dir(bounds%begl:bounds%endl) ! ar_tree_r_dir to opposing (sunlit) wall per unit incident flux
    real(r8) :: ar_tree_r_shadewall_dir(bounds%begl:bounds%endl) ! ar_tree_r_dir to opposing (shaded) wall per unit incident flux
    real(r8) :: ar_tree_r_br_tree_dir(bounds%begl:bounds%endl) ! tree_r_dir to below-roof tree per unit incident flux
    real(r8) :: ar_tree_r_ar_tree_dir(bounds%begl:bounds%endl) ! tree_r_dir to above-roof tree per unit incident flux
    real(r8) :: ar_tree_r_roof_dir(bounds%begl:bounds%endl) ! tree_r_dir to above-roof tree per unit incident flux

    real(r8) :: ar_tree_r_sky_dif(bounds%begl:bounds%endl)     ! ar_tree_r_dif to sky per unit incident flux
    real(r8) :: ar_tree_r_road_dif(bounds%begl:bounds%endl)    ! ar_tree_r_dif to road per unit incident flux
    real(r8) :: ar_tree_r_sunwall_dif(bounds%begl:bounds%endl) ! ar_tree_r_dif to opposing (sunlit) wall per unit incident flux
    real(r8) :: ar_tree_r_shadewall_dif(bounds%begl:bounds%endl) ! ar_tree_r_dif to opposing (shaded) wall per unit incident flux
    real(r8) :: ar_tree_r_br_tree_dif(bounds%begl:bounds%endl) ! tree_r_dif to below-roof tree per unit incident flux
    real(r8) :: ar_tree_r_ar_tree_dif(bounds%begl:bounds%endl) ! tree_r_dif to above-roof tree per unit incident flux
    real(r8) :: ar_tree_r_roof_dif(bounds%begl:bounds%endl) ! tree_r_dif to above-roof tree per unit incident flux

    real(r8) :: br_tree_a_dir(bounds%begl:bounds%endl)         ! absorbed direct solar for br_tree (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: br_tree_a_dif(bounds%begl:bounds%endl)         ! absorbed diffuse solar for br_tree (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: br_tree_r_dir(bounds%begl:bounds%endl)         ! reflected direct solar for br_tree (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: br_tree_r_dif(bounds%begl:bounds%endl)         ! reflected diffuse solar for br_tree (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: br_tree_r_sky_dir(bounds%begl:bounds%endl)     ! br_tree_r_dir to sky per unit incident flux
    real(r8) :: br_tree_r_road_dir(bounds%begl:bounds%endl)    ! br_tree_r_dir to road per unit incident flux
    real(r8) :: br_tree_r_sunwall_dir(bounds%begl:bounds%endl) ! br_tree_r_dir to opposing (sunlit) wall per unit incident flux
    real(r8) :: br_tree_r_shadewall_dir(bounds%begl:bounds%endl) ! br_tree_r_dir to opposing (shaded) wall per unit incident flux
    real(r8) :: br_tree_r_br_tree_dir(bounds%begl:bounds%endl) ! tree_r_dir to below-roof tree per unit incident flux
    real(r8) :: br_tree_r_ar_tree_dir(bounds%begl:bounds%endl) ! tree_r_dir to above-roof tree per unit incident flux
    real(r8) :: br_tree_r_roof_dir(bounds%begl:bounds%endl) ! tree_r_dir to above-roof tree per unit incident flux

    real(r8) :: br_tree_r_sky_dif(bounds%begl:bounds%endl)     ! br_tree_r_dif to sky per unit incident flux
    real(r8) :: br_tree_r_road_dif(bounds%begl:bounds%endl)    ! br_tree_r_dif to road per unit incident flux
    real(r8) :: br_tree_r_sunwall_dif(bounds%begl:bounds%endl) ! br_tree_r_dif to opposing (sunlit) wall per unit incident flux
    real(r8) :: br_tree_r_shadewall_dif(bounds%begl:bounds%endl) ! br_tree_r_dif to opposing (shaded) wall per unit incident flux
    real(r8) :: br_tree_r_br_tree_dif(bounds%begl:bounds%endl) ! tree_r_dif to below-roof tree per unit incident flux
    real(r8) :: br_tree_r_ar_tree_dif(bounds%begl:bounds%endl) ! tree_r_dif to above-roof tree per unit incident flux
    real(r8) :: br_tree_r_roof_dif(bounds%begl:bounds%endl) ! tree_r_dif to above-roof tree per unit incident flux

    real(r8) :: roof_a_dir(bounds%begl:bounds%endl)              ! absorbed direct solar for total roof after "n" reflections per unit incident flux
    real(r8) :: roof_a_dif(bounds%begl:bounds%endl)              ! absorbed diffuse solar for total roof after "n" reflections per unit incident flux
    real(r8) :: roof_r_dir(bounds%begl:bounds%endl)              ! reflected direct solar for total roof after "n" reflections per unit incident flux
    real(r8) :: roof_r_dif(bounds%begl:bounds%endl)              ! reflected diffuse solar for total roof after "n" reflections per unit incident flux
    real(r8) :: roof_r_sky_dir(bounds%begl:bounds%endl)          ! roof_r_dir to sky per unit incident flux
    real(r8) :: roof_r_ar_tree_dir(bounds%begl:bounds%endl) ! roof_r_dir to above-roof tree per unit incident flux
    real(r8) :: roof_r_sky_dif(bounds%begl:bounds%endl)          ! roof_r_dif to sky per unit incident flux
    real(r8) :: roof_r_ar_tree_dif(bounds%begl:bounds%endl) ! roof_r_dif to above-roof tree per unit incident flux

    real(r8) :: canyon_alb_dir(bounds%begl:bounds%endl)          ! direct canyon albedo
    real(r8) :: canyon_alb_dif(bounds%begl:bounds%endl)          ! diffuse canyon albedo

    real(r8) :: stot(bounds%begl:bounds%endl)                    ! sum of radiative terms
    real(r8) :: stot_dir(bounds%begl:bounds%endl)                ! sum of direct radiative terms
    real(r8) :: stot_dif(bounds%begl:bounds%endl)                ! sum of diffuse radiative terms
    !real(r8) :: alb_leaf_dir(bounds%begl:bounds%endl,numrad)                ! direct leaf albedo [landunit, numrad]
    !real(r8) :: alb_leaf_dif(bounds%begl:bounds%endl,numrad)                ! diffuse leaf albedo [landunit, numrad]

    integer  :: l,fl,ib                          ! indices
    integer  :: iter_dir,iter_dif                ! iteration counter
    real(r8) :: crit                             ! convergence criterion
    real(r8) :: err                              ! energy conservation error
    real(r8) :: err_r                          ! energy conservation error (W/m**2)

    integer  :: pass
    integer, parameter :: n = 50                 ! number of interations
    real(r8) :: sabs_road                        ! temporary for absorption over road
    real(r8) :: sref_road                        ! temporary for reflected over road
    real(r8), parameter :: errcrit  = .00001_r8  ! error criteria
    real(r8)  :: A_s(bounds%begl:bounds%endl)                           ! Area of the street canyon (normalized)
    real(r8)  :: A_g(bounds%begl:bounds%endl)                           ! Aarea of the ground (normalized)
    real(r8)  :: A_r(bounds%begl:bounds%endl)                           ! Area of the roof
    real(r8)  :: A_w(bounds%begl:bounds%endl)                           ! Area of the wall (normalized)
    logical  :: debug_write = .false.                  ! true => write out many intermediate variables for debugging 
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(coszen)             == (/bounds%endl/)),         sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(canyon_hwr)         == (/bounds%endl/)),         sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(wtroad_perv)        == (/bounds%endl/)),         sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sdir)               == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sdif)               == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(alb_improad_dir)    == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(alb_perroad_dir)    == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(alb_wall_dir)       == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(alb_roof_dir)       == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(alb_improad_dif)    == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(alb_perroad_dif)    == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(alb_wall_dif)       == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(alb_roof_dif)       == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sdir_road)          == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sdir_sunwall)       == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sdir_shadewall)     == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sdif_road)          == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sdif_sunwall)       == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sdif_shadewall)     == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sref_improad_dir)   == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sref_perroad_dir)   == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sref_improad_dif)   == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sref_perroad_dif)   == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sref_sunwall_dir)   == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sref_sunwall_dif)   == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sref_shadewall_dir) == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sref_shadewall_dif) == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sref_roof_dir)      == (/bounds%endl, numrad/)), sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(sref_roof_dif)      == (/bounds%endl, numrad/)), sourcefile, __LINE__)

    associate(                                                           & 
        kgw1d =>    urbanparams_inst%kgw1d_out  , & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor from ground to  one wall[landunit, nzcanm]    
        kgv1d =>    urbanparams_inst%kgv1d_out   ,& ! Input:  [real(r8) (:,:) ]  Monte carlo view factor from ground to vegetation [landunit, nzcanm]  
        kts1d =>    urbanparams_inst%kts1d_out   ,& ! Input:  [real(r8) (:) ]  Monte carlo view factor from ground to sky[landunit, nzcanm]    

        ksg1d =>    urbanparams_inst%ksg1d_out  , & ! Input:  [real(r8) (:) ]    Monte carlo view factor from sky to ground[landunit]               
        ksw1d =>    urbanparams_inst%ksw1d_out , & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor from sky to wall[landunit, nzcanm]   
        ksr1d =>    urbanparams_inst%ksr1d_out , & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor from sky to roof[landunit, nzcanm]                      
        ksv1d =>    urbanparams_inst%ksv1d_out , & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor from sky to vegetation[landunit, nzcanm]      
                        
        kwg1d =>    urbanparams_inst%kwg1d_out  , & ! Input:  [real(r8) (:,:) ]      Monte carlo view factor of from one wall to ground[landunit]  
        kww1d =>    urbanparams_inst%kww1d_out  , & ! Input:  [real(r8) (:,:,:) ]  Monte carlo view factor from on wall to opposing wall [landunit]  
        kwv1d =>    urbanparams_inst%kwv1d_out  , & ! Input:  [real(r8) (:,:,:) ]  Monte carlo view factor from on wall to vegetation[landunit]  
        kws1d =>    urbanparams_inst%kws1d_out  , & ! Input:  [real(r8) (:,:) ]    Monte carlo view factor from on wall to sky[landunit]  

        kvg1d =>    urbanparams_inst%kvg1d_out  , & ! Input:  [real(r8) (:,:) ]      Monte carlo view factor of from vegetation to ground[landunit]  
        kvw1d =>    urbanparams_inst%kvw1d_out  , & ! Input:  [real(r8) (:,:,:) ]  Monte carlo view factor from vegetation to one wall [landunit]  
        kvv1d =>    urbanparams_inst%kvv1d_out  , & ! Input:  [real(r8) (:,:,:) ]  Monte carlo view factor from vegetation to vegetation[landunit]  
        kvs1d =>    urbanparams_inst%kvs1d_out  , & ! Input:  [real(r8) (:,:) ]    Monte carlo view factor from vegetation to sky[landunit]  
        kvr1d =>    urbanparams_inst%kvr1d_out  , & ! Input:  [real(r8) (:,:,:) ]    Monte carlo view factor from ovegetation to roof[landunit]  
        
        krs1d =>    urbanparams_inst%krs1d_out , & ! Input:  [real(r8) (:,:) ]  Monte carlo view factor from roof to sky [landunit, nzcanm]                 
        krv1d =>    urbanparams_inst%krv1d_out , & ! Input:  [real(r8) (:,:,:) ]  Monte carlo view factor  from roof to vegetation [landunit, nzcanm]              

        sabs_roof_dir   =>    solarabs_inst%sabs_roof_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
        sabs_sunwall_dir   =>    solarabs_inst%sabs_sunwall_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
        sabs_shadewall_dir =>    solarabs_inst%sabs_shadewall_dir_lun , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by shadewall per unit wall area per unit incident flux
        sabs_improad_dir   =>    solarabs_inst%sabs_improad_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux
        sabs_perroad_dir   =>    solarabs_inst%sabs_perroad_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
        sabs_br_tree_dir      =>    solarabs_inst%sabs_br_tree_dir_lun      , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by below-roof treeetation per unit treeetation area per unit incident flux
        sabs_ar_tree_dir      =>    solarabs_inst%sabs_ar_tree_dir_lun      , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by above-roof treeetation per unit treeetation area per unit incident flux
        sabs_tree_dir      =>    solarabs_inst%sabs_tree_dir_lun      , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by above-roof treeetation per unit treeetation area per unit incident flux
        sabs_canyon_dir   =>    solarabs_inst%sabs_canyon_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux

        sabs_roof_dif   =>    solarabs_inst%sabs_roof_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
        sabs_sunwall_dif   =>    solarabs_inst%sabs_sunwall_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
        sabs_shadewall_dif =>    solarabs_inst%sabs_shadewall_dif_lun , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
        sabs_improad_dif   =>    solarabs_inst%sabs_improad_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
        sabs_perroad_dif   =>    solarabs_inst%sabs_perroad_dif_lun    , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux
        sabs_br_tree_dif      =>    solarabs_inst%sabs_br_tree_dif_lun      , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by below-roof treeetation per unit treeetation area per unit incident flux
        sabs_ar_tree_dif      =>    solarabs_inst%sabs_ar_tree_dif_lun      , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by above-roof treeetation per unit treeetation area per unit incident flux
        sabs_tree_dif      =>    solarabs_inst%sabs_tree_dif_lun      , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by above-roof treeetation per unit treeetation area per unit incident flux
        sabs_canyon_dif   =>    solarabs_inst%sabs_canyon_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux

        sref_roof_dir   => solarabs_inst%sref_roof_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar reflected by roof per unit ground area per unit incident flux
        sref_sunwall_dir     => solarabs_inst%sref_sunwall_dir_lun     , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by sunwall per unit wall area per unit incident flux
        sref_shadewall_dir   => solarabs_inst%sref_shadewall_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar reflected by shadewall per unit wall area per unit incident flux
        sref_improad_dir     => solarabs_inst%sref_improad_dir_lun     , & ! Output: [real(r8) (:,:) ]  direct  solar reflected by impervious road per unit ground area per unit incident flux
        sref_perroad_dir     => solarabs_inst%sref_perroad_dir_lun     , & ! Output: [real(r8) (:,:) ]  direct  solar reflected by pervious road per unit ground area per unit incident flux
        sref_br_tree_dir     => solarabs_inst%sref_br_tree_dir_lun     , & ! Output: [real(r8) (:,:) ]  direct  solar reflected by below-roof tree per unit ground area per unit incident flux
        sref_ar_tree_dir     => solarabs_inst%sref_ar_tree_dir_lun     , & ! Output: [real(r8) (:,:) ]  direct  solar reflected by above-roof tree per unit ground area per unit incident flux
        sref_tree_dir     => solarabs_inst%sref_tree_dir_lun     , & ! Output: [real(r8) (:,:) ]  direct  solar reflected by above-roof tree per unit ground area per unit incident flux
        sref_canyon_dir     => solarabs_inst%sref_canyon_dir_lun     , & ! Output: [real(r8) (:,:) ]  direct  solar reflected by above-roof tree per unit ground area per unit incident flux

        sref_roof_dif   => solarabs_inst%sref_roof_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by roof per unit ground area per unit incident flux
        sref_sunwall_dif     => solarabs_inst%sref_sunwall_dif_lun     , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by sunwall per unit wall area per unit incident flux
        sref_shadewall_dif   => solarabs_inst%sref_shadewall_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by shadewall per unit wall area per unit incident flux
        sref_improad_dif     => solarabs_inst%sref_improad_dif_lun     , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by impervious road per unit ground area per unit incident flux
        sref_perroad_dif     => solarabs_inst%sref_perroad_dif_lun     , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by pervious road per unit ground area per unit incident flux
        sref_br_tree_dif     => solarabs_inst%sref_br_tree_dif_lun     , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by below-roof tree per unit ground area per unit incident flux
        sref_ar_tree_dif     => solarabs_inst%sref_ar_tree_dif_lun     , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by above-roof tree per unit ground area per unit incident flux
        sref_tree_dif     => solarabs_inst%sref_tree_dif_lun     , & ! Output: [real(r8) (:,:) ]  diffuse solar reflected by above-roof tree per unit ground area per unit incident flux
        sref_canyon_dif     => solarabs_inst%sref_canyon_dif_lun      & ! Output: [real(r8) (:,:) ]  direct  solar reflected by above-roof tree per unit ground area per unit incident flux
         ) 

      debug_write=.false.    
      
      do fl = 1,num_urbanl 
         l = filter_urbanl(fl)
         wtroad_imperv(l) = 1._r8 - wtroad_perv(l) - wtroad_tree(l)
         A_w(l)=ht_roof(l)
         A_r(l)=ht_roof(l)/(canyon_hwr(l)*(1._r8-wtlunit_roof(l))/wtlunit_roof(l))
         A_g(l)=ht_roof(l)/canyon_hwr(l)
         A_s(l)=A_r(l)+A_g(l)        
      end do
            
      do ib = 1,numrad
         do fl = 1,num_urbanl
            l = filter_urbanl(fl)
            if (coszen(l) > 0._r8) then
            
               if (l==3) then   
                   debug_write=.false. 
               else 
                   debug_write=.false. 
               end if
               ! initial absorption and reflection for road and both walls. 
               ! distribute reflected radiation to sky, road, and walls 
               ! according to appropriate view factor. radiation reflected to 
               ! road and walls will undergo multiple reflections within the canyon. 
               ! do separately for direct beam and diffuse radiation.

               ! direct beam

               road_a_dir(l)              = 0.0_r8
               road_r_dir(l)              = 0.0_r8
               improad_a_dir(l)           = (1._r8-alb_improad_dir(l,ib)) * sdir_road(l,ib) 
               improad_r_dir(l)           =     alb_improad_dir(l,ib)  * sdir_road(l,ib) 
               improad_r_sky_dir(l)       = improad_r_dir(l) * kts1d(l)
               improad_r_sunwall_dir(l)   = improad_r_dir(l) * kgw1d(l,1)
               improad_r_shadewall_dir(l) = improad_r_dir(l) * kgw1d(l,1)
               improad_r_br_tree_dir(l) = improad_r_dir(l) * kgv1d(l,1)
               improad_r_ar_tree_dir(l) = improad_r_dir(l) * kgv1d(l,2)

               err_r=improad_r_dir(l) -improad_r_br_tree_dir(l)*A_v1(l)/A_g(l) -improad_r_ar_tree_dir(l)*A_v2(l)/A_g(l) -improad_r_sky_dir(l)-improad_r_sunwall_dir(l)*A_w(l)/A_g(l)-improad_r_shadewall_dir(l)*A_w(l)/A_g(l)
               if (debug_write) then
                  write(6,*) '-----improad--------l=---', l
                  write(6,*) 'err_r', err_r
               end if   
               
               road_a_dir(l)              = road_a_dir(l) + improad_a_dir(l)*wtroad_imperv(l)
               road_r_dir(l)              = road_r_dir(l) + improad_r_dir(l)*wtroad_imperv(l)

               perroad_a_dir(l)           = (1._r8-alb_perroad_dir(l,ib)) * sdir_road(l,ib) 
               perroad_r_dir(l)           =     alb_perroad_dir(l,ib)  * sdir_road(l,ib) 
               perroad_r_sky_dir(l)       = perroad_r_dir(l) * kts1d(l)
               perroad_r_sunwall_dir(l)   = perroad_r_dir(l) * kgw1d(l,1)
               perroad_r_shadewall_dir(l) = perroad_r_dir(l) * kgw1d(l,1)
               perroad_r_br_tree_dir(l) = perroad_r_dir(l) * kgv1d(l,1)
               perroad_r_ar_tree_dir(l) = perroad_r_dir(l) * kgv1d(l,2)
               
               err_r=perroad_r_dir(l) -perroad_r_br_tree_dir(l)*A_v1(l)/A_g(l) -perroad_r_ar_tree_dir(l)*A_v2(l)/A_g(l) -perroad_r_sky_dir(l)-perroad_r_sunwall_dir(l)*A_w(l)/A_g(l)-perroad_r_shadewall_dir(l)*A_w(l)/A_g (l)
               if (debug_write) then
                  write(6,*) '-----perroad--------l=---', l
                  write(6,*) 'err_r', err_r
               end if   
               ! the ground part of road tree is still pervious road
               road_a_dir(l)              = road_a_dir(l) + perroad_a_dir(l)*(wtroad_perv(l)+wtroad_tree(l))
               road_r_dir(l)              = road_r_dir(l) + perroad_r_dir(l)*(wtroad_perv(l)+wtroad_tree(l))
              
               road_r_sky_dir(l)          = road_r_dir(l) * kts1d(l)
               road_r_sunwall_dir(l)      = road_r_dir(l) * kgw1d(l,1)
               road_r_shadewall_dir(l)    = road_r_dir(l) * kgw1d(l,1)
               road_r_br_tree_dir(l)    = road_r_dir(l) * kgv1d(l,1)
               road_r_ar_tree_dir(l)    = road_r_dir(l) * kgv1d(l,2)

               err_r=road_r_dir(l) -road_r_br_tree_dir(l)*A_v1(l)/A_g(l) -road_r_ar_tree_dir(l)*A_v2(l)/A_g(l) -road_r_sky_dir(l)-road_r_sunwall_dir(l)*A_w(l)/A_g(l)-road_r_shadewall_dir(l)*A_w(l)/A_g(l)
               if (debug_write) then
                  write(6,*) '-----road--------l=---', l
                  write(6,*) 'err_r', err_r
               end if 
               
               sunwall_a_dir(l)           = (1._r8-alb_wall_dir(l,ib)) * sdir_sunwall(l,ib)
               sunwall_r_dir(l)           =     alb_wall_dir(l,ib)  * sdir_sunwall(l,ib)
               sunwall_r_sky_dir(l)       = sunwall_r_dir(l) * kws1d(l,1)
               sunwall_r_road_dir(l)      = sunwall_r_dir(l) * kwg1d(l,1)
               sunwall_r_shadewall_dir(l) = sunwall_r_dir(l) * kww1d(l,1,1)
               sunwall_r_br_tree_dir(l) = sunwall_r_dir(l) * kwv1d(l,1,1)
               sunwall_r_ar_tree_dir(l) = sunwall_r_dir(l) * kwv1d(l,1,2)
               
               err_r=sunwall_r_dir(l) -sunwall_r_br_tree_dir(l)*A_v1(l)/A_w(l) -sunwall_r_ar_tree_dir(l)*A_v2(l)/A_w(l) -sunwall_r_sky_dir(l)-sunwall_r_road_dir(l)*A_g(l)/A_w(l)-sunwall_r_shadewall_dir(l)*A_w(l)/A_w(l)
               if (debug_write) then
                  write(6,*) '------sunwall--------l=---', l
                  write(6,*) 'err_r', err_r
               end if    
               
               shadewall_a_dir(l)         = (1._r8-alb_wall_dir(l,ib)) * sdir_shadewall(l,ib)
               shadewall_r_dir(l)         =     alb_wall_dir(l,ib)  * sdir_shadewall(l,ib)
               shadewall_r_sky_dir(l)     = shadewall_r_dir(l) * kws1d(l,1)
               shadewall_r_road_dir(l)    = shadewall_r_dir(l) * kwg1d(l,1)
               shadewall_r_sunwall_dir(l) = shadewall_r_dir(l) * kww1d(l,1,1)
               shadewall_r_br_tree_dir(l) = shadewall_r_dir(l) * kwv1d(l,1,1)
               shadewall_r_ar_tree_dir(l) = shadewall_r_dir(l) * kwv1d(l,1,2)

               err_r=shadewall_r_dir(l) -shadewall_r_br_tree_dir(l)*A_v1(l)/A_w(l) -shadewall_r_ar_tree_dir(l)*A_v2(l)/A_w(l) -shadewall_r_sky_dir(l)-shadewall_r_road_dir(l)*A_g(l)/A_w(l)-shadewall_r_sunwall_dir(l)*A_w(l)/A_w(l)
               if (debug_write) then
                  write(6,*) '-----shadewall--------l=---', l
                  write(6,*) 'err_r', err_r
               end if   
               
               br_tree_a_dir(l)         = (1._r8-alb_br_tree_dir(l,ib)) * sdir_br_tree(l,ib)
               br_tree_r_dir(l)         =     alb_br_tree_dir(l,ib)  * sdir_br_tree(l,ib)
               br_tree_r_sky_dir(l)     = br_tree_r_dir(l) * kvs1d(l,1)
               br_tree_r_road_dir(l)    = br_tree_r_dir(l) * kvg1d(l,1)
               br_tree_r_sunwall_dir(l) = br_tree_r_dir(l) * kvw1d(l,1,1)
               br_tree_r_shadewall_dir(l) = br_tree_r_dir(l) * kvw1d(l,1,1)   
               br_tree_r_br_tree_dir(l) = br_tree_r_dir(l) * kvv1d(l,1,1)
               br_tree_r_ar_tree_dir(l) = br_tree_r_dir(l) * kvv1d(l,1,2)

               err_r=br_tree_r_dir(l) -br_tree_r_br_tree_dir(l) -br_tree_r_ar_tree_dir(l)*A_v2(l)/A_v1(l) -br_tree_r_sky_dir(l)-br_tree_r_road_dir(l)*A_g(l)/A_v1(l)-br_tree_r_sunwall_dir(l)*A_w(l)/A_v1(l)-br_tree_r_shadewall_dir(l)*A_w(l)/A_v1(l)
               if (debug_write) then
                  write(6,*) '-----br_tree-------l=---', l
                  write(6,*) 'err_r',err_r
               end if  
               
               ar_tree_a_dir(l)         = (1._r8-alb_ar_tree_dir(l,ib)) * sdir_ar_tree(l,ib)
               ar_tree_r_dir(l)         =     alb_ar_tree_dir(l,ib)  * sdir_ar_tree(l,ib)
               ar_tree_r_sky_dir(l)     = ar_tree_r_dir(l) * kvs1d(l,2)
               ar_tree_r_road_dir(l)    = ar_tree_r_dir(l) * kvg1d(l,2)
               ar_tree_r_sunwall_dir(l) = ar_tree_r_dir(l) * kvw1d(l,2,1)
               ar_tree_r_shadewall_dir(l) = ar_tree_r_dir(l) * kvw1d(l,2,1)   
               ar_tree_r_roof_dir(l) = ar_tree_r_dir(l) * kvr1d(l,2,2)
               ar_tree_r_br_tree_dir(l) = ar_tree_r_dir(l) * kvv1d(l,2,1)
               ar_tree_r_ar_tree_dir(l) = ar_tree_r_dir(l) * kvv1d(l,2,2)

               !err_r=ar_tree_r_dir(l) -ar_tree_r_ar_tree_dir(l) -ar_tree_r_br_tree_dir(l)*A_v1(l)/A_v2(l) -ar_tree_r_sky_dir(l)-ar_tree_r_road_dir(l)*A_g(l)/A_v2(l)-ar_tree_r_sunwall_dir(l)*A_w(l)/A_v2(l)-ar_tree_r_shadewall_dir(l)*A_w(l)/A_v2(l)-ar_tree_r_roof_dir(l)*A_r(l)/A_v2(l)
               !if (debug_write) then
                !  write(6,*) '-----ar_tree------l=---', l
                !  write(6,*) 'err_r', err_r  
               !end if   
               roof_a_dir(l)           = (1._r8-alb_roof_dir(l,ib)) * sdir_roof(l,ib)
               roof_r_dir(l)           =     alb_roof_dir(l,ib)  * sdir_roof(l,ib)
               roof_r_sky_dir(l)       = roof_r_dir(l) * krs1d(l,2)
               roof_r_ar_tree_dir(l) = roof_r_dir(l) * krv1d(l,2,2)

               err_r=roof_r_dir(l) -roof_r_ar_tree_dir(l)*A_v2(l)/A_r(l) -roof_r_sky_dir(l)
               if (debug_write) then
                  write(6,*) '-----roof------l=---', l
                  write(6,*) 'err_r', err_r
               end if   

               if (debug_write) then
                  write(6,*) '------------- Direct View Factors and Albedo Debug Output (l =', l, ') -------------'
                  
                  ! Albedo terms
                  write(6,*) 'alb_improad_dir =', alb_improad_dir(l,ib)
                  write(6,*) 'alb_perroad_dir =', alb_perroad_dir(l,ib)
                  write(6,*) 'alb_wall_dir    =', alb_wall_dir(l,ib)
                  write(6,*) 'alb_ar_tree_dir    =', alb_ar_tree_dir(l,ib)
                  write(6,*) 'alb_br_tree_dir    =', alb_br_tree_dir(l,ib)
                  write(6,*) 'alb_roof_dir    =', alb_roof_dir(l,ib)

                  ! Direct solar input
                  write(6,*) 'sdir_road     =', sdir_road(l,ib)
                  write(6,*) 'sdir_sunwall  =', sdir_sunwall(l,ib)
                  write(6,*) 'sdir_shadewall=', sdir_shadewall(l,ib)
                  write(6,*) 'sdir_br_tree  =', sdir_br_tree(l,ib)
                  write(6,*) 'sdir_ar_tree  =', sdir_ar_tree(l,ib)
                  write(6,*) 'sdir_roof=', sdir_roof(l,ib) 

                  ! Direct view factor scalars
                  write(6,*) 'kts1d           =', kts1d(l)
                  write(6,*) 'kgw1d           =', kgw1d(l,1)
                  write(6,*) 'kgv1d_1           =', kgv1d(l,1)
                  write(6,*) 'kgv1d_2           =', kgv1d(l,2)
                  
                  write(6,*) 'kws1d           =', kws1d(l,1)
                  write(6,*) 'kwg1d           =', kwg1d(l,1)
                  write(6,*) 'kww1d           =', kww1d(l,1,1)
                  write(6,*) 'kwv1d_1           =', kwv1d(l,1,1)
                  write(6,*) 'kwv1d_2           =', kwv1d(l,1,2)
                  
                  write(6,*) 'kvs1d_1      =', kvs1d(l,1)
                  write(6,*) 'kvg1d_1      =', kvg1d(l,1)
                  write(6,*) 'kvw1d_1           =', kvw1d(l,1,1)
                  write(6,*) 'kvv1d_1_1           =', kvv1d(l,1,1)
                  write(6,*) 'kvv1d_1_2           =', kvv1d(l,1,2)
                  write(6,*) 'kvg1d_1      =', kvg1d(l,1)
                  
                  write(6,*) 'kvs1d_2      =', kvs1d(l,2)
                  write(6,*) 'kvg1d_2      =', kvg1d(l,2)
                  write(6,*) 'kvw1d_2           =', kvw1d(l,2,1)
                  write(6,*) 'kvv1d_2_1           =', kvv1d(l,2,1)
                  write(6,*) 'kvv1d_2_2           =', kvv1d(l,2,2)
                  write(6,*) 'kvg1d_2      =',  kvg1d(l,2)
                  write(6,*) 'kvr1d_2      =',  kvr1d(l,2,2)
                  
                  write(6,*) 'krs1d           =', krs1d(l,2)
                  write(6,*) 'krv1d_2         =', krv1d(l,2,2)                  

                  ! Summary
                  write(6,*) '-------------------------------------------------------------'
               end if
               
               ! diffuse

               road_a_dif(l)              = 0.0_r8
               road_r_dif(l)              = 0.0_r8
               improad_a_dif(l)           = (1._r8-alb_improad_dif(l,ib)) * sdif_road(l,ib) 
               improad_r_dif(l)           =     alb_improad_dif(l,ib)  * sdif_road(l,ib) 
               improad_r_sky_dif(l)       = improad_r_dif(l) * kts1d(l)
               improad_r_sunwall_dif(l)   = improad_r_dif(l) * kgw1d(l,1)
               improad_r_shadewall_dif(l) = improad_r_dif(l) * kgw1d(l,1)
               improad_r_br_tree_dif(l) = improad_r_dif(l) * kgv1d(l,1)
               improad_r_ar_tree_dif(l) = improad_r_dif(l) * kgv1d(l,2)
               road_a_dif(l)              = road_a_dif(l) + improad_a_dif(l)*wtroad_imperv(l)
               road_r_dif(l)              = road_r_dif(l) + improad_r_dif(l)*wtroad_imperv(l)
               
               err_r=improad_r_dif(l) -improad_r_br_tree_dif(l)*A_v1(l)/A_g(l) -improad_r_ar_tree_dif(l)*A_v2(l)/A_g(l) -improad_r_sky_dif(l)-improad_r_sunwall_dif(l)*A_w(l)/A_g(l)-improad_r_shadewall_dif(l)*A_w(l)/A_g(l)
               if (debug_write) then
                  write(6,*) '-----improad--------l=---', l
                  write(6,*) 'err_r', err_r
               end if   

               perroad_a_dif(l)           = (1._r8-alb_perroad_dif(l,ib)) * sdif_road(l,ib) 
               perroad_r_dif(l)           =     alb_perroad_dif(l,ib)  * sdif_road(l,ib) 
               perroad_r_sky_dif(l)       = perroad_r_dif(l) * kts1d(l)
               perroad_r_sunwall_dif(l)   = perroad_r_dif(l) * kgw1d(l,1)
               perroad_r_shadewall_dif(l) = perroad_r_dif(l) * kgw1d(l,1)
               perroad_r_br_tree_dif(l) = perroad_r_dif(l) * kgv1d(l,1)
               perroad_r_ar_tree_dif(l) = perroad_r_dif(l) * kgv1d(l,2)
               
               err_r=perroad_r_dif(l) -perroad_r_br_tree_dif(l)*A_v1(l)/A_g(l) -perroad_r_ar_tree_dif(l)*A_v2(l)/A_g(l) -perroad_r_sky_dif(l)-perroad_r_sunwall_dif(l)*A_w(l)/A_g(l)-perroad_r_shadewall_dif(l)*A_w(l)/A_g (l)
               if (debug_write) then
                  write(6,*) '-----perroad--------l=---', l
                  write(6,*) 'err_r', err_r
               end if   
               
               road_a_dif(l)              = road_a_dif(l) + perroad_a_dif(l)*(wtroad_perv(l)+wtroad_tree(l))
               road_r_dif(l)              = road_r_dif(l) + perroad_r_dif(l)*(wtroad_perv(l)+wtroad_tree(l))

               road_r_sky_dif(l)          = road_r_dif(l) * kts1d(l)
               road_r_sunwall_dif(l)      = road_r_dif(l) * kgw1d(l,1)
               road_r_shadewall_dif(l)    = road_r_dif(l) * kgw1d(l,1)
               road_r_br_tree_dif(l)    = road_r_dif(l) * kgv1d(l,1)
               road_r_ar_tree_dif(l)    = road_r_dif(l) * kgv1d(l,2)

               err_r=road_r_dif(l) -road_r_br_tree_dif(l)*A_v1(l)/A_g(l) -road_r_ar_tree_dif(l)*A_v2(l)/A_g(l) -road_r_sky_dif(l)-road_r_sunwall_dif(l)*A_w(l)/A_g(l)-road_r_shadewall_dif(l)*A_w(l)/A_g(l)
               if (debug_write) then
                  write(6,*) '-----road--------l=---', l
                  write(6,*) 'err_r', err_r
               end if  

               sunwall_a_dif(l)           = (1._r8-alb_wall_dif(l,ib)) * sdif_sunwall(l,ib)
               sunwall_r_dif(l)           =     alb_wall_dif(l,ib)  * sdif_sunwall(l,ib)
               sunwall_r_sky_dif(l)       = sunwall_r_dif(l) * kws1d(l,1)
               sunwall_r_road_dif(l)      = sunwall_r_dif(l) * kwg1d(l,1)
               sunwall_r_shadewall_dif(l) = sunwall_r_dif(l) * kww1d(l,1,1)
               sunwall_r_br_tree_dif(l) = sunwall_r_dif(l) * kwv1d(l,1,1)
               sunwall_r_ar_tree_dif(l) = sunwall_r_dif(l) * kwv1d(l,1,2)

               err_r=sunwall_r_dif(l) -sunwall_r_br_tree_dif(l)*A_v1(l)/A_w(l) -sunwall_r_ar_tree_dif(l)*A_v2(l)/A_w(l) -sunwall_r_sky_dif(l)-sunwall_r_road_dif(l)*A_g(l)/A_w(l)-sunwall_r_shadewall_dif(l)*A_w(l)/A_w(l)
               if (debug_write) then
                  write(6,*) '------sunwall--------l=---', l
                  write(6,*) 'err_r', err_r
               end if   

               shadewall_a_dif(l)         = (1._r8-alb_wall_dif(l,ib)) * sdif_shadewall(l,ib)
               shadewall_r_dif(l)         =     alb_wall_dif(l,ib)  * sdif_shadewall(l,ib)
               shadewall_r_sky_dif(l)     = shadewall_r_dif(l) * kws1d(l,1)
               shadewall_r_road_dif(l)    = shadewall_r_dif(l) * kwg1d(l,1)
               shadewall_r_sunwall_dif(l) = shadewall_r_dif(l) * kww1d(l,1,1)
               shadewall_r_br_tree_dif(l) = shadewall_r_dif(l) * kwv1d(l,1,1)
               shadewall_r_ar_tree_dif(l) = shadewall_r_dif(l) * kwv1d(l,1,2)

               err_r=shadewall_r_dif(l) -shadewall_r_br_tree_dif(l)*A_v1(l)/A_w(l) -shadewall_r_ar_tree_dif(l)*A_v2(l)/A_w(l) -shadewall_r_sky_dif(l)-shadewall_r_road_dif(l)*A_g(l)/A_w(l)-shadewall_r_sunwall_dif(l)*A_w(l)/A_w(l)
               if (debug_write) then
                  write(6,*) '-----shadewall--------l=---', l
                  write(6,*) 'err_r', err_r
               end if   

               br_tree_a_dif(l)         = (1._r8-alb_br_tree_dif(l,ib)) * sdif_br_tree(l,ib)
               br_tree_r_dif(l)         =     alb_br_tree_dif(l,ib)  * sdif_br_tree(l,ib)
               br_tree_r_sky_dif(l)     = br_tree_r_dif(l) * kvs1d(l,1)
               br_tree_r_road_dif(l)    = br_tree_r_dif(l) * kvg1d(l,1)
               br_tree_r_sunwall_dif(l) = br_tree_r_dif(l) * kvw1d(l,1,1)
               br_tree_r_shadewall_dif(l) = br_tree_r_dif(l) * kvw1d(l,1,1)   
               br_tree_r_br_tree_dif(l) = br_tree_r_dif(l) * kvv1d(l,1,1)
               br_tree_r_ar_tree_dif(l) = br_tree_r_dif(l) * kvv1d(l,1,2)
       
               err_r=br_tree_r_dif(l) -br_tree_r_br_tree_dif(l) -br_tree_r_ar_tree_dif(l)*A_v2(l)/A_v1(l) -br_tree_r_sky_dif(l)-br_tree_r_road_dif(l)*A_g(l)/A_v1(l)-br_tree_r_sunwall_dif(l)*A_w(l)/A_v1(l)-br_tree_r_shadewall_dif(l)*A_w(l)/A_v1(l)
               if (debug_write) then
                  write(6,*) '-----br_tree-------l=---', l
                  write(6,*) 'err_r',err_r
               end if   

               ar_tree_a_dif(l)         = (1._r8-alb_ar_tree_dif(l,ib)) * sdif_ar_tree(l,ib)
               ar_tree_r_dif(l)         =     alb_ar_tree_dif(l,ib)  * sdif_ar_tree(l,ib)
               ar_tree_r_sky_dif(l)     = ar_tree_r_dif(l) * kvs1d(l,2)
               ar_tree_r_road_dif(l)    = ar_tree_r_dif(l) * kvg1d(l,2)
               ar_tree_r_sunwall_dif(l) = ar_tree_r_dif(l) * kvw1d(l,2,1)
               ar_tree_r_shadewall_dif(l) = ar_tree_r_dif(l) * kvw1d(l,2,1)   
               ar_tree_r_roof_dif(l) = ar_tree_r_dif(l) * kvr1d(l,2,2)

               ar_tree_r_br_tree_dif(l) = ar_tree_r_dif(l) * kvv1d(l,2,1)
               ar_tree_r_ar_tree_dif(l) = ar_tree_r_dif(l) * kvv1d(l,2,2)

              ! err_r=ar_tree_r_dif(l) -ar_tree_r_ar_tree_dif(l) -ar_tree_r_br_tree_dif(l)*A_v1(l)/A_v2(l) -ar_tree_r_sky_dif(l)-ar_tree_r_road_dif(l)*A_g(l)/A_v2(l)-ar_tree_r_sunwall_dif(l)*A_w(l)/A_v2(l)-ar_tree_r_shadewall_dif(l)*A_w(l)/A_v2(l)-ar_tree_r_roof_dif(l)*A_r(l)/A_v2(l)
              ! if (debug_write) then
              !    write(6,*) '-----ar_tree------l=---', l
              !    write(6,*) 'err_r', err_r  
               !end if   

               roof_a_dif(l)           = (1._r8-alb_roof_dif(l,ib)) * sdif_roof(l,ib)
               roof_r_dif(l)           =     alb_roof_dif(l,ib)  * sdif_roof(l,ib)
               roof_r_sky_dif(l)       = roof_r_dif(l) * krs1d(l,2)
               roof_r_ar_tree_dif(l) = roof_r_dif(l) * krv1d(l,2,2)
               
               err_r=roof_r_dif(l) -roof_r_ar_tree_dif(l)*A_v2(l)/A_r(l) -roof_r_sky_dif(l)
               if (debug_write) then
                  write(6,*) '-----roof------l=---', l
                  write(6,*) 'err_r', err_r
               end if   

               if (debug_write) then
                  write(6,*) '------------- Diffuse Solar Radiation Debug Output (l =', l, ') -------------'

                  ! Diffuse albedo terms
                  write(6,*) 'alb_improad_dif =', alb_improad_dif(l,ib)
                  write(6,*) 'alb_perroad_dif =', alb_perroad_dif(l,ib)
                  write(6,*) 'alb_wall_dif    =', alb_wall_dif(l,ib)
                  write(6,*) 'alb_ar_tree_dif    =', alb_ar_tree_dif(l,ib)
                  write(6,*) 'alb_br_tree_dif    =', alb_br_tree_dif(l,ib)
                  write(6,*) 'alb_roof_dif    =', alb_roof_dif(l,ib)

                  ! Diffuse solar input
                  write(6,*) 'sdif_road     =', sdif_road(l,ib)
                  write(6,*) 'sdif_sunwall  =', sdif_sunwall(l,ib)
                  write(6,*) 'sdif_shadewall=', sdif_shadewall(l,ib)
                  write(6,*) 'sdif_br_tree  =', sdif_br_tree(l,ib)
                  write(6,*) 'sdif_ar_tree  =', sdif_ar_tree(l,ib)
                  write(6,*) 'sdif_roof=', sdif_roof(l,ib) 

                  write(6,*) '-------------------------------------------------------------'
               end if

               sabs_improad_dir(l,ib)   = improad_a_dir(l)
               sabs_perroad_dir(l,ib)   = perroad_a_dir(l)
               sabs_sunwall_dir(l,ib)   = sunwall_a_dir(l)
               sabs_shadewall_dir(l,ib) = shadewall_a_dir(l)
               sabs_br_tree_dir(l,ib)   = br_tree_a_dir(l)
               sabs_ar_tree_dir(l,ib)   = ar_tree_a_dir(l)
               sabs_roof_dir(l,ib)   = roof_a_dir(l)
               
               sabs_improad_dif(l,ib)   = improad_a_dif(l)
               sabs_perroad_dif(l,ib)   = perroad_a_dif(l)
               sabs_sunwall_dif(l,ib)   = sunwall_a_dif(l)
               sabs_shadewall_dif(l,ib) = shadewall_a_dif(l)             
               sabs_br_tree_dif(l,ib)   = br_tree_a_dif(l)
               sabs_ar_tree_dif(l,ib)   = ar_tree_a_dif(l)
               sabs_roof_dif(l,ib)   = roof_a_dif(l)
               
               sref_improad_dir(l,ib)   = improad_r_sky_dir(l) 
               sref_perroad_dir(l,ib)   = perroad_r_sky_dir(l) 
               sref_sunwall_dir(l,ib)   = sunwall_r_sky_dir(l) 
               sref_shadewall_dir(l,ib) = shadewall_r_sky_dir(l) 
               sref_br_tree_dir(l,ib)   = br_tree_r_sky_dir(l) 
               sref_ar_tree_dir(l,ib)   = ar_tree_r_sky_dir(l) 
               sref_roof_dir(l,ib)   = roof_r_sky_dir(l)
               
               sref_improad_dif(l,ib)   = improad_r_sky_dif(l)
               sref_perroad_dif(l,ib)   = perroad_r_sky_dif(l)
               sref_sunwall_dif(l,ib)   = sunwall_r_sky_dif(l)
               sref_shadewall_dif(l,ib) = shadewall_r_sky_dif(l)
               sref_br_tree_dif(l,ib)   = br_tree_r_sky_dif(l) 
               sref_ar_tree_dif(l,ib)   = ar_tree_r_sky_dif(l) 
               sref_roof_dif(l,ib)   = roof_r_sky_dif(l)

               if (debug_write) then
                  write(6,*) '------------ Absorbed and Reflected SW Flux (sky view component only) ------------'
                  write(6,*) 'l =', l, ', ib =', ib

                  ! Absorbed direct shortwave
                  write(6,*) 'sabs_improad_dir =', sabs_improad_dir(l,ib)
                  write(6,*) 'sabs_perroad_dir =', sabs_perroad_dir(l,ib)
                  write(6,*) 'sabs_sunwall_dir =', sabs_sunwall_dir(l,ib)
                  write(6,*) 'sabs_shadewall_dir =', sabs_shadewall_dir(l,ib)
                  write(6,*) 'sabs_br_tree_dir =', sabs_br_tree_dir(l,ib)
                  write(6,*) 'sabs_ar_tree_dir =', sabs_ar_tree_dir(l,ib)
                  write(6,*) 'sabs_roof_dir =', sabs_roof_dir(l,ib)

                  ! Absorbed diffuse shortwave
                  write(6,*) 'sabs_improad_dif =', sabs_improad_dif(l,ib)
                  write(6,*) 'sabs_perroad_dif =', sabs_perroad_dif(l,ib)
                  write(6,*) 'sabs_sunwall_dif =', sabs_sunwall_dif(l,ib)
                  write(6,*) 'sabs_shadewall_dif =', sabs_shadewall_dif(l,ib)
                  write(6,*) 'sabs_br_tree_dif =', sabs_br_tree_dif(l,ib)
                  write(6,*) 'sabs_ar_tree_dif =', sabs_ar_tree_dif(l,ib)
                  write(6,*) 'sabs_roof_dif =', sabs_roof_dif(l,ib)

                  ! Reflected (to sky) direct shortwave
                  write(6,*) 'sref_improad_dir =', sref_improad_dir(l,ib)
                  write(6,*) 'sref_perroad_dir =', sref_perroad_dir(l,ib)
                  write(6,*) 'sref_sunwall_dir =', sref_sunwall_dir(l,ib)
                  write(6,*) 'sref_shadewall_dir =', sref_shadewall_dir(l,ib)
                  write(6,*) 'sref_br_tree_dir =', sref_br_tree_dir(l,ib)
                  write(6,*) 'sref_ar_tree_dir =', sref_ar_tree_dir(l,ib)
                  write(6,*) 'sref_roof_dir =', sref_roof_dir(l,ib)

                  ! Reflected (to sky) diffuse shortwave
                  write(6,*) 'sref_improad_dif =', sref_improad_dif(l,ib)
                  write(6,*) 'sref_perroad_dif =', sref_perroad_dif(l,ib)
                  write(6,*) 'sref_sunwall_dif =', sref_sunwall_dif(l,ib)
                  write(6,*) 'sref_shadewall_dif =', sref_shadewall_dif(l,ib)
                  write(6,*) 'sref_br_tree_dif =', sref_br_tree_dif(l,ib)
                  write(6,*) 'sref_ar_tree_dif =', sref_ar_tree_dif(l,ib)
                  write(6,*) 'sref_roof_dif =', sref_roof_dif(l,ib) 

                  write(6,*) '----------------------------------------------------------------------------------'
               end if
               
            endif

         end do

         ! absorption and reflection for walls and road with multiple reflections
         ! (i.e., absorb and reflect initial reflection in canyon and allow for 
         ! subsequent scattering)
         !
         ! (1) absorption and reflection of scattered solar radiation
         !     road: reflected fluxes from walls need to be projected to ground area
         !     wall: reflected flux from road needs to be projected to wall area
         !
         ! (2) add absorbed radiation for ith reflection to total absorbed
         !
         ! (3) distribute reflected radiation to sky, road, and walls according to view factors
         !
         ! (4) add solar reflection to sky for ith reflection to total reflection
         !
         ! (5) stop iteration when absorption for ith reflection is less than some nominal amount. 
         !     small convergence criteria is required to ensure solar radiation is conserved
         !
         ! do separately for direct beam and diffuse

         do fl = 1,num_urbanl
            l = filter_urbanl(fl)
            if (l==3) then   
                debug_write=.false. 
            else 
                debug_write=.false. 
            end if 
            if (coszen(l) > 0._r8) then

               ! reflected difect beam
               do iter_dir = 1, n
                  ! step (1)
                  ! The original vf_rw, vf_wr,vf_sr seems to be unweighted view factor
                  stot(l) = (sunwall_r_road_dir(l) + shadewall_r_road_dir(l) + br_tree_r_road_dir(l)&
                            + ar_tree_r_road_dir(l))

                  road_a_dir(l) = 0.0_r8
                  road_r_dir(l) = 0.0_r8
                  improad_a_dir(l) = (1._r8-alb_improad_dir(l,ib)) * stot(l) 
                  improad_r_dir(l) =     alb_improad_dir(l,ib)  * stot(l) 
                  road_a_dir(l)    = road_a_dir(l) + improad_a_dir(l)*wtroad_imperv(l)
                  road_r_dir(l)    = road_r_dir(l) + improad_r_dir(l)*wtroad_imperv(l)
                  perroad_a_dir(l) = (1._r8-alb_perroad_dir(l,ib)) * stot(l) 
                  perroad_r_dir(l) =     alb_perroad_dir(l,ib)  * stot(l) 
                  road_a_dir(l)    = road_a_dir(l) + perroad_a_dir(l)*(wtroad_perv(l)+wtroad_tree(l))
                  road_r_dir(l)    = road_r_dir(l) + perroad_r_dir(l)*(wtroad_perv(l)+wtroad_tree(l))

                  stot(l) = road_r_sunwall_dir(l) + shadewall_r_sunwall_dir(l) + br_tree_r_sunwall_dir(l) &
                            + ar_tree_r_sunwall_dir(l)
                  sunwall_a_dir(l) = (1._r8-alb_wall_dir(l,ib)) * stot(l)
                  sunwall_r_dir(l) =     alb_wall_dir(l,ib)  * stot(l)

                  stot(l) = road_r_shadewall_dir(l) + sunwall_r_shadewall_dir(l) + br_tree_r_shadewall_dir(l) &
                            + ar_tree_r_shadewall_dir(l)
                  shadewall_a_dir(l) = (1._r8-alb_wall_dir(l,ib)) * stot(l)
                  shadewall_r_dir(l) =     alb_wall_dir(l,ib)  * stot(l)

                  stot(l) = road_r_br_tree_dir(l) + sunwall_r_br_tree_dir(l) + shadewall_r_br_tree_dir(l) &
                            + br_tree_r_br_tree_dir(l) +ar_tree_r_br_tree_dir(l)
                  br_tree_a_dir(l) = (1._r8-alb_br_tree_dir(l,ib)) * stot(l)
                  br_tree_r_dir(l) =     alb_br_tree_dir(l,ib)  * stot(l)
    
                  stot(l) = road_r_ar_tree_dir(l) + sunwall_r_ar_tree_dir(l) + shadewall_r_ar_tree_dir(l) &
                        + br_tree_r_ar_tree_dir(l) +ar_tree_r_ar_tree_dir(l)+   roof_r_ar_tree_dir(l)
                  ar_tree_a_dir(l) = (1._r8-alb_ar_tree_dir(l,ib)) * stot(l)
                  ar_tree_r_dir(l) =     alb_ar_tree_dir(l,ib)  * stot(l)

                  stot(l) =  ar_tree_r_roof_dir(l) 
                  roof_a_dir(l) = (1._r8-alb_roof_dir(l,ib)) * stot(l)
                  roof_r_dir(l) =     alb_roof_dir(l,ib)  * stot(l)
                      
                  ! step (2)

                  sabs_improad_dir(l,ib)   = sabs_improad_dir(l,ib)   + improad_a_dir(l)
                  sabs_perroad_dir(l,ib)   = sabs_perroad_dir(l,ib)   + perroad_a_dir(l)
                  sabs_sunwall_dir(l,ib)   = sabs_sunwall_dir(l,ib)   + sunwall_a_dir(l)
                  sabs_shadewall_dir(l,ib)   = sabs_shadewall_dir(l,ib)   + shadewall_a_dir(l)
                  sabs_br_tree_dir(l,ib) = sabs_br_tree_dir(l,ib) + br_tree_a_dir(l)
                  sabs_ar_tree_dir(l,ib) = sabs_ar_tree_dir(l,ib) + ar_tree_a_dir(l)
                  sabs_roof_dir(l,ib) = sabs_roof_dir(l,ib) +roof_a_dir(l)
                  
                  if (debug_write) then
                    write(6,*) '------------ Absorbed and Reflected SW Flux (sky view component only) ------------'
                    write(6,*) 'l =', l, ', ib = n=', ib,n

                    ! Absorbed direct shortwave
                    write(6,*) 'sabs_improad_dir =', sabs_improad_dir(l,ib)
                    write(6,*) 'sabs_perroad_dir =', sabs_perroad_dir(l,ib)
                    write(6,*) 'sabs_sunwall_dir =', sabs_sunwall_dir(l,ib)
                    write(6,*) 'sabs_shadewall_dir =', sabs_shadewall_dir(l,ib)
                    write(6,*) 'sabs_br_tree_dir =', sabs_br_tree_dir(l,ib)
                    write(6,*) 'sabs_ar_tree_dir =', sabs_ar_tree_dir(l,ib)
                    write(6,*) 'sabs_roof_dir =', sabs_roof_dir(l,ib)

                    write(6,*) '----------------------------------------------------------------------------------'
                  end if  
                
                  ! step (3)

                  improad_r_sky_dir(l)       = improad_r_dir(l) * kts1d(l)
                  improad_r_sunwall_dir(l)   = improad_r_dir(l) * kgw1d(l,1)
                  improad_r_shadewall_dir(l) = improad_r_dir(l) * kgw1d(l,1)
                  improad_r_br_tree_dir(l) = improad_r_dir(l) * kgv1d(l,1)
                  improad_r_ar_tree_dir(l) = improad_r_dir(l) * kgv1d(l,2)
                  
                  err_r=improad_r_dir(l) -improad_r_br_tree_dir(l)*A_v1(l)/A_g(l) -improad_r_ar_tree_dir(l)*A_v2(l)/A_g(l) -improad_r_sky_dir(l)-improad_r_sunwall_dir(l)*A_w(l)/A_g(l)-improad_r_shadewall_dir(l)*A_w(l)/A_g(l)
                  if (debug_write) then
                     write(6,*) '-----improad--------l=---n=', l,n
                     write(6,*) 'err_r', err_r
                  end if  
                  
                  perroad_r_sky_dir(l)       = perroad_r_dir(l) * kts1d(l)
                  perroad_r_sunwall_dir(l)   = perroad_r_dir(l) * kgw1d(l,1)
                  perroad_r_shadewall_dir(l) = perroad_r_dir(l) * kgw1d(l,1)
                  perroad_r_br_tree_dir(l) = perroad_r_dir(l) * kgv1d(l,1)
                  perroad_r_ar_tree_dir(l) = perroad_r_dir(l) * kgv1d(l,2)
                  err_r=perroad_r_dir(l) -perroad_r_br_tree_dir(l)*A_v1(l)/A_g(l) -perroad_r_ar_tree_dir(l)*A_v2(l)/A_g(l) -perroad_r_sky_dir(l)-perroad_r_sunwall_dir(l)*A_w(l)/A_g(l)-perroad_r_shadewall_dir(l)*A_w(l)/A_g (l)
                  if (debug_write) then
                     write(6,*) '-----perroad--------l=---n=', l,n
                     write(6,*) 'err_r', err_r
                  end if   
                  road_r_sky_dir(l)          = road_r_dir(l) * kts1d(l)
                  road_r_sunwall_dir(l)      = road_r_dir(l) * kgw1d(l,1)
                  road_r_shadewall_dir(l)    = road_r_dir(l) * kgw1d(l,1)
                  road_r_br_tree_dir(l)    = road_r_dir(l) * kgv1d(l,1)
                  road_r_ar_tree_dir(l)    = road_r_dir(l) * kgv1d(l,2)
                  err_r=road_r_dir(l) -road_r_br_tree_dir(l)*A_v1(l)/A_g(l) -road_r_ar_tree_dir(l)*A_v2(l)/A_g(l) -road_r_sky_dir(l)-road_r_sunwall_dir(l)*A_w(l)/A_g(l)-road_r_shadewall_dir(l)*A_w(l)/A_g(l)
                  if (debug_write) then
                     write(6,*) '-----road--------l=---n=', l,n
                     write(6,*) 'err_r', err_r
                  end if 
                  sunwall_r_sky_dir(l)       = sunwall_r_dir(l) * kws1d(l,1)
                  sunwall_r_road_dir(l)      = sunwall_r_dir(l) * kwg1d(l,1)
                  sunwall_r_shadewall_dir(l) = sunwall_r_dir(l) * kww1d(l,1,1)
                  sunwall_r_br_tree_dir(l) = sunwall_r_dir(l) * kwv1d(l,1,1)
                  sunwall_r_ar_tree_dir(l) = sunwall_r_dir(l) * kwv1d(l,1,2)
                  err_r=sunwall_r_dir(l) -sunwall_r_br_tree_dir(l)*A_v1(l)/A_w(l) -sunwall_r_ar_tree_dir(l)*A_v2(l)/A_w(l) -sunwall_r_sky_dir(l)-sunwall_r_road_dir(l)*A_g(l)/A_w(l)-sunwall_r_shadewall_dir(l)*A_w(l)/A_w(l)
                  if (debug_write) then
                     write(6,*) '------sunwall--------l=---n=', l,n
                     write(6,*) 'err_r', err_r
                  end if                      
                  shadewall_r_sky_dir(l)     = shadewall_r_dir(l) * kws1d(l,1)
                  shadewall_r_road_dir(l)    = shadewall_r_dir(l) * kwg1d(l,1)
                  shadewall_r_sunwall_dir(l) = shadewall_r_dir(l) * kww1d(l,1,1)
                  shadewall_r_br_tree_dir(l) = shadewall_r_dir(l) * kwv1d(l,1,1)
                  shadewall_r_ar_tree_dir(l) = shadewall_r_dir(l) * kwv1d(l,1,2)
                  err_r=shadewall_r_dir(l) -shadewall_r_br_tree_dir(l)*A_v1(l)/A_w(l) -shadewall_r_ar_tree_dir(l)*A_v2(l)/A_w(l) -shadewall_r_sky_dir(l)-shadewall_r_road_dir(l)*A_g(l)/A_w(l)-shadewall_r_sunwall_dir(l)*A_w(l)/A_w(l)
                  if (debug_write) then
                     write(6,*) '-----shadewall--------l=---n=', l,n
                     write(6,*) 'err_r', err_r
                  end if   
                  br_tree_r_sky_dir(l)     = br_tree_r_dir(l) * kvs1d(l,1)
                  br_tree_r_road_dir(l)    = br_tree_r_dir(l) * kvg1d(l,1)
                  br_tree_r_sunwall_dir(l) = br_tree_r_dir(l) * kvw1d(l,1,1)
                  br_tree_r_shadewall_dir(l) = br_tree_r_dir(l) * kvw1d(l,1,1)   
                  br_tree_r_br_tree_dir(l) = br_tree_r_dir(l) * kvv1d(l,1,1)
                  br_tree_r_ar_tree_dir(l) = br_tree_r_dir(l) * kvv1d(l,1,2)
                  err_r=br_tree_r_dir(l) -br_tree_r_br_tree_dir(l) -br_tree_r_ar_tree_dir(l)*A_v2(l)/A_v1(l) -br_tree_r_sky_dir(l)-br_tree_r_road_dir(l)*A_g(l)/A_v1(l)-br_tree_r_sunwall_dir(l)*A_w(l)/A_v1(l)-br_tree_r_shadewall_dir(l)*A_w(l)/A_v1(l)
                  if (debug_write) then
                     write(6,*) '-----br_tree-------l=---n=', l,n
                     write(6,*) 'err_r',err_r
                  end if    
                  ar_tree_r_sky_dir(l)     = ar_tree_r_dir(l) * kvs1d(l,2)
                  ar_tree_r_road_dir(l)    = ar_tree_r_dir(l) * kvg1d(l,2)
                  ar_tree_r_sunwall_dir(l) = ar_tree_r_dir(l) * kvw1d(l,2,1)
                  ar_tree_r_shadewall_dir(l) = ar_tree_r_dir(l) * kvw1d(l,2,1)   
                  ar_tree_r_roof_dir(l) = ar_tree_r_dir(l) * kvr1d(l,2,2)
                  ar_tree_r_br_tree_dir(l) = ar_tree_r_dir(l) * kvv1d(l,2,1)
                  ar_tree_r_ar_tree_dir(l) = ar_tree_r_dir(l) * kvv1d(l,2,2)
                  !err_r=ar_tree_r_dir(l) -ar_tree_r_ar_tree_dir(l) -ar_tree_r_br_tree_dir(l)*A_v1(l)/A_v2(l) -ar_tree_r_sky_dir(l)-ar_tree_r_road_dir(l)*A_g(l)/A_v2(l)-ar_tree_r_sunwall_dir(l)*A_w(l)/A_v2(l)-ar_tree_r_shadewall_dir(l)*A_w(l)/A_v2(l)-ar_tree_r_roof_dir(l)*A_r(l)/A_v2(l)
                  !if (debug_write) then
                  !   write(6,*) '-----ar_tree------l=---n=', l,n
                  !   write(6,*) 'err_r', err_r  
                  !end if   
                  roof_r_sky_dir(l)       = roof_r_dir(l) * krs1d(l,2)
                  roof_r_ar_tree_dir(l) = roof_r_dir(l) * krv1d(l,2,2)
                  
                  err_r=roof_r_dir(l) -roof_r_ar_tree_dir(l)*A_v2(l)/A_r(l) -roof_r_sky_dir(l)
                  if (debug_write) then
                     write(6,*) '-----roof------l=---n=', l,n
                     write(6,*) 'err_r', err_r
                  end if   
                  ! step (4)

                  sref_improad_dir(l,ib)   = sref_improad_dir(l,ib) + improad_r_sky_dir(l)
                  sref_perroad_dir(l,ib)   = sref_perroad_dir(l,ib) + perroad_r_sky_dir(l)
                  sref_sunwall_dir(l,ib)   = sref_sunwall_dir(l,ib) + sunwall_r_sky_dir(l)
                  sref_shadewall_dir(l,ib) = sref_shadewall_dir(l,ib) + shadewall_r_sky_dir(l)
                  sref_br_tree_dir(l,ib)   = sref_br_tree_dir(l,ib) + br_tree_r_sky_dir(l) 
                  sref_ar_tree_dir(l,ib)   = sref_ar_tree_dir(l,ib) + ar_tree_r_sky_dir(l) 
                  sref_roof_dir(l,ib)   = sref_roof_dir(l,ib) + roof_r_sky_dir(l) 
                                
                 if (debug_write) then
                    write(6,*) '------------ Absorbed and Reflected SW Flux (sky view component only) ------------'
                    write(6,*) 'l =', l, ', ib = n=', ib,n

                    ! Reflected (to sky) direct shortwave
                    write(6,*) 'sref_improad_dir =', sref_improad_dir(l,ib)
                    write(6,*) 'sref_perroad_dir =', sref_perroad_dir(l,ib)
                    write(6,*) 'sref_sunwall_dir =', sref_sunwall_dir(l,ib)
                    write(6,*) 'sref_shadewall_dir =', sref_shadewall_dir(l,ib)
                    write(6,*) 'sref_br_tree_dir =', sref_br_tree_dir(l,ib)
                    write(6,*) 'sref_ar_tree_dir =', sref_ar_tree_dir(l,ib)
                    write(6,*) 'sref_roof_dir =', sref_roof_dir(l,ib) 

                    write(6,*) '----------------------------------------------------------------------------------'
                 end if                  
                  ! step (5)
                  if (debug_write) then
                     write(6,*) 'road_a_dir(l), sunwall_a_dir(l), shadewall_a_dir(l),br_tree_a_dir(l),ar_tree_a_dir(l),roof_a_dir(l)',road_a_dir(l), sunwall_a_dir(l), shadewall_a_dir(l),br_tree_a_dir(l),ar_tree_a_dir(l),roof_a_dir(l)
                  end if 
                  crit = max(road_a_dir(l), sunwall_a_dir(l), shadewall_a_dir(l),br_tree_a_dir(l),ar_tree_a_dir(l),roof_a_dir(l))
                  if (crit < errcrit) exit
               end do
               if (iter_dir >= n) then
                  write (iulog,*) 'urban net solar radiation error: no convergence, direct beam'
                  write (iulog,*) 'clm model is stopping'
                  call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
               endif

               ! reflected diffuse

               do iter_dif = 1, n
                  ! step (1)
                  !*canyon_hwr(l) no need to scale?
                  
                  stot(l) = (sunwall_r_road_dif(l) + shadewall_r_road_dif(l) + br_tree_r_road_dif(l)&
                            + ar_tree_r_road_dif(l))

                  road_a_dif(l) = 0.0_r8
                  road_r_dif(l) = 0.0_r8
                  improad_a_dif(l) = (1._r8-alb_improad_dif(l,ib)) * stot(l) 
                  improad_r_dif(l) =     alb_improad_dif(l,ib)  * stot(l) 
                  road_a_dif(l)    = road_a_dif(l) + improad_a_dif(l)*wtroad_imperv(l)
                  road_r_dif(l)    = road_r_dif(l) + improad_r_dif(l)*wtroad_imperv(l)
                  perroad_a_dif(l) = (1._r8-alb_perroad_dif(l,ib)) * stot(l) 
                  perroad_r_dif(l) =     alb_perroad_dif(l,ib)  * stot(l) 
                  road_a_dif(l)    = road_a_dif(l) + perroad_a_dif(l)*(wtroad_perv(l)+wtroad_tree(l))
                  road_r_dif(l)    = road_r_dif(l) + perroad_r_dif(l)*(wtroad_perv(l)+wtroad_tree(l))

                  stot(l) = road_r_sunwall_dif(l) + shadewall_r_sunwall_dif(l) + br_tree_r_sunwall_dif(l) &
                            + ar_tree_r_sunwall_dif(l)
                  sunwall_a_dif(l) = (1._r8-alb_wall_dif(l,ib)) * stot(l)
                  sunwall_r_dif(l) =     alb_wall_dif(l,ib)  * stot(l)

                  stot(l) = road_r_shadewall_dif(l) + sunwall_r_shadewall_dif(l) + br_tree_r_shadewall_dif(l) &
                            + ar_tree_r_shadewall_dif(l)
                  shadewall_a_dif(l) = (1._r8-alb_wall_dif(l,ib)) * stot(l)
                  shadewall_r_dif(l) =     alb_wall_dif(l,ib)  * stot(l)

                  stot(l) = road_r_br_tree_dif(l) + sunwall_r_br_tree_dif(l) + shadewall_r_br_tree_dif(l) &
                            + br_tree_r_br_tree_dif(l) +ar_tree_r_br_tree_dif(l)
                  br_tree_a_dif(l) = (1._r8-alb_br_tree_dif(l,ib)) * stot(l)
                  br_tree_r_dif(l) =     alb_br_tree_dif(l,ib)  * stot(l)
    
                  stot(l) = road_r_ar_tree_dif(l) + sunwall_r_ar_tree_dif(l) + shadewall_r_ar_tree_dif(l) &
                        + br_tree_r_ar_tree_dif(l) +ar_tree_r_ar_tree_dif(l) +  roof_r_ar_tree_dif(l)
                  ar_tree_a_dif(l) = (1._r8-alb_ar_tree_dif(l,ib)) * stot(l)
                  ar_tree_r_dif(l) =     alb_ar_tree_dif(l,ib)  * stot(l)


                  stot(l) =  ar_tree_r_roof_dif(l) 
                  roof_a_dif(l) = (1._r8-alb_roof_dif(l,ib)) * stot(l)
                  roof_r_dif(l) =     alb_roof_dif(l,ib)  * stot(l)
                                
                  ! step (2)

                  sabs_improad_dif(l,ib)   = sabs_improad_dif(l,ib)   + improad_a_dif(l)
                  sabs_perroad_dif(l,ib)   = sabs_perroad_dif(l,ib)   + perroad_a_dif(l)
                  sabs_sunwall_dif(l,ib)   = sabs_sunwall_dif(l,ib)   + sunwall_a_dif(l)
                  sabs_shadewall_dif(l,ib)   = sabs_shadewall_dif(l,ib)   + shadewall_a_dif(l)
                  sabs_br_tree_dif(l,ib) = sabs_br_tree_dif(l,ib) + br_tree_a_dif(l)
                  sabs_ar_tree_dif(l,ib) = sabs_ar_tree_dif(l,ib) + ar_tree_a_dif(l)
                  sabs_roof_dif(l,ib) = sabs_roof_dif(l,ib) + roof_a_dif(l)
                  
                 if (debug_write) then
                    write(6,*) '------------ Absorbed and Reflected SW Flux (sky view component only) ------------'
                    write(6,*) 'l =', l, ', ib = n =', ib,n

                    ! Absorbed diffuse shortwave
                    write(6,*) 'sabs_improad_dif =', sabs_improad_dif(l,ib)
                    write(6,*) 'sabs_perroad_dif =', sabs_perroad_dif(l,ib)
                    write(6,*) 'sabs_sunwall_dif =', sabs_sunwall_dif(l,ib)
                    write(6,*) 'sabs_shadewall_dif =', sabs_shadewall_dif(l,ib)
                    write(6,*) 'sabs_br_tree_dif =', sabs_br_tree_dif(l,ib)
                    write(6,*) 'sabs_ar_tree_dif =', sabs_ar_tree_dif(l,ib)
                    write(6,*) 'sabs_roof_dif =', sabs_roof_dif(l,ib)


                    write(6,*) '----------------------------------------------------------------------------------'
                 end if  
           
                  ! step (3)

                  improad_r_sky_dif(l)       = improad_r_dif(l) * kts1d(l)
                  improad_r_sunwall_dif(l)   = improad_r_dif(l) * kgw1d(l,1)
                  improad_r_shadewall_dif(l) = improad_r_dif(l) * kgw1d(l,1)
                  improad_r_br_tree_dif(l) = improad_r_dif(l) * kgv1d(l,1)
                  improad_r_ar_tree_dif(l) = improad_r_dif(l) * kgv1d(l,2)
                  err_r=improad_r_dif(l) -improad_r_br_tree_dif(l)*A_v1(l)/A_g(l) -improad_r_ar_tree_dif(l)*A_v2(l)/A_g(l) -improad_r_sky_dif(l)-improad_r_sunwall_dif(l)*A_w(l)/A_g(l)-improad_r_shadewall_dif(l)*A_w(l)/A_g(l)
                  if (debug_write) then
                     write(6,*) '-----improad--------l=---n=', l,n
                     write(6,*) 'err_r', err_r
                  end if   
                  perroad_r_sky_dif(l)       = perroad_r_dif(l) * kts1d(l)
                  perroad_r_sunwall_dif(l)   = perroad_r_dif(l) * kgw1d(l,1)
                  perroad_r_shadewall_dif(l) = perroad_r_dif(l) * kgw1d(l,1)
                  perroad_r_br_tree_dif(l) = perroad_r_dif(l) * kgv1d(l,1)
                  perroad_r_ar_tree_dif(l) = perroad_r_dif(l) * kgv1d(l,2)
                  err_r=perroad_r_dif(l) -perroad_r_br_tree_dif(l)*A_v1(l)/A_g(l) -perroad_r_ar_tree_dif(l)*A_v2(l)/A_g(l) -perroad_r_sky_dif(l)-perroad_r_sunwall_dif(l)*A_w(l)/A_g(l)-perroad_r_shadewall_dif(l)*A_w(l)/A_g (l)
                  if (debug_write) then
                     write(6,*) '-----perroad--------l=---n=', l,n
                     write(6,*) 'err_r', err_r
                  end if   
                  road_r_sky_dif(l)          = road_r_dif(l) * kts1d(l)
                  road_r_sunwall_dif(l)      = road_r_dif(l) * kgw1d(l,1)
                  road_r_shadewall_dif(l)    = road_r_dif(l) * kgw1d(l,1)
                  road_r_br_tree_dif(l)    = road_r_dif(l) * kgv1d(l,1)
                  road_r_ar_tree_dif(l)    = road_r_dif(l) * kgv1d(l,2)
                  err_r=road_r_dif(l) -road_r_br_tree_dif(l)*A_v1(l)/A_g(l) -road_r_ar_tree_dif(l)*A_v2(l)/A_g(l) -road_r_sky_dif(l)-road_r_sunwall_dif(l)*A_w(l)/A_g(l)-road_r_shadewall_dif(l)*A_w(l)/A_g(l)
                  if (debug_write) then
                     write(6,*) '-----road--------l=---n=', l,n
                     write(6,*) 'err_r', err_r
                  end if 
                  sunwall_r_sky_dif(l)       = sunwall_r_dif(l) * kws1d(l,1)
                  sunwall_r_road_dif(l)      = sunwall_r_dif(l) * kwg1d(l,1)
                  sunwall_r_shadewall_dif(l) = sunwall_r_dif(l) * kww1d(l,1,1)
                  sunwall_r_br_tree_dif(l) = sunwall_r_dif(l) * kwv1d(l,1,1)
                  sunwall_r_ar_tree_dif(l) = sunwall_r_dif(l) * kwv1d(l,1,2)
                  err_r=sunwall_r_dif(l) -sunwall_r_br_tree_dif(l)*A_v1(l)/A_w(l) -sunwall_r_ar_tree_dif(l)*A_v2(l)/A_w(l) -sunwall_r_sky_dif(l)-sunwall_r_road_dif(l)*A_g(l)/A_w(l)-sunwall_r_shadewall_dif(l)*A_w(l)/A_w(l)
                  if (debug_write) then
                     write(6,*) '------sunwall--------l=---n=', l,n
                     write(6,*) 'err_r', err_r
                  end if                     
                  shadewall_r_sky_dif(l)     = shadewall_r_dif(l) * kws1d(l,1)
                  shadewall_r_road_dif(l)    = shadewall_r_dif(l) * kwg1d(l,1)
                  shadewall_r_sunwall_dif(l) = shadewall_r_dif(l) * kww1d(l,1,1)
                  shadewall_r_br_tree_dif(l) = shadewall_r_dif(l) * kwv1d(l,1,1)
                  shadewall_r_ar_tree_dif(l) = shadewall_r_dif(l) * kwv1d(l,1,2)
                  err_r=shadewall_r_dif(l) -shadewall_r_br_tree_dif(l)*A_v1(l)/A_w(l) -shadewall_r_ar_tree_dif(l)*A_v2(l)/A_w(l) -shadewall_r_sky_dif(l)-shadewall_r_road_dif(l)*A_g(l)/A_w(l)-shadewall_r_sunwall_dif(l)*A_w(l)/A_w(l)
                  if (debug_write) then
                     write(6,*) '-----shadewall--------l=---n=', l,n
                     write(6,*) 'err_r', err_r
                  end if   
                  br_tree_r_sky_dif(l)     = br_tree_r_dif(l) * kvs1d(l,1)
                  br_tree_r_road_dif(l)    = br_tree_r_dif(l) * kvg1d(l,1)
                  br_tree_r_sunwall_dif(l) = br_tree_r_dif(l) * kvw1d(l,1,1)
                  br_tree_r_shadewall_dif(l) = br_tree_r_dif(l) * kvw1d(l,1,1)   
                  br_tree_r_br_tree_dif(l) = br_tree_r_dif(l) * kvv1d(l,1,1)
                  br_tree_r_ar_tree_dif(l) = br_tree_r_dif(l) * kvv1d(l,1,2)
                  err_r=br_tree_r_dif(l) -br_tree_r_br_tree_dif(l) -br_tree_r_ar_tree_dif(l)*A_v2(l)/A_v1(l) -br_tree_r_sky_dif(l)-br_tree_r_road_dif(l)*A_g(l)/A_v1(l)-br_tree_r_sunwall_dif(l)*A_w(l)/A_v1(l)-br_tree_r_shadewall_dif(l)*A_w(l)/A_v1(l)
                  if (debug_write) then
                     write(6,*) '-----br_tree-------l=---n=', l,n
                     write(6,*) 'err_r',err_r
                  end if     
                  ar_tree_r_sky_dif(l)     = ar_tree_r_dif(l) * kvs1d(l,2)
                  ar_tree_r_road_dif(l)    = ar_tree_r_dif(l) * kvg1d(l,2)
                  ar_tree_r_sunwall_dif(l) = ar_tree_r_dif(l) * kvw1d(l,2,1)
                  ar_tree_r_shadewall_dif(l) = ar_tree_r_dif(l) * kvw1d(l,2,1)   
                  ar_tree_r_roof_dif(l) = ar_tree_r_dif(l) * kvr1d(l,2,2)
                  ar_tree_r_br_tree_dif(l) = ar_tree_r_dif(l) * kvv1d(l,2,1)
                  ar_tree_r_ar_tree_dif(l) = ar_tree_r_dif(l) * kvv1d(l,2,2)
                  !err_r=ar_tree_r_dif(l) -ar_tree_r_ar_tree_dif(l) -ar_tree_r_br_tree_dif(l)*A_v1(l)/A_v2(l) -ar_tree_r_sky_dif(l)-ar_tree_r_road_dif(l)*A_g(l)/A_v2(l)-ar_tree_r_sunwall_dif(l)*A_w(l)/A_v2(l)-ar_tree_r_shadewall_dif(l)*A_w(l)/A_v2(l)-ar_tree_r_roof_dif(l)*A_r(l)/A_v2(l)
                  !if (debug_write) then
                  !   write(6,*) '-----ar_tree------l=---n=', l,n
                  !   write(6,*) 'err_r', err_r  
                  !end if  
                  roof_r_sky_dif(l)       = roof_r_dif(l) * krs1d(l,2)
                  roof_r_ar_tree_dif(l) = roof_r_dif(l) * krv1d(l,2,2)
                  err_r=roof_r_dif(l) -roof_r_ar_tree_dif(l)*A_v2(l)/A_r(l) -roof_r_sky_dif(l)
                  if (debug_write) then
                     write(6,*) '-----roof------l=---n=', l,n
                     write(6,*) 'err_r', err_r
                  end if   
                  ! step (4)

                  sref_improad_dif(l,ib)   = sref_improad_dif(l,ib) + improad_r_sky_dif(l)
                  sref_perroad_dif(l,ib)   = sref_perroad_dif(l,ib) + perroad_r_sky_dif(l)
                  sref_sunwall_dif(l,ib)   = sref_sunwall_dif(l,ib) + sunwall_r_sky_dif(l)
                  sref_shadewall_dif(l,ib) = sref_shadewall_dif(l,ib) + shadewall_r_sky_dif(l)
                  sref_br_tree_dif(l,ib)   = sref_br_tree_dif(l,ib) + br_tree_r_sky_dif(l) 
                  sref_ar_tree_dif(l,ib)   = sref_ar_tree_dif(l,ib) + ar_tree_r_sky_dif(l) 
                  sref_roof_dif(l,ib)   = sref_roof_dif(l,ib) + roof_r_sky_dif(l) 

                  
                 if (debug_write) then
                    write(6,*) '------------ Absorbed and Reflected SW Flux (sky view component only) ------------'
                    write(6,*) 'l =', l, ', ib = n=', ib,n

                    ! Reflected (to sky) diffuse shortwave
                    write(6,*) 'sref_improad_dif =', sref_improad_dif(l,ib)
                    write(6,*) 'sref_perroad_dif =', sref_perroad_dif(l,ib)
                    write(6,*) 'sref_sunwall_dif =', sref_sunwall_dif(l,ib)
                    write(6,*) 'sref_shadewall_dif =', sref_shadewall_dif(l,ib)
                    write(6,*) 'sref_br_tree_dif =', sref_br_tree_dif(l,ib)
                    write(6,*) 'sref_ar_tree_dif =', sref_ar_tree_dif(l,ib)
                    write(6,*) 'sref_roof_dif =', sref_roof_dif(l,ib) 

                    write(6,*) '----------------------------------------------------------------------------------'
                 end if
                  ! step (5)

                  crit = max(road_a_dif(l), sunwall_a_dif(l), shadewall_a_dif(l),br_tree_a_dif(l),ar_tree_a_dif(l),roof_a_dif(l))
                  if (crit < errcrit) exit
               end do

               if (iter_dif >= n) then
                  write (iulog,*) 'urban net solar radiation error: no convergence, diffuse'
                  write (iulog,*) 'clm model is stopping'
                  call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
               endif
               
               sref_tree_dir(l,ib)=(sref_br_tree_dir(l,ib)*A_v1(l) + sref_ar_tree_dir(l,ib)*A_v2(l))/(A_v2(l)+A_v1(l))
               sref_tree_dif(l,ib)=(sref_br_tree_dir(l,ib)*A_v1(l) + sref_ar_tree_dir(l,ib)*A_v2(l))/(A_v2(l)+A_v1(l))
               sabs_tree_dir(l,ib)=(sabs_br_tree_dir(l,ib)*A_v1(l) + sabs_ar_tree_dir(l,ib)*A_v2(l))/(A_v2(l)+A_v1(l))
               sabs_tree_dif(l,ib)=(sabs_br_tree_dir(l,ib)*A_v1(l) + sabs_ar_tree_dir(l,ib)*A_v2(l))/(A_v2(l)+A_v1(l))

               ! total reflected by canyon - sum of solar reflection to sky from canyon.
               ! project wall fluxes to horizontal surface
               
                sref_canyon_dir(l,ib) = 0.0_r8
                sref_canyon_dif(l,ib) = 0.0_r8
                
                sref_canyon_dir(l,ib) = sref_canyon_dir(l,ib) + (sref_improad_dir(l,ib)*wtroad_imperv(l)+sref_perroad_dir(l,ib)*(wtroad_perv(l)+wtroad_tree(l)))*A_g(l)/A_s(l) &
                                        +(sref_sunwall_dir(l,ib) + sref_shadewall_dir(l,ib))* A_w(l)/A_s(l) &
                                        + sref_br_tree_dir(l,ib)*A_v1(l)/A_s(l) + sref_ar_tree_dir(l,ib)*A_v2(l)/A_s(l)+ sref_roof_dir(l,ib)*A_r(l)/A_s(l)
                
                sref_canyon_dif(l,ib) = sref_canyon_dif(l,ib) + (sref_improad_dif(l,ib)*wtroad_imperv(l)+sref_perroad_dif(l,ib)*(wtroad_perv(l)+wtroad_tree(l)))*A_g(l)/A_s(l) &
                                        +(sref_sunwall_dif(l,ib) + sref_shadewall_dif(l,ib))* A_w(l)/A_s(l) &
                                        + sref_br_tree_dif(l,ib)*A_v1(l)/A_s(l) + sref_ar_tree_dif(l,ib)*A_v2(l)/A_s(l)+ sref_roof_dif(l,ib)*A_r(l)/A_s(l)
                
                ! total absorbed by canyon. project wall fluxes to horizontal surface
               
                sabs_canyon_dir(l,ib) = 0.0_r8
                sabs_canyon_dif(l,ib) = 0.0_r8
                sabs_canyon_dir(l,ib) = sabs_canyon_dir(l,ib) + (sabs_improad_dir(l,ib)*wtroad_imperv(l)+sabs_perroad_dir(l,ib)*(wtroad_perv(l)+wtroad_tree(l)))*A_g(l)/A_s(l) &
                                       +(sabs_sunwall_dir(l,ib) + sabs_shadewall_dir(l,ib))* A_w(l)/A_s(l) &
                                       + sabs_br_tree_dir(l,ib)*A_v1(l)/A_s(l) + sabs_ar_tree_dir(l,ib)*A_v2(l)/A_s(l)+ sabs_roof_dir(l,ib)*A_r(l)/A_s(l)
               
                sabs_canyon_dif(l,ib) = sabs_canyon_dif(l,ib) + (sabs_improad_dif(l,ib)*wtroad_imperv(l)+sabs_perroad_dif(l,ib)*(wtroad_perv(l)+wtroad_tree(l)))*A_g(l)/A_s(l) &
                                       +(sabs_sunwall_dif(l,ib) + sabs_shadewall_dif(l,ib))* A_w(l)/A_s(l) &
                                       + sabs_br_tree_dif(l,ib)*A_v1(l)/A_s(l) + sabs_ar_tree_dif(l,ib)*A_v2(l)/A_s(l)+ sabs_roof_dif(l,ib)*A_r(l)/A_s(l)
               
                 
                ! conservation check. note: previous conservation checks confirm partioning of total direct
                ! beam and diffuse radiation from atmosphere to road and walls is conserved as
                !    sdir (from atmosphere) = sdir_road + (sdir_sunwall + sdir_shadewall)*canyon_hwr
                !    sdif (from atmosphere) = sdif_road + (sdif_sunwall + sdif_shadewall)*canyon_hwr

                stot_dir(l)=sdir_road(l,ib)*A_g(l)/A_s(l)+(sdir_shadewall(l,ib) + sdir_sunwall(l,ib))* A_w(l)/A_s(l) &
                                    +sdir_roof(l,ib)*A_r(l)/A_s(l) + sdir_br_tree(l,ib)*A_v1(l)/A_s(l)+ sdir_ar_tree(l,ib)*A_v2(l)/A_s(l)
                stot_dif(l)=sdif_road(l,ib)*A_g(l)/A_s(l)+(sdif_shadewall(l,ib) + sdif_sunwall(l,ib))* A_w(l)/A_s(l) &
                                    +sdif_roof(l,ib)*A_r(l)/A_s(l) + sdif_br_tree(l,ib)*A_v1(l)/A_s(l)+ sdif_ar_tree(l,ib)*A_v2(l)/A_s(l) 

                err = stot_dir(l) + stot_dif(l) &
                     - (sabs_canyon_dir(l,ib) + sabs_canyon_dif(l,ib) + sref_canyon_dir(l,ib) + sref_canyon_dif(l,ib))
                if (debug_write) then
                  write(6,*) '--- Solar Absorption Diagnosis for l =', l, ', ib =', ib
                  write(6,*) 'A_g =', A_g(l), ', A_w =', A_w(l), ', A_r =', A_r(l), ', A_s =', A_s(l), ', A_v1 =', A_v1(l), ', A_v2 =', A_v2(l)
                  write(6,*) 'wtroad_imperv =', wtroad_imperv(l)

                  ! Direct
                  write(6,*) 'sabs_improad_dir =', sabs_improad_dir(l,ib)
                  write(6,*) 'sabs_sunwall_dir =', sabs_sunwall_dir(l,ib)
                  write(6,*) 'sabs_shadewall_dir =', sabs_shadewall_dir(l,ib)
                  write(6,*) 'sabs_roof_dir =', sabs_roof_dir(l,ib)
                  write(6,*) 'sabs_br_tree_dir =', sabs_br_tree_dir(l,ib)
                  write(6,*) 'sabs_ar_tree_dir =', sabs_ar_tree_dir(l,ib)

                  ! Diffuse
                  write(6,*) 'sabs_improad_dif =', sabs_improad_dif(l,ib)
                  write(6,*) 'sabs_sunwall_dif =', sabs_sunwall_dif(l,ib)
                  write(6,*) 'sabs_shadewall_dif =', sabs_shadewall_dif(l,ib)
                  write(6,*) 'sabs_roof_dif =', sabs_roof_dif(l,ib)
                  write(6,*) 'sabs_br_tree_dif =', sabs_br_tree_dif(l,ib)
                  write(6,*) 'sabs_ar_tree_dif =', sabs_ar_tree_dif(l,ib)
                  

                  ! Incoming fluxes (sdir/sdif terms)
                  write(6,*) 'sdir_road =', sdir_road(l,ib)
                  write(6,*) 'sdir_sunwall =', sdir_sunwall(l,ib)
                  write(6,*) 'sdir_shadewall =', sdir_shadewall(l,ib)
                  write(6,*) 'sdir_roof =', sdir_roof(l,ib)
                  write(6,*) 'sdir_br_tree =', sdir_br_tree(l,ib)
                  write(6,*) 'sdir_ar_tree =', sdir_ar_tree(l,ib)
                  write(6,*) 'stot_dir =', stot_dir(l)

                  write(6,*) 'sdif_road =', sdif_road(l,ib)
                  write(6,*) 'sdif_sunwall =', sdif_sunwall(l,ib)
                  write(6,*) 'sdif_shadewall =', sdif_shadewall(l,ib)
                  write(6,*) 'sdif_roof =', sdif_roof(l,ib)
                  write(6,*) 'sdif_br_tree =', sdif_br_tree(l,ib)
                  write(6,*) 'sdif_ar_tree =', sdif_ar_tree(l,ib)
                  

                  ! Reflected components
                  write(6,*) 'sabs_canyon_dir =', sabs_canyon_dir(l,ib)
                  write(6,*) 'sref_canyon_dir =', sref_canyon_dir(l,ib)
                  write(6,*) 'stot_dir =', stot_dir(l)
                  write(6,*) 'sabs_canyon_dif =', sabs_canyon_dif(l,ib)
                  write(6,*) 'sref_canyon_dif =', sref_canyon_dif(l,ib)               
                  write(6,*) 'stot_dif =', stot_dif(l)
                  ! Conservation error
                  write(6,*) 'err1 =', stot_dir(l) - (sabs_canyon_dir(l,ib) + sref_canyon_dir(l,ib) )
                  write(6,*) 'err2 =', stot_dif(l) - (sabs_canyon_dif(l,ib) + sref_canyon_dif(l,ib) )
                  write(6,*) 'err =', err
               end if
               !write(6,*) 'err (direct) =', stot_dir(l) - (sabs_canyon_dir(l,ib) + sref_canyon_dir(l,ib) )
               !write(6,*) 'err (diffuse) =', stot_dif(l) - (sabs_canyon_dif(l,ib) + sref_canyon_dif(l,ib) )              
               !write(6,*) 'err =', err
                !if (abs(err) > 0.001_r8 ) then
               !    write(iulog,*)'urban net solar radiation balance error for ib=',ib,' err= ',err
               !    write(iulog,*)' l= ',l,' ib= ',ib 
               !    write(iulog,*)' stot_dir        = ',stot_dir(l)
               !    write(iulog,*)' stot_dif        = ',stot_dif(l)
               !    write(iulog,*)' sabs_canyon_dir = ',sabs_canyon_dir(l,ib)
               !    write(iulog,*)' sabs_canyon_dif = ',sabs_canyon_dif(l,ib)
               !    write(iulog,*)' sref_canyon_dir = ',sref_canyon_dir(l,ib)
               !    write(iulog,*)' sref_canyon_dif = ',sref_canyon_dir(l,ib)
               !    write(iulog,*) 'clm model is stopping'
               !    call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
                !endif

                ! canyon albedo

                canyon_alb_dif(l) = sref_canyon_dif(l,ib) / max(stot_dif(l), 1.e-06_r8)
                canyon_alb_dir(l) = sref_canyon_dir(l,ib) / max(stot_dir(l), 1.e-06_r8)
             end if 


         end do   ! end of landunit loop

      end do   ! end of radiation band loop

    end associate

  end subroutine net_solar

  
end module UrbanAlbedoMod
