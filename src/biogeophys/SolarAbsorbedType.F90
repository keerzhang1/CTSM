module SolarAbsorbedType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_log_mod  , only: errMsg => shr_log_errMsg
  use decompMod    , only : bounds_type
  use clm_varcon   , only : spval
  use clm_varctl   , only : use_luna
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC DATA MEMBERS:
  type, public :: solarabs_type

     ! Solar reflected
     real(r8), pointer :: fsr_patch              (:)   ! patch solar radiation reflected (W/m**2)
     real(r8), pointer :: fsrSF_patch            (:)   ! diagnostic snow-free patch solar radiation reflected (W/m**2)
     real(r8), pointer :: ssre_fsr_patch         (:)   ! snow radiative effect on patch solar radiation reflected (W/m**2)
     ! Solar Absorbed
     real(r8), pointer :: fsa_patch              (:)   ! patch solar radiation absorbed (total) (W/m**2)  
     real(r8), pointer :: fsa_u_patch            (:)   ! patch urban solar radiation absorbed (total) (W/m**2)
     real(r8), pointer :: fsa_r_patch            (:)   ! patch rural solar radiation absorbed (total) (W/m**2)
     real(r8), pointer :: parsun_z_patch         (:,:) ! patch absorbed PAR for sunlit leaves in canopy layer (W/m**2) 
     real(r8), pointer :: parsha_z_patch         (:,:) ! patch absorbed PAR for shaded leaves in canopy layer (W/m**2) 
     real(r8), pointer :: par240d_z_patch        (:,:) ! 10-day running mean of daytime patch absorbed PAR for leaves in canopy layer (W/m**2) 
     real(r8), pointer :: par240x_z_patch        (:,:) ! 10-day running mean of maximum patch absorbed PAR for leaves in canopy layer (W/m**2)
     real(r8), pointer :: par24d_z_patch         (:,:) ! daily accumulated  absorbed PAR for leaves in canopy layer from midnight to current step(J/m**2) 
     real(r8), pointer :: par24x_z_patch         (:,:) ! daily max of patch absorbed PAR for  leaves in canopy layer from midnight to current step(W/m**2)
     real(r8), pointer :: sabg_soil_patch        (:)   ! patch solar radiation absorbed by soil (W/m**2)       
     real(r8), pointer :: sabg_snow_patch        (:)   ! patch solar radiation absorbed by snow (W/m**2)       
     real(r8), pointer :: sabg_patch             (:)   ! patch solar radiation absorbed by ground (W/m**2)     
     real(r8), pointer :: sabg_chk_patch         (:)   ! patch fsno weighted sum (W/m**2)                                   
     real(r8), pointer :: sabg_lyr_patch         (:,:) ! patch absorbed radiation in each snow layer and top soil layer (pft,lyr) [W/m2]
     real(r8), pointer :: sabg_pen_patch         (:)   ! patch (rural) shortwave radiation penetrating top soisno layer [W/m2]

     real(r8), pointer :: sub_surf_abs_SW_patch  (:)   ! patch fraction of solar radiation absorbed below first snow layer
     real(r8), pointer :: sabv_patch             (:)   ! patch solar radiation absorbed by vegetation (W/m**2) 

     real(r8), pointer :: sdir_road_lun   (:,:) ! lun direct beam solar radiation incident on road per unit incident flux  [landunit, numrad]
     real(r8), pointer :: sdir_sunwall_lun   (:,:) ! lun direct beam solar radiation incident on sun-lit wall per unit incident flux  [landunit, numrad]
     real(r8), pointer :: sdir_shadewall_lun   (:,:) ! lun direct beam solar radiation incident on shaded wall per unit incident flux  [landunit, numrad]
     real(r8), pointer :: sdir_roof_lun   (:,:) ! lun direct beam solar radiation incident on roof per unit incident flux  [landunit, numrad]
     real(r8), pointer :: sdir_ar_tree_lun   (:,:) ! lun direct beam solar radiation incident on above-roof tree per unit incident flux  [landunit, numrad]
     real(r8), pointer :: sdir_br_tree_lun   (:,:) ! lun direct beam solar radiation incident on below-roof tree per unit incident flux  [landunit, numrad]
     real(r8), pointer :: sdir_tree_lun   (:,:) ! lun direct beam solar radiation incident on tree per unit incident flux  [landunit, numrad]
     real(r8), pointer :: sdir_force_lun   (:,:) ! lun direct beam radiation  (vis=forc_sols , nir=forc_soll ) for diagnosis

     real(r8), pointer :: sabs_roof_dir_lun       (:,:) ! lun direct solar absorbed by roof per unit ground area per unit incident flux
     real(r8), pointer :: sabs_roof_dif_lun       (:,:) ! lun diffuse solar absorbed by roof per unit ground area per unit incident flux
     real(r8), pointer :: sabs_sunwall_dir_lun    (:,:) ! lun direct solar absorbed by sunwall per unit wall area per unit incident flux
     real(r8), pointer :: sabs_sunwall_dif_lun    (:,:) ! lun diffuse solar absorbed by sunwall per unit wall area per unit incident flux
     real(r8), pointer :: sabs_shadewall_dir_lun  (:,:) ! lun direct solar absorbed by shadewall per unit wall area per unit incident flux
     real(r8), pointer :: sabs_shadewall_dif_lun  (:,:) ! lun diffuse solar absorbed by shadewall per unit wall area per unit incident flux
     real(r8), pointer :: sabs_improad_dir_lun    (:,:) ! lun direct solar absorbed by impervious road per unit ground area per unit incident flux
     real(r8), pointer :: sabs_improad_dif_lun    (:,:) ! lun diffuse solar absorbed by impervious road per unit ground area per unit incident flux
     real(r8), pointer :: sabs_perroad_dir_lun    (:,:) ! lun direct solar absorbed by pervious road per unit ground area per unit incident flux
     real(r8), pointer :: sabs_perroad_dif_lun    (:,:) ! lun diffuse solar absorbed by pervious road per unit ground area per unit incident flux
     real(r8), pointer :: sabs_br_tree_dir_lun      (:,:) ! lun direct solar absorbed by below-roof treeetation per unit treeetation area per unit incident flux
     real(r8), pointer :: sabs_br_tree_dif_lun      (:,:) ! lun diffuse solar absorbed by below-roof treeetation per unit treeetation area per unit incident flux
     real(r8), pointer :: sabs_ar_tree_dir_lun      (:,:) ! lun direct solar absorbed by above-roof treeetation per unit treeetation area per unit incident flux
     real(r8), pointer :: sabs_ar_tree_dif_lun      (:,:) ! lun diffuse solar absorbed by above-roof treeetation per unit treeetation area per unit incident flux
     real(r8), pointer :: sabs_tree_dir_lun   (:,:) ! lun direct  solar absorbed by road tree per unit ground area per unit incident flux
     real(r8), pointer :: sabs_tree_dif_lun   (:,:) ! lun diffuse solar absorbed by road tree per unit ground area per unit incident flux
     real(r8), pointer :: sabs_canyon_dir_lun   (:,:) ! lun direct solar absorbed by impervious road per unit ground area per unit incident flux
     real(r8), pointer :: sabs_canyon_dif_lun   (:,:) ! lun diffuse solar absorbed by impervious road per unit ground area per unit incident flux
     real(r8), pointer :: sref_roof_dir_lun  (:,:) ! lun direct solar reflected by roof per unit ground area per unit incident flux
     real(r8), pointer :: sref_roof_dif_lun  (:,:) ! lun diffuse solar reflected by roof per unit ground area per unit incident flux
     real(r8), pointer :: sref_sunwall_dir_lun    (:,:) ! lun diffuse solar reflected by sunwall per unit wall area per unit incident flux
     real(r8), pointer :: sref_sunwall_dif_lun    (:,:) ! lun diffuse solar reflected by sunwall per unit wall area per unit incident flux
     real(r8), pointer :: sref_shadewall_dir_lun  (:,:) ! lun direct solar reflected by shadewall per unit wall area per unit incident flux
     real(r8), pointer :: sref_shadewall_dif_lun  (:,:) ! lun diffuse solar reflected by shadewall per unit wall area per unit incident flux
     real(r8), pointer :: sref_improad_dir_lun    (:,:) ! lun direct solar reflected by impervious road per unit ground area per unit incident flux
     real(r8), pointer :: sref_improad_dif_lun    (:,:) ! lun diffuse solar reflected by impervious road per unit ground area per unit incident flux
     real(r8), pointer :: sref_perroad_dir_lun    (:,:) ! lun direct solar reflected by pervious road per unit ground area per unit incident flux
     real(r8), pointer :: sref_perroad_dif_lun    (:,:) ! lun diffuse solar reflected by pervious road per unit ground area per unit incident flux
     real(r8), pointer :: sref_br_tree_dir_lun    (:,:) ! lun direct solar reflected by below-roof tree per unit ground area per unit incident flux
     real(r8), pointer :: sref_br_tree_dif_lun    (:,:) ! lun diffuse solar reflected by below-roof tree per unit ground area per unit incident flux
     real(r8), pointer :: sref_ar_tree_dir_lun    (:,:) ! lun direct solar reflected by above-roof tree per unit ground area per unit incident flux
     real(r8), pointer :: sref_ar_tree_dif_lun    (:,:) ! lun diffuse solar reflected by above-roof tree per unit ground area per unit incident flux
     real(r8), pointer :: sref_tree_dir_lun   (:,:) ! lun direct  solar reflected by road tree per unit ground area per unit incident flux
     real(r8), pointer :: sref_tree_dif_lun   (:,:) ! lun diffuse solar reflected by road tree per unit ground area per unit incident flux
     real(r8), pointer :: sref_canyon_dir_lun     (:,:) ! lun direct solar reflected by canyon per unit ground area per unit incident flux
     real(r8), pointer :: sref_canyon_dif_lun     (:,:) ! lun diffuse solar reflected by canyon per unit ground area per unit incident flux

    real(r8), pointer :: albd_perroad_s_lun     (:,:) ! ground direct solar albedo on pervious road with snow effects
    real(r8), pointer :: albi_perroad_s_lun     (:,:) ! ground diffuse solar albedo on pervious road with snow effects

    real(r8), pointer :: lwnet_roof_lun   (:) ! lun net longwave flux at shaded roof
    real(r8), pointer :: lwnet_improad_lun     (:) ! lun net longwave flux at impervious road
    real(r8), pointer :: lw_emi_improad_lun     (:) ! lun emitted longwave flux at impervious road
    real(r8), pointer :: lwnet_perroad_lun     (:) ! lun net longwave flux at pervious road
    real(r8), pointer :: lwnet_tree_perroad_lun     (:) ! lun net longwave flux at pervious road beneath the tree canopy

    real(r8), pointer :: lwnet_sunwall_lun     (:) ! lun net longwave flux at sunlit wall
    real(r8), pointer :: lwnet_shadewall_lun   (:) ! lun net longwave flux at shaded wall
    real(r8), pointer :: lwnet_br_tree_lun     (:) ! lun net longwave flux at below-roof tree
    real(r8), pointer :: lwnet_ar_tree_lun     (:) ! lun net longwave flux at above-roof tree
    real(r8), pointer :: lwnet_tree_lun     (:) ! lun net longwave flux at above-roof tree
    real(r8), pointer :: lwnet_canyon_lun      (:) ! lun net longwave flux at canyon center
    real(r8), pointer :: lwup_roof_lun    (:) ! lun upward longwave flux at shaded roof
    real(r8), pointer :: lwup_improad_lun      (:) ! lun upward longwave flux at impervious road
    real(r8), pointer :: lwup_perroad_lun      (:) ! lun upward longwave flux at pervious road
    real(r8), pointer :: lwup_tree_perroad_lun      (:) ! lun upward longwave flux at pervious road beneath the tree canopy
    real(r8), pointer :: lwup_sunwall_lun      (:) ! lun upward longwave flux at sunlit wall
    real(r8), pointer :: lwup_shadewall_lun    (:) ! lun upward longwave flux at shaded wall
    real(r8), pointer :: lwup_br_tree_lun      (:) ! lun upward longwave flux at below-roof tree
    real(r8), pointer :: lwup_ar_tree_lun      (:) ! lun upward longwave flux at above-roof tree
    real(r8), pointer :: lwup_tree_lun      (:) ! lun upward longwave flux at above-roof tree
    real(r8), pointer :: lwup_canyon_lun       (:) ! lun upward longwave flux at canyon center
    real(r8), pointer :: lwdown_road_lun       (:) ! lun downward longwave flux at road
    real(r8), pointer :: lwdown_sunwall_lun       (:) ! lun downward longwave flux at sunlit wall
    real(r8), pointer :: lwdown_shadewall_lun       (:) ! lun downward longwave flux at shaded wall
    real(r8), pointer :: lwdown_roof_lun       (:) ! lun downward longwave flux at roof
    real(r8), pointer :: lwdown_br_tree_lun       (:) ! lun downward longwave flux at below-roof tree
    real(r8), pointer :: lwdown_ar_tree_lun       (:) ! lun downward longwave flux at above-roof tree

     ! Currently needed by lake code 
     ! TODO (MV 8/20/2014) should be moved in the future
     real(r8), pointer :: fsds_nir_d_patch       (:)   ! patch incident direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: fsds_nir_i_patch       (:)   ! patch incident diffuse nir solar radiation (W/m**2)
     real(r8), pointer :: fsds_nir_d_ln_patch    (:)   ! patch incident direct beam nir solar radiation at local noon (W/m**2)
     real(r8), pointer :: fsr_nir_d_patch        (:)   ! patch reflected direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: fsr_nir_i_patch        (:)   ! patch reflected diffuse nir solar radiation (W/m**2)
     real(r8), pointer :: fsr_nir_d_ln_patch     (:)   ! patch reflected direct beam nir solar radiation at local noon (W/m**2)
     ! optional diagnostic fluxes:
     real(r8), pointer :: fsrSF_nir_d_patch      (:)   ! snow-free patch reflected direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: fsrSF_nir_i_patch      (:)   ! snow-free patch reflected diffuse nir solar radiation (W/m**2)
     real(r8), pointer :: fsrSF_nir_d_ln_patch   (:)   ! snow-free patch reflected direct beam nir solar radiation at local noon (W/m**2)
     real(r8), pointer :: ssre_fsr_nir_d_patch   (:)   ! snow-free patch reflected direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: ssre_fsr_nir_i_patch   (:)   ! snow-free patch reflected diffuse nir solar radiation (W/m**2)
     real(r8), pointer :: ssre_fsr_nir_d_ln_patch(:)   ! snow-free patch reflected direct beam nir solar radiation at local noon (W/m**2)

   contains

     procedure, public  :: Init         
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     
     procedure, public  :: Restart      

  end type solarabs_type
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(solarabs_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar    , only : nlevcan, nlevcan, numrad, nlevsno
    !
    ! !ARGUMENTS:
    class(solarabs_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begl = bounds%begl; endl = bounds%endl

    allocate(this%fsa_patch              (begp:endp))              ; this%fsa_patch              (:)   = nan
    allocate(this%fsa_u_patch            (begp:endp))              ; this%fsa_u_patch            (:)   = nan
    allocate(this%fsa_r_patch            (begp:endp))              ; this%fsa_r_patch            (:)   = nan
    allocate(this%parsun_z_patch         (begp:endp,1:nlevcan))    ; this%parsun_z_patch         (:,:) = nan
    allocate(this%parsha_z_patch         (begp:endp,1:nlevcan))    ; this%parsha_z_patch         (:,:) = nan 
    if(use_luna)then
        allocate(this%par240d_z_patch    (begp:endp,1:nlevcan))    ; this%par240d_z_patch        (:,:) = spval
        allocate(this%par240x_z_patch    (begp:endp,1:nlevcan))    ; this%par240x_z_patch        (:,:) = spval 
        allocate(this%par24d_z_patch     (begp:endp,1:nlevcan))    ; this%par24d_z_patch         (:,:) = spval
        allocate(this%par24x_z_patch     (begp:endp,1:nlevcan))    ; this%par24x_z_patch         (:,:) = spval
    endif    
    allocate(this%sabv_patch             (begp:endp))              ; this%sabv_patch             (:)   = nan
    allocate(this%sabg_patch             (begp:endp))              ; this%sabg_patch             (:)   = nan
    allocate(this%sabg_lyr_patch         (begp:endp,-nlevsno+1:1)) ; this%sabg_lyr_patch         (:,:) = nan
    allocate(this%sabg_pen_patch         (begp:endp))              ; this%sabg_pen_patch         (:)   = nan
    allocate(this%sabg_soil_patch        (begp:endp))              ; this%sabg_soil_patch        (:)   = nan
    allocate(this%sabg_snow_patch        (begp:endp))              ; this%sabg_snow_patch        (:)   = nan
    allocate(this%sabg_chk_patch         (begp:endp))              ; this%sabg_chk_patch         (:)   = nan
    allocate(this%sabs_roof_dir_lun      (begl:endl,1:numrad))     ; this%sabs_roof_dir_lun      (:,:) = nan
    allocate(this%sabs_roof_dif_lun      (begl:endl,1:numrad))     ; this%sabs_roof_dif_lun      (:,:) = nan
    allocate(this%sabs_sunwall_dir_lun   (begl:endl,1:numrad))     ; this%sabs_sunwall_dir_lun   (:,:) = nan
    allocate(this%sabs_sunwall_dif_lun   (begl:endl,1:numrad))     ; this%sabs_sunwall_dif_lun   (:,:) = nan
    allocate(this%sabs_shadewall_dir_lun (begl:endl,1:numrad))     ; this%sabs_shadewall_dir_lun (:,:) = nan
    allocate(this%sabs_shadewall_dif_lun (begl:endl,1:numrad))     ; this%sabs_shadewall_dif_lun (:,:) = nan
    allocate(this%sabs_improad_dir_lun   (begl:endl,1:numrad))     ; this%sabs_improad_dir_lun   (:,:) = nan
    allocate(this%sabs_improad_dif_lun   (begl:endl,1:numrad))     ; this%sabs_improad_dif_lun   (:,:) = nan
    allocate(this%sabs_perroad_dir_lun   (begl:endl,1:numrad))     ; this%sabs_perroad_dir_lun   (:,:) = nan
    allocate(this%sabs_perroad_dif_lun   (begl:endl,1:numrad))     ; this%sabs_perroad_dif_lun   (:,:) = nan 
    allocate(this%sabs_tree_dir_lun   (begl:endl,1:numrad))     ; this%sabs_tree_dir_lun   (:,:) = nan
    allocate(this%sabs_tree_dif_lun   (begl:endl,1:numrad))     ; this%sabs_tree_dif_lun   (:,:) = nan

    allocate(this%sdir_road_lun (begl:endl,1:numrad))     ; this%sdir_road_lun(:,:) = nan
    allocate(this%sdir_sunwall_lun (begl:endl,1:numrad))     ; this%sdir_sunwall_lun (:,:) = nan
    allocate(this%sdir_shadewall_lun (begl:endl,1:numrad))     ; this%sdir_shadewall_lun (:,:) = nan
    allocate(this%sdir_roof_lun (begl:endl,1:numrad))     ; this%sdir_roof_lun (:,:) = nan
    allocate(this%albd_perroad_s_lun     (begl:endl,1:numrad))     ; this%albd_perroad_s_lun     (:,:) = nan
    allocate(this%albi_perroad_s_lun     (begl:endl,1:numrad))     ; this%albi_perroad_s_lun     (:,:) = nan
    allocate(this%sdir_ar_tree_lun (begl:endl,1:numrad))     ; this%sdir_ar_tree_lun (:,:) = nan
    allocate(this%sdir_br_tree_lun (begl:endl,1:numrad))     ; this%sdir_br_tree_lun (:,:) = nan
    allocate(this%sdir_tree_lun (begl:endl,1:numrad))     ; this%sdir_tree_lun (:,:) = nan
    allocate(this%sdir_force_lun (begl:endl,1:numrad))     ; this%sdir_force_lun (:,:) = nan

    allocate(this%sub_surf_abs_SW_patch  (begp:endp))              ; this%sub_surf_abs_SW_patch  (:)   = nan
    allocate(this%fsr_patch              (begp:endp))              ; this%fsr_patch              (:)   = nan
    allocate(this%fsr_nir_d_patch        (begp:endp))              ; this%fsr_nir_d_patch        (:)   = nan
    allocate(this%fsr_nir_i_patch        (begp:endp))              ; this%fsr_nir_i_patch        (:)   = nan
    allocate(this%fsr_nir_d_ln_patch     (begp:endp))              ; this%fsr_nir_d_ln_patch     (:)   = nan
    allocate(this%fsrSF_patch            (begp:endp))              ; this%fsrSF_patch            (:)   = nan
    allocate(this%fsrSF_nir_d_patch      (begp:endp))              ; this%fsrSF_nir_d_patch      (:)   = nan
    allocate(this%fsrSF_nir_i_patch      (begp:endp))              ; this%fsrSF_nir_i_patch      (:)   = nan
    allocate(this%fsrSF_nir_d_ln_patch   (begp:endp))              ; this%fsrSF_nir_d_ln_patch   (:)   = nan
    allocate(this%ssre_fsr_patch         (begp:endp))              ; this%ssre_fsr_patch         (:)   = nan
    allocate(this%ssre_fsr_nir_d_patch   (begp:endp))              ; this%ssre_fsr_nir_d_patch   (:)   = nan
    allocate(this%ssre_fsr_nir_i_patch   (begp:endp))              ; this%ssre_fsr_nir_i_patch   (:)   = nan
    allocate(this%ssre_fsr_nir_d_ln_patch(begp:endp))              ; this%ssre_fsr_nir_d_ln_patch(:)   = nan
    allocate(this%fsds_nir_d_patch       (begp:endp))              ; this%fsds_nir_d_patch       (:)   = nan
    allocate(this%fsds_nir_i_patch       (begp:endp))              ; this%fsds_nir_i_patch       (:)   = nan
    allocate(this%fsds_nir_d_ln_patch    (begp:endp))              ; this%fsds_nir_d_ln_patch    (:)   = nan

    allocate(this%sabs_roof_dir_lun(begl:endl,1:numrad))         ; this%sabs_roof_dir_lun(:,:) = nan
    allocate(this%sabs_roof_dif_lun(begl:endl,1:numrad))         ; this%sabs_roof_dif_lun(:,:) = nan
    allocate(this%sabs_sunwall_dir_lun(begl:endl,1:numrad))      ; this%sabs_sunwall_dir_lun(:,:) = nan
    allocate(this%sabs_sunwall_dif_lun(begl:endl,1:numrad))      ; this%sabs_sunwall_dif_lun(:,:) = nan
    allocate(this%sabs_shadewall_dir_lun(begl:endl,1:numrad))    ; this%sabs_shadewall_dir_lun(:,:) = nan
    allocate(this%sabs_shadewall_dif_lun(begl:endl,1:numrad))    ; this%sabs_shadewall_dif_lun(:,:) = nan
    allocate(this%sabs_improad_dir_lun(begl:endl,1:numrad))      ; this%sabs_improad_dir_lun(:,:) = nan
    allocate(this%sabs_improad_dif_lun(begl:endl,1:numrad))      ; this%sabs_improad_dif_lun(:,:) = nan
    allocate(this%sabs_perroad_dir_lun(begl:endl,1:numrad))      ; this%sabs_perroad_dir_lun(:,:) = nan
    allocate(this%sabs_perroad_dif_lun(begl:endl,1:numrad))      ; this%sabs_perroad_dif_lun(:,:) = nan

    allocate(this%sabs_br_tree_dir_lun(begl:endl,1:numrad))        ; this%sabs_br_tree_dir_lun(:,:) = nan
    allocate(this%sabs_br_tree_dif_lun(begl:endl,1:numrad))        ; this%sabs_br_tree_dif_lun(:,:) = nan
    allocate(this%sabs_ar_tree_dir_lun(begl:endl,1:numrad))        ; this%sabs_ar_tree_dir_lun(:,:) = nan
    allocate(this%sabs_ar_tree_dif_lun(begl:endl,1:numrad))        ; this%sabs_ar_tree_dif_lun(:,:) = nan
    allocate(this%sabs_tree_dir_lun(begl:endl,1:numrad))        ; this%sabs_tree_dir_lun(:,:) = nan
    allocate(this%sabs_tree_dif_lun(begl:endl,1:numrad))        ; this%sabs_tree_dif_lun(:,:) = nan
    allocate(this%sabs_canyon_dir_lun(begl:endl,1:numrad))     ; this%sabs_canyon_dir_lun(:,:) = nan
    allocate(this%sabs_canyon_dif_lun(begl:endl,1:numrad))     ; this%sabs_canyon_dif_lun(:,:) = nan
    allocate(this%sref_roof_dir_lun(begl:endl,1:numrad))    ; this%sref_roof_dir_lun(:,:) = nan
    allocate(this%sref_roof_dif_lun(begl:endl,1:numrad))    ; this%sref_roof_dif_lun(:,:) = nan
    allocate(this%sref_sunwall_dir_lun(begl:endl,1:numrad))      ; this%sref_sunwall_dir_lun(:,:) = nan
    allocate(this%sref_sunwall_dif_lun(begl:endl,1:numrad))      ; this%sref_sunwall_dif_lun(:,:) = nan
    allocate(this%sref_shadewall_dir_lun(begl:endl,1:numrad))    ; this%sref_shadewall_dir_lun(:,:) = nan
    allocate(this%sref_shadewall_dif_lun(begl:endl,1:numrad))    ; this%sref_shadewall_dif_lun(:,:) = nan
    allocate(this%sref_improad_dir_lun(begl:endl,1:numrad))      ; this%sref_improad_dir_lun(:,:) = nan
    allocate(this%sref_improad_dif_lun(begl:endl,1:numrad))      ; this%sref_improad_dif_lun(:,:) = nan
    allocate(this%sref_perroad_dir_lun(begl:endl,1:numrad))      ; this%sref_perroad_dir_lun(:,:) = nan
    allocate(this%sref_perroad_dif_lun(begl:endl,1:numrad))      ; this%sref_perroad_dif_lun(:,:) = nan

    allocate(this%sref_br_tree_dir_lun(begl:endl,1:numrad))      ; this%sref_br_tree_dir_lun(:,:) = nan
    allocate(this%sref_br_tree_dif_lun(begl:endl,1:numrad))      ; this%sref_br_tree_dif_lun(:,:) = nan
    allocate(this%sref_ar_tree_dir_lun(begl:endl,1:numrad))      ; this%sref_ar_tree_dir_lun(:,:) = nan
    allocate(this%sref_ar_tree_dif_lun(begl:endl,1:numrad))      ; this%sref_ar_tree_dif_lun(:,:) = nan
    allocate(this%sref_tree_dir_lun(begl:endl,1:numrad))      ; this%sref_tree_dir_lun(:,:) = nan
    allocate(this%sref_tree_dif_lun(begl:endl,1:numrad))      ; this%sref_tree_dif_lun(:,:) = nan
    allocate(this%sref_canyon_dir_lun(begl:endl,1:numrad))       ; this%sref_canyon_dir_lun(:,:) = nan
    allocate(this%sref_canyon_dif_lun(begl:endl,1:numrad))       ; this%sref_canyon_dif_lun(:,:) = nan
    
    allocate(this%lwnet_roof_lun(begl:endl))       ; this%lwnet_roof_lun(:) = nan
    allocate(this%lwnet_improad_lun(begl:endl))         ; this%lwnet_improad_lun(:) = nan
    allocate(this%lw_emi_improad_lun(begl:endl))         ; this%lw_emi_improad_lun(:) = nan
    allocate(this%lwnet_perroad_lun(begl:endl))         ; this%lwnet_perroad_lun(:) = nan
    allocate(this%lwnet_tree_perroad_lun(begl:endl))         ; this%lwnet_tree_perroad_lun(:) = nan

    allocate(this%lwnet_sunwall_lun(begl:endl))         ; this%lwnet_sunwall_lun(:) = nan
    allocate(this%lwnet_shadewall_lun(begl:endl))       ; this%lwnet_shadewall_lun(:) = nan
    allocate(this%lwnet_br_tree_lun(begl:endl))         ; this%lwnet_br_tree_lun(:) = nan
    allocate(this%lwnet_ar_tree_lun(begl:endl))         ; this%lwnet_ar_tree_lun(:) = nan
    allocate(this%lwnet_tree_lun(begl:endl))         ; this%lwnet_tree_lun(:) = nan
    allocate(this%lwnet_canyon_lun(begl:endl))          ; this%lwnet_canyon_lun(:) = nan
    allocate(this%lwup_roof_lun(begl:endl))        ; this%lwup_roof_lun(:) = nan
    allocate(this%lwup_improad_lun(begl:endl))          ; this%lwup_improad_lun(:) = nan
    allocate(this%lwup_perroad_lun(begl:endl))          ; this%lwup_perroad_lun(:) = nan
    allocate(this%lwup_tree_perroad_lun(begl:endl))          ; this%lwup_tree_perroad_lun(:) = nan
    allocate(this%lwup_sunwall_lun(begl:endl))          ; this%lwup_sunwall_lun(:) = nan
    allocate(this%lwup_shadewall_lun(begl:endl))        ; this%lwup_shadewall_lun(:) = nan
    allocate(this%lwup_br_tree_lun(begl:endl))          ; this%lwup_br_tree_lun(:) = nan
    allocate(this%lwup_ar_tree_lun(begl:endl))          ; this%lwup_ar_tree_lun(:) = nan
    allocate(this%lwup_tree_lun(begl:endl))          ; this%lwup_tree_lun(:) = nan
    allocate(this%lwup_canyon_lun(begl:endl))           ; this%lwup_canyon_lun(:) = nan
    allocate(this%lwup_canyon_lun(begl:endl))           ; this%lwup_canyon_lun(:) = nan
    allocate(this%lwdown_road_lun(begl:endl))           ; this%lwdown_road_lun(:) = nan    
    allocate(this%lwdown_sunwall_lun(begl:endl))           ; this%lwdown_sunwall_lun(:) = nan    
    allocate(this%lwdown_shadewall_lun(begl:endl))           ; this%lwdown_shadewall_lun(:) = nan    
    allocate(this%lwdown_roof_lun(begl:endl))           ; this%lwdown_roof_lun(:) = nan    
    allocate(this%lwdown_br_tree_lun(begl:endl))           ; this%lwdown_br_tree_lun(:) = nan    
    allocate(this%lwdown_ar_tree_lun(begl:endl))           ; this%lwdown_ar_tree_lun(:) = nan    

   end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    use clm_varctl    , only : use_SSRE
    use clm_varpar    , only : nlevsno
    use histFileMod   , only : hist_addfld1d, hist_addfld2d
    use histFileMod   , only : no_snow_normal
    !
    ! !ARGUMENTS:
    class(solarabs_type) :: this
    type(bounds_type), intent(in) :: bounds  

    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: ptr_1d(:)      ! pointer to 1d patch array
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begl = bounds%begl; endl = bounds%endl
    
    ! confirm: if this output is correctly formatted. set_nourb=spval, l2g_scale_type='unity'?
    this%sdir_road_lun(begl:endl,:) = spval 
    call hist_addfld2d (fname='SDIR_ROAD', units='unitless', type2d='numrad', &
      avgflag='A', long_name='Unitless incident direct beam solar on road', &
      ptr_lunit=this%sdir_road_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')

    this%sdir_sunwall_lun(begl:endl,:) = spval 
    call hist_addfld2d (fname='SDIR_SUNWALL', units='unitless', type2d='numrad', &
      avgflag='A', long_name='Unitless incident direct beam solar on sun-lit wall', &
      ptr_lunit=this%sdir_sunwall_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')

    this%sdir_shadewall_lun(begl:endl,:) = spval 
    call hist_addfld2d (fname='SDIR_SHADEWALL', units='unitless', type2d='numrad', &
      avgflag='A', long_name='Unitless incident direct beam solar on shaded wall', &
      ptr_lunit=this%sdir_shadewall_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')   
      
    this%sdir_roof_lun(begl:endl,:) = spval 
    call hist_addfld2d (fname='SDIR_ROOF', units='unitless', type2d='numrad', &
      avgflag='A', long_name='Unitless incident direct beam solar on shaded roof', &
      ptr_lunit=this%sdir_roof_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')   
           
    this%sdir_ar_tree_lun(begl:endl,:) = spval 
    call hist_addfld2d (fname='SDIR_AR_TREE', units='unitless', type2d='numrad', &
      avgflag='A', long_name='Unitless incident direct beam solar on above roof tree', &
      ptr_lunit=this%sdir_ar_tree_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')   

    this%sdir_br_tree_lun(begl:endl,:) = spval 
    call hist_addfld2d (fname='SDIR_BR_TREE', units='unitless', type2d='numrad', &
      avgflag='A', long_name='Unitless incident direct beam solar on below roof tree', &
      ptr_lunit=this%sdir_br_tree_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')   

    this%sdir_tree_lun(begl:endl,:) = spval 
    call hist_addfld2d (fname='SDIR_TREE', units='unitless', type2d='numrad', &
      avgflag='A', long_name='Unitless incident direct beam solar on roof tree', &
      ptr_lunit=this%sdir_tree_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')     
               
    this%sdir_force_lun(begl:endl,:) = spval 
    call hist_addfld2d (fname='SDIR_FORCE', units='W m-2', type2d='numrad', &
      avgflag='A', long_name='Direct beam solar', &
      ptr_lunit=this%sdir_force_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')

    
    this%sabs_roof_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_ROOF_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident direct solar absorbed by roof per unit ground area per unit incident flux', &
       ptr_lunit=this%sabs_roof_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
       
    this%sabs_roof_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_ROOF_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident diffuse solar absorbed by roof per unit ground area per unit incident flux', &
       ptr_lunit=this%sabs_roof_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
       
    this%sabs_sunwall_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_SUNWALL_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident direct solar absorbed by sunwall per unit wall area per unit incident flux', &
       ptr_lunit=this%sabs_sunwall_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
    this%sabs_sunwall_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_SUNWALL_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident diffuse solar absorbed by sunwall per unit wall area per unit incident flux', &
       ptr_lunit=this%sabs_sunwall_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
       
    this%sabs_shadewall_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_SHADEWALL_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident direct solar absorbed by shadewall per unit wall area per unit incident flux', &
       ptr_lunit=this%sabs_shadewall_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
    this%sabs_shadewall_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_SHADEWALL_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident diffuse solar absorbed by shadewall per unit wall area per unit incident flux', &
       ptr_lunit=this%sabs_shadewall_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
       
    this%sabs_improad_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_IMPROAD_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident direct solar absorbed by impervious road per unit ground area per unit incident flux', &
       ptr_lunit=this%sabs_improad_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
    this%sabs_improad_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_IMPROAD_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident diffuse solar absorbed by impervious road per unit ground area per unit incident flux', &
       ptr_lunit=this%sabs_improad_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
       
    this%sabs_perroad_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_PERROAD_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident direct solar absorbed by pervious road per unit ground area per unit incident flux', &
       ptr_lunit=this%sabs_perroad_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
       
    this%sabs_perroad_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_PERROAD_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident diffuse solar absorbed by pervious road per unit ground area per unit incident flux', &
       ptr_lunit=this%sabs_perroad_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')


    this%sabs_br_tree_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_BR_TREE_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident direct solar absorbed by below-roof treeetation per unit treeetation area per unit incident flux', &
       ptr_lunit=this%sabs_br_tree_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
    this%sabs_br_tree_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_BR_TREE_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident diffuse solar absorbed by below-roof treeetation per unit treeetation area per unit incident flux', &
       ptr_lunit=this%sabs_br_tree_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
       
    this%sabs_ar_tree_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_AR_TREE_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident direct solar absorbed by above-roof treeetation per unit treeetation area per unit incident flux', &
       ptr_lunit=this%sabs_ar_tree_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
    this%sabs_ar_tree_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_AR_TREE_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident diffuse solar absorbed by above-roof treeetation per unit treeetation area per unit incident flux', &
       ptr_lunit=this%sabs_ar_tree_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')

   this%sabs_tree_dir_lun(begl:endl,:) = spval
   call hist_addfld2d (fname='SABS_TREE_DIR', units='unitless', type2d='numrad', &
      avgflag='A', long_name='Unitless incident direct solar absorbed by road treeetation per unit treeetation area per unit incident flux', &
      ptr_lunit=this%sabs_tree_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
   this%sabs_tree_dif_lun(begl:endl,:) = spval
   call hist_addfld2d (fname='SABS_TREE_DIF', units='unitless', type2d='numrad', &
      avgflag='A', long_name='Unitless incident diffuse solar absorbed by road treeetation per unit treeetation area per unit incident flux', &
      ptr_lunit=this%sabs_tree_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
      
    this%sabs_canyon_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_CANYON_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident direct solar absorbed by impervious road per unit ground area per unit incident flux', &
       ptr_lunit=this%sabs_canyon_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
    this%sabs_canyon_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SABS_CANYON_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident diffuse solar absorbed by impervious road per unit ground area per unit incident flux', &
       ptr_lunit=this%sabs_canyon_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')

    this%sref_roof_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_ROOF_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident roof direct solar flux lun', &
       ptr_lunit=this%sref_roof_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
    this%sref_roof_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_ROOF_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident roof diffuse solar flux lun', &
       ptr_lunit=this%sref_roof_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
       
    this%sref_sunwall_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_SUNWALL_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident sunwall direct solar flux lun', &
       ptr_lunit=this%sref_sunwall_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')       
    this%sref_sunwall_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_SUNWALL_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident sunwall diffuse solar flux lun', &
       ptr_lunit=this%sref_sunwall_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
       
    this%sref_shadewall_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_SHADEWALL_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident shadewall direct solar flux lun', &
       ptr_lunit=this%sref_shadewall_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
    this%sref_shadewall_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_SHADEWALL_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident shadewall diffuse solar flux lun', &
       ptr_lunit=this%sref_shadewall_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
       
    this%sref_improad_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_IMPROAD_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident improad direct solar flux lun', &
       ptr_lunit=this%sref_improad_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
    this%sref_improad_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_IMPROAD_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident improad diffuse solar flux lun', &
       ptr_lunit=this%sref_improad_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
       
    this%sref_perroad_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_PERROAD_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident perroad direct solar flux lun', &
       ptr_lunit=this%sref_perroad_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
    this%sref_perroad_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_PERROAD_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident perroad diffuse solar flux lun', &
       ptr_lunit=this%sref_perroad_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')

     this%sref_br_tree_dir_lun(begl:endl,:) = spval
     call hist_addfld2d (fname='SREF_BR_TREE_DIR', units='unitless', type2d='numrad', &
        avgflag='A', long_name='Unitless incident brree direct solar flux lun', &
        ptr_lunit=this%sref_br_tree_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
        default='inactive')
     this%sref_br_tree_dif_lun(begl:endl,:) = spval
     call hist_addfld2d (fname='SREF_BR_TREE_DIF', units='unitless', type2d='numrad', &
        avgflag='A', long_name='Unitless incident brree diffuse solar flux lun', &
        ptr_lunit=this%sref_br_tree_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
        default='inactive')
       
    this%sref_ar_tree_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_AR_TREE_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident tree direct solar flux lun', &
       ptr_lunit=this%sref_ar_tree_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
    this%sref_ar_tree_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_AR_TREE_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident tree diffuse solar flux lun', &
       ptr_lunit=this%sref_ar_tree_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')

    this%sref_tree_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_TREE_DIR', units='unitless', type2d='numrad', &
      avgflag='A', long_name='Unitless incident tree direct solar flux lun', &
      ptr_lunit=this%sref_tree_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
    this%sref_tree_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_TREE_DIF', units='unitless', type2d='numrad', &
      avgflag='A', long_name='Unitless incident tree diffuse solar flux lun', &
      ptr_lunit=this%sref_tree_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
      
    this%sref_canyon_dir_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_CANYON_DIR', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident canyon direct solar flux lun', &
       ptr_lunit=this%sref_canyon_dir_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')
    this%sref_canyon_dif_lun(begl:endl,:) = spval
    call hist_addfld2d (fname='SREF_CANYON_DIF', units='unitless', type2d='numrad', &
       avgflag='A', long_name='Unitless incident canyon diffuse solar flux lun', &
       ptr_lunit=this%sref_canyon_dif_lun, set_nourb=spval, l2g_scale_type='unity', &
       default='inactive')

    this%lwnet_roof_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWNET_ROOF', units='W/m^2', &
      avgflag='A', long_name='Incident net longwave flux at shaded roof (W/m^2)', &
      ptr_lunit=this%lwnet_roof_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
    this%lwnet_improad_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWNET_IMPROAD', units='W/m^2', &
      avgflag='A', long_name='Net longwave flux at impervious road (W/m^2)', &
      ptr_lunit=this%lwnet_improad_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')

    this%lw_emi_improad_lun(begl:endl) = spval
    call hist_addfld1d (fname='LW_EMI_IMPROAD', units='W/m^2', &
      avgflag='A', long_name='Emitted longwave flux at impervious road (W/m^2)', &
      ptr_lunit=this%lw_emi_improad_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')

    this%lwnet_perroad_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWNET_PERROAD', units='W/m^2', &
      avgflag='A', long_name='Incident net longwave flux at pervious road (W/m^2)', &
      ptr_lunit=this%lwnet_perroad_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')

    this%lwnet_tree_perroad_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWNET_TREE_PERROAD', units='W/m^2', &
      avgflag='A', long_name='Incident net longwave flux at tree pervious road (W/m^2)', &
      ptr_lunit=this%lwnet_tree_perroad_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')

    this%lwnet_sunwall_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWNET_SUNWALL', units='W/m^2', &
      avgflag='A', long_name='Incident net longwave flux at sunlit wall (W/m^2)', &
      ptr_lunit=this%lwnet_sunwall_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
    this%lwnet_shadewall_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWNET_SHADEWALL', units='W/m^2', &
      avgflag='A', long_name='Incident net longwave flux at shaded wall (W/m^2)', &
      ptr_lunit=this%lwnet_shadewall_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
      
    this%lwnet_br_tree_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWNET_BR_TREE', units='W/m^2', &
      avgflag='A', long_name='Incident net longwave flux at below-roof tree (W/m^2)', &
      ptr_lunit=this%lwnet_br_tree_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
    this%lwnet_ar_tree_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWNET_AR_TREE', units='W/m^2', &
      avgflag='A', long_name='Incident net longwave flux at above-roof tree (W/m^2)', &
      ptr_lunit=this%lwnet_ar_tree_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
      
    this%lwnet_tree_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWNET_TREE', units='W/m^2', &
      avgflag='A', long_name='Incident net longwave flux at road tree (W/m^2)', &
      ptr_lunit=this%lwnet_tree_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
      
    this%lwnet_canyon_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWNET_CANYON', units='W/m^2', &
      avgflag='A', long_name='Incident net longwave flux at canyon center (W/m^2)', &
      ptr_lunit=this%lwnet_canyon_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')

    this%lwdown_road_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWDOWN_ROAD', units='W/m^2', &
      avgflag='A', long_name='Incident downward longwave flux at road (W/m^2)', &
      ptr_lunit=this%lwdown_road_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')   

    this%lwdown_sunwall_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWDOWN_SUNWALL', units='W/m^2', &
      avgflag='A', long_name='Incident downward longwave flux at sunlit wall (W/m^2)', &
      ptr_lunit=this%lwdown_sunwall_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')   

      this%lwdown_shadewall_lun(begl:endl) = spval
      call hist_addfld1d (fname='LWDOWN_SHADEWALL', units='W/m^2', &
        avgflag='A', long_name='Incident downward longwave flux at shaded wall (W/m^2)', &
        ptr_lunit=this%lwdown_shadewall_lun, set_nourb=spval, l2g_scale_type='unity', &
        default='inactive')

      this%lwdown_roof_lun(begl:endl) = spval
      call hist_addfld1d (fname='LWDOWN_ROOF', units='W/m^2', &
        avgflag='A', long_name='Incident downward longwave flux at roof (W/m^2)', &
        ptr_lunit=this%lwdown_roof_lun, set_nourb=spval, l2g_scale_type='unity', &
        default='inactive')

      this%lwdown_br_tree_lun(begl:endl) = spval
      call hist_addfld1d (fname='LWDOWN_BR_TREE', units='W/m^2', &
        avgflag='A', long_name='Incident downward longwave flux at below-roof tree (W/m^2)', &
        ptr_lunit=this%lwdown_br_tree_lun, set_nourb=spval, l2g_scale_type='unity', &
        default='inactive')

      this%lwdown_ar_tree_lun(begl:endl) = spval
      call hist_addfld1d (fname='LWDOWN_AR_TREE', units='W/m^2', &
        avgflag='A', long_name='Incident downward longwave flux at above-roof tree (W/m^2)', &
        ptr_lunit=this%lwdown_ar_tree_lun, set_nourb=spval, l2g_scale_type='unity', &
        default='inactive')
        
    this%lwup_roof_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWUP_ROOF', units='W/m^2', &
      avgflag='A', long_name='Incident upward longwave flux at shaded roof (W/m^2)', &
      ptr_lunit=this%lwup_roof_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
    this%lwup_improad_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWUP_IMPROAD', units='W/m^2', &
      avgflag='A', long_name='Incident upward longwave flux at impervious road (W/m^2)', &
      ptr_lunit=this%lwup_improad_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
    this%lwup_perroad_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWUP_PERROAD', units='W/m^2', &
      avgflag='A', long_name='Incident upward longwave flux at pervious road (W/m^2)', &
      ptr_lunit=this%lwup_perroad_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')

    this%lwup_tree_perroad_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWUP_TREE_PERROAD', units='W/m^2', &
      avgflag='A', long_name='Incident upward longwave flux at tree pervious road (W/m^2)', &
      ptr_lunit=this%lwup_tree_perroad_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')

    this%lwup_sunwall_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWUP_SUNWALL', units='W/m^2', &
      avgflag='A', long_name='Incident upward longwave flux at sunlit wall (W/m^2)', &
      ptr_lunit=this%lwup_sunwall_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
      
    this%lwup_shadewall_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWUP_SHADEWALL', units='W/m^2', &
      avgflag='A', long_name='Incident upward longwave flux at shaded wall (W/m^2)', &
      ptr_lunit=this%lwup_shadewall_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
      
    this%lwup_br_tree_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWUP_BR_TREE', units='W/m^2', &
      avgflag='A', long_name='Incident upward longwave flux at below-roof tree (W/m^2)', &
      ptr_lunit=this%lwup_br_tree_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
      
    this%lwup_ar_tree_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWUP_AR_TREE', units='W/m^2', &
      avgflag='A', long_name='Incident upward longwave flux at above-roof tree (W/m^2)', &
      ptr_lunit=this%lwup_ar_tree_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
      
    this%lwup_tree_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWUP_TREE', units='W/m^2', &
      avgflag='A', long_name='Incident upward longwave flux at road tree (W/m^2)', &
      ptr_lunit=this%lwup_tree_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
      
    this%lwup_canyon_lun(begl:endl) = spval
    call hist_addfld1d (fname='LWUP_CANYON', units='W/m^2', &
      avgflag='A', long_name='Incident upward longwave flux at canyon center (W/m^2)', &
      ptr_lunit=this%lwup_canyon_lun, set_nourb=spval, l2g_scale_type='unity', &
      default='inactive')
     
    this%fsa_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSA', units='W/m^2',  &
         avgflag='A', long_name='absorbed solar radiation', &
         ptr_patch=this%fsa_patch, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSA_ICE', units='W/m^2',  &
         avgflag='A', long_name='absorbed solar radiation (ice landunits only)', &
         ptr_patch=this%fsa_patch, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%fsa_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSA_R', units='W/m^2',  &
         avgflag='A', long_name='Rural absorbed solar radiation', &
         ptr_patch=this%fsa_r_patch, set_spec=spval, default='inactive')

    this%fsa_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSA_U', units='W/m^2',  &
         avgflag='A', long_name='Urban absorbed solar radiation', &
         ptr_patch=this%fsa_u_patch, c2l_scale_type='urbanf', set_nourb=spval, default='inactive')

    this%fsr_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSR', units='W/m^2',  &
         avgflag='A', long_name='reflected solar radiation', &
         ptr_patch=this%fsr_patch, c2l_scale_type='urbanf')
    ! Rename of FSR for Urban intercomparision project
    call hist_addfld1d (fname='SWup', units='W/m^2',  &
         avgflag='A', long_name='upwelling shortwave radiation', &
         ptr_patch=this%fsr_patch, c2l_scale_type='urbanf', default='inactive')
    call hist_addfld1d (fname='FSR_ICE', units='W/m^2',  &
         avgflag='A', long_name='reflected solar radiation (ice landunits only)', &
         ptr_patch=this%fsr_patch, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%sabg_lyr_patch(begp:endp,-nlevsno+1:0) = spval
    data2dptr => this%sabg_lyr_patch(:,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_ABS', units='W/m^2', type2d='levsno',  &
         avgflag='A', long_name='Absorbed solar radiation in each snow layer', &
         ptr_patch=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    call hist_addfld2d (fname='SNO_ABS_ICE', units='W/m^2', type2d='levsno',  &
         avgflag='A', long_name='Absorbed solar radiation in each snow layer (ice landunits only)', &
         ptr_patch=data2dptr, no_snow_behavior=no_snow_normal, &
         l2g_scale_type='ice', default='inactive')

    this%sabv_patch(begp:endp) = spval
    call hist_addfld1d (fname='SABV', units='W/m^2',  &
         avgflag='A', long_name='solar rad absorbed by veg', &
         ptr_patch=this%sabv_patch, c2l_scale_type='urbanf')

    this%sabg_patch(begp:endp) = spval
    call hist_addfld1d (fname='SABG', units='W/m^2',  &
         avgflag='A', long_name='solar rad absorbed by ground', &
         ptr_patch=this%sabg_patch, c2l_scale_type='urbanf')

    this%sabg_pen_patch(begp:endp) = spval
    call hist_addfld1d (fname='SABG_PEN', units='watt/m^2',  &
         avgflag='A', long_name='Rural solar rad penetrating top soil or snow layer', &
         ptr_patch=this%sabg_pen_patch, set_spec=spval)

     ! Currently needed by lake code - TODO should not be here
    this%fsds_nir_d_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSND', units='W/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation', &
         ptr_patch=this%fsds_nir_d_patch)

    this%fsds_nir_i_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir incident solar radiation', &
         ptr_patch=this%fsds_nir_i_patch)

    this%fsds_nir_d_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSNDLN', units='W/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation at local noon', &
         ptr_patch=this%fsds_nir_d_ln_patch)

    this%fsr_nir_d_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRND', units='W/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation', &
         ptr_patch=this%fsr_nir_d_patch, c2l_scale_type='urbanf')

    this%fsr_nir_i_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir reflected solar radiation', &
         ptr_patch=this%fsr_nir_i_patch, c2l_scale_type='urbanf')

    this%fsr_nir_d_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRNDLN', units='W/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation at local noon', &
         ptr_patch=this%fsr_nir_d_ln_patch, c2l_scale_type='urbanf')
    ! diagnostic fluxes for ESM-SnowMIP
    if (use_SSRE) then
       this%fsrSF_patch(begp:endp) = spval
       call hist_addfld1d (fname='FSRSF', units='W/m^2',  &
            avgflag='A', long_name='reflected solar radiation', &
            ptr_patch=this%fsrSF_patch, c2l_scale_type='urbanf')

       this%ssre_fsr_patch(begp:endp) = spval
       call hist_addfld1d (fname='SSRE_FSR', units='W/m^2',  &
            avgflag='A', long_name='surface snow effect on reflected solar radiation', &
            ptr_patch=this%ssre_fsr_patch, c2l_scale_type='urbanf')
       this%fsrSF_nir_d_patch(begp:endp) = spval
       call hist_addfld1d (fname='FSRSFND', units='W/m^2',  &
            avgflag='A', long_name='direct nir reflected solar radiation', &
            ptr_patch=this%fsrSF_nir_d_patch, c2l_scale_type='urbanf')

       this%fsrSF_nir_i_patch(begp:endp) = spval
       call hist_addfld1d (fname='FSRSFNI', units='W/m^2',  &
            avgflag='A', long_name='diffuse nir reflected solar radiation', &
            ptr_patch=this%fsrSF_nir_i_patch, c2l_scale_type='urbanf')

       this%fsrSF_nir_d_ln_patch(begp:endp) = spval
       call hist_addfld1d (fname='FSRSFNDLN', units='W/m^2',  &
            avgflag='A', long_name='direct nir reflected solar radiation at local noon', &
            ptr_patch=this%fsrSF_nir_d_ln_patch, c2l_scale_type='urbanf')

       this%ssre_fsr_nir_d_patch(begp:endp) = spval
       call hist_addfld1d (fname='SSRE_FSRND', units='W/m^2',  &
            avgflag='A', long_name='surface snow effect on direct nir reflected solar radiation', &
            ptr_patch=this%ssre_fsr_nir_d_patch, c2l_scale_type='urbanf')

       this%ssre_fsr_nir_i_patch(begp:endp) = spval
       call hist_addfld1d (fname='SSRE_FSRNI', units='W/m^2',  &
            avgflag='A', long_name='surface snow effect on diffuse nir reflected solar radiation', &
            ptr_patch=this%ssre_fsr_nir_i_patch, c2l_scale_type='urbanf')

       this%ssre_fsr_nir_d_ln_patch(begp:endp) = spval
       call hist_addfld1d (fname='SSRE_FSRNDLN', units='W/m^2',  &
            avgflag='A', long_name='surface snow effect on direct nir reflected solar radiation at local noon', &
            ptr_patch=this%ssre_fsr_nir_d_ln_patch, c2l_scale_type='urbanf')
    end if
    this%sub_surf_abs_SW_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOINTABS', units='-', &
         avgflag='A', long_name='Fraction of incoming solar absorbed by lower snow layers', &
         ptr_patch=this%sub_surf_abs_SW_patch, set_lake=spval, set_urb=spval)

    if(use_luna)then
       ptr_1d => this%par240d_z_patch(:,1)
       call hist_addfld1d (fname='PAR240DZ', units='W/m^2', &
         avgflag='A', long_name='10-day running mean of daytime patch absorbed PAR for leaves for top canopy layer', &
         ptr_patch=ptr_1d, default='inactive')
       ptr_1d => this%par240x_z_patch(:,1)
       call hist_addfld1d (fname='PAR240XZ', units='W/m^2', &
         avgflag='A', long_name='10-day running mean of maximum patch absorbed PAR for leaves for top canopy layer', &
         ptr_patch=ptr_1d, default='inactive')

    endif

  end subroutine InitHistory

  !------------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! Initialize module surface albedos to reasonable values
    !
    use landunit_varcon, only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    class(solarabs_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begl, endl
    !-----------------------------------------------------------------------

    begl = bounds%begl; endl = bounds%endl

    this%sabs_roof_dir_lun      (begl:endl, :) = 0._r8    
    this%sabs_roof_dif_lun      (begl:endl, :) = 0._r8    
    this%sabs_sunwall_dir_lun   (begl:endl, :) = 0._r8
    this%sabs_sunwall_dif_lun   (begl:endl, :) = 0._r8
    this%sabs_shadewall_dir_lun (begl:endl, :) = 0._r8
    this%sabs_shadewall_dif_lun (begl:endl, :) = 0._r8
    this%sabs_improad_dir_lun   (begl:endl, :) = 0._r8
    this%sabs_improad_dif_lun   (begl:endl, :) = 0._r8
    this%sabs_perroad_dir_lun   (begl:endl, :) = 0._r8
    this%sabs_perroad_dif_lun   (begl:endl, :) = 0._r8
    this%sabs_tree_dir_lun   (begl:endl, :) = 0._r8
    this%sabs_tree_dif_lun   (begl:endl, :) = 0._r8

  end subroutine InitCold

  !---------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use shr_infnan_mod , only : shr_infnan_isnan
    use spmdMod        , only : masterproc
    use abortutils     , only : endrun
    use ncdio_pio      , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(solarabs_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    integer :: p
    !---------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='sabs_roof_dir', xtype=ncd_double,  dim1name='landunit',            & 
         dim2name='numrad', switchdim=.true.,                                                                         &
         long_name='direct solar absorbed by roof per unit ground area per unit incident flux', units='',             &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%sabs_roof_dir_lun)

    call restartvar(ncid=ncid, flag=flag, varname='sabs_roof_dif', xtype=ncd_double,  dim1name='landunit',            & 
         dim2name='numrad', switchdim=.true.,                                                                         &
         long_name='diffuse solar absorbed by roof per unit ground area per unit incident flux', units='',            &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%sabs_roof_dif_lun)

    call restartvar(ncid=ncid, flag=flag, varname='sabs_sunwall_dir', xtype=ncd_double,  dim1name='landunit',         & 
         dim2name='numrad', switchdim=.true.,                                                                         &
         long_name='direct solar absorbed by sunwall per unit wall area per unit incident flux', units='',            &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%sabs_sunwall_dir_lun)

    call restartvar(ncid=ncid, flag=flag, varname='sabs_sunwall_dif', xtype=ncd_double,  dim1name='landunit',         & 
         dim2name='numrad', switchdim=.true.,                                                                         &
         long_name='diffuse solar absorbed by sunwall per unit wall area per unit incident flux', units='',           &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%sabs_sunwall_dif_lun)

    call restartvar(ncid=ncid, flag=flag, varname='sabs_shadewall_dir', xtype=ncd_double,  dim1name='landunit',       & 
         dim2name='numrad', switchdim=.true.,                                                                         &
         long_name='direct solar absorbed by shadewall per unit wall area per unit incident flux', units='',          &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%sabs_shadewall_dir_lun)

    call restartvar(ncid=ncid, flag=flag, varname='sabs_shadewall_dif', xtype=ncd_double,  dim1name='landunit',       & 
         dim2name='numrad', switchdim=.true.,                                                                         &
         long_name='diffuse solar absorbed by shadewall per unit wall area per unit incident flux', units='',         &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%sabs_shadewall_dif_lun)

    call restartvar(ncid=ncid, flag=flag, varname='sabs_improad_dir', xtype=ncd_double,  dim1name='landunit',         & 
         dim2name='numrad', switchdim=.true.,                                                                         &
         long_name='direct solar absorbed by impervious road per unit ground area per unit incident flux', units='',  &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%sabs_improad_dir_lun)

    call restartvar(ncid=ncid, flag=flag, varname='sabs_improad_dif', xtype=ncd_double,  dim1name='landunit',         & 
         dim2name='numrad', switchdim=.true.,                                                                         &
         long_name='diffuse solar absorbed by impervious road per unit ground area per unit incident flux', units='', &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%sabs_improad_dif_lun)

    call restartvar(ncid=ncid, flag=flag, varname='sabs_perroad_dir', xtype=ncd_double,  dim1name='landunit',         & 
         dim2name='numrad', switchdim=.true.,                                                                         &
         long_name='direct solar absorbed by pervious road per unit ground area per unit incident flux', units='',    &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%sabs_perroad_dir_lun)

    call restartvar(ncid=ncid, flag=flag, varname='sabs_perroad_dif', xtype=ncd_double,  dim1name='landunit',         & 
         dim2name='numrad', switchdim=.true.,                                                                         &
         long_name='diffuse solar absorbed by pervious road per unit ground area per unit incident flux', units='',   &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%sabs_perroad_dif_lun)

    call restartvar(ncid=ncid, flag=flag, varname='sabs_tree_dir', xtype=ncd_double,  dim1name='landunit',         & 
         dim2name='numrad', switchdim=.true.,                                                                         &
         long_name='direct solar absorbed by road tree per unit ground area per unit incident flux', units='',    &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%sabs_tree_dir_lun)

    call restartvar(ncid=ncid, flag=flag, varname='sabs_tree_dif', xtype=ncd_double,  dim1name='landunit',         & 
         dim2name='numrad', switchdim=.true.,                                                                         &
         long_name='diffuse solar absorbed by road tree per unit ground area per unit incident flux', units='',   &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%sabs_tree_dif_lun)

   if(use_luna)then
      call restartvar(ncid=ncid, flag=flag, varname='par240d', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='10-day running mean of daytime absorbed PAR for leaves in canopy layer', units='W/m**2 leaf', &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%par240d_z_patch )
      call restartvar(ncid=ncid, flag=flag, varname='par24d', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='Accumulative daytime absorbed PAR for leaves in canopy layer for 24 hours', units='J/m**2 leaf', &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%par24d_z_patch )

      call restartvar(ncid=ncid, flag=flag, varname='par240x', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='10-day running mean of maximum absorbed PAR for leaves in canopy layers', units='W/m**2 leaf', &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%par240x_z_patch )
      call restartvar(ncid=ncid, flag=flag, varname='par24x', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='Maximum absorbed PAR for leaves in canopy layer in 24 hours', units='J/m**2 leaf', &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%par24x_z_patch )

      call restartvar(ncid=ncid, flag=flag, varname='parsun', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='Instaneous absorbed PAR for sunlit leaves in canopy layer', units='W/m**2 leaf', &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%parsun_z_patch )
      call restartvar(ncid=ncid, flag=flag, varname='parsha', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='Instaneous absorbed PAR for shaded leaves in canopy layer', units='W/m**2 leaf', &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=this%parsha_z_patch )

   endif

  end subroutine Restart

end module SolarAbsorbedType
