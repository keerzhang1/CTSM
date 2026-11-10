module UrbanParamsType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Urban Constants
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use abortutils   , only : endrun
  use decompMod    , only : bounds_type, subgrid_level_gridcell, subgrid_level_landunit
  use clm_varctl   , only : iulog, fsurdat
  use clm_varcon   , only : grlnd, spval, ispval
  use LandunitType , only : lun

  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: UrbanReadNML      ! Read in the urban namelist items
  public  :: UrbanInput        ! Read in urban input data
  public  :: CheckUrban        ! Check validity of urban points
  public  :: IsSimpleBuildTemp ! If using the simple building temperature method
  public  :: IsProgBuildTemp   ! If using the prognostic building temperature method
  public  :: UrbanPhenology   ! 
  !
  ! PRIVATE MEMBER FUNCTIONS
  private :: montecarlo_view_factors    ! Calculates the view factors between all 'surfaces'  
  private :: ComputeViewFactors        ! Compute Monte Carlo view factors for one urban landunit
  private :: ComputeUrbanRoughness     ! Compute displacement height and roughness length for one urban landunit
  private :: corner_up         ! Determines if a ray strikes the current wall layer or not at the upstream building edge
  private :: corner_dn         ! Determines if a ray strikes a roof or not at the downstream edge of a building
  private :: init_random_seed  ! Initializes the random seed for ray tracing calculations
  private :: ray_up            ! Computes the contribution to view factors of rays traveling upwards
  private :: ray_dn            ! Computes the contribution to view factors of rays traveling downwards
  ! !PRIVATE TYPE
  type urbinp_type
     real(r8), pointer :: canyon_hwr      (:,:)  
     real(r8), pointer :: tree_lai_urb    (:,:,:)
     real(r8), pointer :: tree_sai_urb    (:,:,:)
     integer, pointer :: tree_pft_urb    (:,:) 
     real(r8), pointer :: wtroad_tree    (:,:)      
     real(r8), pointer :: tree_bht_urb     (:,:)  
     real(r8), pointer :: tree_tht_urb     (:,:)      
     real(r8), pointer :: wtlunit_roof    (:,:)  
     real(r8), pointer :: wtroad_perv     (:,:)  
     real(r8), pointer :: wall_to_plan_area_ratio (:,:)
     real(r8), pointer :: em_roof         (:,:)   
     real(r8), pointer :: em_improad      (:,:)  
     real(r8), pointer :: em_perroad      (:,:)  
     real(r8), pointer :: em_tree_urb      (:,:)  
     real(r8), pointer :: em_wall         (:,:)  
     real(r8), pointer :: alb_roof_dir    (:,:,:)  
     real(r8), pointer :: alb_roof_dif    (:,:,:)  
     real(r8), pointer :: alb_improad_dir (:,:,:)  
     real(r8), pointer :: alb_improad_dif (:,:,:)  
     real(r8), pointer :: alb_perroad_dir (:,:,:)  
     real(r8), pointer :: alb_perroad_dif (:,:,:)  
     real(r8), pointer :: alb_tree_urb_dir (:,:,:)  
     real(r8), pointer :: alb_tree_urb_dif (:,:,:)  
     real(r8), pointer :: tran_tree_urb_dir (:,:,:)  
     real(r8), pointer :: tran_tree_urb_dif (:,:,:)  

     real(r8), pointer :: alb_wall_dir    (:,:,:)  
     real(r8), pointer :: alb_wall_dif    (:,:,:)  
     real(r8), pointer :: ht_roof         (:,:)
     real(r8), pointer :: ht_can_eff         (:,:)
     real(r8), pointer :: wind_hgt_canyon (:,:)
     real(r8), pointer :: tk_wall         (:,:,:)
     real(r8), pointer :: tk_roof         (:,:,:)
     real(r8), pointer :: tk_improad      (:,:,:)
     real(r8), pointer :: cv_wall         (:,:,:)
     real(r8), pointer :: cv_roof         (:,:,:)
     real(r8), pointer :: cv_improad      (:,:,:)
     real(r8), pointer :: thick_wall      (:,:)
     real(r8), pointer :: thick_roof      (:,:)
     integer,  pointer :: nlev_improad    (:,:)
     real(r8), pointer :: t_building_min  (:,:)
     real(r8), pointer :: A_v1  (:,:)
     real(r8), pointer :: A_v2  (:,:)

     
  end type urbinp_type
  type (urbinp_type), public :: urbinp   ! urban input derived type

  ! PUBLIC TYPE
  type, public :: urbanparams_type
     ! these output variables are kept for diagnosis purpose for now
     real(r8), pointer :: fww1d_out       (:,:,:) ! Longwave radiation view factor from wall to wall 
     real(r8), pointer :: fvv1d_out       (:,:,:) ! Longwave radiation view factor from vegetation to vegetation 
     real(r8), pointer :: fwv1d_out       (:,:,:) ! Longwave radiation view factor from wall to vegetation 
     real(r8), pointer :: fvw1d_out       (:,:,:) ! Longwave radiation view factor from vegetation to wall 
     real(r8), pointer :: fwr1d_out       (:,:,:) ! Longwave radiation view factor from wall to roof 
     real(r8), pointer :: frw1d_out       (:,:,:) ! Longwave radiation view factor from roof to wall 
     real(r8), pointer :: fvr1d_out       (:,:,:) ! Longwave radiation view factor from vegetation to roof 
     real(r8), pointer :: frv1d_out       (:,:,:) ! Longwave radiation view factor from roof to vegetation 
     
     real(r8), pointer :: fwg1d_out       (:,:)   ! Longwave radiation view factor from wall to ground 
     real(r8), pointer :: fgw1d_out       (:,:)   ! Longwave radiation view factor from ground to wall 
     real(r8), pointer :: fgv1d_out       (:,:)   ! Longwave radiation view factor from ground to vegetation 
     real(r8), pointer :: fsw1d_out       (:,:)   ! Longwave radiation view factor from sky to wall 
     real(r8), pointer :: fvg1d_out       (:,:)   ! Longwave radiation view factor from vegetation to ground 
     real(r8), pointer :: fsr1d_out       (:,:)   ! Longwave radiation view factor from sky to roof 
     real(r8), pointer :: fsv1d_out       (:,:)   ! Longwave radiation view factor from sky to vegetation 
     
     real(r8), pointer :: fws1d_out       (:,:)   ! Longwave radiation view factor from wall to sky 
     real(r8), pointer :: fvs1d_out       (:,:)   ! Longwave radiation view factor from vegetation to sky 
     real(r8), pointer :: frs1d_out       (:,:)   ! Longwave radiation view factor from roof to sky 
     real(r8), pointer :: fsg1d_out       (:)     ! Longwave radiation view factor from sky to ground 
     real(r8), pointer :: fts1d_out       (:)     ! Longwave radiation view factor from ground to sky 

     real(r8), pointer :: kww1d_out       (:,:,:) ! Shortwave radiation view factor from wall to wall 
     real(r8), pointer :: kvv1d_out       (:,:,:) ! Shortwave radiation view factor from vegetation to vegetation 
     real(r8), pointer :: kwv1d_out       (:,:,:) ! Shortwave radiation view factor from wall to vegetation 
     real(r8), pointer :: kvw1d_out       (:,:,:) ! Shortwave radiation view factor from vegetation to wall 
     real(r8), pointer :: kwr1d_out       (:,:,:) ! Shortwave radiation view factor from wall to roof 
     real(r8), pointer :: krw1d_out       (:,:,:) ! Shortwave radiation view factor from roof to wall 
     real(r8), pointer :: kvr1d_out       (:,:,:) ! Shortwave radiation view factor from vegetation to roof 
     real(r8), pointer :: krv1d_out       (:,:,:) ! Shortwave radiation view factor from roof to vegetation 
     
     real(r8), pointer :: kwg1d_out       (:,:)   ! Shortwave radiation view factor from wall to ground 
     real(r8), pointer :: kgw1d_out       (:,:)   ! Shortwave radiation view factor from ground to wall 
     real(r8), pointer :: kgv1d_out       (:,:)   ! Shortwave radiation view factor from ground to vegetation 
     real(r8), pointer :: ksw1d_out       (:,:)   ! Shortwave radiation view factor from sky to wall 
     real(r8), pointer :: kvg1d_out       (:,:)   ! Shortwave radiation view factor from vegetation to ground 
     real(r8), pointer :: ksr1d_out       (:,:)   ! Shortwave radiation view factor from sky to roof 
     real(r8), pointer :: ksv1d_out       (:,:)   ! Shortwave radiation view factor from sky to vegetation 
     
     real(r8), pointer :: kws1d_out       (:,:)   ! Shortwave radiation view factor from wall to sky 
     real(r8), pointer :: kvs1d_out       (:,:)   ! Shortwave radiation view factor from vegetation to sky 
     real(r8), pointer :: krs1d_out       (:,:)   ! Shortwave radiation view factor from roof to sky
     real(r8), pointer :: ksg1d_out       (:)     ! Shortwave radiation view factor from sky to ground 
     real(r8), pointer :: kts1d_out       (:)     ! Shortwave radiation view factor from ground to sky
     
     ! Unweighted longwave radiation view factor outputs
     real(r8), pointer :: vfww_f_out      (:,:,:)  ! Unweighted longwave view factor from wall to wall  
     real(r8), pointer :: vfvv_f_out      (:,:,:)  ! Unweighted longwave view factor from vegetation to vegetation  
     real(r8), pointer :: vfwv_f_out      (:,:,:)  ! Unweighted longwave view factor from wall to vegetation  
     real(r8), pointer :: vfvw_f_out      (:,:,:)  ! Unweighted longwave view factor from vegetation to wall  
     real(r8), pointer :: vfwr_f_out      (:,:,:)  ! Unweighted longwave view factor from wall to roof  
     real(r8), pointer :: vfrw_f_out      (:,:,:)  ! Unweighted longwave view factor from roof to wall  
     real(r8), pointer :: vfvr_f_out      (:,:,:)  ! Unweighted longwave view factor from vegetation to roof  
     real(r8), pointer :: vfrv_f_out      (:,:,:)  ! Unweighted longwave view factor from roof to vegetation  
     real(r8), pointer :: vfwt_f_out      (:,:)    ! Unweighted longwave view factor from wall to ground  
     real(r8), pointer :: vftw_f_out      (:,:)    ! Unweighted longwave view factor from ground to wall  
     real(r8), pointer :: vftv_f_out      (:,:)    ! Unweighted longwave view factor from ground to vegetation  
     real(r8), pointer :: vfvt_f_out      (:,:)    ! Unweighted longwave view factor from vegetation to ground  
     real(r8), pointer :: vfsw_f_out      (:,:)    ! Unweighted longwave view factor from sky to wall  
     real(r8), pointer :: vfsr_f_out      (:,:)    ! Unweighted longwave view factor from sky to roof  
     real(r8), pointer :: vfsv_f_out      (:,:)    ! Unweighted longwave view factor from sky to vegetation  
     ! Unweighted sky view factor outputs (for longwave)
     real(r8), pointer :: svfw_f_out      (:,:)    ! Unweighted longwave sky view factor for wall  
     real(r8), pointer :: svfv_f_out      (:,:)    ! Unweighted longwave sky view factor for vegetation  
     real(r8), pointer :: svfr_f_out      (:,:)    ! Unweighted longwave sky view factor for roof  
     real(r8), pointer :: vfst_f_out      (:)      ! Unweighted longwave view factor from sky to ground  
     real(r8), pointer :: svft_f_out      (:)      ! Unweighted longwave sky view factor of the ground (ground -> sky)

     ! Unweighted shortwave radiation view factor outputs
     real(r8), pointer :: vfww_k_out      (:,:,:)  ! Unweighted shortwave view factor from wall to wall  
     real(r8), pointer :: vfvv_k_out      (:,:,:)  ! Unweighted shortwave view factor from vegetation to vegetation  
     real(r8), pointer :: vfwv_k_out      (:,:,:)  ! Unweighted shortwave view factor from wall to vegetation  
     real(r8), pointer :: vfvw_k_out      (:,:,:)  ! Unweighted shortwave view factor from vegetation to wall  
     real(r8), pointer :: vfwr_k_out      (:,:,:)  ! Unweighted shortwave view factor from wall to roof  
     real(r8), pointer :: vfrw_k_out      (:,:,:)  ! Unweighted shortwave view factor from roof to wall  
     real(r8), pointer :: vfvr_k_out      (:,:,:)  ! Unweighted shortwave view factor from vegetation to roof  
     real(r8), pointer :: vfrv_k_out      (:,:,:)  ! Unweighted shortwave view factor from roof to vegetation  
     real(r8), pointer :: vfwt_k_out      (:,:)    ! Unweighted shortwave view factor from wall to ground  
     real(r8), pointer :: vftw_k_out      (:,:)    ! Unweighted shortwave view factor from ground to wall  
     real(r8), pointer :: vftv_k_out      (:,:)    ! Unweighted shortwave view factor from ground to vegetation  
     real(r8), pointer :: vfvt_k_out      (:,:)    ! Unweighted shortwave view factor from vegetation to ground  
     real(r8), pointer :: vfsw_k_out      (:,:)    ! Unweighted shortwave view factor from sky to wall  
     real(r8), pointer :: vfsr_k_out      (:,:)    ! Unweighted shortwave view factor from sky to roof  
     real(r8), pointer :: vfsv_k_out      (:,:)    ! Unweighted shortwave view factor from sky to vegetation  
     ! Unweighted sky view factor outputs (for shortwave)
     real(r8), pointer :: svfw_k_out      (:,:)    ! Unweighted shortwave sky view factor for wall  
     real(r8), pointer :: svfv_k_out      (:,:)    ! Unweighted shortwave sky view factor for vegetation  
     real(r8), pointer :: svfr_k_out      (:,:)    ! Unweighted shortwave sky view factor for roof  
     real(r8), pointer :: vfst_k_out      (:)      ! Unweighted shortwave view factor from sky to ground  
     real(r8), pointer :: svft_k_out      (:)      ! Unweighted shortwave sky view factor of the ground (ground -> sky)
     
     real(r8), pointer :: frontal_ai_out      (:)      ! frontal area index for output 
     real(r8), pointer :: plan_ai_out      (:)      ! plan area index for output
     real(r8), pointer :: z_d_town_out      (:)      ! displacement height for output
     real(r8), pointer :: z_0_town_out      (:)      ! roughness length for output
     real(r8), allocatable :: wind_hgt_canyon     (:)   ! lun height above road at which wind in canyon is to be computed (m)
     real(r8), allocatable :: em_roof             (:)   ! lun roof emissivity
     real(r8), allocatable :: em_improad          (:)   ! lun impervious road emissivity
     real(r8), allocatable :: em_perroad          (:)   ! lun pervious road emissivity
     real(r8), allocatable :: em_tree_urb          (:)   ! lun road tree emissivity
     real(r8), allocatable :: em_wall             (:)   ! lun wall emissivity
     real(r8), allocatable :: alb_roof_dir        (:,:) ! lun direct  roof albedo
     real(r8), allocatable :: alb_roof_dif        (:,:) ! lun diffuse roof albedo
     real(r8), allocatable :: alb_improad_dir     (:,:) ! lun direct  impervious road albedo
     real(r8), allocatable :: alb_improad_dif     (:,:) ! lun diffuse impervious road albedo
     real(r8), allocatable :: alb_perroad_dir     (:,:) ! lun direct  pervious road albedo
     real(r8), allocatable :: alb_perroad_dif     (:,:) ! lun diffuse pervious road albedo
     real(r8), allocatable :: alb_tree_urb_dir     (:,:) ! lun direct road tree albedo
     real(r8), allocatable :: alb_tree_urb_dif     (:,:) ! lun diffuse road tree albedo
     real(r8), allocatable :: tran_tree_urb_dif     (:,:) ! lun diffuse road tree transmittance
     real(r8), allocatable :: tran_tree_urb_dir     (:,:) ! lun direct road tree transmittance
     real(r8), allocatable :: alb_wall_dir        (:,:) ! lun direct  wall albedo
     real(r8), allocatable :: alb_wall_dif        (:,:) ! lun diffuse wall albedo

     real(r8), allocatable :: tree_elaisai_per_uroad        (:) ! lun road tree exposed leaf + stem area index (elai+esai, snow-adjusted) per unit total road area (m^2/m^2)
     real(r8), allocatable :: tree_lai_per_vega        (:) ! lun road tree LAI per unit tree-covered area, interpolated from monthly surfdata TREE_LAI_URB (m^2/m^2)
     real(r8), allocatable :: tree_sai_per_vega        (:) ! lun road tree SAI per unit tree-covered area, interpolated from monthly surfdata TREE_SAI_URB (m^2/m^2)
     integer, allocatable :: tree_pft_urb        (:) ! lun PFT index of the road tree patch, from surfdata TREE_PFT_URB (-)

     integer , pointer     :: nlev_improad        (:)   ! lun number of impervious road layers (-)
     real(r8), pointer     :: tk_wall             (:,:) ! lun thermal conductivity of urban wall (W/m/K)
     real(r8), pointer     :: tk_roof             (:,:) ! lun thermal conductivity of urban roof (W/m/K)
     real(r8), pointer     :: tk_improad          (:,:) ! lun thermal conductivity of urban impervious road (W/m/K)
     real(r8), pointer     :: cv_wall             (:,:) ! lun heat capacity of urban wall (J/m^3/K)
     real(r8), pointer     :: cv_roof             (:,:) ! lun heat capacity of urban roof (J/m^3/K)
     real(r8), pointer     :: cv_improad          (:,:) ! lun heat capacity of urban impervious road (J/m^3/K)
     real(r8), pointer     :: thick_wall          (:)   ! lun total thickness of urban wall (m)
     real(r8), pointer     :: thick_roof          (:)   ! lun total thickness of urban roof (m)

     real(r8), pointer     :: t_building_min      (:)   ! lun minimum internal building air temperature (K)
     real(r8), pointer     :: eflx_traffic_factor (:)   ! lun multiplicative traffic factor for sensible heat flux from urban traffic (-)
   contains

     procedure, public :: Init

  end type urbanparams_type
  !
  ! !Urban control variables
  character(len= *), parameter, public :: urban_hac_off = 'OFF'                
  character(len= *), parameter, public :: urban_hac_on =  'ON'                 
  character(len= *), parameter, public :: urban_wasteheat_on = 'ON_WASTEHEAT'  
  character(len= 16), public           :: urban_hac = urban_hac_off
  logical, public                      :: urban_explicit_ac = .true.  ! whether to use explicit, time-varying AC adoption rate
  logical, public                      :: urban_traffic = .false.     ! urban traffic fluxes

  ! !PRIVATE MEMBER DATA:
  logical, private    :: ReadNamelist = .false.     ! If namelist was read yet or not
  integer, parameter, private :: BUILDING_TEMP_METHOD_SIMPLE = 0       ! Simple method introduced in CLM4.5
  integer, parameter, private :: BUILDING_TEMP_METHOD_PROG   = 1       ! Prognostic method introduced in CLM5.0
  integer, private :: building_temp_method = BUILDING_TEMP_METHOD_PROG ! Method to calculate the building temperature

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !----------------------------------------------------------------------- 

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar      , only : numrad, nlevurb
    use clm_varctl      , only : use_vancouver, use_mexicocity
    use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon   , only : icol_road_perv, icol_road_imperv, icol_road_perv, icol_road_tree
    use landunit_varcon , only : isturb_MIN
    use clm_varcon        , only : rpi
    !
    ! !ARGUMENTS:
    class(urbanparams_type) :: this
    type(bounds_type)      , intent(in)    :: bounds  
    !
    ! !LOCAL VARIABLES:   
    integer             :: j,l,c,p,g       ! indices
    integer             :: nc,fl,ib, im        ! indices 
    integer             :: dindx           ! urban density type index
    integer             :: ier             ! error status
    real(r8)            :: plan_ai         ! plan area index - ratio building area to plan area (-)
    real(r8)            :: frontal_ai      ! frontal area index of buildings (-)
    real(r8)            :: wall_to_plan_area_ratio  ! provided by PLUMBER site input data
    integer		:: begl, endl
    integer		:: begc, endc
    integer		:: begp, endp
    integer             :: begg, endg
    ! Variables for monthly interpolation of urban tree LAI (mirrors SatellitePhenologyMod)
    integer             :: months(2)             ! the two months to interpolate between (1-12)
    real(r8)            :: timwt(2)              ! interpolation weights (sum to 1)
    integer, parameter  :: nzcanm = 5            ! max vertical levels at urban resolution (for view-factor array allocation)
    !---------------------------------------------------------------------
        
    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begl = bounds%begl; endl = bounds%endl
    begg = bounds%begg; endg = bounds%endg

    ! Allocate urbanparams data structure

    if ( nlevurb > 0 )then
       allocate(this%tk_wall          (begl:endl,nlevurb))  ; this%tk_wall             (:,:) = nan
       allocate(this%tk_roof          (begl:endl,nlevurb))  ; this%tk_roof             (:,:) = nan
       allocate(this%cv_wall          (begl:endl,nlevurb))  ; this%cv_wall             (:,:) = nan
       allocate(this%cv_roof          (begl:endl,nlevurb))  ; this%cv_roof             (:,:) = nan
    end if
    allocate(this%t_building_min      (begl:endl))          ; this%t_building_min      (:)   = nan
    allocate(this%tk_improad          (begl:endl,nlevurb))  ; this%tk_improad          (:,:) = nan
    allocate(this%cv_improad          (begl:endl,nlevurb))  ; this%cv_improad          (:,:) = nan
    allocate(this%thick_wall          (begl:endl))          ; this%thick_wall          (:)   = nan
    allocate(this%thick_roof          (begl:endl))          ; this%thick_roof          (:)   = nan
    allocate(this%nlev_improad        (begl:endl))          ; this%nlev_improad        (:)   = huge(1)
    allocate(this%wind_hgt_canyon     (begl:endl))          ; this%wind_hgt_canyon     (:)   = nan
    allocate(this%em_roof             (begl:endl))          ; this%em_roof             (:)   = nan
    allocate(this%em_improad          (begl:endl))          ; this%em_improad          (:)   = nan
    allocate(this%em_perroad          (begl:endl))          ; this%em_perroad          (:)   = nan
    allocate(this%tree_elaisai_per_uroad          (begl:endl))          ; this%tree_elaisai_per_uroad         (:)   = nan
    allocate(this%tree_pft_urb          (begl:endl))          ; this%tree_pft_urb         (:)   = ispval
    allocate(this%tree_lai_per_vega          (begl:endl))          ; this%tree_lai_per_vega         (:)   = nan
    allocate(this%tree_sai_per_vega          (begl:endl))          ; this%tree_sai_per_vega         (:)   = nan
    allocate(this%em_tree_urb          (begl:endl))          ; this%em_tree_urb          (:)   = nan
    allocate(this%em_wall             (begl:endl))          ; this%em_wall             (:)   = nan
    allocate(this%alb_roof_dir        (begl:endl,numrad))   ; this%alb_roof_dir        (:,:) = nan
    allocate(this%alb_roof_dif        (begl:endl,numrad))   ; this%alb_roof_dif        (:,:) = nan    
    allocate(this%alb_improad_dir     (begl:endl,numrad))   ; this%alb_improad_dir     (:,:) = nan       
    allocate(this%alb_perroad_dir     (begl:endl,numrad))   ; this%alb_perroad_dir     (:,:) = nan       
    allocate(this%alb_tree_urb_dir     (begl:endl,numrad))   ; this%alb_tree_urb_dir     (:,:) = nan       
    allocate(this%alb_improad_dif     (begl:endl,numrad))   ; this%alb_improad_dif     (:,:) = nan       
    allocate(this%alb_perroad_dif     (begl:endl,numrad))   ; this%alb_perroad_dif     (:,:) = nan       
    allocate(this%alb_tree_urb_dif     (begl:endl,numrad))   ; this%alb_tree_urb_dif     (:,:) = nan       
    allocate(this%tran_tree_urb_dif     (begl:endl,numrad))   ; this%tran_tree_urb_dif     (:,:) = nan       
    allocate(this%tran_tree_urb_dir     (begl:endl,numrad))   ; this%tran_tree_urb_dir     (:,:) = nan       
    allocate(this%alb_wall_dir        (begl:endl,numrad))   ; this%alb_wall_dir        (:,:) = nan    
    allocate(this%alb_wall_dif        (begl:endl,numrad))   ; this%alb_wall_dif        (:,:) = nan
    allocate(this%eflx_traffic_factor (begl:endl))          ; this%eflx_traffic_factor (:)   = nan

    allocate(this%fww1d_out           (begl:endl,nzcanm,nzcanm))          ; this%fww1d_out       (:,:,:) = nan
    allocate(this%fvv1d_out           (begl:endl,nzcanm,nzcanm))          ; this%fvv1d_out       (:,:,:) = nan
    allocate(this%fwv1d_out           (begl:endl,nzcanm,nzcanm))          ; this%fwv1d_out       (:,:,:) = nan
    allocate(this%fvw1d_out           (begl:endl,nzcanm,nzcanm))          ; this%fvw1d_out       (:,:,:) = nan
    allocate(this%fwr1d_out           (begl:endl,nzcanm,nzcanm))          ; this%fwr1d_out       (:,:,:) = nan
    allocate(this%frw1d_out           (begl:endl,nzcanm,nzcanm))          ; this%frw1d_out       (:,:,:) = nan
    allocate(this%fvr1d_out           (begl:endl,nzcanm,nzcanm))          ; this%fvr1d_out       (:,:,:) = nan
    allocate(this%frv1d_out           (begl:endl,nzcanm,nzcanm))          ; this%frv1d_out       (:,:,:) = nan

    allocate(this%fwg1d_out           (begl:endl,nzcanm))     ; this%fwg1d_out       (:,:) = nan 
    allocate(this%fgw1d_out           (begl:endl,nzcanm))     ; this%fgw1d_out       (:,:) = nan 
    allocate(this%fgv1d_out           (begl:endl,nzcanm))     ; this%fgv1d_out       (:,:) = nan 
    allocate(this%fsw1d_out           (begl:endl,nzcanm))     ; this%fsw1d_out       (:,:) = nan 
    allocate(this%fvg1d_out           (begl:endl,nzcanm))     ; this%fvg1d_out       (:,:) = nan 
    allocate(this%fsr1d_out           (begl:endl,nzcanm))     ; this%fsr1d_out       (:,:) = nan 
    allocate(this%fsv1d_out           (begl:endl,nzcanm))     ; this%fsv1d_out       (:,:) = nan 
    allocate(this%fws1d_out           (begl:endl,nzcanm))     ; this%fws1d_out       (:,:) = nan
    allocate(this%fvs1d_out           (begl:endl,nzcanm))     ; this%fvs1d_out       (:,:) = nan
    allocate(this%frs1d_out           (begl:endl,nzcanm))     ; this%frs1d_out       (:,:) = nan

    allocate(this%fsg1d_out           (begl:endl))            ; this%fsg1d_out       (:) = nan   
    allocate(this%fts1d_out           (begl:endl))            ; this%fts1d_out       (:) = nan
        
    allocate(this%kww1d_out           (begl:endl,nzcanm,nzcanm))          ; this%kww1d_out       (:,:,:) = nan
    allocate(this%kvv1d_out           (begl:endl,nzcanm,nzcanm))          ; this%kvv1d_out       (:,:,:) = nan
    allocate(this%kwv1d_out           (begl:endl,nzcanm,nzcanm))          ; this%kwv1d_out       (:,:,:) = nan
    allocate(this%kvw1d_out           (begl:endl,nzcanm,nzcanm))          ; this%kvw1d_out       (:,:,:) = nan
    allocate(this%kwr1d_out           (begl:endl,nzcanm,nzcanm))          ; this%kwr1d_out       (:,:,:) = nan
    allocate(this%krw1d_out           (begl:endl,nzcanm,nzcanm))          ; this%krw1d_out       (:,:,:) = nan
    allocate(this%kvr1d_out           (begl:endl,nzcanm,nzcanm))          ; this%kvr1d_out       (:,:,:) = nan
    allocate(this%krv1d_out           (begl:endl,nzcanm,nzcanm))          ; this%krv1d_out       (:,:,:) = nan

    allocate(this%kwg1d_out           (begl:endl,nzcanm))     ; this%kwg1d_out       (:,:) = nan 
    allocate(this%kgw1d_out           (begl:endl,nzcanm))     ; this%kgw1d_out       (:,:) = nan 
    allocate(this%kgv1d_out           (begl:endl,nzcanm))     ; this%kgv1d_out       (:,:) = nan 
    allocate(this%ksw1d_out           (begl:endl,nzcanm))     ; this%ksw1d_out       (:,:) = nan 
    allocate(this%kvg1d_out           (begl:endl,nzcanm))     ; this%kvg1d_out       (:,:) = nan 
    allocate(this%ksr1d_out           (begl:endl,nzcanm))     ; this%ksr1d_out       (:,:) = nan 
    allocate(this%ksv1d_out           (begl:endl,nzcanm))     ; this%ksv1d_out       (:,:) = nan 
    allocate(this%kws1d_out           (begl:endl,nzcanm))     ; this%kws1d_out       (:,:) = nan
    allocate(this%kvs1d_out           (begl:endl,nzcanm))     ; this%kvs1d_out       (:,:) = nan
    allocate(this%krs1d_out           (begl:endl,nzcanm))     ; this%krs1d_out       (:,:) = nan

    allocate(this%ksg1d_out           (begl:endl))            ; this%ksg1d_out       (:) = nan   
    allocate(this%kts1d_out           (begl:endl))            ; this%kts1d_out       (:) = nan

    ! Allocate and initialize unweighted longwave radiation view factor outputs
    allocate(this%vfww_f_out           (begl:endl,nzcanm,nzcanm))          ; this%vfww_f_out       (:,:,:) = nan
    allocate(this%vfvv_f_out           (begl:endl,nzcanm,nzcanm))          ; this%vfvv_f_out       (:,:,:) = nan
    allocate(this%vfwv_f_out           (begl:endl,nzcanm,nzcanm))          ; this%vfwv_f_out       (:,:,:) = nan
    allocate(this%vfvw_f_out           (begl:endl,nzcanm,nzcanm))          ; this%vfvw_f_out       (:,:,:) = nan
    allocate(this%vfwr_f_out           (begl:endl,nzcanm,nzcanm))          ; this%vfwr_f_out       (:,:,:) = nan
    allocate(this%vfrw_f_out           (begl:endl,nzcanm,nzcanm))          ; this%vfrw_f_out       (:,:,:) = nan
    allocate(this%vfvr_f_out           (begl:endl,nzcanm,nzcanm))          ; this%vfvr_f_out       (:,:,:) = nan
    allocate(this%vfrv_f_out           (begl:endl,nzcanm,nzcanm))          ; this%vfrv_f_out       (:,:,:) = nan
    allocate(this%vfwt_f_out           (begl:endl,nzcanm))                 ; this%vfwt_f_out       (:,:)   = nan
    allocate(this%vftw_f_out           (begl:endl,nzcanm))                 ; this%vftw_f_out       (:,:)   = nan
    allocate(this%vftv_f_out           (begl:endl,nzcanm))                 ; this%vftv_f_out       (:,:)   = nan
    allocate(this%vfvt_f_out           (begl:endl,nzcanm))                 ; this%vfvt_f_out       (:,:)   = nan
    allocate(this%vfsw_f_out           (begl:endl,nzcanm))                 ; this%vfsw_f_out       (:,:)   = nan
    allocate(this%vfsr_f_out           (begl:endl,nzcanm))                 ; this%vfsr_f_out       (:,:)   = nan
    allocate(this%vfsv_f_out           (begl:endl,nzcanm))                 ; this%vfsv_f_out       (:,:)   = nan
    allocate(this%vfst_f_out           (begl:endl))                        ; this%vfst_f_out       (:)     = nan
    allocate(this%svft_f_out           (begl:endl))                        ; this%svft_f_out       (:)     = nan

    ! Allocate and initialize unweighted shortwave radiation view factor outputs
    allocate(this%vfww_k_out           (begl:endl,nzcanm,nzcanm))          ; this%vfww_k_out       (:,:,:) = nan
    allocate(this%vfvv_k_out           (begl:endl,nzcanm,nzcanm))          ; this%vfvv_k_out       (:,:,:) = nan
    allocate(this%vfwv_k_out           (begl:endl,nzcanm,nzcanm))          ; this%vfwv_k_out       (:,:,:) = nan
    allocate(this%vfvw_k_out           (begl:endl,nzcanm,nzcanm))          ; this%vfvw_k_out       (:,:,:) = nan
    allocate(this%vfwr_k_out           (begl:endl,nzcanm,nzcanm))          ; this%vfwr_k_out       (:,:,:) = nan
    allocate(this%vfrw_k_out           (begl:endl,nzcanm,nzcanm))          ; this%vfrw_k_out       (:,:,:) = nan
    allocate(this%vfvr_k_out           (begl:endl,nzcanm,nzcanm))          ; this%vfvr_k_out       (:,:,:) = nan
    allocate(this%vfrv_k_out           (begl:endl,nzcanm,nzcanm))          ; this%vfrv_k_out       (:,:,:) = nan
    allocate(this%vfwt_k_out           (begl:endl,nzcanm))                 ; this%vfwt_k_out       (:,:)   = nan
    allocate(this%vftw_k_out           (begl:endl,nzcanm))                 ; this%vftw_k_out       (:,:)   = nan
    allocate(this%vftv_k_out           (begl:endl,nzcanm))                 ; this%vftv_k_out       (:,:)   = nan
    allocate(this%vfvt_k_out           (begl:endl,nzcanm))                 ; this%vfvt_k_out       (:,:)   = nan
    allocate(this%vfsw_k_out           (begl:endl,nzcanm))                 ; this%vfsw_k_out       (:,:)   = nan
    allocate(this%vfsr_k_out           (begl:endl,nzcanm))                 ; this%vfsr_k_out       (:,:)   = nan
    allocate(this%vfsv_k_out           (begl:endl,nzcanm))                 ; this%vfsv_k_out       (:,:)   = nan
    allocate(this%vfst_k_out           (begl:endl))                        ; this%vfst_k_out       (:)     = nan
    allocate(this%svft_k_out           (begl:endl))                        ; this%svft_k_out       (:)     = nan

    ! Allocate and initialize unweighted sky view factor outputs
    allocate(this%svfw_f_out           (begl:endl,nzcanm))                 ; this%svfw_f_out       (:,:)   = nan
    allocate(this%svfv_f_out           (begl:endl,nzcanm))                 ; this%svfv_f_out       (:,:)   = nan
    allocate(this%svfr_f_out           (begl:endl,nzcanm))                 ; this%svfr_f_out       (:,:)   = nan
    allocate(this%svfw_k_out           (begl:endl,nzcanm))                 ; this%svfw_k_out       (:,:)   = nan
    allocate(this%svfv_k_out           (begl:endl,nzcanm))                 ; this%svfv_k_out       (:,:)   = nan
    allocate(this%svfr_k_out           (begl:endl,nzcanm))                 ; this%svfr_k_out       (:,:)   = nan

    allocate(this%z_d_town_out         (begl:endl))                        ; this%z_d_town_out     (:)     = nan
    allocate(this%z_0_town_out         (begl:endl))                        ; this%z_0_town_out     (:)     = nan
    allocate(this%plan_ai_out         (begl:endl))                        ; this%plan_ai_out     (:)     = nan
    allocate(this%frontal_ai_out         (begl:endl))                        ; this%frontal_ai_out     (:)     = nan


    do l = bounds%begl,bounds%endl

       if (lun%urbpoi(l)) then

          g = lun%gridcell(l)
          dindx = lun%itype(l) - isturb_MIN + 1

          this%wind_hgt_canyon(l) = urbinp%wind_hgt_canyon(g,dindx)
          do ib = 1,numrad
             this%alb_roof_dir   (l,ib) = urbinp%alb_roof_dir   (g,dindx,ib)
             this%alb_roof_dif   (l,ib) = urbinp%alb_roof_dif   (g,dindx,ib)
             this%alb_improad_dir(l,ib) = urbinp%alb_improad_dir(g,dindx,ib)
             this%alb_perroad_dir(l,ib) = urbinp%alb_perroad_dir(g,dindx,ib)
             this%alb_tree_urb_dir(l,ib) = urbinp%alb_tree_urb_dir(g,dindx,ib)
             this%alb_improad_dif(l,ib) = urbinp%alb_improad_dif(g,dindx,ib)
             this%alb_perroad_dif(l,ib) = urbinp%alb_perroad_dif(g,dindx,ib)
             this%alb_tree_urb_dif(l,ib) = urbinp%alb_tree_urb_dif(g,dindx,ib)
             this%tran_tree_urb_dir(l,ib) = urbinp%tran_tree_urb_dir(g,dindx,ib)
             this%tran_tree_urb_dif(l,ib) = urbinp%tran_tree_urb_dif(g,dindx,ib)
             this%alb_wall_dir   (l,ib) = urbinp%alb_wall_dir   (g,dindx,ib)
             this%alb_wall_dif   (l,ib) = urbinp%alb_wall_dif   (g,dindx,ib)
          end do
          this%em_roof   (l) = urbinp%em_roof   (g,dindx)
          this%em_improad(l) = urbinp%em_improad(g,dindx)
          this%em_perroad(l) = urbinp%em_perroad(g,dindx)
          this%em_tree_urb(l) = urbinp%em_tree_urb(g,dindx)
          this%em_wall   (l) = urbinp%em_wall   (g,dindx)

         ! Bracketing months and interpolation weights for the current date
         call GetMonthlyInterpWeights(next_time_step=.false., months=months, timwt=timwt)

         this%tree_lai_per_vega(l) = timwt(1)*urbinp%tree_lai_urb(g,dindx,months(1)) &
                                   + timwt(2)*urbinp%tree_lai_urb(g,dindx,months(2))
         this%tree_sai_per_vega(l) = timwt(1)*urbinp%tree_sai_urb(g,dindx,months(1)) &
                                   + timwt(2)*urbinp%tree_sai_urb(g,dindx,months(2))

         ! initialize tree_elaisai_per_uroad using the tree LAI + SAI with no snow effect. This will be updated in the urban phenology module to account for snow effect on tree LAI and SAI.
         this%tree_elaisai_per_uroad(l) = (this%tree_lai_per_vega(l) + this%tree_sai_per_vega(l))*urbinp%wtroad_tree(g,dindx)
         this%tree_pft_urb(l) = urbinp%tree_pft_urb(g,dindx)

          ! Landunit level initialization for urban wall and roof layers and interfaces
          lun%canyon_hwr(l)   = urbinp%canyon_hwr(g,dindx)
          lun%tree_lai_urb_monthly(l,:)          = urbinp%tree_lai_urb(g,dindx,:)
          lun%tree_sai_urb_monthly(l,:)          = urbinp%tree_sai_urb(g,dindx,:)
          lun%wtroad_tree(l)     = urbinp%wtroad_tree(g,dindx)       
          lun%tree_bht_urb(l)          = urbinp%tree_bht_urb(g,dindx)
          lun%tree_tht_urb(l)     = urbinp%tree_tht_urb(g,dindx)       
          lun%wtroad_perv(l)  = urbinp%wtroad_perv(g,dindx)
          lun%ht_roof(l)      = urbinp%ht_roof(g,dindx)
          !confirm: do we really need to assign this value? Or is it just a initialization?
          ! First initialize ht_can_eff= roof height
          lun%ht_can_eff(l)      = urbinp%ht_roof(g,dindx)
          lun%wtlunit_roof(l) = urbinp%wtlunit_roof(g,dindx)
          lun%wall_to_plan_area_ratio(l) = urbinp%wall_to_plan_area_ratio(g,dindx)
          this%tk_wall(l,:)      = urbinp%tk_wall(g,dindx,:)
          this%tk_roof(l,:)      = urbinp%tk_roof(g,dindx,:)
          this%tk_improad(l,:)   = urbinp%tk_improad(g,dindx,:)
          this%cv_wall(l,:)      = urbinp%cv_wall(g,dindx,:)
          this%cv_roof(l,:)      = urbinp%cv_roof(g,dindx,:)
          this%cv_improad(l,:)   = urbinp%cv_improad(g,dindx,:)
          this%thick_wall(l)     = urbinp%thick_wall(g,dindx)
          this%thick_roof(l)     = urbinp%thick_roof(g,dindx)
          this%nlev_improad(l)   = urbinp%nlev_improad(g,dindx)
          this%t_building_min(l) = urbinp%t_building_min(g,dindx)

          ! confirm: what variables belong to lun% and what belongs to this%
          ! when wtroad_tree is small, LAI should be zero
          ! but it is possible to have zero LAI while non-zero wtroad_tree (winter case)
          if ((lun%wtroad_tree(l)  < 1.e-3_r8) .or. (lun%tree_bht_urb(l) < 0.1_r8) .or. (lun%tree_tht_urb(l) < 0.01_r8)) then 
               lun%wtroad_tree(l)  = 0.0_r8
               lun%tree_lai_urb_monthly(l,:) = 0.0_r8
               lun%tree_sai_urb_monthly(l,:) = 0.0_r8
               this%tree_lai_per_vega(l) = 0._r8
               this%tree_sai_per_vega(l) = 0._r8
               this%tree_elaisai_per_uroad(l) = 0._r8
         end if

         ! LAI below the threshold means no tree, so the SAI is dropped as well
         if ((this%tree_lai_per_vega(l)+this%tree_sai_per_vega(l)) < 0.05_r8) then
               this%tree_lai_per_vega(l) = 0._r8
               this%tree_sai_per_vega(l) = 0._r8
               this%tree_elaisai_per_uroad(l) = 0._r8
         end if

          if ((this%tree_lai_per_vega(l)+this%tree_sai_per_vega(l)) > 0.0_r8 .and. (lun%tree_bht_urb(l) < 0.1_r8)) then
              write (iulog,*) 'The bottom of tree crown should be higher than 0.1 when tree exist. tree crown bottom height, tree LAI =',lun%tree_bht_urb(l),this%tree_lai_per_vega(l)
              write (iulog,*) 'clm model is stopping'
              call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
          end if

          if ((this%tree_lai_per_vega(l)+this%tree_sai_per_vega(l)) > 0.0_r8 .and. (lun%tree_tht_urb(l) < 0.01_r8)) then
              write (iulog,*) 'The tree crown depth should be larger than 0.01 when tree exist. tree crown depth, tree LAI =',lun%tree_tht_urb(l),this%tree_lai_per_vega(l)
              write (iulog,*) 'clm model is stopping'
              call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
          end if

          if (lun%tree_bht_urb(l)+lun%tree_tht_urb(l) > 2.0_r8*lun%ht_roof(l)) then
              write (iulog,*) 'The maximum tree height should be lower than two times of building height. tree height, ht_roof =',lun%tree_bht_urb(l)+lun%tree_tht_urb(l),lun%ht_roof(l)
              write (iulog,*) 'clm model is stopping'
              call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
          end if

          if (lun%tree_bht_urb(l) > lun%ht_roof(l)) then
              write (iulog,*) 'The bottom of tree crown should be lower than building height. tree crown bottom height, ht_roof =',lun%tree_bht_urb(l),lun%ht_roof(l)
              write (iulog,*) 'clm model is stopping'
              call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
          end if
          
          ! Compute Monte Carlo view factors for this urban landunit
          call ComputeViewFactors(this, l)
          
          ! Inferred from Sailor and Lu 2004
          if (urban_traffic) then
             this%eflx_traffic_factor(l) = 3.6_r8 * (lun%canyon_hwr(l)-0.5_r8) + 1.0_r8
          else
             this%eflx_traffic_factor(l) = 0.0_r8
          end if

          if (use_vancouver .or. use_mexicocity) then
             ! Freely evolving
             this%t_building_min(l) = 200.00_r8
          else
             if (urban_hac == urban_hac_off) then
                ! Overwrite values read in from urbinp by freely evolving values
                this%t_building_min(l) = 200.00_r8
             end if
          end if

          !----------------------------------------------------------------------------------
          ! Calculate urban land unit aerodynamic constants using Macdonald (1998) as used in
          ! Grimmond and Oke (1999)
          !----------------------------------------------------------------------------------

          ! Set plan area index to roof fraction (see Figure 2 in Grimmond and Oke (1999))
          ! Restrict this value to that of "real cities" as shown in Figure 1 in Grimmond and Ok (1999).
          ! This is to avoid unrealistic combinations of plan and frontal area index that can generate
          ! extremely small values of roughness length that result in stagnant urban canopy air and very
          ! high daytime temperatures.
          ! This is from Keith's note - we need to keep an eye on this as we add more urban land units and more complex urban forms. 
          !  Jackson data occasionally resulted in unrealistic combinations of plan and frontal area index that can generate extremely small values of roughness length ! that result in stagnant urban canopy air and very high daytime temperatures. 
          ! U-surf may have better-behaved data
          ! plan_ai = max(min(lun%wtlunit_roof(l), 0.62_r8), 0.1_r8) 

          ! Derive frontal_ai from equations (7) and (8) in Porson et al. (2010).
          ! See also equation (1) in Masson et al. (2020). The 2/rpi
          ! factor accounts for street orientational averaging.
          ! Restrict this value to that of "real cities" as shown in Figure 1 in Grimmond and Ok (1999).
          ! frontal_ai = max(2._r8 * (1._r8 - plan_ai) * lun%canyon_hwr(l)/ rpi, 0.05_r8)

          ! Use plan_ai = roof fraction. See notes from 6-9-21
         !  plan_ai = lun%wtlunit_roof(l)


          ! Restrict this value to that of "real cities" as shown in Figure 1 in Grimmond and Ok (1999).
          ! This is to avoid unrealistic combinations of plan and frontal area index that can generate
          ! extremely small values of roughness length that result in stagnant urban canopy air and very
          ! high daytime temperatures.

          plan_ai = max(min(lun%wtlunit_roof(l), 0.62_r8), 0.1_r8)


          ! Use relationship derived from Porson (2010) and Masson (2020)
          frontal_ai = lun%wall_to_plan_area_ratio(l)/rpi
          
          this%plan_ai_out(l)=plan_ai
          this%frontal_ai_out(l)=frontal_ai

          ! Calculate aerodynamic displacement height and roughness length
          call ComputeUrbanRoughness(this, l)

       else ! Not urban point 

          this%eflx_traffic_factor(l) = spval
          this%t_building_min(l) = spval

          this%fww1d_out(l,:,:)      = spval
          this%fvv1d_out(l,:,:)      = spval
          this%fwv1d_out(l,:,:)      = spval
          this%fvw1d_out(l,:,:)      = spval
          this%fwr1d_out(l,:,:)      = spval
          this%frw1d_out(l,:,:)      = spval
          this%fvr1d_out(l,:,:)      = spval
          this%frv1d_out(l,:,:)      = spval
          
          this%fwg1d_out(l,:)      = spval
          this%fgw1d_out(l,:)      = spval
          this%fgv1d_out(l,:)      = spval
          this%fsw1d_out(l,:)      = spval
          this%fvg1d_out(l,:)      = spval
          this%fsr1d_out(l,:)      = spval
          this%fsv1d_out(l,:)      = spval
          this%fws1d_out(l,:)      = spval
          this%fvs1d_out(l,:)      = spval
          this%frs1d_out(l,:)      = spval
          this%fts1d_out(l)        = spval 
          this%fsg1d_out(l)        = spval    

          this%kww1d_out(l,:,:)      = spval
          this%kvv1d_out(l,:,:)      = spval
          this%kwv1d_out(l,:,:)      = spval
          this%kvw1d_out(l,:,:)      = spval
          this%kwr1d_out(l,:,:)      = spval
          this%krw1d_out(l,:,:)      = spval
          this%kvr1d_out(l,:,:)      = spval
          this%krv1d_out(l,:,:)      = spval

          this%kwg1d_out(l,:)      = spval
          this%kgw1d_out(l,:)      = spval
          this%kgv1d_out(l,:)      = spval
          this%ksw1d_out(l,:)      = spval
          this%kvg1d_out(l,:)      = spval
          this%ksr1d_out(l,:)      = spval
          this%ksv1d_out(l,:)      = spval
          this%kws1d_out(l,:)      = spval
          this%kvs1d_out(l,:)      = spval
          this%krs1d_out(l,:)      = spval
          this%kts1d_out(l)        = spval
          this%ksg1d_out(l)        = spval  
          
          ! Initialize unweighted longwave radiation view factor outputs
          this%vfww_f_out(l,:,:)      = spval
          this%vfvv_f_out(l,:,:)      = spval
          this%vfwv_f_out(l,:,:)      = spval
          this%vfvw_f_out(l,:,:)      = spval
          this%vfwr_f_out(l,:,:)      = spval
          this%vfrw_f_out(l,:,:)      = spval
          this%vfvr_f_out(l,:,:)      = spval
          this%vfrv_f_out(l,:,:)      = spval
          this%vfwt_f_out(l,:)        = spval
          this%vftw_f_out(l,:)        = spval
          this%vftv_f_out(l,:)        = spval
          this%vfvt_f_out(l,:)        = spval
          this%vfsw_f_out(l,:)        = spval
          this%vfsr_f_out(l,:)        = spval
          this%vfsv_f_out(l,:)        = spval
          this%vfst_f_out(l)          = spval
          this%svft_f_out(l)          = spval

          ! Initialize unweighted shortwave radiation view factor outputs
          this%vfww_k_out(l,:,:)      = spval
          this%vfvv_k_out(l,:,:)      = spval
          this%vfwv_k_out(l,:,:)      = spval
          this%vfvw_k_out(l,:,:)      = spval
          this%vfwr_k_out(l,:,:)      = spval
          this%vfrw_k_out(l,:,:)      = spval
          this%vfvr_k_out(l,:,:)      = spval
          this%vfrv_k_out(l,:,:)      = spval
          this%vfwt_k_out(l,:)        = spval
          this%vftw_k_out(l,:)        = spval
          this%vftv_k_out(l,:)        = spval
          this%vfvt_k_out(l,:)        = spval
          this%vfsw_k_out(l,:)        = spval
          this%vfsr_k_out(l,:)        = spval
          this%vfsv_k_out(l,:)        = spval
          this%vfst_k_out(l)          = spval
          this%svft_k_out(l)          = spval

          ! Initialize unweighted sky view factor outputs
          this%svfw_f_out(l,:)        = spval
          this%svfv_f_out(l,:)        = spval
          this%svfr_f_out(l,:)        = spval
          this%svfw_k_out(l,:)        = spval
          this%svfv_k_out(l,:)        = spval
          this%svfr_k_out(l,:)        = spval

          this%plan_ai_out(l)        = spval
          this%frontal_ai_out(l)        = spval
          this%z_d_town_out(l)        = spval
          this%z_0_town_out(l)        = spval
          
       end if
    end do

    ! Note that we don't deallocate memory for urbinp datatype (call UrbanInput with
    ! mode='finalize') because the arrays are needed for dynamic urban landunits.
    
  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine ComputeViewFactors(this, l)
    !
    ! !DESCRIPTION:
    ! Compute the Monte Carlo radiative view factors for one urban landunit l,
    ! given its canyon geometry and (already interpolated) urban tree LAI
    ! (this%tree_elaisai_per_uroad). Sets up the single-layer urban canyon, derives
    ! the leaf area density and clumping, calls montecarlo_view_factors, and
    ! stores the resulting view factors into this%*_out(l,...).
    !
    ! !ARGUMENTS:
    class(urbanparams_type) , intent(inout) :: this
    integer                 , intent(in)    :: l     ! landunit index
    !
    ! !LOCAL VARIABLES:
    integer, parameter  :: nzcanm = 5        ! Maximum number of vertical levels at urban resolution
    integer             :: maxbhind          ! Maximum vertical layer for highest building or tree
    integer             :: maxind            ! Index of the highest level with roofs or the highest layer with tree foliage, whichever is higher
    real(r8)            :: lad(nzcanm)       ! Canyon-space-average leaf area density [m2 m-3] = LAI per road area / crown
                                             ! depth. NOT the in-crown value (montecarlo forms that as lad/wtroad_tree).
    real(r8)            :: lads(nzcanm)      ! As lad, for shortwave calcs; set equal to lad below
    real(r8)            :: ladl(nzcanm)      ! As lad, for longwave calcs; set equal to lad below
    real(r8)            :: omega(nzcanm)             ! Leaf clumping index 
    real(r8)            :: dzcan             ! Height of buildings [m]
    real(r8)            :: pb_in(nzcanm+1)     ! Probability to have a building with an height equal or higher than each level
    real(r8)            :: ss_in(nzcanm+1)     ! Roof fraction at each level 
    real(r8)            :: wcan              ! Width of the canyons [m]
    real(r8)            :: wbui              ! Width of the buildings [m]
    real(r8)            :: wtroad_tree          ! Canopy cover within the canyon (between-building) space
    real(r8)            :: dray              ! Ray step, normalized by dzcan [-] (NOT metres)
    real(r8)            :: hsky              ! Sky-ray start height as a multiple of the domain top; must be > 1
    integer             :: nrays             ! Points generated on the unit sphere; a foliage layer fires nrays-2 rays.
    integer             :: nsky              ! Number of sky-ray start bins per level (a multiplier, not a ray count)
    real(r8)            :: h1                ! Crown bottom height above ground [m]
    real(r8)            :: h2                ! Crown DEPTH [m], not crown-top height

    !Output of view factor calculation
    !------
    ! Area-weighted longwave view factors
    real(r8)            :: fww1d(nzcanm, nzcanm)  ! Longwave radiation view factor from wall to wall 
    real(r8)            :: fvv1d(nzcanm, nzcanm)  ! Longwave radiation view factor from vegetation to vegetation 
    real(r8)            :: fwv1d(nzcanm, nzcanm)  ! Longwave radiation view factor from wall to vegetation 
    real(r8)            :: fvw1d(nzcanm, nzcanm)  ! Longwave radiation view factor from vegetation to wall 
    real(r8)            :: fwr1d(nzcanm, nzcanm)  ! Longwave radiation view factor from wall to roof 
    real(r8)            :: frw1d(nzcanm, nzcanm)  ! Longwave radiation view factor from roof to wall 
    real(r8)            :: fvr1d(nzcanm, nzcanm)  ! Longwave radiation view factor from vegetation to roof 
    real(r8)            :: frv1d(nzcanm, nzcanm)  ! Longwave radiation view factor from roof to vegetation 
    real(r8)            :: fwg1d(nzcanm)          ! Longwave radiation view factor from wall to ground 
    real(r8)            :: fgw1d(nzcanm)          ! Longwave radiation view factor from ground to wall 
    real(r8)            :: fgv1d(nzcanm)          ! Longwave radiation view factor from ground to vegetation 
    real(r8)            :: fsw1d(nzcanm)          ! Longwave radiation view factor from sky to wall 
    real(r8)            :: fvg1d(nzcanm)          ! Longwave radiation view factor from vegetation to ground 
    real(r8)            :: fsr1d(nzcanm)          ! Longwave radiation view factor from sky to roof 
    real(r8)            :: fsv1d(nzcanm)          ! Longwave radiation view factor from sky to vegetation 
    real(r8)            :: fws1d(nzcanm)          ! Longwave radiation view factor from wall to sky 
    real(r8)            :: fvs1d(nzcanm)          ! Longwave radiation view factor from vegetation to sky 
    real(r8)            :: frs1d(nzcanm)          ! Longwave radiation view factor from roof to sky 
    real(r8)            :: fsg1d                  ! Longwave radiation view factor from sky to ground         
    real(r8)            :: fts1d                  ! Longwave radiation view factor from ground to sky 
    
    ! Area-weighted shortwave view factors
    real(r8)            :: kww1d(nzcanm, nzcanm)  ! Shortwave radiation view factor from wall to wall 
    real(r8)            :: kvv1d(nzcanm, nzcanm)  ! Shortwave radiation view factor from vegetation to vegetation 
    real(r8)            :: kwv1d(nzcanm, nzcanm)  ! Shortwave radiation view factor from wall to vegetation 
    real(r8)            :: kvw1d(nzcanm, nzcanm)  ! Shortwave radiation view factor from vegetation to wall 
    real(r8)            :: kwr1d(nzcanm, nzcanm)  ! Shortwave radiation view factor from wall to roof 
    real(r8)            :: krw1d(nzcanm, nzcanm)  ! Shortwave radiation view factor from roof to wall 
    real(r8)            :: kvr1d(nzcanm, nzcanm)  ! Shortwave radiation view factor from vegetation to roof 
    real(r8)            :: krv1d(nzcanm, nzcanm)  ! Shortwave radiation view factor from roof to vegetation 
    real(r8)            :: kwg1d(nzcanm)          ! Shortwave radiation view factor from wall to ground 
    real(r8)            :: kgw1d(nzcanm)          ! Shortwave radiation view factor from ground to wall 
    real(r8)            :: kgv1d(nzcanm)          ! Shortwave radiation view factor from ground to vegetation 
    real(r8)            :: ksw1d(nzcanm)          ! Shortwave radiation view factor from sky to wall 
    real(r8)            :: kvg1d(nzcanm)          ! Shortwave radiation view factor from vegetation to ground 
    real(r8)            :: ksr1d(nzcanm)          ! Shortwave radiation view factor from sky to roof 
    real(r8)            :: ksv1d(nzcanm)          ! Shortwave radiation view factor from sky to vegetation 
    real(r8)            :: kws1d(nzcanm)          ! Shortwave radiation view factor from wall to sky 
    real(r8)            :: kvs1d(nzcanm)          ! Shortwave radiation view factor from vegetation to sky 
    real(r8)            :: krs1d(nzcanm)          ! Shortwave radiation view factor from roof to sky     
    real(r8)            :: ksg1d                  ! Shortwave radiation view factor from sky to ground         
    real(r8)            :: kts1d                  ! Shortwave radiation view factor from ground to sky   
    
    ! Unweighted longwave radiation view factors (for initialization)
    real(r8)            :: vfww_f(nzcanm,nzcanm)      ! Unweighted longwave radiation view factor from wall to wall
    real(r8)            :: vfvv_f(nzcanm,nzcanm)      ! Unweighted longwave radiation view factor from vegetation to vegetation
    real(r8)            :: vfwv_f(nzcanm,nzcanm)      ! Unweighted longwave radiation view factor from wall to vegetation
    real(r8)            :: vfvw_f(nzcanm,nzcanm)      ! Unweighted longwave radiation view factor from vegetation to wall
    real(r8)            :: vfwr_f(nzcanm,nzcanm)      ! Unweighted longwave radiation view factor from wall to roof
    real(r8)            :: vfrw_f(nzcanm,nzcanm)      ! Unweighted longwave radiation view factor from roof to wall
    real(r8)            :: vfvr_f(nzcanm,nzcanm)      ! Unweighted longwave radiation view factor from vegetation to roof
    real(r8)            :: vfrv_f(nzcanm,nzcanm)      ! Unweighted longwave radiation view factor from roof to vegetation
    real(r8)            :: vfwt_f(nzcanm)             ! Unweighted longwave radiation view factor from wall to ground
    real(r8)            :: vftw_f(nzcanm)             ! Unweighted longwave radiation view factor from ground to wall
    real(r8)            :: vftv_f(nzcanm)             ! Unweighted longwave radiation view factor from ground to vegetation
    real(r8)            :: vfvt_f(nzcanm)             ! Unweighted longwave radiation view factor from vegetation to ground
    real(r8)            :: vfsw_f(nzcanm)             ! Unweighted longwave radiation view factor from sky to wall
    real(r8)            :: vfsr_f(nzcanm)             ! Unweighted longwave radiation view factor from sky to roof
    real(r8)            :: vfsv_f(nzcanm)             ! Unweighted longwave radiation view factor from sky to vegetation
    real(r8)            :: vfst_f                     ! Unweighted longwave radiation view factor from sky to ground    
    real(r8)            :: svft_f                     ! Unweighted longwave sky view factor of the ground (ground -> sky)

    ! Unweighted shortwave radiation view factors (for initialization)
    real(r8)            :: vfww_k(nzcanm,nzcanm)      ! Unweighted shortwave radiation view factor from wall to wall
    real(r8)            :: vfvv_k(nzcanm,nzcanm)      ! Unweighted shortwave radiation view factor from vegetation to vegetation
    real(r8)            :: vfwv_k(nzcanm,nzcanm)      ! Unweighted shortwave radiation view factor from wall to vegetation
    real(r8)            :: vfvw_k(nzcanm,nzcanm)      ! Unweighted shortwave radiation view factor from vegetation to wall
    real(r8)            :: vfwr_k(nzcanm,nzcanm)      ! Unweighted shortwave radiation view factor from wall to roof
    real(r8)            :: vfrw_k(nzcanm,nzcanm)      ! Unweighted shortwave radiation view factor from roof to wall
    real(r8)            :: vfvr_k(nzcanm,nzcanm)      ! Unweighted shortwave radiation view factor from vegetation to roof
    real(r8)            :: vfrv_k(nzcanm,nzcanm)      ! Unweighted shortwave radiation view factor from roof to vegetation
    real(r8)            :: vfwt_k(nzcanm)             ! Unweighted shortwave radiation view factor from wall to ground
    real(r8)            :: vftw_k(nzcanm)             ! Unweighted shortwave radiation view factor from ground to wall
    real(r8)            :: vftv_k(nzcanm)             ! Unweighted shortwave radiation view factor from ground to vegetation
    real(r8)            :: vfvt_k(nzcanm)             ! Unweighted shortwave radiation view factor from vegetation to ground
    real(r8)            :: vfsw_k(nzcanm)             ! Unweighted shortwave radiation view factor from sky to wall
    real(r8)            :: vfsr_k(nzcanm)             ! Unweighted shortwave radiation view factor from sky to roof
    real(r8)            :: vfsv_k(nzcanm)             ! Unweighted shortwave radiation view factor from sky to vegetation
    real(r8)            :: vfst_k                     ! Unweighted shortwave radiation view factor from sky to ground
    real(r8)            :: svft_k                     ! Unweighted shortwave sky view factor of the ground (ground -> sky)
    
    ! Unweighted sky view factors (for initialization)
    real(r8)            :: svfw_f(nzcanm)             ! Unweighted longwave sky view factor for wall
    real(r8)            :: svfv_f(nzcanm)             ! Unweighted longwave sky view factor for vegetation
    real(r8)            :: svfr_f(nzcanm)             ! Unweighted longwave sky view factor for roof
    real(r8)            :: svfw_k(nzcanm)             ! Unweighted shortwave sky view factor for wall
    real(r8)            :: svfv_k(nzcanm)             ! Unweighted shortwave sky view factor for vegetation
    real(r8)            :: svfr_k(nzcanm)             ! Unweighted shortwave sky view factor for roof
    !-----------------------------------------------------------------------
    ! Single-layer urban canyon setup for the view factor calculation
    maxbhind=2
    maxind=2
    ss_in=0._r8
    pb_in=0._r8
    lad=0._r8
    wtroad_tree=0._r8         
    omega=0._r8           
    dray=0._r8 
    h1=0.0_r8
    h2=0.0_r8
    ! ss(k) = pb(k) - pb(k+1)
    ! ss_in(1) = 0 -- no roofs at ground level
    ss_in(2)=1._r8 ! 100% of roofs at z = dzcan
    pb_in(1)=1._r8
    pb_in(2)=1._r8
    
    ! Could change these values if view factors from sky are not accurate, especially nsky
    nsky=1
    hsky=1.5_r8

   dzcan=lun%ht_roof(l)
   wcan=lun%ht_roof(l)/lun%canyon_hwr(l)  
   ! confirmed this expression      
   wbui = lun%ht_roof(l)/(lun%canyon_hwr(l)*(1._r8-lun%wtlunit_roof(l))/lun%wtlunit_roof(l))

   if (this%tree_elaisai_per_uroad(l) > 1e-6_r8) then
      !Eq. 15 in Krayenhoff et al. 2020
      !Use LAI per road area between buildings, which is what used in the Krayenhoff et al. 2020 paper.
      ! confirmed this equation, and it is consistent with the derivation in the Krayenhoff et al. 2020 paper
      ! omega of all layers are the same
      omega=-1.0_r8 / (0.5_r8 * this%tree_elaisai_per_uroad(l) ) * log(1.0_r8 - lun%wtroad_tree(l) * &
            (1.0_r8 - exp(-0.5_r8 * this%tree_elaisai_per_uroad(l)/lun%wtroad_tree(l))))
   else
      omega = 0._r8
   endif     

   ! nrays=50000
   ! dray=0.05_r8*min(min(dzcan,wcan),wbui)/dzcan
   nrays=NAN
   ! dray=NAN*min(min(dzcan,wcan),wbui)/dzcan
   ! Now also add tree crown depth to scale dray
   dray=NAN*min(min(min(dzcan,wcan),wbui),lun%tree_tht_urb(l))/dzcan

   ! nrays=5000
   ! dray=0.05_r8*min(min(dzcan,wcan),wbui)/dzcan
   if (this%tree_elaisai_per_uroad(l) > 1e-6_r8) then
      if ((lun%tree_bht_urb(l)+lun%tree_tht_urb(l))<= lun%ht_roof(l)) then
         lad(1)=max(this%tree_elaisai_per_uroad(l)/lun%tree_tht_urb(l), 1e-6_r8)
         lad(2) = 0._r8
         ! A_v1 unit check: wcan [m] * lad [m2 m-3] * omega [-] * crown depth [m] * 2 = [m], i.e. two-sided,
         ! clumping-weighted leaf area per unit canyon LENGTH. Consistent with A_vs/A_vl inside
         ! montecarlo_view_factors, which are the same quantity divided by dzcan.
         lun%A_v1(l)=max(wcan*lad(1)*omega(1)*lun%tree_tht_urb(l)*2._r8, 1e-6_r8)
         lun%A_v2(l)=0._r8 
      else if ((lun%tree_bht_urb(l)+lun%tree_tht_urb(l)) > lun%ht_roof(l)) then
         ! the lad in the first and second layer should be the same and should be in respect to the realistic tree span
         lad(1)=max(this%tree_elaisai_per_uroad(l)/lun%tree_tht_urb(l), 1e-6_r8)
         lad(2) =max(this%tree_elaisai_per_uroad(l)/lun%tree_tht_urb(l), 1e-6_r8)
         lun%A_v1(l)=max(wcan*lad(1)*omega(1)*(lun%ht_roof(l)-lun%tree_bht_urb(l))*2._r8, 1e-6_r8)
         lun%A_v2(l)=max(wcan*lad(2)*omega(2)*(lun%tree_bht_urb(l)+lun%tree_tht_urb(l)-lun%ht_roof(l))*2._r8, 1e-6_r8)   
      end if  
   else
      lad(1) = 0.0_r8
      lad(2) = 0.0_r8
      lun%A_v1(l) = 0.0_r8   ! Floor only 
      lun%A_v2(l) = 0.0_r8
   end if
   ! These are unused but keep for now:
   lads=lad  ! for shortwave calcs (usually equal to "lad")
   ladl=lad  ! lfor longwave calcs (usually equal to "lad")
   
   wtroad_tree=lun%wtroad_tree(l)
   h1=lun%tree_bht_urb(l)
   h2=lun%tree_tht_urb(l)
   
   ! calculate view factor
   ! checked the input sequence
   call montecarlo_view_factors(nzcanm,dzcan,wcan,wbui,&
         wtroad_tree,lad,lads,ladl,this%tree_elaisai_per_uroad(l),omega,ss_in,pb_in,dray,maxind,&
                  maxbhind,nrays,nsky,hsky,h1,h2,&
                  fww1d,fvv1d,fwv1d,fvw1d,fwr1d,frw1d,fvr1d,&
                  frv1d,fwg1d,fgw1d,fgv1d,fsw1d,fvg1d,fsg1d,fsr1d,&
                  fsv1d,kww1d,kvv1d,kwv1d,kvw1d,kwr1d,krw1d,kvr1d,&
                  krv1d,kwg1d,kgw1d,kgv1d,ksw1d,kvg1d,ksg1d,ksr1d,&
                  ksv1d,kws1d,kvs1d,kts1d,krs1d,fws1d,fvs1d,fts1d,frs1d,&
                  vfww_f,vfvv_f,vfwv_f,vfvw_f,vfwr_f,vfrw_f,vfvr_f,&
                  vfrv_f,vfwt_f,vftw_f,vftv_f,vfsw_f,vfvt_f,vfst_f,vfsr_f,&
                  vfsv_f,vfww_k,vfvv_k,vfwv_k,vfvw_k,vfwr_k,vfrw_k,vfvr_k,&
                  vfrv_k,vfwt_k,vftw_k,vftv_k,vfsw_k,vfvt_k,vfst_k,vfsr_k,&
                  vfsv_k,svfw_k,svfv_k,svft_k,svfr_k,svfw_f,svfv_f,svft_f,svfr_f,l)
   this%fww1d_out(l,:,:)      = fww1d(:,:)
   this%fvv1d_out(l,:,:)      = fvv1d(:,:)
   this%fwv1d_out(l,:,:)      = fwv1d(:,:)
   this%fvw1d_out(l,:,:)      = fvw1d(:,:)
   this%fwr1d_out(l,:,:)      = fwr1d(:,:)
   this%frw1d_out(l,:,:)      = frw1d(:,:)
   this%fvr1d_out(l,:,:)      = fvr1d(:,:)
   this%frv1d_out(l,:,:)      = frv1d(:,:)
   
   this%fwg1d_out(l,:)      = fwg1d(:)
   this%fgw1d_out(l,:)      = fgw1d(:)
   this%fgv1d_out(l,:)      = fgv1d(:)
   this%fsw1d_out(l,:)      = fsw1d(:)
   this%fvg1d_out(l,:)      = fvg1d(:)
   this%fsr1d_out(l,:)      = fsr1d(:)
   this%fsv1d_out(l,:)      = fsv1d(:)
   this%fws1d_out(l,:)      = fws1d(:)
   this%fvs1d_out(l,:)      = fvs1d(:)
   this%frs1d_out(l,:)      = frs1d(:) 
   this%fsg1d_out(l)        = fsg1d
   this%fts1d_out(l)        = fts1d        
   
   this%kww1d_out(l,:,:)      = kww1d(:,:)
   this%kvv1d_out(l,:,:)      = kvv1d(:,:)
   this%kwv1d_out(l,:,:)      = kwv1d(:,:)
   this%kvw1d_out(l,:,:)      = kvw1d(:,:)
   this%kwr1d_out(l,:,:)      = kwr1d(:,:)
   this%krw1d_out(l,:,:)      = krw1d(:,:)
   this%kvr1d_out(l,:,:)      = kvr1d(:,:)
   this%krv1d_out(l,:,:)      = krv1d(:,:)
   
   this%kwg1d_out(l,:)      = kwg1d(:)
   this%kgw1d_out(l,:)      = kgw1d(:)
   this%kgv1d_out(l,:)      = kgv1d(:)
   this%ksw1d_out(l,:)      = ksw1d(:)
   this%kvg1d_out(l,:)      = kvg1d(:)
   this%ksr1d_out(l,:)      = ksr1d(:)
   this%ksv1d_out(l,:)      = ksv1d(:)
   this%kws1d_out(l,:)      = kws1d(:)
   this%kvs1d_out(l,:)      = kvs1d(:)
   this%krs1d_out(l,:)      = krs1d(:) 
   this%ksg1d_out(l)        = ksg1d
   this%kts1d_out(l)        = kts1d   
   
   this%vfww_f_out(l,:,:)      = vfww_f(:,:)
   this%vfvv_f_out(l,:,:)      = vfvv_f(:,:)
   this%vfwv_f_out(l,:,:)      = vfwv_f(:,:)
   this%vfvw_f_out(l,:,:)      = vfvw_f(:,:)
   this%vfwr_f_out(l,:,:)      = vfwr_f(:,:)
   this%vfrw_f_out(l,:,:)      = vfrw_f(:,:)
   this%vfvr_f_out(l,:,:)      = vfvr_f(:,:)
   this%vfrv_f_out(l,:,:)      = vfrv_f(:,:)
   this%vfwt_f_out(l,:)        = vfwt_f(:)
   this%vftw_f_out(l,:)        = vftw_f(:)
   this%vftv_f_out(l,:)        = vftv_f(:)
   this%vfvt_f_out(l,:)        = vfvt_f(:)
   this%vfsw_f_out(l,:)        = vfsw_f(:)
   this%vfsr_f_out(l,:)        = vfsr_f(:)
   this%vfsv_f_out(l,:)        = vfsv_f(:)
   this%vfst_f_out(l)          = vfst_f
   

   this%vfww_k_out(l,:,:)      = vfww_k(:,:)
   this%vfvv_k_out(l,:,:)      = vfvv_k(:,:)
   this%vfwv_k_out(l,:,:)      = vfwv_k(:,:)
   this%vfvw_k_out(l,:,:)      = vfvw_k(:,:)
   this%vfwr_k_out(l,:,:)      = vfwr_k(:,:)
   this%vfrw_k_out(l,:,:)      = vfrw_k(:,:)
   this%vfvr_k_out(l,:,:)      = vfvr_k(:,:)
   this%vfrv_k_out(l,:,:)      = vfrv_k(:,:)
   this%vfwt_k_out(l,:)        = vfwt_k(:)
   this%vftw_k_out(l,:)        = vftw_k(:)
   this%vftv_k_out(l,:)        = vftv_k(:)
   this%vfvt_k_out(l,:)        = vfvt_k(:)
   this%vfsw_k_out(l,:)        = vfsw_k(:)
   this%vfsr_k_out(l,:)        = vfsr_k(:)
   this%vfsv_k_out(l,:)        = vfsv_k(:)
   this%vfst_k_out(l)          = vfst_k

   this%svfw_f_out(l,:)        = svfw_f(:)
   this%svfv_f_out(l,:)        = svfv_f(:)
   this%svfr_f_out(l,:)        = svfr_f(:)
   this%svfw_k_out(l,:)        = svfw_k(:)
   this%svfv_k_out(l,:)        = svfv_k(:)
   this%svfr_k_out(l,:)        = svfr_k(:)

   this%svft_f_out(l)          = svft_f
   this%svft_k_out(l)          = svft_k

  end subroutine ComputeViewFactors

  !-----------------------------------------------------------------------
  subroutine ComputeUrbanRoughness(this, l)
    !
    ! !DESCRIPTION:
    ! Compute the aerodynamic displacement height (z_d_town) and roughness
    ! length (z_0_town) for one urban landunit l, following Macdonald (1998)
    ! as used in Grimmond and Oke (1999), including the tree effect via the
    ! (already interpolated) urban tree LAI (this%tree_elaisai_per_uroad).
    !
    ! !USES:
    use clm_varctl , only : use_vancouver, use_mexicocity
    use clm_varcon , only : vkc
    !
    ! !ARGUMENTS:
    class(urbanparams_type) , intent(inout) :: this
    integer                 , intent(in)    :: l     ! landunit index
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: alpha = 4.43_r8 ! coefficient used to calculate z_d_town
    real(r8), parameter :: beta = 1.0_r8   ! coefficient used to calculate z_d_town
    real(r8), parameter :: C_d = 1.2_r8    ! drag coefficient as used in Grimmond and Oke (1999)
    real(r8)            :: k_opt   ! canopy light extinction parameter
    real(r8)            :: p_2d    ! optical porosity (direct beam transmission through canopy)
    real(r8)            :: p_3d    ! volumetric/aerodynamic porosity
    real(r8)            :: A_pv    ! tree plan area
    real(r8)            :: A_pb    ! building plan area
    real(r8)            :: A_tot    ! total plan area
    real(r8)            :: r_tree  ! tree radius
    real(r8)            :: ht_tree  ! tree height
    real(r8)            :: bv_drag_ratio    ! ratio between vegetation drag CDv and building drag CDb
    ! TODO: confirm the usage of sheltered vs unsheltered frontal area of building and tree
    real(r8)            :: frontal_b_unsh ! unsheltered frontal area of building
    real(r8)            :: frontal_v_unsh ! unsheltered frontal area of tree
    real(r8)            :: frontal_b ! sheltered frontal area of building
    real(r8)            :: frontal_v ! sheltered frontal area of tree
    real(r8)            :: plan_ai_eff ! effective plan area index
    !-----------------------------------------------------------------------
   ! Calculate displacement heightt
   if (use_vancouver) then
      lun%z_d_town(l) = 3.5_r8
   else if (use_mexicocity) then
      lun%z_d_town(l) = 10.9_r8
   else
      ! The same assumption as in UT&C
      r_tree=0.25_r8*lun%wtroad_tree(l)*lun%ht_roof(l)/lun%canyon_hwr(l)
      ht_tree=lun%tree_bht_urb(l)+lun%tree_tht_urb(l)
      k_opt=0.5_r8
      ! doublechek: it is wrong to use tree_elaisai_per_uroad (LAI per road area); for tree porosity, it should be tree_lai_per_vega (LAI per vegetation area)

      ! think: should we use LAI or LAI+SAI?
      p_2d=exp(-k_opt*(this%tree_lai_per_vega(l)+this%tree_sai_per_vega(l)))
      ! p_2d=exp(-k_opt*(this%tree_lai_per_vega(l)))
      p_3d=p_2d**0.4_r8
      A_pv=4.0_r8*r_tree ! essentially it is lun%wtroad_tree(l)*lun%ht_roof(l)/lun%canyon_hwr(l)
      A_pb=lun%ht_roof(l)/lun%canyon_hwr(l)*lun%wtlunit_roof(l)/(1-lun%wtlunit_roof(l))
      A_tot=A_pb+lun%ht_roof(l)/lun%canyon_hwr(l)
      bv_drag_ratio = (-1.251_r8*p_3d**2.0_r8+0.489_r8*p_3d+0.803_r8)/C_d
      lun%ht_can_eff(l) = (lun%ht_roof(l)*A_pb+ht_tree*(1.0_r8-p_3d)*A_pv)/(A_pb+(1.0_r8-p_3d)*A_pv)

      ! The unsheltered and actual frontal area may be swaped in Meili's SI
      ! I *think* unsheltered fromtal area is the part above z-zd, and actual frontal area is ht_roof
      ! there is no need to calculate the unsheltered frontal area but I'll keep them for now
      ! It's best to confirm with Kent and Meili about the definition of unsheltered and actual frontal area.
      plan_ai_eff = (A_pb+(1.0_r8-p_3d)*A_pv)/A_tot
      lun%z_d_town(l) = (1._r8 + alpha**(-plan_ai_eff) * (plan_ai_eff - 1._r8)) * lun%ht_can_eff(l)

      ! I think the frontal area of tree can also be adjusted
      frontal_b = lun%ht_can_eff(l)
      frontal_v = 2_r8 * r_tree ! may use crown depth tree_tht?

      frontal_b_unsh = frontal_b * (lun%ht_can_eff(l) - lun%z_d_town(l))/lun%ht_can_eff(l)
      frontal_v_unsh = frontal_b * (lun%ht_can_eff(l) - lun%z_d_town(l))/lun%ht_can_eff(l)
   end if

   ! Calculate the roughness length
   if (use_vancouver) then
      lun%z_0_town(l) = 0.35_r8
   else if (use_mexicocity) then
      lun%z_0_town(l) = 2.2_r8
   else
      lun%z_0_town(l) = lun%ht_can_eff(l) * (1._r8 - lun%z_d_town(l) / lun%ht_can_eff(l)) * &
         exp(-(1/(vkc**2)*0.5_r8*beta*C_d*(1.0_r8 - lun%z_d_town(l)/lun%ht_can_eff(l))*(frontal_b+bv_drag_ratio*frontal_v)/A_tot)**(-0.5_r8))
   end if

   this%z_d_town_out(l)=lun%z_d_town(l)
   this%z_0_town_out(l)=lun%z_0_town(l)
  end subroutine ComputeUrbanRoughness

  !-----------------------------------------------------------------------
  subroutine UrbanPhenology(bounds, num_urbantreep, filter_urbantreep, urbanparams_inst, &
       canopystate_inst, waterdiagnosticbulk_inst)
    !
    ! !DESCRIPTION:
    ! Time-interpolate the prescribed monthly urban road-tree LAI and SAI
    ! (lun%tree_lai_urb_monthly, lun%tree_sai_urb_monthly) to the current date
    ! and update the per-vegetation-area tree LAI and SAI (tree_lai_per_vega,
    ! tree_sai_per_vega) and the per-road-area tree LAI+SAI (tree_elaisai_per_uroad). Intended to be called every month or every week (decide later) so the urban
    ! tree LAI follows its seasonal cycle, analogous to SatellitePhenology for
    ! natural vegetation. Monthly values are treated as valid at mid-month.
    !
    ! !USES:
    use pftconMod               , only : noveg, nbrdlf_dcd_brl_shrub
    use PatchType               , only : patch
    use CanopyStateType         , only : canopystate_type
    use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_urbantreep       ! number of urban tree patches in filter
    integer                        , intent(in)    :: filter_urbantreep(:) ! urban tree patch filter
    type(urbanparams_type)         , intent(inout) :: urbanparams_inst
    type(canopystate_type)         , intent(inout) :: canopystate_inst
    type(waterdiagnosticbulk_type) , intent(in)    :: waterdiagnosticbulk_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: f                     ! filter index
    integer  :: p, c               ! patch, column indices
    integer  :: l                     ! landunit index
    integer  :: months(2)             ! the two months to interpolate between (1-12)
    real(r8) :: timwt(2)              ! interpolation weights (sum to 1)
    real(r8) :: ol                    ! thickness of canopy layer covered by snow (m)
    real(r8) :: fb                    ! fraction of canopy layer remaining exposed above snow
    !-----------------------------------------------------------------------

    associate(                                             &
         tree_lai_urb_monthly => lun%tree_lai_urb_monthly, & ! Input:  [real(r8) (:,:)] monthly LAI of road tree (m2/m2)
         tree_sai_urb_monthly => lun%tree_sai_urb_monthly, & ! Input:  [real(r8) (:,:)] monthly SAI of road tree (m2/m2)
         wtroad_tree          => lun%wtroad_tree         , & ! Input:  [real(r8) (:)  ] weight of road tree column to total road
         tree_lai_per_vega    => urbanparams_inst%tree_lai_per_vega  , & ! Output: [real(r8) (:)  ] interpolated LAI per vegetation area
         tree_sai_per_vega    => urbanparams_inst%tree_sai_per_vega  , & ! Output: [real(r8) (:)  ] interpolated SAI per vegetation area
         tree_elaisai_per_uroad   => urbanparams_inst%tree_elaisai_per_uroad   ,& ! Output: [real(r8) (:)  ] interpolated LAI+SAI per road area
         h1            => lun%tree_bht_urb                      , & ! Input:  [real(r8) (:)   ]  the height of tree crown base
         h2            => lun%tree_tht_urb                       , & ! Input:  [real(r8) (:)   ]  the depth of tree crown    
         frac_sno           => waterdiagnosticbulk_inst%frac_sno_col   , & ! Input:  [real(r8) (:) ] fraction of ground covered by snow (0 to 1)
         snow_depth         => waterdiagnosticbulk_inst%snow_depth_col , & ! Input:  [real(r8) (:) ] snow height (m)
         tlai               => canopystate_inst%tlai_patch    ,          & ! Output: [real(r8) (:) ] one-sided leaf area index, no burying by snow
         tsai               => canopystate_inst%tsai_patch    ,          & ! Output: [real(r8) (:) ] one-sided stem area index, no burying by snow
         elai               => canopystate_inst%elai_patch    ,          & ! Output: [real(r8) (:) ] one-sided leaf area index with burying by snow
         esai               => canopystate_inst%esai_patch    ,          & ! Output: [real(r8) (:) ] one-sided stem area index with burying by snow
         htop               => canopystate_inst%htop_patch    ,          & ! Output: [real(r8) (:) ] canopy top (m)
         hbot               => canopystate_inst%hbot_patch    ,          & ! Output: [real(r8) (:) ] canopy bottom (m)
         frac_veg_nosno_alb => canopystate_inst%frac_veg_nosno_alb_patch & ! Output: [integer  (:) ] fraction of vegetation not covered by snow (0 OR 1) [-]
         )

      ! Bracketing months and interpolation weights for the upcoming time step.
      ! These depend only on the date, so compute them once for all landunits.
      call GetMonthlyInterpWeights(next_time_step=.true., months=months, timwt=timwt)

      ! another thing worthing doing: call UrbanPhenology once during urbanparamater initialization like ATES-SP already calls SatellitePhenology at init
      ! Think about restart: elai/esai come back from the restart file with SAI and snow burial included. Init recomputes tree_elaisai_per_uroad as LAI only, so the two disagree until the first radiation step.
      ! it is okay to solve this later
      do f = 1, num_urbantreep
         p = filter_urbantreep(f)
         c = patch%column(p)
         l = patch%landunit(p)

         ! urbanparams has a guardian statement to make sure the tree_lai_urb_monthly=0 and
         ! tree_sai_urb_monthly=0 when wtroad_tree(l)/h1(l) is too small
         tree_lai_per_vega(l)  = timwt(1)*tree_lai_urb_monthly(l,months(1)) &
                                 + timwt(2)*tree_lai_urb_monthly(l,months(2))
         tree_sai_per_vega(l)  = timwt(1)*tree_sai_urb_monthly(l,months(1)) &
                                 + timwt(2)*tree_sai_urb_monthly(l,months(2))

         ! Discuss this threshold; whether to apply it for urban trees
         ! LAI below the threshold means no tree, so the SAI is dropped as well
         if ((tree_lai_per_vega(l)+tree_sai_per_vega(l)) < 0.05_r8) then
            tree_lai_per_vega(l) = 0._r8
            tree_sai_per_vega(l) = 0._r8
            tlai(p) = 0.0_r8
            tsai(p) = 0.0_r8
         else
            tsai(p) = tree_sai_per_vega(l) ! SAI per tree area (vegetated area)
            tlai(p) = tree_lai_per_vega(l) ! LAI per tree area (vegetated area)
         end if

         htop(p) =  h1(l) + h2(l)
         hbot(p) =  h1(l)
      
         if ((tree_lai_per_vega(l) > 0.0_r8) .and. (hbot(p)< 0.1_r8)) then
            write (iulog,*) 'The bottom of tree crown should be higher than 0.1 when tree exist. tree crown bottom height, tree LAI =',hbot(p),tree_lai_per_vega(l)
            write (iulog,*) 'clm model is stopping'
            call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
         end if

        ! adjust lai and sai for burying by snow. if exposed lai and sai
         ! are less than 0.05, set equal to zero to prevent numerical
         ! problems associated with very small lai and sai.

         ! snow burial fraction for short vegetation (e.g. grasses, crops) changes with vegetation height
         ! accounts for a 20% bending factor, as used in Lombardozzi et al. (2018) GRL 45(18), 9889-9897

         ! NOTE: The following snow burial code is duplicated in CNVegStructUpdateMod.
         ! Changes in one place should be accompanied by similar changes in the other.

         if (patch%itype(p) > noveg .and. patch%itype(p) <= nbrdlf_dcd_brl_shrub ) then
            ol = min( max(snow_depth(c)-hbot(p), 0._r8), htop(p)-hbot(p))
            fb = 1._r8 - ol / max(1.e-06_r8, htop(p)-hbot(p))
         else
            fb = 1._r8 - (max(min(snow_depth(c),max(0.05_r8,htop(p)*0.8_r8)),0._r8)/(max(0.05_r8,htop(p)*0.8_r8)))
         endif

         ! area weight by snow covered fraction
         ! we do this no matter what use_fates_sp is

         ! if(.not.use_fates_sp)then
         ! Do not set these in FATES_SP mode as they turn on the 'vegsol' filter and also
         ! are duplicated by the FATE variables (in the FATES IFP indexing space)
         
         elai(p) = max(tlai(p)*(1.0_r8 - frac_sno(c)) + tlai(p)*fb*frac_sno(c), 0.0_r8)
         esai(p) = max(tsai(p)*(1.0_r8 - frac_sno(c)) + tsai(p)*fb*frac_sno(c), 0.0_r8)

         ! when tlai < 0.05, elai is always < 0.05 and set to 0 
         if (elai(p) < 0.05_r8) elai(p) = 0._r8
         if (esai(p) < 0.05_r8) esai(p) = 0._r8

         tree_elaisai_per_uroad(l) = (elai(p) + esai(p))*wtroad_tree(l)
         ! Fraction of vegetation free of snow

         if ((elai(p) + esai(p)) >= 0.05_r8) then
            frac_veg_nosno_alb(p) = 1
         else
            frac_veg_nosno_alb(p) = 0
         end if

         ! The urban tree LAI (tree_elaisai_per_uroad) has just been updated, so
         ! recompute the radiative view factors and the aerodynamic parameters
         ! (displacement height, roughness length) for this landunit, both of
         ! which depend on the tree LAI.
         call ComputeViewFactors(urbanparams_inst, l)
         call ComputeUrbanRoughness(urbanparams_inst, l)
      end do

    end associate

  end subroutine UrbanPhenology
  !-----------------------------------------------------------------------
  subroutine GetMonthlyInterpWeights(next_time_step, months, timwt)
    !
    ! !DESCRIPTION:
    ! Compute the two bracketing months and their linear time-interpolation
    ! weights for mid-month-anchored monthly data, following the scheme in
    ! SatellitePhenologyMod. Monthly values are treated as valid at mid-month.
    ! If next_time_step is .true., the weights/months are for the end of the
    ! current time step (nstep+1); otherwise for the current date.
    !
    ! !USES:
    use clm_time_manager , only : get_curr_date, get_step_size_real
    !
    ! !ARGUMENTS:
    logical , intent(in)  :: next_time_step  ! .true. => offset by one time step
    integer , intent(out) :: months(2)       ! the two months to interpolate between (1-12)
    real(r8), intent(out) :: timwt(2)         ! interpolation weights (sum to 1)
    !
    ! !LOCAL VARIABLES:
    integer  :: kyr, kmo, kda, ksec   ! year, month, day, seconds of current date
    integer  :: it(2)                 ! bracketing month offsets (0/1 and 1/2)
    real(r8) :: t                     ! fractional position through current month
    real(r8) :: dtime                 ! land model time step (s)
    integer, parameter :: ndaypm(12) = &
         (/31,28,31,30,31,30,31,31,30,31,30,31/) ! days per month
    !-----------------------------------------------------------------------

    if (next_time_step) then
       dtime = get_step_size_real()
       call get_curr_date(kyr, kmo, kda, ksec, offset=int(dtime))
    else
       call get_curr_date(kyr, kmo, kda, ksec)
    end if

    ! Fractional position through the current month, measured at mid-day
    t = (kda-0.5_r8) / ndaypm(kmo)
    ! if t < 0.5, it(1)=0 (first half of month); otherwise it(1)=1
    it(1) = int(t + 0.5_r8)
    it(2) = it(1) + 1
    ! months(1) and months(2) are the two months to interpolate between
    months(1) = kmo + it(1) - 1
    months(2) = kmo + it(2) - 1
    if (months(1) <  1) months(1) = 12
    if (months(2) > 12) months(2) = 1
    ! two linear interpolation weights (always summing to 1)
    timwt(1) = (it(1)+0.5_r8) - t
    timwt(2) = 1._r8 - timwt(1)

  end subroutine GetMonthlyInterpWeights

  !-----------------------------------------------------------------------
  subroutine UrbanInput(begg, endg, mode)
    !
    ! !DESCRIPTION: 
    ! Allocate memory and read in urban input data
    !
    ! !USES:
    use clm_varpar      , only : numrad, nlevurb, nmonth
    use landunit_varcon , only : numurbl
    use fileutils       , only : getavu, relavu, getfil, opnfil
    use spmdMod         , only : masterproc
    use domainMod       , only : ldomain
    use ncdio_pio       , only : file_desc_t, ncd_io, ncd_inqvdlen, ncd_inqfdims
    use ncdio_pio       , only : ncd_pio_openfile, ncd_pio_closefile, ncd_inqdid, ncd_inqdlen
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: begg, endg
    character(len=*), intent(in) :: mode
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn      ! local file name
    type(file_desc_t)  :: ncid       ! netcdf id
    integer :: dimid                 ! netCDF id
    integer :: nw,n,k,i,j,ni,nj,ns   ! indices
    real(r8), pointer :: tree_lai_urb_1mon(:,:) ! single-month slice of TREE_LAI_URB [gridcell, numurbl]
    real(r8), pointer :: tree_sai_urb_1mon(:,:) ! single-month slice of TREE_SAI_URB [gridcell, numurbl]
                                                ! (pointer, not allocatable: ncd_io's 2d dummy is a pointer)
    integer :: nlevurb_i             ! input grid: number of urban vertical levels
    integer :: numrad_i              ! input grid: number of solar bands (VIS/NIR)
    integer :: numurbl_i             ! input grid: number of urban landunits
    integer :: ier,ret               ! error status
    logical :: isgrid2d              ! true => file is 2d 
    logical :: readvar               ! true => variable is on dataset
    logical :: has_numurbl           ! true => numurbl dimension is on dataset
    character(len=32) :: subname = 'UrbanInput' ! subroutine name
    !-----------------------------------------------------------------------

    if ( nlevurb == 0 ) return

    if (mode == 'initialize') then

       ! Read urban data
       
       if (masterproc) then
          write(iulog,*)' Reading in urban input data from fsurdat file ...'
       end if
       
       call getfil (fsurdat, locfn, 0)
       call ncd_pio_openfile (ncid, locfn, 0)

       if (masterproc) then
          write(iulog,*) subname,trim(fsurdat)
       end if

       ! Check whether this file has new-format urban data
       call ncd_inqdid(ncid, 'numurbl', dimid, dimexist=has_numurbl)

       ! If file doesn't have numurbl, then it is old-format urban;
       ! in this case, set nlevurb to zero
       if (.not. has_numurbl) then
         nlevurb = 0
         if (masterproc) write(iulog,*)'PCT_URBAN is not multi-density, nlevurb set to 0'
       end if

       if ( nlevurb == 0 ) return

       ! Allocate dynamic memory
!-------------------[kz.7]Ray tracing test-------------------------     
       allocate(urbinp%canyon_hwr(begg:endg, numurbl), &  
                urbinp%tree_lai_urb(begg:endg, numurbl, nmonth), &
                urbinp%tree_sai_urb(begg:endg, numurbl, nmonth), &
                urbinp%tree_pft_urb(begg:endg, numurbl), &
                urbinp%wtroad_tree(begg:endg, numurbl), &     
                urbinp%tree_bht_urb(begg:endg, numurbl), &     
                urbinp%tree_tht_urb(begg:endg, numurbl), &     
                !-------------------[kz.7]Ray tracing test-------------------------                                      
                urbinp%wtlunit_roof(begg:endg, numurbl), &  
                urbinp%wtroad_perv(begg:endg, numurbl), &
                urbinp%wall_to_plan_area_ratio(begg:endg, numurbl), &
                urbinp%em_roof(begg:endg, numurbl), &     
                urbinp%em_improad(begg:endg, numurbl), &    
                urbinp%em_perroad(begg:endg, numurbl), &    
                urbinp%em_tree_urb(begg:endg, numurbl), &    
                urbinp%em_wall(begg:endg, numurbl), &    
                urbinp%alb_roof_dir(begg:endg, numurbl, numrad), &    
                urbinp%alb_roof_dif(begg:endg, numurbl, numrad), &    
                urbinp%alb_improad_dir(begg:endg, numurbl, numrad), &    
                urbinp%alb_perroad_dir(begg:endg, numurbl, numrad), &    
                urbinp%alb_tree_urb_dir(begg:endg, numurbl, numrad), &    
                urbinp%alb_improad_dif(begg:endg, numurbl, numrad), &    
                urbinp%alb_perroad_dif(begg:endg, numurbl, numrad), &    
                urbinp%alb_tree_urb_dif(begg:endg, numurbl, numrad), &    
                urbinp%alb_wall_dir(begg:endg, numurbl, numrad), &    
                urbinp%alb_wall_dif(begg:endg, numurbl, numrad), &
                urbinp%tran_tree_urb_dif(begg:endg, numurbl, numrad), &    
                urbinp%tran_tree_urb_dir(begg:endg, numurbl, numrad), &    
                urbinp%ht_roof(begg:endg, numurbl), &
                urbinp%wind_hgt_canyon(begg:endg, numurbl), &
                urbinp%tk_wall(begg:endg, numurbl,nlevurb), &
                urbinp%tk_roof(begg:endg, numurbl,nlevurb), &
                urbinp%tk_improad(begg:endg, numurbl,nlevurb), &
                urbinp%cv_wall(begg:endg, numurbl,nlevurb), &
                urbinp%cv_roof(begg:endg, numurbl,nlevurb), &
                urbinp%cv_improad(begg:endg, numurbl,nlevurb), &
                urbinp%thick_wall(begg:endg, numurbl), &
                urbinp%thick_roof(begg:endg, numurbl), &
                urbinp%nlev_improad(begg:endg, numurbl), &
                urbinp%t_building_min(begg:endg, numurbl), &
                stat=ier)

       if (ier /= 0) then
          call endrun(msg="Allocation error "//errmsg(sourcefile, __LINE__))
       endif

       ! TREE_LAI_URB and TREE_SAI_URB are read one month at a time below (they have a
       ! record dimension on the fsurdat file). Initialize to nan so that a month that never
       ! gets filled shows up immediately rather than being silently used as whatever was in memory.
       urbinp%tree_lai_urb(:,:,:) = nan
       urbinp%tree_sai_urb(:,:,:) = nan

       call ncd_inqfdims (ncid, isgrid2d, ni, nj, ns)
       if (ldomain%ns /= ns .or. ldomain%ni /= ni .or. ldomain%nj /= nj) then
          write(iulog,*)trim(subname), 'ldomain and input file do not match dims '
          write(iulog,*)trim(subname), 'ldomain%ni,ni,= ',ldomain%ni,ni
          write(iulog,*)trim(subname), 'ldomain%nj,nj,= ',ldomain%nj,nj
          write(iulog,*)trim(subname), 'ldomain%ns,ns,= ',ldomain%ns,ns
          call endrun(msg=errmsg(sourcefile, __LINE__))
       end if

       call ncd_inqdid(ncid, 'nlevurb', dimid)
       call ncd_inqdlen(ncid, dimid, nlevurb_i)
       if (nlevurb_i /= nlevurb) then
          write(iulog,*)trim(subname)// ': parameter nlevurb= ',nlevurb, &
               'does not equal input dataset nlevurb= ',nlevurb_i
          call endrun(msg=errmsg(sourcefile, __LINE__))
       endif

       call ncd_inqdid(ncid, 'numrad', dimid)
       call ncd_inqdlen(ncid, dimid, numrad_i)
       if (numrad_i /= numrad) then
          write(iulog,*)trim(subname)// ': parameter numrad= ',numrad, &
               'does not equal input dataset numrad= ',numrad_i
          call endrun(msg=errmsg(sourcefile, __LINE__))
       endif
       call ncd_inqdid(ncid, 'numurbl', dimid)
       call ncd_inqdlen(ncid, dimid, numurbl_i)
       if (numurbl_i /= numurbl) then
          write(iulog,*)trim(subname)// ': parameter numurbl= ',numurbl, &
               'does not equal input dataset numurbl= ',numurbl_i
          call endrun(msg=errmsg(sourcefile, __LINE__))
       endif
!-------------------[kz.8]Ray tracing test-------------------------
       ! TREE_LAI_URB and TREE_SAI_URB are dimensioned (time, numurbl, lat, lon) on the
       ! fsurdat file, with time as the record dimension. Record variables must be read one
       ! slice at a time with an explicit nt index - a single call without nt fills only the
       ! first month and leaves the rest untouched. This follows the MONTHLY_LAI and
       ! MONTHLY_SAI reads in SatellitePhenologyMod.
       allocate(tree_lai_urb_1mon(begg:endg, numurbl), &
                tree_sai_urb_1mon(begg:endg, numurbl), stat=ier)
       if (ier /= 0) then
          call endrun(msg="Allocation error "//errmsg(sourcefile, __LINE__))
       endif

       do k = 1, nmonth
          call ncd_io(ncid=ncid, varname='TREE_LAI_URB', flag='read', data=tree_lai_urb_1mon,&
              dim1name=grlnd, nt=k, readvar=readvar)
          if (.not. readvar) then
            call endrun( msg='ERROR: TREE_LAI_URB NOT on fsurdat file '//errmsg(sourcefile, __LINE__))
          end if
          urbinp%tree_lai_urb(begg:endg, :, k) = tree_lai_urb_1mon(begg:endg, :)

          call ncd_io(ncid=ncid, varname='TREE_SAI_URB', flag='read', data=tree_sai_urb_1mon,&
              dim1name=grlnd, nt=k, readvar=readvar)
          if (.not. readvar) then
            call endrun( msg='ERROR: TREE_SAI_URB NOT on fsurdat file '//errmsg(sourcefile, __LINE__))
          end if
          urbinp%tree_sai_urb(begg:endg, :, k) = tree_sai_urb_1mon(begg:endg, :)
       end do

       deallocate(tree_lai_urb_1mon, tree_sai_urb_1mon)

       call ncd_io(ncid=ncid, varname='TREE_PFT_URB', flag='read', data=urbinp%tree_pft_urb,&
           dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
         call endrun( msg='ERROR: TREE_PFT_URB NOT on fsurdat file '//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='TREE_BHT_URB', flag='read', data=urbinp%tree_bht_urb,&
           dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
         call endrun( msg='ERROR: TREE_BHT_URB NOT on fsurdat file '//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='TREE_THT_URB', flag='read', data=urbinp%tree_tht_urb,&
           dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
         call endrun( msg='ERROR: TREE_THT_URB NOT on fsurdat file '//errmsg(sourcefile, __LINE__))
       end if
!-------------------[kz.8]Ray tracing test-------------------------        
       call ncd_io(ncid=ncid, varname='CANYON_HWR', flag='read', data=urbinp%canyon_hwr,&
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg='ERROR: CANYON_HWR NOT on fsurdat file '//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='WTLUNIT_ROOF', flag='read', data=urbinp%wtlunit_roof, &
            dim1name=grlnd,  readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: WTLUNIT_ROOF NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='WTROAD_PERV', flag='read', data=urbinp%wtroad_perv, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: WTROAD_PERV NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='WALL_TO_PLAN_AREA_RATIO', flag='read', data=urbinp%wall_to_plan_area_ratio, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: WALL_TO_PLAN_AREA_RATIO NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='WTROAD_TREE', flag='read', data=urbinp%wtroad_tree, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: WTROAD_TREE NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='EM_ROOF', flag='read', data=urbinp%em_roof, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: EM_ROOF NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='EM_IMPROAD', flag='read', data=urbinp%em_improad, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: EM_IMPROAD NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='EM_PERROAD', flag='read', data=urbinp%em_perroad, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: EM_PERROAD NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='EM_TREE_URB', flag='read', data=urbinp%em_tree_urb, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: EM_TREE_URB NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='EM_WALL', flag='read', data=urbinp%em_wall, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: EM_WALL NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='HT_ROOF', flag='read', data=urbinp%ht_roof, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: HT_ROOF NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='WIND_HGT_CANYON', flag='read', data=urbinp%wind_hgt_canyon, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: WIND_HGT_CANYON NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='THICK_WALL', flag='read', data=urbinp%thick_wall, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: THICK_WALL NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='THICK_ROOF', flag='read', data=urbinp%thick_roof, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: THICK_ROOF NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='NLEV_IMPROAD', flag='read', data=urbinp%nlev_improad, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: NLEV_IMPROAD NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='T_BUILDING_MIN', flag='read', data=urbinp%t_building_min, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: T_BUILDING_MIN NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='ALB_IMPROAD_DIR', flag='read', data=urbinp%alb_improad_dir, &
            dim1name=grlnd, readvar=readvar)
       if (.not.readvar) then
          call endrun( msg=' ERROR: ALB_IMPROAD_DIR NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='ALB_IMPROAD_DIF', flag='read', data=urbinp%alb_improad_dif, &
            dim1name=grlnd, readvar=readvar)
       if (.not.readvar) then
          call endrun( msg=' ERROR: ALB_IMPROAD_DIF NOT on fsurdat file'//errmsg(sourcefile, __LINE__) )
       end if

       call ncd_io(ncid=ncid, varname='ALB_PERROAD_DIR', flag='read',data=urbinp%alb_perroad_dir, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: ALB_PERROAD_DIR NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='ALB_PERROAD_DIF', flag='read',data=urbinp%alb_perroad_dif, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: ALB_PERROAD_DIF NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='ALB_TREE_URB_DIR', flag='read',data=urbinp%alb_tree_urb_dir, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: ALB_TREE_URB_DIR NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='ALB_TREE_URB_DIF', flag='read',data=urbinp%alb_tree_urb_dif, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: ALB_TREE_URB_DIF NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

              call ncd_io(ncid=ncid, varname='TRAN_TREE_URB_DIR', flag='read',data=urbinp%tran_tree_urb_dir, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: TRAN_TREE_URB_DIR NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='TRAN_TREE_URB_DIF', flag='read',data=urbinp%tran_tree_urb_dif, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: TRAN_TREE_URB_DIF NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='ALB_ROOF_DIR', flag='read', data=urbinp%alb_roof_dir,  &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: ALB_ROOF_DIR NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='ALB_ROOF_DIF', flag='read', data=urbinp%alb_roof_dif,  &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: ALB_ROOF_DIF NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='ALB_WALL_DIR', flag='read', data=urbinp%alb_wall_dir, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: ALB_WALL_DIR NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='ALB_WALL_DIF', flag='read', data=urbinp%alb_wall_dif, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: ALB_WALL_DIF NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='TK_IMPROAD', flag='read', data=urbinp%tk_improad, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: TK_IMPROAD NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='TK_ROOF', flag='read', data=urbinp%tk_roof, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: TK_ROOF NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='TK_WALL', flag='read', data=urbinp%tk_wall, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: TK_WALL NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='CV_IMPROAD', flag='read', data=urbinp%cv_improad, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: CV_IMPROAD NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='CV_ROOF', flag='read', data=urbinp%cv_roof, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: CV_ROOF NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='CV_WALL', flag='read', data=urbinp%cv_wall, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: CV_WALL NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_pio_closefile(ncid)
       if (masterproc) then
          write(iulog,*)' Sucessfully read urban input data' 
          write(iulog,*)
       end if

    else if (mode == 'finalize') then

       if ( nlevurb == 0 ) return       
       deallocate(urbinp%canyon_hwr, &
!-------------------[kz.9]Ray tracing test-------------------------        
                  urbinp%tree_lai_urb, &
                  urbinp%tree_sai_urb, &
                  urbinp%tree_pft_urb, &
                  urbinp%wtroad_tree, & 
                  urbinp%tree_bht_urb, &    
                  urbinp%tree_tht_urb, &     
                  urbinp%wall_to_plan_area_ratio, &
!-------------------[kz.9]Ray tracing test-------------------------                    
                  urbinp%wtlunit_roof, &
                  urbinp%wtroad_perv, &
                  urbinp%em_roof, &
                  urbinp%em_improad, &
                  urbinp%em_perroad, &
                  urbinp%em_tree_urb, &
                  urbinp%em_wall, &
                  urbinp%alb_roof_dir, &
                  urbinp%alb_roof_dif, &
                  urbinp%alb_improad_dir, &
                  urbinp%alb_perroad_dir, &
                  urbinp%alb_tree_urb_dir, &
                  urbinp%alb_improad_dif, &
                  urbinp%alb_perroad_dif, &
                  urbinp%alb_tree_urb_dif, &
                  urbinp%alb_wall_dir, &
                  urbinp%alb_wall_dif, &
                  urbinp%tran_tree_urb_dif, &
                  urbinp%tran_tree_urb_dir, &
                  urbinp%ht_roof, &
                  urbinp%wind_hgt_canyon, &
                  urbinp%tk_wall, &
                  urbinp%tk_roof, &
                  urbinp%tk_improad, &
                  urbinp%cv_wall, &
                  urbinp%cv_roof, &
                  urbinp%cv_improad, &
                  urbinp%thick_wall, &
                  urbinp%thick_roof, &
                  urbinp%nlev_improad, &
                  urbinp%t_building_min, &
                  stat=ier)
       if (ier /= 0) then
          call endrun(msg='initUrbanInput: deallocation error '//errmsg(sourcefile, __LINE__))
       end if
    else
       write(iulog,*)'initUrbanInput error: mode ',trim(mode),' not supported '
       call endrun(msg=errmsg(sourcefile, __LINE__))
    end if

  end subroutine UrbanInput

  !-----------------------------------------------------------------------
  subroutine CheckUrban(begg, endg, pcturb, caller)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Confirm that we have valid urban data for all points with pct urban > 0. If this isn't
    ! true, abort with a message.
    !
    ! !USES:
    use clm_instur      , only : urban_valid
    use landunit_varcon , only : numurbl
    !
    ! !ARGUMENTS:
    implicit none
    integer         , intent(in) :: begg, endg           ! beg & end grid cell indices
    real(r8)        , intent(in) :: pcturb(begg:,:)      ! % urban
    character(len=*), intent(in) :: caller               ! identifier of caller, for more meaningful error messages
    !
    ! !REVISION HISTORY:
    ! Created by Bill Sacks 7/2013, mostly by moving code from surfrd_special
    !
    ! !LOCAL VARIABLES:
    logical :: found
    integer :: nl, n
    integer :: nindx, dindx
    integer :: nlev
    !-----------------------------------------------------------------------

    found = .false.
    do nl = begg,endg
       do n = 1, numurbl
          if ( pcturb(nl,n) > 0.0_r8 ) then
             if ( .not. urban_valid(nl) .or. &
                  urbinp%canyon_hwr(nl,n)            <= 0._r8 .or. &
!-------------------[kz.10]Ray tracing test-------------------------                    
                  any(urbinp%tree_lai_urb(nl,n,:)             < 0._r8) .or. &
                  any(urbinp%tree_sai_urb(nl,n,:)             < 0._r8) .or. &
                  ! note that urban tree pft cannot be 0
                  urbinp%tree_pft_urb(nl,n)              <= 0 .or. &   
                  urbinp%wtroad_tree(nl,n)              < 0._r8 .or. &   
                  urbinp%tree_bht_urb(nl,n)              < 0._r8 .or. &   
                  urbinp%tree_tht_urb(nl,n)              < 0._r8 .or. &   
!-------------------[kz.10]Ray tracing test-------------------------                                     
                  urbinp%em_perroad(nl,n)            <= 0._r8 .or. &
                  urbinp%em_tree_urb(nl,n)            <= 0._r8 .or. &
                  urbinp%em_roof(nl,n)               <= 0._r8 .or. &
                  urbinp%em_wall(nl,n)               <= 0._r8 .or. &
                  urbinp%ht_roof(nl,n)               <= 0._r8 .or. &
                  urbinp%thick_roof(nl,n)            <= 0._r8 .or. &
                  urbinp%thick_wall(nl,n)            <= 0._r8 .or. &
                  urbinp%t_building_min(nl,n)        <= 0._r8 .or. &
                  urbinp%wind_hgt_canyon(nl,n)       <= 0._r8 .or. &
                  urbinp%wtlunit_roof(nl,n)          <= 0._r8 .or. &
                  urbinp%wtroad_perv(nl,n)           < 0._r8 .or. &
                  urbinp%wall_to_plan_area_ratio(nl,n) <= 0._r8 .or. &
                  any(urbinp%alb_improad_dir(nl,n,:) <= 0._r8) .or. &
                  any(urbinp%alb_improad_dif(nl,n,:) <= 0._r8) .or. &
                  any(urbinp%alb_perroad_dir(nl,n,:) <= 0._r8) .or. &
                  any(urbinp%alb_perroad_dif(nl,n,:) <= 0._r8) .or. &
                  any(urbinp%alb_tree_urb_dir(nl,n,:) <= 0._r8) .or. &
                  any(urbinp%alb_tree_urb_dif(nl,n,:) <= 0._r8) .or. &
                  any(urbinp%tran_tree_urb_dir(nl,n,:) <= 0._r8) .or. &
                  any(urbinp%tran_tree_urb_dif(nl,n,:) <= 0._r8) .or. &
                  any(urbinp%alb_roof_dir(nl,n,:)    <= 0._r8) .or. &
                  any(urbinp%alb_roof_dif(nl,n,:)    <= 0._r8) .or. &
                  any(urbinp%alb_wall_dir(nl,n,:)    <= 0._r8) .or. &
                  any(urbinp%alb_wall_dif(nl,n,:)    <= 0._r8) .or. &
                  any(urbinp%tk_roof(nl,n,:)         <= 0._r8) .or. &
                  any(urbinp%tk_wall(nl,n,:)         <= 0._r8) .or. &
                  any(urbinp%cv_roof(nl,n,:)         <= 0._r8) .or. &
                  any(urbinp%cv_wall(nl,n,:)         <= 0._r8)) then
                found = .true.
                nindx = nl
                dindx = n
                exit
             else
                if (urbinp%nlev_improad(nl,n) > 0) then
                   nlev = urbinp%nlev_improad(nl,n)
                   if ( any(urbinp%tk_improad(nl,n,1:nlev) <= 0._r8) .or. &
                        any(urbinp%cv_improad(nl,n,1:nlev) <= 0._r8)) then
                      found = .true.
                      nindx = nl
                      dindx = n
                      exit
                   end if
                end if
             end if
             if (found) exit
          end if
       end do
    end do
    if ( found ) then
       write(iulog,*) trim(caller), ' ERROR: no valid urban data for nl=',nindx
       write(iulog,*)'density type:    ',dindx
       write(iulog,*)'urban_valid:     ',urban_valid(nindx)
       write(iulog,*)'canyon_hwr:      ',urbinp%canyon_hwr(nindx,dindx)
!-------------------[kz.11]Ray tracing test-------------------------         
       write(iulog,*)'tree_lai_urb:             ',urbinp%tree_lai_urb(nindx,dindx,:)
       write(iulog,*)'tree_sai_urb:             ',urbinp%tree_sai_urb(nindx,dindx,:)
        write(iulog,*)'tree_pft_urb:             ',urbinp%tree_pft_urb(nindx,dindx)
       write(iulog,*)'wtroad_tree:        ',urbinp%wtroad_tree(nindx,dindx)   
       write(iulog,*)'tree_bht_urb:        ',urbinp%tree_bht_urb(nindx,dindx)   
       write(iulog,*)'tree_tht_urb:        ',urbinp%tree_tht_urb(nindx,dindx)   
       !-------------------[kz.11]Ray tracing test-------------------------             
       write(iulog,*)'em_improad:      ',urbinp%em_improad(nindx,dindx)
       write(iulog,*)'em_perroad:      ',urbinp%em_perroad(nindx,dindx)
       write(iulog,*)'em_tree_urb:      ',urbinp%em_tree_urb(nindx,dindx)
       write(iulog,*)'em_roof:         ',urbinp%em_roof(nindx,dindx)
       write(iulog,*)'em_wall:         ',urbinp%em_wall(nindx,dindx)
       write(iulog,*)'ht_roof:         ',urbinp%ht_roof(nindx,dindx)
       write(iulog,*)'thick_roof:      ',urbinp%thick_roof(nindx,dindx)
       write(iulog,*)'thick_wall:      ',urbinp%thick_wall(nindx,dindx)
       write(iulog,*)'t_building_min:  ',urbinp%t_building_min(nindx,dindx)
       write(iulog,*)'wind_hgt_canyon: ',urbinp%wind_hgt_canyon(nindx,dindx)
       write(iulog,*)'wtlunit_roof:    ',urbinp%wtlunit_roof(nindx,dindx)
       write(iulog,*)'wtroad_perv:     ',urbinp%wtroad_perv(nindx,dindx)
       write(iulog,*)'wall_to_plan_area_ratio:     ',urbinp%wall_to_plan_area_ratio(nindx,dindx)
       write(iulog,*)'alb_improad_dir: ',urbinp%alb_improad_dir(nindx,dindx,:)
       write(iulog,*)'alb_improad_dif: ',urbinp%alb_improad_dif(nindx,dindx,:)
       write(iulog,*)'alb_perroad_dir: ',urbinp%alb_perroad_dir(nindx,dindx,:)
       write(iulog,*)'alb_perroad_dif: ',urbinp%alb_perroad_dif(nindx,dindx,:)
       write(iulog,*)'alb_tree_urb_dir: ',urbinp%alb_tree_urb_dir(nindx,dindx,:)
       write(iulog,*)'alb_tree_urb_dif: ',urbinp%alb_tree_urb_dif(nindx,dindx,:)
       write(iulog,*)'tran_tree_urb_dir: ',urbinp%tran_tree_urb_dir(nindx,dindx,:)
       write(iulog,*)'tran_tree_urb_dif: ',urbinp%tran_tree_urb_dif(nindx,dindx,:)
       write(iulog,*)'alb_roof_dir:    ',urbinp%alb_roof_dir(nindx,dindx,:)
       write(iulog,*)'alb_roof_dif:    ',urbinp%alb_roof_dif(nindx,dindx,:)
       write(iulog,*)'alb_wall_dir:    ',urbinp%alb_wall_dir(nindx,dindx,:)
       write(iulog,*)'alb_wall_dif:    ',urbinp%alb_wall_dif(nindx,dindx,:)
       write(iulog,*)'tk_roof:         ',urbinp%tk_roof(nindx,dindx,:)
       write(iulog,*)'tk_wall:         ',urbinp%tk_wall(nindx,dindx,:)
       write(iulog,*)'cv_roof:         ',urbinp%cv_roof(nindx,dindx,:)
       write(iulog,*)'cv_wall:         ',urbinp%cv_wall(nindx,dindx,:)
       if (urbinp%nlev_improad(nindx,dindx) > 0) then
          nlev = urbinp%nlev_improad(nindx,dindx)
          write(iulog,*)'tk_improad: ',urbinp%tk_improad(nindx,dindx,1:nlev)
          write(iulog,*)'cv_improad: ',urbinp%cv_improad(nindx,dindx,1:nlev)
       end if
       call endrun(subgrid_index=nindx, subgrid_level=subgrid_level_gridcell, msg=errmsg(sourcefile, __LINE__))
    end if

  end subroutine CheckUrban

  subroutine montecarlo_view_factors(nzcanm,dzcan,wcan,wbui,&
      wtroad_tree,lad,lads,ladl,lai,omega,ss,pb,dray,maxind,maxbhind,n,nsky,hsky,h1,h2,&
      fww1d,fvv1d,fwv1d,fvw1d,fwr1d,frw1d,fvr1d,frv1d,fwg1d,fgw1d,fgv1d,fsw1d,fvg1d,fsg1d,fsr1d,&
      fsv1d,kww1d,kvv1d,kwv1d,kvw1d,kwr1d,krw1d,kvr1d,krv1d,kwg1d,kgw1d,kgv1d,ksw1d,kvg1d,ksg1d,ksr1d,&
      ksv1d,kws1d,kvs1d,kts1d,krs1d,fws1d,fvs1d,fts1d,frs1d,vfww_f,vfvv_f,vfwv_f,vfvw_f,vfwr_f,vfrw_f,vfvr_f,&
      vfrv_f,vfwt_f,vftw_f,vftv_f,vfsw_f,vfvt_f,vfst_f,vfsr_f,vfsv_f,vfww_k,vfvv_k,vfwv_k,vfvw_k,vfwr_k,vfrw_k,vfvr_k,&
      vfrv_k,vfwt_k,vftw_k,vftv_k,vfsw_k,vfvt_k,vfst_k,vfsr_k,vfsv_k,svfw_k,svfv_k,svft_k,svfr_k,svfw_f,svfv_f,svft_f,svfr_f,l)
    !
    ! !DESCRIPTION: 
    !
    ! Calculate the view factors between all 'surfaces' including vegetation layers.
    ! 
    ! !USES:
    !
    ! ARGUMENTS:
    implicit none
    integer         , intent(in) :: nzcanm           ! Maximum number of vertical levels at urban resolution
    real(r8)        , intent(in) :: dzcan            ! Vertical layer thickness [m]; = building height in this single-layer setup.
    real(r8)        , intent(in) :: wcan             ! Width of the canyons [m]
    real(r8)        , intent(in) :: wbui             ! Width of the buildings [m]
    real(r8)        , intent(in) :: wtroad_tree      ! Tree canopy cover as a fraction of the canyon-floor/road area [-]
    real(r8)        , intent(in) :: lad(nzcanm)      ! Canyon-space-average leaf area density of each layer [m2 m-3], i.e. LAI per road
                                                     ! area divided by crown depth. NOT the in-crown value - that is lad_tree below.
    real(r8)        , intent(in) :: lads(nzcanm)     ! As lad, but for shortwave calcs [m2 m-3]. Caller currently sets lads = lad.
    real(r8)        , intent(in) :: lai              ! Exposed tree leaf+stem area index per unit road area [m2 m-2]
                                                     ! (tree_elaisai_per_uroad). Only used as an on/off switch for lad_tree.
    real(r8)        , intent(in) :: ladl(nzcanm)     ! As lad, but for longwave calcs [m2 m-3]. Caller currently sets ladl = lad.
    real(r8)        , intent(in) :: omega(nzcanm)    ! Tree-crown-scale leaf clumping index [-] (Eq. 15 in Krayenhoff et al. 2020)
    real(r8)        , intent(in) :: ss(nzcanm+1)       ! Roof fraction of each level [-]. 
    real(r8)        , intent(in) :: pb(nzcanm+1)       ! Probability to have a building with a height equal to or higher than each level [-].
    real(r8)        , intent(in) :: dray             ! Ray step, normalized by dzcan [-] (NOT metres).
    integer         , intent(in) :: maxind           ! Index of the highest level containing roofs or tree foliage, whichever is higher
    integer         , intent(in) :: maxbhind         ! Index of the highest level containing buildings (highest roof level).
    integer         , intent(in) :: n                ! Number of points generated on the unit sphere. n-2 ray directions are produced;
    integer         , intent(in) :: nsky             ! Number of sky-ray start positions per level (a multiplier, not a ray count);
                                                     ! total rays released from the sky per level = nsky*nk.
    real(r8)        , intent(in) :: hsky             ! Start height of sky rays, as a multiple of the domain-top height (maxind layers). Must be > 1.
    integer         , intent(in) :: l                ! Urban landunit index
    real(r8)        , intent(in) :: h1               ! Crown bottom height above ground [m]; 
    real(r8)        , intent(in) :: h2               ! Crown depth/thickness [m], not crown-top height; 
        
    ! Area-weighted longwave view factors
    real(r8), intent(out) :: fww1d(nzcanm,nzcanm)        ! Longwave view factor from wall to wall 
    real(r8), intent(out) :: fvv1d(nzcanm,nzcanm)        ! Longwave view factor from vegetation to vegetation 
    real(r8), intent(out) :: fwv1d(nzcanm,nzcanm)        ! Longwave view factor from wall to vegetation 
    real(r8), intent(out) :: fvw1d(nzcanm,nzcanm)        ! Longwave view factor from vegetation to wall 
    real(r8), intent(out) :: fwr1d(nzcanm,nzcanm)        ! Longwave view factor from wall to roof 
    real(r8), intent(out) :: frw1d(nzcanm,nzcanm)        ! Longwave view factor from roof to wall 
    real(r8), intent(out) :: fvr1d(nzcanm,nzcanm)        ! Longwave view factor from vegetation to roof 
    real(r8), intent(out) :: frv1d(nzcanm,nzcanm)        ! Longwave view factor from roof to vegetation 
    real(r8), intent(out) :: fwg1d(nzcanm)               ! Longwave view factor from wall to ground 
    real(r8), intent(out) :: fgw1d(nzcanm)               ! Longwave view factor from ground to wall 
    real(r8), intent(out) :: fgv1d(nzcanm)               ! Longwave view factor from ground to vegetation 
    real(r8), intent(out) :: fsw1d(nzcanm)               ! Longwave view factor from sky to wall 
    real(r8), intent(out) :: fvg1d(nzcanm)               ! Longwave view factor from vegetation to ground 
    real(r8), intent(out) :: fsr1d(nzcanm)               ! Longwave view factor from sky to roof 
    real(r8), intent(out) :: fsv1d(nzcanm)               ! Longwave view factor from sky to vegetation 
    real(r8), intent(out) :: fws1d(nzcanm)               ! Longwave view factor from wall to sky
    real(r8), intent(out) :: fvs1d(nzcanm)               ! Longwave view factor from vegetation to sky
    real(r8), intent(out) :: frs1d(nzcanm)               ! Longwave view factor from roof to sky
    real(r8), intent(out) :: fsg1d                       ! Longwave view factor from sky to ground 
    real(r8), intent(out) :: fts1d                       ! Longwave view factor from ground to sky

    ! Area-weighted shortwave view factors
    real(r8), intent(out) :: kww1d(nzcanm,nzcanm)        ! Shortwave view factor from wall to wall 
    real(r8), intent(out) :: kvv1d(nzcanm,nzcanm)        ! Shortwave view factor from vegetation to vegetation 
    real(r8), intent(out) :: kwv1d(nzcanm,nzcanm)        ! Shortwave view factor from wall to vegetation 
    real(r8), intent(out) :: kvw1d(nzcanm,nzcanm)        ! Shortwave view factor from vegetation to wall 
    real(r8), intent(out) :: kwr1d(nzcanm,nzcanm)        ! Shortwave view factor from wall to roof 
    real(r8), intent(out) :: krw1d(nzcanm,nzcanm)        ! Shortwave view factor from roof to wall 
    real(r8), intent(out) :: kvr1d(nzcanm,nzcanm)        ! Shortwave view factor from vegetation to roof 
    real(r8), intent(out) :: krv1d(nzcanm,nzcanm)        ! Shortwave view factor from roof to vegetation 
    real(r8), intent(out) :: kwg1d(nzcanm)               ! Shortwave view factor from wall to ground 
    real(r8), intent(out) :: kgw1d(nzcanm)               ! Shortwave view factor from ground to wall 
    real(r8), intent(out) :: kgv1d(nzcanm)               ! Shortwave view factor from ground to vegetation 
    real(r8), intent(out) :: ksw1d(nzcanm)               ! Shortwave view factor from sky to wall 
    real(r8), intent(out) :: kvg1d(nzcanm)               ! Shortwave view factor from vegetation to ground 
    real(r8), intent(out) :: ksr1d(nzcanm)               ! Shortwave view factor from sky to roof 
    real(r8), intent(out) :: ksv1d(nzcanm)               ! Shortwave view factor from sky to vegetation 
    real(r8), intent(out) :: kws1d(nzcanm)               ! Shortwave view factor from wall to sky
    real(r8), intent(out) :: kvs1d(nzcanm)               ! Shortwave view factor from vegetation to sky
    real(r8), intent(out) :: krs1d(nzcanm)               ! Shortwave view factor from roof to sky
    real(r8), intent(out) :: ksg1d                       ! Shortwave view factor from sky to ground 
    real(r8), intent(out) :: kts1d                       ! Shortwave view factor from ground to sky
    
    ! Unweighted view factors (copied for storage) - longwave (suffix _f) and shortwave (suffix _k)
    real(r8), intent(out) :: vfww_f(nzcanm,nzcanm)       ! Unweighted view factor from wall to wall (longwave)
    real(r8), intent(out) :: vfvv_f(nzcanm,nzcanm)       ! Unweighted view factor from vegetation to vegetation (longwave)
    real(r8), intent(out) :: vfwv_f(nzcanm,nzcanm)       ! Unweighted view factor from wall to vegetation (longwave)
    real(r8), intent(out) :: vfvw_f(nzcanm,nzcanm)       ! Unweighted view factor from vegetation to wall (longwave)
    real(r8), intent(out) :: vfwr_f(nzcanm,nzcanm)       ! Unweighted view factor from wall to roof (longwave)
    real(r8), intent(out) :: vfrw_f(nzcanm,nzcanm)       ! Unweighted view factor from roof to wall (longwave)
    real(r8), intent(out) :: vfvr_f(nzcanm,nzcanm)       ! Unweighted view factor from vegetation to roof (longwave)
    real(r8), intent(out) :: vfrv_f(nzcanm,nzcanm)       ! Unweighted view factor from roof to vegetation (longwave)
    real(r8), intent(out) :: vfwt_f(nzcanm)              ! Unweighted view factor from wall to ground (longwave)
    real(r8), intent(out) :: vftw_f(nzcanm)              ! Unweighted view factor from ground to wall (longwave)
    real(r8), intent(out) :: vftv_f(nzcanm)              ! Unweighted view factor from ground to vegetation (longwave)
    real(r8), intent(out) :: vfvt_f(nzcanm)              ! Unweighted view factor from vegetation to ground (longwave)
    real(r8), intent(out) :: vfsw_f(nzcanm)              ! Unweighted view factor from sky to wall (longwave)
    real(r8), intent(out) :: vfsr_f(nzcanm)              ! Unweighted view factor from sky to roof (longwave)
    real(r8), intent(out) :: vfsv_f(nzcanm)              ! Unweighted view factor from sky to vegetation (longwave)
    real(r8), intent(out) :: vfst_f                      ! Unweighted view factor: sky -> ground (longwave)
    real(r8), intent(out) :: svft_f                      ! Unweighted sky-view factor: ground -> sky (longwave)

    real(r8), intent(out) :: vfww_k(nzcanm,nzcanm)       ! Unweighted view factor from wall to wall (shortwave)
    real(r8), intent(out) :: vfvv_k(nzcanm,nzcanm)       ! Unweighted view factor from vegetation to vegetation (shortwave)
    real(r8), intent(out) :: vfwv_k(nzcanm,nzcanm)       ! Unweighted view factor from wall to vegetation (shortwave)
    real(r8), intent(out) :: vfvw_k(nzcanm,nzcanm)       ! Unweighted view factor from vegetation to wall (shortwave)
    real(r8), intent(out) :: vfwr_k(nzcanm,nzcanm)       ! Unweighted view factor from wall to roof (shortwave)
    real(r8), intent(out) :: vfrw_k(nzcanm,nzcanm)       ! Unweighted view factor from roof to wall (shortwave)
    real(r8), intent(out) :: vfvr_k(nzcanm,nzcanm)       ! Unweighted view factor from vegetation to roof (shortwave)
    real(r8), intent(out) :: vfrv_k(nzcanm,nzcanm)       ! Unweighted view factor from roof to vegetation (shortwave)
    real(r8), intent(out) :: vfwt_k(nzcanm)              ! Unweighted view factor from wall to ground (shortwave)
    real(r8), intent(out) :: vftw_k(nzcanm)              ! Unweighted view factor from ground to wall (shortwave)
    real(r8), intent(out) :: vftv_k(nzcanm)              ! Unweighted view factor from ground to vegetation (shortwave)
    real(r8), intent(out) :: vfvt_k(nzcanm)              ! Unweighted view factor from vegetation to ground (shortwave)
    real(r8), intent(out) :: vfsw_k(nzcanm)              ! Unweighted view factor from sky to wall (shortwave)
    real(r8), intent(out) :: vfsr_k(nzcanm)              ! Unweighted view factor from sky to roof (shortwave)
    real(r8), intent(out) :: vfsv_k(nzcanm)              ! Unweighted view factor from sky to vegetation (shortwave)
    real(r8), intent(out) :: vfst_k                      ! Unweighted view factor from sky to ground (shortwave)
    real(r8), intent(out) :: svft_k                      ! Unweighted sky-view factor: ground -> sky (shortwave)

    real(r8), intent(out) :: svfw_f(nzcanm)              ! Unweighted sky-view factor: wall -> sky (longwave)
    real(r8), intent(out) :: svfv_f(nzcanm)              ! Unweighted sky-view factor: vegetation -> sky (longwave)
    real(r8), intent(out) :: svfr_f(nzcanm)              ! Unweighted sky-view factor: roof -> sky (longwave)
    
    real(r8), intent(out) :: svfw_k(nzcanm)              ! Unweighted sky view factor for wall (shortwave)
    real(r8), intent(out) :: svfv_k(nzcanm)              ! Unweighted sky view factor for vegetation (shortwave)
    real(r8), intent(out) :: svfr_k(nzcanm)              ! Unweighted sky view factor for roof (shortwave)
    
    ! LOCAL VARIABLES:    
    real(r8), parameter :: pi = 3.1415926535897932384626433832795_r8
    real(r8)            :: lad_tree(nzcanm)              ! In-crown (crown-volume average) leaf area density [m2 m-3] = lad/wtroad_tree.
    integer             :: k, izcan, jzcan, kzcan        ! Indices
    real(r8)            :: A_s                           ! Plan area of the urban unit (canyon + building), normalized by dzcan;
                                                         ! this is also the area of the sky 'surface' (= xdom)
    real(r8)            :: A_g                           ! Plan area of the canyon floor / road, normalized by dzcan
    real(r8)            :: A_w(nzcanm)                   ! Area of the wall at each level, normalized by dzcan (= pb(izcan+1))
    real(r8)            :: A_w_max(nzcanm)               ! A_w floored at 1e-6 so it can safely be used as a divisor
    real(r8)            :: A_r(nzcanm)                   ! Area of the roof at each level, normalized by dzcan
    real(r8)            :: A_r_max(nzcanm)               ! A_r floored at 1e-6 so it can safely be used as a divisor
    real(r8)            :: A_vs(nzcanm)                  ! Active two-sided foliage area per layer for shortwave calcs, normalized by
                                                         ! dzcan and already multiplied by the clumping factor omega
    real(r8)            :: A_vs_max(nzcanm)              ! A_vs floored at 1e-6 so it can safely be used as a divisor
    real(r8)            :: A_vl(nzcanm)                  ! Active two-sided foliage area per layer for longwave calcs, normalized by
                                                         ! dzcan and already multiplied by the clumping factor omega
    real(r8)            :: A_vl_max(nzcanm)              ! A_vl floored at 1e-6 so it can safely be used as a divisor
    real(r8)            :: wtot                          ! Total domain width (canyon + building) [m]
    real(r8)            :: xdom                          ! Total domain width normalized by dzcan (= wtot/dzcan)
    real(r8)            :: bldfrac                       ! Fraction of building width in total domain width
    real(r8)            :: kbs_vf                        ! Foliage extinction coefficient used in the ray tracing [-]; 0.5 for a
                                                         ! spherical leaf-angle distribution with rays normal to the target surface
    real(r8)            :: rayy                          ! Ray vertical position, in units of dzcan (so layer izcan spans izcan-1 to izcan)
    real(r8)            :: rayx                          ! Ray horizontal position, in units of dzcan
    real(r8)            :: xfr                           ! rayx folded back into one canyon+building period, still in units of dzcan
    integer             :: kk,nk,ii                      ! Indices
    real(r8)            :: phi, phi1                     ! Azimuthal angle
    real(r8)            :: h                             ! Evenly spaced z-component of a point on the unit sphere (= cos(theta))
    real(r8)            :: theta                         ! Zenith angle calculated as acos(h)
    real(r8)            :: rnum                          ! A random number    
    real(r8)            :: x, y, z2                      ! Cartesian components of a point on the unit sphere (a ray direction).
    real(r8)            :: xxe(n), zze(n)                ! x- and z-components of the ray direction unit vectors, flattened onto the
                                                         ! x-z plane (y removed) and renormalized to unit length
    real(r8)            :: xxe2(n), zze2(n)              ! As xxe/zze, but for the ray set with preferential reflection over transmission
    real(r8)            :: xr, z2r                       ! x- and z-components of a hemispherically (Lambert cosine) distributed ray direction
    real(r8)            :: xxhr(n), zzhr(n)              ! x- and z-components of the hemispherically distributed ray direction unit
                                                         ! vectors; these are the rays released from surfaces (ground, roof, wall, sky)
    real(r8)            :: circ32e(n)                    ! Ratio of 3D ray travel distance to its 2D projection
    real(r8)            :: circ32e2(n)                   ! Ratio of 3D ray travel distance to its 2D projection for preferential reflection
    real(r8)            :: hemi32r(n)                    ! Ratio of 3D ray travel distance to its 2D projection   
    real(r8)            :: rzen                          ! A zenith angle generated using Lambert's cosine law
    real(r8)            :: raz                           ! A random azimuthal angle between 0 and 2*pi
    real(r8)            :: minray                        ! Ray strength below which the ray is no longer tracked
    
    real(r8)            :: vft_tot                      ! Total view factor from ground
    real(r8)            :: vfs_tot                      ! Total view factor from sky
    real(r8)            :: vfr_tot(nzcanm)              ! Total view factors from roof     
    real(r8)            :: vfw_tot(nzcanm)              ! Total view factors from wall
    real(r8)            :: vfv_tot(nzcanm)              ! Total view factors from vegetation 
    integer             :: nsrays                       ! Number of sky-ray start bins across the domain (= maxind*nsky); used only to
                                                        ! spread the sky-ray start positions rayx evenly over the domain
    
    real(r8)            :: svfw(nzcanm)                 ! Sky view factor of each wall level (ray-summed, then divided by the ray count)
    real(r8)            :: svfr(nzcanm)                 ! Sky view factor of each roof level (ray-summed, then divided by the ray count)
    real(r8)            :: svft                         ! Sky view factor of the ground (ray-summed, then divided by the ray count).
    real(r8)            :: svfv(nzcanm)                 ! Sky view factor of each foliage level (ray-summed, then divided by the ray count)
    
    real(r8)            :: vfr(nzcanm)                  ! Single-ray contribution to the view factor to each roof level, from ray_up/ray_dn
    real(r8)            :: vfw(nzcanm)                  ! Single-ray contribution to the view factor to each wall level, from ray_up/ray_dn
    real(r8)            :: vfv(nzcanm)                  ! Single-ray contribution to the view factor to each foliage level, from ray_up/ray_dn
    real(r8)            :: vft                          ! Single-ray contribution to the view factor to the ground, from ray_dn
    
    real(r8)            :: xxt                          ! x-component of the ray direction for rays released from the sky
    real(r8)            :: circ32es                     ! Ratio of 3D ray travel distance to its 2D projection in vegetation vf calculation
    real(r8)            :: zzes                         ! z-direction component in vegetation vf calculation
    real(r8)            :: xxrt                         ! x-component of the ray direction for rays released from roof/road/foliage
    
    real(r8)            :: raystr                       ! Ray strength
    real(r8)            :: bld                          ! 1 if the ray is currently over a building column, 0 if over the canyon
    real(r8)            :: bld_old                      ! Value of bld at the previous ray step (used to detect wall/roof crossings)
    real(r8)            :: pb_old                       ! Pb from the previous step
    
    real(r8)            :: vfww(nzcanm, nzcanm)         ! View factor from wall to wall
    real(r8)            :: vfwr(nzcanm, nzcanm)         ! View factor from wall to roof
    real(r8)            :: vfwv(nzcanm, nzcanm)         ! View factor from wall to vegetation
    real(r8)            :: vfrw(nzcanm, nzcanm)         ! View factor from roof to wall
    real(r8)            :: vfrv(nzcanm, nzcanm)         ! View factor from roof to vegetation
    real(r8)            :: vfvw(nzcanm, nzcanm)         ! View factor from vegetation to wall
    real(r8)            :: vfvr(nzcanm, nzcanm)         ! View factor from vegetation to roof
    real(r8)            :: vfvv(nzcanm, nzcanm)         ! View factor from vegetation to vegetation
    real(r8)            :: vfsr(nzcanm)                 ! View factors from sky to roof
    real(r8)            :: vfsw(nzcanm)                 ! View factors from sky to wall
    real(r8)            :: vfsv(nzcanm)                 ! View factors from sky to vegetation

    real(r8)            :: vfswtmp(nzcanm, nzcanm)      ! Temporary view factors from sky to wall
    real(r8)            :: vfsrtmp(nzcanm, nzcanm)      ! Temporary view factors from sky to roof
    real(r8)            :: vfsvtmp(nzcanm, nzcanm)      ! Temporary view factors from sky to vegetation
    real(r8)            :: vfsttmp(nzcanm)              ! Per-start-bin contribution to vfst from sky rays
    real(r8)            :: vfwt(nzcanm)                 ! View factor from wall to ground
    real(r8)            :: vftw(nzcanm)                 ! View factor from ground to wall
    real(r8)            :: vftv(nzcanm)                 ! View factor from ground to vegetation
    real(r8)            :: vfvt(nzcanm)                 ! View factor from vegetation to ground
    real(r8)            :: vfst                         ! View factor from sky to ground, accumulated over the sky-ray start bins

    logical             :: horiz                        ! .true. if the emitting surface is horizontal (ground, roof, sky);
                                                        ! .false. for a vertical wall
    logical             :: solar                        ! .true. while computing the shortwave (solar) pass, .false. for longwave
    integer             :: start_time, end_time, clock_rate    ! Timekeeping variables
    real(r8)            :: elapsed_time                        ! Elapsed time 
    logical             :: foliage                      ! .true. if the ray originates from a foliage (tree crown) layer
    logical  :: debug_write                  ! true => write out many intermediate variables for debugging
    !----------------------------------------------------------------
    ! Get the clock rate (ticks per second)
    call system_clock(count_rate=clock_rate)
    call system_clock(start_time)
                                    
    debug_write = .false.!.true.
    
    !Initialization
      fww1d=0._r8; fvv1d=0._r8; fwv1d=0._r8; fvw1d=0._r8
      fwr1d=0._r8; frw1d=0._r8; fvr1d=0._r8; frv1d=0._r8
      fwg1d=0._r8; fgw1d=0._r8; fgv1d=0._r8; fsw1d=0._r8
      fvg1d=0._r8; fsr1d=0._r8; fsv1d=0._r8
      fws1d=0._r8; fvs1d=0._r8; frs1d=0._r8
      fsg1d=0._r8; fts1d=0._r8

      kww1d=0._r8; kvv1d=0._r8; kwv1d=0._r8; kvw1d=0._r8
      kwr1d=0._r8; krw1d=0._r8; kvr1d=0._r8; krv1d=0._r8
      kwg1d=0._r8; kgw1d=0._r8; kgv1d=0._r8; ksw1d=0._r8
      kvg1d=0._r8; ksr1d=0._r8; ksv1d=0._r8
      kws1d=0._r8; kvs1d=0._r8; krs1d=0._r8
      ksg1d=0._r8; kts1d=0._r8

      vfww_f=0._r8; vfvv_f=0._r8; vfwv_f=0._r8; vfvw_f=0._r8
      vfwr_f=0._r8; vfrw_f=0._r8; vfvr_f=0._r8; vfrv_f=0._r8
      vfwt_f=0._r8; vftw_f=0._r8; vftv_f=0._r8; vfvt_f=0._r8
      vfsw_f=0._r8; vfsr_f=0._r8; vfsv_f=0._r8
      vfst_f=0._r8; svft_f=0._r8
      svfw_f=0._r8; svfv_f=0._r8; svfr_f=0._r8

      vfww_k=0._r8; vfvv_k=0._r8; vfwv_k=0._r8; vfvw_k=0._r8
      vfwr_k=0._r8; vfrw_k=0._r8; vfvr_k=0._r8; vfrv_k=0._r8
      vfwt_k=0._r8; vftw_k=0._r8; vftv_k=0._r8; vfvt_k=0._r8
      vfsw_k=0._r8; vfsr_k=0._r8; vfsv_k=0._r8
      vfst_k=0._r8; svft_k=0._r8
      svfw_k=0._r8; svfv_k=0._r8; svfr_k=0._r8
    !----------------------------------------------------------------
    ! Foliage extinction coefficient for the ray tracing. Krayenhoff et al. (2014) Sect. 3.3.2 
    ! K_bs = 0.5 is the spherical (random, uniformly distributed) leaf-angle distribution value, and is what
    kbs_vf=1.0_r8/2.0_r8
    
    ! no-tree guard 
    if (lai > 1e-6_r8) then
      lad_tree=lad/wtroad_tree
    else
      lad_tree=0.0_r8
    endif
    
    ! Fixed seed. A deterministic seed is required here so that the
    ! view factors are reproducible run-to-run and bit-for-bit identical on restart.
    call init_random_seed(1234)
    call RANDOM_NUMBER(rnum)
    
    !----------------------------------------------------------------------------
    ! Generate evenly-spaced points on a sphere of radius 1 
    ! Get xxe, zze, circ32e -> the FULL-SPHERE direction set used by rays leaving a foliage layer
    ! (Krayenhoff et al. 2014 Sect. 3.3.1: 'Rays leaving a layer of tree foliage are assumed isotropic and
    !  therefore choose from a full sphere (4 pi solid angle)').
    !----------------------------------------------------------------------------
    ! (This code is originally written by Joseph O'Rourke and Min Xu, June 1997, and was converted to Fortran from C++.)

    ! kk counts the hemispherical (surface) directions accumulated below; phi1 carries the running spiral azimuth
    kk=1
    phi1=0.0_r8
    ! The loop deliberately runs k=2..n-1, skipping the two poles (h=-1 and h=+1) where the spiral azimuth
    ! increment below would divide by sqrt(1-h*h)=0. That is why exactly n-2 foliage directions are produced.
  	do k=2,n-1    
       h=-1.0_r8+2.0_r8*real(k-1,r8)/real(n-1,r8)
       theta=acos(h)

       ! Spiral azimuth step (Rakhmanov-Saff-Zhou)
       phi=phi1+3.6_r8/(sqrt(real(n,r8)*(1._r8-h*h))) 

       ! Wrap into [0,2*pi). Same as modulo(phi,2._r8*pi); 
       phi = phi - (2._r8*pi) * floor(phi / (2._r8*pi))
       phi1=phi
       
       !Convert spherical coordinates (theta, phi) to Cartesian coordinates (x, y, z2)
       x=cos(phi)*sin(theta)
      !  y=sin(phi)*sin(theta)
       z2=h
       
       !----------------------------------------------------------------------------
       ! Flattened into 2-D (remove y-dimension) and increase all 2-D vector lengths to lie on the unit circle
       xxe(k-1)=x/sqrt(x*x+z2*z2)
       zze(k-1)=z2/sqrt(x*x+z2*z2)
       ! 3-D to 2-D ratio of ray travel distance
       circ32e(k-1)=1._r8/sqrt(x*x+z2*z2)		

       !! ADDING IN preferential reflection over transmission (ratio 3:2) for
       ! solar radiation (due to preference of NIR to reflect rather than
       ! transmit through leaves, and direct solar to arrive from above).
       call RANDOM_NUMBER(rnum)
       
       !----------------------------------------------------------------------------
       ! Generate evenly-spaced points on a sphere of radius 1 with preferential refelection
       ! Get xxe2, zze2, circ32e2
       !----------------------------------------------------------------------------
       if (z2 < 0._r8 .and. rnum > 0.8_r8) then
          z2=abs(z2)
       endif
       
       ! Flattened into 2-D (remove y-dimension) and increase all 2-D vector lengths to lie on the unit circle
       xxe2(k-1)=x/sqrt(x*x+z2*z2)
       zze2(k-1)=z2/sqrt(x*x+z2*z2)
       
       ! 3-D to 2-D ratio of ray travel distance
       circ32e2(k-1)=1._r8/sqrt(x*x+z2*z2)
       
       !----------------------------------------------------------------------------
       ! Generate random points with hemispherical distribution
       ! Get xxhr, zzhr, hemi32r
       !----------------------------------------------------------------------------
       ! hemispherical distribution (for rays leaving solid surfaces: ground, roof, wall, sky)
       ! Gate on h, not z2: z2 may have been flipped just above, which would give nk ~ 0.6*(n-2) instead of
       ! (n-2)/2 and break the 2:1 foliage:surface ray ratio required by Krayenhoff et al. (2014) Sect. 3.3.1.
       if (h > 0._r8) then
          ! random hemispherical angles with cosine probability density (Lambert's cosine law)
          call RANDOM_NUMBER(rnum)
          rzen=asin(sqrt(rnum))      ! Actual Lambert cosine law (Siegel and Howell 2002)...
                                     ! equivalent to acos(sqrt(1-rnum)) (Kondo et al. 2001) or...
                                     ! 0.5*acos(1-2*rnum) (Chelle 2006)
!189   continue                                     
          call RANDOM_NUMBER(rnum)
          raz=rnum*2._r8*pi
          xr=sin(rzen)*cos(raz)
          z2r=cos(rzen)

          xxhr(kk)=xr/sqrt(xr*xr+z2r*z2r) 
          zzhr(kk)=z2r/sqrt(xr*xr+z2r*z2r)
          hemi32r(kk)=1._r8/sqrt(xr*xr+z2r*z2r)
                   
          kk=kk+1
       endif
    enddo
    ! Number of hemispherical (surface) ray directions: exactly (n-2)/2, since the gate above tests h
    nk=kk-1
    
    !-----------------------------------------------------------------------
    ! ray strength at which we stop tracking the ray:
    minray=max(0.000001_r8,0.0001_r8*pb(maxbhind))
    
    call system_clock(end_time)        
    ! Calculate elapsed time in seconds
    elapsed_time = real(end_time - start_time,r8) / real(clock_rate,r8)
    
   !  write the elapsed time for each iteration  
   !  write(*, '(A, I0, A, F6.3)') 'Preparing stage elapsed time for l = ', l, ': ', elapsed_time, ' seconds'
        
    !-----------------------------------------------------------------------
    !! VIEW FACTOR CALCULATIONS
    !! Calculate view factors for sfc-sfc, veg-sfc, and veg-veg diffuse exchange using ray tracing
    !! Only need to find view factors for wall on one side due to symmetry
    !-----------------------------------------------------------------------
    ! Width of the two columns (canyon and building)
    wtot=wcan+wbui
    bldfrac=wbui/wtot
    xdom=wtot/dzcan 
    nsrays=(maxind)*nsky 
     
    solar=.false.
    goto 359

358   continue

    solar=.true.
    
    !write(6,*)'----------SOLAR VIEW FACTORS----------'

359   continue

   !  if (.not.solar) then
       !write(6,*)'----------LONGWAVE VIEW FACTORS----------'
   !  endif
    
    ! total view factors (diagnostics to see if they add up to 1.0 for each surface)
    vft_tot=0._r8
    vfs_tot=0._r8
    
    do izcan=1,nzcanm
       vfr_tot(izcan)=0._r8
       vfw_tot(izcan)=0._r8
       vfv_tot(izcan)=0._r8
    enddo
        
    ! MAIN IZCAN LOOP:
    do izcan=1,maxind
       vfst=0._r8
       svft=0._r8
       svfw(izcan)=0._r8
       svfr(izcan)=0._r8
       svfv(izcan)=0._r8

       ! consider run these reset loops to nzcanm instead of maxind.
       do jzcan=1,maxind
         ! w=walls; r=roofs; v=vegetation (trees); t=ground; s=sky
          vfww(izcan,jzcan)=0._r8
          vfwr(izcan,jzcan)=0._r8
          vfwv(izcan,jzcan)=0._r8
          vfrw(izcan,jzcan)=0._r8
          vfrv(izcan,jzcan)=0._r8
          vfvw(izcan,jzcan)=0._r8
          vfvr(izcan,jzcan)=0._r8
          vfvv(izcan,jzcan)=0._r8
          vfsr(jzcan)=0._r8
          vfsw(jzcan)=0._r8
          vfsv(jzcan)=0._r8

          vfswtmp(izcan,jzcan)=0._r8
          vfsrtmp(izcan,jzcan)=0._r8
          vfsvtmp(izcan,jzcan)=0._r8

       enddo
       vfwt(izcan)=0._r8
       vftw(izcan)=0._r8
       vftv(izcan)=0._r8
       vfvt(izcan)=0._r8

       !-----------------------------------------------------
       !write(6,*)'SKY...'
       ! rays starting from SKY.
       horiz=.true.
       do ii=1,nsky
          do kk=1,nk
             raystr=1._r8
             call RANDOM_NUMBER(rnum)
             ! Rays leaving the 'sky' begin at a random horizontal location in equal bins
             rayx=(real((izcan-1)*nsky+ii,r8)-rnum)/real(nsrays,r8)*xdom

             ! Ray's starting vertical position
             rayy=max(real(maxind,r8)*hsky,real(maxind,r8)+2._r8*dray)

             ! make all rays head 'downstream' (this works due to the symmetry of the geometry)
             xxt=xxhr(kk) ! use random hemispherical distribution
             ! mirror reflection so that all rays still head in positive direction
             if (xxt < 0._r8) then
                xxt=-xxt
                rayx=xdom*(2._r8-bldfrac)-rayx
             endif
             
             do while (rayy > real(maxind,r8)+2._r8*dray)
                ! COULD MULTIPLY dray BY A FACTOR PROPORTIONAL TO THE DISTANCE ABOVE REAL(MAXIND), AND WITH A MINIMUM VALUE, TO SPEED THIS UP!
                rayx=rayx+xxt*dray
                rayy=rayy-zzhr(kk)*dray
             enddo
             
             ! test where a ray is horizontally (i.e. building or canyon/tree space):
             ! xfr>bldfrac means canyon, xfr<bldfrac means building (i.e. bld=1)
             xfr = rayx - xdom* floor(rayx /xdom)
             bld=0._r8
             if (xfr/xdom >= 1._r8-bldfrac) bld=1._r8
             bld_old=bld
             pb_old=0._r8
             call ray_dn(minray,rayx,rayy,bld,bld_old,raystr, &
                      xxt,-zzhr(kk),hemi32r(kk),xdom,dray,dzcan,bldfrac,pb, &
                      pb_old,ss,kbs_vf,lad,omega,horiz,h1,h2, vfw,vfr,vfv,vft)

             do kzcan=1,nzcanm
                vfsw(kzcan)=vfsw(kzcan)+vfw(kzcan)
                vfsr(kzcan)=vfsr(kzcan)+vfr(kzcan)
                vfsv(kzcan)=vfsv(kzcan)+vfv(kzcan)
             enddo
             vfst=vfst+vft

          enddo  ! end rays (jj) loop for SKY
       enddo  ! end additional loop (ii) over 'sky locations'
       

       do jzcan=1,nzcanm
          vfswtmp(izcan,jzcan)=vfsw(jzcan)/real(nk,r8)/real(nsky,r8)
          vfsrtmp(izcan,jzcan)=vfsr(jzcan)/real(nk,r8)/real(nsky,r8)
          vfsvtmp(izcan,jzcan)=vfsv(jzcan)/real(nk,r8)/real(nsky,r8)
       enddo

       vfsttmp(izcan)=vfst/real(nk,r8)/real(nsky,r8)
       !-----------------------------------------------------
       ! rays starting from WALLS
       if (izcan <= maxbhind-1) then
          !write(6,*)'WALLS...'
          horiz=.false.
          do kk=1,nk
             ! starting at the downstream building edge (upstream edge of the canyon)
             bld=0._r8
             bld_old=0._r8
             raystr=1._r8
             rayx=0._r8
             ! random:
             call RANDOM_NUMBER(rnum)
             ! ray's starting vertical position is random
             rayy=real(izcan,r8)-rnum
             
             ! use random hemispherical distribution
             if (xxhr(kk) <= 0._r8) goto 545
             ! rays going up
             foliage=.false.
             call ray_up(minray,maxind,rayx,rayy,bld,bld_old,raystr, &
                        xxhr(kk),zzhr(kk),hemi32r(kk),xdom,dray,dzcan,&
                        bldfrac,pb,kbs_vf,lad_tree,lad,omega,horiz,vfw,vfv,h1,h2,foliage)

             do kzcan=1,nzcanm
                vfww(izcan,kzcan)=vfww(izcan,kzcan)+vfw(kzcan)
                vfwv(izcan,kzcan)=vfwv(izcan,kzcan)+vfv(kzcan)
             enddo
             goto 547
545     continue
             ! rays going down
             pb_old=pb(izcan+1)
             call ray_dn(minray,rayx,rayy,bld,bld_old,raystr,&
                        xxhr(kk),zzhr(kk),hemi32r(kk),xdom,dray,dzcan,&
                        bldfrac,pb,pb_old,ss,kbs_vf,lad,omega,horiz,h1,h2, vfw,vfr,vfv,vft)

             do kzcan=1,nzcanm
                vfww(izcan,kzcan)=vfww(izcan,kzcan)+vfw(kzcan)
                vfwr(izcan,kzcan)=vfwr(izcan,kzcan)+vfr(kzcan)
                vfwv(izcan,kzcan)=vfwv(izcan,kzcan)+vfv(kzcan)
             enddo
             vfwt(izcan)=vfwt(izcan)+vft
             goto 548
547    continue
            !! WALLS-SKY
            svfw(izcan)=svfw(izcan)+raystr
548    continue
          enddo  ! end rays (kk) loop for WALLS
          
          svfw(izcan)=svfw(izcan)/real(nk,r8)
          vfwt(izcan)=vfwt(izcan)/real(nk,r8)
          vfw_tot(izcan)=vfw_tot(izcan)+svfw(izcan)+vfwt(izcan)
       endif  ! if izcan has a wall layer

       !-----------------------------------------------------
       ! rays starting from VEGETATION in CANOPY COLUMN
       if (lad(izcan) > 0._r8) then
           !write(6,*)'CANOPY VEGETATION...'
           ! horiz=.true. here does NOT mean 'the emitter is a horizontal surface'. 
           ! it just means that the ray direction is not rotated 
           horiz=.true.
           do k=1,n-2
              ! starting at a random location in the vegetation
              bld=0._r8
              bld_old=0._r8
              raystr=1._r8
              call RANDOM_NUMBER(rnum)
              rayx=rnum*xdom*(1._r8-bldfrac)
              
              call RANDOM_NUMBER(rnum)

              if ((h1+h2)<=dzcan) then 
               ! when (h1+h2)<=dzcan, only izcan=1 is computed
                  rayy=(h1 + h2*rnum)/dzcan
              else if ((h1+h2)>dzcan) then 
                  if (izcan==1) then
                     rayy=(h1 + (dzcan-h1)*rnum)/dzcan
                  else if (izcan==2) then
                     rayy=1.0_r8 + (h1 +h2-dzcan)*rnum/dzcan
                  end if
              end if
              
              if (solar) then
                 xxrt=xxe2(k)
                 zzes=zze2(k)
                 circ32es=circ32e2(k)
              else
                 xxrt=xxe(k)
                 zzes=zze(k)
                 circ32es=circ32e(k)
              endif
              
              ! mirror reflection so that all rays still head in positive direction
              if (xxrt < 0._r8) then
                 xxrt=-xxrt
                 rayx=xdom*(1._r8-bldfrac)-rayx
              endif
	            
              if(zzes <= 0._r8) goto 645
              ! rays going up
              
              foliage=.true.
              call ray_up(minray,maxind,rayx,rayy,bld,bld_old,raystr, &
                         xxrt,zzes,circ32es,xdom,dray,dzcan,bldfrac,pb,&
                         kbs_vf,lad_tree,lad,omega,horiz,vfw,vfv,h1,h2,foliage)    
               
              do kzcan=1,nzcanm
                 vfvw(izcan,kzcan)=vfvw(izcan,kzcan)+vfw(kzcan)
                 vfvv(izcan,kzcan)=vfvv(izcan,kzcan)+vfv(kzcan)
              enddo
              goto 647

645     continue
              ! rays going down
              pb_old=pb(izcan+1)
              call ray_dn(minray,rayx,rayy,bld,bld_old,raystr,&
                        xxrt,zzes,circ32es,xdom,dray,dzcan,bldfrac, &
                        pb,pb_old,ss,kbs_vf,lad,omega,horiz,h1,h2, vfw,vfr,vfv,vft)
                        
              do kzcan=1,nzcanm
                 vfvw(izcan,kzcan)=vfvw(izcan,kzcan)+vfw(kzcan)
                 vfvr(izcan,kzcan)=vfvr(izcan,kzcan)+vfr(kzcan)
                 vfvv(izcan,kzcan)=vfvv(izcan,kzcan)+vfv(kzcan)
              enddo
              vfvt(izcan)=vfvt(izcan)+vft

              goto 648
647    continue
              !! VEGETATION-SKY
              svfv(izcan)=svfv(izcan)+raystr
648    continue
           enddo  ! end rays (k) loop for VEGETATION

           ! n-2 is the number of foliage rays actually fired (the k loop above). 
           svfv(izcan)=svfv(izcan)/real(n-2,r8)
           vfvt(izcan)=vfvt(izcan)/real(n-2,r8)
           vfv_tot(izcan)=vfv_tot(izcan)+vfvt(izcan)+svfv(izcan)
       endif  ! if izcan has a foliage layer

       !-----------------------------------------------------
       ! rays starting from ROOFS
       if (ss(izcan) > 0._r8) then
          !write(6,*)'ROOFS...'
          horiz=.true.
          do kk=1,nk
             ! starting at a random point on the roof
             bld=1._r8
             bld_old=1._r8
             raystr=1._r8
             rayy=real(izcan-1,r8)
             ! start at a random point on the roof
             call RANDOM_NUMBER(rnum)
             rayx=xdom*((1._r8-bldfrac)+rnum*bldfrac)
             
             ! use random hemispherical distribution
             ! make all rays head 'downstream' (this works due to the symmetry of the geometry)
             xxrt=xxhr(kk)
             ! mirror reflection so that all rays still head in positive direction
             if (xxrt < 0._r8) then
                xxrt=-xxrt
                rayx=xdom-(rayx-xdom*(1._r8-bldfrac))
             endif
             foliage=.false.
             call ray_up(minray,maxind,rayx,rayy,bld,bld_old,raystr,&
                       xxrt,zzhr(kk),hemi32r(kk),xdom,dray,dzcan,&
                       bldfrac,pb,kbs_vf,lad_tree,lad,omega,horiz,vfw,vfv,h1,h2,foliage)

             do kzcan=1,nzcanm
                vfrw(izcan,kzcan)=vfrw(izcan,kzcan)+vfw(kzcan)
                vfrv(izcan,kzcan)=vfrv(izcan,kzcan)+vfv(kzcan)
             enddo
             svfr(izcan)=svfr(izcan)+raystr

          enddo  ! end rays (kk) loop for ROOFS
          svfr(izcan)=svfr(izcan)/real(nk,r8)
          vfr_tot(izcan)=vfr_tot(izcan)+svfr(izcan)
       endif

    !-----------------------------------------------------
    enddo  ! end wall/vegetation/roof level (izcan) loop
    !write(6,*)'AFTER IZCAN LOOP'
    !-----------------------------------------------------------------------
    !write(6,*)'ROAD...'
    ! rays starting from ROAD
    horiz=.true.
    do kk=1,nk
       ! starting at a random point on the road
       bld=0._r8
       bld_old=0._r8
       raystr=1._r8
       
       rayy=0._r8
       ! start at a random point on the road
       call RANDOM_NUMBER(rnum)
       rayx=xdom*(1._r8-bldfrac)*rnum
       ! use random hemispherical distribution
       ! make all rays head 'downstream' (this works due to the symmetry of the geometry)
       xxrt=xxhr(kk)
       ! mirror reflection so that all rays still head in positive direction
       if (xxrt < 0._r8) then
          xxrt=-xxrt
          rayx=xdom*(1._r8-bldfrac)-rayx
       endif
       
       foliage=.false.
       call ray_up(minray,maxind,rayx,rayy,bld,bld_old,raystr,&
                 xxrt,zzhr(kk),hemi32r(kk),xdom,dray,dzcan,bldfrac,&
                 pb,kbs_vf,lad_tree,lad,omega,horiz,vfw,vfv,h1,h2,foliage)
   
       do kzcan=1,nzcanm
          vftw(kzcan)=vftw(kzcan)+vfw(kzcan)
          vftv(kzcan)=vftv(kzcan)+vfv(kzcan)
       enddo
       svft=svft+raystr
    enddo  ! end rays (kk) loop for ROAD

    svft=svft/real(nk,r8)

    do izcan=1,nzcanm
       !divide by 2 here because there are two walls at each level, and we want the vf from the road to only one of them

      ! has vf from other surfaces to wall divided by 2 already?
       vftw(izcan)=vftw(izcan)/real(nk,r8)/2._r8
       vftv(izcan)=vftv(izcan)/real(nk,r8)
    enddo   
    !-----------------------------------------------------------------------
    vfst=0._r8
    do izcan=1,maxind
       vfsw(izcan)=0._r8
       vfsr(izcan)=0._r8
       vfsv(izcan)=0._r8
       vfst=vfst+vfsttmp(izcan)
       do jzcan=1,maxind
          vfww(izcan,jzcan)=vfww(izcan,jzcan)/real(nk,r8)
          vfwv(izcan,jzcan)=vfwv(izcan,jzcan)/real(nk,r8)
          vfwr(izcan,jzcan)=vfwr(izcan,jzcan)/real(nk,r8)
          vfw_tot(izcan)=vfw_tot(izcan)+vfww(izcan,jzcan)+ &
                         vfwr(izcan,jzcan)+vfwv(izcan,jzcan)
          ! divide by 2 here because there are two walls at each level, and we want the vf from the roof to only one of them
          vfrw(izcan,jzcan)=vfrw(izcan,jzcan)/real(nk,r8)/2._r8
          vfrv(izcan,jzcan)=vfrv(izcan,jzcan)/real(nk,r8)
          vfr_tot(izcan)=vfr_tot(izcan)+vfrw(izcan,jzcan)*2._r8+vfrv(izcan,jzcan)
         
          vfvr(izcan,jzcan)=vfvr(izcan,jzcan)/real(n-2,r8)
          vfvw(izcan,jzcan)=vfvw(izcan,jzcan)/real(n-2,r8)/2._r8
          vfvv(izcan,jzcan)=vfvv(izcan,jzcan)/real(n-2,r8)          	     
          vfv_tot(izcan)=vfv_tot(izcan)+2._r8*vfvw(izcan,jzcan)+&
                         vfvr(izcan,jzcan)+vfvv(izcan,jzcan)         
          ! the index swap is intended               
          vfsw(izcan)=vfsw(izcan)+vfswtmp(jzcan,izcan)
          vfsr(izcan)=vfsr(izcan)+vfsrtmp(jzcan,izcan)
          vfsv(izcan)=vfsv(izcan)+vfsvtmp(jzcan,izcan)
       enddo
    enddo  
    
    do izcan=1,maxind
      ! they need to be divided by maxind
       vfsw(izcan)=vfsw(izcan)/real(maxind,r8)/2._r8
       vfsr(izcan)=vfsr(izcan)/real(maxind,r8)
       vfsv(izcan)=vfsv(izcan)/real(maxind,r8)
       vfs_tot=vfs_tot+2._r8*vfsw(izcan)+vfsr(izcan)+vfsv(izcan)
    enddo
    ! they need to be divided by maxind
    vfst=vfst/real(maxind,r8)
    vfs_tot=vfs_tot+vfst
 
    vft_tot=svft
    do jzcan=1,maxind	             
       vft_tot=vft_tot+2._r8*vftw(jzcan)+vftv(jzcan)           
    enddo
    
    !Calculate the modified radiation and the streets fluxes 
    !-------------------------------------------------------
    ! areas of all elements (normalized by dzcan, i.e., xdom = domain width/dzcan)
    do izcan=1,maxind
       A_w(izcan)=pb(izcan+1)
       ! two-sided vegetation area 
       if ((h1+h2)<= dzcan) then
          if (izcan==1) then
             A_vs(izcan)=xdom*(1._r8-bldfrac)*lads(izcan)*omega(izcan)*h2*2._r8
             A_vl(izcan)=xdom*(1._r8-bldfrac)*ladl(izcan)*omega(izcan)*h2*2._r8
          else if (izcan==2) then
             A_vs(izcan)=0._r8
             A_vl(izcan)=0._r8
          end if 
       else if ((h1+h2) > dzcan) then
           if (izcan==1) then
              A_vs(izcan)=xdom*(1._r8-bldfrac)*lads(izcan)*omega(izcan)*(dzcan-h1)*2._r8
              A_vl(izcan)=xdom*(1._r8-bldfrac)*ladl(izcan)*omega(izcan)*(dzcan-h1)*2._r8
           else if (izcan==2) then
              A_vs(izcan)=xdom*(1._r8-bldfrac)*lads(izcan)*omega(izcan)*(h1+h2-dzcan)*2._r8
              A_vl(izcan)=xdom*(1._r8-bldfrac)*ladl(izcan)*omega(izcan)*(h1+h2-dzcan)*2._r8
           end if 
       end if
       A_r(izcan)=xdom*bldfrac*ss(izcan)
      
       A_w_max(izcan)=max(1.e-6_r8,A_w(izcan))
       A_vs_max(izcan)=max(1.e-6_r8,A_vs(izcan))
       A_vl_max(izcan)=max(1.e-6_r8,A_vl(izcan))
       A_r_max(izcan)=max(1.e-6_r8,A_r(izcan))
    enddo
   
    A_s=xdom
    A_g=xdom*(1._r8-bldfrac)

    ! calculation of the modified radiation
    if (solar) then
    ! Solar view factors!
    ! view factors multiplied by relative areas (so that correct flux densities are exchanged):  
       do izcan=1,maxind
          do jzcan=1,maxind
             kww1d(izcan,jzcan)=vfww(izcan,jzcan)*A_w(izcan)/A_w_max(jzcan)
             kvv1d(izcan,jzcan)=vfvv(izcan,jzcan)*A_vs(izcan)/A_vs_max(jzcan)
             kwv1d(izcan,jzcan)=vfwv(izcan,jzcan)*A_w(izcan)/A_vs_max(jzcan)
             kvw1d(izcan,jzcan)=vfvw(izcan,jzcan)*A_vs(izcan)/A_w_max(jzcan)
             kwr1d(izcan,jzcan)=vfwr(izcan,jzcan)*A_w(izcan)/A_r_max(jzcan)
             krw1d(izcan,jzcan)=vfrw(izcan,jzcan)*A_r(izcan)/A_w_max(jzcan)
             kvr1d(izcan,jzcan)=vfvr(izcan,jzcan)*A_vs(izcan)/A_r_max(jzcan)
             krv1d(izcan,jzcan)=vfrv(izcan,jzcan)*A_r(izcan)/A_vs_max(jzcan)
             
             vfww_k(izcan,jzcan) = vfww(izcan,jzcan)
             vfvv_k(izcan,jzcan) = vfvv(izcan,jzcan)
             vfwv_k(izcan,jzcan) = vfwv(izcan,jzcan)
             vfvw_k(izcan,jzcan) = vfvw(izcan,jzcan)
             vfwr_k(izcan,jzcan) = vfwr(izcan,jzcan)
             vfrw_k(izcan,jzcan) = vfrw(izcan,jzcan)
             vfvr_k(izcan,jzcan) = vfvr(izcan,jzcan)
             vfrv_k(izcan,jzcan) = vfrv(izcan,jzcan)

          end do !izcan
       end do !jzcan

       do izcan=1,maxind                             
          kwg1d(izcan)=vfwt(izcan)*A_w(izcan)/A_g
          kgw1d(izcan)=vftw(izcan)*A_g/A_w_max(izcan)
          kgv1d(izcan)=vftv(izcan)*A_g/A_vs_max(izcan)
          kvg1d(izcan)=vfvt(izcan)*A_vs(izcan)/A_g

          ksw1d(izcan)=vfsw(izcan)*A_s/A_w_max(izcan)
          ksr1d(izcan)=vfsr(izcan)*A_s/A_r_max(izcan)
          ksv1d(izcan)=vfsv(izcan)*A_s/A_vs_max(izcan)
          kws1d(izcan)=svfw(izcan)
          krs1d(izcan)=svfr(izcan)
          kvs1d(izcan)=svfv(izcan)
          
          vfwt_k(izcan) = vfwt(izcan)
          vftw_k(izcan) = vftw(izcan)
          vftv_k(izcan) = vftv(izcan)
          vfvt_k(izcan) = vfvt(izcan)
          vfsw_k(izcan) = vfsw(izcan)
          vfsr_k(izcan) = vfsr(izcan)
          vfsv_k(izcan) = vfsv(izcan)
          svfw_k(izcan) = svfw(izcan)
          svfr_k(izcan) = svfr(izcan)
          svfv_k(izcan) = svfv(izcan)
          
       end do
       ksg1d=vfst*A_s/A_g
       kts1d=svft
       
       vfst_k=vfst
       svft_k=svft
    else
      ! Longwave view factors!
      ! view factors multiplied by relative areas (so that correct flux densities are exchanged):  
       do izcan=1,maxind
          do jzcan=1,maxind
             fww1d(izcan,jzcan)=vfww(izcan,jzcan)*A_w(izcan)/A_w_max(jzcan)
             fvv1d(izcan,jzcan)=vfvv(izcan,jzcan)*A_vl(izcan)/A_vl_max(jzcan)
             fwv1d(izcan,jzcan)=vfwv(izcan,jzcan)*A_w(izcan)/A_vl_max(jzcan)
             fvw1d(izcan,jzcan)=vfvw(izcan,jzcan)*A_vl(izcan)/A_w_max(jzcan)
             fwr1d(izcan,jzcan)=vfwr(izcan,jzcan)*A_w(izcan)/A_r_max(jzcan)
             frw1d(izcan,jzcan)=vfrw(izcan,jzcan)*A_r(izcan)/A_w_max(jzcan)
             fvr1d(izcan,jzcan)=vfvr(izcan,jzcan)*A_vl(izcan)/A_r_max(jzcan)
             frv1d(izcan,jzcan)=vfrv(izcan,jzcan)*A_r(izcan)/A_vl_max(jzcan)
             
             vfww_f(izcan,jzcan) = vfww(izcan,jzcan)
             vfvv_f(izcan,jzcan) = vfvv(izcan,jzcan)
             vfwv_f(izcan,jzcan) = vfwv(izcan,jzcan)
             vfvw_f(izcan,jzcan) = vfvw(izcan,jzcan)
             vfwr_f(izcan,jzcan) = vfwr(izcan,jzcan)
             vfrw_f(izcan,jzcan) = vfrw(izcan,jzcan)
             vfvr_f(izcan,jzcan) = vfvr(izcan,jzcan)
             vfrv_f(izcan,jzcan) = vfrv(izcan,jzcan)
             
          end do !izcan
       end do !jzcan
      
       do izcan=1,maxind                             
          fwg1d(izcan)=vfwt(izcan)*A_w(izcan)/A_g
          fgw1d(izcan)=vftw(izcan)*A_g/A_w_max(izcan)
          fgv1d(izcan)=vftv(izcan)*A_g/A_vl_max(izcan)
          fvg1d(izcan)=vfvt(izcan)*A_vl(izcan)/A_g

          fsw1d(izcan)=vfsw(izcan)*A_s/A_w_max(izcan)
          fsr1d(izcan)=vfsr(izcan)*A_s/A_r_max(izcan)
          fsv1d(izcan)=vfsv(izcan)*A_s/A_vl_max(izcan)
          fws1d(izcan)=svfw(izcan)
          frs1d(izcan)=svfr(izcan)
          fvs1d(izcan)=svfv(izcan)
          
          vfwt_f(izcan) = vfwt(izcan)
          vftw_f(izcan) = vftw(izcan)
          vftv_f(izcan) = vftv(izcan)
          vfvt_f(izcan) = vfvt(izcan)
          vfsw_f(izcan) = vfsw(izcan)
          vfsr_f(izcan) = vfsr(izcan)
          vfsv_f(izcan) = vfsv(izcan)
          svfw_f(izcan) = svfw(izcan)
          svfr_f(izcan) = svfr(izcan)
          svfv_f(izcan) = svfv(izcan)
          
       enddo
       fts1d=svft
       fsg1d=vfst*A_s/A_g
       
       svft_f=svft
       vfst_f=vfst

    endif
 
    call system_clock(end_time)        
    ! Calculate elapsed time in seconds
    elapsed_time = real(end_time - start_time,r8) / real(clock_rate,r8)
    
    if (debug_write) then
        if (solar) then !print calulation time after solar calculation finishs
           ! write the elapsed time for each iteration  
           write(*, '(A, I0, A, F6.3)') 'Calculation elapsed time for solar calculation l = ', l, ': ', elapsed_time, ' seconds'           
        else ! print longwave 
           ! write the elapsed time for each iteration  
           write(*, '(A, I0, A, F6.3)') 'Calculation elapsed time for longwave calculation l = ', l, ': ', elapsed_time, ' seconds'           
        endif
    end if    
  
       
  
    if (solar) goto 348
    !-----------------------------------------------------------------------
    ! If not solar, go back to calculate view factors for solar
    goto 358
   
348  continue
  end subroutine montecarlo_view_factors 

  !-----------------------------------------------------------------------
  subroutine corner_up(bldfrac,xdom,rayx,rayy,rayyint,zre,strike)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! This routine figures out whether a ray passing a corner at the upstream
    ! edge of a building strikes the roof and the wall of the layer above, or
    ! just the current wall layer (see Fig. 4 of Krayenhoff et al. 2014)
    ! ----------------------------------------------------------------------
    ! !ARGUMENTS:
    implicit none
    real(r8)        , intent(in)  :: bldfrac   ! Fraction of building width in total domain width   
    real(r8)        , intent(in)  :: xdom      ! Total domain width (canyon+building) normalized by dzcan
    real(r8)        , intent(in)  :: rayx      ! Ray horizontal position, in units of dzcan
    real(r8)        , intent(in)  :: rayy      ! Ray vertical position, in units of dzcan (ray_up passes -rayy)
    integer         , intent(in)  :: rayyint   ! Layer boundary just crossed (ray_up passes -rayyint_old)
    real(r8)        , intent(in)  :: zre       ! Zenith angle (atan) of the ray path; compared with the angle to the corner
    logical         , intent(out) :: strike    ! Logical flag indicating strike or not

    strike=.false.
    ! if the following is true, then the ray passing an upstream corner impinged on the roof
    ! (and the wall layer above the current wall layer) over the last ray step
    !if (atan((amod(rayx,xdom)-xdom*(1._r8-bldfrac))/ &
    if (atan(((rayx - xdom * floor(rayx /xdom))-xdom*(1._r8-bldfrac))/ &
                  max(1.e-6_r8,(real(rayyint,r8)-rayy))) > zre) then  
       strike=.true.
    endif
  end subroutine corner_up

  !-----------------------------------------------------------------------
  subroutine corner_dn(xdom,rayx,rayy,rayyint,zre,strike)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! This routine figures out whether a ray passing a corner at the downstream
    ! edge of a building strikes the roof or not (see Fig. 4 of Krayenhoff et al. 2014)
    ! ----------------------------------------------------------------------
    ! !ARGUMENTS:
    implicit none
    real(r8)        , intent(in)  :: xdom      ! Total domain width (canyon+building) normalized by dzcan
    real(r8)        , intent(in)  :: rayx      ! Ray horizontal position, in units of dzcan
    real(r8)        , intent(in)  :: rayy      ! Ray vertical position, in units of dzcan
    integer         , intent(in)  :: rayyint   ! Layer boundary just crossed
    real(r8)        , intent(in)  :: zre       ! Zenith angle (atan) of the ray path; compared with the angle to the corner
    logical         , intent(out) :: strike    ! Logical flag indicating strike or not

    strike=.false.
    ! if the following is true, then the ray passing a downstream corner impinged on the roof
    ! over the last ray step
    !if (atan(amod(rayx,xdom)/max(1.e-6_r8,(real(rayyint,r8)-rayy))).lt.zre) then 
    if (atan((rayx - xdom * floor(rayx /xdom))/ &
                 max(1.e-6_r8,(real(rayyint)-rayy))) < zre) then     
       strike=.true.
    endif
  end subroutine corner_dn    
    
  !-----------------------------------------------------------------------
  subroutine init_random_seed(random_nl)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Initializes the random number generator seed using the system clock.
    !----------------------------------------------------------------------- 
    implicit none   
    integer  , intent(in)              :: random_nl   ! a namelist input
    ! !LOCAL VARIABLES:
    integer                            :: i, n
    integer, dimension(:), allocatable :: seed
      
    ! Get the size of the seed
    call RANDOM_SEED(size = n)

    ! Allocate memory for the seed array
    allocate(seed(n))
      
    ! Compute the seed values based on a user-defined integer
    seed = random_nl + 37 * (/ (i - 1, i = 1, n) /)
    
    !-------if you want the seed to change in every run-------!
    ! integer :: clock
    ! seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    !---------------------------------------------------------!
    
    ! Set the random seed
    call RANDOM_SEED(put = seed)
      
    deallocate(seed)
  end subroutine init_random_seed

  !-----------------------------------------------------------------------
  subroutine ray_up(minray, maxind, rayx, rayy, bld, bld_old, raystr, &
                    xx, zz, hemi32, xdom, dray, dzcan, bldfrac, pb, &
                    kbs, lad_tree, lad, omega, horiz, vfw, vfv,h1,h2,foliage)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Computes the contribution to view factors of rays travelling "upwards"
    ! with a positive vertical component.For vegetation, assumes upward travelling rays leave
    ! leaves clumped within individual tree canopies, seeing a larger lad 
    ! above them than the average canyon lad.
    !-----------------------------------------------------------------------
    ! ARGUMENTS:
    implicit none
    integer,  parameter     :: nzcanm = 5     ! Maximum number of vertical levels at urban resolution
    real(r8), intent(in)    :: minray         ! Minimum ray strength    
    integer,  intent(in)    :: maxind         ! Index of the highest level with roofs or foliage; ray is dropped once
                                              ! rayyint > maxind+1, and its remaining strength becomes the sky view factor
    real(r8), intent(inout) :: rayx, rayy     ! Ray horizontal and VERTICAL position, both in units of dzcan
    real(r8), intent(inout) :: bld, bld_old   ! 1 if over a building column, 0 if over the canyon; this step and the last
    real(r8), intent(inout) :: raystr         ! Ray strength
    real(r8), intent(in)    :: xx, zz         ! x- and z-COMPONENTS of the ray direction unit vector (not coordinates)
    real(r8), intent(in)    :: hemi32         ! Ratio of 3D to 2D ray travel distance (f_i in Krayenhoff et al. 2014 Eq. 9)
    real(r8), intent(in)    :: xdom           ! Total domain width (canyon+building) normalized by dzcan
    real(r8), intent(in)    :: dray           ! Ray step, normalized by dzcan [-] (NOT metres)
    real(r8), intent(in)    :: dzcan          ! Vertical layer thickness [m]; converts normalized lengths back to metres
    real(r8), intent(in)    :: bldfrac        ! Fraction of building width in total domain width
    real(r8), intent(in)    :: pb(nzcanm + 1) ! Probability to have a building with an height equal or higher than each level
    real(r8), intent(in)    :: kbs            ! Foliage extinction coefficient [-]; 0.5 for a spherical leaf-angle distribution
    real(r8), intent(in)    :: lad_tree(nzcanm)! In-crown leaf area density [m2 m-3] = lad/wtroad_tree
    real(r8), intent(in)    :: lad(nzcanm)    ! Canyon-space-average leaf area density [m2 m-3] (not the in-crown value)
    real(r8), intent(in)    :: omega(nzcanm)  ! Leaf clumping index 
    logical,  intent(in)    :: horiz          ! Direction mapping: .true. => (rx,ry)=(xx,zz); .false. rotates onto a wall normal
    real(r8), intent(out)   :: vfw(nzcanm)    ! THIS ray's strength intercepted by each wall level
    real(r8), intent(out)   :: vfv(nzcanm)    ! THIS ray's strength intercepted by each foliage level
    logical,  intent(in)    :: foliage        ! .true. if the ray ORIGINATES in a foliage layer (selects the lad_tree
    real(r8),  intent(in)    :: h1                ! Crown bottom height above ground [m]; divided by dzcan internally
    real(r8),  intent(in)    :: h2                ! Crown DEPTH [m], not crown-top height; divided by dzcan internally
    integer              :: rayyint, rayyint2 ! Current vertical layer, and the wall layer being lit
    integer              :: rayyint_old       ! Previous vertical layer (ceiling value of rayy)
    integer              :: izcan             ! Indices
    real(r8)             :: rx, ry            ! Directions for ray traveling, adjusted based on horiz
    real(r8)             :: xfr               ! rayx folded into one canyon+building period, in units of dzcan
    real(r8)             :: wfact             ! Ray strength intercepted by the wall this step
    real(r8)             :: raystrtmp         ! Ray strength before foliage attenuation, to get the intercepted amount
    logical              :: strike            ! Corner test result: .true. => the lower wall layer was struck
    !-----------------------------------------------------------------------
    ! rays originating from walls
    rx=zz
    ry=xx
    
    ! rays originating from roads/vegetation
    if (horiz) then
       rx=xx
       ry=zz
    endif

    do izcan=1,nzcanm
       vfw(izcan)=0._r8
       vfv(izcan)=0._r8
    enddo

    do while (raystr > minray)
       ! Determine which vertical layer the ray is located in
       rayyint_old=ceiling(rayy)
       rayx=rayx+rx*dray
       rayy=rayy+ry*dray
       rayyint=ceiling(rayy)

       xfr = (rayx - xdom * floor(rayx /xdom))
       
       if (xfr/xdom >= 1._r8-bldfrac) then
          ! interception by walls
          bld=1._r8
          rayyint2=rayyint
          if (rayyint > rayyint_old) then !ray travels across vertical layers
             ! The negated rayy and rayyint_old mirror the upward geometry onto the downward case corner_up was
             ! written for, so that (rayyint - rayy) inside it is the distance past the boundary just crossed in
             ! strick = true means the lower wall layer was struck
             call corner_up(bldfrac,xdom,rayx,-rayy, -rayyint_old,&
                           atan(rx/max(1.e-6_r8,ry)),strike)
             if (strike) then
                ! lower wall layer is being lit
  		          rayyint2=rayyint_old
             endif
          endif
          wfact=raystr*max(bld-bld_old,0._r8)*pb(rayyint2+1)/max(1.e-6_r8,pb(rayyint2+1))*pb(rayyint2+1)
          vfw(rayyint2)=vfw(rayyint2)+wfact
          raystr=raystr-wfact
       else
          bld=0._r8
          ! interception by vegetation in the canopy column
          raystrtmp=raystr
          if ((rayy > h1/dzcan) .and. rayy< (h1+h2)/dzcan) then
             if (foliage) then
                raystr=raystr*exp(-hemi32*dray*dzcan*&
                      (kbs*0.5_r8*(lad(rayyint)+lad_tree(rayyint))*omega(rayyint)))    
             else      
                raystr=raystr*exp(-hemi32*dray*dzcan*(kbs*lad(rayyint)*omega(rayyint)))
             endif
          else
             raystr=raystr ! no attentuation
          end if
          vfv(rayyint)=vfv(rayyint)+(raystrtmp-raystr)
       endif
       if(rayyint > maxind+1) goto 541
       bld_old=bld
    enddo

541   continue
  end subroutine ray_up
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine ray_dn(minray, rayx, rayy, bld, bld_old, raystr, xx, zz, &
                    dist32, xdom, dray, dzcan, bldfrac, pb, pb_old, ss, &
                    kbs, lad, omega, horiz,h1,h2, vfw, vfr, vfv, vft)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Computes the contribution to view factors of rays travelling "downwards"
    ! with a negative vertical component.
    !-----------------------------------------------------------------------
    ! 
    ! !ARGUMENTS:
    implicit none
    integer,  parameter     :: nzcanm = 5     ! Maximum number of vertical levels at urban resolution
    real(r8), intent(in)    :: minray         ! Ray strength below which the ray is no longer tracked
    real(r8), intent(inout) :: rayx, rayy     ! Ray horizontal and vertical position, both in units of dzcan
    real(r8), intent(inout) :: bld, bld_old   ! 1 if over a building column, 0 if over the canyon; 
    real(r8), intent(inout) :: raystr         ! Ray strength
    real(r8), intent(in)    :: xx, zz         ! x- and z-components of the ray direction unit vector (not coordinates)
    real(r8), intent(in)    :: dist32         ! Ratio of 3D to 2D ray travel distance
    real(r8), intent(in)    :: xdom           ! Total domain width (canyon+building) normalized by dzcan
    real(r8), intent(in)    :: dray           ! Ray step, normalized by dzcan [-] 
    real(r8), intent(in)    :: dzcan          ! Vertical layer thickness [m]; converts normalized lengths back to metres
    real(r8), intent(in)    :: bldfrac        ! Fraction of building width in total domain width
    real(r8), intent(in)    :: pb(nzcanm + 1) ! Probability to have a building with an height equal or higher than each level
    real(r8), intent(inout) :: pb_old         ! Pb from the previous step (the pb already accounted for, so pbinc is the
                                              ! newly exposed building fraction); 
    real(r8), intent(in)    :: ss(nzcanm + 1) ! Roof fraction at each level 
    real(r8), intent(in)    :: kbs            ! Foliage extinction coefficient [-]; 0.5 for a spherical leaf-angle distribution
    real(r8), intent(in)    :: lad(nzcanm)    ! Canyon-space-average leaf area density [m2 m-3] (not the in-crown value)
    real(r8), intent(in)    :: omega(nzcanm)  ! Leaf clumping index 
    logical,  intent(in)    :: horiz          ! Direction mapping: .true. => (rx,ry)=(xx,zz); 
    real(r8),  intent(in)    :: h1                ! Crown bottom height above ground [m]; 
    real(r8),  intent(in)    :: h2                ! Crown DEPTH [m], not crown-top height; 
    real(r8), intent(out)   :: vfw(nzcanm)    ! THIS ray's strength intercepted by each wall level
    real(r8), intent(out)   :: vfr(nzcanm)    ! THIS ray's strength intercepted by each roof level
    real(r8), intent(out)   :: vfv(nzcanm)    ! THIS ray's strength intercepted by each foliage level
    real(r8), intent(out)   :: vft            ! THIS ray's strength reaching the ground

    ! !LOCAL VARIABLES:
    ! rayyint = vertical layer holding the ray (= ceiling(rayy)); rayyint2 = wall layer being lit, which is rayyint or rayyint+1 after an upstream-corner strike
    integer              :: rayyint, rayyint2 ! Current vertical layer, and the wall layer being lit
    integer              :: izcan      ! Indices 
    real(r8)             :: rx, ry            ! Directions for ray traveling, adjusted based on horiz
    real(r8)             :: xfr               ! rayx folded into one canyon+building period, in units of dzcan
    real(r8)             :: wfact             ! Ray strength intercepted by the wall this step
    real(r8)             :: rfact             ! Ray strength intercepted by the roof this step
    real(r8)             :: raystrtmp         ! Ray strength before foliage attenuation, to get the intercepted amount
    real(r8)             :: pbinc             ! Increment in building probability.
    real(r8)             :: blde              ! bld for the ROOF test only; forced to 1 after a downstream-corner roof strike
    real(r8)             :: sseff(nzcanm) ! Effective roof probability, rescaled onto the building-free column (Eq. 8)
    ! roof = .false. means this step lit a wall but no roof (upstream corner missed);
    ! strike = corner test result, .true. means the ray hit the roof/upper wall at the corner
    logical              :: roof, strike      ! Ray-tracing state flags
    !-----------------------------------------------------------------------
    ! rays originating from walls
    rx=zz
    ry=xx
    
    ! rays originating from sky/vegetation
    if (horiz) then
       rx=xx
       ry=zz
    endif
    
    do izcan=1,nzcanm
       vfw(izcan)=0._r8
       vfv(izcan)=0._r8
       vfr(izcan)=0._r8
    enddo
    
    vft=0._r8
   
    ! Effective roof probability 
    ! confirm sseff value with Scott
    ! sseff = (0, 1, 0, 0, 0)
    do izcan=1,nzcanm
       if (pb(izcan+1) < 1._r8) then
          sseff(izcan)=ss(izcan)/(1._r8-pb(izcan+1))
       else
          sseff(izcan)=ss(izcan)! effectively ss(izcan)=0
       endif
    enddo

    do while (raystr > minray)
       rayx=rayx+rx*dray
       rayy=rayy+ry*dray
       ! Determine which vertical layer the ray is located in
       rayyint=ceiling(rayy)
       rayyint2=rayyint
       ! the rayyint+1 index translate to the pb index
       ! rayyint counts layers, while pb counts levels
       pbinc=pb(rayyint+1)-pb_old
       ! Test where a ray is horizontally (i.e. building or canyon/tree space):
       ! xfr>bldfrac means canyon, xfr<bldfrac means building (i.e. bld=1)
       xfr = (rayx - xdom * floor(rayx /xdom))
       bld=0._r8
       if(xfr/xdom >= 1._r8-bldfrac) bld=1._r8
       blde=bld

       roof=.true.
 
       ! If the current ray step involves a corner (i.e., crossing both a pb boundary and the building canyon boundary)
       ! if we are crossing a column:
       if (abs(bld-bld_old) > 0.5_r8) then
          if (bld < bld_old) then  ! travel from a building column to a canyon column
             if (pbinc > 0._r8) then !the pb has changed from this step to previous step
                ! building corner at 'downstream' edge of the building
                ! refer to fig.4 of Krayenhoff et al. 2014
                call corner_dn(xdom,rayx,rayy,rayyint, atan(rx/max(1.e-6_r8,(-ry)))&
                               ,strike)
                if (strike) then !ray impinged on the roof
                    blde=1._r8
                    goto 223 !ROOF
                else
                    goto 224 !VEGETATION
                endif
             endif
           ! if we are entering the building column (bld.gt.bld_old)
          else ! travel from a canyon column to a building column
             if (ceiling(rayy-ry*dray)-rayyint > 0) then !the ray traveled across a verctial layer
                 ! building corner or wall layer division at 'upstream' edge of the building
                 call corner_up(bldfrac,xdom,rayx,rayy,rayyint,&
                               atan(rx/max(1.e-6_r8,(-ry))),strike)
                 if (strike) then !the ray impinged on the roof and the wall of the layer above
                    ! next wall layer higher is being lit (and roofs too)
  		               rayyint2=rayyint+1
                 else
                    ! only the present wall layer is being lit (not roofs)
                    roof=.false.
                 endif
             endif		
          endif          
       endif

       if(rayyint2 < 1) goto 223
       !! WALLS - Eq. 4: W(iz) = r_{j-1}*Lambda(iz+1), charged only on a canyon->building crossing, which is what
       !! max(bld-bld_old,0) detects. pb/max(1e-6,pb)*pb is a guarded way of writing pb (0 when pb=0, no 0/0).
       !! Only one wall is tracked since rays always march +x; the caller halves where the two walls must be split.
       wfact=raystr*max(bld-bld_old,0._r8)*pb(rayyint2+1)/max(1.e-6_r8,pb(rayyint2+1))*pb(rayyint2+1)
       vfw(rayyint2)=vfw(rayyint2)+wfact
       raystr=raystr-wfact
223      continue

       if (.not.roof) goto 224
       rfact=blde*raystr*pbinc/max(1.e-6_r8,pbinc)*sseff(rayyint+1)
       ! prevent negative ray strength
       rfact=min(rfact,raystr)        
                          
       vfr(rayyint+1)=vfr(rayyint+1)+rfact                 
       raystr=raystr-rfact
224      continue
       if(rayy <= 0._r8) goto 226
       !! VEGETATION
       ! so far I have not added in the details, e.g. what if a ray crosses vegetation layers or from or into a
       ! building during the ray step? As long as dray is quite small this should not be too important
       if (bld < 0.5_r8) then
          raystrtmp=raystr       
          ! check if the y axis coordinate is within tree crown height
          if ((rayy > h1/dzcan) .and. rayy< ((h1+h2)/dzcan)) then
              ! interception by vegetation in the canopy column
              raystr=raystr*exp(-dist32*dray*dzcan*(kbs*lad(rayyint)*omega(rayyint)))
          else 
          ! if outside of tree crown, the raystrength does not change
              raystr=raystr 
          endif
          vfv(rayyint)=vfv(rayyint)+(raystrtmp-raystr)
       endif
 
       bld_old=bld
      
       if (pbinc > 0._r8) then
          pb_old=pb(rayyint+1)
       endif
    enddo
    
226  continue    

    !! ROADS 
    vft=vft+raystr
    if (raystr > minray .and. ((rayx-rx*dray) - xdom * &
       floor((rayx-rx*dray)/xdom))/xdom > 1._r8-bldfrac .and. pb(2) > 0.999999_r8) then                                        
        write(iulog,*)'PROBLEM (ray_dn), radiation reaching building interior ground,raystr,rayx,rayy=',raystr,rayx,rayy
        write(iulog,*)'xdom,bldfrac',xdom,bldfrac
        write(iulog,*) ((rayx-rx*dray) - xdom * floor((rayx-rx*dray)/xdom))/xdom,1._r8-bldfrac
        write(iulog,*)'izcan',izcan
        write(iulog,*)'xx,zz,dist',xx,zz,dist32
        do izcan=1,nzcanm
           write(iulog,*)'i,ss,sseff',izcan,ss(izcan),sseff(izcan)
        enddo
        write(iulog,*)'minray, pb',minray,pb
        call endrun(msg=errmsg(sourcefile, __LINE__))
    endif
  end subroutine ray_dn
   
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: UrbanReadNML
  !
  ! !INTERFACE:
  !
  subroutine UrbanReadNML ( NLFilename )
    !
    ! !DESCRIPTION:
    !
    ! Read in the urban namelist
    !
    ! !USES:
    use shr_mpi_mod   , only : shr_mpi_bcast
    use abortutils    , only : endrun
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    implicit none
    !
    ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'UrbanReadNML'  ! subroutine name

    namelist / clmu_inparm / urban_hac, urban_explicit_ac, urban_traffic, building_temp_method
    !EOP
    !-----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! Read namelist from input namelist filename
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in clmu_inparm  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'clmu_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clmu_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading clmu_inparm namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR finding clmu_inparm namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )

    end if

    ! Broadcast namelist variables read in
    call shr_mpi_bcast(urban_hac,             mpicom)
    call shr_mpi_bcast(urban_explicit_ac,     mpicom)
    call shr_mpi_bcast(urban_traffic,         mpicom)
    call shr_mpi_bcast(building_temp_method,  mpicom)

    !
    if (urban_traffic) then
       write(iulog,*)'Urban traffic fluxes are not implemented currently'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    !
    if ( masterproc )then
       write(iulog,*) '   urban air conditioning/heating and wasteheat   = ', urban_hac
       write(iulog,*) '   urban explicit air-conditioning adoption rate   = ', urban_explicit_ac
       write(iulog,*) '   urban traffic flux   = ', urban_traffic
    end if

    ReadNamelist = .true.

  end subroutine UrbanReadNML

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: IsSimpleBuildTemp
  !
  ! !INTERFACE:
  !
  logical function IsSimpleBuildTemp( )
    !
    ! !DESCRIPTION:
    !
    ! If the simple building temperature method is being used
    !
    ! !USES:
    implicit none
    !EOP
    !-----------------------------------------------------------------------

    if ( .not. ReadNamelist )then
       write(iulog,*)'Testing on building_temp_method before urban namelist was read in'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    IsSimpleBuildTemp = building_temp_method == BUILDING_TEMP_METHOD_SIMPLE

  end function IsSimpleBuildTemp

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: IsProgBuildTemp
  !
  ! !INTERFACE:
  !
  logical function IsProgBuildTemp( )
    !
    ! !DESCRIPTION:
    !
    ! If the prognostic building temperature method is being used
    !
    ! !USES:
    implicit none
    !EOP
    !-----------------------------------------------------------------------

    if ( .not. ReadNamelist )then
       write(iulog,*)'Testing on building_temp_method before urban namelist was read in'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    IsProgBuildTemp = building_temp_method == BUILDING_TEMP_METHOD_PROG

  end function IsProgBuildTemp

  !-----------------------------------------------------------------------

end module UrbanParamsType
