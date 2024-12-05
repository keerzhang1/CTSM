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
  use clm_varcon   , only : grlnd, spval
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
  !
  !-------------------[kz.1]Ray tracing test-------------------------  
  ! PRIVATE MEMBER FUNCTIONS
  private :: view_factors_v    ! Calculates the view factors between all 'surfaces'  
  private :: corner_up         ! Handle a ray passing a corner at the upstream edge of a building
  private :: corner_dn         ! Handle a ray passing a corner at the downstream edge of a building
  private :: init_random_seed  ! Initialize a random seed
  private :: ray_up            ! Computes the contribution to view factors of rays travelling upwards
  private :: ray_dn            ! Computes the contribution to view factors of rays travelling downwards
  ! confirm:
  private :: ray_up_foliage    ! Consider upward travelling rays leave leaves that are clumped
  !-------------------[kz.1]Ray tracing test-------------------------  
  ! !PRIVATE TYPE
  type urbinp_type
     real(r8), pointer :: canyon_hwr      (:,:)  
     real(r8), pointer :: lai             (:,:)  
     real(r8), pointer :: tree_cov        (:,:)      
     real(r8), pointer :: wtlunit_roof    (:,:)  
     real(r8), pointer :: wtroad_perv     (:,:)  
     real(r8), pointer :: em_roof         (:,:)   
     real(r8), pointer :: em_improad      (:,:)  
     real(r8), pointer :: em_perroad      (:,:)  
     real(r8), pointer :: em_wall         (:,:)  
     real(r8), pointer :: alb_roof_dir    (:,:,:)  
     real(r8), pointer :: alb_roof_dif    (:,:,:)  
     real(r8), pointer :: alb_improad_dir (:,:,:)  
     real(r8), pointer :: alb_improad_dif (:,:,:)  
     real(r8), pointer :: alb_perroad_dir (:,:,:)  
     real(r8), pointer :: alb_perroad_dif (:,:,:)  
     real(r8), pointer :: alb_wall_dir    (:,:,:)  
     real(r8), pointer :: alb_wall_dif    (:,:,:)  
     real(r8), pointer :: ht_roof         (:,:)
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
 
     
  end type urbinp_type
  type (urbinp_type), public :: urbinp   ! urban input derived type

  ! !PUBLIC TYPE
  type, public :: urbanparams_type
!-------------------[kz.2]Ray tracing test-------------------------  
     real(r8), pointer :: fww1d_out       (:,:,:) ! longwave radiation view factor from wall to wall for a given canyon
     real(r8), pointer :: fvv1d_out       (:,:,:) ! view factor from vegetation to vegetation for a given canyon
     real(r8), pointer :: fwv1d_out       (:,:,:) ! view factor from wall to vegetation for a given canyon
     real(r8), pointer :: fvw1d_out       (:,:,:) ! view factor from vegetation to wall for a given canyon
     real(r8), pointer :: fwr1d_out       (:,:,:) ! view factor from wall to roof for a given canyon
     real(r8), pointer :: frw1d_out       (:,:,:) ! view factor from roof to wall for a given canyon
     real(r8), pointer :: fvr1d_out       (:,:,:) ! view factor from vegetation to roof for a given canyon
     real(r8), pointer :: frv1d_out       (:,:,:) ! view factor from roof to vegetation for a given canyon
     
     real(r8), pointer :: fwg1d_out       (:,:)   ! view factor from wall to ground for a given canyon
     real(r8), pointer :: fgw1d_out       (:,:)   ! view factor from ground to wall for a given canyon
     real(r8), pointer :: fgv1d_out       (:,:)   ! view factor from ground to vegetation for a given canyon
     real(r8), pointer :: fsw1d_out       (:,:)   ! view factor from sky to wall for a given canyon
     real(r8), pointer :: fvg1d_out       (:,:)   ! view factor from vegetation to ground for a given canyon
     real(r8), pointer :: fsr1d_out       (:,:)   ! view factor from sky to roof for a given canyon
     real(r8), pointer :: fsv1d_out       (:,:)   ! view factor from sky to vegetation for a given canyon
     
     real(r8), pointer :: fsg1d_out       (:)     ! view factor from sky to ground for a given canyon
     real(r8), pointer :: fws1d_out       (:,:)   ! view factor from wall to sky for a given canyon
     real(r8), pointer :: fvs1d_out       (:,:)   ! view factor from vegetation to sky for a given canyon
     real(r8), pointer :: fts1d_out       (:)     ! view factor from ? to sky for a given canyon
     real(r8), pointer :: frs1d_out       (:,:)   ! view factor from roof to sky for a given canyon
          
     real(r8), pointer :: kww1d_out       (:,:,:) ! shortwave radiation view factor from wall to wall for a given canyon
     real(r8), pointer :: kvv1d_out       (:,:,:) ! view factor from vegetation to vegetation for a given canyon
     real(r8), pointer :: kwv1d_out       (:,:,:) ! view factor from wall to vegetation for a given canyon
     real(r8), pointer :: kvw1d_out       (:,:,:) ! view factor from vegetation to wall for a given canyon
     real(r8), pointer :: kwr1d_out       (:,:,:) ! view factor from wall to roof for a given canyon
     real(r8), pointer :: krw1d_out       (:,:,:) ! view factor from roof to wall for a given canyon
     real(r8), pointer :: kvr1d_out       (:,:,:) ! view factor from vegetation to roof for a given canyon
     real(r8), pointer :: krv1d_out       (:,:,:) ! view factor from roof to vegetation for a given canyon
     
     real(r8), pointer :: kwg1d_out       (:,:)   ! view factor from wall to ground for a given canyon
     real(r8), pointer :: kgw1d_out       (:,:)   ! view factor from ground to wall for a given canyon
     real(r8), pointer :: kgv1d_out       (:,:)   ! view factor from ground to vegetation for a given canyon
     real(r8), pointer :: ksw1d_out       (:,:)   ! view factor from sky to wall for a given canyon
     real(r8), pointer :: kvg1d_out       (:,:)   ! view factor from vegetation to ground for a given canyon
     real(r8), pointer :: ksr1d_out       (:,:)   ! view factor from sky to roof for a given canyon
     real(r8), pointer :: ksv1d_out       (:,:)   ! view factor from sky to vegetation for a given canyon
     
     real(r8), pointer :: ksg1d_out       (:)   ! view factor from sky to ground for a given canyon
     real(r8), pointer :: kws1d_out       (:,:)   ! view factor from wall to sky for a given canyon
     real(r8), pointer :: kvs1d_out       (:,:)   ! view factor from vegetation to sky for a given canyon
     real(r8), pointer :: kts1d_out       (:)     ! view factor from ? to sky for a given canyon
     real(r8), pointer :: krs1d_out       (:,:)   ! view factor from roof to sky for a given canyon
!-------------------[kz.2]Ray tracing test-------------------------  
     real(r8), allocatable :: wind_hgt_canyon     (:)   ! lun height above road at which wind in canyon is to be computed (m)
     real(r8), allocatable :: em_roof             (:)   ! lun roof emissivity
     real(r8), allocatable :: em_improad          (:)   ! lun impervious road emissivity
     real(r8), allocatable :: em_perroad          (:)   ! lun pervious road emissivity
     real(r8), allocatable :: em_wall             (:)   ! lun wall emissivity
     real(r8), allocatable :: alb_roof_dir        (:,:) ! lun direct  roof albedo
     real(r8), allocatable :: alb_roof_dif        (:,:) ! lun diffuse roof albedo
     real(r8), allocatable :: alb_improad_dir     (:,:) ! lun direct  impervious road albedo
     real(r8), allocatable :: alb_improad_dif     (:,:) ! lun diffuse impervious road albedo
     real(r8), allocatable :: alb_perroad_dir     (:,:) ! lun direct  pervious road albedo
     real(r8), allocatable :: alb_perroad_dif     (:,:) ! lun diffuse pervious road albedo
     real(r8), allocatable :: alb_wall_dir        (:,:) ! lun direct  wall albedo
     real(r8), allocatable :: alb_wall_dif        (:,:) ! lun diffuse wall albedo

     integer , pointer     :: nlev_improad        (:)   ! lun number of impervious road layers (-)
     real(r8), pointer     :: tk_wall             (:,:) ! lun thermal conductivity of urban wall (W/m/K)
     real(r8), pointer     :: tk_roof             (:,:) ! lun thermal conductivity of urban roof (W/m/K)
     real(r8), pointer     :: tk_improad          (:,:) ! lun thermal conductivity of urban impervious road (W/m/K)
     real(r8), pointer     :: cv_wall             (:,:) ! lun heat capacity of urban wall (J/m^3/K)
     real(r8), pointer     :: cv_roof             (:,:) ! lun heat capacity of urban roof (J/m^3/K)
     real(r8), pointer     :: cv_improad          (:,:) ! lun heat capacity of urban impervious road (J/m^3/K)
     real(r8), pointer     :: thick_wall          (:)   ! lun total thickness of urban wall (m)
     real(r8), pointer     :: thick_roof          (:)   ! lun total thickness of urban roof (m)

     real(r8), pointer     :: vf_sr               (:)   ! lun view factor of sky for road
     real(r8), pointer     :: vf_wr               (:)   ! lun view factor of one wall for road
     real(r8), pointer     :: vf_sw               (:)   ! lun view factor of sky for one wall
     real(r8), pointer     :: vf_rw               (:)   ! lun view factor of road for one wall
     real(r8), pointer     :: vf_ww               (:)   ! lun view factor of opposing wall for one wall

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
    use clm_varcon      , only : vkc
    use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon   , only : icol_road_perv, icol_road_imperv, icol_road_perv
    use landunit_varcon , only : isturb_MIN
    !
    ! !ARGUMENTS:
    class(urbanparams_type) :: this
    type(bounds_type)      , intent(in)    :: bounds  
    !
    ! !LOCAL VARIABLES:   
    integer             :: j,l,c,p,g       ! indices
    integer             :: nc,fl,ib        ! indices 
    integer             :: dindx           ! urban density type index
    integer             :: ier             ! error status
    real(r8)            :: sumvf           ! sum of view factors for wall or road
    real(r8), parameter :: alpha = 4.43_r8 ! coefficient used to calculate z_d_town
    real(r8), parameter :: beta = 1.0_r8   ! coefficient used to calculate z_d_town
    real(r8), parameter :: C_d = 1.2_r8    ! drag coefficient as used in Grimmond and Oke (1999)
    real(r8)            :: plan_ai         ! plan area index - ratio building area to plan area (-)
    real(r8)            :: frontal_ai      ! frontal area index of buildings (-)
    real(r8)            :: build_lw_ratio  ! building short/long side ratio (-)
    integer		:: begl, endl
    integer		:: begc, endc
    integer		:: begp, endp
    integer             :: begg, endg
    !---------------------------------------------------------------------
    
!-------------------[kz.3]Ray tracing test-------------------------  
    integer, parameter  :: nzcanm = 4       ! Maximum number of vertical levels at urban resolution
    real(r8)            :: lad(nzcanm)           ! Leaf area density in the canyon column [m-1]
    real(r8)            :: lads(nzcanm)          ! Leaf area density in the canyon column [m-1] for shortwave calcs
    real(r8)            :: ladl(nzcanm)          ! Leaf area density in the canyon column [m-1] for longwave calcs
    real(r8)            :: omega(nzcanm)         ! Leaf clumping index      
    real(r8)            :: dzcan            ! height of buildings [m] (originally the thickness of vertical layers in multi-layer model)
    ! confirm:
    real(r8)            :: pb_in(nzcanm)         ! Probability to have a building with equal height? 
    real(r8)            :: ss_in(nzcanm)         ! Probability to have a surface area?
    real(r8)            :: wcan             ! Width of the canyons [m]
    real(r8)            :: wbui             ! Width of the buildings [m]
    ! question: why the tree_cov(11) has a size of 11?
    real(r8)            :: tree_cov(11)     ! Canopy cover within the canyon (between-building) space
    real(r8)            :: dray             ! Ray step
    ! confirm: the ray from sky stars at a location at hsky times the highest tree/building, spread out evenly in the x-direction 
    real(r8)            :: hsky             ! A height scaling factor for sky ray calculations
    integer             :: iurb             ! LCZ
    ! confirm:
    integer             :: maxbhind         ! Maximum building layer index?
    integer             :: maxind           ! Maximum vertical layer index?
    integer             :: nrays            ! Number of rays from foliage layer (# from surfaces will be half);
    integer             :: nsky             ! Number of rays from sky

    !Output of view factor calculation
    !------
    ! Area-weighted longwave view factors
    real(r8)            :: fww1d(nzcanm, nzcanm)  ! View factor from wall to wall for a given canyon
    real(r8)            :: fvv1d(nzcanm, nzcanm)  ! View factor from vegetation to vegetation for a given canyon
    real(r8)            :: fwv1d(nzcanm, nzcanm)  ! View factor from wall to vegetation for a given canyon
    real(r8)            :: fvw1d(nzcanm, nzcanm)  ! View factor from vegetation to wall for a given canyon
    real(r8)            :: fwr1d(nzcanm, nzcanm)  ! View factor from wall to roof for a given canyon
    real(r8)            :: frw1d(nzcanm, nzcanm)  ! View factor from roof to wall for a given canyon
    real(r8)            :: fvr1d(nzcanm, nzcanm)  ! View factor from vegetation to roof for a given canyon
    real(r8)            :: frv1d(nzcanm, nzcanm)  ! View factor from roof to vegetation for a given canyon
    real(r8)            :: fwg1d(nzcanm)          ! View factor from wall to ground for a given canyon
    real(r8)            :: fgw1d(nzcanm)          ! View factor from ground to wall for a given canyon
    real(r8)            :: fgv1d(nzcanm)          ! View factor from ground to vegetation for a given canyon
    real(r8)            :: fsw1d(nzcanm)          ! View factor from sky to wall for a given canyon
    real(r8)            :: fvg1d(nzcanm)          ! View factor from vegetation to ground for a given canyon
    real(r8)            :: fsg1d                  ! View factor from sky to ground for a given canyon
    real(r8)            :: fsr1d(nzcanm)          ! View factor from sky to roof for a given canyon
    real(r8)            :: fsv1d(nzcanm)          ! View factor from sky to vegetation for a given canyon
    real(r8)            :: fws1d(nzcanm)          ! View factor from wall to sky for a given canyon
    real(r8)            :: fvs1d(nzcanm)          ! View factor from vegetation to sky for a given canyon
    ! confirm: 
    real(r8)            :: fts1d                  ! View factor from ground to sky for a given canyon
    real(r8)            :: frs1d(nzcanm)          ! View factor from roof to sky for a given canyon
    ! Area-weighted shortwave view factors
    real(r8)            :: kww1d(nzcanm, nzcanm)  ! View factor from wall to wall for a given canyon
    real(r8)            :: kvv1d(nzcanm, nzcanm)  ! View factor from vegetation to vegetation for a given canyon
    real(r8)            :: kwv1d(nzcanm, nzcanm)  ! View factor from wall to vegetation for a given canyon
    real(r8)            :: kvw1d(nzcanm, nzcanm)  ! View factor from vegetation to wall for a given canyon
    real(r8)            :: kwr1d(nzcanm, nzcanm)  ! View factor from wall to roof for a given canyon
    real(r8)            :: krw1d(nzcanm, nzcanm)  ! View factor from roof to wall for a given canyon
    real(r8)            :: kvr1d(nzcanm, nzcanm)  ! View factor from vegetation to roof for a given canyon
    real(r8)            :: krv1d(nzcanm, nzcanm)  ! View factor from roof to vegetation for a given canyon
    real(r8)            :: kwg1d(nzcanm)          ! View factor from wall to ground for a given canyon
    real(r8)            :: kgw1d(nzcanm)          ! View factor from ground to wall for a given canyon
    real(r8)            :: kgv1d(nzcanm)          ! View factor from ground to vegetation for a given canyon
    real(r8)            :: ksw1d(nzcanm)          ! View factor from sky to wall for a given canyon
    real(r8)            :: kvg1d(nzcanm)          ! View factor from vegetation to ground for a given canyon
    real(r8)            :: ksg1d                  ! View factor from sky to ground for a given canyon
    real(r8)            :: ksr1d(nzcanm)          ! View factor from sky to roof for a given canyon
    real(r8)            :: ksv1d(nzcanm)          ! View factor from sky to vegetation for a given canyon

    real(r8)            :: kws1d(nzcanm)          ! View factor from wall to sky for a given canyon
    real(r8)            :: kvs1d(nzcanm)          ! View factor from vegetation to sky for a given canyon
    ! confirm: 
    real(r8)            :: kts1d                  ! View factor from ground to sky for a given canyon
    real(r8)            :: krs1d(nzcanm)          ! View factor from roof to sky for a given canyon        
!-------------------[kz.3]Ray tracing test-------------------------   
    
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
    allocate(this%vf_sr               (begl:endl))          ; this%vf_sr               (:)   = nan
    allocate(this%vf_wr               (begl:endl))          ; this%vf_wr               (:)   = nan
    allocate(this%vf_sw               (begl:endl))          ; this%vf_sw               (:)   = nan
    allocate(this%vf_rw               (begl:endl))          ; this%vf_rw               (:)   = nan
    allocate(this%vf_ww               (begl:endl))          ; this%vf_ww               (:)   = nan
    allocate(this%wind_hgt_canyon     (begl:endl))          ; this%wind_hgt_canyon     (:)   = nan
    allocate(this%em_roof             (begl:endl))          ; this%em_roof             (:)   = nan
    allocate(this%em_improad          (begl:endl))          ; this%em_improad          (:)   = nan
    allocate(this%em_perroad          (begl:endl))          ; this%em_perroad          (:)   = nan
    allocate(this%em_wall             (begl:endl))          ; this%em_wall             (:)   = nan
    allocate(this%alb_roof_dir        (begl:endl,numrad))   ; this%alb_roof_dir        (:,:) = nan
    allocate(this%alb_roof_dif        (begl:endl,numrad))   ; this%alb_roof_dif        (:,:) = nan    
    allocate(this%alb_improad_dir     (begl:endl,numrad))   ; this%alb_improad_dir     (:,:) = nan       
    allocate(this%alb_perroad_dir     (begl:endl,numrad))   ; this%alb_perroad_dir     (:,:) = nan       
    allocate(this%alb_improad_dif     (begl:endl,numrad))   ; this%alb_improad_dif     (:,:) = nan       
    allocate(this%alb_perroad_dif     (begl:endl,numrad))   ; this%alb_perroad_dif     (:,:) = nan       
    allocate(this%alb_wall_dir        (begl:endl,numrad))   ; this%alb_wall_dir        (:,:) = nan    
    allocate(this%alb_wall_dif        (begl:endl,numrad))   ; this%alb_wall_dif        (:,:) = nan
    allocate(this%eflx_traffic_factor (begl:endl))          ; this%eflx_traffic_factor (:)   = nan
!-------------------[kz.4]Ray tracing test-------------------------  
    allocate(this%fww1d_out          (begl:endl,nzcanm,nzcanm))          ; this%fww1d_out       (:,:,:) = nan
    allocate(this%fvv1d_out          (begl:endl,nzcanm,nzcanm))          ; this%fvv1d_out       (:,:,:) = nan
    allocate(this%fwv1d_out          (begl:endl,nzcanm,nzcanm))          ; this%fwv1d_out       (:,:,:) = nan
    allocate(this%fvw1d_out          (begl:endl,nzcanm,nzcanm))          ; this%fvw1d_out       (:,:,:) = nan
    allocate(this%fwr1d_out          (begl:endl,nzcanm,nzcanm))          ; this%fwr1d_out       (:,:,:) = nan
    allocate(this%frw1d_out          (begl:endl,nzcanm,nzcanm))          ; this%frw1d_out       (:,:,:) = nan
    allocate(this%fvr1d_out          (begl:endl,nzcanm,nzcanm))          ; this%fvr1d_out       (:,:,:) = nan
    allocate(this%frv1d_out          (begl:endl,nzcanm,nzcanm))          ; this%frv1d_out       (:,:,:) = nan

    allocate(this%fwg1d_out          (begl:endl,nzcanm))    ; this%fwg1d_out       (:,:) = nan 
    allocate(this%fgw1d_out          (begl:endl,nzcanm))     ; this%fgw1d_out       (:,:) = nan 
    allocate(this%fgv1d_out          (begl:endl,nzcanm))     ; this%fgv1d_out       (:,:) = nan 
    allocate(this%fsw1d_out          (begl:endl,nzcanm))     ; this%fsw1d_out       (:,:) = nan 
    allocate(this%fvg1d_out          (begl:endl,nzcanm))     ; this%fvg1d_out       (:,:) = nan 
    allocate(this%fsr1d_out          (begl:endl,nzcanm))     ; this%fsr1d_out       (:,:) = nan 
    allocate(this%fsv1d_out          (begl:endl,nzcanm))     ; this%fsv1d_out       (:,:) = nan 

    allocate(this%fsg1d_out          (begl:endl))            ; this%fsg1d_out       (:) = nan   
    allocate(this%fws1d_out          (begl:endl,nzcanm))     ; this%fws1d_out       (:,:) = nan
    allocate(this%fvs1d_out          (begl:endl,nzcanm))     ; this%fvs1d_out       (:,:) = nan
    allocate(this%fts1d_out          (begl:endl))            ; this%fts1d_out       (:) = nan
    allocate(this%frs1d_out          (begl:endl,nzcanm))     ; this%frs1d_out       (:,:) = nan
        
    allocate(this%kww1d_out          (begl:endl,nzcanm,nzcanm))          ; this%kww1d_out       (:,:,:) = nan
    allocate(this%kvv1d_out          (begl:endl,nzcanm,nzcanm))          ; this%kvv1d_out       (:,:,:) = nan
    allocate(this%kwv1d_out          (begl:endl,nzcanm,nzcanm))          ; this%kwv1d_out       (:,:,:) = nan
    allocate(this%kvw1d_out          (begl:endl,nzcanm,nzcanm))          ; this%kvw1d_out       (:,:,:) = nan
    allocate(this%kwr1d_out          (begl:endl,nzcanm,nzcanm))          ; this%kwr1d_out       (:,:,:) = nan
    allocate(this%krw1d_out          (begl:endl,nzcanm,nzcanm))          ; this%krw1d_out       (:,:,:) = nan
    allocate(this%kvr1d_out          (begl:endl,nzcanm,nzcanm))          ; this%kvr1d_out       (:,:,:) = nan
    allocate(this%krv1d_out          (begl:endl,nzcanm,nzcanm))          ; this%krv1d_out       (:,:,:) = nan

    allocate(this%kwg1d_out          (begl:endl,nzcanm))     ; this%kwg1d_out       (:,:) = nan 
    allocate(this%kgw1d_out          (begl:endl,nzcanm))     ; this%kgw1d_out       (:,:) = nan 
    allocate(this%kgv1d_out          (begl:endl,nzcanm))     ; this%kgv1d_out       (:,:) = nan 
    allocate(this%ksw1d_out          (begl:endl,nzcanm))     ; this%ksw1d_out       (:,:) = nan 
    allocate(this%kvg1d_out          (begl:endl,nzcanm))     ; this%kvg1d_out       (:,:) = nan 
    allocate(this%ksr1d_out          (begl:endl,nzcanm))     ; this%ksr1d_out       (:,:) = nan 
    allocate(this%ksv1d_out          (begl:endl,nzcanm))     ; this%ksv1d_out       (:,:) = nan 

    allocate(this%ksg1d_out          (begl:endl))            ; this%ksg1d_out       (:) = nan   
    allocate(this%kws1d_out          (begl:endl,nzcanm))     ; this%kws1d_out       (:,:) = nan
    allocate(this%kvs1d_out          (begl:endl,nzcanm))     ; this%kvs1d_out       (:,:) = nan
    allocate(this%kts1d_out          (begl:endl))            ; this%kts1d_out       (:) = nan
    allocate(this%krs1d_out          (begl:endl,nzcanm))     ; this%krs1d_out       (:,:) = nan

   !-------------------------------
   ! Do not change these values (to maintain single-layer version of code):
   !-------------------------------
    !nzcanm=2  ! Maximum number of layers (one canyon layer, and one above)
    maxbhind=2
    maxind=2
    ! These values give a single-layer canyon:
    ss_in=0._r8
    ss_in(2)=1._r8
    pb_in=0._r8
    pb_in(1)=1._r8
    pb_in(2)=1._r8
    lad=0._r8  ! Average leaf area density within the canyon space
    tree_cov=0._r8  ! Canopy cover within the canyon (between-building) space
    ! If view factors from sky are not accurate, could try changing, especially nsky:
    nsky=1
    hsky=1.5_r8
    !-------------------------------

    !! Hardcoded inputs that can be changed:
    ! building geometry
    !wcan=10.  ! canyon width
    !wbui=10.  ! building width
    !dzcan=10.  ! height of buildings (originally the thickness of vertical layers in multi-layer model)
    ! Change to change tree cover (make sure LAD, tree_cov, and omega all correspond):
    
    !lad(1)=0.5  ! <- Set the LAD (m2/m3) in the canyon here
    !lad(2)=0.  ! <- Set the LAD above the canyon here
    
    !      lad(1)=0.  ! <- Set the LAD in the canyon here
    !      lad(2)=0.5  ! <- Set the LAD above the canyon here
    
    iurb=6  ! LCZ = 6
    !tree_cov(iurb)=0.5  ! <- Set canopy cover for "iurb" LCZ here
    
    omega=0.6_r8  ! Tree crown scale clumping factor (see Eq. 15 in Krayenhoff et al. 2020)
    
    ! Change to change the accuracy of the simulation:
    ! nrays=5000 is often sufficient, whereas 50000 rays should give very accurate results
    !      nrays=5000  ! # of rays from foliage layer (# from surfaces will be half); use more for accuracy - see Krayenhoff et al 2014

    nrays=5000  ! # of rays from foliage layer (# from surfaces will be half); use more for accuracy - see Krayenhoff et al 2014
    ! This value is currently overwritten in the subroutine:
    dray=0.05_r8  ! Ray step; smaller values give better results; for highrise neighbourhoods this value needs to be smaller

    ! Mostly unused but keep for now:
    !lads=lad  ! lad of tree foliage layer for shortwave calcs (usually equal to "lad")
    !ladl=lad  ! lad of tree foliage layer for longwave calcs (usually equal to "lad")
!-------------------[kz.4]Ray tracing test-------------------------  
          
    ! Initialize time constant urban variables
    write(6,*)'nrays = ',nrays
    write(6,*)'dray = ',dray
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
             this%alb_improad_dif(l,ib) = urbinp%alb_improad_dif(g,dindx,ib)
             this%alb_perroad_dif(l,ib) = urbinp%alb_perroad_dif(g,dindx,ib)
             this%alb_wall_dir   (l,ib) = urbinp%alb_wall_dir   (g,dindx,ib)
             this%alb_wall_dif   (l,ib) = urbinp%alb_wall_dif   (g,dindx,ib)
          end do
          this%em_roof   (l) = urbinp%em_roof   (g,dindx)
          this%em_improad(l) = urbinp%em_improad(g,dindx)
          this%em_perroad(l) = urbinp%em_perroad(g,dindx)
          this%em_wall   (l) = urbinp%em_wall   (g,dindx)

          ! Landunit level initialization for urban wall and roof layers and interfaces

          lun%canyon_hwr(l)   = urbinp%canyon_hwr(g,dindx)
          lun%lai(l)          = urbinp%lai(g,dindx)
          lun%tree_cov(l)     = urbinp%tree_cov(g,dindx)          
          lun%wtroad_perv(l)  = urbinp%wtroad_perv(g,dindx)
          lun%ht_roof(l)      = urbinp%ht_roof(g,dindx)
          lun%wtlunit_roof(l) = urbinp%wtlunit_roof(g,dindx)

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
          !-------------------[kz.5]Ray tracing test-------------------------
          ! I tried to add urbanparams_inst as an input in the view_factors_v subroutine
          ! and then save the output view factor as urbanparams_inst%kww1d
          ! This did not work somwhow. Some urbanparams_inst%kww1d is nan (1e36)
          ! After I output urbanparams_inst%kww1d by writing an "InitHistory" with 
          ! hist_addfld2d functions. The output are all nan.

          ! Now, I output the calculated view factor first, and then save then as 
          ! temperature_inst%kww1d_out1 in the TemperatureType.F90. This may not be the neatest way but works..
          
          !-----------------use CLM surface data to calculate view factor -------------------------
          wcan=lun%ht_roof(l)/lun%canyon_hwr(l)  ! canyon width = building height /hwr
          
          ! Building height to building width ratio
          ! building_hwr(l) = canyon_hwr(l)*(1._r8-wtlunit_roof(l))/wtlunit_roof(l)
          ! building width = building height / Building height to building width ratio
          wbui = lun%ht_roof(l)/(lun%canyon_hwr(l)*(1._r8-lun%wtlunit_roof(l))/lun%wtlunit_roof(l))
          dzcan=lun%ht_roof(l)  ! height of buildings (originally the thickness of vertical layers in multi-layer model)

          ! assume LAI = 3 m3/m3
          lad(1)=lun%lai(l)/lun%ht_roof(l)  ! <- Set the LAD (m2/m3) in the canyon here
          lad(2)=0._r8  ! <- Set the LAD above the canyon here
          
          ! Mostly unused but keep for now:
          lads=lad  ! lad of tree foliage layer for shortwave calcs (usually equal to "lad")
          ladl=lad  ! lad of tree foliage layer for longwave calcs (usually equal to "lad")
          
          tree_cov(iurb)=lun%tree_cov(l)
          ! Omega = -1.0 / (0.5 * LAI_can) * log_term
          ! log_term = log(1.0 - tree_cov(iurb) * (1.0 - exp(-0.5 * LAI_tree)))
          ! LAI_tree =lun%lai(l)/lun%tree_cov(l)
          omega=-1.0_r8 / (0.5_r8 * lun%lai(l)) * log(1.0_r8 - tree_cov(iurb) * &
                (1.0_r8 - exp(-0.5_r8 * lun%lai(l)/lun%tree_cov(l))))

          call view_factors_v(nzcanm,dzcan,wcan,wbui,&
                  tree_cov,lad,lads,ladl,omega,ss_in,pb_in,dray,maxind,&
                          maxbhind,nrays,nsky,hsky,&
                          iurb,fww1d,fvv1d,fwv1d,fvw1d,fwr1d,frw1d,fvr1d,&
                          frv1d,fwg1d,fgw1d,fgv1d,fsw1d,fvg1d,fsg1d,fsr1d,&
                          fsv1d,kww1d,kvv1d,kwv1d,kvw1d,kwr1d,krw1d,kvr1d,&
                          krv1d,kwg1d,kgw1d,kgv1d,ksw1d,kvg1d,ksg1d,ksr1d,&
                          ksv1d,kws1d,kvs1d,kts1d,krs1d,fws1d,fvs1d,fts1d,frs1d,l)
                
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

          !write(6,*)'fwr1d(:,:) calculated',fwr1d(:,:)
          !write(6,*)'this%fwr1d_out(l,:,:) saved',this%fwr1d_out(l,:,:)
          !-------------------[kz.5]Ray tracing test------------------------- 
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
          ! View factors for road and one wall in urban canyon (depends only on canyon_hwr)
          ! ---------------------------------------------------------------------------------------
          !                                                        WALL    |
          !                  ROAD                                          |
          !                                                         wall   |
          !          -----\          /-----   -             -  |\----------/
          !              | \  vsr   / |       |         r   |  | \  vww   /   s
          !              |  \      /  |       h         o   w  |  \      /    k
          !        wall  |   \    /   | wall  |         a   |  |   \    /     y
          !              |vwr \  / vwr|       |         d   |  |vrw \  / vsw 
          !              ------\/------       -             -  |-----\/-----
          !                   road                                  wall   |
          !              <----- w ---->                                    |
          !                                                    <---- h --->|
          !
          !    vsr = view factor of sky for road          vrw = view factor of road for wall
          !    vwr = view factor of one wall for road     vww = view factor of opposing wall for wall
          !                                               vsw = view factor of sky for wall
          !    vsr + vwr + vwr = 1                        vrw + vww + vsw = 1
          !
          ! Source: Masson, V. (2000) A physically-based scheme for the urban energy budget in 
          ! atmospheric models. Boundary-Layer Meteorology 94:357-397
          !
          ! - Calculate urban land unit aerodynamic constants using Macdonald (1998) as used in
          ! Grimmond and Oke (1999)
          ! ---------------------------------------------------------------------------------------
          
          ! road -- sky view factor -> 1 as building height -> 0 
          ! and -> 0 as building height -> infinity

          this%vf_sr(l) = sqrt(lun%canyon_hwr(l)**2 + 1._r8) - lun%canyon_hwr(l)
          this%vf_wr(l) = 0.5_r8 * (1._r8 - this%vf_sr(l))

          ! one wall -- sky view factor -> 0.5 as building height -> 0 
          ! and -> 0 as building height -> infinity

          this%vf_sw(l) = 0.5_r8 * (lun%canyon_hwr(l) + 1._r8 - sqrt(lun%canyon_hwr(l)**2+1._r8)) / lun%canyon_hwr(l)
          this%vf_rw(l) = this%vf_sw(l)
          this%vf_ww(l) = 1._r8 - this%vf_sw(l) - this%vf_rw(l)

          ! error check -- make sure view factor sums to one for road and wall
          sumvf = this%vf_sr(l) + 2._r8*this%vf_wr(l)
          if (abs(sumvf-1._r8) > 1.e-06_r8 ) then
             write (iulog,*) 'urban road view factor error',sumvf
             write (iulog,*) 'clm model is stopping'
             call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
          endif
          sumvf = this%vf_sw(l) + this%vf_rw(l) + this%vf_ww(l)
          if (abs(sumvf-1._r8) > 1.e-06_r8 ) then
             write (iulog,*) 'urban wall view factor error',sumvf
             write (iulog,*) 'clm model is stopping'
             call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errmsg(sourcefile, __LINE__))
          endif

          !----------------------------------------------------------------------------------
          ! Calculate urban land unit aerodynamic constants using Macdonald (1998) as used in
          ! Grimmond and Oke (1999)
          !----------------------------------------------------------------------------------

          ! Calculate plan area index 
          plan_ai = lun%canyon_hwr(l)/(lun%canyon_hwr(l) + 1._r8)

          ! Building shape shortside/longside ratio (e.g. 1 = square )
          ! This assumes the building occupies the entire canyon length
          build_lw_ratio = plan_ai

          ! Calculate frontal area index
          frontal_ai = (1._r8 - plan_ai) * lun%canyon_hwr(l)

          ! Adjust frontal area index for different building configuration
          frontal_ai = frontal_ai * sqrt(1/build_lw_ratio) * sqrt(plan_ai)

          ! Calculate displacement height
          if (use_vancouver) then
             lun%z_d_town(l) = 3.5_r8
          else if (use_mexicocity) then
             lun%z_d_town(l) = 10.9_r8
          else
             lun%z_d_town(l) = (1._r8 + alpha**(-plan_ai) * (plan_ai - 1._r8)) * lun%ht_roof(l)
          end if

          ! Calculate the roughness length
          if (use_vancouver) then
             lun%z_0_town(l) = 0.35_r8
          else if (use_mexicocity) then
             lun%z_0_town(l) = 2.2_r8
          else
             lun%z_0_town(l) = lun%ht_roof(l) * (1._r8 - lun%z_d_town(l) / lun%ht_roof(l)) * &
                  exp(-1.0_r8 * (0.5_r8 * beta * C_d / vkc**2 * &
                  (1 - lun%z_d_town(l) / lun%ht_roof(l)) * frontal_ai)**(-0.5_r8))
          end if

       else ! Not urban point 

          this%eflx_traffic_factor(l) = spval
          this%t_building_min(l) = spval

          this%vf_sr(l) = spval
          this%vf_wr(l) = spval
          this%vf_sw(l) = spval
          this%vf_rw(l) = spval
          this%vf_ww(l) = spval
!-------------------[kz.6]Ray tracing test-------------------------     
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
!-------------------[kz.6]Ray tracing test-------------------------     
       end if
    end do

    ! Note that we don't deallocate memory for urbinp datatype (call UrbanInput with
    ! mode='finalize') because the arrays are needed for dynamic urban landunits.
    do l = bounds%begl,bounds%endl
       if (lun%urbpoi(l)) then    
           write(6,*)'do loop saved',this%kwg1d_out(l,:)
       end if 
    end do
    
  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine UrbanInput(begg, endg, mode)
    !
    ! !DESCRIPTION: 
    ! Allocate memory and read in urban input data
    !
    ! !USES:
    use clm_varpar      , only : numrad, nlevurb
    use landunit_varcon , only : numurbl
    use fileutils       , only : getavu, relavu, getfil, opnfil
    use spmdMod         , only : masterproc
    use domainMod       , only : ldomain
    use ncdio_pio       , only : file_desc_t, ncd_io, ncd_inqvdlen, ncd_inqfdims 
    use ncdio_pio       , only : ncd_pio_openfile, ncd_pio_closefile, ncd_inqdid, ncd_inqdlen
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
                urbinp%lai(begg:endg, numurbl), &  
                urbinp%tree_cov(begg:endg, numurbl), &         
                urbinp%wtlunit_roof(begg:endg, numurbl), &  
                urbinp%wtroad_perv(begg:endg, numurbl), &
                urbinp%em_roof(begg:endg, numurbl), &     
                urbinp%em_improad(begg:endg, numurbl), &    
                urbinp%em_perroad(begg:endg, numurbl), &    
                urbinp%em_wall(begg:endg, numurbl), &    
                urbinp%alb_roof_dir(begg:endg, numurbl, numrad), &    
                urbinp%alb_roof_dif(begg:endg, numurbl, numrad), &    
                urbinp%alb_improad_dir(begg:endg, numurbl, numrad), &    
                urbinp%alb_perroad_dir(begg:endg, numurbl, numrad), &    
                urbinp%alb_improad_dif(begg:endg, numurbl, numrad), &    
                urbinp%alb_perroad_dif(begg:endg, numurbl, numrad), &    
                urbinp%alb_wall_dir(begg:endg, numurbl, numrad), &    
                urbinp%alb_wall_dif(begg:endg, numurbl, numrad), &
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
!-------------------[kz.7]Ray tracing test-------------------------                    
       if (ier /= 0) then
          call endrun(msg="Allocation error "//errmsg(sourcefile, __LINE__))
       endif

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
       call ncd_io(ncid=ncid, varname='LAI', flag='read', data=urbinp%lai,&
           dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
         call endrun( msg='ERROR: LAI NOT on fsurdat file '//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='TREE_COV', flag='read', data=urbinp%tree_cov,&
           dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
         call endrun( msg='ERROR: TREE_COV NOT on fsurdat file '//errmsg(sourcefile, __LINE__))
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
!-------------------[kz.9]Ray tracing test-------------------------        
       deallocate(urbinp%canyon_hwr, &
                  urbinp%lai, &     
                  urbinp%tree_cov, &     
                  urbinp%wtlunit_roof, &
                  urbinp%wtroad_perv, &
                  urbinp%em_roof, &
                  urbinp%em_improad, &
                  urbinp%em_perroad, &
                  urbinp%em_wall, &
                  urbinp%alb_roof_dir, &
                  urbinp%alb_roof_dif, &
                  urbinp%alb_improad_dir, &
                  urbinp%alb_perroad_dir, &
                  urbinp%alb_improad_dif, &
                  urbinp%alb_perroad_dif, &
                  urbinp%alb_wall_dir, &
                  urbinp%alb_wall_dif, &
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
!-------------------[kz.9]Ray tracing test-------------------------                    
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
!-------------------[kz.9]Ray tracing test-------------------------  
    found = .false.
    do nl = begg,endg
       do n = 1, numurbl
          if ( pcturb(nl,n) > 0.0_r8 ) then
             if ( .not. urban_valid(nl) .or. &
                  urbinp%canyon_hwr(nl,n)            <= 0._r8 .or. &
                  urbinp%lai(nl,n)                   <= 0._r8 .or. &
                  urbinp%tree_cov(nl,n)              <= 0._r8 .or. &                  
                  urbinp%em_improad(nl,n)            <= 0._r8 .or. &
                  urbinp%em_perroad(nl,n)            <= 0._r8 .or. &
                  urbinp%em_roof(nl,n)               <= 0._r8 .or. &
                  urbinp%em_wall(nl,n)               <= 0._r8 .or. &
                  urbinp%ht_roof(nl,n)               <= 0._r8 .or. &
                  urbinp%thick_roof(nl,n)            <= 0._r8 .or. &
                  urbinp%thick_wall(nl,n)            <= 0._r8 .or. &
                  urbinp%t_building_min(nl,n)        <= 0._r8 .or. &
                  urbinp%wind_hgt_canyon(nl,n)       <= 0._r8 .or. &
                  urbinp%wtlunit_roof(nl,n)          <= 0._r8 .or. &
                  urbinp%wtroad_perv(nl,n)           <= 0._r8 .or. &
                  any(urbinp%alb_improad_dir(nl,n,:) <= 0._r8) .or. &
                  any(urbinp%alb_improad_dif(nl,n,:) <= 0._r8) .or. &
                  any(urbinp%alb_perroad_dir(nl,n,:) <= 0._r8) .or. &
                  any(urbinp%alb_perroad_dif(nl,n,:) <= 0._r8) .or. &
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
!-------------------[kz.9]Ray tracing test-------------------------                  
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
!-------------------[kz.10]Ray tracing test-------------------------         
       write(iulog,*)'lai:             ',urbinp%lai(nindx,dindx)
       write(iulog,*)'tree_cov:        ',urbinp%tree_cov(nindx,dindx)   
!-------------------[kz.10]Ray tracing test-------------------------             
       write(iulog,*)'em_improad:      ',urbinp%em_improad(nindx,dindx)
       write(iulog,*)'em_perroad:      ',urbinp%em_perroad(nindx,dindx)
       write(iulog,*)'em_roof:         ',urbinp%em_roof(nindx,dindx)
       write(iulog,*)'em_wall:         ',urbinp%em_wall(nindx,dindx)
       write(iulog,*)'ht_roof:         ',urbinp%ht_roof(nindx,dindx)
       write(iulog,*)'thick_roof:      ',urbinp%thick_roof(nindx,dindx)
       write(iulog,*)'thick_wall:      ',urbinp%thick_wall(nindx,dindx)
       write(iulog,*)'t_building_min:  ',urbinp%t_building_min(nindx,dindx)
       write(iulog,*)'wind_hgt_canyon: ',urbinp%wind_hgt_canyon(nindx,dindx)
       write(iulog,*)'wtlunit_roof:    ',urbinp%wtlunit_roof(nindx,dindx)
       write(iulog,*)'wtroad_perv:     ',urbinp%wtroad_perv(nindx,dindx)
       write(iulog,*)'alb_improad_dir: ',urbinp%alb_improad_dir(nindx,dindx,:)
       write(iulog,*)'alb_improad_dif: ',urbinp%alb_improad_dif(nindx,dindx,:)
       write(iulog,*)'alb_perroad_dir: ',urbinp%alb_perroad_dir(nindx,dindx,:)
       write(iulog,*)'alb_perroad_dif: ',urbinp%alb_perroad_dif(nindx,dindx,:)
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

!-------------------[kz.11]Ray tracing test------------------------- 
  !----------------------------------------------------------------------- 
  subroutine view_factors_v(nzcanm,dzcan,wcan,wbui,&
      tree_cov,lad,lads,ladl,omega,ss_in,pb_in,dray,maxind,maxbhind,n,nsky,hsky,&
      iurb,fww1d,fvv1d,fwv1d,fvw1d,fwr1d,frw1d,fvr1d,frv1d,fwg1d,fgw1d,fgv1d,fsw1d,&
      fvg1d,fsg1d,fsr1d,fsv1d,kww1d,kvv1d,kwv1d,kvw1d,kwr1d,krw1d,kvr1d,krv1d,kwg1d,&
      kgw1d,kgv1d,ksw1d,kvg1d,ksg1d,ksr1d,ksv1d,kws1d,kvs1d,kts1d,krs1d,fws1d,fvs1d,&
      fts1d,frs1d,l)  
    !
    ! !DESCRIPTION: 
    !
    ! Calculate the view factors between all 'surfaces' including vegetation layers.
    ! 
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    integer         , intent(in) :: nzcanm           ! Maximum number of vertical levels at urban resolution 
    real(r8)        , intent(in) :: dzcan            ! Vertical urban grid resolution [m]
    real(r8)        , intent(in) :: wcan             ! Width of the canyons [m]
    real(r8)        , intent(in) :: wbui             ! Width of the buildings [m]
    real(r8)        , intent(in) :: tree_cov(11)     ! Width of the buildings [m]
    real(r8)        , intent(inout) :: lad(nzcanm)      ! Leaf area density in the canyon column [m-1]
    real(r8)        , intent(inout) :: lads(nzcanm)     ! Leaf area density in the canyon column [m-1]
    real(r8)        , intent(inout) :: ladl(nzcanm)     ! Leaf area density in the canyon column [m-1]
    real(r8)        , intent(in) :: omega(nzcanm)    ! Leaf clumping index      
    real(r8)        , intent(in) :: ss_in(nzcanm)    ! Probability to have a building of a given height
    real(r8)        , intent(in) :: pb_in(nzcanm)    ! Probability to have a building of a given height
    real(r8)        , intent(inout) :: dray             ! Ray step
    real(r8)        , intent(in) :: hsky             ! Ray step
    integer         , intent(in) :: maxind           ! indices
    integer         , intent(in) :: maxbhind         ! indices

    integer         , intent(in) :: l                ! indices
    integer         , intent(in) :: n                ! indices (the number of rays:nrays)
    integer         , intent(in) :: nsky             ! indices
    integer         , intent(in) :: iurb             ! indices
    
    ! Area-weighted longwave view factors
    real(r8), intent(out) :: fww1d(nzcanm,nzcanm)        ! view factor from wall to wall for a given canyon
    real(r8), intent(out) :: fvv1d(nzcanm,nzcanm)        ! view factor from vegetation to vegetation for a given canyon
    real(r8), intent(out) :: fwv1d(nzcanm,nzcanm)        ! view factor from wall to vegetation for a given canyon
    real(r8), intent(out) :: fvw1d(nzcanm,nzcanm)        ! view factor from vegetation to wall for a given canyon
    real(r8), intent(out) :: fwr1d(nzcanm,nzcanm)        ! view factor from wall to roof for a given canyon
    real(r8), intent(out) :: frw1d(nzcanm,nzcanm)        ! view factor from roof to wall for a given canyon
    real(r8), intent(out) :: fvr1d(nzcanm,nzcanm)        ! view factor from vegetation to roof for a given canyon
    real(r8), intent(out) :: frv1d(nzcanm,nzcanm)        ! view factor from roof to vegetation for a given canyon
    real(r8), intent(out) :: fwg1d(nzcanm)               ! view factor from wall to ground for a given canyon
    real(r8), intent(out) :: fgw1d(nzcanm)               ! view factor from ground to wall for a given canyon
    real(r8), intent(out) :: fgv1d(nzcanm)               ! view factor from ground to vegetation for a given canyon
    real(r8), intent(out) :: fsw1d(nzcanm)               ! view factor from sky to wall for a given canyon
    real(r8), intent(out) :: fvg1d(nzcanm)               ! view factor from vegetation to ground for a given canyon
    real(r8), intent(out) :: fsg1d                       ! view factor from sky to ground for a given canyon
    real(r8), intent(out) :: fsr1d(nzcanm)               ! view factor from sky to roof for a given canyon
    real(r8), intent(out) :: fsv1d(nzcanm)               ! view factor from sky to vegetation for a given canyon
    real(r8), intent(out) :: fws1d(nzcanm)               ! longwave view factor from wall to sky
    real(r8), intent(out) :: fvs1d(nzcanm)               ! longwave view factor from vegetation to sky
    real(r8), intent(out) :: fts1d                       ! longwave view factor from ground to sky
    real(r8), intent(out) :: frs1d(nzcanm)               ! longwave view factor from roof to sky
    
    ! Area-weighted shortwave view factors
    real(r8), intent(out) :: kww1d(nzcanm,nzcanm)        ! view factor from wall to wall for a given canyon
    real(r8), intent(out) :: kvv1d(nzcanm,nzcanm)        ! view factor from vegetation to vegetation for a given canyon
    real(r8), intent(out) :: kwv1d(nzcanm,nzcanm)        ! view factor from wall to vegetation for a given canyon
    real(r8), intent(out) :: kvw1d(nzcanm,nzcanm)        ! view factor from vegetation to wall for a given canyon
    real(r8), intent(out) :: kwr1d(nzcanm,nzcanm)        ! view factor from wall to roof for a given canyon
    real(r8), intent(out) :: krw1d(nzcanm,nzcanm)        ! view factor from roof to wall for a given canyon
    real(r8), intent(out) :: kvr1d(nzcanm,nzcanm)        ! view factor from vegetation to roof for a given canyon
    real(r8), intent(out) :: krv1d(nzcanm,nzcanm)        ! view factor from roof to vegetation for a given canyon
    real(r8), intent(out) :: kwg1d(nzcanm)               ! view factor from wall to ground for a given canyon
    real(r8), intent(out) :: kgw1d(nzcanm)               ! view factor from ground to wall for a given canyon
    real(r8), intent(out) :: kgv1d(nzcanm)               ! view factor from ground to vegetation for a given canyon
    real(r8), intent(out) :: ksw1d(nzcanm)               ! view factor from sky to wall for a given canyon
    real(r8), intent(out) :: kvg1d(nzcanm)               ! view factor from vegetation to ground for a given canyon
    real(r8), intent(out) :: ksg1d                       ! view factor from sky to ground for a given canyon
    real(r8), intent(out) :: ksr1d(nzcanm)               ! view factor from sky to roof for a given canyon
    real(r8), intent(out) :: ksv1d(nzcanm)               ! view factor from sky to vegetation for a given canyon
    real(r8), intent(out) :: kws1d(nzcanm)               ! shortwave view factor from wall to sky
    real(r8), intent(out) :: kvs1d(nzcanm)               ! shortwave view factor from vegetation to sky
    real(r8), intent(out) :: kts1d                       ! shortwave view factor from ground to sky
    real(r8), intent(out) :: krs1d(nzcanm)               ! shortwave view factor from roof to sky
    
    !
    ! !LOCAL VARIABLES:    
    real(r8), parameter :: pi = 3.1415926535897932384626433832795_r8
    integer             :: nzcan                         ! Number of levels in the urban grid
    real(r8)            :: lad_tree(nzcanm)              ! sum of view factors for wall or road
    real(r8)            :: lsave(nzcanm)                 ! Leaf area density in the canyon column [m-1]                           
    real(r8)            :: ss(nzcanm)                   ! Probability to have a building of a given height
    ! confirm: 
    real(r8)            :: pb(nzcanm)                   ! Probability to have a building
    integer             :: k, izcan, jzcan, kzcan, count, count2
    real(r8)            :: A_s                          ! Area of the street surface
    real(r8)            :: A_g                          ! Area of the ground
    real(r8)            :: A_w(nzcanm)                  ! Area of the wall
    real(r8)            :: A_v(nzcanm)                  ! Area of the vegetation
    real(r8)            :: A_r(nzcanm)                  ! Area of the roof
    real(r8)            :: A_s_max                      ! Maximum area of the street
    real(r8)            :: A_g_max                      ! Maximum area of the ground
    real(r8)            :: A_w_max(nzcanm)              ! Maximum area of the wall
    real(r8)            :: A_v_max(nzcanm)              ! Maximum area of the vegetation
    real(r8)            :: A_r_max(nzcanm)              ! Maximum area of the roof
    real(r8)            :: h, theta, phi, phi1          ! Angles and heights
    real(r8)            :: x, y, z2                     ! Cartesian coordinates
    real(r8)            :: aztmp                        ! Temporary azimuth angle
    real(r8)            :: xx(n), zz(n)                 ! Arrays for ray tracing
    real(r8)            :: circ32(n), rnum, rand        ! Circular and random numbers
    real(r8)            :: rzen, raz                    ! Zenith and azimuth angles
    real(r8)            :: xr, z2r                      ! Adjusted Cartesian coordinates
    real(r8)            :: xxhr(n), zzhr(n)             ! Hemispherical ray tracing coordinates
    real(r8)            :: hemi32r(n)                   ! Hemispherical ray weights
    real(r8)            :: xxhe(n), zzhe(n)             ! Even distribution for hemisphere
    real(r8)            :: hemi32e(n), circ32e(n)       ! Even distribution weights
    real(r8)            :: xxe(n), zze(n)               ! Even distribution coordinates
    real(r8)            :: xxe2(n), zze2(n)             ! Secondary even distribution
    real(r8)            :: circ32e2(n), circ32es        ! Secondary circular weights
    real(r8)            :: zzes                         ! Adjusted z-direction component
    ! confirm: integer(4)/integer*4= integer?
    integer             :: timeArray(3)                 ! Array for timekeeping
    integer             :: kk                           ! Counter for iterations
    ! confirm: !TJ to stop occasional hangups --> any influence?
    real(r8)            :: rayy                         ! Ray y-coordinate
    real(r8)            :: rayx, xfr                    ! Ray x-coordinate and fractional position
    real(r8)            :: xdom                         ! Domain width
    real(r8)            :: bldfrac                      ! Building fraction
    real(r8)            :: vfwt(nzcanm)                 ! View factor from wall to top
    real(r8)            :: pb_old                       ! Previous building probability
    real(r8)            :: vfww(nzcanm, nzcanm)         ! View factor from wall to wall
    real(r8)            :: vfwr(nzcanm, nzcanm)         ! View factor from wall to roof
    real(r8)            :: vfwv(nzcanm, nzcanm)         ! View factor from wall to vegetation
    real(r8)            :: vfvw(nzcanm, nzcanm)         ! View factor from vegetation to wall
    real(r8)            :: raystr                       ! Ray strength
    real(r8)            :: vfrw(nzcanm, nzcanm)         ! View factor from roof to wall
    real(r8)            :: vfrv(nzcanm, nzcanm)         ! View factor from roof to vegetation
    real(r8)            :: bld, bld_old                 ! Building presence indicators
    real(r8)            :: vfvr(nzcanm, nzcanm)         ! View factor from vegetation to roof
    real(r8)            :: vfvv(nzcanm, nzcanm)         ! View factor from vegetation to vegetation
    real(r8)            :: vftw(nzcanm)                 ! View factor from top to wall
    real(r8)            :: vftv(nzcanm)                 ! View factor from top to vegetation
    real(r8)            :: vfvt(nzcanm)                 ! View factor from vegetation to top
    real(r8)            :: maxtree, maxbh               ! Maximum tree height and building height
    real(r8)            :: minray                       ! Minimum ray strength
    real(r8)            :: wfact, rfact                 ! Weighting and reflection factors
    real(r8)            :: raystrtmp                    ! Temporary ray strength
    real(r8)            :: pbinc                        ! Probability increment
    real(r8)            :: wtot                         ! Total width
    real(r8)            :: svfw(nzcanm), svfr(nzcanm)   ! Summed view factors for walls and roofs
    real(r8)            :: svfv(nzcanm), svft           ! Summed view factors for vegetation and top
    real(r8)            :: vfb(nzcanm), vfv(nzcanm)     ! View factors for vegetation
    real(r8)            :: vft_tot                      ! Total view factor
    real(r8)            :: vfw_tot(nzcanm)              ! Total wall view factors
    real(r8)            :: vfr_tot(nzcanm)              ! Total roof view factors
    real(r8)            :: vfv_tot(nzcanm)              ! Total vegetation view factors
    real(r8)            :: vftw_reci(nzcanm)            ! Reciprocal wall view factors
    real(r8)            :: vfvw_reci(nzcanm, nzcanm)    ! Reciprocal vegetation to wall factors
    real(r8)            :: vfvt_reci(nzcanm)            ! Reciprocal vegetation to top factors
    real(r8)            :: xxtmp, zztmp                 ! Temporary Cartesian coordinates
    real(r8)            :: ang, incr                    ! Angle and increment
    real(r8)            :: vfrw_reci(nzcanm, nzcanm)    ! Reciprocal roof to wall factors
    real(r8)            :: xxrt, xoz, xozmax            ! Ray tracing factors
    real(r8)            :: kbs_vf, rnum2                ! Vegetation factor and random number
    real(r8)            :: vfsw(nzcanm), vfsr(nzcanm)   ! View factors from sky
    real(r8)            :: vfsv(nzcanm)                 ! Sky to vegetation factors
    real(r8)            :: vfst                         ! Total sky factor
    real(r8)            :: xxt                          ! Temporary ray factor
    real(r8)            :: vfs_tot, vft                 ! Total view factors
    real(r8)            :: vfr(nzcanm), vfw(nzcanm)     ! View factors for roofs and walls
    real(r8)            :: vfswtmp(nzcanm, nzcanm)      ! Temporary sky to wall factors
    real(r8)            :: vfsrtmp(nzcanm, nzcanm)      ! Temporary sky to roof factors
    real(r8)            :: vfsttmp(nzcanm)              ! Temporary sky total factors
    real(r8)            :: vfsvtmp(nzcanm, nzcanm)      ! Temporary sky to vegetation factors
    real(r8)            :: rsky                         ! Sky radiation factor
    integer             :: rayyint, nk, maxtrind, j, i  ! Indexes and counters
    integer             :: rayyint2, sse(nzcanm)        ! Secondary ray integer and vegetation
    integer             :: ksave, nsrays, ii, jj, nj    ! Saved indexes and counters
    logical             :: switch, switch2              ! Logical switches
    logical             :: roof, strike, horiz          ! Logical states for ray tracing
    logical             :: vartest, solar               ! Logical flags for tests and solar
    real(r8)            :: vfav(nzcanm, nzcanm)         ! View factors for vegetation
    real(r8)            :: vfva(nzcanm, nzcanm)         ! Vegetation view adjustments
    real(r8)            :: vfwa(nzcanm, nzcanm)         ! Wall view adjustments
    real(r8)            :: vfaw(nzcanm, nzcanm)         ! Wall to vegetation factors
    real(r8)            :: vfra(nzcanm, nzcanm)         ! Roof adjustments
    real(r8)            :: vfaa(nzcanm, nzcanm)         ! All-area adjustments
    real(r8)            :: vfar(nzcanm, nzcanm)         ! Vegetation to roof adjustments
    real(r8)            :: vfsatmp(nzcanm, nzcanm)      ! Temporary sky to all factors
    real(r8)            :: vfta(nzcanm), vfat(nzcanm)   ! Adjusted top factors
    real(r8)            :: svfa(nzcanm)                 ! Summed vegetation factors
    real(r8)            :: vfa(nzcanm)                  ! Vegetation factors
    real(r8)            :: vfa_tot(nzcanm)              ! Total vegetation factors
    real(r8)            :: vfsa(nzcanm)                 ! Sky factors
    real(r8)            :: atotot(nzcanm)               ! Total adjustment factors
    real(r8)            :: atotot2(nzcanm)              ! Secondary total adjustment factors
    real(r8)            :: fact_fgv, fact_fvg           ! Factors for ground to vegetation
    real(r8)            :: fact_kgv, fact_kvg           ! Factors for vegetation to ground
    real(r8)            :: A_vs_max(nzcanm)             ! Maximum vegetation surface area
    real(r8)            :: A_vl_max(nzcanm)             ! Maximum vegetation layer area
    real(r8)            :: A_vs(nzcanm)                 ! Vegetation surface area
    real(r8)            :: A_vl(nzcanm)                 ! Vegetation layer area    
    integer             :: start_time, end_time, clock_rate ! Timekeeping variables
    real(r8)            :: elapsed_time                 ! Elapsed time
                 
    !----------------------------------------------------------------
    ! Get the clock rate (ticks per second)
    call system_clock(count_rate=clock_rate)
    call system_clock(start_time)
                                
    write(6,*)'calculating view factors...',n
    count=0
    count2=0

    do izcan=1,maxind
       lsave(izcan)=lad(izcan)
       write(6,*)'lad,s,l',lad(izcan),lads(izcan),ladl(izcan)
    enddo

    !Initialisation
    !----------------
    nzcan = nzcanm
    ! question: why do we need intemediate variables pb and ss?
    do k = 1,nzcanm
       pb(k) = pb_in(k)
       ss(k) = ss_in(k)
    end do

    write(6,*)'new pb',pb
    write(6,*)'new ss',ss

    !! shortwave attenuation by vegetation
    !  extinction coefficient
    !  spherical leaf angle distribution
    !  since rays for calculating view factors always travel normal to the plan that we care about
    kbs_vf=1.0_r8/2.0_r8
    !! 2-D ray directions derived from evenly-spaced vectors on a (3-D) sphere
    kk=1
    jj=1
     
    phi1=0.0_r8

    call init_random_seed(1234)
    call RANDOM_NUMBER(rnum)

    ! n is the number of rays, or the number of points on the edge of the sphere

    ! Generate evenly-spaced points on a sphere of radius 1.
    !----------------------------------------------------------------------------
    ! (This code is originally written by Joseph O'Rourke and Min Xu, June 1997,
    !  and was converted to Fortran from C++.)
  	do k=2,n-1    
       h=-1.0_r8+2.0_r8*real(k-1,r8)/real(n-1,r8)
       theta=acos(h)

       if (theta.lt.0._r8.or.theta.gt.pi) then
          write(6,*)'Error'
          stop
       endif

       phi=phi1+3.6_r8/(sqrt(real(n,r8)*(1._r8-h*h))) 
       !phi=amod(phi,2._r8*pi)
       phi = phi - (2.*pi) * floor(phi / (2.*pi))
       phi1=phi

       x=cos(phi)*sin(theta)
       y=sin(phi)*sin(theta)
       ! z2=cos(theta); But z2==h, so:
       z2=h
       !----------------------------------------------------------------------------
       ! flattened into 2-D (remove y-dimension)	and increase all 2-D vector lengths to lie on the unit circle
       xxe(k-1)=x/sqrt(x*x+z2*z2)
       zze(k-1)=z2/sqrt(x*x+z2*z2)
       ! 3-D to 2-D ratio of ray travel distance
       circ32e(k-1)=1._r8/sqrt(x*x+z2*z2)		

       !! ADDING IN preferential reflection over transmission (ratio 3:2) for
       ! solar radiation (due to preference of NIR to reflect rather than
       ! transmit through leaves, and direct solar to arrive from above):
       call RANDOM_NUMBER(rnum)
       if (z2.lt.0._r8.and.rnum.gt.0.8_r8) then
          z2=abs(z2)
       endif
       ! flattened into 2-D (remove y-dimension)       and increase all 2-D vector lengths to lie on the unit circle
       xxe2(k-1)=x/sqrt(x*x+z2*z2)
       zze2(k-1)=z2/sqrt(x*x+z2*z2)
       ! 3-D to 2-D ratio of ray travel distance
       circ32e2(k-1)=1._r8/sqrt(x*x+z2*z2)

       !! Not used:
       ! actual even distribution over a sphere:
189   continue
       call RANDOM_NUMBER(rnum)
       rzen=acos(2._r8*rnum-1._r8)
       call RANDOM_NUMBER(rnum)
       raz=rnum*2._r8*pi
       xxtmp=sin(rzen)*cos(raz)
       zztmp=cos(rzen)
       xx(k-1)=xxtmp/sqrt(xxtmp*xxtmp+zztmp*zztmp)
       zz(k-1)=zztmp/sqrt(xxtmp*xxtmp+zztmp*zztmp)
       ! redo if too close to horizontal (it will take too long for the ray to attenuate)
       if (abs(zz(k-1)/xx(k-1)).lt.0.001_r8) goto 189
       ! question: check this if statement
       ! I rewrote it
       if (abs(sqrt(xx(k-1)*xx(k-1)+zz(k-1)*zz(k-1))-1._r8).gt.0.001_r8) then
          write(6,*)'PROB:k,x,z,xz',k-1,xx(k-1),zz(k-1), &
                   sqrt(xx(k-1)*xx(k-1)+zz(k-1)*zz(k-1))
       end if
       circ32(k-1)=1._r8/sqrt(xxtmp*xxtmp+zztmp*zztmp)
293  continue

       ! hemispherical distribution (for sfcs)
       if (z2.gt.0._r8) then
          ! regular hemispherical distribution
          xxhe(kk)=xxe(k-1)
          zzhe(kk)=zze(k-1)
          hemi32e(kk)=1._r8/sqrt(x*x+z2*z2)
          ! random hemispherical angles with cosine probability density (Lambert's cosine law)
          call RANDOM_NUMBER(rnum)

          !        rzen=asin(rnum)   ! Arcsin creates a random zenith angle weighted by cosine
          !        rzen=rnum*pi/2._r8	 ! Evenly distributed amongst all solid angles (not Lambertian)
          rzen=asin(sqrt(rnum))      ! Actual Lambert cosine law (Siegel and Howell 2002)...
                                     ! equivalent to acos(sqrt(1-rnum)) (Kondo et al. 2001) or...
                                     ! 0.5*acos(1-2*rnum) (Chelle 2006)
          call RANDOM_NUMBER(rnum)
          raz=rnum*2._r8*pi
          xr=sin(rzen)*cos(raz)
          z2r=cos(rzen)
          xxhr(kk)=xr/sqrt(xr*xr+z2r*z2r)
          zzhr(kk)=z2r/sqrt(xr*xr+z2r*z2r)
          hemi32r(kk)=1._r8/sqrt(xr*xr+z2r*z2r)
422  continue
          kk=kk+1
       endif

    enddo

    nk=kk-1
    nj=jj-1

    goto 423
    ! rays evenly distributed around a half-circle (2-D) -- not used
    incr=pi/real(nk,r8)
    ang=incr/2._r8
    do kk=1,nk
       xxhr(kk)=cos(ang)
       zzhr(kk)=sin(ang)
       hemi32r(kk)=1._r8
       ang=ang+incr       
    enddo
423  continue

    write(6,*)'nk,nj',nk,nj
    !-----------------------------------------------------------------------
    
    ! find indices with ss>0
    do izcan=1,maxbhind
       sse(izcan)=ss(izcan)
    enddo

    ! ray strength at which we stop tracking the ray:
    minray=max(0.000001_r8,0.0001_r8*pb(maxbhind))

    ! total view factors (diagnostics to see if they add up to 1.0 for each surface)
    vft_tot=0._r8
    vfs_tot=0._r8
    do izcan=1,nzcanm
       vfr_tot(izcan)=0._r8
       vfw_tot(izcan)=0._r8
       vfv_tot(izcan)=0._r8
    enddo

    write(6,*)'before view factor calcs'
    
    call system_clock(end_time)        
    ! Calculate elapsed time in seconds
    elapsed_time = real(end_time - start_time,r8) / real(clock_rate,r8)
    
    ! write the elapsed time for each iteration  
    write(*, '(A, I0, A, F6.3)') 'Preparing stage elapsed time for l = ', l, ': ', elapsed_time, ' seconds'
        
    !-----------------------------------------------------------------------
    !! VIEW FACTOR CALCULATIONS
    !! calculate view factors for sfc-sfc, veg-sfc, and veg-veg diffuse exchange using ray tracing
    ! only need to find view factors for wall on one side due to symmetry

    ! width of the two columns (canyon and building)
    wtot=wcan+wbui
    write(6,*)'wcan,wbui,wtot,dzcan',wcan,wbui,wtot,dzcan
    bldfrac=wbui/wtot
    xdom=wtot/dzcan
    svft=0._r8

    !      dray=dray_in*min(min(dzcan,wcan(ican)),wbui(ican))/dzcan
    ! hard-code the ray step size for view factors:
    dray=0.05_r8*min(min(dzcan,wcan),wbui)/dzcan

    !       nsrays=max(maxind+1,ceiling(wcan(1)/dzcan))*nsky
    nsrays=(maxind)*nsky
    write(6,*)'nsky,hsky',nsky,hsky
    
    solar=.false.
    goto 359

358   continue

    solar=.true.
    write(6,*)'----------SOLAR VIEW FACTORS----------'
    do izcan=1,maxind
       lad(izcan)=lads(izcan)
    enddo

359   continue

    if (.not.solar) then
       write(6,*)'----------LONGWAVE VIEW FACTORS----------'
       do izcan=1,maxind
          lad(izcan)=ladl(izcan)
       enddo
    endif

    vft_tot=0._r8
    vfs_tot=0._r8
    do izcan=1,nzcanm
       vfr_tot(izcan)=0._r8
    enddo
    do izcan=1,nzcanm
       vfw_tot(izcan)=0._r8
       vfv_tot(izcan)=0._r8
    enddo
    
    ! MAIN IZCAN LOOP:
    do izcan=1,maxind
       vfst=0._r8
       svfw(izcan)=0._r8
       svfr(izcan)=0._r8
       svfv(izcan)=0._r8
       do jzcan=1,maxind
         ! w=walls; r=roofs; v=vegetation (trees); t=streets; s=sky
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

       write(6,*)'izcan=',izcan
       !-----------------------------------------------------
       write(6,*)'SKY...'
       ! rays starting from SKY
       write(6,*)'nsrays=',nsrays,xdom,nk,nsky
       horiz=.true.
       do ii=1,nsky
          write(6,*)'ii(sky)=',ii
          do kk=1,nk
             ! starting at a location at hsky times the highest tree/building, spread out evenly in the x-direction 
             raystr=1._r8
             call RANDOM_NUMBER(rnum)
             rayx=(real((izcan-1)*nsky+ii,r8)-rnum)/real(nsrays,r8)*xdom
             if (hsky.lt.1._r8) then
                write(6,*)'hsky (ray start height for sky diffuse view &
                         factors) must be greater than 1; hsky=',hsky
                stop
             endif
             rayy=max(real(maxind,r8)*hsky,real(maxind,r8)+2._r8*dray)

             ! make all rays head 'downstream' (this works due to the symmetry of the geometry)
             xxt=xxhr(kk)
             ! mirror reflection so that all rays still head in positive direction
             if (xxt.lt.0._r8) then
                xxt=-xxt
                rayx=xdom*(2._r8-bldfrac)-rayx
             endif
             do while (rayy.gt.real(maxind,r8)+2._r8*dray)
                ! COULD MULTIPLY dray BY A FACTOR PROPORTIONAL TO THE DISTANCE ABOVE REAL(MAXIND), AND WITH A MINIMUM VALUE, TO SPEED THIS UP!
                rayx=rayx+xxt*dray
                rayy=rayy-zzhr(kk)*dray
             enddo
             ! test where a ray is horizontally (i.e. building or canyon/tree space):
             ! xfr>bldfrac means canyon, xfr<bldfrac means building (i.e. bld=1)
             ! xfr=amod(rayx,xdom)
             xfr = rayx - xdom* floor(rayx /xdom)
             bld=0._r8
             if (xfr/xdom.ge.1._r8-bldfrac) bld=1._r8
             bld_old=bld
             pb_old=0._r8

             call ray_dn(minray,rayx,rayy,bld,bld_old,raystr, &
                      xxt,-zzhr(kk),hemi32r(kk),xdom,dray,dzcan,bldfrac,pb, &
                      pb_old,ss,kbs_vf,lad,omega,horiz,vfw,vfr,vfv,vft)
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
       !       endif

       !-----------------------------------------------------
       ! rays starting from WALLS
       if (izcan.le.maxbhind-1) then
          write(6,*)'WALLS...'
          horiz=.false.
          do kk=1,nk
             ! starting at the downstream building edge (upstream edge of the canyon)
             bld=0._r8
             bld_old=0._r8
             raystr=1._r8
             rayx=0._r8
             ! random:
             call RANDOM_NUMBER(rnum)
             rayy=real(izcan,r8)-rnum
             ! even:
             !         rayy=real(izcan)-(real(kk)-0.5)/real(nk)
             ! centre:
             !         rayy=real(izcan)-0.5
             if (xxhr(kk).le.0._r8) goto 545
             ! rays going up
             call ray_up(minray,maxind,rayx,rayy,bld,bld_old,raystr, &
                        xxhr(kk),zzhr(kk),hemi32r(kk),xdom,dray,dzcan,&
                        bldfrac,pb,kbs_vf,lad,omega,horiz,vfb,vfv)
             do kzcan=1,nzcanm
                vfww(izcan,kzcan)=vfww(izcan,kzcan)+vfb(kzcan)
                vfwv(izcan,kzcan)=vfwv(izcan,kzcan)+vfv(kzcan)
                !vfwup_total=vfwup_total+vfb(kzcan)+vfv(kzcan)
             enddo

             !nrayup(izcan)=nrayup(izcan)+1

             goto 547
545     continue

             ! rays going down
             pb_old=pb(izcan+1)
             call ray_dn(minray,rayx,rayy,bld,bld_old,raystr,&
                        xxhr(kk),zzhr(kk),hemi32r(kk),xdom,dray,dzcan,&
                        bldfrac,pb,pb_old,ss,kbs_vf,lad,omega,horiz,vfw,vfr,vfv,vft)
             do kzcan=1,nzcanm
                vfww(izcan,kzcan)=vfww(izcan,kzcan)+vfw(kzcan)
                vfwr(izcan,kzcan)=vfwr(izcan,kzcan)+vfr(kzcan)
                vfwv(izcan,kzcan)=vfwv(izcan,kzcan)+vfv(kzcan)
                !vfwdn_total=vfwdn_total+vfw(kzcan)+vfv(kzcan)+vfr(kzcan)
             enddo
             vfwt(izcan)=vfwt(izcan)+vft
             !vfwdn_total=vfwdn_total+vft

             !nraydn(izcan)=nraydn(izcan)+1

            goto 548
547    continue
            !! WALLS-SKY
            svfw(izcan)=svfw(izcan)+raystr
            !vfwup_total=vfwup_total+raystr
548    continue

          enddo  ! end rays (kk) loop for WALLS
          
          svfw(izcan)=svfw(izcan)/real(nk,r8)
          vfwt(izcan)=vfwt(izcan)/real(nk,r8)
          vfw_tot(izcan)=vfw_tot(izcan)+svfw(izcan)+vfwt(izcan)
          !vfwdn_total=vfwdn_total/real(nk,r8)*2._r8
          !vfwup_total=vfwup_total/real(nk,r8)*2._r8
       endif  ! if izcan has a wall layer

       !-----------------------------------------------------
       ! rays starting from VEGETATION in CANOPY COLUMN
       if (lad(izcan).gt.0._r8) then
           write(6,*)'CANOPY VEGETATION...'
           horiz=.true.
           do k=1,n-2
              ! starting at a random location in the vegetation
              bld=0._r8
              bld_old=0._r8
              raystr=1._r8
              call RANDOM_NUMBER(rnum)
              rayx=rnum*xdom*(1._r8-bldfrac)
              ! random:
              call RANDOM_NUMBER(rnum)
              rayy=real(izcan,r8)-rnum
              ! make all rays head 'downstream' (this works due to the symmetry of the geometry)
              !          xxrt=xxe(k)
              ! also use different ray directional distributions for solar and
              ! longwave
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
              if (xxrt.lt.0._r8) then
                 xxrt=-xxrt
                 rayx=xdom*(1._r8-bldfrac)-rayx
              endif
	            lad_tree=lad/tree_cov(iurb)

              !if(zze(k).le.0._r8) goto 645
              if(zzes.le.0._r8) goto 645
              ! rays going up
              call ray_up_foliage(minray,maxind,rayx,rayy,bld,bld_old,raystr,&
                       !xxrt,zze(k),circ32e(k),xdom,dray,dzcan,bldfrac,pb,&
                       xxrt,zzes,circ32es,xdom,dray,dzcan,bldfrac,pb, &
                       kbs_vf,lad_tree,lad,omega,horiz,vfb,vfv)
              do kzcan=1,nzcanm
                 vfvw(izcan,kzcan)=vfvw(izcan,kzcan)+vfb(kzcan)
                 vfvv(izcan,kzcan)=vfvv(izcan,kzcan)+vfv(kzcan)
              enddo
              goto 647

645     continue
              ! rays going down
              pb_old=pb(izcan+1)
              call ray_dn(minray,rayx,rayy,bld,bld_old,raystr,&
                        !xxrt,zze(k),circ32e(k),xdom,dray,dzcan,bldfrac,&
                        xxrt,zzes,circ32es,xdom,dray,dzcan,bldfrac, &
                        pb,pb_old,ss,kbs_vf,lad,omega,horiz,vfw,vfr,vfv,vft)
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

           svfv(izcan)=svfv(izcan)/real(n,r8)
           vfvt(izcan)=vfvt(izcan)/real(n,r8)
           vfv_tot(izcan)=vfv_tot(izcan)+vfvt(izcan)+svfv(izcan)
       endif  ! if izcan has a foliage layer
       
       !write(6,*)'AFTER CANOPY FOLIAGE'

       !-----------------------------------------------------
       ! rays starting from ROOFS
       if (ss(izcan).gt.0._r8) then
          write(6,*)'ROOFS...'
          do kk=1,nk
          ! starting at a random point on the roof
             bld=1._r8
             bld_old=1._r8
             raystr=1._r8
             rayy=real(izcan-1,r8)
             ! random:
             ! start at a random point on the roof
             call RANDOM_NUMBER(rnum)
             rayx=xdom*((1._r8-bldfrac)+rnum*bldfrac)
             ! make all rays head 'downstream' (this works due to the symmetry of the geometry)
             xxrt=xxhr(kk)
             ! mirror reflection so that all rays still head in positive direction
             if (xxrt.lt.0._r8) then
                xxrt=-xxrt
                rayx=xdom-(rayx-xdom*(1._r8-bldfrac))
             endif
             horiz=.true.
             call ray_up(minray,maxind,rayx,rayy,bld,bld_old,raystr,&
                       xxrt,zzhr(kk),hemi32r(kk),xdom,dray,dzcan,&
                       bldfrac,pb,kbs_vf,lad,omega,horiz,vfb,vfv)
             do kzcan=1,nzcanm
                vfrw(izcan,kzcan)=vfrw(izcan,kzcan)+vfb(kzcan)
                vfrv(izcan,kzcan)=vfrv(izcan,kzcan)+vfv(kzcan)
             enddo
             svfr(izcan)=svfr(izcan)+raystr

          enddo  ! end rays (kk) loop for ROOFS
          svfr(izcan)=svfr(izcan)/real(nk,r8)
          vfr_tot(izcan)=vfr_tot(izcan)+svfr(izcan)
          write(6,*)'iz,nk,svfr',izcan,nk,svfr(izcan),vfr_tot(1),vfr_tot(2)

       endif

    !-----------------------------------------------------

    enddo  ! end wall/vegetation/roof level (izcan) loop

    write(6,*)'AFTER IZCAN LOOP'

    !-----------------------------------------------------------------------
    write(6,*)'ROAD...'
    ! rays starting from ROAD
    do kk=1,nk
       ! starting at a random point on the road
       bld=0._r8
       bld_old=0._r8
       raystr=1._r8
       
       rayy=0._r8
       ! random:
       ! start at a random point on the road
       call RANDOM_NUMBER(rnum)
       rayx=xdom*(1._r8-bldfrac)*rnum
       ! make all rays head 'downstream' (this works due to the symmetry of the geometry)
       xxrt=xxhr(kk)
       ! mirror reflection so that all rays still head in positive direction
       if (xxrt.lt.0._r8) then
          xxrt=-xxrt
          rayx=xdom*(1._r8-bldfrac)-rayx
       endif
       horiz=.true.
       call ray_up(minray,maxind,rayx,rayy,bld,bld_old,raystr,&
                 xxrt,zzhr(kk),hemi32r(kk),xdom,dray,dzcan,bldfrac,&
                 pb,kbs_vf,lad,omega,horiz,vfb,vfv)
       do kzcan=1,nzcanm
          vftw(kzcan)=vftw(kzcan)+vfb(kzcan)
          vftv(kzcan)=vftv(kzcan)+vfv(kzcan)
       enddo
       svft=svft+raystr
    enddo  ! end rays (kk) loop for ROAD

    svft=svft/real(nk,r8)

    vft_tot=vft_tot+svft
    do izcan=1,nzcanm
       !divide by 2 here because there are two walls at each level, and we want the vf from the road to only one of them
       vftw(izcan)=vftw(izcan)/real(nk,r8)/2._r8
       vftv(izcan)=vftv(izcan)/real(nk,r8)
    enddo   
    !-----------------------------------------------------------------------
    do izcan=1,maxind
       write(6,*)'i,svfw,svfr,svfv',izcan,svfw(izcan),svfr(izcan),svfv(izcan)
    enddo
    write(6,*)'svft',svft

    vfst=0._r8
    do izcan=1,maxind
       vfsw(izcan)=0._r8
       vfsr(izcan)=0._r8
       vfsv(izcan)=0._r8
       do jzcan=1,maxind
          vfww(izcan,jzcan)=vfww(izcan,jzcan)/real(nk,r8)
          vfvv(izcan,jzcan)=vfvv(izcan,jzcan)/real(n,r8)
          vfwv(izcan,jzcan)=vfwv(izcan,jzcan)/real(nk,r8)
          vfwr(izcan,jzcan)=vfwr(izcan,jzcan)/real(nk,r8)
          vfw_tot(izcan)=vfw_tot(izcan)+vfww(izcan,jzcan)+ &
                         vfwr(izcan,jzcan)+vfwv(izcan,jzcan)
          ! divide by 2 here because there are two walls at each level, and we want the vf from the roof to only one of them
          vfrw(izcan,jzcan)=vfrw(izcan,jzcan)/real(nk,r8)/2._r8
          vfrv(izcan,jzcan)=vfrv(izcan,jzcan)/real(nk,r8)
          vfr_tot(izcan)=vfr_tot(izcan)+vfrw(izcan,jzcan)*2._r8+vfrv(izcan,jzcan)
          vfvr(izcan,jzcan)=vfvr(izcan,jzcan)/real(n,r8)
          vfvw(izcan,jzcan)=vfvw(izcan,jzcan)/real(n,r8)/2._r8

          vfrw_reci(jzcan,izcan)=vfwr(izcan,jzcan)*pb(izcan+1)/&
                          (max(ss(jzcan),1.e-6_r8)*xdom*bldfrac)
       enddo
    enddo

    ! Print out view factor values:
    do izcan=1,maxind
       do jzcan=1,maxind	 
          write(6,*)'i,j,vv',izcan,jzcan,vfvv(izcan,jzcan)            
       enddo
    enddo

    vfst=0._r8
    do izcan=1,maxind
       vfst=vfst+vfsttmp(izcan)
       do jzcan=1,maxind
          vfsw(izcan)=vfsw(izcan)+vfswtmp(jzcan,izcan)
          vfsr(izcan)=vfsr(izcan)+vfsrtmp(jzcan,izcan)
          vfsv(izcan)=vfsv(izcan)+vfsvtmp(jzcan,izcan)	 
          write(6,*)'i,j,wr,wv',izcan,jzcan,vfwr(izcan,jzcan),vfwv(izcan,jzcan)            
       enddo
    enddo

    izcan=1
    jzcan=1 
    write(6,*)'i,j,ww',izcan,jzcan,vfww(izcan,jzcan)            

    do izcan=1,maxind
       do jzcan=1,maxind	 
          write(6,*)'i,j,rv,vr',izcan,jzcan,vfrv(izcan,jzcan),vfvr(izcan,jzcan)            
       enddo
    enddo

    do izcan=1,maxind
       jzcan=1
       write(6,*)'i,j,vw',izcan,jzcan,vfvw(izcan,jzcan)	 
       do jzcan=1,maxind	 
          vfv_tot(izcan)=vfv_tot(izcan)+2._r8*vfvw(izcan,jzcan)+&
          vfvr(izcan,jzcan)+vfvv(izcan,jzcan)         
       enddo
    enddo

    do izcan=1,maxind
       jzcan=1 
       write(6,*)'i,j,rw,rwr',izcan,jzcan,vfrw(izcan,jzcan),vfrw_reci(izcan,jzcan)            
    enddo

    izcan=1
    vftw_reci(izcan)=vfwt(izcan)*pb(izcan+1)/(xdom*(1._r8-bldfrac)) 
    write(6,*)'i,wt,tw,twr',izcan,vfwt(izcan),vftw(izcan),vftw_reci(izcan)

    do izcan=1,maxind
       write(6,*)'i,tv,vt',izcan,vftv(izcan),vfvt(izcan)                                     
    enddo

    do izcan=1,maxind
       vfsw(izcan)=vfsw(izcan)/real(maxind,r8)/2._r8
       vfsr(izcan)=vfsr(izcan)/real(maxind,r8)
       vfsv(izcan)=vfsv(izcan)/real(maxind,r8)
       vfs_tot=vfs_tot+2._r8*vfsw(izcan)+vfsr(izcan)+vfsv(izcan)
       write(6,*)'i,sw,sr',izcan,vfsw(izcan),vfsr(izcan)
    enddo
    do izcan=1,maxind
       write(6,*)'i,sv',izcan,vfsv(izcan)
    enddo
    vfst=vfst/real(maxind,r8)
    vfs_tot=vfs_tot+vfst
    write(6,*)'st',vfst

    vft_tot=svft
    do jzcan=1,maxind	             
       vft_tot=vft_tot+2._r8*vftw(jzcan)+vftv(jzcan)           
    enddo
    write(6,*)'vft,vfs',vft_tot,vfs_tot
    do izcan=1,maxind            	 
       write(6,*)'i,vfw,vfr',izcan,vfw_tot(izcan),vfr_tot(izcan)                                     
    enddo
   do izcan=1,maxind            	 
      write(6,*)'i,vfv',izcan,vfv_tot(izcan)                                    
   enddo
   
   !Calculate the modified radiation and the streets fluxes 
   !-------------------------------------------------------
   ! areas of all elements (normalized by dzcan, i.e., xdom = domain width/dzcan)
   do izcan=1,maxind
      A_w(izcan)=pb(izcan+1)
      ! TO GET ACTUAL RADIATION FLUX DENSITIES ON LEAVES, NEED TO MULTIPLY THEM BY OMEGA
      A_v(izcan)=xdom*(1._r8-bldfrac)*lsave(izcan)*omega(izcan)*dzcan*2._r8
      !A_v(izcan)=xdom*(1._r8-bldfrac)*max(lads(izcan),ladl(izcan))*
      !&                      omega(izcan)*dzcan*2._r8
      A_vs(izcan)=xdom*(1._r8-bldfrac)*lads(izcan)*omega(izcan)*dzcan*2._r8
      A_vl(izcan)=xdom*(1._r8-bldfrac)*ladl(izcan)*omega(izcan)*dzcan*2._r8

      A_r(izcan)=xdom*bldfrac*ss(izcan)
      A_w_max(izcan)=max(1.e-6_r8,A_w(izcan))
      A_v_max(izcan)=max(1.e-6_r8,A_v(izcan))
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
         end do !izcan
      end do !jzcan

      do izcan=1,maxind                             
         kwg1d(izcan)=vfwt(izcan)*A_w(izcan)/A_g
         kgw1d(izcan)=vftw(izcan)*A_g/A_w_max(izcan)
         kgv1d(izcan)=vftv(izcan)*A_g/A_vs_max(izcan)
         kvg1d(izcan)=vfvt(izcan)*A_vs(izcan)/A_g
      enddo !izcan

      do izcan=1,maxind
         ksw1d(izcan)=vfsw(izcan)*A_s/A_w_max(izcan)
         ksr1d(izcan)=vfsr(izcan)*A_s/A_r_max(izcan)
         ksv1d(izcan)=vfsv(izcan)*A_s/A_vs_max(izcan)

         kws1d(izcan)=svfw(izcan)
         krs1d(izcan)=svfr(izcan)
         kvs1d(izcan)=svfv(izcan)
      end do
      do izcan=1,maxind
   enddo
   ksg1d=vfst*A_s/A_g
   kts1d=svft

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
         end do !izcan
      end do !jzcan
      
      do izcan=1,maxind                             
         fwg1d(izcan)=vfwt(izcan)*A_w(izcan)/A_g
         fgw1d(izcan)=vftw(izcan)*A_g/A_w_max(izcan)
         fgv1d(izcan)=vftv(izcan)*A_g/A_vl_max(izcan)
         fvg1d(izcan)=vfvt(izcan)*A_vl(izcan)/A_g
      end do !izcan
     
      write(6,*)'printing a fwg1d value...'
      write(6,*)'dimension 1',fwg1d(1)
      !stop

      do izcan=1,maxind

         fsw1d(izcan)=vfsw(izcan)*A_s/A_w_max(izcan)
         fsr1d(izcan)=vfsr(izcan)*A_s/A_r_max(izcan)
         fsv1d(izcan)=vfsv(izcan)*A_s/A_vl_max(izcan)
         fws1d(izcan)=svfw(izcan)
         frs1d(izcan)=svfr(izcan)
         fvs1d(izcan)=svfv(izcan)

         goto 237
         fsw1d(izcan)=svfw(izcan)
         fsr1d(izcan)=svfr(izcan)
         fsv1d(izcan)=svfv(izcan)
237  continue
      enddo

      fts1d=svft
      fsg1d=vfst*A_s/A_g
   endif

   call system_clock(end_time)        
   ! Calculate elapsed time in seconds
   elapsed_time = real(end_time - start_time,r8) / real(clock_rate,r8)
  
   if (solar) then !print calulation time after solar calculation finishs
      ! write the elapsed time for each iteration  
      write(*, '(A, I0, A, F6.3)') 'Calculation elapsed time for solar calculation l = ', l, ': ', elapsed_time, ' seconds'           
   else ! print longwave 
      ! write the elapsed time for each iteration  
      write(*, '(A, I0, A, F6.3)') 'Calculation elapsed time for longwave calculation l = ', l, ': ', elapsed_time, ' seconds'           
   endif           

       
   if (solar) goto 348

   !-----------------------------------------------------------------------
   ! Now go back and compute view factors for solar!!
   
   goto 358
   
348  continue
 

   do izcan=1,maxind
      do jzcan=1,maxind
      enddo
   enddo

   do izcan=1,maxind
      ! debug
      !lad(izcan)=lsave(izcan)
   enddo
   
   call system_clock(end_time)        
   ! Calculate elapsed time in seconds
   elapsed_time = real(end_time - start_time,r8) / real(clock_rate,r8)
    
   ! write the elapsed time for each iteration  
   write(*, '(A, I0, A, F6.3)') 'Finish elapsed time for l = ', l, ': ', elapsed_time, ' seconds'
   
  end subroutine view_factors_v 

  !-----------------------------------------------------------------------
  subroutine corner_up(bldfrac,xdom,rayx,rayy,rayyint,zre,strike)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! This routine figures out whether a ray passing a corner at the upstream
    ! edge of a building strikes the roof and the wall of the layer above, or
    ! just the current wall layer     !
    ! ----------------------------------------------------------------------
    ! !ARGUMENTS:
    implicit none
    integer         , intent(in)  :: rayyint          ! 
    real(r8)        , intent(in)  :: bldfrac          !
    real(r8)        , intent(in)  :: xdom             !
    real(r8)        , intent(in)  :: rayx             !
    real(r8)        , intent(in)  :: zre              !
    real(r8)        , intent(in)  :: rayy             !
    logical         , intent(out) :: strike           ! 

    strike=.false.
    ! if the following is true, then the ray passing an upstream corner impinged on the roof
    ! (and the wall layer above the current wall layer) over the last ray step
    !if (atan((amod(rayx,xdom)-xdom*(1._r8-bldfrac))/ &
    if (atan(((rayx - xdom * floor(rayx /xdom))-xdom*(1._r8-bldfrac))/ &
                  max(1.e-6_r8,(real(rayyint,r8)-rayy))).gt.zre) then  
       strike=.true.
    endif
  end subroutine corner_up

  
  !-----------------------------------------------------------------------
  subroutine corner_dn(xdom,rayx,rayy,rayyint,zre,strike)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! This routine figures out whether a ray passing a corner at the downstream
    ! edge of a building strikes the roof or not
    ! ----------------------------------------------------------------------
    ! !ARGUMENTS:
    implicit none
    integer         , intent(in)  :: rayyint          ! 
    real(r8)        , intent(in)  :: xdom             !
    real(r8)        , intent(in)  :: rayx             !
    real(r8)        , intent(in)  :: zre              !
    real(r8)        , intent(in)  :: rayy             !
    logical         , intent(out) :: strike           ! 

    strike=.false.
    ! if the following is true, then the ray passing a downstream corner impinged on the roof
    ! over the last ray step
    !if (atan(amod(rayx,xdom)/max(1.e-6_r8,(real(rayyint,r8)-rayy))).lt.zre) then 
    if (atan((rayx - xdom * floor(rayx /xdom))/ &
                 max(1.e-6_r8,(real(rayyint)-rayy))).lt.zre) then     
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
    integer         , intent(in)  :: random_nl         ! 
    ! !LOCAL VARIABLES:
    integer                       :: i, n, clock
    integer, dimension(:), allocatable :: seed
      
    ! Get the size of the seed
    call RANDOM_SEED(size = n)

    ! Allocate memory for the seed array
    allocate(seed(n))
      
    ! Get the current time from the system clock
    call SYSTEM_CLOCK(count = clock)
    ! Compute the seed values based on system clock
    !seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    ! Compute the seed values based on a user-defined integer
    seed = random_nl + 37 * (/ (i - 1, i = 1, n) /)
    
    ! Set the random seed
    call RANDOM_SEED(put = seed)
      
    deallocate(seed)
  end subroutine init_random_seed

  !-----------------------------------------------------------------------
  subroutine ray_up(minray, maxind, rayx, rayy, bld, bld_old, raystr, &
                    xx, zz, hemi32, xdom, dray, dzcan, bldfrac, pb, &
                    kbs, lad, omega, horiz, vfb, vfv)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Computes the contribution to view factors of rays travelling "upwards"
    ! with a positive vertical component.
    !-----------------------------------------------------------------------
    !!note: define a constant variable layer numbers
    ! !ARGUMENTS:
    implicit none
    integer,  parameter  :: nzcanm = 4       ! Maximum number of vertical levels at urban resolution
    real(r8), intent(in) :: minray
    real(r8), intent(inout) :: rayx, rayy, bld, bld_old, raystr
    integer,  intent(in) :: maxind
    real(r8), intent(in) :: xx, zz, hemi32, xdom, dray, dzcan, bldfrac
    real(r8), intent(in) :: kbs
    real(r8), intent(in) :: omega(nzcanm), lad(nzcanm)
    real(r8), intent(in) :: pb(nzcanm+1)
    logical,  intent(in) :: horiz

    real(r8), intent(out) :: vfb(nzcanm), vfv(nzcanm)
    !
    ! !LOCAL VARIABLES:
    integer              :: rayyint, izcan, rayyint_old, rayyint2
    real(r8)             :: rx, ry, xfr, wfact, raystrtmp, raystra
    logical              :: strike
    !-----------------------------------------------------------------------
    !        raystra=raystr
    ! rays originating from walls
    rx=zz
    ry=xx
    ! rays originating from roads/vegetation
    if (horiz) then
       rx=xx
       ry=zz
    endif

    do izcan=1,nzcanm
       vfb(izcan)=0._r8
       vfv(izcan)=0._r8
    enddo

    do while (raystr.gt.minray)
       rayyint_old=ceiling(rayy)
       rayx=rayx+rx*dray
       rayy=rayy+ry*dray
       rayyint=ceiling(rayy)
       !xfr=amod(rayx,xdom)
       xfr = (rayx - xdom * floor(rayx /xdom))
       
       if (xfr/xdom.ge.1._r8-bldfrac) then
          ! interception by walls
          bld=1._r8
          ! if(bld.gt.bld_old) raystra=raystr
          rayyint2=rayyint
          if (rayyint.gt.rayyint_old) then
             ! the inputs to corner_up are slightly different than in ray_dn
             call corner_up(bldfrac,xdom,rayx,-rayy, -rayyint_old,&
                           atan(rx/max(1.e-6_r8,ry)),strike)
             if (strike) then
                ! lower wall layer is being lit
  		          rayyint2=rayyint_old
             endif
          endif
          wfact=raystr*max(bld-bld_old,0._r8)*pb(rayyint2+1)/max(1.e-6_r8,pb(rayyint2+1))*pb(rayyint2+1)
          vfb(rayyint2)=vfb(rayyint2)+wfact
          raystr=raystr-wfact
       else
          bld=0._r8
          ! interception by vegetation in the canopy column
          raystrtmp=raystr
          raystr=raystr*exp(-hemi32*dray*dzcan*(kbs*lad(rayyint)*omega(rayyint)))
          vfv(rayyint)=vfv(rayyint)+(raystrtmp-raystr)
        endif
        if(rayyint.gt.maxind+1) goto 541
        bld_old=bld
    enddo

541   continue

  end subroutine ray_up
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine ray_dn(minray, rayx, rayy, bld, bld_old, raystr, xx, zz, &
                    dist32, xdom, dray, dzcan, bldfrac, pb, pb_old, ss, &
                    kbs, lad, omega, horiz, vfw, vfr, vfv, vft)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Computes the contribution to view factors of rays travelling "downwards"
    ! with a negative vertical component.
    !-----------------------------------------------------------------------
    ! !note: run the original code (turn on the debug mode)
    ! !ARGUMENTS:
    implicit none
    integer,  parameter  :: nzcanm = 4       ! Maximum number of vertical levels at urban resolution
    real(r8), intent(in) :: minray
    real(r8), intent(inout) :: rayx, rayy, bld, bld_old, raystr, pb_old
    real(r8), intent(in) :: xx, zz, dist32, xdom, dray, dzcan, bldfrac
    real(r8), intent(in) ::  pb(nzcanm + 1), ss(nzcanm + 1)
    real(r8), intent(in) :: kbs, lad(nzcanm), omega(nzcanm)
    logical,  intent(in)  :: horiz

    real(r8), intent(out) :: vfw(nzcanm), vfr(nzcanm), vfv(nzcanm), vft

    !
    ! !LOCAL VARIABLES:
    integer              :: rayyint, izcan, kzcan, rayyint2
    real(r8)             :: rx, ry, xfr, wfact, raystrtmp, pbinc
    real(r8)             :: rfact, blde, sseff(nzcanm + 1),raystra
    logical              :: roof, strike
    !-----------------------------------------------------------------------
    !         raystra=raystr
    !         raystrb=raystr

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
    ! rays going down
    rayyint=ceiling(rayy)

    do izcan=1,nzcanm
       if (pb(izcan+1).lt.1._r8) then
          sseff(izcan)=ss(izcan)/(1._r8-pb(izcan+1))
       else
          sseff(izcan)=ss(izcan)
       endif
    enddo

    do while (raystr.gt.minray)
       rayx=rayx+rx*dray
       rayy=rayy+ry*dray
       rayyint=ceiling(rayy)
       rayyint2=rayyint
       pbinc=pb(rayyint+1)-pb_old
       ! test where a ray is horizontally (i.e. building or canyon/tree space):
       ! xfr>bldfrac means canyon, xfr<bldfrac means building (i.e. bld=1)
       !xfr=amod(rayx,xdom)
       xfr = (rayx - xdom * floor(rayx /xdom))
       bld=0._r8
       if(xfr/xdom.ge.1._r8-bldfrac) bld=1._r8
       blde=bld

       roof=.true.
 
       ! If the current ray step involves a corner (i.e., crossing both a pb boundary and the building canyon boundary)
       ! if we are crossing a column:
       if (abs(bld-bld_old).gt.0.5_r8) then
          if (bld.lt.bld_old) then 
             if (pbinc.gt.0._r8) then
                ! building corner at 'downstream' edge of the building
                call corner_dn(xdom,rayx,rayy,rayyint, atan(rx/max(1.e-6_r8,(-ry)))&
                               ,strike)
                if (strike) then
                    blde=1._r8
                    goto 223
                else
                    goto 224
                endif
             endif
           ! if we are entering the building column (bld.gt.bld_old)
          else
             if (ceiling(rayy-ry*dray)-rayyint.gt.0) then
                 ! if(rayy.le.0.) goto 225
                 ! building corner or wall layer division at 'upstream' edge of the building
                 call corner_up(bldfrac,xdom,rayx,rayy,rayyint,&
                               atan(rx/max(1.e-6,(-ry))),strike)
                 if (strike) then
                    ! next wall layer higher is being lit (and roofs too)
  		               rayyint2=rayyint+1
                 else
                    ! only the present wall layer is being lit (not roofs)
                    roof=.false.
                 endif
             endif		
          endif          
       endif

       if(rayyint2.lt.1) goto 223
       !! WALLS
       wfact=raystr*max(bld-bld_old,0._r8)*pb(rayyint2+1)/max(1.e-6_r8,pb(rayyint2+1))*pb(rayyint2+1)
       vfw(rayyint2)=vfw(rayyint2)+wfact
       raystr=raystr-wfact

223      continue
       if(.not.roof) goto 224
       !! ROOFS
       rfact=blde*raystr*pbinc/max(1.e-6_r8,pbinc)*sseff(rayyint+1)
       ! important, otherwise could end up with negative ray strength!
       !          rfact=min(rfact,raystr)                           
       vfr(rayyint+1)=vfr(rayyint+1)+rfact                 
       raystr=raystr-rfact
224      continue
       if(rayy.le.0._r8) goto 226
       !! VEGETATION
       ! so far I have not added in the details, e.g. what if a ray crosses vegetation layers or from or into a
       ! building during the ray step? As long as dray is quite small this should not be too important
       if (bld.lt.0.5_r8) then
       ! interception by vegetation in the canopy column
          raystrtmp=raystr
          raystr=raystr*exp(-dist32*dray*dzcan*(kbs*lad(rayyint)*omega(rayyint)))
          vfv(rayyint)=vfv(rayyint)+(raystrtmp-raystr)
       endif
 
       if (rayy.lt.0._r8) goto 226
       bld_old=bld
      
       if (pbinc.gt.0._r8) then
          pb_old=pb(rayyint+1)
       endif

    enddo

    goto 226

225  continue

    !         if(xdom*(1.-bldfrac)) then
    call corner_up(bldfrac,xdom,rayx,rayy,rayyint, atan(rx/max(1.e-6_r8,(-ry))),strike)
    if (strike) then
       ! next wall layer higher is being lit
       rayyint2=rayyint+1
       !! WALLS (at bottom, near road)
       wfact=raystr*max(bld-bld_old,0._r8)*pb(rayyint2+1)/max(1.e-6_r8,pb(rayyint2+1))*pb(rayyint2+1)
       vfw(rayyint2)=vfw(rayyint2)+wfact
       raystr=raystr-wfact
       !! GROUND-LEVEL 'ROOFS'
       rfact=raystr*sseff(rayyint+1)
       ! important, otherwise could end up with negative ray strength!                           
       vfr(rayyint+1)=vfr(rayyint+1)+rfact                 
       raystr=raystr-rfact
    endif
    !endif 

226  continue

    !! ROADS
    vft=vft+raystr
    !if (raystr.gt.minray.and.amod(rayx-rx*dray,xdom)/xdom.gt.1._r8-bldfrac.and.pb(2).gt.0.999999_r8) then   
    if (raystr.gt.minray.and.((rayx-rx*dray) - xdom * &
       floor((rayx-rx*dray)/xdom))/xdom.gt.1._r8-bldfrac.and.pb(2).gt.0.999999_r8) then                                        
        write(6,*)'PROBLEM (ray_dn), radiation reaching building interior ground,raystr,rayx,rayy=',raystr,rayx,rayy
        write(6,*)'xdom,bldfrac',xdom,bldfrac
        write(6,*)amod(rayx-rx*dray,xdom)/xdom,1._r8-bldfrac
        write(6,*)'izcan',izcan
        write(6,*)'xx,zz,dist',xx,zz,dist32
        do izcan=1,nzcanm
           write(6,*)'i,ss,sseff',izcan,ss(izcan),sseff(izcan)
        enddo
        write(6,*)'minray, pb',minray,pb
        stop
    endif
  end subroutine ray_dn

  !-----------------------------------------------------------------------
  subroutine ray_up_foliage(minray, maxind, rayx, rayy, bld, bld_old, raystr, &
                            xx, zz, hemi32, xdom, dray, dzcan, bldfrac, pb, &
                            kbs, lad_tree, lad, omega, horiz, vfb, vfv)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Computes the contribution to view factors of rays travelling "upwards"
    ! with a positive vertical component. Assumes upward travelling rays leave
    ! leaves clumped within individual tree canopies, seeing a larger lad 
    ! above them than the average canyon lad.
    !-----------------------------------------------------------------------
    ! !ARGUMENTS:
    implicit none
    integer,  parameter  :: nzcanm = 4       ! Maximum number of vertical levels at urban resolution     
    real(r8), intent(in) :: minray
    real(r8), intent(inout) :: rayx, rayy, bld, bld_old, raystr
    integer,  intent(in) :: maxind
    real(r8), intent(in) :: xx, zz, hemi32, xdom, dray, dzcan, bldfrac
    real(r8), intent(in) :: kbs
    real(r8), intent(in) :: omega(nzcanm), lad(nzcanm), lad_tree(nzcanm)
    real(r8), intent(in) :: pb(nzcanm+1)
    logical,  intent(in) :: horiz

    real(r8), intent(out) :: vfb(nzcanm), vfv(nzcanm)
    
    ! confirm: these output are not needed?
    !real(r8), intent(out) :: vfv2(nzcanm),vfa2(nzcanm+1),vfa(nzcanm+1)

    ! !LOCAL VARIABLES:
    integer              :: rayyint, izcan, rayyint_old, rayyint2
    real(r8)             :: rx, ry, xfr, wfact, raystrtmp
    real(r8)             :: raystra
    logical              :: strike
    ! confirm: these output are not needed?
    ! real(r8)             :: atotot(nzcanm+1),atotot2(nzcanm+1)
    !-----------------------------------------------------------------------
    !        raystra=raystr
    ! rays originating from walls
    rx=zz
    ry=xx
    ! rays originating from roads/vegetation
    if (horiz) then
       rx = xx
       ry = zz
    endif

    do izcan=1,nzcanm
       vfb(izcan)=0._r8
       vfv(izcan)=0._r8
    enddo

    do while (raystr.gt.minray)
       rayyint_old=ceiling(rayy)
       rayx=rayx+rx*dray
       rayy=rayy+ry*dray
       rayyint=ceiling(rayy)
       !xfr=amod(rayx,xdom)
       xfr = rayx - xdom* floor(rayx /xdom)
       if (xfr/xdom.ge.1._r8-bldfrac) then
          ! interception by walls
          bld=1._r8
          !if(bld.gt.bld_old) raystra=raystr
          rayyint2=rayyint
          if (rayyint.gt.rayyint_old) then
             ! the inputs to corner_up are slightly different than in ray_dn
             call corner_up(bldfrac,xdom,rayx,-rayy, -rayyint_old,&
                            atan(rx/max(1.e-6_r8,ry)),strike)
             if (strike) then
                ! lower wall layer is being lit
  		          rayyint2=rayyint_old
             endif
          endif
          wfact=raystr*max(bld-bld_old,0._r8)*pb(rayyint2+1)/max(1.e-6_r8,pb(rayyint2+1))*pb(rayyint2+1)
          vfb(rayyint2)=vfb(rayyint2)+wfact
          raystr=raystr-wfact
       else
          bld=0._r8
          ! interception by vegetation and air in the canopy column
          raystrtmp=raystr
          raystr=raystr*exp(-hemi32*dray*dzcan*&
          !(kbs*0.5*(lad(rayyint)+lad_tree*omega(rayyint)&
                 (kbs*0.5_r8*(lad(rayyint)+lad_tree(rayyint))*omega(rayyint)))

          vfv(rayyint)=vfv(rayyint)+(raystrtmp-raystr)
       endif
       
       if(rayyint.gt.maxind+1) goto 541
       bld_old=bld
    enddo

541   continue
  end subroutine ray_up_foliage

!-------------------[kz.11]Ray tracing test-------------------------   
   
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




