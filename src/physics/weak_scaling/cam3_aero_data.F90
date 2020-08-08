module cam3_aero_data
!-----------------------------------------------------------------------
!
! Purposes:
!       read, store, interpolate, and return fields
!         of aerosols to CAM.  The initialization
!         file (mass.nc) is assumed to be a monthly climatology
!         of aerosols from MATCH (on a sigma pressure
!         coordinate system).
!       also provide a "background" aerosol field to correct
!         for any deficiencies in the physical parameterizations
!         This fields is a "tuning" parameter.
!       Public methods:
!       (1) - initialization
!          read aerosol masses from external file
!             also pressure coordinates
!          convert from monthly average values to mid-month values
!       (2) - interpolation (time and vertical)
!          interpolate onto pressure levels of CAM
!          interpolate to time step of CAM
!          return mass of aerosols
!
!-----------------------------------------------------------------------

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_scam_mod,   only: shr_scam_GetCloseLatLon
  use spmd_utils,     only: masterproc
  use ppgrid,         only: pcols, pver, pverp, begchunk, endchunk
  use phys_grid,      only: get_ncols_p
  use time_manager,   only: get_curr_calday
  use infnan,         only: nan, assignment(=)
  use cam_abortutils, only: endrun
  use scamMod,        only: scmlon,scmlat,single_column
  use error_messages, only: handle_ncerr
  use physics_types,  only: physics_state
  use boundarydata,   only: boundarydata_init, boundarydata_type
  use perf_mod,       only: t_startf, t_stopf
  use cam_logfile,    only: iulog
  use netcdf

  implicit none
  private
  save

  public :: &
     cam3_aero_data_readnl,       & ! read namelist
     cam3_aero_data_register,     & ! register these aerosols with pbuf2d
     cam3_aero_data_init,         & ! read from file, interpolate onto horiz grid
     cam3_aero_data_timestep_init   ! update data-aerosols to this timestep

  ! namelist variables
  logical, public :: cam3_aero_data_on = .false.
  character(len=256) :: bndtvaer = 'bndtvaer'   ! full pathname for time-variant aerosol mass climatology dataset

  ! naer is number of species in climatology
  integer, parameter :: naer = 11

  real(r8), parameter :: wgt_sscm = 6.0_r8 / 7.0_r8 ! Fraction of total seasalt mass in coarse mode

  ! indices to aerosol array (species portion)
  integer, parameter :: &
      idxSUL   =  1, &
      idxSSLTA =  2, & ! accumulation mode
      idxSSLTC =  3, & ! coarse mode
      idxOCPHO =  8, &
      idxBCPHO =  9, &
      idxOCPHI =  10, &
      idxBCPHI = 11

  ! indices to sections of array that represent
  ! groups of aerosols
  integer, parameter :: &
      idxSSLTfirst    = 2, numSSLT  = 2, &
      idxDUSTfirst    = 4, &
      numDUST         = 4, &
      idxCARBONfirst = 8, &
      numCARBON      = 4

  ! names of aerosols are they are represented in
  ! the climatology file.
  ! Appended '_V' indicates field has been vertically summed.
  character(len=8), parameter :: aerosol_name(naer) =  &
     (/"MSUL_V  "&
      ,"MSSLTA_V"&
      ,"MSSLTC_V"&
      ,"MDUST1_V"&
      ,"MDUST2_V"&
      ,"MDUST3_V"&
      ,"MDUST4_V"&
      ,"MOCPHO_V"&
      ,"MBCPHO_V"&
      ,"MOCPHI_V"&
      ,"MBCPHI_V"/)

  ! number of different "groups" of aerosols
  integer, parameter :: num_aer_groups=4

  ! which group does each bin belong to?
  integer, dimension(naer), parameter ::  &
      group =(/1,2,2,3,3,3,3,4,4,4,4/)

  ! name of each group
  character(len=10), dimension(num_aer_groups), parameter :: &
      aerosol_names = (/'sul  ','sslt ','dust ','car  '/)

  ! this boundarydata_type is used for datasets in the ncols format only.
  type(boundarydata_type) :: aerosol_datan

  integer :: aernid = -1           ! netcdf id for aerosol file (init to invalid)
  integer :: species_id(naer) = -1 ! netcdf_id of each aerosol species (init to invalid)
  integer :: Mpsid                 ! netcdf id for MATCH PS
  integer :: nm = 1                ! index to prv month in array. init to 1 and toggle between 1 and 2
  integer :: np = 2                ! index to nxt month in array. init to 2 and toggle between 1 and 2
  integer :: mo_nxt = huge(1)      ! index to nxt month in file

  real(r8) :: cdaym                ! calendar day of prv month
  real(r8) :: cdayp                ! calendar day of next month

  ! aerosol mass
  real(r8), allocatable :: aer_mass(:, :, :, :)

  ! Days into year for mid month date
  ! This variable is dumb, the dates are in the dataset to be read in but they are
  ! slightly different than this so getting rid of it causes a change which
  ! exceeds roundoff.
  real(r8) :: Mid(12) = (/16.5_r8,  46.0_r8,  75.5_r8, 106.0_r8, 136.5_r8, 167.0_r8, &
                         197.5_r8, 228.5_r8, 259.0_r8, 289.5_r8, 320.0_r8, 350.5_r8 /)

  !  values read from file and temporary values used for interpolation
  !
  !  aerosolc is:
  !  Cumulative Mass at midpoint of each month
  !    on CAM's horizontal grid (col)
  !    on MATCH's levels (lev)
  !  aerosolc
  integer, parameter :: paerlev = 28       ! number of levels for aerosol fields (MUST = naerlev)
  integer :: naerlev                       ! size of level dimension in MATCH data
  integer :: naerlon
  integer :: naerlat

  ! indices for fields in the physics buffer
  integer :: cam3_sul_idx, cam3_ssam_idx, cam3_sscm_idx, &
      cam3_dust1_idx, cam3_dust2_idx, cam3_dust3_idx, cam3_dust4_idx,&
      cam3_ocpho_idx, cam3_bcpho_idx, cam3_ocphi_idx, cam3_bcphi_idx

!================================================================================================
contains
!================================================================================================

subroutine cam3_aero_data_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'cam3_aero_data_readnl'

   namelist /cam3_aero_data_nl/ cam3_aero_data_on, bndtvaer
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cam3_aero_data_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, cam3_aero_data_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(cam3_aero_data_on, 1, mpilog, 0, mpicom)
   call mpibcast(bndtvaer, len(bndtvaer), mpichar, 0, mpicom)
#endif

   ! Prevent using these before they are set.
   cdaym = nan
   cdayp = nan

end subroutine cam3_aero_data_readnl

!================================================================================================

subroutine cam3_aero_data_register

   ! register old prescribed aerosols with physics buffer

   use physics_buffer, only: pbuf_add_field, dtype_r8

   call pbuf_add_field('cam3_sul',  'physpkg',dtype_r8,(/pcols,pver/),cam3_sul_idx)
   call pbuf_add_field('cam3_ssam', 'physpkg',dtype_r8,(/pcols,pver/),cam3_ssam_idx)
   call pbuf_add_field('cam3_sscm', 'physpkg',dtype_r8,(/pcols,pver/),cam3_sscm_idx)
   call pbuf_add_field('cam3_dust1','physpkg',dtype_r8,(/pcols,pver/),cam3_dust1_idx)
   call pbuf_add_field('cam3_dust2','physpkg',dtype_r8,(/pcols,pver/),cam3_dust2_idx)
   call pbuf_add_field('cam3_dust3','physpkg',dtype_r8,(/pcols,pver/),cam3_dust3_idx)
   call pbuf_add_field('cam3_dust4','physpkg',dtype_r8,(/pcols,pver/),cam3_dust4_idx)
   call pbuf_add_field('cam3_ocpho','physpkg',dtype_r8,(/pcols,pver/),cam3_ocpho_idx)
   call pbuf_add_field('cam3_bcpho','physpkg',dtype_r8,(/pcols,pver/),cam3_bcpho_idx)
   call pbuf_add_field('cam3_ocphi','physpkg',dtype_r8,(/pcols,pver/),cam3_ocphi_idx)
   call pbuf_add_field('cam3_bcphi','physpkg',dtype_r8,(/pcols,pver/),cam3_bcphi_idx)

end subroutine cam3_aero_data_register

!================================================================================================

subroutine cam3_aero_data_init(phys_state)
!------------------------------------------------------------------
!  Reads in:
!     file from which to read aerosol Masses on CAM grid. Currently
!        assumed to be MATCH ncep runs, averaged by month.
!     NOTE (Data have been externally interpolated onto CAM grid
!        and backsolved to provide Mid-month values)
!
!  Populates:
!     module variables:
!       aerosolc(pcols,paerlev+1,begchunk:endchunk,naer,2))
!       aerosolc(  column_index
!                , level_index (match levels)
!                , chunk_index
!                , species_index
!                , month = 1:2 )
!       M_hybi(level_index = Lev_MATCH) = pressure at mid-level.
!       M_ps_cam_col(column,chunk,month) ! PS from MATCH on Cam Columns
!
!  Method:
!    read data from file
!    allocate memory for storage of aerosol data on CAM horizontal grid
!    distribute data to remote nodes
!    populates the module variables
!
!------------------------------------------------------------------
   use ioFileMod,    only: getfil

#if ( defined SPMD )
   use mpishorthand
#endif
   type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

! local variables

   integer :: naerlev

   integer dateid                       ! netcdf id for date variable
   integer secid                        ! netcdf id for seconds variable
   integer londimid                     ! netcdf id for longitude dimension
   integer latdimid                     ! netcdf id for latitude dimension
   integer levdimid                     ! netcdf id for level dimension

   integer timesiz                      ! number of time samples (=12) in netcdf file
   integer latid                        ! netcdf id for latitude variable
   integer Mhybiid                      ! netcdf id for MATCH hybi
   integer timeid                       ! netcdf id for time variable
   integer dimids(nf90_max_var_dims)      ! variable shape
   integer :: start(4)                  ! start vector for netcdf calls
   integer :: kount(4)                  ! count vector for netcdf calls
   integer mo                           ! month index
   integer m                            ! constituent index
   integer :: n                         ! loop index
   integer :: i,j,k                     ! spatial indices
   integer :: date_aer(12)              ! Date on aerosol dataset (YYYYMMDD)
   integer :: attnum                    ! attribute number
   integer :: ierr                      ! netcdf return code
   real(r8) ::  coldata(paerlev)    ! aerosol field read in from dataset
   integer :: ret
   integer mo_prv                       ! index to previous month
   integer latidx,lonidx

   character(len=8) :: aname                   ! temporary aerosol name
   character(len=8) :: tmp_aero_name(naer) ! name for input to boundary data

   character(len=256) :: locfn          ! netcdf local filename to open
!
! aerosol_data will be read in from the aerosol boundary dataset, then scattered to chunks
! after filling in the bottom level with zeros
!
   real(r8), allocatable :: aerosol_data(:,:,:)    ! aerosol field read in from dataset
   real(r8), allocatable :: aerosol_field(:,:,:)   ! (plon,paerlev+1,plat)  aerosol field to be scattered
   real(r8) :: caldayloc                           ! calendar day of current timestep
   real(r8) :: closelat,closelon

   character(len=*), parameter :: subname = 'cam3_aero_data_init'
   !------------------------------------------------------------------

   call endrun(subname//': not supported with weak scaling fix')
end subroutine cam3_aero_data_init

!================================================================================================

subroutine cam3_aero_data_timestep_init(pbuf2d,  phys_state)
!------------------------------------------------------------------
!
!  Input:
!     time at which aerosol masses are needed (get_curr_calday())
!     chunk index
!     CAM's vertical grid (pint)
!
!  Output:
!     values for Aerosol Mass at time specified by get_curr_calday
!     on vertical grid specified by pint (aer_mass) :: aerosol at time t
!
!  Method:
!     first determine which indexs of aerosols are the bounding data sets
!     interpolate both onto vertical grid aerm(),aerp().
!     from those two, interpolate in time.
!
!------------------------------------------------------------------

   use interpolate_data, only: get_timeinterp_factors

   use physics_buffer, only: physics_buffer_desc, dtype_r8, pbuf_set_field, pbuf_get_chunk
   use cam_logfile,     only: iulog
   use ppgrid,          only: begchunk,endchunk
   use physconst,       only: gravit

!
! aerosol fields interpolated to current time step
!   on pressure levels of this time step.
! these should be made read-only for other modules
! Is allocation done correctly here?
!

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   type(physics_state), intent(in), dimension(begchunk:endchunk) :: phys_state

!
! Local workspace
!
   type(physics_buffer_desc), pointer :: phys_buffer_chunk(:)
   real(r8) :: pint(pcols,pverp)  ! interface pres.
   integer :: c                           ! chunk index
   real(r8) caldayloc                     ! calendar day of current timestep
   real(r8) fact1, fact2                  ! time interpolation factors

   integer i, k, j                        ! spatial indices
   integer m                              ! constituent index
   integer lats(pcols),lons(pcols)        ! latitude and longitudes of column
   integer ncol                           ! number of columns
   integer lchnk                          ! chunk index

   real(r8) speciesmin(naer)              ! minimal value for each species
!
! values before current time step "the minus month"
! aerosolm(pcols,pver) is value of preceeding month's aerosol masses
! aerosolp(pcols,pver) is value of next month's aerosol masses
!  (think minus and plus or values to left and right of point to be interpolated)
!
   real(r8) aerosolm(pcols,pver,naer,begchunk:endchunk) ! aerosol mass from MATCH in column,level at previous (minus) month
!
! values beyond (or at) current time step "the plus month"
!
   real(r8) aerosolp(pcols,pver,naer,begchunk:endchunk) ! aerosol mass from MATCH in column,level at next (plus) month
   real(r8) :: mass_to_mmr(pcols,pver)

   character(len=*), parameter :: subname = 'cam3_aero_data_timestep_init'

   logical error_found
   !------------------------------------------------------------------
   call endrun(subname//': not supported with weak scaling fix')

end subroutine cam3_aero_data_timestep_init

!================================================================================================

subroutine vert_interpolate (Match_ps, pint, n, aerosol_mass, ncol, c)
!--------------------------------------------------------------------
! Input: match surface pressure, cam interface pressure,
!        month index, number of columns, chunk index
!
! Output: Aerosol mass mixing ratio (aerosol_mass)
!
! Method:
!         interpolate column mass (cumulative) from match onto
!           cam's vertical grid (pressure coordinate)
!         convert back to mass mixing ratio
!
!--------------------------------------------------------------------

   real(r8), intent(out) :: aerosol_mass(pcols,pver,naer)  ! aerosol mass from MATCH
   real(r8), intent(in) :: Match_ps(pcols)                ! surface pressure at a particular month
   real(r8), intent(in) :: pint(pcols,pverp)              ! interface pressure from CAM

   integer, intent(in) :: ncol,c                          ! chunk index and number of columns
   integer, intent(in) :: n                               ! prv or nxt month index
!
! Local workspace
!
   integer m                           ! index to aerosol species
   integer kupper(pcols)               ! last upper bound for interpolation
   integer i, k, kk, kkstart, kount    ! loop vars for interpolation
   integer isv, ksv, msv               ! loop indices to save

   logical bad                         ! indicates a bad point found
   logical lev_interp_comp             ! interpolation completed for a level
   logical error_found

   real(r8) aerosol(pcols,pverp,naer)  ! cumulative mass of aerosol in column beneath upper
                                       ! interface of level in column at particular month
   real(r8) dpl, dpu                   ! lower and upper intepolation factors
   real(r8) v_coord                    ! vertical coordinate
   real(r8) AER_diff                   ! temp var for difference between aerosol masses

   character(len=*), parameter :: subname = 'cam3_aero_data.vert_interpolate'
   !-----------------------------------------------------------------------
   call endrun(subname//': not supported with weak scaling fix')
end subroutine vert_interpolate

!================================================================================================

subroutine aerint (phys_state)

   type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

   integer :: ntmp                                ! used in index swapping
   integer :: start(4)                            ! start vector for netcdf calls
   integer :: kount(4)                            ! count vector for netcdf calls
   integer :: i,j,k                               ! spatial indices
   integer :: m                                   ! constituent index
   integer :: cols, cole
   integer :: lchnk, ncol
   real(r8) :: caldayloc                          ! calendar day of current timestep
   real(r8) :: aerosol_data(naerlon,naerlat,paerlev)    ! aerosol field read in from dataset
   real(r8) :: aerosol_field(naerlon,paerlev+1,naerlat) ! aerosol field to be scattered
   integer latidx,lonidx
   real(r8) closelat,closelon

   character(len=*), parameter :: subname = 'cam3_aero_data.aerint'
   !-----------------------------------------------------------------------
   call endrun(subname//': not supported with weak scaling fix')
end subroutine aerint

end module cam3_aero_data
