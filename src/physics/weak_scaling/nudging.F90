module nudging
!=====================================================================
!
! Purpose: Implement Nudging of the model state of U,V,T,Q, and/or PS
!          toward specified values from analyses.
!
! Author: Patrick Callaghan
!
! Description:
!
!    This module assumes that the user has {U,V,T,Q,PS} values from analyses
!    which have been preprocessed onto the current model grid and adjusted
!    for differences in topography. It is also assumed that these resulting
!    values and are stored in individual files which are indexed with respect
!    to year, month, day, and second of the day. When the model is inbetween
!    the given begining and ending times, a relaxation forcing is added to
!    nudge the model toward the analyses values determined from the forcing
!    option specified. After the model passes the ending analyses time, the
!    forcing discontinues.
!
!    Some analyses products can have gaps in the available data, where values
!    are missing for some interval of time. When files are missing, the nudging
!    force is switched off for that interval of time, so we effectively 'coast'
!    thru the gap.
!
!    Currently, the nudging module is set up to accomodate nudging of PS
!    values, however that functionality requires forcing that is applied in
!    the selected dycore and is not yet implemented.
!
!    The nudging of the model toward the analyses data is controlled by
!    the 'nudging_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which nudging is applied, the strength of the nudging
!    tendencies, and its spatial distribution.
!
!    FORCING:
!    --------
!    Nudging tendencies are applied as a relaxation force between the current
!    model state values and target state values derived from the avalilable
!    analyses. The form of the target values is selected by the 'Nudge_Force_Opt'
!    option, the timescale of the forcing is determined from the given
!    'Nudge_TimeScale_Opt', and the nudging strength Alpha=[0.,1.] for each
!    variable is specified by the 'Nudge_Xcoef' values. Where X={U,V,T,Q,PS}
!
!           F_nudge = Alpha*((Target-Model(t_curr))/TimeScale
!
!
!    WINDOWING:
!    ----------
!    The region of applied nudging can be limited using Horizontal/Vertical
!    window functions that are constructed using a parameterization of the
!    Heaviside step function.
!
!    The Heaviside window function is the product of separate horizonal and vertical
!    windows that are controled via 12 parameters:
!
!        Nudge_Hwin_lat0:     Specify the horizontal center of the window in degrees.
!        Nudge_Hwin_lon0:     The longitude must be in the range [0,360] and the
!                             latitude should be [-90,+90].
!        Nudge_Hwin_latWidth: Specify the lat and lon widths of the window as positive
!        Nudge_Hwin_lonWidth: values in degrees.Setting a width to a large value (e.g. 999)
!                             renders the window a constant in that direction.
!        Nudge_Hwin_latDelta: Controls the sharpness of the window transition with a
!        Nudge_Hwin_lonDelta: length in degrees. Small non-zero values yeild a step
!                             function while a large value yeilds a smoother transition.
!        Nudge_Hwin_Invert  : A logical flag used to invert the horizontal window function
!                             to get its compliment.(e.g. to nudge outside a given window).
!
!        Nudge_Vwin_Lindex:   In the vertical, the window is specified in terms of model
!        Nudge_Vwin_Ldelta:   level indcies. The High and Low transition levels should
!        Nudge_Vwin_Hindex:   range from [0,(NLEV+1)]. The transition lengths are also
!        Nudge_Vwin_Hdelta:   specified in terms of model indices. For a window function
!                             constant in the vertical, the Low index should be set to 0,
!                             the High index should be set to (NLEV+1), and the transition
!                             lengths should be set to 0.001
!        Nudge_Vwin_Invert  : A logical flag used to invert the vertical window function
!                             to get its compliment.
!
!        EXAMPLE: For a channel window function centered at the equator and independent
!                 of the vertical (30 levels):
!                        Nudge_Hwin_lat0     = 0.         Nudge_Vwin_Lindex = 0.
!                        Nudge_Hwin_latWidth = 30.        Nudge_Vwin_Ldelta = 0.001
!                        Nudge_Hwin_latDelta = 5.0        Nudge_Vwin_Hindex = 31.
!                        Nudge_Hwin_lon0     = 180.       Nudge_Vwin_Hdelta = 0.001
!                        Nudge_Hwin_lonWidth = 999.       Nudge_Vwin_Invert = .false.
!                        Nudge_Hwin_lonDelta = 1.0
!                        Nudge_Hwin_Invert   = .false.
!
!                 If on the other hand one wanted to apply nudging at the poles and
!                 not at the equator, the settings would be similar but with:
!                        Nudge_Hwin_Invert = .true.
!
!    A user can preview the window resulting from a given set of namelist values before
!    running the model. Lookat_NudgeWindow.ncl is a script avalable in the tools directory
!    which will read in the values for a given namelist and display the resulting window.
!
!    The module is currently configured for only 1 window function. It can readily be
!    extended for multiple windows if the need arises.
!
!
! Input/Output Values:
!    Forcing contributions are available for history file output by
!    the names:    {'Nudge_U','Nudge_V','Nudge_T',and 'Nudge_Q'}
!    The target values that the model state is nudged toward are available for history
!    file output via the variables:  {'Target_U','Target_V','Target_T',and 'Target_Q'}
!
!    &nudging_nl
!      Nudge_Model         - LOGICAL toggle to activate nudging.
!                              TRUE  -> Nudging is on.
!                              FALSE -> Nudging is off.                            [DEFAULT]
!
!      Nudge_Path          - CHAR path to the analyses files.
!                              (e.g. '/glade/scratch/USER/inputdata/nudging/ERAI-Data/')
!
!      Nudge_File_Template - CHAR Analyses filename with year, month, day, and second
!                                 values replaced by %y, %m, %d, and %s respectively.
!                              (e.g. '%y/ERAI_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc')
!
!      Nudge_Times_Per_Day - INT Number of analyses files available per day.
!                              1 --> daily analyses.
!                              4 --> 6 hourly analyses.
!                              8 --> 3 hourly.
!
!      Model_Times_Per_Day - INT Number of times to update the model state (used for nudging)
!                                each day. The value is restricted to be longer than the
!                                current model timestep and shorter than the analyses
!                                timestep. As this number is increased, the nudging
!                                force has the form of newtonian cooling.
!                              48 --> 1800 Second timestep.
!                              96 -->  900 Second timestep.
!
!      Nudge_Beg_Year      - INT nudging begining year.  [1979- ]
!      Nudge_Beg_Month     - INT nudging begining month. [1-12]
!      Nudge_Beg_Day       - INT nudging begining day.   [1-31]
!      Nudge_End_Year      - INT nudging ending year.    [1979-]
!      Nudge_End_Month     - INT nudging ending month.   [1-12]
!      Nudge_End_Day       - INT nudging ending day.     [1-31]
!
!      Nudge_Force_Opt     - INT Index to select the nudging Target for a relaxation
!                                forcing of the form:
!                                where (t'==Analysis times ; t==Model Times)
!
!                              0 -> NEXT-OBS: Target=Anal(t'_next)                 [DEFAULT]
!                              1 -> LINEAR:   Target=(F*Anal(t'_curr) +(1-F)*Anal(t'_next))
!                                                 F =(t'_next - t_curr )/Tdlt_Anal
!
!      Nudge_TimeScale_Opt - INT Index to select the timescale for nudging.
!                                where (t'==Analysis times ; t==Model Times)
!
!                              0 -->  TimeScale = 1/Tdlt_Anal                      [DEFAULT]
!                              1 -->  TimeScale = 1/(t'_next - t_curr )
!
!      Nudge_Uprof         - INT index of profile structure to use for U.  [0,1,2]
!      Nudge_Vprof         - INT index of profile structure to use for V.  [0,1,2]
!      Nudge_Tprof         - INT index of profile structure to use for T.  [0,1,2]
!      Nudge_Qprof         - INT index of profile structure to use for Q.  [0,1,2]
!      Nudge_PSprof        - INT index of profile structure to use for PS. [0,N/A]
!
!                                The spatial distribution is specified with a profile index.
!                                 Where:  0 == OFF      (No Nudging of this variable)
!                                         1 == CONSTANT (Spatially Uniform Nudging)
!                                         2 == HEAVISIDE WINDOW FUNCTION
!
!      Nudge_Ucoef         - REAL fractional nudging coeffcient for U.
!      Nudge_Vcoef         - REAL fractional nudging coeffcient for V.
!      Nudge_Tcoef         - REAL fractional nudging coeffcient for T.
!      Nudge_Qcoef         - REAL fractional nudging coeffcient for Q.
!      Nudge_PScoef        - REAL fractional nudging coeffcient for PS.
!
!                                 The strength of the nudging is specified as a fractional
!                                 coeffcient between [0,1].
!
!      Nudge_Hwin_lat0     - REAL latitudinal center of window in degrees.
!      Nudge_Hwin_lon0     - REAL longitudinal center of window in degrees.
!      Nudge_Hwin_latWidth - REAL latitudinal width of window in degrees.
!      Nudge_Hwin_lonWidth - REAL longitudinal width of window in degrees.
!      Nudge_Hwin_latDelta - REAL latitudinal transition length of window in degrees.
!      Nudge_Hwin_lonDelta - REAL longitudinal transition length of window in degrees.
!      Nudge_Hwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!      Nudge_Vwin_Lindex   - REAL LO model index of transition
!      Nudge_Vwin_Hindex   - REAL HI model index of transition
!      Nudge_Vwin_Ldelta   - REAL LO transition length
!      Nudge_Vwin_Hdelta   - REAL HI transition length
!      Nudge_Vwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!    /
!
!================
!
! TO DO:
! -----------
!    ** Implement Ps Nudging????
!
!=====================================================================
  ! Useful modules
  !------------------
  use shr_kind_mod,   only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,   only:timemgr_time_ge,timemgr_time_inc,get_curr_date,get_step_size
  use cam_abortutils, only:endrun
  use spmd_utils  ,   only:masterproc
  use cam_logfile ,   only:iulog
#ifdef SPMD
  use mpishorthand
#endif

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: Nudge_Model,Nudge_ON
  public:: nudging_readnl
  public:: nudging_init
  public:: nudging_timestep_init
  public:: nudging_timestep_tend
  private::nudging_update_analyses_se
  private::nudging_update_analyses_eul
  private::nudging_update_analyses_fv
  private::nudging_set_PSprofile
  private::nudging_set_profile
  private::calc_DryStaticEnergy

  ! Nudging Parameters
  !--------------------
  logical          :: Nudge_Model       =.false.
  logical          :: Nudge_ON          =.false.
  logical          :: Nudge_Initialized =.false.
  character(len=cl):: Nudge_Path
  character(len=cs):: Nudge_File,Nudge_File_Template
  integer          :: Nudge_Force_Opt
  integer          :: Nudge_TimeScale_Opt
  integer          :: Nudge_TSmode
  integer          :: Nudge_Times_Per_Day
  integer          :: Model_Times_Per_Day
  real(r8)         :: Nudge_Ucoef,Nudge_Vcoef
  integer          :: Nudge_Uprof,Nudge_Vprof
  real(r8)         :: Nudge_Qcoef,Nudge_Tcoef
  integer          :: Nudge_Qprof,Nudge_Tprof
  real(r8)         :: Nudge_PScoef
  integer          :: Nudge_PSprof
  integer          :: Nudge_Beg_Year ,Nudge_Beg_Month
  integer          :: Nudge_Beg_Day  ,Nudge_Beg_Sec
  integer          :: Nudge_End_Year ,Nudge_End_Month
  integer          :: Nudge_End_Day  ,Nudge_End_Sec
  integer          :: Nudge_Curr_Year,Nudge_Curr_Month
  integer          :: Nudge_Curr_Day ,Nudge_Curr_Sec
  integer          :: Nudge_Next_Year,Nudge_Next_Month
  integer          :: Nudge_Next_Day ,Nudge_Next_Sec
  integer          :: Nudge_Step
  integer          :: Model_Curr_Year,Model_Curr_Month
  integer          :: Model_Curr_Day ,Model_Curr_Sec
  integer          :: Model_Next_Year,Model_Next_Month
  integer          :: Model_Next_Day ,Model_Next_Sec
  integer          :: Model_Step
  real(r8)         :: Nudge_Hwin_lat0
  real(r8)         :: Nudge_Hwin_latWidth
  real(r8)         :: Nudge_Hwin_latDelta
  real(r8)         :: Nudge_Hwin_lon0
  real(r8)         :: Nudge_Hwin_lonWidth
  real(r8)         :: Nudge_Hwin_lonDelta
  logical          :: Nudge_Hwin_Invert = .false.
  real(r8)         :: Nudge_Hwin_lo
  real(r8)         :: Nudge_Hwin_hi
  real(r8)         :: Nudge_Vwin_Hindex
  real(r8)         :: Nudge_Vwin_Hdelta
  real(r8)         :: Nudge_Vwin_Lindex
  real(r8)         :: Nudge_Vwin_Ldelta
  logical          :: Nudge_Vwin_Invert =.false.
  real(r8)         :: Nudge_Vwin_lo
  real(r8)         :: Nudge_Vwin_hi
  real(r8)         :: Nudge_Hwin_latWidthH
  real(r8)         :: Nudge_Hwin_lonWidthH
  real(r8)         :: Nudge_Hwin_max
  real(r8)         :: Nudge_Hwin_min

  ! Nudging State Arrays
  !-----------------------
  integer Nudge_nlon,Nudge_nlat,Nudge_ncol,Nudge_nlev
  real(r8),allocatable::Target_U     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_V     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_T     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_S     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_Q     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_PS    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: Model_U     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_V     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_T     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_S     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_Q     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_PS    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Utau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Vtau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Stau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Qtau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_PStau (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Ustep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Vstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Sstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Qstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_PSstep(:,:)    !(pcols,begchunk:endchunk)

  ! Nudging Observation Arrays
  !-----------------------------
  integer               Nudge_NumObs
  integer,allocatable:: Nudge_ObsInd(:)
  logical ,allocatable::Nudge_File_Present(:)
  real(r8),allocatable::Nobs_U (:,:,:,:) !(pcols,pver,begchunk:endchunk,Nudge_NumObs)
  real(r8),allocatable::Nobs_V (:,:,:,:) !(pcols,pver,begchunk:endchunk,Nudge_NumObs)
  real(r8),allocatable::Nobs_T (:,:,:,:) !(pcols,pver,begchunk:endchunk,Nudge_NumObs)
  real(r8),allocatable::Nobs_Q (:,:,:,:) !(pcols,pver,begchunk:endchunk,Nudge_NumObs)
  real(r8),allocatable::Nobs_PS(:,:,:)   !(pcols,begchunk:endchunk,Nudge_NumObs)

contains
  !================================================================
  subroutine nudging_readnl(nlfile)
   !
   ! NUDGING_READNL: Initialize default values controlling the Nudging
   !                 process. Then read namelist values to override
   !                 them.
   !===============================================================
   use ppgrid        ,only: pver
   use namelist_utils,only:find_group_name
   use units         ,only:getunit,freeunit
   !
   ! Arguments
   !-------------
   character(len=*),intent(in)::nlfile
   !
   ! Local Values
   !---------------
   integer ierr,unitn

   namelist /nudging_nl/ Nudge_Model,Nudge_Path,                       &
                         Nudge_File_Template,Nudge_Force_Opt,          &
                         Nudge_TimeScale_Opt,                          &
                         Nudge_Times_Per_Day,Model_Times_Per_Day,      &
                         Nudge_Ucoef ,Nudge_Uprof,                     &
                         Nudge_Vcoef ,Nudge_Vprof,                     &
                         Nudge_Qcoef ,Nudge_Qprof,                     &
                         Nudge_Tcoef ,Nudge_Tprof,                     &
                         Nudge_PScoef,Nudge_PSprof,                    &
                         Nudge_Beg_Year,Nudge_Beg_Month,Nudge_Beg_Day, &
                         Nudge_End_Year,Nudge_End_Month,Nudge_End_Day, &
                         Nudge_Hwin_lat0,Nudge_Hwin_lon0,              &
                         Nudge_Hwin_latWidth,Nudge_Hwin_lonWidth,      &
                         Nudge_Hwin_latDelta,Nudge_Hwin_lonDelta,      &
                         Nudge_Hwin_Invert,                            &
                         Nudge_Vwin_Lindex,Nudge_Vwin_Hindex,          &
                         Nudge_Vwin_Ldelta,Nudge_Vwin_Hdelta,          &
                         Nudge_Vwin_Invert

   ! Nudging is NOT initialized yet, For now
   ! Nudging will always begin/end at midnight.
   !--------------------------------------------
   Nudge_Initialized =.false.
   Nudge_ON          =.false.
   Nudge_Beg_Sec=0
   Nudge_End_Sec=0

   ! Set Default Namelist values
   !-----------------------------
   Nudge_Model         = .false.
   Nudge_Path          = './Data/YOTC_ne30np4_001/'
   Nudge_File_Template = 'YOTC_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc'
   Nudge_Force_Opt     = 0
   Nudge_TimeScale_Opt = 0
   Nudge_TSmode        = 0
   Nudge_Times_Per_Day = 4
   Model_Times_Per_Day = 4
   Nudge_Ucoef         = 0._r8
   Nudge_Vcoef         = 0._r8
   Nudge_Qcoef         = 0._r8
   Nudge_Tcoef         = 0._r8
   Nudge_PScoef        = 0._r8
   Nudge_Uprof         = 0
   Nudge_Vprof         = 0
   Nudge_Qprof         = 0
   Nudge_Tprof         = 0
   Nudge_PSprof        = 0
   Nudge_Beg_Year      = 2008
   Nudge_Beg_Month     = 5
   Nudge_Beg_Day       = 1
   Nudge_End_Year      = 2008
   Nudge_End_Month     = 9
   Nudge_End_Day       = 1
   Nudge_Hwin_lat0     = 0._r8
   Nudge_Hwin_latWidth = 9999._r8
   Nudge_Hwin_latDelta = 1.0_r8
   Nudge_Hwin_lon0     = 180._r8
   Nudge_Hwin_lonWidth = 9999._r8
   Nudge_Hwin_lonDelta = 1.0_r8
   Nudge_Hwin_Invert   = .false.
   Nudge_Hwin_lo       = 0.0_r8
   Nudge_Hwin_hi       = 1.0_r8
   Nudge_Vwin_Hindex   = float(pver+1)
   Nudge_Vwin_Hdelta   = 0.001_r8
   Nudge_Vwin_Lindex   = 0.0_r8
   Nudge_Vwin_Ldelta   = 0.001_r8
   Nudge_Vwin_Invert   = .false.
   Nudge_Vwin_lo       = 0.0_r8
   Nudge_Vwin_hi       = 1.0_r8

   ! Read in namelist values
   !------------------------
   if(masterproc) then
     unitn = getunit()
     open(unitn,file=trim(nlfile),status='old')
     call find_group_name(unitn,'nudging_nl',status=ierr)
     if(ierr.eq.0) then
       read(unitn,nudging_nl,iostat=ierr)
       if(ierr.ne.0) then
         call endrun('nudging_readnl:: ERROR reading namelist')
       endif
     endif
     close(unitn)
     call freeunit(unitn)
   endif

   ! Set hi/lo values according to the given '_Invert' parameters
   !--------------------------------------------------------------
   if(Nudge_Hwin_Invert) then
     Nudge_Hwin_lo = 1.0_r8
     Nudge_Hwin_hi = 0.0_r8
   else
     Nudge_Hwin_lo = 0.0_r8
     Nudge_Hwin_hi = 1.0_r8
   endif

   if(Nudge_Vwin_Invert) then
     Nudge_Vwin_lo = 1.0_r8
     Nudge_Vwin_hi = 0.0_r8
   else
     Nudge_Vwin_lo = 0.0_r8
     Nudge_Vwin_hi = 1.0_r8
   endif

   ! Check for valid namelist values
   !----------------------------------
   if((Nudge_Hwin_lat0.lt.-90._r8).or.(Nudge_Hwin_lat0.gt.+90._r8)) then
     write(iulog,*) 'NUDGING: Window lat0 must be in [-90,+90]'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lat0=',Nudge_Hwin_lat0
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_lon0.lt.0._r8).or.(Nudge_Hwin_lon0.ge.360._r8)) then
     write(iulog,*) 'NUDGING: Window lon0 must be in [0,+360)'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lon0=',Nudge_Hwin_lon0
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Vwin_Lindex.gt.Nudge_Vwin_Hindex)                         .or. &
      (Nudge_Vwin_Hindex.gt.float(pver+1)).or.(Nudge_Vwin_Hindex.lt.0._r8).or. &
      (Nudge_Vwin_Lindex.gt.float(pver+1)).or.(Nudge_Vwin_Lindex.lt.0._r8)   ) then
     write(iulog,*) 'NUDGING: Window Lindex must be in [0,pver+1]'
     write(iulog,*) 'NUDGING: Window Hindex must be in [0,pver+1]'
     write(iulog,*) 'NUDGING: Lindex must be LE than Hindex'
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Lindex=',Nudge_Vwin_Lindex
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Hindex=',Nudge_Vwin_Hindex
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_latDelta.le.0._r8).or.(Nudge_Hwin_lonDelta.le.0._r8).or. &
      (Nudge_Vwin_Hdelta  .le.0._r8).or.(Nudge_Vwin_Ldelta  .le.0._r8)    ) then
     write(iulog,*) 'NUDGING: Window Deltas must be positive'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_latDelta=',Nudge_Hwin_latDelta
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lonDelta=',Nudge_Hwin_lonDelta
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Hdelta=',Nudge_Vwin_Hdelta
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Ldelta=',Nudge_Vwin_Ldelta
     call endrun('nudging_readnl:: ERROR in namelist')

   endif

   if((Nudge_Hwin_latWidth.le.0._r8).or.(Nudge_Hwin_lonWidth.le.0._r8)) then
     write(iulog,*) 'NUDGING: Window widths must be positive'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_latWidth=',Nudge_Hwin_latWidth
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lonWidth=',Nudge_Hwin_lonWidth
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   ! Broadcast namelist variables
   !------------------------------
#ifdef SPMD
   call mpibcast(Nudge_Path         ,len(Nudge_Path)         ,mpichar,0,mpicom)
   call mpibcast(Nudge_File_Template,len(Nudge_File_Template),mpichar,0,mpicom)
   call mpibcast(Nudge_Model        , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Initialized  , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_ON           , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Force_Opt    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_TimeScale_Opt, 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_TSmode       , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Model_Times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Ucoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Tcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Qcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_PScoef       , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Uprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Vprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Tprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Qprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_PSprof       , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Hwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lat0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_latWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_latDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lon0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lonWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lonDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_Invert,   1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Vwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Hindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Hdelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Lindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Ldelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Invert,   1, mpilog, 0, mpicom)
#endif

   ! End Routine
   !------------
   return
  end subroutine ! nudging_readnl
  !================================================================


  !================================================================
  subroutine nudging_init
   !
   ! NUDGING_INIT: Allocate space and initialize Nudging values
   !===============================================================
   use ppgrid        ,only: pver,pcols,begchunk,endchunk
   use error_messages,only: alloc_err
   use dycore        ,only: dycore_is
   use dyn_grid      ,only: get_horiz_grid_dim_d
   use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p
   use cam_history   ,only: addfld
   use shr_const_mod ,only: SHR_CONST_PI
   use filenames     ,only: interpret_filename_spec

   ! Local values
   !----------------
   integer  Year,Month,Day,Sec
   integer  YMD1,YMD
   logical  After_Beg,Before_End
   integer  istat,lchnk,ncol,icol,ilev
   integer  hdim1_d,hdim2_d
   integer  dtime
   real(r8) rlat,rlon
   real(r8) Wprof(pver)
   real(r8) lonp,lon0,lonn,latp,lat0,latn
   real(r8) Val1_p,Val2_p,Val3_p,Val4_p
   real(r8) Val1_0,Val2_0,Val3_0,Val4_0
   real(r8) Val1_n,Val2_n,Val3_n,Val4_n
   integer               nn

   ! Get the time step size
   !------------------------
   dtime = get_step_size()

   call endrun('nudging_init: Nudging not supported with weak scaling fix')
  end subroutine ! nudging_init
  !================================================================


  !================================================================
  subroutine nudging_timestep_init(phys_state)
   !
   ! NUDGING_TIMESTEP_INIT:
   !                 Check the current time and update Model/Nudging
   !                 arrays when necessary. Toggle the Nudging flag
   !                 when the time is withing the nudging window.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state
   use constituents ,only: cnst_get_ind
   use dycore       ,only: dycore_is
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use filenames    ,only: interpret_filename_spec
   use ESMF

   ! Arguments
   !-----------
   type(physics_state),intent(in):: phys_state(begchunk:endchunk)

   call endrun('nudging_timestep_init: Nudging not supported with weak scaling fix')
  end subroutine ! nudging_timestep_init
  !================================================================


  !================================================================
  subroutine nudging_timestep_tend(phys_state,phys_tend)
   !
   ! NUDGING_TIMESTEP_TEND:
   !                If Nudging is ON, return the Nudging contributions
   !                to forcing using the current contents of the Nudge
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use cam_history  ,only: outfld

   ! Arguments
   !-------------
   type(physics_state), intent(in) :: phys_state
   type(physics_ptend), intent(out):: phys_tend
   call endrun('nudging_timestep_tend: Nudging not supported with weak scaling fix')
  end subroutine ! nudging_timestep_tend
  !================================================================


  !================================================================
  subroutine nudging_update_analyses_se(anal_file)
   !
   ! NUDGING_UPDATE_ANALYSES_SE:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer ncol,plev,istat
   integer ncid,varid
   real(r8) Xanal(Nudge_ncol,Nudge_nlev)
   real(r8) PSanal(Nudge_ncol)
   real(r8) Lat_anal(Nudge_ncol)
   real(r8) Lon_anal(Nudge_ncol)
   integer  nn,Nindex

   ! Rotate Nudge_ObsInd() indices, then check the existence of the analyses
   ! file; broadcast the updated indices and file status to all the other MPI nodes.
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   call endrun('nudging_update_analyses_se: Nudging not supported with weak scaling fix')
  end subroutine ! nudging_update_analyses_se
  !================================================================


  !================================================================
  subroutine nudging_update_analyses_eul(anal_file)
   !
   ! NUDGING_UPDATE_ANALYSES_EUL:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer nlon,nlat,plev,istat
   integer ncid,varid
   integer ilat,ilon,ilev
   real(r8) Xanal(Nudge_nlon,Nudge_nlat,Nudge_nlev)
   real(r8) PSanal(Nudge_nlon,Nudge_nlat)
   real(r8) Lat_anal(Nudge_nlat)
   real(r8) Lon_anal(Nudge_nlon)
   real(r8) Xtrans(Nudge_nlon,Nudge_nlev,Nudge_nlat)
   integer  nn,Nindex

   call endrun('nudging_update_analyses_eul: Nudging not supported with weak scaling fix')
  end subroutine ! nudging_update_analyses_eul
  !================================================================


  !================================================================
  subroutine nudging_update_analyses_fv(anal_file)
   !
   ! NUDGING_UPDATE_ANALYSES_FV:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer nlon,nlat,plev,istat
   integer ncid,varid
   integer ilat,ilon,ilev
   real(r8) Xanal(Nudge_nlon,Nudge_nlat,Nudge_nlev)
   real(r8) PSanal(Nudge_nlon,Nudge_nlat)
   real(r8) Lat_anal(Nudge_nlat)
   real(r8) Lon_anal(Nudge_nlon)
   real(r8) Xtrans(Nudge_nlon,Nudge_nlev,Nudge_nlat)
   integer  nn,Nindex

   call endrun('nudging_update_analyses_fv: Nudging not supported with weak scaling fix')
  end subroutine ! nudging_update_analyses_fv
  !================================================================


  !================================================================
  subroutine nudging_set_profile(rlat,rlon,Nudge_prof,Wprof,nlev)
   !
   ! NUDGING_SET_PROFILE: for the given lat,lon, and Nudging_prof, set
   !                      the verical profile of window coeffcients.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   integer  nlev,Nudge_prof
   real(r8) rlat,rlon
   real(r8) Wprof(nlev)

   call endrun('nudging_set_profile: Nudging not supported with weak scaling fix')
  end subroutine ! nudging_set_profile
  !================================================================


  !================================================================
  real(r8) function nudging_set_PSprofile(rlat,rlon,Nudge_PSprof)
   !
   ! NUDGING_SET_PSPROFILE: for the given lat and lon set the surface
   !                      pressure profile value for the specified index.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   real(r8) rlat,rlon
   integer  Nudge_PSprof

   call endrun('nudging_set_PSprofile: Nudging not supported with weak scaling fix')
  end function ! nudging_set_PSprofile
  !================================================================


  !================================================================
  subroutine calc_DryStaticEnergy(t, q, phis, ps, dse, ncol)
   !
   ! calc_DryStaticEnergy: Given the temperature, specific humidity, surface pressure,
   !                       and surface geopotential for a chunk containing 'ncol' columns,
   !                       calculate and return the corresponding dry static energy values.
   !--------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pver, pverp
   use dycore,       only: dycore_is
   use hycoef,       only: hyai, hybi, ps0, hyam, hybm
   use physconst,    only: zvir, gravit, cpair, rair
   !
   ! Input/Output arguments
   !-----------------------
   integer , intent(in) :: ncol      ! Number of columns in chunk
   real(r8), intent(in) :: t(:,:)    ! (pcols,pver) - temperature
   real(r8), intent(in) :: q(:,:)    ! (pcols,pver) - specific humidity
   real(r8), intent(in) :: ps(:)     ! (pcols)      - surface pressure
   real(r8), intent(in) :: phis(:)   ! (pcols)      - surface geopotential
   real(r8), intent(out):: dse(:,:)  ! (pcols,pver)  - dry static energy
   !
   call endrun('calc_DryStaticEnergy: Nudging not supported with weak scaling fix')
  end subroutine calc_DryStaticEnergy
  !================================================================

end module nudging
