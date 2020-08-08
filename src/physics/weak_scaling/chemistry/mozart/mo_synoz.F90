
module mo_synoz

  !--------------------------------------------------------------------
  !	... synthetic stratospheric ozone emission source
  ! Shadow version for weak scaling patch
  !--------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog

  implicit none

  save

  real(r8), allocatable :: po3(:,:,:)

  private
  public :: synoz_inti
  public :: po3

contains

  subroutine synoz_inti( )
    !-----------------------------------------------------------------------
    ! 	... initialize synoz emissions
    !	    note: the emissions are in in units of molecules/cm**3/s
    !-----------------------------------------------------------------------

    use cam_abortutils, only: endrun

    call endrun('synoz_inti: synoz not supported with weak scaling')

  end subroutine synoz_inti

end module mo_synoz
