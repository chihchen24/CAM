module polar_avg
!-----------------------------------------------------------------------
!
! Purpose:
!  These routines are used by the fv dycore to set the collocated
!  pole points at the limits of the latitude dimension to the same
!  value.
!
! Methods:
!  The shadow file exists to prevent having to compile physics/cam
!   version which contains calls to unsupported phys_grid interfaces
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
   use shr_kind_mod,  only: r8 => shr_kind_r8

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private
   save

!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public :: &
      polar_average           ! support for LR dycore polar averaging

   interface polar_average
      module procedure polar_average2d, polar_average3d
   end interface

   CONTAINS
!
!========================================================================
!
   subroutine polar_average2d(field)
      use cam_abortutils, only: endrun
      use phys_grid,      only: begchunk, endchunk, pcols
!-----------------------------------------------------------------------
! Purpose: Set the collocated pole points at the limits of the latitude
!          dimension to the same value.
! Author: J. Edwards
!-----------------------------------------------------------------------
!
! Arguments
!
     real(r8), intent(inout) :: field(pcols,begchunk:endchunk)
!
     call endrun('polar_average2d: not supported for weak scaling')

   end subroutine polar_average2d

!
!========================================================================
!

   subroutine polar_average3d(nlev, field)
      use cam_abortutils, only: endrun
      use phys_grid,      only: begchunk, endchunk, pcols
!-----------------------------------------------------------------------
! Purpose: Set the collocated pole points at the limits of the latitude
!          dimension to the same value.
! Author: J. Edwards
!-----------------------------------------------------------------------
!
! Arguments
!
     integer, intent(in) :: nlev
     real(r8), intent(inout) :: field(pcols,nlev,begchunk:endchunk)

     call endrun('polar_average3d: not supported for weak scaling')

   end subroutine polar_average3d

end module polar_avg
