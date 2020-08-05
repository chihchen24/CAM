module gmean_mod
   !-----------------------------------------------------------------------
   !
   ! Purpose:
   ! Perform mixed layer global calculations for energy conservation checks.
   !
   ! Method:
   ! Reproducible (semi-scalable):
   !    Gather to a master processor which sums all columns in order, or,
   !    gather to a few processors, each of which sums its block in order
   !    then sends its sum to the master processor for a final sum.
   !    Chunking is controlled by the namelist variable, <global_max_fields>.
   !    Sums should be reproducible as long as the value of this variable
   !    does not change.
   !
   !-----------------------------------------------------------------------
   use shr_kind_mod,  only: r8 => shr_kind_r8
   use perf_mod,      only: t_startf, t_stopf
   use spmd_utils,    only: MPI_COMM_NULL, column_redist_t
   use cam_logfile,   only: iulog

   implicit none
   private
   save

   public :: gmean_init ! Initialize reproducible sum algorithm (if enabled)
   public :: gmean ! compute global mean of 2D fields on physics decomposition
   public :: gmean_finalize

   interface gmean
      module procedure gmean_arr
      module procedure gmean_scl
   end interface gmean

   ! Private data
   logical                       :: initialized = .false.
   logical                       :: use_repro_sum = .false.
   integer                       :: max_nflds = -1 ! # fields for disp. decomp
   type(column_redist_t)         :: column_reorder

CONTAINS

   !
   !========================================================================
   !

   subroutine gmean_init(max_nflds_in)
      use phys_control,   only: phys_getopts
      use cam_abortutils, only: endrun
      use spmd_utils,     only: masterproc, mpicom, npes, iam
      use spmd_utils,     only: MPI_INTEGER, MPI_REAL8, MPI_SUM, MPI_MAX
      use ppgrid,         only: pcols, begchunk, endchunk
      use phys_grid,      only: global_max_fields => phys_global_max_fields
      use phys_grid,      only: num_global_phys_cols, columns_on_task
      use phys_grid,      only: init_col_assem_p
      !-----------------------------------------------------------------------
      !
      ! Purpose:
      ! Pre-compute layout and communication information for efficient
      !  run-time compuatation of global sums
      !
      !-----------------------------------------------------------------------
      !
      ! Dummy arguments
      !
      integer,          intent(in)  :: max_nflds_in ! maximum number of fields
      ! Local variables
      !
      integer                       :: total_size
      integer                       :: num_blocks
      integer                       :: num_sum_cols
      integer                       :: max_block_size
      integer                       :: findex
      integer                       :: dest_spacing
      integer                       :: col_start, col_end
      integer                       :: color
      integer                       :: ierr
      integer                       :: block_sum_comm
      integer,          allocatable :: send_cnts(:)
      integer,          allocatable :: col_info(:)
      character(len=256)            :: errmsg
      integer,          parameter   :: max_recv = 16 ! max for gather
      integer,          parameter   :: max_pcnt = 3  ! max ratio for gather
      integer,          parameter   :: default_max_fields = 100000
      character(len=*), parameter   :: subname = 'gmean_init: '
      !
      !-----------------------------------------------------------------------
      !
      call phys_getopts(cam_repro_sum_out=use_repro_sum)
      if (.not. use_repro_sum) then
         call gmean_finalize()
         return
      end if
      if (initialized) then
         if (max_nflds /= max_nflds_in) then
            call gmean_finalize()
         else
            return ! Should not need to do anything
         end if
      end if
      ! Perform reproducible sum. Note, this depends on the value
      ! of <global_max_fields> not changing between runs
      total_size = num_global_phys_cols * max_nflds_in
      ! First try
      if (global_max_fields <= 0) then
         ! max global size not configured, go for approximately tot size / 10^5
         num_blocks = MAX(total_size / default_max_fields, 1)
      else
         num_blocks = MAX(total_size / global_max_fields, 1)
      end if
      ! We need to be able to process at least one entire field at a time
      !   even if this violates the <global_max_fields> limit. This algorithm
      !   requires enough memory for at least one additional field.
      do
         num_sum_cols = total_size / (num_blocks * max_nflds_in)
         max_block_size = MIN(num_sum_cols * max_nflds_in, global_max_fields)
         if ((num_blocks * max_block_size) < total_size) then
            num_blocks = num_blocks + 1
         else
            exit
         end if
      end do
      column_reorder%max_nflds = max_nflds_in
      if (num_blocks > npes) then
         ! Will need to stage sums in chunks if nflds gets too high
         num_blocks = npes
         do
            column_reorder%max_nflds = MAX(column_reorder%max_nflds - 1, 1)
            total_size = num_global_phys_cols * column_reorder%max_nflds
            num_sum_cols = total_size / (num_blocks * column_reorder%max_nflds)
            max_block_size = num_sum_cols * column_reorder%max_nflds
            if (max_block_size <= global_max_fields) then
               exit
            else if (column_reorder%max_nflds == 1) then
               exit
            end if
         end do
      end if
      ! If we get here, we should not have any allocated arrays, check.
      if (allocated(column_reorder%dest_tasks)) then
         call endrun(subname//'dest_tasks allocated (should not be)')
      end if
      if (associated(column_reorder%col_starts)) then
         call endrun(subname//'col_starts allocated (should not be)')
      end if
      if (allocated(column_reorder%task_sizes)) then
         call endrun(subname//'task_sizes allocated (should not be)')
      end if
      if (allocated(column_reorder%recv_cnts)) then
         call endrun(subname//'recv_cnts allocated (should not be)')
      end if
      if (allocated(column_reorder%recv_disps)) then
         call endrun(subname//'recv_disps allocated (should not be)')
      end if
      if (allocated(column_reorder%recv_reorder)) then
         call endrun(subname//'recv_reorder allocated (should not be)')
      end if
      allocate(column_reorder%dest_tasks(0:num_blocks-1))
      allocate(column_reorder%col_starts(0:num_blocks-1))
      dest_spacing = MAX(npes / num_blocks, 1)
      do findex = 0, num_blocks - 1
         column_reorder%dest_tasks(findex) = dest_spacing * findex ! < npes
         column_reorder%col_starts(findex) = (num_sum_cols * findex) + 1
      end do
      ! Get redistribution info from physics grid
      call init_col_assem_p(column_reorder)
      ! Create communicator for receiving tasks
      block_sum_comm = mpicom
      if (ANY(column_reorder%dest_tasks == iam)) then
         color = 1
      else
         color = 2
      end if
      if (size(column_reorder, 1) < npes) then
         call MPI_Comm_split(mpicom, color, iam, block_sum_comm)
      end if
      if (color /= 1) then
         call mpi_comm_free(block_sum_comm, ierr)
         block_sum_comm = MPI_COMM_NULL
      end if
      column_reorder%mpi_comm = block_sum_comm
      ! Retrieve number of columns (column_reorder%task_sizes) for each block
      allocate(column_reorder%recv_cnts(0:npes-1))
      column_reorder%recv_cnts = 0
      allocate(send_cnts(0:npes-1))
      ! Initialize send_cnts to the send count for each destination task
      send_cnts = 0
      send_cnts(column_reorder%dest_tasks(:)) = column_reorder%task_sizes(:)
      if ((num_blocks >= max_recv) .or. ((npes / num_blocks) < max_pcnt)) then
         call MPI_Alltoall(send_cnts, 1, MPI_INTEGER,                         &
              column_reorder%recv_cnts, 1, MPI_INTEGER, mpicom, ierr)
      else
         do findex = 0, num_blocks - 1
            call MPI_gather(column_reorder%task_sizes(findex), 1,             &
                 MPI_INTEGER, column_reorder%recv_cnts, 1, MPI_INTEGER,       &
                 column_reorder%dest_tasks(findex), mpicom, ierr)
         end do
      end if
      ! Compute displacements from counts
      allocate(column_reorder%send_disps(0:npes - 1))
      column_reorder%send_disps(0) = 0
      allocate(column_reorder%recv_disps(0:npes - 1))
      column_reorder%recv_disps(0) = 0
      do findex = 1, npes - 1
         column_reorder%send_disps(findex) =                                  &
              column_reorder%send_disps(findex - 1) + send_cnts(findex - 1)
         column_reorder%recv_disps(findex) =                                  &
              column_reorder%recv_disps(findex - 1) +                         &
              column_reorder%recv_cnts(findex - 1)
      end do
      ! Compute the summation order for each task with a block of data
      allocate(col_info(SUM(column_reorder%recv_cnts)))
      allocate(column_reorder%recv_reorder(size(col_info)))
      col_info = 0
      call MPI_alltoallv(column_reorder%task_indices, send_cnts,              &
           column_reorder%send_disps, MPI_INTEGER, col_info,                  &
           column_reorder%recv_cnts, column_reorder%recv_disps,               &
           MPI_INTEGER, mpicom, ierr)
      ! col_info now has the global column index for each received quantity
      do findex = 1, size(col_info)
         ! color is the location where this column belongs in summation order
         color = col_info(findex) - column_reorder%col_starts(iam) + 1
         if (color < 1) then
            write(errmsg, *) 'recv_order underflow on ',findex,', ',color
            call endrun(subname//errmsg)
         else if (color > size(column_reorder%recv_reorder)) then
            write(errmsg, *) 'recv_order overflow on ',findex,', ',color
            call endrun(subname//errmsg)
         end if
         column_reorder%recv_reorder(color) = findex
      end do
      deallocate(col_info)
      ! Done !
      max_nflds = max_nflds_in
      initialized = .true.
      ! Cleanup
      deallocate(send_cnts)
      deallocate(col_info)
   end subroutine gmean_init

   subroutine gmean_arr(arr, arr_gmean, nflds)
      use cam_abortutils, only: endrun
      use spmd_utils,     only: masterprocid, mpicom, npes, iam
      use spmd_utils,     only: MPI_REAL8, MPI_SUM
      use ppgrid,         only: pcols, begchunk, endchunk
      use phys_grid,      only: global_max_fields => phys_global_max_fields
      use phys_grid,      only: num_global_phys_cols, columns_on_task
      use phys_grid,      only: weighted_sum_p, weighted_field_p
      !-----------------------------------------------------------------------
      !
      ! Purpose:
      ! Compute the global mean of each field in "arr" in the physics
      ! chunked decomposition
      !
      !-----------------------------------------------------------------------
      !
      ! Dummy arguments
      !
      integer,           intent(in)  :: nflds            ! number of fields
      real(r8),          intent(in)  :: arr(pcols, begchunk:endchunk, nflds)
      real(r8),          intent(out) :: arr_gmean(nflds) ! global means
      !
      ! Local variables
      !
      integer                       :: fld_start, fld_end, num_flds
      integer                       :: fld_ind, col_ind
      integer                       :: windex
      integer                       :: ierr
      integer,          allocatable :: send_cnts(:), send_disps(:)
      integer,          allocatable :: recv_cnts(:), recv_disps(:)
      real(r8)                      :: wsum(nflds) ! Weighted sums
      real(r8),         allocatable :: weighted_field(:,:)
      real(r8),         allocatable :: sum_arr(:,:)
      character(len=*), parameter   :: subname = 'gmean_arr: '
      !
      !-----------------------------------------------------------------------
      !
      if (.not. initialized) then
         call endrun(subname//'gmean not initialized')
      end if

      arr_gmean = 0.0_r8
      if (use_repro_sum) then
         ! Perform reproducible sum. Note, this depends on the value
         ! of global_max_fields not changing between runs
         allocate(send_cnts(npes))
         send_cnts = 0
         allocate(send_disps(npes))
         allocate(recv_cnts(npes))
         allocate(recv_disps(npes))
         do fld_start = 1, nflds, column_reorder%max_nflds
            fld_end = MIN(nflds, fld_start + column_reorder%max_nflds - 1)
            num_flds = fld_end - fld_start + 1
            send_cnts(column_reorder%dest_tasks(:)) =                         &
                 column_reorder%task_sizes(:) * nflds
            send_disps = column_reorder%send_disps * nflds
            recv_cnts = column_reorder%recv_cnts * nflds
            recv_disps = column_reorder%recv_disps * nflds
            ! Retrieve the raw data
            call weighted_field_p(arr(:,:,fld_start:fld_end), weighted_field)
            call MPI_alltoallv(weighted_field, send_cnts, send_disps,         &
                 MPI_REAL8, sum_arr, recv_cnts, recv_disps,                   &
                 MPI_REAL8, mpicom, ierr)
            ! Sum in global column order
            do col_ind = 1, size(column_reorder%recv_reorder, 1)
               windex = column_reorder%recv_reorder(col_ind)
               do fld_ind = fld_start, fld_end
                  wsum(fld_ind) = wsum(fld_ind) +                             &
                       weighted_field(fld_ind-fld_start+1, windex)
               end do
            end do
         end do
         ! Calculate a global sum (only receiving tasks)
         arr_gmean = wsum
         call MPI_gather(arr_gmean, nflds, MPI_REAL8, wsum, nflds,            &
              MPI_REAL8, masterprocid, column_reorder%mpi_comm, ierr)
         ! Finally, broadcast the result
         if (npes > 1) then
            call mpi_Bcast(arr_gmean, nflds, MPI_REAL8, masterprocid,         &
                 mpicom, ierr)
         end if
         deallocate(send_cnts)
         deallocate(send_disps)
         deallocate(recv_cnts)
         deallocate(recv_disps)
      else
         do fld_ind = 1, nflds
            ! Get weighted sum from physics
            call weighted_sum_p(arr(:, :, fld_ind), wsum(fld_ind))
         end do
         call MPI_allreduce(wsum, arr_gmean, nflds, MPI_REAL8, MPI_SUM,       &
              mpicom, ierr)
      end if
   end subroutine gmean_arr

   !
   !========================================================================
   !

   subroutine gmean_scl (arr, gmean)
      use ppgrid,        only: pcols, begchunk, endchunk
      !-----------------------------------------------------------------------
      !
      ! Purpose:
      ! Compute the global mean of each field in "arr" in the physics
      ! chunked decomposition
      !
      !-----------------------------------------------------------------------
      !
      ! Arguments
      !
      real(r8), intent(in)  :: arr(pcols,begchunk:endchunk)
      ! Input array, chunked
      real(r8), intent(out) :: gmean ! global means
      !
      ! Local workspace
      !
      integer, parameter    :: nflds = 1
      real(r8)              :: gmean_array(nflds)

      call gmean_arr(reshape(arr, (/ pcols, endchunk-begchunk+1, nflds /)),   &
           gmean_array, nflds)
      gmean = gmean_array(1)

   end subroutine gmean_scl

   !
   !========================================================================
   !

   subroutine gmean_finalize()
      !-----------------------------------------------------------------------
      !
      ! Purpose: Deallocate arrays and reset to uninitialized state
      !
      !-----------------------------------------------------------------------

      ! Local variable
      integer :: ierr

      initialized = .false.
      use_repro_sum = .false.
      if (column_reorder%mpi_comm /= MPI_COMM_NULL) then
         call MPI_comm_free(column_reorder%mpi_comm, ierr)
         column_reorder%mpi_comm = MPI_COMM_NULL
      end if
      if (associated(column_reorder%dest_tasks)) then
         deallocate(column_reorder%dest_tasks)
         nullify(column_reorder%dest_tasks)
      end if
      if (associated(column_reorder%col_starts)) then
         deallocate(column_reorder%col_starts)
         nullify(column_reorder%col_starts)
      end if
      if (associated(column_reorder%recv_cnts)) then
         deallocate(column_reorder%recv_cnts)
         nullify(column_reorder%recv_cnts)
      end if
      if (associated(column_reorder%recv_disps)) then
         deallocate(column_reorder%recv_disps)
         nullify(column_reorder%recv_disps)
      end if
      if (associated(column_reorder%recv_reorder)) then
         deallocate(column_reorder%recv_reorder)
         nullify(column_reorder%recv_reorder)
      end if
      if (associated(column_reorder%task_sizes)) then
         deallocate(column_reorder%task_sizes)
         nullify(column_reorder%task_sizes)
      end if
      if (associated(column_reorder%task_indices)) then
         deallocate(column_reorder%task_indices)
         nullify(column_reorder%task_indices)
      end if
      if (associated(column_reorder%send_disps)) then
         deallocate(column_reorder%send_disps)
         nullify(column_reorder%send_disps)
      end if
      if (associated(column_reorder%send_reorder)) then
         deallocate(column_reorder%send_reorder)
         nullify(column_reorder%send_reorder)
      end if
   end subroutine gmean_finalize

end module gmean_mod
