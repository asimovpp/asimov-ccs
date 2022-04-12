program test_mat_set_entries

  use testing_lib

  use kinds
  use constants, only : insert_mode, add_mode
  use types, only : ccs_matrix, ccs_mesh
  use mesh_utils, only : global_start, local_count, build_square_mesh
  
  implicit none

  class(ccs_matrix), allocatable :: M
  integer(ccs_int) :: n
  type(ccs_mesh) :: mesh
  real(ccs_real), parameter :: elt_val = 1.0_ccs_real
  
  call init()

  do n = 1, 100
    mesh = build_square_mesh(par_env, n, 1.0_ccs_real)
    
    call init_matrix()

    call set_matrix(n, add_mode)
    call set_matrix(n, insert_mode)
    
    call clean_matrix()
  end do
  
  call fin()

contains

  subroutine init_matrix()

    use utils, only : set_size, initialise
    use mat, only : create_matrix

    type(matrix_spec) :: mat_sizes

    call initialise(mat_sizes)
    call set_size(par_env, mesh, mat_sizes)

    call create_matrix(mat_sizes, M)
    
  end subroutine init_matrix

  subroutine set_matrix(n, mode)

    use utils, only : set_mode, set_row, set_entry, set_values, clear_entries, update
    use mat, only : create_matrix_values
    
    integer(ccs_int), intent(in) :: n
    integer(ccs_int), intent(in) :: mode

    integer(ccs_int) :: nblocks !> How many blocks should I split my elements into?
    integer(ccs_int) :: nlocal  !> How many elements do I own?
    integer(ccs_int) :: nrows   !> How many rows to set simultaneously?

    type(matrix_values) :: val_dat
    
    integer(ccs_int) :: i
    integer(ccs_int) :: j

    integer(ccs_int) :: idxl
    integer(ccs_int) :: idxg
    integer(ccs_int) :: offset

    offset = global_start(n, par_env%proc_id, par_env%num_procs) - 1
    nlocal = local_count(n, par_env%proc_id, par_env%num_procs)
    
    nrows = 1_ccs_int
    nblocks = nlocal / nrows
    
    ! print*,"On rank ", par_env%proc_id," offset is ", offset
    ! print*,"On rank ", par_env%proc_id," nlocal is ", nlocal
    ! print*,"On rank ", par_env%proc_id," nrows is ", nrows
    ! print*,"On rank ", par_env%proc_id," nblocks is ", nblocks

    call create_matrix_values(nrows, val_dat)
    call set_mode(mode, val_dat)
    
    do i = 1_ccs_int, nblocks
      call clear_entries(val_dat)

      do j = 1_ccs_int, nrows
        ! Compute local and global indices
        idxl = j + (i - 1) * nrows
        idxg = idxl + offset

        ! print*,"On rank ", par_env%proc_id," idxl is ", idxl
        ! print*,"On rank ", par_env%proc_id," idxg is ", idxg
        ! print*,"On rank ", par_env%proc_id," setting row"
        call set_row(idxl, val_dat) ! TODO: this should work on local indices...
        ! print*,"On rank ", par_env%proc_id," setting entry"
        call set_entry(elt_val, val_dat)
      end do
      
      call set_values(val_dat, M) ! TODO: this should support setting multiple value simultaneously
    end do
    ! ! TODO: remainder loop (required for blocksize > 1)...

    call update(M)
    
  end subroutine set_matrix
  
  subroutine clean_matrix()

    deallocate(M)
    
  end subroutine clean_matrix
  
end program test_mat_set_entries
