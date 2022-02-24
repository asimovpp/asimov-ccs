program test_mat_set_entries

  use testing_lib

  use kinds
  use constants, only : insert_mode, add_mode
  use types, only : matrix
  use mesh_utils, only : global_start, local_count
  
  implicit none

  class(matrix), allocatable :: M
  integer(accs_int) :: n
  real(accs_real), parameter :: elt_val = 1.0_accs_real
  
  call init()

  do n = 1, 100
    call init_matrix(n)

    call set_matrix(n, add_mode)
    call set_matrix(n, insert_mode)
    
    call clean_matrix()
  end do
  
  call fin()

contains

  subroutine init_matrix(n)

    use utils, only : set_global_size, initialise
    use mat, only : create_matrix
    
    integer(accs_int), intent(in) :: n

    type(matrix_init_data) :: mat_sizes

    call initialise(mat_sizes)
    call set_global_size(mat_sizes, n, n, par_env)

    call create_matrix(mat_sizes, M)
    
  end subroutine init_matrix

  subroutine set_matrix(n, mode)

    use utils, only : set_mode, set_row, set_entry, set_values, clear_entries, update
    use mat, only : create_matrix_values
    
    integer(accs_int), intent(in) :: n
    integer(accs_int), intent(in) :: mode

    integer(accs_int) :: nblocks !> How many blocks should I split my elements into?
    integer(accs_int) :: nlocal  !> How many elements do I own?
    integer(accs_int) :: nrows   !> How many rows to set simultaneously?

    type(matrix_values) :: val_dat
    
    integer(accs_int) :: i
    integer(accs_int) :: j

    integer(accs_int) :: idxl
    integer(accs_int) :: idxg
    integer(accs_int) :: offset

    offset = global_start(n, par_env%proc_id, par_env%num_procs) - 1
    nlocal = local_count(n, par_env%proc_id, par_env%num_procs)
    
    nrows = 1_accs_int
    nblocks = nlocal / nrows
    
    call create_matrix_values(nrows, val_dat)
    call set_mode(mode, val_dat)
    
    do i = 1_accs_int, nblocks
      call clear_entries(val_dat)

      do j = 1_accs_int, nrows
        ! Compute local and global indices
        idxl = j + (i - 1) * nrows
        idxg = idxl + offset

        call set_row(idxg, val_dat) ! TODO: this should work on local indices...
        call set_entry(elt_val, val_dat)
      end do
      
      call set_values(val_dat, M) ! TODO: this should support setting multiple value simultaneously
    end do
    ! TODO: remainder loop (required for blocksize > 1)...

    call update(M)
    
  end subroutine set_matrix
  
  subroutine clean_matrix()

    deallocate(M)
    
  end subroutine clean_matrix
  
end program test_mat_set_entries
