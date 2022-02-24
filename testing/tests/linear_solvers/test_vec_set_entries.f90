!> @brief Test the cell/face centres of a square mesh.
!
!> @description The cell/face centres of a mesh should all fall within the meshed domain, for a
!!              square mesh \f$x\in[0,1]^d\f$.
program test_vec_set_entries

  use testing_lib

  use kinds
  use constants, only : insert_mode, add_mode
  use types, only : vector
  
  implicit none

  class(vector), allocatable :: v
  integer(accs_int) :: n
  real(accs_real), parameter :: elt_val = 1.0_accs_real
  
  call init()

  do n = 1, 100
    call init_vector(n)

    call set_vector(n, add_mode)
    call test_vector(n)
    call set_vector(n, insert_mode)
    call test_vector(n)
    
    call clean_vector()
  end do
  
  call fin()

contains

  subroutine init_vector(n)

    use utils, only : set_global_size, initialise
    use vec, only : create_vector
    
    integer(accs_int), intent(in) :: n

    type(vector_init_data) :: vec_sizes

    call initialise(vec_sizes)
    call set_global_size(vec_sizes, n, par_env)

    call create_vector(vec_sizes, v)
    
  end subroutine init_vector

  subroutine set_vector(n, mode)

    use utils, only : set_mode, set_row, set_entry, set_values, clear_entries, update
    use vec, only : create_vector_values
    
    integer(accs_int), intent(in) :: n
    integer(accs_int), intent(in) :: mode

    integer(accs_int) :: nblocks !> How many blocks should I split my elements into?
    integer(accs_int) :: nlocal  !> How many elements do I own?
    integer(accs_int) :: nrows   !> How many rows to set simultaneously?

    type(vector_values) :: val_dat
    
    integer(accs_int) :: i
    integer(accs_int) :: j

    integer(accs_int) :: idxl
    integer(accs_int) :: idxg
    integer(accs_int) :: offset

    offset = global_start(n, par_env%proc_id, par_env%num_procs) - 1
    nlocal = local_count(n, par_env%proc_id, par_env%num_procs)
    
    nrows = 1_accs_int
    nblocks = nlocal / nrows
    
    call create_vector_values(nrows, val_dat)
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
      
      call set_values(val_dat, v) ! TODO: this should support setting multiple value simultaneously
    end do
    ! TODO: remainder loop (required for blocksize > 1)...

    call update(v)
    
  end subroutine set_vector

  subroutine test_vector(n)

    use vec, only : norm
    
    integer(accs_int), intent(in) :: n

    real(accs_real) :: expectation

    expectation = sqrt(real(n, accs_real))
    if (norm(v, 2) /= expectation) then
      stop
    end if
    
  end subroutine test_vector
  
  subroutine clean_vector()

    deallocate(v)
    
  end subroutine clean_vector
  
  integer function global_start(n, procid, nproc)

    integer(accs_int), intent(in) :: n
    integer(accs_int), intent(in) :: procid
    integer(accs_int), intent(in) :: nproc

    !! Each PE gets an equal split of the problem with any remainder split equally between the lower
    !! PEs.
    global_start = procid * (n / nproc) + min(procid, modulo(n, nproc))

    !! Fortran indexing
    global_start = global_start + 1
    
  end function global_start

  integer function local_count(n, procid, nproc)

    integer(accs_int), intent(in) :: n
    integer(accs_int), intent(in) :: procid
    integer(accs_int), intent(in) :: nproc

    if (procid < n) then
      local_count = global_start(n, procid, nproc)
      if (procid < (nproc - 1)) then
        local_count = global_start(n, procid + 1, nproc) - local_count
      else
        local_count = n - (local_count - 1)
      end if
    else
      local_count = 0
    end if
    
  end function local_count
  
end program test_vec_set_entries
