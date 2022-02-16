!> @brief Test the cell/face centres of a square mesh.
!
!> @description The cell/face centres of a mesh should all fall within the meshed domain, for a
!!              square mesh \f$x\in[0,1]^d\f$.
program test_square_mesh_centres

  use testing_lib

  use kinds
  use types, only : vector
  
  implicit none

  class(vector), allocatable :: v
  integer(accs_int) :: n
  real(accs_real), parameter :: elt_val = 1.0_accs_real
  
  call init()

  do n = 1, 100
    call init_vector(n)
    call set_vector(n)
    call test_vetor(n)
    call clean_vector()
  end do
  
  call fin()

contains

  subroutine init_vector(n)

    use utils, only : set_global_size
    use vec, only : create_vector
    
    integer(accs_int), intent(in) :: n

    type(vector_init_data) :: vec_sizes

    call initialise(vec_sizes)
    call set_global_size(vec_sizes, n, par_env)

    call create_vector(vec_sizes, v)
    
  end subroutine init_vector

  subroutine set_vector(n)

    integer(accs_int), intent(in) :: n

    integer(accs_int) :: nblocks !> How many blocks should I split my elements into?
    integer(accs_int) :: nlocal  !> How many elements do I own?
    integer(accs_int) :: nrows   !> How many rows to set simultaneously?

    integer(accs_int) :: i
    integer(accs_int) :: j
    
    nrows = 1_accs_int
    nblocks = nlocal / nrows
    
    call create_vector_values(nrows, val_dat)
    call set_mode(insert_mode, val_dat)
    
    do i = 1_accs_int, nblocks
      do j = 1_accs_int, nrows
        idxg = i + offset

        call set_row(idxg, val_dat) ! TODO: this should work on local indices...
        call set_val(elt_val, val_dat)

        call set_values(val_dat, v) ! TODO: this should support setting multiple value simultaneously
      end do
    end do

    deallocate(val_dat)
    
  end subroutine set_vector

  subroutine test_vector(n)

    use vec, only : norm
    
    integer(accs_int), intent(in) :: n

    real(accs_real) :: expectation

    expectation = real(n, accs_real)
    if (norm(v, 2) /= expectation) then
      stop
    end if
    
  end subroutine test_vector
  
  subroutine clean_vector()

    deallocate(v)
    
  end subroutine clean_vector
  
end program test_square_mesh_centres
