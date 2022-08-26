!> @brief Test the cell/face centres of a square mesh.
!
!> @description The cell/face centres of a mesh should all fall within the meshed domain, for a
!!              square mesh \f$x\in[0,1]^d\f$.
program test_vec_set_entries

  use testing_lib

  use kinds
  use constants, only : insert_mode, add_mode
  use types, only : ccs_vector, ccs_mesh
  use mesh_utils, only : global_start, local_count, build_square_mesh
  use utils, only : zero
  
  implicit none

  class(ccs_vector), allocatable :: v
  integer(ccs_int) :: n
  type(ccs_mesh) :: mesh
  real(ccs_real), parameter :: elt_val = 1.0_ccs_real
  
  call init()

  do n = 1, 100
    mesh = build_square_mesh(par_env, n, 1.0_ccs_real)

    call init_vector()

    call zero(v)
    call set_vector(add_mode)
    call test_vector()
    call set_vector(insert_mode)
    call test_vector()
    
    call clean_vector()
  end do
  
  call fin()

contains

  subroutine init_vector()

    use utils, only : set_size, initialise
    use vec, only : create_vector

    type(vector_spec) :: vec_sizes

    call initialise(vec_sizes)
    call set_size(par_env, mesh, vec_sizes)

    call create_vector(vec_sizes, v)
  end subroutine init_vector

  subroutine set_vector(mode)

    use types, only : cell_locator
    use utils, only : set_mode, set_row, set_entry, set_values, clear_entries, update
    use vec, only : create_vector_values
    use meshing, only : set_cell_location, get_global_index
    
    integer(ccs_int), intent(in) :: mode

    integer(ccs_int) :: nblocks !> How many blocks should I split my elements into?
    integer(ccs_int) :: nrows   !> How many rows to set simultaneously?

    type(vector_values) :: val_dat
    type(cell_locator) :: loc_p
    
    integer(ccs_int) :: i
    integer(ccs_int) :: j

    integer(ccs_int) :: index_p
    integer(ccs_int) :: global_index_p

    associate(nlocal => mesh%topo%local_num_cells)
    
      nrows = 1_ccs_int
      nblocks = nlocal / nrows

      call create_vector_values(nrows, val_dat)
      call set_mode(mode, val_dat)
    
      do i = 1_ccs_int, nblocks
        call clear_entries(val_dat)
        
        do j = 1_ccs_int, nrows
          ! Compute local and global indices
          index_p = j + (i - 1) * nrows

          call set_cell_location(mesh, index_p, loc_p)
          call get_global_index(loc_p, global_index_p)

          call set_row(global_index_p, val_dat) ! TODO: this should work on local indices...
          call set_entry(elt_val, val_dat)
        end do
      
        call set_values(val_dat, v) ! TODO: this should support setting multiple value simultaneously
      end do
      ! TODO: remainder loop (required for blocksize > 1)...

    end associate
    
    call update(v)
    
  end subroutine set_vector

  subroutine test_vector()

    use solver, only : norm

    real(ccs_real) :: expectation
    real(ccs_real) :: test_value
    
    expectation = sqrt(real(mesh%topo%global_num_cells, ccs_real))
    test_value = norm(v, 2)
    if (test_value /= expectation) then
      print *, "AAARGH, FAIL: ", test_value, expectation
      stop
    end if
    
  end subroutine test_vector
  
  subroutine clean_vector()

    deallocate(v)
    
  end subroutine clean_vector
  
end program test_vec_set_entries
