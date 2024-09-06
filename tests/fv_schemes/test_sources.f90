!> Test case for source terms
!

program test_sources

  use testing_lib

  use ccs_base, only: bnd_names_default
  use types, only: ccs_matrix, ccs_mesh, ccs_vector

  use meshing, only: get_centre, get_local_num_cells, &
       create_cell_locator, get_volume
  use vec, only: get_vector_data, restore_vector_data
  use meshing, only: set_mesh_object, nullify_mesh_object

  use utils, only: zero

  implicit none

  ! Mesh / geometry information
  integer(ccs_int), parameter :: n = 5
  real(ccs_real), parameter :: l = 1.0_ccs_real

  ! Linear system
  class(ccs_vector), allocatable :: rhs ! Right hand side vector
  class(ccs_vector), allocatable :: x   ! Solution vector
  class(ccs_vector), allocatable :: S   ! Source vector
  class(ccs_matrix), allocatable :: M   ! System matrix
  
  call init()

  call init_case()
  call test_fixed_source()
  call test_linear_source()

  call nullify_mesh_object()
  
  call fin()

contains

  !v Tests the addition of a fixed source term to the righthand side (RHS) vector.
  !
  !  Using the finite volume method, the expectation is that the RHS vector will contain the source
  !  term multiplied by the cell volumes after this operation.
  subroutine test_fixed_source()

    use fv, only: add_fixed_source
    
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: index_p
    type(cell_locator) :: loc_p
    real(ccs_real), dimension(3) :: x_p
    real(ccs_real) :: V_p

    real(ccs_real), dimension(:), pointer :: rhs_data
    real(ccs_real) :: rhs_expect
    
    call zero(rhs) ! Just a precaution / to simplify the error check.

    call add_fixed_source(S, rhs)

    call get_vector_data(rhs, rhs_data)
    call get_local_num_cells(local_num_cells)
    do index_p = 1, local_num_cells
       call create_cell_locator(index_p, loc_p)
       call get_centre(loc_p, x_p)
       call get_volume(loc_p, V_p)

       rhs_expect = compute_source(x_p) * V_p
       call assert_eq(rhs_data(index_p), rhs_expect, "RHS does not correspond to adding a fixed source term!")
    end do
    call restore_vector_data(rhs, rhs_data)
    
  end subroutine test_fixed_source

  !v Tests the addition of a linear source term to the system matrix.
  !
  !  Using the finite volume method, the expectation is that the matrix diagonal will contain the
  !  source term multiplied by the cell volumes after this operation.
  !  This is tested by performing a matvec with a vector of 1s, the result vector should then
  !  contain the diagonal coeficient which can be compared with the expected value.
  subroutine test_linear_source()

    use mat, only: mat_vec_product, finalise_matrix

    use fv, only: add_linear_source
    
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: index_p
    type(cell_locator) :: loc_p
    real(ccs_real), dimension(3) :: x_p
    real(ccs_real) :: V_p

    real(ccs_real), dimension(:), pointer :: rhs_data
    real(ccs_real) :: rhs_expect

    call zero(rhs) ! Just a precaution / to simplify the error check.
    call zero(M)

    call add_linear_source(S, M)

    call finalise_matrix(M)
    
    !! Compute RHS
    call mat_vec_product(M, x, rhs)
    
    call get_vector_data(rhs, rhs_data)
    call get_local_num_cells(local_num_cells)
    do index_p = 1, local_num_cells
       call create_cell_locator(index_p, loc_p)
       call get_centre(loc_p, x_p)
       call get_volume(loc_p, V_p)

       rhs_expect = compute_source(x_p) * set_solution(x_p) * V_p
       call assert_eq(rhs_data(index_p), rhs_expect, "RHS does not correspond to mat-vec with linear source term!")
    end do
    call restore_vector_data(rhs, rhs_data)

  end subroutine test_linear_source
  
  !v Creates the system matrix and vectors, setting them to some initial values.
  subroutine init_case()

    use constants, only: cell
    use types, only: matrix_spec, vector_spec

    use vec, only: create_vector, set_vector_location
    use mat, only: create_matrix, set_nnz
    
    use utils, only: set_size, initialise, update
    use mesh_utils, only: build_mesh
    
    type(matrix_spec) :: mat_sizes
    type(vector_spec) :: vec_sizes

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: index_p
    type(cell_locator) :: loc_p
    real(ccs_real), dimension(3) :: x_p
    real(ccs_real) :: V_P

    real(ccs_real), dimension(:), pointer :: x_data
    real(ccs_real), dimension(:), pointer :: S_data
    
    ! Initialise mesh
    mesh = build_mesh(par_env, shared_env, n, n, n, l, &
         bnd_names_default)
    call set_mesh_object(mesh)

    ! Initialise vectors
    call initialise(vec_sizes)
    call set_vector_location(cell, vec_sizes)
    call set_size(par_env, mesh, vec_sizes)

    call create_vector(vec_sizes, rhs)
    call create_vector(vec_sizes, x)
    call create_vector(vec_sizes, S)

    call zero(rhs)
    call zero(x)
    call zero(S)

    call update(rhs)
    call update(x)
    call update(S)
    
    call get_vector_data(x, x_data)
    call get_vector_data(S, S_data)
    call get_local_num_cells(local_num_cells)
    do index_p = 1, local_num_cells
       call create_cell_locator(index_p, loc_p)
       call get_centre(loc_p, x_p)
       call get_volume(loc_p, V_P)

       x_data(index_p) = set_solution(x_p)
       S_data(index_p) = compute_source(x_p) * V_P ! We need to pass the integrated source
    end do
    call restore_vector_data(x, x_data)
    call restore_vector_data(S, S_data)

    call update(rhs)
    call update(x)
    call update(S)
    
    ! Initialise matrix
    call initialise(mat_sizes)
    call set_size(par_env, mesh, mat_sizes)
    call set_nnz(1, mat_sizes) ! Only need room to store the diagonal
    call create_matrix(mat_sizes, M)

    call zero(M)
    
  end subroutine init_case

  real(ccs_real) function compute_source(x)

    real(ccs_real), dimension(3), intent(in) :: x

    compute_source = 1.0_ccs_real + sqrt(sum(x))

  end function compute_source

  real(ccs_real) function set_solution(x)

    real(ccs_real), dimension(3), intent(in) :: x

    set_solution = 1.0_ccs_real + 0 * sqrt(sum(x))

  end function set_solution
  
end program test_sources
