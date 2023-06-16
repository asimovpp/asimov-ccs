!> Test case for source terms
!

program test_sources

  use testing_lib

  use types, only: ccs_matrix, ccs_mesh, ccs_vector

  use meshing, only: get_centre, get_local_num_cells, &
       create_cell_locator
  use vec, only: get_vector_data, restore_vector_data

  use utils, only: zero

  implicit none

  ! Mesh / geometry information
  integer(ccs_int), parameter :: n = 5
  real(ccs_real), parameter :: l = 1.0_ccs_real
  type(ccs_mesh) :: mesh

  ! Linear system
  class(ccs_vector), allocatable :: rhs ! Right hand side vector
  class(ccs_vector), allocatable :: x   ! Solution vector
  class(ccs_vector), allocatable :: S   ! Source vector
  class(ccs_matrix), allocatable :: M   ! System matrix
  
  call init()

  call init_case()
  call test_fixed_source()
  
  call fin()

contains

  !v Tests the addition of a fixed source term to the righthand side (RHS) vector.
  !
  !  Using the finite volume method, the expectation is that the RHS vector will contain the source
  !  term multiplied by the cell volumes after this operation.
  subroutine test_fixed_source()

    use fv, only: add_fixed_source
    
    use meshing, only: get_volume
    
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: index_p
    type(cell_locator) :: loc_p
    real(ccs_real), dimension(3) :: x_p
    real(ccs_real) :: V_p

    real(ccs_real), dimension(:), pointer :: rhs_data
    real(ccs_real) :: rhs_expect
    
    call zero(rhs) ! Just a precaution / to simplify the error check.

    call add_fixed_source(mesh, S, rhs)

    call get_vector_data(rhs, rhs_data)
    call get_local_num_cells(mesh, local_num_cells)
    do index_p = 1, local_num_cells
       call create_cell_locator(mesh, index_p, loc_p)
       call get_centre(loc_p, x_p)
       call get_volume(loc_p, V_p)

       rhs_expect = compute_source(x_p) * V_p
       call assert_eq(rhs_data(index_p), rhs_expect, "BLAH")
    end do
    call restore_vector_data(rhs, rhs_data)
    
  end subroutine test_fixed_source
  
  !v Creates the system matrix and vectors, setting them to some initial values.
  subroutine init_case()

    use types, only: matrix_spec, vector_spec

    use vec, only: create_vector
    use mat, only: create_matrix
    
    use utils, only: set_size, initialise, update
    use mesh_utils, only: build_mesh
    
    type(matrix_spec) :: mat_sizes
    type(vector_spec) :: vec_sizes

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: index_p
    type(cell_locator) :: loc_p
    real(ccs_real), dimension(3) :: x_p

    real(ccs_real), dimension(:), pointer :: x_data
    real(ccs_real), dimension(:), pointer :: S_data
    
    ! Initialise mesh
    mesh = build_mesh(par_env, n, n, n, l)

    ! Initialise vectors
    call initialise(vec_sizes)
    call set_size(par_env, mesh, vec_sizes)

    call create_vector(vec_sizes, rhs)
    call create_vector(vec_sizes, x)
    call create_vector(vec_sizes, S)

    call zero(rhs)
    call zero(x)
    call zero(S)

    call get_vector_data(x, x_data)
    call get_vector_data(S, S_data)
    call get_local_num_cells(mesh, local_num_cells)
    do index_p = 1, local_num_cells
       call create_cell_locator(mesh, index_p, loc_p)
       call get_centre(loc_p, x_p)

       x_data(index_p) = set_solution(x_p)
       S_data(index_p) = compute_source(x_p)
    end do
    call restore_vector_data(x, x_data)
    call restore_vector_data(S, S_data)

    call update(rhs)
    call update(x)
    call update(S)
    
    ! Initialise matrix
    call initialise(mat_sizes)
    call set_size(par_env, mesh, mat_sizes)

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
