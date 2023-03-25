!v Program file for Poisson case
!
!  Based on prototype/ex3 a port of PETSc ksp/tutorial/ex3.c to ASiMoV-CCS style code.
!  This case demonstrates setting up a linear system and solving it with ASiMoV-CCS, note
!  the code is independent of PETSc.
!  The example case solves the equation
!  \[
!    {\nabla^2} p = f
!  \]
!  in the unit square with Dirichlet boundary conditions
!  \[
!    p\left(\boldsymbol{x}\right) = y,\ \boldsymbol{x}\in\partial\Omega
!  \]

module problem_setup

  use constants, only : ndim
  use kinds, only : ccs_int, ccs_real
  use types, only : ccs_mesh, cell_locator, face_locator

  use meshing, only : set_face_location, set_cell_location, get_centre
  
  implicit none

  private

  public :: eval_solution
  public :: eval_cell_rhs

  !> Interface to evaluate exact solution.
  interface eval_solution
    module procedure eval_solution_cell
    module procedure eval_solution_face
  end interface eval_solution
  
contains

  !v Evaluate the exact solution.
  !
  !  Used to set the Dirichlet BCs and also the reference solution for testing the numerical
  !  solution. Thus this should reflect changes to the forcing function.
  function eval_solution_coordinates(x) result(r)

    real(ccs_real), dimension(:), intent(in) :: x
    real(ccs_real) :: r

    associate (y => x(2))
      r = y
    end associate

  end function eval_solution_coordinates

  function eval_solution_cell(loc_p) result(r)

    type(cell_locator), intent(in) :: loc_p
    real(ccs_real) :: r

    real(ccs_real), dimension(ndim) :: x

    call get_centre(loc_p, x)
    r = eval_solution_coordinates(x)
    
  end function eval_solution_cell

  function eval_solution_face(loc_f) result(r)

    type(face_locator), intent(in) :: loc_f
    real(ccs_real) :: r

    real(ccs_real), dimension(ndim) :: x

    call get_centre(loc_f, x)
    r = eval_solution_coordinates(x)
    
  end function eval_solution_face

  !> Apply forcing function
  pure subroutine eval_cell_rhs(x, y, H, r)

    real(ccs_real), intent(in) :: x, y, H
    real(ccs_real), intent(out) :: r

    r = 0.0_ccs_real &
        + 0.0_ccs_real * (x + y + H) ! Silence unused dummy argument error

  end subroutine eval_cell_rhs

end module problem_setup

module poisson_discretisation

  use constants, only : ndim, add_mode, insert_mode
  use kinds, only : ccs_int, ccs_real
  use types, only : ccs_vector, vector_values, &
       ccs_matrix, matrix_values_spec, matrix_values, &
       ccs_mesh, cell_locator, neighbour_locator, face_locator

  use mat, only : create_matrix_values, set_matrix_values_spec_ncols, set_matrix_values_spec_nrows
  use meshing, only : get_local_num_cells, set_cell_location, get_centre, get_volume, &
       get_global_index, count_neighbours, set_neighbour_location, get_boundary_status, &
       set_face_location, get_face_area
  use utils, only : clear_entries, set_mode, set_col, set_row, set_entry, set_values
  use vec, only : create_vector_values

  use problem_setup, only : eval_solution
  
  implicit none

  private

  public :: discretise_poisson
  public :: apply_dirichlet_bcs

  interface
    module subroutine discretise_poisson(mesh, M)
      type(ccs_mesh), intent(in) :: mesh    !< The mesh the problem is defined upon
      class(ccs_matrix), intent(inout) :: M !< The system matrix
    end subroutine discretise_poisson

    module subroutine apply_dirichlet_bcs(mesh, M, b)
      type(ccs_mesh), intent(in) :: mesh    !< The mesh the problem is defined upon
      class(ccs_matrix), intent(inout) :: M !< The system matrix
      class(ccs_vector), intent(inout) :: b !< The system righthand side vector
    end subroutine apply_dirichlet_bcs
  end interface
  
end module poisson_discretisation

program poisson

  use poisson_discretisation
  use problem_setup
  
  ! ASiMoV-CCS uses
  use constants, only: ndim, add_mode, insert_mode
  use kinds, only: ccs_real, ccs_int
  use case_config, only: velocity_solver_method_name, velocity_solver_precon_name, &
                         pressure_solver_method_name, pressure_solver_precon_name
  use types, only: vector_spec, ccs_vector, matrix_spec, ccs_matrix, &
                   equation_system, linear_solver, ccs_mesh, cell_locator, face_locator, &
                   neighbour_locator, vector_values, matrix_values, matrix_values_spec, field, central_field
  use meshing, only: set_cell_location, set_face_location, set_neighbour_location, get_local_num_cells
  use vec, only: create_vector
  use mat, only: create_matrix, set_nnz, create_matrix_values, set_matrix_values_spec_nrows, &
                 set_matrix_values_spec_ncols
  use solver, only: create_solver, solve, set_equation_system, axpy, norm, &
       set_solver_method, set_solver_precon
  use utils, only: update, begin_update, end_update, finalise, initialise, &
                   set_size, zero, &
                   set_values, clear_entries, set_values, set_row, set_col, set_entry, set_mode
  use vec, only: create_vector_values
  use mesh_utils, only: build_square_mesh
  use meshing, only: get_face_area, get_centre, get_volume, get_global_index, &
                     count_neighbours, get_boundary_status
  use parallel_types, only: parallel_environment
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, &
                      read_command_line_arguments, &
                      timer, sync
  use timestepping, only: update_old_values, finalise_timestep, activate_timestepping, &
                          set_timestep, initialise_old_values, apply_timestep

  implicit none

  class(parallel_environment), allocatable, target :: par_env
  class(ccs_vector), allocatable, target :: b, diag
  class(field), allocatable, target :: u
  class(ccs_vector), allocatable :: u_exact
  class(ccs_matrix), allocatable, target :: M
  class(linear_solver), allocatable :: poisson_solver

  type(vector_spec) :: vec_properties
  type(matrix_spec) :: mat_properties
  type(equation_system) :: poisson_eq
  type(ccs_mesh) :: mesh

  integer(ccs_int) :: cps = 10 !< Default value for cells per side
  integer(ccs_int) :: num_steps = 20 !< Default number of timesteps
  integer(ccs_int) :: t
  real(ccs_real) :: dt = 1e-8

  real(ccs_real) :: err_norm

  double precision :: start_time
  double precision :: end_time

  call initialise_parallel_environment(par_env)
  call read_command_line_arguments(par_env, cps=cps)

  allocate (central_field :: u)
  
  ! set solver and preconditioner info
  velocity_solver_method_name = "gmres"
  velocity_solver_precon_name = "bjacobi"
  pressure_solver_method_name = "cg"
  pressure_solver_precon_name = "gamg"

  call sync(par_env)
  call timer(start_time)

  call initialise_poisson(par_env)

  ! Initialise with default values
  call initialise(vec_properties)
  call initialise(mat_properties)
  call initialise(poisson_eq)

  call set_size(par_env, mesh, mat_properties)
  call set_nnz(5, mat_properties)

  ! Create right-hand-side and solution vectors
  call set_size(par_env, mesh, vec_properties)
  call create_vector(vec_properties, u_exact)
  call create_vector(vec_properties, diag)
  call create_vector(vec_properties, u%values)

  ! Initialise with exact solution
  call set_exact_sol(u_exact)

  call update(u%values)
  call initialise_old_values(vec_properties, u)

  call activate_timestepping()
  call set_timestep(dt)

  do t = 1, num_steps

    call update_old_values(u)

    ! Create stiffness matrix
    call create_matrix(mat_properties, M)
    call create_vector(vec_properties, b)
    
    call zero(M)
    call zero(b)

    call discretise_poisson(mesh, M)
    ! Evaluate right-hand-side vector
    call eval_rhs(mesh, b)

    call apply_timestep(mesh, u, diag, M, b)

    call update(M) ! Start the parallel assembly for M
    call update(u%values) ! Start the parallel assembly for u
    call update(b) ! Start the parallel assembly for b

    ! Modify matrix and right-hand-side vector to apply Dirichlet boundary conditions
    call apply_dirichlet_bcs(mesh, M, b)
    call update(b) ! Start the parallel assembly for b
    call finalise(M)


    ! Create linear solver & set options
    call set_equation_system(par_env, b, u%values, M, poisson_eq)
    if (t .eq. 1) then
      call create_solver(poisson_eq, poisson_solver)
      call set_solver_method(pressure_solver_method_name, poisson_solver)
      call set_solver_precon(pressure_solver_precon_name, poisson_solver)
    end if
    call solve(poisson_solver)


    call finalise_timestep()

    call axpy(-1.0_ccs_real, u_exact, u%values)
    err_norm = norm(u%values, 2) * mesh%geo%h

    if (par_env%proc_id == par_env%root) then
      print *, " TIME = ", t, "Error norm = ", err_norm
    end if
  end do
 

  ! Clean up
  deallocate (u)
  deallocate (b)
  deallocate (u_exact)
  deallocate (M)
  deallocate (poisson_solver)

  call timer(end_time)

  if (par_env%proc_id == par_env%root) then
    print *, "Elapsed time = ", (end_time - start_time)
  end if

  call cleanup_parallel_environment(par_env)

contains

  subroutine set_exact_sol(u_exact)

    class(ccs_vector), intent(inout) :: u_exact

    type(vector_values) :: vec_values
    integer(ccs_int) :: i, local_num_cells

    type(cell_locator) :: loc_p
    integer(ccs_int) :: global_index_p
    integer(ccs_int) :: nrows_working_set

    nrows_working_set = 1_ccs_int
    call create_vector_values(nrows_working_set, vec_values)
    call set_mode(insert_mode, vec_values)

    call get_local_num_cells(mesh, local_num_cells)
    do i = 1, local_num_cells
      call clear_entries(vec_values)
      call set_cell_location(mesh, i, loc_p)
      call get_global_index(loc_p, global_index_p)

      call set_row(global_index_p, vec_values)
      call set_entry(eval_solution(loc_p), vec_values)
      call set_values(vec_values, u_exact)
    end do
    deallocate (vec_values%global_indices)
    deallocate (vec_values%values)

    call update(u_exact)
  end subroutine set_exact_sol

  subroutine initialise_poisson(par_env)

    class(parallel_environment) :: par_env

    mesh = build_square_mesh(par_env, cps, 1.0_ccs_real)

  end subroutine initialise_poisson

  !> Forcing function
  subroutine eval_rhs(mesh, b)

    type(ccs_mesh), intent(in) :: mesh
    class(ccs_vector), intent(inout) :: b

    integer(ccs_int) :: nloc
    integer(ccs_int) :: i
    real(ccs_real) :: r

    type(vector_values) :: val_dat

    type(cell_locator) :: loc_p
    real(ccs_real), dimension(ndim) :: cc
    real(ccs_real) :: V
    integer(ccs_int) :: global_index_p
    integer(ccs_int) :: nrows_working_set

    nrows_working_set = 1_ccs_int
    call create_vector_values(nrows_working_set, val_dat)
    call set_mode(add_mode, val_dat)

    call get_local_num_cells(mesh, nloc)
    associate (h => mesh%geo%h)
      ! this is currently setting 1 vector value at a time
      ! consider changing to doing all the updates in one go
      ! to do only 1 call to eval_cell_rhs and set_values
      do i = 1, nloc
        call clear_entries(val_dat)

        call set_cell_location(mesh, i, loc_p)
        call get_centre(loc_p, cc)
        call get_volume(loc_p, V)
        call get_global_index(loc_p, global_index_p)
        associate (x => cc(1), y => cc(2))
          call eval_cell_rhs(x, y, h**2, r)
          r = V * r
          call set_row(global_index_p, val_dat)
          call set_entry(r, val_dat)
          call set_values(val_dat, b)
        end associate
      end do
    end associate

    deallocate (val_dat%global_indices)
    deallocate (val_dat%values)

  end subroutine eval_rhs

end program poisson
