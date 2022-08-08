!v  Program file for Poisson case
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

program poisson

  ! ASiMoV-CCS uses
  use constants, only: ndim, add_mode, insert_mode
  use kinds, only: ccs_real, ccs_int
  use types, only: vector_spec, ccs_vector, matrix_spec, ccs_matrix, &
                   equation_system, linear_solver, ccs_mesh, cell_locator, face_locator, &
                   neighbour_locator, vector_values, matrix_values, matrix_values_spec
  use meshing, only: set_cell_location, set_face_location, set_neighbour_location
  use vec, only: create_vector
  use mat, only: create_matrix, set_nnz, create_matrix_values, set_matrix_values_spec_nrows, &
                 set_matrix_values_spec_ncols
  use solver, only: create_solver, solve, set_equation_system, axpy, norm
  use utils, only: update, begin_update, end_update, finalise, initialise, &
                   set_size, &
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

  implicit none

  class(parallel_environment), allocatable, target :: par_env
  class(ccs_vector), allocatable, target :: u, b
  class(ccs_vector), allocatable :: u_exact
  class(ccs_matrix), allocatable, target :: M
  class(linear_solver), allocatable :: poisson_solver

  type(vector_spec) :: vec_properties
  type(matrix_spec) :: mat_properties
  type(equation_system) :: poisson_eq
  type(ccs_mesh) :: mesh

  integer(ccs_int) :: cps = 10 !< Default value for cells per side

  real(ccs_real) :: err_norm

  double precision :: start_time
  double precision :: end_time

  call initialise_parallel_environment(par_env)
  call read_command_line_arguments(par_env, cps=cps)

  call sync(par_env)
  call timer(start_time)

  call initialise_poisson(par_env)

  ! Initialise with default values
  call initialise(vec_properties)
  call initialise(mat_properties)
  call initialise(poisson_eq)

  ! Create stiffness matrix
  call set_size(par_env, mesh, mat_properties)
  call set_nnz(5, mat_properties)
  call create_matrix(mat_properties, M)

  call discretise_poisson(M)

  call begin_update(M) ! Start the parallel assembly for M

  ! Create right-hand-side and solution vectors
  call set_size(par_env, mesh, vec_properties)
  call create_vector(vec_properties, b)
  call create_vector(vec_properties, u_exact)
  call create_vector(vec_properties, u)

  call begin_update(u) ! Start the parallel assembly for u

  ! Evaluate right-hand-side vector
  call eval_rhs(b)

  call begin_update(b) ! Start the parallel assembly for b
  call end_update(M) ! Complete the parallel assembly for M
  call end_update(b) ! Complete the parallel assembly for b

  ! Modify matrix and right-hand-side vector to apply Dirichlet boundary conditions
  call apply_dirichlet_bcs(M, b)
  call begin_update(b) ! Start the parallel assembly for b
  call finalise(M)

  call end_update(u) ! Complete the parallel assembly for u
  call end_update(b) ! Complete the parallel assembly for b

  ! Create linear solver & set options
  call set_equation_system(par_env, b, u, M, poisson_eq)
  call create_solver(poisson_eq, poisson_solver)
  call solve(poisson_solver)

  ! Check solution
  call set_exact_sol(u_exact)
  call axpy(-1.0_ccs_real, u_exact, u)

  err_norm = norm(u, 2) * mesh%h
  if (par_env%proc_id == par_env%root) then
    print *, "Norm of error = ", err_norm
  end if

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

  subroutine eval_rhs(b)

    class(ccs_vector), intent(inout) :: b

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

    associate (nloc => mesh%nlocal, &
               h => mesh%h)
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

  !> Apply forcing function
  pure subroutine eval_cell_rhs(x, y, H, r)

    real(ccs_real), intent(in) :: x, y, H
    real(ccs_real), intent(out) :: r

    r = 0.0_ccs_real &
        + 0.0_ccs_real * (x + y + H) ! Silence unused dummy argument error

  end subroutine eval_cell_rhs

  subroutine discretise_poisson(M)

    class(ccs_matrix), intent(inout) :: M

    type(matrix_values_spec) :: mat_val_spec
    type(matrix_values) :: mat_coeffs
    integer(ccs_int) :: i, j

    integer(ccs_int) :: row, col
    real(ccs_real) :: coeff_f, coeff_p, coeff_nb

    type(face_locator) :: loc_f
    real(ccs_real) :: A

    integer(ccs_int) :: global_index_p
    type(cell_locator) :: loc_p
    integer(ccs_int) :: nnb

    type(neighbour_locator) :: loc_nb
    logical :: is_boundary
    integer(ccs_int) :: global_index_nb

    ! Loop over cells
    do i = 1, mesh%nlocal
      !^ @todo Doing this in a loop is awful code - malloc maximum coefficients per row once,
      !        filling from front, and pass the number of coefficients to be set, requires
      !        modifying the matrix_values type and the implementation of set_values applied to
      !        matrices.
      call set_cell_location(mesh, i, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)

      call set_matrix_values_spec_nrows(1_ccs_int, mat_val_spec)
      call set_matrix_values_spec_ncols((1_ccs_int + nnb), mat_val_spec)
      call create_matrix_values(mat_val_spec, mat_coeffs)
      call set_mode(insert_mode, mat_coeffs)

      row = global_index_p
      coeff_p = 0.0_ccs_real

      ! Loop over faces
      do j = 1, nnb
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)

        if (.not. is_boundary) then
          ! Interior face

          call set_face_location(mesh, i, j, loc_f)
          call get_face_area(loc_f, A)
          coeff_f = (1.0 / mesh%h) * A

          call get_global_index(loc_nb, global_index_nb)

          coeff_p = coeff_p - coeff_f
          coeff_nb = coeff_f
          col = global_index_nb
        else
          col = -1
          coeff_nb = 0.0_ccs_real
        end if

        call set_row(row, mat_coeffs)
        call set_col(col, mat_coeffs)
        call set_entry(coeff_nb, mat_coeffs)

      end do

      ! Add the diagonal entry
      col = row
      call set_row(row, mat_coeffs)
      call set_col(col, mat_coeffs)
      call set_entry(coeff_p, mat_coeffs)

      ! Set the values
      call set_values(mat_coeffs, M)

    end do

  end subroutine discretise_poisson

  subroutine apply_dirichlet_bcs(M, b)

    implicit none

    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b

    integer(ccs_int) :: i, j
    real(ccs_real) :: boundary_coeff, boundary_val

    integer(ccs_int) :: idx, row, col
    real(ccs_real) :: r, coeff

    type(vector_values) :: vec_values
    type(matrix_values_spec) :: mat_val_spec
    type(matrix_values) :: mat_coeffs

    type(face_locator) :: loc_f
    real(ccs_real) :: A

    type(cell_locator) :: loc_p
    integer(ccs_int) :: global_index_p
    type(neighbour_locator) :: loc_nb

    integer(ccs_int) :: nnb
    logical :: is_boundary

    integer(ccs_int) :: nrows_working_set

    call set_matrix_values_spec_nrows(1_ccs_int, mat_val_spec)
    call set_matrix_values_spec_ncols(1_ccs_int, mat_val_spec)
    call create_matrix_values(mat_val_spec, mat_coeffs)
    call set_mode(add_mode, mat_coeffs)

    nrows_working_set = 1_ccs_int
    call create_vector_values(nrows_working_set, vec_values)
    call set_mode(add_mode, vec_values)

    do i = 1, mesh%nlocal
      if (minval(mesh%neighbour_indices(:, i)) < 0) then
        call clear_entries(mat_coeffs)
        call clear_entries(vec_values)
        call set_cell_location(mesh, i, loc_p)
        call get_global_index(loc_p, global_index_p)
        coeff = 0.0_ccs_real
        r = 0.0_ccs_real

        row = global_index_p
        col = global_index_p
        idx = global_index_p

        call count_neighbours(loc_p, nnb)
        do j = 1, nnb

          call set_neighbour_location(loc_p, j, loc_nb)
          call get_boundary_status(loc_nb, is_boundary)

          if (is_boundary) then
            call set_face_location(mesh, i, j, loc_f)
            call get_face_area(loc_f, A)
            boundary_coeff = (2.0 / mesh%h) * A
            boundary_val = rhs_val(i, j)

            ! Coefficient
            coeff = coeff - boundary_coeff

            ! RHS vector
            r = r - boundary_val * boundary_coeff
          end if

        end do

        call set_row(row, mat_coeffs)
        call set_col(col, mat_coeffs)
        call set_entry(coeff, mat_coeffs)

        call set_row(row, vec_values)
        call set_entry(r, vec_values)

        call set_values(mat_coeffs, M)
        call set_values(vec_values, b)

      end if
    end do

    deallocate (vec_values%global_indices)
    deallocate (vec_values%values)

  end subroutine apply_dirichlet_bcs

  subroutine set_exact_sol(u_exact)

    class(ccs_vector), intent(inout) :: u_exact

    type(vector_values) :: vec_values
    integer(ccs_int) :: i

    type(cell_locator) :: loc_p
    integer(ccs_int) :: global_index_p
    integer(ccs_int) :: nrows_working_set

    nrows_working_set = 1_ccs_int
    call create_vector_values(nrows_working_set, vec_values)
    call set_mode(insert_mode, vec_values)

    do i = 1, mesh%nlocal
      call clear_entries(vec_values)
      call set_cell_location(mesh, i, loc_p)
      call get_global_index(loc_p, global_index_p)

      call set_row(global_index_p, vec_values)
      call set_entry(rhs_val(i), vec_values)
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

  function rhs_val(i, f) result(r)

    integer(ccs_int), intent(in) :: i !< Cell index
    integer(ccs_int), intent(in), optional :: f !< Face index (local wrt cell)

    type(cell_locator) :: loc_p
    type(face_locator) :: loc_f

    real(ccs_real), dimension(ndim) :: x
    real(ccs_real) :: r

    if (present(f)) then
      ! Face-centred value
      call set_face_location(mesh, i, f, loc_f)
      call get_centre(loc_f, x)
      associate (y => x(2))
        r = y
      end associate
    else
      ! Cell-centred value
      call set_cell_location(mesh, i, loc_p)
      call get_centre(loc_p, x)
      associate (y => x(2))
        r = y
      end associate
    end if

  end function rhs_val

end program poisson
