!> @brief Submodule file pv_coupling_simple.smod
!
!> @details Implementation of the SIMPLE algorithm for pressure-velocity coupling.

submodule (pv_coupling) pv_coupling_simple

  use kinds, only: accs_real, accs_int
  use types, only: vector_init_data, vector, matrix_init_data, matrix, &
                   linear_system, linear_solver, mesh, &
                   field, upwind_field, central_field, bc_config, &
                   vector_values, cell_locator, face_locator, neighbour_locator, &
                   matrix_values
  use fv, only: compute_fluxes, calc_mass_flux, update_gradient_component
  use vec, only: create_vector, vec_reciprocal, get_vector_data, restore_vector_data
  use mat, only: create_matrix, set_nnz, get_matrix_diagonal
  use utils, only: update, initialise, finalise, set_global_size, &
                   set_values, pack_entries, mult
  use solver, only: create_solver, solve, set_linear_system, axpy
  use parallel_types, only: parallel_environment
  use constants, only: insert_mode, add_mode, ndim
  use meshing, only: get_face_area, get_global_index, count_neighbours, &
                     get_boundary_status, get_face_normal, &
                     set_neighbour_location, &
                     set_face_location, set_cell_location, &
                     get_volume

  implicit none

  contains

  !> @ brief Solve Navier-Stokes equations using the SIMPLE algorithm
  !
  !> @param[in]  par_env   - parallel environment
  !> @param[in]  cell_mesh - the mesh
  !> @param[in,out] u, v   - arrays containing velocity fields in x, y directions
  !> @param[in,out] p      - array containing pressure values
  !> @param[in,out] pp     - array containing pressure-correction values
  module subroutine solve_nonlinear(par_env, cell_mesh, cps, it_start, it_end, u, v, p, pp)

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env
    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: cps, it_start, it_end
    class(field), intent(inout) :: u, v, p, pp
    
    ! Local variables
    integer(accs_int) :: i
    class(vector), allocatable :: source
    class(matrix), allocatable :: M
    class(vector), allocatable :: invAu, invAv

    ! This could be done better, but Fortran arrays of allocatable things won't let me...
    class(vector), allocatable :: pgradx
    class(vector), allocatable :: pgrady
    
    type(vector_init_data) :: vec_sizes
    type(matrix_init_data) :: mat_sizes
    type(linear_system)    :: lin_system
    type(bc_config) :: bcs

    logical :: converged
    
    ! Initialise linear system
    call initialise(mat_sizes)
    call initialise(vec_sizes)
    call initialise(lin_system)

    ! Create coefficient matrix
    call set_global_size(mat_sizes, cell_mesh, par_env)
    call set_nnz(mat_sizes, 5)
    call create_matrix(mat_sizes, M)

    ! Create RHS vector
    call set_global_size(vec_sizes, cell_mesh, par_env)
    call create_vector(vec_sizes, source)

    ! Create vectors for storing inverse of velocity central coefficients
    call create_vector(vec_sizes, invAu)
    call create_vector(vec_sizes, invAv)

    ! Create vectors for pressure gradients
    call create_vector(vec_sizes, pgradx)
    call create_vector(vec_sizes, pgrady)

    outerloop: do i = it_start, it_end

      ! Get pressure gradient
      call update_gradient_component(cell_mesh, 1, p, pgradx)
      call update_gradient_component(cell_mesh, 2, p, pgrady)
      
      ! Solve momentum equation with guessed pressure and velocity fields (eq. 4)
      call calculate_velocity(par_env,cell_mesh, bcs, cps, M, source, lin_system, u, v, pgradx, pgrady, invAu, invAv)

      ! Calculate pressure correction from mass imbalance (sub. eq. 11 into eq. 8)
      call calculate_pressure_correction(par_env, cell_mesh, M, source, lin_system, u, v, pp)

      ! Update pressure field with pressure correction
      call update_pressure(pp, p)
      
      ! Update velocity with velocity correction (eq. 6)
      call update_velocity(cell_mesh, invAu, invAv, pp, pgradx, pgrady, u, v)

      ! Update face velocity (need data type for faces) (eq. 9)

      ! Todo:
      !call calculate_scalars()

      call check_convergence(converged)
      if (converged) then
        exit outerloop
      endif

    end do outerloop

  end subroutine solve_nonlinear

  subroutine calculate_velocity(par_env, cell_mesh, bcs, cps, M, vec, lin_sys, u, v, pgradx, pgrady, invAu, invAv)

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env
    type(mesh), intent(in)        :: cell_mesh
    type(bc_config), intent(in)   :: bcs
    integer(accs_int), intent(in) :: cps
    class(matrix), allocatable, intent(inout)  :: M
    class(vector), allocatable, intent(inout)  :: vec
    type(linear_system), intent(inout) :: lin_sys
    type(field), intent(inout)    :: u, v
    class(vector), intent(in) :: pgradx, pgrady
    class(vector), intent(inout)    :: invAu, invAv

    ! Local variables
    class(linear_solver), allocatable :: lin_solver
    real(accs_real) :: alpha !> Underrelaxation factor


    ! Set underrelaxation factor
    ! TODO: read from input
    alpha = 0.7_accs_real
    
    ! u-velocity
    ! ----------
    ! Calculate fluxes and populate coefficient matrix
    call compute_fluxes(u, u, v, cell_mesh, bcs, cps, M, vec)

    ! Calculate pressure source term and populate RHS vector
    call calculate_pressure_source(cell_mesh, pgradx, vec)

    ! Underrelax the equations
    call underrelax(cell_mesh, alpha, u, invAu, M, vec)
    
    ! Store reciprocal of central coefficient
    call get_matrix_diagonal(M, invAu)
    call vec_reciprocal(invAu)
    
    ! Assembly of coefficient matrix and source vector
    call update(M)
    call update(vec)

    ! Create linear solver
    call set_linear_system(lin_sys, vec, u%vec, M, par_env)
    call create_solver(lin_sys, lin_solver)

    ! Solve the linear system
    call solve(lin_solver)

    ! v-velocity
    ! ----------
    ! Calculate fluxes and populate coefficient matrix
    call compute_fluxes(v, u, v, cell_mesh, bcs, cps, M, vec)

    ! Calculate pressure source term and populate RHS vector
    call calculate_pressure_source(cell_mesh, pgrady, vec)

    ! Store reciprocal of central coefficient
    call get_matrix_diagonal(M, invAv)
    call vec_reciprocal(invAv)

    ! Underrelax the equations
    call underrelax(cell_mesh, alpha, v, invAv, M, vec)

    ! Assembly of coefficient matrix and source vector
    call update(M)
    call update(vec)

    ! Create linear solver
    call set_linear_system(lin_sys, vec, v%vec, M, par_env)
    call create_solver(lin_sys, lin_solver)

    ! Solve the linear system
    call solve(lin_solver)

    ! Clean up
    deallocate(lin_solver)

  end subroutine calculate_velocity

  subroutine calculate_pressure_source(cell_mesh, pgrad, vec)

    ! Arguments
    class(mesh), intent(in) :: cell_mesh
    class(vector), intent(in) :: pgrad
    class(vector), intent(inout) :: vec

    ! Local variables
    type(vector_values) :: vec_values
    type(cell_locator) :: self_loc
    integer(accs_int) :: self_idx, local_idx
    real(accs_real) :: r
    real(accs_real), dimension(:), pointer :: pgrad_data

    real(accs_real) :: V
    
    allocate(vec_values%idx(1))
    allocate(vec_values%val(1))

    vec_values%mode = add_mode

    ! Temporary storage for p values
    call get_vector_data(pgrad, pgrad_data)
    
    ! Loop over cells
    do local_idx = 1, cell_mesh%nlocal
      call set_cell_location(self_loc, cell_mesh, local_idx)
      call get_global_index(self_loc, self_idx)

      call get_volume(self_loc, V)
      
      r = -pgrad_data(local_idx) * V
      call pack_entries(vec_values, 1, self_idx, r)
      call set_values(vec_values, vec)

    end do

    deallocate(vec_values%idx)
    deallocate(vec_values%val)

    call restore_vector_data(pgrad, pgrad_data)
    
  end subroutine calculate_pressure_source


  subroutine calculate_pressure_correction(par_env, cell_mesh, M, vec, lin_sys, u, v, pp)

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env
    class(mesh), intent(in) :: cell_mesh
    class(matrix), allocatable, intent(inout)  :: M
    class(vector), allocatable, intent(inout)  :: vec
    type(linear_system), intent(inout) :: lin_sys
    class(field), intent(in) :: u, v
    class(field), intent(out) :: pp

    ! Local variables
    type(matrix_values) :: mat_coeffs
    type(vector_values) :: vec_values
    type(cell_locator) :: self_loc
    type(neighbour_locator) :: ngb_loc
    type(face_locator) :: face_loc
    class(linear_solver), allocatable :: lin_solver
    integer(accs_int) :: self_idx, ngb_idx, local_idx
    integer(accs_int) :: j
    integer(accs_int) :: n_ngb
    integer(accs_int) :: row, col
    integer(accs_int) :: bc_flag
    real(accs_real) :: face_area
    real(accs_real), dimension(ndim) :: face_normal
    real(accs_real) :: r
    real(accs_real) :: coeff_f, coeff_p, coeff_nb
    real(accs_real), dimension(:), pointer :: u_data, v_data
    logical :: is_boundary

    allocate(vec_values%idx(1))
    allocate(vec_values%val(1))

    mat_coeffs%mode = insert_mode
    vec_values%mode = insert_mode

    ! Temporary storage
    call get_vector_data(u%vec, u_data)
    call get_vector_data(v%vec, v_data)
    
    ! Loop over cells
    do local_idx = 1, cell_mesh%nlocal
      call set_cell_location(self_loc, cell_mesh, local_idx)
      call get_global_index(self_loc, self_idx)
      call count_neighbours(self_loc, n_ngb)

      allocate(mat_coeffs%rglob(1))
      allocate(mat_coeffs%cglob(1 + n_ngb))
      allocate(mat_coeffs%val(1 + n_ngb))

      row = self_idx
      coeff_p = 0.0_accs_real
      r = 0.0_accs_real

      ! Loop over faces
      do j = 1, n_ngb
        call set_neighbour_location(ngb_loc, self_loc, j)
        call get_boundary_status(ngb_loc, is_boundary)
        call get_global_index(ngb_loc, ngb_idx)

        if (.not. is_boundary) then
          ! Interior face
          call set_face_location(face_loc, cell_mesh, local_idx, j)
          call get_face_area(face_loc, face_area)
          call get_face_normal(face_loc, face_normal)
          coeff_f = (1.0 / cell_mesh%h) * face_area  ! multiply by -d/(1+gam.d)

          coeff_p = coeff_p - coeff_f
          coeff_nb = coeff_f
          col = ngb_idx

          ! RHS vector
          bc_flag = 0
          r = r + calc_mass_flux(u_data, v_data, ngb_idx, self_idx, face_normal, bc_flag)

        else
          col = -1
          coeff_nb = 0.0_accs_real
        endif
        call pack_entries(mat_coeffs, 1, j+1, row, col, coeff_nb)
        call pack_entries(vec_values, 1, self_idx, r)

      end do

      ! Add the diagonal entry
      col = row
      call pack_entries(mat_coeffs, 1, 1, row, col, coeff_p)

      ! Set the values
      call set_values(mat_coeffs, M)
      call set_values(vec_values, vec)

    end do

    ! Assembly of coefficient matrix and source vector
    call update(M)
    call update(vec)

    ! Create linear solver
    call set_linear_system(lin_sys, vec, pp%vec, M, par_env)
    call create_solver(lin_sys, lin_solver)

    ! Solve the linear system
    call solve(lin_solver)

    ! Clean up
    deallocate(lin_solver)
    deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)

    call restore_vector_data(u%vec, u_data)
    call restore_vector_data(v%vec, v_data)
    
  end subroutine calculate_pressure_correction


  subroutine update_pressure(pp, p)

    ! Arguments
    class(field), intent(in) :: pp
    class(field), intent(inout) :: p

    ! Local variables
    real(accs_real) :: alpha   !< Under-relaxation factor

    ! Set under-relaxation factor (todo: read this from input file)
    alpha = 0.3_accs_real

    call axpy(alpha, pp%vec, p%vec)
    
  end subroutine update_pressure

  subroutine update_velocity(cell_mesh, invAu, invAv, pp, ppgradx, ppgrady, u, v)

    ! Arguments
    class(mesh), intent(in) :: cell_mesh
    class(vector), intent(in) :: invAu, invAv
    class(field), intent(in) :: pp
    class(vector), intent(inout)  :: ppgradx, ppgrady
    class(field), intent(inout) :: u, v

    ! First update gradients
    call update_gradient_component(cell_mesh, 1, pp, ppgradx)
    call update_gradient_component(cell_mesh, 2, pp, ppgrady)

    ! Multiply gradients by inverse diagonal coefficients
    call mult(invAu, ppgradx)
    call mult(invAv, ppgrady)

    ! Compute correction source on velocity
    call calculate_pressure_source(cell_mesh, ppgradx, u%vec)
    call calculate_pressure_source(cell_mesh, ppgrady, v%vec)
    
  end subroutine update_velocity


  subroutine calculate_scalars()


  end subroutine calculate_scalars


  subroutine check_convergence(converged)

    ! Arguments
    logical, intent(inout) :: converged

    converged = .true. ! XXX: temporary
    
  end subroutine check_convergence

  subroutine underrelax(cell_mesh, alpha, phi, diag, M, b)

    use mat, only : set_matrix_diagonal
    
    type(mesh), intent(in) :: cell_mesh
    real(accs_real), intent(in) :: alpha
    class(field), intent(in) :: phi
    class(vector), intent(inout) :: diag
    class(matrix), intent(inout) :: M
    class(vector), intent(inout) :: b

    real(accs_real), dimension(:), pointer :: diag_data
    real(accs_real), dimension(:), pointer :: phi_data
    real(accs_real), dimension(:), pointer :: b_data

    integer(accs_int) :: i
    
    call get_matrix_diagonal(M, diag)

    call get_vector_data(phi%vec, phi_data)
    call get_vector_data(diag, diag_data)
    call get_vector_data(b, b_data)
    
    do i = 1, cell_mesh%nlocal
      diag_data(i) = diag_data(i) / alpha

      b_data(i) = b_data(i) + (1.0_accs_real - alpha) * diag_data(i) * phi_data(i)
    end do

    call restore_vector_data(phi%vec, phi_data)
    call restore_vector_data(diag, diag_data)
    call restore_vector_data(b, b_data)
    
    call set_matrix_diagonal(diag, M)
    
  end subroutine underrelax
  
end submodule pv_coupling_simple
