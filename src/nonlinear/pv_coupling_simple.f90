!> @brief Submodule file pv_coupling_simple.smod
!
!> @details Implementation of the SIMPLE algorithm for pressure-velocity coupling.

submodule (pv_coupling) pv_coupling_simple

  use kinds, only: accs_real, accs_int
  use types, only: vector_init_data, vector, matrix_init_data, matrix, &
                   linear_system, linear_solver, mesh, set_global_matrix_size, &
                  field, upwind_field, central_field
  use fv, only: compute_fluxes
  use vec, only: create_vector
  use mat, only: create_matrix, set_nnz
  use utils, only: update, initialise, finalise
  use solver, only: create_solver, solve, set_linear_system
  use parallel_types, only: parallel_environment
!  use petsctypes, only: matrix_petsc, vector_petsc

  implicit none

  contains

  !> @ brief Solve Navier-Stokes equations using the SIMPLE algorithm
  !
  !> @param[in]  par_env   - parallel environment
  !> @param[in]  cell_mesh - the mesh
  !> @param[in,out] u, v   - arrays containing velocity fields in x, y directions
  !> @param[in,out] p      - array containing pressure values
  !> @param[in,out] pp     - array containing pressure-correction values
  module subroutine solve_nonlinear(par_env, cell_mesh, it_start, it_end, u, v, p, pp)

    ! Arguments
    class(parallel_environment), intent(in) :: par_env
    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: it_start, it_end
    class(field), intent(inout) :: u, v, p, pp
    
    ! Local variables
    integer(accs_int) :: i
    class(vector), allocatable :: source
    class(matrix), allocatable :: M
    class(vector), allocatable :: invAu, invAv
    
    type(vector_init_data) :: vec_sizes
    type(matrix_init_data) :: mat_sizes
    type(linear_system)    :: lin_system

    ! Initialise linear system
    call initialise(mat_sizes)
    call initialise(vec_sizes)
    call initialise(lin_system)

    ! Create coefficient matrix
    call set_global_size(mat_sizes, cell_mesh%nglobal, cell_mesh%nglobal, par_env)
    call set_nnz(mat_sizes, 5)
    call create_matrix(mat_sizes, M)

    ! Create RHS and solution vectors
    call set_global_size(vec_sizes, cell_mesh%nglobal, par_env)
    call create_vector(vec_sizes, source)
    call create_vector(vec_sizes, sol)

    ! Create vectors for storing inverse of velocity central coefficients
    call create_vector(vec_sizes, invAu)
    call create_vector(vec_sizes, invAv)

    outerloop: do i = it_start, it_end
      ! Solve momentum equation with guessed pressure and velocity fields (eq. 4)
      call calculate_velocity(cell_mesh, bcs, cps, M, source, lin_system u, v, invAu, invAv)

      ! Calculate pressure correction from mass imbalance (sub. eq. 11 into eq. 8)
      call calculate_pressure_correction(cell_mesh, M, source, lin_system, u, v, pp)

      ! Update pressure field with pressure correction
      call update_pressure(cell_mesh, pp, p)

      ! Update velocity with velocity correction (eq. 6)
      call update_velocity(u, v, pp, invAu, invAv)

      ! Update face velocity (need data type for faces) (eq. 9)

      ! Todo:
      !call calculate_scalars()

      call check_convergence()
      if (converged) then
        exit outerloop
      endif

    end do outerloop

  end subroutine solve_nonlinear

  subroutine calculate_velocity(cell_mesh, bcs, cps, M, vec, lin_sys, u, v, invAu, invAv)

    ! Arguments
    type(mesh), intent(in)        :: cell_mesh
    type(bc_config), intent(in)   :: bcs
    integer(accs_int), intent(in) :: cps
    class(matrix), intent(inout)  :: M
    class(vector), intent(inout)  :: vec
    type(linear_system), intent(inout) :: lin_sys
    type(field), intent(inout)    :: u, v
    class(vector), intent(inout)    :: invAu, invAv

    ! Local variables
    class(linear_solver), allocatable :: lin_solver

    ! u-velocity
    ! ----------
    ! Calculate fluxes and populate coefficient matrix
    call compute_fluxes(u, u, v, cell_mesh, bcs, cps, M, vec)

    ! Calculate pressure source term and populate RHS vector
    call compute_pressure_source(cell_mesh, p, vec)

    ! Store reciprocal of central coefficient
    call get_matrix_diagonal(M, invAu)
    do i = 1, cell_mesh%nlocal
      invAu(i) = 1.0_accs_real / invAu(i)
    end do

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

    ! Store reciprocal of central coefficient
    call get_matrix_diagonal(M, invAv)
    do i = 1, cell_mesh_nlocal
      invAv(i) = 1.0_accs_real / invAv(i)
    end do

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

  subroutine calculate_pressure_source(cell_mesh, p, vec)

    ! Arguments
    class(mesh), intent(in) :: cell_mesh
    class(field), intent(in) :: p
    class(vector), intent(inout) :: vec

    ! Local variables
    type(vector_values) :: vec_values
    type(cell_locator) :: self_loc
    type(neighbour_locator) :: ngb_loc
    type(face_locator) :: face_loc
    integer(accs_int) :: self_idx, ngb_idx, local_idx
    integer(accs_int) :: j
    integer(accs_int) :: n_ngb
    real(accs_real) :: face_area
    real(accs_real) :: r
    logical :: is_boundary

    allocate(vec_values%idx(1))
    allocate(vec_values%val(1))

    vec_values%mode = insert_mode

    ! Loop over cells
    do local_idx = 1, cell_mesh%nlocal
      call set_cell_location(self_loc, cell_mesh, local_idx)
      call get_global_index(self_loc, self_idx)
      call count_neighbours(self_loc, n_ngb)

      r = 0.0_accs_real

      ! Loop over faces
      do j = 1, n_ngb
        call set_neighbour_location(ngb_loc, self_loc, j)
        call get_global_index(ngb_loc, ngb_idx)
        call get_boundary_status(ngb_loc, is_boundary)

        if (.not. is_boundary) then
          ! Interior face
          call set_face_location(face_loc, cell_mesh, local_idx, j)
          call get_face_area(face_loc, face_area)

          ! Calculate p at face using CDS
          r = r - 0.5_accs_real * (p(self_idx) + p(ngb_idx)) * face_area
        else
          ! Boundary face (zero gradient?)
          r = r - p(self_idx) * face_area
        endif
      end do
          
      call pack_entries(vec_values, 1, self_idx, r)
      call set_values(vec_values, vec)

    end do

    deallocate(vec_values%idx)
    deallocate(vec_values%val)

  end subroutine calculate_pressure_source


  subroutine calculate_pressure_correction(cell_mesh, M, vec, lin_sys, u, v, pp)

    ! Arguments
    class(mesh), intent(in) :: cell_mesh
    class(matrix), intent(inout)  :: M
    class(vector), intent(inout)  :: vec
    type(linear_system), intent(inout) :: lin_sys
    class(field), intent(in) :: u, v
    class(field), intent(out) :: pp

    ! Local variables
    type(matrix_values) :: mat_coeffs
    type(vector_values) :: vec_values
    type(neighbour_locator) :: ngb_loc
    type(face_locator) :: face_loc
    class(linear_solver), allocatable :: lin_solver
    integer(accs_int) :: self_idx, ngb_idx, local_idx
    integer(accs_int) :: j
    integer(accs_int) :: n_ngb
    real(accs_real) :: face_area
    real(accs_real), dimension(ndim) :: face_normal
    real(accs_real) :: r
    logical :: is_boundary

    allocate(vec_values%idx(1))
    allocate(vec_values%val(1))

    mat_coeffs%mode = insert_mode
    vec_values%mode = insert_mode
    
    ! Loop over cells
    do local_idx = 1, cell_mesh%nlocal
      call set_cell_location(self_loc, cell_mesh, local_idx)
      call get_global_index(self_loc, self_idx)
      call count_neighbours(self_loc, n_ngb)

      allocate(mat_coeffs%rglob(1))
      allocate(mat_coeffs%cglob(1 + n_ngb))
      allocate(mat_coeffs%val(1 + n_ngb))

      row = self_loc
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
          coeff_f = (1.0 / cell_mesh%h) * A  ! multiply by -d/(1+gam.d)

          coeff_p = coeff_p - coeff_f
          coeff_nb = coeff_f
          col = ngb_idx

          ! RHS vector
          r = r + calc_mass_flux(u, v, ngb_idx, self_idx, face_area, face_normal, bc_flag)

        else
          col = -1
          coeff_nb = 0.0_accs_real
        endif
        call pack_entries(mat_coeffs, 1, j+1, row, col, coeff_nb)
        call pack_entries(vec_values, 1, idx, r)

      end do

      ! Add the diagonal entry
      col = row
      call pack_entries(mat_coeffs, 1, 1, row, col, coeff_p)

      ! Set the values
      call set_values(mat_coeffs, M)
      call set_values(vec_values, vec)

      deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)

    end do

    ! Assembly of coefficient matrix and source vector
    call update(M)
    call update(vec)

    ! Create linear solver
    call set_linear_system(lin_system, vec, pp, M, par_env)
    call create_solver(lin_system, lin_solver)

    ! Solve the linear system
    call solve(lin_solver)

    ! Clean up
    deallocate(lin_solver)

  end subroutine calculate_pressure_correction


  subroutine update_pressure(cell_mesh, pp, p)

    ! Arguments
    class(mesh), intent(in) :: cell_mesh
    class(field), intent(in) :: pp
    class(field), intent(inout) :: p

    ! Local variables
    integer(accs_int) :: icell
    real(accs_real) :: alpha   !< Under-relaxation factor

    ! Loop over cells
    do icell = 1, cell_mesh%nlocal
      p%vec(icell) = p%vec(icell) + alpha*pp%vec(icell)
    end do
         
  end subroutine update_pressure

  subroutine update_velocity(cell_mesh, invAu, invAv, pp, u, v)

    ! Arguments
    class(mesh), intent(in) :: cell_mesh
    class(vector), intent(in) :: invAu, invAv
    class(field), intent(in)  :: pp
    class(field), intent(inout) :: u, v

    ! Local variables
    type(vector_values) :: vec_values
    type(cell_locator) :: self_loc
    type(neighbour_locator) :: ngb_loc
    type(face_locator) :: face_loc
    integer(accs_int) :: self_idx, ngb_idx, local_idx
    integer(accs_int) :: j
    integer(accs_int) :: n_ngb
    real(accs_real) :: face_area
    real(accs_real) :: up, vp
    logical :: is_boundary

    allocate(vec_values%idx(1))
    allocate(vec_values%val(1))

    vec_values%mode = add_mode
     
    ! Loop over cells
    do local_idx = 1, cell_mesh%nlocal
      call set_cell_location(self_loc, cell_mesh, local_idx)
      call get_global_index(self_loc, self_idx)
      call count_neighbours(cell_location, n_ngb)

      up = 0.0_accs_real
      vp = 0.0_accs_real

      ! Loop over faces to calculate pressure correction gradient
      do j = 1, n_ngb
        call set_neighbour_location(ngb_loc, self_loc, j)
        call get_global_index(ngb_loc, ngb_idx)
        call get_boundary_status(ngb_loc, is_boundary)

        if (.not. is_boundary) then
          ! Interior face
          call set_face_location(face_loc, cell_mesh, local_idx, j)
          call get_face_area(face_loc, face_area)
          
          up = up - 0.5_accs_real * (pp(self_idx) + pp(ngb_idx)) * face_area

        else
          ! Boundary face (zero gradient?)
          up = up - pp(self_idx) * face_area
        endif

      end do

      vp = up

      up = up * invAu(self_idx)
      vp = vp * invAv(self_idx)

      ! Update u and v vectors with velocity correction
      call pack_entries(vec_values, 1, self_idx, up)
      call set_values(vec_values, u)

      call pack_entries(vec_values, 1, self_idx, vp)
      call set_values(vec_values, v)

    end do

  end subroutine update_velocity


  subroutine calculate_scalars()


  end subroutine calculate_scalars


  subroutine check_convergence(converged)

    ! Arguments
    logical, intent(inout) :: converged

  end subroutine check_convergence

end submodule pv_coupling_simple
