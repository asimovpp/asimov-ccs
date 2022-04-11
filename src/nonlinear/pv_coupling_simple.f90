!> @brief Submodule file pv_coupling_simple.smod
!
!> @details Implementation of the SIMPLE algorithm for pressure-velocity coupling.

submodule (pv_coupling) pv_coupling_simple

  use kinds, only: accs_real, accs_int
  use types, only: vector_init_data, vector, matrix_init_data, matrix, linear_system, &
                   linear_solver, mesh, field, bc_config, vector_values, cell_locator, &
                   face_locator, neighbour_locator, matrix_values
  use fv, only: compute_fluxes, calc_mass_flux, update_gradient
  use vec, only: create_vector, vec_reciprocal, get_vector_data, restore_vector_data
  use mat, only: create_matrix, set_nnz, get_matrix_diagonal
  use utils, only: update, initialise, finalise, set_global_size, set_values, pack_entries, &
                   mult, scale, zero
  use solver, only: create_solver, solve, set_linear_system, axpy, norm
  use parallel_types, only: parallel_environment
  use constants, only: insert_mode, add_mode, ndim
  use meshing, only: get_face_area, get_global_index, get_local_index, count_neighbours, &
                     get_boundary_status, get_face_normal, set_neighbour_location, set_face_location, &
                     set_cell_location, get_volume
  
  implicit none

contains

  !> @ brief Solve Navier-Stokes equations using the SIMPLE algorithm
  !
  !> @param[in]  par_env   - parallel environment
  !> @param[in]  cell_mesh - the mesh
  !> @param[in,out] u, v   - fields containing velocity fields in x, y directions
  !> @param[in,out] p      - field containing pressure values
  !> @param[in,out] pp     - field containing pressure-correction values
  !> @param[in,out] mf     - field containing the face-centred velocity flux 
  module subroutine solve_nonlinear(par_env, cell_mesh, cps, it_start, it_end, u, v, p, pp, mf)

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env
    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: cps, it_start, it_end
    class(field), intent(inout) :: u, v, p, pp, mf
    
    ! Local variables
    integer(accs_int) :: i
    class(vector), allocatable :: source
    class(matrix), allocatable :: M
    class(vector), allocatable :: invAu, invAv
    
    type(vector_init_data) :: vec_sizes
    type(matrix_init_data) :: mat_sizes
    type(linear_system)    :: lin_system
    type(bc_config) :: bcs

    logical :: converged
    
    ! Initialise linear system
    print *, "NONLINEAR: init"
    call initialise(mat_sizes)
    call initialise(vec_sizes)
    call initialise(lin_system)

    ! Create coefficient matrix
    print *, "NONLINEAR: setup matrix"
    call set_global_size(mat_sizes, cell_mesh, par_env)
    call set_nnz(mat_sizes, 5)
    call create_matrix(mat_sizes, M)

    ! Create RHS vector
    print *, "NONLINEAR: setup RHS"
    call set_global_size(vec_sizes, cell_mesh, par_env)
    call create_vector(vec_sizes, source)

    ! Create vectors for storing inverse of velocity central coefficients
    print *, "NONLINEAR: setup ind coeff"
    call create_vector(vec_sizes, invAu)
    call create_vector(vec_sizes, invAv)
    
    ! Get pressure gradient
    print *, "NONLINEAR: compute grad p"
    call update_gradient(cell_mesh, p)

    outerloop: do i = it_start, it_end

      print *, "NONLINEAR: iteration ", i

      ! Solve momentum equation with guessed pressure and velocity fields (eq. 4)
      print *, "NONLINEAR: guess velocity"
      call calculate_velocity(par_env, cell_mesh, bcs, cps, M, source, lin_system, u, v, mf, p, invAu, invAv)

      ! Calculate pressure correction from mass imbalance (sub. eq. 11 into eq. 8)
      print *, "NONLINEAR: mass imbalance"
      call compute_mass_imbalance(par_env, cell_mesh, u, v, p, invAu, invAv, mf, source)
      print *, "NONLINEAR: compute p'"
      call calculate_pressure_correction(par_env, cell_mesh, M, source, lin_system, invAu, invAv, pp)
      
      ! Update velocity with velocity correction (eq. 6)
      print *, "NONLINEAR: correct face velocity"
      call update_face_velocity(cell_mesh, pp, invAu, invAv, mf)
      print *, "NONLINEAR: correct velocity"
      call update_velocity(cell_mesh, invAu, invAv, pp, u, v)

      ! Update pressure field with pressure correction
      print *, "NONLINEAR: correct pressure"
      call update_pressure(pp, p)
      print *, "NONLINEAR: compute gradp"
      call update_gradient(cell_mesh, p)

      ! Todo:
      !call calculate_scalars()

      call check_convergence(converged)
      if (converged) then
        exit outerloop
      endif

    end do outerloop

  end subroutine solve_nonlinear

  !> @brief Computes the guessed velocity fields based on a frozen pressure field
  !
  !> @param[in]    par_env      - the parallel environment
  !> @param[in]    cell_mesh    - the mesh
  !> @param[in]    bcs          - boundary conditions
  !> @param[in]    cps          - cells per side of a square mesh (remove)
  !> @param[inout] M            - matrix object
  !> @param[inout] vec          - vector object
  !> @param[inout] lin_sys      - linear system object
  !> @param[inout] u, v         - the x and y velocity fields
  !> @param[in]    mf           - the face velocity flux
  !> @param[in]    p            - the pressure field
  !> @param[inout] invAu, invAv - vectors containing the inverse momentum coefficients
  !
  !> @description Given an initial guess of a pressure field form the momentum equations (as scalar
  !!              equations) and solve to obtain an intermediate velocity field u* that will not
  !!              satisfy continuity.
  subroutine calculate_velocity(par_env, cell_mesh, bcs, cps, M, vec, lin_sys, u, v, mf, p, invAu, invAv)

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env
    type(mesh), intent(in)         :: cell_mesh
    type(bc_config), intent(inout) :: bcs
    integer(accs_int), intent(in)  :: cps
    class(matrix), allocatable, intent(inout)  :: M
    class(vector), allocatable, intent(inout)  :: vec
    type(linear_system), intent(inout) :: lin_sys
    class(field), intent(inout)    :: u, v
    class(field), intent(in) :: mf
    class(field), intent(in) :: p
    class(vector), intent(inout)    :: invAu, invAv

    
    ! u-velocity
    ! ----------

    ! TODO: Do boundaries properly
    bcs%bc_type(:) = 0 !> Fixed zero BC
    bcs%bc_type(4) = 1 !> Fixed one BC at lid
    call calculate_velocity_component(par_env, cell_mesh, bcs, cps, M, vec, lin_sys, u, 1, mf, p, invAu)
    
    ! v-velocity
    ! ----------
    
    ! TODO: Do boundaries properly
    bcs%bc_type(:) = 0 !> Fixed zero BC
    call calculate_velocity_component(par_env, cell_mesh, bcs, cps, M, vec, lin_sys, v, 2, mf, p, invAv)

  end subroutine calculate_velocity

  subroutine calculate_velocity_component(par_env, cell_mesh, bcs, cps, M, vec, lin_sys, u, component, mf, p, invAu)

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env
    type(mesh), intent(in)         :: cell_mesh
    type(bc_config), intent(inout) :: bcs
    integer(accs_int), intent(in)  :: cps
    class(matrix), allocatable, intent(inout)  :: M
    class(vector), allocatable, intent(inout)  :: vec
    type(linear_system), intent(inout) :: lin_sys
    class(field), intent(inout)    :: u
    class(field), intent(in) :: mf
    class(field), intent(in) :: p
    class(vector), intent(inout)    :: invAu
    integer(accs_int), intent(in) :: component
    
    ! Local variables
    class(linear_solver), allocatable :: lin_solver
    real(accs_real) :: alpha !> Underrelaxation factor

    ! Set underrelaxation factor
    ! TODO: read from input
    alpha = 0.9_accs_real

    ! First zero matrix/RHS
    call zero(vec)
    call zero(M)
    
    ! Calculate fluxes and populate coefficient matrix
    print *, "GV: compute u flux"
    call compute_fluxes(u, mf, cell_mesh, bcs, cps, M, vec)

    ! Calculate pressure source term and populate RHS vector
    print *, "GV: compute u gradp"
    if (component == 1) then
      call calculate_momentum_pressure_source(cell_mesh, p%gradx, vec)
    else if (component == 2) then
      call calculate_momentum_pressure_source(cell_mesh, p%grady, vec)
    else
      print *, "Unsupported vector component: ", component
      stop 1
    end if
    
    ! Underrelax the equations
    print *, "GV: underrelax u"
    call underrelax(cell_mesh, alpha, u, invAu, M, vec)
    
    ! Store reciprocal of central coefficient
    print *, "GV: get u diag"
    call get_matrix_diagonal(M, invAu)
    call vec_reciprocal(invAu)
    
    ! Assembly of coefficient matrix and source vector
    print *, "GV: build u lin sys"
    call update(M)
    call update(vec)
    call update(invAu)
    call finalise(M)

    ! Create linear solver
    call set_linear_system(lin_sys, vec, u%vec, M, par_env)
    call create_solver(lin_sys, lin_solver)

    ! Solve the linear system
    print *, "GV: solve u"
    call solve(lin_solver)

    ! Clean up
    deallocate(lin_solver)
    
  end subroutine calculate_velocity_component
  
  !> @brief Adds the momentum source due to pressure gradient
  !
  !> @param[in]    cell_mesh - the mesh
  !> @param[in]    pgrad     - the pressure gradient
  !> @param[inout] vec       - the momentum equation RHS vector
  subroutine calculate_momentum_pressure_source(cell_mesh, pgrad, vec)

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
    
  end subroutine calculate_momentum_pressure_source

  !> @brief Solves the pressure correction equation
  !
  !> @param[in]    par_env      - the parallel environment
  !> @param[in]    cell_mesh    - the mesh
  !> @param[inout] M            - matrix object
  !> @param[inout] vec          - the RHS vector
  !> @param[inout] lin_sys      - linear system object
  !> @param[in]    invAu, invAv - inverse diagonal momentum coefficients
  !> @param[inout] pp           - the pressure correction field
  !
  !> @description Solves the pressure correction equation formed by the mass-imbalance.
  subroutine calculate_pressure_correction(par_env, cell_mesh, M, vec, lin_sys, invAu, invAv, pp)

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env
    class(mesh), intent(in) :: cell_mesh
    class(matrix), allocatable, intent(inout)  :: M
    class(vector), allocatable, intent(inout)  :: vec
    type(linear_system), intent(inout) :: lin_sys
    class(vector), intent(in) :: invAu, invAv
    class(field), intent(inout) :: pp

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
    real(accs_real) :: face_area
    real(accs_real), dimension(ndim) :: face_normal
    real(accs_real) :: r
    real(accs_real) :: coeff_f, coeff_p, coeff_nb
    logical :: is_boundary

    real(accs_real), dimension(:), pointer :: invAu_data
    real(accs_real), dimension(:), pointer :: invAv_data
    
    real(accs_real) :: Vp
    real(accs_real) :: Vnb
    real(accs_real) :: Vf
    real(accs_real) :: invAp
    real(accs_real) :: invAnb
    real(accs_real) :: invAf

    integer(accs_int) :: idxnb

    integer(accs_int) :: cps   ! Cells per side
    integer(accs_int) :: rcrit ! Global index of approximate central cell
    
    ! First zero matrix
    call zero(M)

    ! The computed mass imbalance is +ve, to have a +ve diagonal coefficient we need to negate this.
    print *, "P': negate RHS"
    call scale(-1.0_accs_real, vec)
    call update(vec)
    
    allocate(vec_values%idx(1))
    allocate(vec_values%val(1))

    mat_coeffs%mode = insert_mode
    vec_values%mode = add_mode    ! We already have a mass-imbalance vector, BCs get ADDED

    call update(M)
    
    print *, "P': get invA"
    call get_vector_data(invAu, invAu_data)
    call get_vector_data(invAv, invAv_data)
    
    ! Loop over cells
    print *, "P': cell loop"
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
        call set_face_location(face_loc, cell_mesh, local_idx, j)
        call get_face_area(face_loc, face_area)
        call get_face_normal(face_loc, face_normal)

        call get_boundary_status(face_loc, is_boundary)
        
        if (.not. is_boundary) then
          ! Interior face
          call set_neighbour_location(ngb_loc, self_loc, j)
          call get_global_index(ngb_loc, ngb_idx)
          call get_local_index(ngb_loc, idxnb)
          coeff_f = (1.0 / cell_mesh%h) * face_area

          call get_volume(self_loc, Vp)
          call get_volume(ngb_loc, Vnb)
          Vf = 0.5_accs_real * (Vp + Vnb)

          invAp = 0.5_accs_real * (invAu_data(local_idx) + invAv_data(local_idx))
          invAnb = 0.5_accs_real * (invAu_data(idxnb) + invAv_data(idxnb))
          invAf = 0.5_accs_real * (invAp + invAnb)

          coeff_f = -(Vf * invAf) * coeff_f
          
          coeff_p = coeff_p - coeff_f
          coeff_nb = coeff_f
          col = ngb_idx
        else
          ! XXX: Fixed velocity BC - no pressure correction
          col = -1
          coeff_nb = 0.0_accs_real
        endif
        call pack_entries(mat_coeffs, 1, j+1, row, col, coeff_nb)

      end do

      ! XXX: Need to fix pressure somewhere
      !!     Row is the global index - should be unique
      !!     Locate approximate centre of mesh (assuming a square)
      cps = int(sqrt(real(cell_mesh%nglobal)), accs_int)
      rcrit = (cps / 2) * (1 + cps)
      if (row == rcrit) then
        coeff_p = coeff_p + 1.0e30 ! Force diagonal to be huge -> zero solution (approximately).
        print *, "Fixed coeff_p", coeff_p, " at ", row
      end if
      
      ! Add the diagonal entry
      col = row
      call pack_entries(mat_coeffs, 1, 1, row, col, coeff_p)

      call pack_entries(vec_values, 1, self_idx, r)

      ! Set the values
      call set_values(mat_coeffs, M)
      call set_values(vec_values, vec)

      deallocate(mat_coeffs%rglob)
      deallocate(mat_coeffs%cglob)
      deallocate(mat_coeffs%val)
    end do

    print *, "P': restore invA"
    call restore_vector_data(invAu, invAu_data)
    call restore_vector_data(invAv, invAv_data)

    ! Assembly of coefficient matrix and source vector
    print *, "P': assemble matrix, RHS"
    call update(M)
    call update(vec)
    call finalise(M)
    
    ! Create linear solver
    print *, "P': create lin sys"
    call set_linear_system(lin_sys, vec, pp%vec, M, par_env)
    call create_solver(lin_sys, lin_solver)

    ! Solve the linear system
    print *, "P': solve"
    call solve(lin_solver)

    ! Clean up
    deallocate(lin_solver)
    
  end subroutine calculate_pressure_correction

  !> @brief Computes the per-cell mass imbalance, updating the face velocity flux as it does so.
  subroutine compute_mass_imbalance(par_env, cell_mesh, u, v, p, invAu, invAv, mf, b)

    class(parallel_environment), intent(in) :: par_env
    type(mesh), intent(in) :: cell_mesh !> The mesh object
    class(field), intent(inout) :: u       !> The x velocity component
    class(field), intent(inout) :: v       !> The y velocity component
    class(field), intent(inout) :: p       !> The pressure field
    class(vector), intent(in) :: invAu  !> The inverse x momentum equation diagonal coefficient
    class(vector), intent(in) :: invAv  !> The inverse y momentum equation diagonal coefficient
    class(field), intent(inout) :: mf   !> The face velocity flux
    class(vector), intent(inout) :: b   !> The per-cell mass imbalance

    type(vector_values) :: vec_values
    integer(accs_int) :: i !> Cell counter
    integer(accs_int) :: j !> Cell-face counter

    type(cell_locator) :: loc_p !> Central cell locator object
    type(face_locator) :: loc_f !> Face locator object

    integer(accs_int) :: idxp_g  !> Central cell global index
    real(accs_real) :: face_area !> Face area
    integer(accs_int) :: idxf    !> Face index
    integer(accs_int) :: nnb     !> Cell neighbour count

    real(accs_real), dimension(:), pointer :: mf_data     !> Data array for the mass flux
    real(accs_real), dimension(:), pointer :: u_data      !> Data array for x velocity component
    real(accs_real), dimension(:), pointer :: v_data      !> Data array for y velocity component
    real(accs_real), dimension(:), pointer :: p_data      !> Data array for pressure
    real(accs_real), dimension(:), pointer :: pgradx_data !> Data array for pressure x gradient
    real(accs_real), dimension(:), pointer :: pgrady_data !> Data array for pressure y gradient
    real(accs_real), dimension(:), pointer :: invAu_data  !> Data array for inverse x momentum
                                                          !! diagonal coefficient
    real(accs_real), dimension(:), pointer :: invAv_data  !> Data array for inverse y momentum
                                                          !! diagonal coefficient

    logical :: is_boundary            !> Boundary indicator
    type(neighbour_locator) :: loc_nb !> Neighbour cell locator object
    integer(accs_int) :: idxnb        !> Neighbour cell index
    
    real(accs_real) :: mib !> Cell mass imbalance

    allocate(vec_values%idx(1))
    allocate(vec_values%val(1))
    vec_values%mode = insert_mode

    ! First zero RHS
    call zero(b)

    ! Update vectors to make sure all data is up to date
    call update(u%vec)
    call update(v%vec)
    call update(p%vec)
    call update(p%gradx)
    call update(p%grady)
    call get_vector_data(mf%vec, mf_data)
    call get_vector_data(u%vec, u_data)
    call get_vector_data(v%vec, v_data)
    call get_vector_data(p%vec, p_data)
    call get_vector_data(p%gradx, pgradx_data)
    call get_vector_data(p%grady, pgrady_data)
    call get_vector_data(invAu, invAu_data)
    call get_vector_data(invAv, invAv_data)
    
    do i = 1, cell_mesh%nlocal
      call set_cell_location(loc_p, cell_mesh, i)
      call get_global_index(loc_p, idxp_g)
      call count_neighbours(loc_p, nnb)

      mib = 0.0_accs_real

      do j = 1, nnb
        call set_face_location(loc_f, cell_mesh, i, j)
        call get_face_area(loc_f, face_area)
        call get_local_index(loc_f, idxf)

        ! Check face orientation
        call get_boundary_status(loc_f, is_boundary)
        if (.not. is_boundary) then
          call set_neighbour_location(loc_nb, loc_p, j)
          call get_local_index(loc_nb, idxnb)
          if (idxnb < i) then
            face_area = -face_area
          else
            ! Compute mass flux through face
            mf_data(idxf) = calc_mass_flux(u_data, v_data, &
                 p_data, pgradx_data, pgrady_data, &
                 invAu_data, invAv_data, &
                 loc_f)
          end if
        end if
        
        mib = mib + mf_data(idxf) * face_area
      end do
      
      call pack_entries(vec_values, 1, idxp_g, mib)
      call set_values(vec_values, b)
    end do

    call restore_vector_data(mf%vec, mf_data)
    call restore_vector_data(u%vec, u_data)
    call restore_vector_data(v%vec, v_data)
    call restore_vector_data(p%vec, p_data)
    call restore_vector_data(p%gradx, pgradx_data)
    call restore_vector_data(p%grady, pgrady_data)
    call restore_vector_data(invAu, invAu_data)
    call restore_vector_data(invAv, invAv_data)
    ! Update vectors on exit (just in case)
    call update(u%vec)
    call update(v%vec)
    call update(p%vec)
    call update(p%gradx)
    call update(p%grady)

    call update(b)
    call update(mf%vec)

    mib = norm(b, 2)
    if (par_env%proc_id == par_env%root) then
      print *, "SIMPLE intermediate mass imbalance: ", mib
    end if
    
  end subroutine compute_mass_imbalance

  !> @brief Corrects the pressure field, using explicit underrelaxation
  subroutine update_pressure(pp, p)

    ! Arguments
    class(field), intent(in) :: pp
    class(field), intent(inout) :: p

    ! Local variables
    real(accs_real) :: alpha   !< Under-relaxation factor

    ! Set under-relaxation factor (todo: read this from input file)
    alpha = 0.1_accs_real

    call axpy(alpha, pp%vec, p%vec)
    
  end subroutine update_pressure

  !> @brief Corrects the velocity field using the pressure correction gradient
  subroutine update_velocity(cell_mesh, invAu, invAv, pp, u, v)

    use vec, only : zero_vector
    
    ! Arguments
    class(mesh), intent(in) :: cell_mesh
    class(vector), intent(in) :: invAu, invAv
    class(field), intent(inout) :: pp
    class(field), intent(inout) :: u, v

    ! First update gradients
    call zero_vector(pp%gradx)
    call zero_vector(pp%grady)
    call update_gradient(cell_mesh, pp)

    ! Multiply gradients by inverse diagonal coefficients
    call mult(invAu, pp%gradx)
    call mult(invAv, pp%grady)

    ! Compute correction source on velocity
    call calculate_momentum_pressure_source(cell_mesh, pp%gradx, u%vec)
    call calculate_momentum_pressure_source(cell_mesh, pp%grady, v%vec)

    call update(u%vec)
    call update(v%vec)
    
  end subroutine update_velocity

  !> @brief Corrects the face velocity flux using the pressure correction
  subroutine update_face_velocity(cell_mesh, pp, invAu, invAv, mf)

    type(mesh), intent(in) :: cell_mesh
    class(field), intent(inout) :: pp
    class(vector), intent(in) :: invAu
    class(vector), intent(in) :: invAv
    class(field), intent(inout) :: mf
    
    integer(accs_int) :: i

    real(accs_real) :: mf_prime
    real(accs_real), dimension(:), allocatable :: zero_arr
    real(accs_real), dimension(:), pointer :: mf_data
    real(accs_real), dimension(:), pointer :: pp_data
    real(accs_real), dimension(:), pointer :: invAu_data
    real(accs_real), dimension(:), pointer :: invAv_data

    type(cell_locator) :: loc_p
    integer(accs_int) :: nnb
    integer(accs_int) :: j
    type(face_locator) :: loc_f
    integer(accs_int) :: idxf

    logical :: is_boundary
    type(neighbour_locator) :: loc_nb
    integer(accs_int) :: idxnb
    
    ! Update vector to make sure data is up to date
    call update(pp%vec)
    call get_vector_data(pp%vec, pp_data)
    call get_vector_data(invAu, invAu_data)
    call get_vector_data(invAv, invAv_data)
    call get_vector_data(mf%vec, mf_data)
    
    allocate(zero_arr(size(pp_data)))
    zero_arr(:) = 0.0_accs_real

    ! XXX: This should really be a face loop
    do i = 1, cell_mesh%nlocal
      call set_cell_location(loc_p, cell_mesh, i)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb
        call set_face_location(loc_f, cell_mesh, i, j)
        call get_boundary_status(loc_f, is_boundary)
        if (.not. is_boundary) then
          call set_neighbour_location(loc_nb, loc_p, j)
          call get_local_index(loc_nb, idxnb)
          if (i < idxnb) then
            mf_prime = calc_mass_flux(zero_arr, zero_arr, &
                 pp_data, zero_arr, zero_arr, &
                 invAu_data, invAv_data, &
                 loc_f)

            call get_local_index(loc_f, idxf)
            mf_data(idxf) = mf_data(idxf) + mf_prime
          end if
        end if
      end do
    end do

    deallocate(zero_arr)

    call restore_vector_data(pp%vec, pp_data)
    call restore_vector_data(invAu, invAu_data)
    call restore_vector_data(invAv, invAv_data)
    call restore_vector_data(mf%vec, mf_data)

    call update(mf%vec)
    ! Update vector on exit (just in case)
    call update(pp%vec)
    
  end subroutine update_face_velocity

  subroutine calculate_scalars()


  end subroutine calculate_scalars


  subroutine check_convergence(converged)

    ! Arguments
    logical, intent(inout) :: converged

    converged = .false. ! XXX: temporary - force run for maximum iterations
    
  end subroutine check_convergence

  !> @brief Applies implicit underrelaxation to an equation
  !
  !> @description Extracts the diagonal coefficient of a matrix and divides by the URF, adding a
  !!              proportional explicit term to the RHS vector.
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

    print *, "UR: get diagonal vec"
    call finalise(M)
    call get_matrix_diagonal(M, diag)

    print *, "UR: get phi, diag, b"
    call get_vector_data(phi%vec, phi_data)
    call get_vector_data(diag, diag_data)
    call update(b)
    call get_vector_data(b, b_data)

    print *, "UR: apply UR"
    do i = 1, cell_mesh%nlocal
      diag_data(i) = diag_data(i) / alpha

      b_data(i) = b_data(i) + (1.0_accs_real - alpha) * diag_data(i) * phi_data(i)
    end do

    print *, "UR: Restore data"
    call restore_vector_data(phi%vec, phi_data)
    call restore_vector_data(diag, diag_data)
    call restore_vector_data(b, b_data)

    print *, "UR: Set matrix diagonal"
    call set_matrix_diagonal(diag, M)
    
  end subroutine underrelax
  
end submodule pv_coupling_simple
