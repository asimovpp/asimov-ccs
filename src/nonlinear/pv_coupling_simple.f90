!>  Submodule file pv_coupling_simple.smod
!
!>  Implementation of the SIMPLE algorithm for pressure-velocity coupling.

submodule (pv_coupling) pv_coupling_simple
#include "ccs_macros.inc"
  use kinds, only: ccs_real, ccs_int
  use types, only: vector_spec, ccs_vector, matrix_spec, ccs_matrix, equation_system, &
                   linear_solver, ccs_mesh, field, bc_config, vector_values, cell_locator, &
                   face_locator, neighbour_locator, matrix_values
  use fv, only: compute_fluxes, calc_mass_flux, update_gradient
  use vec, only: create_vector, vec_reciprocal, get_vector_data, restore_vector_data, scale_vec, &
       create_vector_values
  use mat, only: create_matrix, set_nnz, get_matrix_diagonal
  use utils, only: update, initialise, finalise, set_size, set_values, pack_entries, &
                   mult, zero, clear_entries, set_entry, set_row, set_mode, &
                   str, exit_print
  use utils, only: debug_print
  use solver, only: create_solver, solve, set_equation_system, axpy, norm
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
  !> @param[in]  mesh - the mesh
  !> @param[in,out] u, v   - fields containing velocity fields in x, y directions
  !> @param[in,out] p      - field containing pressure values
  !> @param[in,out] p_prime     - field containing pressure-correction values
  !> @param[in,out] mf     - field containing the face-centred velocity flux 
  module subroutine solve_nonlinear(par_env, mesh, cps, it_start, it_end, u, v, p, p_prime, mf)

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env
    type(ccs_mesh), intent(in) :: mesh
    integer(ccs_int), intent(in) :: cps, it_start, it_end
    class(field), intent(inout) :: u, v, p, p_prime, mf
    
    ! Local variables
    integer(ccs_int) :: i
    class(ccs_vector), allocatable :: source
    class(ccs_matrix), allocatable :: M
    class(ccs_vector), allocatable :: invAu, invAv
    
    type(vector_spec) :: vec_properties
    type(matrix_spec) :: mat_properties
    type(equation_system)    :: lin_system

    logical :: converged
    
    ! Initialise linear system
    call dprint("NONLINEAR: init")
    call initialise(vec_properties)
    call initialise(mat_properties)
    call initialise(lin_system)

    ! Create coefficient matrix
    call dprint("NONLINEAR: setup matrix")
    call set_size(par_env, mesh, mat_properties)
    call set_nnz(5, mat_properties)
    call create_matrix(mat_properties, M)

    ! Create RHS vector
    call dprint("NONLINEAR: setup RHS")
    call set_size(par_env, mesh, vec_properties)
    call create_vector(vec_properties, source)

    ! Create vectors for storing inverse of velocity central coefficients
    call dprint("NONLINEAR: setup ind coeff")
    call create_vector(vec_properties, invAu)
    call create_vector(vec_properties, invAv)
    
    ! Get pressure gradient
    call dprint("NONLINEAR: compute grad p")
    call update_gradient(mesh, p)

    outerloop: do i = it_start, it_end

      call dprint("NONLINEAR: iteration " // str(i))

      ! Solve momentum equation with guessed pressure and velocity fields (eq. 4)
      call dprint("NONLINEAR: guess velocity")
      call calculate_velocity(par_env, mesh, cps, mf, p, M, source, lin_system, u, v, invAu, invAv)

      ! Calculate pressure correction from mass imbalance (sub. eq. 11 into eq. 8)
      call dprint("NONLINEAR: mass imbalance")
      call compute_mass_imbalance(par_env, mesh, invAu, invAv, u, v, p, mf, source)
      call dprint("NONLINEAR: compute p'")
      call calculate_pressure_correction(par_env, mesh, invAu, invAv, M, source, lin_system, p_prime)
      
      ! Update velocity with velocity correction (eq. 6)
      call dprint("NONLINEAR: correct face velocity")
      call update_face_velocity(mesh, invAu, invAv, p_prime, mf)
      call dprint("NONLINEAR: correct velocity")
      call update_velocity(mesh, invAu, invAv, p_prime, u, v)

      ! Update pressure field with pressure correction
      call dprint("NONLINEAR: correct pressure")
      call update_pressure(p_prime, p)
      call dprint("NONLINEAR: compute gradp")
      call update_gradient(mesh, p)

      ! Todo:
      !call calculate_scalars()

      call check_convergence(converged)
      if (converged) then
        exit outerloop
      endif

    end do outerloop

  end subroutine solve_nonlinear

  !>  Computes the guessed velocity fields based on a frozen pressure field
  !
  !> @param[in]    par_env      - the parallel environment
  !> @param[in]    mesh         - the mesh
  !> @param[in]    cps          - cells per side of a square mesh (remove)
  !> @param[in]    mf           - the face velocity flux
  !> @param[in]    p            - the pressure field
  !> @param[inout] M            - matrix object
  !> @param[inout] vec          - vector object
  !> @param[inout] lin_sys      - linear system object
  !> @param[inout] u, v         - the x and y velocity fields
  !> @param[inout] invAu, invAv - vectors containing the inverse momentum coefficients
  !
  !> @description Given an initial guess of a pressure field form the momentum equations (as scalar
  !!              equations) and solve to obtain an intermediate velocity field u* that will not
  !!              satisfy continuity.
  subroutine calculate_velocity(par_env, mesh, cps, mf, p, M, vec, lin_sys, u, v, invAu, invAv)

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env
    type(ccs_mesh), intent(in)         :: mesh
    integer(ccs_int), intent(in)  :: cps
    class(field), intent(in) :: mf
    class(field), intent(in) :: p
    class(ccs_matrix), allocatable, intent(inout)  :: M
    class(ccs_vector), allocatable, intent(inout)  :: vec
    type(equation_system), intent(inout) :: lin_sys
    class(field), intent(inout)    :: u, v
    class(ccs_vector), intent(inout)    :: invAu, invAv

    
    ! u-velocity
    ! ----------

    ! TODO: Do boundaries properly
    !u%bcs%bc_type(:) = 0 !< Fixed zero BC
    !u%bcs%bc_type(4) = 1 !< Fixed one BC at lid
    call calculate_velocity_component(par_env, mesh, cps, mf, p, 1, M, vec, lin_sys, u, invAu)
    
    ! v-velocity
    ! ----------
    
    ! TODO: Do boundaries properly
    !v%bcs%bc_type(:) = 0 !< Fixed zero BC
    call calculate_velocity_component(par_env, mesh, cps, mf, p, 2, M, vec, lin_sys, v, invAv)

  end subroutine calculate_velocity

  subroutine calculate_velocity_component(par_env, mesh, cps, mf, p, component, M, vec, lin_sys, u, invAu)

    use case_config, only: velocity_relax

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env
    type(ccs_mesh), intent(in)         :: mesh
    integer(ccs_int), intent(in)  :: cps
    class(field), intent(in) :: mf
    class(field), intent(in) :: p
    integer(ccs_int), intent(in) :: component
    class(ccs_matrix), allocatable, intent(inout)  :: M
    class(ccs_vector), allocatable, intent(inout)  :: vec
    type(equation_system), intent(inout) :: lin_sys
    class(field), intent(inout)    :: u
    class(ccs_vector), intent(inout)    :: invAu
    
    ! Local variables
    class(linear_solver), allocatable :: lin_solver

    ! First zero matrix/RHS
    call zero(vec)
    call zero(M)
    
    ! Calculate fluxes and populate coefficient matrix
    call dprint("GV: compute u flux")
    call compute_fluxes(u, mf, mesh, cps, M, vec)

    ! Calculate pressure source term and populate RHS vector
    call dprint("GV: compute u gradp")
    if (component == 1) then
      call calculate_momentum_pressure_source(mesh, p%x_gradients, vec)
    else if (component == 2) then
      call calculate_momentum_pressure_source(mesh, p%y_gradients, vec)
    else
      call error_abort("Unsupported vector component: " // str(component))
    end if
    
    ! Underrelax the equations
    call dprint("GV: underrelax u")
    call underrelax(mesh, velocity_relax, u, invAu, M, vec)
    
    ! Store reciprocal of central coefficient
    call dprint("GV: get u diag")
    call get_matrix_diagonal(M, invAu)
    call vec_reciprocal(invAu)
    
    ! Assembly of coefficient matrix and source vector
    call dprint("GV: build u lin sys")
    call update(M)
    call update(vec)
    call update(invAu)
    call finalise(M)

    ! Create linear solver
    call set_equation_system(par_env, vec, u%values, M, lin_sys)
    call create_solver(lin_sys, lin_solver)

    ! Solve the linear system
    call dprint("GV: solve u")
    call solve(lin_solver)

    ! Clean up
    deallocate(lin_solver)
    
  end subroutine calculate_velocity_component
  
  !v Adds the momentum source due to pressure gradient
  subroutine calculate_momentum_pressure_source(mesh, p_gradients, vec)

    ! Arguments
    class(ccs_mesh), intent(in) :: mesh           !< the mesh
    class(ccs_vector), intent(in) :: p_gradients  !< the pressure gradient
    class(ccs_vector), intent(inout) :: vec       !< the momentum equation RHS vector

    ! Local variables
    type(vector_values) :: vec_values
    type(cell_locator) :: loc_p
    integer(ccs_int) :: global_index_p, index_p
    real(ccs_real) :: r
    real(ccs_real), dimension(:), pointer :: p_gradient_data

    real(ccs_real) :: V
    
    call create_vector_values(1_ccs_int, vec_values)
    call set_mode(add_mode, vec_values)

    ! Temporary storage for p values
    call get_vector_data(p_gradients, p_gradient_data)
    
    ! Loop over cells
    do index_p = 1, mesh%nlocal
      call clear_entries(vec_values)
      
      call set_cell_location(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)

      call get_volume(loc_p, V)
      
      r = -p_gradient_data(index_p) * V
      call set_row(global_index_p, vec_values)
      call set_entry(r, vec_values)
      call set_values(vec_values, vec)
    end do

    deallocate(vec_values%global_indices)
    deallocate(vec_values%values)

    call restore_vector_data(p_gradients, p_gradient_data)
    
  end subroutine calculate_momentum_pressure_source

  !>  Solves the pressure correction equation
  !
  !> @description Solves the pressure correction equation formed by the mass-imbalance.
  subroutine calculate_pressure_correction(par_env, mesh, invAu, invAv, M, vec, lin_sys, p_prime)

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env  !< the parallel environment
    class(ccs_mesh), intent(in) :: mesh                              !< the mesh
    class(ccs_vector), intent(in) :: invAu, invAv                    !< inverse diagonal momentum coefficients
    class(ccs_matrix), allocatable, intent(inout)  :: M              !< matrix object
    class(ccs_vector), allocatable, intent(inout)  :: vec            !< the RHS vector
    type(equation_system), intent(inout) :: lin_sys                  !< linear system object
    class(field), intent(inout) :: p_prime                           !< the pressure correction field

    ! Local variables
    type(matrix_values) :: mat_coeffs
    type(vector_values) :: vec_values
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f
    class(linear_solver), allocatable :: lin_solver
    integer(ccs_int) :: global_index_p, global_index_nb, index_p
    integer(ccs_int) :: j
    integer(ccs_int) :: nnb
    integer(ccs_int) :: row, col
    real(ccs_real) :: face_area
    real(ccs_real), dimension(ndim) :: face_normal
    real(ccs_real) :: r
    real(ccs_real) :: coeff_f, coeff_p, coeff_nb
    logical :: is_boundary

    real(ccs_real), dimension(:), pointer :: invAu_data
    real(ccs_real), dimension(:), pointer :: invAv_data
    
    real(ccs_real) :: Vp
    real(ccs_real) :: V_nb
    real(ccs_real) :: Vf
    real(ccs_real) :: invA_p
    real(ccs_real) :: invA_nb
    real(ccs_real) :: invA_f

    integer(ccs_int) :: index_nb

    integer(ccs_int) :: cps   ! Cells per side
    integer(ccs_int) :: rcrit ! Global index of approximate central cell
    
    ! First zero matrix
    call zero(M)

    ! The computed mass imbalance is +ve, to have a +ve diagonal coefficient we need to negate this.
    call dprint("P': negate RHS")
    call scale_vec(-1.0_ccs_real, vec)
    call update(vec)
    
    call create_vector_values(1_ccs_int, vec_values)
    call set_mode(add_mode, vec_values)

    mat_coeffs%setter_mode = insert_mode

    call update(M)
    
    call dprint("P': get invA")
    call get_vector_data(invAu, invAu_data)
    call get_vector_data(invAv, invAv_data)
    
    ! Loop over cells
    call dprint("P': cell loop")
    do index_p = 1, mesh%nlocal
      call clear_entries(vec_values)
      
      call set_cell_location(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)

      allocate(mat_coeffs%global_row_indices(1))
      allocate(mat_coeffs%global_col_indices(1 + nnb))
      allocate(mat_coeffs%values(1 + nnb))

      row = global_index_p
      coeff_p = 0.0_ccs_real
      r = 0.0_ccs_real

      ! Loop over faces
      do j = 1, nnb
        call set_face_location(mesh, index_p, j, loc_f)
        call get_face_area(loc_f, face_area)
        call get_face_normal(loc_f, face_normal)

        call get_boundary_status(loc_f, is_boundary)
        
        if (.not. is_boundary) then
          ! Interior face
          call set_neighbour_location(loc_p, j, loc_nb)
          call get_global_index(loc_nb, global_index_nb)
          call get_local_index(loc_nb, index_nb)
          coeff_f = (1.0 / mesh%h) * face_area

          call get_volume(loc_p, Vp)
          call get_volume(loc_nb, V_nb)
          Vf = 0.5_ccs_real * (Vp + V_nb)

          invA_p = 0.5_ccs_real * (invAu_data(index_p) + invAv_data(index_p))
          invA_nb = 0.5_ccs_real * (invAu_data(index_nb) + invAv_data(index_nb))
          invA_f = 0.5_ccs_real * (invA_p + invA_nb)

          coeff_f = -(Vf * invA_f) * coeff_f
          
          coeff_p = coeff_p - coeff_f
          coeff_nb = coeff_f
          col = global_index_nb
        else
          ! XXX: Fixed velocity BC - no pressure correction
          col = -1
          coeff_nb = 0.0_ccs_real
        endif
        call pack_entries(1, j+1, row, col, coeff_nb, mat_coeffs)

      end do

      ! XXX: Need to fix pressure somewhere
      !!     Row is the global index - should be unique
      !!     Locate approximate centre of mesh (assuming a square)
      cps = int(sqrt(real(mesh%nglobal)), ccs_int)
      rcrit = (cps / 2) * (1 + cps)
      if (row == rcrit) then
        coeff_p = coeff_p + 1.0e30 ! Force diagonal to be huge -> zero solution (approximately).
        call dprint("Fixed coeff_p" // str(coeff_p) // " at " // str(row))
      end if
      
      ! Add the diagonal entry
      col = row
      call pack_entries(1, 1, row, col, coeff_p, mat_coeffs)

      call set_row(global_index_p, vec_values)
      call set_entry(r, vec_values)

      ! Set the values
      call set_values(mat_coeffs, M)
      call set_values(vec_values, vec)

      deallocate(mat_coeffs%global_row_indices)
      deallocate(mat_coeffs%global_col_indices)
      deallocate(mat_coeffs%values)
    end do

    call dprint("P': restore invA")
    call restore_vector_data(invAu, invAu_data)
    call restore_vector_data(invAv, invAv_data)

    ! Assembly of coefficient matrix and source vector
    call dprint("P': assemble matrix, RHS")
    call update(M)
    call update(vec)
    call finalise(M)
    
    ! Create linear solver
    call dprint("P': create lin sys")
    call set_equation_system(par_env, vec, p_prime%values, M, lin_sys)
    call create_solver(lin_sys, lin_solver)

    ! Solve the linear system
    call dprint("P': solve")
    call solve(lin_solver)

    ! Clean up
    deallocate(lin_solver)
    
  end subroutine calculate_pressure_correction

  !>  Computes the per-cell mass imbalance, updating the face velocity flux as it does so.
  subroutine compute_mass_imbalance(par_env, mesh, invAu, invAv, u, v, p, mf, b)

    class(parallel_environment), intent(in) :: par_env
    type(ccs_mesh), intent(in) :: mesh      !< The mesh object
    class(ccs_vector), intent(in) :: invAu  !< The inverse x momentum equation diagonal coefficient
    class(ccs_vector), intent(in) :: invAv  !< The inverse y momentum equation diagonal coefficient
    class(field), intent(inout) :: u        !< The x velocity component
    class(field), intent(inout) :: v        !< The y velocity component
    class(field), intent(inout) :: p        !< The pressure field
    class(field), intent(inout) :: mf       !< The face velocity flux
    class(ccs_vector), intent(inout) :: b   !< The per-cell mass imbalance

    type(vector_values) :: vec_values
    integer(ccs_int) :: i !< Cell counter
    integer(ccs_int) :: j !< Cell-face counter

    type(cell_locator) :: loc_p !< Central cell locator object
    type(face_locator) :: loc_f !< Face locator object

    integer(ccs_int) :: global_index_p  !< Central cell global index
    real(ccs_real) :: face_area !< Face area
    integer(ccs_int) :: index_f    !< Face index
    integer(ccs_int) :: nnb     !< Cell neighbour count

    real(ccs_real), dimension(:), pointer :: mf_data     !< Data array for the mass flux
    real(ccs_real), dimension(:), pointer :: u_data      !< Data array for x velocity component
    real(ccs_real), dimension(:), pointer :: v_data      !< Data array for y velocity component
    real(ccs_real), dimension(:), pointer :: p_data      !< Data array for pressure
    real(ccs_real), dimension(:), pointer :: dpdx_data !< Data array for pressure x gradient
    real(ccs_real), dimension(:), pointer :: dpdy_data !< Data array for pressure y gradient
    real(ccs_real), dimension(:), pointer :: invAu_data  !< Data array for inverse x momentum
                                                          !! diagonal coefficient
    real(ccs_real), dimension(:), pointer :: invAv_data  !< Data array for inverse y momentum
                                                          !! diagonal coefficient

    logical :: is_boundary            !< Boundary indicator
    type(neighbour_locator) :: loc_nb !< Neighbour cell locator object
    integer(ccs_int) :: index_nb      !< Neighbour cell index
    
    real(ccs_real) :: mib !< Cell mass imbalance

    call create_vector_values(1_ccs_int, vec_values)
    call set_mode(insert_mode, vec_values)
    
    ! First zero RHS
    call zero(b)

    ! Update vectors to make sure all data is up to date
    call update(u%values)
    call update(v%values)
    call update(p%values)
    call update(p%x_gradients)
    call update(p%y_gradients)
    call get_vector_data(mf%values, mf_data)
    call get_vector_data(u%values, u_data)
    call get_vector_data(v%values, v_data)
    call get_vector_data(p%values, p_data)
    call get_vector_data(p%x_gradients, dpdx_data)
    call get_vector_data(p%y_gradients, dpdy_data)
    call get_vector_data(invAu, invAu_data)
    call get_vector_data(invAv, invAv_data)
    
    do i = 1, mesh%nlocal
      call clear_entries(vec_values)
      
      call set_cell_location(mesh, i, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)

      mib = 0.0_ccs_real

      do j = 1, nnb
        call set_face_location(mesh, i, j, loc_f)
        call get_face_area(loc_f, face_area)
        call get_local_index(loc_f, index_f)

        ! Check face orientation
        call get_boundary_status(loc_f, is_boundary)
        if (.not. is_boundary) then
          call set_neighbour_location(loc_p, j, loc_nb)
          call get_local_index(loc_nb, index_nb)
          if (index_nb < i) then
            face_area = -face_area
          else
            ! Compute mass flux through face
            mf_data(index_f) = calc_mass_flux(u_data, v_data, &
                 p_data, dpdx_data, dpdy_data, &
                 invAu_data, invAv_data, &
                 loc_f)
          end if
        end if
        
        mib = mib + mf_data(index_f) * face_area
      end do

      call set_row(global_index_p, vec_values)
      call set_entry(mib, vec_values)
      call set_values(vec_values, b)
    end do

    call restore_vector_data(mf%values, mf_data)
    call restore_vector_data(u%values, u_data)
    call restore_vector_data(v%values, v_data)
    call restore_vector_data(p%values, p_data)
    call restore_vector_data(p%x_gradients, dpdx_data)
    call restore_vector_data(p%y_gradients, dpdy_data)
    call restore_vector_data(invAu, invAu_data)
    call restore_vector_data(invAv, invAv_data)
    ! Update vectors on exit (just in case)
    call update(u%values)
    call update(v%values)
    call update(p%values)
    call update(p%x_gradients)
    call update(p%y_gradients)

    call update(b)
    call update(mf%values)

    mib = norm(b, 2)
    if (par_env%proc_id == par_env%root) then
      print *, "SIMPLE intermediate mass imbalance: " // str(mib)
    end if
    
  end subroutine compute_mass_imbalance

  !>  Corrects the pressure field, using explicit underrelaxation
  subroutine update_pressure(p_prime, p)

    use case_config, only: pressure_relax

    ! Arguments
    class(field), intent(in) :: p_prime   !< pressure correction
    class(field), intent(inout) :: p      !< the pressure field being corrected

    call axpy(pressure_relax, p_prime%values, p%values)
    
  end subroutine update_pressure

  !>  Corrects the velocity field using the pressure correction gradient
  subroutine update_velocity(mesh, invAu, invAv, p_prime, u, v)

    use vec, only : zero_vector
    
    ! Arguments
    class(ccs_mesh), intent(in) :: mesh             !< The mesh
    class(ccs_vector), intent(in) :: invAu, invAv   !< The inverse x, y momentum equation diagonal coefficients
    class(field), intent(inout) :: p_prime          !< The pressure correction
    class(field), intent(inout) :: u, v             !< The x, y velocities being corrected

    ! First update gradients
    call zero_vector(p_prime%x_gradients)
    call zero_vector(p_prime%y_gradients)
    call update_gradient(mesh, p_prime)

    ! Multiply gradients by inverse diagonal coefficients
    call mult(invAu, p_prime%x_gradients)
    call mult(invAv, p_prime%y_gradients)

    ! Compute correction source on velocity
    call calculate_momentum_pressure_source(mesh, p_prime%x_gradients, u%values)
    call calculate_momentum_pressure_source(mesh, p_prime%y_gradients, v%values)

    call update(u%values)
    call update(v%values)
    
  end subroutine update_velocity

  !>  Corrects the face velocity flux using the pressure correction
  subroutine update_face_velocity(mesh, invAu, invAv, p_prime, mf)

    type(ccs_mesh), intent(in) :: mesh              !< The mesh
    class(ccs_vector), intent(in) :: invAu, invAv   !< The inverse x, y momentum equation diagonal coefficients
    class(field), intent(inout) :: p_prime          !< The pressure correction
    class(field), intent(inout) :: mf               !< The face velocity being corrected
    
    integer(ccs_int) :: i

    real(ccs_real) :: mf_prime
    real(ccs_real), dimension(:), allocatable :: zero_arr
    real(ccs_real), dimension(:), pointer :: mf_data
    real(ccs_real), dimension(:), pointer :: pp_data
    real(ccs_real), dimension(:), pointer :: invAu_data
    real(ccs_real), dimension(:), pointer :: invAv_data

    type(cell_locator) :: loc_p
    integer(ccs_int) :: nnb
    integer(ccs_int) :: j
    type(face_locator) :: loc_f
    integer(ccs_int) :: index_f

    logical :: is_boundary
    type(neighbour_locator) :: loc_nb
    integer(ccs_int) :: index_nb
    
    ! Update vector to make sure data is up to date
    call update(p_prime%values)
    call get_vector_data(p_prime%values, pp_data)
    call get_vector_data(invAu, invAu_data)
    call get_vector_data(invAv, invAv_data)
    call get_vector_data(mf%values, mf_data)
    
    allocate(zero_arr(size(pp_data)))
    zero_arr(:) = 0.0_ccs_real

    ! XXX: This should really be a face loop
    do i = 1, mesh%nlocal
      call set_cell_location(mesh, i, loc_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb
        call set_face_location(mesh, i, j, loc_f)
        call get_boundary_status(loc_f, is_boundary)
        if (.not. is_boundary) then
          call set_neighbour_location(loc_p, j, loc_nb)
          call get_local_index(loc_nb, index_nb)
          if (i < index_nb) then
            mf_prime = calc_mass_flux(zero_arr, zero_arr, &
                 pp_data, zero_arr, zero_arr, &
                 invAu_data, invAv_data, &
                 loc_f)

            call get_local_index(loc_f, index_f)
            mf_data(index_f) = mf_data(index_f) + mf_prime
          end if
        end if
      end do
    end do

    deallocate(zero_arr)

    call restore_vector_data(p_prime%values, pp_data)
    call restore_vector_data(invAu, invAu_data)
    call restore_vector_data(invAv, invAv_data)
    call restore_vector_data(mf%values, mf_data)

    call update(mf%values)
    ! Update vector on exit (just in case)
    call update(p_prime%values)
    
  end subroutine update_face_velocity

  subroutine check_convergence(converged)

    ! Arguments
    logical, intent(inout) :: converged

    converged = .false. ! XXX: temporary - force run for maximum iterations
    
  end subroutine check_convergence

  !>  Applies implicit underrelaxation to an equation
  !
  !> @description Extracts the diagonal coefficient of a matrix and divides by the URF, adding a
  !!              proportional explicit term to the RHS vector.
  subroutine underrelax(mesh, alpha, phi, diag, M, b)

    use mat, only : set_matrix_diagonal
    
    type(ccs_mesh), intent(in) :: mesh
    real(ccs_real), intent(in) :: alpha
    class(field), intent(in) :: phi
    class(ccs_vector), intent(inout) :: diag
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b

    real(ccs_real), dimension(:), pointer :: diag_data
    real(ccs_real), dimension(:), pointer :: phi_data
    real(ccs_real), dimension(:), pointer :: b_data

    integer(ccs_int) :: i

    call dprint("UR: get diagonal vec")
    call finalise(M)
    call get_matrix_diagonal(M, diag)

    call dprint("UR: get phi, diag, b")
    call get_vector_data(phi%values, phi_data)
    call get_vector_data(diag, diag_data)
    call update(b)
    call get_vector_data(b, b_data)

    call dprint("UR: apply UR")
    do i = 1, mesh%nlocal
      diag_data(i) = diag_data(i) / alpha

      b_data(i) = b_data(i) + (1.0_ccs_real - alpha) * diag_data(i) * phi_data(i)
    end do

    call dprint("UR: Restore data")
    call restore_vector_data(phi%values, phi_data)
    call restore_vector_data(diag, diag_data)
    call restore_vector_data(b, b_data)

    call dprint("UR: Set matrix diagonal")
    call set_matrix_diagonal(diag, M)
    
  end subroutine underrelax
  
end submodule pv_coupling_simple
