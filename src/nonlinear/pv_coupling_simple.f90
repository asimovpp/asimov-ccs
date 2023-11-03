!v Submodule file pv_coupling_simple.smod
!
!  Implementation of the SIMPLE algorithm for pressure-velocity coupling.

submodule(pv_coupling) pv_coupling_simple
#include "ccs_macros.inc"
  use case_config, only: velocity_solver_method_name, velocity_solver_precon_name, &
                         pressure_solver_method_name, pressure_solver_precon_name
  use types, only: vector_spec, ccs_vector, matrix_spec, ccs_matrix, equation_system, &
                   linear_solver, bc_config, vector_values, cell_locator, &
                   face_locator, neighbour_locator, matrix_values, matrix_values_spec, upwind_field
  use fv, only: compute_fluxes, calc_mass_flux, update_gradient
  use vec, only: create_vector, vec_reciprocal, get_vector_data, restore_vector_data, scale_vec, &
                 create_vector_values, set_vector_location, zero_vector, vec_aypx, &
                 mult_vec_vec
  use mat, only: create_matrix, set_nnz, get_matrix_diagonal, set_matrix_values_spec_nrows, &
                 set_matrix_values_spec_ncols, create_matrix_values, mat_vec_product
  use utils, only: update, initialise, finalise, set_size, set_values, &
                   mult, zero, clear_entries, set_entry, set_row, set_col, set_mode, &
                   str, exit_print

  use utils, only: debug_print, get_field, get_fluid_solver_selector
  use solver, only: create_solver, solve, set_equation_system, axpy, norm, set_solver_method, set_solver_precon
  use constants, only: insert_mode, add_mode, ndim, cell, field_u, field_v, field_w, field_p, field_p_prime, &
                        field_mf, field_viscosity, field_density         
  use meshing, only: get_face_area, get_global_index, get_local_index, count_neighbours, &
                     get_boundary_status, get_face_normal, create_neighbour_locator, create_face_locator, &
                     create_cell_locator, get_volume, get_distance, &
                     get_local_num_cells, get_face_interpolation, &
                     get_global_num_cells, &
                     get_max_faces
  use scalars, only: update_scalars
  use timestepping, only: update_old_values, finalise_timestep, get_current_step, get_current_time
  use bc_constants, only: bc_type_dirichlet

  implicit none

  integer(ccs_int), save :: varp = 0

contains

  !> Solve Navier-Stokes equations using the SIMPLE algorithm
  module subroutine solve_nonlinear(par_env, mesh, it_start, it_end, res_target, &
                                    flow_solver_selector, flow)

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env   !< parallel environment
    type(ccs_mesh), intent(in) :: mesh                                !< the mesh
    integer(ccs_int), intent(in) :: it_start
    integer(ccs_int), intent(in) :: it_end
    real(ccs_real), intent(in) :: res_target                          !< Target residual
    type(fluid_solver_selector), intent(in) :: flow_solver_selector   !< determines which fluid fields need to be solved for
    type(fluid), intent(inout) :: flow                                !< The structure containting all the fluid fields

    ! Local variables
    integer(ccs_int) :: i
    class(ccs_vector), allocatable :: source
    class(ccs_matrix), allocatable :: M
    class(ccs_vector), allocatable :: invAu, invAv, invAw
    class(ccs_vector), allocatable :: res
    real(ccs_real), dimension(:), allocatable :: residuals
    integer(ccs_int) :: max_faces ! The maximum number of faces per cell

    type(vector_spec) :: vec_properties
    type(matrix_spec) :: mat_properties
    type(equation_system) :: lin_system
    class(linear_solver), allocatable :: lin_solverP !< Pressure correction linear solver

    logical :: converged

    integer(ccs_int) :: nvar ! Number of flow variables to solve
    integer(ccs_int) :: ivar ! Counter for flow variables

    logical :: u_sol !< solve u velocity field
    logical :: v_sol !< solve v velocity field
    logical :: w_sol !< solve w velocity field
    logical :: p_sol !< solve pressure field
    class(field), pointer :: u       !< velocity fields in x direction
    class(field), pointer :: v       !< velocity fields in y direction
    class(field), pointer :: w       !< velocity field in z direction
    class(field), pointer :: p       !< field containing pressure values
    class(field), pointer :: p_prime !< field containing pressure-correction values
    class(field), pointer :: mf      !< field containing the face-centred velocity flux
    class(field), pointer :: viscosity !< field containing the viscosity
    class(field), pointer :: density !< field containing the density

    call get_field(flow, field_u, u)
    call get_field(flow, field_v, v)
    call get_field(flow, field_w, w)
    call get_field(flow, field_p, p)
    call get_field(flow, field_p_prime, p_prime)
    call get_field(flow, field_mf, mf)
    call get_field(flow, field_viscosity, viscosity)
    call get_field(flow, field_density, density)
    call get_fluid_solver_selector(flow_solver_selector, field_u, u_sol)
    call get_fluid_solver_selector(flow_solver_selector, field_v, v_sol)
    call get_fluid_solver_selector(flow_solver_selector, field_w, w_sol)
    call get_fluid_solver_selector(flow_solver_selector, field_p, p_sol)

    ! Initialising SIMPLE solver
    nvar = 0
    ivar = 0
    converged = .false.

    call update_old_values(u)
    call update_old_values(v)
    call update_old_values(w)

    ! Initialise linear system
    call dprint("NONLINEAR: init")
    call initialise(vec_properties)
    call initialise(mat_properties)
    call initialise(lin_system)

    ! Create coefficient matrix
    call dprint("NONLINEAR: setup matrix")
    call get_max_faces(mesh, max_faces)
    call set_size(par_env, mesh, mat_properties)
    call set_nnz(max_faces + 1, mat_properties)
    call create_matrix(mat_properties, M)

    ! Create RHS vector
    call dprint("NONLINEAR: setup RHS")
    call set_size(par_env, mesh, vec_properties)
    call create_vector(vec_properties, source)

    ! Create vectors for storing inverse of velocity central coefficients
    call dprint("NONLINEAR: setup ind coeff")
    call create_vector(vec_properties, invAu)
    call create_vector(vec_properties, invAv)
    call create_vector(vec_properties, invAw)

    ! Create vectors for storing residuals
    call dprint("NONLINEAR: setup residuals")
    call create_vector(vec_properties, res)
    if (u_sol) nvar = nvar + 1
    if (v_sol) nvar = nvar + 1
    if (w_sol) nvar = nvar + 1
    if (p_sol) nvar = nvar + 2 ! (Pressure residual & mass imbalance)
    allocate (residuals(2 * nvar))
    residuals(:) = 0.0_ccs_real

    ! Get pressure gradient
    call dprint("NONLINEAR: compute gradients")
    call update_gradient(mesh, p)
    if (u_sol) call update_gradient(mesh, u)
    if (v_sol) call update_gradient(mesh, v)
    if (w_sol) call update_gradient(mesh, w)

    outerloop: do i = it_start, it_end
      call dprint("NONLINEAR: iteration " // str(i))

      ! Solve momentum equation with guessed pressure and velocity fields (eq. 4)
      call dprint("NONLINEAR: guess velocity")
      call calculate_velocity(par_env, mesh, flow, flow_solver_selector, ivar, M, source, &
                              lin_system, invAu, invAv, invAw, res, residuals)

            ! Calculate pressure correction from mass imbalance (sub. eq. 11 into eq. 8)
      call dprint("NONLINEAR: mass imbalance")
      call compute_mass_imbalance(mesh, invAu, invAv, invAw, ivar, flow, source, residuals)
      call dprint("NONLINEAR: compute p'")
      call calculate_pressure_correction(par_env, mesh, invAu, invAv, invAw, M, source, lin_system, p_prime, lin_solverP)

      ! Update velocity with velocity correction (eq. 6)
      call dprint("NONLINEAR: correct face velocity")
      call update_face_velocity(mesh, invAu, invAv, invAw, p_prime, mf, res, residuals)
      call dprint("NONLINEAR: correct velocity")
      call update_velocity(mesh, invAu, invAv, invAw, flow)

      ! Update pressure field with pressure correction
      call dprint("NONLINEAR: correct pressure")
      call update_pressure(mesh, p_prime, p)

      !< density values are in single digits (same as i/p)

      ! Transport scalars
      ! XXX: Should we distinguish active scalars (update in non-linear loop) and passive scalars
      !      (single update per timestep)?
      call update_scalars(par_env, mesh, flow)

      !< density values change to exponential here after update

      call check_convergence(par_env, i, residuals, res_target, &
                             flow_solver_selector, converged)
      if (converged) then
        call dprint("NONLINEAR: converged!")
        if (par_env%proc_id == par_env%root) then
          write (*, *)
          write (*, '(a)') 'Converged!'
          write (*, *)
        end if
        exit outerloop
      end if

      ! density values are in exponential
    end do outerloop

    deallocate (lin_solverP)
    call finalise_timestep()

    ! Free up memory
    deallocate (residuals)

  end subroutine solve_nonlinear

  !v Computes the guessed velocity fields based on a frozen pressure field
  !
  !  Given an initial guess of a pressure field form the momentum equations (as scalar
  !  equations) and solve to obtain an intermediate velocity field u* that will not
  !  satisfy continuity.
  subroutine calculate_velocity(par_env, mesh, flow, flow_solver_selector, ivar, M, vec, &
                                lin_sys, invAu, invAv, invAw, res, residuals)

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env !< the parallel environment
    type(ccs_mesh), intent(in) :: mesh                   !< the mesh
    type(fluid_solver_selector), intent(in) :: flow_solver_selector
    type(fluid), intent(inout) :: flow
    integer(ccs_int), intent(inout) :: ivar              !< flow variable counter
    class(ccs_matrix), allocatable, intent(inout) :: M   !< matrix object
    class(ccs_vector), allocatable, intent(inout) :: vec !< vector object
    type(equation_system), intent(inout) :: lin_sys      !< linear system object
    class(ccs_vector), intent(inout) :: invAu            !< vector containing the inverse x momentum coefficients
    class(ccs_vector), intent(inout) :: invAv            !< vector containing the inverse y momentum coefficients
    class(ccs_vector), intent(inout) :: invAw            !< vector containing the inverse z momentum coefficients
    class(ccs_vector), intent(inout) :: res              !< residual field
    real(ccs_real), dimension(:), intent(inout) :: residuals !< RMS and L-inf of residuals for each flow variable

    ! Local variables
    logical, save :: first_time = .true.
    integer(ccs_int), save :: varu = 0
    integer(ccs_int), save :: varv = 0
    integer(ccs_int), save :: varw = 0

    logical :: u_sol
    logical :: v_sol
    logical :: w_sol
    class(field), pointer :: u
    class(field), pointer :: v
    class(field), pointer :: w
    class(field), pointer :: p

    call get_field(flow, field_u, u)
    call get_field(flow, field_v, v)
    call get_field(flow, field_w, w)
    call get_field(flow, field_p, p)
    call get_fluid_solver_selector(flow_solver_selector, field_u, u_sol)
    call get_fluid_solver_selector(flow_solver_selector, field_v, v_sol)
    call get_fluid_solver_selector(flow_solver_selector, field_w, w_sol)

    ! Set flow variable identifiers (for residuals)
    if (first_time) then
      if (u_sol) then
        ivar = ivar + 1
        varu = ivar
      end if
      if (v_sol) then
        ivar = ivar + 1
        varv = ivar
      end if
      if (w_sol) then
        ivar = ivar + 1
        varw = ivar
      end if
      first_time = .false.
    end if

    ! u-velocity
    ! ----------
    if (u_sol) then
      call calculate_velocity_component(flow, par_env, varu, mesh, p, 1, M, vec, lin_sys, u, invAu, res, residuals)
    end if

    ! v-velocity
    ! ----------
    if (v_sol) then
      call calculate_velocity_component(flow, par_env, varv, mesh, p, 2, M, vec, lin_sys, v, invAv, res, residuals)
    end if

    ! w-velocity
    ! ----------
    if (w_sol) then
      call calculate_velocity_component(flow, par_env, varw, mesh, p, 3, M, vec, lin_sys, w, invAw, res, residuals)
    end if

  end subroutine calculate_velocity

  subroutine calculate_velocity_component(flow, par_env, ivar, mesh, p, component, M, vec, lin_sys, u, invAu, input_res, residuals)

    use case_config, only: velocity_relax
    use timestepping, only: apply_timestep

    ! Arguments
    type(fluid), intent(inout) :: flow
    class(parallel_environment), allocatable, intent(in) :: par_env
    integer(ccs_int), intent(in) :: ivar
    type(ccs_mesh), intent(in) :: mesh
    class(field), pointer :: mf
    class(field), pointer :: viscosity
    class(field), pointer :: density
    class(field), intent(inout) :: p
    integer(ccs_int), intent(in) :: component
    class(ccs_matrix), allocatable, intent(inout) :: M
    class(ccs_vector), allocatable, intent(inout) :: vec
    type(equation_system), intent(inout) :: lin_sys
    class(field), target, intent(inout) :: u
    class(ccs_vector), intent(inout) :: invAu
    class(ccs_vector), target, intent(inout) :: input_res
    class(ccs_vector), pointer :: res
    real(ccs_real), dimension(:), intent(inout) :: residuals

    ! Local variables
    class(linear_solver), allocatable :: lin_solver
    integer(ccs_int) :: nvar ! Number of flow variables to solve
    integer(ccs_int) :: global_num_cells

    ! First zero matrix/RHS
    call zero(vec)
    call zero(M)
    
    ! Select either field residuals or reuse 'input_res'
    if (allocated(u%residuals)) then
      res => u%residuals
    else
      res => input_res
    end if

    ! Zero residual vector
    call zero(res)

    ! Calculate fluxes and populate coefficient matrix
    if (component == 1) then
      call dprint("GV: compute u flux")
    else if (component == 2) then
      call dprint("GV: compute v flux")
    else if (component == 3) then
      call dprint("GV: compute w flux")
    else
      call error_abort("Unsupported vector component: " // str(component))
    end if

    call get_field(flow, field_mf, mf)
    call get_field(flow, field_viscosity, viscosity)
    call get_field(flow, field_density, density)
    
    call compute_fluxes(u, mf, viscosity, density, mesh, component, M, vec)

    call apply_timestep(mesh, u, invAu, M, vec)

    ! Calculate pressure source term and populate RHS vector
    call dprint("GV: compute u gradp")
    if (component == 1) then
      call calculate_momentum_pressure_source(mesh, p%x_gradients, vec)
    else if (component == 2) then
      call calculate_momentum_pressure_source(mesh, p%y_gradients, vec)
    else if (component == 3) then
      call calculate_momentum_pressure_source(mesh, p%z_gradients, vec)
    end if

    !calculate viscous source term and populate RHS vector 
    call dprint("compute viscous souce term")
    call calculate_momentum_viscous_source (mesh, flow, component, vec)

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

    ! Compute residual
    call mat_vec_product(M, u%values, res)
    call vec_aypx(vec, -1.0_ccs_real, res)
    ! Stores RMS of residuals
    call get_global_num_cells(mesh, global_num_cells)
    residuals(ivar) = norm(res, 2) / sqrt(real(global_num_cells))
    ! Stores Linf norm of residuals
    nvar = int(size(residuals) / 2_ccs_int)
    residuals(ivar + nvar) = norm(res, 0)

    ! Create linear solver
    if (allocated(u%values%name)) then
      call set_equation_system(par_env, vec, u%values, M, lin_sys, u%values%name)
    else
      call set_equation_system(par_env, vec, u%values, M, lin_sys)
    end if
    call create_solver(lin_sys, lin_solver)

    ! Customise linear solver
    call set_solver_method(velocity_solver_method_name, lin_solver)
    call set_solver_precon(velocity_solver_precon_name, lin_solver)

    ! Solve the linear system
    call dprint("GV: solve u")
    call solve(lin_solver)

    ! Clean up
    deallocate (lin_solver)

  end subroutine calculate_velocity_component

  !v Adds the momentum source due to pressure gradient
  subroutine calculate_momentum_pressure_source(mesh, p_gradients, vec)

    ! Arguments
    class(ccs_mesh), intent(in) :: mesh             !< the mesh
    class(ccs_vector), intent(inout) :: p_gradients !< the pressure gradient
    class(ccs_vector), intent(inout) :: vec         !< the momentum equation RHS vector

    ! Local variables
    type(vector_values) :: vec_values
    type(cell_locator) :: loc_p
    integer(ccs_int) :: global_index_p, index_p
    real(ccs_real) :: r
    real(ccs_real), dimension(:), pointer :: p_gradient_data
    integer(ccs_int) :: local_num_cells
    real(ccs_real) :: V

    call create_vector_values(1_ccs_int, vec_values)
    call set_mode(add_mode, vec_values)

    ! Temporary storage for p values
    call get_vector_data(p_gradients, p_gradient_data)

    ! Loop over cells
    call get_local_num_cells(mesh, local_num_cells)
    do index_p = 1, local_num_cells
      call clear_entries(vec_values)

      call create_cell_locator(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)

      call get_volume(loc_p, V)

      r = -p_gradient_data(index_p) * V
      call set_row(global_index_p, vec_values)
      call set_entry(r, vec_values)
      call set_values(vec_values, vec)
    end do

    deallocate (vec_values%global_indices)
    deallocate (vec_values%values)

    call restore_vector_data(p_gradients, p_gradient_data)

  end subroutine calculate_momentum_pressure_source

  !v Adds the momentum source due to variation in viscosity
  subroutine calculate_momentum_viscous_source (mesh, flow, component, vec)
    type(ccs_mesh), intent(in) :: mesh  !< the mesh
    type(fluid), intent(inout) :: flow
    integer(ccs_int), intent(in) :: component   !< integer indicating direction of velocity field component
    class(ccs_vector), allocatable, intent(inout) :: vec !< the momentum equation RHS vector
    class(field), pointer :: u  ! x-component of velocity 
    class(field), pointer :: v  ! y-component of velocity
    class(field), pointer :: w  ! z-component of velocity
    class(field), pointer :: viscosity
 
    ! Local variables
    type(vector_values) :: vec_values
    type(cell_locator) :: loc_p
    integer(ccs_int) :: global_index_p, index_p
    real(ccs_real) :: r1, r2, r3  ! variables involved in calculating the source term
    real(ccs_real), dimension(:), pointer :: dux_data, dvx_data, dwx_data
    real(ccs_real), dimension(:), pointer :: duy_data, dvy_data, dwy_data
    real(ccs_real), dimension(:), pointer :: duz_data, dvz_data, dwz_data
    real(ccs_real), dimension(3) :: duvw, duvwp, duvwf
    integer(ccs_int) :: local_num_cells
    real(ccs_real) :: Vol
    integer(ccs_int) :: nnb
    integer(ccs_int) :: j
    integer(ccs_int) :: index_nb
    type(neighbour_locator) :: loc_nb
    logical :: is_boundary
    type(face_locator) :: loc_f
    real(ccs_real), dimension(ndim) :: face_normal
    real(ccs_real) :: interpolation_factor   
    real(ccs_real) :: viscosity_face
    real(ccs_real) :: face_area
    real(ccs_real), dimension(:), pointer :: viscosity_data

    call get_field(flow, field_u, u)
    call get_field(flow, field_v, v)
    call get_field(flow, field_w, w)
    call get_field(flow, field_viscosity, viscosity)

    call create_vector_values(1_ccs_int, vec_values)
    call set_mode(add_mode, vec_values)

    ! Extracting the necessary data    
    call get_vector_data(viscosity%values, viscosity_data)
    ! x-component velocity gradient values
    call get_vector_data(u%x_gradients, dux_data)
    call get_vector_data(v%x_gradients, dvx_data)
    call get_vector_data(w%x_gradients, dwx_data)
    ! y-component velocity gradient values
    call get_vector_data(u%y_gradients, duy_data)
    call get_vector_data(v%y_gradients, dvy_data)
    call get_vector_data(w%y_gradients, dwy_data)
    ! z-component velocity gradient values
    call get_vector_data(u%z_gradients, duz_data)
    call get_vector_data(v%z_gradients, dvz_data)
    call get_vector_data(w%z_gradients, dwz_data)

    ! Loop over cells
    call get_local_num_cells(mesh, local_num_cells)
    do index_p = 1, local_num_cells

      call clear_entries(vec_values)
      call create_cell_locator(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call get_volume(loc_p, Vol)
      call count_neighbours(loc_p, nnb)

      r1=0.0_ccs_real
      r2=0.0_ccs_real

      do j=1,nnb
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_local_index(loc_nb, index_nb)
        call get_boundary_status(loc_nb, is_boundary)
        call create_face_locator(mesh, index_p, j, loc_f)
        call get_face_normal(loc_f, face_normal)
        call get_face_area(loc_f, face_area)
        call get_face_interpolation(loc_f, interpolation_factor)

        !evaluating gradients for neighbouring cell
        if(component==1) then ! x-component of velocity
          ! present cell gradients
          duvwp(1)=dux_data(index_p)
          duvwp(2)=dvx_data(index_p)
          duvwp(3)=dwx_data(index_p)
          ! neighbouring cell gradients
          duvwf(1)=dux_data(index_nb)
          duvwf(2)=dvx_data(index_nb)
          duvwf(3)=dwx_data(index_nb)

          if(.not.is_boundary) then ! no boundary face
            duvw(1)=(interpolation_factor*duvwp(1))+((1.0_ccs_real-interpolation_factor)*duvwf(1))
            duvw(2)=(interpolation_factor*duvwp(2))+((1.0_ccs_real-interpolation_factor)*duvwf(2))
            duvw(3)=(interpolation_factor*duvwp(3))+((1.0_ccs_real-interpolation_factor)*duvwf(3))
            viscosity_face=(interpolation_factor*viscosity_data(index_p))+((1.0_ccs_real-interpolation_factor)*viscosity_data(index_nb))
            r1=face_area*viscosity_face*dot_product(duvw,face_normal)
          else ! boundary face
            r1=face_area*viscosity_data(index_p)*dot_product(duvwp,face_normal)
          end if
          r2=r1+r2
        else if (component == 2) then ! y-component of velocity
          ! present cell gradients
          duvwp(1)=duy_data(index_p)
          duvwp(2)=dvy_data(index_p)
          duvwp(3)=dwy_data(index_p)
          ! neighbouring cell gradients
          duvwf(1)=duy_data(index_nb)
          duvwf(2)=dvy_data(index_nb)
          duvwf(3)=dwy_data(index_nb)
          
          if(.not.is_boundary) then ! no boundary face
            duvw(1)=(interpolation_factor*duvwp(1))+((1.0_ccs_real-interpolation_factor)*duvwf(1))
            duvw(2)=(interpolation_factor*duvwp(2))+((1.0_ccs_real-interpolation_factor)*duvwf(2))
            duvw(3)=(interpolation_factor*duvwp(3))+((1.0_ccs_real-interpolation_factor)*duvwf(3))
            viscosity_face=(interpolation_factor*viscosity_data(index_p))+((1.0_ccs_real-interpolation_factor)*viscosity_data(index_nb))
            r1=face_area*viscosity_face*dot_product(duvw,face_normal)
          else ! boundary face
            r1=face_area*viscosity_data(index_p)*dot_product(duvwp,face_normal)
          end if
          r2=r1+r2
        else if(component == 3) then ! z-component of velocity
          ! present cell gradients
          duvwp(1)=duz_data(index_p)
          duvwp(2)=dvz_data(index_p)
          duvwp(3)=dwz_data(index_p)
          ! neighbouring cell gradients
          duvwf(1)=duz_data(index_nb)
          duvwf(2)=dvz_data(index_nb)
          duvwf(3)=dwz_data(index_nb)
          
          if(.not.is_boundary) then ! no boundary face
            duvw(1)=(interpolation_factor*duvwp(1))+((1.0_ccs_real-interpolation_factor)*duvwf(1))
            duvw(2)=(interpolation_factor*duvwp(2))+((1.0_ccs_real-interpolation_factor)*duvwf(2))
            duvw(3)=(interpolation_factor*duvwp(3))+((1.0_ccs_real-interpolation_factor)*duvwf(3))
            viscosity_face=(interpolation_factor*viscosity_data(index_p))+((1.0_ccs_real-interpolation_factor)*viscosity_data(index_nb))
            r1=face_area*viscosity_face*dot_product(duvw,face_normal)
          else ! boundary face
            r1=face_area*viscosity_data(index_p)*dot_product(duvwp,face_normal)
          end if
          r2=r1+r2
        end if 
      end do

      r3=r2*Vol
      call set_row(global_index_p, vec_values)
      call set_entry(r3, vec_values)
      call set_values(vec_values, vec)
    end do

    deallocate (vec_values%global_indices)
    deallocate (vec_values%values)

    ! Restoring the necessary data    
    call restore_vector_data(viscosity%values, viscosity_data)
    ! x-component velocity gradient values
    call restore_vector_data(u%x_gradients, dux_data)
    call restore_vector_data(v%x_gradients, dvx_data)
    call restore_vector_data(w%x_gradients, dwx_data)
    ! y-component velocity gradient values
    call restore_vector_data(u%y_gradients, duy_data)
    call restore_vector_data(v%y_gradients, dvy_data)
    call restore_vector_data(w%y_gradients, dwy_data)
    ! z-component velocity gradient values
    call restore_vector_data(u%z_gradients, duz_data)
    call restore_vector_data(v%z_gradients, dvz_data)
    call restore_vector_data(w%z_gradients, dwz_data)
    
  end subroutine calculate_momentum_viscous_source



  !v Solves the pressure correction equation
  !
  !  Solves the pressure correction equation formed by the mass-imbalance.
  subroutine calculate_pressure_correction(par_env, mesh, invAu, invAv, invAw, M, vec, lin_sys, p_prime, lin_solver)

    use fv, only: compute_boundary_coeffs

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env !< the parallel environment
    class(ccs_mesh), intent(in) :: mesh                             !< the mesh
    class(ccs_vector), intent(inout) :: invAu, invAv, invAw            !< inverse diagonal momentum coefficients
    class(ccs_matrix), allocatable, intent(inout) :: M              !< matrix object
    class(ccs_vector), allocatable, intent(inout) :: vec            !< the RHS vector
    type(equation_system), intent(inout) :: lin_sys                 !< linear system object
    class(field), intent(inout) :: p_prime                          !< the pressure correction field
    class(linear_solver), allocatable, intent(inout) :: lin_solver  !< Linear solver that is being reused

    ! Local variables
    type(matrix_values) :: mat_coeffs
    type(vector_values) :: vec_values
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: global_index_p, global_index_nb, index_p
    integer(ccs_int) :: j
    integer(ccs_int) :: nnb
    integer(ccs_int) :: row, col
    real(ccs_real) :: face_area
    real(ccs_real), dimension(ndim) :: face_normal
    real(ccs_real) :: r
    real(ccs_real) :: coeff_f, coeff_p, coeff_nb
    real(ccs_real) :: aPb, bP
    logical :: is_boundary

    real(ccs_real), dimension(:), pointer :: invAu_data
    real(ccs_real), dimension(:), pointer :: invAv_data
    real(ccs_real), dimension(:), pointer :: invAw_data

    real(ccs_real) :: Vp
    real(ccs_real) :: V_nb
    real(ccs_real) :: Vf
    real(ccs_real) :: invA_p
    real(ccs_real) :: invA_nb
    real(ccs_real) :: invA_f

    integer(ccs_int) :: index_nb

    integer(ccs_int) :: cps   ! Cells per side
    integer(ccs_int) :: rcrit ! Global index of approximate central cell

    ! Specify block size (how many elements to set at once?)
    integer(ccs_int) :: block_nrows
    integer(ccs_int) :: block_ncols

    type(matrix_values_spec) :: mat_val_spec

    real(ccs_real) :: uSwitch, vSwitch, wSwitch
    real(ccs_real) :: problem_dim
    real(ccs_real) :: interpol_factor

    real(ccs_real), dimension(ndim) :: dx
    real(ccs_real) :: dxmag

    integer(ccs_int) :: global_num_cells

    ! First zero matrix
    call zero(M)

    ! The computed mass imbalance is +ve, to have a +ve diagonal coefficient we need to negate this.
    call dprint("P': negate RHS")
    call scale_vec(-1.0_ccs_real, vec)
    call update(vec)

    call create_vector_values(1_ccs_int, vec_values)
    call set_mode(add_mode, vec_values)

    call update(M)

    call dprint("P': get invA")
    call get_vector_data(invAu, invAu_data)
    call get_vector_data(invAv, invAv_data)
    call get_vector_data(invAw, invAw_data)

    uSwitch = 1.0_ccs_real
    vSwitch = 1.0_ccs_real
    wSwitch = 1.0_ccs_real
    problem_dim = uSwitch + vSwitch + wSwitch

    ! Loop over cells
    call dprint("P': cell loop")
    call get_local_num_cells(mesh, local_num_cells)
    do index_p = 1, local_num_cells
      call clear_entries(vec_values)

      call create_cell_locator(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)

      allocate (mat_coeffs%global_row_indices(1))
      allocate (mat_coeffs%global_col_indices(1 + nnb))
      allocate (mat_coeffs%values(1 + nnb))

      block_nrows = 1_ccs_int
      block_ncols = 1_ccs_int + nnb
      call set_matrix_values_spec_nrows(block_nrows, mat_val_spec)
      call set_matrix_values_spec_ncols(block_ncols, mat_val_spec)
      call create_matrix_values(mat_val_spec, mat_coeffs)
      call set_mode(insert_mode, mat_coeffs)

      row = global_index_p
      coeff_p = 0.0_ccs_real
      r = 0.0_ccs_real

      call get_volume(loc_p, Vp)
      invA_p = (uSwitch * invAu_data(index_p) + vSwitch * invAv_data(index_p) + wSwitch * invAw_data(index_p)) / problem_dim

      ! Loop over faces
      do j = 1, nnb
        call create_face_locator(mesh, index_p, j, loc_f)
        call get_face_area(loc_f, face_area)
        call get_face_normal(loc_f, face_normal)

        call get_boundary_status(loc_f, is_boundary)

        if (.not. is_boundary) then
          ! Interior face
          call create_neighbour_locator(loc_p, j, loc_nb)
          call get_global_index(loc_nb, global_index_nb)
          call get_local_index(loc_nb, index_nb)
          call get_face_interpolation(loc_f, interpol_factor)

          call get_distance(loc_p, loc_nb, dx)
          dxmag = sqrt(sum(dx**2))
          coeff_f = (1.0 / dxmag) * face_area

          call get_volume(loc_nb, V_nb)
          Vf = interpol_factor * Vp + (1.0_ccs_real - interpol_factor) * V_nb

          invA_nb = (uSwitch * invAu_data(index_nb) + vSwitch * invAv_data(index_nb) + wSwitch * invAw_data(index_nb)) / problem_dim
          invA_f = interpol_factor * invA_p + (1.0_ccs_real - interpol_factor) * invA_nb

          coeff_f = -(Vf * invA_f) * coeff_f

          coeff_nb = coeff_f
          col = global_index_nb
        else
          call get_distance(loc_p, loc_f, dx)
          dxmag = sqrt(sum(dx**2))

          coeff_f = (1.0 / (2 * dxmag)) * face_area
          coeff_f = -(Vp * invA_p) * coeff_f

          call compute_boundary_coeffs(p_prime, 0, loc_p, loc_f, face_normal, aPb, bP)
          coeff_p = coeff_p + coeff_f * aPb
          r = r - coeff_f * bP
          col = -1 ! Don't attempt to set neighbour coefficients
          coeff_nb = 0.0
        end if
        coeff_p = coeff_p - coeff_f

        call set_row(row, mat_coeffs)
        call set_col(col, mat_coeffs)
        call set_entry(coeff_nb, mat_coeffs)
        ! call clear_entries(mat_coeffs)

      end do

      ! XXX: Need to fix pressure somewhere
      !      Row is the global index - should be unique
      !      Locate approximate centre of mesh (assuming a square)
      if (.not. any(p_prime%bcs%bc_types(:) == bc_type_dirichlet)) then
        call get_global_num_cells(mesh, global_num_cells)
        cps = int(sqrt(real(global_num_cells)), ccs_int)
        rcrit = (cps / 2) * (1 + cps)
        if (row == rcrit) then
          coeff_p = coeff_p + 1.0e30 ! Force diagonal to be huge -> zero solution (approximately).
          call dprint("Fixed coeff_p" // str(coeff_p) // " at " // str(row))
        end if
      end if

      ! Add the diagonal entry
      col = row
      call set_row(row, mat_coeffs)
      call set_col(col, mat_coeffs)
      call set_entry(coeff_p, mat_coeffs)

      call set_row(global_index_p, vec_values)
      call set_entry(r, vec_values)

      ! Set the values
      call set_values(mat_coeffs, M)
      call set_values(vec_values, vec)
      call clear_entries(mat_coeffs)

      deallocate (mat_coeffs%global_row_indices)
      deallocate (mat_coeffs%global_col_indices)
      deallocate (mat_coeffs%values)
    end do

    call dprint("P': restore invA")
    call restore_vector_data(invAu, invAu_data)
    call restore_vector_data(invAv, invAv_data)
    call restore_vector_data(invAw, invAw_data)

    ! Assembly of coefficient matrix and source vector
    call dprint("P': assemble matrix, RHS")
    call update(M)
    call update(vec)
    call finalise(M)

    ! Create linear solver
    call dprint("P': create lin sys")
    if (allocated(p_prime%values%name)) then
      call set_equation_system(par_env, vec, p_prime%values, M, lin_sys, p_prime%values%name)
    else
      call set_equation_system(par_env, vec, p_prime%values, M, lin_sys)
    end if
    call create_solver(lin_sys, lin_solver)

    ! Customise linear solver
    call set_solver_method(pressure_solver_method_name, lin_solver)
    call set_solver_precon(pressure_solver_precon_name, lin_solver)

    ! Solve the linear system
    call dprint("P': solve")
    call solve(lin_solver)

  end subroutine calculate_pressure_correction

  !> Computes the per-cell mass imbalance, updating the face velocity flux as it does so.
  subroutine compute_mass_imbalance(mesh, invAu, invAv, invAw, ivar, flow, input_b, residuals)

    type(ccs_mesh), intent(in) :: mesh      !< The mesh object
    class(ccs_vector), intent(inout) :: invAu  !< The inverse x momentum equation diagonal coefficient
    class(ccs_vector), intent(inout) :: invAv  !< The inverse y momentum equation diagonal coefficient
    class(ccs_vector), intent(inout) :: invAw  !< The inverse z momentum equation diagonal coefficient
    integer(ccs_int), intent(inout) :: ivar !< Counter for flow variables
    type(fluid), intent(inout) :: flow
    class(ccs_vector), target, intent(inout) :: input_b   !< The per-cell mass imbalance
    real(ccs_real), dimension(:), intent(inout) :: residuals !< Residual for each equation

    class(ccs_vector), pointer :: b   !< The per-cell mass imbalance
    type(vector_values) :: vec_values
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i   ! Cell counter
    integer(ccs_int) :: j   ! Cell-face counter

    type(cell_locator) :: loc_p   ! Central cell locator object
    type(face_locator) :: loc_f   ! Face locator object

    integer(ccs_int) :: global_index_p  ! Central cell global index
    real(ccs_real) :: face_area         ! Face area
    integer(ccs_int) :: index_f         ! Face index
    integer(ccs_int) :: nnb             ! Cell neighbour count

    real(ccs_real), dimension(:), pointer :: mf_data      ! Data array for the mass flux
    real(ccs_real), dimension(:), pointer :: p_data       ! Data array for pressure
    real(ccs_real), dimension(:), pointer :: dpdx_data    ! Data array for pressure x gradient
    real(ccs_real), dimension(:), pointer :: dpdy_data    ! Data array for pressure y gradient
    real(ccs_real), dimension(:), pointer :: dpdz_data    ! Data array for pressure z gradient
    real(ccs_real), dimension(:), pointer :: invAu_data   ! Data array for inverse x momentum
    ! diagonal coefficient
    real(ccs_real), dimension(:), pointer :: invAv_data   ! Data array for inverse y momentum
    ! diagonal coefficient
    real(ccs_real), dimension(:), pointer :: invAw_data ! Data array for inverse z momentum
    ! diagonal coefficient

    logical :: is_boundary            ! Boundary indicator
    type(neighbour_locator) :: loc_nb ! Neighbour cell locator object
    integer(ccs_int) :: index_nb      ! Neighbour cell index

    real(ccs_real) :: mib ! Cell mass imbalance
    integer(ccs_int) :: nvar ! Number of flow variables to solve
    integer(ccs_int) :: global_num_cells

    logical, save :: first_time = .true.

    class(field), pointer :: u        !< The x velocity component
    class(field), pointer :: v        !< The y velocity component
    class(field), pointer :: w        !< The z velocity component
    class(field), pointer :: p        !< The pressure field
    class(field), pointer :: mf       !< The face velocity flux

    call get_field(flow, field_u, u)
    call get_field(flow, field_v, v)
    call get_field(flow, field_w, w)
    call get_field(flow, field_p, p)
    call get_field(flow, field_mf, mf)

    ! Set variable index for pressure
    if (first_time) then
      ivar = ivar + 1
      varp = ivar
      first_time = .false.
    end if

    call create_vector_values(1_ccs_int, vec_values)
    call set_mode(insert_mode, vec_values)

    if (allocated(p%residuals)) then
      b => p%residuals
    else
      b => input_b
    end if

    ! First zero RHS
    call zero(b)

    ! Update vectors to make sure all data is up to date
    call update(u%values)
    call update(v%values)
    call update(w%values)
    call update(p%values)
    call update(p%x_gradients)
    call update(p%y_gradients)
    call update(p%z_gradients)
    call get_vector_data(mf%values, mf_data)
    call get_vector_data(p%values, p_data)
    call get_vector_data(p%x_gradients, dpdx_data)
    call get_vector_data(p%y_gradients, dpdy_data)
    call get_vector_data(p%z_gradients, dpdz_data)
    call get_vector_data(invAu, invAu_data)
    call get_vector_data(invAv, invAv_data)
    call get_vector_data(invAw, invAw_data)

    call get_local_num_cells(mesh, local_num_cells)
    do i = 1, local_num_cells
      call clear_entries(vec_values)

      call create_cell_locator(mesh, i, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)

      mib = 0.0_ccs_real

      do j = 1, nnb
        call create_face_locator(mesh, i, j, loc_f)
        call get_face_area(loc_f, face_area)
        call get_local_index(loc_f, index_f)

        ! Check face orientation
        call get_boundary_status(loc_f, is_boundary)
        if (.not. is_boundary) then
          call create_neighbour_locator(loc_p, j, loc_nb)
          call get_local_index(loc_nb, index_nb)
          if (index_nb < i) then
            face_area = -face_area
          else
            ! Compute mass flux through face
            mf_data(index_f) = calc_mass_flux(u, v, w, &
                                              p_data, dpdx_data, dpdy_data, dpdz_data, &
                                              invAu_data, invAv_data, invAw_data, &
                                              loc_f, p%enable_cell_corrections)
          end if
        else
          ! Compute mass flux through face
          mf_data(index_f) = calc_mass_flux(u, v, w, &
                                            p_data, dpdx_data, dpdy_data, dpdz_data, &
                                            invAu_data, invAv_data, invAw_data, &
                                            loc_f, .false.)
        end if

        mib = mib + mf_data(index_f) * face_area
      end do

      call set_row(global_index_p, vec_values)
      call set_entry(mib, vec_values)
      call set_values(vec_values, b)
    end do

    call restore_vector_data(mf%values, mf_data)
    call restore_vector_data(p%values, p_data)
    call restore_vector_data(p%x_gradients, dpdx_data)
    call restore_vector_data(p%y_gradients, dpdy_data)
    call restore_vector_data(p%z_gradients, dpdz_data)
    call restore_vector_data(invAu, invAu_data)
    call restore_vector_data(invAv, invAv_data)
    call restore_vector_data(invAw, invAw_data)
    ! Update vectors on exit (just in case)
    call update(u%values)
    call update(v%values)
    call update(w%values)
    call update(p%values)
    call update(p%x_gradients)
    call update(p%y_gradients)
    call update(p%z_gradients)

    call update(b)
    call update(mf%values)

    ! Pressure residual
    ! Stores RMS of residuals
    call get_global_num_cells(mesh, global_num_cells)
    residuals(varp) = norm(b, 2) / sqrt(real(global_num_cells))
    ! Stores Linf norm of residuals
    nvar = int(size(residuals) / 2_ccs_int)
    residuals(varp + nvar) = norm(b, 0)

  end subroutine compute_mass_imbalance

  !> Corrects the pressure field, using explicit underrelaxation
  subroutine update_pressure(mesh, p_prime, p)

    use case_config, only: pressure_relax

    ! Arguments
    class(ccs_mesh), intent(in) :: mesh    !< The mesh
    class(field), intent(in) :: p_prime !< pressure correction
    class(field), intent(inout) :: p    !< the pressure field being corrected

    call axpy(pressure_relax, p_prime%values, p%values)

    call update_gradient(mesh, p)

  end subroutine update_pressure

  !> Corrects the velocity field using the pressure correction gradient
  subroutine update_velocity(mesh, invAu, invAv, invAw, flow)

    use vec, only: zero_vector

    ! Arguments
    class(ccs_mesh), intent(in) :: mesh    !< The mesh
    class(ccs_vector), intent(in) :: invAu !< The inverse x momentum equation diagonal coefficient
    class(ccs_vector), intent(in) :: invAv !< The inverse y momentum equation diagonal coefficient
    class(ccs_vector), intent(in) :: invAw !< The inverse z momentum equation diagonal coefficient
    type(fluid), intent(inout) :: flow

    class(field), pointer :: p_prime !< The pressure correction
    class(field), pointer :: u       !< The x velocities being corrected
    class(field), pointer :: v       !< The y velocities being corrected
    class(field), pointer :: w       !< The z velocities being corrected

    call get_field(flow, field_u, u)
    call get_field(flow, field_v, v)
    call get_field(flow, field_w, w)
    call get_field(flow, field_p_prime, p_prime)

    ! First update gradients
    call zero_vector(p_prime%x_gradients)
    call zero_vector(p_prime%y_gradients)
    call zero_vector(p_prime%z_gradients)
    call update_gradient(mesh, p_prime)

    ! Multiply gradients by inverse diagonal coefficients
    call mult(invAu, p_prime%x_gradients)
    call mult(invAv, p_prime%y_gradients)
    call mult(invAw, p_prime%z_gradients)

    ! Compute correction source on velocity
    call calculate_momentum_pressure_source(mesh, p_prime%x_gradients, u%values)
    call calculate_momentum_pressure_source(mesh, p_prime%y_gradients, v%values)
    call calculate_momentum_pressure_source(mesh, p_prime%z_gradients, w%values)

    call update(u%values)
    call update(v%values)
    call update(w%values)

    call update_gradient(mesh, u)
    call update_gradient(mesh, v)
    call update_gradient(mesh, w)

  end subroutine update_velocity

  !> Corrects the face velocity flux using the pressure correction
  subroutine update_face_velocity(mesh, invAu, invAv, invAw, p_prime, mf, b, residuals)

    type(ccs_mesh), intent(in) :: mesh                               !< The mesh
    class(ccs_vector), intent(inout) :: invAu                        !< The inverse x momentum equation diagonal coefficient
    class(ccs_vector), intent(inout) :: invAv                        !< The inverse y momentum equation diagonal coefficient
    class(ccs_vector), intent(inout) :: invAw                        !< The inverse z momentum equation diagonal coefficient
    class(field), intent(inout) :: p_prime                           !< The pressure correction
    class(field), intent(inout) :: mf                                !< The face velocity being corrected
    class(ccs_vector), intent(inout) :: b   !< The per-cell mass imbalance
    real(ccs_real), dimension(:), intent(inout) :: residuals !< Residual for each equation

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i

    real(ccs_real) :: mf_prime
    real(ccs_real), dimension(:), allocatable :: zero_arr
    real(ccs_real), dimension(:), pointer :: mf_data
    real(ccs_real), dimension(:), pointer :: pp_data
    real(ccs_real), dimension(:), pointer :: invAu_data
    real(ccs_real), dimension(:), pointer :: invAv_data
    real(ccs_real), dimension(:), pointer :: invAw_data

    type(cell_locator) :: loc_p
    integer(ccs_int) :: nnb
    integer(ccs_int) :: j
    type(face_locator) :: loc_f
    integer(ccs_int) :: index_f

    logical :: is_boundary
    type(neighbour_locator) :: loc_nb
    integer(ccs_int) :: index_nb

    integer(ccs_int) :: global_index_p  ! Central cell global index
    real(ccs_real) :: face_area         ! Face area
    real(ccs_real) :: mib
    integer(ccs_int) :: nvar ! Number of flow variables to solve
    type(vector_values) :: vec_values

    integer(ccs_int) :: global_num_cells

    call create_vector_values(1_ccs_int, vec_values)
    call set_mode(insert_mode, vec_values)
    call zero(b)

    ! Update vector to make sure data is up to date
    call update(p_prime%values)
    call get_vector_data(p_prime%values, pp_data)
    call get_vector_data(invAu, invAu_data)
    call get_vector_data(invAv, invAv_data)
    call get_vector_data(invAw, invAw_data)
    call get_vector_data(mf%values, mf_data)

    allocate (zero_arr(size(pp_data)))
    zero_arr(:) = 0.0_ccs_real

    ! XXX: This should really be a face loop
    call get_local_num_cells(mesh, local_num_cells)
    do i = 1, local_num_cells
      call clear_entries(vec_values)
      mib = 0.0_ccs_real

      call create_cell_locator(mesh, i, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb
        call create_face_locator(mesh, i, j, loc_f)
        call get_local_index(loc_f, index_f)
        call get_face_area(loc_f, face_area)
        call get_boundary_status(loc_f, is_boundary)
        if (.not. is_boundary) then
          call create_neighbour_locator(loc_p, j, loc_nb)
          call get_local_index(loc_nb, index_nb)
          if (i < index_nb) then
            mf_prime = calc_mass_flux(pp_data, zero_arr, zero_arr, zero_arr, &
                                      invAu_data, invAv_data, invAw_data, loc_f, p_prime%enable_cell_corrections)

            mf_data(index_f) = mf_data(index_f) + mf_prime
          else
            face_area = -face_area
          end if
        end if

        mib = mib + mf_data(index_f) * face_area
      end do

      call set_row(global_index_p, vec_values)
      call set_entry(mib, vec_values)
      call set_values(vec_values, b)
    end do

    deallocate (zero_arr)

    call restore_vector_data(p_prime%values, pp_data)
    call restore_vector_data(invAu, invAu_data)
    call restore_vector_data(invAv, invAv_data)
    call restore_vector_data(invAw, invAw_data)
    call restore_vector_data(mf%values, mf_data)

    call update(mf%values)
    ! Update vector on exit (just in case)
    call update(p_prime%values)

    !! Get corrected mass-imbalance
    call update(b)

    ! Stores RMS of residuals
    call get_global_num_cells(mesh, global_num_cells)
    residuals(varp + 1) = norm(b, 2) / sqrt(real(global_num_cells))
    ! Stores Linf norm of residuals
    nvar = int(size(residuals) / 2_ccs_int)
    residuals(varp + 1 + nvar) = norm(b, 0)

  end subroutine update_face_velocity

  subroutine check_convergence(par_env, itr, residuals, res_target, &
                               flow_solver_selector, converged)

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env !< The parallel environment
    integer(ccs_int), intent(in) :: itr                             !< Iteration count
    real(ccs_real), dimension(:), intent(in) :: residuals           !< RMS and Linf of residuals for each equation
    real(ccs_real), intent(in) :: res_target                        !< Target residual
    type(fluid_solver_selector), intent(in) :: flow_solver_selector
    logical, intent(inout) :: converged                             !< Has solution converged (true/false)

    ! Local variables
    integer :: io_unit
    integer(ccs_int) :: step                            !< The current time-step
    real(ccs_real) :: time                              !< The current time
    integer(ccs_int) :: nvar              ! Number of variables (u,v,w,p,etc)
    integer(ccs_int) :: i
    character(len=60) :: fmt              ! Format string for writing out residuals
    character(len=60) :: prefix           ! prefix for residual norms
    logical, save :: first_time = .true.  ! Whether first time this subroutine is called

    logical :: u_sol                                    !< Is x-velocity being solved (true/false)
    logical :: v_sol                                    !< Is y-velocity being solved (true/false)
    logical :: w_sol                                    !< Is z-velocity being solved (true/false)
    logical :: p_sol                                    !< Is pressure field being solved (true/false)

    call get_fluid_solver_selector(flow_solver_selector, field_u, u_sol)
    call get_fluid_solver_selector(flow_solver_selector, field_v, v_sol)
    call get_fluid_solver_selector(flow_solver_selector, field_w, w_sol)
    call get_fluid_solver_selector(flow_solver_selector, field_p, p_sol)
    call get_current_step(step)
    call get_current_time(time)

    nvar = int(size(residuals) / 2_ccs_int)

    ! Print residuals
    if (par_env%proc_id == par_env%root) then
      if (first_time) then
        ! Write header
        open (newunit=io_unit, file="residuals.log", status="replace", form="formatted")

        write (*, *)
        if (step >= 0) then
          write (*, '(a6, 1x, a6)', advance='no') 'Step', 'Iter'

          write (io_unit, '(a6, 1x, a12, 1x, a6)', advance='no') '#step', 'time', 'iter'
        else
          write (*, '(a6)', advance='no') 'Iter'

          write (io_unit, '(a6)', advance='no') '#iter'
        end if
        do i = 1, 2
          if (u_sol) write (*, '(1x,a12)', advance='no') 'u'
          if (v_sol) write (*, '(1x,a12)', advance='no') 'v'
          if (w_sol) write (*, '(1x,a12)', advance='no') 'w'
          if (p_sol) write (*, '(1x,a12)', advance='no') 'p'
          if (p_sol) write (*, '(1x,a12)', advance='no') '|div(u)|'

          if (i == 1) then
            prefix = "L2_"
          else
            prefix = "Linf_"
          end if
          if (u_sol) write (io_unit, '(1x,a12)', advance='no') trim(prefix) // 'u'
          if (v_sol) write (io_unit, '(1x,a12)', advance='no') trim(prefix) // 'v'
          if (w_sol) write (io_unit, '(1x,a12)', advance='no') trim(prefix) // 'w'
          if (p_sol) write (io_unit, '(1x,a12)', advance='no') trim(prefix) // 'p'
          if (p_sol) write (io_unit, '(1x,a12)', advance='no') trim(prefix) // 'div(u)'
        end do
        write (*, *)
        write (io_unit, *)
        first_time = .false.
      else
        open (newunit=io_unit, file="residuals.log", status="old", form="formatted", position="append")
      end if

      ! Write step, iteration and residuals
      if (step >= 0) then
        fmt = '(i6,1x,i6,' // str(2 * nvar) // '(1x,e12.4))'
        write (*, fmt) step, itr, residuals(1:2 * nvar)

        fmt = '(i6,1x,e12.4,1x,i6,' // str(2 * nvar) // '(1x,e12.4))'
        write (io_unit, fmt) step, time, itr, residuals(1:2 * nvar)
      else
        fmt = '(i6,' // str(2 * nvar) // '(1x,e12.4))'
        write (*, fmt) itr, residuals(1:2 * nvar)

        write (io_unit, fmt) itr, residuals(1:2 * nvar)
      end if
      close (io_unit)
    end if

    ! checks if RMS of residuals is below target
    if (maxval(residuals(1:nvar)) < res_target) converged = .true.

    if (maxval(residuals) > huge(1.0_ccs_real)) then
      call error_abort("Computation diverging")
    end if

  end subroutine check_convergence

  !v Applies implicit underrelaxation to an equation
  !
  !  Extracts the diagonal coefficient of a matrix and divides by the URF, adding a
  !  proportional explicit term to the RHS vector.
  subroutine underrelax(mesh, alpha, phi, diag, M, b)

    use mat, only: set_matrix_diagonal

    type(ccs_mesh), intent(in) :: mesh
    real(ccs_real), intent(in) :: alpha
    class(field), intent(inout) :: phi
    class(ccs_vector), intent(inout) :: diag
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b

    real(ccs_real), dimension(:), pointer :: diag_data
    real(ccs_real), dimension(:), pointer :: phi_data
    real(ccs_real), dimension(:), pointer :: b_data

    integer(ccs_int) :: local_num_cells
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
    call get_local_num_cells(mesh, local_num_cells)
    do i = 1, local_num_cells
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
