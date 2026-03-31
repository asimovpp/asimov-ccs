!> Program file for Poiseuille case
module poiseuille_core
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use core
  use ccs_base, only: mesh
  use constants, only: ccs_string_len, ndim
  use kinds, only: ccs_real, ccs_int, ccs_long
  use kinds, only: CCS_MPI_PRECISION
  use types, only: field, fluid, ccs_mesh, bc_profile
  use parallel, only: timer, is_root
  use parallel_types, only: parallel_environment
  use utils, only:  calc_kinetic_energy, calc_enstrophy
  use fields, only: get_field, dealloc_fluid_fields
  use boundary_conditions, only: set_bc_profile
  use timestepping, only: reset_timestepping
  use meshing, only: get_total_num_cells, set_mesh_object, get_global_num_cells, nullify_mesh_object, get_local_num_cells
  use io_visualisation, only: reset_io_visualisation
  use utils, only: str, exit_print, reset_outputlist_counter
  use timers, only: timer_init, timer_register_start, timer_register, timer_start, timer_stop, &
                    timer_print, timer_get_time, timer_print_all, timer_reset
  use logging, only: log_unit_out

  implicit none

  public :: run_poiseuille

  integer(ccs_int), dimension(:), allocatable :: variable_types              ! cell centred upwind, central, etc.

  ! Global variables to pass errors to/from postprocessing
  real(ccs_real), dimension(3) :: pois_error_L2_global, pois_error_Linf_global
  
  contains

  subroutine run_poiseuille(par_env, shared_env, error_L2, error_Linf, input_mesh)
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
    real(ccs_real), dimension(3), intent(out) :: error_L2 !< L2 norm of the error for the U, V and P fields respectively
    real(ccs_real), dimension(3), intent(out) :: error_Linf !< Linf norm of the error for the U, V and P fields respectively
    type(ccs_mesh), intent(inout), optional :: input_mesh !< mesh object to use, if not provided, the build_square_mesh is used

    class(field), pointer :: u
    type(bc_profile), allocatable :: profile

    integer(ccs_int) :: timer_index_total
    integer(ccs_int) :: timer_index_init

    type(fluid) :: flow_fields

    type(ccs_options) :: run_options

    call timer_register_start("Elapsed time", timer_index_total)
    call timer_register_start("Init time", timer_index_init)

    ! Read case name and runtime parameters from configuration file
    call get_config(par_env, run_options)

    if (is_root(par_env)) write(log_unit_out,*) "Starting ", run_options%paths%case_name, " case!"

    if (present(input_mesh)) then
      mesh = input_mesh
      call set_mesh_object(mesh)
    else
      call initialise_mesh(par_env, shared_env, run_options)
    end if

    ! Initialise fields
    if (is_root(par_env)) write(log_unit_out,*) "Initialise fields"
    call initialise_fields(par_env, run_options, flow_fields)

    ! Create and initialise field vectors
    if (is_root(par_env)) write(log_unit_out,*) "Initialise field vectors"
    
    
    ! Set to 1st boundary condition (inlet)
    call get_field(flow_fields, "u", u)
    call get_inlet_profile(run_options, profile)
    call set_bc_profile(u, profile, 1)
    nullify(u)

    ! Initialise velocity field
    if (is_root(par_env)) write(log_unit_out,*) "Initialise velocity field"
    call initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)

    ! Solve using SIMPLE algorithm
    if (is_root(par_env)) write(log_unit_out,*) "Start SIMPLE"

    call timer_stop(timer_index_init)

    pois_error_L2_global = 0.0_ccs_real
    pois_error_Linf_global = 0.0_ccs_real
    call run_solver(par_env, run_options, eval_sources, postproc_poiseuille, flow_fields)
    error_L2 = pois_error_L2_global
    error_Linf = pois_error_Linf_global
    
    ! Clean-up

    call timer_stop(timer_index_total)

    call timer_print_all(par_env)

    call reset_timestepping()
    call reset_outputlist_counter()
    call reset_io_visualisation()
    call timer_reset()
    call dealloc_fluid_fields(flow_fields)
    call nullify_mesh_object()

  end subroutine

  subroutine postproc_poiseuille(par_env, flow_fields)

    class(parallel_environment), allocatable, intent(in) :: par_env
    type(fluid), intent(in) :: flow_fields
    class(field), pointer :: u, v, w, p

    call get_field(flow_fields, "u", u)
    call get_field(flow_fields, "v", v)
    call get_field(flow_fields, "w", w)
    call get_field(flow_fields, "p", p)
    call calc_kinetic_energy(par_env, u, v, w)
    call calc_enstrophy(par_env, u, v, w)

    call calc_error(par_env, u, v, p, pois_error_L2_global, pois_error_Linf_global)
    nullify(u)
    nullify(v)
    nullify(w)
    nullify(p)
    
  end subroutine postproc_poiseuille

  pure subroutine get_init_flow(loc_p, field_name, init_val)
    use types, only: cell_locator
    use meshing, only: get_centre

    type(cell_locator), intent(in) :: loc_p
    character(len=*), intent(in) :: field_name
    real(ccs_real), intent(inout) :: init_val

    ! Silence compiler warnings
    associate(foo => loc_p, bar => field_name, baz => init_val)
    end associate

  end subroutine
  
  pure subroutine get_init_mass_flux(loc_f, init_val)
    use types, only: face_locator
    type(face_locator), intent(in) :: loc_f
    real(ccs_real), intent(inout) :: init_val

    associate (foo => loc_f, bar => init_val)
    end associate

  end subroutine

  subroutine get_inlet_profile(run_options, profile)

    type(ccs_options), intent(in) :: run_options
    type(bc_profile), allocatable, intent(out) :: profile
    integer(ccs_int) :: n, i
    real(ccs_real) :: y, h, mu, P

    integer(ccs_int) :: cps

    cps = run_options%mesh%cps
    n = 200*3*cps

    allocate(profile)

    allocate(profile%centre(3))
    allocate(profile%values(n))
    allocate(profile%coordinates(n))
    h = 1.0_ccs_real
    mu = 0.01_ccs_real
    P = 8*mu 

    profile%centre(:) = [ 0, 0, 0 ]

    do i=1, n
      y =  real(i-1, ccs_real)*h/real(n-1, ccs_real)
      profile%coordinates(i) = y      
      profile%values(i) = P*y*(h-y)/ (2.0_ccs_real*mu)
    end do

  end subroutine
  
  subroutine calc_error(par_env, u, v, p, error_L2, error_Linf)

    use constants, only: ndim
    use types, only: cell_locator
    use utils, only: str

    use vec, only: get_vector_data, restore_vector_data

    use meshing, only: get_centre, create_cell_locator, get_local_num_cells

    use parallel, only: allreduce
    use parallel_types_mpi, only: parallel_environment_mpi
    use timestepping, only: get_current_step

    class(parallel_environment), intent(in) :: par_env !< The parallel environment
    class(field), intent(inout) :: u, v, p
    real(ccs_real), dimension(3), intent(out) :: error_L2
    real(ccs_real), dimension(3), intent(out) :: error_Linf

    real(ccs_real), dimension(3) :: error_L2_local
    real(ccs_real), dimension(3) :: error_Linf_local

    real(ccs_real) :: u_an, v_an, p_an
    real(ccs_real), dimension(:), pointer :: u_data, v_data, p_data

    real(ccs_real) :: mu, rho, nu, x, y

    logical, save :: first_time = .true.

    type(cell_locator) :: loc_p
    real(ccs_real), dimension(ndim) :: x_p
    integer(ccs_int) :: index_p, local_num_cells

    character(len=ccs_string_len) :: fmt
    integer(ccs_int) :: step

    integer(ccs_int) :: global_num_cells

    integer :: io_unit

    integer :: ierr

    mu = 0.01_ccs_real ! XXX: currently hardcoded somewhere
    rho = 1.0_ccs_real ! XXX: implicitly 1 throughout
    nu = mu / rho

    error_Linf_local(:) = 0.0_ccs_real
    error_L2_local(:) = 0.0_ccs_real

    call get_vector_data(u%values, u_data)
    call get_vector_data(v%values, v_data)
    call get_vector_data(p%values, p_data)
    call get_current_step(step)

    call get_local_num_cells(local_num_cells)
    do index_p = 1, local_num_cells

      call create_cell_locator(index_p, loc_p)
      call get_centre(loc_p, x_p)

      ! Compute analytical solution
      x = x_p(1)
      y = x_p(2)
      u_an = 8*mu*y*(1-y)/(2*mu)
      v_an = 0.0_ccs_real
      p_an = -8*mu*(x-1)

      error_L2_local(1) = error_L2_local(1) + (u_an - u_data(index_p))**2
      error_L2_local(2) = error_L2_local(2) + (v_an - v_data(index_p))**2
      error_L2_local(3) = error_L2_local(3) + (p_an - p_data(index_p))**2

      error_Linf_local(1) = max(error_Linf_local(1), abs(u_an - u_data(index_p)))
      error_Linf_local(2) = max(error_Linf_local(2), abs(v_an - v_data(index_p)))
      error_Linf_local(3) = max(error_Linf_local(3), abs(p_an - p_data(index_p)))

    end do
    call restore_vector_data(u%values, u_data)
    call restore_vector_data(v%values, v_data)
    call restore_vector_data(p%values, p_data)

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_AllReduce(error_L2_local, error_L2, size(error_L2), CCS_MPI_PRECISION, MPI_SUM, par_env%comm, ierr)
      call MPI_AllReduce(error_Linf_local, error_Linf, size(error_Linf), CCS_MPI_PRECISION, MPI_MAX, par_env%comm, ierr)
    class default
      call error_abort("ERROR: Unknown type")
    end select

    call get_global_num_cells(global_num_cells)
    error_L2(:) = sqrt(error_L2(:) / global_num_cells)

    if (is_root(par_env)) then
      if (first_time) then
        first_time = .false.
        open (newunit=io_unit, file="err.log", status="replace", form="formatted")
      else
        open (newunit=io_unit, file="err.log", status="old", form="formatted", position="append")
      end if
      fmt = '(I0,' // str(2 * size(error_L2)) // '(1x,e12.4))'
      write (io_unit, fmt) step, error_L2, error_Linf
      close (io_unit)
    end if

  end subroutine calc_error

  !> Case-specific source terms
  subroutine eval_sources(flow, phi, R, S)
    use types, only: fluid, field, ccs_vector
    use fv, only: zero_sources

    type(fluid), intent(in) :: flow !< Provides access to full flow field
    class(field), intent(in) :: phi !< Field being transported
    class(ccs_vector), intent(inout) :: R !< Work vector (for evaluating linear/implicit sources)
    class(ccs_vector), intent(inout) :: S !< Work vector (for evaluating fixed/explicit sources)
    
    ! Dummy implementation - just zeros the sources, see sero_sources for example implementation
    call zero_sources(flow, phi, R, S)
    
  end subroutine eval_sources

end module poiseuille_core
