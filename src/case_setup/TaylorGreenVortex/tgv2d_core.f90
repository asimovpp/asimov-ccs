module tgv2d_core
#include "ccs_macros.inc"

  use mpi
  
  use core
  use ccs_base, only: mesh
  use constants, only: ccs_string_len
  use kinds, only: ccs_real, ccs_int
  use types, only: fluid, field, ccs_mesh
  use parallel, only: is_root
  use parallel_types, only: parallel_environment
  use meshing, only: get_global_num_cells, set_mesh_object, nullify_mesh_object
  use mesh_utils, only: build_square_mesh
  use utils, only: exit_print, calc_kinetic_energy, calc_enstrophy, &
                   add_field_to_outputlist, reset_outputlist_counter
  use fields, only: get_field, add_field, dealloc_fluid_fields, set_is_field_solved
  use timestepping, only: reset_timestepping
  use io_visualisation, only: reset_io_visualisation
  use profiler, only: profiler_init, profiler_shutdown, profiler_begin_region, profiler_end_region

  implicit none

  public :: run_tgv2d

  integer(ccs_int), dimension(:), allocatable :: variable_types              ! cell centred upwind, central, etc.

  ! Global variables to pass error calculations to the postprocessing subroutine
  real(ccs_real), dimension(3) :: tgv2d_error_L2_global, tgv2d_error_Linf_global
  
contains

  subroutine run_tgv2d(par_env, shared_env, error_L2, error_Linf, input_mesh, input_dt, input_num_steps)
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    real(ccs_real), dimension(3), intent(out) :: error_L2 !< L2 norm of the error for the U, V and P fields respectively
    real(ccs_real), dimension(3), intent(out) :: error_Linf !< Linf norm of the error for the U, V and P fields respectively
    type(ccs_mesh), intent(inout), optional :: input_mesh !< mesh object to use, if not provided, the build_square_mesh is used
    real(ccs_real), intent(in), optional :: input_dt !< timestep, if not provided, the yaml config option is used
    integer(ccs_int), intent(in), optional :: input_num_steps !< number of timesteps, if not provided, the yaml config option is used

    integer(ccs_int) :: irank ! MPI rank ID
    integer(ccs_int) :: isize ! Size of MPI world

    type(fluid) :: flow_fields

    type(ccs_options) :: run_options
    
    irank = par_env%proc_id
    isize = par_env%num_procs

    call profiler_begin_region('Total elapsed time')
    call profiler_begin_region('Total initialisation')

    call get_config(par_env, run_options)
    ! call configure_parallelism(run_options, par_env, shared_env)

    ! Create a square mesh
    if (present(input_mesh)) then
      mesh = input_mesh
      call set_mesh_object(mesh)
    else
      call initialise_mesh(par_env, shared_env, run_options)
    end if

    if (present(input_dt)) then
      run_options%solve%dt = input_dt
    end if
    if (present(input_num_steps)) then
      run_options%solve%num_steps = input_num_steps
    end if

    ! Initialise fields
    if (is_root(par_env)) print *, "Initialise fields"
    call initialise_fields(par_env, run_options, flow_fields)

    ! Initialise velocity field
    if (is_root(par_env)) print *, "Initialise velocity field"
    call initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)

    ! Solve using SIMPLE algorithm
    if (is_root(par_env)) print *, "Start SIMPLE"
    
    tgv2d_error_L2_global = 0.0_ccs_real
    tgv2d_error_Linf_global = 0.0_ccs_real
    call profiler_end_region('Total initialisation')
    call run_solver(par_env, run_options, eval_sources, postproc_tgv, flow_fields)
    call profiler_end_region('Total elapsed time')
    error_L2 = tgv2d_error_L2_global
    error_Linf = tgv2d_error_Linf_global

    ! Clean-up
    call reset_timestepping()
    call reset_outputlist_counter()
    call reset_io_visualisation()
    call dealloc_fluid_fields(flow_fields)
    call nullify_mesh_object()

  end subroutine run_tgv2d

  pure subroutine get_init_flow(loc_p, field_name, init_val)

    use constants, only: ndim
    use types, only: cell_locator
    use meshing, only: get_centre
    
    type(cell_locator), intent(in) :: loc_p
    character(len=*), intent(in) :: field_name
    real(ccs_real), intent(inout) :: init_val
    real(ccs_real), dimension(ndim) :: x_p
    real(ccs_real) :: rho

    select case (field_name)
    case ("u")
      call get_centre(loc_p, x_p)
      init_val = sin(x_p(1)) * cos(x_p(2))
    case ("v")
      call get_centre(loc_p, x_p)
      init_val = -cos(x_p(1)) * sin(x_p(2))
    case ("p")
      call get_centre(loc_p, x_p)
      rho = 1.0_ccs_real ! XXX: implicitly 1 throughout
      init_val = +(rho / 4.0_ccs_real) * (cos(2 * x_p(1)) + cos(2 * x_p(2)))
    case default
      init_val = init_val + 0.0_ccs_real ! Keep the default value
    end select
    
  end subroutine get_init_flow
  
  pure subroutine get_init_mass_flux(loc_f, init_val)
    
    use constants, only: ndim
    use types, only: face_locator
    use meshing, only: get_face_normal, get_centre

    type(face_locator), intent(in) :: loc_f
    real(ccs_real), intent(inout) :: init_val
    real(ccs_real), dimension(ndim) :: x_f
    real(ccs_real), dimension(ndim) :: face_normal

    call get_face_normal(loc_f, face_normal)
    call get_centre(loc_f, x_f)
  
    init_val = sin(x_f(1)) * cos(x_f(2)) * face_normal(1) &
             - cos(x_f(1)) * sin(x_f(2)) * face_normal(2)

  end subroutine get_init_mass_flux

  subroutine calc_tgv2d_error(par_env, u, v, w, p, error_L2, error_Linf)

    use constants, only: ndim
    use types, only: cell_locator
    use utils, only: str

    use vec, only: get_vector_data, restore_vector_data

    use meshing, only: get_centre, create_cell_locator, get_local_num_cells

    use parallel, only: allreduce
    use parallel_types_mpi, only: parallel_environment_mpi
    use timestepping, only: get_current_time, get_current_step

    class(parallel_environment), intent(in) :: par_env !< The parallel environment
    class(field), intent(inout) :: u, v, w, p
    real(ccs_real), dimension(3), intent(out) :: error_L2
    real(ccs_real), dimension(3), intent(out) :: error_Linf

    real(ccs_real), dimension(3) :: error_L2_local
    real(ccs_real), dimension(3) :: error_Linf_local

    real(ccs_real) :: ft
    real(ccs_real) :: u_an, v_an, w_an, p_an
    real(ccs_real), dimension(:), pointer :: u_data, v_data, w_data, p_data

    real(ccs_real) :: mu, rho, nu

    logical, save :: first_time = .true.

    type(cell_locator) :: loc_p
    real(ccs_real), dimension(ndim) :: x_p
    integer(ccs_int) :: index_p, local_num_cells

    character(len=ccs_string_len) :: fmt
    real(ccs_real) :: time
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
    call get_vector_data(w%values, w_data)
    call get_vector_data(p%values, p_data)

    call get_current_time(time)
    call get_current_step(step)

    call get_local_num_cells(local_num_cells)
    do index_p = 1, local_num_cells

      call create_cell_locator(index_p, loc_p)
      call get_centre(loc_p, x_p)

      ! Compute analytical solution
      ft = exp(-2 * nu * time)
      ! u_an = cos(x_p(1)) * sin(x_p(2)) * ft
      ! v_an = -sin(x_p(1)) * cos(x_p(2)) * ft
      u_an = sin(x_p(1)) * cos(x_p(2)) * ft
      v_an = -cos(x_p(1)) * sin(x_p(2)) * ft
      w_an = 0.0_ccs_real
      p_an = +(rho / 4.0_ccs_real) * (cos(2 * x_p(1)) + cos(2 * x_p(2))) * (ft**2)

      !print *, x_p(1), x_p(2), u_an, u_data(index_p), v_an, v_data(index_p)
      error_L2_local(1) = error_L2_local(1) + (u_an - u_data(index_p))**2
      error_L2_local(2) = error_L2_local(2) + (v_an - v_data(index_p))**2
      !error_L2_local(3) = error_L2_local(3) + (w_an - w_data(index_p))**2
      error_L2_local(3) = error_L2_local(3) + (p_an - p_data(index_p))**2

      error_Linf_local(1) = max(error_Linf_local(1), abs(u_an - u_data(index_p)))
      error_Linf_local(2) = max(error_Linf_local(2), abs(v_an - v_data(index_p)))
      !error_Linf_local(3) = max(error_Linf_local(3), abs(w_an - w_data(index_p)))
      error_Linf_local(3) = max(error_Linf_local(3), abs(p_an - p_data(index_p)))

    end do
    call restore_vector_data(u%values, u_data)
    call restore_vector_data(v%values, v_data)
    call restore_vector_data(w%values, w_data)
    call restore_vector_data(p%values, p_data)

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_AllReduce(error_L2_local, error_L2, size(error_L2), MPI_DOUBLE_PRECISION, MPI_SUM, par_env%comm, ierr)
      call MPI_AllReduce(error_Linf_local, error_Linf, size(error_Linf), MPI_DOUBLE_PRECISION, MPI_MAX, par_env%comm, ierr)
    class default
      call error_abort("ERROR: Unknown type")
    end select

    call get_global_num_cells(global_num_cells)
    error_L2(:) = sqrt(error_L2(:) / global_num_cells)

    if (is_root(par_env)) then
      if (first_time) then
        first_time = .false.
        open (newunit=io_unit, file="tgv2d-err.log", status="replace", form="formatted")
      else
        open (newunit=io_unit, file="tgv2d-err.log", status="old", form="formatted", position="append")
      end if
      fmt = '(I0,' // str(2 * size(error_L2)) // '(1x,e12.4))'
      write (io_unit, fmt) step, error_L2, error_Linf
      close (io_unit)
    end if

  end subroutine calc_tgv2d_error

  subroutine postproc_tgv(par_env, flow_fields)

    class(parallel_environment), allocatable, intent(in) :: par_env
    type(fluid), intent(in) :: flow_fields

    class(field), pointer:: u, v, w, p

    call get_field(flow_fields, "u", u)
    call get_field(flow_fields, "v", v)
    call get_field(flow_fields, "w", w)
    call get_field(flow_fields, "p", p)
    call calc_kinetic_energy(par_env, u, v, w)
    call calc_enstrophy(par_env, u, v, w)
    call calc_tgv2d_error(par_env, u, v, w, p, tgv2d_error_L2_global, tgv2d_error_Linf_global)
    nullify(u)
    nullify(v)
    nullify(w)
    nullify(p)

  end subroutine postproc_tgv

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

end module tgv2d_core

