!> Program file for TaylorGreenVortex case
program tgv
#include "ccs_macros.inc"

  use petscsys
  use petscvec

  use core
  use constants, only: ndim
  use meshing, only: nullify_mesh_object
  use kinds, only: ccs_real, ccs_int, ccs_long
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      is_root
  use parallel_types, only: parallel_environment
  use types, only: fluid, field
  use utils, only: exit_print, &
                   calc_kinetic_energy, calc_enstrophy, &
                   add_field_to_outputlist, get_field, add_field, &
                   set_is_field_solved, &
                   str, debug_print, &
                   dealloc_fluid_fields
  use timers, only: timer_init, timer_register_start, timer_register, timer_start, timer_stop, timer_print, &
                    timer_get_time, timer_print_all, timer_export_csv, timer_get_index

  implicit none

  class(parallel_environment), allocatable :: par_env
  class(parallel_environment), allocatable :: shared_env

  type(ccs_options) :: run_options

  integer(ccs_int) :: irank  ! MPI rank ID
  integer(ccs_int) :: isize  ! Size of MPI world

  integer(ccs_int) :: timer_index_total
  integer(ccs_int) :: timer_index_init
  integer(ccs_int) :: timer_index_io_sol
  integer(ccs_int) :: timer_index_sol

  double precision :: sol_time, io_time

  type(fluid) :: flow_fields

  ! Launch MPI
  call initialise_parallel_environment(par_env)
  call timer_init()

  call get_config(par_env, run_options)
  call configure_parallelism(run_options, par_env, shared_env)
  irank = par_env%proc_id
  isize = par_env%num_procs

  call timer_register_start("Elapsed time", timer_index_total, is_total_time=.true.)

  call timer_register_start("Init time", timer_index_init)

  if (is_root(par_env)) print *, "Starting ", run_options%paths%case_name, " case!"

  call initialise_mesh(par_env, shared_env, run_options)

  ! Initialise fields
  if (is_root(par_env)) print *, "Initialise fields"

  call initialise_fields(par_env, run_options, flow_fields)
  
  ! Initialise velocity field
  if (is_root(par_env)) print *, "Initialise velocity field"
  call initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)

  call timer_stop(timer_index_init)
  call run_solver(par_env, run_options, eval_sources, postproc_tgv, flow_fields)
  call timer_stop(timer_index_total)

  call timer_print_all(par_env)
  call timer_export_csv(par_env)

  call timer_get_index("Solver time inc I/O", timer_index_sol)
  call timer_get_index("I/O time for solution", timer_index_io_sol)
  call timer_get_time(timer_index_sol, sol_time)
  call timer_get_time(timer_index_io_sol, io_time)
  if (is_root(par_env)) then
    write(*,'(A30, F10.4, A)') "Solver time no I/O:", sol_time - io_time, " s"
    write(*,'(A30, F10.4, A)') "Average time/step (no I/O):", (sol_time - io_time) / run_options%solve%num_steps, " s"
  end if

  call dealloc_fluid_fields(flow_fields)
  call nullify_mesh_object()
  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

contains

  pure subroutine get_init_flow(loc_p, field_name, init_val)
    use types, only: cell_locator
    use meshing, only: get_centre
    type(cell_locator), intent(in) :: loc_p
    character(len=*), intent(in) :: field_name
    real(ccs_real), intent(inout) :: init_val
    real(ccs_real), dimension(ndim) :: x_p

    select case (field_name)
    case ("u")
      call get_centre(loc_p, x_p)
      init_val = sin(x_p(1)) * cos(x_p(2)) * cos(x_p(3))
    case ("v")
      call get_centre(loc_p, x_p)
      init_val = -cos(x_p(1)) * sin(x_p(2)) * cos(x_p(3))
    case default
      init_val = init_val + 0.0_ccs_real ! Keep the default value
    end select
    
  end subroutine get_init_flow
  
  pure subroutine get_init_mass_flux(loc_f, init_val)
    use types, only: face_locator
    use meshing, only: get_face_normal, get_centre
    type(face_locator), intent(in) :: loc_f
    real(ccs_real), intent(inout) :: init_val
    real(ccs_real), dimension(ndim) :: x_f
    real(ccs_real), dimension(ndim) :: face_normal

    call get_face_normal(loc_f, face_normal)
    call get_centre(loc_f, x_f)
  
    init_val = sin(x_f(1)) * cos(x_f(2)) * cos(x_f(3)) * face_normal(1) &
             - cos(x_f(1)) * sin(x_f(2)) * cos(x_f(3)) * face_normal(2)

  end subroutine get_init_mass_flux

  subroutine postproc_tgv(par_env, flow_fields)
    use timestepping, only: get_current_step

    class(parallel_environment), allocatable, intent(in) :: par_env
    type(fluid), intent(in) :: flow_fields
    integer(ccs_int) :: step
    logical :: first_time

    class(field), pointer:: u, v, w

    call get_field(flow_fields, "u", u)
    call get_field(flow_fields, "v", v)
    call get_field(flow_fields, "w", w)
    call calc_kinetic_energy(par_env, u, v, w)
    call calc_enstrophy(par_env, u, v, w)

    call get_current_step(step)

    if (modulo(step, 20) == 0) then
      first_time = (step == 20)
    ! Centre line
      call export_line(par_env, [1.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real], &
                                [0.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real], 40, u, run_options%mesh%cps, first_time)
  
      call export_line(par_env, [1.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real], &
                                [0.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real], 40, v, run_options%mesh%cps, first_time)
  
      call export_line(par_env, [1.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real], &
                                [0.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real], 40, w, run_options%mesh%cps, first_time)

    ! 1/4 line                            
      call export_line(par_env, [1.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real], &
                                [0.0_ccs_real, run_options%mesh%domain_size / 4.0_ccs_real, 0.0_ccs_real], 40, u, run_options%mesh%cps, .false.)
  
      call export_line(par_env, [1.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real], &
                                [0.0_ccs_real, run_options%mesh%domain_size / 4.0_ccs_real, 0.0_ccs_real], 40, v, run_options%mesh%cps, .false.)
  
      call export_line(par_env, [1.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real], &
                                [0.0_ccs_real, run_options%mesh%domain_size / 4.0_ccs_real, 0.0_ccs_real], 40, w, run_options%mesh%cps, .false.)
    end if

    nullify(u)
    nullify(v)
    nullify(w)

  end subroutine postproc_tgv


  ! Get natural index from x_loc and cps, only works for 3D cartesian meshes
  pure subroutine get_natural_index(x_loc, cps, ii)
    real(ccs_real), intent(in), dimension(3) :: x_loc
    integer(ccs_int), intent(in) :: cps
    integer(ccs_int), intent(out) :: ii
    integer(ccs_int) :: nx, ny, nz
    real(ccs_real) :: x, y, z

    ii = 0
    nx = cps
    ny = cps
    nz = cps
    x = x_loc(1)/run_options%mesh%domain_size - (1.0_ccs_real/nx)/2
    y = x_loc(2)/run_options%mesh%domain_size - (1.0_ccs_real/ny)/2
    z = x_loc(3)/run_options%mesh%domain_size - (1.0_ccs_real/nz)/2

    ii = floor(x*nx + 0.5, ccs_int) 
    ii = ii + nx* floor(y*ny + 0.5, ccs_int)
    ii = ii + nx*ny* floor(z*nz + 0.5, ccs_int)
    ii = ii +1

  end subroutine


  !> Get the 8 natural indice of cells part of the trilineal stencil around x_loc
  pure subroutine get_stencil(x_loc, cps, natural_ids)
    real(ccs_real), intent(in), dimension(3) :: x_loc
    integer(ccs_int), intent(in) :: cps
    integer(ccs_int), intent(out), dimension(8) :: natural_ids
    real(ccs_real), dimension(3) :: s_loc
    real(ccs_real) :: h

    h = run_options%mesh%domain_size / cps /2
    s_loc(:) = x_loc(:) + [ h, h, h]
    call get_natural_index(s_loc, cps, natural_ids(1))

    s_loc(:) = x_loc(:) + [ -h, -h, h]
    call get_natural_index(s_loc, cps, natural_ids(2))

    s_loc(:) = x_loc(:) + [ -h, h, h]
    call get_natural_index(s_loc, cps, natural_ids(3))
    
    s_loc(:) = x_loc(:) + [ h, -h, h]
    call get_natural_index(s_loc, cps, natural_ids(4))

    s_loc(:) = x_loc(:) + [ h, h, -h]
    call get_natural_index(s_loc, cps, natural_ids(5))

    s_loc(:) = x_loc(:) + [ -h, -h, -h]
    call get_natural_index(s_loc, cps, natural_ids(6))

    s_loc(:) = x_loc(:) + [ -h, h, -h]
    call get_natural_index(s_loc, cps, natural_ids(7))
    
    s_loc(:) = x_loc(:) + [ h, -h, -h]
    call get_natural_index(s_loc, cps, natural_ids(8))

  end subroutine

  ! get array of local ids from natural indices, cells not on local rank will 
  ! have a "0" local id
  pure subroutine get_local_ids(natural_ids, local_ids)
    use ccs_base, only: mesh
    use meshing, only: get_local_num_cells
    integer(ccs_int), intent(in), dimension(8) :: natural_ids
    integer(ccs_int), intent(out), dimension(8) :: local_ids
    integer(ccs_int) :: local_num_cells
    integer(ccs_int), dimension(1) :: loc 
    integer(ccs_int) :: i

    local_ids(:) = 0
    call get_local_num_cells(local_num_cells)
    do i=1, 8
      loc = findloc(mesh%topo%natural_indices(1:local_num_cells), natural_ids(i))
      if (loc(1) /= 0) then
        local_ids(i) = loc(1)
      end if
    end do

  end subroutine

  ! Gather field values of stencil points for location x_loc. Non local stencil value will be filled with 0
  subroutine gather_local_stencil_values(x_loc, cps, phi, stencil_values)
    use vec, only: get_vector_data_readonly, restore_vector_data_readonly
    real(ccs_real), intent(in), dimension(3) :: x_loc
    integer(ccs_int), intent(in) :: cps
    class(field), pointer, intent(in) :: phi
    real(ccs_real), dimension(8), intent(out) :: stencil_values
    integer(ccs_int), dimension(8) :: natural_ids, local_ids
    real(ccs_real), dimension(:), pointer :: phi_data
    integer(ccs_int) :: i

    stencil_values(:) = 0.0_ccs_real

    call get_stencil(x_loc, cps, natural_ids)
    call get_local_ids(natural_ids, local_ids)

    call get_vector_data_readonly(phi%values, phi_data)

    do i=1, 8
      if (local_ids(i) /= 0) then
        stencil_values(i) = phi_data(local_ids(i))
      end if
    end do

    call restore_vector_data_readonly(phi%values, phi_data)

  end subroutine


  ! Stores in `value` the interpolated value for location x_loc
  ! Handles parallel communication, every rank gets the interpolated value
  ! -> lots of ways to improve performance: - caching natural index -> local index map
  !                                         - not using all reduce, but point to point if needed
  !                                         - doing parallel communication for every single interpolation at once
  subroutine get_interpolated_value(par_env, x_loc, cps, phi, interpolated_value)
    use parallel_types_mpi, only: parallel_environment_mpi

    class(parallel_environment), allocatable, intent(in) :: par_env
    real(ccs_real), intent(in), dimension(3) :: x_loc
    integer(ccs_int), intent(in) :: cps
    class(field), pointer, intent(in) :: phi
    real(ccs_real), intent(out) :: interpolated_value
    real(ccs_real), dimension(8) :: stencil_values
    integer :: ierr

    call gather_local_stencil_values(x_loc, cps, phi, stencil_values)

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Allreduce(MPI_IN_PLACE, stencil_values, 8, MPI_DOUBLE_PRECISION, MPI_SUM, par_env%comm, ierr)
    class default
      call error_abort("invalid parallel environment")
    end select

    interpolated_value = sum(stencil_values)/8.0_ccs_real

  end subroutine


  ! Export phi on a line defined by its direction and offset from the centre of the domain
  subroutine export_line(par_env, direction, offset, num_values, phi, cps, first_time)
    use timestepping, only: get_current_time, get_current_step

    class(parallel_environment), allocatable, intent(in) :: par_env
    real(ccs_real), dimension(3), intent(in) :: direction
    real(ccs_real), dimension(3), intent(in) :: offset
    integer(ccs_int), intent(in) :: num_values
    class(field), pointer, intent(in) :: phi
    integer(ccs_int), intent(in) :: cps
    logical, intent(in) :: first_time

    real(ccs_real), dimension(:,:), allocatable :: x_loc
    real(ccs_real), dimension(:), allocatable :: values
    real(ccs_real) :: time, domain_size
    integer(ccs_int) :: step
    integer(ccs_int) :: i
    integer :: io_unit
    integer(ccs_int) :: timer_id

    call timer_register_start("export_line", timer_id)

    allocate(values(num_values))
    allocate(x_loc(3, num_values))

    domain_size = run_options%mesh%domain_size
    call get_current_time(time)
    call get_current_step(step)

    do i=1, num_values

      x_loc(:, i) = [domain_size/2, domain_size/2, domain_size/2]*(1.0_ccs_real - direction) + offset &
      + direction*i*domain_size/(num_values+1)

      call get_interpolated_value(par_env, x_loc(:, i), cps, phi, values(i))
    end do

    if (is_root(par_env)) then

      if (first_time) then
        open (newunit=io_unit, file="line_export_"// trim(phi%name) //".dat", status="replace", form="formatted")
        write (io_unit, *) "# step, time, x, y, z, phi"
      else
        open (newunit=io_unit, file="line_export_"// trim(phi%name) //".dat", status="old", form="formatted", position="append")
      end if

      do i=1, num_values
        write (io_unit, "(i0, 1x, 5(e23.16, 1x))") step, time, x_loc(:, i), values(i)
      end do

      close (io_unit)
    end if

    call timer_stop(timer_id)

  end subroutine



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

end program tgv
