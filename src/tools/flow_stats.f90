!v Set of tools to measure flow statistics
module flow_stats

  use mpi

  use kinds, only: ccs_real, ccs_int
  use types, only: fluid, cell_locator, field
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi

  use fields, only: get_field
  use meshing, only: get_local_num_cells, get_global_num_cells, create_cell_locator, get_volume
  use parallel, only: is_root
  use timestepping, only: get_timestep, timestepping_is_active

  implicit none

  private
  public :: report_cfl

contains

  subroutine report_cfl(par_env, flow)
    use logging, only: log_unit_out
    class(parallel_environment), intent(in) :: par_env
    type(fluid), intent(in) :: flow

    class(field), pointer :: u, v, w

    real(ccs_real) :: cfl_max, cfl_avg, cfl_i
    real(ccs_real) :: dt
    real(ccs_real) :: dx, V_p
    type(cell_locator) :: loc_p

    real(ccs_real) :: vel

    integer(ccs_int) :: i, nlocal, nglobal

    integer :: ierr

    if (.not. timestepping_is_active()) then
      return
    end if

    dt = get_timestep()

    cfl_max = 0.0_ccs_real
    cfl_avg = 0.0_ccs_real

    call get_field(flow, "u", u)
    call get_field(flow, "v", v)
    call get_field(flow, "w", w)

    call get_local_num_cells(nlocal)
    call get_global_num_cells(nglobal)
    do i = 1, nlocal
      vel = norm2([u%values_ro(i), v%values_ro(i), w%values_ro(i)])

      call create_cell_locator(i, loc_p)
      call get_volume(loc_p, V_p)
      dx = get_lscale(V_p)

      cfl_i = cfl(vel, dt, dx)
      cfl_max = max(cfl_i, cfl_max)
      cfl_avg = cfl_avg + cfl_i / nglobal
    end do

    nullify (u, v, w)

    !! Reduce the MAX/AVG CFL numbers
    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Allreduce(MPI_IN_PLACE, cfl_max, 1, MPI_DOUBLE, MPI_MAX, par_env%comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, cfl_avg, 1, MPI_DOUBLE, MPI_SUM, par_env%comm, ierr)
    class default
      error stop "Unsupported parallel environment"
    end select

    if (is_root(par_env)) then
      write (log_unit_out, *) "CFL Max: ", cfl_max
      write (log_unit_out, *) "CFL Avg: ", cfl_avg
    end if

  end subroutine report_cfl

  pure real(ccs_real) function get_lscale(vol) result(l)
    real(ccs_real), intent(in) :: vol

    l = vol**(1.0 / 3.0)
  end function get_lscale

  pure real(ccs_real) function cfl(u, dt, dx)
    real(ccs_real), intent(in) :: u
    real(ccs_real), intent(in) :: dt
    real(ccs_real), intent(in) :: dx

    cfl = u * dt / dx

  end function cfl

end module flow_stats
