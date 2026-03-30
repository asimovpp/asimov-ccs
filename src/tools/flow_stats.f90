!v Set of tools to measure flow statistics
module flow_stats

  use kinds, only: ccs_real
  use types, only: parallel_environment, fluid

  implicit none
  
  private

contains

  subroutine report_cfl(par_env, flow)
    class(parallel_environment), intent(in) :: par_env
    type(fluid), intent(in) :: flow

    if (.not. transient) then
       return
    end if

    cfl_max = 0.0_ccs_real
    cfl_avg = 0.0_ccs_real

    do i = 1, nlocal
       vel = get_uscale(u, v, w)
       dx = get_lscale(V_p)
       
       cfl_i = cfl(vel, dt, dx)
       cfl_max = max(cfl_i, cfl_max)
       cfl_avg = cfl_avg + cfl_i / nglobal
    end do

    if (is_root(par_env)) then
       print *, "CFL Max: ", cfl_max
       print *, "CFL Avg: ", cfl_avg
    end if
    
  end subroutine report_cfl

  pure real(ccs_real) function get_uscale(u, v, w) result(vel)
    real(ccs_real), intent(in) :: u, v, w

    vel = sqrt(u**2 + v**2 + w**2)
  end function get_uscale
  
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
