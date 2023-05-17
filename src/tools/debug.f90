!v Set of tools used for debugging purposes
module debug
#include "ccs_macros.inc"

  use kinds
  use types

contains


!v Export fields as well as mpi rank and local cell ID to text file
  subroutine dump_flow_tofile(par_env, mesh, u, v, w, p)

    use constants, only: ndim
    use types, only: cell_locator
    use vec, only: get_vector_data, restore_vector_data
    use meshing, only: set_cell_location, get_centre, &
                       get_local_num_cells
    use parallel, only: sync
    use parallel_types, only: parallel_environment

    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env !< The parallel environment
    class(ccs_mesh), intent(in) :: mesh
    class(field), intent(inout) :: u, v, w, p

    ! Local variables
    integer :: io_unit, irank
    integer(ccs_int) :: n_local
    integer(ccs_int) :: index_p
    type(cell_locator) :: loc_p

    real(ccs_real), dimension(ndim) :: x_p
    real(ccs_real), dimension(:), pointer :: u_data, v_data, w_data, p_data

    call get_local_num_cells(mesh, n_local)

    call get_vector_data(u%values, u_data)
    call get_vector_data(v%values, v_data)
    call get_vector_data(w%values, w_data)

    call get_vector_data(p%values, p_data)

    do irank = 0, par_env%num_procs -1

      call sync(par_env)

      if (irank == par_env%proc_id) then
        if (irank == 0) then
          open (newunit=io_unit, file="flow_dump.dat", status="replace", form="formatted")
        else
          open (newunit=io_unit, file="flow_dump.dat", status="old", form="formatted", position="append")
        end if

        do index_p = 1, n_local
          call set_cell_location(mesh, index_p, loc_p)
          call get_centre(loc_p, x_p)

          write (io_unit, '(7(1x,e12.4), 1x,i12, 1x,i12)') x_p(:), u_data(index_p), v_data(index_p), w_data(index_p), p_data(index_p), irank, index_p

        end do
        close (io_unit)
      end if
    end do

    call restore_vector_data(u%values, u_data)
    call restore_vector_data(v%values, v_data)
    call restore_vector_data(w%values, w_data)
    call restore_vector_data(p%values, p_data)

  end subroutine

end module debug
