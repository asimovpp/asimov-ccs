!v Submodule file io_visualisation_adios2.smod
!
!  Implementation (using MPI and ADIOS2) of parallel IO visualisation-realted routines
submodule(io_visualisation) io_visualisation_adios2
#include "ccs_macros.inc"

  use io, only: initialise_io, cleanup_io, configure_io, open_file, close_file, &
                write_array
  use adios2
  use adios2_types, only: adios2_io_process
  use utils, only: exit_print

  implicit none

contains

  !> Write the field data to file
  module subroutine write_fields(par_env, case_name, mesh, output_list, step, maxstep)

    use kinds, only: ccs_long
    use constants, only: ndim, adiosconfig
    use vec, only: get_vector_data, restore_vector_data
    use types, only: field_ptr, cell_locator
    use case_config, only: write_gradients
    use meshing, only: get_local_num_cells, get_global_num_cells, &
                       set_cell_location, &
                       get_global_index

    ! Arguments
    class(parallel_environment), allocatable, target, intent(in) :: par_env  !< The parallel environment
    character(len=:), allocatable, intent(in) :: case_name                   !< The case name
    type(ccs_mesh), intent(in) :: mesh                                       !< The mesh
    type(field_ptr), dimension(:), intent(inout) :: output_list              !< List of fields to output
    integer(ccs_int), optional, intent(in) :: step                           !< The current time-step count
    integer(ccs_int), optional, intent(in) :: maxstep                        !< The maximum time-step count

    ! Local variables
    character(len=:), allocatable :: sol_file     ! Solution file name
    character(len=:), allocatable :: adios2_file  ! ADIOS2 config file name
    character(len=:), allocatable :: data_name    ! String for storing data path in file

    class(io_environment), allocatable, save :: io_env
    class(io_process), allocatable, save :: sol_writer

    integer(ccs_long), dimension(1) :: sel_shape
    integer(ccs_long), dimension(1) :: sel_start
    integer(ccs_long), dimension(1) :: sel_count

    integer(ccs_long), dimension(2) :: sel2_shape
    integer(ccs_long), dimension(2) :: sel2_start
    integer(ccs_long), dimension(2) :: sel2_count

    real(ccs_real), dimension(:), pointer :: data

    integer(ccs_int) :: i

    integer(ccs_int) :: global_num_cells

    type(cell_locator) :: loc_p
    integer(ccs_int) :: index_global
    
    sol_file = case_name // '.sol.h5'
    adios2_file = case_name // adiosconfig

    if (present(step)) then
      ! Unsteady case
      if (step == 1) then
        call initialise_io(par_env, adios2_file, io_env)
        call configure_io(io_env, "sol_writer", sol_writer)
        call open_file(sol_file, "write", sol_writer)
      end if
    else
      ! Steady case
      call initialise_io(par_env, adios2_file, io_env)
      call configure_io(io_env, "sol_writer", sol_writer)
      call open_file(sol_file, "write", sol_writer)
    end if

    call get_global_num_cells(mesh, global_num_cells)

    ! Need to get data relating to first cell
    call set_cell_location(mesh, 1, loc_p)

    call get_global_index(loc_p, index_global)
    
    ! 1D data
    sel_shape(1) = global_num_cells
    sel_start(1) = index_global - 1
    call get_local_num_cells(mesh, sel_count(1))

    ! 2D data
    sel2_shape(1) = ndim
    sel2_shape(2) = global_num_cells
    sel2_start(1) = 0
    sel2_start(2) = index_global - 1
    sel2_count(1) = ndim
    call get_local_num_cells(mesh, sel2_count(2))

    ! Begin step
    call begin_step(sol_writer)

    ! Loop over output list and write out
    do i = 1, size(output_list)
      ! Check whether pointer is associated with a field
      if (.not. associated(output_list(i)%ptr)) exit

      call get_vector_data(output_list(i)%ptr%values, data)
      data_name = "/" // trim(output_list(i)%name)
      call write_array(sol_writer, data_name, sel_shape, sel_start, sel_count, data)
      call restore_vector_data(output_list(i)%ptr%values, data)
    end do

    ! Write out gradients, if required (e.g. for calculating enstrophy)
    if (write_gradients) then
      do i = 1, size(output_list)
        ! Check whether pointer is associated with a field
        if (.not. associated(output_list(i)%ptr)) exit

        ! x-gradient
        call get_vector_data(output_list(i)%ptr%x_gradients, data)
        data_name = "/d" // trim(output_list(i)%name) // "dx"
        call write_array(sol_writer, data_name, sel_shape, sel_start, sel_count, data)
        call restore_vector_data(output_list(i)%ptr%x_gradients, data)

        ! y-gradient
        call get_vector_data(output_list(i)%ptr%y_gradients, data)
        data_name = "/d" // trim(output_list(i)%name) // "dy"
        call write_array(sol_writer, data_name, sel_shape, sel_start, sel_count, data)
        call restore_vector_data(output_list(i)%ptr%y_gradients, data)

        ! z-gradient
        call get_vector_data(output_list(i)%ptr%z_gradients, data)
        data_name = "/d" // trim(output_list(i)%name) // "dz"
        call write_array(sol_writer, data_name, sel_shape, sel_start, sel_count, data)
        call restore_vector_data(output_list(i)%ptr%z_gradients, data)
      end do
    end if

    ! End step
    call end_step(sol_writer)

    ! Close the file and finalise ADIOS2 IO environment
    if (present(step)) then
      ! Unsteady case
      if (step == maxstep) then
        call close_file(sol_writer)
        call cleanup_io(io_env)
      end if
    else
      ! Steady case
      call close_file(sol_writer)
      call cleanup_io(io_env)
    end if

  end subroutine write_fields

  !> Begin a new step in the data file
  subroutine begin_step(io_proc)
    class(io_process), intent(inout) :: io_proc !< The IO process

    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_begin_step(io_proc%engine, ierr)

    class default
      call error_abort("Unknown IO process handler type")

    end select
  end subroutine begin_step

  !> End current step in the data file and flush IO
  subroutine end_step(io_proc)
    class(io_process), intent(inout) :: io_proc !< The IO process

    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_end_step(io_proc%engine, ierr)

    class default
      call error_abort("Unknown IO process handler type")

    end select
  end subroutine end_step

end submodule
