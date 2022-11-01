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
  module subroutine write_fields(par_env, case_name, mesh, output_list, step, maxstep, dt)

    use kinds, only: ccs_long
    use constants, only: ndim, adiosconfig
    use vec, only : get_vector_data, restore_vector_data
    use types, only: field_ptr
    use case_config, only: write_gradients

    ! Arguments
    class(parallel_environment), allocatable, target, intent(in) :: par_env  !< The parallel environment
    character(len=:), allocatable, intent(in) :: case_name                   !< The case name
    type(ccs_mesh), intent(in) :: mesh                                       !< The mesh
    type(field_ptr), dimension(:), intent(inout) :: output_list              !< List of fields to output
    integer(ccs_int), optional, intent(in) :: step                           !< The current time-step count
    integer(ccs_int), optional, intent(in) :: maxstep                        !< The maximum time-step count
    real(ccs_real), optional, intent(in) :: dt                               !< The time-step size

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

    integer(ccs_int) :: ierr
    integer(ccs_int) :: i
    integer(ccs_int), save :: step_counter = 0

    sol_file = case_name//'.sol.h5'
    adios2_file = case_name//adiosconfig

    if (present(step)) then
      ! Unsteady case
      if (step == 1) then
        call initialise_io(par_env, adios2_file, io_env)
        call configure_io(io_env, "sol_writer", sol_writer)
        call open_file(sol_file, "write", sol_writer)
      endif
    else
      ! Steady case
      call initialise_io(par_env, adios2_file, io_env)
      call configure_io(io_env, "sol_writer", sol_writer)
      call open_file(sol_file, "write", sol_writer)
    endif

    ! 1D data
    sel_shape(1) = mesh%topo%global_num_cells
    sel_start(1) = mesh%topo%global_indices(1) - 1
    sel_count(1) = mesh%topo%local_num_cells

    ! 2D data
    sel2_shape(1) = ndim
    sel2_shape(2) = mesh%topo%global_num_cells
    sel2_start(1) = 0
    sel2_start(2) = mesh%topo%global_indices(1) - 1
    sel2_count(1) = ndim
    sel2_count(2) = mesh%topo%local_num_cells

    ! Begin step
    call begin_step(sol_writer)

    ! Loop over output list and write out
    do i = 1, size(output_list)
      call get_vector_data(output_list(i)%ptr%values, data)
      data_name = "/" // trim(output_list(i)%name)
      call write_array(sol_writer, data_name, sel_shape, sel_start, sel_count, data)
      call restore_vector_data(output_list(i)%ptr%values, data)
    end do

    ! Write out gradients, if required (e.g. for calculating enstrophy)
    if (write_gradients) then
      do i = 1, size(output_list)
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
    endif

    ! End step
    call end_step(sol_writer)

    ! Close the file and finalise ADIOS2 IO environment
    if (present(step)) then
      ! Unsteady case
      if (step == maxstep) then
        call close_file(sol_writer)
        call cleanup_io(io_env)
      endif
    else
      ! Steady case
      call close_file(sol_writer)
      call cleanup_io(io_env)
    endif

  end subroutine

  !< Begin a new step in the data file
  subroutine begin_step(io_proc)
    class(io_process), intent(inout) :: io_proc

    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_begin_step(io_proc%engine, ierr)

    class default
      call error_abort("Unknown IO process handler type")

    end select
  end subroutine

  !< End current step in the data file and flush IO
  subroutine end_step(io_proc)
    class(io_process), intent(inout) :: io_proc

    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_end_step(io_proc%engine, ierr)

    class default
      call error_abort("Unknown IO process handler type")

    end select
  end subroutine

end submodule