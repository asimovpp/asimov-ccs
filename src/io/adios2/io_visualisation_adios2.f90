!v Submodule file io_visualisation_adios2.smod
!
!  Implementation (using MPI and ADIOS2) of parallel IO visualisation-realted routines
submodule(io_visualisation) io_visualisation_adios2
#include "ccs_macros.inc"

  use io, only: initialise_io, cleanup_io, configure_io, open_file, close_file, &
                write_array, read_array
  use adios2
  use adios2_types, only: adios2_io_process
  use utils, only: exit_print, get_field, update, set_values
  use timers, only: timer_register, timer_start, timer_stop
  use types, only: field
  
  implicit none

  logical, save :: initial_step = .true.

contains

  module subroutine reset_io_visualisation_module()

    initial_step = .true.

  end subroutine

  !> Read the field data from file
  module subroutine read_fields(par_env, case_name, mesh, flow, step, maxstep)

    use kinds, only: ccs_long
    use constants, only: ndim, adiosconfig
    use types, only: cell_locator
    use meshing, only: get_global_num_cells, create_cell_locator, get_global_index, get_natural_index, get_local_num_cells
    use utils, only: get_natural_data
    use vec, only: get_vector_data, restore_vector_data
    use io, only: get_num_steps
    use vec, only: get_global_data_vec

    ! Arguments
    class(parallel_environment), allocatable, target, intent(in) :: par_env  !< The parallel environment
    character(len=:), allocatable, intent(in) :: case_name                   !< The case name
    type(ccs_mesh), intent(in) :: mesh                                       !< The mesh
    type(fluid), intent(inout) :: flow                                       !< The flow variables
    integer(ccs_int), optional, intent(in) :: step                           !< The current time-step count
    integer(ccs_int), optional, intent(in) :: maxstep                        !< The maximum time-step count

    ! Local variables
    character(len=:), allocatable :: sol_file     ! Solution file name
    character(len=:), allocatable :: adios2_file  ! ADIOS2 config file name
    character(len=:), allocatable :: data_name    ! String for storing data path in file
    integer(ccs_int) :: global_num_cells

    class(io_environment), allocatable, save :: io_env
    class(io_process), allocatable, save :: sol_reader

    integer(ccs_long), dimension(1) :: sel_shape
    integer(ccs_long), dimension(1) :: sel_start
    integer(ccs_long), dimension(1) :: sel_count

    integer(ccs_long), dimension(2) :: sel2_shape
    integer(ccs_long), dimension(2) :: sel2_start
    integer(ccs_long), dimension(2) :: sel2_count

    real(ccs_real), dimension(:), allocatable, target :: data
    real(ccs_real), dimension(:), allocatable :: re_order_data

    integer(ccs_int) :: timer_index_nat_data_output
    integer(ccs_int) :: timer_index_nat_data
    integer(ccs_int) :: timer_index_output
    integer(ccs_int) :: timer_index_grad

    integer(ccs_int) :: i

    integer(ccs_long) :: steps

    type(cell_locator) :: loc_p
    integer(ccs_int) :: index_global

    real(ccs_long), dimension(:), pointer :: output_data
    class(field), pointer :: phi

    sol_file = case_name // '.sol.h5'
    adios2_file = case_name // adiosconfig

    call timer_register("Get natural data (output)", timer_index_nat_data_output)
    call timer_register("Get natural data (grads)", timer_index_nat_data)
    call timer_register("Read output time", timer_index_output)
    call timer_register("Read gradients time", timer_index_grad)

    if (present(step)) then
      ! Unsteady case
      if (initial_step) then
        call initialise_io(par_env, adios2_file, io_env)
        call configure_io(io_env, "sol_reader", sol_reader)
        call open_file(sol_file, "read", sol_reader)

        initial_step = .false.
      end if
    else
      ! Steady case
      call initialise_io(par_env, adios2_file, io_env)
      call configure_io(io_env, "sol_reader", sol_reader)
      call open_file(sol_file, "read", sol_reader)
    end if

    call get_global_num_cells(global_num_cells)

    ! Need to get data relating to first cell
    call create_cell_locator(1, loc_p)

    call get_global_index(loc_p, index_global)

    ! 1D data
    sel_shape(1) = global_num_cells
    sel_start(1) = index_global - 1
    call get_local_num_cells(sel_count(1))

    ! 2D data
    sel2_shape(1) = ndim
    sel2_shape(2) = global_num_cells
    sel2_start(1) = 0
    sel2_start(2) = index_global - 1
    sel2_count(1) = ndim
    call get_local_num_cells(sel2_count(2))

    ! Get number of steps in solution file
    call get_num_steps(sol_reader, steps)
    steps = steps - 1  ! Set to max. step (count starts from 0)

    ! Loop over output list and write out
    call timer_start(timer_index_output)
    do i = 1, size(flow%fields)
      call get_field(flow, i, phi)
      if (phi%output) then
        ! XXX: This seems unnecessary?
        call timer_start(timer_index_nat_data_output)
        call get_natural_data(par_env, mesh, phi%values, data)
        call timer_stop(timer_index_nat_data_output)

        data_name = "/" // trim(phi%name)

        call read_array(sol_reader, data_name, sel_start, sel_count, data, steps)
        call get_vector_data (phi%values, output_data)
        output_data = data
        call restore_vector_data (phi%values, output_data)
        call get_global_data_vec(par_env, mesh, phi%values, re_order_data)
        ! re-ordering here. 

        call get_vector_data(phi%values, output_data)
        output_data = re_order_data

        ! XXX: This doesn't appear to do anything
        ! call get_local_num_cells(n_local)
        ! do index_p = 1, n_local
        !   call create_cell_locator(index_p, loc_p) 
        !   call get_global_index(loc_p, global_index_p)
        !   call get_natural_index(loc_p, natural_index_p)
        ! end do

        call restore_vector_data(phi%values, output_data)
        call update(phi%values)
        
      end if 
      nullify(phi)
    end do
    
    call timer_stop(timer_index_output)

    if (allocated(data)) then
      deallocate (data)
    end if

    ! Close the file and finalise ADIOS2 IO environment
    if (present(step)) then
      ! Unsteady case
      if (step == maxstep) then
        call close_file(sol_reader)
        call cleanup_io(io_env)
      end if
    else
      ! Steady case
      call close_file(sol_reader)
      call cleanup_io(io_env)
    end if
  
  end subroutine read_fields


  !> Write the field data to file
  module subroutine write_fields(par_env, case_name, mesh, flow, step, maxstep)

    use kinds, only: ccs_long
    use constants, only: ndim, adiosconfig
    use vec, only: get_vector_data, restore_vector_data
    use types, only: field_ptr, cell_locator
    use case_config, only: write_gradients
    use meshing, only: get_local_num_cells, get_global_num_cells, &
                       create_cell_locator, &
                       get_global_index
    use utils, only: get_natural_data
    use parallel, only: timer

    ! Arguments
    class(parallel_environment), allocatable, target, intent(in) :: par_env  !< The parallel environment
    character(len=:), allocatable, intent(in) :: case_name                   !< The case name
    type(ccs_mesh), intent(in) :: mesh                                       !< The mesh
    type(fluid), intent(inout) :: flow                                       !< The flow variables
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

    real(ccs_real), dimension(:), allocatable :: data

    integer(ccs_int) :: timer_index_nat_data_output
    integer(ccs_int) :: timer_index_nat_data
    integer(ccs_int) :: timer_index_output
    integer(ccs_int) :: timer_index_grad

    integer(ccs_int) :: i

    integer(ccs_int) :: global_num_cells

    type(cell_locator) :: loc_p
    integer(ccs_int) :: index_global

    class(field), pointer :: phi
    
    sol_file = case_name // '.sol.h5'
    adios2_file = case_name // adiosconfig

    call timer_register("Get natural data (output)", timer_index_nat_data_output)
    call timer_register("Get natural data (grads)", timer_index_nat_data)
    call timer_register("Write output time", timer_index_output)
    call timer_register("Write gradients time", timer_index_grad)

    if (present(step)) then
      ! Unsteady case
      if (initial_step) then
        call initialise_io(par_env, adios2_file, io_env)
        call configure_io(io_env, "sol_writer", sol_writer)
        call open_file(sol_file, "write", sol_writer)

        initial_step = .false.
      end if
    else
      ! Steady case
      call initialise_io(par_env, adios2_file, io_env)
      call configure_io(io_env, "sol_writer", sol_writer)
      call open_file(sol_file, "write", sol_writer)
    end if

    call get_global_num_cells(global_num_cells)

    ! Need to get data relating to first cell
    call create_cell_locator(1, loc_p)

    call get_global_index(loc_p, index_global)

    ! 1D data
    sel_shape(1) = global_num_cells
    sel_start(1) = index_global - 1
    call get_local_num_cells(sel_count(1))

    ! 2D data
    sel2_shape(1) = ndim
    sel2_shape(2) = global_num_cells
    sel2_start(1) = 0
    sel2_start(2) = index_global - 1
    sel2_count(1) = ndim
    call get_local_num_cells(sel2_count(2))

    ! Begin step
    call begin_step(sol_writer)

    ! Loop over output list and write out
    call timer_start(timer_index_output)
    do i = 1, size(flow%fields)
      call get_field(flow, i, phi)
      if (phi%output) then
      
        call timer_start(timer_index_nat_data_output)
        call get_natural_data(par_env, mesh, phi%values, data)
        call timer_stop(timer_index_nat_data_output)
        data_name = "/" // trim(phi%name)
        call write_array(sol_writer, data_name, sel_shape, sel_start, sel_count, data)

        ! Store residuals if available
        if (allocated(phi%residuals)) then
          call timer_start(timer_index_nat_data_output)
          call get_natural_data(par_env, mesh, phi%residuals, data)
          call timer_stop(timer_index_nat_data_output)
          data_name = "/" // trim(phi%name // "_res")
          call write_array(sol_writer, data_name, sel_shape, sel_start, sel_count, data)
        end if
      end if
    end do
    call timer_stop(timer_index_output)
   

    ! Write out gradients, if required (e.g. for calculating enstrophy)
    if (write_gradients) then
      call timer_start(timer_index_grad)
      do i = 1, size(flow%fields)
        call get_field(flow, i, phi)

        if (phi%output) then
          ! x-gradient
          call timer_start(timer_index_nat_data)
          call get_natural_data(par_env, mesh, phi%x_gradients, data)
          call timer_stop(timer_index_nat_data)
          data_name = "/d" // trim(phi%name) // "dx"
          call write_array(sol_writer, data_name, sel_shape, sel_start, sel_count, data)

          ! y-gradient
          call timer_start(timer_index_nat_data)
          call get_natural_data(par_env, mesh, phi%y_gradients, data)
          call timer_stop(timer_index_nat_data)
          data_name = "/d" // trim(phi%name) // "dy"
          call write_array(sol_writer, data_name, sel_shape, sel_start, sel_count, data)

          ! z-gradient
          call timer_start(timer_index_nat_data)
          call get_natural_data(par_env, mesh, phi%z_gradients, data)
          call timer_stop(timer_index_nat_data)
          data_name = "/d" // trim(phi%name) // "dz"
          call write_array(sol_writer, data_name, sel_shape, sel_start, sel_count, data)
        end if
      end do
      call timer_stop(timer_index_grad)
    end if

    if (allocated(data)) then
      deallocate (data)
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
