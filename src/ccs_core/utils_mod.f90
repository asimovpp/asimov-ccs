!v Module file utils.mod
!
!  Provides utility functions for ASiMoV-CCS, these should be polymorphic on their input
!  and call type-specific implementations of the interface in other modules.

module utils
#include "ccs_macros.inc"

  use iso_c_binding

  use vec, only: set_vector_values, update_vector, begin_update_vector, end_update_vector, &
                 initialise_vector, set_vector_size, &
                 set_vector_values_mode, set_vector_values_row, set_vector_values_entry, &
                 clear_vector_values_entries, &
                 mult_vec_vec, scale_vec, zero_vector
  use mat, only: set_matrix_values, update_matrix, begin_update_matrix, end_update_matrix, &
                 initialise_matrix, finalise_matrix, set_matrix_size, &
                 set_matrix_values_mode, set_matrix_values_row, set_matrix_values_col, set_matrix_values_entry, &
                 clear_matrix_values_entries, zero_matrix
  use solver, only: initialise_equation_system
  use kinds, only: ccs_int, ccs_real

  implicit none

  private

  public :: set_values
  public :: set_entry
  public :: clear_entries
  public :: begin_update
  public :: end_update
  public :: update
  public :: finalise
  public :: initialise
  public :: set_size
  public :: mult
  public :: zero
  public :: set_mode
  public :: set_row
  public :: set_col
  public :: str
  public :: debug_print
  public :: exit_print
  public :: calc_kinetic_energy
  public :: calc_enstrophy

  !> Generic interface to set values on an object.
  interface set_values
    module procedure set_vector_values
    module procedure set_matrix_values
  end interface set_values

  !> Generic interface to store a value for later setting.
  interface set_entry
    module procedure set_vector_values_entry
    module procedure set_matrix_values_entry
  end interface set_entry

  !> Generic interface to set the storage mode (ADD/INSERT) of an object.
  interface set_mode
    module procedure set_vector_values_mode
    module procedure set_matrix_values_mode
  end interface set_mode

  !> Generic interface to set the working row.
  interface set_row
    module procedure set_vector_values_row
    module procedure set_matrix_values_row
  end interface set_row
  !> Generic interface to set the working column.
  interface set_col
    module procedure set_matrix_values_col
  end interface set_col

  !> Generic interface to indicate an object is ready for use.
  interface finalise
    module procedure finalise_matrix
  end interface finalise

  !> Generic interface to clear stored entries.
  interface clear_entries
    module procedure clear_vector_values_entries
    module procedure clear_matrix_values_entries
  end interface clear_entries

  !> Generic interface to perform parallel update of an object.
  interface update
    module procedure update_vector
    module procedure update_matrix
  end interface update

  !v Generic interface to begin parallel update of an object.
  !
  !  This is to allow overlapping comms and compute.
  interface begin_update
    module procedure begin_update_vector
    module procedure begin_update_matrix
  end interface begin_update

  !v Generic interface to end parallel update of an object.
  !
  !  This is to allow overlapping comms and compute.
  interface end_update
    module procedure end_update_vector
    module procedure end_update_matrix
  end interface end_update

  !> Generic interface to initialse vectors, matrices and linear systems
  interface initialise
    module procedure initialise_vector
    module procedure initialise_matrix
    module procedure initialise_equation_system
  end interface initialise

  !> Generic interface to set vector and matrix sizes
  interface set_size
    module procedure set_vector_size
    module procedure set_matrix_size
  end interface set_size

  !>  Generic interface to perform multiplications
  interface mult
    module procedure mult_vec_vec
  end interface mult

  !> Generic interface to zero an object
  interface zero
    module procedure zero_vector
    module procedure zero_matrix
  end interface zero

  !> Generic interface to converting numbers to strings
  interface str
    module procedure int2str
    module procedure real2str
  end interface str

  !> Generic interface to debug printer
  interface debug_print
    module procedure debug_print_actual
    module procedure noop
  end interface debug_print

contains

  !> Print a message, along with with its location.
  subroutine debug_print_actual(msg, filepath, line)
    use mpi

    character(*), intent(in) :: msg      !< text to be printed
    character(*), intent(in) :: filepath !< calling file, added by __FILE__ macro
    integer, intent(in) :: line          !< line number in calling file, added by __LINE__ macro

    character(len=:), allocatable :: filename
    integer :: slash_position
    integer :: rank, ierror
    logical :: init_flag

    slash_position = scan(trim(filepath), "/", back=.true.)
    if (slash_position .gt. 0) then
      filename = filepath(slash_position + 1:len_trim(filepath))
    else
      filename = filepath
    end if

    call mpi_initialized(init_flag, ierror)
    if (init_flag) then
      call mpi_comm_rank(MPI_COMM_WORLD, rank, ierror)
      print *, trim(filename), "(", int2str(line), ")[", int2str(rank), "] : ", msg
    else
      print *, trim(filename), "(", int2str(line), ") : ", msg
    end if
  end subroutine

  !> No-op routine, does nothing.
  subroutine noop()
  end subroutine

  !> Convert integer to string.
  function int2str(in_int, format_str) result(out_string)
    integer(ccs_int), intent(in) :: in_int           !< integer to convert
    character(*), optional, intent(in) :: format_str !< format string to use
    character(:), allocatable :: out_string          !< formatted string from input integer

    character(32) :: tmp_string

    if (present(format_str)) then
      write (tmp_string, format_str) in_int
    else
      write (tmp_string, *) in_int
    end if
    out_string = trim(adjustl(tmp_string))
  end function

  !> Convert real to string.
  function real2str(in_real, format_str) result(out_string)
    real(ccs_real), intent(in) :: in_real            !< real number to convert
    character(*), optional, intent(in) :: format_str !< format string to use
    character(:), allocatable :: out_string          !< formatted string from input real

    character(32) :: tmp_string

    if (present(format_str)) then
      write (tmp_string, format_str) in_real
    else
      write (tmp_string, *) in_real
    end if
    out_string = trim(adjustl(tmp_string))
  end function

  !> Print a message and stop program execution.
  subroutine exit_print(msg, filepath, line)
    character(*), intent(in) :: msg      !< text to be printed
    character(*), intent(in) :: filepath !< calling file, added by __FILE__ macro
    integer, intent(in) :: line          !< line number in calling file, added by __LINE__ macro

    call debug_print(msg, filepath, line)
    stop 1
  end subroutine exit_print

  !TODO: move this subroutine to more appropriate module
  !> Calculate kinetic energy over density
  subroutine calc_kinetic_energy(par_env, mesh, t, u, v, w)

    use constants, only: ndim, ccs_string_len
    use types, only: field, ccs_mesh
    use vec, only: get_vector_data, restore_vector_data
    use parallel, only: allreduce, error_handling
    use parallel_types_mpi, only: parallel_environment_mpi
    use parallel_types, only: parallel_environment
    use mpi

    class(parallel_environment), allocatable, intent(in) :: par_env !< parallel environment
    type(ccs_mesh), intent(in) :: mesh !< the mesh
    integer(ccs_int), intent(in) :: t !< timestep
    class(field), intent(inout) :: u !< solve x velocity field
    class(field), intent(inout) :: v !< solve y velocity field
    class(field), intent(inout) :: w !< solve z velocity field

    real(ccs_real) :: ek_local, ek_global, volume_local, volume_global
    real(ccs_real), dimension(:), pointer :: u_data, v_data, w_data
    real(ccs_real) :: rho
    integer(ccs_int) :: index_p
    character(len=ccs_string_len) :: fmt
    integer(ccs_int) :: ierr

    logical, save :: first_time = .true.
    integer :: io_unit

    rho = 1.0_ccs_real

    ek_local = 0.0_ccs_real
    ek_global = 0.0_ccs_real
    volume_local = 0.0_ccs_real
    volume_global = 0.0_ccs_real

    call get_vector_data(u%values, u_data)
    call get_vector_data(v%values, v_data)
    call get_vector_data(w%values, w_data)

    do index_p = 1, mesh%topo%local_num_cells

      ek_local = ek_local + 0.5 * rho * mesh%geo%volumes(index_p) * &
                 (u_data(index_p)**2 + v_data(index_p)**2 + w_data(index_p)**2)

      volume_local = volume_local + mesh%geo%volumes(index_p)

    end do

    call restore_vector_data(u%values, u_data)
    call restore_vector_data(v%values, v_data)
    call restore_vector_data(w%values, w_data)

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_AllReduce(ek_local, ek_global, 1, MPI_DOUBLE, MPI_SUM, par_env%comm, ierr)
      call error_handling(ierr, "mpi", par_env)
      call MPI_AllReduce(volume_local, volume_global, 1, MPI_DOUBLE, MPI_SUM, par_env%comm, ierr)
      call error_handling(ierr, "mpi", par_env)
    class default
      call error_abort("ERROR: Unknown type")
    end select

    ek_global = ek_global / volume_global

    if (par_env%proc_id == par_env%root) then
      if (first_time) then
        first_time = .false.
        open (newunit=io_unit, file="tgv2d-ek.log", status="replace", form="formatted")
      else
        open (newunit=io_unit, file="tgv2d-ek.log", status="old", form="formatted", position="append")
      end if
      fmt = '(I0,1(1x,e12.4))'
      write (io_unit, fmt) t, ek_global
      close (io_unit)
    end if

  end subroutine calc_kinetic_energy

  !TODO: move this subroutine to more appropriate module
  !> Calculate enstrophy
  subroutine calc_enstrophy(par_env, mesh, t, u, v, w)

    use constants, only: ndim, ccs_string_len
    use types, only: field, ccs_mesh
    use vec, only: get_vector_data, restore_vector_data
    use parallel, only: allreduce, error_handling
    use parallel_types_mpi, only: parallel_environment_mpi
    use parallel_types, only: parallel_environment
    use mpi

    class(parallel_environment), allocatable, intent(in) :: par_env !< parallel environment
    type(ccs_mesh), intent(in) :: mesh !< the mesh
    integer(ccs_int), intent(in) :: t !< timestep
    class(field), intent(inout) :: u !< solve x velocity field
    class(field), intent(inout) :: v !< solve y velocity field
    class(field), intent(inout) :: w !< solve z velocity field

    real(ccs_real) :: ens_local, ens_global
    real(ccs_real), dimension(:), pointer :: dudy, dudz, dvdx, dvdz, dwdx, dwdy
    integer(ccs_int) :: index_p
    character(len=ccs_string_len) :: fmt
    integer(ccs_int) :: ierr

    logical, save :: first_time = .true.
    integer :: io_unit

    ens_local = 0.0_ccs_real
    ens_global = 0.0_ccs_real

    call get_vector_data(u%y_gradients, dudy)
    call get_vector_data(u%z_gradients, dudz)
    call get_vector_data(v%x_gradients, dvdx)
    call get_vector_data(v%z_gradients, dvdz)
    call get_vector_data(w%x_gradients, dwdx)
    call get_vector_data(w%y_gradients, dwdy)

    do index_p = 1, mesh%topo%local_num_cells

      ens_local = ens_local + (dwdy(index_p) - dvdz(index_p))**2 + &
                  (dudz(index_p) - dwdx(index_p))**2 + &
                  (dvdx(index_p) - dudy(index_p))**2

    end do
    ens_local = 0.5 * ens_local

    call restore_vector_data(u%y_gradients, dudy)
    call restore_vector_data(u%z_gradients, dudz)
    call restore_vector_data(v%x_gradients, dvdx)
    call restore_vector_data(v%z_gradients, dvdz)
    call restore_vector_data(w%x_gradients, dwdx)
    call restore_vector_data(w%y_gradients, dwdy)

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_AllReduce(ens_local, ens_global, 1, MPI_DOUBLE, MPI_SUM, par_env%comm, ierr)
      call error_handling(ierr, "mpi", par_env)
    class default
      call error_abort("ERROR: Unknown type")
    end select

    if (par_env%proc_id == par_env%root) then
      if (first_time) then
        first_time = .false.
        open (newunit=io_unit, file="tgv2d-ens.log", status="replace", form="formatted")
      else
        open (newunit=io_unit, file="tgv2d-ens.log", status="old", form="formatted", position="append")
      end if
      fmt = '(I0,1(1x,e12.4))'
      write (io_unit, fmt) t, ens_global
      close (io_unit)
    end if

  end subroutine calc_enstrophy

end module utils
