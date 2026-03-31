!v Logging utility to specify where to redirect output.
module logging
#include "ccs_macros.inc"

  use iso_fortran_env, only : output_unit, error_unit

  implicit none

  logical :: initialised = .false.

  integer :: log_unit_out = output_unit
  integer :: log_unit_err = error_unit

  character(256) :: message_buf

  contains

  !> Initialise logging module
  subroutine initialise_logging(par_env, file_name)
    use parallel_types, only: parallel_environment
    use parallel, only: is_root

    implicit none

    class(parallel_environment), allocatable, intent(in) :: par_env !< parallel environment
    character(len=*), intent(in) :: file_name
    integer :: ierr

    if (is_root(par_env)) then
      open(newunit=log_unit_out, file=file_name, status="replace", iostat=ierr)
      if (ierr /= 0) then
        write(log_unit_out,*) "ERROR: Could not open log file " // file_name
        stop 1
      end if
      initialised = .true.
    end if

  end subroutine initialise_logging

  !> Finalise logging module
  subroutine finalise_logging(par_env)
    use parallel_types, only: parallel_environment
    use parallel, only: is_root

    implicit none

    class(parallel_environment), allocatable, intent(in) :: par_env !< parallel environment
    integer :: ierr

    if (is_root(par_env)) then
      close(log_unit_out, iostat=ierr)
      if (ierr /= 0) then
        write(log_unit_out, *) "ERROR: Could not close log file."
      end if
      initialised = .false.
    end if

  end subroutine finalise_logging

end module logging
