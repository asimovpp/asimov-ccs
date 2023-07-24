!> Testing library
module testing_lib
#include "ccs_macros.inc"

  use MPI

  use kinds
  use types
  use parallel
  use parallel_types
  use parallel_types_mpi
  use utils, only: str, exit_print
  use constants, only: ccs_split_type_low_high, ccs_split_type_shared

  implicit none

  public :: assert_eq, assert_neq, assert_lt, assert_gt, assert_bool, assert_le, assert_ge

  !> Assert equality
  interface assert_eq
    procedure assert_eq_integer_rank0
    procedure assert_eq_integer_rank1
    procedure assert_eq_real_rank0
    procedure assert_eq_real_rank1
    procedure assert_eq_string
  end interface

  !> Assert first argument is less than the second argument
  interface assert_lt
    procedure assert_lt_integer
    procedure assert_lt_real
  end interface

  !> Assert first argument is less than or equal to the second argument
  interface assert_le
    procedure assert_le_integer
    procedure assert_le_real
  end interface

  !> Assert first argument is greater than the second argument
  interface assert_gt
    procedure assert_gt_integer
    procedure assert_gt_real
  end interface

  !> Assert first argument is greater than or equal to the second argument
  interface assert_ge
    procedure assert_ge_integer
    procedure assert_ge_real
  end interface

  !> Assert that input is True
  interface assert_bool
    procedure assert_bool_rank0
    procedure assert_bool_rank1
  end interface

  !> Assert inequality
  interface assert_neq
    procedure assert_neq_integer
    procedure assert_neq_real
    procedure assert_neq_string
  end interface

  !> Core procedure for comparing numbers
  interface a_eq
    procedure a_eq_integer
    procedure a_eq_real
  end interface

  !> Printing utilities for multi-dimensional test results
  interface print_failed
    procedure print_failed_integer
    procedure print_failed_real
    procedure print_failed_bool
  end interface

  class(parallel_environment), allocatable, target :: par_env
  class(parallel_environment), allocatable, target :: shared_env
  class(parallel_environment), allocatable, target :: roots_env
  integer(ccs_err) :: ierr
  integer :: real_type
  character(1024) :: message

  real(ccs_real), parameter :: eps = epsilon(0.0_ccs_real)

contains

  !v Test initialisation
  !
  !  Performs initialisation for the test (setting up parallel environment, etc.)
  subroutine init()

    integer, allocatable :: seed(:)
    integer :: n
    logical :: use_mpi_splitting

    if (kind(0.0_ccs_real) == kind(0.0d0)) then
      real_type = MPI_DOUBLE
    else
      real_type = MPI_FLOAT
    end if

    call initialise_parallel_environment(par_env)
    use_mpi_splitting = .false.
    call create_new_par_env(par_env, ccs_split_type_low_high, use_mpi_splitting, shared_env)
    !call create_new_par_env(par_env, ccs_split_type_shared, use_mpi_splitting, shared_env)
    call create_shared_roots_comm(par_env, shared_env, roots_env)

    ! XXX: This would be a good candidate for a testing library
    call random_seed(size=n)
    allocate (seed(n))
    call random_seed(get=seed)
    if (par_env%proc_id == par_env%root) then
      print *, "Using seed: ", seed
      print *, "----------------------------------"
    end if
    deallocate (seed)

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Barrier(par_env%comm, ierr)
    class default
      call error_abort("ERROR: Unknown parallel environment!")
    end select

  end subroutine init

  !v Test finalisation
  !
  !  Performs finalisation for the test (tearing down parallel environment, etc.)
  subroutine fin()

    call cleanup_parallel_environment(par_env)

  end subroutine fin

  !v Helper function to get a random number in parallel
  !
  !  Generates a random number and broadcasts the value on the root of the parallel
  !  environment, ensuring a uniform value is used.
  !
  !  @note Does this belong in the parallel module?
  real(ccs_real) function parallel_random(par_env)

    class(parallel_environment), intent(in) :: par_env !< parallel environment

    call random_number(parallel_random)

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Bcast(parallel_random, 1, real_type, par_env%root, par_env%comm, ierr)
    class default
      call error_abort("ERROR: Unknown parallel environment!")
    end select

  end function parallel_random

  !v Test failure stop
  !
  !  Stop a test, provide an error message, do cleanup etc.
  subroutine stop_test(message)

    character(*), intent(in) :: message !< Error message

    print *, "(rank "   //   str(par_env%proc_id)   //   ") ", trim(message)

    ! other PEs might not have encountered a test failure
    ! fin()

    stop 1

  end subroutine stop_test

  !v Stop test or return test result
  !
  !  If outval is present, this simply returns the value of res.
  !  If outval is not present, the test will be stopped if res is False and message will be printed.
  subroutine return_or_stop(res, message, outval)

    logical, intent(in) :: res               !< Evaluation result
    character(*), intent(in) :: message      !< Error message
    logical, optional, intent(out) :: outval !< Output value to replace stopping the test

    if (present(outval)) then
      outval = res
    else
      if (.not. res) then
        call stop_test(message)
      end if
    end if

  end subroutine return_or_stop

!==========================Numeric comparison operators
  !> Compare two integer values for equality
  elemental logical function a_eq_integer(a, b) result(comparison)

    integer(ccs_int), intent(in) :: a
    integer(ccs_int), intent(in) :: b

    comparison = (a == b)

  end function a_eq_integer

  !> Compare two real values for equality within a tolerance
  elemental logical function a_eq_real(a, b) result(comparison)

    real(ccs_real), intent(in) :: a
    real(ccs_real), intent(in) :: b

    ! TODO: double check we are happy with this evaluation
    comparison = (abs(a - b) <= 10 * epsilon(b) * abs(b))

  end function a_eq_real
!==========================

!==========================Integer equality
  !> Single integer comparison
  subroutine assert_eq_integer_rank0(received, expected, message, outval)

    integer(ccs_int), intent(in) :: received !< Test value
    integer(ccs_int), intent(in) :: expected !< Reference value
    character(*), intent(in) :: message      !< Error message
    logical, optional, intent(out) :: outval !< Output value to replace stopping the test

    call return_or_stop(a_eq(received, expected), &
                        message   //   " Expected: "   //   str(expected)   //   " Received: "   //   str(received), &
                        outval)

  end subroutine assert_eq_integer_rank0

  !> Integer array comparison
  subroutine assert_eq_integer_rank1(received, expected, message, outval)

    integer(ccs_int), dimension(:), intent(in) :: received !< Test values
    integer(ccs_int), dimension(:), intent(in) :: expected !< Reference values
    character(*), intent(in) :: message                    !< Error message
    logical, optional, intent(out) :: outval               !< Output value to replace stopping the test

    call return_or_stop(all(a_eq(received, expected)), &
                        message   //   print_failed(received, expected), &
                        outval)

  end subroutine assert_eq_integer_rank1
!==========================

!==========================Real equality
  !> Single real comparison
  subroutine assert_eq_real_rank0(received, expected, message, outval)

    real(ccs_real), intent(in) :: received   !< Test value
    real(ccs_real), intent(in) :: expected   !< Reference value
    character(*), intent(in) :: message      !< Error message
    logical, optional, intent(out) :: outval !< Output value to replace stopping the test

    call return_or_stop(a_eq(received, expected), &
                        message   //   " Expected: "   //   str(expected)   //   " Received: "   //   str(received), &
                        outval)

  end subroutine assert_eq_real_rank0

  !> Real array comparison
  subroutine assert_eq_real_rank1(received, expected, message, outval)

    real(ccs_real), dimension(:), intent(in) :: received !< Test values
    real(ccs_real), dimension(:), intent(in) :: expected !< Reference values
    character(*), intent(in) :: message                  !< Error message
    logical, optional, intent(out) :: outval             !< Output value to replace stopping the test

    call return_or_stop(all(a_eq(received, expected)), &
                        message   //   print_failed(received, expected), &
                        outval)

  end subroutine assert_eq_real_rank1
!==========================

!==========================String equality
  !> String comparison
  subroutine assert_eq_string(received, expected, message, outval)

    character(*), intent(in) :: received     !< Test value
    character(*), intent(in) :: expected     !< Reference value
    character(*), intent(in) :: message      !< Error message
    logical, optional, intent(out) :: outval !< Output value to replace stopping the test

    call return_or_stop(received == expected, &
                        message   //   " Expected: "   //   expected   //   " Received: "   //   received, &
                        outval)

  end subroutine assert_eq_string
!==========================

!==========================Inequality
  !> Integer comparison
  subroutine assert_neq_integer(received, notexpected, message, outval)

    integer(ccs_int), intent(in) :: received    !< Test value
    integer(ccs_int), intent(in) :: notexpected !< Reference value
    character(*), intent(in) :: message         !< Error message
    logical, optional, intent(out) :: outval    !< Output value to replace stopping the test

    call return_or_stop(.not. a_eq(received, notexpected), &
                        message   //   " Not Expected: "   //   str(notexpected)   //   " Received: "   //   str(received), &
                        outval)

  end subroutine assert_neq_integer

  !> Real comparison
  subroutine assert_neq_real(received, notexpected, message, outval)

    real(ccs_real), intent(in) :: received    !< Test value
    real(ccs_real), intent(in) :: notexpected !< Reference value
    character(*), intent(in) :: message       !< Error message
    logical, optional, intent(out) :: outval  !< Output value to replace stopping the test

    call return_or_stop(.not. a_eq(received, notexpected), &
                        message   //   " Not Expected: "   //   str(notexpected)   //   " Received: "   //   str(received), &
                        outval)

  end subroutine assert_neq_real

  !> String comparison
  subroutine assert_neq_string(received, notexpected, message, outval)

    character(*), intent(in) :: received     !< Test value
    character(*), intent(in) :: notexpected  !< Reference value
    character(*), intent(in) :: message      !< Error message
    logical, optional, intent(out) :: outval !< Output value to replace stopping the test

    call return_or_stop(.not. received == notexpected, &
                        message   //   " Not Expected: "   //   notexpected   //   " Received: "   //   received, &
                        outval)

  end subroutine assert_neq_string
!==========================

!==========================Less-than
  !> Integer comparison
  subroutine assert_lt_integer(received, upper_limit, message, outval)

    integer(ccs_int), intent(in) :: received    !< Test value
    integer(ccs_int), intent(in) :: upper_limit !< Reference value
    character(*), intent(in) :: message         !< Error message
    logical, optional, intent(out) :: outval    !< Output value to replace stopping the test

    call return_or_stop(received < upper_limit, &
                        message   //   "Upper limit allowed: "   //   str(upper_limit)   //   " Received: "   //   str(received), &
                        outval)

  end subroutine assert_lt_integer

  !> Real comparison
  subroutine assert_lt_real(received, upper_limit, message, outval)

    real(ccs_real), intent(in) :: received    !< Test value
    real(ccs_real), intent(in) :: upper_limit !< Reference value
    character(*), intent(in) :: message       !< Error message
    logical, optional, intent(out) :: outval  !< Output value to replace stopping the test

    call return_or_stop(received < upper_limit, &
                        message   //   "Upper limit allowed: "   //   str(upper_limit)   //   " Received: "   //   str(received), &
                        outval)

  end subroutine assert_lt_real
!==========================

!==========================Less-than-or-equal
  !> Integer comparison
  subroutine assert_le_integer(received, upper_limit, message, outval)

    integer(ccs_int), intent(in) :: received    !< Test value
    integer(ccs_int), intent(in) :: upper_limit !< Reference value
    character(*), intent(in) :: message         !< Error message
    logical, optional, intent(out) :: outval    !< Output value to replace stopping the test

    call return_or_stop(received <= upper_limit, &
                        message   //   "Upper limit allowed: "   //   str(upper_limit)   //   " Received: "   //   str(received), &
                        outval)

  end subroutine assert_le_integer

  !> Real comparison
  subroutine assert_le_real(received, upper_limit, message, outval)

    real(ccs_real), intent(in) :: received    !< Test value
    real(ccs_real), intent(in) :: upper_limit !< Reference value
    character(*), intent(in) :: message       !< Error message
    logical, optional, intent(out) :: outval  !< Output value to replace stopping the test

    call return_or_stop(received <= upper_limit, &
                        message   //   "Upper limit allowed: "   //   str(upper_limit)   //   " Received: "   //   str(received), &
                        outval)

  end subroutine assert_le_real
!==========================

!==========================Greater-than
  !> Integer comparison
  subroutine assert_gt_integer(received, lower_limit, message, outval)

    integer(ccs_int), intent(in) :: received    !< Test value
    integer(ccs_int), intent(in) :: lower_limit !< Reference value
    character(*), intent(in) :: message         !< Error message
    logical, optional, intent(out) :: outval    !< Output value to replace stopping the test

    call return_or_stop(received > lower_limit, &
                        message   //   "Lower limit allowed: "   //   str(lower_limit)   //   " Received: "   //   str(received), &
                        outval)

  end subroutine assert_gt_integer

  !> Real comparison
  subroutine assert_gt_real(received, lower_limit, message, outval)

    real(ccs_real), intent(in) :: received    !< Test value
    real(ccs_real), intent(in) :: lower_limit !< Reference value
    character(*), intent(in) :: message       !< Error message
    logical, optional, intent(out) :: outval  !< Output value to replace stopping the test

    call return_or_stop(received > lower_limit, &
                        message   //   "Lower limit allowed: "   //   str(lower_limit)   //   " Received: "   //   str(received), &
                        outval)

  end subroutine assert_gt_real
!==========================

!==========================Greater-than-or-equal
  !> Integer comparison
  subroutine assert_ge_integer(received, lower_limit, message, outval)

    integer(ccs_int), intent(in) :: received    !< Test value
    integer(ccs_int), intent(in) :: lower_limit !< Reference value
    character(*), intent(in) :: message         !< Error message
    logical, optional, intent(out) :: outval    !< Output value to replace stopping the test

    call return_or_stop(received >= lower_limit, &
                        message   //   "Lower limit allowed: "   //   str(lower_limit)   //   " Received: "   //   str(received), &
                        outval)

  end subroutine assert_ge_integer

  !> Real comparison
  subroutine assert_ge_real(received, lower_limit, message, outval)

    real(ccs_real), intent(in) :: received    !< Test value
    real(ccs_real), intent(in) :: lower_limit !< Reference value
    character(*), intent(in) :: message       !< Error message
    logical, optional, intent(out) :: outval  !< Output value to replace stopping the test

    call return_or_stop(received >= lower_limit, &
                        message   //   "Lower limit allowed: "   //   str(lower_limit)   //   " Received: "   //   str(received), &
                        outval)

  end subroutine assert_ge_real
!==========================

!==========================Bool tests
  !> Single bool test
  subroutine assert_bool_rank0(received, message, outval)

    logical, intent(in) :: received          !< Test value
    character(*), intent(in) :: message      !< Error message
    logical, optional, intent(out) :: outval !< Output value to replace stopping the test

    call return_or_stop(received, &
                        message   //   " Expected: T (true) Received: "   //   str(received), &
                        outval)

  end subroutine assert_bool_rank0

  !> Bool array test
  subroutine assert_bool_rank1(received, message, outval)

    logical, dimension(:), intent(in) :: received !< Test values
    character(*), intent(in) :: message           !< Error message
    logical, optional, intent(out) :: outval      !< Output value to replace stopping the test

    call return_or_stop(all(received), &
                        message   //   print_failed(received), &
                        outval)

  end subroutine assert_bool_rank1
!==========================

!==========================Printing functions
  !> Print integers
  function print_failed_integer(received, expected) result(msg)

    integer(ccs_int), dimension(:), intent(in) :: expected !< Test values
    integer(ccs_int), dimension(:), intent(in) :: received !< Reference values
    character(len=:), allocatable :: msg                   !< Constructed message

    integer :: i
    logical, dimension(:), allocatable :: mask
    allocate (mask(size(received)))
    mask = a_eq(received, expected)

    msg = new_line('a')   //   "Index Expected Received"   //   new_line('a')
    do i = 1, size(mask)
      if (.not. mask(i)) then
        msg = msg   //   str(i)   //   achar(9)   //   str(expected(i))   //   achar(9)   //   str(received(i))   //   new_line('a')
      end if
    end do

  end function print_failed_integer

  !> Print reals
  function print_failed_real(received, expected) result(msg)

    real(ccs_real), dimension(:), intent(in) :: expected !< Test values
    real(ccs_real), dimension(:), intent(in) :: received !< Reference values
    character(len=:), allocatable :: msg                 !< Constructed message

    integer :: i
    logical, dimension(:), allocatable :: mask
    allocate (mask(size(received)))
    mask = a_eq(received, expected)

    msg = new_line('a')   //   "Index Expected Received"   //   new_line('a')
    do i = 1, size(mask)
      if (.not. mask(i)) then
        msg = msg   //   str(i)   //   achar(9)   //   str(expected(i))   //   achar(9)   //   str(received(i))   //   new_line('a')
      end if
    end do

  end function print_failed_real

  !> Print bools
  function print_failed_bool(received) result(msg)

    logical, dimension(:), intent(in) :: received !< Test values
    character(len=:), allocatable :: msg          !< Constructed message
    integer :: i

    msg = new_line('a')   //   "Index Received"   //   new_line('a')
    do i = 1, size(received)
      if (.not. received(i)) then
        msg = msg   //   str(i)   //   achar(9)   //   "FALSE"   //   new_line('a')
      end if
    end do

  end function print_failed_bool
!==========================

end module testing_lib
