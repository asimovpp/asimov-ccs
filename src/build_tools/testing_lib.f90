!> @brief Testing library
module testing_lib

  use MPI

  use kinds
  use types
  use parallel
  use parallel_types
  use parallel_types_mpi

  implicit none

  class(parallel_environment), allocatable, target :: par_env
  integer(accs_err) :: ierr
  integer :: real_type
  character(1024) :: message

  real(accs_real), parameter :: eps = epsilon(0.0_accs_real)
  
contains

  !> @brief Test initialisation
  !
  !> @description Performs initialisation for the test (setting up parallel environment, etc.)
  subroutine init()

    integer, allocatable :: seed(:)
    integer :: n

    if (kind(0.0_accs_real) == kind(0.0d0)) then
      real_type = MPI_DOUBLE
    else
      real_type = MPI_FLOAT
    end if

    call initialise_parallel_environment(par_env)

    ! XXX: This would be a good candidate for a testing library
    call random_seed(size=n)
    allocate(seed(n))
    call random_seed(get=seed)
    if (par_env%proc_id == par_env%root) then
      print *, "Using seed: ", seed
      print *, "----------------------------------"
    end if
    deallocate(seed)

    select type(par_env)
    type is(parallel_environment_mpi)
      call MPI_Barrier(par_env%comm, ierr)
    class default
      print *, "ERROR: Unknown parallel environment!"
      stop
    end select
    
  end subroutine init

  !> @brief Test finalisation
  !
  !> @description Performs finalisation for the test (tearing down parallel environment, etc.)
  subroutine fin()

    call cleanup_parallel_environment(par_env)

  end subroutine fin

  !> @brief Helper function to get a random number in parallel
  !
  !> @description Generates a random number and broadcasts the value on the root of the parallel
  !! environment, ensuring a uniform value is used.
  !
  !> @note Does this belong in the parallel module?
  real(accs_real) function parallel_random(par_env)

    class(parallel_environment), intent(in) :: par_env

    call random_number(parallel_random)

    select type(par_env)
    type is(parallel_environment_mpi)
      call MPI_Bcast(parallel_random, 1, real_type, par_env%root, par_env%comm, ierr)
    class default
      print *, "ERROR: Unknown parallel environment!"
      stop
    end select
    
  end function parallel_random

  !> @brief Test failure stop
  !
  !> @description Stop a test, provide an error message, do cleanup etc.
  subroutine stop_test(message)

    character(*), intent(in) :: message
    character(len=32) :: id_str

    write (id_str, "(I0)") par_env%proc_id
    print *, "("//trim(id_str)//") ", trim(message)

    ! other PEs might not have encountered a test failure
    ! fin()

    stop 1
  end subroutine stop_test
  
  !> @brief Assertion for integer equality
  !
  !> @description Check whether input integers are equal. If not, construct message, print and stop.
  subroutine assert_equal(a, b, msg_format)

    integer(accs_int), intent(in) :: a
    integer(accs_int), intent(in) :: b
    character(*), intent(in) :: msg_format
    character(1024) :: message

    if (a /= b) then
      write (message, msg_format) a, b
      call stop_test(message)
    end if

  end subroutine assert_equal
  
  
end module testing_lib
