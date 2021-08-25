!> @brief Submodule file parallel_types_mpi.smod
!>
!> @details Implementation of parallel type set/get routines that
!> use MPI, e.g. MPI constants
submodule (parallel_types) parallel_types_mpi

  use mpi
  use iso_fortran_env

  implicit none

  contains

  !> @brief Set the reduction operators to their MPI values
  !>
  !> @param[out] reduction_operator_mpi rop - Variable that holds
  !> reduction operator values for MPI
  module subroutine set_reduction_operators(rop)

    class(reduction_operator), intent(out) :: rop

    select type (rop)
      type is (reduction_operator_mpi)   
      rop%sum_op = int(MPI_SUM,kind=4)
      rop%min_op = int(MPI_MIN,kind=4)
      rop%max_op = int(MPI_MAX,kind=4)
      rop%prod_op = int(MPI_PROD,kind=4)
      rop%land_op = int(MPI_LAND,kind=4)
      rop%lor_op = int(MPI_LOR,kind=4)
      rop%band_op = int(MPI_BAND,kind=4)
      rop%bor_op = int(MPI_BOR,kind=4)
      rop%maxloc_op = int(MPI_MAXLOC,kind=4)
      rop%minloc_op = int(MPI_MINLOC,kind=4)
    end select

  end subroutine

end submodule parallel_types_mpi