!> @brief Submodule file parallel_types_mpi.smod
!>
!> @details Implementation of parallel type set/get routines that
!> use MPI, e.g. MPI constants
submodule (parallel_types) parallel_types_mpi

  implicit none

  contains

  !> @brief Set the reduction operators to their MPI values
  !>
  !> @param[out] reduction_operator_mpi rop - Variable that holds
  !> reduction operator values for MPI
  module subroutine set_reduction_operators(rop)

    class(reduction_operator), intent(inout) :: rop

    select type (rop)

    type is (reduction_operator_mpi)   
      rop%sum_op = MPI_SUM
      rop%min_op = MPI_MIN
      rop%max_op = MPI_MAX
      rop%prod_op = MPI_PROD
      rop%land_op = MPI_LAND
      rop%lor_op = MPI_LOR
      rop%band_op = MPI_BAND
      rop%bor_op = MPI_BOR
      rop%maxloc_op = MPI_MAXLOC
      rop%minloc_op = MPI_MINLOC

    class default
      write(*,*) "Unsupported parallel environment"    

    end select

  end subroutine 

end submodule parallel_types_mpi