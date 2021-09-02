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
  module subroutine set_mpi_reduction_operator(this)

    class(parallel_environment), intent(inout) :: this

    select type (this)

    type is (parallel_environment_mpi)   
      this%sum_op%op = MPI_SUM
      this%min_op%op = MPI_MIN
      this%max_op%op = MPI_MAX
      this%prod_op%op = MPI_PROD
      this%land_op%op = MPI_LAND
      this%lor_op%op = MPI_LOR
      this%band_op%op = MPI_BAND
      this%bor_op%op = MPI_BOR
      this%maxloc_op%op = MPI_MAXLOC
      this%minloc_op%op = MPI_MINLOC

    class default
      write(*,*) "Unsupported parallel environment"    

    end select

  end subroutine 

end submodule parallel_types_mpi