!v Module file parallel_types.mod
!
!  Module that defines the parallel environment types for ASiMoV-CCS
!
!  @build mpi
module parallel_types

  use kinds, only: ccs_int

  implicit none

  private

  type, public :: reduction_operator !< placeholder reduction operator type
  end type reduction_operator

  !v parallel environment type with common parameters
  type, public :: parallel_environment
    integer(ccs_int) :: proc_id = huge(0_ccs_int)   !< process id
    integer(ccs_int) :: num_procs = huge(0_ccs_int) !< number of processes
    integer(ccs_int) :: root = huge(0_ccs_int)      !< root process
  end type parallel_environment

end module parallel_types
