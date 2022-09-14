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
    integer(ccs_int) :: proc_id !< process id
    integer(ccs_int) :: num_procs !< number of processes
    integer(ccs_int) :: root !< root process
  end type parallel_environment

end module parallel_types
