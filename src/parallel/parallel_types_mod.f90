!>  Module file parallel_types.mod
!
!> @build mpi
!
!>  Module that defines the parallel environment types for ASiMoV-CCS
module parallel_types

  use kinds, only: ccs_int
  
  implicit none

  private 

  !>  placeholder reduction operator type
  type, public :: reduction_operator
  end type reduction_operator

  !v  parallel environment type with common parameters
  !        process id, number of processes and root process
  type, public :: parallel_environment
    integer(ccs_int) :: proc_id
    integer(ccs_int) :: num_procs
    integer(ccs_int) :: root
  end type parallel_environment

end module parallel_types
