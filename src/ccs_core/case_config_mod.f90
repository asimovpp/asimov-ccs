!>  Module file case_config.mod
!
!> @details Provides concrete types and bases of extensible types.

module case_config

    use kinds, only : ccs_int, ccs_real
    use parallel_types, only: parallel_environment
    
    implicit none
  
    private

    integer(ccs_int), public :: num_steps = huge(0)

    real(ccs_real), public :: velocity_relax = huge(0.0)
    real(ccs_real), public :: pressure_relax = huge(0.0)

end module case_config