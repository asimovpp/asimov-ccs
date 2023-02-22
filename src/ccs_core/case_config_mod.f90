!v Module file case_config.mod
!
!  Provides concrete types and bases of extensible types.

module case_config

  use kinds, only: ccs_int, ccs_real
  use parallel_types, only: parallel_environment

  implicit none

  private

  ! Number of timesteps and iterations
  integer(ccs_int), public :: num_steps = huge(0)
  integer(ccs_int), public :: num_iters = huge(0)

  ! Number of cells per side (used for mesh generation)
  integer(ccs_int), public :: cps = huge(0)

  ! Frequency (in terms of timesteps) of writing solution to file
  integer(ccs_int), public :: write_frequency = huge(0)

  real(ccs_real), public :: velocity_relax = huge(0.0)
  real(ccs_real), public :: pressure_relax = huge(0.0)

  real(ccs_real), public :: domain_size = huge(0.0)
  real(ccs_real), public :: cfl = huge(0.0)
  real(ccs_real), public :: dt = huge(0.0)

  real(ccs_real), public :: res_target = huge(0.0)

  character(len=:), allocatable, save, public :: case_name
  character(len=:), allocatable, save, public :: velocity_solver_method_name
  character(len=:), allocatable, save, public :: velocity_solver_precon_name

  character(len=:), allocatable, save, public :: pressure_solver_method_name
  character(len=:), allocatable, save, public :: pressure_solver_precon_name

  logical, public :: write_gradients = .false.

end module case_config
