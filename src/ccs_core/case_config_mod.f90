!v Module file case_config.mod
!
!  Provides concrete types and bases of extensible types.

module case_config

  use kinds, only: ccs_int, ccs_real
  use parallel_types, only: parallel_environment

  implicit none

  private

  ! Restart a simulation
  logical, save, public :: restart = .false.
  logical, save, public :: unsteady = .false.

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

  ! Fluid properties
  real(ccs_real), public :: viscocity = huge(0.0)

  character(len=:), allocatable, save, public :: case_name
  character(len=:), allocatable, save, public :: velocity_solver_method_name
  character(len=:), allocatable, save, public :: velocity_solver_precon_name

  character(len=:), allocatable, save, public :: pressure_solver_method_name
  character(len=:), allocatable, save, public :: pressure_solver_precon_name

  logical, save, public :: write_gradients = .false.

  ! Logical to toggle whether the partition quality and matrix bandwidth is computed or not
  logical, save, public :: compute_partqual = .true.
  logical, save, public :: compute_bwidth = .true.

end module case_config
