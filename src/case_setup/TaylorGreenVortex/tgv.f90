
program tgv

  use, intrinsic :: iso_fortran_env, only:  output_unit
  use yaml, only:       parse, error_length
  use read_config
  use kinds, only : accs_real

  implicit none

  class(*), pointer :: root
  character(len=error_length) :: error

  ! Case title
  character(len=:), allocatable :: title

  ! Number of steps
  integer :: steps

  ! Initialisation
  character(len=:), allocatable :: init

  ! Reference numbers
  integer :: pref_at_cell
  real(accs_real) :: pressure
  real(accs_real) :: temperature
  real(accs_real) :: density
  real(accs_real) :: viscosity

  ! Solve
  character(len=:), allocatable :: w_sol

  ! Solvers
  integer :: u_solver
  integer :: v_solver

  ! Unsteady solution
  character(len=:), allocatable :: transient_type
  real(accs_real) :: dt
  real(accs_real) :: gamma
  integer :: max_sub_steps

  ! Target residual
  real(accs_real) :: residual

  ! Monitor cell
  integer :: monitor_cell

  ! Convection schemes
  integer :: u_conv
  integer :: v_conv
  integer :: w_conv

  ! Blending factors
  real(accs_real) :: u_blend
  real(accs_real) :: v_blend
  real(accs_real) :: w_blend

  ! Relaxation factors
  real(accs_real) :: u_relax
  real(accs_real) :: v_relax
  real(accs_real) :: p_relax

  ! Output frequency & iteration
  integer :: output_freq
  integer :: output_iter

  ! Plot format
  character(len=:), allocatable :: plot_format

  ! Post type & vars
  character(len=:), allocatable :: post_type
  character(len=2), dimension(10) :: post_vars = "  "

  ! Boundardies
  character(len=16), dimension(:), allocatable :: bnd_region
  character(len=16), dimension(:), allocatable :: bnd_type
  real(accs_real), dimension(:,:), allocatable :: bnd_vector

  ! Read TGV configuration file
  call read_configuration()

  contains

  subroutine read_configuration()

    root => parse("./tgv_config.yaml", error=error)
    if (error/='') then
      print*,trim(error)
      stop 1
    endif
    
    ! Get title
    call get_case_name(root, title)
    
    ! Get reference numbers
    call get_reference_number(root, p_ref=pressure, temp_ref=temperature, dens_ref=density, &
                               visc_ref=viscosity, pref_at_cell=pref_at_cell)

    ! Get steps
    call get_steps(root, steps)

    ! Variables to solve
    call get_solve(root, w_sol=w_sol)

    ! Solvers
    call get_solver(root, u_solver=u_solver, v_solver=v_solver)  
    
    ! Get initilisation
    call get_init(root, init)

    ! Get unsteady solution parameters
    call get_transient(root, transient_type, dt, gamma, max_sub_steps)

    ! Get target resdiual
    call get_target_residual(root, residual)

    ! Get monitor cell ID
    call get_monitor_cell(root, monitor_cell)

    ! Get convection schemes
    call get_convection_scheme(root, u_conv, v_conv, w_conv)

    ! Get blending factors
    call get_blending_factor(root, u_blend, v_blend, w_blend)

    ! Get relaxation factors
    call get_relaxation_factor(root, u_relax=u_relax, v_relax=v_relax, p_relax=p_relax)

    ! Get output frequency and iteration
    call get_output_frequency(root, output_freq, output_iter)

    ! Get plot format
    call get_plot_format(root, plot_format)

    ! Get output type and variables
    call get_output_type(root, post_type, post_vars)

    ! Get boundaries
    call get_boundaries(root, bnd_region, bnd_type, bnd_vector)
     
    deallocate(root)

  end subroutine

end program tgv
