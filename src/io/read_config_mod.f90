!> @brief Module file read_config.mod
!>
!> @details Module defining interface to read YAML config file

module read_config

  use kinds, only : accs_real
    
  implicit none
  
  private  

  public :: get_case_name
  public :: get_steps
  public :: get_init
  public :: get_reference_numbers
  public :: get_solve
  public :: get_solvers
  public :: get_transient
  public :: get_target_residual
  public :: get_monitor_cell
  public :: get_convection_scheme
  public :: get_blending_factor
  public :: get_output_frequency
  public :: get_plot_format
  public :: get_output_type
  public :: get_boundaries

  interface

    !> @brief Get the name of the test case
    !
    !> @details Get the case name for the configuration file and 
    !! store it in a string.
    !
    !> @param[in] root - the YAML root node    
    !> @param[in,out] title - the case name string    
    module subroutine get_case_name(root, title)
      class(*), pointer, intent(in) :: root
      character(len=:), allocatable, intent(inout) :: title
    end subroutine
      
    !> @brief Get the number of steps
    !
    !> @details Get the maximum number of iterations 
    !! to be preformed in the current run 
    !
    !> @param[in] root - the YAML root node    
    !> @param[in,out] steps - the maximum number of iterations    
    module  subroutine get_steps(root, steps)
      class(*), pointer, intent(in) :: root
      integer, intent(inout) :: steps
    end subroutine
      
    !> @brief Get source of initial values
    !
    !> @details Get the source of the initial values - accepted
    !! values are "user", "field" or "step" 
    !
    !> @param[in] root - the YAML root node    
    !> @param[in,out] init - the source of the initial values    
    module  subroutine get_init(root, init)
      class(*), pointer, intent(in) :: root
      character(len=5), intent(inout) :: init
    end subroutine

    !> @brief Get reference numbers
    !
    !> @details Get the reference numbers, the fluid properties 
    !! and the operating condition 
    !
    !> @param[in] root - the YAML root node    
    !> @param[in,out] pressure - reference pressure      
    !> @param[in,out] temperature - reference temperature      
    !> @param[in,out] density - reference density      
    !> @param[in,out] viscosity - laminar viscosity      
    !> @param[in,out] pref_at_cell - cell at which the reference pressure is set      
    module subroutine get_reference_numbers(root, pressure, temperature, density, viscosity, pref_at_cell)
      class(*), pointer, intent(in) :: root
      integer, intent(inout) :: pressure
      integer, intent(inout) :: temperature
      integer, intent(inout) :: density
      integer, intent(inout) :: pref_at_cell
      real(accs_real), intent(inout) :: viscosity
    end subroutine

    !> @brief Get variables to be solved
    !
    !> @details By default, all variables will be solved. Using this 
    !! "solve" keyword, the user can specifically request that 
    !! certain variables will not be solved by setting in to "off"
    !
    !> @param[in] root - the YAML root node    
    !> @param[in,out] u_sol - solve u on/off
    !> @param[in,out] v_sol - solve v on/off
    !> @param[in,out] w_sol - solve w on/off
    !> @param[in,out] p_sol - solve p on/off
    !
    !> @todo extend list of variables 
    module subroutine get_solve(root, u_sol, v_sol, w_sol, p_sol)
      class(*), pointer, intent(in) :: root
      character(len=3), intent(inout) :: u_sol
      character(len=3), intent(inout) :: v_sol
      character(len=3), intent(inout) :: w_sol
      character(len=3), intent(inout) :: p_sol
    end subroutine

    !> @brief Get solvers to be used
    !
    !> @details Get the solvers that are to be used for each of
    !! the variables. Solver types are defined by integer values
    !
    !> @param[in] root - the YAML root node    
    !> @param[in,out] u_solver - solver to be used for u
    !> @param[in,out] v_solver - solver to be used for v
    !> @param[in,out] w_solver - solver to be used for w
    !> @param[in,out] p_solver - solver to be used for p
    !
    !> @todo extend list of variables   
    module subroutine get_solvers(root, u_solver, v_solver, w_solver, p_solver)
      class(*), pointer, intent(in) :: root
      integer, intent(inout) :: u_solver
      integer, intent(inout) :: v_solver
      integer, intent(inout) :: w_solver
      integer, intent(inout) :: p_solver
    end subroutine

    !> @brief Get solvers to be used
    !
    !> @details Get the solvers that are to be used for each of
    !! the variables. Solver types are defined by integer values
    !
    !> @param[in] root - the YAML root node    
    !> @param[in,out] u_solver - solver to be used for u
    !> @param[in,out] v_solver - solver to be used for v
    !> @param[in,out] w_solver - solver to be used for w
    !> @param[in,out] p_solver - solver to be used for p
    !
    !> @todo extend list of variables   
    module subroutine get_transient(root, transient_type, dt, gamma, max_sub_steps)
      class(*), pointer, intent(in) :: root
      character(len=5), intent(inout) :: transient_type
      real(accs_real), intent(inout) :: dt
      real(accs_real), intent(inout) :: gamma
      integer, intent(inout) :: max_sub_steps
    end subroutine

    module  subroutine get_target_residual(root, residual)
      class(*), pointer, intent(in) :: root
      real(accs_real), intent(inout) :: residual
    end subroutine

    module  subroutine get_monitor_cell(root, monitor_cell)
      class(*), pointer, intent(in) :: root
      integer, intent(inout) :: monitor_cell
    end subroutine

    module  subroutine get_convection_scheme(root, u_conv, v_conv, w_conv)
      class(*), pointer, intent(in) :: root
      integer, intent(inout) :: u_conv
      integer, intent(inout) :: v_conv
      integer, intent(inout) :: w_conv
    end subroutine

    module subroutine get_blending_factor(root, u_blend, v_blend, p_blend)
      class(*), pointer, intent(in) :: root
      real(accs_real), intent(inout) :: u_blend
      real(accs_real), intent(inout) :: v_blend
      real(accs_real), intent(inout) :: p_blend
    end subroutine

    module subroutine get_output_frequency(root, output_freq, output_iter)
      class(*), pointer, intent(in) :: root
      integer, intent(inout) :: output_freq
      integer, intent(inout) :: output_iter
    end subroutine

    module subroutine get_plot_format(root, plot_format)
      class(*), pointer, intent(in) :: root
      character(len=:), allocatable, intent(inout) :: plot_format
    end subroutine

    module subroutine get_output_type(root, post_type, post_vars)
      class(*), pointer, intent(in) :: root
      character(len=:), allocatable, intent(inout) :: post_type
      character(len=2), dimension(10), intent(inout) :: post_vars    
    end subroutine

    module subroutine get_boundaries(root, bnd_region, bnd_type, bnd_vector)
      class(*), pointer, intent(in) :: root
      character(len=16), dimension(:), allocatable, intent(inout) :: bnd_region
      character(len=16), dimension(:), allocatable, intent(inout) :: bnd_type
      real(accs_real), dimension(:,:), allocatable, intent(inout) :: bnd_vector
    end subroutine

  end interface
end module read_config
  