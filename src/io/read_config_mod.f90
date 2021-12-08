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
  public :: get_reference_number
  public :: get_solve
  public :: get_solver
  public :: get_transient
  public :: get_target_residual
  public :: get_monitor_cell
  public :: get_convection_scheme
  public :: get_blending_factor
  public :: get_relaxation_factor
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
    !> @param[in] root - the entry point to the config file    
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
    !> @param[in] root - the entry point to the config file    
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
    !> @param[in] root - the entry point to the config file    
    !> @param[in,out] init - the source of the initial values    
    module  subroutine get_init(root, init, u_init, v_init, w_init, te_init, ed_init)
      class(*), pointer, intent(in) :: root
      character(len=:), allocatable, intent(inout) :: init
      integer, optional, intent(inout) :: u_init
      integer, optional, intent(inout) :: v_init
      integer, optional, intent(inout) :: w_init
      integer, optional, intent(inout) :: te_init
      integer, optional, intent(inout) :: ed_init

    end subroutine

    !> @brief Get reference numbers
    !
    !> @details Get the reference numbers, the fluid properties 
    !! and the operating condition 
    !
    !> @param[in] root - the entry point to the config file    
    !> @param[in,out] p_ref - reference pressure 
    !> @param[in,out] p_total - total pressure 
    !> @param[in,out] temp_ref - reference temperature      
    !> @param[in,out] dens_ref - reference density      
    !> @param[in,out] visc_ref - laminar viscosity      
    !> @param[in,out] velo_ref - reference velocity      
    !> @param[in,out] len_ref - reference length, used to define the Reynolds number of the flow      
    !> @param[in,out] pref_at_cell - cell at which the reference pressure is set      
    module subroutine get_reference_number(root, p_ref, p_total, temp_ref, &
                                            dens_ref, visc_ref, velo_ref, len_ref, pref_at_cell)
      class(*), pointer, intent(in) :: root
      real(accs_real), optional, intent(inout) :: p_ref
      real(accs_real), optional, intent(inout) :: p_total
      real(accs_real), optional, intent(inout) :: temp_ref
      real(accs_real), optional, intent(inout) :: dens_ref
      real(accs_real), optional, intent(inout) :: visc_ref
      real(accs_real), optional, intent(inout) :: velo_ref
      real(accs_real), optional, intent(inout) :: len_ref      
      integer, optional, intent(inout) :: pref_at_cell
    end subroutine

    !> @brief Get variables to be solved
    !
    !> @details By default, all variables will be solved. Using this 
    !! "solve" keyword, the user can specifically request that 
    !! certain variables will not be solved by setting in to "off"
    !
    !> @param[in] root - the entry point to the config file    
    !> @param[in,out] u_sol - solve u on/off
    !> @param[in,out] v_sol - solve v on/off
    !> @param[in,out] w_sol - solve w on/off
    !> @param[in,out] p_sol - solve p on/off
    !
    !> @todo extend list of variables 
    module subroutine get_solve(root, u_sol, v_sol, w_sol, p_sol)
      class(*), pointer, intent(in) :: root
      character(len=:), allocatable, optional, intent(inout) :: u_sol
      character(len=:), allocatable, optional, intent(inout) :: v_sol
      character(len=:), allocatable, optional, intent(inout) :: w_sol
      character(len=:), allocatable, optional, intent(inout) :: p_sol
    end subroutine

    !> @brief Get solvers to be used
    !
    !> @details Get the solvers that are to be used for each of
    !! the variables. Solver types are defined by integer values
    !
    !> @param[in] root - the entry point to the config file    
    !> @param[in,out] u_solver - solver to be used for u
    !> @param[in,out] v_solver - solver to be used for v
    !> @param[in,out] w_solver - solver to be used for w
    !> @param[in,out] p_solver - solver to be used for p
    !
    !> @todo extend list of variables   
    module subroutine get_solver(root, u_solver, v_solver, w_solver, p_solver, te_solver, ed_solver)
      class(*), pointer, intent(in) :: root
      integer, optional, intent(inout) :: u_solver
      integer, optional, intent(inout) :: v_solver
      integer, optional, intent(inout) :: w_solver
      integer, optional, intent(inout) :: p_solver
      integer, optional, intent(inout) :: te_solver
      integer, optional, intent(inout) :: ed_solver
    end subroutine

    !> @brief Get transient status
    !
    !> @details Enables/disables unsteady solution algorithm
    !
    !> @param[in] root - the entry point to the config file   
    !> @param[in,out] transient_type - "euler" (first order) or "quad" (second order)
    !> @param[in,out] dt - time interval (seconds) between two consecutive time steps
    !> @param[in,out] gamma - euler blending factor which blends quad
    !> @param[in,out] max_sub_step - maximum number of sub-iterations at each time step
    module subroutine get_transient(root, transient_type, dt, gamma, max_sub_steps)
      class(*), pointer, intent(in) :: root
      character(len=:), allocatable, intent(inout) :: transient_type
      real(accs_real), intent(inout) :: dt
      real(accs_real), intent(inout) :: gamma
      integer, intent(inout) :: max_sub_steps
    end subroutine

    !> @brief Get target residual
    !
    !> @details Get the convergence criterion. 
    !! The calculation will stop when the residuals (L2-norm) of the 
    !! governing equations are less than the target residual.
    !
    !> @param[in] root - the entry point to the config file   
    !> @param[in,out] residual - convergence criterion
    module  subroutine get_target_residual(root, residual)
      class(*), pointer, intent(in) :: root
      real(accs_real), intent(inout) :: residual
    end subroutine

    !> @brief Get grid cell to monitor
    !
    !> @details Get the grid cell at which to monitor the values
    !! of the flow variables (U,V,W,P,TE,ED and T)
    !
    !> @param[in] root - the entry point to the config file   
    !> @param[in,out] monitor_cell - grid cell ID
    module  subroutine get_monitor_cell(root, monitor_cell)
      class(*), pointer, intent(in) :: root
      integer, intent(inout) :: monitor_cell
    end subroutine

    !> @brief Get convection schemes 
    !
    !> @details Get convection schemes to be used for the 
    !! different variables. The convection schemes are defined
    !! by integer values.
    !
    !> @param[in] root - the entry point to the config file   
    !> @param[in,out] u_conv - convection scheme for u
    !> @param[in,out] v_conv - convection scheme for v
    !> @param[in,out] w_conv - convection scheme for w
    !> @param[in,out] te_conv - convection scheme for te
    !> @param[in,out] ed_conv - convection scheme for ed
    module  subroutine get_convection_scheme(root, u_conv, v_conv, w_conv, te_conv, ed_conv)
      class(*), pointer, intent(in) :: root
      integer, optional, intent(inout) :: u_conv
      integer, optional, intent(inout) :: v_conv
      integer, optional, intent(inout) :: w_conv
      integer, optional, intent(inout) :: te_conv
      integer, optional, intent(inout) :: ed_conv
    end subroutine

    !> @brief Get blending factor values 
    !
    !> @details Get blending factors
    !
    !> @param[in] root - the entry point to the config file
    !> @param[in,out] u_blend - blending factor for u
    !> @param[in,out] v_blend - blending factor for v
    !> @param[in,out] w_blend - blending factor for w
    !> @param[in,out] te_blend - blending factor for te
    !> @param[in,out] ed_blend - blending factor for ed
    module subroutine get_blending_factor(root, u_blend, v_blend, w_blend, te_blend, ed_blend)
      class(*), pointer, intent(in) :: root
      real(accs_real), optional, intent(inout) :: u_blend
      real(accs_real), optional, intent(inout) :: v_blend
      real(accs_real), optional, intent(inout) :: w_blend
      real(accs_real), optional, intent(inout) :: te_blend
      real(accs_real), optional, intent(inout) :: ed_blend
    end subroutine

    !> @brief Get relaxation factor values 
    !
    !> @details Get relaxation factors
    !
    !> @param[in] root - the entry point to the config file
    !> @param[in,out] u_relax - relaxation factor for u
    !> @param[in,out] v_relax - relaxation factor for v
    !> @param[in,out] p_relax - relaxation factor for p
    !> @param[in,out] te_relax - relaxation factor for te
    !> @param[in,out] ed_relax - relaxation factor for ed
    module subroutine get_relaxation_factor(root, u_relax, v_relax, p_relax, te_relax, ed_relax)
      class(*), pointer, intent(in) :: root
      real(accs_real), optional, intent(inout) :: u_relax
      real(accs_real), optional, intent(inout) :: v_relax
      real(accs_real), optional, intent(inout) :: p_relax
      real(accs_real), optional, intent(inout) :: te_relax
      real(accs_real), optional, intent(inout) :: ed_relax
    end subroutine

    !> @brief Get output frequency 
    !
    !> @details Get output frequency, set with keywords "every"
    !! "iter" or both.
    !
    !> @param[in] root - the entry point to the config file
    !> @param[inout] output_freq - the frequency of writing output files
    !> @param[inout] output iter - output files are written every output_iter/steps
    module subroutine get_output_frequency(root, output_freq, output_iter)
      class(*), pointer, intent(in) :: root
      integer, intent(inout) :: output_freq
      integer, intent(inout) :: output_iter
    end subroutine

    !> @brief Get output file format 
    !
    !> @param[in] root - the entry point to the config file
    !> @param[inout] plot_format - output format (e.g. vtk)
    module subroutine get_plot_format(root, plot_format)
      class(*), pointer, intent(in) :: root
      character(len=:), allocatable, intent(inout) :: plot_format
    end subroutine

    !> @brief Get output type variables 
    !
    !> @param[in] root - the entry point to the config file
    !> @param[inout] post_type - values at cell centres or cell vertices?
    !> @param[inout] post_vars - variables to be written out
    module subroutine get_output_type(root, post_type, post_vars)
      class(*), pointer, intent(in) :: root
      character(len=:), allocatable, intent(inout) :: post_type
      character(len=2), dimension(10), intent(inout) :: post_vars    
    end subroutine


    !> @brief Get boundary condntions 
    !
    !> @param[in] root - the entry point to the config file
    !> @param[inout] bnd_region - array of boundary region names
    !> @param[inout] bnd_type - array of boundary types (e.g. periodic, symmetric, ...)
    !> @param[inout] bnd_vector - array of boundary vectors

    module subroutine get_boundaries(root, bnd_region, bnd_type, bnd_vector)
      class(*), pointer, intent(in) :: root
      character(len=16), dimension(:), allocatable, intent(inout) :: bnd_region
      character(len=16), dimension(:), allocatable, intent(inout) :: bnd_type
      real(accs_real), dimension(:,:), allocatable, intent(inout) :: bnd_vector
    end subroutine

  end interface
end module read_config
  