!>  Module file read_config.mod
!
!>  Module defining interface to read YAML config file

module read_config

  use kinds, only : ccs_real, ccs_int
  use types, only : bc_config
    
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
  public :: get_bc_field_data   

  interface

    !>  Get the name of the test case
    !
    !v  Get the case name for the configuration file and 
    !   store it in a string.
    module subroutine get_case_name(config_file, title)
      class(*), pointer, intent(in) :: config_file            !< the entry point to the config file    
      character(len=:), allocatable, intent(inout) :: title   !< the case name string    
    end subroutine
      
    !>  Get the number of steps
    !
    !v  Get the maximum number of iterations 
    !   to be preformed in the current run 
    module  subroutine get_steps(config_file, steps)
      class(*), pointer, intent(in) :: config_file  !< the entry point to the config file   
      integer, intent(inout) :: steps               !< the maximum number of iterations    
    end subroutine
      
    !>  Get source of initial values
    !
    !v  Get the source of the initial values - accepted
    !   values are "user", "field" or "step" 
    module  subroutine get_init(config_file, init, u_init, v_init, w_init, te_init, ed_init)
      class(*), pointer, intent(in) :: config_file
      character(len=:), allocatable, intent(inout) :: init
      integer, optional, intent(inout) :: u_init
      integer, optional, intent(inout) :: v_init
      integer, optional, intent(inout) :: w_init
      integer, optional, intent(inout) :: te_init
      integer, optional, intent(inout) :: ed_init

    end subroutine

    !>  Get reference numbers
    !
    !>  Get the reference numbers, the fluid properties and the operating condition 
    !
    !> @param[in] config_file - the entry point to the config file    
    !> @param[in,out] p_ref - reference pressure 
    !> @param[in,out] p_total - total pressure 
    !> @param[in,out] temp_ref - reference temperature      
    !> @param[in,out] dens_ref - reference density      
    !> @param[in,out] visc_ref - laminar viscosity      
    !> @param[in,out] velo_ref - reference velocity      
    !> @param[in,out] len_ref - reference length, used to define the Reynolds number of the flow      
    !> @param[in,out] pref_at_cell - cell at which the reference pressure is set      
    module subroutine get_reference_number(config_file, p_ref, p_total, temp_ref, &
                                            dens_ref, visc_ref, velo_ref, len_ref, pref_at_cell)
      class(*), pointer, intent(in) :: config_file
      real(ccs_real), optional, intent(inout) :: p_ref
      real(ccs_real), optional, intent(inout) :: p_total
      real(ccs_real), optional, intent(inout) :: temp_ref
      real(ccs_real), optional, intent(inout) :: dens_ref
      real(ccs_real), optional, intent(inout) :: visc_ref
      real(ccs_real), optional, intent(inout) :: velo_ref
      real(ccs_real), optional, intent(inout) :: len_ref      
      integer, optional, intent(inout) :: pref_at_cell
    end subroutine

    !>  Get variables to be solved
    !
    !v  By default, all variables will be solved. Using this 
    ! "solve" keyword, the user can specifically request that 
    ! certain variables will not be solved by setting in to "off"
    !
    !> @todo extend list of variables
    module subroutine get_solve(config_file, u_sol, v_sol, w_sol, p_sol)
      class(*), pointer, intent(in) :: config_file                      !< the entry point to the config file    
      character(len=:), allocatable, optional, intent(inout) :: u_sol   !< solve u on/off
      character(len=:), allocatable, optional, intent(inout) :: v_sol   !< solve v on/off
      character(len=:), allocatable, optional, intent(inout) :: w_sol   !< solve w on/off
      character(len=:), allocatable, optional, intent(inout) :: p_sol   !< solve p on/off
    end subroutine

    !>  Get solvers to be used
    !
    !v  Get the solvers that are to be used for each of
    ! the variables. Solver types are defined by integer values
    !
    !> @param[in] config_file  
    !> @param[in,out] u_solver 
    !> @param[in,out] v_solver 
    !> @param[in,out] w_solver 
    !> @param[in,out] p_solver 
    !
    !> @todo extend list of variables   
    module subroutine get_solver(config_file, u_solver, v_solver, w_solver, p_solver, te_solver, ed_solver)
      class(*), pointer, intent(in) :: config_file    !< the entry point to the config file    
      integer, optional, intent(inout) :: u_solver    !< solver to be used for u
      integer, optional, intent(inout) :: v_solver    !< solver to be used for v
      integer, optional, intent(inout) :: w_solver    !< solver to be used for w
      integer, optional, intent(inout) :: p_solver    !< solver to be used for p
      integer, optional, intent(inout) :: te_solver   !< solver to be used for te
      integer, optional, intent(inout) :: ed_solver   !< solver to be used for ed
    end subroutine

    !>  Get transient status
    !
    !>  Enables/disables unsteady solution algorithm
    !
    !> @param[in] config_file - the entry point to the config file   
    !> @param[in,out] transient_type - "euler" (first order) or "quad" (second order)
    !> @param[in,out] dt - time interval (seconds) between two consecutive time steps
    !> @param[in,out] euler_blend - gamma, euler blending factor which blends quad
    !> @param[in,out] max_sub_step - maximum number of sub-iterations at each time step
    module subroutine get_transient(config_file, transient_type, dt, euler_blend, max_sub_steps)
      class(*), pointer, intent(in) :: config_file
      character(len=:), allocatable, intent(inout) :: transient_type
      real(ccs_real), intent(inout) :: dt
      real(ccs_real), intent(inout) :: euler_blend
      integer, intent(inout) :: max_sub_steps
    end subroutine

    !>  Get target residual
    !
    !v Get the convergence criterion. 
    !  The calculation will stop when the residuals (L2-norm) of the 
    !  governing equations are less than the target residual.
    module  subroutine get_target_residual(config_file, residual)
      class(*), pointer, intent(in) :: config_file  !< the entry point to the config file   
      real(ccs_real), intent(inout) :: residual     !< convergence criterion
    end subroutine

    !>  Get grid cell to monitor
    !
    !v Get the grid cell at which to monitor the values
    !  of the flow variables (U,V,W,P,TE,ED and T)
    !
    !> @param[in] config_file      
    !> @param[in,out] monitor_cell 
    module  subroutine get_monitor_cell(config_file, monitor_cell)
      class(*), pointer, intent(in) :: config_file  !< the entry point to the config file   
      integer, intent(inout) :: monitor_cell        !< grid cell ID
    end subroutine

    !>  Get convection schemes 
    !
    !v  Get convection schemes to be used for the 
    !   different variables. The convection schemes are defined
    !   by integer values.
    module  subroutine get_convection_scheme(config_file, u_conv, v_conv, w_conv, te_conv, ed_conv)
      class(*), pointer, intent(in) :: config_file  !< the entry point to the config file   
      integer, optional, intent(inout) :: u_conv    !< convection scheme for u
      integer, optional, intent(inout) :: v_conv    !< convection scheme for v
      integer, optional, intent(inout) :: w_conv    !< convection scheme for w
      integer, optional, intent(inout) :: te_conv   !< convection scheme for te
      integer, optional, intent(inout) :: ed_conv   !< convection scheme for ed
    end subroutine

    !>  Get blending factor values 
    !
    !>  Get blending factors
    !
    !> @param[in] config_file - the entry point to the config file
    !> @param[in,out] u_blend - blending factor for u
    !> @param[in,out] v_blend - blending factor for v
    !> @param[in,out] w_blend - blending factor for w
    !> @param[in,out] te_blend - blending factor for te
    !> @param[in,out] ed_blend - blending factor for ed
    module subroutine get_blending_factor(config_file, u_blend, v_blend, w_blend, te_blend, ed_blend)
      class(*), pointer, intent(in) :: config_file
      real(ccs_real), optional, intent(inout) :: u_blend
      real(ccs_real), optional, intent(inout) :: v_blend
      real(ccs_real), optional, intent(inout) :: w_blend
      real(ccs_real), optional, intent(inout) :: te_blend
      real(ccs_real), optional, intent(inout) :: ed_blend
    end subroutine

    !>  Get relaxation factor values 
    !
    !> Get relaxation factors
    !
    !> @param[in] config_file - the entry point to the config file
    !> @param[in,out] u_relax - relaxation factor for u
    !> @param[in,out] v_relax - relaxation factor for v
    !> @param[in,out] p_relax - relaxation factor for p
    !> @param[in,out] te_relax - relaxation factor for te
    !> @param[in,out] ed_relax - relaxation factor for ed
    module subroutine get_relaxation_factor(config_file, u_relax, v_relax, p_relax, te_relax, ed_relax)
      class(*), pointer, intent(in) :: config_file
      real(ccs_real), optional, intent(inout) :: u_relax
      real(ccs_real), optional, intent(inout) :: v_relax
      real(ccs_real), optional, intent(inout) :: p_relax
      real(ccs_real), optional, intent(inout) :: te_relax
      real(ccs_real), optional, intent(inout) :: ed_relax
    end subroutine

    !>  Get output frequency 
    !
    !>  Get output frequency, set with keywords "every", "iter", or both.
    !
    !> @param[in] config_file - the entry point to the config file
    !> @param[inout] output_freq - the frequency of writing output files
    !> @param[inout] output iter - output files are written every output_iter/steps
    module subroutine get_output_frequency(config_file, output_freq, output_iter)
      class(*), pointer, intent(in) :: config_file
      integer, intent(inout) :: output_freq
      integer, intent(inout) :: output_iter
    end subroutine

    !>  Get output file format 
    !
    !> @param[in] config_file - the entry point to the config file
    !> @param[inout] plot_format - output format (e.g. vtk)
    module subroutine get_plot_format(config_file, plot_format)
      class(*), pointer, intent(in) :: config_file
      character(len=:), allocatable, intent(inout) :: plot_format
    end subroutine

    !>  Get output type variables 
    !
    !> @param[in] config_file - the entry point to the config file
    !> @param[inout] post_type - values at cell centres or cell vertices?
    !> @param[inout] post_vars - variables to be written out
    module subroutine get_output_type(config_file, post_type, post_vars)
      class(*), pointer, intent(in) :: config_file
      character(len=:), allocatable, intent(inout) :: post_type
      character(len=2), dimension(10), intent(inout) :: post_vars    
    end subroutine
    
    module subroutine get_bc_field_data(config_file, bc_field, bcs)
      class(*), pointer, intent(in) :: config_file
      character(len=*), intent(in) :: bc_field      
      type(bc_config), intent(inout) :: bcs           
    end subroutine
  end interface
end module read_config
  
