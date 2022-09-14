!v Module file read_config.mod
!
!  Module defining interface to read YAML config file

module read_config

  use kinds, only: ccs_real, ccs_int
  use types, only: bc_config, field
  use constants, only: ccs_string_len

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
  public :: get_bc_variables
  public :: get_bc_field
  public :: get_boundary_count

  interface

    !v Get the name of the test case
    !
    !  Get the case name for the configuration file and
    !  store it in a string.
    module subroutine get_case_name(config_file, title)
      class(*), pointer, intent(in) :: config_file            !< the entry point to the config file
      character(len=:), allocatable, intent(inout) :: title   !< the case name string
    end subroutine

    !v Get the number of steps
    !
    !  Get the maximum number of iterations
    !  to be preformed in the current run
    module subroutine get_steps(config_file, steps)
      class(*), pointer, intent(in) :: config_file  !< the entry point to the config file
      integer, intent(inout) :: steps               !< the maximum number of iterations
    end subroutine

    !v Get source of initial values
    !
    !  Get the source of the initial values - accepted
    !  values are "user", "field" or "step"
    module subroutine get_init(config_file, init, u_init, v_init, w_init, te_init, ed_init)
      class(*), pointer, intent(in) :: config_file
      character(len=:), allocatable, intent(inout) :: init
      integer, optional, intent(inout) :: u_init
      integer, optional, intent(inout) :: v_init
      integer, optional, intent(inout) :: w_init
      integer, optional, intent(inout) :: te_init
      integer, optional, intent(inout) :: ed_init

    end subroutine

    !v Get reference numbers
    !
    !  Get the reference numbers, the fluid properties and the operating condition
    module subroutine get_reference_number(config_file, p_ref, p_total, temp_ref, &
                                           dens_ref, visc_ref, velo_ref, len_ref, pref_at_cell)
      class(*), pointer, intent(in) :: config_file        !< the entry point to the config file
      real(ccs_real), optional, intent(inout) :: p_ref    !< reference pressure
      real(ccs_real), optional, intent(inout) :: p_total  !< total pressure
      real(ccs_real), optional, intent(inout) :: temp_ref !< reference temperature
      real(ccs_real), optional, intent(inout) :: dens_ref !< reference density
      real(ccs_real), optional, intent(inout) :: visc_ref !< laminar viscosity
      real(ccs_real), optional, intent(inout) :: velo_ref !< reference velocity
      real(ccs_real), optional, intent(inout) :: len_ref  !< reference length, used to define the Reynolds number of the flow
      integer, optional, intent(inout) :: pref_at_cell    !< cell at which the reference pressure is set
    end subroutine

    !v Get variables to be solved
    !
    !  By default, all variables will be solved. Using this
    !  "solve" keyword, the user can specifically request that
    !  certain variables will not be solved by setting in to "off"
    !
    !  @todo extend list of variables
    module subroutine get_solve(config_file, u_sol, v_sol, w_sol, p_sol)
      class(*), pointer, intent(in) :: config_file                      !< the entry point to the config file
      character(len=:), allocatable, optional, intent(inout) :: u_sol   !< solve u on/off
      character(len=:), allocatable, optional, intent(inout) :: v_sol   !< solve v on/off
      character(len=:), allocatable, optional, intent(inout) :: w_sol   !< solve w on/off
      character(len=:), allocatable, optional, intent(inout) :: p_sol   !< solve p on/off
    end subroutine

    !v Get solvers to be used
    !
    !  Get the solvers that are to be used for each of
    !  the variables. Solver types are defined by integer values
    !  @todo extend list of variables
    module subroutine get_solver(config_file, u_solver, v_solver, w_solver, p_solver, te_solver, ed_solver)
      class(*), pointer, intent(in) :: config_file    !< the entry point to the config file
      integer, optional, intent(inout) :: u_solver    !< solver to be used for u
      integer, optional, intent(inout) :: v_solver    !< solver to be used for v
      integer, optional, intent(inout) :: w_solver    !< solver to be used for w
      integer, optional, intent(inout) :: p_solver    !< solver to be used for p
      integer, optional, intent(inout) :: te_solver   !< solver to be used for te
      integer, optional, intent(inout) :: ed_solver   !< solver to be used for ed
    end subroutine

    !v Get transient status
    !
    !  Enables/disables unsteady solution algorithm
    module subroutine get_transient(config_file, transient_type, dt, euler_blend, max_sub_steps)
      class(*), pointer, intent(in) :: config_file                   !< the entry point to the config file
      character(len=:), allocatable, intent(inout) :: transient_type !< "euler" (first order) or "quad" (second order)
      real(ccs_real), intent(inout) :: dt                            !< time interval (seconds) between two consecutive time steps
      real(ccs_real), intent(inout) :: euler_blend                   !< gamma, euler blending factor which blends quad
      integer, intent(inout) :: max_sub_steps                        !< maximum number of sub-iterations at each time step
    end subroutine

    !v Get target residual
    !
    !  Get the convergence criterion.
    !  The calculation will stop when the residuals (L2-norm) of the
    !  governing equations are less than the target residual.
    module subroutine get_target_residual(config_file, residual)
      class(*), pointer, intent(in) :: config_file  !< the entry point to the config file
      real(ccs_real), intent(inout) :: residual     !< convergence criterion
    end subroutine

    !v Get grid cell to monitor
    !
    !  Get the grid cell at which to monitor the values
    !  of the flow variables (U,V,W,P,TE,ED and T)
    module subroutine get_monitor_cell(config_file, monitor_cell)
      class(*), pointer, intent(in) :: config_file  !< the entry point to the config file
      integer, intent(inout) :: monitor_cell        !< grid cell ID
    end subroutine

    !v Get convection schemes
    !
    !  Get convection schemes to be used for the
    !  different variables. The convection schemes are defined
    !  by integer values.
    module subroutine get_convection_scheme(config_file, u_conv, v_conv, w_conv, te_conv, ed_conv)
      class(*), pointer, intent(in) :: config_file  !< the entry point to the config file
      integer, optional, intent(inout) :: u_conv    !< convection scheme for u
      integer, optional, intent(inout) :: v_conv    !< convection scheme for v
      integer, optional, intent(inout) :: w_conv    !< convection scheme for w
      integer, optional, intent(inout) :: te_conv   !< convection scheme for te
      integer, optional, intent(inout) :: ed_conv   !< convection scheme for ed
    end subroutine

    !v Get blending factor values
    !
    !  Get blending factors
    module subroutine get_blending_factor(config_file, u_blend, v_blend, w_blend, te_blend, ed_blend)
      class(*), pointer, intent(in) :: config_file        !< the entry point to the config file
      real(ccs_real), optional, intent(inout) :: u_blend  !< blending factor for u
      real(ccs_real), optional, intent(inout) :: v_blend  !< blending factor for v
      real(ccs_real), optional, intent(inout) :: w_blend  !< blending factor for w
      real(ccs_real), optional, intent(inout) :: te_blend !< blending factor for te
      real(ccs_real), optional, intent(inout) :: ed_blend !< blending factor for ed
    end subroutine

    !v Get relaxation factor values
    !
    !  Get relaxation factors
    module subroutine get_relaxation_factor(config_file, u_relax, v_relax, p_relax, te_relax, ed_relax)
      class(*), pointer, intent(in) :: config_file        !< the entry point to the config file
      real(ccs_real), optional, intent(inout) :: u_relax  !< relaxation factor for u
      real(ccs_real), optional, intent(inout) :: v_relax  !< relaxation factor for v
      real(ccs_real), optional, intent(inout) :: p_relax  !< relaxation factor for p
      real(ccs_real), optional, intent(inout) :: te_relax !< relaxation factor for te
      real(ccs_real), optional, intent(inout) :: ed_relax !< relaxation factor for ed
    end subroutine

    !v Get output frequency
    !
    !  Get output frequency, set with keywords "every", "iter", or both.
    module subroutine get_output_frequency(config_file, output_freq, output_iter)
      class(*), pointer, intent(in) :: config_file !< the entry point to the config file
      integer, intent(inout) :: output_freq        !< the frequency of writing output files
      integer, intent(inout) :: output_iter        !< output files are written every output_iter/steps
    end subroutine

    !> Get output file format
    module subroutine get_plot_format(config_file, plot_format)
      class(*), pointer, intent(in) :: config_file !< the entry point to the config file
      character(len=:), allocatable, intent(inout) :: plot_format !< output format (e.g. vtk)
    end subroutine

    !> Get output type variables
    module subroutine get_output_type(config_file, post_type, post_vars)
      class(*), pointer, intent(in) :: config_file                !< the entry point to the config file
      character(len=:), allocatable, intent(inout) :: post_type   !< values at cell centres or cell vertices?
      character(len=2), dimension(10), intent(inout) :: post_vars !< variables to be written out
    end subroutine

    !> Gets the specified field value from the config file and writes to given bcs struct
    module subroutine get_bc_field(config_file, bc_field, phi, required)
      class(*), pointer, intent(in) :: config_file  !< pointer to configuration file
      character(len=*), intent(in) :: bc_field      !< string indicating which field to read from BCs
      class(field), intent(inout) :: phi            !< field structure
      logical, optional, intent(in) :: required     !< flag indicating whether field is required
    end subroutine

    !> Gets variables that bcs are defined for
    module subroutine get_bc_variables(filename, variables)
      character(len=*), intent(in) :: filename                                            !< name of the config file
      character(len=ccs_string_len), dimension(:), allocatable, intent(out) :: variables  !< string array indicating variables used in BCs
    end subroutine

    !> Gets the number of boundaries
    module subroutine get_boundary_count(filename, n_boundaries)
      character(len=*), intent(in) :: filename      !< name of the config file
      integer(ccs_int), intent(out) :: n_boundaries !< number of boundaries
    end subroutine
  end interface
end module read_config
