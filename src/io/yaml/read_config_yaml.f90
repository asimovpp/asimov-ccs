!v Submodule implementation file read_config_utils.mod
!
!  Module implementing the interface to read YAML config file
!
!  @build yaml
submodule(read_config) read_config_utils
#include "ccs_macros.inc"

  use utils, only: exit_print, debug_print, str
  use yaml_types, only: type_dictionary, &
                        type_error, &
                        type_list, &
                        type_list_item, &
                        type_scalar
  use yaml, only: parse, error_length
  use boundary_conditions, only: set_bc_attribute, allocate_bc_arrays

  implicit none

  interface get_value
    module procedure get_integer_value
    module procedure get_real_value
    module procedure get_string_value
  end interface

contains

  !> Gets the integer value associated with the keyword from dict
  subroutine get_integer_value(dict, keyword, int_val)
    class(*), pointer, intent(in) :: dict     !< The dictionary
    character(len=*), intent(in) :: keyword   !< The key 
    integer, intent(out) :: int_val           !< The corresponding value

    type(type_error), pointer :: io_err

    select type (dict)
    type is (type_dictionary)

      int_val = dict%get_integer(keyword, error=io_err)
      call error_handler(io_err)

    class default
      call error_abort("Unknown type")
    end select

    if (associated(io_err) .eqv. .true.) then
      call error_abort("Error reading " // keyword)
    end if

  end subroutine
    
  !v Gets the real value specified by the keyword from the dictionary. Returns a flag indicating 
  !  whether the key-value pair is present in the dictionary. Takes a flag indicating whether the 
  !  value is required.
  subroutine get_real_value(dict, keyword, real_val, value_present, required)
    class (*), pointer, intent(in) :: dict            !< The dictionary to read from
    character (len=*), intent(in) :: keyword          !< The key to read
    real(ccs_real), intent(out)  :: real_val          !< The value read from the dictionary
    logical, intent(inout), optional :: value_present !< Indicates whether the key-value pair is present in the dictionary
    logical, intent(in), optional :: required         !< Flag indicating whether the value is required

    type(type_error), pointer :: io_err

    select type (dict)
    type is (type_dictionary)

      real_val = dict%get_real(keyword,error=io_err)
      if (present(value_present)) then
        if (associated(io_err)) then 
          value_present = .false.
        else
          value_present = .true.
        end if
      end if
      if (present(required)) then
        if (required .eqv. .true.) then
          call error_handler(io_err)  
        end if
      end if
      
    class default
      call error_abort("Unknown type")
    end select

    if((associated(io_err) .eqv. .true.) .and. present(required)) then 
      if (required .eqv. .true.) then
        call error_abort("Error reading " // keyword)
      end if
    end if

  end subroutine

  !> Gets the string associated with the keyword from dict
  subroutine get_string_value(dict, keyword, string_val)
    class(*), pointer, intent(in) :: dict                       !< The dictionary
    character(len=*), intent(in) :: keyword                     !< The key
    character(len=:), allocatable, intent(inout) :: string_val  !< The corresponding value

    type(type_error), pointer :: io_err

    select type (dict)
    type is (type_dictionary)

      string_val = trim(dict%get_string(keyword, error=io_err))
      call error_handler(io_err)

    class default
      call error_abort("Unknown type")
    end select

    if (associated(io_err) .eqv. .true.) then
      call error_abort("Error reading " // keyword)
    end if

  end subroutine

  subroutine error_handler(io_err)
    type(type_error), pointer, intent(inout) :: io_err

    if (associated(io_err)) then
      print *, trim(io_err%message)
    end if

  end subroutine

  !v Get the name of the test case
  !
  !  Get the case name for the configuration file and store it in a string.
  module subroutine get_case_name(config_file, title)
    class(*), pointer, intent(in) :: config_file          !< the entry point to the config file
    character(len=:), allocatable, intent(inout) :: title !< the case name string

    call get_value(config_file, "title", title)

  end subroutine

  !v Get the number of steps
  !
  !  Get the maximum number of iterations to be preformed in the current run
  module subroutine get_steps(config_file, steps)
    class(*), pointer, intent(in) :: config_file !< the entry point to the config file
    integer, intent(inout) :: steps              !< the maximum number of iterations

    call get_value(config_file, 'steps', steps)

  end subroutine

  !v Get source of initial values
  !
  !  Get the source of the initial values - accepted values are "user", "field" or "step"
  module subroutine get_init(config_file, init, u_init, v_init, w_init, te_init, ed_init)
    class(*), pointer, intent(in) :: config_file         !< the entry point to the config file
    character(len=:), allocatable, intent(inout) :: init !< the source of the initial values (user or field)
    integer, optional, intent(inout) :: u_init           !< initial value for u
    integer, optional, intent(inout) :: v_init           !< initial value for v
    integer, optional, intent(inout) :: w_init           !< initial value for w
    integer, optional, intent(inout) :: te_init          !< initial value for te
    integer, optional, intent(inout) :: ed_init          !< initial value for ed

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type (config_file)
    type is (type_dictionary)

      dict => config_file%get_dictionary('init', required=.true., error=io_err)

      call get_value(dict, "type", init)

      if (present(u_init)) then
        call get_value(dict, "u", u_init)
      end if

      if (present(v_init)) then
        call get_value(dict, "v", v_init)
      end if

      if (present(w_init)) then
        call get_value(dict, "w", w_init)
      end if

      if (present(te_init)) then
        call get_value(dict, "te", te_init)
      end if

      if (present(ed_init)) then
        call get_value(dict, "ed", ed_init)
      end if

    class default
      call error_abort("Unknown type")
    end select

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

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type (config_file)
    type is (type_dictionary)

      dict => config_file%get_dictionary('reference_numbers', required=.true., error=io_err)

      ! Pressure
      if (present(p_ref)) then
        call get_value(dict, "pressure", p_ref)
      end if

      ! Pressure_total
      if (present(p_total)) then
        call get_value(dict, "pressure_total", p_total)
      end if

      ! Temperature
      if (present(temp_ref)) then
        call get_value(dict, "temperature", temp_ref)
      end if

      ! Density
      if (present(dens_ref)) then
        call get_value(dict, "density", dens_ref)
      end if

      ! Viscosity
      if (present(visc_ref)) then
        call get_value(dict, "viscosity", visc_ref)
      end if

      ! Velocity
      if (present(velo_ref)) then
        call get_value(dict, "velocity", velo_ref)
      end if

      ! Length
      if (present(len_ref)) then
        call get_value(dict, "length", len_ref)
      end if

      ! Pref_at_cell
      if (present(pref_at_cell)) then
        call get_value(dict, "pref_at_cell", pref_at_cell)
      end if

    class default
      call error_abort("Unknown type")
    end select

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

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type (config_file)
    type is (type_dictionary)

      dict => config_file%get_dictionary('solve', required=.true., error=io_err)

      ! Solve u?
      if (present(u_sol)) then
        call get_value(dict, "u", u_sol)
      end if

      ! Solve v?
      if (present(v_sol)) then
        call get_value(dict, "v", v_sol)
      end if

      ! Solve w?
      if (present(w_sol)) then
        call get_value(dict, "w", w_sol)
      end if

      ! Solve p?
      if (present(p_sol)) then
        call get_value(dict, "p", p_sol)
      end if

    class default
      call error_abort("Unknown type")
    end select

  end subroutine

  !v Get solvers to be used
  !
  !  Get the solvers that are to be used for each of the variables. Solver types are defined by integer values
  !
  !  @todo extend list of variables
  module subroutine get_solver(config_file, u_solver, v_solver, w_solver, p_solver, te_solver, ed_solver)

    class(*), pointer, intent(in) :: config_file  !< the entry point to the config file
    integer, optional, intent(inout) :: u_solver  !< solver to be used for u
    integer, optional, intent(inout) :: v_solver  !< solver to be used for v
    integer, optional, intent(inout) :: w_solver  !< solver to be used for w
    integer, optional, intent(inout) :: p_solver  !< solver to be used for p
    integer, optional, intent(inout) :: te_solver !< solver to be used for te
    integer, optional, intent(inout) :: ed_solver !< solver to be used for ed

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type (config_file)
    type is (type_dictionary)

      dict => config_file%get_dictionary('solver', required=.true., error=io_err)

      ! Get u_solver
      if (present(u_solver)) then
        call get_value(dict, "u", u_solver)
      end if

      ! Get v_solver
      if (present(v_solver)) then
        call get_value(dict, "v", v_solver)
      end if

      ! Get w_solver
      if (present(w_solver)) then
        call get_value(dict, "w", w_solver)
      end if

      ! Get p_solver
      if (present(p_solver)) then
        call get_value(dict, "p", p_solver)
      end if

      ! Get te_solver
      if (present(te_solver)) then
        call get_value(dict, "te", te_solver)
      end if

      ! Get ed_solver
      if (present(ed_solver)) then
        call get_value(dict, "ed", ed_solver)
      end if

    class default
      call error_abort("Unknown type")
    end select

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

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type (config_file)
    type is (type_dictionary)

      dict => config_file%get_dictionary('transient', required=.false., error=io_err)

      ! Transient type (euler/quad)
      call get_value(dict, "type", transient_type)

      ! Dt
      call get_value(dict, "dt", dt)

      ! Gamma
      call get_value(dict, "gamma", euler_blend)

      ! Maximum number of sub steps
      call get_value(dict, "max_sub_steps", max_sub_steps)

    class default
      call error_abort("Unknown type")
    end select

  end subroutine

  !v Get target residual
  !
  !  Get the convergence criterion.
  !  The calculation will stop when the residuals (L2-norm) of the
  !  governing equations are less than the target residual.
  module subroutine get_target_residual(config_file, residual)
    class(*), pointer, intent(in) :: config_file  !< the entry point to the config file
    real(ccs_real), intent(inout) :: residual     !< convergence criterion

    call get_value(config_file, "target_residual", residual)

  end subroutine

  !v Get grid cell to monitor
  !
  !  Get the grid cell at which to monitor the values of the flow variables (U,V,W,P,TE,ED and T)
  module subroutine get_monitor_cell(config_file, monitor_cell)
    class(*), pointer, intent(in) :: config_file !< the entry point to the config file
    integer, intent(inout) :: monitor_cell       !< grid cell ID

    call get_value(config_file, "monitor_cell", monitor_cell)

  end subroutine

  !v Get convection schemes
  !
  !  Get convection schemes to be used for the different variables. The convection schemes are defined
  !  by integer values.
  module subroutine get_convection_scheme(config_file, u_conv, v_conv, w_conv, te_conv, ed_conv)
    class(*), pointer, intent(in) :: config_file !< the entry point to the config file
    integer, optional, intent(inout) :: u_conv   !< convection scheme for u
    integer, optional, intent(inout) :: v_conv   !< convection scheme for v
    integer, optional, intent(inout) :: w_conv   !< convection scheme for w
    integer, optional, intent(inout) :: te_conv  !< convection scheme for te
    integer, optional, intent(inout) :: ed_conv  !< convection scheme for ed

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type (config_file)
    type is (type_dictionary)

      dict => config_file%get_dictionary('convection_scheme', required=.false., error=io_err)

      if (present(u_conv)) then
        call get_value(dict, "u", u_conv)
      end if

      if (present(v_conv)) then
        call get_value(dict, "v", v_conv)
      end if

      if (present(w_conv)) then
        call get_value(dict, "w", w_conv)
      end if

      if (present(te_conv)) then
        call get_value(dict, "te", te_conv)
      end if

      if (present(ed_conv)) then
        call get_value(dict, "ed", ed_conv)
      end if

    class default
      call error_abort("Unknown type")
    end select

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

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type (config_file)
    type is (type_dictionary)

      dict => config_file%get_dictionary('blending_factor', required=.false., error=io_err)

      if (present(u_blend)) then
        call get_value(dict, "u", u_blend)
      end if

      if (present(v_blend)) then
        call get_value(dict, "v", v_blend)
      end if

      if (present(w_blend)) then
        call get_value(dict, "w", w_blend)
      end if

      if (present(te_blend)) then
        call get_value(dict, "te", te_blend)
      end if

      if (present(ed_blend)) then
        call get_value(dict, "ed", ed_blend)
      end if

    class default
      call error_abort("Unknown type")
    end select

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

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type (config_file)
    type is (type_dictionary)

      dict => config_file%get_dictionary('relaxation_factor', required=.false., error=io_err)

      if (present(u_relax)) then
        call get_value(dict, "u", u_relax)
      end if

      if (present(v_relax)) then
        call get_value(dict, "v", v_relax)
      end if

      if (present(p_relax)) then
        call get_value(dict, "p", p_relax)
      end if

      if (present(te_relax)) then
        call get_value(dict, "te", te_relax)
      end if

      if (present(ed_relax)) then
        call get_value(dict, "ed", ed_relax)
      end if

    class default
      call error_abort("Unknown type")
    end select

  end subroutine

  !v Get output frequency
  !
  !  Get output frequency, set with keywords "every", "iter", or both.
  module subroutine get_output_frequency(config_file, output_freq, output_iter)
    class(*), pointer, intent(in) :: config_file !< the entry point to the config file
    integer, intent(inout) :: output_freq        !< the frequency of writing output files
    integer, intent(inout) :: output_iter        !< output files are written every output_iter/steps

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type (config_file)
    type is (type_dictionary)

      dict => config_file%get_dictionary('output', required=.false., error=io_err)

      call get_value(dict, "every", output_freq)
      call get_value(dict, "iter", output_iter)

    class default
      call error_abort("Unknown type")
    end select

  end subroutine

  !> Get output file format
  module subroutine get_plot_format(config_file, plot_format)
    class(*), pointer, intent(in) :: config_file                !< the entry point to the config file
    character(len=:), allocatable, intent(inout) :: plot_format !< output format (e.g. vtk)

    call get_value(config_file, "plot_format", plot_format)

  end subroutine

  !> Get output type and variables
  module subroutine get_output_type(config_file, post_type, post_vars)
    class(*), pointer, intent(in) :: config_file                !< the entry point to the config file
    character(len=:), allocatable, intent(inout) :: post_type   !< values at cell centres or cell vertices?
    character(len=2), dimension(10), intent(inout) :: post_vars !< variables to be written out

    class(*), pointer :: dict
    class(type_list), pointer :: list
    class(type_list_item), pointer :: item
    type(type_error), pointer :: io_err
    integer :: idx

    select type (config_file)
    type is (type_dictionary)

      dict => config_file%get_dictionary('post', required=.false., error=io_err)

      call get_value(dict, "type", post_type)

      select type (dict)
      type is (type_dictionary)

        list => dict%get_list('variables', required=.false., error=io_err)
        call error_handler(io_err)

        item => list%first
        idx = 1
        do while (associated(item))
          select type (element => item%node)
          class is (type_scalar)
            post_vars(idx) = trim(element%string)
            print *, post_vars(idx)
            item => item%next
            idx = idx + 1
          end select
        end do

      class default
        call error_abort("Unknown type")
      end select

    class default
      call error_abort("Unknown type")
    end select

  end subroutine

  module subroutine get_n_boundaries(filename, n_boundaries)
    character(len=*), intent(in) :: filename      !< name of the config file
    integer(ccs_int), intent(out) :: n_boundaries !< number of boundaries

    class(*), pointer :: config_file
    class(*), pointer :: dict
    character(len=error_length) :: error
    type(type_error), pointer :: io_err

    config_file => parse(filename, error=error)
    if (error/='') then
      call error_abort(trim(error))
    endif

    select type (config_file)
    type is (type_dictionary)
      dict => config_file%get_dictionary("boundaries", required=.true., error=io_err)
      call error_handler(io_err)

      call get_value(dict, "n_boundaries", n_boundaries)
    class default
       call error_abort("type unhandled")
    end select
  end subroutine get_n_boundaries

  module subroutine get_bc_variables(filename, variables)
    character(len=*), intent(in) :: filename                              !< name of the config file
    character(len=6), dimension(:), allocatable, intent(out) :: variables !< string array indicating variables used in BCs
    
    class(*), pointer :: config_file
    class(*), pointer :: dict
    class(*), pointer :: dict_var
    type(type_error), pointer :: io_err
    integer(ccs_int) :: i
    integer(ccs_int) :: n_var
    character(len=25) :: key
    character(len=:), allocatable :: variable
    character(len=error_length) :: error

    config_file => parse(filename, error=error)
    if (error/='') then
      call error_abort(trim(error))
    endif

    select type (config_file)
    type is (type_dictionary)
      dict => config_file%get_dictionary("variables", required=.true., error=io_err)
      call error_handler(io_err)
      
      call get_value(dict, "n_variables", n_var)
      allocate(variables(n_var))

      do i = 1, n_var
        write(key, '(A, I0)') "variable_", i 
        select type (dict)
        type is (type_dictionary)
           dict_var => dict%get_dictionary(key, required=.true., error=io_err)
           call error_handler(io_err)
           call get_value(dict_var, "name", variable)
           write(variables(i), '(A)') trim(variable)
        class default
          call error_abort("type unhandled")
        end select
      end do
    class default
      call error_abort("type unhandled")
    end select
  end subroutine get_bc_variables

  !> Gets the specified field value from the config file and writes to given bcs struct
  module subroutine get_bc_field(config_file, bc_field, phi) 
    class(*), pointer, intent(in) :: config_file  !< pointer to configuration file
    character(len=*), intent(in) :: bc_field      !< string indicating which field to read from BCs
    class(field), intent(inout) :: phi            !< field structure
    
    ! local variables
    class(*), pointer :: dict
    class(*), pointer :: dict2
    integer(ccs_int) :: i
    integer(ccs_int) :: n_boundaries
    type(type_error), pointer :: io_err
    character(len=:), allocatable :: bc_field_string
    character(len=25) :: boundary_index
    integer(ccs_int) :: bc_id

    class(*), pointer :: variable_dict
    character(len=25) :: variable
    character(len=:), allocatable :: bc_type
    real(ccs_real) :: bc_value
    logical :: field_exists
    
    select type (config_file)
    type is (type_dictionary)
      dict => config_file%get_dictionary("boundaries", required=.true., error=io_err)
      call error_handler(io_err)

      i = 1
      n_boundaries = size(phi%bcs%ids)
      do while (i <= n_boundaries)
        write(boundary_index, '(A, I0)') "boundary_", i
        select type (dict)
        type is (type_dictionary)
          dict2 => dict%get_dictionary(boundary_index, required=.true., error=io_err)
          call error_handler(io_err)

          select case (bc_field)
          case ("name")
            call dprint("starting reading name")
            call get_value(dict2, bc_field, bc_field_string)
            call set_bc_attribute(i, bc_field, bc_field_string, phi%bcs)
            call dprint("finished reading name")
          case ("id")
            call dprint("starting reading id")
            call get_value(dict2, bc_field, bc_id)
            call set_bc_attribute(i, bc_id, phi%bcs)
            call dprint("finished reading id")
          case default
            call dprint("starting reading " // bc_field)
            !call get_bc_field_data(dict2, i, bc_field, phi)
            select type (dict2)
            type is (type_dictionary)
              write (variable, '(A, A)') "variable_", bc_field
              variable_dict => dict2%get_dictionary(variable, required=.false., error=io_err)
              !call error_handler(io_err)
              call dprint("read dictionary " // bc_field)

              call get_value(variable_dict, "type", bc_type)
              call set_bc_attribute(i, "type", bc_type, phi%bcs)
              call dprint("read type " // bc_type)

              call get_value(variable_dict, "value", bc_value, field_exists)
              call dprint("read value " // bc_field // " " // str(bc_value))
              if (field_exists) then
                call set_bc_attribute(i, bc_value, phi%bcs)
              end if
              call dprint("finished reading " // bc_field)
            end select
          end select
        end select
        i = i+1
      end do
    class default
      call error_abort("type unhandled")
    end select
  end subroutine get_bc_field

  subroutine get_bc_field_data(dict, boundary_index, bc_field, phi)
    class(*), pointer, intent(in) :: dict
    integer(ccs_int), intent(in) :: boundary_index
    character(len=*), intent(in) :: bc_field       !< string indicating which field to read from BCs
    class(field), intent(inout) :: phi

    class(*), pointer :: variable_dict
    character(len=25) :: variable
    !character(len=:), allocatable :: bc_type
    real(ccs_real) :: bc_value
    logical :: field_exists
    type(type_error), pointer :: io_err

    select type (dict)
    type is (type_dictionary)
      write (variable, '(A, A)') "variable_", bc_field
      variable_dict => dict%get_dictionary(variable, required=.false., error=io_err)
      !call error_handler(io_err)

      !call get_value(variable_dict, "type", bc_type)
      !call set_bc_attribute(boundary_index, bc_field, bc_type, phi%bcs)

      call get_value(variable_dict, "value", bc_value, field_exists)
      if (field_exists) then
        call set_bc_attribute(boundary_index, bc_value, phi%bcs)
      end if
    class default
      call error_abort("type unhandled")
    end select
  end subroutine get_bc_field_data

end submodule read_config_utils
