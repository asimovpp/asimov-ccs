!v Submodule implementation file read_config_utils.mod
!
!  Module implementing the interface to read YAML config file
!
!  @build yaml
submodule(read_config) read_config_utils
#include "ccs_macros.inc"

  use constants, only: cell_centred_central, cell_centred_upwind
  use utils, only: exit_print, debug_print, str, get_scheme_id
  use fortran_yaml_c_interface, only: parse
  use fortran_yaml_c, only: type_dictionary, &
                            type_error, &
                            type_list, &
                            type_list_item, &
                            type_scalar
  use boundary_conditions, only: set_bc_real_value, set_bc_id, set_bc_type, allocate_bc_arrays

  implicit none

contains

  !> Gets the integer value associated with the keyword from dict
  module subroutine get_integer_value(dict, keyword, int_val, value_present, required)
    class(*), pointer, intent(in) :: dict     !< The dictionary
    character(len=*), intent(in) :: keyword   !< The key
    integer, intent(inout) :: int_val           !< The corresponding value
    logical, intent(out), optional :: value_present !< Indicates whether the key-value pair is present in the dictionary
    logical, intent(in), optional :: required         !< Flag indicating whether the value is required. Absence implies not required

    type(type_error), allocatable :: io_err

    select type (dict)
    type is (type_dictionary)

      int_val = dict%get_integer(keyword, error=io_err)

      if (allocated(io_err)) then
        if (present(value_present)) then
          value_present = .false.
        end if
        if(present(required) .and. required .eqv. .true.) then
          call error_abort("Error reading keyword " // keyword // ". Possibly missing keyword in yaml file.")
        end if
      else
        if (present(value_present)) then
          value_present = .true.
        end if
      end if

    class default
      call error_abort("Unknown type")
    end select
   
  end subroutine

  !v Gets the real value specified by the keyword from the dictionary. Returns a flag indicating
  !  whether the key-value pair is present in the dictionary. Takes a flag indicating whether the
  !  value is required.
  module subroutine get_real_value(dict, keyword, real_val, value_present, required)
    class(*), pointer, intent(in) :: dict            !< The dictionary to read from
    character(len=*), intent(in) :: keyword          !< The key to read
    real(ccs_real), intent(inout) :: real_val          !< The value read from the dictionary
    logical, intent(out), optional :: value_present !< Indicates whether the key-value pair is present in the dictionary
    logical, intent(in), optional :: required         !< Flag indicating whether the value is required. Absence implies not required

    type(type_error), allocatable :: io_err

    select type (dict)
    type is (type_dictionary)

      real_val = dict%get_real(keyword, error=io_err)

      if (allocated(io_err)) then
        if (present(value_present)) then
          value_present = .false.
        end if
        if(present(required) .and. required .eqv. .true.) then
          call error_abort("Error reading keyword " // keyword // ". Possibly missing keyword in yaml file.")
        end if
      else
        if (present(value_present)) then
          value_present = .true.
        end if
      end if

    class default
      call error_abort("Unknown type")
    end select

  end subroutine

  !> Gets the string associated with the keyword from dict
  module subroutine get_string_value(dict, keyword, string_val, value_present, required)
    class(*), pointer, intent(in) :: dict                       !< The dictionary
    character(len=*), intent(in) :: keyword                     !< The key
    character(len=:), allocatable, intent(out) :: string_val  !< The corresponding value
    logical, intent(out), optional :: value_present           !< Indicates whether the key-value pair is present in the dictionary
    logical, optional, intent(in) :: required                   !< Flag indicating whether result is required. Absence implies not required.

    type(type_error), allocatable :: io_err

    select type (dict)
    type is (type_dictionary)

      string_val = ""
      string_val = dict%get_string(keyword, error=io_err)

      if (allocated(io_err)) then
        ! print *, "Keyword ", keyword, " not in YAML file"
        if (present(value_present)) then
          value_present = .false.
        end if
        if(present(required) .and. required .eqv. .true.) then
          call error_abort("Error reading keyword " // keyword // ". Possibly missing keyword in yaml file.")
        end if
      else
        string_val = trim(string_val)
        if (present(value_present)) then
          value_present = .true.
        end if
      end if
      
    class default
      call error_abort("Unknown type")
    end select

  end subroutine

  module subroutine get_logical_value(dict, keyword, logical_val, value_present, required)
    class(*), pointer, intent(in) :: dict                       !< The dictionary
    character(len=*), intent(in) :: keyword                     !< The key
    logical, intent(inout) :: logical_val                       !< The corresponding value
    logical, intent(out), optional :: value_present           !< Indicates whether the key-value pair is present in the dictionary
    logical, intent(in), optional :: required                   !< Flag indicating whether result is required. Absence implies not required.

    type(type_error), allocatable :: io_err

    select type (dict)
    type is (type_dictionary)

      logical_val = dict%get_logical(keyword, error=io_err)
      if (allocated(io_err)) then
        ! print *, "Keyword ", keyword, " not in YAML file"
        if (present(value_present)) then
          value_present = .false.
        end if
        if(present(required) .and. required .eqv. .true.) then
          call error_abort("Error reading keyword " // keyword // ". Possibly missing keyword in yaml file.")
        end if
      else
        if (present(value_present)) then
          value_present = .true.
        end if
      end if

    class default
      call error_abort("Unknown type")
    end select

  end subroutine get_logical_value

  logical function parse_logic_string(str)

    character(len=*), intent(in) :: str

    if(str == 'true') then
      parse_logic_string = .true.
    else if (str == 'false') then
      parse_logic_string = .false.
    else 
      parse_logic_string = .false.
      call error_abort("The string " // str // " cannot be parsed as a logic value")
    end if

  end function parse_logic_string
 
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
    type(type_error), allocatable :: io_err

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
    type(type_error), allocatable :: io_err

    select type (config_file)
    type is (type_dictionary)

      dict => config_file%get_dictionary('reference_numbers', required=.true., error=io_err)
      if (.not. allocated(io_err)) then

        ! Pressure
        if (present(p_ref)) then
          call get_value(dict, "pressure", p_ref, required=.false.)
        end if

        ! Pressure_total
        if (present(p_total)) then
          call get_value(dict, "pressure_total", p_total, required = .false.)
        end if

        ! Temperature
        if (present(temp_ref)) then
          call get_value(dict, "temperature", temp_ref, required = .false.)
        end if

        ! Density
        if (present(dens_ref)) then
          call get_value(dict, "density", dens_ref, required = .false.)
        end if

        ! Viscosity
        if (present(visc_ref)) then
          call get_value(dict, "viscosity", visc_ref, required = .false.)
        end if

        ! Velocity
        if (present(velo_ref)) then
          call get_value(dict, "velocity", velo_ref, required = .false.)
        end if

        ! Length
        if (present(len_ref)) then
          call get_value(dict, "length", len_ref, required = .false.)
        end if

        ! Pref_at_cell
        if (present(pref_at_cell)) then
          call get_value(dict, "pref_at_cell", pref_at_cell, required = .false.)
        end if
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
  module subroutine get_solve(config_file, solved_variables)
    
    class(*), pointer, intent(in) :: config_file                      !< the entry point to the config file
    character(len=ccs_string_len), dimension(:), allocatable, intent(out) :: solved_variables

    class(*), pointer :: dict
    class(*), pointer :: dict_var
    type(type_error), allocatable :: io_err
    integer(ccs_int) :: n_var
    integer :: i
    character(len=25) :: key
    character(len=:), allocatable :: solved
    character(len=:), allocatable :: variable
    character(len=ccs_string_len) :: var
    logical :: val_present

   allocate(solved_variables(0))
    
    select type (config_file)
    type is (type_dictionary)

      dict => config_file%get_dictionary('variables', required=.true., error=io_err)
      call get_value(dict, "n_variables", n_var)

      select type (dict)
      type is(type_dictionary)
        do i = 1, n_var
          write (key, '(A, I0)') "variable_", i
          print *, i
          dict_var => dict%get_dictionary(key, required=.true., error=io_err)
          call get_value(dict_var, 'solve', solved, value_present=val_present, required=.false.)
          if (val_present) then
            if (trim(solved) == "on") then
              call get_value(dict_var, "name", variable)
              var = adjustl(variable)
              solved_variables = [solved_variables, var]
            end if
          end if
        end do
        class default
          call error_abort("Unknown type")
      end select

    class default
      call error_abort("Unknown type")
    end select

    print *, "Solved variables = ", solved_variables

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
    type(type_error), allocatable :: io_err

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
    type(type_error), allocatable :: io_err

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
    type(type_error), allocatable :: io_err

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
    type(type_error), allocatable :: io_err

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
  module subroutine get_relaxation_factors(config_file, u_relax, v_relax, p_relax, te_relax, ed_relax)
    class(*), pointer, intent(in) :: config_file        !< the entry point to the config file
    real(ccs_real), optional, intent(inout) :: u_relax  !< relaxation factor for u
    real(ccs_real), optional, intent(inout) :: v_relax  !< relaxation factor for v
    real(ccs_real), optional, intent(inout) :: p_relax  !< relaxation factor for p
    real(ccs_real), optional, intent(inout) :: te_relax !< relaxation factor for te
    real(ccs_real), optional, intent(inout) :: ed_relax !< relaxation factor for ed

    class(*), pointer :: dict
    type(type_error), allocatable :: io_err

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

  !> Get output file format
  module subroutine get_plot_format(config_file, plot_format)
    class(*), pointer, intent(in) :: config_file                !< the entry point to the config file
    character(len=:), allocatable, intent(inout) :: plot_format !< output format (e.g. vtk)

    call get_value(config_file, "plot_format", plot_format)

  end subroutine

  !> Get output type and variables
  module subroutine get_output_type(config_file, post_type, post_vars)
    class(*), pointer, intent(in) :: config_file            !< the entry point to the config file
    character(len=:), allocatable, intent(out) :: post_type !< values at cell centres or cell vertices?
    character(len=ccs_string_len), dimension(:), allocatable, intent(out) :: post_vars !< variables to be written out

    class(*), pointer :: dict
    class(type_list), pointer :: list
    class(type_list_item), pointer :: item
    type(type_error), allocatable :: io_err
    character(len=ccs_string_len) :: elt

    allocate(post_vars(0)) ! Set initial size to zero
    
    select type (config_file)
    type is (type_dictionary)

      dict => config_file%get_dictionary('post', required=.false., error=io_err)

      print *, "+++++ 1"
      if (.not. allocated(io_err)) then
        ! print *, "Reading post type"
        ! call get_value(dict, "type", post_type)
        post_type = "center"
        select type (dict)
        type is (type_dictionary)

          print *, "+++++ 2"
          list => dict%get_list('variables', required=.false., error=io_err)
          ! call error_handler(io_err)

          print *, "+++++ 2.1"
          if (.not. allocated(io_err)) then
            print *, "+++++ 2.2"
            item => list%first
            print *, "+++++ 3"
            do while (associated(item))
              select type (element => item%node)
              class is (type_scalar)
                print *, "+++++ 4"
                elt = ""
                elt = element%string
                print *, elt
                post_vars = [post_vars, elt]
                item => item%next
              end select
            end do
          else
            print *, "COULDN'T FIND POST VARIABLES"
          end if

        class default
          call error_abort("Unknown type")
        end select
      end if

    class default
      call error_abort("Unknown type")
    end select

  end subroutine get_output_type

  module subroutine get_boundary_count(filename, n_boundaries)
    character(len=*), intent(in) :: filename
    integer(ccs_int), intent(out) :: n_boundaries

    class(*), pointer :: config_file
    class(*), pointer :: dict
    character(:), allocatable :: error
    type(type_error), allocatable :: io_err

    config_file => parse(filename, error)
    if (allocated(error)) then
      call error_abort(trim(error))
    end if

    select type (config_file)
    type is (type_dictionary)
      dict => config_file%get_dictionary("boundaries", required=.true., error=io_err)
      !call error_handler(io_err)

      call get_value(dict, "n_boundaries", n_boundaries)
    class default
      call error_abort("type unhandled")
    end select
  end subroutine get_boundary_count

  module subroutine get_store_residuals(filename, store_residuals)
    character(len=*), intent(in) :: filename
    logical, intent(out) :: store_residuals

    class(*), pointer :: config_file
    class(*), pointer :: dict
    character(:), allocatable :: error
    type(type_error), allocatable :: io_err
    logical :: value_present

    config_file => parse(filename, error)
    if (allocated(error)) then
      call error_abort(trim(error))
    end if

    select type (config_file)
    type is (type_dictionary)
      dict => config_file%get_dictionary("variables", required=.true., error=io_err)
      !call error_handler(io_err)

      call get_value(dict, "store_residuals", store_residuals, value_present)
      ! do not store residuals by default
      if (.not. value_present) then
        store_residuals = .false.
      end if
    class default
      call error_abort("type unhandled")
    end select
  end subroutine get_store_residuals

  module subroutine get_enable_cell_corrections(filename, enable_cell_corrections)
    character(len=*), intent(in) :: filename
    logical, intent(out) :: enable_cell_corrections

    class(*), pointer :: config_file
    class(*), pointer :: dict
    character(:), allocatable :: error
    type(type_error), allocatable :: io_err
    logical :: value_present

    config_file => parse(filename, error)
    if (allocated(error)) then
      call error_abort(trim(error))
    end if

    select type (config_file)
    type is (type_dictionary)
      dict => config_file%get_dictionary("variables", required=.true., error=io_err)
      !call error_handler(io_err)

      call get_value(dict, "enable_cell_corrections", enable_cell_corrections, value_present)
      ! do not store residuals by default
      if (.not. value_present) then
        enable_cell_corrections = .true.
      end if
    class default
      call error_abort("type unhandled")
    end select
  end subroutine get_enable_cell_corrections

  module subroutine get_variables(config_file, variables)
    class(*), pointer, intent(in) :: config_file
    character(len=ccs_string_len), dimension(:), allocatable, intent(out) :: variables

    class(*), pointer :: dict
    class(*), pointer :: dict_var
    type(type_error), allocatable :: io_err
    integer(ccs_int) :: i
    integer(ccs_int) :: n_var
    character(len=25) :: key
    character(len=:), allocatable :: variable

    select type (config_file)
    type is (type_dictionary)
      dict => config_file%get_dictionary("variables", required=.true., error=io_err)
      ! call error_handler(io_err)

      call get_value(dict, "n_variables", n_var)
      if (allocated(variables)) then
        deallocate(variables)
      end if
      allocate (variables(n_var))

      do i = 1, n_var
        write (key, '(A, I0)') "variable_", i
        select type (dict)
        type is (type_dictionary)
          dict_var => dict%get_dictionary(key, required=.true., error=io_err)
          ! call error_handler(io_err)
          call get_value(dict_var, "name", variable)
          print *, "Found variable ", variable
          write (variables(i), '(A)') trim(variable)
        class default
          call error_abort("type unhandled")
        end select
      end do
    class default
      call error_abort("type unhandled")
    end select
  end subroutine get_variables

  module subroutine get_variable_types(config_file, variable_types)
    class(*), pointer, intent(in) :: config_file
    integer(ccs_int), dimension(:), allocatable, intent(out) :: variable_types

    ! class(*), pointer :: config_file
    class(*), pointer :: dict
    class(*), pointer :: dict_var
    type(type_error), allocatable :: io_err
    integer(ccs_int) :: i
    integer(ccs_int) :: n_var
    character(len=25) :: key
    character(len=:), allocatable :: scheme

    select type (config_file)
    type is (type_dictionary)
      dict => config_file%get_dictionary("variables", required=.true., error=io_err)
      ! call error_handler(io_err)

      call get_value(dict, "n_variables", n_var)
      allocate (variable_types(n_var))

      do i = 1, n_var
        write (key, '(A, I0)') "variable_", i
        select type (dict)
        type is (type_dictionary)
          dict_var => dict%get_dictionary(key, required=.true., error=io_err)
          ! call error_handler(io_err)
          call get_value(dict_var, "type", scheme)
          print *, scheme
          variable_types(i) = get_scheme_id(trim(scheme))
        class default
          call error_abort("type unhandled")
        end select
      end do
    class default
      call error_abort("type unhandled")
    end select
  end subroutine get_variable_types

  module subroutine get_bc_field(config_file, bc_field, phi, required)
    class(*), pointer, intent(in) :: config_file
    character(len=*), intent(in) :: bc_field
    class(field), intent(inout) :: phi
    logical, optional, intent(in) :: required

    ! local variables
    class(*), pointer :: dict
    class(*), pointer :: dict2
    integer(ccs_int) :: i
    integer(ccs_int) :: n_boundaries
    type(type_error), allocatable :: io_err
    character(len=:), allocatable :: bc_field_string
    character(len=25) :: boundary_index

    class(*), pointer :: variable_dict
    character(len=128) :: variable
    character(len=:), allocatable :: bc_type
    real(ccs_real) :: bc_value
    logical :: field_exists

    field_exists = .true.
    select type (config_file)
    type is (type_dictionary)
      dict => config_file%get_dictionary("boundaries", required=.true., error=io_err)
      ! call error_handler(io_err)

      i = 1
      n_boundaries = size(phi%bcs%ids)
      do while (i <= n_boundaries)
        write (boundary_index, '(A, I0)') "boundary_", i
        select type (dict)
        type is (type_dictionary)
          dict2 => dict%get_dictionary(boundary_index, required=.true., error=io_err)
          ! call error_handler(io_err)

          select case (bc_field)
          case ("name")
            call get_value(dict2, bc_field, bc_field_string)
            call set_bc_id(i, bc_field_string, phi%bcs)
          case ("type")
            call get_value(dict2, bc_field, bc_type, field_exists, required=required)
            if (field_exists) then
              call set_bc_type(i, bc_type, phi%bcs)
            end if
          case ("value")
            call get_value(dict2, bc_field, bc_value, field_exists, required=required)
            if (field_exists) then
              call set_bc_real_value(i, bc_value, phi%bcs)
            end if
          case default
            select type (dict2)
            type is (type_dictionary)
              write (variable, '(A, A)') "variable_", trim(bc_field)
              variable_dict => dict2%get_dictionary(trim(variable), required=.false., error=io_err)
              ! call error_handler(io_err)

              if (associated(variable_dict)) then
                call get_value(variable_dict, "type", bc_type, field_exists)
                if (field_exists) then
                  call set_bc_type(i, bc_type, phi%bcs)
                end if

                call get_value(variable_dict, "value", bc_value, field_exists)
                if (field_exists) then
                  call set_bc_real_value(i, bc_value, phi%bcs)
                end if
              end if
            end select
          end select
        end select
        i = i + 1
      end do
    class default
      call error_abort("type unhandled")
    end select
  end subroutine get_bc_field

end submodule read_config_utils
