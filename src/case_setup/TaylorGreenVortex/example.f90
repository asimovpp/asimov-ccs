
program example

  use yaml, only: parse, error_length
  use yaml_types, only: type_node, type_dictionary, type_error, real_kind, &
                        type_list, type_list_item
  
  class(type_node), pointer :: root
  character(len=error_length) :: error
  
  class(type_dictionary), pointer :: dict
  class(type_list), pointer :: list
  class(type_list_item), pointer :: item
  type(type_error), pointer :: io_err
  
  character(len=:), allocatable :: string

  character(len=:), allocatable :: title

  ! Number of steps
  integer :: steps

  ! Initialisation
  character(len=5) :: init = "     "

  ! Reference numbers
  integer :: pressure
  integer :: temperature
  integer :: density
  integer :: pref_at_cell
  real(real_kind)    :: viscosity

  ! Solve
  character(len=3) :: u_sol = "   "
  character(len=3) :: v_sol = "   "
  character(len=3) :: w_sol = "   "
  character(len=3) :: p_sol = "   "

  ! Solvers
  integer :: u_solver
  integer :: v_solver
  integer :: w_solver
  integer :: p_solver

  ! Unsteady solution
  character(len=5) :: transient_type = "     "  
  real(real_kind) :: dt
  real(real_kind) :: gamma
  integer :: max_sub_steps

  ! Target residual
  real(real_kind) :: residual

  ! Monitor cell
  integer :: monitor_cell

  root => parse("./tgv_config.yaml", error = error)
  if (error/='') then
    print*,trim(error)
    stop 1
  endif
  
  select type (root)
  class is (type_dictionary)

    ! Get title
    call get_case_name(root, title)
    
    ! Get reference numbers
    call get_reference_numbers(root, pressure, temperature, density, viscosity, pref_at_cell)

    ! Get steps
    call get_steps(root, steps)

    ! Variables to solve
    call get_solve(root, u_sol, v_sol, w_sol, p_sol)

    ! Solvers
    call get_solvers(root, u_solver, v_solver, w_solver, p_solver)  
    
    ! Get initilisation
    call get_init(root, init)

    ! Get unsteady solution parameters
    call get_transient(root, transient_type, dt, gamma, max_sub_steps)

    ! Get target resdiual
    call get_target_residual(root, residual)

    ! Get monitor cell ID
    call get_monitor_cell(root, monitor_cell)

  end select

  call root%finalize()
  
  deallocate(root)


contains

  subroutine error_handler(io_err)
    type (type_error), pointer :: io_err
    if (associated(io_err)) then
      print*,trim(io_err%message)
      !stop 1
    endif
  end subroutine

  subroutine get_case_name(root, title)
    class(type_dictionary), pointer, intent(in) :: root
    character(len=:), allocatable, intent(inout) :: title

    title = trim(root%get_string('title',error=io_err))
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'title=',title

  end subroutine

  subroutine get_steps(root, steps)
    class(type_dictionary), pointer, intent(in) :: root
    integer :: steps

    steps = root%get_integer('steps',error=io_err)
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'steps=',steps

  end subroutine

  subroutine get_init(root, init)
    class(type_dictionary), pointer, intent(in) :: root
    character(len=5), intent(inout) :: init

    init = trim(root%get_string('init',error=io_err))
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'init=',init

  end subroutine


  subroutine get_reference_numbers(root, p, temp, density, viscosity, pref)
    class(type_dictionary), pointer, intent(in) :: root
    integer, intent(inout) :: p
    integer, intent(inout) :: temp
    integer, intent(inout) :: density
    integer, intent(inout) :: pref
    real(real_kind), intent(inout) :: viscosity

    class(type_dictionary), pointer :: dict

    dict => root%get_dictionary('reference_numbers',required=.true.,error=io_err)

    print*,"Reference numbers: "
    
    ! Pressure
    pressure = dict%get_integer("pressure",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"pressure= ",pressure
    
    ! Temperature
    temperature = dict%get_integer("temperature",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"temperature= ",temperature

    ! Density
    density = dict%get_integer("density",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"density= ",density
    
    ! Pref_at_cell
    pref_at_cell = dict%get_integer("pref_at_cell",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"pref_at_cell= ",pref_at_cell

    ! Viscosity
    viscosity = dict%get_real("viscosity",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"viscosity= ",viscosity
  
  end subroutine

  subroutine get_solve(root, u_sol, v_sol, w_sol, p_sol)

    class(type_dictionary), pointer, intent(in) :: root
    character(len=3), intent(inout) :: u_sol
    character(len=3), intent(inout) :: v_sol
    character(len=3), intent(inout) :: w_sol
    character(len=3), intent(inout) :: p_sol

    class(type_dictionary), pointer :: dict

    dict => root%get_dictionary('solve',required=.true.,error=io_err)

    print*,"Solve (on/off): "

    ! Solve u?
    u_sol = trim(dict%get_string('u',error=io_err))
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,'u=',u_sol
    
    ! Solve v?
    v_sol = trim(dict%get_string('v',error=io_err))
    call error_handler(io_err)    
    if(associated(io_err) == .false.) print*,'v=',v_sol

    ! Solve w?
    w_sol = trim(dict%get_string('w',error=io_err))
    call error_handler(io_err)   
    if(associated(io_err) == .false.) print*,'w=',w_sol
    
    ! Solve p?
    p_sol = trim(dict%get_string('p',error=io_err))
    call error_handler(io_err)    
    if(associated(io_err) == .false.) print*,'p=',p_sol

  end subroutine

  subroutine get_solvers(root, u_solver, v_solver, w_solver, p_solver)

    class(type_dictionary), pointer, intent(in) :: root
    integer, intent(inout) :: u_solver
    integer, intent(inout) :: v_solver
    integer, intent(inout) :: w_solver
    integer, intent(inout) :: p_solver

    class(type_dictionary), pointer :: dict

    dict => root%get_dictionary('solver',required=.true.,error=io_err)

    print*,"Solver: "

    ! Get u_solver
    u_solver = dict%get_integer('u',error=io_err)
    call error_handler(io_err)   
    if(associated(io_err) == .false.) print*,'u=',u_solver
    
    ! Get v_solver
    v_solver = dict%get_integer('u',error=io_err)
    call error_handler(io_err)    
    if(associated(io_err) == .false.) print*,'v=',v_solver

    ! Get w_solver
    w_solver = dict%get_integer('w',error=io_err)
    call error_handler(io_err)   
    if(associated(io_err) == .false.) print*,'w=',w_solver
    
    ! Get p_solver
    p_solver = dict%get_integer('p',error=io_err)
    call error_handler(io_err)    
    if(associated(io_err) == .false.) print*,'p=',p_solver

  end subroutine

  subroutine get_transient(root, transient_type, dt, gamma, max_sub_steps)
    class(type_dictionary), pointer, intent(in) :: root
    character(len=5), intent(inout) :: transient_type
    real(real_kind), intent(inout) :: dt
    real(real_kind), intent(inout) :: gamma
    integer, intent(inout) :: max_sub_steps

    class(type_dictionary), pointer :: dict

    dict => root%get_dictionary('transient',required=.false.,error=io_err)

    print*,"Transient: "
    
    ! Transient type (euler/quad)
    transient_type = dict%get_string("type",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"type= ",transient_type
    
    ! Dt
    dt = dict%get_real("dt",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"dt= ",dt

    ! Gamma
    gamma = dict%get_integer("gamma",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"gamma= ",gamma
    
    ! Maximum number of sub steps
    max_sub_steps = dict%get_integer("max_sub_steps",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"max_sub_steps= ",max_sub_steps

  end subroutine

  subroutine get_target_residual(root, residual)
    class(type_dictionary), pointer, intent(in) :: root
    real(real_kind), intent(inout) :: residual

    residual = root%get_real('target_residual',error=io_err)
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'target_residual=',residual

  end subroutine

  subroutine get_monitor_cell(root, monitor_cell)
    class(type_dictionary), pointer, intent(in) :: root
    integer, intent(inout) :: monitor_cell

    monitor_cell = root%get_integer('monitor_cell',error=io_err)
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'monitor_cell=',monitor_cell

  end subroutine

end program


