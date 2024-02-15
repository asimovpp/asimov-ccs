!v Submodule file io_visualisation_common.smod
!
!  An implementation of the visualisation-related IO routines
!
!  @dont_fail_linter

submodule(io_visualisation) io_visualisation_common
#include "ccs_macros.inc"

  use constants, only: ndim
  use timers, only: timer_init, timer_register_start, timer_register, timer_start, timer_stop, timer_print, timer_get_time, timer_print_all

  use types, only: field
  use utils, only: get_field
  
  implicit none

  character(len=1), dimension(4), parameter :: skip_fields = &
       [ "u", "v", "w", &
          "p" ]

  logical, save :: initial_step = .true.
    

contains

  module subroutine reset_io_visualisation

    initial_step = .true.
    call reset_io_visualisation_module()

  end subroutine

  !> Write the flow solution for the current time-step to file
  module subroutine write_solution(par_env, case_name, mesh, flow, step, maxstep, dt)

    use parallel, only: timer

    ! Arguments
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    character(len=:), allocatable, intent(in) :: case_name                  !< The case name
    type(ccs_mesh), intent(in) :: mesh                                      !< The mesh
    type(fluid), intent(inout) :: flow                                       !< The flow variables
    integer(ccs_int), optional, intent(in) :: step                          !< The current time-step count
    integer(ccs_int), optional, intent(in) :: maxstep                       !< The maximum time-step count
    real(ccs_real), optional, intent(in) :: dt                              !< The time-step size
    integer(ccs_int) :: timer_index_write_field
    integer(ccs_int) :: timer_index_write_xdmf

    call timer_register("Write fields time", timer_index_write_field)
    call timer_register("Write xdmf time", timer_index_write_xdmf)

    ! Write the required fields ('heavy' data)
    if (present(step) .and. present(maxstep)) then
      ! Unsteady case
      call timer_start(timer_index_write_field)
      call write_fields(par_env, case_name, mesh, flow, step, maxstep)
      call timer_stop(timer_index_write_field)
    else
      ! Steady case
      call timer_start(timer_index_write_field)
      call write_fields(par_env, case_name, mesh, flow)
      call timer_stop(timer_index_write_field)
    end if

    ! Write the XML descriptor ('light' data)
    if (present(step) .and. present(maxstep) .and. present(dt)) then
      ! Unsteady case
      call timer_start(timer_index_write_xdmf)
      call write_xdmf(par_env, case_name, flow, step, maxstep, dt)
      call timer_stop(timer_index_write_xdmf)
    else
      ! Steady case
      call timer_start(timer_index_write_xdmf)
      call write_xdmf(par_env, case_name, flow)
      call timer_stop(timer_index_write_xdmf)
    end if

  end subroutine

  !> Write the XML descriptor file, which describes the grid and flow data in the 'heavy' data files
  module subroutine write_xdmf(par_env, case_name, flow, step, maxstep, dt)

    use case_config, only: write_gradients
    use meshing, only: get_global_num_cells, get_vert_per_cell, get_global_num_vertices, get_mesh_generated

    ! Arguments
    class(parallel_environment), allocatable, target, intent(in) :: par_env  !< The parallel environment
    character(len=:), allocatable, intent(in) :: case_name                   !< The case name
    type(fluid), intent(inout) :: flow                                       !< The flow variables
    integer(ccs_int), optional, intent(in) :: step                           !< The current time-step count
    integer(ccs_int), optional, intent(in) :: maxstep                        !< The maximum time-step count
    real(ccs_real), optional, intent(in) :: dt                               !< The time-step size

    ! Local variables
    character(len=:), allocatable :: xdmf_file   ! Name of the XDMF (XML) file
    character(len=:), allocatable :: sol_file    ! Name of the solution file
    character(len=:), allocatable :: geo_file    ! Name of the mesh file
    character(len=50) :: fmt                     ! Format string
    integer(ccs_int), save :: ioxdmf             ! IO unit of the XDMF file
    integer(ccs_int), save :: step_counter = 0   ! ADIOS2 step counter
    integer(ccs_int) :: num_vel_cmp             ! Number of velocity components in output field list
    integer(ccs_int) :: i                       ! Loop counter

    character(len=2), parameter :: l1 = '  '           ! Indentation level 1
    character(len=4), parameter :: l2 = '    '         ! Indentation level 2
    character(len=6), parameter :: l3 = '      '       ! Indentation level 3
    character(len=8), parameter :: l4 = '        '     ! Indentation level 4
    character(len=10), parameter :: l5 = '          '   ! Indentation level 5
    character(len=12), parameter :: l6 = '            ' ! Indentation level 6

    integer(ccs_int) :: ncel
    integer(ccs_int) :: vert_per_cell
    integer(ccs_int) :: nvrt

    logical :: is_generated
    character(len=:), allocatable :: mesh_data_root ! Where in the mesh data path is topology/geometry stored?

    class(field), pointer :: phi
    
    xdmf_file = case_name // '.sol.xmf'
    sol_file = case_name // '.sol.h5'
    geo_file = case_name // '.geo'

    ! On first call, write the header of the XML file
    if (par_env%proc_id == par_env%root) then
      if (present(step)) then
        ! Unsteady case
        if (initial_step) then
          ! Open file
          open (newunit=ioxdmf, file=xdmf_file, status='unknown')

          ! Write file contents
          write (ioxdmf, '(a)') '<?xml version = "1.0"?>'
          write (ioxdmf, '(a)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
          write (ioxdmf, '(a)') '<Xdmf Version = "2.0">'
          write (ioxdmf, '(a,a)') l1, '<Domain>'
          write (ioxdmf, '(a,a)') l2, '<Grid Name = "Unsteady" GridType = "Collection" CollectionType = "Temporal">'

          initial_step = .false.
        end if
      else
        ! Steady case
        ! Open file
        open (newunit=ioxdmf, file=xdmf_file, status='unknown')

        ! Write file contents
        write (ioxdmf, '(a)') '<?xml version = "1.0"?>'
        write (ioxdmf, '(a)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
        write (ioxdmf, '(a)') '<Xdmf Version = "2.0">'
        write (ioxdmf, '(a,a)') l1, '<Domain>'
      end if

      call get_global_num_cells(ncel)
      call get_vert_per_cell(vert_per_cell)
      call get_global_num_vertices(nvrt)

      write (ioxdmf, '(a,a)') l3, '<Grid Name = "Mesh">'

      if (present(step)) then
        write (ioxdmf, '(a,a,f10.7,a)') l4, '<Time Value = "', step * dt, '" />'
      end if

      ! Check whether mesh was read or generated and set data path root appropriately
      call get_mesh_generated(is_generated)
      if (is_generated) then
        mesh_data_root = "/Step0"
      else
        mesh_data_root = ""
      end if

      ! Topology
      if (vert_per_cell == 4) then
        write (ioxdmf, '(a,a,i0,a)') l4, '<Topology Type = "Quadrilateral" NumberOfElements = "', ncel, '" BaseOffset = "1">'
      else
        write (ioxdmf, '(a,a,i0,a)') l4, '<Topology Type = "Hexahedron" NumberOfElements = "', ncel, '" BaseOffset = "1">'
      end if

      fmt = '(a,a,i0,1x,i0,3(a))'
      write (ioxdmf, fmt) l5, '<DataItem Dimensions = "', ncel, vert_per_cell, '" Format = "HDF">', &
        trim(geo_file), ':' // mesh_data_root // '/cell/vertices</DataItem>'
      write (ioxdmf, '(a,a)') l4, '</Topology>'

      ! Geometry
      write (ioxdmf, '(a,a)') l4, '<Geometry Type = "XYZ">'
      write (ioxdmf, fmt) l5, '<DataItem Dimensions = "', nvrt, ndim, '" Format = "HDF">', trim(geo_file), &
        ':' // mesh_data_root // '/vert</DataItem>'
      write (ioxdmf, '(a,a)') l4, '</Geometry>'

      ! Velocity vector
      ! Count number of velocity components in list of fields to be written out
      num_vel_cmp = 0
      do i = 1, size(flow%fields)
        call get_field(flow, i, phi)
        if (phi%output) then
          if ((trim(phi%name) == 'u') .or. (trim(phi%name) == 'v') .or. (trim(phi%name) == 'w')) then
            num_vel_cmp = num_vel_cmp + 1
          end if
        end if
      end do

      if (num_vel_cmp > 0) then
        write (ioxdmf, '(a,a)') l4, '<Attribute Name = "velocity" AttributeType = "Vector" Center = "Cell">'

        fmt = '(a,a,i0,1x,i0,a)'
        if (num_vel_cmp == 1) then
          write (ioxdmf, fmt) l5, '<DataItem Dimensions = "', ncel, 1, '" ItemType = "Function" Function = "JOIN($0)">'
        else if (num_vel_cmp == 2) then
          write (ioxdmf, fmt) l5, '<DataItem Dimensions = "', ncel, 2, '" ItemType = "Function" Function = "JOIN($0, $1)">'
        else if (num_vel_cmp == 3) then
          write (ioxdmf, fmt) l5, '<DataItem Dimensions = "', ncel, 3, '" ItemType = "Function" Function = "JOIN($0, $1, $2)">'
        end if

        fmt = '(a,a,i0,3(a),i0,a)'

        do i = 1, size(flow%fields)
          call get_field(flow, i, phi)
          if (phi%output) then
            if ((trim(phi%name) == 'u') .or. (trim(phi%name) == 'v') .or. (trim(phi%name) == 'w')) then
              write (ioxdmf, fmt) l6, '<DataItem Format = "HDF" Dimensions = "', ncel, '">', trim(sol_file), ':/Step', &
                   step_counter, '/'//trim(phi%name)//'</DataItem>'
            end if
          end if
        end do
        write (ioxdmf, '(a,a)') l5, '</DataItem>'
        write (ioxdmf, '(a,a)') l4, '</Attribute>'
      end if

      ! Pressure
      fmt = '(a,a,i0,3(a),i0,a)'
      do i = 1, size(flow%fields)
        call get_field(flow, i, phi)
        if (trim(phi%name) == 'p') then
          write (ioxdmf, '(a,a)') l4, '<Attribute Name = "pressure" AttributeType = "Scalar" Center = "Cell">'
          write (ioxdmf, fmt) l5, '<DataItem Dimensions = "', ncel, '" Format = "HDF">', trim(sol_file), ':/Step', &
            step_counter, '/p</DataItem>'
          write (ioxdmf, '(a,a)') l4, '</Attribute>'
        end if
      end do

      ! Scalars
      fmt = '(a,a,i0,3(a),i0,a)'
      do i = 1, size(flow%fields)
        call get_field(flow, i, phi)
        if (phi%output) then
          if (.not. any(skip_fields == trim(phi%name))) then
            write (ioxdmf, '(a,a)') l4, '<Attribute Name = "'//phi%name//'" AttributeType = "Scalar" Center = "Cell">'
            write (ioxdmf, fmt) l5, '<DataItem Dimensions = "', ncel, '" Format = "HDF">', trim(sol_file), ':/Step', &
                 step_counter, '/'//trim(phi%name)//'</DataItem>'
            write (ioxdmf, '(a,a)') l4, '</Attribute>'
          end if
        end if
      end do

      ! Kinetic Energy
      if (num_vel_cmp > 0) then
        write (ioxdmf, '(a,a)') l4, '<Attribute Name = "kinetic energy" AttributeType = "Scalar" Center = "Cell">'

        fmt = '(a,a,i0,a,a)'

        if (num_vel_cmp == 1) then
          write (ioxdmf, fmt) l5, '<DataItem Dimensions = "', ncel, '" ItemType = "Function"', &
            ' Function = "0.5 * ($0*$0)">'
        else if (num_vel_cmp == 2) then
          write (ioxdmf, fmt) l5, '<DataItem Dimensions = "', ncel, '" ItemType = "Function"', &
            ' Function = "0.5 * ($0*$0 + $1*$1)">'
        else if (num_vel_cmp == 3) then
          write (ioxdmf, fmt) l5, '<DataItem Dimensions = "', ncel, '" ItemType = "Function"', &
            ' Function = "0.5 * ($0*$0 + $1*$1 + $2*$2)">'
        end if

        fmt = '(a,a,i0,3(a),i0,a)'

        do i = 1, size(flow%fields)
          call get_field(flow, i, phi)
          if (phi%output) then
            if ((trim(phi%name) == 'u') .or. (trim(phi%name) == 'v') .or. (trim(phi%name) == 'w')) then
              write (ioxdmf, fmt) l6, '<DataItem Format = "HDF" Dimensions = "', ncel, '">', trim(sol_file), ':/Step', &
                   step_counter, '/'//trim(phi%name)//'</DataItem>'
            end if
          end if
        end do
        write (ioxdmf, '(a,a)') l5, '</DataItem>'
        write (ioxdmf, '(a,a)') l4, '</Attribute>'
      end if

      ! Enstrophy
      if (write_gradients .and. (num_vel_cmp == 3)) then
        write (ioxdmf, '(a,a)') l4, '<Attribute Name = "enstrophy" AttributeType = "Scalar" Center = "Cell">'

        fmt = '(a,a,i0,a,a)'
        write (ioxdmf, fmt) l5, '<DataItem Dimensions = "', ncel, '" ItemType = "Function"', &
          ' Function = "0.5 * (($5-$3)*($5-$3) + ($1-$4)*($1-$4) + ($2-$0)*($2-$0))">'

        fmt = '(a,a,i0,3(a),i0,a)'
        write (ioxdmf, fmt) l6, '<DataItem Format = "HDF" Dimensions = "', ncel, '">', trim(sol_file), ':/Step', &
          step_counter, '/dudy</DataItem>'
        write (ioxdmf, fmt) l6, '<DataItem Format = "HDF" Dimensions = "', ncel, '">', trim(sol_file), ':/Step', &
          step_counter, '/dudz</DataItem>'
        write (ioxdmf, fmt) l6, '<DataItem Format = "HDF" Dimensions = "', ncel, '">', trim(sol_file), ':/Step', &
          step_counter, '/dvdx</DataItem>'
        write (ioxdmf, fmt) l6, '<DataItem Format = "HDF" Dimensions = "', ncel, '">', trim(sol_file), ':/Step', &
          step_counter, '/dvdz</DataItem>'
        write (ioxdmf, fmt) l6, '<DataItem Format = "HDF" Dimensions = "', ncel, '">', trim(sol_file), ':/Step', &
          step_counter, '/dwdx</DataItem>'
        write (ioxdmf, fmt) l6, '<DataItem Format = "HDF" Dimensions = "', ncel, '">', trim(sol_file), ':/Step', &
          step_counter, '/dwdy</DataItem>'
        write (ioxdmf, '(a,a)') l5, '</DataItem>'
        write (ioxdmf, '(a,a)') l4, '</Attribute>'
      end if

      ! Residuals
      do i = 1, size(flow%fields)
        call get_field(flow, i, phi)
        if (allocated(phi%residuals)) then
          write (ioxdmf, '(a,a)') l4, '<Attribute Name = "residuals_' // trim(phi%name) // '" AttributeType = "Scalar" Center = "Cell">'
          write (ioxdmf, fmt) l5, '<DataItem Dimensions = "', ncel, '" Format = "HDF">', trim(sol_file), ':/Step', &
            step_counter, '/' // trim(phi%name) // '_res</DataItem>'
          write (ioxdmf, '(a,a)') l4, '</Attribute>'
        end if
      end do

      write (ioxdmf, '(a,a)') l3, '</Grid>'

      flush (ioxdmf)

      ! On final call, write the closing tags and close the XML file
      if (present(step)) then
        ! Unsteady case
        if (step == maxstep) then
          write (ioxdmf, '(a,a)') l2, '</Grid>'
          write (ioxdmf, '(a,a)') l1, '</Domain>'
          write (ioxdmf, '(a)') '</Xdmf>'
          close (ioxdmf)
        end if
      else
        ! Steady case
        write (ioxdmf, '(a,a)') l1, '</Domain>'
        write (ioxdmf, '(a)') '</Xdmf>'
        close (ioxdmf)
      end if
    end if ! root

    ! Increment ADIOS2 step counter
    step_counter = step_counter + 1

  end subroutine write_xdmf

end submodule
