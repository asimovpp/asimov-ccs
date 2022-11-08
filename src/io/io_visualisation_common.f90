!v Submodule file io_visualisation_common.smod
!
!  An implementation of the visualisation-related IO routines

submodule(io_visualisation) io_visualisation_common
#include "ccs_macros.inc"

  use constants, only: ndim

  implicit none

contains

  !> Write the flow solution for the current time-step to file
  module subroutine write_solution(par_env, case_name, mesh, output_list, step, maxstep, dt)

    ! Arguments
    class(parallel_environment), allocatable, target, intent(in) :: par_env  !< The parallel environment
    character(len=:), allocatable, intent(in) :: case_name                   !< The case name
    type(ccs_mesh), intent(in) :: mesh                                       !< The mesh
    type(field_ptr), dimension(:), intent(inout) :: output_list              !< List of fields to output
    integer(ccs_int), optional, intent(in) :: step                           !< The current time-step count
    integer(ccs_int), optional, intent(in) :: maxstep                        !< The maximum time-step count
    real(ccs_real), optional, intent(in) :: dt                               !< The time-step size

    ! Write the required fields ('heavy' data)
    if (present(step) .and. present(maxstep) .and. present(dt)) then
      ! Unsteady case
      call write_fields(par_env, case_name, mesh, output_list, step, maxstep, dt)
    else
      ! Steady case
      call write_fields(par_env, case_name, mesh, output_list)
    end if

    ! Write the XML descriptor ('light' data)
    if (present(step) .and. present(maxstep) .and. present(dt)) then
      ! Unsteady case
      call write_xdmf(par_env, case_name, mesh, output_list, step, maxstep, dt)
    else
      ! Steady case
      call write_xdmf(par_env, case_name, mesh, output_list)
    end if

  end subroutine

  !> Write the XML descriptor file, which describes the grid and flow data in the 'heavy' data files
  module subroutine write_xdmf(par_env, case_name, mesh, output_list, step, maxstep, dt)

    use case_config, only: write_gradients

    ! Arguments
    class(parallel_environment), allocatable, target, intent(in) :: par_env  !< The parallel environment
    character(len=:), allocatable, intent(in) :: case_name                   !< The case name
    type(ccs_mesh), intent(in) :: mesh                                       !< The mesh
    type(field_ptr), dimension(:), intent(inout) :: output_list              !< List of fields to output
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

    xdmf_file = case_name // '.sol.xmf'
    sol_file = case_name // '.sol.h5'
    geo_file = case_name // '.geo'

    ! On first call, write the header of the XML file
    if (par_env%proc_id == par_env%root) then
      if (present(step)) then
        ! Unsteady case
        if (step == 1) then
          ! Open file
          open (newunit=ioxdmf, file=xdmf_file, status='unknown')

          ! Write file contents
          write (ioxdmf, '(a)') '<?xml version="1.0"?>'
          write (ioxdmf, '(a)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
          write (ioxdmf, '(a)') '<Xdmf Version="2.0">'
          write (ioxdmf, '(a)') '  <Domain>'
          write (ioxdmf, '(a)') '    <Grid Name="Unsteady" GridType="Collection" CollectionType="Temporal">'
        end if
      else
        ! Steady case
        ! Open file
        open (newunit=ioxdmf, file=xdmf_file, status='unknown')

        ! Write file contents
        write (ioxdmf, '(a)') '<?xml version="1.0"?>'
        write (ioxdmf, '(a)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
        write (ioxdmf, '(a)') '<Xdmf Version="2.0">'
        write (ioxdmf, '(a)') '  <Domain>'
      end if

      associate (ncel => mesh%topo%global_num_cells, &
                 nvrt => mesh%topo%global_num_vertices)

        write (ioxdmf, '(a)') '      <Grid Name="Mesh">'

        if (present(step)) then
          write (ioxdmf, '(a,f10.7,a)') '        <Time Value = "', step * dt, '" />'
        end if

        ! Topology
        if (mesh%topo%vert_per_cell == 4) then
          write (ioxdmf, '(a,i0,a)') '        <Topology Type="Quadrilateral" NumberOfElements="', ncel, '" BaseOffset="1">'
        else
          write (ioxdmf, '(a,i0,a)') '        <Topology Type="Hexahedron" NumberOfElements="', ncel, '" BaseOffset="1">'
        end if

        fmt = '(a,i0,1x,i0,3(a))'
        write (ioxdmf, fmt) '          <DataItem Dimensions="', ncel, mesh%topo%vert_per_cell, '" Format="HDF">', &
          trim(geo_file), ':/Step0/cell/vertices</DataItem>'
        write (ioxdmf, '(a)') '        </Topology>'

        ! Geometry
        write (ioxdmf, '(a)') '        <Geometry Type="XYZ">'
        write (ioxdmf, fmt) '          <DataItem Dimensions="', nvrt, ndim, '" Format="HDF">', trim(geo_file), &
          ':/Step0/vert</DataItem>'
        write (ioxdmf, '(a)') '        </Geometry>'

        ! Velocity vector
        write (ioxdmf, '(a)') '        <Attribute Name="velocity" AttributeType="Vector" Center="Cell">'

        fmt = '(a,i0,1x,i0,a)'
        write (ioxdmf, fmt) '          <DataItem Dimensions="', ncel, ndim, '" ItemType="Function" Function="JOIN($0, $1, $2)">'

        fmt = '(a,i0,3(a),i0,a)'
        write (ioxdmf, fmt) '            <DataItem Format="HDF" Dimensions="', ncel, '">', trim(sol_file), ':/Step', &
          step_counter, '/u</DataItem>'
        write (ioxdmf, fmt) '            <DataItem Format="HDF" Dimensions="', ncel, '">', trim(sol_file), ':/Step', &
          step_counter, '/v</DataItem>'
        write (ioxdmf, fmt) '            <DataItem Format="HDF" Dimensions="', ncel, '">', trim(sol_file), ':/Step', &
          step_counter, '/w</DataItem>'
        write (ioxdmf, '(a)') '          </DataItem>'
        write (ioxdmf, '(a)') '        </Attribute>'

        ! Pressure
        write (ioxdmf, '(a)') '        <Attribute Name="pressure" AttributeType="Scalar" Center="Cell">'
        write (ioxdmf, fmt) '          <DataItem Dimensions="', ncel, '" Format="HDF">', trim(sol_file), ':/Step', &
          step_counter, '/p</DataItem>'
        write (ioxdmf, '(a)') '        </Attribute>'

        ! Kinetic Energy
        write (ioxdmf, '(a)') '        <Attribute Name="kinetic energy" AttributeType="Scalar" Center="Cell">'

        fmt = '(a,i0,a,a)'
        write (ioxdmf, fmt) '          <DataItem Dimensions="', ncel, '" ItemType="Function"', &
          ' Function="0.5 * ($0*$0 + $1*$1 + $2*$2)">'

        fmt = '(a,i0,3(a),i0,a)'
        write (ioxdmf, fmt) '            <DataItem Format="HDF" Dimensions="', ncel, '">', trim(sol_file), ':/Step', &
          step_counter, '/u</DataItem>'
        write (ioxdmf, fmt) '            <DataItem Format="HDF" Dimensions="', ncel, '">', trim(sol_file), ':/Step', &
          step_counter, '/v</DataItem>'
        write (ioxdmf, fmt) '            <DataItem Format="HDF" Dimensions="', ncel, '">', trim(sol_file), ':/Step', &
          step_counter, '/w</DataItem>'
        write (ioxdmf, '(a)') '          </DataItem>'
        write (ioxdmf, '(a)') '        </Attribute>'

        ! Enstrophy
        if (write_gradients) then
          write (ioxdmf, '(a)') '        <Attribute Name="enstrophy" AttributeType="Scalar" Center="Cell">'

          fmt = '(a,i0,a,a)'
          write (ioxdmf, fmt) '          <DataItem Dimensions="', ncel, '" ItemType="Function"', &
            ' Function="0.5 * (($5-$3)*($5-$3) + ($1-$4)*($1-$4) + ($2-$0)*($2-$0))">'

          fmt = '(a,i0,3(a),i0,a)'
          write (ioxdmf, fmt) '            <DataItem Format="HDF" Dimensions="', ncel, '">', trim(sol_file), ':/Step', &
            step_counter, '/dudy</DataItem>'
          write (ioxdmf, fmt) '            <DataItem Format="HDF" Dimensions="', ncel, '">', trim(sol_file), ':/Step', &
            step_counter, '/dudz</DataItem>'
          write (ioxdmf, fmt) '            <DataItem Format="HDF" Dimensions="', ncel, '">', trim(sol_file), ':/Step', &
            step_counter, '/dvdx</DataItem>'
          write (ioxdmf, fmt) '            <DataItem Format="HDF" Dimensions="', ncel, '">', trim(sol_file), ':/Step', &
            step_counter, '/dvdz</DataItem>'
          write (ioxdmf, fmt) '            <DataItem Format="HDF" Dimensions="', ncel, '">', trim(sol_file), ':/Step', &
            step_counter, '/dwdx</DataItem>'
          write (ioxdmf, fmt) '            <DataItem Format="HDF" Dimensions="', ncel, '">', trim(sol_file), ':/Step', &
            step_counter, '/dwdy</DataItem>'
          write (ioxdmf, '(a)') '          </DataItem>'
          write (ioxdmf, '(a)') '        </Attribute>'
        end if

        write (ioxdmf, '(a)') '      </Grid>'

      end associate

      flush (ioxdmf)

      ! On final call, write the closing tags and close the XML file
      if (present(step)) then
        ! Unsteady case
        if (step == maxstep) then
          write (ioxdmf, '(a)') '    </Grid>'
          write (ioxdmf, '(a)') '  </Domain>'
          write (ioxdmf, '(a)') '</Xdmf>'
          close (ioxdmf)
        end if
      else
        ! Steady case
        write (ioxdmf, '(a)') '  </Domain>'
        write (ioxdmf, '(a)') '</Xdmf>'
        close (ioxdmf)
      end if
    end if ! root

    ! Increment ADIOS2 step counter
    step_counter = step_counter + 1

  end subroutine write_xdmf

end submodule
