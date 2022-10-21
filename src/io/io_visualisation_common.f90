submodule(io_visualisation) io_visualisation_common
#include "ccs_macros.inc"

  use constants, only: ndim

  implicit none

contains

  module subroutine write_solution(par_env, case_name, step, maxstep, dt, mesh, output_field_list)

    ! Arguments
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    character(len=:), allocatable, intent(in) :: case_name
    integer(ccs_int), intent(in) :: step
    integer(ccs_int), intent(in) :: maxstep
    real(ccs_real), intent(in) :: dt
    type(ccs_mesh), intent(in) :: mesh
    type(output_list), dimension(:), intent(inout) :: output_field_list

    call write_fields(par_env, case_name, step, maxstep, dt, mesh, output_field_list)

    call write_xdmf(par_env, case_name, step, maxstep, dt, mesh, output_field_list)

  end subroutine

  module subroutine write_xdmf(par_env, case_name, step, maxstep, dt, mesh, output_field_list)

    use case_config, only: write_gradients

    ! Arguments
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    character(len=:), allocatable, intent(in) :: case_name
    integer(ccs_int), intent(in) :: step
    integer(ccs_int), intent(in) :: maxstep
    real(ccs_real), intent(in) :: dt
    type(ccs_mesh), intent(in) :: mesh
    type(output_list), dimension(:), intent(inout) :: output_field_list

    ! Local variables
    character(len=:), allocatable :: xdmf_file
    character(len=:), allocatable :: sol_file
    character(len=:), allocatable :: geo_file
    integer(ccs_int), save :: ioxdmf
    integer(ccs_int), save :: step_counter = 0

    xdmf_file = case_name//'.sol.xmf'
    sol_file = case_name//'.sol.h5'
    geo_file = case_name//'.geo'

    if (step == 1) then
      if (par_env%proc_id == par_env%root) then
        ! Open file
        open(newunit=ioxdmf, file=xdmf_file, status='unknown')
  
        ! Write file contents
        write(ioxdmf, '(a)')        '<?xml version="1.0"?>'
        write(ioxdmf, '(a)')        '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
        write(ioxdmf, '(a)')        '<Xdmf Version="2.0">'
        write(ioxdmf, '(a)')        '  <Domain>'
        write(ioxdmf, '(a)')        '    <Grid Name="Unsteady" GridType="Collection" CollectionType="Temporal">'
      endif
    endif

    associate (ncel => mesh%topo%global_num_cells, &
        nvrt => mesh%topo%global_num_vertices)

      write(ioxdmf, '(a)')            '      <Grid Name="Mesh">'
      write(ioxdmf, '(a,f10.7,a)')    '        <Time Value = "',step*dt,'" />'

      ! Topology
      if (mesh%topo%vert_per_cell == 4) then
        write(ioxdmf, '(a,i0,a)')   '        <Topology Type="Quadrilateral" NumberOfElements="',ncel,'" BaseOffset="1">'
      else
        write(ioxdmf, '(a,i0,a)')   '        <Topology Type="Hexahedron" NumberOfElements="',ncel,'" BaseOffset="1">'
      endif
      write(ioxdmf, '(a,i0,1x,i0,3(a))') '          <DataItem Dimensions="',ncel,mesh%topo%vert_per_cell,'" Format="HDF">',trim(geo_file),':/Step0/cell/vertices</DataItem>'          
      write(ioxdmf, '(a)')            '        </Topology>'

      ! Geometry
      write(ioxdmf, '(a)')            '        <Geometry Type="XYZ">'
      write(ioxdmf, '(a,i0,1x,i0,3(a))') '          <DataItem Dimensions="',nvrt,ndim,'" Format="HDF">',trim(geo_file),':/Step0/vert</DataItem>'        
      write(ioxdmf, '(a)')            '        </Geometry>'

      ! Velocity vector
      write(ioxdmf, '(a)')            '        <Attribute Name="velocity" AttributeType="Vector" Center="Cell">'
      write(ioxdmf, '(a,i0,1x,i0,a,a)') '          <DataItem Dimensions="',ncel,ndim,'" ItemType="Function"', &
                                                   ' Function="JOIN($0, $1, $2)">'
      write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/u</DataItem>'
      write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/v</DataItem>'
      write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/w</DataItem>'
      write(ioxdmf, '(a)')            '          </DataItem>'
      write(ioxdmf, '(a)')            '        </Attribute>'

      ! Pressure
      write(ioxdmf, '(a)')            '        <Attribute Name="pressure" AttributeType="Scalar" Center="Cell">'
      write(ioxdmf, '(a,i0,3(a),i0,a)') '          <DataItem Dimensions="',ncel,'" Format="HDF">',trim(sol_file),':/Step',step_counter,'/p</DataItem>'          
      write(ioxdmf, '(a)')            '        </Attribute>'

      ! Kinetic Energy
      write(ioxdmf, '(a)')            '        <Attribute Name="kinetic energy" AttributeType="Scalar" Center="Cell">'
      write(ioxdmf, '(a,i0,a,a)')     '          <DataItem Dimensions="',ncel,'" ItemType="Function"', &
                                                 ' Function="0.5 * ($0*$0 + $1*$1 + $2*$2)">'
      write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/u</DataItem>'
      write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/v</DataItem>'
      write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/w</DataItem>'
      write(ioxdmf, '(a)')            '          </DataItem>'
      write(ioxdmf, '(a)')            '        </Attribute>'

      ! Enstrophy
      if (write_gradients) then
        write(ioxdmf, '(a)')            '        <Attribute Name="enstrophy" AttributeType="Scalar" Center="Cell">'
        write(ioxdmf, '(a,i0,a,a)')     '          <DataItem Dimensions="',ncel,'" ItemType="Function"', &
                                                   ' Function="0.5 * (($5-$3)*($5-$3) + ($1-$4)*($1-$4) + ($2-$0)*($2-$0))">'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/dudy</DataItem>'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/dudz</DataItem>'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/dvdx</DataItem>'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/dvdz</DataItem>'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/dwdx</DataItem>'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/dwdy</DataItem>'
        write(ioxdmf, '(a)')            '          </DataItem>'
        write(ioxdmf, '(a)')            '        </Attribute>'
      endif

      write(ioxdmf, '(a)')            '      </Grid>'

    end associate

    ! Close file
    if (step == maxstep) then
      write(ioxdmf, '(a)')          '    </Grid>'
      write(ioxdmf, '(a)')          '  </Domain>'
      write(ioxdmf, '(a)')          '</Xdmf>'
      close(ioxdmf)
    endif

    step_counter = step_counter + 1

  end subroutine write_xdmf

end submodule