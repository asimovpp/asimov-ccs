submodule (partitioning) partitioning_common

implicit none

  contains

  module subroutine read_topology(par_env, case_name, topo)

    use constants, only: geoext, adiosconfig
    use io, only: initialise_io, cleanup_io, configure_io, &
                  open_file, close_file, &
                  read_scalar, read_array
    use types, only: io_environment, io_process
    use parallel, only: read_command_line_arguments

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    character(len=:), allocatable :: case_name
    type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition

    ! Local variables
    character(len=:), allocatable :: geo_file    ! Geo file name
    character(len=:), allocatable :: adios2_file ! ADIOS2 config file name

    class(io_environment), allocatable :: io_env
    class(io_process), allocatable :: geo_reader

    integer(ccs_long), dimension(1) :: sel_start
    integer(ccs_long), dimension(1) :: sel_count

    geo_file = case_name//geoext
    adios2_file = case_name//adiosconfig

    call initialise_io(par_env, adios2_file, io_env)
    call configure_io(io_env, "geo_reader", geo_reader)  
  
    call open_file(geo_file, "read", geo_reader)
  
    ! Read attribute "ncel" - the total number of cells
    call read_scalar(geo_reader, "ncel", topo%global_num_cells)
    ! Read attribute "nfac" - the total number of faces
    call read_scalar(geo_reader, "nfac", topo%global_num_faces)
    ! Read attribute "maxfaces" - the maximum number of faces per cell
    call read_scalar(geo_reader, "maxfaces", topo%max_faces)

    allocate(topo%face_edge_end1(topo%global_num_faces))
    allocate(topo%face_edge_end2(topo%global_num_faces))

    sel_start(1) = 0 ! Global index to start reading from
    sel_count(1) = topo%global_num_faces ! How many elements to read in total

    ! Read arrays face/cell1 and face/cell2
    call read_array(geo_reader, "/face/cell1", sel_start, sel_count, topo%face_edge_end1)
    call read_array(geo_reader, "/face/cell2", sel_start, sel_count, topo%face_edge_end2)

    ! Close the file and ADIOS2 engine
    call close_file(geo_reader)

    ! Finalise the ADIOS2 IO environment
    call cleanup_io(io_env)

  end subroutine read_topology

end submodule