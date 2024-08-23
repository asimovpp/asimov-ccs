program bc_test
#include "ccs_macros.inc"

  ! ASiMoV-CCS uses
  use ccs_base, only: mesh
  use testing_lib
  use kinds, only: ccs_real, ccs_int
  use types, only: field, central_field
  use utils, only: debug_print, exit_print, str
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays
  use read_config, only: get_variables, get_boundary_count, get_boundary_names
  use constants, only: ccs_string_len
  use bc_constants

  use fortran_yaml_c_interface, only: parse

  implicit none

  class(field), allocatable :: u, v, w, p
  type(central_field), dimension(:), allocatable :: phi
  integer(ccs_int) :: i, j
  integer(ccs_int) :: n_boundaries
  character(len=*), parameter :: config_file = "test_read_bc_config.in"
  character(len=ccs_string_len), dimension(:), allocatable :: variable_names
  character(len=6) :: variable_name
  integer(ccs_int), dimension(:), allocatable :: bc_names, bc_ids, bc_types
  real(ccs_real), dimension(:), allocatable :: bc_values

  class(*), pointer :: parsed_config_file
  character(:), allocatable :: error

  character(len=128), dimension(:), allocatable :: bnd_names

    call init()

  ! Init velocities and scalar
  allocate (central_field :: u)
  allocate (central_field :: v)
  allocate (central_field :: w)
  allocate (central_field :: p)

  ! Read bc configuration
  ! First get the number of boundaries and variables
  parsed_config_file => parse(config_file, error)
  if (allocated(error)) then
    call error_abort(trim(error))
  end if
  call get_boundary_names(config_file, bnd_names)
  mesh%bnd_names = bnd_names ! Need to mock a mesh with bounday names
  print *, mesh%bnd_names
  call get_boundary_count(config_file, n_boundaries)
  call get_variables(parsed_config_file, variable_names)
  allocate (phi(size(variable_names) - 4))

  ! Check that the number of boundaries and the variables are read in correctly
  call assert_eq(n_boundaries, 4, "n_boundaries doesn't match expected value")
  call assert_eq(variable_names(1), "u", "variable name doesn't match expected value")
  call assert_eq(variable_names(2), "v", "variable name doesn't match expected value")
  call assert_eq(variable_names(3), "w", "variable name doesn't match expected value")
  call assert_eq(variable_names(4), "p", "variable name doesn't match expected value")
  do i = 5, size(variable_names)
    write (variable_name, "(A, I0)") "phi_", i - 4
    call assert_eq(variable_names(i), variable_name, "variable name doesn't match expected value")
  end do

  ! Allocate arrays for the BCs on a per-variable basis
  call allocate_bc_arrays(n_boundaries, u%bcs)
  call allocate_bc_arrays(n_boundaries, v%bcs)
  call allocate_bc_arrays(n_boundaries, w%bcs)
  call allocate_bc_arrays(n_boundaries, p%bcs)
  do i = 1, size(variable_names) - 4
    call allocate_bc_arrays(n_boundaries, phi(i)%bcs)
  end do

  ! Finally read in the BCs
  call read_bc_config(config_file, "u", u)
  call read_bc_config(config_file, "v", v)
  call read_bc_config(config_file, "w", w)
  call read_bc_config(config_file, "p", p)
  do i = 1, size(variable_names) - 4
    write (variable_name, "(A, I0)") "phi_", i
    call read_bc_config(config_file, trim(variable_name), phi(i))
  end do

  ! Check that they're read in correctly
  allocate (bc_names(n_boundaries))
  allocate (bc_ids(n_boundaries))
  allocate (bc_types(n_boundaries))
  allocate (bc_values(n_boundaries))
  bc_ids = [ 1, 2, 3, 4 ]
  bc_types = bc_type_dirichlet
  bc_values = 1
  call check_bcs(u%bcs, bc_ids, bc_types, bc_values)
  bc_types = bc_type_neumann
  bc_values = 0
  call check_bcs(v%bcs, bc_ids, bc_types, bc_values)
  bc_types(1) = bc_type_wall
  call check_bcs(w%bcs, bc_ids, bc_types, bc_values)
  bc_types = bc_type_extrapolate
  call check_bcs(p%bcs, bc_ids, bc_types)
  do j = 1, size(variable_names) - 4
    call check_bcs(phi(j)%bcs, bc_ids, bc_types)
  end do

  call fin()

contains

  subroutine check_bcs(bcs, ids, types, values)
    type(bc_config), intent(in) :: bcs
    integer(ccs_int), dimension(:), intent(in) :: ids
    integer(ccs_int), dimension(:), intent(in) :: types
    real(ccs_real), dimension(:), intent(in), optional :: values

    do i = 1, size(bcs%ids)
      call assert_eq(bcs%ids(i), ids(i), "bc id does not match.")
      call assert_eq(bcs%bc_types(i), types(i), "bc type does not match.")
      if (present(values)) then
        call assert_eq(bcs%values(i), values(i), "bc value does not match.")
      end if
    end do
  end subroutine check_bcs
end program bc_test
