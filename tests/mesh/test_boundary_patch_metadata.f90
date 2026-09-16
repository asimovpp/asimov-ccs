!> @brief Verify durable boundary-patch metadata from configuration through BC translation.
!> The input uses Slip_Wall, In-let, OUT let, and inlet to check name normalisation, fallback
!> from absent or unrecognised patch types, and precedence of explicit wall and generic types.
!> Direct classification checks also cover compound names, inlet/outlet aliases, ASCII whitespace,
!> and matching slipwall before wall. A generated mesh must retain this metadata after per-field
!> numerical boundary conditions are translated.
program test_boundary_patch_metadata

  use testing_lib
  use boundary_conditions, only: allocate_bc_arrays, translate_bcs_phi
  use bc_constants, only: bc_type_neumann, bc_type_wall
  use core, only: ccs_options
  use kinds, only: ccs_real
  use mesh_utils, only: build_square_mesh
  use read_config, only: get_boundary_names, get_boundary_types
  use types, only: boundary_patch, ccs_mesh, field, &
                   classify_boundary_patch_type, patch_type_generic, &
                   patch_type_inflow, patch_type_outflow, patch_type_slipwall, &
                   patch_type_wall
  use fortran_yaml_c_interface, only: parse

  implicit none

  character(len=128), dimension(:), allocatable :: bnd_names
  character(len=128), dimension(:), allocatable :: bnd_types
  character(len=:), allocatable :: error
  class(*), pointer :: config_file
  real(ccs_real), dimension(3, 4) :: bnd_normals
  type(boundary_patch), dimension(:), allocatable :: patches
  type(ccs_mesh) :: generated_mesh
  type(ccs_options) :: run_options
  type(field) :: phi

  call init()

  config_file => parse('test_boundary_patch_metadata.in', error)
  if (allocated(error)) then
    error stop trim(error)
  end if
  call get_boundary_names(config_file, bnd_names)
  call get_boundary_types(config_file, bnd_types)

  call assert_eq(bnd_names(1), 'Slip_Wall', 'Expected boundary name')
  call assert_eq(bnd_types(1), 'wall', 'Expected top-level boundary type')
  call assert_eq(bnd_types(2), 'neumann', 'Expected numerical boundary type')
  call assert_eq(bnd_types(3), '', 'Expected absent boundary type to be empty')

  call assert_eq(classify_boundary_patch_type('bottomwall'), 'wall', &
                 'Expected compound bottomwall name to classify as wall')
  call assert_eq(classify_boundary_patch_type('walls_pilot'), 'wall', &
                 'Expected compound walls_pilot name to classify as wall')
  call assert_eq(classify_boundary_patch_type('inlet'), 'inflow', &
                 'Expected inlet alias to classify as inflow')
  call assert_eq(classify_boundary_patch_type('outlet'), 'outflow', &
                 'Expected outlet alias to classify as outflow')
  call assert_eq(classify_boundary_patch_type('inlet', 'neumann'), 'inflow', &
                 'Expected invalid configured type to retain name fallback')
  call assert_eq(classify_boundary_patch_type('inlet', 'wall'), 'wall', &
                 'Expected configured type to take precedence')
  call assert_eq(classify_boundary_patch_type('SLIP-WALL', 'dirichlet'), 'slipwall', &
                 'Expected slipwall to precede wall after normalization')
  call assert_eq(classify_boundary_patch_type('slip' // new_line('a') // achar(9) &
                                               // achar(11) // achar(12) // achar(13) // 'wall', 'dirichlet'), &
                 'slipwall', 'Expected all ASCII whitespace to be normalized')
  call assert_eq(classify_boundary_patch_type('unknown', 'slipwall'), &
                 'slipwall', 'Expected configured type to take precedence')

  run_options%mesh%bnd_names = bnd_names
  run_options%mesh%bnd_types = bnd_types
  run_options%mesh%cps = 2
  run_options%mesh%domain_size = 1.0_ccs_real
  generated_mesh = build_square_mesh(par_env, shared_env, run_options)
  patches = generated_mesh%boundary_patches

  call assert_eq(patches(1)%name, 'Slip_Wall', 'Expected stable patch name')
  call assert_eq(patches(1)%patch_type, 'wall', &
                 'Expected configured type to override name')
  call assert_eq(patches(2)%patch_type, 'inflow', &
                 'Expected inlet alias to classify as inflow')
  call assert_eq(patches(3)%patch_type, 'outflow', &
                 'Expected outlet alias to classify as outflow')
  call assert_eq(patches(4)%patch_type, 'generic', &
                 'Expected explicit generic type to override inlet alias')

  call allocate_bc_arrays(4, phi%bcs)
  phi%name = 'u'
  phi%bcs%bc_types = bc_type_wall
  phi%bcs%values = 0.0_ccs_real
  bnd_normals = 0.0_ccs_real
  bnd_normals(1, :) = 1.0_ccs_real
  call translate_bcs_phi(bnd_normals, phi)

  call assert_eq(phi%bcs%bc_types(1), bc_type_neumann, &
                 'Expected field BC translation to run')
  call assert_eq(generated_mesh%boundary_patches(1)%patch_type, 'wall', &
                 'Expected mesh patch metadata to survive BC translation')

  call fin()

end program test_boundary_patch_metadata
