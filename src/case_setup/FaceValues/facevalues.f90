!> @brief Program file for testing face-centred value functionality

program facevalues

    use kinds, only: ccs_int, ccs_real
    use parallel_types, only: parallel_environment
    use parallel, only: initialise_parallel_environment, &
                        cleanup_parallel_environment
    use types, only: vector_spec, ccs_mesh, field, face_field
    use utils, only: set_size, initialise
    use mesh_utils, only: build_square_mesh
    use vec, only: create_vector, set_vector_location
    use constants, only: face

    implicit none

    class(parallel_environment), allocatable, target :: par_env
    class(field), allocatable :: mf

    type(vector_spec) :: vec_sizes
    type(ccs_mesh) :: mesh

    ! integer(ccs_int) :: nfaces
    integer(ccs_int) :: cps = 3 ! Cells per side of the mesh

    call initialise_parallel_environment(par_env)

    ! Create a square mesh
    mesh = build_square_mesh(par_env, cps, 1.0_ccs_real)

    allocate(face_field :: mf)

    call initialise(vec_sizes)

    ! Setup vector size to store face-centred values (rather than cell-centred values)
    call set_vector_location(face, vec_sizes)

    call set_size(par_env, mesh, vec_sizes)
    call create_vector(vec_sizes, mf%values)

    ! ! View the contents of the vector
    ! call vec_view(vec_sizes, mf%values)

    call cleanup_parallel_environment(par_env)

end program facevalues
