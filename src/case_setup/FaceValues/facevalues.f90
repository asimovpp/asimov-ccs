!> @brief Program file for testing face-centred value functionality

program facevalues

    use kinds, only: accs_int, accs_real
    use parallel_types, only: parallel_environment
    use parallel, only: initialise_parallel_environment, &
                        cleanup_parallel_environment
    use types, only: vector_init_data, mesh, field, face_field
    use utils, only: set_global_size, initialise
    use mesh_utils, only: build_square_mesh
    use vec, only: create_vector, set_vector_location
    use constants, only: face

    implicit none

    class(parallel_environment), allocatable, target :: par_env
    class(field), allocatable :: mf

    type(vector_init_data) :: vec_sizes
    type(mesh) :: square_mesh

    ! integer(accs_int) :: nfaces
    integer(accs_int) :: cps = 3 ! Cells per side of the mesh

    call initialise_parallel_environment(par_env)

    ! Create a square mesh
    square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)

    allocate(face_field :: mf)

    call initialise(vec_sizes)

    ! Setup vector size to store face-centred values (rather than cell-centred values)
    call set_vector_location(vec_sizes, face)

    call set_global_size(vec_sizes, square_mesh, par_env)
    call create_vector(vec_sizes, mf%vec)

    ! ! View the contents of the vector
    ! call vec_view(vec_sizes, mf%vec)

    call cleanup_parallel_environment(par_env)

end program facevalues
