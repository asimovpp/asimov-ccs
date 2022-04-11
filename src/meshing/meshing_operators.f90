submodule (meshing) meshing_operators

  implicit none

contains

  !> @brief Returns the distance between two cell centres
  !
  !> @param[in]  cell_locator      loc_p  - The cell distance is measured from.
  !> @param[in]  neighbour_locator loc_nb - The cell distance is measured to.
  !> @param[out] accs_real         dx     - ndim-array of the distance
  module subroutine get_neighbour_distance(loc_p, loc_nb, dx)
    type(cell_locator), intent(in) :: loc_p
    type(neighbour_locator), intent(in) :: loc_nb
    real(accs_real), dimension(ndim), intent(out) :: dx

    real(accs_real), dimension(ndim) :: xp
    real(accs_real), dimension(ndim) :: xnb

    call get_centre(loc_p, xp)
    call get_centre(loc_nb, xnb)

    dx(:) = xnb(:) - xp(:)  
    
  end subroutine get_neighbour_distance

  !> @brief Returns the distance from cell to face centres
  !
  !> @param[in]  cell_locator loc_p - The cell distance is measured from.
  !> @param[in]  face_locator loc_f - The face distance is measured to.
  !> @param[out] accs_real    dx    - ndim-array of the distance
  module subroutine get_face_distance(loc_p, loc_f, dx)
    type(cell_locator), intent(in) :: loc_p
    type(face_locator), intent(in) :: loc_f
    real(accs_real), dimension(ndim), intent(out) :: dx

    real(accs_real), dimension(ndim) :: xp
    real(accs_real), dimension(ndim) :: xf

    call get_centre(loc_p, xp)
    call get_centre(loc_f, xf)

    dx(:) = xf(:) - xp(:)
    
  end subroutine get_face_distance
    
end submodule meshing_operators
