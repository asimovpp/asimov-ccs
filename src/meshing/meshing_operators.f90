submodule (meshing) meshing_operators

  implicit none

contains

  module subroutine get_neighbour_distance(loc_p, loc_nb, dx)
    type(cell_locator), intent(in) :: loc_p
    type(neighbour_locator), intent(in) :: loc_nb
    real(accs_real), intent(out) :: dx

    real(accs_real), dimension(ndim) :: xp
    real(accs_real), dimension(ndim) :: xnb

    call get_centre(loc_p, xp)
    call get_centre(loc_nb, xnb)

    dx = sqrt(sum((xp - xnb)**2))
    
  end subroutine get_neighbour_distance
    
end submodule meshing_operators
