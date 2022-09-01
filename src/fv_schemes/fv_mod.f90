!v Module file fv.mod
!
!  An interface to finite volume implementations (CDS, UDS, etc.)

module fv

  use kinds, only: ccs_real, ccs_int
  use types, only: ccs_matrix, ccs_vector, ccs_mesh, field, upwind_field, central_field, bc_config, face_locator

  implicit none

  private

  public :: compute_fluxes
  public :: calc_advection_coeff
  public :: calc_diffusion_coeff
  public :: calc_mass_flux
  public :: calc_cell_coords
  public :: update_gradient

  interface calc_advection_coeff
    module procedure calc_advection_coeff_cds
    module procedure calc_advection_coeff_uds
  end interface calc_advection_coeff

  interface

    !> Calculates advection coefficient for neighbouring cell using CDS discretisation
    module subroutine calc_advection_coeff_cds(phi, mf, bc, coeff)
      type(central_field), intent(in) :: phi !< scalar (central) field
      real(ccs_real), intent(in) :: mf       !< mass flux at the face
      integer(ccs_int), intent(in) :: bc     !< flag indicating whether cell is on boundary
      real(ccs_real), intent(out) :: coeff   !< advection coefficient to be calculated
    end subroutine calc_advection_coeff_cds

    !> Calculates advection coefficient for neighbouring cell using UDS discretisation
    module subroutine calc_advection_coeff_uds(phi, mf, bc, coeff)
      type(upwind_field), intent(in) :: phi !< scalar (upwind) field
      real(ccs_real), intent(in) :: mf      !< mass flux at the face
      integer(ccs_int), intent(in) :: bc    !< flag indicating whether cell is on boundary
      real(ccs_real), intent(out) :: coeff  !< advection coefficient to be calculated
    end subroutine calc_advection_coeff_uds

    !> Sets the diffusion coefficient
    ! XXX: why is this a function when the equivalent advection ones are subroutines?
    module function calc_diffusion_coeff(index_p, index_nb, mesh) result(coeff)
      integer(ccs_int), intent(in) :: index_p  !< the local cell index
      integer(ccs_int), intent(in) :: index_nb !< the local neigbouring cell index
      type(ccs_mesh), intent(in) :: mesh       !< the mesh structure
      real(ccs_real) :: coeff                  !< the diffusion coefficient
    end function calc_diffusion_coeff

    !> Computes fluxes and assign to matrix and RHS
    module subroutine compute_fluxes(phi, mf, mesh, cps, M, vec)
      class(field), intent(in) :: phi         !< scalar field structure
      class(field), intent(in) :: mf          !< mass flux field structure (defined at faces)
      type(ccs_mesh), intent(in) :: mesh      !< the mesh being used
      integer(ccs_int), intent(in) :: cps     !< the number of cells per side in the (square) mesh
      class(ccs_matrix), intent(inout) :: M   !< Data structure containing matrix to be filled
      class(ccs_vector), intent(inout) :: vec !< Data structure containing RHS vector to be filled
    end subroutine

    !> Calculates mass flux across given face. Note: assumes rho = 1 and uniform grid
    module function calc_mass_flux(u, v, w, p, dpdx, dpdy, dpdz, invAu, invAv, invAw, loc_f) result(flux)
      real(ccs_real), dimension(:), intent(in) :: u     !< x velocities
      real(ccs_real), dimension(:), intent(in) :: v     !< y velocities
      real(ccs_real), dimension(:), intent(in) :: w     !< z velocities
      real(ccs_real), dimension(:), intent(in) :: p     !< array containing pressure
      real(ccs_real), dimension(:), intent(in) :: dpdx  !< pressure gradients in x
      real(ccs_real), dimension(:), intent(in) :: dpdy  !< pressure gradients in y
      real(ccs_real), dimension(:), intent(in) :: dpdz  !< pressure gradients in z
      real(ccs_real), dimension(:), intent(in) :: invAu !< inverse momentum diagonal in x
      real(ccs_real), dimension(:), intent(in) :: invAv !< inverse momentum diagonal in y
      real(ccs_real), dimension(:), intent(in) :: invAw !< inverse momentum diagonal in z
      type(face_locator), intent(in) :: loc_f           !< face locator
      real(ccs_real) :: flux                            !< the flux across the boundary
    end function calc_mass_flux

    !> Calculates the row and column indices from flattened vector index. Assumes square mesh
    module subroutine calc_cell_coords(index, cps, row, col)
      integer(ccs_int), intent(in) :: index !< cell index
      integer(ccs_int), intent(in) :: cps   !< number of cells per side
      integer(ccs_int), intent(out) :: row  !< cell row within mesh
      integer(ccs_int), intent(out) :: col  !< cell column within mesh
    end subroutine calc_cell_coords

    !v Performs an update of the gradients of a field.
    !  @note This will perform a parallel update of the gradient fields to ensure halo cells are
    !  correctly updated on other PEs. @endnote
    module subroutine update_gradient(mesh, phi)
      type(ccs_mesh), intent(in) :: mesh !< the mesh
      class(field), intent(inout) :: phi !< the field whose gradients we want to update
    end subroutine update_gradient

  end interface

end module fv
