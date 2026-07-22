!v Submodule file fv_equations_momentum.smod
submodule(fv_equations) fv_equations_momentum

  use meshing, only: get_global_index, get_local_index, get_volume

  implicit none

contains

  module subroutine momentum_init(self, max_faces, component)
    class(momentum_equation), intent(inout) :: self
    integer(ccs_int), intent(in) :: max_faces
    integer(ccs_int), intent(in), optional :: component

    if (.not. present(component)) then
      error stop "momentum_equation%init: missing velocity component"
    end if

    call ensure_momentum_payload(self)
    call scalar_transport_init(self, max_faces, component)

    select type (payload => self%payload)
    type is (momentum_payload)
      payload%volume_p = 0.0_ccs_real
      payload%pressure_gradient_p = 0.0_ccs_real
    class default
      error stop "momentum_equation%init: invalid payload type"
    end select

  end subroutine momentum_init

  module subroutine momentum_gather_pressure_source(self, p_gradient, loc_p)
    class(momentum_equation), intent(inout) :: self
    real(ccs_real), dimension(:), intent(in) :: p_gradient
    type(cell_locator), intent(in) :: loc_p

    integer(ccs_int) :: index_p

    call ensure_momentum_payload(self)
    call get_local_index(loc_p, index_p)

    select type (payload => self%payload)
    type is (momentum_payload)
      payload%index_p = index_p
      payload%pressure_gradient_p = p_gradient(index_p)

      call get_global_index(loc_p, payload%global_index_p)
      call get_volume(loc_p, payload%volume_p)
    class default
      error stop "momentum_equation%gather_pressure_source: invalid payload type"
    end select

  end subroutine momentum_gather_pressure_source

  module subroutine momentum_apply_pressure_source(self)
    class(momentum_equation), intent(inout) :: self

    real(ccs_real) :: source

    call ensure_momentum_payload(self)

    select type (payload => self%payload)
    type is (momentum_payload)
      source = -payload%pressure_gradient_p * payload%volume_p

      call self%prepare_row(payload%global_index_p)
      call self%set_rhs(source)
    class default
      error stop "momentum_equation%apply_pressure_source: invalid payload type"
    end select

  end subroutine momentum_apply_pressure_source

  module subroutine ensure_momentum_payload(self)
    class(momentum_equation), intent(inout) :: self

    if (allocated(self%payload)) then
      select type (payload => self%payload)
      type is (momentum_payload)
        return
      class default
        deallocate (self%payload)
        self%capacity = 0_ccs_int
      end select
    end if

    allocate (momentum_payload :: self%payload)

  end subroutine ensure_momentum_payload

end submodule fv_equations_momentum
