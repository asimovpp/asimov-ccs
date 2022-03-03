!> @brief boundary conditions module
!
!> @details Various BC related functionality. Need to expand.

module boundary_conditions
  use types, only: bc_config
  use kinds, only: accs_int, accs_real
  
  implicit none

  private
  public :: read_bc_config

  contains

  !> @brief Wrapper for reading config file and assigning data to BC structure
  !
  !> @param[in] filename - name of config file
  !> @param[out] bcs     - boundary conditions structure
  subroutine read_bc_config(filename, bcs) 
    use yaml, only: parse, error_length
    character(len=*), intent(in) :: filename
    type(bc_config), intent(out) :: bcs

    class(*), pointer :: config_file
    character(len=error_length) :: error

    config_file => parse(filename, error=error)
    if (error/='') then
      print*,trim(error)
      stop 1
    endif

    call get_bcs(config_file, bcs)
  end subroutine read_bc_config
  
  !> @brief Assigns bc data to structure
  !
  !> @param[in] config_file - pointer to configuration file
  !> @param[out] bcs        - boundary conditions structure
  subroutine get_bcs(config_file, bcs)
    use bc_constants
    use read_config, only: get_boundaries

    class(*), pointer, intent(in) :: config_file
    type(bc_config), intent(out) :: bcs
    character(len=16), dimension(:), allocatable :: region
    character(len=16), dimension(:), allocatable :: bc_type
    real(accs_real), dimension(:,:), allocatable :: bc_data
    integer(accs_int) :: i

    call get_boundaries(config_file, region, bc_type, bc_data)

    do i = 1, size(region)
      select case(region(i))
        case("left")
          bcs%region(i) = bc_region_left
        case("right")
          bcs%region(i) = bc_region_right
        case("top")
          bcs%region(i) = bc_region_top
        case("bottom")
          bcs%region(i) = bc_region_bottom
        case default
          print *, 'invalid BC region selected'
          stop 
      end select

      select case(bc_type(i))
        case("periodic")
          bcs%bc_type(i) = bc_type_periodic
        case("sym")
          bcs%bc_type(i) = bc_type_sym
        case("dirichlet")
          bcs%bc_type(i) = bc_type_dirichlet
        case("const_grad")
          bcs%bc_type(i) = bc_type_const_grad
        case default
          print *, 'invalid BC type selected'
          stop 
      end select
    end do
    ! ALEXEI: since we are working in 2D for now we only need a subset of the vectors specified
    bcs%endpoints(:,:) = bc_data(2:,:2)
  end subroutine get_bcs
end module boundary_conditions
