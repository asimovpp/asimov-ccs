!> @brief Module file initialise_fields
!
!> @details Performs initialisation of the fields describing a problem.

module initialise_fields

  implicit none

  private
  public :: initialise_field

contains
  
  !> @brief Initialise a field for the current problem
  !
  !> @details Takes an object containing the fields for the problem and some way of traversing them -
  !!          i.e. the mesh - and the name of the field to be set.
  !!          Applies corresponding initialisation function for named field, setting the initial
  !!          values.
  subroutine initialise_field()
  end subroutine initialise_field

end module initialise_fields
