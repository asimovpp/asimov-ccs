!> @brief Module file petsctypes.mod
!> @build petsc
!
!> @details Provides petsc-extended types.

module petsctypes

  use petscksp, only : tKSP
  use petscvec, only : tVec
  use petscmat, only : tMat

  use kinds, only : accs_err, accs_int
  use types, only : ccs_vector, matrix, linear_solver
  
  implicit none

  private

  !> @brief Implements the vector class backed by a PETSc vector
  type, public, extends(ccs_vector) :: vector_petsc
    type(tVec) :: v      !> The PETSc vector
    type(tVec) :: vl     !> The "local" PETSc vector (inc. ghost points)
    logical :: allocated !> Indicates whether the PETSc vector has been allocated
    logical :: ghosted   !> Does this vector have ghost points?
    integer(accs_int) :: mode !> Current mode for setting values
    logical :: modeset        !> Is the current mode still valid? i.e. does vector need updated before switching modes?
  contains
    final :: free_vector_petsc
  end type vector_petsc

  !> @brief Implements the matrix class backed by a PETSc matrix
  type, public, extends(matrix) :: matrix_petsc
     type(tMat) :: M      !> The PETSc matrix
     logical :: allocated !> Indicates whether the PETSc matrix has been allocated
    integer(accs_int) :: mode !> Current mode for setting values
    logical :: modeset        !> Is the current mode still valid? i.e. does matrix need updated before switching modes?
   contains
     final :: free_matrix_petsc
  end type matrix_petsc

  type, public, extends(linear_solver) :: linear_solver_petsc
     type(tKSP) :: KSP !> The PETSc solver
     logical :: allocated
   contains
     final :: free_linear_solver_petsc
  end type linear_solver_petsc

  interface
    module subroutine free_vector_petsc(v)
      type(vector_petsc), intent(inout) :: v
    end subroutine

    module subroutine free_matrix_petsc(M)
      type(matrix_petsc), intent(inout) :: M
    end subroutine

    module subroutine free_linear_solver_petsc(solver)
      type(linear_solver_petsc), intent(inout) :: solver
    end subroutine

   end interface

contains
  
  !> @brief Destroys a PETSc-backed vector.
  !
  !> @details Destructor called by deallocating a vector_petsc - confirms the PETSc vector object is
  !!          allocated and calls the necessary destructor on the wrapped PETSc vector object, sets
  !!          the allocated flag to .false. to prevent double free's.
  !> @param[in/out] vector v - the vector to be destroyed.
  module subroutine free_vector_petsc(v)
    
    use petscvec, only: VecDestroy
    
    type(vector_petsc), intent(inout) :: v

    integer(accs_err) :: ierr !> Error code

    if (v%allocated) then
       call VecDestroy(v%v, ierr)
       v%allocated = .false.
    else
       print *, "WARNING: attempted double free of vector"
    end if
    
  end subroutine
  
  !> @brief Destroys a PETSc-backed matrix.
  !
  !> @details Destructor called by deallocating a matrix_petsc - confirms the PETSc matrix object is
  !!          allocated and calls the necessary destructor on the wrapped PETSc matrix object, sets
  !!          the allocated flag to .false. to prevent double free's.
  !> @param[in/out] matrix M - the matrix to be destroyed.
  module subroutine free_matrix_petsc(M)
    
   use petscmat, only: MatDestroy

    type(matrix_petsc), intent(inout) :: M

    integer(accs_err) :: ierr !> Error code

    if (M%allocated) then
       call MatDestroy(M%M, ierr)
       M%allocated = .false.
    else
       print *, "WARNING: attempted double free of matrix"
    end if

  end subroutine

  !> @brief Destroys a PETSc-backed linear solver.
  !
  !> @details Destructor called by deallocating a linear_solver_petsc - confirms the PETSc vector
  !!          object is allocated and calls the necessary destructor on the wrapped
  !!          PETSc linear_solver object, sets the allocated flag to .false. to prevent double
  !!          free's.
  !> @param[in/out] linear_solver solver - the linear solver to be destroyed.
  module subroutine free_linear_solver_petsc(solver)
    
    use petscksp, only: KSPDestroy

    type(linear_solver_petsc), intent(inout) :: solver

    integer(accs_err) :: ierr !> Error code

    if (solver%allocated) then
       call KSPDestroy(solver%KSP, ierr)
       solver%allocated = .false.
    else
       print *, "WARNING: attempted double free of linear solver"
    end if

  end subroutine 

end module petsctypes
