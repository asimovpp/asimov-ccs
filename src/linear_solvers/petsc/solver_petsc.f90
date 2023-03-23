!v Submodule file solver_petsc.smod
!
!  An implementation of a PETSc solver
!
!  @build petsc
submodule(solver) solver_petsc
#include "ccs_macros.inc"

  use kinds, only: ccs_err
  use petsctypes, only: linear_solver_petsc, matrix_petsc, vector_petsc
  use parallel_types_mpi, only: parallel_environment_mpi
  use utils, only: update, exit_print

  implicit none

contains

  !> Create a new PETSc solver object.
  module subroutine create_solver(linear_system, solver)

    use petsc, only: PETSC_TRUE
    use petscksp, only: KSPCreate, KSPSetOperators, KSPSetFromOptions, KSPSetInitialGuessNonzero
    

    type(equation_system), intent(in) :: linear_system        !< Data structure containing equation system to be solved.
    class(linear_solver), allocatable, intent(inout) :: solver  !< The linear solver returned allocated.

    integer(ccs_err) :: ierr ! Error code

    if (allocated(solver)) then
      return
    end if

    allocate (linear_solver_petsc :: solver)

    select type (solver)
    type is (linear_solver_petsc)

      select type (par_env => linear_system%par_env)
      type is (parallel_environment_mpi)

        solver%linear_system = linear_system

        associate (comm => par_env%comm, &
                   ksp => solver%KSP, &
                   M => solver%linear_system%matrix)

          select type (M)
          type is (matrix_petsc)

            call KSPCreate(comm, ksp, ierr)
            if (ierr /= 0) then
              call error_abort("Error in creating solver KSP")
            end if
            call KSPSetOperators(ksp, M%M, M%M, ierr)
            call KSPSetFromOptions(ksp, ierr)
            call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)

          class default
            call error_abort("ERROR: Trying to use non-PETSc matrix with PETSc solver.")

          end select

        end associate

      class default
        call error_abort("Unknown parallel environment")

      end select

      solver%allocated = .true.

    class default
      call error_abort("Unknown solver type")

    end select

  end subroutine

  !> Solve the linear system in a PETSc solver.
  module subroutine solve(solver)

    use petscksp, only: KSPSolve

    class(linear_solver), intent(inout) :: solver   !< The linear solver object.

    integer(ccs_err) :: ierr ! Error code

    select type (solver)
    type is (linear_solver_petsc)

      associate (ksp => solver%KSP)

        select type (b => solver%linear_system%rhs)
        type is (vector_petsc)

          select type (u => solver%linear_system%solution)
          type is (vector_petsc)
            call KSPSolve(ksp, b%v, u%v, ierr)
            call update(u)
            if (ierr /= 0) then
              call error_abort("ERROR in linear solve.")
            end if

          class default
            call error_abort("ERROR: Trying to use non-PETSc vector for solution with PETSc solver.")
          end select

        class default
          call error_abort("ERROR: Trying to use non-PETSc vector for RHS with PETSc solver.")
        end select

      end associate

    class default
      call error_abort("Unknown solver type")

    end select

  end subroutine

  !> Interface to set the primary method of a linear solver
  module subroutine set_solver_method(method_name, solver)

    use petscksp, only: KSPSetType, KSPSetFromOptions

    ! Arguments
    character(len=*), intent(in) :: method_name   !< String naming the linear solver to be used.
    class(linear_solver), intent(inout) :: solver !< The linear solver object

    ! Local
    integer(ccs_err) :: ierr ! Error code

    select type (solver)
    type is (linear_solver_petsc)
      associate (ksp => solver%KSP)
        ! Set linear solver type directly from method name
        call KSPSetType(ksp, method_name, ierr)

        if (allocated(solver%linear_system%name)) then
          call KSPSetOptionsPrefix(ksp, solver%linear_system%name//':', ierr)
        endif
        call KSPSetFromOptions(ksp, ierr)
        
      end associate
    class default
      call error_abort("ERROR: Unknown solver type")
    end select

  end subroutine set_solver_method

  !> Interface to set the preconditioner of a linear solver
  module subroutine set_solver_precon(precon_name, solver)

    use petscksp, only: KSPGetPC
    use petscpc, only: tPC, PCSetType
    use petsc, only: PETSC_TRUE

    ! Arguments
    character(len=*), intent(in) :: precon_name   !< String naming the preconditioner to be used.
    class(linear_solver), intent(inout) :: solver !< The linear solver object

    ! Local
    type(tPC) :: pc          ! PETSc preconditioner object
    integer(ccs_err) :: ierr ! Error code

    select type (solver)
    type is (linear_solver_petsc)
      associate (ksp => solver%KSP)
        call KSPGetPC(ksp, pc, ierr)

        ! Set preconditioner type directly using precon_name
        call PCSetType(pc, precon_name, ierr)
        call PCSetReusePreconditioner(pc, PETSC_TRUE, ierr)
        
        ! Allow command-line options to override settings in source or config file
        if (allocated(solver%linear_system%name)) then
          call PCSetOptionsPrefix(pc, solver%linear_system%name//':', ierr)
        endif
        call PCSetFromOptions(pc, ierr)

      end associate
    class default
      call error_abort("ERROR: Unknown solver type")
    end select

  end subroutine set_solver_precon

end submodule solver_petsc
