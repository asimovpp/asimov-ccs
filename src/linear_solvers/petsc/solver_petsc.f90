!v Submodule file solver_petsc.smod
!
!  An implementation of a PETSc solver
!
!  @build petsc
submodule(solver) solver_petsc
#include "ccs_macros.inc"
#include <petscversion.h>

  use kinds, only: ccs_err
  use petsctypes, only: linear_solver_petsc, matrix_petsc, vector_petsc
  use parallel_types_mpi, only: parallel_environment_mpi
  use utils, only: update, exit_print
  use logging, only: log_unit_out
  use constants, only: ccs_string_len

  implicit none

contains

  !> Create a new PETSc solver object.
  module subroutine create_solver(linear_system, solver)

    use petsc, only: PETSC_TRUE
    use petscksp, only: KSPCreate, KSPSetOperators, KSPSetFromOptions, KSPSetInitialGuessNonzero

    type(equation_system), intent(in) :: linear_system        !< Data structure containing equation system to be solved.
    class(linear_solver), allocatable, intent(inout) :: solver  !< The linear solver returned allocated.

    integer(ccs_err) :: ierr ! Error code
    logical :: first_creation

    first_creation = .false.
    if (.not. allocated(solver)) then
      allocate (linear_solver_petsc :: solver)
      first_creation = .true.
    end if

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

            if (first_creation) then
              call KSPCreate(comm, ksp, ierr)
              if (ierr /= 0) then
                call error_abort("Error in creating solver KSP")
              end if
              call KSPSetOperators(ksp, M%M, M%M, ierr)
              call KSPSetFromOptions(ksp, ierr)
              call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
            end if

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
              write(log_unit_out,*) "ERROR in linear solve."
              call error_abort("ERROR in linear solve.")
            end if

          class default
            write(log_unit_out,*)"ERROR: Trying to use non-PETSc vector for solution with PETSc solver."
            call error_abort("ERROR: Trying to use non-PETSc vector for solution with PETSc solver.")
          end select

        class default
          write(log_unit_out,*) "ERROR: Trying to use non-PETSc vector for RHS with PETSc solver."
          call error_abort("ERROR: Trying to use non-PETSc vector for RHS with PETSc solver.")
        end select

      end associate

    class default
      call error_abort("Unknown solver type")

    end select

  end subroutine

  !> Interface to set the primary method of a linear solver
  module subroutine set_solver_method(method_name, solver)

#if PETSC_VERSION_GE(3,23,0)
    use petscksp, only: KSPSetType, KSPGetType, KSPSetFromOptions, KSPSetOptionsPrefix
#else
    use petscksp, only: KSPSetType, KSPSetFromOptions
#endif
    ! Arguments
    character(len=*), intent(in) :: method_name   !< String naming the linear solver to be used.
    class(linear_solver), intent(inout) :: solver !< The linear solver object

    ! Local
    character(len=ccs_string_len) :: petsc_method_name
    integer(ccs_err) :: ierr ! Error code

    select type (solver)
    type is (linear_solver_petsc)
      associate (ksp => solver%KSP)
        ! Set linear solver type directly from method name
        call KSPSetType(ksp, trim(method_name), ierr)

        if (ierr /= 0) then
          call error_abort("ERROR: setting solver method failed, " // trim(method_name) // " solver likely unsuported.")
        end if

#if PETSC_VERSION_LT(3,23,0)
        call KSPGetType(ksp, petsc_method_name, ierr)
        if (trim(petsc_method_name) /= trim(method_name)) then
          call error_abort("ERROR: petsc solver method ("//trim(petsc_method_name)//") doesn't match requested method ("// trim(method_name)//")")
        end if
#endif

        if (allocated(solver%linear_system%name)) then
          call KSPSetOptionsPrefix(ksp, solver%linear_system%name // ':', ierr)
        end if
        call KSPSetFromOptions(ksp, ierr)

      end associate
    class default
      call error_abort("ERROR: Unknown solver type")
    end select

  end subroutine set_solver_method

  !> Interface to set the preconditioner of a linear solver
  module subroutine set_solver_precon(precon_name, solver)

#if PETSC_VERSION_GE(3,23,0)
    use petscksp, only: KSPGetPC, tPC, PCSetType, PCGetType, PCSetReusePreconditioner, PCSetFromOptions, &
                        PCSetOptionsPrefix, KSPSetOptionsPrefix
#else
    use petscpc, only: tPC, PCSetType, PCSetReusePreconditioner, PCSetFromOptions
#endif

    use petsc, only: PETSC_TRUE

    ! Arguments
    character(len=*), intent(in) :: precon_name   !< String naming the preconditioner to be used.
    class(linear_solver), intent(inout) :: solver !< The linear solver object

    ! Local
    type(tPC) :: pc          ! PETSc preconditioner object
    character(len=ccs_string_len) :: petsc_precon_name
    integer(ccs_err) :: ierr ! Error code

    select type (solver)
    type is (linear_solver_petsc)
      associate (ksp => solver%KSP)
        call KSPGetPC(ksp, pc, ierr)

        ! Set preconditioner type directly using precon_name
        call PCSetType(pc, trim(precon_name), ierr)

        if (ierr /= 0) then
          call error_abort("ERROR: setting preconditioner method failed, " // trim(precon_name) // " preconditioner likely unsuported.")
        end if

#if PETSC_VERSION_LT(3,23,0)
        call PCGetType(pc, petsc_precon_name, ierr)
        if (trim(petsc_precon_name) /= trim(precon_name)) then
          call error_abort("ERROR: petsc precon method ("//trim(petsc_precon_name)//") doesn't match requested precon ("// trim(precon_name)//")")
        end if
#endif

        call PCSetReusePreconditioner(pc, PETSC_TRUE, ierr)

        ! Allow command-line options to override settings in source or config file
        if (allocated(solver%linear_system%name)) then
          call PCSetOptionsPrefix(pc, solver%linear_system%name // ':', ierr)
        end if
        call PCSetFromOptions(pc, ierr)

      end associate
    class default
      call error_abort("ERROR: Unknown solver type")
    end select

  end subroutine set_solver_precon

end submodule solver_petsc
