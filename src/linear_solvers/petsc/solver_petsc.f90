!>  Submodule file solver_petsc.smod
!
!>  @build petsc
!
!>  An implementation of a PETSc solver
submodule (solver) solver_petsc

  use kinds, only : ccs_err
  use petsctypes, only : linear_solver_petsc, matrix_petsc, vector_petsc
  use parallel_types_mpi, only: parallel_environment_mpi
  use utils, only: update
  
  implicit none

contains

  !>  Create a new PETSc solver object.
  module subroutine create_solver(linear_system, solver)

    use petsc, only : PETSC_TRUE
    use petscksp, only : KSPCreate, KSPSetOperators, KSPSetFromOptions, KSPSetInitialGuessNonzero

    type(equation_system), intent(in) :: linear_system        !< Data structure containing equation system to be solved.
    class(linear_solver), allocatable, intent(out) :: solver  !< The linear solver returned allocated.

    integer(ccs_err) :: ierr !< Error code
    
    allocate(linear_solver_petsc :: solver)
    
    select type(solver)
      type is(linear_solver_petsc)

        select type (par_env => linear_system%par_env)
          type is(parallel_environment_mpi)
            
            solver%linear_system = linear_system
       
            associate(comm => par_env%comm, &
                      ksp  => solver%KSP, &
                      M    => solver%linear_system%matrix)

              select type(M)
                type is(matrix_petsc)

                  call KSPCreate(comm, ksp, ierr)
                  if (ierr /= 0) then
                    print *, "Error in creating solver KSP"
                    stop
                  end if
                  call KSPSetOperators(ksp, M%M, M%M, ierr)
                  call KSPSetFromOptions(ksp, ierr)
                  call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)

                class default
                  print *, "ERROR: Trying to use non-PETSc matrix with PETSc solver!"
                  stop

              end select

            end associate

          class default
            print *, "Unknown parallel environment"
    
        end select         

        solver%allocated = .true.

      class default
        print *, "Unknown solver type"
        stop

    end select
    
  end subroutine

  !>  Solve the linear system in a PETSc solver.
  module subroutine solve(solver)

    use petscksp, only : KSPSolve
    
    class(linear_solver), intent(inout) :: solver   !< The linear solver object.

    integer(ccs_err) :: ierr !< Error code
    
    select type(solver)
      type is(linear_solver_petsc)

        associate(ksp => solver%KSP)

          select type (b => solver%linear_system%rhs)
            type is(vector_petsc)

              select type(u => solver%linear_system%solution)
                type is(vector_petsc)
                  call KSPSolve(ksp, b%v, u%v, ierr)
                  call update(u)
                  if (ierr /= 0) then
                    print *, "ERROR in linear solve!"
                  end if

              class default
                print *, "ERROR: Trying to use non-PETSc vector for solution with PETSc solver!"
                stop
              end select

            class default
              print *, "ERROR: Trying to use non-PETSc vector for RHS with PETSc solver!"
              stop
          end select

        end associate

      class default
        print *, "Unknown solver type"
        stop

    end select
       
  end subroutine
  
end submodule solver_petsc
