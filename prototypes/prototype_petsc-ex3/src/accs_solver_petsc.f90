submodule (accs_solver) accs_solver_petsc

  use accs_kinds, only : accs_int, accs_err
  use accs_petsctypes, only : linear_solver_petsc, matrix_petsc, vector_petsc
  
  implicit none

contains

  module subroutine create_solver(eqsys, solver)

    use petsc, only : PETSC_TRUE
    use petscksp, only : KSPCreate, KSPSetOperators, KSPSetFromOptions, KSPSetInitialGuessNonzero

    type(linear_system), intent(in) :: eqsys
    class(linear_solver), allocatable, intent(out) :: solver

    integer(accs_err) :: ierr
    
    allocate(linear_solver_petsc :: solver)
    
    select type(solver)
    type is(linear_solver_petsc)
       solver%eqsys = eqsys
       
       associate(comm => solver%eqsys%comm, &
            ksp => solver%KSP, &
            M => solver%eqsys%M)
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
       solver%allocated = .true.
    class default
       print *, "Unknown solver type"
       stop
    end select
    
  end subroutine

  module subroutine solve(solver)

    use petscksp, only : KSPSolve
    
    class(linear_solver), intent(inout) :: solver

    integer(accs_err) :: ierr
    
    select type(solver)
    type is(linear_solver_petsc)

       associate(ksp => solver%KSP, &
            b => solver%eqsys%rhs, &
            u => solver%eqsys%sol)
         select type (b)
         type is(vector_petsc)
            select type(u)
            type is(vector_petsc)
               call KSPSolve(ksp, b%v, u%v, ierr)
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
  
end submodule accs_solver_petsc
