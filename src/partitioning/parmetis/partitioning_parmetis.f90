submodule(partitioning) partitioning_parmetis
#include "ccs_macros.inc"

  use kinds, only: ccs_int, ccs_real, ccs_long
  use utils, only: str, debug_print
  use parallel_types_mpi, only: parallel_environment_mpi
  use meshing, only: set_local_num_cells, get_local_num_cells

  implicit none

  interface
    subroutine partition_parmetiskway(vtxdist, xadj, adjncy, vwgt, adjwgt, &
                                    wgtflag, numflag, ncon, num_procs, &
                                    tpwgts, ubvec, options, &
                                    edgecuts, local_partition, comm) bind(c)
      use iso_c_binding

      integer(c_long), dimension(*) :: vtxdist
      integer(c_long), dimension(*) :: xadj
      integer(c_long), dimension(*) :: adjncy
      integer(c_long), dimension(*) :: vwgt
      integer(c_long), dimension(*) :: adjwgt
      integer(c_long) :: wgtflag ! Set to 0 for "no weights"
      integer(c_long) :: numflag ! Numbering scheme - 1 means Fortran style
      integer(c_long) :: ncon
      integer(c_long) :: num_procs
      real(c_float), dimension(*) :: tpwgts
      real(c_float), dimension(*) :: ubvec
      integer(c_long), dimension(*) :: options
      integer(c_long) :: edgecuts
      integer(c_long), dimension(*) :: local_partition
      integer(c_int) :: comm
    end subroutine
  end interface

contains

  !v Partition the mesh
  !
  ! Use Parmetis library to compute a k-way vertex separator given a k-way partition of the graph.
  ! The graph can be weighted or unweighted.
  module subroutine partition_kway(par_env, mesh)

    use mpi
    use iso_c_binding

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    ! Local variables
    integer(ccs_long), dimension(:), allocatable :: tmp_partition
    integer(ccs_int) :: local_part_size
    integer(ccs_int) :: irank
    integer(ccs_int) :: ierr
    integer(ccs_int) :: i

    integer(c_long), dimension(:), allocatable :: vtxdist
    integer(c_long), dimension(:), allocatable :: xadj
    integer(c_long), dimension(:), allocatable :: adjncy
    integer(c_long), dimension(:), allocatable :: vwgt
    integer(c_long), dimension(:), allocatable :: adjwgt
    integer(c_long) :: wgtflag ! Set to 0 for "no weights"
    integer(c_long) :: numflag ! Numbering scheme - 1 means Fortran style
    integer(c_long) :: ncon
    integer(c_long) :: num_procs
    real(c_float), dimension(:), allocatable :: tpwgts
    real(c_float), dimension(:), allocatable :: ubvec
    integer(c_long), dimension(:), allocatable :: options
    integer(c_long) :: edgecuts
    integer(c_long), dimension(:), allocatable :: local_partition
    integer(c_int) :: comm

    ! Values mostly hardcoded for now
    wgtflag = 0 ! No weights
    numflag = 0 ! Use C-style indexing for now
    ncon = 3
    num_procs = par_env%num_procs
    edgecuts = -1 ! XXX: silence unused variable warning

    allocate(ubvec(ncon))
    allocate(tpwgts(ncon * num_procs))
    allocate(options(0:2))

    options(0) = 0 ! 0 = default values, 1 = values specified in (1) and (2)
    options(1) = 1 ! Output verbosity - 1 gives timing information
    options(2) = 2023 ! Random number seed

    ubvec(:) = 1.05 ! Imbalance tolerance for each vertex weight, 1.05 is recommended value
    tpwgts(:) = 1.0 / real(num_procs, c_float) ! XXX: Not quite correct, though probably does not matter as
                              ! we are not using weights
                              ! Fraction of vertex weight that should be distributed 
                              ! to each sub-domain. Sum of tpwgts(:) should be 1.

    allocate (tmp_partition(mesh%topo%global_num_cells)) ! Temporary partition array
    tmp_partition = 0

    irank = par_env%proc_id ! Current rank

    vtxdist = mesh%topo%vtxdist - 1
    xadj = mesh%topo%xadj - 1
    adjncy = mesh%topo%adjncy - 1

    adjwgt = mesh%topo%adjwgt
    vwgt = mesh%topo%vwgt

    ! Set weights to 1
    adjwgt = 1
    vwgt = 1

    ! Number of elements in local partition array
    ! Needed for gathering loca partitions into global partition array
    local_part_size = size(mesh%topo%local_partition)

    allocate (local_partition(local_part_size))

    ! Partitioning an unweighted graph
    select type (par_env)
    type is (parallel_environment_mpi)

      comm = par_env%comm

      call partition_parmetiskway(vtxdist, xadj, adjncy, vwgt, adjwgt, &
                                wgtflag, numflag, ncon, num_procs, &
                                tpwgts, ubvec, options, &
                                edgecuts, local_partition, comm)

      mesh%topo%local_partition(:) = local_partition(:)

      do i = 1, local_part_size
        tmp_partition(i + vtxdist(irank + 1)) = mesh%topo%local_partition(i)
      end do

      call MPI_AllReduce(tmp_partition, mesh%topo%global_partition, mesh%topo%global_num_cells, &
                         MPI_LONG, MPI_SUM, par_env%comm, ierr)

    class default
      print *, "ERROR: Unknown parallel environment!"
    end select

    call dprint("Number of edgecuts: " // str(int(edgecuts)))

    deallocate (tmp_partition)

  end subroutine partition_kway

  !v Compute the input arrays for the partitioner
  !
  ! Using the topology object, compute the input arrays for the Parmetis partitioner
  ! Input arrays for the partitioner are: vtxdist, xadj and adjncy
  module subroutine compute_partitioner_input(par_env, mesh)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    call compute_partitioner_input_generic(par_env, mesh)

  end subroutine compute_partitioner_input

end submodule
