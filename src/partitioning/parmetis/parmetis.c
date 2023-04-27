#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <mpi.h>
#include <parmetis.h>

void partition_parmetiskway(idx_t *vtxdist,idx_t *xadj,idx_t *adjncy,
                              idx_t *vwgt,idx_t *adjwgt,
                              idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts,
                              float *tpwgts, float *ubvec, idx_t *options, 
                              idx_t *edgecut, idx_t *part, MPI_Fint *ccomm)

{
  int mype, npes;
  MPI_Comm comm;

  comm=MPI_Comm_f2c(*ccomm);
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, 
                       wgtflag, numflag, ncon, nparts, 
	                     tpwgts, ubvec, options, 
                       edgecut, part, &comm);

}