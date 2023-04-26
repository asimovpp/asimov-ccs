#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <mpi.h>
#include <parmetis.h>

void partition_parmetiskway(idxtype *vtxdist,idxtype *xadj,idxtype *adjncy,
                              idxtype *vwgt,idxtype *adjwgt,
                              int *wgtflag, int *numflag, int *ncon, int *nparts,
                              float *tpwgts, float *ubvec, int *options, 
                              int *edgecut, idxtype *part, MPI_Fint *ccomm)

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