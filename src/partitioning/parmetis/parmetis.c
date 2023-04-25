#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <mpi.h>
// #include <parmetis.h>

#ifndef _MSC_VER
#define __cdecl
#endif

void partition_parmetiskway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy,
                            idxtype *vwgt, idxtype *adjwgt,
                            int *nparts, double *imbalance, int *suppress_output, 
                            int *seed, int *mode, int *edgecut, idxtype *part, 
                            MPI_Fint *ccomm)
{
    // int mype, npes, i;
    // MPI_Comm comm;

    // comm=MPI_Comm_f2c(*ccomm);
    // MPI_Comm_size(comm, &npes);
    // MPI_Comm_rank(comm, &mype);

    // bool suppress = (*suppress_output == 1) ? true : false;

/* 
Parmetis Interface definition:

int __cdecl ParMETIS_V3_PartKway(
             idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
	     idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts, 
	     real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, 
	     MPI_Comm *comm);
*/

    // ParMETIS_V3_PartKWay(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ncon,
    //                     nparts, tpwgts, ubvec, options, edgecut, part, &comm);

    printf("Nothing to see here - yet.");
}