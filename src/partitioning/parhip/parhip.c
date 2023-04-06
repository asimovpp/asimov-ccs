#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <mpi.h>
#include <parhip_interface.h>

void partition_parhipkway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy,
                            idxtype *vwgt, idxtype *adjwgt,
                            int *nparts, double *imbalance, int *suppress_output, 
                            int *seed, int *mode, int *edgecut, idxtype *part, 
                            MPI_Fint *ccomm)
{
    int mype, npes, i;
    MPI_Comm comm;

    comm=MPI_Comm_f2c(*ccomm);
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &mype);

    bool suppress = (*suppress_output == 1) ? true : false;

/* 
ParHIP Interface definition:

void ParHIPPartitionKWay(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt,
                         int *nparts, double* imbalance, bool suppress_output, int seed, int mode, int *edgecut, idxtype *part, 
                         MPI_Comm *comm);

*/

    ParHIPPartitionKWay(vtxdist, xadj, adjncy, vwgt, adjwgt,
                        nparts, imbalance, suppress, *seed, *mode, edgecut, part, &comm);

}