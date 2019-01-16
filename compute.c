#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>
#include "compute.h"

int main(int argc, char const *argv[])
{
    int ret;
    
    const size_t N1 = 20, N2 = 20;
    const double *const restrict x1 = calloc(N1, sizeof(double));
    const double *const restrict x2 = calloc(N2, sizeof(double));

    const double h1 = 1/N1;
    const double h2 = 1/N2;
    const double eps = 1e-5;

    int np, mp;

    MPI_Comm grid_comm;
    
    int dims[2];
    bool periodic[2];
    int coordinates[2];
    MPI_Comm row_comm;
    dims[0] = N1; dims[1] = N2;
    periodic[0] = periodic[1] = false;
    coords[0] = 0; coords[1] = 1;
    
    ret = MPI_Init(&argc, &argv);
    if (!ret)
        goto on_error;
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    MPI_Comm_rank(MPI_COMM_WORLD,&mp);
    
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, true, &grid_comm);
    MPI_Cart_coords(grid_comm, mp, 2, coordinates);


    MPI_Finalize();
    free(x1);
    free(x2);
    return 0;
on_error:
    free(x1);
    free(x2);
    return -1;
}
