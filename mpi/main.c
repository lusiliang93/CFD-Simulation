#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "mpi.h"
#include "compute.h"

int main(int argc,char* argv[]){
	int opt=0;
	char *inputname=NULL;
	double startTime,endTime;

    //Initialize MPI
    MPI_Init(&argc,&argv);

	/*Read command line arguments*/
	do {
		opt = getopt(argc,argv,"f:"); /*why type -f?*/
		switch(opt){
			case 'f':
                inputname=optarg;
			break;

			default:
            break;
		}
	}while (opt!=-1);

	/* read inputs*/
	if (inputname==NULL){
		printf("Inputname is incorrect.\n");
		MPI_Finalize();
		return -1;
	}

	// Get process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);

    // Get total number of processes specificed at start of run
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    //Run computation
    startTime = MPI_Wtime();
    compute(procID,nproc,inputname);
    endTime = MPI_Wtime();

    printf("elapsed time for proc %d: %f\n", procID, endTime - startTime);
    return 0;
}
