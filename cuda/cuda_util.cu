#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "init.h"
#include "uvp.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>

__global__ void setbound_kernel_x(double** u, double** v, int imax, int jmax){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx>=1&&idx<jmax+1){
        u[idx][0] = 0;
        u[idx][imax] = 0;
        v[idx][0] = -v[idx][1];
        v[idx][imax+1] = -v[idx][imax];
    }
}

__global__ void setbound_kernel_y(double** u, double** v, int imax, int jmax){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int us = 1;
    if(idx>=1&&idx<imax+1){
        v[0][idx] = 0;
        v[jmax][idx] = 0;
        u[0][idx] = -u[1][idx];
        u[jmax+1][idx] = 2*us - u[jmax][idx];
    }
}

void setbound(double **u,double **v,int imax,int jmax,int wW, int wE,int wN,int wS){
    int nBlocks  = (jmax+1 + THREADSPB-1)/THREADSPB;
    setbound_kernel_x<<<nBlocks, THREADSPB>>>(u,v,imax,jmax);
    nBlocks  = (imax+1 + THREADSPB-1)/THREADSPB;
    setbound_kernel_y<<<nBlocks, THREADSPB>>>(u,v,imax,jmax);
    return;
}