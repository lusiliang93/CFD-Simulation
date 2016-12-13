#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#include "cuda_util.h"

#define THREADSPB 256

void cuda_init(int imax, int jmax){
}

#define get_index(i,j,jmax) ((jmax+2)*i+j)

__global__ void setbound_kernel_x(double* cudaDevice_u, double* cudaDevice_v, double* cudaDevice_u2, double* cudaDevice_v2, int imax, int jmax){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx>=1&&idx<jmax+1){
        cudaDevice_u2[get_index(idx, 0, jmax)] = 0;
        cudaDevice_u2[get_index(idx, imax, jmax)] = 0;
        cudaDevice_v2[get_index(idx, 0, jmax)] = -cudaDevice_v[get_index(idx, 1, jmax)];
        cudaDevice_v2[get_index(idx, imax+1, jmax)] = -cudaDevice_v[get_index(idx, imax, jmax)];
    }
}

__global__ void setbound_kernel_y(double* cudaDevice_u, double* cudaDevice_v, double* cudaDevice_u2, double* cudaDevice_v2, int imax, int jmax){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int us = 1;
    if(idx>=1&&idx<imax+1){
        cudaDevice_v2[get_index(0, idx, jmax)] = 0;
        cudaDevice_v2[get_index(jmax, idx, jmax)] = 0;
        cudaDevice_u2[get_index(0, idx, jmax)] = -cudaDevice_u[get_index(1, idx, jmax)];
        cudaDevice_u2[get_index(jmax+1, idx, jmax)] = 2*us - cudaDevice_u[get_index(jmax, idx, jmax)];
    }
}

void setbound(double **u,double **v,int imax,int jmax,int wW, int wE,int wN,int wS){
    double *cudaDevice_u, *cudaDevice_v, *cudaDevice_p, *cudaDevice_f, *cudaDevice_g, *cudaDevice_rhs;
    double *cudaDevice_u2, *cudaDevice_v2, *cudaDevice_p2, *cudaDevice_f2, *cudaDevice_g2, *cudaDevice_rhs2;
    cudaMalloc(&cudaDevice_u, (imax+2)*(jmax+2) *sizeof(double));
    cudaMalloc(&cudaDevice_v, (imax+2)*(jmax+2)*sizeof(double));
    cudaMalloc(&cudaDevice_p, (imax+2)*(jmax+2)*sizeof(double));
    cudaMalloc(&cudaDevice_f, (imax+2)*(jmax+2)*sizeof(double));
    cudaMalloc(&cudaDevice_g, (imax+2)*(jmax+2)*sizeof(double));
    cudaMalloc(&cudaDevice_rhs, (imax+2)*(jmax+2)*sizeof(double));
    cudaMalloc(&cudaDevice_u2, (imax+2)*(jmax+2)*sizeof(double));
    cudaMalloc(&cudaDevice_v2, (imax+2)*(jmax+2)*sizeof(double));
    cudaMalloc(&cudaDevice_p2, (imax+2)*(jmax+2)*sizeof(double));
    cudaMalloc(&cudaDevice_f2, (imax+2)*(jmax+2)*sizeof(double));
    cudaMalloc(&cudaDevice_g2, (imax+2)*(jmax+2)*sizeof(double));
    cudaMalloc(&cudaDevice_rhs2, (imax+2)*(jmax+2)*sizeof(double));

    cudaMemcpy(cudaDevice_u, u, sizeof(double)*(imax+2)*(jmax+2), cudaMemcpyHostToDevice);
    cudaMemcpy(cudaDevice_v, v, sizeof(double)*(imax+2)*(jmax+2), cudaMemcpyHostToDevice);
    printf("ijmax: %d %d\n", imax, jmax);
    int nBlocks = (jmax+1 + THREADSPB-1)/THREADSPB;
    setbound_kernel_x<<<nBlocks, THREADSPB>>>(cudaDevice_u, cudaDevice_v, cudaDevice_u2, cudaDevice_v2,imax,jmax);
    nBlocks = (imax+1 + THREADSPB-1)/THREADSPB;
    setbound_kernel_y<<<nBlocks, THREADSPB>>>(cudaDevice_u, cudaDevice_v, cudaDevice_u2, cudaDevice_v2,imax,jmax);
    cudaMemcpy(u, cudaDevice_u2, sizeof(double)*(imax+2)*(jmax+2), cudaMemcpyDeviceToHost);
    cudaMemcpy(v, cudaDevice_v2, sizeof(double)*(imax+2)*(jmax+2), cudaMemcpyDeviceToHost);
    return;
}