#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#include "cuda_util.h"
#include <thrust/scan.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <thrust/extrema.h>
#include <thrust/reduce.h>
#include "cublas_v2.h"

#define THREADSPB 256
#define get_index(i,j) ((jmax+2)*i+j)

double *cudaDevice_u, *cudaDevice_v, *cudaDevice_p, *cudaDevice_f, *cudaDevice_g, *cudaDevice_rhs;
double *cudaDevice_u2, *cudaDevice_v2, *cudaDevice_p2, *cudaDevice_f2, *cudaDevice_g2, *cudaDevice_rhs2;

__global__ void fill_val(double* p, int length, int val){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    p[idx] = val;
}

/* sum all the value in vector device_p */
double sum_vector(double* device_p, int length){
    thrust::device_ptr<double> d_ptr = thrust::device_pointer_cast(device_p);
    double sum = thrust::reduce(d_ptr, d_ptr+length, (double)0.0, thrust::plus<double>());
    return sum;
}

/* return the max value in vector device_p */
double max_vector(double* device_p, int length){
    cublasHandle_t handle;
    cublasCreate(&handle);
    double mymax = 0.0;
    int max_idx = 0;
    cublasIdamax(handle, length, device_p, 1, &max_idx);
    mymax = device_p[max_idx];
    return mymax;
}

void cuda_init(int imax, int jmax){
    cudaMalloc(&cudaDevice_u, (imax+2)*(jmax+2)*sizeof(double));
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

    int nBlocks = ((imax+2)*(jmax+2) + THREADSPB-1)/THREADSPB;
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_u, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_v, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_p, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_f, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_g, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_rhs, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_u2, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_v2, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_p2, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_f2, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_g2, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_rhs2, (imax+2)*(jmax+2), 0);
}

/* copy from matrix to matrix 2(matrix 2 is the stale data from last iteration and is  read-only) */
void copy_matrix(int imax, int jmax){
    double* tmp_u = cudaDevice_u2;
    double* tmp_v = cudaDevice_v2;
    double* tmp_p = cudaDevice_p2;
    double* tmp_f = cudaDevice_f2;
    double* tmp_g = cudaDevice_g2;
    double* tmp_rhs = cudaDevice_rhs2;
    cudaDevice_u2 = cudaDevice_u;
    cudaDevice_v2 = cudaDevice_v;
    cudaDevice_p2 = cudaDevice_p;
    cudaDevice_f2 = cudaDevice_f;
    cudaDevice_g2 = cudaDevice_g;
    cudaDevice_rhs2 = cudaDevice_rhs;
    cudaDevice_u = tmp_u;
    cudaDevice_v = tmp_v;
    cudaDevice_p = tmp_p;
    cudaDevice_f = tmp_f;
    cudaDevice_g = tmp_g;
    cudaDevice_rhs = tmp_rhs;
    int nBlocks = ((imax+2)*(jmax+2) + THREADSPB-1)/THREADSPB;
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_u, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_v, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_p, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_f, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_g, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_rhs, (imax+2)*(jmax+2), 0);
}

double comp_delt(int imax, int jmax,double delx,double dely,double Re,double tau){
    double first,second,third,min;
    double delta = 1/(delx*delx)+1/(dely*dely);
    first = Re/2/delta;
    min=first;
    int length = (imax+2)*(jmax+2);
    second = delx/abs(max_vector(cudaDevice_u2, length));
    third= dely/abs(max_vector(cudaDevice_v2, length));
    if(min>second){
        min=second;
        if(min>third)
            min=third;
    }
    else{
        if(min>third)
            min=third;
    }
    double ret = tau*min;
    return ret;
}

__global__ void setbound_kernel_x(double* cudaDevice_u, double* cudaDevice_v, double* cudaDevice_u2, double* cudaDevice_v2, int imax, int jmax){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int j = idx;
    if(idx>=1&&idx<jmax+1){
        cudaDevice_u[get_index(j, 0)] = 0;
        cudaDevice_u[get_index(j, imax)] = 0;
        cudaDevice_v[get_index(j, 0)] = -cudaDevice_v2[get_index(j, 1)];
        cudaDevice_v[get_index(j, imax+1)] = -cudaDevice_v2[get_index(j, imax)];
    }
}

__global__ void setbound_kernel_y(double* cudaDevice_u, double* cudaDevice_v, double* cudaDevice_u2, double* cudaDevice_v2, int imax, int jmax){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i = idx;
    int us = 1;
    if(idx>=1&&idx<imax+1){
        cudaDevice_v[get_index(0, i)] = 0;
        cudaDevice_v[get_index(jmax, i)] = 0;
        cudaDevice_u[get_index(0, i)] = -cudaDevice_u2[get_index(1, i)];
        cudaDevice_u[get_index(jmax+1, i)] = 2*us - cudaDevice_u2[get_index(jmax, i)];
    }
}

void setbound(int imax,int jmax,int wW, int wE,int wN,int wS){
    int nBlocks = (jmax+1 + THREADSPB-1)/THREADSPB;
    setbound_kernel_x<<<nBlocks, THREADSPB>>>(cudaDevice_u, cudaDevice_v, cudaDevice_u2, cudaDevice_v2,imax,jmax);
    nBlocks = (imax+1 + THREADSPB-1)/THREADSPB;
    setbound_kernel_y<<<nBlocks, THREADSPB>>>(cudaDevice_u, cudaDevice_v, cudaDevice_u2, cudaDevice_v2,imax,jmax);
    return;
}

__global__ void init_uvp_kernel(double* cudaDevice_u, double* cudaDevice_v, double* cudaDevice_p, int imax, int jmax, int UI, int VI, int PI){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int j = idx/(jmax+2);
    int i = idx%(jmax+2);
    if(idx<(imax+2)*(jmax+2)){
        cudaDevice_u[get_index(j,i)] = UI;
        cudaDevice_v[get_index(j,i)] = VI;
        cudaDevice_p[get_index(j,i)] = PI;
    }
}

void init_uvp(int imax, int jmax,int UI, int VI, int PI){
    int nBlocks = ((jmax+2)*(imax+2) + THREADSPB-1)/THREADSPB;
    init_uvp_kernel<<<nBlocks, THREADSPB>>>(cudaDevice_u, cudaDevice_v, cudaDevice_p, imax,jmax,UI,VI,PI);
}

__global__ void comp_fg_kernel_1(double* cudaDevice_u2, double* cudaDevice_v2, double* cudaDevice_f, double* cudaDevice_g, int imax, int jmax){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i = idx;
    int j = idx;
    if(j>=1&&j<jmax+1){
        cudaDevice_f[get_index(j,0)] = cudaDevice_u2[get_index(j,0)];
        cudaDevice_f[get_index(j,imax)] = cudaDevice_u2[get_index(j,imax)];
    }
    if(i>=1&&i<imax+1){
        cudaDevice_g[get_index(0,i)] = cudaDevice_v2[get_index(0,i)];
        cudaDevice_g[get_index(jmax,i)] = cudaDevice_v2[get_index(jmax,i)];
    }
}

__global__ void comp_fg_kernel_2(double* cudaDevice_u2, double* cudaDevice_v2, double* cudaDevice_f, double* cudaDevice_g, int imax, int jmax, double delt,double delx,double dely,double gx,double gy,double gamma,double Re){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int j = idx/(jmax+2);
    int i = idx%(jmax+2);
    double a,b,c,d,e,ff,gg,h,va,vb,u2x,uvy,u2x2,u2y2;
    double ua,ub,uvx,v2y,v2x2,v2y2;
    if(j>=1&&j<jmax+1){
        if(i>=1&&i<imax){
            a = cudaDevice_u2[get_index(j,i)] + cudaDevice_u2[get_index(j,i+1)];
            b = cudaDevice_u2[get_index(j,i-1)] + cudaDevice_u2[get_index(j,i)];
            c = cudaDevice_u2[get_index(j,i)] - cudaDevice_u2[get_index(j,i+1)];
            d = cudaDevice_u2[get_index(j,i-1)] - cudaDevice_u2[get_index(j,i)];
            e = cudaDevice_u2[get_index(j,i)] + cudaDevice_u2[get_index(j+1,i)];
            ff = cudaDevice_u2[get_index(j-1,i)] + cudaDevice_u2[get_index(j,i)];
            gg = cudaDevice_u2[get_index(j,i)] - cudaDevice_u2[get_index(j+1,i)];
            h = cudaDevice_u2[get_index(j-1,i)] - cudaDevice_u2[get_index(j,i)];
            va = cudaDevice_v2[get_index(j,i)] + cudaDevice_v2[get_index(j,i+1)];
            vb = cudaDevice_v2[get_index(j-1,i)] + cudaDevice_v2[get_index(j-1,i+1)];
            u2x = 1/delx*((a/2)*(a/2)-(b/2)*(b/2))+gamma*1/delx*(abs(a)/2*c/2-abs(b)/2*d/2);
            uvy = 1/dely*(va/2*e/2-vb/2*ff/2)+gamma*1/dely*(abs(va)/2*gg/2-abs(vb)/2*h/2);
            u2x2 = (cudaDevice_u2[get_index(j,i+1)] - 2*cudaDevice_u2[get_index(j,i)] + cudaDevice_u2[get_index(j,i-1)])/(delx*delx);
            u2y2 = (cudaDevice_u2[get_index(j+1,i)] - 2*cudaDevice_u2[get_index(j,i)] + cudaDevice_u2[get_index(j-1,i)])/(dely*dely);
            cudaDevice_f[get_index(j,i)] = cudaDevice_u2[get_index(j,i)] + delt*(1/Re*(u2x2+u2y2)-u2x-uvy+gx);
        }
    }

    if(j>=1&&j<jmax){
        if(i>=1&&i<imax+1){
            a = cudaDevice_v2[get_index(j,i)] + cudaDevice_v2[get_index(j,i+1)];
            b = cudaDevice_v2[get_index(j,i-1)] + cudaDevice_v2[get_index(j,i)];
            c = cudaDevice_v2[get_index(j,i)] - cudaDevice_v2[get_index(j,i+1)];
            d = cudaDevice_v2[get_index(j,i-1)] - cudaDevice_v2[get_index(j,i)];
            e = cudaDevice_v2[get_index(j,i)] + cudaDevice_v2[get_index(j+1,i)];
            ff = cudaDevice_v2[get_index(j-1,i)] + cudaDevice_v2[get_index(j,i)];
            gg = cudaDevice_v2[get_index(j,i)] - cudaDevice_v2[get_index(j+1,i)];
            h = cudaDevice_v2[get_index(j-1,i)] - cudaDevice_v2[get_index(j,i)];
            ua = cudaDevice_u2[get_index(j,i)] + cudaDevice_u2[get_index(j+1,i)];
            ub = cudaDevice_u2[get_index(j,i-1)] + cudaDevice_u2[get_index(j+1,i-1)];
            uvx = 1/delx*(ua/2*a/2-ub/2*b/2)+gamma*1/delx*(abs(ua)/2*c/2-abs(ub)/2*d/2);
            v2y = 1/dely*((e/2)*(e/2)-(ff/2)*(ff/2))+gamma*1/dely*(abs(e)/2*gg/2-abs(ff)/2*h/2);
            v2x2 = (cudaDevice_v2[get_index(j,i+1)] - 2*cudaDevice_v2[get_index(j,i)] + cudaDevice_v2[get_index(j,i-1)])/(delx*delx);
            v2y2 = (cudaDevice_v2[get_index(j+1,i)] - 2*cudaDevice_v2[get_index(j,i)] + cudaDevice_v2[get_index(j-1,i)])/(dely*dely);
            cudaDevice_g[get_index(j,i)] = cudaDevice_v2[get_index(j,i)] + delt*(1/Re*(v2x2+v2y2)-uvx-v2y+gy);
        }
    }
}

void comp_fg(int imax, int jmax,double delt,double delx,double dely,double gx,double gy,double gamma,double Re){
    int nBlocks = ((imax+2)*(jmax+2) + THREADSPB-1)/THREADSPB;
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_f, (imax+2)*(jmax+2), 0);
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_g, (imax+2)*(jmax+2), 0);

    nBlocks = (max(imax,jmax)+2 + THREADSPB-1)/THREADSPB;
    comp_fg_kernel_1<<<nBlocks, THREADSPB>>>(cudaDevice_u2, cudaDevice_v2, cudaDevice_f, cudaDevice_g, imax, jmax);

    nBlocks = ((imax+2)*(jmax+2) + THREADSPB-1)/THREADSPB;
    comp_fg_kernel_2<<<nBlocks, THREADSPB>>>(cudaDevice_u2, cudaDevice_v2, cudaDevice_f, cudaDevice_g, imax, jmax, delt, delx, dely, gx, gy, gamma, Re);
}

__global__ void comp_rhs_kernel(double* cudaDevice_f2, double* cudaDevice_g2, double* cudaDevice_rhs, int imax, int jmax, double delx, double dely, double delt){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int j = idx/(jmax+2);
    int i = idx%(jmax+2);
    if(j>=1&&j<jmax+1){
        if(i>=1&&i<imax+1){
            int tmp = (cudaDevice_f2[get_index(j,i)]-cudaDevice_f2[get_index(j,i-1)])/delx + (cudaDevice_g2[get_index(j,i)]-cudaDevice_g2[get_index(j-1,i)])/dely;
            cudaDevice_rhs[get_index(j,i)] = 1/delt * tmp;
        }
    }
}

void comp_rhs(int imax, int jmax,double delt,double delx,double dely){
    int nBlocks = ((imax+2)*(jmax+2) + THREADSPB-1)/THREADSPB;
    fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_rhs, (imax+2)*(jmax+2), 0);

    nBlocks = ((imax+2)*(jmax+2) + THREADSPB-1)/THREADSPB;
    comp_rhs_kernel<<<nBlocks, THREADSPB>>>(cudaDevice_f2, cudaDevice_g2, cudaDevice_rhs, imax, jmax, delx, dely, delt);
}

__global__ void poisson_kernel_1(double* cudaDevice_p, double* cudaDevice_p2, double* cudaDevice_r, int imax, int jmax){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int j = idx;
    int i = idx;
    if(j>=1&&j<jmax+1){
        cudaDevice_p[get_index(j,0)] = cudaDevice_p2[get_index(j,1)];
        cudaDevice_p[get_index(j,imax+1)] = cudaDevice_p2[get_index(j,imax)];
    }
    if(i>=1&&i<imax+1){
        cudaDevice_p[get_index(0,i)] = cudaDevice_p2[get_index(1,i)];
        cudaDevice_p[get_index(jmax+1,i)] = cudaDevice_p2[get_index(jmax,i)];
    }
}

__global__ void poisson_kernel_2(double* cudaDevice_r, double* cudaDevice_p, double* cudaDevice_p2, double* cudaDevice_rhs2, int imax, int jmax, double delx, double dely, double omg){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int j = idx/(jmax+2);
    int i = idx%(jmax+2);
    int eiw,eie,ejs,ejn;
    if(j>=1&&j<jmax+1){
        if(i>=1&&i<imax+1){
            eiw=1;eie=1;ejs=1;ejn=1;
            cudaDevice_p[get_index(j,i)] = (1-omg)*cudaDevice_p2[get_index(j,i)]
            +
            omg/((eie+eiw)/(delx*delx)+(ejn+ejs)/(dely*dely)) * (
                (eie*cudaDevice_p2[get_index(j,i+1)]+eiw*cudaDevice_p2[get_index(j,i-1)])/(delx*delx)
                +(ejn*cudaDevice_p2[get_index(j+1,i)]+ejs*cudaDevice_p2[get_index(j-1,i)])/(dely*dely)
                -cudaDevice_rhs2[get_index(j,i)]
            );

            cudaDevice_r[get_index(j,i)] = (
                eie*(cudaDevice_p2[get_index(j,i+1)]-cudaDevice_p2[get_index(j,i)])
                -eiw*(cudaDevice_p2[get_index(j,i)]-cudaDevice_p2[get_index(j,i-1)])
                )/(delx*delx)
            +    (
                ejn*(cudaDevice_p2[get_index(j+1,i)]-cudaDevice_p2[get_index(j,i)])
                -ejs*(cudaDevice_p2[get_index(j,i)]-cudaDevice_p2[get_index(j-1,i)])
                )/(dely*dely)
            - cudaDevice_rhs2[get_index(j,i)];

            cudaDevice_r[get_index(j,i)] = cudaDevice_r[get_index(j,i)]*cudaDevice_r[get_index(j,i)];
        }
    }
}

int poisson(int imax, int jmax,double delx,double dely,double eps,int itermax,double omg){
    int it;
    double sum;
    double res;
    double* cudaDevice_r;
    cudaMalloc(&cudaDevice_r, (imax+2)*(jmax+2) *sizeof(double));
    for(it=0;it<itermax;it++){
        int nBlocks = ((imax+2)*(jmax+2) + THREADSPB-1)/THREADSPB;
        /* Init of r to 0 can be moved out of the loop */
        fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_r, (imax+2)*(jmax+2), 0);

        nBlocks = (max(imax,jmax)+2 + THREADSPB-1)/THREADSPB;
        poisson_kernel_1<<<nBlocks, THREADSPB>>>(cudaDevice_p, cudaDevice_p2, cudaDevice_r, imax, jmax);

        nBlocks = ((imax+2)*(jmax+2) + THREADSPB-1)/THREADSPB;
        poisson_kernel_2<<<nBlocks, THREADSPB>>>(cudaDevice_r, cudaDevice_p, cudaDevice_p2, cudaDevice_rhs2, imax, jmax, delx, dely, omg);

        sum = sum_vector(cudaDevice_r, (imax+2)*(jmax+2));
        res=sqrt(sum/(imax*jmax));
        if(res<eps){
            break;
        }
        /* copy p to stale p(p2) */
        double* tmp_p = cudaDevice_p2;
        cudaDevice_p2 = cudaDevice_p;
        cudaDevice_p = tmp_p;
    }
    cudaFree(cudaDevice_r);
    return it;
}

__global__ void adap_uv_kernel(double* cudaDevice_u, double* cudaDevice_v, double* cudaDevice_f2, double* cudaDevice_g2, double* cudaDevice_p2, int imax,int jmax,double delt,double delx,double dely){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int j = idx/(jmax+2);
    int i = idx%(jmax+2);
    if(j>=1&&j<jmax+1){
        if(i>=1&&i<imax){
            cudaDevice_u[get_index(j,i)] = cudaDevice_f2[get_index(j,i)]-delt/delx*(cudaDevice_p2[get_index(j,i+1)]-cudaDevice_p2[get_index(j,i)]);
        }
    }
    if(j>=1&&j<jmax){
        if(i>=1&&i<imax+1){
            cudaDevice_v[get_index(j,i)] = cudaDevice_g2[get_index(j,i)]-delt/dely*(cudaDevice_p2[get_index(j+1,i)]-cudaDevice_p2[get_index(j,i)]);
        }
    }
}

void adap_uv(int imax, int jmax, double delt, double delx, double dely){
    int nBlocks = nBlocks = ((imax+2)*(jmax+2) + THREADSPB-1)/THREADSPB;
    adap_uv_kernel<<<nBlocks, THREADSPB>>>(cudaDevice_u, cudaDevice_v, cudaDevice_f2, cudaDevice_g2, cudaDevice_p2, imax, jmax, delt, delx, dely);
}
