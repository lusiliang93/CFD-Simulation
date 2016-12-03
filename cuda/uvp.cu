#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "init.h"
#include "uvp.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>

extern float toBW(int bytes, float sec);

void printCudaInfo() {
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    printf("---------------------------------------------------------\n");
    printf("Found %d CUDA devices\n", deviceCount);
    for (int i=0; i<deviceCount; i++) {
        cudaDeviceProp deviceProps;
        cudaGetDeviceProperties(&deviceProps, i);
        printf("Device %d: %s\n", i, deviceProps.name);
        printf("   SMs:        %d\n", deviceProps.multiProcessorCount);
        printf("   Global mem: %.0f MB\n",
               static_cast<float>(deviceProps.totalGlobalMem) / (1024 * 1024));
        printf("   CUDA Cap:   %d.%d\n", deviceProps.major, deviceProps.minor);
    }
    printf("---------------------------------------------------------\n");
}

template<typename T>
__global__ 
void kernelSum(const T* __restrict__ input, 
        const size_t lda, // pitch of input in words of sizeof(T)
        T* __restrict__ per_block_results, 
                const size_t n)
{
    extern __shared__ T sdata[];

    T x = 0.0;
    const T * p = &input[blockIdx.x * lda];
    // Accumulate per thread partial sum
    for(int i=threadIdx.x; i < n; i += blockDim.x) {
        x += p[i];
    }

    // load thread partial sum into shared memory
    sdata[threadIdx.x] = x;
    __syncthreads();

    for(int offset = blockDim.x / 2; offset > 0; offset >>= 1) {
        if(threadIdx.x < offset) {
            sdata[threadIdx.x] += sdata[threadIdx.x + offset];
        }
        __syncthreads();
    }

    // thread 0 writes the final result
    if(threadIdx.x == 0) {
        per_block_results[blockIdx.x] = sdata[0];
    }
}

void CDFCuda(int N, float alpha, float* xarray, float* yarray, float* resultarray) {

    int totalBytes = sizeof(float) * 3 * N;
    size_t size = N * sizeof(float);

    // compute number of blocks and threads per block
    const int threadsPerBlock = 512;
    const int blocks = (N + threadsPerBlock - 1) / threadsPerBlock;

    float* device_x;
    float* device_y;
    float* device_result;

    // allocate device memory buffers on the GPU using cudaMalloc
    cudaMalloc(&device_x, size);
    cudaMalloc(&device_y, size);
    cudaMalloc(&device_result, size);

    // start timing after allocation of device memory
    double startTime = CycleTimer::currentSeconds();

    // copy input arrays to the GPU using cudaMemcpy
    cudaMemcpy(device_x, xarray, size, cudaMemcpyHostToDevice);
    cudaMemcpy(device_y, yarray, size, cudaMemcpyHostToDevice);

    // run kernel
    double kernelStartTime = CycleTimer::currentSeconds();
    saxpy_kernel<<<blocks, threadsPerBlock>>>(N, alpha, device_x, device_y, device_result);
    cudaThreadSynchronize();
    double kernelEndTime = CycleTimer::currentSeconds();

    // copy result from GPU using cudaMemcpy
    cudaMemcpy(resultarray, device_result, size, cudaMemcpyDeviceToHost);

    // end timing after result has been copied back into host memory
    double endTime = CycleTimer::currentSeconds();

    cudaError_t errCode = cudaPeekAtLastError();
    if (errCode != cudaSuccess) {
        fprintf(stderr, "WARNING: A CUDA error occured: code=%d, %s\n", errCode, cudaGetErrorString(errCode));
    }

    double overallDuration = endTime - startTime;
    printf("Overall: %.3f ms\t\t[%.3f GB/s]\n", 1000.f * overallDuration, toBW(totalBytes, overallDuration));
    printf("Kernel run time: %.3f ms\n", 1000.f * (kernelEndTime-kernelStartTime));
    
    // free memory buffers on the GPU
    cudaFree(device_x);
    cudaFree(device_y);
    cudaFree(device_result);
}

__global__ void fg_kernel(double** u, double** v, double** f, double** g, int imax, int jmax, double delt, double delx, double dely, double gx, double gy, double gamma, double Re){
	int j = 
	int i = 
	// Note: loop start and stop conditions are different
	a=u[j][i]+u[j][i+1];
	b=u[j][i-1]+u[j][i];
	c=u[j][i]-u[j][i+1];
	d=u[j][i-1]-u[j][i];
	e=u[j][i]+u[j+1][i];
	ff=u[j-1][i]+u[j][i];
	gg=u[j][i]-u[j+1][i];
	h=u[j-1][i]-u[j][i];
	va=v[j][i]+v[j][i+1];
	vb=v[j-1][i]+v[j-1][i+1];
	u2x=1/delx*((a/2)*(a/2)-(b/2)*(b/2))
	+gamma*1/delx*(abs(a)/2*c/2-abs(b)/2*d/2);
	uvy=1/dely*(va/2*e/2-vb/2*ff/2)
	+gamma*1/dely*(abs(va)/2*gg/2-abs(vb)/2*h/2);
	u2x2=(u[j][i+1]-2*u[j][i]+u[j][i-1])/(delx*delx);
	u2y2=(u[j+1][i]-2*u[j][i]+u[j-1][i])/(dely*dely);
	f[j][i]=u[j][i]+delt*(1/Re*(u2x2+u2y2)-u2x-uvy+gx);

	a=v[j][i]+v[j][i+1];
	b=v[j][i-1]+v[j][i];
	c=v[j][i]-v[j][i+1];
	d=v[j][i-1]-v[j][i];
	e=v[j][i]+v[j+1][i];
	ff=v[j-1][i]+v[j][i];
	gg=v[j][i]-v[j+1][i];
	h=v[j-1][i]-v[j][i];
	ua=u[j][i]+u[j+1][i];
	ub=u[j][i-1]+u[j+1][i-1];
	uvx=1/delx*(ua/2*a/2-ub/2*b/2)
	+gamma*1/delx*(abs(ua)/2*c/2-abs(ub)/2*d/2);
	v2y=1/dely*((e/2)*(e/2)-(ff/2)*(ff/2))
	+gamma*1/dely*(abs(e)/2*gg/2-abs(ff)/2*h/2);
	v2x2=(v[j][i+1]-2*v[j][i]+v[j][i-1])/(delx*delx);
	v2y2=(v[j+1][i]-2*v[j][i]+v[j-1][i])/(dely*dely);
	g[j][i]=v[j][i]+delt*(1/Re*(v2x2+v2y2)-uvx-v2y+gy);
}

void comp_fg(double **u,double **v,double **f,double **g, int imax,int jmax,double delt,double delx,double dely,double gx,double gy,double gamma,double Re){
	int j,i;
	double a,b,c,d,e,ff,gg,h,va,vb,u2x,uvy,u2x2,u2y2;
	double ua,ub,uvx,v2y,v2x2,v2y2;
    for(j=0;j<jmax+2;j++){
        for(i=0;i<imax+2;i++){
            f[j][i]=0;
            g[j][i]=0;
        }
    }
	for(j=1;j<jmax+1;j++){
		f[j][0]=u[j][0];
		f[j][imax]=u[j][imax];
	}
	for(i=1;i<imax+1;i++){
		g[0][i]=v[0][i];
		g[jmax][i]=v[jmax][i];
	}
	fg_kernel<<<>>>(u,v,f,g,imax,jmax,delt,delx,dely,gx,gy,gamma,Re);
	return;
}

__global__ void rhs_kernel(double** f, double** g, double** rhs, int imax, int jmax, double delt, double delx, double dely){
	int j = 
	int i = 
	rhs[j][i]=1/delt*((f[j][i]-f[j][i-1])/delx+(g[j][i]-g[j-1][i])/dely);
}

void comp_rhs(double **f, double **g,double **rhs,int imax,int jmax,double delt,double delx,double dely){
	int j,i;
    for(j=0;j<jmax+2;j++){
        for(i=0;i<imax+2;i++){
            rhs[j][i]=0;
        }
    }
	rhs_kernel<<<>>>(f,g,rhs,imax,jmax,delt,delx,dely);
	return;
}

__global__ void poisson_kernel(double** p, double** rhs, double delx, double dely, double omg){
	int j = 
	int i = 
	eiw=1;eie=1;ejs=1;ejn=1;
	p[j][i]=(1-omg)*p[j][i]
	+omg/((eie+eiw)/(delx*delx)+(ejn+ejs)/(dely*dely))
	*((eie*p[j][i+1]+eiw*p[j][i-1])/(delx*delx)
	+(ejn*p[j+1][i]+ejs*p[j-1][i])/(dely*dely)-rhs[j][i]);

	int tmp =(eie*(p[j][i+1]-p[j][i])-eiw*(p[j][i]-p[j][i-1]))/(delx*delx)
	+(ejn*(p[j+1][i]-p[j][i])-ejs*(p[j][i]-p[j-1][i]))/(dely*dely)-rhs[j][i];
	sum[j][i] = tmp*tmp;
}

int poisson(double **p,double **rhs,int imax,int jmax,double delx,double dely,double eps,int itermax,double omg){
	int it,j,i,eiw,eie,ejs,ejn;
	double **r;
    double res;
    double** sum = RMATRIX(0, jmax+1,0,imax+1);
	for(it=0;it<itermax;it++){
		for(j=1;j<jmax+1;j++){
			p[j][0]=p[j][1];
			p[j][imax+1]=p[j][imax];
		}
		for(i=1;i<imax+1;i++){
			p[0][i]=p[1][i];
			p[jmax+1][i]=p[jmax][i];
		}
		poisson_kernel<<<>>>(p, rhs, imax, jmax, delx, dely, omg);
		int reduce_sum = kernelSum(sum);
        res=sqrt(reduce_sum/(imax*jmax));
        if(res<eps){
            printf("Converged...%f\n",res);
        	break;
        }
	}
	return it;
}

__global__ void adap_kernel(double** u, double** v, double** f, double** g, double** p, double delx, double dely){
	int j = 
	int i = 
	u[j][i]=f[j][i]-delt/delx*(p[j][i+1]-p[j][i]);
	v[j][i]=g[j][i]-delt/dely*(p[j+1][i]-p[j][i]);
}

void adap_uv(double **u,double **v,double **f,double **g,double **p,int imax,int jmax,double delt,double delx,double dely){
    adap_kernel<<<>>>(u,v,f,g,p,imax,jmax,delx,dely);
	return;
}
