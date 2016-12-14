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
fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_u2, (imax+2)*(jmax+2), 0);
fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_v2, (imax+2)*(jmax+2), 0);
fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_p2, (imax+2)*(jmax+2), 0);
fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_f2, (imax+2)*(jmax+2), 0);
fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_g2, (imax+2)*(jmax+2), 0);
fill_val<<<nBlocks, THREADSPB>>>(cudaDevice_rhs2, (imax+2)*(jmax+2), 0);

cublasHandle_t handle;
    cublasCreate(&handle);
    double mymax = 0.0;
    int max_idx = 0;
    cublasIdamax(handle, length, device_p, 1, &max_idx);
    double* tmp = (double*)malloc(sizeof(double));
    cudaMemcpy(tmp,&(device_p[max_idx]),sizeof(double),cudaMemcpyDeviceToHost);
    mymax = *tmp;
    free(tmp);
    return mymax;

void CFDCuda(int N, float alpha, float* xarray, float* yarray, float* resultarray) {

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

