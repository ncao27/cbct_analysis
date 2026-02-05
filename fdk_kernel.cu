#include "mex.h"
#include "cuda_runtime.h"
#include <cmath>
#include <vector>
#include <omp.h> 

void weigh_dist()
/*
This function takes in the forward projection and weighs each pixel by its distance. 
*/
{

}

void ramp_filter()
/*
Performs a ramp filter by using the fft on the forward projections
*/
{

}

__global__ fdk_backproj_kernel()
/*
The kernel that executes fdk. 
*/
{

}

void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    float* d_proj  = (float*) mxGPUGetDataReadOnly(prhs[0]);
    float* d_recon = (float*) mxGPUGetData(prhs[1]);

    int nx = (int) mxGetScalar(prhs[2]);
    int ny = (int) mxGetScalar(prhs[3]);
    int nz = (int) mxGetScalar(prhs[4]);

    // this is CUDA specific, but basically, you launch a grid and each grid has a bunch of blocks
    // where each block has it's own threads and the threads do parallel computation
    dim3 block(8,8,8);
    dim3 grid(
        (nx+7)/8,
        (ny+7)/8,
        (nz+7)/8
    );

    fdk_backproj_kernel<<<grid, block>>>(
        d_proj, d_recon,
        nx, ny, nz,
    );

    cudaDeviceSynchronize();

}