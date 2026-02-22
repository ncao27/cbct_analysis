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

/*
% do ramp filtering
    for v = 1:nv
        
        % first do fft to frequency domain, then center projections on zero
        P = fftshift(fft(projections(:,v,view)));
        P = P .* ramp;

        % ifft the projects back
        projections(:,v,view) = real(ifft(ifftshift(P)));
    end*/

__global__ fdk_backproj_kernel(
    const float* proj,   // [nu * nv]
    float* recon,        // [nx * ny * nz]
    int nx, int ny, int nz,
    int nu, int nv,
    float dx, float dy, float dz,
    float du, float dv,
    float sid, float sdd,
    float theta
)
/*
The kernel that executes fdk. 
*/
{
    // Explanatio: in the mexFunction, we declared the grid dimensions (how 
    // many blocks there are) and the number of threads per block. So here 
    // we are basically setting the thread that is being used. This is equivalent
    // to a for loop except without explicity defining a for loop. The whole idea
    // behind CUDA is to define things voxel by voxel and not use for loops.
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    // safeguard to stop executing this thread if ix, iy, iz go out of bounds
    if (ix >= nx || iy >= ny || iz >= nz) return;

    // voxel position
    float x = (ix - nx/2.0f) * dx;
    float y = (iy - ny/2.0f) * dy;
    float z = (iz - nz/2.0f) * dz;

    // get the angles 
    float cosT = cosf(theta);
    float sinT = sinf(theta);

}

void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    mexPrintf("Hello from C++ Siddon!\n");

    float* d_proj  = (float*) mxGPUGetDataReadOnly(prhs[0]);
    float* d_recon = (float*) mxGPUGetData(prhs[1]);

    // get the number of grid points
    int nx = (int) mxGetScalar(prhs[2]);
    int ny = (int) mxGetScalar(prhs[3]);
    int nz = (int) mxGetScalar(prhs[4]);

    // get the number of detector plane grid points
    int nu = (int) mxGetScalar(prhs[5]);
    int nv = (int) mxGetScalar(prhs[6]);

    // get the grid point sizes of the volume
    float dx = (float) mxGetScalar(prhs[7]);
    float dy = (float) mxGetScalar(prhs[8]);
    float dz = (float) mxGetScalar(prhs[9]);

    // get the detector pixel sizes
    float du = (float) mxGetScalar(prhs[10]);
    float dv = (float) mxGetScalar(prhs[11]);

    // get the distance to center of vol and dist to detector
    float sid = (float) mxGetScalar(prhs[12]);
    float sdd = (float) mxGetScalar(prhs[13]);

    // get the angle associated with the view
    float theta = (float) mxGetScalar(prhs[14]);

    // this is CUDA specific, but basically, you launch a grid and each grid has a bunch of blocks
    // where each block has it's own threads and the threads do parallel computation

    // basically here 8*8*8 = 512, so each block has 512 threads
    dim3 block(8,8,8);
    dim3 grid(
        (nx+7)/8,
        (ny+7)/8,
        (nz+7)/8
    );

    fdk_backproj_kernel<<<grid, block>>>(
        d_proj, d_recon,
        nx, ny, nz,
        nu, nv,
        dx, dy, dz,
        du, dv,
        sid, sdd,
        theta
    );

    cudaDeviceSynchronize();

}