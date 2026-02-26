#include "mex.h"
#include "cuda_runtime.h"
#include "gpu/mxGPUArray.h"
#include <cmath>
#include <vector>
#include <omp.h> 


__global__ void fdk_backproj_kernel(
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
    // Explanation: in the mexFunction, we declared the grid dimensions (how 
    // many blocks there are) and the number of threads per block. So here 
    // we are basically setting the thread that is being used. This is equivalent
    // to a for loop except without explicity defining a for loop. The whole idea
    // behind CUDA is to define things voxel by voxel and not use for loops.
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    // safeguard to stop executing this thread if ix, iy, iz go out of bounds
    if (ix >= nx || iy >= ny || iz >= nz) return;

    // voxel position, 2.0f means float literal
    float x = (ix - (nx - 1)/2.0f) * dx;
    float y = (iy - (ny - 1)/2.0f) * dy;
    float z = (iz - (nz - 1)/2.0f) * dz;

    // get the angles, we do cosf to do 32-bit math which is faster
    float cosT = cosf(theta);
    float sinT = sinf(theta);

    // the denominator for CBCT, if it's less than 0 exit
    float denom = sid - x*cosT - y*sinT;
    if (denom <= 0.0f) return;

    // the u and v detector pixels that corresponds to the voxel
    float u = (sdd / denom) * (x*sinT - y*cosT) / du + (nu - 1)/2.0f;
    float v = (sdd / denom) * z / dv + (nv - 1)/2.0f;

    // round it because the above gives a floating point number
    int iu = (int) roundf(u);
    int iv = (int) roundf(v);
    /*
    int iu0 = floorf(u);
int iv0 = floorf(v);
float du = u - iu0;
float dv = v - iv0;*/

    // check if the calculations are actually within the bounds of the detector panel or not
    if (iu < 0 || iu >= nu || iv < 0 || iv >= nv) return;

    // find the index number corresponding to the projection and the volume 
    int projIdx = iu + iv * nu;
    int volIdx  = ix + iy * nx + iz * nx * ny;

    // define the weighting factor of the specific pixel
    float w = (sid * sid) / (denom * denom);

    // add the contribution
    recon[volIdx] += proj[projIdx] * w;
}

void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{

    mxInitGPU();

    // Convert MATLAB gpuArray inputs to mxGPUArray
    const mxGPUArray* projGPU  = mxGPUCreateFromMxArray(prhs[0]);
    mxGPUArray* reconGPU = (mxGPUArray*) mxGPUCreateFromMxArray(prhs[1]);

    // Get device pointers
    const float* d_proj = (const float*) mxGPUGetDataReadOnly(projGPU);
    float* d_recon = (float*) mxGPUGetData(reconGPU);

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