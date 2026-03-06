#include "mex.h"
#include "cuda_runtime.h"
#include "gpu/mxGPUArray.h"
#include <cmath>
#include <vector>
#include <omp.h> 

// call the siddon_kernel
__global__ void siddon_kernel(
    const float* vol,    // [nx * ny * nz] volume on GPU
    float* proj,         // [nu * nv] projection output on GPU
    int nx, int ny, int nz,
    int nu, int nv,
    float dx, float dy, float dz,
    const double* x_plane,
    const double* y_plane,
    const double* z_plane,
    double xs,
    double ys,
    double zs,
    const double* xd,
    const double* yd,
    const double* zd
)
/*
Siddon kernel where we will write Siddon reconstruction-related code
*/
{   
    // lauch thread from blocks; confusing but threads literally enumerate pixels on the detector panel
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    const double x0 = x_plane[0];
    const double y0 = y_plane[0];
    const double z0 = z_plane[0];

    // safeguard to stop executing this thread if ix, iy, iz go out of bounds
    if (i >= nu || j >= nv) return;

    // find the distance from specific pixel to source
    double dxr = xd[i] - xs;
    double dyr = yd[i] - ys;
    double dzr = zd[j] - zs;
    double L = sqrt(dxr*dxr + dyr*dyr + dzr*dzr);

    // check to see if we are out of range or not
    if (fabs(dxr) < 1e-12 || fabs(dyr) < 1e-12 || fabs(dzr) < 1e-12) return;

    // calculate the entry and exit of every thread
    double ax0 = (x_plane[0] - xs) / dxr;
    double ax1 = (x_plane[nx] - xs) / dxr;
    double ay0 = (y_plane[0] - ys) / dyr;
    double ay1 = (y_plane[ny] - ys) / dyr;
    double az0 = (z_plane[0] - zs) / dzr;
    double az1 = (z_plane[nz] - zs) / dzr;

    double a_min = fmax(fmax(fmin(ax0, ax1), fmin(ay0, ay1)), fmin(az0, az1));
    double a_max = fmin(fmin(fmax(ax0, ax1), fmax(ay0, ay1)), fmax(az0, az1));

    if (a_min >= a_max) return;

    // -----------------------------
    // Initial voxel indices
    // -----------------------------
    double xm = xs + a_min * dxr;
    double ym = ys + a_min * dyr;
    double zm = zs + a_min * dzr;

    int ix = (int)floor((xm - x0) / dx);
    int iy = (int)floor((ym - y0) / dy);
    int iz = (int)floor((zm - z0) / dz);

    // Clamp to volume
    ix = max(0, min(ix, nx - 1));
    iy = max(0, min(iy, ny - 1));
    iz = max(0, min(iz, nz - 1));

    // -----------------------------
    // Step directions
    // -----------------------------
    int sx = (dxr > 0) ? 1 : -1;
    int sy = (dyr > 0) ? 1 : -1;
    int sz = (dzr > 0) ? 1 : -1;

    // -----------------------------
    // Next boundary distances
    // -----------------------------
    double ax = (x_plane[ix + (sx > 0)] - xs) / dxr;
    double ay = (y_plane[iy + (sy > 0)] - ys) / dyr;
    double az = (z_plane[iz + (sz > 0)] - zs) / dzr;

    double dax = dx / fabs(dxr);
    double day = dy / fabs(dyr);
    double daz = dz / fabs(dzr);

    // -----------------------------
    // Siddon traversal
    // -----------------------------
    double a = a_min;
    double sum = 0.0;

    while (a < a_max) {
        double a_next;

        if (ax <= ay && ax <= az) {
            a_next = ax;
            ax += dax;
            ix += sx;
        } else if (ay <= az) {
            a_next = ay;
            ay += day;
            iy += sy;
        } else {
            a_next = az;
            az += daz;
            iz += sz;
        }

        if (ix >= 0 && ix < (int)nx &&
            iy >= 0 && iy < (int)ny &&
            iz >= 0 && iz < (int)nz) {

            double len = (a_next - a) * L;
            sum += vol[ix + iy*nx + iz*nx*ny] * len;
        }

        a = a_next;
    }

    proj[(int) (i + j*nu)] = sum;
}


void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
    mxInitGPU();

    // Convert MATLAB volume into GPU array (volume, projection array with zeros)
    const mxGPUArray* volGPU  = mxGPUCreateFromMxArray(prhs[0]);
    mxGPUArray* projGPU = (mxGPUArray*) mxGPUCreateFromMxArray(prhs[1]);

    // Get GPU device pointers so that CUDA kernel an access the data
    const float* d_vol = (const float*) mxGPUGetDataReadOnly(volGPU);
    float* d_proj = (float*) mxGPUGetData(projGPU);

    // get dimensions of the volume
    int nx = (int)mxGetDimensions(prhs[0])[0];
    int ny = (int)mxGetDimensions(prhs[0])[1];
    int nz = (int)mxGetDimensions(prhs[0])[2];

    // get the x_plane, y_plane, and z_plane
    const double* x_plane = mxGetPr(prhs[2]);
    const double* y_plane = mxGetPr(prhs[3]);
    const double* z_plane = mxGetPr(prhs[4]);

    // get the source coordinates
    double xs = static_cast<double>(mxGetScalar(prhs[5]));
    double ys = static_cast<double>(mxGetScalar(prhs[6]));
    double zs = static_cast<double>(mxGetScalar(prhs[7]));

    //get the detector coordinates
    const double* xd = mxGetPr(prhs[8]);
    const double* yd = mxGetPr(prhs[9]);
    const double* zd = mxGetPr(prhs[10]);

    // get the number of detector plane grid points
    int nu = (int) mxGetScalar(prhs[11]);
    int nv = (int) mxGetScalar(prhs[12]);

    // get the grid point sizes of the volume
    const double dx = mxGetScalar(prhs[13]);
    const double dy = mxGetScalar(prhs[14]);
    const double dz = mxGetScalar(prhs[15]);

    // we have to make the matlab cpu arrays into gpu arrays
    double *d_x_plane, *d_y_plane, *d_z_plane;
    double *d_xd, *d_yd, *d_zd;

    cudaMalloc(&d_x_plane, (nx+1)*sizeof(double));
    cudaMalloc(&d_y_plane, (ny+1)*sizeof(double));
    cudaMalloc(&d_z_plane, (nz+1)*sizeof(double));

    cudaMemcpy(d_x_plane, x_plane, (nx+1)*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y_plane, y_plane, (ny+1)*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_z_plane, z_plane, (nz+1)*sizeof(double), cudaMemcpyHostToDevice);


    cudaMalloc(&d_xd, nu*sizeof(double));
    cudaMalloc(&d_yd, nu*sizeof(double));
    cudaMalloc(&d_zd, nv*sizeof(double));

    cudaMemcpy(d_xd, xd, nu*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_yd, yd, nu*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_zd, zd, nv*sizeof(double), cudaMemcpyHostToDevice);



    // we now define the number of threads and number of blocks used
    // matlab goes though each view, so realistically we only worry about one panel
    // 400 < 1024 threads per block, still valid
    dim3 block(16,16);
    dim3 grid(
        (nu + 15)/16,
        (nv + 15)/16
    );

    // call the siddon_kernel
    siddon_kernel<<<grid, block>>>(
        d_vol, d_proj,
        nx, ny, nz,
        nu, nv,
        dx, dy, dz,
        d_x_plane, d_y_plane, d_z_plane,
        xs, ys, zs,
        d_xd, d_yd, d_zd
    );

    // get real time feedback of where the error comes from
    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();

    if (err != cudaSuccess) {
        mexErrMsgIdAndTxt("CUDA:siddon_kernel",
            cudaGetErrorString(err));
    }
}