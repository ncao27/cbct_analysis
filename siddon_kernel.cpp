#include "mex.h"
#include <cmath>
#include <vector>
#include <omp.h> 

double* MakeArray(
    const double num, const double dim
){
    /*
    Make an array of size dim filled with num and return.  
    */
    double* arr = new double[(int)dim];
    for(int i = 0; i < dim; i++){
        arr[i] = num;
    }
    return arr;
}

void SiddonProject(
    const double* vol, mwSize nx, mwSize ny, mwSize nz,
    const double* x_plane, const double* y_plane, const double* z_plane,
    const double xnplane, const double ynplane, const double znplane,
    const double dx, const double dy, const double dz,
    const double xs, const double ys, const double zs,
    const double* xd, const double* yd, const double* zd,
    double nu, double nv,
    double* prj
){

    const double x0 = x_plane[0];
    const double y0 = y_plane[0];
    const double z0 = z_plane[0];

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nu; i++){
        for(int j = 0; j < nv; j++){
            double dxr = xd[i] - xs;
            double dyr = yd[i] - ys;
            double dzr = zd[j] - zs;

            double L = sqrt(dxr*dxr + dyr*dyr + dzr*dzr);

            if (fabs(dxr) < 1e-12 || fabs(dyr) < 1e-12 || fabs(dzr) < 1e-12) {
                continue;
            }

            // 0 is entry, 1 is exit
            double ax0 = (x_plane[0] - xs) / dxr;
            double ax1 = (x_plane[nx] - xs) / dxr;
            double ay0 = (y_plane[0] - ys) / dyr;
            double ay1 = (y_plane[ny] - ys) / dyr;
            double az0 = (z_plane[0] - zs) / dzr;
            double az1 = (z_plane[nz] - zs) / dzr;

            double a_min = fmax(fmax(fmin(ax0, ax1), fmin(ay0, ay1)), fmin(az0, az1));
            double a_max = fmin(fmin(fmax(ax0, ax1), fmax(ay0, ay1)), fmax(az0, az1));

            if (a_min >= a_max) {
                continue;
            }

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
            ix = std::max(0, std::min(ix, (int)nx - 1));
            iy = std::max(0, std::min(iy, (int)ny - 1));
            iz = std::max(0, std::min(iz, (int)nz - 1));
            
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

            prj[(int) (i + j*nu)] = sum;
        }
    }
}

// MATLAB always looks for a function named mexFunction. It’s like main() in C.
// Also, this entire function is written in mex syntax while SiddonProject is written in cpp
void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    mexPrintf("Hello from C++ Siddon!\n");

    // check if we really have the right number of inputs and outputs (12 inputs)
    if (nrhs != 15)
        mexErrMsgTxt("15 inputs required.");

    // check if outputs requirement is fulfilled
    if (nlhs != 1)
        mexErrMsgTxt("1 output required.");

    /* READ IN ALL INPUTS*/

    // read in the CT volume
    const mxArray* vol_mx = prhs[0];
    const double* vol = mxGetPr(vol_mx); // NOTE: newer matlab versions prefer mxGetDoubles, but mine is 2018 so we use getPr
    const mwSize* volDims = mxGetDimensions(vol_mx);
    mwSize nx = volDims[0];
    mwSize ny = volDims[1];
    mwSize nz = volDims[2];

    // get the x, y, and z planes
    const double* x_plane = mxGetPr(prhs[1]);
    const double* y_plane = mxGetPr(prhs[2]);
    const double* z_plane = mxGetPr(prhs[3]);

    // get dimensions of x, y, ad z planes
    double xnplane = static_cast<double>(mxGetNumberOfElements(prhs[1]));
    double ynplane = static_cast<double>(mxGetNumberOfElements(prhs[2]));
    double znplane = static_cast<double>(mxGetNumberOfElements(prhs[3]));

    // get the detector centers
    double xs = static_cast<double>(mxGetScalar(prhs[4]));
    double ys = static_cast<double>(mxGetScalar(prhs[5]));
    double zs = static_cast<double>(mxGetScalar(prhs[6]));

    // get the global coordinates of detector pixels
    const double* xd = mxGetPr(prhs[7]);
    const double* yd = mxGetPr(prhs[8]);
    const double* zd = mxGetPr(prhs[9]);

    // get the number of detector pixels in u and v directions
    double nu = static_cast<double>(mxGetScalar(prhs[10]));
    double nv = static_cast<double>(mxGetScalar(prhs[11]));

    // the spacing of pixels
    const double dx = *mxGetPr(prhs[12]);
    const double dy = *mxGetPr(prhs[13]);
    const double dz = *mxGetPr(prhs[14]);

    // define the output
    plhs[0] = mxCreateDoubleMatrix(nu, nv, mxREAL);
    double* prj = mxGetPr(plhs[0]);

    // call the siddon project helper function
    SiddonProject(
    vol, nx, ny, nz,
    x_plane, y_plane, z_plane,
    xnplane, ynplane, znplane,
    dx, dy, dz,
    xs, ys, zs,
    xd, yd, zd,
    nu, nv,
    prj
    );
}