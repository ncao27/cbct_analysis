#include "mex.h"
#include <cmath>
#include <vector>
#include <omp.h> 
#include <algorithm> // for std::max, std::min

void SiddonProject2D(
    const double* vol, mwSize nx, mwSize ny,
    const double* x_plane, const double* y_plane,
    const double dx, const double dy,
    const double xs, const double ys,
    const double* xd, const double* yd,
    double nu,
    double* prj
){
    const double x0 = x_plane[0];
    const double y0 = y_plane[0];

    // Removed collapse(2) since we only have one loop over nu now
    #pragma omp parallel for 
    for(int i = 0; i < (int)nu; i++){
        
        double dxr = xd[i] - xs;
        double dyr = yd[i] - ys;

        double L = sqrt(dxr*dxr + dyr*dyr);

        // Skip if source and detector are the same point
        if (fabs(dxr) < 1e-12 && fabs(dyr) < 1e-12) {
            prj[i] = 0.0;
            continue;
        }

        // Safeguard against perfectly vertical/horizontal parallel rays
        if (fabs(dxr) < 1e-12) dxr = 1e-12;
        if (fabs(dyr) < 1e-12) dyr = 1e-12;

        // Parametric boundaries
        double ax0 = (x_plane[0] - xs) / dxr;
        double ax1 = (x_plane[nx] - xs) / dxr;
        double ay0 = (y_plane[0] - ys) / dyr;
        double ay1 = (y_plane[ny] - ys) / dyr;

        double a_min = fmax(fmin(ax0, ax1), fmin(ay0, ay1));
        double a_max = fmin(fmax(ax0, ax1), fmax(ay0, ay1));

        if (a_min >= a_max) {
            prj[i] = 0.0;
            continue;
        }

        // -----------------------------
        // Initial pixel indices
        // -----------------------------
        double xm = xs + a_min * dxr;
        double ym = ys + a_min * dyr;

        int ix = (int)floor((xm - x0) / dx);
        int iy = (int)floor((ym - y0) / dy);

        // Clamp to 2D image boundaries
        ix = std::max(0, std::min(ix, (int)nx - 1));
        iy = std::max(0, std::min(iy, (int)ny - 1));
        
        // -----------------------------
        // Step directions
        // -----------------------------
        int sx = (dxr > 0) ? 1 : -1;
        int sy = (dyr > 0) ? 1 : -1;

        // -----------------------------
        // Next boundary distances
        // -----------------------------
        double ax = (x_plane[ix + (sx > 0)] - xs) / dxr;
        double ay = (y_plane[iy + (sy > 0)] - ys) / dyr;

        double dax = dx / fabs(dxr);
        double day = dy / fabs(dyr);

        // -----------------------------
        // Siddon traversal (2D)
        // -----------------------------
        double a = a_min;
        double sum = 0.0;

        while (a < a_max) {
            double a_next;

            if (ax <= ay) {
                a_next = ax;
                ax += dax;
                ix += sx;
            } else {
                a_next = ay;
                ay += day;
                iy += sy;
            }

            // Check if we are inside the 2D grid
            if (ix >= 0 && ix < (int)nx &&
                iy >= 0 && iy < (int)ny) {

                double len = (a_next - a) * L;
                sum += vol[ix + iy*nx] * len; // 2D flat indexing
            }

            a = a_next;
        }

        prj[i] = sum;
    }
}

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    // We now only need 10 inputs for 2D.
    if (nrhs != 10)
        mexErrMsgTxt("10 inputs required for 2D Siddon (vol, x_plane, y_plane, xs, ys, xd, yd, nu, dx, dy).");

    if (nlhs != 1)
        mexErrMsgTxt("1 output required.");

    // 0: Image (2D)
    const mxArray* vol_mx = prhs[0];
    const double* vol = mxGetPr(vol_mx); 
    const mwSize* volDims = mxGetDimensions(vol_mx);
    mwSize nx = volDims[0];
    mwSize ny = volDims[1];

    // 1 & 2: Grid Planes
    const double* x_plane = mxGetPr(prhs[1]);
    const double* y_plane = mxGetPr(prhs[2]);

    // 3 & 4: Source Position
    double xs = static_cast<double>(mxGetScalar(prhs[3]));
    double ys = static_cast<double>(mxGetScalar(prhs[4]));

    // 5 & 6: Detector pixel coordinates
    const double* xd = mxGetPr(prhs[5]);
    const double* yd = mxGetPr(prhs[6]);

    // 7: Number of detector pixels
    double nu = static_cast<double>(mxGetScalar(prhs[7]));

    // 8 & 9: Pixel spacings
    const double dx = *mxGetPr(prhs[8]);
    const double dy = *mxGetPr(prhs[9]);

    // Create 1D Output array (size: nu x 1)
    plhs[0] = mxCreateDoubleMatrix(nu, 1, mxREAL);
    double* prj = mxGetPr(plhs[0]);

    // Call the 2D projection function
    SiddonProject2D(
        vol, nx, ny,
        x_plane, y_plane,
        dx, dy,
        xs, ys,
        xd, yd,
        nu,
        prj
    );
}