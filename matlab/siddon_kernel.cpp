#include "mex.h"

void SiddonProject(
    const double* vol, mwSize nx, mwSize ny, mwSize nz,
    const double* x_plane, const double* y_plane, const double* z_plane,
    const double xs, const double ys, const double zs,
    const double* xd, const double* yd, const double* zd,
    double nu, double nv,
    double* prj
){

}

// MATLAB always looks for a function named mexFunction. It’s like main() in C.
void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    mexPrintf("Hello from C++ Siddon!\n");

    // check if we really have the right number of inputs and outputs (12 inputs)
    if (nrhs != 12)
        mexErrMsgTxt("12 inputs required.");

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
    mwSize xnplane = mxGetNumberOfElements(prhs[1]);
    mwSize ynplane = mxGetNumberOfElements(prhs[2]);
    mwSize znplane = mxGetNumberOfElements(prhs[3]);

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

    // define the output
    plhs[0] = mxCreateDoubleMatrix(nu, nv, mxREAL);
    double* prj = mxGetPr(plhs[0]);

    // call the siddon project helper function
    SiddonProject(
    vol, nx, ny, nz,
    x_plane, y_plane, z_plane,
    xs, ys, zs,
    xd, yd, zd,
    nu, nv,
    prj
    );
}
