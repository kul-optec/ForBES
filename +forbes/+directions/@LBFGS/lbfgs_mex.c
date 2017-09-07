#include "mex.h"
#include "../libLBFGS.h"

#define IS_REAL_SPARSE_MAT(P) (mxGetNumberOfDimensions(P) == 2 && \
    mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_DENSE_MAT(P) (mxGetNumberOfDimensions(P) == 2 && \
    !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_DENSE_VEC(P) ((mxGetNumberOfDimensions(P) == 1 || \
    (mxGetNumberOfDimensions(P) == 2 && (mxGetN(P) == 1 || mxGetM(P) == 1))) && \
    !mxIsSparse(P) && mxIsDouble(P))
#define IS_INT32_DENSE_VEC(P) ((mxGetNumberOfDimensions(P) == 1 || \
    (mxGetNumberOfDimensions(P) == 2 && (mxGetN(P) == 1 || mxGetM(P) == 1))) && \
    !mxIsSparse(P) && mxIsInt32(P))
#define IS_REAL_SCALAR(P) (IS_REAL_DENSE_VEC(P) && mxGetNumberOfElements(P) == 1)
#define IS_INT32_SCALAR(P) (IS_INT32_DENSE_VEC(P) && mxGetNumberOfElements(P) == 1)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int n, mem, curridx, currmem, dir_dims[2];
  double * dir, * s, * y, * ys, H, * g, * alpha;

  if (nrhs != 7) {
    mexErrMsgTxt("lbfgs_mex: you should provide exactly 7 arguments.");
    return;
  }
  if (nlhs > 1) {
    mexErrMsgTxt("lbfgs_mex: too many output arguments.");
    return;
  }
  if (!IS_REAL_DENSE_MAT(prhs[0])) {
    mexErrMsgTxt("lbfgs_mex: 1st argument must be a double, dense matrix.");
    return;
  }
  if (!IS_REAL_DENSE_MAT(prhs[1])) {
    mexErrMsgTxt("lbfgs_mex: 2nd argument must be a double, dense matrix.");
    return;
  }
  if (!IS_REAL_DENSE_VEC(prhs[2])) {
    mexErrMsgTxt("lbfgs_mex: 3rd argument must be a double, dense vector.");
    return;
  }
  if (!IS_REAL_SCALAR(prhs[3])) {
    mexErrMsgTxt("lbfgs_mex: 4rd argument must be a double scalar.");
    return;
  }
  if (!IS_REAL_DENSE_VEC(prhs[4])) {
    mexErrMsgTxt("lbfgs_mex: 5th argument must be a double, dense vector.");
    return;
  }
  if (!IS_INT32_SCALAR(prhs[5])) {
    mexErrMsgTxt("lbfgs_mex: 6th argument must be a 32-bit integer.");
    return;
  }
  if (!IS_INT32_SCALAR(prhs[6])) {
    mexErrMsgTxt("lbfgs_mex: 7th argument must be a 32-bit integer.");
    return;
  }

  s = mxGetPr(prhs[0]);
  y = mxGetPr(prhs[1]);
  ys = mxGetPr(prhs[2]);
  H = mxGetScalar(prhs[3]);
  g = mxGetPr(prhs[4]);
  curridx = (int)mxGetScalar(prhs[5])-1;
	currmem = (int)mxGetScalar(prhs[6]);

  n = mxGetDimensions(prhs[0])[0];
  mem = mxGetDimensions(prhs[0])[1];
  dir_dims[0] = n;
  dir_dims[1] = 1;

	alpha = mxCalloc(mem, sizeof(double));

  dir_dims[0] = n;
  dir_dims[1] = 1;
  plhs[0] = mxCreateNumericArray(2, dir_dims, mxDOUBLE_CLASS, mxREAL);
  dir = mxGetPr(plhs[0]);

  libLBFGS_buffer b;
  b.n = n;
  b.mem = mem;
  b.currmem = currmem;
  b.curridx = curridx;
  b.s_n_m = s;
  b.y_n_m = y;
  b.ys_m = ys;
  b.alpha_m = alpha;
  libLBFGS_matvec(&b, H, g, dir);

	mxFree(alpha);
}
