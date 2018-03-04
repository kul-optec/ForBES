/*
    Copyright (C) 2015, Lorenzo Stella and Panagiotis Patrinos

    This file is part of ForBES.

    ForBES is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ForBES is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with ForBES. If not, see <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "mex.h"

#define IS_REAL_SPARSE_MAT(P) (mxGetNumberOfDimensions(P) == 2 && \
    mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_DENSE_MAT(P) (mxGetNumberOfDimensions(P) == 2 && \
    !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_DENSE_3MAT(P) (mxGetNumberOfDimensions(P) == 3 && \
    !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_DENSE_VEC(P) ((mxGetNumberOfDimensions(P) == 1 || \
    (mxGetNumberOfDimensions(P) == 2 && (mxGetN(P) == 1 || mxGetM(P) == 1))) && \
    !mxIsSparse(P) && mxIsDouble(P))
#define IS_INT32_DENSE_VEC(P) ((mxGetNumberOfDimensions(P) == 1 || \
    (mxGetNumberOfDimensions(P) == 2 && (mxGetN(P) == 1 || mxGetM(P) == 1))) && \
    !mxIsSparse(P) && mxIsInt32(P))
#define IS_REAL_SCALAR(P) (IS_REAL_DENSE_VEC(P) && mxGetNumberOfElements(P) == 1)
#define IS_INT32_SCALAR(P) (IS_INT32_DENSE_VEC(P) && mxGetNumberOfElements(P) == 1)

void ZERO(int n, double * x)
{
    int i;
    for (i=0; i<n; i++) {
        x[i] = 0.0;
    }
}

void AXPY(int n, double * y, double a, double * x)
{
    int i;
    for (i=0; i<n; i++) {
        y[i] += a*x[i];
    }
}

void DENSE_MATVEC_INC(int m, int n, double * r, double * A, double * x)
{
    int i, j;
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            r[i] += A[j*m+i]*x[j];
        }
    }
}

void DENSE_TMATVEC_INC(int m, int n, double * r, double alpha, double * A, double * x)
{
    int i, j;
    double inc;
    for (j = 0; j < m; j++) {
        inc = 0.0;
        for (i = 0; i < n; i++) {
            inc += A[j*n+i]*x[i];
        }
        r[j] += alpha*inc;
    }
}

void L_SOLVE(double * L, double * x, int n)
{
    int i, j;
    double sum;
    *x = *x / *L;
    for (i = 1; i < n; i++) {
        sum = x[i];
        for (j = 0; j < i; j++) {
            sum -= (x[j] * L[j * n + i]);
        }
        x[i] = sum / L[i * (n + 1)];
    }
}

void Lt_SOLVE(double * L, double * x, int n)
{
    int i, j;
    double sum;    
    x[n - 1] = x[n - 1] / L[n * n - 1];
    for (i = 1; i < n; i++) {
        sum = x[n - 1 - i];
        for (j = 0; j < i; j++) {
            sum -= (L[n * (n - 1 - i) + (n - i + j)] * x[n - i + j]);
        }
        x[n - 1 - i] = sum / L[(n - 1 - i)*(1 + n)];
    }
}

/*** input  w has dimension (N+1)*n+N*m ***/
/*** output x has dimension (N+1)*n+N*m ***/
void RICCATI_SOLVE(int n, int m, int N, double * x, double * w, double * x0_n,
    double * A_n_n, double * B_n_m, double * LRs_m_m_N, double * Ks_m_n_N, 
    double * Ms_m_n_N, double * Ls_n_n_N, double * es_n_N, double * temp_n)
{
    double * ptr1, * ptr2, * ptr3; /* convenient */
    int i, j, mm = m*m, mn = m*n, nn = n*n, mpn = m+n;
    int size_x_w = (N+1)*n+N*m;
    ZERO(n*N, es_n_N);
    ZERO(N*(n+m)+n, x);
    for (i = 1; i <= n; i++) {
        es_n_N[n*N-i] = -w[size_x_w-i];
    }
    for (j = N-2; j >= 0; j--) {
        ptr1 = &(es_n_N[j*n]);
        AXPY(n, ptr1, -1.0, &(w[(j+1)*mpn]));
        DENSE_MATVEC_INC(n, n, ptr1, &(Ls_n_n_N[(j+1)*nn]), &(es_n_N[(j+1)*n]));
        DENSE_TMATVEC_INC(n, m, ptr1, -1.0, &(Ks_m_n_N[(j+1)*mn]), &(w[(j+1)*(mpn)+n]));
    }
    memcpy(x, x0_n, n*sizeof(double)); /* copy initial point */
    for (j = 0; j < N; j++) {
        ptr1 = &(x[j*(mpn)]);       /* points to x(j) */
        ptr2 = &(x[j*(mpn)+n]);     /* points to u(j) */
        ptr3 = &(x[(j+1)*(mpn)]);   /* points to x(j+1) */
        memcpy(ptr2, &(w[j*(mpn)+n]), m*sizeof(double)); /* put rhs in place */
        L_SOLVE(&(LRs_m_m_N[j*mm]), ptr2, m);        /* forward solve */
        Lt_SOLVE(&(LRs_m_m_N[j*mm]), ptr2, m);       /* backward solve*/
        DENSE_MATVEC_INC(m, n, ptr2, &(Ks_m_n_N[j*mn]), ptr1);
        DENSE_MATVEC_INC(m, n, ptr2, &(Ms_m_n_N[j*mn]), &(es_n_N[j*n]));
        DENSE_MATVEC_INC(n, n, ptr3, A_n_n, ptr1);
        DENSE_MATVEC_INC(n, m, ptr3, B_n_m, ptr2);
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n, m, N;
    double * x, * w, * x0, * A, * B, * LRs, * Ks, * Ms, * Ls;
    double * es, * temp_n;

    if (nrhs != 11) {
        mexErrMsgTxt("RiccatiSolve: you should provide exactly 11 arguments.");
        return;
    }
    if (nlhs > 2) {
        mexErrMsgTxt("RiccatiSolve: too many output arguments.");
        return;
    }
    if (!IS_REAL_DENSE_VEC(prhs[0])) {
        mexErrMsgTxt("RiccatiSolve: 1st argument must be a double, dense vector.");
        return;
    }
    if (!IS_REAL_DENSE_VEC(prhs[1])) {
        mexErrMsgTxt("RiccatiSolve: 2nd argument must be a double, dense vector.");
        return;
    }
    if (!IS_REAL_DENSE_MAT(prhs[2])) {
        mexErrMsgTxt("RiccatiSolve: 3rd argument must be a double, dense matrix.");
        return;
    }
    if (!IS_REAL_DENSE_MAT(prhs[3])) {
        mexErrMsgTxt("RiccatiSolve: 4th argument must be a double, dense matrix.");
        return;
    }
    if (!IS_REAL_DENSE_3MAT(prhs[4])) {
        mexErrMsgTxt("RiccatiSolve: 5th argument must be a double, dense 3D matrix.");
        return;
    }
    if (!IS_REAL_DENSE_3MAT(prhs[5])) {
        mexErrMsgTxt("RiccatiSolve: 6th argument must be a double, dense 3D matrix.");
        return;
    }
    if (!IS_REAL_DENSE_3MAT(prhs[6])) {
        mexErrMsgTxt("RiccatiSolve: 7th argument must be a double, dense 3D matrix.");
        return;
    }
    if (!IS_REAL_DENSE_3MAT(prhs[7])) {
        mexErrMsgTxt("RiccatiSolve: 8th argument must be a double, dense 3D matrix.");
        return;
    }
    if (!IS_INT32_SCALAR(prhs[8])) {
        mexErrMsgTxt("RiccatiSolve: 9th argument must be a 32-bit integer.");
        return;
    }
    if (!IS_INT32_SCALAR(prhs[9])) {
        mexErrMsgTxt("RiccatiSolve: 10th argument must be a 32-bit integer.");
        return;
    }
    if (!IS_INT32_SCALAR(prhs[10])) {
        mexErrMsgTxt("RiccatiSolve: 11th argument must be a 32-bit integer.");
        return;
    }

    w = mxGetPr(prhs[0]);
    x0 = mxGetPr(prhs[1]);
    A = mxGetPr(prhs[2]);
    B = mxGetPr(prhs[3]);
    LRs = mxGetPr(prhs[4]);
    Ks = mxGetPr(prhs[5]);
    Ms = mxGetPr(prhs[6]);
    Ls = mxGetPr(prhs[7]);
    n = mxGetScalar(prhs[8]);
    m = mxGetScalar(prhs[9]);
    N = mxGetScalar(prhs[10]);
    const mwSize x_dims[2] = { n ,1};
    plhs[0] = mxCreateDoubleScalar(0);
    plhs[1] = mxCreateNumericArray(2, x_dims, mxDOUBLE_CLASS, mxREAL);
    x = mxGetPr(plhs[1]);
    
    es = mxCalloc(n*N, sizeof(double));
    temp_n = mxCalloc(n, sizeof(double));
    
    RICCATI_SOLVE(n, m, N, x, w, x0, A, B, LRs, Ks, Ms, Ls, es, temp_n);

    mxFree(es);
    mxFree(temp_n);
}
