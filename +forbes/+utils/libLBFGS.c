#include <stdlib.h>
#include "libLBFGS.h"

libLBFGS_buffer * libLBFGS_init(int n, int mem) {
    libLBFGS_buffer * b = (libLBFGS_buffer *)malloc(sizeof(libLBFGS_buffer));
    b->n = n;
    b->mem = mem;
    b->currmem = 0;
    b->curridx = -1;
    b->s_n_m = (double *)malloc(mem*n*sizeof(double));
    b->y_n_m = (double *)malloc(mem*n*sizeof(double));
    b->ys_m = (double *)malloc(mem*sizeof(double));
    b->alpha_m = (double *)malloc(mem*sizeof(double));
    return b;
}

int libLBFGS_push(libLBFGS_buffer * b, double * s, double * y) {
    int i, base;
    b->curridx += 1;
    if (b->curridx >= b->mem) b->curridx = 0;
    b->currmem += 1;
    if (b->currmem > b->mem) b->currmem = b->mem;
    base = b->curridx*b->n;
    for (i=0; i<b->n; i++) {
        b->s_n_m[base+i] = s[i];
        b->y_n_m[base+i] = y[i];
        b->ys_m[b->curridx] += s[i]*y[i];
    }
    return 0;
}

int libLBFGS_matvec(libLBFGS_buffer * b, double H, double * g_n, double * dir_n) {
    double beta;
    int i, j, k;
    int n = b->n;
    int mem = b->mem;
    int currmem = b->currmem;
    int curridx = b->curridx;
    double * s_n_m = b->s_n_m;
    double * y_n_m = b->y_n_m;
    double * ys_m = b->ys_m;
    double * alpha_m = b->alpha_m;

    for (j=0; j<n; j++) dir_n[j] = g_n[j];

    i = curridx;
    for (k=0; k<currmem; k++) {
        alpha_m[i] = 0;
        for (j=0; j<n; j++) {
            alpha_m[i] += s_n_m[i*n+j]*dir_n[j];
        }
        alpha_m[i] /= ys_m[i];
        for (j=0; j<n; j++) {
            dir_n[j] = dir_n[j] - alpha_m[i]*y_n_m[i*n+j];
        }
        i = i-1;
        if (i<0) i = mem-1;
    }

    for (j=0; j<n; j++) {
        dir_n[j] = H*dir_n[j];
    }

    i = i+1;
    if (i>=mem) i = 0;
    for (k=0; k<currmem; k++) {
        beta = 0;
        for (j=0; j<n; j++) {
            beta += y_n_m[i*n+j]*dir_n[j];
        }
        beta /= ys_m[i];
        for (j=0; j<n; j++) {
            dir_n[j] = dir_n[j] + (alpha_m[i]-beta)*s_n_m[i*n+j];
        }
        i = i+1;
        if (i>=mem) i = 0;
    }
    return 0;
}

int libLBFGS_reset(libLBFGS_buffer * b) {
    b->currmem = 0;
    b->curridx = -1;
    return 0;
}

int libLBFGS_clear(libLBFGS_buffer * b) {
    free(b->s_n_m);
    free(b->y_n_m);
    free(b->ys_m);
    free(b->alpha_m);
    free(b);
    return 0;
}
