#include <stdio.h>
#include <stdlib.h>
#include "forbes_helpers.h"

/**
 * C_INDBOX calculates blah blah
 *
 *Arguments:
 * double *prox        : The value of the proximal operator at x
 * double *val         :
 * const double *x     :
 * const double *lower :
 * const double *upper :
 * const int dim       :
 *
 *Note: prox, x, upper and lower need to allocate at least `dim*sizeof(double)`
 *      bits in memory space.
 */
void c_indBox(
        double *prox,
        double *val, 
        const double *x, 
        const double *lower, 
        const double *upper,
        const int dim)
{    
    unsigned int i;    
    for (i=0; i<dim; i++){
        prox[i] = MIN(upper[i], MAX(lower[i], x[i]));
    }
    *val = 0.0;
    
}

