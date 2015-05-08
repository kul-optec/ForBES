#include "../c_indBox.c"

#define DIM 20

int main(void)
{
    double  prox[DIM],
            x[DIM],
            upper[DIM],
            lower[DIM],
            val;
    unsigned int i;
    
    for (i = 0; i < DIM; i++)
    {
        x[i] = (i-10)/2;
        upper[i] = 6;
        lower[i] = 1;
    }
    
    c_indBox(prox, &val, x, lower, upper, DIM);
    
    for (i = 0; i < DIM; i++)
    {
        printf("prox[%d] = %g\n", i, prox[i]);
    }
    printf("OK!\n");
}