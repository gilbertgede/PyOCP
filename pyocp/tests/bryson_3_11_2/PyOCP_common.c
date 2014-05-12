
#include "PyOCP_common.h"

/* A numerical differentiation function */
Number num_diff(Number (*f)(Number *x), Number *x0, int var, int len_x0)
{
    Number *dx, *xll, *xl, *xu, *xuu;
    Number h0, temp;

    dx = (Number*)malloc(sizeof(Number) * len_x0);
    dx[var] = 1e-4;

    xll = (Number*)malloc(sizeof(Number) * len_x0);
    xl = (Number*)malloc(sizeof(Number) * len_x0);
    xu = (Number*)malloc(sizeof(Number) * len_x0);
    xuu = (Number*)malloc(sizeof(Number) * len_x0);
    for (int i = 0; i < len_x0; i++)
        xll[i] = x0[i];
    for (int i = 0; i < len_x0; i++)
        xl[i] = x0[i];
    for (int i = 0; i < len_x0; i++)
        xu[i] = x0[i];
    for (int i = 0; i < len_x0; i++)
        xuu[i] = x0[i];

    for (int i = 0; i < 5; i++)
    {
        xll[var] = x0[var] - 2 * dx[var];
        xl[var] = x0[var] - dx[var];
        xu[var] = x0[var] + dx[var];
        xuu[var] = x0[var] + 2 * dx[var];
        h0 = dx[var];
        temp = 0.5  / pow(h0, 3) * (-f(xll) + 2 * f(xl) - 2 * f(xu) + f(xuu));
        if (temp < DBL_EPSILON)
        {
            dx[var] = 1;
            break;
        }
        dx[var] = pow(3 * DBL_EPSILON / temp, 1.0 / 3.0);
        if (abs(h0 - dx[var]) < DBL_EPSILON)
            break;
    }
    h0 = dx[var];

    xu[var] = x0[var] + h0;
    xl[var] = x0[var] - h0;

    return 0.5 / h0 * (f(xu) - f(xl));
}

