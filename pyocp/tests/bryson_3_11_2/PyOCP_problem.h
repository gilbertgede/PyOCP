
#ifndef PYOCP_COMMON_H
  #include "PyOCP_common.h"
#endif
#ifndef PYOCP_PROBLEM_H
  #define PYOCP_PROBLEM_H
#endif

;


#include "Double_Integrator__0.h"
#include "Double_Integrator__1.h"
#include "Double_Integrator__2.h"


/* Function Declarations */
Bool eval_f(Index n, Number* x, Bool new_x,
            Number* obj_value, UserDataPtr user_data);

Bool eval_grad_f(Index n, Number* x, Bool new_x,
                 Number* grad_f, UserDataPtr user_data);

Bool eval_g(Index n, Number* x, Bool new_x,
            Index m, Number* g, UserDataPtr user_data);

Bool eval_jac_g(Index n, Number *x, Bool new_x,
                Index m, Index nele_jac,
                Index *iRow, Index *jCol, Number *values,
                UserDataPtr user_data);

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
            Index m, Number *lambda, Bool new_lambda,
            Index nele_hess, Index *iRow, Index *jCol,
            Number *values, UserDataPtr user_data);

/* User data structure definition (for user arguments) */
struct MyUserData
{
  Number l;
  };


