
#include "Double_Integrator__1.h"






/* Main Program */
/*
int main()
{
  Index n=-1;                          // number of variables
  Index m=-1;                          // number of constraints
  Number* x_L = NULL;                  // lower bounds on x
  Number* x_U = NULL;                  // upper bounds on x
  Number* g_L = NULL;                  // lower bounds on g
  Number* g_U = NULL;                  // upper bounds on g
  IpoptProblem nlp = NULL;             // IpoptProblem
  enum ApplicationReturnStatus status; // Solve return code
  Number* x = NULL;                    // starting point and solution vector
  Number* mult_g = NULL;               // constraint multipliers at the solution
  Number* mult_x_L = NULL;             // lower bound multipliers at the solution
  Number* mult_x_U = NULL;             // upper bound multipliers at the solution
  Number obj;                          // objective value
  Index i;                             // generic counter

  // Number of nonzeros in the Jacobian of the constraints
  Index nele_jac = 338;
  // Number of nonzeros in the Hessian of the Lagrangian (lower or
  // upper triangual part only)
  //Index nele_hess = 0; // To be updated once hessian is
  //written
  Index nele_hess = 0;
  // indexing style for matrices
  Index index_style = 0; // C-style; start counting of rows and column indices at 0

  // our user data for the function evalutions.
  struct MyUserData user_data;
  // Initialize the user data
    Number l = 0.111111111111111;
  user_data.l = l;


  // set the number of variables and allocate space for the bounds
  n=42;
  x_L = (Number*)malloc(sizeof(Number)*n);
  x_U = (Number*)malloc(sizeof(Number)*n);
  // set the values for the variable bounds
    x_L[0 + 0] = -HUGE_VAL;
  x_U[0 + 0] = HUGE_VAL;
  x_L[0 + 1] = -HUGE_VAL;
  x_U[0 + 1] = HUGE_VAL;
  x_L[0 + 2] = -HUGE_VAL;
  x_U[0 + 2] = HUGE_VAL;
  x_L[0 + 3] = -HUGE_VAL;
  x_U[0 + 3] = HUGE_VAL;
  x_L[0 + 4] = -HUGE_VAL;
  x_U[0 + 4] = HUGE_VAL;
  x_L[0 + 5] = -HUGE_VAL;
  x_U[0 + 5] = HUGE_VAL;
  x_L[0 + 6] = -HUGE_VAL;
  x_U[0 + 6] = HUGE_VAL;
  x_L[0 + 7] = -HUGE_VAL;
  x_U[0 + 7] = HUGE_VAL;
  x_L[0 + 8] = -HUGE_VAL;
  x_U[0 + 8] = HUGE_VAL;
  x_L[0 + 9] = -HUGE_VAL;
  x_U[0 + 9] = HUGE_VAL;
  x_L[0 + 10] = -HUGE_VAL;
  x_U[0 + 10] = HUGE_VAL;
  x_L[0 + 11] = -HUGE_VAL;
  x_U[0 + 11] = HUGE_VAL;
  x_L[0 + 12] = -HUGE_VAL;
  x_U[0 + 12] = HUGE_VAL;
  x_L[0 + 13] = -HUGE_VAL;
  x_U[0 + 13] = HUGE_VAL;
  x_L[0 + 14] = -HUGE_VAL;
  x_U[0 + 14] = HUGE_VAL;
  x_L[0 + 15] = -HUGE_VAL;
  x_U[0 + 15] = HUGE_VAL;
  x_L[0 + 16] = -HUGE_VAL;
  x_U[0 + 16] = HUGE_VAL;
  x_L[0 + 17] = -HUGE_VAL;
  x_U[0 + 17] = HUGE_VAL;
  x_L[0 + 18] = -HUGE_VAL;
  x_U[0 + 18] = HUGE_VAL;
  x_L[0 + 19] = -HUGE_VAL;
  x_U[0 + 19] = HUGE_VAL;
  x_L[0 + 20] = -HUGE_VAL;
  x_U[0 + 20] = HUGE_VAL;
  x_L[0 + 21] = -HUGE_VAL;
  x_U[0 + 21] = HUGE_VAL;
  x_L[0 + 22] = -HUGE_VAL;
  x_U[0 + 22] = HUGE_VAL;
  x_L[0 + 23] = -HUGE_VAL;
  x_U[0 + 23] = HUGE_VAL;
  x_L[0 + 24] = -HUGE_VAL;
  x_U[0 + 24] = HUGE_VAL;
  x_L[0 + 25] = -HUGE_VAL;
  x_U[0 + 25] = HUGE_VAL;
  x_L[0 + 26] = -HUGE_VAL;
  x_U[0 + 26] = HUGE_VAL;
  x_L[0 + 27] = -HUGE_VAL;
  x_U[0 + 27] = HUGE_VAL;
  x_L[0 + 28] = -HUGE_VAL;
  x_U[0 + 28] = HUGE_VAL;
  x_L[0 + 29] = -HUGE_VAL;
  x_U[0 + 29] = HUGE_VAL;
  x_L[0 + 30] = -HUGE_VAL;
  x_U[0 + 30] = HUGE_VAL;
  x_L[0 + 31] = -HUGE_VAL;
  x_U[0 + 31] = HUGE_VAL;
  x_L[0 + 32] = -HUGE_VAL;
  x_U[0 + 32] = HUGE_VAL;
  x_L[0 + 33] = -HUGE_VAL;
  x_U[0 + 33] = HUGE_VAL;
  x_L[0 + 34] = -HUGE_VAL;
  x_U[0 + 34] = HUGE_VAL;
  x_L[0 + 35] = -HUGE_VAL;
  x_U[0 + 35] = HUGE_VAL;
  x_L[0 + 36] = -HUGE_VAL;
  x_U[0 + 36] = HUGE_VAL;
  x_L[0 + 37] = -HUGE_VAL;
  x_U[0 + 37] = HUGE_VAL;
  x_L[0 + 38] = -HUGE_VAL;
  x_U[0 + 38] = HUGE_VAL;
  x_L[0 + 39] = -HUGE_VAL;
  x_U[0 + 39] = HUGE_VAL;
  x_L[0 + 40] = -HUGE_VAL;
  x_U[0 + 40] = HUGE_VAL;
  x_L[0 + 41] = -HUGE_VAL;
  x_U[0 + 41] = HUGE_VAL;


  // set the number of constraints and allocate space for the bounds
  m=45;
  g_L = (Number*)malloc(sizeof(Number)*m);
  g_U = (Number*)malloc(sizeof(Number)*m);
  // set the values of the constraint bounds
    g_L[0 + 0] = 0;
  g_U[0 + 0] = 0;
  g_L[0 + 1] = 0;
  g_U[0 + 1] = 0;
  g_L[0 + 2] = 0;
  g_U[0 + 2] = 0;
  g_L[0 + 3] = 0;
  g_U[0 + 3] = 0;
  g_L[0 + 4] = 0;
  g_U[0 + 4] = HUGE_VAL;
  g_L[0 + 5] = 0;
  g_U[0 + 5] = 0;
  g_L[0 + 6] = 0;
  g_U[0 + 6] = 0;
  g_L[0 + 7] = 0;
  g_U[0 + 7] = 0;
  g_L[0 + 8] = 0;
  g_U[0 + 8] = HUGE_VAL;
  g_L[0 + 9] = 0;
  g_U[0 + 9] = 0;
  g_L[0 + 10] = 0;
  g_U[0 + 10] = 0;
  g_L[0 + 11] = 0;
  g_U[0 + 11] = 0;
  g_L[0 + 12] = 0;
  g_U[0 + 12] = HUGE_VAL;
  g_L[0 + 13] = 0;
  g_U[0 + 13] = 0;
  g_L[0 + 14] = 0;
  g_U[0 + 14] = 0;
  g_L[0 + 15] = 0;
  g_U[0 + 15] = 0;
  g_L[0 + 16] = 0;
  g_U[0 + 16] = HUGE_VAL;
  g_L[0 + 17] = 0;
  g_U[0 + 17] = 0;
  g_L[0 + 18] = 0;
  g_U[0 + 18] = 0;
  g_L[0 + 19] = 0;
  g_U[0 + 19] = 0;
  g_L[0 + 20] = 0;
  g_U[0 + 20] = HUGE_VAL;
  g_L[0 + 21] = 0;
  g_U[0 + 21] = 0;
  g_L[0 + 22] = 0;
  g_U[0 + 22] = 0;
  g_L[0 + 23] = 0;
  g_U[0 + 23] = 0;
  g_L[0 + 24] = 0;
  g_U[0 + 24] = HUGE_VAL;
  g_L[0 + 25] = 0;
  g_U[0 + 25] = 0;
  g_L[0 + 26] = 0;
  g_U[0 + 26] = 0;
  g_L[0 + 27] = 0;
  g_U[0 + 27] = 0;
  g_L[0 + 28] = 0;
  g_U[0 + 28] = HUGE_VAL;
  g_L[0 + 29] = 0;
  g_U[0 + 29] = 0;
  g_L[0 + 30] = 0;
  g_U[0 + 30] = 0;
  g_L[0 + 31] = 0;
  g_U[0 + 31] = 0;
  g_L[0 + 32] = 0;
  g_U[0 + 32] = HUGE_VAL;
  g_L[0 + 33] = 0;
  g_U[0 + 33] = 0;
  g_L[0 + 34] = 0;
  g_U[0 + 34] = 0;
  g_L[0 + 35] = 0;
  g_U[0 + 35] = 0;
  g_L[0 + 36] = 0;
  g_U[0 + 36] = HUGE_VAL;
  g_L[0 + 37] = 0;
  g_U[0 + 37] = 0;
  g_L[0 + 38] = 0;
  g_U[0 + 38] = 0;
  g_L[0 + 39] = 0;
  g_U[0 + 39] = 0;
  g_L[0 + 40] = 0;
  g_U[0 + 40] = HUGE_VAL;
  g_L[0 + 41] = 0;
  g_U[0 + 41] = 0;
  g_L[0 + 42] = 0;
  g_U[0 + 42] = 0;
  g_L[0 + 43] = 0;
  g_U[0 + 43] = 0;
  g_L[0 + 44] = 0;
  g_U[0 + 44] = 0;


  // create the IpoptProblem
  nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
                           index_style, &Double_Integrator__1_eval_f,
                           &Double_Integrator__1_eval_g,
                           &Double_Integrator__1_eval_grad_f,
                           &Double_Integrator__1_eval_jac_g,
                           &Double_Integrator__1_eval_h);

  // We can free the memory now - the values for the bounds have been
  // copied internally in CreateIpoptProblem
  free(x_L);
  free(x_U);
  free(g_L);
  free(g_U);

  // Set some options.  Note the following ones are only examples,
  // they might not be suitable for your problem.
  AddIpoptNumOption(nlp, "tol", 1e-7);
  AddIpoptStrOption(nlp, "mu_strategy", "adaptive");
  AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");
  //AddIpoptStrOption(nlp, "derivative_test", "first-order");
  //AddIpoptStrOption(nlp, "output_file", "ipopt.out");

  // allocate space for the initial point and set the values
  x = (Number*)malloc(sizeof(Number)*n);
    x[0] = 0.111111111111111;
  x[1] = 0.0;
  x[2] = 2.00000000000000;
  x[3] = 0.0;
  x[4] = 0.111111111111111;
  x[5] = 0.0;
  x[6] = 2.00000000000000;
  x[7] = 0.0;
  x[8] = 0.111111111111111;
  x[9] = 0.0;
  x[10] = 2.00000000000000;
  x[11] = 0.0;
  x[12] = 0.111111111111111;
  x[13] = 0.0;
  x[14] = 2.00000000000000;
  x[15] = 0.0;
  x[16] = 0.111111111111111;
  x[17] = 0.0;
  x[18] = 2.00000000000000;
  x[19] = 0.0;
  x[20] = 0.111111111111111;
  x[21] = 0.0;
  x[22] = 2.00000000000000;
  x[23] = 0.0;
  x[24] = 0.111111111111111;
  x[25] = 0.0;
  x[26] = 2.00000000000000;
  x[27] = 0.0;
  x[28] = 0.111111111111111;
  x[29] = 0.0;
  x[30] = 2.00000000000000;
  x[31] = 0.0;
  x[32] = 0.111111111111111;
  x[33] = 0.0;
  x[34] = 2.00000000000000;
  x[35] = 0.0;
  x[36] = 0.111111111111111;
  x[37] = 0.0;
  x[38] = 2.00000000000000;
  x[39] = 0.0;
  x[40] = 0.300000000000000;
  x[41] = 0.700000000000000;


  // allocate space to store the bound multipliers at the solution
  mult_g = (Number*)malloc(sizeof(Number)*m);
  mult_x_L = (Number*)malloc(sizeof(Number)*n);
  mult_x_U = (Number*)malloc(sizeof(Number)*n);

  // solve the problem
  status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, &user_data);

  if (status == Solve_Succeeded) {
    printf("\n\nSolution of the primal variables, x\n");
    for (i=0; i<n; i++) {
      printf("x[%d] = %e\n", i, x[i]);
    }

    printf("\n\nSolution of the ccnstraint multipliers, lambda\n");
    for (i=0; i<m; i++) {
      printf("lambda[%d] = %e\n", i, mult_g[i]);
    }
    printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
    for (i=0; i<n; i++) {
      printf("z_L[%d] = %e\n", i, mult_x_L[i]);
    }
    for (i=0; i<n; i++) {
      printf("z_U[%d] = %e\n", i, mult_x_U[i]);
    }

    printf("\n\nObjective value\n");
    printf("f(x*) = %e\n", obj);

    FILE *fs;
    fs = fopen("Double_Integrator__1_results.csv", "w");
    if(fs == NULL){
        printf("Couldn't open an output file\n");
        return -1;
    }
    for (int iii = 0; iii < n; iii++)
        fprintf(fs, "%g, ", x[iii]);
    fprintf(fs, "\n");
    for (int iii = 0; iii < n; iii++)
        fprintf(fs, "%g, ", mult_x_L[iii]);
    fprintf(fs, "\n");
    for (int iii = 0; iii < n; iii++)
        fprintf(fs, "%g, ", mult_x_U[iii]);
    fprintf(fs, "\n");
    for (int iii = 0; iii < m; iii++)
        fprintf(fs, "%g, ", mult_g[iii]);
    fprintf(fs, "\n");
    fclose(fs);

  }
  else {
    printf("\n\nERROR OCCURRED DURING IPOPT OPTIMIZATION.\n");
  }

  // free allocated memory
  FreeIpoptProblem(nlp);
  free(x);
  free(mult_g);
  free(mult_x_L);
  free(mult_x_U);

  return (int)status;
}
*/

/* Function Implementations */
Bool Double_Integrator__1_eval_f(Index n, Number* X, Bool new_x,
            Number* obj_value, UserDataPtr user_data)
{
  struct MyUserData* my_data = user_data;

  assert(n == 42);

    Number l = my_data->l;


  *obj_value = 0;


  return TRUE;
}

Bool Double_Integrator__1_eval_grad_f(Index n, Number* X, Bool new_x,
                 Number* grad_f, UserDataPtr user_data)
{
  struct MyUserData* my_data = user_data;

  assert(n == 42);

    Number l = my_data->l;


    grad_f[0] = 0.0;
  grad_f[1] = 0.0;
  grad_f[2] = 0.0;
  grad_f[3] = 0.0;
  grad_f[4] = 0.0;
  grad_f[5] = 0.0;
  grad_f[6] = 0.0;
  grad_f[7] = 0.0;
  grad_f[8] = 0.0;
  grad_f[9] = 0.0;
  grad_f[10] = 0.0;
  grad_f[11] = 0.0;
  grad_f[12] = 0.0;
  grad_f[13] = 0.0;
  grad_f[14] = 0.0;
  grad_f[15] = 0.0;
  grad_f[16] = 0.0;
  grad_f[17] = 0.0;
  grad_f[18] = 0.0;
  grad_f[19] = 0.0;
  grad_f[20] = 0.0;
  grad_f[21] = 0.0;
  grad_f[22] = 0.0;
  grad_f[23] = 0.0;
  grad_f[24] = 0.0;
  grad_f[25] = 0.0;
  grad_f[26] = 0.0;
  grad_f[27] = 0.0;
  grad_f[28] = 0.0;
  grad_f[29] = 0.0;
  grad_f[30] = 0.0;
  grad_f[31] = 0.0;
  grad_f[32] = 0.0;
  grad_f[33] = 0.0;
  grad_f[34] = 0.0;
  grad_f[35] = 0.0;
  grad_f[36] = 0.0;
  grad_f[37] = 0.0;
  grad_f[38] = 0.0;
  grad_f[39] = 0.0;
  grad_f[40] = 0.0;
  grad_f[41] = 0.0;


  return TRUE;
}

Bool Double_Integrator__1_eval_g(Index n, Number* X, Bool new_x,
            Index m, Number* g, UserDataPtr user_data)
{
  struct MyUserData* my_data = user_data;

  assert(42 == n);
  assert(45 == m);

    Number l = my_data->l;


      Number t[10];
    Number grid[10];
    t[0] = X[40];
    t[1] = 0.888888888888889*X[40] + 0.111111111111111*X[41];
    t[2] = 0.777777777777778*X[40] + 0.222222222222222*X[41];
    t[3] = 0.666666666666667*X[40] + 0.333333333333333*X[41];
    t[4] = 0.555555555555556*X[40] + 0.444444444444444*X[41];
    t[5] = 0.444444444444444*X[40] + 0.555555555555556*X[41];
    t[6] = 0.333333333333333*X[40] + 0.666666666666667*X[41];
    t[7] = 0.222222222222222*X[40] + 0.777777777777778*X[41];
    t[8] = 0.111111111111111*X[40] + 0.888888888888889*X[41];
    t[9] = 1.0*X[41];
    grid[0] = 0.0;
    grid[1] = 0.111111111111111;
    grid[2] = 0.222222222222222;
    grid[3] = 0.333333333333333;
    grid[4] = 0.444444444444444;
    grid[5] = 0.555555555555556;
    grid[6] = 0.666666666666667;
    grid[7] = 0.777777777777778;
    grid[8] = 0.888888888888889;
    grid[9] = 1.00000000000000;


    g[0] += -X[0];
  g[1] += -X[1];
  g[2] += -X[2];
  g[3] += -X[40];
  g[41] += X[36];
  g[42] += X[37];
  g[43] += X[38];
  g[44] += X[41];
  Number local[7];
  for (int iii = 0; iii < 10; iii++)
  {
    local[0] = X[iii * 4 + 0];
    local[1] = X[iii * 4 + 1];
    local[2] = X[iii * 4 + 2];
    local[3] = X[iii * 4 + 3];
    local[4] = X[40];
    local[5] = X[41];
    local[6] = t[iii];
    g[iii * 4 + 0 + 4] = -local[0] + l;
  }
  Number d_local[12];
  for (int iii = 0; iii < 9; iii++)
  {
    d_local[0] = X[iii * 4 + 0];
    d_local[1] = X[iii * 4 + 1];
    d_local[2] = X[iii * 4 + 2];
    d_local[3] = X[iii * 4 + 3];
    d_local[4] = X[iii * 4 + 4];
    d_local[5] = X[iii * 4 + 5];
    d_local[6] = X[iii * 4 + 6];
    d_local[7] = X[iii * 4 + 7];
    d_local[8] = X[40];
    d_local[9] = X[41];
    d_local[10] = grid[iii];
    d_local[11] = grid[iii + 1];
    g[iii * 4 + 0 + 5] = -d_local[0] + d_local[4] - (-1.0L/6.0L*d_local[10]*(-d_local[8] + d_local[9]) + (1.0L/6.0L)*d_local[11]*(-d_local[8] + d_local[9]))*(3*d_local[1] + 3*d_local[5] + 4*(d_local[3] - d_local[7])*(-1.0L/8.0L*d_local[10]*(-d_local[8] + d_local[9]) + (1.0L/8.0L)*d_local[11]*(-d_local[8] + d_local[9])));
    g[iii * 4 + 1 + 5] = -d_local[1] + d_local[5] - (3*d_local[3] + 3*d_local[7])*(-1.0L/6.0L*d_local[10]*(-d_local[8] + d_local[9]) + (1.0L/6.0L)*d_local[11]*(-d_local[8] + d_local[9]));
    g[iii * 4 + 2 + 5] = -d_local[2] + d_local[6] - (-1.0L/6.0L*d_local[10]*(-d_local[8] + d_local[9]) + (1.0L/6.0L)*d_local[11]*(-d_local[8] + d_local[9]))*((1.0L/2.0L)*pow(d_local[3], 2) + (1.0L/2.0L)*pow(d_local[7], 2) + 2*pow((1.0L/2.0L)*d_local[3] + (1.0L/2.0L)*d_local[7], 2));
  }


  return TRUE;
}

Bool Double_Integrator__1_eval_jac_g(Index n, Number *X, Bool new_x,
                Index m, Index nele_jac,
                Index *iRow, Index *jCol, Number *values,
                UserDataPtr user_data)
{
  struct MyUserData* my_data = user_data;
  if (values == NULL) {
    /* return the structure of the jacobian */

        iRow[0] = 0;
    jCol[0] = 0;
    iRow[1] = 1;
    jCol[1] = 1;
    iRow[2] = 2;
    jCol[2] = 2;
    iRow[3] = 3;
    jCol[3] = 40;
    iRow[4] = 41;
    jCol[4] = 36;
    iRow[5] = 42;
    jCol[5] = 37;
    iRow[6] = 43;
    jCol[6] = 38;
    iRow[7] = 44;
    jCol[7] = 41;
    iRow[8] = 4;
    jCol[8] = 0;
    iRow[9] = 4;
    jCol[9] = 1;
    iRow[10] = 4;
    jCol[10] = 2;
    iRow[11] = 4;
    jCol[11] = 3;
    iRow[12] = 4;
    jCol[12] = 40;
    iRow[13] = 4;
    jCol[13] = 41;
    iRow[14] = 8;
    jCol[14] = 4;
    iRow[15] = 8;
    jCol[15] = 5;
    iRow[16] = 8;
    jCol[16] = 6;
    iRow[17] = 8;
    jCol[17] = 7;
    iRow[18] = 8;
    jCol[18] = 40;
    iRow[19] = 8;
    jCol[19] = 41;
    iRow[20] = 12;
    jCol[20] = 8;
    iRow[21] = 12;
    jCol[21] = 9;
    iRow[22] = 12;
    jCol[22] = 10;
    iRow[23] = 12;
    jCol[23] = 11;
    iRow[24] = 12;
    jCol[24] = 40;
    iRow[25] = 12;
    jCol[25] = 41;
    iRow[26] = 16;
    jCol[26] = 12;
    iRow[27] = 16;
    jCol[27] = 13;
    iRow[28] = 16;
    jCol[28] = 14;
    iRow[29] = 16;
    jCol[29] = 15;
    iRow[30] = 16;
    jCol[30] = 40;
    iRow[31] = 16;
    jCol[31] = 41;
    iRow[32] = 20;
    jCol[32] = 16;
    iRow[33] = 20;
    jCol[33] = 17;
    iRow[34] = 20;
    jCol[34] = 18;
    iRow[35] = 20;
    jCol[35] = 19;
    iRow[36] = 20;
    jCol[36] = 40;
    iRow[37] = 20;
    jCol[37] = 41;
    iRow[38] = 24;
    jCol[38] = 20;
    iRow[39] = 24;
    jCol[39] = 21;
    iRow[40] = 24;
    jCol[40] = 22;
    iRow[41] = 24;
    jCol[41] = 23;
    iRow[42] = 24;
    jCol[42] = 40;
    iRow[43] = 24;
    jCol[43] = 41;
    iRow[44] = 28;
    jCol[44] = 24;
    iRow[45] = 28;
    jCol[45] = 25;
    iRow[46] = 28;
    jCol[46] = 26;
    iRow[47] = 28;
    jCol[47] = 27;
    iRow[48] = 28;
    jCol[48] = 40;
    iRow[49] = 28;
    jCol[49] = 41;
    iRow[50] = 32;
    jCol[50] = 28;
    iRow[51] = 32;
    jCol[51] = 29;
    iRow[52] = 32;
    jCol[52] = 30;
    iRow[53] = 32;
    jCol[53] = 31;
    iRow[54] = 32;
    jCol[54] = 40;
    iRow[55] = 32;
    jCol[55] = 41;
    iRow[56] = 36;
    jCol[56] = 32;
    iRow[57] = 36;
    jCol[57] = 33;
    iRow[58] = 36;
    jCol[58] = 34;
    iRow[59] = 36;
    jCol[59] = 35;
    iRow[60] = 36;
    jCol[60] = 40;
    iRow[61] = 36;
    jCol[61] = 41;
    iRow[62] = 40;
    jCol[62] = 36;
    iRow[63] = 40;
    jCol[63] = 37;
    iRow[64] = 40;
    jCol[64] = 38;
    iRow[65] = 40;
    jCol[65] = 39;
    iRow[66] = 40;
    jCol[66] = 40;
    iRow[67] = 40;
    jCol[67] = 41;
    iRow[68] = 5;
    jCol[68] = 0;
    iRow[69] = 5;
    jCol[69] = 1;
    iRow[70] = 5;
    jCol[70] = 2;
    iRow[71] = 5;
    jCol[71] = 3;
    iRow[72] = 5;
    jCol[72] = 4;
    iRow[73] = 5;
    jCol[73] = 5;
    iRow[74] = 5;
    jCol[74] = 6;
    iRow[75] = 5;
    jCol[75] = 7;
    iRow[76] = 5;
    jCol[76] = 40;
    iRow[77] = 5;
    jCol[77] = 41;
    iRow[78] = 6;
    jCol[78] = 0;
    iRow[79] = 6;
    jCol[79] = 1;
    iRow[80] = 6;
    jCol[80] = 2;
    iRow[81] = 6;
    jCol[81] = 3;
    iRow[82] = 6;
    jCol[82] = 4;
    iRow[83] = 6;
    jCol[83] = 5;
    iRow[84] = 6;
    jCol[84] = 6;
    iRow[85] = 6;
    jCol[85] = 7;
    iRow[86] = 6;
    jCol[86] = 40;
    iRow[87] = 6;
    jCol[87] = 41;
    iRow[88] = 7;
    jCol[88] = 0;
    iRow[89] = 7;
    jCol[89] = 1;
    iRow[90] = 7;
    jCol[90] = 2;
    iRow[91] = 7;
    jCol[91] = 3;
    iRow[92] = 7;
    jCol[92] = 4;
    iRow[93] = 7;
    jCol[93] = 5;
    iRow[94] = 7;
    jCol[94] = 6;
    iRow[95] = 7;
    jCol[95] = 7;
    iRow[96] = 7;
    jCol[96] = 40;
    iRow[97] = 7;
    jCol[97] = 41;
    iRow[98] = 9;
    jCol[98] = 4;
    iRow[99] = 9;
    jCol[99] = 5;
    iRow[100] = 9;
    jCol[100] = 6;
    iRow[101] = 9;
    jCol[101] = 7;
    iRow[102] = 9;
    jCol[102] = 8;
    iRow[103] = 9;
    jCol[103] = 9;
    iRow[104] = 9;
    jCol[104] = 10;
    iRow[105] = 9;
    jCol[105] = 11;
    iRow[106] = 9;
    jCol[106] = 40;
    iRow[107] = 9;
    jCol[107] = 41;
    iRow[108] = 10;
    jCol[108] = 4;
    iRow[109] = 10;
    jCol[109] = 5;
    iRow[110] = 10;
    jCol[110] = 6;
    iRow[111] = 10;
    jCol[111] = 7;
    iRow[112] = 10;
    jCol[112] = 8;
    iRow[113] = 10;
    jCol[113] = 9;
    iRow[114] = 10;
    jCol[114] = 10;
    iRow[115] = 10;
    jCol[115] = 11;
    iRow[116] = 10;
    jCol[116] = 40;
    iRow[117] = 10;
    jCol[117] = 41;
    iRow[118] = 11;
    jCol[118] = 4;
    iRow[119] = 11;
    jCol[119] = 5;
    iRow[120] = 11;
    jCol[120] = 6;
    iRow[121] = 11;
    jCol[121] = 7;
    iRow[122] = 11;
    jCol[122] = 8;
    iRow[123] = 11;
    jCol[123] = 9;
    iRow[124] = 11;
    jCol[124] = 10;
    iRow[125] = 11;
    jCol[125] = 11;
    iRow[126] = 11;
    jCol[126] = 40;
    iRow[127] = 11;
    jCol[127] = 41;
    iRow[128] = 13;
    jCol[128] = 8;
    iRow[129] = 13;
    jCol[129] = 9;
    iRow[130] = 13;
    jCol[130] = 10;
    iRow[131] = 13;
    jCol[131] = 11;
    iRow[132] = 13;
    jCol[132] = 12;
    iRow[133] = 13;
    jCol[133] = 13;
    iRow[134] = 13;
    jCol[134] = 14;
    iRow[135] = 13;
    jCol[135] = 15;
    iRow[136] = 13;
    jCol[136] = 40;
    iRow[137] = 13;
    jCol[137] = 41;
    iRow[138] = 14;
    jCol[138] = 8;
    iRow[139] = 14;
    jCol[139] = 9;
    iRow[140] = 14;
    jCol[140] = 10;
    iRow[141] = 14;
    jCol[141] = 11;
    iRow[142] = 14;
    jCol[142] = 12;
    iRow[143] = 14;
    jCol[143] = 13;
    iRow[144] = 14;
    jCol[144] = 14;
    iRow[145] = 14;
    jCol[145] = 15;
    iRow[146] = 14;
    jCol[146] = 40;
    iRow[147] = 14;
    jCol[147] = 41;
    iRow[148] = 15;
    jCol[148] = 8;
    iRow[149] = 15;
    jCol[149] = 9;
    iRow[150] = 15;
    jCol[150] = 10;
    iRow[151] = 15;
    jCol[151] = 11;
    iRow[152] = 15;
    jCol[152] = 12;
    iRow[153] = 15;
    jCol[153] = 13;
    iRow[154] = 15;
    jCol[154] = 14;
    iRow[155] = 15;
    jCol[155] = 15;
    iRow[156] = 15;
    jCol[156] = 40;
    iRow[157] = 15;
    jCol[157] = 41;
    iRow[158] = 17;
    jCol[158] = 12;
    iRow[159] = 17;
    jCol[159] = 13;
    iRow[160] = 17;
    jCol[160] = 14;
    iRow[161] = 17;
    jCol[161] = 15;
    iRow[162] = 17;
    jCol[162] = 16;
    iRow[163] = 17;
    jCol[163] = 17;
    iRow[164] = 17;
    jCol[164] = 18;
    iRow[165] = 17;
    jCol[165] = 19;
    iRow[166] = 17;
    jCol[166] = 40;
    iRow[167] = 17;
    jCol[167] = 41;
    iRow[168] = 18;
    jCol[168] = 12;
    iRow[169] = 18;
    jCol[169] = 13;
    iRow[170] = 18;
    jCol[170] = 14;
    iRow[171] = 18;
    jCol[171] = 15;
    iRow[172] = 18;
    jCol[172] = 16;
    iRow[173] = 18;
    jCol[173] = 17;
    iRow[174] = 18;
    jCol[174] = 18;
    iRow[175] = 18;
    jCol[175] = 19;
    iRow[176] = 18;
    jCol[176] = 40;
    iRow[177] = 18;
    jCol[177] = 41;
    iRow[178] = 19;
    jCol[178] = 12;
    iRow[179] = 19;
    jCol[179] = 13;
    iRow[180] = 19;
    jCol[180] = 14;
    iRow[181] = 19;
    jCol[181] = 15;
    iRow[182] = 19;
    jCol[182] = 16;
    iRow[183] = 19;
    jCol[183] = 17;
    iRow[184] = 19;
    jCol[184] = 18;
    iRow[185] = 19;
    jCol[185] = 19;
    iRow[186] = 19;
    jCol[186] = 40;
    iRow[187] = 19;
    jCol[187] = 41;
    iRow[188] = 21;
    jCol[188] = 16;
    iRow[189] = 21;
    jCol[189] = 17;
    iRow[190] = 21;
    jCol[190] = 18;
    iRow[191] = 21;
    jCol[191] = 19;
    iRow[192] = 21;
    jCol[192] = 20;
    iRow[193] = 21;
    jCol[193] = 21;
    iRow[194] = 21;
    jCol[194] = 22;
    iRow[195] = 21;
    jCol[195] = 23;
    iRow[196] = 21;
    jCol[196] = 40;
    iRow[197] = 21;
    jCol[197] = 41;
    iRow[198] = 22;
    jCol[198] = 16;
    iRow[199] = 22;
    jCol[199] = 17;
    iRow[200] = 22;
    jCol[200] = 18;
    iRow[201] = 22;
    jCol[201] = 19;
    iRow[202] = 22;
    jCol[202] = 20;
    iRow[203] = 22;
    jCol[203] = 21;
    iRow[204] = 22;
    jCol[204] = 22;
    iRow[205] = 22;
    jCol[205] = 23;
    iRow[206] = 22;
    jCol[206] = 40;
    iRow[207] = 22;
    jCol[207] = 41;
    iRow[208] = 23;
    jCol[208] = 16;
    iRow[209] = 23;
    jCol[209] = 17;
    iRow[210] = 23;
    jCol[210] = 18;
    iRow[211] = 23;
    jCol[211] = 19;
    iRow[212] = 23;
    jCol[212] = 20;
    iRow[213] = 23;
    jCol[213] = 21;
    iRow[214] = 23;
    jCol[214] = 22;
    iRow[215] = 23;
    jCol[215] = 23;
    iRow[216] = 23;
    jCol[216] = 40;
    iRow[217] = 23;
    jCol[217] = 41;
    iRow[218] = 25;
    jCol[218] = 20;
    iRow[219] = 25;
    jCol[219] = 21;
    iRow[220] = 25;
    jCol[220] = 22;
    iRow[221] = 25;
    jCol[221] = 23;
    iRow[222] = 25;
    jCol[222] = 24;
    iRow[223] = 25;
    jCol[223] = 25;
    iRow[224] = 25;
    jCol[224] = 26;
    iRow[225] = 25;
    jCol[225] = 27;
    iRow[226] = 25;
    jCol[226] = 40;
    iRow[227] = 25;
    jCol[227] = 41;
    iRow[228] = 26;
    jCol[228] = 20;
    iRow[229] = 26;
    jCol[229] = 21;
    iRow[230] = 26;
    jCol[230] = 22;
    iRow[231] = 26;
    jCol[231] = 23;
    iRow[232] = 26;
    jCol[232] = 24;
    iRow[233] = 26;
    jCol[233] = 25;
    iRow[234] = 26;
    jCol[234] = 26;
    iRow[235] = 26;
    jCol[235] = 27;
    iRow[236] = 26;
    jCol[236] = 40;
    iRow[237] = 26;
    jCol[237] = 41;
    iRow[238] = 27;
    jCol[238] = 20;
    iRow[239] = 27;
    jCol[239] = 21;
    iRow[240] = 27;
    jCol[240] = 22;
    iRow[241] = 27;
    jCol[241] = 23;
    iRow[242] = 27;
    jCol[242] = 24;
    iRow[243] = 27;
    jCol[243] = 25;
    iRow[244] = 27;
    jCol[244] = 26;
    iRow[245] = 27;
    jCol[245] = 27;
    iRow[246] = 27;
    jCol[246] = 40;
    iRow[247] = 27;
    jCol[247] = 41;
    iRow[248] = 29;
    jCol[248] = 24;
    iRow[249] = 29;
    jCol[249] = 25;
    iRow[250] = 29;
    jCol[250] = 26;
    iRow[251] = 29;
    jCol[251] = 27;
    iRow[252] = 29;
    jCol[252] = 28;
    iRow[253] = 29;
    jCol[253] = 29;
    iRow[254] = 29;
    jCol[254] = 30;
    iRow[255] = 29;
    jCol[255] = 31;
    iRow[256] = 29;
    jCol[256] = 40;
    iRow[257] = 29;
    jCol[257] = 41;
    iRow[258] = 30;
    jCol[258] = 24;
    iRow[259] = 30;
    jCol[259] = 25;
    iRow[260] = 30;
    jCol[260] = 26;
    iRow[261] = 30;
    jCol[261] = 27;
    iRow[262] = 30;
    jCol[262] = 28;
    iRow[263] = 30;
    jCol[263] = 29;
    iRow[264] = 30;
    jCol[264] = 30;
    iRow[265] = 30;
    jCol[265] = 31;
    iRow[266] = 30;
    jCol[266] = 40;
    iRow[267] = 30;
    jCol[267] = 41;
    iRow[268] = 31;
    jCol[268] = 24;
    iRow[269] = 31;
    jCol[269] = 25;
    iRow[270] = 31;
    jCol[270] = 26;
    iRow[271] = 31;
    jCol[271] = 27;
    iRow[272] = 31;
    jCol[272] = 28;
    iRow[273] = 31;
    jCol[273] = 29;
    iRow[274] = 31;
    jCol[274] = 30;
    iRow[275] = 31;
    jCol[275] = 31;
    iRow[276] = 31;
    jCol[276] = 40;
    iRow[277] = 31;
    jCol[277] = 41;
    iRow[278] = 33;
    jCol[278] = 28;
    iRow[279] = 33;
    jCol[279] = 29;
    iRow[280] = 33;
    jCol[280] = 30;
    iRow[281] = 33;
    jCol[281] = 31;
    iRow[282] = 33;
    jCol[282] = 32;
    iRow[283] = 33;
    jCol[283] = 33;
    iRow[284] = 33;
    jCol[284] = 34;
    iRow[285] = 33;
    jCol[285] = 35;
    iRow[286] = 33;
    jCol[286] = 40;
    iRow[287] = 33;
    jCol[287] = 41;
    iRow[288] = 34;
    jCol[288] = 28;
    iRow[289] = 34;
    jCol[289] = 29;
    iRow[290] = 34;
    jCol[290] = 30;
    iRow[291] = 34;
    jCol[291] = 31;
    iRow[292] = 34;
    jCol[292] = 32;
    iRow[293] = 34;
    jCol[293] = 33;
    iRow[294] = 34;
    jCol[294] = 34;
    iRow[295] = 34;
    jCol[295] = 35;
    iRow[296] = 34;
    jCol[296] = 40;
    iRow[297] = 34;
    jCol[297] = 41;
    iRow[298] = 35;
    jCol[298] = 28;
    iRow[299] = 35;
    jCol[299] = 29;
    iRow[300] = 35;
    jCol[300] = 30;
    iRow[301] = 35;
    jCol[301] = 31;
    iRow[302] = 35;
    jCol[302] = 32;
    iRow[303] = 35;
    jCol[303] = 33;
    iRow[304] = 35;
    jCol[304] = 34;
    iRow[305] = 35;
    jCol[305] = 35;
    iRow[306] = 35;
    jCol[306] = 40;
    iRow[307] = 35;
    jCol[307] = 41;
    iRow[308] = 37;
    jCol[308] = 32;
    iRow[309] = 37;
    jCol[309] = 33;
    iRow[310] = 37;
    jCol[310] = 34;
    iRow[311] = 37;
    jCol[311] = 35;
    iRow[312] = 37;
    jCol[312] = 36;
    iRow[313] = 37;
    jCol[313] = 37;
    iRow[314] = 37;
    jCol[314] = 38;
    iRow[315] = 37;
    jCol[315] = 39;
    iRow[316] = 37;
    jCol[316] = 40;
    iRow[317] = 37;
    jCol[317] = 41;
    iRow[318] = 38;
    jCol[318] = 32;
    iRow[319] = 38;
    jCol[319] = 33;
    iRow[320] = 38;
    jCol[320] = 34;
    iRow[321] = 38;
    jCol[321] = 35;
    iRow[322] = 38;
    jCol[322] = 36;
    iRow[323] = 38;
    jCol[323] = 37;
    iRow[324] = 38;
    jCol[324] = 38;
    iRow[325] = 38;
    jCol[325] = 39;
    iRow[326] = 38;
    jCol[326] = 40;
    iRow[327] = 38;
    jCol[327] = 41;
    iRow[328] = 39;
    jCol[328] = 32;
    iRow[329] = 39;
    jCol[329] = 33;
    iRow[330] = 39;
    jCol[330] = 34;
    iRow[331] = 39;
    jCol[331] = 35;
    iRow[332] = 39;
    jCol[332] = 36;
    iRow[333] = 39;
    jCol[333] = 37;
    iRow[334] = 39;
    jCol[334] = 38;
    iRow[335] = 39;
    jCol[335] = 39;
    iRow[336] = 39;
    jCol[336] = 40;
    iRow[337] = 39;
    jCol[337] = 41;

  }
  else {
    /* return the values of the jacobian of the constraints */

      Number l = my_data->l;


        Number t[10];
    Number grid[10];
    t[0] = X[40];
    t[1] = 0.888888888888889*X[40] + 0.111111111111111*X[41];
    t[2] = 0.777777777777778*X[40] + 0.222222222222222*X[41];
    t[3] = 0.666666666666667*X[40] + 0.333333333333333*X[41];
    t[4] = 0.555555555555556*X[40] + 0.444444444444444*X[41];
    t[5] = 0.444444444444444*X[40] + 0.555555555555556*X[41];
    t[6] = 0.333333333333333*X[40] + 0.666666666666667*X[41];
    t[7] = 0.222222222222222*X[40] + 0.777777777777778*X[41];
    t[8] = 0.111111111111111*X[40] + 0.888888888888889*X[41];
    t[9] = 1.0*X[41];
    grid[0] = 0.0;
    grid[1] = 0.111111111111111;
    grid[2] = 0.222222222222222;
    grid[3] = 0.333333333333333;
    grid[4] = 0.444444444444444;
    grid[5] = 0.555555555555556;
    grid[6] = 0.666666666666667;
    grid[7] = 0.777777777777778;
    grid[8] = 0.888888888888889;
    grid[9] = 1.00000000000000;


        Number l_diffs[0];
    Number d_diffs[0];
    Number xtemp[0];
    Number local[7];
    Number d_local[12];
    Number diffs[0];



    values[0] = -1;
    values[1] = -1;
    values[2] = -1;
    values[3] = -1;
    values[4] = 1;
    values[5] = 1;
    values[6] = 1;
    values[7] = 1;
    for (int iii = 0; iii < 10; iii++)
    {
      local[0] = X[iii * 4 + 0];
      local[1] = X[iii * 4 + 1];
      local[2] = X[iii * 4 + 2];
      local[3] = X[iii * 4 + 3];
      local[4] = X[40];
      local[5] = X[41];
      local[6] = t[iii];
      values[8 + iii * 6] = -1;
      values[9 + iii * 6] = 0;
      values[10 + iii * 6] = 0;
      values[11 + iii * 6] = 0;
      values[12 + iii * 6] = 0;
      values[13 + iii * 6] = 0;
    }
    for (int iii = 0; iii < 9; iii++)
    {
      d_local[0] = X[iii * 4 + 0];
      d_local[1] = X[iii * 4 + 1];
      d_local[2] = X[iii * 4 + 2];
      d_local[3] = X[iii * 4 + 3];
      d_local[4] = X[iii * 4 + 4];
      d_local[5] = X[iii * 4 + 5];
      d_local[6] = X[iii * 4 + 6];
      d_local[7] = X[iii * 4 + 7];
      d_local[8] = X[40];
      d_local[9] = X[41];
      d_local[10] = grid[iii];
      d_local[11] = grid[iii + 1];
      values[68 + iii * 30] = -1;
      values[69 + iii * 30] = (1.0L/2.0L)*d_local[10]*(-d_local[8] + d_local[9]) - 1.0L/2.0L*d_local[11]*(-d_local[8] + d_local[9]);
      values[70 + iii * 30] = 0;
      values[71 + iii * 30] = -(-1.0L/2.0L*d_local[10]*(-d_local[8] + d_local[9]) + (1.0L/2.0L)*d_local[11]*(-d_local[8] + d_local[9]))*(-1.0L/6.0L*d_local[10]*(-d_local[8] + d_local[9]) + (1.0L/6.0L)*d_local[11]*(-d_local[8] + d_local[9]));
      values[72 + iii * 30] = 1;
      values[73 + iii * 30] = (1.0L/2.0L)*d_local[10]*(-d_local[8] + d_local[9]) - 1.0L/2.0L*d_local[11]*(-d_local[8] + d_local[9]);
      values[74 + iii * 30] = 0;
      values[75 + iii * 30] = -(-1.0L/6.0L*d_local[10]*(-d_local[8] + d_local[9]) + (1.0L/6.0L)*d_local[11]*(-d_local[8] + d_local[9]))*((1.0L/2.0L)*d_local[10]*(-d_local[8] + d_local[9]) - 1.0L/2.0L*d_local[11]*(-d_local[8] + d_local[9]));
      values[76 + iii * 30] = -4*((1.0L/8.0L)*d_local[10] - 1.0L/8.0L*d_local[11])*(d_local[3] - d_local[7])*(-1.0L/6.0L*d_local[10]*(-d_local[8] + d_local[9]) + (1.0L/6.0L)*d_local[11]*(-d_local[8] + d_local[9])) - ((1.0L/6.0L)*d_local[10] - 1.0L/6.0L*d_local[11])*(3*d_local[1] + 3*d_local[5] + 4*(d_local[3] - d_local[7])*(-1.0L/8.0L*d_local[10]*(-d_local[8] + d_local[9]) + (1.0L/8.0L)*d_local[11]*(-d_local[8] + d_local[9])));
      values[77 + iii * 30] = -(-1.0L/6.0L*d_local[10] + (1.0L/6.0L)*d_local[11])*(3*d_local[1] + 3*d_local[5] + 4*(d_local[3] - d_local[7])*(-1.0L/8.0L*d_local[10]*(-d_local[8] + d_local[9]) + (1.0L/8.0L)*d_local[11]*(-d_local[8] + d_local[9]))) - 4*(-1.0L/8.0L*d_local[10] + (1.0L/8.0L)*d_local[11])*(d_local[3] - d_local[7])*(-1.0L/6.0L*d_local[10]*(-d_local[8] + d_local[9]) + (1.0L/6.0L)*d_local[11]*(-d_local[8] + d_local[9]));
      values[78 + iii * 30] = 0;
      values[79 + iii * 30] = -1;
      values[80 + iii * 30] = 0;
      values[81 + iii * 30] = (1.0L/2.0L)*d_local[10]*(-d_local[8] + d_local[9]) - 1.0L/2.0L*d_local[11]*(-d_local[8] + d_local[9]);
      values[82 + iii * 30] = 0;
      values[83 + iii * 30] = 1;
      values[84 + iii * 30] = 0;
      values[85 + iii * 30] = (1.0L/2.0L)*d_local[10]*(-d_local[8] + d_local[9]) - 1.0L/2.0L*d_local[11]*(-d_local[8] + d_local[9]);
      values[86 + iii * 30] = -((1.0L/6.0L)*d_local[10] - 1.0L/6.0L*d_local[11])*(3*d_local[3] + 3*d_local[7]);
      values[87 + iii * 30] = -(-1.0L/6.0L*d_local[10] + (1.0L/6.0L)*d_local[11])*(3*d_local[3] + 3*d_local[7]);
      values[88 + iii * 30] = 0;
      values[89 + iii * 30] = 0;
      values[90 + iii * 30] = -1;
      values[91 + iii * 30] = -(2*d_local[3] + d_local[7])*(-1.0L/6.0L*d_local[10]*(-d_local[8] + d_local[9]) + (1.0L/6.0L)*d_local[11]*(-d_local[8] + d_local[9]));
      values[92 + iii * 30] = 0;
      values[93 + iii * 30] = 0;
      values[94 + iii * 30] = 1;
      values[95 + iii * 30] = -(d_local[3] + 2*d_local[7])*(-1.0L/6.0L*d_local[10]*(-d_local[8] + d_local[9]) + (1.0L/6.0L)*d_local[11]*(-d_local[8] + d_local[9]));
      values[96 + iii * 30] = -((1.0L/6.0L)*d_local[10] - 1.0L/6.0L*d_local[11])*((1.0L/2.0L)*pow(d_local[3], 2) + (1.0L/2.0L)*pow(d_local[7], 2) + 2*pow((1.0L/2.0L)*d_local[3] + (1.0L/2.0L)*d_local[7], 2));
      values[97 + iii * 30] = -(-1.0L/6.0L*d_local[10] + (1.0L/6.0L)*d_local[11])*((1.0L/2.0L)*pow(d_local[3], 2) + (1.0L/2.0L)*pow(d_local[7], 2) + 2*pow((1.0L/2.0L)*d_local[3] + (1.0L/2.0L)*d_local[7], 2));
    }


  }

  return TRUE;
}


Bool Double_Integrator__1_eval_h(Index n, Number *X, Bool new_x, Number obj_factor,
            Index m, Number *lambda, Bool new_lambda,
            Index nele_hess, Index *iRow, Index *jCol,
            Number *values, UserDataPtr user_data)
{
  return 0;
}

//
//Bool eval_h(Index n, Number *X, Bool new_x, Number obj_factor,
//            Index m, Number *lambda, Bool new_lambda,
//            Index nele_hess, Index *iRow, Index *jCol,
//            Number *values, UserDataPtr user_data)
//{
//  Index idx = 0; /* nonzero element counter */
//  Index row = 0; /* row counter for loop */
//  Index col = 0; /* col counter for loop */
//  if (values == NULL) {
//    /* return the structure. This is a symmetric matrix, fill the lower left
//     * triangle only. */
//
//    /* the hessian for this problem is actually dense */
//    idx=0;
//    for (row = 0; row < 4; row++) {
//      for (col = 0; col <= row; col++) {
//        iRow[idx] = row;
//        jCol[idx] = col;
//        idx++;
//      }
//    }
//
//    assert(idx == nele_hess);
//  }
//  else {
//    /* return the values. This is a symmetric matrix, fill the lower left
//     * triangle only */
//
//    /* fill the objective portion */
//    values[0] = obj_factor * (2*x[3]);               /* 0,0 */
//
//    values[1] = obj_factor * (x[3]);                 /* 1,0 */
//    values[2] = 0;                                   /* 1,1 */
//
//    values[3] = obj_factor * (x[3]);                 /* 2,0 */
//    values[4] = 0;                                   /* 2,1 */
//    values[5] = 0;                                   /* 2,2 */
//
//    values[6] = obj_factor * (2*x[0] + x[1] + x[2]); /* 3,0 */
//    values[7] = obj_factor * (x[0]);                 /* 3,1 */
//    values[8] = obj_factor * (x[0]);                 /* 3,2 */
//    values[9] = 0;                                   /* 3,3 */
//
//
//    /* add the portion for the first constraint */
//    values[1] += lambda[0] * (x[2] * x[3]);          /* 1,0 */
//
//    values[3] += lambda[0] * (x[1] * x[3]);          /* 2,0 */
//    values[4] += lambda[0] * (x[0] * x[3]);          /* 2,1 */
//
//    values[6] += lambda[0] * (x[1] * x[2]);          /* 3,0 */
//    values[7] += lambda[0] * (x[0] * x[2]);          /* 3,1 */
//    values[8] += lambda[0] * (x[0] * x[1]);          /* 3,2 */
//
//    /* add the portion for the second constraint */
//    values[0] += lambda[1] * 2;                      /* 0,0 */
//
//    values[2] += lambda[1] * 2;                      /* 1,1 */
//
//    values[5] += lambda[1] * 2;                      /* 2,2 */
//
//    values[9] += lambda[1] * 2;                      /* 3,3 */
//  }
//
//  return TRUE;
//}


