
#include "Double_Integrator__0.h"






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
  Index nele_jac = 300;
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
  n=41;
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


  // set the number of constraints and allocate space for the bounds
  m=44;
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
  g_U[0 + 3] = HUGE_VAL;
  g_L[0 + 4] = 0;
  g_U[0 + 4] = 0;
  g_L[0 + 5] = 0;
  g_U[0 + 5] = 0;
  g_L[0 + 6] = 0;
  g_U[0 + 6] = 0;
  g_L[0 + 7] = 0;
  g_U[0 + 7] = HUGE_VAL;
  g_L[0 + 8] = 0;
  g_U[0 + 8] = 0;
  g_L[0 + 9] = 0;
  g_U[0 + 9] = 0;
  g_L[0 + 10] = 0;
  g_U[0 + 10] = 0;
  g_L[0 + 11] = 0;
  g_U[0 + 11] = HUGE_VAL;
  g_L[0 + 12] = 0;
  g_U[0 + 12] = 0;
  g_L[0 + 13] = 0;
  g_U[0 + 13] = 0;
  g_L[0 + 14] = 0;
  g_U[0 + 14] = 0;
  g_L[0 + 15] = 0;
  g_U[0 + 15] = HUGE_VAL;
  g_L[0 + 16] = 0;
  g_U[0 + 16] = 0;
  g_L[0 + 17] = 0;
  g_U[0 + 17] = 0;
  g_L[0 + 18] = 0;
  g_U[0 + 18] = 0;
  g_L[0 + 19] = 0;
  g_U[0 + 19] = HUGE_VAL;
  g_L[0 + 20] = 0;
  g_U[0 + 20] = 0;
  g_L[0 + 21] = 0;
  g_U[0 + 21] = 0;
  g_L[0 + 22] = 0;
  g_U[0 + 22] = 0;
  g_L[0 + 23] = 0;
  g_U[0 + 23] = HUGE_VAL;
  g_L[0 + 24] = 0;
  g_U[0 + 24] = 0;
  g_L[0 + 25] = 0;
  g_U[0 + 25] = 0;
  g_L[0 + 26] = 0;
  g_U[0 + 26] = 0;
  g_L[0 + 27] = 0;
  g_U[0 + 27] = HUGE_VAL;
  g_L[0 + 28] = 0;
  g_U[0 + 28] = 0;
  g_L[0 + 29] = 0;
  g_U[0 + 29] = 0;
  g_L[0 + 30] = 0;
  g_U[0 + 30] = 0;
  g_L[0 + 31] = 0;
  g_U[0 + 31] = HUGE_VAL;
  g_L[0 + 32] = 0;
  g_U[0 + 32] = 0;
  g_L[0 + 33] = 0;
  g_U[0 + 33] = 0;
  g_L[0 + 34] = 0;
  g_U[0 + 34] = 0;
  g_L[0 + 35] = 0;
  g_U[0 + 35] = HUGE_VAL;
  g_L[0 + 36] = 0;
  g_U[0 + 36] = 0;
  g_L[0 + 37] = 0;
  g_U[0 + 37] = 0;
  g_L[0 + 38] = 0;
  g_U[0 + 38] = 0;
  g_L[0 + 39] = 0;
  g_U[0 + 39] = HUGE_VAL;
  g_L[0 + 40] = 0;
  g_U[0 + 40] = 0;
  g_L[0 + 41] = 0;
  g_U[0 + 41] = 0;
  g_L[0 + 42] = 0;
  g_U[0 + 42] = 0;
  g_L[0 + 43] = 0;
  g_U[0 + 43] = 0;


  // create the IpoptProblem
  nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
                           index_style, &Double_Integrator__0_eval_f,
                           &Double_Integrator__0_eval_g,
                           &Double_Integrator__0_eval_grad_f,
                           &Double_Integrator__0_eval_jac_g,
                           &Double_Integrator__0_eval_h);

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
    x[0] = 0.0;
  x[1] = 1.00000000000000;
  x[2] = 0.0;
  x[3] = -6.00000000000000;
  x[4] = 0.0123456790123457;
  x[5] = 0.888888888888889;
  x[6] = 0.222222222222222;
  x[7] = -5.33333333333333;
  x[8] = 0.0246913580246914;
  x[9] = 0.777777777777778;
  x[10] = 0.444444444444444;
  x[11] = -4.66666666666667;
  x[12] = 0.0370370370370370;
  x[13] = 0.666666666666667;
  x[14] = 0.666666666666667;
  x[15] = -4.00000000000000;
  x[16] = 0.0493827160493827;
  x[17] = 0.555555555555556;
  x[18] = 0.888888888888889;
  x[19] = -3.33333333333333;
  x[20] = 0.0617283950617284;
  x[21] = 0.444444444444444;
  x[22] = 1.11111111111111;
  x[23] = -2.66666666666667;
  x[24] = 0.0740740740740741;
  x[25] = 0.333333333333333;
  x[26] = 1.33333333333333;
  x[27] = -2.00000000000000;
  x[28] = 0.0864197530864197;
  x[29] = 0.222222222222222;
  x[30] = 1.55555555555556;
  x[31] = -1.33333333333333;
  x[32] = 0.0987654320987654;
  x[33] = 0.111111111111111;
  x[34] = 1.77777777777778;
  x[35] = -0.666666666666667;
  x[36] = 0.111111111111111;
  x[37] = 0.0;
  x[38] = 2.00000000000000;
  x[39] = 0.0;
  x[40] = 0.300000000000000;


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
    fs = fopen("Double_Integrator__0_results.csv", "w");
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
Bool Double_Integrator__0_eval_f(Index n, Number* X, Bool new_x,
            Number* obj_value, UserDataPtr user_data)
{
  struct MyUserData* my_data = user_data;

  assert(n == 41);

    Number l = my_data->l;


  *obj_value = 0;


  return TRUE;
}

Bool Double_Integrator__0_eval_grad_f(Index n, Number* X, Bool new_x,
                 Number* grad_f, UserDataPtr user_data)
{
  struct MyUserData* my_data = user_data;

  assert(n == 41);

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


  return TRUE;
}

Bool Double_Integrator__0_eval_g(Index n, Number* X, Bool new_x,
            Index m, Number* g, UserDataPtr user_data)
{
  struct MyUserData* my_data = user_data;

  assert(41 == n);
  assert(44 == m);

    Number l = my_data->l;


      Number t[10];
    Number grid[10];
    t[0] = 0;
    t[1] = 0.111111111111111*X[40];
    t[2] = 0.222222222222222*X[40];
    t[3] = 0.333333333333333*X[40];
    t[4] = 0.444444444444444*X[40];
    t[5] = 0.555555555555556*X[40];
    t[6] = 0.666666666666667*X[40];
    t[7] = 0.777777777777778*X[40];
    t[8] = 0.888888888888889*X[40];
    t[9] = 1.0*X[40];
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


    g[0] = -X[0];
  g[1] = -X[1] + 1;
  g[2] = -X[2];
  g[40] += X[36];
  g[41] += X[37];
  g[42] += X[38];
  g[43] += X[40];
  Number local[6];
  for (int iii = 0; iii < 10; iii++)
  {
    local[0] = X[iii * 4 + 0];
    local[1] = X[iii * 4 + 1];
    local[2] = X[iii * 4 + 2];
    local[3] = X[iii * 4 + 3];
    local[4] = X[40];
    local[5] = t[iii];
    g[iii * 4 + 0 + 3] = -local[0] + l;
  }
  Number d_local[11];
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
    d_local[9] = grid[iii];
    d_local[10] = grid[iii + 1];
    g[iii * 4 + 0 + 4] = -d_local[0] + d_local[4] - ((1.0L/6.0L)*d_local[10]*d_local[8] - 1.0L/6.0L*d_local[8]*d_local[9])*(3*d_local[1] + 3*d_local[5] + 4*(d_local[3] - d_local[7])*((1.0L/8.0L)*d_local[10]*d_local[8] - 1.0L/8.0L*d_local[8]*d_local[9]));
    g[iii * 4 + 1 + 4] = -d_local[1] + d_local[5] - (3*d_local[3] + 3*d_local[7])*((1.0L/6.0L)*d_local[10]*d_local[8] - 1.0L/6.0L*d_local[8]*d_local[9]);
    g[iii * 4 + 2 + 4] = -d_local[2] + d_local[6] - ((1.0L/6.0L)*d_local[10]*d_local[8] - 1.0L/6.0L*d_local[8]*d_local[9])*((1.0L/2.0L)*pow(d_local[3], 2) + (1.0L/2.0L)*pow(d_local[7], 2) + 2*pow((1.0L/2.0L)*d_local[3] + (1.0L/2.0L)*d_local[7], 2));
  }


  return TRUE;
}

Bool Double_Integrator__0_eval_jac_g(Index n, Number *X, Bool new_x,
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
    iRow[3] = 40;
    jCol[3] = 36;
    iRow[4] = 41;
    jCol[4] = 37;
    iRow[5] = 42;
    jCol[5] = 38;
    iRow[6] = 43;
    jCol[6] = 40;
    iRow[7] = 3;
    jCol[7] = 0;
    iRow[8] = 3;
    jCol[8] = 1;
    iRow[9] = 3;
    jCol[9] = 2;
    iRow[10] = 3;
    jCol[10] = 3;
    iRow[11] = 3;
    jCol[11] = 40;
    iRow[12] = 7;
    jCol[12] = 4;
    iRow[13] = 7;
    jCol[13] = 5;
    iRow[14] = 7;
    jCol[14] = 6;
    iRow[15] = 7;
    jCol[15] = 7;
    iRow[16] = 7;
    jCol[16] = 40;
    iRow[17] = 11;
    jCol[17] = 8;
    iRow[18] = 11;
    jCol[18] = 9;
    iRow[19] = 11;
    jCol[19] = 10;
    iRow[20] = 11;
    jCol[20] = 11;
    iRow[21] = 11;
    jCol[21] = 40;
    iRow[22] = 15;
    jCol[22] = 12;
    iRow[23] = 15;
    jCol[23] = 13;
    iRow[24] = 15;
    jCol[24] = 14;
    iRow[25] = 15;
    jCol[25] = 15;
    iRow[26] = 15;
    jCol[26] = 40;
    iRow[27] = 19;
    jCol[27] = 16;
    iRow[28] = 19;
    jCol[28] = 17;
    iRow[29] = 19;
    jCol[29] = 18;
    iRow[30] = 19;
    jCol[30] = 19;
    iRow[31] = 19;
    jCol[31] = 40;
    iRow[32] = 23;
    jCol[32] = 20;
    iRow[33] = 23;
    jCol[33] = 21;
    iRow[34] = 23;
    jCol[34] = 22;
    iRow[35] = 23;
    jCol[35] = 23;
    iRow[36] = 23;
    jCol[36] = 40;
    iRow[37] = 27;
    jCol[37] = 24;
    iRow[38] = 27;
    jCol[38] = 25;
    iRow[39] = 27;
    jCol[39] = 26;
    iRow[40] = 27;
    jCol[40] = 27;
    iRow[41] = 27;
    jCol[41] = 40;
    iRow[42] = 31;
    jCol[42] = 28;
    iRow[43] = 31;
    jCol[43] = 29;
    iRow[44] = 31;
    jCol[44] = 30;
    iRow[45] = 31;
    jCol[45] = 31;
    iRow[46] = 31;
    jCol[46] = 40;
    iRow[47] = 35;
    jCol[47] = 32;
    iRow[48] = 35;
    jCol[48] = 33;
    iRow[49] = 35;
    jCol[49] = 34;
    iRow[50] = 35;
    jCol[50] = 35;
    iRow[51] = 35;
    jCol[51] = 40;
    iRow[52] = 39;
    jCol[52] = 36;
    iRow[53] = 39;
    jCol[53] = 37;
    iRow[54] = 39;
    jCol[54] = 38;
    iRow[55] = 39;
    jCol[55] = 39;
    iRow[56] = 39;
    jCol[56] = 40;
    iRow[57] = 4;
    jCol[57] = 0;
    iRow[58] = 4;
    jCol[58] = 1;
    iRow[59] = 4;
    jCol[59] = 2;
    iRow[60] = 4;
    jCol[60] = 3;
    iRow[61] = 4;
    jCol[61] = 4;
    iRow[62] = 4;
    jCol[62] = 5;
    iRow[63] = 4;
    jCol[63] = 6;
    iRow[64] = 4;
    jCol[64] = 7;
    iRow[65] = 4;
    jCol[65] = 40;
    iRow[66] = 5;
    jCol[66] = 0;
    iRow[67] = 5;
    jCol[67] = 1;
    iRow[68] = 5;
    jCol[68] = 2;
    iRow[69] = 5;
    jCol[69] = 3;
    iRow[70] = 5;
    jCol[70] = 4;
    iRow[71] = 5;
    jCol[71] = 5;
    iRow[72] = 5;
    jCol[72] = 6;
    iRow[73] = 5;
    jCol[73] = 7;
    iRow[74] = 5;
    jCol[74] = 40;
    iRow[75] = 6;
    jCol[75] = 0;
    iRow[76] = 6;
    jCol[76] = 1;
    iRow[77] = 6;
    jCol[77] = 2;
    iRow[78] = 6;
    jCol[78] = 3;
    iRow[79] = 6;
    jCol[79] = 4;
    iRow[80] = 6;
    jCol[80] = 5;
    iRow[81] = 6;
    jCol[81] = 6;
    iRow[82] = 6;
    jCol[82] = 7;
    iRow[83] = 6;
    jCol[83] = 40;
    iRow[84] = 8;
    jCol[84] = 4;
    iRow[85] = 8;
    jCol[85] = 5;
    iRow[86] = 8;
    jCol[86] = 6;
    iRow[87] = 8;
    jCol[87] = 7;
    iRow[88] = 8;
    jCol[88] = 8;
    iRow[89] = 8;
    jCol[89] = 9;
    iRow[90] = 8;
    jCol[90] = 10;
    iRow[91] = 8;
    jCol[91] = 11;
    iRow[92] = 8;
    jCol[92] = 40;
    iRow[93] = 9;
    jCol[93] = 4;
    iRow[94] = 9;
    jCol[94] = 5;
    iRow[95] = 9;
    jCol[95] = 6;
    iRow[96] = 9;
    jCol[96] = 7;
    iRow[97] = 9;
    jCol[97] = 8;
    iRow[98] = 9;
    jCol[98] = 9;
    iRow[99] = 9;
    jCol[99] = 10;
    iRow[100] = 9;
    jCol[100] = 11;
    iRow[101] = 9;
    jCol[101] = 40;
    iRow[102] = 10;
    jCol[102] = 4;
    iRow[103] = 10;
    jCol[103] = 5;
    iRow[104] = 10;
    jCol[104] = 6;
    iRow[105] = 10;
    jCol[105] = 7;
    iRow[106] = 10;
    jCol[106] = 8;
    iRow[107] = 10;
    jCol[107] = 9;
    iRow[108] = 10;
    jCol[108] = 10;
    iRow[109] = 10;
    jCol[109] = 11;
    iRow[110] = 10;
    jCol[110] = 40;
    iRow[111] = 12;
    jCol[111] = 8;
    iRow[112] = 12;
    jCol[112] = 9;
    iRow[113] = 12;
    jCol[113] = 10;
    iRow[114] = 12;
    jCol[114] = 11;
    iRow[115] = 12;
    jCol[115] = 12;
    iRow[116] = 12;
    jCol[116] = 13;
    iRow[117] = 12;
    jCol[117] = 14;
    iRow[118] = 12;
    jCol[118] = 15;
    iRow[119] = 12;
    jCol[119] = 40;
    iRow[120] = 13;
    jCol[120] = 8;
    iRow[121] = 13;
    jCol[121] = 9;
    iRow[122] = 13;
    jCol[122] = 10;
    iRow[123] = 13;
    jCol[123] = 11;
    iRow[124] = 13;
    jCol[124] = 12;
    iRow[125] = 13;
    jCol[125] = 13;
    iRow[126] = 13;
    jCol[126] = 14;
    iRow[127] = 13;
    jCol[127] = 15;
    iRow[128] = 13;
    jCol[128] = 40;
    iRow[129] = 14;
    jCol[129] = 8;
    iRow[130] = 14;
    jCol[130] = 9;
    iRow[131] = 14;
    jCol[131] = 10;
    iRow[132] = 14;
    jCol[132] = 11;
    iRow[133] = 14;
    jCol[133] = 12;
    iRow[134] = 14;
    jCol[134] = 13;
    iRow[135] = 14;
    jCol[135] = 14;
    iRow[136] = 14;
    jCol[136] = 15;
    iRow[137] = 14;
    jCol[137] = 40;
    iRow[138] = 16;
    jCol[138] = 12;
    iRow[139] = 16;
    jCol[139] = 13;
    iRow[140] = 16;
    jCol[140] = 14;
    iRow[141] = 16;
    jCol[141] = 15;
    iRow[142] = 16;
    jCol[142] = 16;
    iRow[143] = 16;
    jCol[143] = 17;
    iRow[144] = 16;
    jCol[144] = 18;
    iRow[145] = 16;
    jCol[145] = 19;
    iRow[146] = 16;
    jCol[146] = 40;
    iRow[147] = 17;
    jCol[147] = 12;
    iRow[148] = 17;
    jCol[148] = 13;
    iRow[149] = 17;
    jCol[149] = 14;
    iRow[150] = 17;
    jCol[150] = 15;
    iRow[151] = 17;
    jCol[151] = 16;
    iRow[152] = 17;
    jCol[152] = 17;
    iRow[153] = 17;
    jCol[153] = 18;
    iRow[154] = 17;
    jCol[154] = 19;
    iRow[155] = 17;
    jCol[155] = 40;
    iRow[156] = 18;
    jCol[156] = 12;
    iRow[157] = 18;
    jCol[157] = 13;
    iRow[158] = 18;
    jCol[158] = 14;
    iRow[159] = 18;
    jCol[159] = 15;
    iRow[160] = 18;
    jCol[160] = 16;
    iRow[161] = 18;
    jCol[161] = 17;
    iRow[162] = 18;
    jCol[162] = 18;
    iRow[163] = 18;
    jCol[163] = 19;
    iRow[164] = 18;
    jCol[164] = 40;
    iRow[165] = 20;
    jCol[165] = 16;
    iRow[166] = 20;
    jCol[166] = 17;
    iRow[167] = 20;
    jCol[167] = 18;
    iRow[168] = 20;
    jCol[168] = 19;
    iRow[169] = 20;
    jCol[169] = 20;
    iRow[170] = 20;
    jCol[170] = 21;
    iRow[171] = 20;
    jCol[171] = 22;
    iRow[172] = 20;
    jCol[172] = 23;
    iRow[173] = 20;
    jCol[173] = 40;
    iRow[174] = 21;
    jCol[174] = 16;
    iRow[175] = 21;
    jCol[175] = 17;
    iRow[176] = 21;
    jCol[176] = 18;
    iRow[177] = 21;
    jCol[177] = 19;
    iRow[178] = 21;
    jCol[178] = 20;
    iRow[179] = 21;
    jCol[179] = 21;
    iRow[180] = 21;
    jCol[180] = 22;
    iRow[181] = 21;
    jCol[181] = 23;
    iRow[182] = 21;
    jCol[182] = 40;
    iRow[183] = 22;
    jCol[183] = 16;
    iRow[184] = 22;
    jCol[184] = 17;
    iRow[185] = 22;
    jCol[185] = 18;
    iRow[186] = 22;
    jCol[186] = 19;
    iRow[187] = 22;
    jCol[187] = 20;
    iRow[188] = 22;
    jCol[188] = 21;
    iRow[189] = 22;
    jCol[189] = 22;
    iRow[190] = 22;
    jCol[190] = 23;
    iRow[191] = 22;
    jCol[191] = 40;
    iRow[192] = 24;
    jCol[192] = 20;
    iRow[193] = 24;
    jCol[193] = 21;
    iRow[194] = 24;
    jCol[194] = 22;
    iRow[195] = 24;
    jCol[195] = 23;
    iRow[196] = 24;
    jCol[196] = 24;
    iRow[197] = 24;
    jCol[197] = 25;
    iRow[198] = 24;
    jCol[198] = 26;
    iRow[199] = 24;
    jCol[199] = 27;
    iRow[200] = 24;
    jCol[200] = 40;
    iRow[201] = 25;
    jCol[201] = 20;
    iRow[202] = 25;
    jCol[202] = 21;
    iRow[203] = 25;
    jCol[203] = 22;
    iRow[204] = 25;
    jCol[204] = 23;
    iRow[205] = 25;
    jCol[205] = 24;
    iRow[206] = 25;
    jCol[206] = 25;
    iRow[207] = 25;
    jCol[207] = 26;
    iRow[208] = 25;
    jCol[208] = 27;
    iRow[209] = 25;
    jCol[209] = 40;
    iRow[210] = 26;
    jCol[210] = 20;
    iRow[211] = 26;
    jCol[211] = 21;
    iRow[212] = 26;
    jCol[212] = 22;
    iRow[213] = 26;
    jCol[213] = 23;
    iRow[214] = 26;
    jCol[214] = 24;
    iRow[215] = 26;
    jCol[215] = 25;
    iRow[216] = 26;
    jCol[216] = 26;
    iRow[217] = 26;
    jCol[217] = 27;
    iRow[218] = 26;
    jCol[218] = 40;
    iRow[219] = 28;
    jCol[219] = 24;
    iRow[220] = 28;
    jCol[220] = 25;
    iRow[221] = 28;
    jCol[221] = 26;
    iRow[222] = 28;
    jCol[222] = 27;
    iRow[223] = 28;
    jCol[223] = 28;
    iRow[224] = 28;
    jCol[224] = 29;
    iRow[225] = 28;
    jCol[225] = 30;
    iRow[226] = 28;
    jCol[226] = 31;
    iRow[227] = 28;
    jCol[227] = 40;
    iRow[228] = 29;
    jCol[228] = 24;
    iRow[229] = 29;
    jCol[229] = 25;
    iRow[230] = 29;
    jCol[230] = 26;
    iRow[231] = 29;
    jCol[231] = 27;
    iRow[232] = 29;
    jCol[232] = 28;
    iRow[233] = 29;
    jCol[233] = 29;
    iRow[234] = 29;
    jCol[234] = 30;
    iRow[235] = 29;
    jCol[235] = 31;
    iRow[236] = 29;
    jCol[236] = 40;
    iRow[237] = 30;
    jCol[237] = 24;
    iRow[238] = 30;
    jCol[238] = 25;
    iRow[239] = 30;
    jCol[239] = 26;
    iRow[240] = 30;
    jCol[240] = 27;
    iRow[241] = 30;
    jCol[241] = 28;
    iRow[242] = 30;
    jCol[242] = 29;
    iRow[243] = 30;
    jCol[243] = 30;
    iRow[244] = 30;
    jCol[244] = 31;
    iRow[245] = 30;
    jCol[245] = 40;
    iRow[246] = 32;
    jCol[246] = 28;
    iRow[247] = 32;
    jCol[247] = 29;
    iRow[248] = 32;
    jCol[248] = 30;
    iRow[249] = 32;
    jCol[249] = 31;
    iRow[250] = 32;
    jCol[250] = 32;
    iRow[251] = 32;
    jCol[251] = 33;
    iRow[252] = 32;
    jCol[252] = 34;
    iRow[253] = 32;
    jCol[253] = 35;
    iRow[254] = 32;
    jCol[254] = 40;
    iRow[255] = 33;
    jCol[255] = 28;
    iRow[256] = 33;
    jCol[256] = 29;
    iRow[257] = 33;
    jCol[257] = 30;
    iRow[258] = 33;
    jCol[258] = 31;
    iRow[259] = 33;
    jCol[259] = 32;
    iRow[260] = 33;
    jCol[260] = 33;
    iRow[261] = 33;
    jCol[261] = 34;
    iRow[262] = 33;
    jCol[262] = 35;
    iRow[263] = 33;
    jCol[263] = 40;
    iRow[264] = 34;
    jCol[264] = 28;
    iRow[265] = 34;
    jCol[265] = 29;
    iRow[266] = 34;
    jCol[266] = 30;
    iRow[267] = 34;
    jCol[267] = 31;
    iRow[268] = 34;
    jCol[268] = 32;
    iRow[269] = 34;
    jCol[269] = 33;
    iRow[270] = 34;
    jCol[270] = 34;
    iRow[271] = 34;
    jCol[271] = 35;
    iRow[272] = 34;
    jCol[272] = 40;
    iRow[273] = 36;
    jCol[273] = 32;
    iRow[274] = 36;
    jCol[274] = 33;
    iRow[275] = 36;
    jCol[275] = 34;
    iRow[276] = 36;
    jCol[276] = 35;
    iRow[277] = 36;
    jCol[277] = 36;
    iRow[278] = 36;
    jCol[278] = 37;
    iRow[279] = 36;
    jCol[279] = 38;
    iRow[280] = 36;
    jCol[280] = 39;
    iRow[281] = 36;
    jCol[281] = 40;
    iRow[282] = 37;
    jCol[282] = 32;
    iRow[283] = 37;
    jCol[283] = 33;
    iRow[284] = 37;
    jCol[284] = 34;
    iRow[285] = 37;
    jCol[285] = 35;
    iRow[286] = 37;
    jCol[286] = 36;
    iRow[287] = 37;
    jCol[287] = 37;
    iRow[288] = 37;
    jCol[288] = 38;
    iRow[289] = 37;
    jCol[289] = 39;
    iRow[290] = 37;
    jCol[290] = 40;
    iRow[291] = 38;
    jCol[291] = 32;
    iRow[292] = 38;
    jCol[292] = 33;
    iRow[293] = 38;
    jCol[293] = 34;
    iRow[294] = 38;
    jCol[294] = 35;
    iRow[295] = 38;
    jCol[295] = 36;
    iRow[296] = 38;
    jCol[296] = 37;
    iRow[297] = 38;
    jCol[297] = 38;
    iRow[298] = 38;
    jCol[298] = 39;
    iRow[299] = 38;
    jCol[299] = 40;

  }
  else {
    /* return the values of the jacobian of the constraints */

      Number l = my_data->l;


        Number t[10];
    Number grid[10];
    t[0] = 0;
    t[1] = 0.111111111111111*X[40];
    t[2] = 0.222222222222222*X[40];
    t[3] = 0.333333333333333*X[40];
    t[4] = 0.444444444444444*X[40];
    t[5] = 0.555555555555556*X[40];
    t[6] = 0.666666666666667*X[40];
    t[7] = 0.777777777777778*X[40];
    t[8] = 0.888888888888889*X[40];
    t[9] = 1.0*X[40];
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
    Number local[6];
    Number d_local[11];
    Number diffs[0];



    values[0] = -1;
    values[1] = -1;
    values[2] = -1;
    values[3] = 1;
    values[4] = 1;
    values[5] = 1;
    values[6] = 1;
    for (int iii = 0; iii < 10; iii++)
    {
      local[0] = X[iii * 4 + 0];
      local[1] = X[iii * 4 + 1];
      local[2] = X[iii * 4 + 2];
      local[3] = X[iii * 4 + 3];
      local[4] = X[40];
      local[5] = t[iii];
      values[7 + iii * 5] = -1;
      values[8 + iii * 5] = 0;
      values[9 + iii * 5] = 0;
      values[10 + iii * 5] = 0;
      values[11 + iii * 5] = 0;
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
      d_local[9] = grid[iii];
      d_local[10] = grid[iii + 1];
      values[57 + iii * 27] = -1;
      values[58 + iii * 27] = -1.0L/2.0L*d_local[10]*d_local[8] + (1.0L/2.0L)*d_local[8]*d_local[9];
      values[59 + iii * 27] = 0;
      values[60 + iii * 27] = -((1.0L/6.0L)*d_local[10]*d_local[8] - 1.0L/6.0L*d_local[8]*d_local[9])*((1.0L/2.0L)*d_local[10]*d_local[8] - 1.0L/2.0L*d_local[8]*d_local[9]);
      values[61 + iii * 27] = 1;
      values[62 + iii * 27] = -1.0L/2.0L*d_local[10]*d_local[8] + (1.0L/2.0L)*d_local[8]*d_local[9];
      values[63 + iii * 27] = 0;
      values[64 + iii * 27] = -(-1.0L/2.0L*d_local[10]*d_local[8] + (1.0L/2.0L)*d_local[8]*d_local[9])*((1.0L/6.0L)*d_local[10]*d_local[8] - 1.0L/6.0L*d_local[8]*d_local[9]);
      values[65 + iii * 27] = -4*((1.0L/8.0L)*d_local[10] - 1.0L/8.0L*d_local[9])*(d_local[3] - d_local[7])*((1.0L/6.0L)*d_local[10]*d_local[8] - 1.0L/6.0L*d_local[8]*d_local[9]) - ((1.0L/6.0L)*d_local[10] - 1.0L/6.0L*d_local[9])*(3*d_local[1] + 3*d_local[5] + 4*(d_local[3] - d_local[7])*((1.0L/8.0L)*d_local[10]*d_local[8] - 1.0L/8.0L*d_local[8]*d_local[9]));
      values[66 + iii * 27] = 0;
      values[67 + iii * 27] = -1;
      values[68 + iii * 27] = 0;
      values[69 + iii * 27] = -1.0L/2.0L*d_local[10]*d_local[8] + (1.0L/2.0L)*d_local[8]*d_local[9];
      values[70 + iii * 27] = 0;
      values[71 + iii * 27] = 1;
      values[72 + iii * 27] = 0;
      values[73 + iii * 27] = -1.0L/2.0L*d_local[10]*d_local[8] + (1.0L/2.0L)*d_local[8]*d_local[9];
      values[74 + iii * 27] = -((1.0L/6.0L)*d_local[10] - 1.0L/6.0L*d_local[9])*(3*d_local[3] + 3*d_local[7]);
      values[75 + iii * 27] = 0;
      values[76 + iii * 27] = 0;
      values[77 + iii * 27] = -1;
      values[78 + iii * 27] = -(2*d_local[3] + d_local[7])*((1.0L/6.0L)*d_local[10]*d_local[8] - 1.0L/6.0L*d_local[8]*d_local[9]);
      values[79 + iii * 27] = 0;
      values[80 + iii * 27] = 0;
      values[81 + iii * 27] = 1;
      values[82 + iii * 27] = -(d_local[3] + 2*d_local[7])*((1.0L/6.0L)*d_local[10]*d_local[8] - 1.0L/6.0L*d_local[8]*d_local[9]);
      values[83 + iii * 27] = -((1.0L/6.0L)*d_local[10] - 1.0L/6.0L*d_local[9])*((1.0L/2.0L)*pow(d_local[3], 2) + (1.0L/2.0L)*pow(d_local[7], 2) + 2*pow((1.0L/2.0L)*d_local[3] + (1.0L/2.0L)*d_local[7], 2));
    }


  }

  return TRUE;
}


Bool Double_Integrator__0_eval_h(Index n, Number *X, Bool new_x, Number obj_factor,
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


