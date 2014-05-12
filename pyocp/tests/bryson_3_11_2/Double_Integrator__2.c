
#include "Double_Integrator__2.h"






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
  Index nele_jac = 299;
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
  m=43;
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


  // create the IpoptProblem
  nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
                           index_style, &Double_Integrator__2_eval_f,
                           &Double_Integrator__2_eval_g,
                           &Double_Integrator__2_eval_grad_f,
                           &Double_Integrator__2_eval_jac_g,
                           &Double_Integrator__2_eval_h);

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
  x[4] = 0.0987654320987654;
  x[5] = -0.111111111111111;
  x[6] = 2.22222222222222;
  x[7] = -0.666666666666667;
  x[8] = 0.0864197530864197;
  x[9] = -0.222222222222222;
  x[10] = 2.44444444444444;
  x[11] = -1.33333333333333;
  x[12] = 0.0740740740740741;
  x[13] = -0.333333333333333;
  x[14] = 2.66666666666667;
  x[15] = -2.00000000000000;
  x[16] = 0.0617283950617284;
  x[17] = -0.444444444444444;
  x[18] = 2.88888888888889;
  x[19] = -2.66666666666667;
  x[20] = 0.0493827160493827;
  x[21] = -0.555555555555556;
  x[22] = 3.11111111111111;
  x[23] = -3.33333333333333;
  x[24] = 0.0370370370370370;
  x[25] = -0.666666666666667;
  x[26] = 3.33333333333333;
  x[27] = -4.00000000000000;
  x[28] = 0.0246913580246914;
  x[29] = -0.777777777777778;
  x[30] = 3.55555555555556;
  x[31] = -4.66666666666667;
  x[32] = 0.0123456790123457;
  x[33] = -0.888888888888889;
  x[34] = 3.77777777777778;
  x[35] = -5.33333333333333;
  x[36] = 0.0;
  x[37] = -1.00000000000000;
  x[38] = 4.00000000000000;
  x[39] = -6.00000000000000;
  x[40] = 0.700000000000000;


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
    fs = fopen("Double_Integrator__2_results.csv", "w");
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
Bool Double_Integrator__2_eval_f(Index n, Number* X, Bool new_x,
            Number* obj_value, UserDataPtr user_data)
{
  struct MyUserData* my_data = user_data;

  assert(n == 41);

    Number l = my_data->l;


  *obj_value = X[38];


  return TRUE;
}

Bool Double_Integrator__2_eval_grad_f(Index n, Number* X, Bool new_x,
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
  grad_f[39] = 0.0;
  grad_f[40] = 0.0;
  grad_f[38] = 1;


  return TRUE;
}

Bool Double_Integrator__2_eval_g(Index n, Number* X, Bool new_x,
            Index m, Number* g, UserDataPtr user_data)
{
  struct MyUserData* my_data = user_data;

  assert(41 == n);
  assert(43 == m);

    Number l = my_data->l;


      Number t[10];
    Number grid[10];
    t[0] = X[40];
    t[1] = 0.888888888888889*X[40] + 0.111111111111111;
    t[2] = 0.777777777777778*X[40] + 0.222222222222222;
    t[3] = 0.666666666666667*X[40] + 0.333333333333333;
    t[4] = 0.555555555555556*X[40] + 0.444444444444444;
    t[5] = 0.444444444444444*X[40] + 0.555555555555556;
    t[6] = 0.333333333333333*X[40] + 0.666666666666667;
    t[7] = 0.222222222222222*X[40] + 0.777777777777778;
    t[8] = 0.111111111111111*X[40] + 0.888888888888889;
    t[9] = 1.00000000000000;
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
  g[42] += X[37] + 1;
  Number local[6];
  for (int iii = 0; iii < 10; iii++)
  {
    local[0] = X[iii * 4 + 0];
    local[1] = X[iii * 4 + 1];
    local[2] = X[iii * 4 + 2];
    local[3] = X[iii * 4 + 3];
    local[4] = X[40];
    local[5] = t[iii];
    g[iii * 4 + 0 + 4] = -local[0] + l;
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
    g[iii * 4 + 0 + 5] = -d_local[0] + d_local[4] - ((1.0L/6.0L)*d_local[10]*(-d_local[8] + 1) - 1.0L/6.0L*d_local[9]*(-d_local[8] + 1))*(3*d_local[1] + 3*d_local[5] + 4*(d_local[3] - d_local[7])*((1.0L/8.0L)*d_local[10]*(-d_local[8] + 1) - 1.0L/8.0L*d_local[9]*(-d_local[8] + 1)));
    g[iii * 4 + 1 + 5] = -d_local[1] + d_local[5] - (3*d_local[3] + 3*d_local[7])*((1.0L/6.0L)*d_local[10]*(-d_local[8] + 1) - 1.0L/6.0L*d_local[9]*(-d_local[8] + 1));
    g[iii * 4 + 2 + 5] = -d_local[2] + d_local[6] - ((1.0L/6.0L)*d_local[10]*(-d_local[8] + 1) - 1.0L/6.0L*d_local[9]*(-d_local[8] + 1))*((1.0L/2.0L)*pow(d_local[3], 2) + (1.0L/2.0L)*pow(d_local[7], 2) + 2*pow((1.0L/2.0L)*d_local[3] + (1.0L/2.0L)*d_local[7], 2));
  }


  return TRUE;
}

Bool Double_Integrator__2_eval_jac_g(Index n, Number *X, Bool new_x,
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
    iRow[6] = 4;
    jCol[6] = 0;
    iRow[7] = 4;
    jCol[7] = 1;
    iRow[8] = 4;
    jCol[8] = 2;
    iRow[9] = 4;
    jCol[9] = 3;
    iRow[10] = 4;
    jCol[10] = 40;
    iRow[11] = 8;
    jCol[11] = 4;
    iRow[12] = 8;
    jCol[12] = 5;
    iRow[13] = 8;
    jCol[13] = 6;
    iRow[14] = 8;
    jCol[14] = 7;
    iRow[15] = 8;
    jCol[15] = 40;
    iRow[16] = 12;
    jCol[16] = 8;
    iRow[17] = 12;
    jCol[17] = 9;
    iRow[18] = 12;
    jCol[18] = 10;
    iRow[19] = 12;
    jCol[19] = 11;
    iRow[20] = 12;
    jCol[20] = 40;
    iRow[21] = 16;
    jCol[21] = 12;
    iRow[22] = 16;
    jCol[22] = 13;
    iRow[23] = 16;
    jCol[23] = 14;
    iRow[24] = 16;
    jCol[24] = 15;
    iRow[25] = 16;
    jCol[25] = 40;
    iRow[26] = 20;
    jCol[26] = 16;
    iRow[27] = 20;
    jCol[27] = 17;
    iRow[28] = 20;
    jCol[28] = 18;
    iRow[29] = 20;
    jCol[29] = 19;
    iRow[30] = 20;
    jCol[30] = 40;
    iRow[31] = 24;
    jCol[31] = 20;
    iRow[32] = 24;
    jCol[32] = 21;
    iRow[33] = 24;
    jCol[33] = 22;
    iRow[34] = 24;
    jCol[34] = 23;
    iRow[35] = 24;
    jCol[35] = 40;
    iRow[36] = 28;
    jCol[36] = 24;
    iRow[37] = 28;
    jCol[37] = 25;
    iRow[38] = 28;
    jCol[38] = 26;
    iRow[39] = 28;
    jCol[39] = 27;
    iRow[40] = 28;
    jCol[40] = 40;
    iRow[41] = 32;
    jCol[41] = 28;
    iRow[42] = 32;
    jCol[42] = 29;
    iRow[43] = 32;
    jCol[43] = 30;
    iRow[44] = 32;
    jCol[44] = 31;
    iRow[45] = 32;
    jCol[45] = 40;
    iRow[46] = 36;
    jCol[46] = 32;
    iRow[47] = 36;
    jCol[47] = 33;
    iRow[48] = 36;
    jCol[48] = 34;
    iRow[49] = 36;
    jCol[49] = 35;
    iRow[50] = 36;
    jCol[50] = 40;
    iRow[51] = 40;
    jCol[51] = 36;
    iRow[52] = 40;
    jCol[52] = 37;
    iRow[53] = 40;
    jCol[53] = 38;
    iRow[54] = 40;
    jCol[54] = 39;
    iRow[55] = 40;
    jCol[55] = 40;
    iRow[56] = 5;
    jCol[56] = 0;
    iRow[57] = 5;
    jCol[57] = 1;
    iRow[58] = 5;
    jCol[58] = 2;
    iRow[59] = 5;
    jCol[59] = 3;
    iRow[60] = 5;
    jCol[60] = 4;
    iRow[61] = 5;
    jCol[61] = 5;
    iRow[62] = 5;
    jCol[62] = 6;
    iRow[63] = 5;
    jCol[63] = 7;
    iRow[64] = 5;
    jCol[64] = 40;
    iRow[65] = 6;
    jCol[65] = 0;
    iRow[66] = 6;
    jCol[66] = 1;
    iRow[67] = 6;
    jCol[67] = 2;
    iRow[68] = 6;
    jCol[68] = 3;
    iRow[69] = 6;
    jCol[69] = 4;
    iRow[70] = 6;
    jCol[70] = 5;
    iRow[71] = 6;
    jCol[71] = 6;
    iRow[72] = 6;
    jCol[72] = 7;
    iRow[73] = 6;
    jCol[73] = 40;
    iRow[74] = 7;
    jCol[74] = 0;
    iRow[75] = 7;
    jCol[75] = 1;
    iRow[76] = 7;
    jCol[76] = 2;
    iRow[77] = 7;
    jCol[77] = 3;
    iRow[78] = 7;
    jCol[78] = 4;
    iRow[79] = 7;
    jCol[79] = 5;
    iRow[80] = 7;
    jCol[80] = 6;
    iRow[81] = 7;
    jCol[81] = 7;
    iRow[82] = 7;
    jCol[82] = 40;
    iRow[83] = 9;
    jCol[83] = 4;
    iRow[84] = 9;
    jCol[84] = 5;
    iRow[85] = 9;
    jCol[85] = 6;
    iRow[86] = 9;
    jCol[86] = 7;
    iRow[87] = 9;
    jCol[87] = 8;
    iRow[88] = 9;
    jCol[88] = 9;
    iRow[89] = 9;
    jCol[89] = 10;
    iRow[90] = 9;
    jCol[90] = 11;
    iRow[91] = 9;
    jCol[91] = 40;
    iRow[92] = 10;
    jCol[92] = 4;
    iRow[93] = 10;
    jCol[93] = 5;
    iRow[94] = 10;
    jCol[94] = 6;
    iRow[95] = 10;
    jCol[95] = 7;
    iRow[96] = 10;
    jCol[96] = 8;
    iRow[97] = 10;
    jCol[97] = 9;
    iRow[98] = 10;
    jCol[98] = 10;
    iRow[99] = 10;
    jCol[99] = 11;
    iRow[100] = 10;
    jCol[100] = 40;
    iRow[101] = 11;
    jCol[101] = 4;
    iRow[102] = 11;
    jCol[102] = 5;
    iRow[103] = 11;
    jCol[103] = 6;
    iRow[104] = 11;
    jCol[104] = 7;
    iRow[105] = 11;
    jCol[105] = 8;
    iRow[106] = 11;
    jCol[106] = 9;
    iRow[107] = 11;
    jCol[107] = 10;
    iRow[108] = 11;
    jCol[108] = 11;
    iRow[109] = 11;
    jCol[109] = 40;
    iRow[110] = 13;
    jCol[110] = 8;
    iRow[111] = 13;
    jCol[111] = 9;
    iRow[112] = 13;
    jCol[112] = 10;
    iRow[113] = 13;
    jCol[113] = 11;
    iRow[114] = 13;
    jCol[114] = 12;
    iRow[115] = 13;
    jCol[115] = 13;
    iRow[116] = 13;
    jCol[116] = 14;
    iRow[117] = 13;
    jCol[117] = 15;
    iRow[118] = 13;
    jCol[118] = 40;
    iRow[119] = 14;
    jCol[119] = 8;
    iRow[120] = 14;
    jCol[120] = 9;
    iRow[121] = 14;
    jCol[121] = 10;
    iRow[122] = 14;
    jCol[122] = 11;
    iRow[123] = 14;
    jCol[123] = 12;
    iRow[124] = 14;
    jCol[124] = 13;
    iRow[125] = 14;
    jCol[125] = 14;
    iRow[126] = 14;
    jCol[126] = 15;
    iRow[127] = 14;
    jCol[127] = 40;
    iRow[128] = 15;
    jCol[128] = 8;
    iRow[129] = 15;
    jCol[129] = 9;
    iRow[130] = 15;
    jCol[130] = 10;
    iRow[131] = 15;
    jCol[131] = 11;
    iRow[132] = 15;
    jCol[132] = 12;
    iRow[133] = 15;
    jCol[133] = 13;
    iRow[134] = 15;
    jCol[134] = 14;
    iRow[135] = 15;
    jCol[135] = 15;
    iRow[136] = 15;
    jCol[136] = 40;
    iRow[137] = 17;
    jCol[137] = 12;
    iRow[138] = 17;
    jCol[138] = 13;
    iRow[139] = 17;
    jCol[139] = 14;
    iRow[140] = 17;
    jCol[140] = 15;
    iRow[141] = 17;
    jCol[141] = 16;
    iRow[142] = 17;
    jCol[142] = 17;
    iRow[143] = 17;
    jCol[143] = 18;
    iRow[144] = 17;
    jCol[144] = 19;
    iRow[145] = 17;
    jCol[145] = 40;
    iRow[146] = 18;
    jCol[146] = 12;
    iRow[147] = 18;
    jCol[147] = 13;
    iRow[148] = 18;
    jCol[148] = 14;
    iRow[149] = 18;
    jCol[149] = 15;
    iRow[150] = 18;
    jCol[150] = 16;
    iRow[151] = 18;
    jCol[151] = 17;
    iRow[152] = 18;
    jCol[152] = 18;
    iRow[153] = 18;
    jCol[153] = 19;
    iRow[154] = 18;
    jCol[154] = 40;
    iRow[155] = 19;
    jCol[155] = 12;
    iRow[156] = 19;
    jCol[156] = 13;
    iRow[157] = 19;
    jCol[157] = 14;
    iRow[158] = 19;
    jCol[158] = 15;
    iRow[159] = 19;
    jCol[159] = 16;
    iRow[160] = 19;
    jCol[160] = 17;
    iRow[161] = 19;
    jCol[161] = 18;
    iRow[162] = 19;
    jCol[162] = 19;
    iRow[163] = 19;
    jCol[163] = 40;
    iRow[164] = 21;
    jCol[164] = 16;
    iRow[165] = 21;
    jCol[165] = 17;
    iRow[166] = 21;
    jCol[166] = 18;
    iRow[167] = 21;
    jCol[167] = 19;
    iRow[168] = 21;
    jCol[168] = 20;
    iRow[169] = 21;
    jCol[169] = 21;
    iRow[170] = 21;
    jCol[170] = 22;
    iRow[171] = 21;
    jCol[171] = 23;
    iRow[172] = 21;
    jCol[172] = 40;
    iRow[173] = 22;
    jCol[173] = 16;
    iRow[174] = 22;
    jCol[174] = 17;
    iRow[175] = 22;
    jCol[175] = 18;
    iRow[176] = 22;
    jCol[176] = 19;
    iRow[177] = 22;
    jCol[177] = 20;
    iRow[178] = 22;
    jCol[178] = 21;
    iRow[179] = 22;
    jCol[179] = 22;
    iRow[180] = 22;
    jCol[180] = 23;
    iRow[181] = 22;
    jCol[181] = 40;
    iRow[182] = 23;
    jCol[182] = 16;
    iRow[183] = 23;
    jCol[183] = 17;
    iRow[184] = 23;
    jCol[184] = 18;
    iRow[185] = 23;
    jCol[185] = 19;
    iRow[186] = 23;
    jCol[186] = 20;
    iRow[187] = 23;
    jCol[187] = 21;
    iRow[188] = 23;
    jCol[188] = 22;
    iRow[189] = 23;
    jCol[189] = 23;
    iRow[190] = 23;
    jCol[190] = 40;
    iRow[191] = 25;
    jCol[191] = 20;
    iRow[192] = 25;
    jCol[192] = 21;
    iRow[193] = 25;
    jCol[193] = 22;
    iRow[194] = 25;
    jCol[194] = 23;
    iRow[195] = 25;
    jCol[195] = 24;
    iRow[196] = 25;
    jCol[196] = 25;
    iRow[197] = 25;
    jCol[197] = 26;
    iRow[198] = 25;
    jCol[198] = 27;
    iRow[199] = 25;
    jCol[199] = 40;
    iRow[200] = 26;
    jCol[200] = 20;
    iRow[201] = 26;
    jCol[201] = 21;
    iRow[202] = 26;
    jCol[202] = 22;
    iRow[203] = 26;
    jCol[203] = 23;
    iRow[204] = 26;
    jCol[204] = 24;
    iRow[205] = 26;
    jCol[205] = 25;
    iRow[206] = 26;
    jCol[206] = 26;
    iRow[207] = 26;
    jCol[207] = 27;
    iRow[208] = 26;
    jCol[208] = 40;
    iRow[209] = 27;
    jCol[209] = 20;
    iRow[210] = 27;
    jCol[210] = 21;
    iRow[211] = 27;
    jCol[211] = 22;
    iRow[212] = 27;
    jCol[212] = 23;
    iRow[213] = 27;
    jCol[213] = 24;
    iRow[214] = 27;
    jCol[214] = 25;
    iRow[215] = 27;
    jCol[215] = 26;
    iRow[216] = 27;
    jCol[216] = 27;
    iRow[217] = 27;
    jCol[217] = 40;
    iRow[218] = 29;
    jCol[218] = 24;
    iRow[219] = 29;
    jCol[219] = 25;
    iRow[220] = 29;
    jCol[220] = 26;
    iRow[221] = 29;
    jCol[221] = 27;
    iRow[222] = 29;
    jCol[222] = 28;
    iRow[223] = 29;
    jCol[223] = 29;
    iRow[224] = 29;
    jCol[224] = 30;
    iRow[225] = 29;
    jCol[225] = 31;
    iRow[226] = 29;
    jCol[226] = 40;
    iRow[227] = 30;
    jCol[227] = 24;
    iRow[228] = 30;
    jCol[228] = 25;
    iRow[229] = 30;
    jCol[229] = 26;
    iRow[230] = 30;
    jCol[230] = 27;
    iRow[231] = 30;
    jCol[231] = 28;
    iRow[232] = 30;
    jCol[232] = 29;
    iRow[233] = 30;
    jCol[233] = 30;
    iRow[234] = 30;
    jCol[234] = 31;
    iRow[235] = 30;
    jCol[235] = 40;
    iRow[236] = 31;
    jCol[236] = 24;
    iRow[237] = 31;
    jCol[237] = 25;
    iRow[238] = 31;
    jCol[238] = 26;
    iRow[239] = 31;
    jCol[239] = 27;
    iRow[240] = 31;
    jCol[240] = 28;
    iRow[241] = 31;
    jCol[241] = 29;
    iRow[242] = 31;
    jCol[242] = 30;
    iRow[243] = 31;
    jCol[243] = 31;
    iRow[244] = 31;
    jCol[244] = 40;
    iRow[245] = 33;
    jCol[245] = 28;
    iRow[246] = 33;
    jCol[246] = 29;
    iRow[247] = 33;
    jCol[247] = 30;
    iRow[248] = 33;
    jCol[248] = 31;
    iRow[249] = 33;
    jCol[249] = 32;
    iRow[250] = 33;
    jCol[250] = 33;
    iRow[251] = 33;
    jCol[251] = 34;
    iRow[252] = 33;
    jCol[252] = 35;
    iRow[253] = 33;
    jCol[253] = 40;
    iRow[254] = 34;
    jCol[254] = 28;
    iRow[255] = 34;
    jCol[255] = 29;
    iRow[256] = 34;
    jCol[256] = 30;
    iRow[257] = 34;
    jCol[257] = 31;
    iRow[258] = 34;
    jCol[258] = 32;
    iRow[259] = 34;
    jCol[259] = 33;
    iRow[260] = 34;
    jCol[260] = 34;
    iRow[261] = 34;
    jCol[261] = 35;
    iRow[262] = 34;
    jCol[262] = 40;
    iRow[263] = 35;
    jCol[263] = 28;
    iRow[264] = 35;
    jCol[264] = 29;
    iRow[265] = 35;
    jCol[265] = 30;
    iRow[266] = 35;
    jCol[266] = 31;
    iRow[267] = 35;
    jCol[267] = 32;
    iRow[268] = 35;
    jCol[268] = 33;
    iRow[269] = 35;
    jCol[269] = 34;
    iRow[270] = 35;
    jCol[270] = 35;
    iRow[271] = 35;
    jCol[271] = 40;
    iRow[272] = 37;
    jCol[272] = 32;
    iRow[273] = 37;
    jCol[273] = 33;
    iRow[274] = 37;
    jCol[274] = 34;
    iRow[275] = 37;
    jCol[275] = 35;
    iRow[276] = 37;
    jCol[276] = 36;
    iRow[277] = 37;
    jCol[277] = 37;
    iRow[278] = 37;
    jCol[278] = 38;
    iRow[279] = 37;
    jCol[279] = 39;
    iRow[280] = 37;
    jCol[280] = 40;
    iRow[281] = 38;
    jCol[281] = 32;
    iRow[282] = 38;
    jCol[282] = 33;
    iRow[283] = 38;
    jCol[283] = 34;
    iRow[284] = 38;
    jCol[284] = 35;
    iRow[285] = 38;
    jCol[285] = 36;
    iRow[286] = 38;
    jCol[286] = 37;
    iRow[287] = 38;
    jCol[287] = 38;
    iRow[288] = 38;
    jCol[288] = 39;
    iRow[289] = 38;
    jCol[289] = 40;
    iRow[290] = 39;
    jCol[290] = 32;
    iRow[291] = 39;
    jCol[291] = 33;
    iRow[292] = 39;
    jCol[292] = 34;
    iRow[293] = 39;
    jCol[293] = 35;
    iRow[294] = 39;
    jCol[294] = 36;
    iRow[295] = 39;
    jCol[295] = 37;
    iRow[296] = 39;
    jCol[296] = 38;
    iRow[297] = 39;
    jCol[297] = 39;
    iRow[298] = 39;
    jCol[298] = 40;

  }
  else {
    /* return the values of the jacobian of the constraints */

      Number l = my_data->l;


        Number t[10];
    Number grid[10];
    t[0] = X[40];
    t[1] = 0.888888888888889*X[40] + 0.111111111111111;
    t[2] = 0.777777777777778*X[40] + 0.222222222222222;
    t[3] = 0.666666666666667*X[40] + 0.333333333333333;
    t[4] = 0.555555555555556*X[40] + 0.444444444444444;
    t[5] = 0.444444444444444*X[40] + 0.555555555555556;
    t[6] = 0.333333333333333*X[40] + 0.666666666666667;
    t[7] = 0.222222222222222*X[40] + 0.777777777777778;
    t[8] = 0.111111111111111*X[40] + 0.888888888888889;
    t[9] = 1.00000000000000;
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
    values[3] = -1;
    values[4] = 1;
    values[5] = 1;
    for (int iii = 0; iii < 10; iii++)
    {
      local[0] = X[iii * 4 + 0];
      local[1] = X[iii * 4 + 1];
      local[2] = X[iii * 4 + 2];
      local[3] = X[iii * 4 + 3];
      local[4] = X[40];
      local[5] = t[iii];
      values[6 + iii * 5] = -1;
      values[7 + iii * 5] = 0;
      values[8 + iii * 5] = 0;
      values[9 + iii * 5] = 0;
      values[10 + iii * 5] = 0;
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
      values[56 + iii * 27] = -1;
      values[57 + iii * 27] = -1.0L/2.0L*d_local[10]*(-d_local[8] + 1) + (1.0L/2.0L)*d_local[9]*(-d_local[8] + 1);
      values[58 + iii * 27] = 0;
      values[59 + iii * 27] = -((1.0L/6.0L)*d_local[10]*(-d_local[8] + 1) - 1.0L/6.0L*d_local[9]*(-d_local[8] + 1))*((1.0L/2.0L)*d_local[10]*(-d_local[8] + 1) - 1.0L/2.0L*d_local[9]*(-d_local[8] + 1));
      values[60 + iii * 27] = 1;
      values[61 + iii * 27] = -1.0L/2.0L*d_local[10]*(-d_local[8] + 1) + (1.0L/2.0L)*d_local[9]*(-d_local[8] + 1);
      values[62 + iii * 27] = 0;
      values[63 + iii * 27] = -(-1.0L/2.0L*d_local[10]*(-d_local[8] + 1) + (1.0L/2.0L)*d_local[9]*(-d_local[8] + 1))*((1.0L/6.0L)*d_local[10]*(-d_local[8] + 1) - 1.0L/6.0L*d_local[9]*(-d_local[8] + 1));
      values[64 + iii * 27] = -(-1.0L/6.0L*d_local[10] + (1.0L/6.0L)*d_local[9])*(3*d_local[1] + 3*d_local[5] + 4*(d_local[3] - d_local[7])*((1.0L/8.0L)*d_local[10]*(-d_local[8] + 1) - 1.0L/8.0L*d_local[9]*(-d_local[8] + 1))) - 4*(-1.0L/8.0L*d_local[10] + (1.0L/8.0L)*d_local[9])*(d_local[3] - d_local[7])*((1.0L/6.0L)*d_local[10]*(-d_local[8] + 1) - 1.0L/6.0L*d_local[9]*(-d_local[8] + 1));
      values[65 + iii * 27] = 0;
      values[66 + iii * 27] = -1;
      values[67 + iii * 27] = 0;
      values[68 + iii * 27] = -1.0L/2.0L*d_local[10]*(-d_local[8] + 1) + (1.0L/2.0L)*d_local[9]*(-d_local[8] + 1);
      values[69 + iii * 27] = 0;
      values[70 + iii * 27] = 1;
      values[71 + iii * 27] = 0;
      values[72 + iii * 27] = -1.0L/2.0L*d_local[10]*(-d_local[8] + 1) + (1.0L/2.0L)*d_local[9]*(-d_local[8] + 1);
      values[73 + iii * 27] = -(-1.0L/6.0L*d_local[10] + (1.0L/6.0L)*d_local[9])*(3*d_local[3] + 3*d_local[7]);
      values[74 + iii * 27] = 0;
      values[75 + iii * 27] = 0;
      values[76 + iii * 27] = -1;
      values[77 + iii * 27] = -(2*d_local[3] + d_local[7])*((1.0L/6.0L)*d_local[10]*(-d_local[8] + 1) - 1.0L/6.0L*d_local[9]*(-d_local[8] + 1));
      values[78 + iii * 27] = 0;
      values[79 + iii * 27] = 0;
      values[80 + iii * 27] = 1;
      values[81 + iii * 27] = -(d_local[3] + 2*d_local[7])*((1.0L/6.0L)*d_local[10]*(-d_local[8] + 1) - 1.0L/6.0L*d_local[9]*(-d_local[8] + 1));
      values[82 + iii * 27] = -(-1.0L/6.0L*d_local[10] + (1.0L/6.0L)*d_local[9])*((1.0L/2.0L)*pow(d_local[3], 2) + (1.0L/2.0L)*pow(d_local[7], 2) + 2*pow((1.0L/2.0L)*d_local[3] + (1.0L/2.0L)*d_local[7], 2));
    }


  }

  return TRUE;
}


Bool Double_Integrator__2_eval_h(Index n, Number *X, Bool new_x, Number obj_factor,
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


