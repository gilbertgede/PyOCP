IPOPT_c_template = \
"""\
#include "{base_filename}.h"

/* Main Program */
int main()
{{
  Index n=-1;                          /* number of variables */
  Index m=-1;                          /* number of constraints */
  Number* x_L = NULL;                  /* lower bounds on x */
  Number* x_U = NULL;                  /* upper bounds on x */
  Number* g_L = NULL;                  /* lower bounds on g */
  Number* g_U = NULL;                  /* upper bounds on g */
  IpoptProblem nlp = NULL;             /* IpoptProblem */
  enum ApplicationReturnStatus status; /* Solve return code */
  Number* x = NULL;                    /* starting point and solution vector */
  Number* mult_g = NULL;               /* constraint multipliers
             at the solution */
  Number* mult_x_L = NULL;             /* lower bound multipliers
             at the solution */
  Number* mult_x_U = NULL;             /* upper bound multipliers
             at the solution */
  Number obj;                          /* objective value */
  Index i;                             /* generic counter */

  /* Number of nonzeros in the Jacobian of the constraints */
  Index nele_jac = {num_nonzero_jac};
  /* Number of nonzeros in the Hessian of the Lagrangian (lower or
     upper triangual part only) */
  //Index nele_hess = {num_nonzero_hess}; // To be updated once hessian is
  //written
  Index nele_hess = 0;
  /* indexing style for matrices */
  Index index_style = 0; /* C-style; start counting of rows and column
             indices at 0 */

  /* our user data for the function evalutions. */
  struct MyUserData user_data;
  /* Initialize the user data */
  {arg_struct_fill}

  /* set the number of variables and allocate space for the bounds */
  n={num_nlp_vars};
  x_L = (Number*)malloc(sizeof(Number)*n);
  x_U = (Number*)malloc(sizeof(Number)*n);
  /* set the values for the variable bounds */
  {nlp_bounds_definition}

  /* set the number of constraints and allocate space for the bounds */
  m={num_nlp_constraints};
  g_L = (Number*)malloc(sizeof(Number)*m);
  g_U = (Number*)malloc(sizeof(Number)*m);
  /* set the values of the constraint bounds */
  {constraints_bound_definition}

  /* create the IpoptProblem */
  nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
                           index_style, &{base_filename}_eval_f,
                           &{base_filename}_eval_g,
                           &{base_filename}_eval_grad_f,
                           &{base_filename}_eval_jac_g,
                           &{base_filename}_eval_h);

  /* We can free the memory now - the values for the bounds have been
     copied internally in CreateIpoptProblem */
  free(x_L);
  free(x_U);
  free(g_L);
  free(g_U);

  /* Set some options.  Note the following ones are only examples,
     they might not be suitable for your problem. */
  AddIpoptNumOption(nlp, "tol", 1e-7);
  AddIpoptStrOption(nlp, "mu_strategy", "adaptive");
  AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");
  //AddIpoptStrOption(nlp, "output_file", "ipopt.out");

  /* allocate space for the initial point and set the values */
  x = (Number*)malloc(sizeof(Number)*n);
  {initial_nlp_values}

  /* allocate space to store the bound multipliers at the solution */
  mult_g = (Number*)malloc(sizeof(Number)*m);
  mult_x_L = (Number*)malloc(sizeof(Number)*n);
  mult_x_U = (Number*)malloc(sizeof(Number)*n);

  /* solve the problem */
  status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, &user_data);

  if (status == Solve_Succeeded) {{
    printf("\\n\\nSolution of the primal variables, x\\n");
    for (i=0; i<n; i++) {{
      printf("x[%d] = %e\\n", i, x[i]);
    }}

    printf("\\n\\nSolution of the ccnstraint multipliers, lambda\\n");
    for (i=0; i<m; i++) {{
      printf("lambda[%d] = %e\\n", i, mult_g[i]);
    }}
    printf("\\n\\nSolution of the bound multipliers, z_L and z_U\\n");
    for (i=0; i<n; i++) {{
      printf("z_L[%d] = %e\\n", i, mult_x_L[i]);
    }}
    for (i=0; i<n; i++) {{
      printf("z_U[%d] = %e\\n", i, mult_x_U[i]);
    }}

    printf("\\n\\nObjective value\\n");
    printf("f(x*) = %e\\n", obj);

    fs = fopen("{base_filename}_results.csv", "w");
    if(fs == NULL){{
        printf("Couldn't open an output file\n");
        return;
    }}
    for (int iii = 0; iii < n; iii++)
        fprintf(fs, "%g, ", x[iii]);
    fprintf(fs, "\\n");
    for (int iii = 0; iii < n; iii++)
        fprintf(fs, "%g, ", mult_x_L[iii]);
    fprintf(fs, "\\n");
    for (int iii = 0; iii < n; iii++)
        fprintf(fs, "%g, ", mult_x_U[iii]);
    fprintf(fs, "\\n");
    for (int iii = 0; iii < m; iii++)
        fprintf(fs, "%g, ", mult_g[iii]);
    fprintf(fs, "\\n");
    fclose(fs);

  }}
  else {{
    printf("\\n\\nERROR OCCURRED DURING IPOPT OPTIMIZATION.\\n");
  }}

  /* free allocated memory */
  FreeIpoptProblem(nlp);
  free(x);
  free(mult_g);
  free(mult_x_L);
  free(mult_x_U);

  return (int)status;
}}


/* Function Implementations */
Bool {base_filename}_eval_f(Index n, Number* X, Bool new_x,
            Number* obj_value, UserDataPtr user_data)
{{
  struct MyUserData* my_data = user_data;

  assert(n == {num_nlp_vars});

  {arg_struct_unpack}

  {objective_definition}

  return TRUE;
}}

Bool {base_filename}_eval_grad_f(Index n, Number* X, Bool new_x,
                 Number* grad_f, UserDataPtr user_data)
{{
  struct MyUserData* my_data = user_data;

  assert(n == {num_nlp_vars});

  {arg_struct_unpack}

  {grad_f_definition}

  return TRUE;
}}

Bool {base_filename}_eval_g(Index n, Number* X, Bool new_x,
            Index m, Number* g, UserDataPtr user_data)
{{
  struct MyUserData* my_data = user_data;

  assert({num_nlp_vars} == n);
  assert({num_nlp_constraints} == m);

  {arg_struct_unpack}

  {t_ccode}

  {g_definition}

  return TRUE;
}}

Bool {base_filename}_eval_jac_g(Index n, Number *X, Bool new_x,
                Index m, Index nele_jac,
                Index *iRow, Index *jCol, Number *values,
                UserDataPtr user_data)
{{
  struct MyUserData* my_data = user_data;
  if (values == NULL) {{
    /* return the structure of the jacobian */

    {jac_g_rows_cols}
  }}
  else {{
    /* return the values of the jacobian of the constraints */

    {arg_struct_unpack}

    {t_ccode}

    {jac_g_values_only}

  }}

  return TRUE;
}}


Bool {base_filename}_eval_h(Index n, Number *X, Bool new_x, Number obj_factor,
            Index m, Number *lambda, Bool new_lambda,
            Index nele_hess, Index *iRow, Index *jCol,
            Number *values, UserDataPtr user_data)
{{
  return 0;
}}

//
//Bool eval_h(Index n, Number *X, Bool new_x, Number obj_factor,
//            Index m, Number *lambda, Bool new_lambda,
//            Index nele_hess, Index *iRow, Index *jCol,
//            Number *values, UserDataPtr user_data)
//{{
//  Index idx = 0; /* nonzero element counter */
//  Index row = 0; /* row counter for loop */
//  Index col = 0; /* col counter for loop */
//  if (values == NULL) {{
//    /* return the structure. This is a symmetric matrix, fill the lower left
//     * triangle only. */
//
//    /* the hessian for this problem is actually dense */
//    idx=0;
//    for (row = 0; row < 4; row++) {{
//      for (col = 0; col <= row; col++) {{
//        iRow[idx] = row;
//        jCol[idx] = col;
//        idx++;
//      }}
//    }}
//
//    assert(idx == nele_hess);
//  }}
//  else {{
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
//  }}
//
//  return TRUE;
//}}


"""






IPOPT_h_template = \
"""\
#include "IpStdCInterface.h"
#include <float.h>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
{user_h_include}

/* Function Declarations */
Bool {base_filename}_eval_f(Index n, Number* x, Bool new_x,
            Number* obj_value, UserDataPtr user_data);

Bool {base_filename}_eval_grad_f(Index n, Number* x, Bool new_x,
                 Number* grad_f, UserDataPtr user_data);

Bool {base_filename}_eval_g(Index n, Number* x, Bool new_x,
            Index m, Number* g, UserDataPtr user_data);

Bool {base_filename}_eval_jac_g(Index n, Number *x, Bool new_x,
                Index m, Index nele_jac,
                Index *iRow, Index *jCol, Number *values,
                UserDataPtr user_data);

Bool {base_filename}_eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
            Index m, Number *lambda, Bool new_lambda,
            Index nele_hess, Index *iRow, Index *jCol,
            Number *values, UserDataPtr user_data);

/* This is an example how user_data can be used. */
{arg_struct_def}


/* Wrapping any user funcs which are multiarg to single vector arg. */
{wrap_user_funcs}


/* A numerical differentiation function */
Number num_diff(Number (*f)(Number *x), Number *x0, int var, int len_x0)
{{
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
    {{
        xll[var] = x0[var] - 2 * dx[var];
        xl[var] = x0[var] - dx[var];
        xu[var] = x0[var] + dx[var];
        xuu[var] = x0[var] + 2 * dx[var];
        h0 = dx[var];
        temp = 0.5  / pow(h0, 3) * (-f(xll) + 2 * f(xl) - 2 * f(xu) + f(xuu));
        if (temp < DBL_EPSILON)
        {{
            dx[var] = 1;
            break;
        }}
        dx[var] = pow(3 * DBL_EPSILON / temp, 1.0 / 3.0);
        if (abs(h0 - dx[var]) < DBL_EPSILON)
            break;
    }}
    h0 = dx[var];

    xu[var] = x0[var] + h0;
    xl[var] = x0[var] - h0;

    return 0.5 / h0 * (f(xu) - f(xl));
}};

"""


IPOPT_makefile_template = \
"""\
EXE = {base_filename}
OBJS = {base_filename}.o
ADDLIBS =
ADDINCFLAGS =

CC = cc
CFLAGS = -std=c99 -O3 -pipe -DNDEBUG -Wimplicit -Wparentheses -Wsequence-point -Wreturn-type -Wcast-qual -Wall -Wno-unknown-pragmas -Wno-long-long   -DIPOPT_BUILD
CLINKFLAGS =

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
INCL = `PKG_CONFIG_PATH=/opt/local/lib64/pkgconfig:/opt/local/lib/pkgconfig:/opt/local/share/pkgconfig: pkg-config --cflags ipopt` $(ADDINCFLAGS)
#INCL = -I`$(CYGPATH_W) /opt/local/include/coin`  $(ADDINCFLAGS)

# Linker flags
LIBS = `PKG_CONFIG_PATH=/opt/local/lib64/pkgconfig:/opt/local/lib/pkgconfig:/opt/local/share/pkgconfig: pkg-config --libs ipopt` -lstdc++ -lm
##LIBS = -link -libpath:`$(CYGPATH_W) /opt/local/lib` libipopt.lib -framework vecLib -framework vecLib -lm  -ldl -lstdc++ -lm
#LIBS = -L/opt/local/lib -lipopt -framework vecLib -framework vecLib -lm  -ldl -lstdc++ -lm

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

all: $(EXE)

.SUFFIXES: .c .o .obj

$(EXE): $(OBJS)
	bla=;\\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \\
	$(CC) $(CFLAGS) $(CLINKFLAGS) -o $@ $$bla $(LIBS) $(ADDLIBS)

clean:
	rm -f $(EXE) $(OBJS) ipopt.out

.c.o:
	$(CC) $(CFLAGS) $(INCL) -c -o $@ $<


.c.obj:
	$(CC) $(CFLAGS) $(INCL) -c -o $@ `$(CYGPATH_W) '$<'`

"""

