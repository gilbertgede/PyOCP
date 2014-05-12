
#ifndef PYOCP_COMMON_H
  #define PYOCP_COMMON_H
#endif

#ifndef IPSTDCINT
  #define IPSTDCINT
  #include "IpStdCInterface.h"
#endif

#ifndef FLOAT_H
  #define FLOAT_H
  #include <float.h>
#endif

#ifndef STDARG_H
  #define STDARG_H
  #include <stdarg.h>
#endif

#ifndef STDLIB_H
  #define STDLIB_H
  #include <stdlib.h>
#endif

#ifndef ASSERT_H
  #define ASSERT_H
  #include <assert.h>
#endif

#ifndef STDIO_H
  #define STDIO_H
  #include <stdio.h>
#endif

#ifndef MATH_H
  #define MATH_H
  #include <math.h>
#endif

#ifndef NUM_DIFF_DEF
  #define NUM_DIFF_DEF
  Number num_diff(Number (*f)(Number *x), Number *x0, int var, int len_x0);
#endif

