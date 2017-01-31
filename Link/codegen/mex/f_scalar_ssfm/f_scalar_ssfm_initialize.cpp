/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * f_scalar_ssfm_initialize.cpp
 *
 * Code generation for function 'f_scalar_ssfm_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "f_scalar_ssfm.h"
#include "f_scalar_ssfm_initialize.h"
#include "_coder_f_scalar_ssfm_mex.h"
#include "f_scalar_ssfm_data.h"

/* Function Definitions */
void f_scalar_ssfm_initialize()
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (f_scalar_ssfm_initialize.cpp) */
