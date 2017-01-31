/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * f_scalar_ssfm_terminate.cpp
 *
 * Code generation for function 'f_scalar_ssfm_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "f_scalar_ssfm.h"
#include "f_scalar_ssfm_terminate.h"
#include "_coder_f_scalar_ssfm_mex.h"
#include "f_scalar_ssfm_data.h"

/* Function Definitions */
void f_scalar_ssfm_atexit()
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void f_scalar_ssfm_terminate()
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (f_scalar_ssfm_terminate.cpp) */
