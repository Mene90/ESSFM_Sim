/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * f_scalar_ssfm_mexutil.cpp
 *
 * Code generation for function 'f_scalar_ssfm_mexutil'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "f_scalar_ssfm.h"
#include "f_scalar_ssfm_mexutil.h"
#include "f_scalar_ssfm_data.h"

/* Function Definitions */
emlrtCTX emlrtGetRootTLSGlobal()
{
  return emlrtRootTLSGlobal;
}

void emlrtLockerFunction(EmlrtLockeeFunction aLockee, const emlrtConstCTX aTLS,
  void *aData)
{
  omp_set_lock(&emlrtLockGlobal);
  emlrtCallLockeeFunction(aLockee, aTLS, aData);
  omp_unset_lock(&emlrtLockGlobal);
}

/* End of code generation (f_scalar_ssfm_mexutil.cpp) */
