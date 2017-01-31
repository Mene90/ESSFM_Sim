/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * f_scalar_ssfm.h
 *
 * Code generation for function 'f_scalar_ssfm'
 *
 */

#ifndef F_SCALAR_SSFM_H
#define F_SCALAR_SSFM_H

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "omp.h"
#include "f_scalar_ssfm_types.h"

/* Function Declarations */
extern void b_r2br_r2dit_trig(const emlrtStack *sp, const emxArray_creal_T *x,
  int32_T n1_unsigned, const emxArray_real_T *costab, const emxArray_real_T
  *sintab, emxArray_creal_T *y);
extern void f_scalar_ssfm(const emlrtStack *sp, const Channel *ch, real_T Pavg,
  Signal *sig);
extern void get_algo_sizes(const emlrtStack *sp, int32_T n1, boolean_T useRadix2,
  int32_T *N2blue, int32_T *nRows);
extern void r2br_r2dit_trig(const emlrtStack *sp, const emxArray_creal_T *x,
  int32_T n1_unsigned, const emxArray_real_T *costab, const emxArray_real_T
  *sintab, emxArray_creal_T *y);
extern void r2br_r2dit_trig_impl(const emlrtStack *sp, const emxArray_creal_T *x,
  int32_T unsigned_nRows, const emxArray_real_T *costab, const emxArray_real_T
  *sintab, emxArray_creal_T *y);

#endif

/* End of code generation (f_scalar_ssfm.h) */
