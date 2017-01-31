/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * power.cpp
 *
 * Code generation for function 'power'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "f_scalar_ssfm.h"
#include "power.h"
#include "error.h"
#include "eml_int_forloop_overflow_check.h"
#include "f_scalar_ssfm_emxutil.h"
#include "scalexpAlloc.h"
#include "f_scalar_ssfm_data.h"

/* Variable Definitions */
static emlrtRTEInfo b_emlrtRTEI = { 1, /* lineNo */
  14,                                  /* colNo */
  "power",                             /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\ops\\power.m"/* pName */
};

/* Function Definitions */
boolean_T fltpower_domain_error(const emxArray_real_T *a, real_T b)
{
  boolean_T p;
  int32_T k;
  p = false;
  if (a->size[1] == 1) {
    if (a->data[0] < 0.0) {
      p = !(b == b);
    }
  } else {
    if (!(b == b)) {
      for (k = 0; k < a->size[1]; k++) {
        if (p || (a->data[k] < 0.0)) {
          p = true;
        } else {
          p = false;
        }
      }
    }
  }

  return p;
}

void power(const emlrtStack *sp, const emxArray_real_T *a, emxArray_real_T *y)
{
  emxArray_real_T *x;
  int32_T ub_loop;
  int32_T loop_ub;
  boolean_T overflow;
  int32_T k;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T(sp, &x, 2, &b_emlrtRTEI, true);
  st.site = &h_emlrtRSI;
  b_st.site = &i_emlrtRSI;
  ub_loop = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = a->size[1];
  emxEnsureCapacity(&b_st, (emxArray__common *)x, ub_loop, (int32_T)sizeof
                    (real_T), &b_emlrtRTEI);
  loop_ub = a->size[0] * a->size[1];
  for (ub_loop = 0; ub_loop < loop_ub; ub_loop++) {
    x->data[ub_loop] = a->data[ub_loop];
  }

  c_st.site = &k_emlrtRSI;
  ub_loop = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = a->size[1];
  emxEnsureCapacity(&c_st, (emxArray__common *)y, ub_loop, (int32_T)sizeof
                    (real_T), &c_emlrtRTEI);
  if (!dimagree(y, a)) {
    emlrtErrorWithMessageIdR2012b(&c_st, &ab_emlrtRTEI, "MATLAB:dimagree", 0);
  }

  ub_loop = a->size[1];
  c_st.site = &l_emlrtRSI;
  overflow = ((!(1 > a->size[1])) && (a->size[1] > 2147483646));
  if (overflow) {
    d_st.site = &m_emlrtRSI;
    check_forloop_overflow_error(&d_st);
  }

  emlrtEnterParallelRegion(&b_st, omp_in_parallel());

#pragma omp parallel for \
 num_threads(emlrtAllocRegionTLSs(b_st.tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs()))

  for (k = 1; k <= ub_loop; k++) {
    y->data[k - 1] = x->data[k - 1] * x->data[k - 1];
  }

  emlrtExitParallelRegion(&b_st, omp_in_parallel());
  emxFree_real_T(&x);
  if (fltpower_domain_error(a, 2.0)) {
    b_st.site = &j_emlrtRSI;
    error(&b_st);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (power.cpp) */
