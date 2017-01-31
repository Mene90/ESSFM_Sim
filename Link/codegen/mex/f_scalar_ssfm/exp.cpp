/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * exp.cpp
 *
 * Code generation for function 'exp'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "f_scalar_ssfm.h"
#include "exp.h"
#include "eml_int_forloop_overflow_check.h"
#include "f_scalar_ssfm_data.h"

/* Function Definitions */
void b_exp(const emlrtStack *sp, emxArray_creal_T *x)
{
  int32_T nx;
  boolean_T overflow;
  int32_T k;
  real_T x_re;
  real_T r;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &n_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  nx = x->size[1];
  b_st.site = &o_emlrtRSI;
  overflow = ((!(1 > x->size[1])) && (x->size[1] > 2147483646));
  if (overflow) {
    c_st.site = &m_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }

  for (k = 0; k + 1 <= nx; k++) {
    if (x->data[k].im == 0.0) {
      x_re = muDoubleScalarExp(x->data[k].re);
      r = 0.0;
    } else if (muDoubleScalarIsInf(x->data[k].im) && muDoubleScalarIsInf(x->
                data[k].re) && (x->data[k].re < 0.0)) {
      x_re = 0.0;
      r = 0.0;
    } else {
      r = muDoubleScalarExp(x->data[k].re / 2.0);
      x_re = r * (r * muDoubleScalarCos(x->data[k].im));
      r *= r * muDoubleScalarSin(x->data[k].im);
    }

    x->data[k].re = x_re;
    x->data[k].im = r;
  }
}

/* End of code generation (exp.cpp) */
