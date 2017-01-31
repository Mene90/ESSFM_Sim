/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fft.cpp
 *
 * Code generation for function 'fft'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "f_scalar_ssfm.h"
#include "fft.h"
#include "f_scalar_ssfm_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "f_scalar_ssfm_data.h"

/* Variable Definitions */
static emlrtRSInfo r_emlrtRSI = { 14,  /* lineNo */
  "fft",                               /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\fft.m"/* pathName */
};

static emlrtRSInfo ub_emlrtRSI = { 251,/* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

static emlrtRTEInfo f_emlrtRTEI = { 1, /* lineNo */
  14,                                  /* colNo */
  "fft",                               /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\fft.m"/* pName */
};

/* Function Declarations */
static void generate_twiddle_tables(const emlrtStack *sp, int32_T nRows,
  boolean_T useRadix2, emxArray_real_T *costab, emxArray_real_T *sintab,
  emxArray_real_T *sintabinv);

/* Function Definitions */
static void generate_twiddle_tables(const emlrtStack *sp, int32_T nRows,
  boolean_T useRadix2, emxArray_real_T *costab, emxArray_real_T *sintab,
  emxArray_real_T *sintabinv)
{
  emxArray_real_T *costab1q;
  real_T e;
  int32_T nRowsD4;
  int32_T nd2;
  int32_T k;
  int32_T n2;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T(sp, &costab1q, 2, &m_emlrtRTEI, true);
  e = 6.2831853071795862 / (real_T)nRows;
  nRowsD4 = nRows / 2 / 2;
  st.site = &x_emlrtRSI;
  nd2 = costab1q->size[0] * costab1q->size[1];
  costab1q->size[0] = 1;
  costab1q->size[1] = nRowsD4 + 1;
  emxEnsureCapacity(&st, (emxArray__common *)costab1q, nd2, (int32_T)sizeof
                    (real_T), &l_emlrtRTEI);
  costab1q->data[0] = 1.0;
  nd2 = nRowsD4 / 2;
  for (k = 1; k <= nd2; k++) {
    costab1q->data[k] = muDoubleScalarCos(e * (real_T)k);
  }

  for (k = nd2 + 1; k < nRowsD4; k++) {
    costab1q->data[k] = muDoubleScalarSin(e * (real_T)(nRowsD4 - k));
  }

  costab1q->data[nRowsD4] = 0.0;
  if (!useRadix2) {
    st.site = &y_emlrtRSI;
    nRowsD4 = costab1q->size[1] - 1;
    n2 = (costab1q->size[1] - 1) << 1;
    nd2 = costab->size[0] * costab->size[1];
    costab->size[0] = 1;
    costab->size[1] = n2 + 1;
    emxEnsureCapacity(&st, (emxArray__common *)costab, nd2, (int32_T)sizeof
                      (real_T), &l_emlrtRTEI);
    nd2 = sintab->size[0] * sintab->size[1];
    sintab->size[0] = 1;
    sintab->size[1] = n2 + 1;
    emxEnsureCapacity(&st, (emxArray__common *)sintab, nd2, (int32_T)sizeof
                      (real_T), &l_emlrtRTEI);
    costab->data[0] = 1.0;
    sintab->data[0] = 0.0;
    nd2 = sintabinv->size[0] * sintabinv->size[1];
    sintabinv->size[0] = 1;
    sintabinv->size[1] = n2 + 1;
    emxEnsureCapacity(&st, (emxArray__common *)sintabinv, nd2, (int32_T)sizeof
                      (real_T), &l_emlrtRTEI);
    for (k = 1; k <= nRowsD4; k++) {
      sintabinv->data[k] = costab1q->data[nRowsD4 - k];
    }

    for (k = costab1q->size[1]; k <= n2; k++) {
      sintabinv->data[k] = costab1q->data[k - nRowsD4];
    }

    for (k = 1; k <= nRowsD4; k++) {
      costab->data[k] = costab1q->data[k];
      sintab->data[k] = -costab1q->data[nRowsD4 - k];
    }

    for (k = costab1q->size[1]; k <= n2; k++) {
      costab->data[k] = -costab1q->data[n2 - k];
      sintab->data[k] = -costab1q->data[k - nRowsD4];
    }
  } else {
    st.site = &ab_emlrtRSI;
    nRowsD4 = costab1q->size[1] - 1;
    n2 = (costab1q->size[1] - 1) << 1;
    nd2 = costab->size[0] * costab->size[1];
    costab->size[0] = 1;
    costab->size[1] = n2 + 1;
    emxEnsureCapacity(&st, (emxArray__common *)costab, nd2, (int32_T)sizeof
                      (real_T), &l_emlrtRTEI);
    nd2 = sintab->size[0] * sintab->size[1];
    sintab->size[0] = 1;
    sintab->size[1] = n2 + 1;
    emxEnsureCapacity(&st, (emxArray__common *)sintab, nd2, (int32_T)sizeof
                      (real_T), &l_emlrtRTEI);
    costab->data[0] = 1.0;
    sintab->data[0] = 0.0;
    for (k = 1; k <= nRowsD4; k++) {
      costab->data[k] = costab1q->data[k];
      sintab->data[k] = -costab1q->data[nRowsD4 - k];
    }

    for (k = costab1q->size[1]; k <= n2; k++) {
      costab->data[k] = -costab1q->data[n2 - k];
      sintab->data[k] = -costab1q->data[k - nRowsD4];
    }

    nd2 = sintabinv->size[0] * sintabinv->size[1];
    sintabinv->size[0] = 1;
    sintabinv->size[1] = 0;
    emxEnsureCapacity(sp, (emxArray__common *)sintabinv, nd2, (int32_T)sizeof
                      (real_T), &l_emlrtRTEI);
  }

  emxFree_real_T(&costab1q);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void fft(const emlrtStack *sp, const emxArray_creal_T *x, emxArray_creal_T *y)
{
  boolean_T useRadix2;
  int32_T n1;
  int32_T idx;
  emxArray_creal_T *b_y1;
  emxArray_real_T *costab;
  emxArray_real_T *sintab;
  emxArray_real_T *sintabinv;
  emxArray_creal_T *wwc;
  emxArray_creal_T *fy;
  emxArray_creal_T *fv;
  int32_T N2blue;
  int32_T rt;
  int32_T b_x[1];
  int32_T nInt2m1;
  emxArray_creal_T c_x;
  int32_T nInt2;
  int32_T k;
  int32_T b_y;
  real_T nt_im;
  real_T nt_re;
  real_T wwc_re;
  real_T fv_im;
  real_T wwc_im;
  real_T fv_re;
  real_T b_fv_re;
  real_T b_fv_im;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  st.site = &r_emlrtRSI;
  if ((x->size[1] == 1) || (x->size[1] != 1)) {
    useRadix2 = true;
  } else {
    useRadix2 = false;
  }

  if (!useRadix2) {
    emlrtErrorWithMessageIdR2012b(&st, &db_emlrtRTEI,
      "Coder:toolbox:autoDimIncompatibility", 0);
  }

  n1 = x->size[1];
  b_st.site = &s_emlrtRSI;
  idx = x->size[1];
  emxInit_creal_T(&b_st, &b_y1, 1, &h_emlrtRTEI, true);
  emxInit_real_T(&b_st, &costab, 2, &f_emlrtRTEI, true);
  emxInit_real_T(&b_st, &sintab, 2, &f_emlrtRTEI, true);
  emxInit_real_T(&b_st, &sintabinv, 2, &f_emlrtRTEI, true);
  emxInit_creal_T(&b_st, &wwc, 1, &i_emlrtRTEI, true);
  emxInit_creal_T(&b_st, &fy, 1, &j_emlrtRTEI, true);
  emxInit_creal_T(&b_st, &fv, 1, &k_emlrtRTEI, true);
  if ((idx == 0) || (x->size[1] == 0)) {
    idx = x->size[1];
    rt = b_y1->size[0];
    b_y1->size[0] = idx;
    emxEnsureCapacity(&b_st, (emxArray__common *)b_y1, rt, (int32_T)sizeof
                      (creal_T), &e_emlrtRTEI);
    idx = x->size[1];
    if (x->size[1] > idx) {
      idx = b_y1->size[0];
      rt = b_y1->size[0];
      b_y1->size[0] = idx;
      emxEnsureCapacity(&b_st, (emxArray__common *)b_y1, rt, (int32_T)sizeof
                        (creal_T), &f_emlrtRTEI);
      for (rt = 0; rt < idx; rt++) {
        b_y1->data[rt].re = 0.0;
        b_y1->data[rt].im = 0.0;
      }
    }
  } else {
    useRadix2 = ((x->size[1] & (x->size[1] - 1)) == 0);
    c_st.site = &t_emlrtRSI;
    get_algo_sizes(&c_st, x->size[1], useRadix2, &N2blue, &idx);
    c_st.site = &u_emlrtRSI;
    generate_twiddle_tables(&c_st, idx, useRadix2, costab, sintab, sintabinv);
    if (useRadix2) {
      idx = x->size[1];
      b_x[0] = idx;
      c_x = *x;
      c_x.size = (int32_T *)&b_x;
      c_x.numDimensions = 1;
      c_st.site = &v_emlrtRSI;
      r2br_r2dit_trig(&c_st, &c_x, x->size[1], costab, sintab, b_y1);
    } else {
      c_st.site = &w_emlrtRSI;
      d_st.site = &lb_emlrtRSI;
      e_st.site = &nb_emlrtRSI;
      nInt2m1 = (x->size[1] + x->size[1]) - 1;
      rt = wwc->size[0];
      wwc->size[0] = nInt2m1;
      emxEnsureCapacity(&e_st, (emxArray__common *)wwc, rt, (int32_T)sizeof
                        (creal_T), &f_emlrtRTEI);
      idx = x->size[1];
      rt = 0;
      wwc->data[x->size[1] - 1].re = 1.0;
      wwc->data[x->size[1] - 1].im = 0.0;
      nInt2 = x->size[1] << 1;
      for (k = 1; k < n1; k++) {
        b_y = (k << 1) - 1;
        if (nInt2 - rt <= b_y) {
          rt += b_y - nInt2;
        } else {
          rt += b_y;
        }

        nt_im = -3.1415926535897931 * (real_T)rt / (real_T)x->size[1];
        if (nt_im == 0.0) {
          nt_re = 1.0;
          nt_im = 0.0;
        } else {
          nt_re = muDoubleScalarCos(nt_im);
          nt_im = muDoubleScalarSin(nt_im);
        }

        wwc->data[idx - 2].re = nt_re;
        wwc->data[idx - 2].im = -nt_im;
        idx--;
      }

      idx = 0;
      for (k = nInt2m1 - 1; k >= n1; k--) {
        wwc->data[k] = wwc->data[idx];
        idx++;
      }

      d_st.site = &mb_emlrtRSI;
      idx = x->size[1];
      nInt2 = muIntScalarMin_sint32(n1, idx);
      idx = x->size[1];
      rt = b_y1->size[0];
      b_y1->size[0] = idx;
      emxEnsureCapacity(&d_st, (emxArray__common *)b_y1, rt, (int32_T)sizeof
                        (creal_T), &g_emlrtRTEI);
      idx = x->size[1];
      if (x->size[1] > idx) {
        idx = b_y1->size[0];
        rt = b_y1->size[0];
        b_y1->size[0] = idx;
        emxEnsureCapacity(&d_st, (emxArray__common *)b_y1, rt, (int32_T)sizeof
                          (creal_T), &f_emlrtRTEI);
        for (rt = 0; rt < idx; rt++) {
          b_y1->data[rt].re = 0.0;
          b_y1->data[rt].im = 0.0;
        }
      }

      idx = 0;
      e_st.site = &pb_emlrtRSI;
      if ((!(1 > nInt2)) && (nInt2 > 2147483646)) {
        f_st.site = &m_emlrtRSI;
        check_forloop_overflow_error(&f_st);
      }

      for (k = 0; k + 1 <= nInt2; k++) {
        nt_re = wwc->data[(n1 + k) - 1].re;
        nt_im = wwc->data[(n1 + k) - 1].im;
        wwc_re = x->data[idx].re;
        fv_im = x->data[idx].im;
        wwc_im = x->data[idx].im;
        fv_re = x->data[idx].re;
        b_y1->data[k].re = nt_re * wwc_re + nt_im * fv_im;
        b_y1->data[k].im = nt_re * wwc_im - nt_im * fv_re;
        idx++;
      }

      e_st.site = &qb_emlrtRSI;
      useRadix2 = ((!(nInt2 + 1 > x->size[1])) && (x->size[1] > 2147483646));
      if (useRadix2) {
        f_st.site = &m_emlrtRSI;
        check_forloop_overflow_error(&f_st);
      }

      while (nInt2 + 1 <= n1) {
        b_y1->data[nInt2].re = 0.0;
        b_y1->data[nInt2].im = 0.0;
        nInt2++;
      }

      e_st.site = &rb_emlrtRSI;
      r2br_r2dit_trig_impl(&e_st, b_y1, N2blue, costab, sintab, fy);
      e_st.site = &sb_emlrtRSI;
      r2br_r2dit_trig(&e_st, wwc, N2blue, costab, sintab, fv);
      rt = fy->size[0];
      emxEnsureCapacity(&d_st, (emxArray__common *)fy, rt, (int32_T)sizeof
                        (creal_T), &f_emlrtRTEI);
      idx = fy->size[0];
      for (rt = 0; rt < idx; rt++) {
        nt_re = fy->data[rt].re;
        nt_im = fy->data[rt].im;
        b_fv_re = fv->data[rt].re;
        b_fv_im = fv->data[rt].im;
        fy->data[rt].re = nt_re * b_fv_re - nt_im * b_fv_im;
        fy->data[rt].im = nt_re * b_fv_im + nt_im * b_fv_re;
      }

      e_st.site = &tb_emlrtRSI;
      b_r2br_r2dit_trig(&e_st, fy, N2blue, costab, sintabinv, fv);
      idx = 0;
      e_st.site = &ub_emlrtRSI;
      useRadix2 = ((!(x->size[1] > wwc->size[0])) && (wwc->size[0] > 2147483646));
      if (useRadix2) {
        f_st.site = &m_emlrtRSI;
        check_forloop_overflow_error(&f_st);
      }

      for (k = x->size[1] - 1; k + 1 <= wwc->size[0]; k++) {
        nt_re = wwc->data[k].re;
        b_fv_re = fv->data[k].re;
        nt_im = wwc->data[k].im;
        b_fv_im = fv->data[k].im;
        wwc_re = wwc->data[k].re;
        fv_im = fv->data[k].im;
        wwc_im = wwc->data[k].im;
        fv_re = fv->data[k].re;
        b_y1->data[idx].re = nt_re * b_fv_re + nt_im * b_fv_im;
        b_y1->data[idx].im = wwc_re * fv_im - wwc_im * fv_re;
        idx++;
      }
    }
  }

  emxFree_creal_T(&fv);
  emxFree_creal_T(&fy);
  emxFree_creal_T(&wwc);
  emxFree_real_T(&sintabinv);
  emxFree_real_T(&sintab);
  emxFree_real_T(&costab);
  rt = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = x->size[1];
  emxEnsureCapacity(&st, (emxArray__common *)y, rt, (int32_T)sizeof(creal_T),
                    &f_emlrtRTEI);
  idx = x->size[1];
  for (rt = 0; rt < idx; rt++) {
    y->data[rt] = b_y1->data[rt];
  }

  emxFree_creal_T(&b_y1);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (fft.cpp) */
