/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * f_scalar_ssfm.cpp
 *
 * Code generation for function 'f_scalar_ssfm'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "f_scalar_ssfm.h"
#include "error.h"
#include "eml_int_forloop_overflow_check.h"
#include "f_scalar_ssfm_emxutil.h"
#include "exp.h"
#include "power.h"
#include "fft.h"
#include "scalexpAlloc.h"
#include "f_scalar_ssfm_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 6,     /* lineNo */
  "f_scalar_ssfm",                     /* fcnName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 18,  /* lineNo */
  "f_scalar_ssfm",                     /* fcnName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 26,  /* lineNo */
  "f_scalar_ssfm",                     /* fcnName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 27,  /* lineNo */
  "f_scalar_ssfm",                     /* fcnName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 36,  /* lineNo */
  "f_scalar_ssfm",                     /* fcnName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 37,  /* lineNo */
  "f_scalar_ssfm",                     /* fcnName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 45,  /* lineNo */
  "f_scalar_ssfm",                     /* fcnName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 55,  /* lineNo */
  "f_scalar_ssfm",                     /* fcnName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 56,  /* lineNo */
  "f_scalar_ssfm",                     /* fcnName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pathName */
};

static emlrtRSInfo jb_emlrtRSI = { 430,/* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

static emlrtRSInfo kb_emlrtRSI = { 308,/* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

static emlrtRSInfo ob_emlrtRSI = { 40, /* lineNo */
  "bluestein_setup",                   /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\bluestein_setup.m"/* pathName */
};

static emlrtRSInfo vb_emlrtRSI = { 25, /* lineNo */
  "ifft",                              /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\ifft.m"/* pathName */
};

static emlrtRSInfo yb_emlrtRSI = { 244,/* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

static emlrtRSInfo ac_emlrtRSI = { 62, /* lineNo */
  "f_scalar_ssfm",                     /* fcnName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pathName */
};

static emlrtRSInfo bc_emlrtRSI = { 63, /* lineNo */
  "f_scalar_ssfm",                     /* fcnName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pathName */
};

static emlrtRTEInfo emlrtRTEI = { 1,   /* lineNo */
  20,                                  /* colNo */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pName */
};

static emlrtRTEInfo d_emlrtRTEI = { 53,/* lineNo */
  16,                                  /* colNo */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pName */
};

static emlrtRTEInfo n_emlrtRTEI = { 296,/* lineNo */
  9,                                   /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

static emlrtRTEInfo o_emlrtRTEI = { 417,/* lineNo */
  14,                                  /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

static emlrtRTEInfo p_emlrtRTEI = { 293,/* lineNo */
  9,                                   /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

static emlrtRTEInfo q_emlrtRTEI = { 261,/* lineNo */
  14,                                  /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

static emlrtRTEInfo r_emlrtRTEI = { 1, /* lineNo */
  14,                                  /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

static emlrtRTEInfo s_emlrtRTEI = { 60,/* lineNo */
  16,                                  /* colNo */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pName */
};

static emlrtRTEInfo u_emlrtRTEI = { 5, /* lineNo */
  5,                                   /* colNo */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pName */
};

static emlrtRTEInfo v_emlrtRTEI = { 6, /* lineNo */
  5,                                   /* colNo */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pName */
};

static emlrtRTEInfo w_emlrtRTEI = { 8, /* lineNo */
  5,                                   /* colNo */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pName */
};

static emlrtRTEInfo x_emlrtRTEI = { 18,/* lineNo */
  5,                                   /* colNo */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pName */
};

static emlrtRTEInfo y_emlrtRTEI = { 55,/* lineNo */
  5,                                   /* colNo */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pName */
};

static emlrtRTEInfo eb_emlrtRTEI = { 109,/* lineNo */
  5,                                   /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

static emlrtECInfo emlrtECI = { 2,     /* nDims */
  6,                                   /* lineNo */
  13,                                  /* colNo */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  27,                                  /* lineNo */
  28,                                  /* colNo */
  "xi",                                /* aName */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo fb_emlrtRTEI = { 30,/* lineNo */
  11,                                  /* colNo */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pName */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  37,                                  /* lineNo */
  32,                                  /* colNo */
  "xi",                                /* aName */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m",/* pName */
  0                                    /* checkKind */
};

static emlrtECInfo b_emlrtECI = { 2,   /* nDims */
  56,                                  /* lineNo */
  16,                                  /* colNo */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pName */
};

static emlrtECInfo c_emlrtECI = { 2,   /* nDims */
  62,                                  /* lineNo */
  11,                                  /* colNo */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pName */
};

static emlrtECInfo d_emlrtECI = { 2,   /* nDims */
  63,                                  /* lineNo */
  11,                                  /* colNo */
  "f_scalar_ssfm",                     /* fName */
  "C:\\Users\\mene9\\Documents\\MATLAB\\My Simulator\\Link\\f_scalar_ssfm.m"/* pName */
};

/* Function Declarations */
static void b_generate_twiddle_tables(const emlrtStack *sp, int32_T nRows,
  boolean_T useRadix2, emxArray_real_T *costab, emxArray_real_T *sintab,
  emxArray_real_T *sintabinv);
static void eml_fft(const emlrtStack *sp, const emxArray_creal_T *x, int32_T n,
                    emxArray_creal_T *y);
static void scalar_lin_step(const emlrtStack *sp, const emxArray_real_T *betaxdz,
  emxArray_creal_T *ux);
static void scalar_nl_step(const emlrtStack *sp, emxArray_creal_T *ux, real_T xi);

/* Function Definitions */
static void b_generate_twiddle_tables(const emlrtStack *sp, int32_T nRows,
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
      sintab->data[k] = costab1q->data[nRowsD4 - k];
    }

    for (k = costab1q->size[1]; k <= n2; k++) {
      costab->data[k] = -costab1q->data[n2 - k];
      sintab->data[k] = costab1q->data[k - nRowsD4];
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

static void eml_fft(const emlrtStack *sp, const emxArray_creal_T *x, int32_T n,
                    emxArray_creal_T *y)
{
  int32_T idx;
  emxArray_real_T *costab;
  emxArray_real_T *sintab;
  emxArray_real_T *sintabinv;
  emxArray_creal_T *wwc;
  emxArray_creal_T *fy;
  emxArray_creal_T *fv;
  boolean_T useRadix2;
  int32_T rt;
  int32_T N2blue;
  int32_T b_x[1];
  int32_T b_y;
  emxArray_creal_T c_x;
  int32_T nInt2m1;
  int32_T nInt2;
  int32_T k;
  real_T denom_im;
  real_T denom_re;
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
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  idx = x->size[0];
  emxInit_real_T(sp, &costab, 2, &r_emlrtRTEI, true);
  emxInit_real_T(sp, &sintab, 2, &r_emlrtRTEI, true);
  emxInit_real_T(sp, &sintabinv, 2, &r_emlrtRTEI, true);
  emxInit_creal_T(sp, &wwc, 1, &i_emlrtRTEI, true);
  emxInit_creal_T(sp, &fy, 1, &j_emlrtRTEI, true);
  emxInit_creal_T(sp, &fv, 1, &k_emlrtRTEI, true);
  if ((idx == 0) || (n == 0)) {
    rt = y->size[0];
    y->size[0] = n;
    emxEnsureCapacity(sp, (emxArray__common *)y, rt, (int32_T)sizeof(creal_T),
                      &e_emlrtRTEI);
    idx = x->size[0];
    if (n > idx) {
      b_y = y->size[0];
      rt = y->size[0];
      y->size[0] = b_y;
      emxEnsureCapacity(sp, (emxArray__common *)y, rt, (int32_T)sizeof(creal_T),
                        &r_emlrtRTEI);
      for (rt = 0; rt < b_y; rt++) {
        y->data[rt].re = 0.0;
        y->data[rt].im = 0.0;
      }
    }
  } else {
    useRadix2 = ((n & (n - 1)) == 0);
    st.site = &t_emlrtRSI;
    get_algo_sizes(&st, n, useRadix2, &N2blue, &idx);
    st.site = &u_emlrtRSI;
    b_generate_twiddle_tables(&st, idx, useRadix2, costab, sintab, sintabinv);
    if (useRadix2) {
      b_x[0] = x->size[0];
      c_x = *x;
      c_x.size = (int32_T *)&b_x;
      c_x.numDimensions = 1;
      st.site = &v_emlrtRSI;
      b_r2br_r2dit_trig(&st, &c_x, n, costab, sintab, y);
    } else {
      st.site = &w_emlrtRSI;
      b_st.site = &lb_emlrtRSI;
      c_st.site = &nb_emlrtRSI;
      nInt2m1 = (n + n) - 1;
      rt = wwc->size[0];
      wwc->size[0] = nInt2m1;
      emxEnsureCapacity(&c_st, (emxArray__common *)wwc, rt, (int32_T)sizeof
                        (creal_T), &r_emlrtRTEI);
      idx = n;
      rt = 0;
      wwc->data[n - 1].re = 1.0;
      wwc->data[n - 1].im = 0.0;
      nInt2 = n << 1;
      d_st.site = &ob_emlrtRSI;
      if ((!(1 > n - 1)) && (n - 1 > 2147483646)) {
        e_st.site = &m_emlrtRSI;
        check_forloop_overflow_error(&e_st);
      }

      for (k = 1; k < n; k++) {
        b_y = (k << 1) - 1;
        if (nInt2 - rt <= b_y) {
          rt += b_y - nInt2;
        } else {
          rt += b_y;
        }

        denom_im = 3.1415926535897931 * (real_T)rt / (real_T)n;
        if (denom_im == 0.0) {
          denom_re = 1.0;
          denom_im = 0.0;
        } else {
          denom_re = muDoubleScalarCos(denom_im);
          denom_im = muDoubleScalarSin(denom_im);
        }

        wwc->data[idx - 2].re = denom_re;
        wwc->data[idx - 2].im = -denom_im;
        idx--;
      }

      idx = 0;
      for (k = nInt2m1 - 1; k >= n; k--) {
        wwc->data[k] = wwc->data[idx];
        idx++;
      }

      b_st.site = &mb_emlrtRSI;
      idx = x->size[0];
      nInt2 = muIntScalarMin_sint32(n, idx);
      rt = y->size[0];
      y->size[0] = n;
      emxEnsureCapacity(&b_st, (emxArray__common *)y, rt, (int32_T)sizeof
                        (creal_T), &g_emlrtRTEI);
      idx = x->size[0];
      if (n > idx) {
        b_y = y->size[0];
        rt = y->size[0];
        y->size[0] = b_y;
        emxEnsureCapacity(&b_st, (emxArray__common *)y, rt, (int32_T)sizeof
                          (creal_T), &r_emlrtRTEI);
        for (rt = 0; rt < b_y; rt++) {
          y->data[rt].re = 0.0;
          y->data[rt].im = 0.0;
        }
      }

      idx = 0;
      c_st.site = &pb_emlrtRSI;
      if ((!(1 > nInt2)) && (nInt2 > 2147483646)) {
        d_st.site = &m_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }

      for (k = 0; k + 1 <= nInt2; k++) {
        denom_re = wwc->data[(n + k) - 1].re;
        denom_im = wwc->data[(n + k) - 1].im;
        wwc_re = x->data[idx].re;
        fv_im = x->data[idx].im;
        wwc_im = x->data[idx].im;
        fv_re = x->data[idx].re;
        y->data[k].re = denom_re * wwc_re + denom_im * fv_im;
        y->data[k].im = denom_re * wwc_im - denom_im * fv_re;
        idx++;
      }

      c_st.site = &qb_emlrtRSI;
      if ((!(nInt2 + 1 > n)) && (n > 2147483646)) {
        d_st.site = &m_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }

      while (nInt2 + 1 <= n) {
        y->data[nInt2].re = 0.0;
        y->data[nInt2].im = 0.0;
        nInt2++;
      }

      c_st.site = &rb_emlrtRSI;
      r2br_r2dit_trig_impl(&c_st, y, N2blue, costab, sintab, fy);
      c_st.site = &sb_emlrtRSI;
      r2br_r2dit_trig(&c_st, wwc, N2blue, costab, sintab, fv);
      rt = fy->size[0];
      emxEnsureCapacity(&b_st, (emxArray__common *)fy, rt, (int32_T)sizeof
                        (creal_T), &r_emlrtRTEI);
      idx = fy->size[0];
      for (rt = 0; rt < idx; rt++) {
        denom_re = fy->data[rt].re;
        denom_im = fy->data[rt].im;
        b_fv_re = fv->data[rt].re;
        b_fv_im = fv->data[rt].im;
        fy->data[rt].re = denom_re * b_fv_re - denom_im * b_fv_im;
        fy->data[rt].im = denom_re * b_fv_im + denom_im * b_fv_re;
      }

      c_st.site = &tb_emlrtRSI;
      b_r2br_r2dit_trig(&c_st, fy, N2blue, costab, sintabinv, fv);
      idx = 0;
      c_st.site = &yb_emlrtRSI;
      useRadix2 = ((!(n > wwc->size[0])) && (wwc->size[0] > 2147483646));
      if (useRadix2) {
        d_st.site = &m_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }

      for (k = n - 1; k + 1 <= wwc->size[0]; k++) {
        denom_re = wwc->data[k].re;
        b_fv_re = fv->data[k].re;
        denom_im = wwc->data[k].im;
        b_fv_im = fv->data[k].im;
        wwc_re = wwc->data[k].re;
        fv_im = fv->data[k].im;
        wwc_im = wwc->data[k].im;
        fv_re = fv->data[k].re;
        y->data[idx].re = denom_re * b_fv_re + denom_im * b_fv_im;
        y->data[idx].im = wwc_re * fv_im - wwc_im * fv_re;
        denom_re = wwc->data[k].re;
        b_fv_re = fv->data[k].re;
        denom_im = wwc->data[k].im;
        b_fv_im = fv->data[k].im;
        wwc_re = wwc->data[k].re;
        fv_im = fv->data[k].im;
        wwc_im = wwc->data[k].im;
        fv_re = fv->data[k].re;
        y->data[idx].re = denom_re * b_fv_re + denom_im * b_fv_im;
        y->data[idx].im = wwc_re * fv_im - wwc_im * fv_re;
        denom_re = y->data[idx].re;
        denom_im = y->data[idx].im;
        if (denom_im == 0.0) {
          y->data[idx].re = denom_re / (real_T)n;
          y->data[idx].im = 0.0;
        } else if (denom_re == 0.0) {
          y->data[idx].re = 0.0;
          y->data[idx].im = denom_im / (real_T)n;
        } else {
          y->data[idx].re = denom_re / (real_T)n;
          y->data[idx].im = denom_im / (real_T)n;
        }

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
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

static void scalar_lin_step(const emlrtStack *sp, const emxArray_real_T *betaxdz,
  emxArray_creal_T *ux)
{
  emxArray_creal_T *Hf;
  int32_T x;
  int32_T loop_ub;
  emxArray_creal_T *b_x;
  int32_T c_x[2];
  int32_T b_Hf[2];
  real_T x_re;
  real_T x_im;
  real_T Hf_re;
  boolean_T b0;
  real_T Hf_im;
  emxArray_creal_T *b_y1;
  int32_T d_x[2];
  emxArray_creal_T e_x;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_creal_T1(sp, &Hf, 2, &y_emlrtRTEI, true);
  x = Hf->size[0] * Hf->size[1];
  Hf->size[0] = 1;
  Hf->size[1] = betaxdz->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)Hf, x, (int32_T)sizeof(creal_T),
                    &d_emlrtRTEI);
  loop_ub = betaxdz->size[0] * betaxdz->size[1];
  for (x = 0; x < loop_ub; x++) {
    Hf->data[x].re = betaxdz->data[x] * 0.0;
    Hf->data[x].im = -betaxdz->data[x];
  }

  emxInit_creal_T1(sp, &b_x, 2, &d_emlrtRTEI, true);
  st.site = &p_emlrtRSI;
  b_exp(&st, Hf);
  st.site = &q_emlrtRSI;
  fft(&st, ux, b_x);
  for (x = 0; x < 2; x++) {
    c_x[x] = b_x->size[x];
  }

  for (x = 0; x < 2; x++) {
    b_Hf[x] = Hf->size[x];
  }

  if ((c_x[0] != b_Hf[0]) || (c_x[1] != b_Hf[1])) {
    emlrtSizeEqCheckNDR2012b(&c_x[0], &b_Hf[0], (emlrtECInfo *)&b_emlrtECI, sp);
  }

  st.site = &q_emlrtRSI;
  x = b_x->size[0] * b_x->size[1];
  b_x->size[0] = 1;
  emxEnsureCapacity(&st, (emxArray__common *)b_x, x, (int32_T)sizeof(creal_T),
                    &d_emlrtRTEI);
  loop_ub = b_x->size[0];
  x = b_x->size[1];
  loop_ub *= x;
  for (x = 0; x < loop_ub; x++) {
    x_re = b_x->data[x].re;
    x_im = b_x->data[x].im;
    Hf_re = Hf->data[x].re;
    Hf_im = Hf->data[x].im;
    b_x->data[x].re = x_re * Hf_re - x_im * Hf_im;
    b_x->data[x].im = x_re * Hf_im + x_im * Hf_re;
  }

  emxFree_creal_T(&Hf);
  b_st.site = &vb_emlrtRSI;
  if ((b_x->size[1] == 1) || (b_x->size[1] != 1)) {
    b0 = true;
  } else {
    b0 = false;
  }

  if (!b0) {
    emlrtErrorWithMessageIdR2012b(&b_st, &db_emlrtRTEI,
      "Coder:toolbox:autoDimIncompatibility", 0);
  }

  emxInit_creal_T(&b_st, &b_y1, 1, &h_emlrtRTEI, true);
  d_x[0] = b_x->size[1];
  d_x[1] = 1;
  e_x = *b_x;
  e_x.size = (int32_T *)&d_x;
  e_x.numDimensions = 1;
  c_st.site = &s_emlrtRSI;
  eml_fft(&c_st, &e_x, b_x->size[1], b_y1);
  loop_ub = b_x->size[1];
  x = ux->size[0] * ux->size[1];
  ux->size[0] = 1;
  ux->size[1] = loop_ub;
  emxEnsureCapacity(&b_st, (emxArray__common *)ux, x, (int32_T)sizeof(creal_T),
                    &d_emlrtRTEI);
  emxFree_creal_T(&b_x);
  for (x = 0; x < loop_ub; x++) {
    ux->data[ux->size[0] * x] = b_y1->data[x];
  }

  emxFree_creal_T(&b_y1);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

static void scalar_nl_step(const emlrtStack *sp, emxArray_creal_T *ux, real_T xi)
{
  emxArray_real_T *b_ux;
  int32_T i4;
  int32_T loop_ub;
  emxArray_real_T *r0;
  emxArray_real_T *c_ux;
  emxArray_real_T *r1;
  int32_T iv2[2];
  int32_T iv3[2];
  emxArray_creal_T *r2;
  real_T y_re;
  real_T y_im;
  int32_T d_ux[2];
  int32_T iv4[2];
  real_T re;
  real_T im;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T(sp, &b_ux, 2, &s_emlrtRTEI, true);
  i4 = b_ux->size[0] * b_ux->size[1];
  b_ux->size[0] = 1;
  b_ux->size[1] = ux->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)b_ux, i4, (int32_T)sizeof(real_T),
                    &s_emlrtRTEI);
  loop_ub = ux->size[0] * ux->size[1];
  for (i4 = 0; i4 < loop_ub; i4++) {
    b_ux->data[i4] = ux->data[i4].re;
  }

  emxInit_real_T(sp, &r0, 2, &s_emlrtRTEI, true);
  emxInit_real_T(sp, &c_ux, 2, &s_emlrtRTEI, true);
  st.site = &ac_emlrtRSI;
  power(&st, b_ux, r0);
  i4 = c_ux->size[0] * c_ux->size[1];
  c_ux->size[0] = 1;
  c_ux->size[1] = ux->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)c_ux, i4, (int32_T)sizeof(real_T),
                    &s_emlrtRTEI);
  loop_ub = ux->size[0] * ux->size[1];
  emxFree_real_T(&b_ux);
  for (i4 = 0; i4 < loop_ub; i4++) {
    c_ux->data[i4] = ux->data[i4].im;
  }

  emxInit_real_T(sp, &r1, 2, &s_emlrtRTEI, true);
  st.site = &ac_emlrtRSI;
  power(&st, c_ux, r1);
  emxFree_real_T(&c_ux);
  for (i4 = 0; i4 < 2; i4++) {
    iv2[i4] = r0->size[i4];
  }

  for (i4 = 0; i4 < 2; i4++) {
    iv3[i4] = r1->size[i4];
  }

  emxInit_creal_T1(sp, &r2, 2, &s_emlrtRTEI, true);
  if ((iv2[0] != iv3[0]) || (iv2[1] != iv3[1])) {
    emlrtSizeEqCheckNDR2012b(&iv2[0], &iv3[0], (emlrtECInfo *)&c_emlrtECI, sp);
  }

  y_re = xi * 0.0;
  y_im = -xi;
  i4 = r2->size[0] * r2->size[1];
  r2->size[0] = 1;
  r2->size[1] = r0->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)r2, i4, (int32_T)sizeof(creal_T),
                    &s_emlrtRTEI);
  loop_ub = r0->size[0] * r0->size[1];
  for (i4 = 0; i4 < loop_ub; i4++) {
    r2->data[i4].re = (r0->data[i4] + r1->data[i4]) * y_re;
    r2->data[i4].im = (r0->data[i4] + r1->data[i4]) * y_im;
  }

  emxFree_real_T(&r1);
  emxFree_real_T(&r0);
  st.site = &bc_emlrtRSI;
  b_exp(&st, r2);
  for (i4 = 0; i4 < 2; i4++) {
    d_ux[i4] = ux->size[i4];
  }

  for (i4 = 0; i4 < 2; i4++) {
    iv4[i4] = r2->size[i4];
  }

  if ((d_ux[0] != iv4[0]) || (d_ux[1] != iv4[1])) {
    emlrtSizeEqCheckNDR2012b(&d_ux[0], &iv4[0], (emlrtECInfo *)&d_emlrtECI, sp);
  }

  i4 = ux->size[0] * ux->size[1];
  ux->size[0] = 1;
  emxEnsureCapacity(sp, (emxArray__common *)ux, i4, (int32_T)sizeof(creal_T),
                    &s_emlrtRTEI);
  loop_ub = ux->size[1];
  for (i4 = 0; i4 < loop_ub; i4++) {
    y_re = ux->data[ux->size[0] * i4].re;
    y_im = ux->data[ux->size[0] * i4].im;
    re = r2->data[r2->size[0] * i4].re;
    im = r2->data[r2->size[0] * i4].im;
    ux->data[ux->size[0] * i4].re = y_re * re - y_im * im;
    ux->data[ux->size[0] * i4].im = y_re * im + y_im * re;
  }

  emxFree_creal_T(&r2);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void b_r2br_r2dit_trig(const emlrtStack *sp, const emxArray_creal_T *x, int32_T
  n1_unsigned, const emxArray_real_T *costab, const emxArray_real_T *sintab,
  emxArray_creal_T *y)
{
  int32_T SZ1;
  int32_T istart;
  int32_T nRowsD2;
  int32_T nRowsD4;
  int32_T ix;
  int32_T ju;
  int32_T i;
  boolean_T tst;
  real_T temp_re;
  real_T temp_im;
  real_T r;
  int32_T j;
  real_T twid_im;
  int32_T ihi;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &jb_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  SZ1 = x->size[0];
  istart = muIntScalarMin_sint32(SZ1, n1_unsigned) - 1;
  nRowsD2 = n1_unsigned / 2;
  nRowsD4 = nRowsD2 / 2;
  ix = y->size[0];
  y->size[0] = n1_unsigned;
  emxEnsureCapacity(&st, (emxArray__common *)y, ix, (int32_T)sizeof(creal_T),
                    &n_emlrtRTEI);
  if (n1_unsigned > x->size[0]) {
    SZ1 = y->size[0];
    ix = y->size[0];
    y->size[0] = SZ1;
    emxEnsureCapacity(&st, (emxArray__common *)y, ix, (int32_T)sizeof(creal_T),
                      &o_emlrtRTEI);
    for (ix = 0; ix < SZ1; ix++) {
      y->data[ix].re = 0.0;
      y->data[ix].im = 0.0;
    }
  }

  ix = 0;
  ju = 0;
  SZ1 = 0;
  b_st.site = &kb_emlrtRSI;
  if ((!(1 > istart)) && (istart > 2147483646)) {
    c_st.site = &m_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }

  for (i = 1; i <= istart; i++) {
    y->data[SZ1] = x->data[ix];
    SZ1 = n1_unsigned;
    tst = true;
    while (tst) {
      SZ1 >>= 1;
      ju ^= SZ1;
      tst = ((ju & SZ1) == 0);
    }

    SZ1 = ju;
    ix++;
  }

  y->data[SZ1] = x->data[ix];
  if (n1_unsigned > 1) {
    for (i = 0; i <= n1_unsigned - 2; i += 2) {
      temp_re = y->data[i + 1].re;
      temp_im = y->data[i + 1].im;
      y->data[i + 1].re = y->data[i].re - y->data[i + 1].re;
      y->data[i + 1].im = y->data[i].im - y->data[i + 1].im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }
  }

  SZ1 = 2;
  ix = 4;
  ju = 1 + ((nRowsD4 - 1) << 2);
  while (nRowsD4 > 0) {
    for (i = 0; i < ju; i += ix) {
      temp_re = y->data[i + SZ1].re;
      temp_im = y->data[i + SZ1].im;
      y->data[i + SZ1].re = y->data[i].re - temp_re;
      y->data[i + SZ1].im = y->data[i].im - temp_im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }

    istart = 1;
    for (j = nRowsD4; j < nRowsD2; j += nRowsD4) {
      r = costab->data[j];
      twid_im = sintab->data[j];
      i = istart;
      ihi = istart + ju;
      while (i < ihi) {
        temp_re = r * y->data[i + SZ1].re - twid_im * y->data[i + SZ1].im;
        temp_im = r * y->data[i + SZ1].im + twid_im * y->data[i + SZ1].re;
        y->data[i + SZ1].re = y->data[i].re - temp_re;
        y->data[i + SZ1].im = y->data[i].im - temp_im;
        y->data[i].re += temp_re;
        y->data[i].im += temp_im;
        i += ix;
      }

      istart++;
    }

    nRowsD4 /= 2;
    SZ1 = ix;
    ix <<= 1;
    ju -= SZ1;
  }

  if (y->size[0] > 1) {
    r = 1.0 / (real_T)y->size[0];
    ix = y->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)y, ix, (int32_T)sizeof(creal_T),
                      &o_emlrtRTEI);
    SZ1 = y->size[0];
    for (ix = 0; ix < SZ1; ix++) {
      y->data[ix].re *= r;
      y->data[ix].im *= r;
    }
  }
}

void f_scalar_ssfm(const emlrtStack *sp, const Channel *ch, real_T Pavg, Signal *
                   sig)
{
  emxArray_real_T *omega;
  real_T Leff;
  int32_T i3;
  int32_T loop_ub;
  emxArray_real_T *beta;
  int32_T b_beta;
  emxArray_real_T *xi;
  emxArray_real_T *y;
  boolean_T overflow;
  int32_T k;
  int32_T c_beta[2];
  int32_T b_y[2];
  emxArray_creal_T *ux;
  emxArray_real_T *d_beta;
  int32_T i;
  emxArray_real_T *e_beta;
  emxArray_real_T *f_beta;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
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
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T(sp, &omega, 2, &u_emlrtRTEI, true);

  /* SCALAR_SSFM Summary of this function goes here */
  /*    Detailed explanation goes here */
  Leff = 6.2831853071795862 * sig->SYMBOLRATE;
  i3 = omega->size[0] * omega->size[1];
  omega->size[0] = 1;
  omega->size[1] = sig->FN->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)omega, i3, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = sig->FN->size[0] * sig->FN->size[1];
  for (i3 = 0; i3 < loop_ub; i3++) {
    omega->data[i3] = Leff * sig->FN->data[i3] * 1.0E+9;
  }

  emxInit_real_T(sp, &beta, 2, &v_emlrtRTEI, true);

  /*  [rad/s] */
  st.site = &emlrtRSI;
  power(&st, omega, beta);
  i3 = beta->size[0] * beta->size[1];
  beta->size[0] = 1;
  emxEnsureCapacity(sp, (emxArray__common *)beta, i3, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = beta->size[0];
  b_beta = beta->size[1];
  loop_ub *= b_beta;
  for (i3 = 0; i3 < loop_ub; i3++) {
    beta->data[i3] = 0.5 * beta->data[i3] * ch->b2;
  }

  emxInit_real_T(sp, &xi, 2, &x_emlrtRTEI, true);
  st.site = &emlrtRSI;
  b_st.site = &h_emlrtRSI;
  c_st.site = &i_emlrtRSI;
  i3 = xi->size[0] * xi->size[1];
  xi->size[0] = 1;
  xi->size[1] = omega->size[1];
  emxEnsureCapacity(&c_st, (emxArray__common *)xi, i3, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = omega->size[0] * omega->size[1];
  for (i3 = 0; i3 < loop_ub; i3++) {
    xi->data[i3] = omega->data[i3];
  }

  emxInit_real_T(&c_st, &y, 2, &emlrtRTEI, true);
  d_st.site = &k_emlrtRSI;
  i3 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = omega->size[1];
  emxEnsureCapacity(&d_st, (emxArray__common *)y, i3, (int32_T)sizeof(real_T),
                    &c_emlrtRTEI);
  if (!dimagree(y, omega)) {
    emlrtErrorWithMessageIdR2012b(&d_st, &ab_emlrtRTEI, "MATLAB:dimagree", 0);
  }

  loop_ub = omega->size[1];
  d_st.site = &l_emlrtRSI;
  overflow = ((!(1 > omega->size[1])) && (omega->size[1] > 2147483646));
  if (overflow) {
    e_st.site = &m_emlrtRSI;
    check_forloop_overflow_error(&e_st);
  }

  emlrtEnterParallelRegion(&c_st, omp_in_parallel());

#pragma omp parallel for \
 num_threads(emlrtAllocRegionTLSs(c_st.tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs()))

  for (k = 1; k <= loop_ub; k++) {
    y->data[k - 1] = muDoubleScalarPower(xi->data[k - 1], 3.0);
  }

  emlrtExitParallelRegion(&c_st, omp_in_parallel());
  if (fltpower_domain_error(omega, 3.0)) {
    c_st.site = &j_emlrtRSI;
    error(&c_st);
  }

  i3 = y->size[0] * y->size[1];
  y->size[0] = 1;
  emxEnsureCapacity(sp, (emxArray__common *)y, i3, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = y->size[0];
  b_beta = y->size[1];
  loop_ub *= b_beta;
  for (i3 = 0; i3 < loop_ub; i3++) {
    y->data[i3] = y->data[i3] * ch->b3 / 6.0;
  }

  for (i3 = 0; i3 < 2; i3++) {
    c_beta[i3] = beta->size[i3];
  }

  for (i3 = 0; i3 < 2; i3++) {
    b_y[i3] = y->size[i3];
  }

  if ((c_beta[0] != b_y[0]) || (c_beta[1] != b_y[1])) {
    emlrtSizeEqCheckNDR2012b(&c_beta[0], &b_y[0], (emlrtECInfo *)&emlrtECI, sp);
  }

  i3 = beta->size[0] * beta->size[1];
  beta->size[0] = 1;
  emxEnsureCapacity(sp, (emxArray__common *)beta, i3, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = beta->size[0];
  b_beta = beta->size[1];
  loop_ub *= b_beta;
  for (i3 = 0; i3 < loop_ub; i3++) {
    beta->data[i3] += y->data[i3];
  }

  emxFree_real_T(&y);

  /* ;get(sig,'FIELDX'); */
  if (muDoubleScalarAbs(ch->alphalin * ch->dz) > 1.0E-6) {
    Leff = (1.0 - muDoubleScalarExp(-ch->alphalin * ch->dz)) / ch->alphalin;
  } else {
    Leff = ch->dz;
  }

  if (muDoubleScalarIsNaN(ch->nstep - 1.0)) {
    i3 = omega->size[0] * omega->size[1];
    omega->size[0] = 1;
    omega->size[1] = 1;
    emxEnsureCapacity(sp, (emxArray__common *)omega, i3, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    omega->data[0] = rtNaN;
  } else if (ch->nstep - 1.0 < 0.0) {
    i3 = omega->size[0] * omega->size[1];
    omega->size[0] = 1;
    omega->size[1] = 0;
    emxEnsureCapacity(sp, (emxArray__common *)omega, i3, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
  } else if (muDoubleScalarIsInf(ch->nstep - 1.0) && (0.0 == ch->nstep - 1.0)) {
    i3 = omega->size[0] * omega->size[1];
    omega->size[0] = 1;
    omega->size[1] = 1;
    emxEnsureCapacity(sp, (emxArray__common *)omega, i3, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    omega->data[0] = rtNaN;
  } else {
    i3 = omega->size[0] * omega->size[1];
    omega->size[0] = 1;
    omega->size[1] = (int32_T)muDoubleScalarFloor(ch->nstep - 1.0) + 1;
    emxEnsureCapacity(sp, (emxArray__common *)omega, i3, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    loop_ub = (int32_T)muDoubleScalarFloor(ch->nstep - 1.0);
    for (i3 = 0; i3 <= loop_ub; i3++) {
      omega->data[omega->size[0] * i3] = i3;
    }
  }

  i3 = omega->size[0] * omega->size[1];
  omega->size[0] = 1;
  emxEnsureCapacity(sp, (emxArray__common *)omega, i3, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = omega->size[0];
  b_beta = omega->size[1];
  loop_ub *= b_beta;
  for (i3 = 0; i3 < loop_ub; i3++) {
    omega->data[i3] = -ch->alphalin * (ch->dz * omega->data[i3]);
  }

  st.site = &b_emlrtRSI;
  b_st.site = &n_emlrtRSI;
  i3 = xi->size[0] * xi->size[1];
  xi->size[0] = 1;
  xi->size[1] = omega->size[1];
  emxEnsureCapacity(&b_st, (emxArray__common *)xi, i3, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = omega->size[0] * omega->size[1];
  for (i3 = 0; i3 < loop_ub; i3++) {
    xi->data[i3] = omega->data[i3];
  }

  c_st.site = &o_emlrtRSI;
  overflow = ((!(1 > omega->size[1])) && (omega->size[1] > 2147483646));
  if (overflow) {
    d_st.site = &m_emlrtRSI;
    check_forloop_overflow_error(&d_st);
  }

  for (loop_ub = 0; loop_ub + 1 <= omega->size[1]; loop_ub++) {
    xi->data[loop_ub] = muDoubleScalarExp(xi->data[loop_ub]);
  }

  emxFree_real_T(&omega);
  Leff *= ch->gamma;
  i3 = xi->size[0] * xi->size[1];
  xi->size[0] = 1;
  emxEnsureCapacity(sp, (emxArray__common *)xi, i3, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = xi->size[0];
  b_beta = xi->size[1];
  loop_ub *= b_beta;
  for (i3 = 0; i3 < loop_ub; i3++) {
    xi->data[i3] = Leff * xi->data[i3] * Pavg;
  }

  emxInit_creal_T1(sp, &ux, 2, &w_emlrtRTEI, true);
  Leff = ch->dz / 2.0;

  /*                            HALF DZ GVD                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    % */
  /*                                 DZ SPM                      % */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  i3 = ux->size[0] * ux->size[1];
  ux->size[0] = 1;
  ux->size[1] = sig->FIELDX->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)ux, i3, (int32_T)sizeof(creal_T),
                    &emlrtRTEI);
  loop_ub = sig->FIELDX->size[0] * sig->FIELDX->size[1];
  for (i3 = 0; i3 < loop_ub; i3++) {
    ux->data[i3] = sig->FIELDX->data[i3];
  }

  emxInit_real_T(sp, &d_beta, 2, &emlrtRTEI, true);
  i3 = d_beta->size[0] * d_beta->size[1];
  d_beta->size[0] = 1;
  d_beta->size[1] = beta->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)d_beta, i3, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = beta->size[0] * beta->size[1];
  for (i3 = 0; i3 < loop_ub; i3++) {
    d_beta->data[i3] = beta->data[i3] * Leff;
  }

  st.site = &c_emlrtRSI;
  scalar_lin_step(&st, d_beta, ux);
  i3 = xi->size[1];
  if (!(1 <= i3)) {
    emlrtDynamicBoundsCheckR2012b(1, 1, i3, (emlrtBCInfo *)&emlrtBCI, sp);
  }

  st.site = &d_emlrtRSI;
  scalar_nl_step(&st, ux, xi->data[0]);

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  emlrtForLoopVectorCheckR2012b(2.0, 1.0, ch->nstep, mxDOUBLE_CLASS, (int32_T)
    (ch->nstep + -1.0), (emlrtRTEInfo *)&fb_emlrtRTEI, sp);
  i = 0;
  emxFree_real_T(&d_beta);
  emxInit_real_T(sp, &e_beta, 2, &emlrtRTEI, true);
  while (i <= (int32_T)(ch->nstep + -1.0) - 1) {
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /*                            DZ GVD                       % */
    /*                            DZ SPM                       % */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    i3 = e_beta->size[0] * e_beta->size[1];
    e_beta->size[0] = 1;
    e_beta->size[1] = beta->size[1];
    emxEnsureCapacity(sp, (emxArray__common *)e_beta, i3, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    loop_ub = beta->size[0] * beta->size[1];
    for (i3 = 0; i3 < loop_ub; i3++) {
      e_beta->data[i3] = beta->data[i3] * ch->dz;
    }

    st.site = &e_emlrtRSI;
    scalar_lin_step(&st, e_beta, ux);
    i3 = xi->size[1];
    b_beta = (int32_T)(2.0 + (real_T)i);
    if (!((b_beta >= 1) && (b_beta <= i3))) {
      emlrtDynamicBoundsCheckR2012b(b_beta, 1, i3, (emlrtBCInfo *)&b_emlrtBCI,
        sp);
    }

    st.site = &f_emlrtRSI;
    scalar_nl_step(&st, ux, xi->data[b_beta - 1]);

    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    i++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  emxFree_real_T(&e_beta);
  emxFree_real_T(&xi);
  emxInit_real_T(sp, &f_beta, 2, &emlrtRTEI, true);

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*                       LAST HALF DZ GVD                      % */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  i3 = f_beta->size[0] * f_beta->size[1];
  f_beta->size[0] = 1;
  f_beta->size[1] = beta->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)f_beta, i3, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = beta->size[0] * beta->size[1];
  for (i3 = 0; i3 < loop_ub; i3++) {
    f_beta->data[i3] = beta->data[i3] * Leff;
  }

  emxFree_real_T(&beta);
  st.site = &g_emlrtRSI;
  scalar_lin_step(&st, f_beta, ux);

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  i3 = sig->FIELDX->size[0] * sig->FIELDX->size[1];
  sig->FIELDX->size[0] = 1;
  sig->FIELDX->size[1] = ux->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)sig->FIELDX, i3, (int32_T)sizeof
                    (creal_T), &emlrtRTEI);
  loop_ub = ux->size[0] * ux->size[1];
  emxFree_real_T(&f_beta);
  for (i3 = 0; i3 < loop_ub; i3++) {
    sig->FIELDX->data[i3] = ux->data[i3];
  }

  emxFree_creal_T(&ux);

  /* set(sig,'FIELDX',ux); */
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void get_algo_sizes(const emlrtStack *sp, int32_T n1, boolean_T useRadix2,
                    int32_T *N2blue, int32_T *nRows)
{
  int32_T nn1m1;
  int32_T pmax;
  int32_T pmin;
  boolean_T exitg1;
  int32_T p;
  int32_T pow2p;
  *N2blue = 1;
  if (useRadix2) {
    *nRows = n1;
  } else {
    nn1m1 = (n1 + n1) - 1;
    pmax = 31;
    if (nn1m1 <= 1) {
      pmax = 0;
    } else {
      pmin = 0;
      exitg1 = false;
      while ((!exitg1) && (pmax - pmin > 1)) {
        p = (pmin + pmax) >> 1;
        pow2p = 1 << p;
        if (pow2p == nn1m1) {
          pmax = p;
          exitg1 = true;
        } else if (pow2p > nn1m1) {
          pmax = p;
        } else {
          pmin = p;
        }
      }
    }

    *N2blue = 1 << pmax;
    if (!(*N2blue <= (n1 << 2))) {
      emlrtErrorWithMessageIdR2012b(sp, &eb_emlrtRTEI,
        "Coder:builtins:AssertionFailed", 0);
    }

    *nRows = *N2blue;
  }
}

void r2br_r2dit_trig(const emlrtStack *sp, const emxArray_creal_T *x, int32_T
                     n1_unsigned, const emxArray_real_T *costab, const
                     emxArray_real_T *sintab, emxArray_creal_T *y)
{
  int32_T SZ1;
  int32_T istart;
  int32_T nRowsD2;
  int32_T nRowsD4;
  int32_T ix;
  int32_T ju;
  int32_T i;
  boolean_T tst;
  real_T temp_re;
  real_T temp_im;
  int32_T j;
  real_T twid_re;
  real_T twid_im;
  int32_T ihi;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &jb_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  SZ1 = x->size[0];
  istart = muIntScalarMin_sint32(SZ1, n1_unsigned) - 1;
  nRowsD2 = n1_unsigned / 2;
  nRowsD4 = nRowsD2 / 2;
  SZ1 = y->size[0];
  y->size[0] = n1_unsigned;
  emxEnsureCapacity(&st, (emxArray__common *)y, SZ1, (int32_T)sizeof(creal_T),
                    &n_emlrtRTEI);
  if (n1_unsigned > x->size[0]) {
    ix = y->size[0];
    SZ1 = y->size[0];
    y->size[0] = ix;
    emxEnsureCapacity(&st, (emxArray__common *)y, SZ1, (int32_T)sizeof(creal_T),
                      &o_emlrtRTEI);
    for (SZ1 = 0; SZ1 < ix; SZ1++) {
      y->data[SZ1].re = 0.0;
      y->data[SZ1].im = 0.0;
    }
  }

  ix = 0;
  ju = 0;
  SZ1 = 0;
  b_st.site = &kb_emlrtRSI;
  if ((!(1 > istart)) && (istart > 2147483646)) {
    c_st.site = &m_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }

  for (i = 1; i <= istart; i++) {
    y->data[SZ1] = x->data[ix];
    SZ1 = n1_unsigned;
    tst = true;
    while (tst) {
      SZ1 >>= 1;
      ju ^= SZ1;
      tst = ((ju & SZ1) == 0);
    }

    SZ1 = ju;
    ix++;
  }

  y->data[SZ1] = x->data[ix];
  if (n1_unsigned > 1) {
    for (i = 0; i <= n1_unsigned - 2; i += 2) {
      temp_re = y->data[i + 1].re;
      temp_im = y->data[i + 1].im;
      y->data[i + 1].re = y->data[i].re - y->data[i + 1].re;
      y->data[i + 1].im = y->data[i].im - y->data[i + 1].im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }
  }

  SZ1 = 2;
  ix = 4;
  ju = 1 + ((nRowsD4 - 1) << 2);
  while (nRowsD4 > 0) {
    for (i = 0; i < ju; i += ix) {
      temp_re = y->data[i + SZ1].re;
      temp_im = y->data[i + SZ1].im;
      y->data[i + SZ1].re = y->data[i].re - temp_re;
      y->data[i + SZ1].im = y->data[i].im - temp_im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }

    istart = 1;
    for (j = nRowsD4; j < nRowsD2; j += nRowsD4) {
      twid_re = costab->data[j];
      twid_im = sintab->data[j];
      i = istart;
      ihi = istart + ju;
      while (i < ihi) {
        temp_re = twid_re * y->data[i + SZ1].re - twid_im * y->data[i + SZ1].im;
        temp_im = twid_re * y->data[i + SZ1].im + twid_im * y->data[i + SZ1].re;
        y->data[i + SZ1].re = y->data[i].re - temp_re;
        y->data[i + SZ1].im = y->data[i].im - temp_im;
        y->data[i].re += temp_re;
        y->data[i].im += temp_im;
        i += ix;
      }

      istart++;
    }

    nRowsD4 /= 2;
    SZ1 = ix;
    ix <<= 1;
    ju -= SZ1;
  }
}

void r2br_r2dit_trig_impl(const emlrtStack *sp, const emxArray_creal_T *x,
  int32_T unsigned_nRows, const emxArray_real_T *costab, const emxArray_real_T
  *sintab, emxArray_creal_T *y)
{
  int32_T SZ1;
  int32_T istart;
  int32_T nRowsD2;
  int32_T nRowsD4;
  int32_T ix;
  int32_T ju;
  int32_T i;
  boolean_T tst;
  real_T temp_re;
  real_T temp_im;
  int32_T j;
  real_T twid_re;
  real_T twid_im;
  int32_T ihi;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  SZ1 = x->size[0];
  istart = muIntScalarMin_sint32(SZ1, unsigned_nRows) - 1;
  nRowsD2 = unsigned_nRows / 2;
  nRowsD4 = nRowsD2 / 2;
  SZ1 = y->size[0];
  y->size[0] = unsigned_nRows;
  emxEnsureCapacity(sp, (emxArray__common *)y, SZ1, (int32_T)sizeof(creal_T),
                    &p_emlrtRTEI);
  SZ1 = x->size[0];
  if (unsigned_nRows > SZ1) {
    ix = y->size[0];
    SZ1 = y->size[0];
    y->size[0] = ix;
    emxEnsureCapacity(sp, (emxArray__common *)y, SZ1, (int32_T)sizeof(creal_T),
                      &q_emlrtRTEI);
    for (SZ1 = 0; SZ1 < ix; SZ1++) {
      y->data[SZ1].re = 0.0;
      y->data[SZ1].im = 0.0;
    }
  }

  ix = 0;
  ju = 0;
  SZ1 = 0;
  st.site = &kb_emlrtRSI;
  if ((!(1 > istart)) && (istart > 2147483646)) {
    b_st.site = &m_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }

  for (i = 1; i <= istart; i++) {
    y->data[SZ1] = x->data[ix];
    SZ1 = unsigned_nRows;
    tst = true;
    while (tst) {
      SZ1 >>= 1;
      ju ^= SZ1;
      tst = ((ju & SZ1) == 0);
    }

    SZ1 = ju;
    ix++;
  }

  y->data[SZ1] = x->data[ix];
  if (unsigned_nRows > 1) {
    for (i = 0; i <= unsigned_nRows - 2; i += 2) {
      temp_re = y->data[i + 1].re;
      temp_im = y->data[i + 1].im;
      y->data[i + 1].re = y->data[i].re - y->data[i + 1].re;
      y->data[i + 1].im = y->data[i].im - y->data[i + 1].im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }
  }

  SZ1 = 2;
  ix = 4;
  ju = 1 + ((nRowsD4 - 1) << 2);
  while (nRowsD4 > 0) {
    for (i = 0; i < ju; i += ix) {
      temp_re = y->data[i + SZ1].re;
      temp_im = y->data[i + SZ1].im;
      y->data[i + SZ1].re = y->data[i].re - temp_re;
      y->data[i + SZ1].im = y->data[i].im - temp_im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }

    istart = 1;
    for (j = nRowsD4; j < nRowsD2; j += nRowsD4) {
      twid_re = costab->data[j];
      twid_im = sintab->data[j];
      i = istart;
      ihi = istart + ju;
      while (i < ihi) {
        temp_re = twid_re * y->data[i + SZ1].re - twid_im * y->data[i + SZ1].im;
        temp_im = twid_re * y->data[i + SZ1].im + twid_im * y->data[i + SZ1].re;
        y->data[i + SZ1].re = y->data[i].re - temp_re;
        y->data[i + SZ1].im = y->data[i].im - temp_im;
        y->data[i].re += temp_re;
        y->data[i].im += temp_im;
        i += ix;
      }

      istart++;
    }

    nRowsD4 /= 2;
    SZ1 = ix;
    ix <<= 1;
    ju -= SZ1;
  }
}

/* End of code generation (f_scalar_ssfm.cpp) */
