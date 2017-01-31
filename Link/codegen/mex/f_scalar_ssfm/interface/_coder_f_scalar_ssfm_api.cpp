/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_f_scalar_ssfm_api.cpp
 *
 * Code generation for function '_coder_f_scalar_ssfm_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "f_scalar_ssfm.h"
#include "_coder_f_scalar_ssfm_api.h"
#include "f_scalar_ssfm_emxutil.h"
#include "f_scalar_ssfm_data.h"

/* Variable Definitions */
static emlrtRTEInfo t_emlrtRTEI = { 1, /* lineNo */
  1,                                   /* colNo */
  "_coder_f_scalar_ssfm_api",          /* fName */
  ""                                   /* pName */
};

/* Function Declarations */
static Channel b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *Pavg,
  const char_T *identifier);
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *sig, const
  char_T *identifier, Signal *y);
static Channel emlrt_marshallIn(const emlrtStack *sp, const mxArray *ch, const
  char_T *identifier);
static const mxArray *emlrt_marshallOut(const emlrtStack *sp, const Signal u);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, Signal *y);
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_creal_T *y);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_creal_T *ret);
static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);

/* Function Definitions */
static Channel b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  Channel y;
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[6] = { "b2", "b3", "dz", "alphalin", "gamma",
    "nstep" };

  static const int32_T dims = 0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 6, fieldNames, 0U, &dims);
  thisId.fIdentifier = "b2";
  y.b2 = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "b2")),
    &thisId);
  thisId.fIdentifier = "b3";
  y.b3 = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "b3")),
    &thisId);
  thisId.fIdentifier = "dz";
  y.dz = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "dz")),
    &thisId);
  thisId.fIdentifier = "alphalin";
  y.alphalin = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "alphalin")), &thisId);
  thisId.fIdentifier = "gamma";
  y.gamma = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "gamma")), &thisId);
  thisId.fIdentifier = "nstep";
  y.nstep = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "nstep")), &thisId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *Pavg,
  const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = c_emlrt_marshallIn(sp, emlrtAlias(Pavg), &thisId);
  emlrtDestroyArray(&Pavg);
  return y;
}

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *sig, const
  char_T *identifier, Signal *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  f_emlrt_marshallIn(sp, emlrtAlias(sig), &thisId, y);
  emlrtDestroyArray(&sig);
}

static Channel emlrt_marshallIn(const emlrtStack *sp, const mxArray *ch, const
  char_T *identifier)
{
  Channel y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(ch), &thisId);
  emlrtDestroyArray(&ch);
  return y;
}

static const mxArray *emlrt_marshallOut(const emlrtStack *sp, const Signal u)
{
  const mxArray *y;
  emxArray_creal_T *b_u;
  emxArray_real_T *c_u;
  int32_T i0;
  int32_T loop_ub;
  const mxArray *b_y;
  const mxArray *m0;
  const mxArray *c_y;
  const mxArray *d_y;
  real_T *pData;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_creal_T1(sp, &b_u, 2, (emlrtRTEInfo *)NULL, true);
  emxInit_real_T(sp, &c_u, 2, (emlrtRTEInfo *)NULL, true);
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 0, NULL));
  i0 = b_u->size[0] * b_u->size[1];
  b_u->size[0] = 1;
  b_u->size[1] = u.FIELDX->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)b_u, i0, (int32_T)sizeof(creal_T),
                    (emlrtRTEInfo *)NULL);
  loop_ub = u.FIELDX->size[0] * u.FIELDX->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_u->data[i0] = u.FIELDX->data[i0];
  }

  b_y = NULL;
  m0 = emlrtCreateNumericArray(2, *(int32_T (*)[2])b_u->size, mxDOUBLE_CLASS,
    mxCOMPLEX);
  emlrtExportNumericArrayR2013b(sp, m0, (void *)&b_u->data[0], 8);
  emlrtAssign(&b_y, m0);
  emlrtAddField(y, b_y, "FIELDX", 0);
  c_y = NULL;
  m0 = emlrtCreateDoubleScalar(u.SYMBOLRATE);
  emlrtAssign(&c_y, m0);
  emlrtAddField(y, c_y, "SYMBOLRATE", 0);
  i0 = c_u->size[0] * c_u->size[1];
  c_u->size[0] = 1;
  c_u->size[1] = u.FN->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)c_u, i0, (int32_T)sizeof(real_T),
                    (emlrtRTEInfo *)NULL);
  loop_ub = u.FN->size[0] * u.FN->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    c_u->data[i0] = u.FN->data[i0];
  }

  d_y = NULL;
  m0 = emlrtCreateNumericArray(2, *(int32_T (*)[2])c_u->size, mxDOUBLE_CLASS,
    mxREAL);
  pData = (real_T *)mxGetPr(m0);
  i0 = 0;
  for (loop_ub = 0; loop_ub < c_u->size[1]; loop_ub++) {
    pData[i0] = c_u->data[c_u->size[0] * loop_ub];
    i0++;
  }

  emlrtAssign(&d_y, m0);
  emlrtAddField(y, d_y, "FN", 0);
  emxFree_real_T(&c_u);
  emxFree_creal_T(&b_u);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
  return y;
}

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, Signal *y)
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[3] = { "FIELDX", "SYMBOLRATE", "FN" };

  static const int32_T dims = 0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 3, fieldNames, 0U, &dims);
  thisId.fIdentifier = "FIELDX";
  g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "FIELDX")),
                     &thisId, y->FIELDX);
  thisId.fIdentifier = "SYMBOLRATE";
  y->SYMBOLRATE = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0,
    "SYMBOLRATE")), &thisId);
  thisId.fIdentifier = "FN";
  h_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "FN")),
                     &thisId, y->FN);
  emlrtDestroyArray(&u);
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_creal_T *y)
{
  j_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  k_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_creal_T *ret)
{
  static const int32_T dims[2] = { 1, -1 };

  boolean_T bv0[2] = { false, true };

  int32_T iv0[2];
  int32_T i1;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", true, 2U, dims, &bv0[0],
    iv0);
  i1 = ret->size[0] * ret->size[1];
  ret->size[0] = iv0[0];
  ret->size[1] = iv0[1];
  emxEnsureCapacity(sp, (emxArray__common *)ret, i1, (int32_T)sizeof(creal_T),
                    (emlrtRTEInfo *)NULL);
  emlrtImportArrayR2015b(sp, src, ret->data, 8, true);
  emlrtDestroyArray(&src);
}

static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  static const int32_T dims[2] = { 1, -1 };

  boolean_T bv1[2] = { false, true };

  int32_T iv1[2];
  int32_T i2;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims, &bv1[0],
    iv1);
  i2 = ret->size[0] * ret->size[1];
  ret->size[0] = iv1[0];
  ret->size[1] = iv1[1];
  emxEnsureCapacity(sp, (emxArray__common *)ret, i2, (int32_T)sizeof(real_T),
                    (emlrtRTEInfo *)NULL);
  emlrtImportArrayR2015b(sp, src, ret->data, 8, false);
  emlrtDestroyArray(&src);
}

void f_scalar_ssfm_api(const mxArray * const prhs[3], const mxArray *plhs[1])
{
  Signal sig;
  Channel ch;
  real_T Pavg;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInitStruct_Signal(&st, &sig, &t_emlrtRTEI, true);

  /* Marshall function inputs */
  ch = emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[0]), "ch");
  Pavg = d_emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[1]), "Pavg");
  e_emlrt_marshallIn(&st, emlrtAliasP((const mxArray *)prhs[2]), "sig", &sig);

  /* Invoke the target function */
  f_scalar_ssfm(&st, &ch, Pavg, &sig);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(&st, sig);
  emxFreeStruct_Signal(&sig);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_f_scalar_ssfm_api.cpp) */
