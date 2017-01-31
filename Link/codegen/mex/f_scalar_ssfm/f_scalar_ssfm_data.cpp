/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * f_scalar_ssfm_data.cpp
 *
 * Code generation for function 'f_scalar_ssfm_data'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "f_scalar_ssfm.h"
#include "f_scalar_ssfm_data.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;
omp_lock_t emlrtLockGlobal;
omp_nest_lock_t emlrtNestLockGlobal;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131435U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "f_scalar_ssfm",                     /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

emlrtRSInfo h_emlrtRSI = { 49,         /* lineNo */
  "power",                             /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\ops\\power.m"/* pathName */
};

emlrtRSInfo i_emlrtRSI = { 58,         /* lineNo */
  "power",                             /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\ops\\power.m"/* pathName */
};

emlrtRSInfo j_emlrtRSI = { 61,         /* lineNo */
  "power",                             /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\ops\\power.m"/* pathName */
};

emlrtRSInfo k_emlrtRSI = { 73,         /* lineNo */
  "applyScalarFunction",               /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\eml\\+coder\\+internal\\applyScalarFunction.m"/* pathName */
};

emlrtRSInfo l_emlrtRSI = { 132,        /* lineNo */
  "applyScalarFunction",               /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\eml\\+coder\\+internal\\applyScalarFunction.m"/* pathName */
};

emlrtRSInfo m_emlrtRSI = { 20,         /* lineNo */
  "eml_int_forloop_overflow_check",    /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\eml\\eml_int_forloop_overflow_check.m"/* pathName */
};

emlrtRSInfo n_emlrtRSI = { 9,          /* lineNo */
  "exp",                               /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\elfun\\exp.m"/* pathName */
};

emlrtRSInfo o_emlrtRSI = { 24,         /* lineNo */
  "applyScalarFunctionInPlace",        /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\eml\\+coder\\+internal\\applyScalarFunctionInPlace.m"/* pathName */
};

emlrtRSInfo s_emlrtRSI = { 45,         /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo t_emlrtRSI = { 73,         /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo u_emlrtRSI = { 74,         /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo v_emlrtRSI = { 76,         /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo w_emlrtRSI = { 78,         /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo x_emlrtRSI = { 124,        /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo y_emlrtRSI = { 128,        /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo ab_emlrtRSI = { 135,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo bb_emlrtRSI = { 482,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo cb_emlrtRSI = { 485,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo db_emlrtRSI = { 511,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo eb_emlrtRSI = { 514,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo fb_emlrtRSI = { 517,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo gb_emlrtRSI = { 521,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo hb_emlrtRSI = { 536,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo ib_emlrtRSI = { 540,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo lb_emlrtRSI = { 160,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo mb_emlrtRSI = { 163,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo nb_emlrtRSI = { 9,         /* lineNo */
  "bluestein_setup",                   /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\bluestein_setup.m"/* pathName */
};

emlrtRSInfo pb_emlrtRSI = { 214,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo qb_emlrtRSI = { 219,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo rb_emlrtRSI = { 229,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo sb_emlrtRSI = { 232,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo tb_emlrtRSI = { 237,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo wb_emlrtRSI = { 527,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRSInfo xb_emlrtRSI = { 531,       /* lineNo */
  "eml_fft",                           /* fcnName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pathName */
};

emlrtRTEInfo c_emlrtRTEI = { 16,       /* lineNo */
  9,                                   /* colNo */
  "scalexpAlloc",                      /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\eml\\+coder\\+internal\\scalexpAlloc.m"/* pName */
};

emlrtRTEInfo e_emlrtRTEI = { 67,       /* lineNo */
  9,                                   /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

emlrtRTEInfo g_emlrtRTEI = { 208,      /* lineNo */
  9,                                   /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

emlrtRTEInfo h_emlrtRTEI = { 45,       /* lineNo */
  9,                                   /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

emlrtRTEInfo i_emlrtRTEI = { 160,      /* lineNo */
  1,                                   /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

emlrtRTEInfo j_emlrtRTEI = { 228,      /* lineNo */
  5,                                   /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

emlrtRTEInfo k_emlrtRTEI = { 231,      /* lineNo */
  5,                                   /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

emlrtRTEInfo l_emlrtRTEI = { 115,      /* lineNo */
  40,                                  /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

emlrtRTEInfo m_emlrtRTEI = { 124,      /* lineNo */
  5,                                   /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

emlrtRTEInfo ab_emlrtRTEI = { 17,      /* lineNo */
  19,                                  /* colNo */
  "scalexpAlloc",                      /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\eml\\+coder\\+internal\\scalexpAlloc.m"/* pName */
};

emlrtRTEInfo db_emlrtRTEI = { 18,      /* lineNo */
  19,                                  /* colNo */
  "eml_fft",                           /* fName */
  "C:\\Program Files\\MATLAB\\R2016b\\toolbox\\eml\\lib\\matlab\\datafun\\private\\eml_fft.m"/* pName */
};

/* End of code generation (f_scalar_ssfm_data.cpp) */
