/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * f_scalar_ssfm_types.h
 *
 * Code generation for function 'f_scalar_ssfm'
 *
 */

#ifndef F_SCALAR_SSFM_TYPES_H
#define F_SCALAR_SSFM_TYPES_H

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
typedef struct {
  real_T b2;
  real_T b3;
  real_T dz;
  real_T alphalin;
  real_T gamma;
  real_T nstep;
} Channel;

#ifndef struct_emxArray_creal_T
#define struct_emxArray_creal_T

struct emxArray_creal_T
{
  creal_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_creal_T*/

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

typedef struct {
  emxArray_creal_T *FIELDX;
  real_T SYMBOLRATE;
  emxArray_real_T *FN;
} Signal;

#ifndef struct_emxArray__common
#define struct_emxArray__common

struct emxArray__common
{
  void *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray__common*/
#endif

/* End of code generation (f_scalar_ssfm_types.h) */
