@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2016b
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2016b\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=f_scalar_ssfm_mex
set MEX_NAME=f_scalar_ssfm_mex
set MEX_EXT=.mexw64
call setEnv.bat
echo # Make settings for f_scalar_ssfm > f_scalar_ssfm_mex.mki
echo COMPILER=%COMPILER%>> f_scalar_ssfm_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> f_scalar_ssfm_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> f_scalar_ssfm_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> f_scalar_ssfm_mex.mki
echo LINKER=%LINKER%>> f_scalar_ssfm_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> f_scalar_ssfm_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> f_scalar_ssfm_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> f_scalar_ssfm_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> f_scalar_ssfm_mex.mki
echo BORLAND=%BORLAND%>> f_scalar_ssfm_mex.mki
echo OMPFLAGS=/openmp >> f_scalar_ssfm_mex.mki
echo OMPLINKFLAGS=/nodefaultlib:vcomp /LIBPATH:"C:\PROGRA~1\MATLAB\R2016b\bin\win64" >> f_scalar_ssfm_mex.mki
echo EMC_COMPILER=msvcpp140>> f_scalar_ssfm_mex.mki
echo EMC_CONFIG=optim>> f_scalar_ssfm_mex.mki
"C:\Program Files\MATLAB\R2016b\bin\win64\gmake" -B -f f_scalar_ssfm_mex.mk
