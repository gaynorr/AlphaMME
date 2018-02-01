#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _AlphaMME_calcD(SEXP);
extern SEXP _AlphaMME_calcG(SEXP);
extern SEXP _AlphaMME_calcGIbs(SEXP);
extern SEXP _AlphaMME_fastDist(SEXP);
extern SEXP _AlphaMME_fastPairDist(SEXP, SEXP);
extern SEXP _AlphaMME_gaussKernel(SEXP, SEXP);
extern SEXP _AlphaMME_readMat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaMME_solveAniModel(SEXP, SEXP, SEXP);
extern SEXP _AlphaMME_solveMKM(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaMME_solveMVM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaMME_solveRRBLUP(SEXP, SEXP, SEXP);
extern SEXP _AlphaMME_solveRRBLUPMK(SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaMME_solveRRBLUPMV(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaMME_solveUVM(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_AlphaMME_calcD",         (DL_FUNC) &_AlphaMME_calcD,         1},
    {"_AlphaMME_calcG",         (DL_FUNC) &_AlphaMME_calcG,         1},
    {"_AlphaMME_calcGIbs",      (DL_FUNC) &_AlphaMME_calcGIbs,      1},
    {"_AlphaMME_fastDist",      (DL_FUNC) &_AlphaMME_fastDist,      1},
    {"_AlphaMME_fastPairDist",  (DL_FUNC) &_AlphaMME_fastPairDist,  2},
    {"_AlphaMME_gaussKernel",   (DL_FUNC) &_AlphaMME_gaussKernel,   2},
    {"_AlphaMME_readMat",       (DL_FUNC) &_AlphaMME_readMat,       6},
    {"_AlphaMME_solveAniModel", (DL_FUNC) &_AlphaMME_solveAniModel, 3},
    {"_AlphaMME_solveMKM",      (DL_FUNC) &_AlphaMME_solveMKM,      5},
    {"_AlphaMME_solveMVM",      (DL_FUNC) &_AlphaMME_solveMVM,      6},
    {"_AlphaMME_solveRRBLUP",   (DL_FUNC) &_AlphaMME_solveRRBLUP,   3},
    {"_AlphaMME_solveRRBLUPMK", (DL_FUNC) &_AlphaMME_solveRRBLUPMK, 4},
    {"_AlphaMME_solveRRBLUPMV", (DL_FUNC) &_AlphaMME_solveRRBLUPMV, 5},
    {"_AlphaMME_solveUVM",      (DL_FUNC) &_AlphaMME_solveUVM,      4},
    {NULL, NULL, 0}
};

void R_init_AlphaMME(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
