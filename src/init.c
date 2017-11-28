#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void aclambda(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void engen2008Ir(void *, void *, void *, void *, void *, void *, void *);
extern void engen2008r(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void fctdomain(void *, void *, void *, void *, void *, void *, void *);
extern void permutksel(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void randenfar(void *, void *, void *, void *, void *, void *);
extern void rankma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"aclambda",    (DL_FUNC) &aclambda,    11},
    {"engen2008Ir", (DL_FUNC) &engen2008Ir,  7},
    {"engen2008r",  (DL_FUNC) &engen2008r,  10},
    {"fctdomain",   (DL_FUNC) &fctdomain,    7},
    {"permutksel",  (DL_FUNC) &permutksel,  19},
    {"randenfar",   (DL_FUNC) &randenfar,    6},
    {"rankma",      (DL_FUNC) &rankma,      10},
    {NULL, NULL, 0}
};

void R_init_adehabitatHS(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
