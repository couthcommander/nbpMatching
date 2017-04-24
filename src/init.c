#include <R_ext/RS.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern void F77_NAME(mwrap)(int *n, int *wt, int *nmatch, int *prcn);

static const R_FortranMethodDef FortranEntries[] = {
    {"mwrap", (DL_FUNC) &F77_NAME(mwrap), 4},
    {NULL, NULL, 0}
};

void R_init_nbpMatching(DllInfo *dll) {
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
