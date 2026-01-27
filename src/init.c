#include <R.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "fcros.h"

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

const static R_CallMethodDef R_CallDef[] = {
      CALLDEF(rmat, 4),
      CALLDEF(rmat2, 5),
      CALLDEF(moyStd, 3),
      CALLDEF(tproba, 5),
      CALLDEF(fc2, 5),
      CALLDEF(merge, 12),
      {NULL, NULL, 0}
};

void R_init_fcros(DllInfo *dll)  {
     R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);   
     R_useDynamicSymbols(dll, FALSE);
     R_forceSymbols(dll, TRUE);
}
