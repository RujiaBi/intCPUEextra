#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" void R_init_intCPUEextra(DllInfo *dll) {
  R_registerRoutines(dll, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
