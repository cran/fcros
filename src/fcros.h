#include <Rinternals.h> // for SEXP

extern SEXP rmat(SEXP fvectC, SEXP nC, SEXP m1C, SEXP m2C);
extern SEXP rmat2(SEXP rvectC, SEXP nC, SEXP mC, SEXP idxC, SEXP m2C);
extern SEXP tproba(SEXP moyC, SEXP stdC, SEXP nC, SEXP dlC, SEXP emC);
extern SEXP fc2(SEXP rvectC, SEXP nC, SEXP mC, SEXP idxC, SEXP m2C);
extern SEXP moyStd(SEXP rvectC, SEXP nC, SEXP mC);
extern SEXP merge(SEXP nSegC, SEXP segIdSC, SEXP segIdEC, SEXP segLBC, 
       SEXP segUBC, SEXP segValC, SEXP segProbaC, SEXP fcallC, SEXP L2RC, 
       SEXP ndC,  SEXP dmC, SEXP sigmaC);
