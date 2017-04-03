
#include<R.h>
#include<Rmath.h>
#include <Rinternals.h>

#include "fcros.h"

/* this function is called in fcros() and in pfco() */
SEXP rmat(SEXP fvectC, SEXP nC, SEXP m1C, SEXP m2C) {

     int i, j, k, nn, mm1, mm2, n;
     double deno, nume, *fvect, *rvect, *fc;
     SEXP rvectC, FCC, results;
     static const char *resultNames[]={"rvectC", "FCC", ""};

     fvect = REAL(fvectC);
     nn = asInteger(nC);
     mm1 = asInteger(m1C);
     mm2 = asInteger(m2C);
     n = nn*mm1*mm2;

    PROTECT(results = mkNamed(VECSXP, resultNames));
    rvectC = SET_VECTOR_ELT(results, 0, allocVector(REALSXP, n));
    rvect = REAL(rvectC);
    FCC = SET_VECTOR_ELT(results, 1, allocVector(REALSXP, nn));
    fc = REAL(FCC);

    /* calculation of rvectC  */
     for (i=0; i<mm1; i++) {
         for (j=0; j<mm2; j++) {
             for (k=0; k<nn; k++) {
                 rvect[i*mm2*nn+j*nn+k] = fvect[nn*mm1+j*nn+k] - fvect[i*nn+k];
             }
             k++;
         }
     }

     /* calculation of FCC  */
     for (i=0; i<nn; i++) {
         deno = nume = 0.0;
         for (j=0; j<mm1; j++) deno += pow(2, fvect[i+j*nn]);
         for (j=0; j<mm2; j++) nume += pow(2, fvect[nn*mm1+i+j*nn]);
         if (deno != 0.0) 
            fc[i] = (nume*mm1)/(deno*mm2);
         else fc[i] = 1000;
     }

     UNPROTECT(1);
     return results;
}/* end of function rmat() */

/* This function is used in the functions rmat2() and fc2() */
void qSort(double *table, int p, int r) {
     int i, j, itmp;
     double x, tmp;

     i=p;   j=r;     itmp = (i+j+1) / 2;
     x = table[itmp];
     do	{
        while (table[i] < x) 	i++;
        while (table[j] > x) 	j--;
        if (i<=j)       {
           tmp = table[i];    
           table[i] = table[j];
           table[j] = tmp;
           i++; j--;
        }
     } while (i<=j);
     if (p<j)	qSort(table,p,j);
     if (r>i)	qSort(table,i,r);
}/* end of function qSort()  */

SEXP rmat2(SEXP rvectC, SEXP nC, SEXP mC, SEXP idxC, SEXP m2C) {

     int i, j, n, m, m2, *idx, nm2;
     double *vect, *rvect, *rvect2;
     SEXP rvect2C;

     rvect = REAL(rvectC);
     n = asInteger(nC);
     m = asInteger(mC);
     idx = INTEGER(idxC);
     m2 = asInteger(m2C);
     nm2 = n*m2;
     vect = (double *)R_alloc(m, sizeof(double));

     PROTECT(rvect2C = allocVector(REALSXP, nm2));
     rvect2 = REAL(rvect2C);

     /* calculation of rvect2C */
     for (i=0; i<n; i++) {
         for (j=0; j<m; j++)  vect[j] = rvect[j*n+i];
         qSort(vect, 0, m-1);
         for (j=0; j<m2; j++) rvect2[j*n+i] = vect[idx[j]-1];
     }
     UNPROTECT(1);
     return rvect2C;
}/* end of function rmat2C object  */

SEXP tproba(SEXP moyC, SEXP stdC, SEXP nC, SEXP dlC, SEXP emC) {

     int i, n, dl;
     double tval, em, *moy, *std, *proba; 
     SEXP probaC;

     moy = REAL(moyC);
     std = REAL(stdC);
     n = asInteger(nC);
     dl = asInteger(dlC);
     em = asReal(emC);
     PROTECT(probaC = allocVector(REALSXP, n));
     proba = REAL(probaC);

     /* calculation of proba */
     for (i=0; i<n; i++) {
         tval = fabs(sqrt(dl+1.0)*((moy[i]-em)/std[i]));
         GetRNGstate();
         proba[i] = 2.0*(pt(tval, dl, 0, 0));
         PutRNGstate();
     }
     UNPROTECT(1);
     return probaC;
}/* end of function tproba() */

SEXP fc2(SEXP rvectC, SEXP nC, SEXP mC, SEXP idxC, SEXP m2C) {

     int i, j, n, m, m2, *idx;
     double *vect, tmp, *rvect, *fc2;
     SEXP fc2C;

     rvect = REAL(rvectC);
     n = asInteger(nC);
     m = asInteger(mC);
     idx = INTEGER(idxC);
     m2 = asInteger(m2C);
     PROTECT(fc2C = allocVector(REALSXP, n));
     fc2 = REAL(fc2C);
     vect = (double *)R_alloc(m, sizeof(double));

     /* calculation of fc2C  */
     for (i=0; i<n; i++) {
         for (j=0; j<m; j++)  vect[j] = rvect[j*n+i];
         qSort(vect, 0, m-1);
         tmp = 0.0;
         for (j=0; j<m2; j++) tmp += pow(2.0, vect[idx[j]-1]);
         fc2[i] = tmp / (double)m2;
     }
     UNPROTECT(1);
     return fc2C;
}/* end of function fc2() */


SEXP moyStd(SEXP rvectC, SEXP nC, SEXP mC) {

     int i, j, n, m;
     double tmp, *vect, *rvect, *moy2, *std2;
     SEXP moyC, stdC, results;
     static const char *resultNames[]={"moyC", "stdC", ""};

     rvect = REAL(rvectC);
     n = asInteger(nC);
     m = asInteger(mC);
     vect = (double *)R_alloc(m, sizeof(double));

     PROTECT(results = mkNamed(VECSXP, resultNames));
     moyC = SET_VECTOR_ELT(results, 0, allocVector(REALSXP, n));
     moy2 = REAL(moyC);
     stdC = SET_VECTOR_ELT(results, 1, allocVector(REALSXP, n));
     std2 = REAL(stdC);

      /* calculation of moyC and stdC */
     for (i=0; i<n; i++) {
         tmp = 0.0;
         for (j=0; j<m; j++)  {
             vect[j] = rvect[j*n+i];
             tmp += vect[j];
         }
         moy2[i] = tmp/(double)m;
         tmp = 0.0;
         for (j=0; j<m; j++) tmp += (vect[j]-moy2[i]) * (vect[j]-moy2[i]);
         std2[i] = sqrt(tmp/(double)(m-1));
     }
     UNPROTECT(1);
     return results;
}/* end of function moyStd()  */

SEXP merge(SEXP nSegC, SEXP segIdSC, SEXP segIdEC, SEXP segLBC,
                SEXP segUBC, SEXP segValC, SEXP segProbaC, SEXP fcallC,
                SEXP L2RC, SEXP ndC, SEXP dmC, SEXP sigmaC) {

     int i, j, k, i1, i2, j1, j2, nb, v1, v2, seg_d, seg_v, ns, nbSeg;
     int *nSeg, *segIdS2, *segIdS, *segIdE2, *segIdE, *fcall;
     double dmean, sig, d1, dnb, ndd, xbar, stat;
     double  *segVal2, *segVal, *segProba2, *segProba, *L2R,  *segLB2,
             *segLB, *segUB2, *segUB;
     SEXP r_nSeg, r_segLB, r_segUB, r_segVal, r_segProba, r_segIdS, 
          r_segIdE, results;
     static const char *resultNames[]={"nSeg", "segLB", "segUB", "segVal",
            "segProba", "segIdS", "segIdE", ""};

     nbSeg = asInteger(nSegC);
     segIdS2 = INTEGER(segIdSC);
     segIdE2 = INTEGER(segIdEC);
     segLB2 = REAL(segLBC);
     segUB2 = REAL(segUBC);
     segVal2 = REAL(segValC);
     segProba2 = REAL(segProbaC);
     fcall = INTEGER(fcallC);
     L2R = REAL(L2RC);
     ndd = asInteger(ndC);
     dmean = asReal(dmC);
     sig = asReal(sigmaC);

     PROTECT(results = mkNamed(VECSXP, resultNames));
     r_nSeg = SET_VECTOR_ELT(results, 0, allocVector(INTSXP, 1));
     nSeg = INTEGER(r_nSeg);
     r_segLB = SET_VECTOR_ELT(results, 1, allocVector(REALSXP, nbSeg));
     segLB = REAL(r_segLB);
     r_segUB = SET_VECTOR_ELT(results, 2, allocVector(REALSXP, nbSeg));
     segUB = REAL(r_segUB);
     r_segVal = SET_VECTOR_ELT(results, 3, allocVector(REALSXP, nbSeg));
     segVal = REAL(r_segVal);
     r_segProba = SET_VECTOR_ELT(results, 4, allocVector(REALSXP, nbSeg));
     segProba = REAL(r_segProba);
     r_segIdS = SET_VECTOR_ELT(results, 5, allocVector(INTSXP, nbSeg));
     segIdS = INTEGER(r_segIdS);
     r_segIdE = SET_VECTOR_ELT(results, 6, allocVector(INTSXP, nbSeg));
     segIdE = INTEGER(r_segIdE);

     for (i=0; i<nbSeg; i++) {
         segLB[i] = segLB2[i];
         segUB[i] = segUB2[i];
         segVal[i] = segVal2[i];
         segProba[i] = segProba2[i];
         segIdS[i] = segIdS2[i];
         segIdE[i] = segIdE2[i];
     }

     /* start merging of segments data  */
     j = 0;
     while (j < nbSeg-1) {
           j1 = segIdE[j];
           j2 = segIdS[j+1];
           nb = j2-j1+1;
           dnb = nb*dmean; /* average dist between segment j and segment j+1  */
           d1 = dnb;
           seg_d = seg_v = v1 = v2 = 0;
           if (nb < ndd) d1 = segLB[j+1] - segUB[j];
           if (d1 < dnb) seg_d = 1;
           if ((fcall[j1-1] > 0) && (fcall[j2-1] > 0)) v1 = 1;  /* amplification  */
           if ((fcall[j1-1] < 0) && (fcall[j2-1] < 0)) v2 = 1;  /* deletion       */
           if (v1 || v2) seg_v = 1;

           /* check for merging  */
           if (seg_d && seg_v) { /* merging and stay at the same segment  */
              segUB[j] = segUB[j+1];
              segIdE[j] = segIdE[j+1];
              ns = segIdE[j]-segIdS[j]+1;
              xbar = 0; i1 = segIdS[j]; i2 = segIdE[j];
              for (i=i1; i<=i2; i++) xbar += L2R[i-1];
              xbar /= (double) ns;
              segVal[j] = xbar;
              stat = sqrt((double) ns)*xbar/sig;
              GetRNGstate();
              segProba[j] = pnorm(stat, 0, 1, 1, 0);
              PutRNGstate();
              for (k=j+1; k<nbSeg-1; k++) {
                  segLB[k] = segLB[k+1];
                  segUB[k] = segUB[k+1];
                  segVal[k] = segVal[k+1];
                  segProba[k] = segProba[k+1];
                  segIdS[k] = segIdS[k+1];
                  segIdE[k] = segIdE[k+1];
              }
              nbSeg--;
           } else { /* no merging, go to next segment  */
                  j++;
           }
     }
     UNPROTECT(1);
     (*nSeg) = nbSeg;
     return results;
}/* end of function merge() */
