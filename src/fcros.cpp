#include<R.h>
#include<Rmath.h>

extern "C" {
void rmat(double *fvectC, int *nC, int *m1C, int *m2C,
                 double *rvectC, double *FCC) {
     int i, j, k, nn, mm1, mm2;
     double deno, nume;

     nn = (*nC);
     mm1 = (*m1C);
     mm2 = (*m2C);

     // calculation of rvectC
     for (i=0; i<mm1; i++) {
         for (j=0; j<mm2; j++) {
             for (k=0; k<nn; k++) {
                 rvectC[i*mm2*nn+j*nn+k] = fvectC[nn*mm1+j*nn+k] - fvectC[i*nn+k];
             }
             k++;
         }
     }
     // calculation of FCC
     for (i=0; i<nn; i++) {
         deno = nume = 0.0;
         for (j=0; j<mm1; j++) deno += pow(2, fvectC[i+j*nn]);
         for (j=0; j<mm2; j++) nume += pow(2, fvectC[nn*mm1+i+j*nn]);
         if (deno != 0.0) 
            FCC[i] = (nume*mm1)/(deno*mm2);
         else FCC[i] = 1000;
     }
}// end of function rmat()

void merge(int *nSeg, int *segIdS, int *segIdE, double *segLB, double *segUB,
           double *segVal, double *segProba, int *fcall, double *L2R, int *nd,
           double *dm, double *sigma) {

     int i, j, k, nbSeg, i1, i2, j1, j2, nb, v1, v2, seg_d, seg_v, ns;
     double dmean, sig, d1, dnb, ndd, xbar, stat;

     nbSeg = (*nSeg);
     dmean = (*dm);
     sig = (*sigma);
     ndd = (*nd);

     // start merging of segments data
     j = 0;
     while (j < nbSeg-1) {
           j1 = segIdE[j];
           j2 = segIdS[j+1];
           nb = j2-j1+1;
           dnb = nb*dmean; // average dist between segment j and segment j+1
           d1 = dnb;
           seg_d = seg_v = v1 = v2 = 0;
           if (nb < ndd) d1 = segLB[j+1] - segUB[j];
           if (d1 < dnb) seg_d = 1;
           if ((fcall[j1-1] > 0) && (fcall[j2-1] > 0)) v1 = 1;  // amplification
           if ((fcall[j1-1] < 0) && (fcall[j2-1] < 0)) v2 = 1;  // deletion
           if (v1 || v2) seg_v = 1;
           
           // check for merging
           if (seg_d && seg_v) { // merging and stay at the same segment
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
           } else { // no merging, go to next segment
                  j++;
           }
     }
     (*nSeg) = nbSeg;
}// end of function merge()

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
}// end of function qSort()

void rmat2(double *rvectC, int *nC, int *mC,
                  int *idxC, int *m2C, double *rvect2C) {
     int i, j, n, m, m2;
     double *vect;

     n = (*nC);
     m = (*mC);
     m2 = (*m2C);
     vect = new double[m];

     // calculation of rvect2C
     for (i=0; i<n; i++) {
         for (j=0; j<m; j++)  vect[j] = rvectC[j*n+i];
         qSort(vect, 0, m-1);
         for (j=0; j<m2; j++) rvect2C[j*n+i] = vect[idxC[j]-1];
     }
     delete[] vect;
}// end of function rmat2()

void moyStd(double *rvectC, int *nC, int *mC, double *moyC, double *stdC) {
     int i, j, n, m;
     double tmp, *vect;

     n = (*nC);
     m = (*mC);
     vect = new double[m];

     // calculation of moyC and stdC
     for (i=0; i<n; i++) {
         tmp = 0.0;
         for (j=0; j<m; j++)  {
             vect[j] = rvectC[j*n+i];
             tmp += vect[j];
         }
         moyC[i] = tmp/(double)m;
         tmp = 0.0;
         for (j=0; j<m; j++) tmp += (vect[j]-moyC[i]) * (vect[j]-moyC[i]);
         stdC[i] = sqrt(tmp/(double)(m-1));
     }
     delete[] vect;
}// end of function moyStd()

void tproba(double *moyC, double *stdC, int *nC, int *dlC, double *emC, 
                   double *probaC) {
     int i, n, dl;
     double tval, em;

     n = (*nC);
     em = (*emC);
     dl = (*dlC);

     // calculation of probaC
     for (i=0; i<n; i++) {
         tval = fabs(sqrt(dl+1.0)*((moyC[i]-em)/stdC[i]));
         GetRNGstate();
         probaC[i] = 2.0*(pt(tval, dl, 0, 0));
         PutRNGstate();
     }
}// end of function tproba()

void fc2(double *rvectC, int *nC, int *mC,
                  int *idxC, int *m2C, double *fc2C) {
     int i, j, n, m, m2;
     double *vect, tmp;

     n = (*nC);
     m = (*mC);
     m2 = (*m2C);
     vect = new double[m];

     // calculation of fc2C
     for (i=0; i<n; i++) {
         for (j=0; j<m; j++)  vect[j] = rvectC[j*n+i];
         qSort(vect, 0, m-1);
         tmp = 0.0;
         for (j=0; j<m2; j++) tmp += pow(2.0, vect[idxC[j]-1]);
         fc2C[i] = tmp / (double)m2;
     }
     delete[] vect;
}// end of function fc2()

}  // extern "C"
