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
}

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
}
}  // extern "C"
