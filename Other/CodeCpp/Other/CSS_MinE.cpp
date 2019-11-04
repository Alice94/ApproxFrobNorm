#include <math.h>
#include <time.h>
#include "Poly.cpp"
#include "BidiagonalizeCSS.cpp"
using namespace std;

// Selects the t-th column of B, knowing that the current iteration is t. Returns the chosen index (1-based!)
long int SelectNextColumn_MinE(long int* m, long int* n, double* B, long int* k, long int* t);

long int SelectNextColumn_MinE(long int* m, long int* n, double* B, long int* k, long int* t) {
  
  // Initialization for time checks
  clock_t t1, t2;
  clock_t clockInitialSVD = 0, clockBidiagonalize = 0, clockPoly = 0, clockBidiagonalSVD = 0;
   
  // Inizialization for SVD
  double* s = (double*)malloc(min(*n,*m)*sizeof(double));
  long int lwork_local = 7*min(*m,*n)+4*max(*n,*m);
  double* work_local = (double*)malloc(lwork_local*sizeof(double));
  long int info;
  long int ldvt = 1;
  char uplo_SVD = 'U';
  long int p = min(*m, *n);
  
  // Inizialization for QR factorization and multiplication
  char uplo_cpy = 'A';
  char side = 'L';
  char trans = 'T';
  
  // Initialization for Bidiagonalize
  double* e = (double*)malloc((p-1)*sizeof(double));
  double* d = (double*)malloc(p*sizeof(double));
   
  // Initialization for Poly
  double* res_poly = (double*)malloc(((*m)+2)*sizeof(double));
  double* res_poly1 = (double*)malloc(((*m)+2)*sizeof(double));
  double* res_poly2 = (double*)malloc(((*m)+2)*sizeof(double));
  
  // General useful stuff
  long int one = 1;
  long int zero = 0;
  char notrans = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  
  double* C = (double*)malloc((*m)*(*n)*sizeof(double));
  double* U = (double*)malloc((*m)*(p)*sizeof(double));
 
  double* v;
  v = (double*)malloc((*m)*sizeof(double));
  double* w = (double*)malloc((*m)*sizeof(double));
  
  long int bestCol = 0;
  double MinRatio = -1;
    
  // Do the SVD of B, result is saved in d and C
  char jobu = 'S';
  char jobvt = 'N';
  dlacpy(&uplo_cpy, m, n, B, m, C, m);
  t1 = clock();
  dgesvd(&jobu, &jobvt, m, n, C, m, s, U, m, NULL, &one, work_local, &lwork_local, &info);
  t2 = clock();
  clockInitialSVD = (t2 - t1);
  
  // Test one column at a time
  for (long int i=0; i<*n; i++) {
    // Check the norm
    if (dnrm2_(m, B + (*m)*i, &one) == 0) continue;
    
    // Compute singular values of the updated matrix
    dlacpy_(&uplo_cpy, m, &one, B + (*m)*i, m, v, m);  
    double norm2 = dnrm2_(m, v, &one);
    for (long int j=0; j<*m; j++) w[j] = v[j]/norm2;
    dgemm_(&trans, &notrans, &p, &one, m, &alpha, U, m, w, m, &beta, v, m);
    for (long int j=0; j<p; j++) d[j] = s[j];
    
    t1 = clock();
    BidiagonalizeCSS(&p, v, d, e);
    t2 = clock();
    clockBidiagonalize += (t2-t1);
    d[0] = d[0]*(1-v[0]*v[0]);
    e[0] = e[0]*(1-v[0]*v[0]);
    
    t1 = clock();
    dbdsqr_(&uplo_SVD, &p, &zero, &zero, &zero, d, e, NULL, &one, NULL, &one, NULL, &one, work_local, &info); // SVD of bidiagonal matrix
    t2 = clock();
    clockBidiagonalSVD += (t2-t1);
    for (long int j=0; j<p; j++) d[j] = d[j]*d[j]; // square the singular values
    
    // Summation algorithm to get coefficients of the characteristic polynomial
    t1 = clock();
    Poly(&p, d, k, t, res_poly1, res_poly2, work_local, &lwork_local);
    t2 = clock();
    clockPoly += (t2-t1);
    
    // Check the ratio of the two interesting coefficients
    double ratio = abs(res_poly1[*k-*t])/abs(res_poly1[*k-*t-1])*exp(res_poly2[*k-*t]-res_poly2[*k-*t-1]);
    if (MinRatio < 0 || ratio<MinRatio) {
      MinRatio = ratio;
      bestCol = i;
    }
  }
  
  // Free all the space
  free(U);
  free(work_local);
  free(e);
  free(d);
  free(C);
  free(res_poly);
  free(res_poly1);
  free(res_poly2);
  free(v); 
  free(w);
  
  // Return
  return(bestCol+1);
}
