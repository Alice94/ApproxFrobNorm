#include <math.h>
#include <time.h>
#include <algorithm>
#include "Poly.cpp"
#include "BidiagonalizeCA.cpp"
#include "Tridiagonalize.cpp"
using namespace std;

// Selects the t-th cross of B, knowing that the current iteration is t. Writes the chosen indices (1-based!) in it, jt
void SelectNextCross_EarlyStop(long int* m, long int* n, double* B, long int* k, long int* t, double* it, double* jt, double* upperBound, double* numIndices);

void SelectNextCross_EarlyStop(long int* m, long int* n, double* B, long int* k, long int* t, double* it, double* jt, double* upperBound, double* numIndices) {
  *numIndices = 0;
  
  // Initialization for time checks
  clock_t t1, t2;
  clock_t clockInitialSVD = 0, clockBidiagonalize = 0, clockPoly = 0, clockBidiagonalSVD = 0;
  
  // Inizialization for SVD
  char jobu, jobvt;
  long int lwork_local = 7*min(*m,*n)+4*max(*n,*m);
  double* work_local = (double*)malloc(lwork_local*sizeof(double));
  long int info;
  long int ldvt = 1;
  char uplo_SVD = 'U';
  
  char uplo_cpy = 'A';
  char side = 'L';
  char trans = 'T';
  char notrans = 'N';
  char charN = 'N';
  
  double alpha;
  
  // Useful numbers
  long int zero = 0;
  long int one = 1;
  long int two = 2;
  long int three = 3;
  
  double oneDouble = 1.0;
  double zeroDouble = 0.0;
  
  long int p = min(*m, *n);
  long int p1 = p-1;
  long int p2 = p-2;
  double* C = (double*)malloc((*m)*(*n)*sizeof(double)); // for storing a copy of B
  double* U = (double*)malloc((*m)*(p)*sizeof(double)); // for initial SVD
  double* VT = (double*)malloc(p*(*n)*sizeof(double)); // for initial SVD
  double* s = (double*)malloc(p*sizeof(double)); // for initial SVD
  double* a = (double*)malloc(p*sizeof(double)); // column vector for update
  double* b = (double*)malloc(p*sizeof(double)); // row vector for update
  double* e = (double*)malloc(p1*sizeof(double)); // for bidiagonalization
  double* d = (double*)malloc(p*sizeof(double)); // for bidiagonalization
  double* f = (double*)malloc(p2*sizeof(double)); // for "tridiagonalization"
  double* AB = (double*)malloc(3*p*sizeof(double)); // for "tridiagonalization"
  double* res_poly = (double*)malloc((p+2)*sizeof(double)); // for characteristic polynomial
  double* res_poly1 = (double*)malloc(((*m)+2)*sizeof(double));
  double* res_poly2 = (double*)malloc(((*m)+2)*sizeof(double));
  
  
  long int bestCol = 0;
  double MinRatio = -1;
  
  // Preprocessing: Compute thin SVD of B    
  // Do the SVD of B, result is saved in s, U and VT
  dlacpy(&uplo_cpy, m, n, B, m, C, m);
  t1 = clock();
  jobu = 'S';
  jobvt = 'S';
  dgesvd(&jobu, &jobvt, m, n, C, m, s, U, m, VT, &p, work_local, &lwork_local, &info);
  t2 = clock();
  clockInitialSVD = (t2 - t1);
  
  // Sort the elements of the matrix B in descending order
  pair<double, pair<long int, long int> >* sortedIndices;
  sortedIndices = (pair<double, pair<long int, long int> >*)malloc((*n)*(*m)*(sizeof(pair<double, pair<long int, long int> >)));
  for (long int i=0; i<*m; i++) {
    for (long int j=0; j<*n; j++) {
      sortedIndices[i + j*(*m)].first = abs(B[i + j*(*m)]);
      sortedIndices[i + j*(*m)].second = make_pair(i,j);
    }
  }
  sort(sortedIndices, sortedIndices + (*m)*(*n));
  
  // Test one pair of indices at a time
  bool found = false;
  for (long int ii=(*m)*(*n)-1; ii>=0 && !found; ii--) {
      *numIndices = *numIndices + 1;
      long int i = sortedIndices[ii].second.first;
      long int j = sortedIndices[ii].second.second;
      // Check if it is already exactly zero
      if (B[j*(*m) + i] == 0) continue;
      
      t1 = clock();
      // Compute the two new vectors a & b (they have length p)
      dgemv(&trans, m, &p, &oneDouble, U, m, B + j*(*m), &one, &zeroDouble, a, &one);
      alpha = -1/B[j*(*m) + i];
      dgemv(&notrans, &p, n, &alpha, VT, &p, B + i, m, &zeroDouble, b, &one);   
      
      // Now we want to transform (Sigma - ab^T) into bidiagonal form
      // Step 1: reduce to bidiagonal plus first row
      dcopy(&p, s, &one, d, &one);
      BidiagonalizeCA(&p, a, b, d, e);
      
      // Step 2: reduce to upper-tridiagonal
      for (long int h=0; h<p; h++) {
	b[h] *= a[0];
      }
      d[0] += b[0];
      e[0] += b[1];
      
      Tridiagonalize(&p, d, e, f, b);      
      
      // Step 3: reduce to bidiagonal, Lapack does it :)
      AB[0] = AB[1] = AB[3] = 0;
      dcopy(&p, d, &one, AB+2, &three);
      dcopy(&p1, e, &one, AB+4, &three);
      dcopy(&p2, f, &one, AB+6, &three);
      
      dgbbrd(&charN, &p, &p, &zero, &zero, &two, AB, &three, d, e, NULL, &one, NULL, &one, NULL, &one, work_local, &info);

      t2 = clock();
      clockBidiagonalize += (t2-t1);
      
      t1 = clock();
      dbdsqr_(&uplo_SVD, &p, &zero, &zero, &zero, d, e, NULL, &one, NULL, &one, NULL, &one, work_local, &info); // SVD of bidiagonal matrix
      t2 = clock();
      clockBidiagonalSVD += (t2-t1);
      
      // Summation algorithm to get coefficients of the characteristic polynomial
      t1 = clock();
      for (long int h=0; h<p; h++) {
	d[h] = d[h]*d[h]; // square the singular values
      }
      Poly(m, d, k, t, res_poly1, res_poly2, work_local, &lwork_local);
      t2 = clock();
      clockPoly += (t2-t1);
      
      // Check the ratio of the two interesting coefficients
      double ratio = abs(res_poly1[*k-*t])/abs(res_poly1[*k-*t-1])*exp(res_poly2[*k-*t]-res_poly2[*k-*t-1]);
      if (MinRatio < 0 || ratio<MinRatio) {
	MinRatio = ratio;
	*it = double(i+1);
	*jt = double(j+1);
	if (ratio*(*k-*t)*(*k-*t)<=*upperBound) {
	  found = true;
	}
      }
  }
  
  // Free all the space
  free(work_local);
  free(res_poly);
  free(res_poly2);
  free(res_poly1);
  free(U);
  free(C);
  free(VT);
  free(AB);
  free(s);
  free(a);
  free(b);
  free(d);
  free(e);
  free(f);
}
