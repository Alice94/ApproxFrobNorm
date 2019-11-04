#include <math.h>
#include <stdio.h>
using namespace std;

void Poly(long int* m, double* a, long int* k, long int* t, double* res1, double* res2, double* work, long int* lwork) {
  double* work1 = work;
  double* work2 = work+(*m)+1;
  double x;
  double minExponent;
  
  for (long int i=0; i<*m; i++) {
    res1[i] = 0;
    res2[i] = 0;
  }
  res1[0] = 1;
  
  work1[0] = 1;
  work2[0] = 0;
  
  work2[1] = floor(log(a[0]));
  work1[1] = a[0]/exp(work2[1]);
  for (long int i=2; i<*m+1; i++) {
    work1[i] = 0;
    work2[i] = 0;
  }
  
  for (long int i=0; i<*m-1; i++) {
    for (long int j=1; j<=(*k)-(*t); j++) {
      minExponent = min(work2[j], work2[j-1]);
      if (work1[j] == 0) {
	x = work1[j-1]*a[i+1];
	minExponent = work2[j-1];
      }
      else if (work1[j-1] == 0) {
	x = work1[j];
	minExponent = work2[j];
      }
      else {
	x = work1[j-1]*a[i+1]*exp(work2[j-1]-minExponent) + work1[j]*exp(work2[j]-minExponent);
      }
      if (x == 0) {
	res1[j] = 0;
	res2[j] = 0;
      }
      else {
	res2[j] = floor(log(x));
	res1[j] = x/exp(res2[j]);
	res2[j] += minExponent;
      }
    }
    for (long int j=0; j<=(*k)-(*t); j++) {
      work1[j] = res1[j];
      work2[j] = res2[j];
    }
  }
  
  for (long int i=0; i<*m+1; i++) {
    res1[i] = work1[i];
    res2[i] = work2[i];
  }
}
