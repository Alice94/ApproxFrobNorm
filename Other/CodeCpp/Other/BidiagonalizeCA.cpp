#include <math.h>
#include <stdio.h>
using namespace std;

void BidiagonalizeCA(long int* P, double* a, double* b, double* d, double* e) {
  double s, c;
  double bulgeUp = 0, bulgeDown = 0;
  double tmp;
  long int n = *P;
  for (int i=0; i<n-1; i++) e[i] = 0;
  double tol = 1e-50;
  
  for (long int i=n-1; i>0; i--) {
    // Blue
    if (max(abs(a[i-1]), abs(a[i])) > tol) {
      tmp = sqrt(a[i-1]*a[i-1] + a[i]*a[i]);
      c = a[i-1]/tmp;
      s = a[i]/tmp;

      // update a
      a[i-1] = c*a[i-1]+s*a[i];
      a[i] = 0;
      // update B
      bulgeDown = -s*d[i-1];
      d[i-1] = c*d[i-1];
      e[i-1] = s*d[i];
      d[i] = c*d[i];
      if (i+1 != n) {
	bulgeUp = s*e[i];
	e[i] = c*e[i];
      }
    }
    
    // Yellow
    if (max(abs(d[i]), abs(bulgeDown)) > tol) {
      // update B
      tmp = sqrt(d[i]*d[i] + bulgeDown*bulgeDown);
      c = d[i]/tmp;
      s = bulgeDown/tmp;
      
      d[i] = s*bulgeDown + c*d[i];
      bulgeDown = 0;
      tmp = c*d[i-1] - s*e[i-1];
      e[i-1] = s*d[i-1] + c*e[i-1];
      d[i-1] = tmp;
      
      // update b
      tmp = c*b[i-1] - s*b[i];
      b[i] = s*b[i-1] + c*b[i];
      b[i-1] = tmp;
    }
    
    for (long int j=i; j<n-1; j++) {
      // Orange
      if (max(abs(e[j-1]), abs(bulgeUp)) > tol) {
	// update B
	tmp = sqrt(e[j-1]*e[j-1] + bulgeUp*bulgeUp);
	c = e[j-1]/tmp;
	s = bulgeUp/tmp;
	
	e[j-1] = c*e[j-1] + s*bulgeUp;
	bulgeUp = 0;
	tmp = c*d[j] + s*e[j];
	e[j] = -s*d[j] + c*e[j];
	d[j] = tmp;
	bulgeDown = s*d[j+1];
	d[j+1] = c*d[j+1];
	
	// update b
	tmp = c*b[j] + s*b[j+1];
	b[j+1] = -s*b[j] + c*b[j+1];
	b[j] = tmp;
      }
      
      // Red - only update B because the corresponding entries of a are already 0
      if (max(abs(d[j]), abs(bulgeDown)) > tol) {
	tmp = sqrt(d[j]*d[j] + bulgeDown*bulgeDown);
	c = d[j]/tmp;
	s = bulgeDown/tmp;
	
	d[j] = c*d[j] + s*bulgeDown;
	bulgeDown = 0;
	tmp = c*e[j] + s*d[j+1];
	d[j+1] = -s*e[j] + c*d[j+1];
	e[j] = tmp;
	if (j+2 != n) {
	  bulgeUp = s*e[j+1];
	  e[j+1] = c*e[j+1];
	}
      }
    }
  }  
}
