#include <math.h>
#include <stdio.h>
using namespace std;

void BidiagonalizeCSS(long int* M, double* v, double* d, double* e) {
  double s, c;
  double bulgeUp = 0, bulgeDown = 0;
  double tmp;
  long int m = *M;
  for (int i=0; i<m-1; i++) e[i] = 0;
  double tol = 1e-50;
  
  for (long int i=m-1; i>0; i--) {
    // Blue
    if (max(abs(v[i-1]), abs(v[i])) > tol) {
      tmp = sqrt(v[i-1]*v[i-1] + v[i]*v[i]);
      c = v[i-1]/tmp;
      s = v[i]/tmp;

      // update v
      v[i-1] = c*v[i-1]+s*v[i];
      v[i] = 0;
      // update A
      bulgeDown = -s*d[i-1];
      d[i-1] = c*d[i-1];
      e[i-1] = s*d[i];
      d[i] = c*d[i];
      if (i+1 != m) {
	bulgeUp = s*e[i];
	e[i] = c*e[i];
      }
    }
    
    // Yellow - update A
    if (max(abs(d[i]), abs(bulgeDown)) > tol) {
      tmp = sqrt(d[i]*d[i] + bulgeDown*bulgeDown);
      c = d[i]/tmp;
      s = bulgeDown/tmp;
      
      d[i] = s*bulgeDown + c*d[i];
      bulgeDown = 0;
      tmp = c*d[i-1] - s*e[i-1];
      e[i-1] = s*d[i-1] + c*e[i-1];
      d[i-1] = tmp;
    }
    
    for (long int j=i; j<m-1; j++) {
      // Orange - update A
      if (max(abs(e[j-1]), abs(bulgeUp)) > tol) {
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
      }
      
      // Red - update A
      if (max(abs(d[j]), abs(bulgeDown)) > tol) {
	tmp = sqrt(d[j]*d[j] + bulgeDown*bulgeDown);
	c = d[j]/tmp;
	s = bulgeDown/tmp;
	
	d[j] = c*d[j] + s*bulgeDown;
	bulgeDown = 0;
	tmp = c*e[j] + s*d[j+1];
	d[j+1] = -s*e[j] + c*d[j+1];
	e[j] = tmp;
	if (j+2 != m) {
	  bulgeUp = s*e[j+1];
	  e[j+1] = c*e[j+1];
	}
      }
    }
  }  
}
