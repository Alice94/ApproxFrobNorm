#include <math.h>
#include <stdio.h>
using namespace std;

void Tridiagonalize(long int* p, double* d, double* e, double* f, double* b) {
  double tol = 1e-50;
  double c, s, tmp;
  double bulgeUp = 0, bulgeDown = 0;
  long int row, col;
  
  for (long int i=*p-3; i>0; i--) {
    // First transformation, zero out the (i+2)th element of b
    // Introduces bulgeDown and f[i]
    if (max(abs(b[i+1]), abs(b[i+2])) > tol) {
      tmp = sqrt(b[i+1]*b[i+1] + b[i+2]*b[i+2]);
      c = b[i+1]/tmp;
      s = b[i+2]/tmp;
      
      b[i+1] = tmp;
      b[i+2] = 0;
      
      f[i] = -s*e[i];
      e[i] = c*e[i];
      tmp = c*d[i+1] + s*e[i+1];
      e[i+1] = -s*d[i+1] + c*e[i+1];
      d[i+1] = tmp;
      bulgeDown = s*d[i+2];
      d[i+2] = c*d[i+2];
      
    }
    row = i+2; 
    // Second transformation, zeroes the bulgeDown and maybe introduces a bulgeUp
    if (max(abs(bulgeDown), abs(d[i+1])) > tol) {
      tmp = sqrt(bulgeDown*bulgeDown + d[i+1]*d[i+1]);
      s = bulgeDown/tmp;
      c = d[i+1]/tmp;
      
      d[i+1] = tmp;
      bulgeDown = 0;
      tmp = c*e[i+1] + s*d[i+2];
      d[i+2] = -s*e[i+1] + c*d[i+2];
      e[i+1] = tmp;
      if (i < *p-3) {
	tmp = c*f[i+1] + s*e[i+2];
	e[i+2] = -s*f[i+1] + c*e[i+2];
	f[i+1] = tmp;
      }
      if (i+2 < (*p)-2) {
	bulgeUp = s*f[i+2];
	f[i+2] = c*f[i+2];
      }
    }  
    row = i+2;
    
    while (row < (*p)-2) {
      // there is a bulgeUp to chase!
      col = row+2;
      if (max(abs(bulgeUp), abs(f[col-3])) > tol) {
	tmp = sqrt(bulgeUp*bulgeUp + f[col-3]*f[col-3]);
	c = f[col-3]/tmp;
	s = bulgeUp/tmp;
	
	bulgeUp = 0;
	f[col-3] = tmp;
	
	tmp = e[col-2]*c + f[col-2]*s;
	f[col-2] = -e[col-2]*s + f[col-2]*c;
	e[col-2] = tmp;
	
	tmp = d[col-1]*c + e[col-1]*s;
	e[col-1] = -d[col-1]*s + e[col-1]*c;
	d[col-1] = tmp;
	
	bulgeDown = s*d[col];
	d[col] = c*d[col];
      }
      
      // now there is a bulgeDown to chase...
      row = col;
      if (max(abs(bulgeDown), abs(d[row-1])) > tol) {
	tmp = sqrt(bulgeDown*bulgeDown + d[row-1]*d[row-1]);
	c = d[row-1]/tmp;
	s = bulgeDown/tmp;
	
	bulgeDown = 0;
	d[row-1] = tmp;
	
	tmp = c*e[row-1] + s*d[row];
	d[row] = -s*e[row-1] + c*d[row];
	e[row-1] = tmp;
	
	if (row <= (*p)-2) {
	  tmp = c*f[row-1] + s*e[row];
	  e[row] = -s*f[row-1] + c*e[row];
	  f[row-1] = tmp;
	}
	
	if (row <= (*p)-3) {
	  bulgeUp = s*f[row];
	  f[row] = c*f[row];
	}
      }
    }
  }
  
  f[0] = b[2];
}
