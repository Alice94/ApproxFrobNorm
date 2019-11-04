/*==========================================================
 * mexCSS_MinE.cpp
 * 
 * the function takes
 * B = double*, matrix MxN obtained at step t
 * K = integer, number of columns to be selected
 * T = integer, index of the current column to select

 * 
 * and returns
 * res = integer, index of the returned column
 * 
 * Compile with:
 * mex -lblas -llapack mexCSS_MinE.cpp
 * 
 *========================================================*/

#include <iostream>
#include "mex.h"
#include "blas.h"
#include "lapack.h"
#include "../Other/CSS_MinE.cpp"
using namespace std;

/* The gateway function. */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

    /* Check for proper number of arguments */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", "MEXCPP requires 3 input arguments: B, k, t");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout", "MEXCPP requires one output argument: index of the next selected column");
    }
    
    /* Check if the input is of proper type */
    if (!mxIsDouble(prhs[0])) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin", "First argument has to be matrix of doubles.");
    }
    if (!mxIsDouble(prhs[1]) || !mxIsScalar(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsScalar(prhs[2])) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin", "Second and third argument have to be double scalar.");
    }
    
    /* Acquire pointers to the input data */
    
    double* B = mxGetPr(prhs[0]);
    long int K = round(mxGetScalar(prhs[1]));
    long int T = round(mxGetScalar(prhs[2]))-1; // -1 because 0-based
    
    long int M = (long int)mxGetM(prhs[0]);
    long int N = (long int)mxGetN(prhs[0]);
  
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    double* p = mxGetPr(plhs[0]);
    *p = (double)SelectNextColumn_MinE(&M, &N, B, &K, &T);
}
