/*==========================================================
 * mexCA_EarlyStop.cpp
 * 
 * the function takes
 * B = double*, matrix MxN obtained at step t
 * K = integer, target rank
 * T = integer, index of the current pair of indices to select
 * upperBound = double, the value (k+1)*(sigma_{k+1}^2(A) + ... + sigma_m(A)^2)
 * 
 * and returns
 * it, jt = integer, indices of the new row and column
 * numIndices = number of indices it had to check
 * 
 * Compile with:
 * mex -lblas - llapack mexCA_EarlyStop.cpp
 * 
 *========================================================*/

#include <iostream>
#include "mex.h"
#include "blas.h"
#include "lapack.h"
#include "../Other/CA_EarlyStop.cpp"
using namespace std;

/* The gateway function. */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

    /* Check for proper number of arguments */
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", "MEXCPP requires 4 input arguments: B, k, t, upperBound");
    }
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout", "MEXCPP requires three output argument: index of the next selected row and column, number of analyzed indices");
    }
    
    /* Check if the input is of proper type */
    if (!mxIsDouble(prhs[0])) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin", "First argument has to be matrix of doubles.");
    }
    if (!mxIsDouble(prhs[1]) || !mxIsScalar(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsScalar(prhs[2]) || !mxIsDouble(prhs[3]) || !mxIsScalar(prhs[3])) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin", "Second, third and fourth argument have to be double scalar.");
    }
    
    /* Acquire pointers to the input data */
    
    double* B = mxGetPr(prhs[0]);
    long int K = round(mxGetScalar(prhs[1]));
    long int T = round(mxGetScalar(prhs[2]))-1; // -1 because 0-based
    double upperBound = mxGetScalar(prhs[3]);
    
    long int M = (long int)mxGetM(prhs[0]);
    long int N = (long int)mxGetN(prhs[0]);
  
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    double* it = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    double* jt = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    double* numIndices = mxGetPr(plhs[2]);
    SelectNextCross_EarlyStop(&M, &N, B, &K, &T, it, jt, &upperBound, numIndices);
}













