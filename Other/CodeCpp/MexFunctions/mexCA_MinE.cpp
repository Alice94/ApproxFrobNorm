/*==========================================================
 * mexCA_MinE.cpp
 * 
 * the function takes
 * B = double*, matrix MxN obtained at step t
 * K = integer, target rank
 * T = integer, index of the current pair of indices to select
 * 
 * and returns
 * it, jt = integer, indices of the new row and column
 * 
 * Compile with:
 * mex -lblas -llapack mexCA_MinE.cpp
 * 
 *========================================================*/

#include <iostream>
#include "mex.h"
#include "blas.h"
#include "lapack.h"
#include "../Other/CA_MinE.cpp"
using namespace std;

/* The gateway function. */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

    /* Check for proper number of arguments */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", "MEXCPP requires 3 input arguments: B, k, t");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout", "MEXCPP requires two output argument: index of the next selected row and column");
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
    double* it = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    double* jt = mxGetPr(plhs[1]);
    SelectNextCross_MinE(&M, &N, B, &K, &T, it, jt);
}













