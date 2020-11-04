# Low-rank approximation in the Frobenius norm

This contains the implementations of the algorithms presented in [1], as well as scripts that reproduce all the numerical experiments contained there.

----------------------------------------------------------------

### Installation

Compile the mex files running the script "compileMex" in Matlab (staying in this folder).

----------------------------------------------------------------

### Numerical experiments
The scripts to reproduce the numerical experiments presented in [1] are in the folder "NumericalExperiments". 

Here is the correspondence between the Figures in the paper and the scripts in the NumericalExperiments folder

- Figure 1 --> TestCSS1.m

- Figure 2 --> TestCSS2.m

- Figure 3 --> TestCSS3.m

- Figure 4 --> TestCUR1.m, TestCUR2.m, TestCUR3.m

- Figure 6 --> TestCA1.m

- Figure 7 --> TestCA2.m

- Figure 8 --> TestCA3.m

- Figure 9 --> PlotResidualNorm.m

- Figure 10 --> TestTensor1.m (left), TestTensor2.m (right)

- Example at the end of Section 3.1.1 --> TestBadMatrixCSS.m

- Other examples in Section 3.2.3 --> IntermediateResGrowth.m, TestBadMatrixCA.m, TestSPDmatrixCA.m

- Example at the end of page 5 --> CharPolyUpdateFails.m

----------------------------------------------------------------

### Functions for matrix/tensor approximation
List of Matlab functions in folder "ApproxFunctions"

- Column subset selection (Algorithm 1): CSS_EarlyStop.m, CSS_minE.m

- Matrix CUR approximation (Algorithm 2): CUR_EarlyStop.m, CUR_minE.m

- Cross approximation (Algorithm 3): CA_EarlyStop.m, CA_minE.m

- Tensor approximation (Algorithm 4): Tensor_EarlyStop.m, Tensor_minE.m

----------------------------------------------------------------
### References

[1] A. Cortinovis and D. Kressner. "Low-rank approximation in the Frobenius norm by column and row subset selection". arXiv preprint arXiv:1908.06059, 2019.
