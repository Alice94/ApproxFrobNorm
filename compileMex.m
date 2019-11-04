% This compiles the mex functions
% Execute it in the folder where the script is placed.

cd Other/CodeCpp/MexFunctions;
mex -lblas -llapack mexCSS_MinE.cpp
mex -lblas -llapack mexCSS_EarlyStop.cpp
mex -lblas -llapack mexCA_MinE.cpp
mex -lblas -llapack mexCA_EarlyStop.cpp
cd ../../..
