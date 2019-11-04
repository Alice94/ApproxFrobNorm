% Bad matrix for cross approximation
clear; clc;
n = 6;
theta = 0.1;
k = 5;

L = -tril(ones(n),-1)*cos(theta) + eye(n);
d = sin(theta).^(0:2:(2*n-2));
A = L*diag(d)*L';

s = svd(A);

[I1, J1] = CA_MinE(A, k);
A1 = A - A(:,J1)/A(I1,J1)*A(I1,:);
errACA1 = norm(A1, 'fro');

upperBound = (k+1)*sqrt(sum(s(k+1:end).^2));
errBest = sqrt(sum(s(k+1:end).^2));

A3 = A - A(:,1:k)/A(1:k,1:k)*A(1:k,:);
errFirstK = norm(A3,'fro');

disp('Matrix A')
disp(A)
disp('Selected indices')
disp([I1', J1']);
disp('Error of the cross approximation')
disp(errACA1);
disp('Error of best rank-k approximation')
disp(errBest)
disp('Upper bound given by theorem')
disp(upperBound)
disp('Error of cross approximation using first k indices')
disp(errFirstK)
