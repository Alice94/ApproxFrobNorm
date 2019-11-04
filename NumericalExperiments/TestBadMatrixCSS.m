% Bad matrix for row/column selection

clear;
clc
n = 6;
alpha = 0.1;
L = -tril(ones(n),-1) + eye(n);
[L, ~] = qr(L);
d = alpha.^(0:2:(2*n-2));
A = L*diag(d)*L';
k = n-1;
s = svd(A);

[A1, I1, I2] = CUR_MinE(A, k);
err1 = norm(A-A1, 'fro');

upperBound = sqrt(2*(k+1))*sqrt(sum(s(k+1:end).^2));
errBest = sqrt(sum(s(k+1:end).^2));

[Q1, ~] = qr(A(:,1:k), 0);
A3 = A - Q1*Q1'*A*Q1*Q1';
errFirstK = norm(A3, 'fro');

disp('Matrix A')
disp(A)
disp('Selected indices')
disp([I1', I2']);
disp('Error of the CUR approximation')
disp(err1);
disp('Error of best rank-k approximation')
disp(errBest)
disp('Upper bound given by theorem')
disp(upperBound)
disp('Error of CUR approximation using first k indices')
disp(errFirstK)


