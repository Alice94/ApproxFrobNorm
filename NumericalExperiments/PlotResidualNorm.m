% This script produces a plot describing the intermediate growth of the
% Frobenius norms of the residuals 
% A - A(:,[j1...jt])/A([i1...it],[j1...jt])*A([i1...it],:),
% which is responsible for the poor accuracy achieved by Algorithm 3 for
% some choices of matrices and ranks.

clear;

% Construct the matrix A_poly
m = 50;
n = 100;
A = zeros(m,n);
p = 10;
for i=1:m
    for j=1:n
        A(i,j) = 1/((i/n)^p + (j/n)^p)^(1/p);
    end
end

% Get the indices from Algorithm 3
k = 43;
[I, J] = CA_MinE(A, k);

B = A;
resNorm1 = zeros(1, k);
condExp1 = zeros(1, k);
for t = 1:k
    % residual norm
    B = B - B(:,J(t))*B(I(t),:)/B(I(t), J(t));
    resNorm1(t) = norm(B, 'fro');
    % conditional expectation
    s = svd(B);
    p = myPoly(s.^2);
    r = p(k-t+2)/p(k-t+1);
    condExp1(t) = r*(k-t+1)^2;
end

% Get indices from Algorithm 3, early stop
[I, J] = CA_EarlyStop(A, k);

B = A;
resNorm2 = zeros(1, k);
condExp2 = zeros(1, k);
for t = 1:k
    % residual norm
    B = B - B(:,J(t))*B(I(t),:)/B(I(t), J(t));
    resNorm2(t) = norm(B, 'fro');
    % conditional expectation
    s = svd(B);
    p = myPoly(s.^2);
    r = p(k-t+2)/p(k-t+1);
    condExp2(t) = r*(k-t+1)^2;
end

semilogy(resNorm1)
hold on
semilogy(resNorm2)
xlabel('iteration')
ylabel('Frobenius norm of residual')
legend('Algorithm 3', 'Alg. 3, early stop')
hold off



