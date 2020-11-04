function [I, J, numIndices] = CA_EarlyStop(A, k)
% function [I, J] = CA_MinE(A, k)
% Given a matrix A and rank k, outputs index sets I, J given by Algorithm 3
% in [Cortinovis/Kressner'2019] with early stop.

I = [];
J = [];
B = A;

numIndices = 0;

s = svd(A);
upperBound = (k+1)^2*sum(s(k+1:end).^2);

for t=1:k
    [it, jt, res] = mexCA_EarlyStop(B, k, t, upperBound);
    numIndices = numIndices + res;
    
    I = [I, it];
    J = [J, jt];
    B = B - B(:,jt)*(B(it,:)/B(it,jt));
    B(it,:) = zeros(size(B(it,:)));
    B(:,jt) = zeros(size(B(:,jt)));
end
