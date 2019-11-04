function [I, J] = CA_MinE(A, k)
% function [I, J] = CA_MinE(A, k)
% Given a matrix A and rank k, outputs index sets I, J given by Algorithm 3
% in [Cortinovis/Kressner'2019], minimizing the conditional expectation.

I = [];
J = [];
B = A;

for t=1:k
    [it, jt] = mexCA_MinE(B, k, t);
    I = [I, it];
    J = [J, jt];
    B = B - B(:,jt)*B(it,:)/B(it,jt);
    B(it,:) = zeros(size(B(it,:)));
    B(:,jt) = zeros(size(B(:,jt)));
end
