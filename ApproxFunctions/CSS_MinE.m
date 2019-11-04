function S = CSS_MinE(A, k)
% function S = CSS_MinE(A, k)
% 
% Given a matrix A and rank k, selects columns S using Algorithm 1 in
% [Cortinovis/Kressner'2019], minimizing the conditional expectation.

[m, n] = size(A);

% Initialize set of columns
S = [];
B = A;

for t = 1:k
    x = mexCSS_MinE(B, k, t);
    S = [S, x];
    B1 = A(:,S);
    [Q, ~] = qr(B1, 0);
    B = A - Q*(Q'*A);   
    B(:,S) = zeros(m, length(S));
end
