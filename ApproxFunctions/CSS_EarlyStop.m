function [S, colNum] = CSS_EarlyStop(A, k)
% function [S, colNum] = CSS_EarlyStop(A, k)
% 
% Given a matrix A and rank k, selects columns S using Algorithm 1 in
% [Cortinovis/Kressner'2019]with early stop. colNum is the number of
% columns that the algorithm had to consider (it will be >= k)

[m, n] = size(A);

% Initialize set of columns
S = [];
B = A;
svA = svd(A);
upperBound = (k+1)*sum(svA(k+1:end).^2);
colNum = [];

for t = 1:k
    [x, colNum(t)] = mexCSS_EarlyStop(B, k, t, upperBound);
    S = [S, x];
    B1 = A(:,S);
    [Q, ~] = qr(B1, 0);
    B = A - Q*(Q'*A);    
    B(:,S) = zeros(m, length(S));
end

colNum = sum(colNum);
