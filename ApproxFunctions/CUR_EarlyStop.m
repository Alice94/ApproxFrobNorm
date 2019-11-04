function [B, S1, S2] = CUR_EarlyStop(A, k)
% function [B, S1, S2] = CUR_EarlyStop(A, k)
% Given matrix A and rank k outputs a CUR approximation B using rows S1 and
% columns S2. Uses Algorithm 2 in [Cortinovis/Kressner'2019] with early
% stop.

[n, m] = size(A);

[S1, ~] = CSS_EarlyStop(A', k);
[S2, ~] = CSS_EarlyStop(A,k);

B1 = A(S1,:);
B2 = A(:,S2);

[Q1, ~] = qr(B1', 0);
[Q2, ~] = qr(B2, 0);
B = Q2*(Q2'*A*Q1)*Q1';
