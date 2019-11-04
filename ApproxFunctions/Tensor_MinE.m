function [X, S1, S2, S3] = Tensor_MinE(A, k1, k2, k3)
% function [X, S1, S2, S3] = Tensor_MinE(A, k1, k2, k3)
% 
% Given a 3-dimensional tensor A and ranks k1, k2, k3
% returns approximation X and fibser indices S1, S2, S3 of cardinality
% k1, k2, k3 respectively.
% Uses Algorithm 4 in [Cortinovis/Kressner'2019] with minimization of
% conditional expectation.

A1 = tensor_matricization(A, 1);
S1 = CSS_MinE(A1, k1);

A2 = tensor_matricization(A, 2);
S2 = CSS_MinE(A2, k2);

A3 = tensor_matricization(A, 3);
S3 = CSS_MinE(A3, k3);

% Create core tensor
[Q1, ~] = qr(A1(:,S1), 0);
[Q2, ~] = qr(A2(:,S2), 0);
[Q3, ~] = qr(A3(:,S3), 0);

C = tensor_mu_mode_multiplication(A,1,Q1');
C = tensor_mu_mode_multiplication(C,2,Q2');
C = tensor_mu_mode_multiplication(C,3,Q3');

% Construct approximated tensor
X = tensor_mu_mode_multiplication(C,1,Q1);
X = tensor_mu_mode_multiplication(X,2,Q2);
X = tensor_mu_mode_multiplication(X,3,Q3);