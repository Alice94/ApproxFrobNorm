% Symmetric positive definite matrix for which the best cross approximation is not symmetrix

clear; clc;
A = [1.87, -1.82, -2.11; -1.82, 1.87, 2.11; -2.11, 2.11, 2.54];
n = 3;
bestErrDiag = 100000;
bestErr = 100000;
disp('Matrix A for which we want a rank-1 cross approximation:')
disp(A)
disp('It is symmetric positive semidefinite, its eigenvalues are')
disp(eig(A))
disp('Best rank-1 approximation error in Frobenius norm:')
s = svd(A);
s = sqrt(sum(s(2:end).^2));
disp(s);
disp('Bound given by theorem:')
disp(s*2);

disp('[i, j, error of cross approx in Frob norm]');
for i=1:n
    for j=1:n
        B = A - A(:,j)*A(i,:)/A(i,j);
        err = norm(B, 'fro');
        disp([i, j, err]);
        if (i==j)
            bestErrDiag = min(bestErrDiag, err);
        end
        bestErr = min(bestErr, err);
    end
end


disp('Best rank-1 cross approximation error in Frobenius norm:')
disp(bestErr)
disp('Best diagonal rank-1 cross approximation error in Frobenius norm:')
disp(bestErrDiag)


