% This is an example of small matrix for which we can observe intermediate
% residual growth during the derandomized cross approximation algorithm

A = [-10^(-4), 3, -4; 4, 1, 2; 8, -1, 1];
k = 2;
disp('Matrix A')
disp(A)
disp('Aiming for a cross approximation  of rank')
disp(k)

disp('We look at the first step (t=1) of the algorithm and for each pair of indices (i,j) we compute')
disp('E[ | A - A(:,Y)A(X,:)/A(X,Y) |_F^2 \mid X1 = i, Y1 = j]')
for i=1:3
    for j=1:3
        B = A - A(:,j)*A(i,:)/A(i,j);
        s = svd(B);
        E = 4*s(2);
        disp([i, j, E])
    end
end


disp('So the pair that minimizes E is (1,1)')
disp('The residual A - A(:,1)A(1,:)/A(1,1) is')
disp(A - A(:,1)*A(1,:)/A(1,1))
disp('and its Frobenius norm is ')
disp(norm(A - A(:,1)*A(1,:)/A(1,1), 'fro'))
disp('which is')
disp(norm(A - A(:,1)*A(1,:)/A(1,1), 'fro')/norm(A,'fro'))
disp('times bigger than the norm of A.')