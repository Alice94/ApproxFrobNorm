A = [6.583644e-7,  8.113362e-3; 8.113362e-3, 100]

c1 = CSS_MinE(A, 1);
c2 = CSS_CharPolyUpdate(A, 1);
c3 = CSS_CharPolyUpdate(vpa(A),1);

disp('Column subset selection with update of singular values chooses column')
disp(c1)
disp('and attains approximation error')
disp(norm(A - A(:,c1)/norm(A(:,c1)) * (A(:,c1)/norm(A(:,c1)))' * A, 'fro'))

disp('Column subset selection with update of characteristic polynomial chooses column')
disp(c2)
disp('and attains approximation error')
disp(norm(A - A(:,c2)/norm(A(:,c2)) * (A(:,c2)/norm(A(:,c2)))' * A, 'fro'))

disp('Column subset selection with update of characteristic polynomial, applied on vpa(A), chooses column')
disp(c3)
disp('and attains approximation error')
disp(norm(A - A(:,c3)/norm(A(:,c3)) * (A(:,c3)/norm(A(:,c3)))' * A, 'fro'))
disp('(This is to check that the algorithm which updated the characteristic polynomial is mathematically correct).')


