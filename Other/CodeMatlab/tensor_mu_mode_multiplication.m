function Y = tensor_mu_mode_multiplication(X, mu, A)
    Y = tensor_matricization(X, mu);
    Y = A*Y;
    sz = size(X);
    sz2 = [size(A,1), sz(1:(mu-1)), sz((mu+1):end)];
    Y = reshape(Y, sz2);
    perm = [2:mu, 1, (mu+1):length(sz)];
    Y = permute(Y, perm);
end