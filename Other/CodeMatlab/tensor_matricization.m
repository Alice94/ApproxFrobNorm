function Y = tensor_matricization(X, mu)
    % function Y = tensor_matricization(X, mu) builds the mu-mode
    % matricization of the tensor X
    D = 1;
    for i=[1:3]
        D = D*size(X,i);
    end
    
    perm = [mu, 1:(mu-1), (mu+1):3];
    X = permute(X, perm);
    
    Y = reshape(X, size(X,1), D/size(X,1));
end