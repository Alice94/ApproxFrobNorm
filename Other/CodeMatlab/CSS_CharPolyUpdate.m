function S = CS_CharPolyUpdate(A, k)
% function S = CS_CharPolyUpdate(A, k)
% 
% Derandomized column-subset selection
% Algorithm 4 from [Deshpande / Rademacher 2010],
% outputs a subset of k columns such that
% || A - \Pi_S(A) ||_F <= sqrt(k+1) || A - A_k ||_F
%
% This uses Algorithm 3 in [Deshpande / Rademacher 2010] to 
% update the coefficients of the characteristic polynomial 

[m, n] = size(A);

% Initialize set of columns
S = [];
B = A;

for t = 1:k
    % Initialization
    BestI = 1;
    MinRatio = -1;
    [~, Sigma, V] = svd(B, 'econ');
    sigma = diag(Sigma).^2;
    f = myPoly(sigma);
    G = zeros(m,m+1);
    for i=1:m
       G(i,2:end) = myPoly([sigma(1:i-1)', sigma(i+1:m)']);
    end
    
    % Try out the different new columns
    for i = 1:n
        % Pick the i-th column and normalize it
        bi = B(:,i);
        if (norm(bi) == 0)
            continue;
        end
        
        cmkt1 = 0;
        cmkt = 0;
        for j=1:m
           cmkt1 = cmkt1 + sigma(j)^2*V(i,j)^2*G(j,k-t+2);
           cmkt = cmkt + sigma(j)^2*V(i,j)^2*G(j,k-t+1);
        end
        cmkt1 = cmkt1/norm(bi)^2 - f(k-t+2);
        cmkt = cmkt/norm(bi)^2 - f(k-t+1);
        Ratio = abs(cmkt1/cmkt);
        
        if ((MinRatio < 0) || (Ratio < MinRatio))
            BestI = i;
            MinRatio = Ratio;
        end
    end
    
    if (MinRatio == -1)
        disp('CharPolyUpdate -- None of the rows is acceptable')
        disp(t);
    end
        
    % Select the best one
    S = [S, BestI];
    % Recompute B to be more precise
    B1 = A(:,S);
    [Q, ~] = qr(B1, 0);
    B = A - Q*(Q'*A);   
    B(:,S) = zeros(m, length(S));

end
