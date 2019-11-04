% Hilbert tensor
% Construct tensor A
clear; clc;
n = 50;
A = zeros(n,n,n);
for i=1:n
    for j=1:n
        for h=1:n
            A(i,j,h) = 1/(i+j+h-1);
        end
    end
end

maxK = 20;
numericalRank = rank(tensor_matricization(A, 1));

% Initialize errors
err1 = zeros(1, maxK);
err2 = zeros(1, maxK);

% Initialize numCol
numCol = zeros(1, maxK);

% Initialize times
time1 = zeros(1, maxK);
time2 = zeros(1, maxK);

for k=1:maxK
    tic;
    [X1, ~, ~, ~] = Tensor_MinE(A, k, k, k);
    time1(k) = toc;
    err1(k) = norm(tensor_matricization(A-X1,1), 'fro');
    
    tic;
    [X2, ~, ~, ~, numCol(k)] = Tensor_EarlyStop(A, k, k, k);
    time2(k) = toc;
    err2(k) = norm(tensor_matricization(A-X2, 1), 'fro');
end

hold off
semilogy(err1)
hold on
semilogy(err2)
SVD = svd(tensor_matricization(A,1));
best1 = zeros(maxK);
for i=1:maxK
    best1(i) = sqrt(3*sum(SVD(i+1:end).^2));
end
semilogy(best1, 'k:')
hleg1 = legend( 'Algorithm 1', 'Alg.1, early stop', 'quasi-best approx. error');
set(hleg1,'Location','best')
xlabel('k');
ylabel('log of approx. error')