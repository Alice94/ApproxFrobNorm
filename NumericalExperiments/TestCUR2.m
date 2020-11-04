% exp matrix 100 x 200
% Construct matrix A
clear; clc;
m = 100;
n = 200;
A = zeros(m,n);
for i=1:m
    for j=1:n
        A(i,j) = exp(-0.3*abs(i-j)/n);
    end
end
maxK = min(rank(A)+10, m-1);
numericalRank = rank(A);
s = svd(A);

% Initialize errors
err1 = zeros(1, maxK);
err2 = zeros(1, maxK);
upperBound = zeros(1, maxK);
errBest = zeros(1, maxK);

% Initialize chosen indices arrays
chosenCol1 = zeros(maxK, maxK);
chosenRow1 = zeros(maxK, maxK);
chosenCol2 = zeros(maxK, maxK);
chosenRow2 = zeros(maxK, maxK);

% Initialize times
time1 = zeros(1, maxK);
time2 = zeros(1, maxK);

for k=1:maxK
    tic;
    [A1, I1, J1] = CUR_MinE(A, k);
    time1(k) = toc;
    err1(k) = norm(A-A1, 'fro');
    chosenRow1(k, 1:k) = I1;
    chosenCol1(k, 1:k) = J1;
    
    tic;
    [A2, I2, J2] = CUR_EarlyStop(A, k);
    time2(k) = toc;
    err2(k) = norm(A-A2, 'fro');
    chosenRow2(k, 1:k) = I2;
    chosenCol2(k, 1:k) = J2;
    
    upperBound(k) = sqrt(2*(k+1))*sqrt(sum(s(k+1:end).^2));
    errBest(k) = sqrt(sum(s(k+1:end).^2));
end

hold off
semilogy(err1)
hold on
semilogy(err2)
semilogy(upperBound, ':')
semilogy(errBest, ':');
hleg1 = legend( 'Algorithm 1', 'Alg.1, early stop', 'upper bound', 'best rank-k approx. error');
set(hleg1,'Location','best')
xlabel('k');
ylabel('approximation error')
