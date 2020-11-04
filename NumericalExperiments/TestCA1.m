% Hilbert matrix 100 x 100
% Construct matrix A
clear; clc;
n = 100;
A = hilb(n);
maxK = rank(A)+10;
numericalRank = rank(A);
s = svd(A);

% Initialize errors
errACA1 = zeros(1, maxK);
errACA2 = zeros(1, maxK);
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
numIndices = zeros(1, maxK);

for k=1:maxK
    tic;
    [I1, J1] = CA_MinE(A, k);
    time1(k) = toc;
    A1 = A - A(:,J1)/A(I1,J1)*A(I1,:);
    errACA1(k) = norm(A1, 'fro');
    chosenRow1(k, 1:k) = I1;
    chosenCol1(k, 1:k) = J1;
    
    tic;
    [I2, J2, numIndices(k)] = CA_EarlyStop(A, k);
    time2(k) = toc;
    A2 = A - A(:,J2)/A(I2,J2)*A(I2,:);
    errACA2(k) = norm(A2, 'fro');
    chosenRow2(k, 1:k) = I2;
    chosenCol2(k, 1:k) = J2;
    
    upperBound(k) = (k+1)*sqrt(sum(s(k+1:end).^2));
    errBest(k) = sqrt(sum(s(k+1:end).^2));
end

set(gcf, 'Position',  [100, 100, 1200, 400])

subplot(2,2,[1 3])
hold off
semilogy(errACA1)
hold on
semilogy(errACA2)
semilogy(upperBound, ':')
semilogy(errBest, ':');
hleg1 = legend( 'Algorithm 3', 'Alg 3, early stop', 'upper bound', 'best rank-k approx. error');
set(hleg1,'Location','best')
xlabel('k');
ylabel('approximation error')
title('Matrix approximation error')

subplot(2,2,2)
hold off
plot(time1./time2)
title('Speed-up by stopping early')
xlabel('k')
hleg1 = legend('ratio time(Alg. 3)/time(Alg. 3, early stop)');
set(hleg1,'Location','best')

subplot(2,2,4)
hold off
plot(numIndices);
xlabel('k');
title('Number of considered indices')

