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
errGivens = zeros(1, maxK);
errLazy = zeros(1, maxK);
upperBound = zeros(1, maxK);
errBest = zeros(1, maxK);

% Initialize chosen column arrays
chosenColGivens = zeros(maxK, maxK);
chosenColLazy = zeros(maxK, maxK);

% Initialize times
timeGivens = zeros(1, maxK);
timeLazy = zeros(1, maxK);

numCol = [];

for k=1:maxK
    tic;
    S1 = CSS_MinE(A, k);
    timeGivens(k) = toc;
    [Q, ~] = qr(A(:,S1), 0);
    A1 = A - Q*Q'*A;
    errGivens(k) = norm(A1, 'fro');
    chosenColGivens(k, 1:k) = S1;
    
    tic;
    [S2, numCol(k)] = CSS_EarlyStop(A, k);
    timeLazy(k) = toc;
    [Q, ~] = qr(A(:,S2), 0);
    A2 = A - Q*Q'*A;
    errLazy(k) = norm(A2, 'fro');
    chosenColLazy(k, 1:k) = S2;
    
    upperBound(k) = sqrt((k+1)*sum(s(k+1:end).^2));
    errBest(k) = sqrt(sum(s(k+1:end).^2));
end

set(gcf, 'Position',  [100, 100, 1200, 400])

subplot(2,2,[1 3])
hold off
semilogy(errGivens)
hold on
semilogy(errLazy)
semilogy(upperBound, ':')
semilogy(errBest, ':');
hleg1 = legend( 'Algorithm 1', 'Alg.1, early stop', 'upper bound', 'best rank-k approx. error');
set(hleg1,'Location','best')
xlabel('k');
ylabel('approximation error')
title('Approximation error')

subplot(2,2,2)
hold off
plot(timeGivens./timeLazy)
title('Speed-up by stopping early')
xlabel('k')
hleg1 = legend('ratio time(Algorithm 1)/time(Alg.1, early stop)');
set(hleg1,'Location','best')

subplot(2,2,4)
hold off
plot(numCol)
xlabel('k')
title('Number of columns considered')
