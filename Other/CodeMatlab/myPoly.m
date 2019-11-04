function p = myPoly(v)

n = length(v);
p = zeros(1,n+1);
p(1) = 1;

q = [1 v(1) zeros(1,n-1)];
for i=1:n-1
    for j=2:n+1
        p(j) = v(i+1)*q(j-1) + q(j);
    end
    q = p;
end
p = q;

