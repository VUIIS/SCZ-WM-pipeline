function [c,R2] = multpolyfit(x,y,n)

m = size(x,2);
c = zeros(n+1,m);
r = zeros(size(y));

for k = 1:m
    M = repmat(x(:,k),1,n+1);
    M = bsxfun(@power,M,0:n);
    c(:,k) = M\y(:,k);
    r(:,k) = M*c(:,k)-y(:,k);
end

sserr = sum(r.^2);
sstot = sum(bsxfun(@minus,y,mean(y)).^2);
R2 = 1 - sserr./sstot;
