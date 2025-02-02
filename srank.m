function r = srank(A,tol)
[~,s,~] = svds(A);
s=diag(s);
if nargin==1
   tol = max(size(A)) * eps(max(s));
end
r = sum(s > tol);