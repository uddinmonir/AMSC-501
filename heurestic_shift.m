function p=heurestic_shift(A,E,l0,kl, ks)
n = size(A,1);                   
if kl >= n, error('kl must be smaller than n!'); end
if ks >= n, error('ks must be smaller than n!'); end
if (2*(l0) >= kl+ks), error('2*l0 must be smaller than kl+ks!'); 
end
rw = [];
if kl > 0 
 [Hl,Vl] = ritz_large(E,A,kl);
  rwp = eig(Hl(1:kl,1:kl));                  % =: R_+
  rw = [rw; rwp];
end
if ks > 0
   [Hs,Vs]=ritz_small(E,A,ks);
   ds = ones(ks,1)./eig(Hs(1:ks,1:ks));
     rw = [rw; ds];
end
  rw0 = rw; rw = [];
   for j = 1:length(rw0)
    if real(rw0(j))<0
      rw = [rw; rw0(j)];
    end
   end
  p=lp_mnmx(rw,l0);
return
%% large magnitude Ritz value
function [H,V] = ritz_large(E,A,k,rv) 
na = nargin;
n = size(E,1);                 
if k >= n-1, error('k must be smaller than the order of A!');
end
if na<7, rv = ones(n,1); end 

V = zeros(n,k+1);
H = zeros(k+1,k);

V(:,1) = (1.0/norm(rv))*rv;

beta = 0;

for j = 1:k
    if j > 1
    H(j,j-1) = beta;
    V(:,j) = (1.0/beta)*rv;
    end
  x = V(:,j);
  w=E\(A*x);
  rv = w;
  
  for i = 1:j
    H(i,j) = V(:,i)'*w;
    rv = rv-H(i,j)*V(:,i);
  end

  beta = norm(rv);
  H(j+1,j) = beta;
 
end  

V(:,k+1) = (1.0/beta)*rv;
%% small magnetude ritzvalue
function [H,V] = ritz_small(E,A,k,rv) 
na = nargin;
n = size(A,1);               
if k >= n-1, error('k must be smaller than the order of A!'); 
end
if na<7, rv = randn(n,1); end 
V = zeros(n,k+1);
H = zeros(k+1,k);
V(:,1) = (1.0/norm(rv))*rv;
beta = 0;
for j = 1:k
   if j > 1
    H(j,j-1) = beta;
    V(:,j) = (1.0/beta)*rv;
  end
   w = A\(E*V(:,j));
   rv = w;
   for i = 1:j
    H(i,j) = V(:,i)'*w;
    rv = rv-H(i,j)*V(:,i);
  end
  beta = norm(rv);
  H(j+1,j) = beta;
 end  
V(:,k+1) = (1.0/beta)*rv;

%% solve the minimax problem
function p = lp_mnmx(rw,l0)
if length(rw)<l0
  error('length(rw) must be at least l0.');
end
max_rr = +Inf;                      
for i = 1:length(rw)
  max_r = lp_s(rw(i),rw); 
  if max_r < max_rr
    p0 = rw(i);
    max_rr = max_r;
  end
end  

if imag(p0)
  p = [ p0; conj(p0) ];
else
  p = p0;                            
end  
[max_r,i] = lp_s(p,rw);         

while size(p,1) < l0 
  p0 = rw(i);
  if imag(p0)
    p = [ p; p0; conj(p0) ];
  else
    p = [ p; p0]; 
  end
  [max_r,i] = lp_s(p,rw);
end
%% Choose further parameters
function [max_r,ind] = lp_s(p,set)
max_r = -1;
ind = 0;
for i = 1:length(set)
  x = set(i);
  rr = 1;
  for j = 1:length(p)
    rr = rr*abs(p(j)-x)/abs(p(j)+x);
  end  
  if rr > max_r
    max_r = rr;
    ind = i;
  end
end  
