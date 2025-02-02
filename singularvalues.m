function [s,SV] =singularvalues(E,A,B,C,D,Wmin,Wup,N)
   s=logspace(Wmin,Wup, N);
    for k=1:N
     G=C*((1j*s(k)*E-A)\B)+D;
     SV(k)=max(svds(G));
    end
loglog(s,SV)
end  
       