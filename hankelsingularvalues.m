function sigma =hankelsingularvalues(E,A,B,C)
   R=lyapchol(A,B,E);
     L=lyap(A',C',E');
     sigma=svd(R'*L);
    semilogy(sigma)
end  
       