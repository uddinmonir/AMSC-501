function h2 =h2norm(E,A,B,C)
   P=lyap(A,B*B',[],E);
   h2=sqrt(trace(C*P*C'));
end  
    