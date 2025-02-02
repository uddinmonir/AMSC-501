y = my_step(E,A,B,C,Tmax,N);
h=(Tmax-0)/N;
 t = linspace(0,Tmax,N);
 yout = zeros(size(C,1),n);
 [L,U,P,Q] = lu(E-h*A);
  for i=1:n
    x=Q*(U\(L\(P*(E*x+h*B))));
    y(:,i)=C*x;
  end
   plot(t,y)