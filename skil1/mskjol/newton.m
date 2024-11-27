function x=newton(x0,tol)
%skilgreina hér fallið f(x)
%f=@(x) ...
%og f'(x)
%Df=@(x) ...
x=x0;oldx=x-100;
%count=0;
while abs(x-oldx)>tol
   %count=count+1;
   oldx=x;
   x=x-f(x)/Df(x); 
end
%fprintf("Counter newton %0.1f", count)
end