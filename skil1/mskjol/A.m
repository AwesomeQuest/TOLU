function y = A(x,k,L)
    y = cosh(k*x) - cos(k*x) + ((cos(k*L) + cosh(k*L))/(sin(k*L) + sinh(k*L)))*(sin(k*x) - sinh(k*x));
end