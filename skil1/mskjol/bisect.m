function xc = bisect(a,b,tol)
counter=0;
if sign(f(a))*sign(f(b)) >= 0
  error('f(a)f(b)<0 not satisfied!') %ceases execution
end
fa=f(a);
while (b-a)/2>tol
  counter= counter+1; %Adding to the counter
  c=(a+b)/2;
  fc=f(c);
  if fc == 0              %c is a solution, done
    break
  end
  if sign(fc)*sign(fa)<0  %a and c make the new interval
    b=c;
  else                    %c and b make the new interval
    a=c;fa=fc;
  end
end
xc=(a+b)/2;               %new midpoint is best estimate
fprintf("counter  %1.0f", counter) %Prenta counterinn