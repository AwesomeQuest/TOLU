Td = ((2:10) .*10 + 273.15)';
mud = [1.003, 0.799, 0.657, 0.548, 0.467, 0.405, 0.355, 0.316, 0.283]'*1e-3;

V0 = 0.01;
T0 = 20 + 273.15;
T1 = 80 + 273.15;
L  = 0.1;




simpson(@(x) x^3, 0,1,4)


function I=simpson(f,a,b,n)
    h=(b-a)/n;
    x=a;y=f(a);
    for i=1:n-1
        x=x+h;
        y=y+2*f(x);
    end
    x=a;
    for i=1:n
        y=y+4*f((x+x+h)/2);
        x=x+h;
    end
    y=y+f(b);
    I=h/6*y;
end