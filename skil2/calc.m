L  = 0.1;
V(L/2,0.5e-10)

function I = V(y,n)
    beta1 = [
         -8.94378370228994
         -839.223323491023
          421125.826871677
    ]; % error  2.68698517476343e-06

    V0 = 0.01;
    T0 = 20 + 273.15;
    T1 = 80 + 273.15;
    L  = 0.1;

    mu12 = @(T,c) exp(c(1) + c(2) ./T' + c(3)./T'.^2);

    T = @(x) T0 + x*(T1-T0)/L;
    mu = @(xi) mu12(T(xi),beta1);
    K = adapquadcount(@(xi) 1/mu(xi), 0,L,n);
    U = adapquadcount(@(xi) 1/mu(xi), 0,y,n);
    disp(K(2)+U(2))
    I = V0/K(1) * U(1);
end

function I=adapquad(f,a,b,tol)
    I0 = simpson(f,a,b,1);
    I1 = simpson(f,a,b,2);
    if abs(I1-I0) < 10*tol
        I = I1;
    else
        I = adapquad(f,a,(a+b)/2,tol/2) + adapquad(f,(a+b)/2,b,tol/2);
    end
end
function I=adapquadcount(f,a,b,tol)
    I0 = simpson(f,a,b,1);
    I1 = simpson(f,a,b,2);
    if abs(I1-I0) < 10*tol
        I = [I1, 1];
    else
        I = adapquadcount(f,a,(a+b)/2,tol/2) + adapquadcount(f,(a+b)/2,b,tol/2);
    end
end

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