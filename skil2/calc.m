T = ((2:10) .*10 + 273.15)';
mu = [1.003, 0.799, 0.657, 0.548, 0.467, 0.405, 0.355, 0.316, 0.283]'*1e-3;


f = @(c) exp(c(1)./T + c(2)*T + c(3)*T.^2) - mu;
Df = @(c) [
    1./T T T .^2
    ] .* exp(c(1)./T + c(2)*T + c(3)*T.^2);

betahat = newton(f,Df,zeros(3,1),1e-5)
norm(f(betahat))

function u=newton(f,Df,x0,tol)
    prev = x0;
    curr = x0 - Df(x0)\f(x0);

    while norm(curr - prev) > tol
        prev = curr;
        curr = curr - Df(curr)\f(curr);
    end
    u = curr;
end
