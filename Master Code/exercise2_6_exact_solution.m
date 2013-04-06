function Phi = exercise2_6_exact_solution(v, D, alpha, phi0, X, T)

gamma = sqrt(1 + 4*alpha*D/v^2);
Phi = (exp((v*(1-gamma)*X)/(2*D)) .* erfc((X-v*gamma*T)./(2*sqrt(D*T))) ...
     + exp((v*(1+gamma)*X)/(2*D)) .* erfc((X+v*gamma*T)./(2*sqrt(D*T)))) * phi0/2;
Phi(isnan(Phi)) = 1;