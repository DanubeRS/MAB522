function Phi = assignment1_analytical(u, D, lambda, c0, X, T)
% Solution to dc/dt + u*dc/dx = D*d^2c/dx^2 - lambda*c, x > 0
% subject to c(0,t) = c0 and lim x->inf c(x,t) = 0 and
% initially c(x,0) = 0, x > 0.

gamma = sqrt(1 + 4*lambda*D/u^2);
Phi = (exp((u*(1-gamma)*X)/(2*D)) .* erfc((X-u*gamma*T)./(2*sqrt(D*T))) ...
     + exp((u*(1+gamma)*X)/(2*D)) .* erfc((X+u*gamma*T)./(2*sqrt(D*T)))) * c0/2;
Phi(isnan(Phi)) = 1;