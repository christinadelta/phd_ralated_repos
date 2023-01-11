function [Xfit, LL, BIC] = model1_fit(a, r)

obFunc = @(x) model1_lik(a, r, x);

X0 = rand;
LB = 0;
UB = 1;
[Xfit, NegLL] = fmincon(obFunc, X0, [], [], [], [], LB, UB);

LL = -NegLL;
BIC = length(X0) * log(length(a)) + 2*NegLL;

return