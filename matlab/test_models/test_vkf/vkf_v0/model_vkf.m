function [loglik, tx, c] = model_vkf(params,data)


ux = @(x)(1./(1+exp(-x)));

y = data;
lambda = ux(params(1));
v0 = 10*ux(params(2));
omega = exp(params(3));        
tx = [lambda, v0, omega];

[~,k,v,m] = vkf_bin(y,lambda,v0,omega);
m = [.5*ones(1,size(m,2)); m];

c = corr(v,k,'type','spearman');

mu = m(1:end-1,:);
p = mu.*y + (1-mu).*(1-y);
loglik = sum(sum(log(p+eps)));


end