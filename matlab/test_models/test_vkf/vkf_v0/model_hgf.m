function [loglik, tx, c] = model_hgf(params,data)
ux = @(x)(1./(1+exp(-x)));

y = data;
nu = ux(params(1));
kappa = ux(params(2));        
omega = params(3) -1;
tx = [nu kappa omega];

[~,bad_traj,m, mu3, sigma2] = hgf_bin(y,nu,kappa,omega);
m = 1./(1+exp(-m));

mu3 = mu3(1:end-1,:);
sigma2 = sigma2(2:end,:);
c = corr(mu3,sigma2,'type','spearman'); 
    
mu1 = m(1:end-1,:);
p = mu1.*y + (1-mu1).*(1-y);        
loglik = sum(sum(log(p+eps)));
if bad_traj>0
    loglik = -10^16;
end   


end
