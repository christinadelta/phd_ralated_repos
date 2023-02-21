function nckv = supp(outc)

o       = outc;

rng(0);
nsim    = 1000;
ckv     = nan(nsim,2);

i = 1;
while i <= nsim

    % run vkf
    lambda      = rand;
    v0          = exp(randn);
    omega       = exp(randn);

    initpars.lambda         = lambda;
    initpars.init_vol       = v0;
    initpars.omega          = omega;

    [predictions, signals]  = vkf_v1(o,initpars);
    % [~, k, v]   = vkf_bin(o,lambda,v0,omega);

    % extract output
    m                       = signals.predictions; % extract only predicted-state only for probability of vertical 
    v                       = signals.volatility;
    k                       = signals.learning_rate;
    ckv(i,1)                = corr(v,k,'type','spearman');
    
    % run hgf
    nu          = rand;
    kappa       = rand;
    omega       = randn - 3;
    
    [~, bad_traj, ~, mu3, sigma2] = hgf_bin(o,nu,kappa,omega);    
    v = mu3(1:end-1);
    sigma2 = sigma2(2:end);
    ckv(i,2) = corr(sigma2,v,'type','spearman');
    
    i = i+1;
    if bad_traj, i = i-1; end
end

nckv = mean(ckv<0);

end