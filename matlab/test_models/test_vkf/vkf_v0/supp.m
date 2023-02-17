function nckv = supp
simcat = 'basic';
pipedir = getdefaults('pipedir');
fname = fullfile(pipedir,simcat,'bin.mat');
data = load(fname);
o = data.o';

rng(0);
nsim = 1000;
ckv = nan(nsim,2);

i = 1;
while i<=nsim
    lambda = rand;
    v0 = exp(randn);
    omega = exp(randn);

    [~, k, v] = vkf_bin(o,lambda,v0,omega);
    ckv(i,1) = corr(v,k,'type','spearman');

    nu = rand;
    kappa = rand;
    omega = randn - 3;
    
    [~, bad_traj, ~, mu3, sigma2] = hgf_bin(o,nu,kappa,omega);    
    v = mu3(1:end-1);
    sigma2 = sigma2(2:end);
    ckv(i,2) = corr(sigma2,v,'type','spearman');
    
    i = i+1;
    if bad_traj, i = i-1; end
end

nckv = mean(ckv<0);

end