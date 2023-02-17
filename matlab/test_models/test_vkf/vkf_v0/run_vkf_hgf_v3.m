% testing vkf and HGF (figure 6) with simulated data

clear all
clc

simcat      = 'basic';
pipedir     = getdefaults('pipedir');
fname       = fullfile(pipedir,simcat,'bin.mat');
data        = load(fname);
o           = data.o';
x           = data.x';

[tx_vkf, tx_hgf] = sim_lr_vol;

% get vkf parameters
lambda      = tx_vkf(1);
v0          = tx_vkf(2);
omega       = tx_vkf(3);

% run vkf model
[m, k, v]   = vkf_bin(o,lambda,v0,omega);
val1        = m;
vol1        = v;
lr1         = k;

% get hgf parameters
nu          = tx_hgf(1);
kappa       = tx_hgf(2);
omega       = tx_hgf(3);

% run hgf model
[~, ~, mu2, mu3, sigma2] = hgf_bin(o,nu,kappa,omega);    

m           = mu2(1:end-1);
v           = (mu3(1:end-1));
sigma2      = sigma2(2:end);    
val2        = m;
vol2        = v;
lr2         = sigma2;

% plot results
fig_plot3(x,o,val1,vol1,lr1,val2,vol2,lr2);

