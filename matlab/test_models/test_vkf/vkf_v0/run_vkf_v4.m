% run VKF for different combinations of parameter values

clc
clear all

simcat = 'basic';
pipedir = getdefaults('pipedir');
fname = fullfile(pipedir,simcat,'bin.mat');
data = load(fname);
o = data.o(1:80)';
x = data.x(1:80)';

% 4 different combs of parameter values 
tx = [.1 .5 .1;.5 .5 .1;.1 .2 .1;.1 .5 .2];

T   = length(o); % number of trials
val = nan(T,size(tx,1));
vol = nan(T,size(tx,1));

% run model for different combs 
for j = 1:size(tx,1)

    [m1, ~, v1] = vkf_bin(o,tx(j,1),tx(j,2),tx(j,3)); % input (order): outcomes, lambda, v0, omega
    
    % no need of k1 (learning rate for now)
    % could also output learning rate
    val(:,j)    = m1; % state predictions 
    vol(:,j)    = v1; % volatility 
end

% plot figure 
fig_plot(x,vol,val,tx);