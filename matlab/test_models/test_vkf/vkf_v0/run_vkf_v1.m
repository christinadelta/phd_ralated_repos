% testing vkf with simulated data

[y,x1]      = timeseries_bin; % simulate data (outcome [y] and feedbackprobs [x])
lambda      = .1; % volatility learning rate
v0          = .1; % initial volatility 
omega       = .1; % noise parameter 

% run vkf model for binary choices
[m1, k1, v1] = vkf_bin(y,lambda,v0,omega);