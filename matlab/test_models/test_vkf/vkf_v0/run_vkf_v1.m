% testing vkf (figure 2) with simulated data

clc
clear all

[y,x1]          = timeseries_bin; % simulate data (outcome [y] and feedbackprobs [x])
lambda          = .1; % volatility update rate
v0              = .1; % initial volatility (or uncertainty?)
omega           = .1; % noise parameter (the higher noise parameter, the higher the speed of volatility update is)

% run vkf model for binary choices
[m1, k1, v1]    = vkf_bin(y,lambda,v0,omega);

% do stuff that I don't understand :/
m1              = 1./(1+exp(-m1)); % prediction errors
val(:,1)        = m1; 
vol(:,1)        = v1; % volatility 
kal(:,1)        = k1; % learning rate

% not required (will comment out) 
x   = x1;
v   = vol;
lr  = kal;
m   = val;

fig_plot(x1,vol,kal,val,y); % plot results