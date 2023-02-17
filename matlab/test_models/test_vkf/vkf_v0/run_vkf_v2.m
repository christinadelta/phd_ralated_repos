% testing vkf (figure 3) with simulated data


clc
clear all

% simulate data
[y, x1] = timeseries_bin2;

lambda = .1;
v0 = .1;
omega = .1;

% run model
[m1, k1, v1] = vkf_bin(y,lambda,v0,omega);


m1 = 1./(1+exp(-m1));
val(:,1) = m1;
vol(:,1) = v1;
kal(:,1) = k1;

% not required (will comment out) 
x   = x1;
v   = vol;
lr  = kal;
m   = val;

fig_plot2(x1,vol,kal,val,y); % plot results


