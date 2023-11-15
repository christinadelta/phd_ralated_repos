% alt stc-vol model based on the stc-vol model by Piray & Daw, 2021 
% VERSION 1 -- November 2023
% @christinadelta 

% comments will go here 

%% clear sruff

clear all
clc

%% define paths 
% set paths
simpath         = fullfile(pwd, 'datasim');     addpath(simpath);
modelpath       = fullfile(pwd, 'model');       addpath(modelpath);
plotpath        = fullfile(pwd, 'plotting');    addpath(plotpath);

%% generate some data

data            = generateData();

%% define parameters for the model (version 1)
% we'll now run a model that estimates voaltility only 

% we need:
% 1. outcomes 
% 2. lambda (volatility update term)
% 3. sigma (constant/noise parameter)
params                      = struct('lambda_v',.1, 'v0',.1,'sigma',.1, 'y', data.o);
[predict, vol, lr, sigm]    = kf_v1(params);

% we need to bring predited signal to positive values only (between 0 and 1)
m1 = 1./(1+exp(-predict));

%% plot trial-by-trial 

col = 1;
row = 3;
x   = data.x;
y   = data.o;
f   = plotsignal(col,row,m1,vol,lr,x,y);

%% plot learning rates for stable vs volatile

h = plotlrsV1(lr, data.volind);

%%
signals = kf_v1_1(params);

% unpack signals struct
posterior_m = signals.predictions;
vol2 = signals.volatility;
lr2 = signals.learning_rate;

m2 = 1./(1+exp(-posterior_m));

%% plot 2nd version of the kalman filter with volatility

col = 1;
row = 3;
x   = data.x;
y   = data.o;
f   = plotsignal(col,row,m2,vol2,lr2,x,y);

%% plot 2nd version learning rates for stable vs volatile

h = plotlrsV1(lr2, data.volind);

%% build the second part of the model (add second noise type)

% include stc and lambda_s
params                      = struct('lambda_v',.1, 'lambda_s', .1, 'v0',.1,'s0',.1, 'y', data.o);
[signals]                   = kf_v2(params);

m   = signals.predictions;
lr  = signals.learning_rate;
vol = signals.volatility;
stc = signals.stochasticity;

m = 1./(1+exp(-m));

%% plot the kalman filter with volatility + stochasticty 

col = 1;
row = 4;
x   = data.x;
y   = data.o;
f   = plotsignal2(col,row,m,vol,stc,lr,x,y);

%%










