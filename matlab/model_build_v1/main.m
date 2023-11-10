% alt stc-vol model based on the stc-vol model by Piray & Daw, 2021 
% VERSION 1 -- November 2023
% @christinadelta 

% comments will go here 

%% clear sruff

clear all
clc

%% define paths 
% set paths
simpath         = fullfile(pwd, 'datasim'); addpath(simpath);
modelpath       = fullfile(pwd, 'model');    addpath(modelpath);
% plotpath        = fullfile(pwd, 'plotting');    addpath(plotpath);

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

% we need to bring predited signal to positive values only (between 0 and
% 1)
m1 = 1./(1+exp(-predict));

%% plot trial-by-trial 



