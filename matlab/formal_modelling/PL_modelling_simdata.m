% main script for modelling perceptual learning simulations & empirical data 
 % created September 2023

 % 1. simulate data
 % 2. run healthy model with simulated data ( all three stc levels)
 % 3. run healthy model and lession models with simulated for each stc
 % level seperately
 % 4. for simulated and for empirical data, estimate learning rates, vol
 % stc and performance and plot them trial-by-trial like in figure 2 in the
 % Piray, 2021 paper

 %--------------------------

 clear all
clc

%% set paths
addpath(fullfile(pwd,'functions'))

outpath         = fullfile(pwd, 'output');      addpath(outpath);
figpath         = fullfile(pwd, 'figures');     addpath(figpath);
modelpath       = fullfile(pwd, 'stcvol');      addpath(modelpath);
plotpath        = fullfile(pwd, 'plotting');    addpath(plotpath);
modelfitpath    = fullfile(pwd, 'modelfit');    addpath(modelfitpath);

%% define initial variables 

% initialise variables 
subjects        = 1;
condition       = 6;                        % stable & volatile / small, medium & large stochasticity
task            = 1;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.80 .20;                 % small stochasticity probabilities
    .70 .30;                                % medium stoch probabilities
    .60 .40];                               % large stochasticity probabilities (either 60:40 or 64:36)
trials          = 140;                      % total trials per stc level
condtrials      = {70,[20,10,10,30]};       % 50 per stochasticity condition in stable env and 10 trials in volatile condition;
outtype         = 2;                        % if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes
nCues           = 2;
nOut            = 2;

%% simulate data for the model

data            = PLsimdata_v1(probabilities, trials, condtrials, outtype);

%%