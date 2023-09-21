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

 %% 

 clear all
 clc

 %% define paths 
 % set paths
addpath(fullfile(pwd,'functions'))

outpath         = fullfile(pwd, 'output');      addpath(outpath);
figpath         = fullfile(pwd, 'figures');     addpath(figpath);
modelpath       = fullfile(pwd, 'stcvol');    addpath(modelpath);
plotpath        = fullfile(pwd, 'plotting');    addpath(plotpath);

%% define initial variables 

% initialise variables 
subjects        = 1;
condition       = 6;                        % stable & volatile / small, medium & large stochasticity
task            = 2;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.80 .20;                 % small stochasticity probabilities
    .30 .70;                                % medium stochasticity
    .60 .40];                               % large stochasticity probabilities (either 60:40 or 64:36)
trials          = 140;                      % total trials
condtrials      = {40,[20,10,10,20]};       % 60 per stochasticity condition in stable env and 16 trials in volatile condition;
outtype         = 2;                        % if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes
nCues           = 2;

%% simulate data

data            = ALsimdata_v1(probabilities, trials, condtrials, outtype);

%% run model with simulated data 

% run model wiith simulated data first
mout            = simodel_step_one(data,probabilities, trials,condtrials, outtype);



