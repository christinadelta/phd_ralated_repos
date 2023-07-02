% second version of the aversive learning task with the
% stochasticity-volatility model 

% STEPS:

% created May July 2023
% changes introduced:
% changed the environment of the experiment
% changed the number of trials

clear all
clc

%% TESTING STCS:
% For large stochasticity will test contigences: [68% 32%] & [60% 40%]
% Stable condition will include 40 trials. The environment looks like this:
% Low stc trials: [40 25 20 35 40] -- first 40 & last 40 constitute stable env
% High stch: [40 25 20 35 40]

%% set paths and init variables

% set paths
addpath(fullfile(pwd,'functions'))

outpath         = fullfile(pwd, 'output');      addpath(outpath);
figpath         = fullfile(pwd, 'figures');     addpath(figpath);
modelpath       = fullfile(pwd, 'volstoch');    addpath(modelpath);
plotpath        = fullfile(pwd, 'plotting');    addpath(plotpath);

% initialise variables 
subjects        = 1;
condition       = 4;                        % stable & volatile / small & large stochasticity
task            = 2;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.80 .20;                 % small stochasticity probabilities
    .68 .32];                               % large stochasticity probabilities (either 60:40 or 64:36)
trials          = 160;                      % total trials
condtrials      = [40 16];                  % 80 per stochasticity condition in stable env and 16 trials in volatile condition;
outtype         = 2;                        % if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes
nCues           = 2;

% define initial learning rate inverse temperature parameters (don't think
% that this will be used)
% params          = [.2 .6;
%     3 6]; %alpha and beta parameteres 

% models to run:
mtorun = 3; % (first test the PF using the exact same params as Piray (2021), then run healthy model and stochasticity model

%% simulate dataset

for sub = 1:subjects

    % simulate dataset(s)
    data            = action_simdataV1(condition, probabilities, trials,condtrials, outpath, outtype, task);
    alldata{sub,1}  = data;


end % end of subjects loop




