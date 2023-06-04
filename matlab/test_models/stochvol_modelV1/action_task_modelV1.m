% first version of the aversive learning task with the
% stochasticity-volatility model 

% STEPS:

% created May 2023 

clear all
clc

%% set paths and init variables

% set paths
addpath(fullfile(pwd,'functions'))

outpath         = fullfile(pwd, 'output');      addpath(outpath);
figpath         = fullfile(pwd, 'figures');     addpath(figpath);
modelpath       = fullfile(pwd, 'volstoch');   addpath(modelpath);

% initialise variables 
subjects        = 1;
condition       = 4;                        % stable & volatile / small & large stochasticity
% task            = 1;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.88 .12;                 % small stochasticity probabilities
    .60 .40];                               % large stochasticity probabilities (either 60:40 or 64:36)
trials          = 200;                      % total trials
condtrials      = [100 25];                 % 100 per stochasticity condition in stable env and 25 trials in volatile condition;
outtype         = 2;                        % if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes

% define initial learning rate inverse temperature parameters 
params          = [.2 .6;
    3 6]; %alpha and beta parameteres 

%% simulate dataset

for sub = 1:subjects

    % simulate dataset(s)
    data            = action_simdataV1(condition, probabilities, trials,condtrials, outpath, outtype);
    alldata{sub,1}  = data;

    % Simulate responses using the softmax function.
    % Im still a bit sceptic about whether I should use this one or not. 
    softmaxOutput   = responseModel_v1(params, data, outpath);

end % end of subjects loop

