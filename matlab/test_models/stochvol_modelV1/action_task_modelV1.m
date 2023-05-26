% first version of the aversive learning task with the
% stochasticity-volatility model 

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
    .62 .38];                               % large stochasticity probabilities
trials          = 200;                      % total trials
condtrials      = [100 25];                 % 100 per stochasticity condition in stable env and 25 trials in volatile condition;

%% simulate dataset

for sub = 1:subjects

    % simulate dataset
    data = action_simdataV1(condition, probabilities, trials,condtrials, outpath);
    


end % end of subjects loop

