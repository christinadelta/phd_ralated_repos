% the script tests different stochasticity levels 

% some code should go here 

%% clear stuff

clear all
clc

%% define paths 
% set paths
addpath(fullfile(pwd,'functions'))

outpath         = fullfile(pwd, 'output');      addpath(outpath);
figpath         = fullfile(pwd, 'figures');     addpath(figpath);
modelpath       = fullfile(pwd, 'stcvol');      addpath(modelpath);
plotpath        = fullfile(pwd, 'plotting');    addpath(plotpath);

%% define initial variables 

% initialise variables 
subjects        = 1;
condition       = 6;                        % stable & volatile / small, medium & large stochasticity
task            = 2;                        % stable without switch (if task = 2 then stable with one switch)
trials          = 140;                      % total trials
condtrials      = {70,[30,10,10,20]};
outtype         = 2;                        % if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes
nCues           = 2;

% create probabiltiies cell
all_probabilities = {[.95 .05; .90 .10; .80 .20], [.90 .10; .80 .20;.70 .30],...
    [.85 .15; .75 .25; .65 .35], [.80 .20; .70 .30; .60 .40]};

%% simulate data and models 

for i = 1:size(all_probabilities,2)

    probabilities   = all_probabilities{1,i};
    
    % simulate some data
    data            = ALsimdata_v2(probabilities, trials, condtrials, outtype);

    % simulate model 
    mout            = simodel_step_one(data,probabilities, trials,condtrials, outtype);
    
    % extract mean and all learning rates
    allLRs_mean{1,i}    = mout.mean_data.mean_alphas{1,1};
    all_lrs{1,i}        = mout.sim_data.lrs;

end % end of for loop 

%% plot learning-rates for different stochasticity levels 



