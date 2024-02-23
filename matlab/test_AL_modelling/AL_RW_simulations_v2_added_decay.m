% Action Learning RW simulations version 2

% added a decay/forgeting parameter to model the natural forgetting process
% or decrease in sakience of unchosen options over time 

%% housekeeping commands

clear all
clc

%% set figure-docking as default 

set(0,'DefaultFigureWindowStyle','docked')

%% define paths etc 
% paths for subjects data 
startpath       = pwd;
datadir         = fullfile(startpath,'data');
subs            = dir(fullfile(datadir, '*sub-*'));
nsubs           = length(subs);

figpath         = fullfile(pwd, 'figures');     addpath(figpath);

%% simulate some data to extract underlying loss rate stc and vol indecies

% initialise variables 
subjects        = 1;
condition       = 6;                        % stable & volatile / small, medium & large stochasticity
task            = 2;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.90 .10;                 % small stochasticity probabilities
    .80 .20;                                % medium stochasticity
    .70 .30];                               % large stochasticity probabilities (either 70:30 or 60:40)
total_trials    = 140;                      % total trials
condtrials      = {70,[30,10,10,20]};
nCues           = 2;
params          = [0.3 1.5 0.3];

data            = ALsimdata_v2(probabilities, total_trials,condtrials);

%% simulate model 

modelout = modelRW_v1(params, data);

%% extract model output

% for plotting we we need Qvalues, choice probabilities and choices both conditions, so, concatinate the
% trials of the two conditions
Qvals   = modelout.Qvals;
Ps      = modelout.allPs; %
choices = modelout.a;
a       = 2 - choices; % convert to be 1 and 0
x       = data.x;
correct = modelout.correct;

% get score
score   = mean(correct);

%% plot RL model

h               = plotRW_v1(x,Ps,Qvals,a);

%%




