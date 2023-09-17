% AVERSIVE LEARNING TASK ONLINE VERSION

% Created: @christinadelta -- August 2023

% The script creates a spreadsheet (list of trials) of the aversive learning task to be used
% in the ASD_PAL online study. 

% ---------------------------------------------
%% clear window and workspace

clear all
clc

%% define initial variables 

% initialise variables 
subjects        = 1;
condition       = 6;                        % stable & volatile / small, medium & large stochasticity
task            = 2;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.80 .20;                 % small stochasticity probabilities
    .30 .70;                                % medium stochasticity
    .60 .40];                               % large stochasticity probabilities (either 60:40 or 64:36)
trials          = 140;                      % total trials
condtrials      = {70,[30,20,10,10]};       % 60 per stochasticity condition in stable env and 16 trials in volatile condition;
outtype         = 1;                        % if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes
nCues           = 2;

%% create trials list

% make trial list
% data            = avlearn_trial_list(condition, probabilities, trials, condtrials, outtype, task);
data            = avlearn_trial_list_v2(condition, probabilities, trials, condtrials, outtype, task);



%% make tables


