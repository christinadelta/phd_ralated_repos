% ASSOCIATIVE LEARNING TASK ONLINE VERSION

% Created: @christinadelta -- August 2023

% The script creates a spreadsheet (list of trials) of the associative learning task to be used
% in the ASD_PAL online study. 

%-----------------------
%%

clear all
clc

%% define initial variables 

% initialise variables 
subjects        = 1;
condition       = 6;                        % stable & volatile / small, medium & large stochasticity
task            = 1;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.80 .20;                 % small stochasticity probabilities
    .70 .30;                                % medium stoch probabilities
    .60 .40];                               % large stochasticity probabilities (either 60:40 or 64:36)
trials          = 140;                      % total trials
condtrials      = {70,[30,20,10,10]};       % 70 per stochasticity condition in stable env and 10 trials in volatile condition;
outtype         = 1;                        % if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes
nCues           = 2;
nOut            = 2;

%% create trials list

% data = perclearn_trial_list(condition, probabilities, trials,condtrials, outtype, task);
data = perclearn_trial_list_v2(condition, probabilities, trials,condtrials, outtype, task);



