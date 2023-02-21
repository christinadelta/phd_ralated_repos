% Bayesian Learning - Version 2 -- Bayes with & without volatility 

% Created: 12/2/2023 -- version 1
% Modified 21/2/2023 -- version 2

% 2nd version of a simple bayesian learning model:
% The script goes through Bayesian modeling without with and volatility 

% Now let's start building a simple Bayesian model for the aversive
% learning task.

% This version includes:
% 1. simulating datasets for volatile/stable conditions
% 2. Running simple Bayesian model where the main goal is to estimate the
% the q values (predictions of the true underlying probabilities)
% 3. adding volatility component
% 4. plotting the results 

clear all  
clc

%% initialise variables

% set paths
addpath(fullfile(pwd,'functions'))
outpath = fullfile(pwd, 'output'); addpath(outpath);
figpath = fullfile(pwd, 'figures'); addpath(figpath);

% initialise variables 
subjects            = 1;
params              = [.25 4];                % alpha and beta values 
condition           = 2;                      % only stable condition for now (if 2 = stable and volatile)
task                = 1;                      % stable without switch (if task = 2 then stable with one switch)

if condition == 1
    probs           = [.75 .25];              % probabilities of the stable condition
else
    probs           = [.75 .25; .80 .20];     % probabilities of the stable + volatile condition
end

labels              = {'alpha', 'beta'};      % for ploting
condstring          = {'stable', 'volatile'}; % for ploting 
trials              = 100;                    % per volatility condition

H                   = 1/25;                   % used to calculate the prior (changing the value affects model predictions)

%% simulate dataset(s)

% for now simulate one dataset per condition 
for cond = 1:condition
    data{1,cond}                = avlearn_simulate_v1(cond, probs, trials, outpath, task); % data is a structure containaining the simulation output
end

% for plotting we we need outcomes and underlying probabilities both conditions, so, concatinate the
% trials of the two conditions
feedback           = cat(1,data{1,1}.feedback, data{1,2}.feedback);
feedbackprob       = cat(1,data{1,1}.feedbackprob, data{1,2}.feedbackprob); % for ploting

% what data do we need for running the model?
q                   = feedbackprob;    % true underlying probabilities
y                   = feedback(:,1);   % feedback -- 1 = vertical gabor, 0 = horizontal gabor, is the highly rewarded

%% 










