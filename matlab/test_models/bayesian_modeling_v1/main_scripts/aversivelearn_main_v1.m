%  AVERSIVE LEARNING MODELING USING THE VOLATILE KALMAN FILTER -- VERSION 1. 

% Date created: 16/02/2023
% modified 1: 

% The script simulates data the aversive learning task, runs vkf model and
% plots outpus. 
% All functions are called using this script.
% Outputs are also stored in tables for easy access and reading
% figures are saved in .svg and .fig formats
% 

% Information about the task: Aversive learning task Version 1:
% A simple 2-armed bandit choice. More details in the simulation function:
% aversivelearn_sim_v1.m

clear all  
clc

%% set paths & initialise variables

% set paths
addpath(fullfile(pwd,'functions'))
addpath(fullfile(pwd,'main_scripts'))
outpath = fullfile(pwd, 'output'); addpath(outpath);
figpath = fullfile(pwd, 'figures'); addpath(figpath);

% initialise variables 
subjects        = 1;
condition       = 1;                      % only stable condition for now (if 2 = stable and volatile)
task            = 1;                      % stable without switch (if task = 2 then stable with one switch)

if condition == 1
    probs           = [.75 .25];          % probabilities of the stable condition
else
    probs           = [.75 .25; .80 .20]; % probabilities of the stable + volatile condition

end

model_params    = {}; % init struct to store all model parameters 

% init parameters fo the VKF model 
model_params.lambda             = 0.1;
model_params.init_vol           = 0.1;
model_params.variance           = 0.1;
model_params.omega              = 0.1;

labels                          = {'lambda', 'initial volatility'...
    'variance', 'omega'};                                   % for ploting
condstring                      = {'stable', 'volatile'};   % for ploting 
trials                          = 100;                      % per volatility condition

%% simulate dataset 

for sub = 1:subjects
    % simulate dataset
    data                = aversivelearn_sim_v2(condition, probs, trials, outpath, task); % data is a structure containaining the simulation output
    allsub_data{1,sub}  = data; % if many subjects add them to cell

end
    
%% model the dataset(s)

for sub = 1:subjects

    this_data               = allsub_data{1,sub};
    


end % end of subjects loop






